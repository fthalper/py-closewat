#!/usr/bin/env python3
import sys
import math
import os
from typing import List, Dict, Tuple, Optional, Callable, Any, Union

# Constants from the original C code
MAXHBONDSQ = 3.2 * 3.2  # default max H-bond dist squared
MINHBONDSQ = 2.5 * 2.5  # default min H-bond dist squared
CLOSEHSQ = 1.5 * 1.5    # closest approach to hydrogen
CLOSEHET = 3.9 * 3.9    # default min dist to heteroatom
NUMGAP = 10             # gap in residue number
HBARELY = 0.0891
OBARELY = 0.1984

# Diagnostic strings
DIAGSTR = ["faraway", " ok now", "  metal", "  leave", "altconf", "replace", "   bump", "   edit"]

class PDBRecord:
    """Equivalent to the PDBRECORD struct in C"""
    def __init__(self):
        self.p_rtype = ""        # Record type (e.g., "ATOM  ", "HETATM")
        self.p_attype = ""       # Atom type (e.g., " CA ", " O  ")
        self.p_resname = ""      # Residue name (e.g., "ALA", "HOH")
        self.p_atomid = ""       # Atom ID
        self.p_conf = ' '        # Conformer ID (e.g., 'A', 'B', etc.)
        self.p_chainid = ' '     # Chain ID
        self.p_nconfs = 1        # Number of conformers
        self.p_chaino = ' '      # Original chain ID
        self.p_confo = ' '       # Original conformer ID
        self.p_nconfo = 1        # Original number of conformers
        self.p_diag = 0          # Diagnostic code
        self.p_rff2 = 0          # Reserved for future use
        self.p_xc = 0.0          # X coordinate
        self.p_yc = 0.0          # Y coordinate
        self.p_zc = 0.0          # Z coordinate
        self.p_occ = 0.0         # Occupancy
        self.p_bval = 0.0        # B value (temperature factor)
        self.p_occo = 0.0        # Original occupancy
        self.p_bvo = 0.0         # Original B value
        self.p_dfc = 0.0         # Distance from closest
        self.p_atnum = 0         # Atom number
        self.p_rff3 = 0          # Reserved for future use
        self.p_resnum = 0        # Residue number
        self.p_resno = 0         # Original residue number
        self.p_nbr = None        # Pointer to neighbor (will be index in Python)

class Chain:
    """Equivalent to the CHAIN struct in C"""
    def __init__(self):
        self.c_chainid = ' '     # Chain ID
        self.c_crffu = 0         # Reserved for future use
        self.c_minwat = 0        # Minimum water number
        self.c_curwat = 0        # Current water number
        self.c_wat0 = 0          # First water index
        self.c_watl1 = 0         # Last water index + 1
        self.c_watm0 = 0         # First multiple-conformer water index
        self.c_watl = 0          # Last water index

class TotalSt:
    """Equivalent to the TOTALST struct in C"""
    def __init__(self):
        self.tnhet = 0           # Number of heteroatoms
        self.tnwatat = 0         # Number of water atoms
        self.tnwathet = 0        # Number of water heteroatoms
        self.tnch = 0            # Number of chains
        self.tnatoms = 0         # Total number of atoms
        self.tnpat = 0           # Number of protein atoms
        self.tnpwat = 0          # Number of protein waters
        self.tnwaters = 0        # Number of waters
        self.t1ch_s = 0          # Flag for single chain
        self.thighbq = 0         # Flag for high B quality
        self.tnclose = 0         # Number of close contacts
        self.tbump = 0           # Flag for bump
        self.talert = [0] * 8    # Alert counters
        self.tmeanbw = 0.0       # Mean B value for waters
        self.tmaxhbsq = MAXHBONDSQ  # Max H-bond distance squared
        self.tminhbsq = MINHBONDSQ  # Min H-bond distance squared
        self.tminhsq = CLOSEHSQ     # Min hydrogen distance squared
        self.tminhet = CLOSEHET     # Min heteroatom distance
        self.tpat = []           # List of protein atoms
        self.tpatp = 0           # Index into tpat
        self.tpwa = []           # List of water atoms
        self.tpwap = 0           # Index into tpwa
        self.tpate = 0           # End of protein atoms
        self.tchs = []           # List of chains
        self.tfpi = None         # Input file
        self.tfpwat = None       # Output file
        self.tfpl = None         # Log file

def outrec(prec: PDBRecord, fp) -> None:
    """Writes out a PDBRecord structure as a string"""
    # Format the atom ID for output
    atomid = prec.p_atomid
    if len(atomid) == 1:
        atomid = f" {atomid}"
    
    # Format the PDB record line
    line = f"{prec.p_rtype:<6}{prec.p_atnum:5d} {prec.p_attype:<4}{prec.p_conf}{prec.p_resname:<3} {prec.p_chainid}{prec.p_resnum:4d}    {prec.p_xc:8.3f}{prec.p_yc:8.3f}{prec.p_zc:8.3f}{prec.p_occ:6.2f}{prec.p_bval:6.2f}          {atomid}  "
    
    # Write to file
    fp.write(line + "\n")

def pdbdist(p0: PDBRecord, p1: PDBRecord) -> float:
    """Calculates the distance squared between two atoms"""
    dr = (p0.p_xc - p1.p_xc) ** 2 + (p0.p_yc - p1.p_yc) ** 2 + (p0.p_zc - p1.p_zc) ** 2
    return dr

def isconformer(top: TotalSt, p0: PDBRecord, p1: PDBRecord, checkdist: int) -> int:
    """Returns true only if p1 is a later conformer than p0"""
    if p0 is None or p1 is None:
        return 0
    if p0.p_chainid != p1.p_chainid:
        return 0
    if p0.p_resnum != p1.p_resnum:
        return 0
    if p0.p_conf not in ['A', 'B', 'C']:
        return 0
    if p1.p_conf not in ['B', 'C', 'D']:
        return 0
    
    dif = ord(p1.p_conf) - ord(p0.p_conf)
    if dif < 1 or dif > 3:
        return 0
    if abs(p0.p_bval - p1.p_bval) > 0.005:
        return 0
    
    if checkdist:
        return 0 if pdbdist(p0, p1) > top.tminhbsq else 1
    else:
        return 1

def occsort(pr0: PDBRecord, pr1: PDBRecord) -> int:
    """Compare routine for sorting on the basis of occupancy"""
    return 1 if pr0.p_occ < pr1.p_occ else -1

def occbsort(pr0: PDBRecord, pr1: PDBRecord) -> int:
    """Compare routine for sorting, on conformer, then on occupancy, and finally on B value"""
    if pr1.p_conf != ' ':
        if pr0.p_conf != ' ':
            if pr0.p_resnum == pr1.p_resnum and pr0.p_chainid == pr1.p_chainid:
                dif = ord(pr1.p_conf) - ord(pr0.p_conf)
                return -1 if dif > 0 else (0 if dif == 0 else 1)
            else:
                return -1 if pr0.p_bval < pr1.p_bval else 1
        else:
            return -1
    
    if pr0.p_conf != ' ':
        return 1
    
    if pr0.p_conf == ' ' and pr1.p_conf != ' ':
        return -1
    if pr1.p_conf == ' ' and pr0.p_conf != ' ':
        return 1
    
    if pr0.p_occ > pr1.p_occ:
        return -1
    if pr0.p_occ < pr1.p_occ:
        return 1
    
    return -1 if pr0.p_bval < pr1.p_bval else 1

def chrsort(pr0: PDBRecord, pr1: PDBRecord) -> int:
    """Compare routine for sorting, on chain ID and secondarily on residue #"""
    if pr0.p_chainid < pr1.p_chainid:
        return -1
    elif pr0.p_chainid > pr1.p_chainid:
        return 1
    
    if pr0.p_resnum < pr1.p_resnum:
        return -1
    elif pr0.p_resnum > pr1.p_resnum:
        return 1
    
    return -1 if pr0.p_conf < pr1.p_conf else 1

def bonly(pr0: PDBRecord, pr1: PDBRecord) -> int:
    """Sort routine to operate only based on conformer letter and B value"""
    if pr0.p_resnum == pr1.p_resnum and pr0.p_conf != pr1.p_conf:
        dif = ord(pr0.p_conf) - ord(pr1.p_conf)
        return -1 if dif < 0 else 1
    
    return -1 if pr0.p_bval <= pr1.p_bval else 1

def strtorec(line: str, ptp: PDBRecord) -> None:
    """Converts a PDB text string into a PDB record"""
    legal_conf = " ABCD"
    
    # Initialize the record
    ptp.p_nbr = None
    
    # Extract record type (columns 1-6)
    ptp.p_rtype = line[0:6].strip()
    
    # Extract atom number (columns 7-11)
    ptp.p_atnum = int(line[6:11].strip())
    
    # Extract atom type (columns 13-16)
    ptp.p_attype = line[12:16].strip()
    
    # Extract conformer ID (column 17)
    ptp.p_confo = ptp.p_conf = line[16] if 16 < len(line) else ' '
    ptp.p_nconfo = ptp.p_nconfs = 1  # this can be modified later
    
    # Only legal conformers are A,B,C,D
    if ptp.p_conf not in legal_conf:
        ptp.p_conf = ' '
    
    if ptp.p_conf != ' ':
        ptp.p_nconfo = 2  # first stab
    
    # Extract residue name (columns 18-20)
    ptp.p_resname = line[17:20].strip() if 20 <= len(line) else ""
    
    # Extract chain ID (column 22)
    ptp.p_chaino = ptp.p_chainid = line[21] if 21 < len(line) else ' '
    
    # Extract residue number (columns 23-26)
    ptp.p_resno = ptp.p_resnum = int(line[22:26].strip()) if 26 <= len(line) else 0
    
    # Extract coordinates (columns 27-54)
    ptp.p_xc = float(line[26:38].strip()) if 38 <= len(line) else 0.0
    ptp.p_yc = float(line[38:46].strip()) if 46 <= len(line) else 0.0
    ptp.p_zc = float(line[46:54].strip()) if 54 <= len(line) else 0.0
    
    # Extract occupancy and B value (columns 55-66)
    ptp.p_occo = ptp.p_occ = float(line[54:60].strip()) if 60 <= len(line) else 0.0
    ptp.p_bvo = ptp.p_bval = float(line[60:66].strip()) if 66 <= len(line) else 0.0
    
    # Extract atom ID (columns 77-78)
    if len(line) >= 78:
        ptp.p_atomid = line[76:78].strip()
    else:
        ptp.p_atomid = ""

def findchain(line: str, top: TotalSt) -> None:
    """Find chain information in the PDB line"""
    if line.startswith("ATOM  ") or line.startswith("HETATM"):
        chain_id = line[21] if len(line) > 21 else ' '
        
        # Check if this chain already exists
        for i, chain in enumerate(top.tchs):
            if chain.c_chainid == chain_id:
                return
        
        # If not, add a new chain
        if len(top.tchs) < 50:  # Original C code allocated 50 chains
            new_chain = Chain()
            new_chain.c_chainid = chain_id
            top.tchs.append(new_chain)
            top.tnch = len(top.tchs)

def getpdblin(line: str, top: TotalSt) -> None:
    """Process a line from the PDB file during the first pass"""
    if line.startswith("ATOM  ") or line.startswith("HETATM"):
        findchain(line, top)
        
        if line.startswith("ATOM  "):
            top.tnatoms += 1
            top.tnpat += 1
        elif line.startswith("HETATM"):
            top.tnatoms += 1
            
            # Check if it's a water molecule
            resname = line[17:20].strip() if len(line) > 20 else ""
            if resname == "HOH" or resname == "WAT":
                top.tnwatat += 1
            else:
                top.tnhet += 1

def procargs(ac: int, av: List[str], top: TotalSt) -> None:
    """Process command line arguments"""
    # Initialize top structure
    top.tmeanbw = top.tnhet = top.tnwatat = 0
    top.tnatoms = top.tnpat = top.tnpwat = top.tnwaters = 0
    top.tfpl = None
    top.tnch = 5000
    top.tchs = [Chain() for _ in range(50)]  # Allocate 50 chains
    top.t1ch_s = top.thighbq = 0
    top.tmaxhbsq = MAXHBONDSQ
    top.tminhbsq = MINHBONDSQ
    top.tminhsq = CLOSEHSQ
    top.tminhet = CLOSEHET
    top.tbump = 0
    
    thismaxhb = thisminhb = thishyd = thishet = 0.0
    
    progname = av[0]
    av = av[1:]
    ac -= 1
    
    # Process command line options
    while ac > 0 and av[0].startswith('-'):
        for cp in av[0][1:]:
            cp = cp.upper()
            if cp == 'S':
                top.t1ch_s = 1
            elif cp == 'H':
                top.thighbq = 1
            elif cp == 'X':
                thismaxhb = float(av[0][av[0].index('X')+1:])
                break
            elif cp == 'M':
                thisminhb = float(av[0][av[0].index('M')+1:])
                break
            elif cp == 'L':
                thishyd = float(av[0][av[0].index('L')+1:])
                break
            elif cp == 'O':
                thishet = float(av[0][av[0].index('O')+1:])
                break
            elif cp == 'B':
                top.tbump = 1
        
        av = av[1:]
        ac -= 1
    
    # Set parameters based on command line options
    if 0.05 < thismaxhb < 10.0:
        top.tmaxhbsq = thismaxhb * thismaxhb
    if 0.05 < thisminhb < 10.0:
        top.tminhbsq = thisminhb * thisminhb
    if 0.05 < thishyd < 10.0:
        top.tminhsq = thishyd * thishyd
    if 0.05 < thishet < 10.0:
        top.tminhet = thishet * thishet
    
    # Set up input and output files
    top.tfpi = sys.stdin
    top.tfpwat = sys.stdout
    top.tfpl = open("closewat.log", "w")
    
    if ac > 0:
        try:
            top.tfpi = open(av[0], "r")
        except:
            top.tfpi = sys.stdin
        
        if ac > 1:
            try:
                top.tfpwat = open(av[1], "w")
            except:
                top.tfpwat = sys.stdout
    
    top.tfpl.write(f"Run of {progname}\n")

def ready(top: TotalSt) -> int:
    """Prepare for processing by allocating memory and rewinding input file"""
    # Allocate memory for protein atoms
    top.tpat = [PDBRecord() for _ in range(top.tnatoms + 10)]
    top.tpatp = 0
    
    # Allocate memory for water atoms
    top.tpwa = [PDBRecord() for _ in range(top.tnwatat + 10)]
    top.tpwap = 0
    
    # Rewind input file if it's not stdin
    if top.tfpi != sys.stdin:
        top.tfpi.seek(0)
    
    return 0

def procpdblin(line: str, top: TotalSt) -> None:
    """Process a line from the PDB file during the second pass"""
    if not (line.startswith("ATOM  ") or line.startswith("HETATM")):
        return
    
    # Check if it's a water molecule
    resname = line[17:20].strip() if len(line) > 20 else ""
    is_water = resname == "HOH" or resname == "WAT"
    
    if line.startswith("ATOM  "):
        # Process protein atom
        if top.tpatp < len(top.tpat):
            strtorec(line, top.tpat[top.tpatp])
            top.tpatp += 1
    elif line.startswith("HETATM"):
        if is_water:
            # Process water atom
            if top.tpwap < len(top.tpwa):
                strtorec(line, top.tpwa[top.tpwap])
                top.tpwap += 1
        else:
            # Process other heteroatom
            if top.tpatp < len(top.tpat):
                strtorec(line, top.tpat[top.tpatp])
                top.tpatp += 1

def ismetal(atstr: str) -> int:
    """Check if atom is a metal"""
    metals = ["LI", "NA", "K", "RB", "CS", "FR", "BE", "MG", "CA", "SR", "BA", "RA",
              "AL", "GA", "IN", "SN", "TL", "PB", "BI", "SC", "TI", "V", "CR", "MN",
              "FE", "CO", "NI", "CU", "ZN", "Y", "ZR", "NB", "MO", "TC", "RU", "RH",
              "PD", "AG", "CD", "LA", "HF", "TA", "W", "RE", "OS", "IR", "PT", "AU",
              "HG", "CE", "PR", "ND", "PM", "SM", "EU", "GD", "TB", "DY", "HO", "ER",
              "TM", "YB", "LU", "TH", "PA", "U", "NP", "PU"]
    
    atstr = atstr.strip()
    for metal in metals:
        if atstr == metal:
            return 1
    return 0

def noconformers(top: TotalSt) -> int:
    """Get rid of input AHOH, BHOH, etc.
    Since we're going to identify conformers ourselves,
    we remove the input ones. We rescale the occupancy and B values."""
    maxres = -1000
    ncleared = 0
    
    # First determine the highest input residue number
    for i in range(top.tpwap):
        pwap = top.tpwa[i]
        if pwap.p_resnum > maxres:
            maxres = pwap.p_resnum
    
    if maxres < -999:
        return 0
    
    for i in range(top.tpwap):
        pwap = top.tpwa[i]
        if pwap.p_conf != ' ':
            if pwap.p_conf != 'A':
                maxres += 1
                pwap.p_resnum = maxres
            
            corocc = 0.9 if pwap.p_occ < 0.81 else 1.0 / math.sqrt(pwap.p_occ)
            pwap.p_bval *= corocc
            
            if pwap.p_bval > 80.0:
                print(f"noconformers:B={pwap.p_bval:8.2f} to {pwap.p_conf}{pwap.p_resname} {pwap.p_chainid}{pwap.p_resnum:4d} @[{pwap.p_xc:5.1f},{pwap.p_yc:5.1f},{pwap.p_zc:5.1f}]: corocc={corocc:7.4f}")
            
            if pwap.p_occ < 0.5:
                pwap.p_occ *= 2.0
            else:
                pwap.p_occ = 1.0
            
            pwap.p_conf = ' '
            ncleared += 1
    
    if ncleared:
        print(f" {ncleared:6d} waters had their conformation symbols cleared")
        if top.tfpl is not None:
            top.tfpl.write(f" {ncleared:6d} waters had their conformation symbols cleared\n")
    
    return ncleared

def closestwat(top: TotalSt) -> None:
    """For each water, find the nearest-neighbor water.
    If it's closer than 2.5 A, we report it as a possible collision.
    If it's between 2.5 and 3.2 A, we report it as a hydrogen bond.
    We check all waters, not just the ones that occur later in the list than the current one."""
    for i in range(top.tpwap - 1):
        p0 = top.tpwa[i]
        mindist = 1.0e10
        p2 = None
        
        for j in range(top.tpwap):
            p1 = top.tpwa[j]
            if i != j:  # p1 != p0
                dist = pdbdist(p0, p1)
                if dist < mindist:
                    mindist = dist
                    p2 = p1
        
        if mindist < top.tmaxhbsq and p2 is not None:
            sum_water(top, p0, p2, mindist)
            if mindist < top.tminhbsq:
                relateem(top, p0, p2)

def thisthird(top: TotalSt, p0: PDBRecord) -> None:
    """Analyze neighbors of a water molecule"""
    if p0 is None or p0.p_nbr is None:
        return
    
    # In the C code, p_nbr is a pointer to another PDBRecord
    # In our Python implementation, it's a direct reference to another PDBRecord object
    p1 = p0.p_nbr
    
    # Check if p1 is already a conformer
    if p1.p_conf != ' ':
        return
    
    # Calculate distance between p0 and p1
    dist = pdbdist(p0, p1)
    
    # If they're too far apart, return
    if dist > top.tmaxhbsq:
        return
    
    # If they're close enough, relate them
    if dist < top.tminhbsq:
        relateem(top, p0, p1)

def meanbw(top: TotalSt) -> None:
    """Calculate mean B value for waters"""
    if top.tpwap <= 0:
        top.tmeanbw = 0.0
        return
    
    sum_b = 0.0
    count = 0
    
    for i in range(top.tpwap):
        pwap = top.tpwa[i]
        sum_b += pwap.p_bval
        count += 1
    
    if count > 0:
        top.tmeanbw = sum_b / count
        if top.tfpl is not None:
            top.tfpl.write(f"Mean B value for waters = {top.tmeanbw:.2f}\n")

def reduceocc(top: TotalSt) -> None:
    """Reduce occupancy of waters with B > 1.2*<B>"""
    if top.tmeanbw <= 0.0:
        return
    
    bthresh = 1.2 * top.tmeanbw
    nreduced = 0
    
    for i in range(top.tpwap):
        pwap = top.tpwa[i]
        if pwap.p_bval > bthresh and pwap.p_occ > 0.5:
            # Calculate reduction factor based on B value
            factor = top.tmeanbw / pwap.p_bval
            if factor < 0.5:
                factor = 0.5
            
            # Reduce occupancy
            pwap.p_occ *= factor
            nreduced += 1
    
    if nreduced > 0 and top.tfpl is not None:
        top.tfpl.write(f"Reduced occupancy of {nreduced} waters with B > {bthresh:.2f}\n")

def makechains(top: TotalSt) -> None:
    """Figure out which chain each water belongs to"""
    # Initialize chain information
    for i in range(len(top.tchs)):
        chain = top.tchs[i]
        chain.c_wat0 = 0
        chain.c_watl1 = 0
        chain.c_watm0 = 0
        chain.c_watl = 0
        chain.c_minwat = 1  # Default starting water number
    
    # Assign each water to the nearest protein chain
    for i in range(top.tpwap):
        pwap = top.tpwa[i]
        nearest_chain = None
        min_dist = float('inf')
        
        # Find the nearest protein chain
        for j in range(top.tpatp):
            pat = top.tpat[j]
            if pat.p_chainid != ' ':  # Skip atoms without a chain
                dist = pdbdist(pwap, pat)
                if dist < min_dist:
                    min_dist = dist
                    nearest_chain = pat.p_chainid
        
        # Assign the water to the nearest chain
        if nearest_chain is not None:
            pwap.p_chainid = nearest_chain
            
            # Update the chain's water count
            for j in range(len(top.tchs)):
                if top.tchs[j].c_chainid == nearest_chain:
                    top.tchs[j].c_watl += 1
                    break

def adjustmult(top: TotalSt) -> None:
    """Count and modify multiple conformers"""
    # Count the number of waters in each chain
    chain_counts = {}
    for i in range(top.tpwap):
        pwap = top.tpwa[i]
        if pwap.p_chainid not in chain_counts:
            chain_counts[pwap.p_chainid] = 0
        chain_counts[pwap.p_chainid] += 1
    
    # Initialize residue numbers for each chain
    for i in range(len(top.tchs)):
        chain = top.tchs[i]
        chain.c_curwat = chain.c_minwat
    
    # Assign residue numbers to waters based on their chain
    for i in range(top.tpwap):
        pwap = top.tpwa[i]
        for j in range(len(top.tchs)):
            chain = top.tchs[j]
            if chain.c_chainid == pwap.p_chainid:
                # If this is a new conformer (A or B), assign it the same residue number
                if pwap.p_conf in ['A', 'B'] and i > 0 and top.tpwa[i-1].p_conf in ['A', 'B'] and \
                   top.tpwa[i-1].p_chainid == pwap.p_chainid and top.tpwa[i-1].p_resnum == pwap.p_resnum:
                    continue  # Keep the same residue number
                else:
                    # Assign a new residue number
                    pwap.p_resnum = chain.c_curwat
                    chain.c_curwat += 1
                break

def sortmults(top: TotalSt) -> None:
    """Reorder multiple-conformer waters within each chain"""
    # Group waters by chain and residue number
    chain_res_groups = {}
    
    for i in range(top.tpwap):
        pwap = top.tpwa[i]
        key = (pwap.p_chainid, pwap.p_resnum)
        if key not in chain_res_groups:
            chain_res_groups[key] = []
        chain_res_groups[key].append(pwap)
    
    # Sort each group by conformer ID
    for key, group in chain_res_groups.items():
        if len(group) > 1:
            # Sort by conformer ID (A before B before C, etc.)
            group.sort(key=lambda x: x.p_conf)
            
            # Ensure that if there's an 'A' conformer, it has higher occupancy than 'B'
            if len(group) >= 2 and group[0].p_conf == 'A' and group[1].p_conf == 'B':
                if group[0].p_occ < group[1].p_occ:
                    # Swap occupancies
                    group[0].p_occ, group[1].p_occ = group[1].p_occ, group[0].p_occ

def proximity(top: TotalSt) -> None:
    """For all waters, identify close contacts"""
    # Check each water against all protein atoms
    for i in range(top.tpwap):
        pwap = top.tpwa[i]
        
        # Skip waters that are already marked for special handling
        if pwap.p_diag != 0:
            continue
        
        # Check against all protein atoms
        for j in range(top.tpatp):
            pat = top.tpat[j]
            
            # Calculate distance
            dist = pdbdist(pwap, pat)
            
            # Check for close contacts
            if dist < top.tminhsq:
                # Too close to a hydrogen - mark for attention
                pwap.p_diag = 6  # Bump
                top.talert[6] += 1
                
                if top.tfpl is not None:
                    top.tfpl.write(f"Water {pwap.p_chainid}{pwap.p_resnum} too close to {pat.p_chainid}{pat.p_resnum} {pat.p_attype}: {math.sqrt(dist):.2f} Å\n")
                
                break
            elif dist < top.tminhet and ismetal(pat.p_attype):
                # Close to a metal - mark as metal-coordinated
                pwap.p_diag = 2  # Metal
                top.talert[2] += 1
                
                if top.tfpl is not None:
                    top.tfpl.write(f"Water {pwap.p_chainid}{pwap.p_resnum} coordinated to metal {pat.p_chainid}{pat.p_resnum} {pat.p_attype}: {math.sqrt(dist):.2f} Å\n")
                
                break

def insert_singles(top: TotalSt) -> None:
    """Move single-conformation but conformation-marked waters into the empty zone
    between the unconformation-marked waters and the multiple-conformation waters"""
    # Count waters by type
    unmarked_count = 0
    single_conf_marked = []
    multi_conf = []
    
    for i in range(top.tpwap):
        pwap = top.tpwa[i]
        if pwap.p_conf == ' ':
            unmarked_count += 1
        elif pwap.p_nconfs == 1:
            single_conf_marked.append(pwap)
        else:
            multi_conf.append(pwap)
    
    # If there are no single-conformation marked waters, nothing to do
    if not single_conf_marked:
        return
    
    # Create a new ordered list of waters
    new_waters = []
    
    # Add unmarked waters first
    for i in range(top.tpwap):
        if top.tpwa[i].p_conf == ' ':
            new_waters.append(top.tpwa[i])
    
    # Add single-conformation marked waters
    new_waters.extend(single_conf_marked)
    
    # Add multiple-conformation waters
    new_waters.extend(multi_conf)
    
    # Replace the original water list with the reordered one
    for i in range(min(len(new_waters), top.tpwap)):
        top.tpwa[i] = new_waters[i]

def finalmult(top: TotalSt) -> None:
    """Final count of multiple conformers"""
    # Count waters by conformer type
    unmarked_count = 0
    a_conf_count = 0
    b_conf_count = 0
    other_conf_count = 0
    
    for i in range(top.tpwap):
        pwap = top.tpwa[i]
        if pwap.p_conf == ' ':
            unmarked_count += 1
        elif pwap.p_conf == 'A':
            a_conf_count += 1
        elif pwap.p_conf == 'B':
            b_conf_count += 1
        else:
            other_conf_count += 1
    
    # Calculate total number of waters
    total_waters = unmarked_count + a_conf_count + b_conf_count + other_conf_count
    
    # Calculate number of water positions (A and B conformers of the same residue count as one position)
    water_positions = unmarked_count + a_conf_count  # Each A conformer represents one position
    
    if top.tfpl is not None:
        top.tfpl.write(f"\nFinal water statistics:\n")
        top.tfpl.write(f"  Total water molecules: {total_waters}\n")
        top.tfpl.write(f"  Unmarked waters: {unmarked_count}\n")
        top.tfpl.write(f"  'A' conformers: {a_conf_count}\n")
        top.tfpl.write(f"  'B' conformers: {b_conf_count}\n")
        top.tfpl.write(f"  Other conformers: {other_conf_count}\n")
        top.tfpl.write(f"  Total water positions: {water_positions}\n")
        
        # Report on diagnostic codes
        for i in range(8):
            if top.talert[i] > 0:
                top.tfpl.write(f"  {DIAGSTR[i]}: {top.talert[i]}\n")

def sum_water(top: TotalSt, pwap: PDBRecord, thisp: PDBRecord, mindist: float) -> None:
    """Sum water information"""
    warmsg = "clash"
    hbmsg = "Hbond"
    
    if top.tfpl is None:
        return
    
    top.tfpl.write(
        f"{pwap.p_chainid}{pwap.p_resnum:4d} @{pwap.p_xc:7.3f} {pwap.p_yc:7.3f} {pwap.p_zc:7.3f} {warmsg if mindist < top.tminhbsq else hbmsg} {math.sqrt(mindist):5.2f} A to"
    )
    top.tfpl.write(
        f" {thisp.p_chainid}{thisp.p_resnum:4d} @{thisp.p_xc:7.3f} {thisp.p_yc:7.3f} {thisp.p_zc:7.3f}\n"
    )

def relateem(top: TotalSt, pwap: PDBRecord, thisp: PDBRecord) -> None:
    """Relate two water molecules.
    This defines AHOH and BHOH records to replace the input records.
    If they're already AHOH and BHOH relative to one another, we skip them."""
    if pwap is None or thisp is None:
        return
    
    # For now, we'll only merge waters that are currently unmerged; but we'll remember the others
    if pwap.p_conf != ' ' or thisp.p_conf != ' ':
        if (pwap.p_chainid == thisp.p_chainid and 
            pwap.p_resnum == thisp.p_resnum and
            ((pwap.p_conf == 'A' and thisp.p_conf == 'B') or
             (pwap.p_conf == 'B' and thisp.p_conf == 'A'))):
            return
        
        if pwap.p_conf != ' ':
            thisp.p_nbr = pwap  # In Python, we can directly reference the object
        if thisp.p_conf != ' ':
            pwap.p_nbr = thisp
        return
    
    occ00 = pwap.p_occ
    occ01 = thisp.p_occ
    b00 = pwap.p_bval
    b01 = thisp.p_bval
    
    # Avoid division by zero
    total_occ = occ00 + occ01
    if total_occ <= 0.0 or b00 <= 0.0 or b01 <= 0.0:
        occ0 = 0.5
        occ1 = 0.5
    else:
        occ0 = (occ00 / total_occ) * (b00 + b01) / (2 * b00)
        occ1 = (occ01 / total_occ) * (b00 + b01) / (2 * b01)
    
    if occ0 > occ1:
        if occ0 < 0.5:
            occ0 = 0.5
        elif occ0 > 0.75:
            occ0 = 0.75
        occ1 = 1.0 - occ0
    else:
        if occ1 < 0.5:
            occ1 = 0.5
        elif occ1 > 0.75:
            occ1 = 0.75
        occ0 = 1.0 - occ1
    
    wt = abs(b00 - b01) / (b00 + b01)
    if wt < 0.07:
        bscal = 0.8
    else:
        bscal = 0.9 if wt < 0.3 else 0.96
    
    b0 = bscal * (b00 if b00 < b01 else b01)
    
    pwap.p_nconfs = thisp.p_nconfs = 2
    pwap.p_occ = occ0
    thisp.p_occ = occ1
    
    if b0 > 90.0:
        print(f"Assigning B={b0:8.2f} to {pwap.p_conf}{pwap.p_resname} {pwap.p_chainid}{pwap.p_resnum:4d} @[{pwap.p_xc:5.1f},{pwap.p_yc:5.1f},{pwap.p_zc:5.1f}] in relateem")
    
    pwap.p_bval = thisp.p_bval = b0
    
    if occ0 >= occ1:
        # pwap is more than 50% occupied: make it primary
        thisp.p_resnum = pwap.p_resnum
        thisp.p_chainid = pwap.p_chainid
        pwap.p_conf = 'A'
        thisp.p_conf = 'B'
    else:
        # thisp is more than 50% occupied: make it primary
        pwap.p_resnum = thisp.p_resnum
        pwap.p_chainid = thisp.p_chainid
        pwap.p_conf = 'B'
        thisp.p_conf = 'A'

def main() -> int:
    """Main function"""
    # Create and initialize the total structure
    tos = TotalSt()
    
    # Process command line arguments
    procargs(len(sys.argv), sys.argv, tos)
    
    # First pass: read the PDB file to count atoms and chains
    for line in tos.tfpi:
        getpdblin(line, tos)
    
    # Prepare for second pass
    if ready(tos) != 0:
        return 1
    
    # Second pass: process the PDB file
    for line in tos.tfpi:
        procpdblin(line, tos)
    
    # Calculate the number of waters
    tos.tnwaters = tos.tpwap
    
    # Get rid of input AHOH, BHOH, etc.
    noconformers(tos)
    
    # Set the end of protein atoms
    tos.tpate = tos.tpatp
    
    # For each water, find the nearest-neighbor water
    closestwat(tos)
    
    # Analyze neighbors of waters
    for i in range(tos.tnwaters):
        pwap = tos.tpwa[i]
        if pwap.p_nbr is not None and pwap.p_conf == ' ':
            # This water has a close neighbor but isn't labeled yet
            thisthird(tos, pwap)
    
    # Calculate mean B value for waters
    meanbw(tos)
    
    # If user requested it, reduce occupancy of waters with B > 1.2*<B>
    if tos.thighbq:
        reduceocc(tos)
    
    # Sort water records on occupancy and B value
    tos.tpwa.sort(key=lambda x: (x.p_conf != ' ', -x.p_occ, x.p_bval))
    
    # Figure out which chain each water belongs to
    makechains(tos)
    
    # Count and modify multiple conformers
    adjustmult(tos)
    
    # Sort waters on chain # and residue #
    tos.tpwa.sort(key=lambda x: (x.p_chainid, x.p_resnum, x.p_conf))
    
    # Reorder multiple-conformer waters within each chain
    sortmults(tos)
    
    # Initialize alert counters
    tos.talert = [0] * 8
    
    # If there's only one chain and the user has specified the 'S' flag,
    # then just write to a chain called S
    if len(tos.tchs) == 1 and tos.t1ch_s and tos.tchs[0].c_chainid != 'S':
        tos.tchs[0].c_curwat = tos.tchs[0].c_minwat = 1
        i = tos.tpwa[0].p_resnum - 1
        for j in range(tos.tnwaters):
            tos.tpwa[j].p_chainid = 'S'
            tos.tpwa[j].p_resnum -= i
    
    # For all waters, identify close contacts
    proximity(tos)
    
    # Move single-conformation but conformation-marked waters
    insert_singles(tos)
    
    # Final count of multiple conformers
    finalmult(tos)
    
    # Write out the waters
    for i in range(tos.tnwaters):
        outrec(tos.tpwa[i], tos.tfpwat)
    
    # Clean up
    if tos.tfpl != None:
        tos.tfpl.close()
    
    if tos.tfpi != sys.stdin:
        tos.tfpi.close()
    
    if tos.tfpwat != sys.stdout:
        tos.tfpwat.close()
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
