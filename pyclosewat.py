#!/usr/bin/env python3
import sys
import math
import os
import argparse
from typing import List, Dict, Tuple, Optional, Callable, Any, Union
from contextlib import contextmanager

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

def cpypdb(pin: PDBRecord, pou: PDBRecord) -> None:
    """Copies the PDBRecord pin into pou"""
    pou.p_rtype = pin.p_rtype
    pou.p_attype = pin.p_attype
    pou.p_resname = pin.p_resname
    pou.p_atomid = pin.p_atomid
    pou.p_conf = pin.p_conf
    pou.p_chainid = pin.p_chainid
    pou.p_nconfs = pin.p_nconfs
    pou.p_chaino = pin.p_chaino
    pou.p_confo = pin.p_confo
    pou.p_nconfo = pin.p_nconfo
    pou.p_diag = pin.p_diag
    pou.p_rff2 = pin.p_rff2
    pou.p_xc = pin.p_xc
    pou.p_yc = pin.p_yc
    pou.p_zc = pin.p_zc
    pou.p_occ = pin.p_occ
    pou.p_bval = pin.p_bval
    pou.p_occo = pin.p_occo
    pou.p_bvo = pin.p_bvo
    pou.p_dfc = pin.p_dfc
    pou.p_atnum = pin.p_atnum
    pou.p_rff3 = pin.p_rff3
    pou.p_resnum = pin.p_resnum
    pou.p_resno = pin.p_resno
    pou.p_nbr = pin.p_nbr

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

def get_sort_key_occ(pr: PDBRecord) -> Tuple[float, float]:
    """Get sorting key for occupancy-based sorting."""
    return (-pr.p_occ, pr.p_bval)  # Negative for descending order

def get_sort_key_occb(pr: PDBRecord) -> Tuple[bool, float, float]:
    """Get sorting key for conformer, occupancy, and B-value based sorting."""
    has_conf = pr.p_conf != ' '
    return (has_conf, -pr.p_occ, pr.p_bval)

def get_sort_key_chain(pr: PDBRecord) -> Tuple[str, int, str]:
    """Get sorting key for chain-based sorting."""
    return (pr.p_chainid, pr.p_resnum, pr.p_conf)

def get_sort_key_b(pr: PDBRecord) -> Tuple[int, float]:
    """Get sorting key for B-value based sorting."""
    return (pr.p_resnum, pr.p_bval)

def strtorec(line: str, ptp: PDBRecord) -> None:
    """Converts a PDB text string into a PDB record with proper error handling."""
    legal_conf = " ABCD"
    
    # Initialize the record
    ptp.p_nbr = None
    
    # Ensure minimum line length for PDB format
    if len(line) < 66:  # Minimum length for required fields
        raise ValueError(f"PDB line too short: {len(line)} chars")
    
    try:
        # Extract record type (columns 1-6)
        ptp.p_rtype = line[0:6].strip()
        
        # Extract atom number (columns 7-11)
        ptp.p_atnum = int(line[6:11].strip())
        
        # Extract atom type (columns 13-16)
        ptp.p_attype = line[12:16].strip()
        
        # Extract conformer ID (column 17)
        ptp.p_confo = ptp.p_conf = line[16] if len(line) > 16 else ' '
        ptp.p_nconfo = ptp.p_nconfs = 1  # this can be modified later
        
        # Only legal conformers are A,B,C,D
        if ptp.p_conf not in legal_conf:
            ptp.p_conf = ' '
        
        if ptp.p_conf != ' ':
            ptp.p_nconfo = 2  # first stab
        
        # Extract residue name (columns 18-20)
        ptp.p_resname = line[17:20].strip() if len(line) > 20 else ""
        
        # Extract chain ID (column 22)
        ptp.p_chaino = ptp.p_chainid = line[21] if len(line) > 21 else ' '
        
        # Extract residue number (columns 23-26)
        ptp.p_resno = ptp.p_resnum = int(line[22:26].strip()) if len(line) > 26 else 0
        
        # Extract coordinates (columns 27-54)
        ptp.p_xc = float(line[26:38].strip()) if len(line) > 38 else 0.0
        ptp.p_yc = float(line[38:46].strip()) if len(line) > 46 else 0.0
        ptp.p_zc = float(line[46:54].strip()) if len(line) > 54 else 0.0
        
        # Extract occupancy and B value (columns 55-66)
        ptp.p_occo = ptp.p_occ = float(line[54:60].strip()) if len(line) > 60 else 0.0
        ptp.p_bvo = ptp.p_bval = float(line[60:66].strip()) if len(line) > 66 else 0.0
        
        # Extract atom ID (columns 77-78)
        ptp.p_atomid = line[76:78].strip() if len(line) > 78 else ""
        
    except (ValueError, IndexError) as e:
        raise ValueError(f"Error parsing PDB line: {e}\nLine: {line}")

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

def procargs(args, top: TotalSt, tfpl) -> None:
    """Process command line arguments using argparse. File objects are passed in from main."""
    # Initialize top structure
    top.tmeanbw = top.tnhet = top.tnwatat = 0
    top.tnatoms = top.tnpat = top.tnpwat = top.tnwaters = 0
    top.tfpl = tfpl
    top.tnch = 0  # Initialize chain count to 0
    top.tchs = [Chain() for _ in range(50)]  # Allocate 50 chains
    top.t1ch_s = args.S
    top.thighbq = args.H
    top.tmaxhbsq = MAXHBONDSQ
    top.tminhbsq = MINHBONDSQ
    top.tminhsq = CLOSEHSQ
    top.tminhet = CLOSEHET
    top.tbump = args.B

    # Set parameters based on command line options
    if args.X and 0.05 < args.X < 10.0:
        top.tmaxhbsq = args.X * args.X
    if args.M and 0.05 < args.M < 10.0:
        top.tminhbsq = args.M * args.M
    if args.L and 0.05 < args.L < 10.0:
        top.tminhsq = args.L * args.L
    if args.O and 0.05 < args.O < 10.0:
        top.tminhet = args.O * args.O

    if tfpl is not None:
        tfpl.write(f"Run of {sys.argv[0]}\n")

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
        if hasattr(top.tfpi, 'seek') and callable(top.tfpi.seek):
            try:
                top.tfpi.seek(0)
            except Exception as e:
                if top.tfpl and hasattr(top.tfpl, 'write'):
                    top.tfpl.write(f"Error seeking input file in ready(): {e}\n")
                # If seek fails, try to reopen as a fallback if possible
                file_name_to_reopen = None
                if hasattr(top.tfpi, 'name') and isinstance(top.tfpi.name, str) and top.tfpi.name and top.tfpi.name != '<stdin>':
                    file_name_to_reopen = top.tfpi.name
                
                try:
                    if hasattr(top.tfpi, 'close') and callable(top.tfpi.close):
                        top.tfpi.close()
                except Exception:
                    pass

                if file_name_to_reopen:
                    try:
                        top.tfpi = open(file_name_to_reopen, 'r')
                    except Exception as e_reopen:
                        if top.tfpl and hasattr(top.tfpl, 'write'):
                            top.tfpl.write(f"Error reopening file '{file_name_to_reopen}' in ready() after seek failed: {e_reopen}\n")
                        return 1 # Indicate error
                else:
                    if top.tfpl and hasattr(top.tfpl, 'write'):
                        top.tfpl.write(f"Error: Input stream is not sys.stdin, seek failed, and no valid name to reopen in ready().\n")
                    return 1 # Indicate error
        else:
            # If we can't seek (e.g. not supported, or missing attribute), we need to try to reopen the file by name
            file_name_to_reopen = None
            if hasattr(top.tfpi, 'name') and isinstance(top.tfpi.name, str) and top.tfpi.name and top.tfpi.name != '<stdin>':
                file_name_to_reopen = top.tfpi.name
            
            try:
                if hasattr(top.tfpi, 'close') and callable(top.tfpi.close):
                    top.tfpi.close()
            except Exception:
                pass # Ignore errors on close, as we are trying to reopen anyway

            if file_name_to_reopen:
                try:
                    top.tfpi = open(file_name_to_reopen, 'r')
                except Exception as e:
                    if top.tfpl and hasattr(top.tfpl, 'write'):
                        top.tfpl.write(f"Error reopening file '{file_name_to_reopen}' in ready(): {e}\n")
                    return 1 # Indicate error
            else:
                # Cannot reopen if it's not a named file (e.g. some other stream)
                if top.tfpl and hasattr(top.tfpl, 'write'):
                    top.tfpl.write(f"Error: Input stream is not sys.stdin, not seekable, and has no valid name to reopen in ready().\n")
                return 1 # Indicate error
    
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
    for pwap in top.tpwa[:top.tpwap]: # Iterate over active waters
        if pwap.p_resnum > maxres:
            maxres = pwap.p_resnum
    
    if maxres < -999:
        return 0
    
    for pwap in top.tpwa[:top.tpwap]: # Iterate over active waters
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
        
        for j in range(top.tpwap): # This inner loop index 'j' is necessary for comparison
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
    """Analyze neighbors of a water molecule and detect quartets"""
    if p0 is None or p0.p_nbr is None:
        return
    
    # Get the chain of water molecules
    p1 = p0.p_nbr
    if p1 is None:
        return
    
    # Check if p1 is already a conformer
    if p1.p_conf != ' ':
        return
    
    # Calculate distance between p0 and p1
    dist = pdbdist(p0, p1)
    
    # If they're too far apart, return
    if dist > top.tmaxhbsq:
        return
    
    # Get p2, the neighbor of p1 that's not p0
    p2 = p1.p_nbr
    if p2 is None or p2 == p0:
        # No third water, relate p0 and p1 as a pair
        if dist < top.tminhbsq:
            relateem(top, p0, p1)
        return
    
    # Check if p2 is already a conformer
    if p2.p_conf != ' ':
        # Relate p0 and p1 as a pair
        if dist < top.tminhbsq:
            relateem(top, p0, p1)
        return
    
    # Get p3, the neighbor of p2 that's not p1 and not p0
    p3 = p2.p_nbr
    if p3 is None or p3 == p1 or p3 == p0:
        # No fourth water, handle as a trio
        # Calculate distances
        d01 = dist
        d12 = pdbdist(p1, p2)
        
        # If all distances are small enough, relate all three
        if d01 < top.tminhbsq and d12 < top.tminhbsq:
            # Apply conformer grouping to trio
            p0.p_conf = 'A'
            p1.p_conf = 'B'
            p2.p_conf = 'C'
            p0.p_nconfs = p1.p_nconfs = p2.p_nconfs = 3
            
            # Adjust occupancies - for a trio, use 0.33 each
            p0.p_occ = p1.p_occ = p2.p_occ = 0.33
            
            # Average and scale B values
            avg_b = (p0.p_bval + p1.p_bval + p2.p_bval) / 3.0
            p0.p_bval = p1.p_bval = p2.p_bval = 0.9 * avg_b
            
            # Report the triplet
            reportmult(top, p0, 0)
        else:
            # If only two are close enough, relate those two
            if d01 < top.tminhbsq:
                relateem(top, p0, p1)
            if d12 < top.tminhbsq:
                relateem(top, p1, p2)
        return
    
    # Check if p3 is already a conformer
    if p3.p_conf != ' ':
        # Handle as a trio (p0, p1, p2)
        d01 = dist
        d12 = pdbdist(p1, p2)
        
        if d01 < top.tminhbsq and d12 < top.tminhbsq:
            p0.p_conf = 'A'
            p1.p_conf = 'B'
            p2.p_conf = 'C'
            p0.p_nconfs = p1.p_nconfs = p2.p_nconfs = 3
            p0.p_occ = p1.p_occ = p2.p_occ = 0.33
            avg_b = (p0.p_bval + p1.p_bval + p2.p_bval) / 3.0
            p0.p_bval = p1.p_bval = p2.p_bval = 0.9 * avg_b
            reportmult(top, p0, 0)
        else:
            if d01 < top.tminhbsq:
                relateem(top, p0, p1)
            if d12 < top.tminhbsq:
                relateem(top, p1, p2)
        return
    
    # Check if this is a closed quartet (p3 points back to p0)
    if p3.p_nbr == p0:
        # We have a quartet - analyze and adjust
        # Calculate all distances
        d01 = dist
        d12 = pdbdist(p1, p2)
        d23 = pdbdist(p2, p3)
        d30 = pdbdist(p3, p0)
        
        # If all distances are small, mark as a quartet
        if d01 < top.tmaxhbsq and d12 < top.tmaxhbsq and d23 < top.tmaxhbsq and d30 < top.tmaxhbsq:
            # Check if we should split into 2+2
            split4(top, p0)
            
            # Determine if we should adjust as quartet or reorganize
            if (d01 < top.tminhbsq and d12 < top.tminhbsq and 
                d23 < top.tminhbsq and d30 < top.tminhbsq):
                # All are close enough - adjust as a full quartet
                adjustqb(top, p0)
            else:
                # Not all are close enough - reorganize based on distances
                reorg4(top, p0, p1, p2, p3)
        else:
            # Not all are in range, handle as pairs
            if d01 < top.tminhbsq:
                relateem(top, p0, p1)
            if d12 < top.tminhbsq:
                relateem(top, p1, p2)
            if d23 < top.tminhbsq:
                relateem(top, p2, p3)
            if d30 < top.tminhbsq:
                relateem(top, p3, p0)
    else:
        # Not a closed quartet - handle as separate parts
        d01 = dist
        d12 = pdbdist(p1, p2)
        d23 = pdbdist(p2, p3)
        
        # Handle as pairs or triplet + single
        if d01 < top.tminhbsq and d12 < top.tminhbsq and d23 < top.tminhbsq:
            # All three pairs are close - handle as triplet and single
            p0.p_conf = 'A'
            p1.p_conf = 'B'
            p2.p_conf = 'C'
            p0.p_nconfs = p1.p_nconfs = p2.p_nconfs = 3
            p0.p_occ = p1.p_occ = p2.p_occ = 0.33
            avg_b = (p0.p_bval + p1.p_bval + p2.p_bval) / 3.0
            p0.p_bval = p1.p_bval = p2.p_bval = 0.9 * avg_b
            reportmult(top, p0, 0)
            
            # p3 remains separate
        else:
            # Handle as pairs
            if d01 < top.tminhbsq:
                relateem(top, p0, p1)
            if d12 < top.tminhbsq:
                relateem(top, p1, p2)
            if d23 < top.tminhbsq:
                relateem(top, p2, p3)

def meanbw(top: TotalSt) -> None:
    """Calculate mean B value for waters"""
    if top.tpwap <= 0:
        top.tmeanbw = 0.0
        return
    
    sum_b = 0.0
    count = 0
    
    for pwap in top.tpwa[:top.tpwap]: # Iterate over active waters
        sum_b += pwap.p_bval
        count += 1
    
    if count > 0:
        top.tmeanbw = sum_b / count
        if top.tfpl is not None:
            top.tfpl.write(f"Mean B value for waters = {top.tmeanbw:.2f}\n")

def adjustqb(top: TotalSt, p0: PDBRecord) -> None:
    """Adjusts quality (occupancy) and B values for multiple conformers"""
    if p0 is None or p0.p_nbr is None:
        return
    
    # Get the neighboring waters
    p1 = p0.p_nbr
    if p1 is None:
        return
    
    p2 = p1.p_nbr if p1 is not None and p1.p_nbr is not None and p1.p_nbr != p0 else None
    p3 = p2.p_nbr if p2 is not None and p2.p_nbr is not None and p2.p_nbr != p1 and p2.p_nbr != p0 else None
    
    # Array to hold the records for processing
    records = [p0]
    if p1 is not None:
        records.append(p1)
    if p2 is not None:
        records.append(p2)
    if p3 is not None:
        records.append(p3)
    
    # Count total number of records
    nrecs = len(records)
    if nrecs <= 1:
        return
    
    # Calculate average B value
    total_b = sum(p.p_bval for p in records)
    avg_b = total_b / nrecs
    
    # Calculate average occupancy
    total_occ = sum(p.p_occ for p in records)
    avg_occ = total_occ / nrecs
    
    # Adjust B values and occupancies
    scale_factor = 0.9
    for i, p in enumerate(records):
        # Set conformer ID based on position
        p.p_conf = chr(ord('A') + i)
        p.p_nconfs = nrecs
        
        # Adjust B value - use a scaled average
        new_b = scale_factor * avg_b
        if new_b > 80.0:
            print(f"adjustqb: B={new_b:.2f} to {p.p_conf}{p.p_resname} {p.p_chainid}{p.p_resnum:4d}")
            new_b = 80.0
        p.p_bval = new_b
        
        # Adjust occupancy - for pairs, use 0.5/0.5, for triplets or more, distribute evenly
        if nrecs == 2:
            p.p_occ = 0.5
        elif nrecs == 3:
            p.p_occ = 0.33
        elif nrecs == 4:
            p.p_occ = 0.25
        else:
            p.p_occ = 1.0 / nrecs
    
    # Report the adjustment
    reportmult(top, p0, 2)  # Report as "Adjusted"

def split4(top: TotalSt, pwap: PDBRecord) -> None:
    """Splits a group of four waters into two groups of 2 if possible"""
    if pwap is None or pwap.p_nbr is None:
        return
    
    # Get the 4 waters
    p0 = pwap
    p1 = p0.p_nbr
    if p1 is None:
        return
    
    p2 = p1.p_nbr
    if p2 is None or p2 == p0:
        return
    
    p3 = p2.p_nbr
    if p3 is None or p3 == p0 or p3 == p1:
        return
    
    # Check if this is a cycle of 4
    if p3.p_nbr != p0:
        return
    
    # Calculate pairwise distances
    d01 = pdbdist(p0, p1)
    d12 = pdbdist(p1, p2)
    d23 = pdbdist(p2, p3)
    d30 = pdbdist(p3, p0)
    d02 = pdbdist(p0, p2)
    d13 = pdbdist(p1, p3)
    
    # Find the min and max distances
    all_dists = [d01, d12, d23, d30, d02, d13]
    min_dist = min(all_dists)
    max_dist = max(all_dists)
    
    # If the ratio of max to min is too high, split the group
    if max_dist / min_dist > 1.5:
        # Determine which pair to split
        if d02 == max_dist:
            # Split between 0-2
            p0.p_nbr = p3
            p2.p_nbr = p1
            p1.p_nbr = p2
            p3.p_nbr = p0
        elif d13 == max_dist:
            # Split between 1-3
            p0.p_nbr = p1
            p1.p_nbr = p0
            p2.p_nbr = p3
            p3.p_nbr = p2
        elif d01 == max_dist:
            # Split between 0-1
            p0.p_nbr = p3
            p3.p_nbr = p2
            p2.p_nbr = p3
            p1.p_nbr = None  # Break the cycle
        elif d12 == max_dist:
            # Split between 1-2
            p0.p_nbr = p1
            p1.p_nbr = p0
            p2.p_nbr = None
            p3.p_nbr = None
        elif d23 == max_dist:
            # Split between 2-3
            p0.p_nbr = p1
            p1.p_nbr = p2
            p2.p_nbr = p1
            p3.p_nbr = None
        elif d30 == max_dist:
            # Split between 3-0
            p1.p_nbr = p0
            p0.p_nbr = p1
            p2.p_nbr = p3
            p3.p_nbr = p2
        
        # Report the split
        top.tfpl.write(f"Split four waters: {p0.p_chainid}{p0.p_resnum}, {p1.p_chainid}{p1.p_resnum}, "
                       f"{p2.p_chainid}{p2.p_resnum}, {p3.p_chainid}{p3.p_resnum}\n")
        top.tfpl.write(f"  Max distance: {math.sqrt(max_dist):.2f} Å, Min distance: {math.sqrt(min_dist):.2f} Å\n")

def reorg4(top: TotalSt, p0: PDBRecord, p1: PDBRecord, p2: PDBRecord, p3: PDBRecord) -> None:
    """Reorganizes four consecutive PDBRECORDs that point to each other"""
    if p0 is None or p1 is None or p2 is None or p3 is None:
        return
    
    # Calculate all pairwise distances
    d01 = pdbdist(p0, p1)
    d02 = pdbdist(p0, p2)
    d03 = pdbdist(p0, p3)
    d12 = pdbdist(p1, p2)
    d13 = pdbdist(p1, p3)
    d23 = pdbdist(p2, p3)
    
    # Create a list of all distances and their indices
    distances = [
        (d01, 0, 1), (d02, 0, 2), (d03, 0, 3),
        (d12, 1, 2), (d13, 1, 3), (d23, 2, 3)
    ]
    
    # Sort by distance
    distances.sort(key=lambda x: x[0])
    
    # Check if this is already a 2+2 split pattern
    if (p0.p_nbr == p1 and p1.p_nbr == p0 and p2.p_nbr == p3 and p3.p_nbr == p2) or \
       (p0.p_nbr == p2 and p2.p_nbr == p0 and p1.p_nbr == p3 and p3.p_nbr == p1) or \
       (p0.p_nbr == p3 and p3.p_nbr == p0 and p1.p_nbr == p2 and p2.p_nbr == p1):
        return
    
    # Otherwise organize into the 3 closest pairs
    # Start with the first (closest) pair
    i, j = distances[0][1], distances[0][2]
    waters = [p0, p1, p2, p3]
    waters[i].p_nbr = waters[j]
    waters[j].p_nbr = waters[i]
    
    # Find the next closest non-overlapping pair
    for dist, idx1, idx2 in distances[1:]:
        if idx1 != i and idx1 != j and idx2 != i and idx2 != j:
            waters[idx1].p_nbr = waters[idx2]
            waters[idx2].p_nbr = waters[idx1]
            break
    
    top.tfpl.write(f"Reorganized four waters: {p0.p_chainid}{p0.p_resnum}, {p1.p_chainid}{p1.p_resnum}, "
                   f"{p2.p_chainid}{p2.p_resnum}, {p3.p_chainid}{p3.p_resnum}\n")
    top.tfpl.write(f"  Paired: ({waters[i].p_chainid}{waters[i].p_resnum}-{waters[j].p_chainid}{waters[j].p_resnum})"
                   f" and ({waters[idx1].p_chainid}{waters[idx1].p_resnum}-{waters[idx2].p_chainid}{waters[idx2].p_resnum})\n")
    
    # Report the reorganization
    reportmult(top, p0, 2)

def reduceocc(top: TotalSt) -> None:
    """Reduce occupancy of waters with B > 1.2*<B>"""
    if top.tmeanbw <= 0.0:
        return
    
    bthresh = 1.2 * top.tmeanbw
    nreduced = 0
    
    for pwap in top.tpwa[:top.tpwap]: # Iterate over active waters
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
    for chain_idx in range(top.tnch): # Iterate over active chains
        chain = top.tchs[chain_idx]
        chain.c_wat0 = 0
        chain.c_watl1 = 0
        chain.c_watm0 = 0
        chain.c_watl = 0
        chain.c_minwat = 1  # Default starting water number
    
    # Assign each water to the nearest protein chain
    for pwap in top.tpwa[:top.tpwap]: # Iterate over active waters
        nearest_chain_id = None # Store ID, not object, to avoid issues if chain list changes
        min_dist = float('inf')
        
        # Find the nearest protein chain
        for pat in top.tpat[:top.tpatp]: # Iterate over active protein atoms
            if pat.p_chainid != ' ':  # Skip atoms without a chain
                dist = pdbdist(pwap, pat)
                if dist < min_dist:
                    min_dist = dist
                    nearest_chain_id = pat.p_chainid
        
        # Assign the water to the nearest chain
        if nearest_chain_id is not None:
            pwap.p_chainid = nearest_chain_id
            
            # Update the chain's water count
            for chain_idx in range(top.tnch): # Iterate over active chains
                if top.tchs[chain_idx].c_chainid == nearest_chain_id:
                    top.tchs[chain_idx].c_watl += 1
                    break

def adjustmult(top: TotalSt) -> None:
    """Count and modify multiple conformers"""
    # Count the number of waters in each chain
    chain_counts: Dict[str, int] = {} # Added type hint for clarity
    for pwap in top.tpwa[:top.tpwap]: # Iterate over active waters
        if pwap.p_chainid not in chain_counts:
            chain_counts[pwap.p_chainid] = 0
        chain_counts[pwap.p_chainid] += 1
    
    # Initialize residue numbers for each chain
    for chain_idx in range(top.tnch): # Iterate over active chains
        chain = top.tchs[chain_idx]
        chain.c_curwat = chain.c_minwat
    
    # Assign residue numbers to waters based on their chain
    # This loop needs index access for comparison with top.tpwa[i-1]
    for i in range(top.tpwap):
        pwap = top.tpwa[i]
        for chain_idx in range(top.tnch): # Iterate over active chains
            chain = top.tchs[chain_idx]
            if chain.c_chainid == pwap.p_chainid:
                # If this is a new conformer (A or B), assign it the same residue number
                if pwap.p_conf in ['A', 'B'] and i > 0 and \
                   top.tpwa[i-1].p_conf in ['A', 'B'] and \
                   top.tpwa[i-1].p_chainid == pwap.p_chainid and \
                   top.tpwa[i-1].p_resnum == pwap.p_resnum:
                    continue  # Keep the same residue number
                else:
                    # Assign a new residue number
                    pwap.p_resnum = chain.c_curwat
                    chain.c_curwat += 1
                break

def sortmults(top: TotalSt) -> None:
    """Reorder multiple-conformer waters within each chain"""
    # Group waters by chain and residue number
    chain_res_groups: Dict[Tuple[str, int], List[PDBRecord]] = {} # Added type hint
    
    for pwap in top.tpwa[:top.tpwap]: # Iterate over active waters
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
    for pwap in top.tpwa[:top.tpwap]: # Iterate over active waters
        
        # Skip waters that are already marked for special handling
        if pwap.p_diag != 0:
            continue
        
        # Check against all protein atoms
        for pat in top.tpat[:top.tpatp]: # Iterate over active protein atoms
            
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
    unmarked_waters = []
    single_conf_marked = []
    multi_conf = []
    
    for pwap in top.tpwa[:top.tpwap]: # Iterate over active waters
        if pwap.p_conf == ' ':
            unmarked_waters.append(pwap)
        elif pwap.p_nconfs == 1:
            single_conf_marked.append(pwap)
        else:
            multi_conf.append(pwap)
    
    # If there are no single-conformation marked waters, nothing to do
    if not single_conf_marked:
        return
    
    # Create a new ordered list of waters
    new_waters = unmarked_waters + single_conf_marked + multi_conf
        
    # Replace the original water list with the reordered one, only up to tpwap
    for i in range(min(len(new_waters), top.tpwap)):
        top.tpwa[i] = new_waters[i]

def finalmult(top: TotalSt) -> None:
    """Final count of multiple conformers"""
    # Count waters by conformer type
    unmarked_count = 0
    a_conf_count = 0
    b_conf_count = 0
    other_conf_count = 0
    
    for pwap in top.tpwa[:top.tpwap]: # Iterate over active waters
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

def warn4(p0: PDBRecord, p1: PDBRecord, whic: int) -> None:
    """Prints a warning about water molecules"""
    if whic == 0:  # More than 2 conformers
        print(f"More than 2 {p0.p_resname} conformers at {p0.p_chainid}{p0.p_resnum}")
    elif whic == 1:  # Overlapping waters
        print(f"Water overlap: {p0.p_chainid}{p0.p_resnum} and {p1.p_chainid}{p1.p_resnum}")
    elif whic == 2:  # Multiple close waters
        print(f"Multiple close waters at {p0.p_chainid}{p0.p_resnum}")
    elif whic == 3:  # B-value disparity
        print(f"B-value ratio > 2: {p0.p_chainid}{p0.p_resnum}({p0.p_bval:.1f}) and {p1.p_chainid}{p1.p_resnum}({p1.p_bval:.1f})")

def reportmult(top: TotalSt, p0: PDBRecord, reporttype: int) -> None:
    """Describes a duo, trio, or quartet of water molecules"""
    if top.tfpl is None:
        return
    
    # Get water positions
    p1 = p0.p_nbr if p0.p_nbr is not None else None
    p2 = p1.p_nbr if p1 is not None and p1.p_nbr is not None and p1.p_nbr != p0 else None
    p3 = p2.p_nbr if p2 is not None and p2.p_nbr is not None and p2.p_nbr != p1 and p2.p_nbr != p0 else None
    
    # Check for cyclic relationships
    if p3 is not None and p3.p_nbr is not None and p3.p_nbr != p2 and p3.p_nbr != p1 and p3.p_nbr != p0:
        # More than 4 in a cycle - report only first 4
        p3.p_nbr = None
    
    # Format the report string based on the number of linked waters and report type
    report_prefix = ""
    if reporttype == 0:    # Initial/Diagnostic
        report_prefix = "Detected"
    elif reporttype == 1:  # Already labeled
        report_prefix = "Already labeled"
    elif reporttype == 2:  # Adjusted/Modified
        report_prefix = "Adjusted"
    
    if p3 is not None:
        # Quartet
        top.tfpl.write(f"{report_prefix} quartet of waters: {p0.p_chainid}{p0.p_resnum}, {p1.p_chainid}{p1.p_resnum}, "
                       f"{p2.p_chainid}{p2.p_resnum}, {p3.p_chainid}{p3.p_resnum}\n")
    elif p2 is not None:
        # Triplet
        top.tfpl.write(f"{report_prefix} triplet of waters: {p0.p_chainid}{p0.p_resnum}, {p1.p_chainid}{p1.p_resnum}, "
                       f"{p2.p_chainid}{p2.p_resnum}\n")
    elif p1 is not None:
        # Duo
        top.tfpl.write(f"{report_prefix} pair of waters: {p0.p_chainid}{p0.p_resnum}, {p1.p_chainid}{p1.p_resnum}\n")
    else:
        # Single
        top.tfpl.write(f"{report_prefix} single water: {p0.p_chainid}{p0.p_resnum}\n")
    
    # Report distances if it's a diagnostic report
    if reporttype == 0:
        if p1 is not None:
            dist = math.sqrt(pdbdist(p0, p1))
            top.tfpl.write(f"  Distance {p0.p_chainid}{p0.p_resnum}-{p1.p_chainid}{p1.p_resnum}: {dist:.2f} Å\n")
        if p2 is not None:
            dist = math.sqrt(pdbdist(p1, p2))
            top.tfpl.write(f"  Distance {p1.p_chainid}{p1.p_resnum}-{p2.p_chainid}{p2.p_resnum}: {dist:.2f} Å\n")
        if p3 is not None:
            dist = math.sqrt(pdbdist(p2, p3))
            top.tfpl.write(f"  Distance {p2.p_chainid}{p2.p_resnum}-{p3.p_chainid}{p3.p_resnum}: {dist:.2f} Å\n")

def main() -> int:
    """Main function with improved error handling and file management."""
    parser = argparse.ArgumentParser(description='Process PDB files for water analysis.')
    parser.add_argument('-S', action='store_true', help='Single chain mode')
    parser.add_argument('-H', action='store_true', help='High B quality mode')
    parser.add_argument('-B', action='store_true', help='Bump mode')
    parser.add_argument('-X', type=float, help='Maximum H-bond distance')
    parser.add_argument('-M', type=float, help='Minimum H-bond distance')
    parser.add_argument('-L', type=float, help='Minimum hydrogen distance')
    parser.add_argument('-O', type=float, help='Minimum heteroatom distance')
    parser.add_argument('input_file', nargs='?', default='-', help='Input PDB file (default: stdin)')
    parser.add_argument('output_file', nargs='?', default='-', help='Output PDB file (default: stdout)')
    args = parser.parse_args()

    tfpi_obj = None
    tfpwat_obj = None
    tfpl_obj = None

    try:
        # Open files
        if args.input_file == '-' or args.input_file == '':
            tfpi_obj = sys.stdin
        else:
            tfpi_obj = open(args.input_file, 'r')
        
        if args.output_file == '-' or args.output_file == '':
            tfpwat_obj = sys.stdout
        else:
            tfpwat_obj = open(args.output_file, 'w')
        
        tfpl_obj = open("closewat.log", "w")

        tos = TotalSt()
        tos.tfpi = tfpi_obj
        tos.tfpwat = tfpwat_obj
        tos.tfpl = tfpl_obj
        
        procargs(args, tos, tfpl_obj)

        input_lines = []
        is_stdin = (tos.tfpi == sys.stdin)

        if is_stdin:
            input_lines = list(tos.tfpi) # Read all lines from stdin at once

        # First pass: read the PDB file to count atoms and chains
        iterator1 = input_lines if is_stdin else tos.tfpi
        for line in iterator1:
            getpdblin(line, tos)

        # Prepare for second pass
        if ready(tos) != 0:
            # Error in ready(), files will be closed in finally
            return 1

        # Second pass: process the PDB file
        iterator2 = input_lines if is_stdin else tos.tfpi
        for line in iterator2:
            try:
                procpdblin(line, tos)
            except ValueError as e:
                print(f"Warning: {e}", file=sys.stderr)
                if tos.tfpl: 
                    tos.tfpl.write(f"Warning processing line: {e}\n")
                continue

        tos.tnwaters = tos.tpwap
        noconformers(tos)
        tos.tpate = tos.tpatp
        closestwat(tos)
        
        # Loop for thisthird - original indexing seems important for p_nbr logic.
        for i in range(tos.tnwaters): 
            pwap = tos.tpwa[i]
            if pwap.p_nbr is not None and pwap.p_conf == ' ':
                thisthird(tos, pwap)
        
        meanbw(tos)
        if tos.thighbq:
            reduceocc(tos)
        
        # Sort the active part of tpwa
        if tos.tnwaters > 0: # Only sort if there are waters
            active_waters = tos.tpwa[:tos.tnwaters]
            active_waters.sort(key=get_sort_key_occb)
            tos.tpwa[:tos.tnwaters] = active_waters

        makechains(tos)
        adjustmult(tos)
        
        if tos.tnwaters > 0: # Only sort if there are waters
            active_waters = tos.tpwa[:tos.tnwaters]
            active_waters.sort(key=get_sort_key_chain)
            tos.tpwa[:tos.tnwaters] = active_waters
        
        sortmults(tos)
        
        tos.talert = [0] * 8 # Reset alerts
        
        if len(tos.tchs) == 1 and tos.t1ch_s and tos.tchs[0].c_chainid != 'S':
            # Ensure tnwaters is positive before accessing tpwa[0]
            if tos.tnwaters > 0:
                first_resnum_offset = tos.tpwa[0].p_resnum - 1 
                for k_idx in range(tos.tnwaters): 
                    tos.tpwa[k_idx].p_chainid = 'S'
                    tos.tpwa[k_idx].p_resnum -= first_resnum_offset
        
        proximity(tos)
        insert_singles(tos) # Reorders tos.tpwa up to tpwap (which is tnwaters)
        finalmult(tos)
        
        for k_idx in range(tos.tnwaters): 
            outrec(tos.tpwa[k_idx], tos.tfpwat)
        
        return 0
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc(file=sys.stderr)
        return 1
    finally:
        if tfpl_obj is not None and not tfpl_obj.closed:
            tfpl_obj.close()
        # Check if tfpi_obj was opened and is not stdin before closing
        if tfpi_obj is not None and tfpi_obj != sys.stdin and not tfpi_obj.closed:
            tfpi_obj.close()
        # Check if tfpwat_obj was opened and is not stdout before closing
        if tfpwat_obj is not None and tfpwat_obj != sys.stdout and not tfpwat_obj.closed:
            tfpwat_obj.close()

if __name__ == "__main__":
    sys.exit(main())
