#include	<stdio.h>
#include	<stdlib.h>
#include	<math.h>
#include	<string.h>
#define	MAXHBONDSQ	(3.2*3.2)	/* default max H-bond dist squared */
#define	MINHBONDSQ	(2.5*2.5)	/* default min H-bond dist squared */
#define	CLOSEHSQ	(1.5 * 1.5)	/* closest approach to hydrogen */
#define	CLOSEHET	(3.9 * 3.9)	/* default min dist to heteroatom */
#define	NUMGAP		10		/* gap in residue number */
#define	HBARELY		0.0891
#define	OBARELY		0.1984

/* This reads the a PDB and finds the water molecules therein.
 It sorts them on occupancy and B value, and then analyzes them in terms
 of which macromolecular chain the water molecules are closest to.
 It then renumbers the waters so that they're associated with their closest
 macromolecular chain. It also recalculates occupancies and B values
 for waters that have been bundled or unbundled. */

typedef unsigned short	zushort;
typedef int (*QPE)(const void *, const void *);

typedef struct {
		char		p_rtype[8];
		char		p_attype[8];
		char		p_resname[8];
		char		p_atomid[4];
		char		p_conf, p_chainid, p_nconfs, p_chaino;
		char		p_confo, p_nconfo, p_diag, p_rff2;
		float		p_xc, p_yc, p_zc, p_occ, p_bval;
		float		p_occo, p_bvo, p_dfc;
		unsigned short	p_atnum, p_rff3;
		short		p_resnum, p_resno;
		long		*p_nbr;
		} PDBRECORD;

typedef struct {
		char		c_chainid, c_crffu;
		short		c_minwat;
		short		c_curwat;
		short		c_wat0, c_watl1, c_watm0, c_watl;
		} CHAIN;

typedef struct {
		zushort		tnhet, tnwatat, tnwathet, tnch;
		zushort 	tnatoms, tnpat, tnpwat, tnwaters;
		zushort		t1ch_s, thighbq, tnclose, tbump;
		zushort		talert[8];
		float		tmeanbw, tmaxhbsq, tminhbsq, tminhsq, tminhet;
		PDBRECORD	*tpat, *tpatp, *tpwa, *tpwap, *tpate;
		CHAIN		*tchs;
		FILE		*tfpi, *tfpwat, *tfpl;
		} TOTALST;

extern void exit(int exval);

static char	*diagstr[] = { "faraway", " ok now", "  metal",
		"  leave", "altconf", "replace", "   bump", "   edit"};

static void outrec(PDBRECORD *prec, FILE *fp)
{ /* writes out a PDBRECORD structure as a string */
	fprintf(fp,
	 "%6s%5d %4s%c%3s %c%4hd    %8.3f%8.3f%8.3f%6.2f%6.2f",
		prec->p_rtype, (int)(prec->p_atnum), prec->p_attype,
		prec->p_conf, prec->p_resname, prec->p_chainid, prec->p_resnum,
		prec->p_xc, prec->p_yc, prec->p_zc,
		prec->p_occ, prec->p_bval);
	if (prec->p_atomid[1] == '\0')
	  {
		prec->p_atomid[1] = prec->p_atomid[0];
		prec->p_atomid[0] = ' ';
	  }
	fprintf(fp, "          %c%c  \n",
		prec->p_atomid[0], prec->p_atomid[1]);
}

static double pdbdist(PDBRECORD *p0, PDBRECORD *p1)
{ /* This calculates the distance squared between two atoms */
	double	dr;

	dr =	(p0->p_xc - p1->p_xc) * (p0->p_xc - p1->p_xc) +
		(p0->p_yc - p1->p_yc) * (p0->p_yc - p1->p_yc) +
		(p0->p_zc - p1->p_zc) * (p0->p_zc - p1->p_zc);
	return dr;
}

static int isconformer(TOTALST *top,
	PDBRECORD *p0, PDBRECORD *p1, int checkdist)
{ /* returns true only if p1 is a later conformer than p0.
   If checkdist is nonzero, confirm that the p0 and p1 are near each other */
	int	dif;

	if ((NULL == p0) || (NULL == p1)) return 0;
	if (p0->p_chainid != p1->p_chainid) return 0;
	if (p0->p_resnum != p1->p_resnum) return 0;
	if ((p0->p_conf != 'A') && (p0->p_conf != 'B') && (p0->p_conf != 'C'))
		return 0;
	if ((p1->p_conf != 'B') && (p1->p_conf != 'C') && (p1->p_conf != 'D'))
		return 0;
	dif = p1->p_conf - p0->p_conf;
	if ((dif < 1) || (dif > 3)) return 0;
	if (fabs(p0->p_bval - p1->p_bval) > 0.005) return 0;
	if (checkdist) return (pdbdist(p0, p1) > top->tminhbsq) ? 0 : 1;
	else	return 1;
}

static int occsort(PDBRECORD *pr0, PDBRECORD *pr1)
{ /* compare routine for sorting on the basis of occupancy */
	return (pr0->p_occ < pr1->p_occ) ? 1 : -1;
}

static int occbsort(PDBRECORD *pr0, PDBRECORD *pr1)
{ /* compare routine for sorting, on conformer, then
	on occupancy, and finally on B value. */
	int	dif;

	if (pr1->p_conf != ' ')
	  {
		if (pr0->p_conf != ' ')
		  {
			if ((pr0->p_resnum == pr1->p_resnum) &&
			 (pr0->p_chainid == pr1->p_chainid))
			  {
				dif = ((int)(pr1->p_conf)) -
					((int)(pr0->p_conf));
				return (dif > 0) ? -1 : ((dif == 0) ? 0 : 1);
			  }
			else	return (pr0->p_bval < pr1->p_bval) ? -1 : 1;
		  }
		else return -1;
	  }
	if (pr0->p_conf != ' ') return 1;
	if ((pr0->p_conf == ' ') && (pr1->p_conf != ' ')) return -1;
	if ((pr1->p_conf == ' ') && (pr0->p_conf != ' ')) return 1;
	if (pr0->p_occ > pr1->p_occ) return -1;
	if (pr0->p_occ < pr1->p_occ) return 1;
	return (pr0->p_bval < pr1->p_bval) ? -1 : 1;
}

static int chrsort(PDBRECORD *pr0, PDBRECORD *pr1)
{ /* compare routine for sorting, on chain ID and secondarily on residue # */
	if (pr0->p_chainid < pr1->p_chainid) return -1;
	else if (pr0->p_chainid > pr1->p_chainid) return 1;
	if (pr0->p_resnum < pr1->p_resnum) return -1;
	else if (pr0->p_resnum > pr1->p_resnum) return 1;
	return (pr0->p_conf < pr1->p_conf) ? -1 : 1;
}

static int bonly(PDBRECORD *pr0, PDBRECORD *pr1)
{ /* sort routine to operate only based on conformer letter and B value */
	int	dif;

	if ((pr0->p_resnum == pr1->p_resnum) && (pr0->p_conf != pr1->p_conf))
	  {
		dif = ((int)(pr0->p_conf)) - ((int)(pr1->p_conf));
		return (dif < 0) ? -1 : 1;
	  }
	return (pr0->p_bval <= pr1->p_bval) ? -1 : 1;
}

static void strtorec(char *line, PDBRECORD *ptp)
{ /* Converts a PDB text string into a PDB record */
	int	i, j;
	char	tchar[40];
	char	legalconf[6] = " ABCD";

	ptp->p_nbr = (long *)(PDBRECORD *)NULL;
	for (i = 0; i < 6; i++) ptp->p_rtype[i] = line[i];
	ptp->p_rtype[6] = '\0';
	for (i = 0; i < 5; i++) tchar[i] = line[i+6];
	tchar[5] = '\0';	ptp->p_atnum = atoi(tchar);
	for (i = 0; i < 4; i++) ptp->p_attype[i] = line[i+12];
	ptp->p_attype[4] = '\0';
	ptp->p_confo = ptp->p_conf = line[16];
	ptp->p_nconfo = ptp->p_nconfs = 1; /* this can be modified later */
	for (j = -1, i = 0; i < 5; i++) /* only legal conformers are A,B,C,D */
	   {
		if (ptp->p_conf == legalconf[i])
		  {
			j = i;	break;
		  }
	   }
	if (j == -1) ptp->p_conf = ' ';
	if (ptp->p_conf != ' ') ptp->p_nconfo = 2; /* first stab */
	for (i = 0; i < 3; i++) ptp->p_resname[i] = line[i+17];
	ptp->p_resname[3] = ptp->p_resname[4] = '\0';
	ptp->p_chaino = ptp->p_chainid = line[21];
	for (i = 0; i < 4; i++) tchar[i] = line[i+22];
	tchar[4] = '\0';
	ptp->p_resno = ptp->p_resnum = (short)atoi(tchar);
	for (i = 0; i < 12; i++) tchar[i] = line[i+26];
	tchar[12] = '\0';	ptp->p_xc = atof(tchar);
	for (i = 0; i < 8; i++) tchar[i] = line[i+38];
	tchar[8] = '\0';	ptp->p_yc = atof(tchar);
	for (i = 0; i < 8; i++) tchar[i] = line[i+46];
	tchar[8] = '\0';	ptp->p_zc = atof(tchar);
	for (i = 0; i < 6; i++) tchar[i] = line[i+54];
	tchar[6] = '\0';	ptp->p_occo = ptp->p_occ = atof(tchar);
	for (i = 0; i < 6; i++) tchar[i] = line[i+60];
	tchar[6] = '\0';	ptp->p_bvo = ptp->p_bval = atof(tchar);
	for (j = 0, i = 66; i < 80; i++)
	   if ((' ' != line[i]) && ('\n' != line[i]))
		ptp->p_atomid[j++] = line[i];
	ptp->p_atomid[j] = '\0';
}

static void findchain(char *li, TOTALST *top)
{ /* This looks to see whether chain tc has already been recorded;
   if it hasn't, it remembers it */
	int	i, nch;
	char	tc;
	CHAIN	*chp;

	nch = top->tnch;	/* number of chains already found */
	tc = *(li + 21);	/* this record's chain ID */
	if (nch > 4000) /* no chains found yet */
	  {
		chp = top->tchs;	nch = 0;
		chp->c_chainid = tc;
		chp->c_minwat = (short)atoi(li + 22);
		nch++;
	  }
	else /* see if this is a new chain */
	  {
		for (chp = top->tchs, i = 0; i < nch; i++, chp++)
			if (tc == chp->c_chainid) break;
		if (i >= nch) /* this is a new one: remember it */
		  {
			chp = top->tchs + nch;
			chp->c_chainid = tc;
			chp->c_minwat = (short)atoi(li + 22);
			nch++;
		  }
	  }
	top->tnch = nch;
}

static void getpdblin(char *li, TOTALST *top)
{ /* This reads a PDB line to see what it contains and aggregates
  counts of ATOM and HETATM records, plus the number of HOH entries */
	int	iswat;
	char	*aty;

	if ((li[0] != 'A') && (li[0] != 'H')) return; /* not ATOM or HETATM */
	aty = li + 17;
	iswat = (*aty == 'H') && (*(aty+1) == 'O') && (*(aty+2) == 'H');
	if (!strncmp(li, "ATOM  ", 6))
	  {
		top->tnatoms += 1;
		if (iswat) top->tnwatat += 1;
		else findchain(li, top);
	  }
	else if (!strncmp(li, "HETATM", 6))
	  {
		top->tnhet += 1;
		if (iswat) top->tnwathet += 1;
		else findchain(li, top);
	  }
}

static void procargs(int ac, char *av[], TOTALST *top)
{ /* argument processor for this program */
	char		*progname, *cp;
	double		thismaxhb, thisminhb, thishyd, thishet;

	top->tmeanbw = top->tnhet = top->tnwatat = top->tnwathet =
	 top->tnatoms = top->tnpat = top->tnpwat = top->tnwaters = 0;
	top->tfpl = NULL;
	top->tnch = 5000;
	if (NULL == (top->tchs = (CHAIN *)calloc(50, sizeof (CHAIN))))
	  {
		fprintf(stderr, "Cannot allocate memory for chain records\n");
		exit(1);
	  }
	top->t1ch_s = top->thighbq = 0;
	top->tmaxhbsq = MAXHBONDSQ;
	top->tminhbsq = MINHBONDSQ;
	top->tminhsq = CLOSEHSQ;
	top->tminhet = CLOSEHET;
	top->tbump = 0;
	thismaxhb = thisminhb = thishyd = thishet = 0.;
	progname = *av;	av++;	ac--;
	while ((ac > 0) && (**av == '-'))
	  {
		for (cp = (*av) + 1; (NULL != cp) && ('\0' != *cp); cp++)
		   {
			if (('a' <= *cp) && (*cp <= 'z')) *cp -= 'a' - 'A';
			if (*cp == 'S') top->t1ch_s = 1;
			else if (*cp == 'H') top->thighbq = 1;
			else if (*cp == 'X')
			  {
				thismaxhb = atof(cp + 1);	break;
			  }
			else if (*cp == 'M')
			  {
				thisminhb = atof(cp + 1);	break;
			  }
			else if (*cp == 'L')
			  {
				thishyd = atof(cp + 1);	break;
			  }
			else if (*cp == 'O')
			  {
				thishet = atof(cp + 1);	break;
			  }
			else if (*cp == 'B') top->tbump = 1;
		   }
		av++;	ac--;
	  }
	if ((0.05 < thismaxhb) && (thismaxhb < 10.)) top->tmaxhbsq = thismaxhb;
	if ((0.05 < thisminhb) && (thisminhb < 10.)) top->tminhbsq = thisminhb;
	if ((0.05 < thishyd) && (thishyd < 10.)) top->tminhsq = thishyd;
	if ((0.05 < thishet) && (thishet < 10.)) top->tminhet = thishet;
	top->tfpi = stdin;	top->tfpwat = stdout;
	top->tfpl = fopen("closewat.log", "w");
	if (ac)
	  {
		if (NULL == (top->tfpi = fopen(*av, "r")))
			top->tfpi = stdin;
		if (ac > 1)
		  {
			if (NULL == (top->tfpwat = fopen(*(av+1), "w")))
				top->tfpwat = stdout;
		  }
	  }
	fprintf(top->tfpl, "Run of %s\n", progname);
}

static void sum_avail(TOTALST *top, FILE *fpl)
{
	int	atnotwat, hetnotwat;

	atnotwat = top->tnatoms - top->tnwatat;
	hetnotwat = top->tnhet - top->tnwathet;
	fprintf(fpl,
 "  ATOM records    HETATM records     Total        # of\n");
	fprintf(fpl,
 " non-water water non-water water non-water water Chains\n");
	fprintf(fpl, " %9d %5d %9d %5d %9d %5d %6d\n",
	 atnotwat, top->tnwatat, hetnotwat, top->tnwathet,
	 atnotwat + hetnotwat, top->tnwatat + top->tnwathet, top->tnch);
}

static int ready(TOTALST *top)
{ /* cleans up some counts, rewinds the input, reserves memory */
	top->tnatoms += (top->tnhet - top->tnwathet);
	if (top->tnatoms < 1) /* give up if no atoms */
	  {
		fprintf(stderr, "No atom records detected in input\n");
		fprintf(stderr, "tnhet, tnwatat, tnwathet, tnch:\n");
		fprintf(stderr, "%5d, %7d, %8d, %4d\n",
		 top->tnhet, top->tnwatat, top->tnwathet, top->tnch);
		fprintf(stderr, "tnatoms, tnpat, tnpwat, tnwaters\n");
		fprintf(stderr, "%7d, %5d, %6d, %8d\n",
		 top->tnatoms, top->tnpat, top->tnpwat, top->tnwaters);
		return 1;
	  }
	if ((top->tnwaters = top->tnwatat + top->tnwathet) < 1)
	  { /* nothing to do if no waters */
		fprintf(stderr, "No water records detected in input\n");
		return 1;
	  }
	if (NULL == (top->tpat = (PDBRECORD *)calloc(top->tnatoms,
	 sizeof (PDBRECORD)))) /* reserve memory for PDB atom records */
	  {
		fprintf(stderr, "Error allocating memory for PDB records\n");
		return -1;
	  }
	if (NULL == (top->tpwa = (PDBRECORD *)calloc(top->tnwaters,
	 sizeof (PDBRECORD)))) /* reserve memory for waters */
	  {
		fprintf(stderr, "Error allocating memory for PDB waters\n");
		return -2;
	  }
	if (-1 == (int)fseek(top->tfpi, 0L, 0)) /* rewind the input file */
	  {
		fprintf(stderr, "Error rewinding input\n");
		return -3;
	  }
	sum_avail(top, stderr);
	if (NULL != top->tfpl) sum_avail(top, top->tfpl);
	top->tpatp = top->tpat;	top->tpwap = top->tpwa;
	return 0;
}

static void procpdblin(char *li, TOTALST *top)
{ /* populate the appropriate PDBRECORD structures */
	char		tc;
	int		ich, nch;
	PDBRECORD	*patp;
	CHAIN		*chp;

	if ((strncmp(li, "ATOM  ", 6)) && (strncmp(li, "HETATM", 6))) return;
	patp = top->tpatp;	nch = top->tnch;
	if (strncmp(&(li[17]), "HOH", 3)) /* not a water record */
	  {
		strtorec(li, patp); /* keep PDBRECORD information */
		tc = patp->p_chainid;
		for (chp = top->tchs, ich = 0; ich < nch; ich++, chp++)
		   {
			if (chp->c_chainid == tc)
			  {
				if (chp->c_minwat < patp->p_resnum)
				  {
					chp->c_minwat = patp->p_resnum;
				  }
				break;
			  }
		   }
		(top->tpatp)++;
	  }
	else { /* this is a water record: keep the PDB record */
		strtorec(li, top->tpwap);	(top->tpwap)++;
	     }
}

static int ismetal(char *atstr)
{ /* returns true if atstr is the name of a metal.
  For now we'll assume anything other than C,N,O,S,P, and H are metals.*/
	int		i;
	char		*cp, *cpo;
	static char	instr[80];
	static char	notmet[] = "CNOSPH";

	for (i = 0, cp = atstr, cpo = &(instr[0]); i < 4; i++, cp++)
	   {
		if ('\0' == *cp) break;
		if (' ' != *cp)
		  {	*cpo = *cp;
			if (('a' <= *cpo) && (*cpo <= 'z')) *cpo -= 'a' - 'A';
			cpo++;
		  }
	   }
	*cpo = '\0';
	if (strlen(instr) > 1) return 1;
	for (i = 0, cp = &(notmet[0]); i < 6; i++, cp++)
		if (instr[0] == notmet[i]) return 0;
	return 1;
}

static int seeneighbor(PDBRECORD *pwap, PDBRECORD *pwaq,
	double dsq, TOTALST *top)
{ /* returns a code value describing the relationship between
 the water record pwap and the neighbor record pwaq.
	Code values:
	0:	pwap and pwaq are distant from one another
	1:	pwap and pwaq are already independent of one another
		(e.g., the input conformer codes are different)
	2:	pwaq is a metal and therefore it's okay to be close
	3:	pwap and pwaq are independent in the input but
		 not in the proposed output
	4:	pwap and pwaq will become independent if
		 the conformer code for pwap is altered
	5:	pwap needs to be made a non-blank conformer
		 and its occupancy must be 0.5 or lower
	6:	pwap is just barely too close to pwaq and must be edited
	7:	pwap is considerably too close to pwaq
		 an must be manually edited.
	It also remembers that code in cases where we might want to use it
	for some other purpose */
	int	diag;

	if (ismetal(pwaq->p_atomid)) diag = 2;
	else if (pwaq->p_conf == ' ')
	  {
			if ((pwaq->p_attype[0] == 'H') ||
			 (pwaq->p_attype[1] == 'H'))
				diag = (dsq < top->tminhsq - HBARELY) ? 7 : 6;
			else
				diag = (dsq < top->tminhbsq - OBARELY)? 7 : 6;
	  }
	else if (pwap->p_conf == ' ')
	  {
		if ((pwap->p_confo != ' ') && (pwap->p_confo != pwaq->p_confo))
			diag = 3;
		else	diag = 5;
	  }
	else if (pwap->p_conf != pwaq->p_confo)
	  {
		if (pwap->p_confo == ' ') diag = 4;
		else if (pwap->p_confo != pwaq->p_conf) diag = 1;
		else diag = 4;
	  }
	else diag = 1;
	pwap->p_diag = diag;
	return diag;
}

static void closetitle(TOTALST *top)
{ /* prints a suitable header atop the list of waters near other atoms */
	fprintf(top->tfpl,
	 " Water molecules within %7.2f A of neighbors (or %7.2f A of H's)\n",
	 sqrt(top->tminhbsq), sqrt(top->tminhsq));
	fprintf(top->tfpl, "    Coordinates of Water Distance");
	fprintf(top->tfpl, " Renumbered   Original   Nearest Neighbor Atom\n");
	fprintf(top->tfpl, "       x       y       z               ");
	fprintf(top->tfpl,
	 "Water      Water  AtomID Renumbered   Original Diagnosis\n");
}

static void printclose(PDBRECORD *pwap, PDBRECORD *pwaq, double dsq, FILE *fp)
{
 
	fprintf(fp, " %7.3f %7.3f %7.3f %8.3f",
	 pwap->p_xc, pwap->p_yc, pwap->p_zc, sqrt(dsq));
	fprintf(fp, " %c%3s %c%4hd %c%3s %c%4hd ",
	 pwap->p_conf, pwap->p_resname, pwap->p_chainid, pwap->p_resnum,
	 pwap->p_confo, pwap->p_resname, pwap->p_chaino, pwap->p_resno);
	fprintf(fp, "   %4s %c%3s %c%4hd %c%3s %c%4hd",
	 pwaq->p_attype,
	 pwaq->p_conf, pwaq->p_resname, pwaq->p_chainid, pwaq->p_resnum,
	 pwaq->p_confo, pwaq->p_resname, pwaq->p_chaino, pwaq->p_resno);
	if ((0 < pwap->p_diag) && (pwap->p_diag <= 7))
		fprintf(fp, " %7s\n", diagstr[(int)(pwap->p_diag)]);
	else fprintf(fp, "\n");
}

static void swap01(PDBRECORD *pw0, PDBRECORD *pw1)
{ /* Swap the contents of records pw0 and pw1 */
	char		*cp0, *cp1;
	size_t		i, reclen;
	PDBRECORD	tem;

	reclen = sizeof(PDBRECORD);
	
	for (cp1 = (char *)&tem, cp0 = (char *)pw0, i = 0; i < reclen; i++)
	   {
		*cp1 = *cp0;	cp1++;	cp0++;
	   }
	for (cp1 = (char *)pw0,cp0 = (char *)pw1, i = 0; i < reclen; i++)
	   {
		*cp1 = *cp0;	cp1++;	cp0++;
	   }
	for (cp1 = (char *)pw1, cp0 = (char *)&tem, i = 0; i < reclen; i++)
	   {
		*cp1 = *cp0;	cp1++;	cp0++;
	   }
}

static char confchange(PDBRECORD *pwp, TOTALST *top)
{
	PDBRECORD 	*oth;

	if (pwp->p_conf == 'A')
	  {
		if (pwp >= top->tpwap - 1)
		  {
			pwp->p_conf = 'B';	pwp->p_diag = 4; return 'B';
		  }
		oth = pwp + 1;
		if ((oth->p_conf != 'B') || (oth->p_resnum != pwp->p_resnum))
		  {
			pwp->p_conf = 'B';	pwp->p_diag = 4;
		  }
		else
		  {
			swap01(pwp, oth);
			pwp->p_conf = 'A';	oth->p_conf = 'B';
			pwp->p_diag = oth->p_diag = 4;
		  }
		return pwp->p_conf;
	  }
	else if (pwp->p_conf == 'B')
	  {
		if (pwp <= top->tpwa)
		  {
			pwp->p_conf = 'A';	pwp->p_diag = 4;
			return pwp->p_conf;
		  }
		oth = pwp - 1;
		if ((oth->p_conf != 'A') || (oth->p_resnum != pwp->p_resnum))
		  {
			pwp->p_conf = 'A'; pwp->p_diag = 4; return pwp->p_conf;
		  }
		else
		  {
			swap01(pwp, oth);
			pwp->p_conf = 'B';	oth->p_conf = 'A';
			pwp->p_diag = oth->p_diag = 4;
		  }
		return pwp->p_conf;
	  }
	return '\0';
}

static void diagclose(PDBRECORD *pwap, PDBRECORD *pwaq,
 double dsq, TOTALST *top)
{ /* this diagnoses a close contact between water pwap and other record pwaq */
	int	diagnosis;
	char	avail;
	double	alpha, alphap;

	if ((pwaq->p_attype[0] == 'H') || (pwaq->p_attype[1] == 'H'))
	  { /* if it's a hydrogen, you need to be closer */
		if (dsq >= top->tminhsq) return;
	  }
	/* if it's already okay, don't count it as an error */
	if ((diagnosis = seeneighbor(pwap, pwaq, dsq, top)) < 2) return;
	if ((pwap->p_resnum == 796) || (pwap->p_resnum == 797))
		printclose(pwap, pwaq, dsq, stderr);
	/* if the user wants the software to do the bump, we do it here */
	if ((diagnosis == 6) && (top->tbump > 0))
	  {
		if ((pwaq->p_attype[0] == 'H') || (pwaq->p_attype[1] == 'H'))
			alpha = 1. / sqrt(dsq / top->tminhsq);
		else	alpha = 1. / sqrt(dsq / top->tminhbsq);
		/* we nudge alpha just slightly higher so that the
		ending distance between neighbors is a tiny bit farther */
		alpha *= 1.004;
		alphap = 1. - alpha;
		pwap->p_xc = pwap->p_xc * alpha + pwaq->p_xc * alphap;
		pwap->p_yc = pwap->p_yc * alpha + pwaq->p_yc * alphap;
		pwap->p_zc = pwap->p_zc * alpha + pwaq->p_zc * alphap;
		return; /* we're finished, in this case */
	  }
	/* resolve situations where pwaq has a nonblank conformer */
	if (pwaq->p_conf != ' ')
	  {
		if (pwap->p_conf == ' ')
		  { /* if pwap->p_conf is blank, make it either A or B */
			pwap->p_conf = (pwaq->p_conf == 'A') ? 'B' : 'A';
			if (pwap->p_occ > 0.5) pwap->p_occ = 0.5;
			return;
		  }
		/* if pwap is already nonblank and happens to be
		 a different conformer from pwaq, we're done */
		else if (pwap->p_conf != pwaq->p_conf) return;
		/* otherwise we have to get clever. I hope this works. */
		else if ('\0' != (avail = confchange(pwap, top)))
			pwap->p_conf = avail;
		else if ('\0' != (avail = confchange(pwaq, top)))
			pwaq->p_conf = avail;
		else
		  {
			fprintf(stderr,
		 "Warning: unresolved conformers for %c%4s %c%4hd\n",
	 	 pwap->p_conf, pwap->p_resname,
		 pwap->p_chainid, pwap->p_resnum);
			printclose(pwap, pwaq, dsq, stderr);
		  }
	  }
	if (top->tnclose == 0) closetitle(top);
	printclose(pwap, pwaq, dsq, top->tfpl);
	if ((0 < pwap->p_diag) && (pwap->p_diag <= 7))
		(top->talert[(int)pwap->p_diag])++;
	(top->tnclose)++;
	return;
}

static void relateem(TOTALST *top, PDBRECORD *pwap, PDBRECORD *thisp)
{
 /* This defines AHOH and BHOH records to replace the input records.
  If they're already AHOH and BHOH relative to one another, we skip them. */
	double		bscal, wt, occ0, occ1, occ00, occ01, b0, b00, b01;

	if ((NULL == pwap) || (NULL == thisp)) return;
	/* for now, we'll only merge waters that are
	 currently unmerged; but we'll remember the others. */
	if ((pwap->p_conf != ' ') || (thisp->p_conf != ' '))
	  {
		if ((pwap->p_chainid == thisp->p_chainid) &&
		 (pwap->p_resnum == thisp->p_resnum) &
		 (((pwap->p_conf == 'A') && (thisp->p_conf == 'B')) ||
		  ((pwap->p_conf == 'B') && (thisp->p_conf == 'A')))) return;
		if (pwap->p_conf != ' ') thisp->p_nbr = (long *)pwap;
		if (thisp->p_conf != ' ') pwap->p_nbr = (long *)thisp;
		return;
	  }
	occ00 = pwap->p_occ;	occ01 = thisp->p_occ;
	b00 = pwap->p_bval;	b01 = thisp->p_bval;
	occ0 = (occ00 / (occ00 + occ01)) * (b00 + b01) / (2 * b00);
	occ1 = (occ01 / (occ00 + occ01)) * (b00 + b01) / (2 * b01);
	if (occ0 > occ1)
	  {
		if (occ0 < 0.5) occ0 = 0.5;
		else if (occ0 > 0.75) occ0 = 0.75;
		occ1 = 1. - occ0;
	  }
	else {
		if (occ1 < 0.5) occ1 = 0.5;
		else if (occ1 > 0.75) occ1 = 0.75;
		occ0 = 1. - occ1;
	     }
	wt = fabs(b00 - b01) / (b00 + b01);
	if (wt < 0.07) bscal = 0.8;
	else bscal = (wt < 0.3) ? 0.9 : 0.96;
	b0 = bscal * ((b00 < b01) ? b00 : b01);
	pwap->p_nconfs = thisp->p_nconfs = 2;
	pwap->p_occ = occ0;	thisp->p_occ = occ1;
	if (b0 > 90.) fprintf(stderr,
 "Assigning B=%8.2f to %c%s %c%4hd @[%5.1f,%5.1f,%5.1f] in relateem\n",
	 b0, pwap->p_conf, pwap->p_resname,
	 pwap->p_chainid, pwap->p_resnum, pwap->p_xc, pwap->p_yc, pwap->p_zc);
	pwap->p_bval = thisp->p_bval = b0;
	if (occ0 >= occ1)
	  { /* pwap is more than 50% occupied: make it primary */
		thisp->p_resnum = pwap->p_resnum;
		thisp->p_chainid = pwap->p_chainid;
		pwap->p_conf = 'A';	thisp->p_conf = 'B';
	  }
	else
	  { /* thisp is more than 50% occupied: make it primary */
		pwap->p_resnum = thisp->p_resnum;
		pwap->p_chainid = thisp->p_chainid;
		pwap->p_conf = 'B';	thisp->p_conf = 'A';
	  }
}

static void sum_water(TOTALST *top,
 PDBRECORD *pwap, PDBRECORD *thisp, double mindist)
{
	static char	warmsg[] = "clash";
	static char	hbmsg[] = "Hbond";

	if (NULL == top->tfpl) return;
	fprintf(top->tfpl,
	 "%c%4d @%7.3f %7.3f %7.3f %s %5.2f A to",
	 pwap->p_chainid, pwap->p_resnum, pwap->p_xc, pwap->p_yc, pwap->p_zc,
	 ((mindist < top->tminhbsq) ? warmsg : hbmsg), sqrt(mindist));
	fprintf(top->tfpl,
	 " %c%4d @%7.3f %7.3f %7.3f\n",
	 thisp->p_chainid, thisp->p_resnum,
	 thisp->p_xc, thisp->p_yc, thisp->p_zc);
}

static int noconformers(TOTALST *top)
{ /* Since we're going to identify conformers ourselves,
   we remove the input ones. We rescale the occupancy and B values. */
	int		maxres, ncleared;
	double		corocc;
	PDBRECORD	*pwap;

	/* first determine the highest input residue number */
	for (maxres = -1000, pwap = top->tpwa; pwap < top->tpwap; pwap++)
		if (pwap->p_resnum > maxres) maxres = pwap->p_resnum;
	if (maxres < -999) return 0;
	for (ncleared = 0, pwap = top->tpwa; pwap < top->tpwap; pwap++)
		if (pwap->p_conf != ' ')
		  {
			if (pwap->p_conf != 'A')
			  {
				maxres++;	pwap->p_resnum = maxres;
			  }
			corocc = (pwap->p_occ < 0.81) ? 0.9 :
			 1. / sqrt(pwap->p_occ);
			pwap->p_bval *= corocc;
			if (pwap->p_bval > 80.) fprintf(stderr,
 "noconformers:B=%8.2f to %c%s %c%4hd @[%5.1f,%5.1f,%5.1f]: corocc=%7.4f\n",
			 pwap->p_bval, pwap->p_conf, pwap->p_resname,
			 pwap->p_chainid, pwap->p_resnum, pwap->p_xc,
			 pwap->p_yc, pwap->p_zc, corocc);
			if (pwap->p_occ < 0.5) pwap->p_occ *= 2.;
			else pwap->p_occ = 1.;
			pwap->p_conf = ' ';
			ncleared++;
		  }
	if (ncleared)
	  {
		fprintf(stderr,
	 " %6d waters had their conformation symbols cleared\n", ncleared);
		if (NULL != top->tfpl) fprintf(top->tfpl,
	 " %6d waters had their conformation symbols cleared\n", ncleared);
	  }
	return ncleared;
}

static void closestwat(TOTALST *top)
{ /* For each water, this finds the closest water to it.
	If it's closer than 2.5 A, we report it as a possible collision.
	If it's between 2.5 and 3.2 A, we report it as a hydrogen bond.
	In the current version we check all waters, not just the ones
	that occur later in the list than the current one. */

	PDBRECORD	*p0, *p1, *p2;
	double		dist, mindist;

	for (p0 = top->tpwa; p0 < top->tpwap - 1; p0++)
	   {
		mindist = 1.e10;
		for (p1 = top->tpwa; p1 < top->tpwap; p1++)
		   if (p1 != p0)
		     {
			if ((dist = pdbdist(p0, p1)) < mindist)
			  {
				mindist = dist;	p2 = p1;
			  }
		     }
		if (mindist < top->tmaxhbsq)
		  {
			sum_water(top, p0, p2, mindist);
			if (mindist < top->tminhbsq) relateem(top, p0, p2);
		  }
	   }
}

static void warn4(PDBRECORD *p0, PDBRECORD *p1, int whic)
{
	fprintf(stderr,
	 "Warn: %c%s %c%hd %7.2f %7.2f %7.2f %6.2f %7.2f near %s ",
		p0->p_conf, p0->p_resname, p0->p_chainid, p0->p_resnum,
		p0->p_xc, p0->p_yc, p0->p_zc, p0->p_occ, p0->p_bval,
		(whic == 1)? "thisp" : "other");
	fprintf(stderr, "%c%s %c%hd\n",
		p1->p_conf, p1->p_resname, p1->p_chainid, p1->p_resnum);
}

static void thisthird(TOTALST *top, PDBRECORD *p0)
{
	double		occsum;
	PDBRECORD	*p1, *p2, *p3, *wasa, *wasc, *small, *large;

	if (NULL == (p2 = (PDBRECORD *)p0->p_nbr)) return;
	if (p2->p_nconfs > 3) /* p2 already has 3 neighbors */
	  {	warn4(p0, p2, 0);
		return; /* already ABCD present */
	  }
	/* find the conformer associated with p2 */
	for (p1 = top->tpwa; p1 < top->tpwap; p1++)
		if (isconformer(top, p2, p1, 1) ||
			isconformer(top, p1, p2, 1)) break;
	if (p1 >= top->tpwap)
	  { /* failed to find the pair: give up */
		fprintf(stderr, "Couldn't find pair for %c%s %c%hd\n",
		 p2->p_conf, p2->p_resname,
		 p2->p_chainid, p2->p_resnum);
		return;
	  }
	if (p1->p_nconfs > 3)
	  { /* p1 aloready has 3 neighbors */
		warn4(p0, p1, 1);
		return;
	  }
	p2->p_nconfs = p1->p_nconfs = p0->p_nconfs = 3;
	p0->p_chainid = p2->p_chainid;
	p0->p_resnum = p2->p_resnum;
	if (p2->p_bval > 90.) fprintf(stderr,
 "Assigning B=%8.2f to %c%s %c%4hd @[%5.1f,%5.1f,%5.1f] in thisthird\n",
	 p2->p_bval, p0->p_conf, p0->p_resname,
	 p0->p_chainid, p0->p_resnum, p0->p_xc, p0->p_yc, p0->p_zc);
	p0->p_bval = p2->p_bval;
	if ((p2->p_conf == 'C') || (p1->p_conf == 'C'))
	  { /* if either p1 or p2 used to be C, then we now have 4 */
		if (p2->p_conf == 'C')
		  {
			wasc = p2;	wasa = p1;
		  }
		else {
			wasc = p1;	wasa = p2;
		     }
		for (p3 = top->tpwa; p3 < top->tpwap; p3++)
		   if ((p3 != p0) && (p3 != wasa) && (p3 != wasc) &&
		    isconformer(top, p3, wasc, 0) &&
		    isconformer(top, wasa, p3, 0))
		     {
			occsum = p1->p_occ + p2->p_occ + p3->p_occ;
			if (occsum > 0.001)
			  {
				p1->p_occ *= 0.75 / occsum;
				p2->p_occ *= 0.75 / occsum;
				p3->p_occ *= 0.75 / occsum;
			  }
			p0->p_chainid = p1->p_chainid;
			p0->p_resnum = p1->p_resnum;
			p0->p_occ = 0.25;
			if (p1->p_bval > 90.) fprintf(stderr,
 "Assigning B=%8.2f to %c%s %c%4hd @[%5.1f,%5.1f,%5.1f] in thisthird, b\n",
	 p1->p_bval, p0->p_conf, p0->p_resname,
	 p0->p_chainid, p0->p_resnum, p0->p_xc, p0->p_yc, p0->p_zc);
			p0->p_bval = p1->p_bval;
			p0->p_conf = 'D';
			p0->p_nconfs = p2->p_nconfs = p1->p_nconfs
			 = p3->p_nconfs = 4;
			return;
		     }
		return;
	  }
	p0->p_conf = 'C'; /* here we're just going from 2 to 3 */
	occsum = p1->p_occ + p2->p_occ;
	/* figure out which should be which */
	if (p0->p_bval < p1->p_bval)
	  { /* B value was lower than for the pair */
		if (p2->p_conf == 'B')
		  {
			p2->p_conf = 'C';
			p1->p_conf = 'B';
		  }
		else {
			p2->p_conf = 'B';
			p1->p_conf = 'C';
		     }
		p0->p_conf = 'A';
		p2->p_occ *= 0.64 / occsum;
		p1->p_occ *= 0.64 / occsum;
	  }
	else {
		p2->p_occ *= 0.68 / occsum;
		p1->p_occ *= 0.68 / occsum;
		if (p2->p_occ > p1->p_occ)
		  {
			small = p1;	large = p2;
		  }
		else {
			large = p1;	small = p2;
		     }
		large->p_conf = 'A';
		if (p0->p_occ > small->p_occ)
		  {
			p0->p_conf = 'B';
			small->p_conf = 'C';
		  }
		else {
			small->p_conf = 'B';
			p0->p_conf = 'C';
		     }
	     }
	/* occupancies of the three conformers should add to one */
	p0->p_occ = 1. - (p2->p_occ + p1->p_occ);
}

static void meanbw(TOTALST *top)
{ /* calculate the mean B value for all waters in this file */
	PDBRECORD	*p0;

	top->tmeanbw = 0.;
	for (p0 = top->tpwa; p0 < top->tpwap; p0++) top->tmeanbw += p0->p_bval;
	top->tmeanbw /= (double)(int)(top->tpwap - top->tpwa);
	fprintf(stderr, " Mean B value for waters: %8.3f\n", top->tmeanbw);
	if (NULL != top->tfpl) fprintf(top->tfpl,
		" Mean B value for waters: %8.3f\n", top->tmeanbw);
}

static void reduceocc(TOTALST *top)
{ /* reduce the occupancy of waters with big B values */
	int		i;
	double		occscale, rescale;
	PDBRECORD	*p0;

	for (i = 0, p0 = top->tpwa; p0 < top->tpwap; p0++)
	   if ((p0->p_occ > 0.99) && (p0->p_bval > 1.2 * top->tmeanbw))
	     {
		occscale = p0->p_bval / top->tmeanbw;
		rescale = pow(occscale, 0.33333);
		if (occscale >= 1.2)
		  {
			if (rescale > 1.25) rescale = 1.25;
			p0->p_bval /= rescale;
			p0->p_occ /= occscale;	i++;
		  }
	     }
	fprintf(stderr, " B and Q values adjusted for %6d %s\n",
	 i, (i == 1) ? "water" : "waters");
	if (NULL != top->tfpl) fprintf(top->tfpl,
	 " B and Q values adjusted for %6d %s\n",
	 i, (i == 1) ? "water" : "waters");
}

static void reportmult(TOTALST *top, PDBRECORD *p0, int reporttype)
{ /* Describes a duo, trio, or quartet of water molecules */
	PDBRECORD	*p1, *p2, *p3;
	char		*label;
	static char	labelquar[] = "Quad: ";
	static char	labeltrio[] = "Trio: ";
	static char	labelduo[] =  " Duo: ";

	if (NULL == top->tfpl) return;
	p1 = p0 + 1;	p2 = p1 + 1;	p3 = p2 + 1;
	if (reporttype == 4) label = labelquar;
	else if (reporttype == 3) label = labeltrio;
	else label = labelduo;
	fprintf(top->tfpl,
	 " %s %c%4hd: %8.3f %8.3f %8.3f %5.2f %6.2f (was %c%hd%c)\n",
	 label, p0->p_chainid, p0->p_resnum,
	 p0->p_xc, p0->p_yc, p0->p_zc, p0->p_occ, p0->p_bval,
	 p0->p_chaino, p0->p_resno, p0->p_confo);
	fprintf(top->tfpl,
	 "               %8.3f %8.3f %8.3f %5.2f %6.2f (was %c%hd%c)\n",
	 p1->p_xc, p1->p_yc, p1->p_zc, p1->p_occ, p1->p_bval,
	 p1->p_chaino, p1->p_resno, p1->p_confo);
	if (reporttype <= 2) return;
	fprintf(top->tfpl,
	 "               %8.3f %8.3f %8.3f %5.2f %6.2f (was %c%hd%c)\n",
	 p2->p_xc, p2->p_yc, p2->p_zc, p2->p_occ, p2->p_bval,
	 p2->p_chaino, p2->p_resno, p2->p_confo);
	if (reporttype <= 3) return;
	fprintf(top->tfpl,
	 "               %8.3f %8.3f %8.3f %5.2f %6.2f (was %c%hd%c)\n",
	 p3->p_xc, p3->p_yc, p3->p_zc, p3->p_occ, p3->p_bval,
	 p3->p_chaino, p3->p_resno, p3->p_confo);
}

static PDBRECORD *nearchain(PDBRECORD *pwap, TOTALST *top)
{ /* this finds the closest non-water chain to a particular water.
  It returns the next record we should be examining after this one. */
	int		ich;
	char		origid;
	unsigned short	orignum;
	double		dist, lx, ly, lz, mindist, rx, ry, rz;
	PDBRECORD	*pwe, *patp, *thisp;
	PDBRECORD	*next, *nextp1, *nextp2;
	CHAIN		*chp;

	/* the calculation will only be performed for the predominant
	 conformer on multi-conformer waters, which will be 'A' conformer. */
	next = NULL;	if (pwap < top->tpwap - 1) next = pwap + 1;
	if (pwap->p_conf == 'C') return next;
	if (pwap->p_conf == 'B')
	  {
		if (isconformer(top, pwap, next, 1)) return next + 1;
		else return next;
	  }
	lx = pwap->p_xc; ly = pwap->p_yc; lz = pwap->p_zc;
	thisp = NULL;	mindist = 1.e10;
	for (patp = top->tpat; patp < top->tpate; patp++)
	   { /* go through the non-water records */
		rx = patp->p_xc; ry = patp->p_yc; rz = patp->p_zc;
		/* calculate distance of this record from water */
		dist = (rx - lx) * (rx - lx) + (ry - ly) * (ry - ly) +
			(rz - lz) * (rz - lz);
		if (dist < mindist) /* if closest, remember it */
		  {
			mindist = dist;	thisp = patp;
		  }
	   }
	if (NULL == thisp) /* assign chain to this water */
	  {
		fprintf(stderr,
 "Warning: no chain assigned to H2O %c %d @ [%6.2f %6.2f %6.2f]\n",
		 pwap->p_chainid, pwap->p_resnum, lx, ly, lz);
		return pwap + 1;
	  }
	origid = pwap->p_chainid;	orignum = pwap->p_resnum;
	pwe = top->tpwap;	next = pwap + 1;
	for (ich = 0, chp = top->tchs; ich < top->tnch; ich++, chp++)
	   {
		if (chp->c_chainid == thisp->p_chainid)
		  {
			pwap->p_chainid = thisp->p_chainid;
			pwap->p_resnum = chp->c_curwat;
			if (pwap->p_conf == ' ')
			  {
				(chp->c_curwat)++;	return pwap+1;
			  }
			(chp->c_curwat)++;
			if ((pwap->p_conf == 'A') && (next < pwe) &&
			 (next->p_chainid == origid) &&
			 (next->p_resnum == orignum) &&
			 (next->p_conf == 'B'))
			  {
				nextp1 = next + 1;
				if (((nextp1) < pwe) &&
				 (nextp1->p_chainid == origid) &&
				 (nextp1->p_resnum == orignum) &&
				 (nextp1->p_conf == 'C'))
				  {
					next->p_chainid =
					 nextp1->p_chainid = pwap->p_chainid;
					next->p_resnum =
					 nextp1->p_resnum = pwap->p_resnum;
					nextp2 = nextp1 + 1;
					if ((nextp2 < pwe) &&
					 (nextp2->p_chainid == origid) &&
					 (nextp2->p_resnum == orignum) &&
					 (nextp2->p_conf == 'D'))
					  {
						nextp2->p_chainid =
							pwap->p_chainid;
						nextp2->p_resnum =
							pwap->p_resnum;
						reportmult(top, pwap, 4);
						return nextp2 + 1;
					  }
					else {
						reportmult(top, pwap, 3);
						return nextp1 + 1;
					     }
				  }
				else {
					next->p_chainid = pwap->p_chainid;
					next->p_resnum = pwap->p_resnum;
					reportmult(top, pwap, 2);
					return nextp1;
				     }
			  }
			else return next;
		  }
	   }
	if (ich >= top->tnch) fprintf(stderr,
 "Warning: %c %4d unassigned: near %c %4d %c%s @ [%7.2f %7.2f %7.2f]\n",
		pwap->p_chainid, pwap->p_resnum,
		thisp->p_chainid, thisp->p_resnum, thisp->p_conf,
		thisp->p_resname, thisp->p_xc, thisp->p_yc, thisp->p_zc);
	return next;
}

static void makechains(TOTALST *top)
{
	int		ich;
	CHAIN		*chp;
	PDBRECORD	*pwap;

	for (ich = 0, chp = top->tchs; ich < top->tnch; ich++, chp++)
	   { /* define minimum water residue number for each chain */
		chp->c_minwat += 100;
		chp->c_minwat = (short)(100 * (int)(chp->c_minwat / 100) + 1);
		chp->c_curwat = chp->c_minwat;
	   }
	/* for each water, find the chain closest to it */
	for (pwap = top->tpwa; pwap < top->tpwap; )
		pwap = nearchain(pwap, top);
}

static void adjustqb(TOTALST *top, PDBRECORD *p0)
{
 /* This tweaks the occupancy and Debye-Waller factors for a
  group of water conformers. In the process, we correct roundoff
  errors for the occupancy values */
	PDBRECORD	*p1, *p2, *p3;
	double		qp0, qp1, qp2, qp3;
	double		scq, sumq, bbar, sumqp, bmin, varn, dsumocc;
	int		occ1h0, occ1h1, occ1h2, occ1h3, sumocc;

	p1 = p0 + 1;
	if (p0->p_bvo < 2.) p0->p_bvo = 2.;
	if (p1->p_bvo < 2.) p1->p_bvo = 2.;
	bmin = p0->p_bvo;	if (p1->p_bvo < p0->p_bvo) bmin = p1->p_bvo;
	sumq = p0->p_occo + p1->p_occo;
	bbar = p0->p_occo * p0->p_bvo + p1->p_occo * p1->p_bvo;
	if (p0->p_nconfs > 2)
	  {
		p2 = p1 + 1;
		if (p2->p_bvo < 2.) p2->p_bvo = 2.;
		if (p2->p_bvo < bmin) bmin = p2->p_bvo;
		sumq += p2->p_occo;	bbar += p2->p_occo * p2->p_bvo;
		if (p2->p_bvo < bmin) bmin = p2->p_bvo;
		if (p0->p_nconfs > 3)
		  {
			p3 = p2 + 1;
			if (p3->p_bvo < 2.) p3->p_bvo = 2.;
			if (p3->p_bvo < bmin) bmin = p3->p_bvo;
			sumq += p3->p_occo; bbar += p3->p_occo * p3->p_bvo;
			if (p3->p_bvo < bmin) bmin = p3->p_bvo;
		  }
	  }
	if (sumq < 0.0001) return;
	bbar /= sumq;
	sumqp = qp0 = p0->p_occo * bbar / p0->p_bvo;
	qp1 = p1->p_occo * bbar / p1->p_bvo;	sumqp += qp1;
	varn = (p0->p_bvo - bbar) * (p0->p_bvo - bbar) * p0->p_occo +
		(p1->p_bvo - bbar) * (p1->p_bvo - bbar) * p1->p_occo;
	if (p0->p_nconfs > 2)
	  {
		qp2 = p2->p_occo * bbar / p2->p_bvo;	sumqp += qp2;
		varn += (p2->p_bvo - bbar) * (p2->p_bvo - bbar) * p2->p_occo;
		if (p0->p_nconfs > 3)
		  {
			qp3 = p3->p_occo * bbar / p3->p_bvo;	sumqp += qp3;
			varn += (p3->p_bvo - bbar) *
				(p3->p_bvo - bbar) * p3->p_occo;
		  }
		else qp3 = 0.;
	  }
	else qp2 = 0.;
	if (sumqp < 0.0001) return;
	varn = (varn > 0.0001) ? sqrt(varn / sumq) / bbar : 0.;
	 /* if original occupancies were low, don't scale B's down as much */
	if ((scq = sumq / (double)(p0->p_nconfs)) < 0.9)
	  {
		if (scq < 0.49) scq = 0.49;	bmin /= sqrt(scq);
	  }
	if (varn < 0.1) bmin = 0.8 * bmin;
	else bmin = (varn < 0.2) ? 0.9 * bmin : 0.95 * bmin;
	p0->p_occ = qp0 / sumqp;	p1->p_occ = qp1 / sumqp;
	if (bmin > 90.) fprintf(stderr,
 "Assigning B=%8.2f to %c%s %c%4hd @[%5.1f,%5.1f,%5.1f] in adjustqb\n",
	 bmin, p0->p_conf, p0->p_resname,
	 p0->p_chainid, p0->p_resnum, p0->p_xc, p0->p_yc, p0->p_zc);
	p1->p_bval = p0->p_bval = bmin;
	if (p0->p_nconfs > 2)
	  {
		p2->p_occ = qp2 / sumqp;
		p2->p_bval = bmin;
		if (p0->p_nconfs > 3)
		  {
			p3->p_occ = qp3 / sumqp; p3->p_bval = bmin;
		  }
	  }
	/* having adjusted the B and Q values, organize them by
	 decreasing value of Q */
	if (p0->p_nconfs == 2) /* if they're backwards, just reverse A for B */
	  {
		if (p0->p_occ < p1->p_occ)
		  {
			p0->p_conf = 'B';	p1->p_conf = 'A';
		  }
	  }
	else { /* for > 2 conformers, sort them by occupancy and redefine */
		qsort((void *)p0, (size_t)(p0->p_nconfs), sizeof (PDBRECORD),
	 	 (QPE)occsort); /* sort water records on occupancy */
		p0->p_conf = 'A';	p1->p_conf = 'B';
		p2->p_conf = 'C';
		if (p0->p_nconfs > 3) (p2+1)->p_conf = 'D';
	    }
	occ1h0 = p0->p_occ * 100. + 0.4999;
	occ1h1 = p1->p_occ * 100. + 0.4999;
	occ1h2 = occ1h3 = 0;
	if (p0->p_nconfs > 2)
	  {
		occ1h2 = p2->p_occ * 100. + 0.4999;
		if (p0->p_nconfs > 3)
			occ1h3 = p3->p_occ * 100. + 0.4999;
	  }
	sumocc = occ1h0 + occ1h1 + occ1h2 + occ1h3;
	dsumocc = ((double)sumocc) / 100.;
	if (dsumocc < 0.995)
	  {
		p0->p_occ = ((double)(occ1h0 + 1)) / 100.;
		dsumocc = p0->p_occ + p1->p_occ;
		if (p0->p_nconfs > 2)
		  {
			dsumocc += p2->p_occ;
			if (p0->p_nconfs > 3) dsumocc += (p2+1)->p_occ;
		  }
		if (dsumocc > 1.0) p0->p_occ -= (dsumocc - 0.9995);
	  }
	else if (dsumocc > 1.005)
		p1->p_occ = ((double)(occ1h1 - 1)) / 100.;
}

static void adjustmult(TOTALST *top)
{ /* Counts how many (AHOH,BHOH) pairs, (AHOH,BHOH,CHOH) triples,
	and (AHOH,BHOH,CHOH,DHOH) quads there are in the list.
	At the same time, adjust the B and Q values for the group */
	int		npairs, ntriples, nquads;
	PDBRECORD	*p0, *p1, *p2, *p3, *firstmult;

	npairs = ntriples = nquads = 0;
	/* We count the number of pairs, trios, and quads. Meanwhile,
	 because occupancy values are only kept to 2 significant figures,
	 we may occasionally get roundoff errors. Fix that. */
	firstmult = NULL;
	for (p0 = top->tpwa; p0 < top->tpwap - 1; p0++)
	   if (' ' != p0->p_conf)
	     {
		if (p0->p_conf == 'A')
		  {
			if (NULL == firstmult) firstmult = p0;
			p1 = p0 + 1;
			if (p1->p_conf == 'B')
			  {
				p2 = p1 + 1;
				if ((p2 < top->tpwap) && (p2->p_conf == 'C'))
				  {
					p3 = p2 + 1;
					if ((p3 < top->tpwap) &&
						(p3->p_conf == 'D'))
					  {
						nquads++;
						adjustqb(top, p0);
						p0 += 2; /* skip B, C, D */
					  }
					else {
						ntriples++;
						adjustqb(top, p0);
						p0 += 1; /* skip B, C */
					     }
				  }
				else {
					adjustqb(top, p0);
					npairs++;
				     }
			  }
		  }
	     }
	fprintf(stderr,
	 "Multiple-conformation waters: %5d pairs, %5d triples, %5d quads\n",
		npairs, ntriples, nquads);
	if (NULL != top->tfpl) fprintf(top->tfpl,
	 "Multiple-conformation waters: %5d pairs, %5d triples, %5d quads\n",
		npairs, ntriples, nquads);
}

static void finalmult(TOTALST *top)
{ /* Counts how many (AHOH,BHOH) pairs, (AHOH,BHOH,CHOH) triples,
	and (AHOH,BHOH,CHOH,DHOH) quads there are in the list.
	This time we don't fiddle with them: we leave them as they are. */
	int		npairs, ntriples, nquads, nonconf;
	int		lonea, loneb, lonec, loned, totwat;
	short		p0num, p1num, p2num;
	PDBRECORD	*p0, *p1, *p2, *p3;

	npairs = ntriples = nquads = nonconf = 0;
	lonea = loneb = lonec = loned = 0;
	/* We count the number of pairs, trios, and quads */
	for (p0 = top->tpwa; p0 < top->tpwap - 1; p0++)
	 {
	   if (' ' != p0->p_conf)
	     {
		p0num = p0->p_resnum;
		if (p0->p_conf == 'A')
		  {
			p1 = p0 + 1;	p1num = p1->p_resnum;
			if ((p1 < top->tpwap) && (p1->p_conf == 'B') &&
			 (p0num == p1num))
			  {
				p2 = p1 + 1;	p2num = p2->p_resnum;
				if ((p2 < top->tpwap) && (p2->p_conf == 'C')
				 && (p2num == p0num))
				  {
					p3 = p2 + 1;
					if ((p3 < top->tpwap) &&
						(p3->p_conf == 'D'))
					  {
						nquads++;
						p0 += 2; /* skip B, C, D */
					  }
					else {
						ntriples++;
						p0 += 1; /* skip B, C */
					     }
				  }
				else npairs++;
				p0++;
			  }
			else lonea++;
		  }
		else 
		  {
			if (p0->p_conf == 'B') loneb++;
			else if (p0->p_conf == 'C') lonec++;
			else if (p0->p_conf == 'D') loned++;
		  }
	     }
	    else nonconf++;
	 }
	totwat = top->tpwap - top->tpwa;
	fprintf(stderr, "%6d total waters: %5d unconformed,", totwat, nonconf);
	fprintf(stderr, "lone conformers:%3d A,%3d B,%3d C,%3d D\n",
		lonea, loneb, lonec, loned);
	fprintf(stderr, "%5d pairs, %5d triples, %5d quads\n",
		npairs, ntriples, nquads);
	if (NULL == top->tfpl) return;
	fprintf(top->tfpl, "%6d total waters: %5d unconformed, ",
		totwat, nonconf);
	fprintf(top->tfpl, "lone conformers:%3d A,%3d B,%3d C,%3d D\n",
		lonea, loneb, lonec, loned);
	fprintf(top->tfpl, "%5d pairs, %5d triples, %5d quads\n",
		npairs, ntriples, nquads);
}

static void cpypdb(PDBRECORD *pin, PDBRECORD *pou)
{ /* This copies the PDBRECORD pin into the space for pou */
	char	*cpo, *cpi;
	size_t	i, sizped;

	sizped = sizeof (PDBRECORD);
	for (cpo = (char *)pou, cpi = (char *)pin, i = 0;
	 i < sizped; i++, cpo++, cpi++) *cpo = *cpi;
}

static void incr_reswats(PDBRECORD *p3, int resnoff, TOTALST *top)
{
	int		i;
	CHAIN		*chp;
	PDBRECORD	*pwaq;

	/* find the last water in this chain
	 and remember its residue number */
	if (p3->p_chainid == 'S') chp = top->tchs;
	else {
		for (i = 0, chp = top->tchs; i < top->tnch; i++, chp++)
		   if (chp->c_chainid == p3->p_chainid) break;
		if (i >= top->tnch) return;
	     }
	/* increment residue numbers for remaining waters in this chain*/
	for (pwaq = p3 + 1; (NULL != pwaq) && (pwaq < top->tpwap) &&
	 (pwaq->p_chainid == p3->p_chainid); pwaq++) pwaq->p_resnum += resnoff;
}

static void oneplus3(PDBRECORD *p0, PDBRECORD *p1, PDBRECORD *p2,
 PDBRECORD *p3, TOTALST *top)
{ /* handles a case where the 1st water is alone and the other is a triplet */
	double	dsumocc;

	p3->p_resnum = p2->p_resnum = p1->p_resnum = p0->p_resnum + 1;
	p0->p_conf = 'A';	p1->p_conf = 'C';
	p2->p_conf = 'A';	p3->p_conf = 'B';
	swap01(p1, p3);		swap01(p1, p2);
	p0->p_occ = 0.5;	p0->p_bval *= 1.2;
	dsumocc = p1->p_occ + p2->p_occ + p3->p_occ;
	if (dsumocc < 0.01) dsumocc = 1.;
	else if (dsumocc > 1.) dsumocc = 1.;
	p1->p_occ /= dsumocc; p2->p_occ /= dsumocc; p3->p_occ /= dsumocc;
	p1->p_bval *= 1.05; p2->p_bval *= 1.05; p3->p_bval *= 1.05;
	incr_reswats(p3, 1, top);	return;
}

static void reorg4(PDBRECORD *p0, PDBRECORD *p1,
 PDBRECORD *p2, PDBRECORD *p3, TOTALST *top)
{ /* This reorganizes four consecutive PDBRECORDs that point to
	conformers of the same water into some other arrangement,
	if it can get away with doing so. */
	int		i, ngtcut, isgtcut, minpair, resnoff;
	int		isgtcutp, ngtcutp;
	double		occsum, minpairdsq, dsq[6], aved12, aved34;

	dsq[0] = pdbdist(p0, p1);	dsq[1] = pdbdist(p0, p2);
	dsq[2] = pdbdist(p0, p3);	dsq[3] = pdbdist(p1, p2);
	dsq[4] = pdbdist(p1, p3);	dsq[5] = pdbdist(p2, p3);
	for (isgtcut = ngtcut = i = 0; i < 6; i++)
	   {
		if (dsq[i] > MINHBONDSQ)
		  {
			isgtcut |= (1 << i);
			ngtcut++;
		  }
	   }
	fprintf(stderr, "Split 4-conformer water %c%4d into ",
	 p0->p_chainid, p0->p_resnum);
	if (ngtcut < 1)
	  {
		fprintf(stderr,
		 "No squared distances less than %7.3f\n",
		  MINHBONDSQ);
		return;
	  }
	if (ngtcut == 1)
	  {
		fprintf(stderr, "A + ABC\n");
		if (isgtcut == 02) swap01(p1, p2);
		else if (isgtcut == 04) swap01(p1, p3);
		else if (isgtcut == 010) swap01(p0, p2);
		else if (isgtcut == 020) swap01(p0, p3);
		else if (isgtcut == 040)
		  {
			swap01(p0, p2);	swap01(p1, p3);
		  }
		oneplus3(p0, p1, p2, p3, top);	return;
	  }
	else if ((ngtcut == 2) &&
	 ((isgtcut != 014) && (isgtcut != 022) && (isgtcut != 041)))
	  {
		fprintf(stderr, "ABC + A\n");
		if (isgtcut == 05) swap01(p1, p2);
		else if (isgtcut == 03) swap01(p1, p3);
		else if (isgtcut == 030) swap01(p0, p1);
		else if (isgtcut == 021)
		  {
			swap01(p0, p1);	swap01(p1, p2);
		  }
		else if (isgtcut == 011)
		  {
			swap01(p0, p1); swap01(p1, p3);
		  }
		else if (isgtcut == 050)
		  {
			swap01(p0, p2); swap01(p1, p2);
		  }
		else if (isgtcut == 042) swap01(p0, p2);
		else if (isgtcut == 012)
		  {
			swap01(p0, p2); swap01(p1, p3);
		  }
		else if (isgtcut == 060)
		  {
			swap01(p0, p3); swap01(p1, p3);
		  }
		else if (isgtcut == 044) swap01(p0, p3);
		else if (isgtcut == 024)
		  {
			swap01(p0, p3); swap01(p1, p2);
		  }
		else if (isgtcut != 06)
		  {
			fprintf(stderr, "Failure\n");
			return;
		  }
		oneplus3(p0, p1, p2, p3, top);	return;
	  }
	if ((022 == (isgtcut & 022)) || (041 == (isgtcut & 041)) ||
	 (014 == (isgtcut & 014)))
	  {
		fprintf(stderr, "2 AB pairs\n");
		if (041 == (isgtcut & 041)) swap01(p1, p2);
		else if (014 == (isgtcut & 014)) swap01(p3, p2);
		dsq[0] = pdbdist(p0, p1);	dsq[1] = pdbdist(p0, p2);
		dsq[2] = pdbdist(p0, p3);	dsq[3] = pdbdist(p1, p2);
		dsq[4] = pdbdist(p1, p3);	dsq[5] = pdbdist(p2, p3);
		for (isgtcutp = ngtcutp = i = 0; i < 6; i++)
		   {
			if (dsq[i] > MINHBONDSQ)
			  {
				isgtcutp |= (1 << i);	ngtcutp++;
			  }
		   }
		if (022 != (isgtcutp & 022))
		  {
			if (041 == (isgtcut & 041)) swap01(p1, p2);
			else if (014 == (isgtcut & 014)) swap01(p3, p2);
			dsq[0] = pdbdist(p0, p1); dsq[1] = pdbdist(p0, p2);
			dsq[2] = pdbdist(p0, p3); dsq[3] = pdbdist(p1, p2);
			dsq[4] = pdbdist(p1, p3); dsq[5] = pdbdist(p2, p3);
			for (isgtcutp = ngtcutp = i = 0; i < 6; i++)
			   {
				if (dsq[i] > MINHBONDSQ)
				  {
					isgtcutp |= (1 << i);
					ngtcutp++;
				  }
			   }
			if (022 != (isgtcutp & 022))
			  {
				fprintf(stderr, "Unable to complete split\n");
				return;
			  }
		  }
		/* p0&p1 are 1 residue, p2&p3 are another */
		p0->p_conf = p2->p_conf = 'A';
		p1->p_conf = p3->p_conf = 'B';
		occsum = p2->p_occ + p3->p_occ;	if (occsum > 1.) occsum = 1.;
		if (occsum < 0.01) p2->p_occ = p3->p_occ = 0.5;
		else {
			p2->p_occ /= occsum;	p3->p_occ /= occsum;
		     }
		p2->p_bval *= 1.1;	p3->p_bval *= 1.1;
		p3->p_resnum = p2->p_resnum = p0->p_resnum + 1;
		resnoff = 1;
		/* adjust the occupancies and B values for p0 and p1 */
		occsum = p0->p_occ + p1->p_occ;	if (occsum > 1.) occsum = 1.;
		if (occsum < 0.01) p0->p_occ = p1->p_occ = 0.5;
		else {
			p0->p_occ /= occsum;	p1->p_occ /= occsum;
		     }
		p0->p_bval *= 1.1;	p1->p_bval *= 1.1;
		incr_reswats(p3, resnoff, top);	return;
	  }
	if ((isgtcut == 07) || (isgtcut == 021) || (isgtcut == 052) ||
	 (isgtcut == 064))
	  {
		fprintf(stderr, "Blank + ABC\n");
		if (isgtcut == 021) swap01(p0, p1);
		if (isgtcut == 052) swap01(p0, p2);
		if (isgtcut == 064) swap01(p0, p3);
		p0->p_conf = ' ';
		p1->p_conf = 'A'; p2->p_conf = 'B'; p3->p_conf = 'C';
		p3->p_resnum = p2->p_resnum = p1->p_resnum = p0->p_resnum + 1;
		p0->p_occ = 0.5;
		p3->p_occ = p2->p_occ = p1->p_occ = 0.3333;
		p0->p_bval *= 1.1;
		resnoff = 1;
		incr_reswats(p3, resnoff, top);	return;
	  }
	else if ((isgtcut == 013) || (isgtcut == 025) || (isgtcut == 046)
	 || (isgtcut == 070))
	    {
		fprintf(stderr, "one pair+2 singles\n");
		for  (minpair = -1, i = 0, minpairdsq = 10000.; i < 6; i++)
		   {
			if (dsq[i] < minpairdsq)
			  {
				minpairdsq = dsq[i];	minpair = i;
			  }
		   }
		if (minpair < 0)
		  {
			fprintf(stderr, "Nonsensical pairing\n");
			return;
		  }
		if (minpair == 1) swap01(p1, p2);
		else if (minpair == 2) swap01(p1, p3);
		else if (minpair == 3) swap01(p0, p2);
		else if (minpair == 4) swap01(p0, p3);
		else if (minpair == 5)
		  {
			swap01(p0, p2);	swap01(p3, p1);
		  }
		dsq[0] = pdbdist(p0, p1);	dsq[1] = pdbdist(p0, p2);
		dsq[2] = pdbdist(p0, p3);	dsq[3] = pdbdist(p1, p2);
		dsq[4] = pdbdist(p1, p3);	dsq[5] = pdbdist(p2, p3);
		/* so now the minimum pair is [0,1]. */
		aved12 = 0.5 * (dsq[1] + dsq[2]);
		aved34 = 0.5 * (dsq[3] + dsq[4]);
		if (aved12 < aved34)
		  {
			swap01(p0, p1);
			dsq[0] = pdbdist(p0, p1); dsq[1] = pdbdist(p0, p2);
			dsq[2] = pdbdist(p0, p3); dsq[3] = pdbdist(p1, p2);
			dsq[4] = pdbdist(p1, p3); dsq[5] = pdbdist(p2, p3);
		  }
		p0->p_conf = 'A';	p1->p_conf = 'B';
		p2->p_conf = 'A';	p3->p_conf = 'A';
		p1->p_resnum = p0->p_resnum;
		p2->p_resnum = p1->p_resnum + 1;
		p3->p_resnum = p2->p_resnum + 1;
		p1->p_occ = p0->p_occ = 0.5;
		p3->p_occ = p2->p_occ = 0.5;
		p0->p_bval *= 1.1;
		p3->p_occ = p2->p_occ = 0.5;
		p2->p_bval *= 1.2;	p3->p_bval *= 1.2;
		p2->p_resnum = p0->p_resnum + 1;
		p3->p_resnum = p2->p_resnum + 1;
		resnoff = 2;
		incr_reswats(p3, resnoff, top);	return;
	  }
	else {
		fprintf(stderr,
 "No separable pairs exist: isgtcut = 0%o\n", isgtcut);
		fprintf(stderr, "dsq = %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f\n",
		 dsq[0], dsq[1], dsq[2], dsq[3], dsq[4], dsq[5]);
		return;
	     }
}

static void split4(PDBRECORD *pwap, TOTALST *top)
{ /* This splits a group of four waters into two groups of 2 if possible */
	PDBRECORD	*pm1, *pm2, *pm3;

	pm1 = pwap - 1; pm2 = pm1 - 1; pm3 = pm2 - 1;
	if (pm3 < top->tpwa) return;
	if ((pm3->p_conf != 'A') || (pm2->p_conf != 'B') || 
	 (pm1->p_conf != 'C'))
	  {
		fprintf(stderr, "%c%4s %c%4hd isn't preceded by A,B,C\n",
		 pwap->p_conf, pwap->p_resname,
		 pwap->p_chainid, pwap->p_resnum);
		fprintf(stderr, "%c%4s %c%4hd, %c%4s %c%4hd, %c%4s %c%4hd\n",
		 pm3->p_conf, pm3->p_resname, pm3->p_chainid, pm3->p_resnum,
		 pm2->p_conf, pm2->p_resname, pm2->p_chainid, pm2->p_resnum,
		 pm1->p_conf, pm1->p_resname, pm1->p_chainid, pm1->p_resnum);
		return;
	  }
	if (strcmp(pm3->p_resname, "HOH") ||
	 strcmp(pm2->p_resname, "HOH") || strcmp(pm1->p_resname, "HOH"))
	  {
		fprintf(stderr,
		 "Residues preceding %c%4s aren't all water\n",
			pwap->p_conf, pwap->p_resname);
		return;
	  }
	if ((pm3->p_chainid != pwap->p_chainid) ||
	 (pm2->p_chainid != pwap->p_chainid) ||
	 (pm1->p_chainid != pwap->p_chainid))
	  {
		fprintf(stderr,
		 "Residues preceding %c%4s have different chain IDs\n",
			pwap->p_conf, pwap->p_resname);
		return;
	  }
	if ((pm3->p_resnum != pwap->p_resnum) ||
	 (pm2->p_resnum != pwap->p_resnum) || (pm1->p_resnum != pwap->p_resnum))
	  {
		fprintf(stderr,
	 "Alert: changing residue number for %c%4s %c%4hd to %4hd\n",
		 pwap->p_conf, pwap->p_resname, pwap->p_chainid,
		 pwap->p_resnum, pm1->p_resnum);
		pwap->p_resnum = pm1->p_resnum;
	  }
	/* Phew. That was boring. Now try to split it */
	reorg4(pm3, pm2, pm1, pwap, top);
	return;
}

static void proximity(TOTALST *top)
{ /* for each water we identify close contacts */
	int		i, i0;
	PDBRECORD	*pwap, *pat, *pclose;
	double		dsq, closestsq;

	fprintf(top->tfpl, " %6d waters, %6d non-waters\n",
	 (int)(top->tpwap - top->tpwa), (int)(top->tpatp - top->tpat));
	fprintf(stderr, " %6d waters, %6d non-waters\n",
	 (int)(top->tpwap - top->tpwa), (int)(top->tpatp - top->tpat));
	/* look for waters that have four conformations (A,B,C,D).
	 In some cases we can create two conformations of two each */
	for (pwap = top->tpwa; pwap < top->tpwap; pwap++)
	   if (pwap->p_conf == 'D') split4(pwap, top);
	for (top->tnclose = 0, pwap = top->tpwa; pwap < top->tpwap; pwap++)
	   {
		for (pat = top->tpat; pat < top->tpatp; pat++)
		   if ((dsq = pdbdist(pwap, pat)) < top->tminhbsq)
			diagclose(pwap, pat, dsq, top);
		/* for the water proximities, we only look for waters earlier
		 in the list rather than later so we don't double-count */
		for (pat = top->tpwa; pat < pwap; pat++)
		   if ((pwap->p_resnum != pat->p_resnum) &&
		     ((dsq = pdbdist(pwap, pat)) < top->tminhbsq))
			diagclose(pwap, pat, dsq, top);
	   }
	fprintf(top->tfpl, " %6d close %s noted", top->tnclose,
	 (top->tnclose == 1) ? "contact": "contacts");
	fprintf(stderr, " %6d close %s noted", top->tnclose,
	 (top->tnclose == 1) ? "contact": "contacts");
	for (i0 = -1, i = 0; i < 8; i++)
	   if (top->talert[i] > 0)
	     {
		if (i0 == -1) i0 = i;
		fprintf(top->tfpl, "%c %3hd %s",
			(i0 == i) ? ':' : ',', top->talert[i], diagstr[i]);
		fprintf(stderr, "%c %3hd %s",
			(i0 == i) ? ':' : ',', top->talert[i], diagstr[i]);
	     }
	fprintf(top->tfpl, "\n");
	fprintf(stderr, "\n");
	for (top->tnclose = 0, pwap = top->tpwa; pwap < top->tpwap; pwap++)
	   {
		closestsq = 1.e9;	pclose = (PDBRECORD *)NULL;
		for (pat = top->tpat; pat < top->tpatp; pat++)
		   if ((pat->p_attype[0] != 'C') && (pat->p_attype[0] != 'H') &&
		    (pat->p_attype[1] != 'C') && (pat->p_attype[1] != 'H'))
		     {
			dsq = pdbdist(pwap, pat);
			if (dsq < closestsq)
			  {
				closestsq = dsq;	pclose = pat;
			  }
		     }
		for (pat = top->tpwa; pat < top->tpwap; pat++)
		   if (pat != pwap)
		     {
			dsq = pdbdist(pwap, pat);
			if (dsq < closestsq)
			  {
				closestsq = dsq;	pclose = pat;
			  }
		     }
		if ((NULL != pclose) && (closestsq > top->tminhet))
		  {
			if (!(top->tnclose)) fprintf(top->tfpl,
 "Waters that are more than %8.2f A from the nearest heteroatom neighbor:\n",
			 sqrt(top->tminhet));
			pwap->p_diag = 0;
			printclose(pwap, pclose, closestsq, top->tfpl);
			(top->tnclose)++;
		  }
	   }
	fprintf(top->tfpl, " %6d %s too far from all neighbors\n",
	 top->tnclose, (top->tnclose == 1) ? "water is" : "waters are");
	fprintf(stderr, " %6d %s too far from all neighbors\n",
	 top->tnclose, (top->tnclose == 1) ? "water is" : "waters are");
}

static void sortmults(TOTALST *top)
{
	char		tc;
	short		minrnum, trn, prevn;
	int		ich;
	PDBRECORD	*p0, *p00, *pne, *pq;
	CHAIN		*tch;

	tc = ' ';	pne = NULL;	p0 = top->tpwa;
	while (1)
	     {
		p00 = p0;
		for ( ; p0 < top->tpwap; p0++)
		   if (p0->p_conf != ' ')
		     {	/* skip past non-conformers at beginning */
			pne = p0;	tc = p0->p_chainid;	break;
		     }
		tch = (CHAIN *)NULL;
		for (tch = top->tchs, ich = 0; ich < top->tnch; ich++, tch++)
			if (tch->c_chainid == tc) break;
		if ((p0 >= top->tpwap) || (' ' == tc)) return;
		if (NULL != tch)
		  {
			tch->c_wat0 = p00 - top->tpwa; /* 0th water in chain */
			tch->c_watl1 = p0 - top->tpwa; /* final nonconf H2O */
			tch->c_watm0 = tch->c_watl1 + 1; /* 1st H2O with conf*/
		  }
		/* work past all the multiple conformers in this chain */
		for ( ; (p0 < top->tpwap) && (p0->p_chainid == tc) &&
		 (p0->p_conf != ' '); p0++) ;
		if (p0 > pne)
		  {
			qsort(pne, (size_t)(p0 - pne), sizeof (PDBRECORD),
				(QPE)bonly);
			/* having sorted them, we need to renumber them */
			minrnum = 10000;
			for (pq = pne; pq < p0; pq++)
			  if (pq->p_resnum < minrnum) minrnum = pq->p_resnum;
			trn = minrnum;	prevn = pne->p_resnum;
			for (pq = pne; pq < p0; pq++)
			   {
				if (pq->p_resnum != prevn)
				  {
					prevn = pq->p_resnum;	trn++;
				  }
				pq->p_resnum = trn;
			   }
			/* add NUMGAP to each of these so there's a gap
			 between the last single-conformation water
			 and the first multiple-conformation water */
			for (pq = pne; pq < p0; pq++) pq->p_resnum += NUMGAP;
			/* remember the final water for this chain */
			if (NULL != tch) tch->c_watl = p0 - top->tpwa;
		   }
		else  fprintf(stderr,
		 "Warning: Water %d, %d appear out of order\n",
		 (int)(p0 - top->tpwa), (int)(pne - top->tpwa));
	     }
}

static void insert_singles(TOTALST *top)
{
	int		i;
	CHAIN		*tch;
	PDBRECORD	*pwa0, *pwap, *pwaq, *pwae, *pwlas0, tempwa;

	pwa0 = top->tpwa;
	for (tch = top->tchs, i = 0; i  < top->tnch; i++, tch++)
	   {
		pwap = pwa0 + (int)(tch->c_wat0);
		pwlas0 = pwa0 + (int)(tch->c_watl1);
		pwae = pwlas0 - 1;
		for ( ; pwap < pwlas0; )
		   {
			if (pwap->p_conf != ' ')
			  {
				cpypdb(pwap, &tempwa);
				for (pwaq = pwap; pwaq < pwae; pwaq++)
					cpypdb(pwaq + 1, pwaq);
				cpypdb(&tempwa, pwaq);
				pwlas0--;	(tch->c_watl1)--;
			  }
			else pwap++;
		   }
		/* renumber the residues in this chain */
		pwaq = pwa0 + (int)(tch->c_wat0);
		for (pwap = pwaq; pwap < pwa0 + tch->c_watm0 - 1; pwap++)
			pwap->p_resnum = tch->c_minwat + (pwap - pwaq);
	   }
}

int main(int argc, char *argv[])
{
	int		i;
	char		lin[180];
	PDBRECORD	*pwap;
	TOTALST		tos;

	procargs(argc, argv, &tos); /* get arguments, reserve memory */
	while (NULL != fgets(lin, 120, tos.tfpi))
		getpdblin(lin, &tos); /* read */
	if (0 != (i = ready(&tos))) /* clean up counts and rewind input */
	  {
		free((void *)(tos.tchs));
		exit(i);
	  }
	while (NULL != fgets(lin, 80, tos.tfpi)) procpdblin(lin, &tos);
	tos.tnwaters = tos.tpwap - tos.tpwa;
	(void)noconformers(&tos); /* get rid of input AHOH,BHOH,etc. */
	tos.tpate = tos.tpatp;
	closestwat(&tos); /* for each water:find nearest-neighbor water*/
	for (pwap = tos.tpwa; pwap < tos.tpwap; pwap++)
	   if ((NULL != pwap->p_nbr) && (pwap->p_conf == ' '))
	     { /* this water has a close neighbor but isn't labeled yet */
		thisthird(&tos, pwap); /* analyze neighbors of pwap */
	     }
	meanbw(&tos);
	/* if user requested it, reduce occupancy of waters with B > 1.2*<B> */
	if (tos.thighbq) reduceocc(&tos);
	qsort((void *)tos.tpwa, (size_t)(tos.tnwaters), sizeof (PDBRECORD),
	 (QPE)occbsort); /* sort water records on occupancy and B value */
	makechains(&tos); /* figure out which chain each water belongs to */
	adjustmult(&tos); /* count and modify multiple conformers */
	qsort((void *)tos.tpwa, (size_t)(tos.tnwaters), sizeof (PDBRECORD),
		(QPE)chrsort); /* sort waters on chain # and residue # */
	/* Reorder multiple-conformer waters within each chain */
	sortmults(&tos);
	for (i = 0; i < 8; i++) tos.talert[i] = 0;
	/* if there's only one chain and the user has specified the 'S' flag,
	 then just write to a chain called S */
	if ((1 == tos.tnch) && (tos.t1ch_s) && ('S' != (tos.tchs)->c_chainid))
	  {
		(tos.tchs)->c_curwat = (tos.tchs)->c_minwat = 1;
		i = (tos.tpwa)->p_resnum - 1;
		for (pwap = tos.tpwa; pwap < tos.tpwap; pwap++)
		   {
			pwap->p_chainid = 'S';	pwap->p_resnum -= i;
		   }
	  }
	proximity(&tos); /* for all waters, identify close contacts */
	/* move the single-conformation but nonetheless conformation-marked
	 waters into the empty zone between the unconformation-marked
	 waters and the multiple-conformation waters */
	insert_singles(&tos);
	finalmult(&tos); /* final count of multiple conformers */
	for (pwap = tos.tpwa; pwap < tos.tpwap; pwap++)
		outrec(pwap, tos.tfpwat); /* write out the waters */
	free((void *)(tos.tpat));
	free((void *)(tos.tpwa));
	free((void *)(tos.tchs));
	exit(0);
}

/* end of closewat.c */
