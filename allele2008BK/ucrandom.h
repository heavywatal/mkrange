


short randc(void);  /* random integer 0〜32767 */
void srand( unsigned short seed);
unsigned long lsrandc(unsigned long  seed); 
unsigned long randomizec(void);/* initialize randc() */
short rnd( short x); /*  random integer 0 <= rnd <= x  */
double drnd( double x); /*  random double 0 <= drnd <= x  */
short zrnd( short x); /* random double -x <= zrnd <= x  */
double nrnd( double x); /*   Normal RaNDom */
short prnd( double n);
double urnd(void);/* 0≦urnd()<1 */
double normalselect (double mean, double sd);
short randombit();
void randomove (long x, long y, long d,  long *ntx, long *nty);
void ArrangeFreq (short n, double freq, short *genn1, short *genn2);
short rndfrom1( short x);
void GerateRandomperm (short n, short *a);
unsigned long longrndc();
double longurnd();
unsigned long llsrandc(unsigned long  seed); 
unsigned long lrandomizec(void);/* initialize randc() */
short rounds(double d);