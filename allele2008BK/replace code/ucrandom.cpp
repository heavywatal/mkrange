
#include "uglobal.h"
#include "ucrandom.h"
#include <math.h> 
#include <time.h>
#define RAND_MAX 32767 
#define ULONG_MAX 4294967295

static unsigned long next=1; /* n๑ */
static unsigned long lnext=1; /* n๑ */

unsigned long longrndc() 
{
	lnext = lnext*1566083941UL + 1;
return lnext;
	
}

double longurnd() 
{
	return (1.0 / (ULONG_MAX + 1.0) ) * longrndc(); 
}


short randc(void) /* 0`32767 */
{
  next=next*1103515245+12345;
  return (unsigned short) (next/65536) % 32768 ;
                   
} /* rand */


unsigned long lsrandc(unsigned long  seed) /* ํ(0`ULONG_MAX) */
{
  return next=seed;
} /* lsrand */

unsigned long llsrandc(unsigned long  seed) /* ํ(0`ULONG_MAX) */
{
  return lnext=seed;
} /* lsrand */

unsigned long randomizec(void)
{
return lsrandc(time(NULL));
}

unsigned long lrandomizec(void)
{
return llsrandc(time(NULL));
}


/* ฎึ  integral RaNDom */
/* 0`x */
short rnd( short x) /* 0`RAND_MAX */
{
  return randc() % (x+1);
} /* rnd */

/* ฎึ  integral RaNDom */
short rndfrom1( short x) /* 1`RAND_MAX */
{
short s;
  do {
 s= randc() % (x+1);
 }while(s<=0);
  return s;
} /* rnd */






/* urnd.f */

/* Pสึ  Unit RaNDom */
/* 0urnd()<1 */

double urnd(void) 
{
  return randc()/(RAND_MAX+1.);
} /* urnd */


/* drnd.f */


/* ภึ  Double RaNDom */
/* 0`x */
double drnd( double x) /* ลๅl */
{
  return x/RAND_MAX*randc();
} /* drnd */

/* zrnd.f */
/* lฎึ  Zero RaNDom */
/* -x`+x */
short zrnd( short x) /* โฮl 0`RAND_MAX */
{
  return randc() % (2*x+1)-x;
} /* zrnd */


/* nrnd.f */
/* ณKึ  Normal RaNDom */
/* WณKชz*x */

double nrnd( double x) /* {ฆ */
{
  static int f=0;
  static double k,q;
  
  if(f)
  {
    f=0;
    return x*k*sin(q);
  }

  f=1;
  k=sqrt(-2*log(1-urnd()));
  q=2*3.141592653589793*urnd();
  return x*k*cos(q);
} /* nrnd */

double normalselect (double mean, double sd)

	{
		return sd * nrnd(1) + mean;
	}



/* prnd.f */
/* |A\  Poisson Random */
/* ๑^Pสิ */

short prnd( double n)  /* ฝฯ๑^Pสิ */
{
  short k=0;
  n=exp(n)*urnd();
  while(n>1)n*=urnd(),k++;
  return k;
} /* prnd */

short randombit()
{
	short			k, g;
	
		k = rnd(1000);
		g = k % 2;
		if (g == 0)
		{
			return 0;
		}
		else
		{
			return 1;
			}
	}



void randomove (long x, long y, long d,  long *ntx, long *nty)
	{
		short xspeed, yspeed;
		double angle;
		//ntx=&xx;
		//nty=&yy;
		angle = urnd()*2*3.141592653589;
		xspeed = rounds(cos(angle) * d);
		yspeed = rounds(sin(angle) * d);
		*ntx = x + xspeed;
		*nty = y + yspeed;
	}


void ArrangeFreq (short n, double freq, short *genn1, short *genn2)
{
	double 		A, B, FA, FR1, FR2;
	short		i,  AA, AB, BB;
	
		
				A = n * 2 * freq;
				B = n * 2 * (1 - freq);

				AA = 0;
				AB = 0;
				BB = 0;
				for( i = 1; i<=n ; i++)
					{
						FA = A / (A + B);
						FR1 = drnd(1);
						FR2 = drnd(1);
						if (FR1 <= FA) 
							{
								A = A - 1;
								FA = A / (A + B);
								if (FR2 <= FA) 
									{
										genn1[i] = 1;
										genn2[i] = 1;
										AA = AA + 1;
										A = A - 1;
									}
								else
									{
										genn1[ i] = 1;
										genn2[i] = 0;
										AB = AB + 1;
										B = B - 1;
									}
							}
						else
							{
								B = B - 1;
								FA = A / (A + B);
								if (FR2 <= FA) 
									{
										genn1[ i] = 0;
										genn2[ i] = 1;
										AB = AB + 1;
										A = A - 1;
									}
								else
									{
										genn1[ i] = 0;
										genn2[i] = 0;
										BB = BB + 1;
										B = B - 1;
									}
							}
					}
			

}


void GerateRandomperm (short n, short *a)
{
	short	i, j, t;
	
		for (i = 1;i<=n;i++)
			{
			a[i]= i;
			}
		for( i= n; i>=2; i--)
			{
				j = floor(urnd() * i) + 1;
				t = a[i];
				a[i] = a[j];
				a[j] = t;
			}
	}

short rounds(double d)
{
double a, b;

a = floor(d );

b=d-a;
if (b>=0.5)
	{ 
	return (short)ceil(d);
	}
	else
	{
		return a;
	}
}












