#include "uglobal.h"
#include <stdio.h>
#include "list.h"

class Cball {				/* Class Declaration				*/

public:
		short xp,yp, ix,iy, homerangesize, sexi,mx,my;
		short	ages, encounter,mother,nomating,mateP,nocandiate;
		double	ResouceP, MatingP,ChoiceP, ReprodIsolP,fitness,dfitness,mdistance;
		gene	gene1, gene2,genotype;
	      list<Cball>::iterator candidatemate[3000];


		void Iball(short x, short y, short s);

		double ResourceM();
		short choice();
		short choiceM();
		short select();
		void Igene(short nthgene, short gf, short gl);
		void nreproduction (list<Cball> *ablist,  short nogene, double mdis, double fdis, double mr, double nmr);
		void measurefitness(list<Cball>  *clist,  double RR, double gradient,double Vs,double K,double Range, double MS);
		};

/*******************function********************/


void Newball(short n, short male, list<Cball>  *list1, double fr);
void matingcount(list<Cball>  *clist, list<Cball>::iterator focalindiv,list<Cball>::iterator *matp, long noi,  short matingsize,  short *dens);
void SaveF(list<Cball>  *clist,short g,short gg,double md, long n);
short Resource(short a, short b);
void SaveE(short g,short gg);
void SaveA(list<Cball>  *clist,short g,short gg,double md, long n, short clas );
void AssignBucket(list<Cball>  *list1);//Individuals is stored in each bucket according to the position.

