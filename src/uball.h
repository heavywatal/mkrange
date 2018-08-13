#include "uglobal.h"
#include <cstdio>
#include <list>

class Cball {				/* Class Declaration				*/

public:
		short xp,yp, ix,iy, homerangesize, sexi,mx,my;
		short	encounter,nomating,mateP,nocandiate;//mother,
		double	ResouceP, MatingP,ChoiceP, ReprodIsolP,fitness,dfitness,mdistance;
		gene	gene1, gene2,genotype;
	      std::list<Cball>::iterator candidatemate[1500];
		void Iball(short x, short y, short s);
		double ResourceM();
		//short choice();
		//short choiceM();
		//short select();
		void Igene(short nthgene, short gf, short gl);
		void nreproduction (std::list<Cball> *ablist,  short nogene, double mdis, double fdis, double mr, double nmr);
		void measurefitness(std::list<Cball>  *clist,  double RR, double gradient,double Vs,double K,double Range,double MS);
		};

/*******************function********************/


void Newball(short n, short male, std::list<Cball>  *list1, double fr);
void Newball2008(short n, short male, std::list<Cball>  *list1, double fr);
void matingcount(std::list<Cball>  *clist, std::list<Cball>::iterator focalindiv,std::list<Cball>::iterator *matp, long noi,  short matingsize,  short *dens);
void SaveF(std::list<Cball>  *clist,short g,short gg,double md, long n);
short Resource(short a, short b);
void SaveE(short g,short gg);
void SaveA(std::list<Cball>  *clist,short g,short gg,double md, long n, short clas );
void AssignBucket(std::list<Cball>  *list1);//Individuals is stored in each bucket according to the position.

