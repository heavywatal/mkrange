#include <cstdio>
#include <list>

class Cball {
    typedef short gene[140];
  public:
    short xp, yp, ix, iy, sexi;
    short nomating, nocandiate;
    double fitness, dfitness, mdistance;
    gene gene1, gene2, genotype;
    std::list<Cball>::iterator candidatemate[1500];
    void Iball(short x, short y, short s);
    double ResourceM();
    void Igene(short nthgene, short gf, short gl);
    void nreproduction (std::list<Cball>* ablist, short nogene, double mdis, double fdis, double mr, double nmr);
    void measurefitness(double RR, double gradient, double Vs, double K, double Range, double MS);
};

void Newball(std::list<Cball>* list1);
void Newball2008(short n, short male, std::list<Cball>* list1, double fr);
void matingcount(std::list<Cball>::iterator focalindiv, std::list<Cball>::iterator* matp, short matingsize, short* dens);
void SaveF(std::list<Cball>* clist, short g, short gg, long n, short nogene);
short Resource(short a, short b);
void SaveE(short g, short gg);
void SaveA(std::list<Cball>* clist, short g, short gg, long n, short clas );
void AssignBucket(std::list<Cball>* list1);
//Individuals is stored in each bucket according to the position.
