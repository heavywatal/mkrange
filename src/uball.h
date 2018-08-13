#include <cstdio>
#include <list>

class Cball {
    typedef short gene[140];
  public:
    Cball (short x, short y, short s)
    : xp(x), yp(y), ix(x), iy(y), sexi(s),
      nomating(0), nocandiate(0),
      fitness(0), dfitness(0), mdistance(0) {}
    // initialize genes
    void Igene(short nthgene, short gf, short gl) {
         gene1[nthgene] = gf;
         gene2[nthgene] = gl;
         genotype[nthgene] = gf + gl;
    }
    double ResourceM() const;
    void nreproduction(std::list<Cball>* ablist, short nogene, double mdis, double fdis, double mr, double nmr);
    void measurefitness(double RR, double gradient, double Vs, double K, double Range, double MS);

    short xp, yp, ix, iy, sexi;
    short nomating, nocandiate;
    double fitness, dfitness, mdistance;
    gene gene1, gene2, genotype;
    std::list<Cball>::iterator candidatemate[1500];
};

void Newball(std::list<Cball>* list1);
void Newball2008(short n, short male, std::list<Cball>* list1, double fr);
short matingcount(std::list<Cball>::iterator focalindiv, std::list<Cball>::iterator* matp, short matingsize);
void SaveF(const std::list<Cball>& clist, short g, short gg, long n, short nogene);
short Resource(short a, short b);
void SaveE(short g, short gg);
void SaveA(const std::list<Cball>& clist, short g, short gg, long n, short clas );
void AssignBucket(std::list<Cball>* list1);
//Individuals is stored in each bucket according to the position.
