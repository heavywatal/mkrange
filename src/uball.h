#include <cstdio>
#include <vector>
#include <list>

class Cball {
    typedef short gene[140];
  public:
    Cball(int x, int y, int s)
    : xp(x), yp(y), sexi(s),
      nomating(0),
      fitness(0), dfitness(0), mdistance(0) {}
    Cball(const std::vector<int>&);

    // initialize genes
    void Igene(int nthgene, short gf, short gl) {
         gene1[nthgene] = gf;
         gene2[nthgene] = gl;
         genotype[nthgene] = gf + gl;
    }
    double ResourceM() const;
    void nreproduction(const Cball& male, std::list<Cball>* ablist, int nogene, double mdis, double fdis, double mr, double nmr);
    void measurefitness(double RR, double gradient, double Vs, double K, double Range, double MS);
    size_t matingcount(int matingsize) const;
    double distance(const Cball&) const;

    int xp, yp, sexi;
    int nomating;
    double fitness, dfitness, mdistance;
    gene gene1, gene2, genotype;
    std::vector<std::list<Cball>::iterator> candidatemate;
};

void Newball(const char* infile, std::list<Cball>* list1);
void Newball2008(int n, int male, std::list<Cball>* list1, double fr);
void SaveF(const std::list<Cball>& clist, int g, int gg, size_t n, int nogene);
void SaveE(int g, int gg);
void SaveA(const std::list<Cball>& clist, int g, int gg, size_t n, int clas);
void AssignBucket(std::list<Cball>* list1);
//Individuals is stored in each bucket according to the position.
