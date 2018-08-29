#include <vector>
#include <list>

class Cball {
  public:
    Cball(int x, int y, int s,
          const std::vector<short>& g1, const std::vector<short>& g2)
    : xp(x), yp(y), sexi(s),
      gene1(g1), gene2(g2),
      nomating(0),
      dfitness(0.0), nooffspring(0),
      resource_(0.0) {
        set_resource();
    }
    Cball(const std::vector<int>&);

    std::vector<short> make_gamete(double mr, double nmr) const;
    void set_resource();
    void nreproduction(const Cball& male, std::list<Cball>* ablist, double mdis, double fdis, double mr, double nmr);
    void measurefitness(double RR, double gradient, double Vs, double K, double Range, double MS);
    unsigned matingcount(int matingsize) const;
    double distance(const Cball&) const;
    double resource() const {return resource_;};

    int xp, yp, sexi;
    std::vector<short> gene1, gene2;
    unsigned nomating;
    double dfitness;
    unsigned nooffspring;
    std::vector<std::list<Cball>::iterator> candidatemate;
  private:
    double resource_;
};

void Newball(const char* infile, std::list<Cball>* list1);
void Newball2008(unsigned n, unsigned male, std::list<Cball>* list1, double fr);
void SaveF(const std::list<Cball>& clist, unsigned g, unsigned gg, unsigned n);
void SaveE(unsigned g, unsigned gg);
void SaveA(const std::list<Cball>& clist, unsigned g, unsigned gg, unsigned n, unsigned clas);
void AssignBucket(std::list<Cball>* list1);
//Individuals is stored in each bucket according to the position.
