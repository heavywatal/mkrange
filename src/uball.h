#pragma once
#ifndef MKRANGE_BALL_HPP_
#define MKRANGE_BALL_HPP_

#include <vector>
#include <list>

class Cball {
  public:
    Cball(int x, int y, int s,
          const std::vector<short>& g1, const std::vector<short>& g2)
    : xp(x), yp(y), sexi(s),
      gene1(g1), gene2(g2),
      dfitness(0.0), nooffspring(0),
      nomating(0), resource_(calc_resource()) {}
    Cball(const std::vector<int>&);

    void nreproduction(const Cball& male, std::list<Cball>* ablist, double mdis, double fdis, double mr, double nmr);
    void measurefitness(double RR, double gradient, double Vs, double K, double Range, double MS);
    unsigned matingcount(int matingsize) const;
    double resource() const {return resource_;};

    int xp, yp, sexi;
    std::vector<short> gene1, gene2;
    double dfitness;
    unsigned nooffspring;
    std::vector<std::list<Cball>::const_iterator> candidatemate;

  private:
    double calc_resource() const;
    double distance(const Cball&) const;
    std::vector<short> make_gamete(double mr, double nmr) const;

    unsigned nomating;
    double resource_;
};

void Newball(unsigned n, double fr, std::list<Cball>* list1);
void Newball(const char* infile, std::list<Cball>* list1);
void AssignBucket(const std::list<Cball>& list1);
//Individuals is stored in each bucket according to the position.
void SaveF(const std::list<Cball>& clist, unsigned g, unsigned gg, unsigned n);
void SaveA(const std::list<Cball>& clist, unsigned g, unsigned gg, unsigned n, unsigned clas);
void SaveE(unsigned g, unsigned gg);

#endif// MKRANGE_BALL_HPP_
