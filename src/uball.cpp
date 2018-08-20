#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <sstream>
#include <iostream>
#include "ucrandom.h"
#include "uball.h"
#include "read_array.hpp"

extern int xrange;
extern int minxrange;
extern int yrange;
extern int xmi;
extern int xma;
extern double Vp;
extern int nloci;
extern int nopoly;
extern double allele_effect;

std::vector<std::list<Cball>::iterator> gridindiv[162][7];

// calculate resouce use phenotypes from the genotype
double Cball::ResourceM() const {
    double sum = 0;
    for (int i = 11; i <= 10 + nloci; ++i) {
        sum += allele_effect * (gene1[i] + gene2[i]);
    }
    sum += nrnd(std::sqrt(Vp));
    return sum;
}

double Cball::distance(const Cball& other) const {
    const int y2 = (yp >= other.yp) ? (other.yp + yrange) : (other.yp - yrange);
    const double dist1 = std::sqrt((xp - other.xp) * (xp - other.xp) + (yp - other.yp) * (yp - other.yp));
    const double dist2 = std::sqrt((xp - other.xp) * (xp - other.xp) + (yp - y2) * (yp - y2));
    return std::min(dist1, dist2);
}

// reproduction
void Cball::nreproduction (const Cball& male, std::list<Cball>* ablist, int nogene, double mdis, double fdis, double mr, double nmr) {
    ++nomating;
    const int nooffspring = static_cast<int>(fitness);
    std::vector<short> og1(nogene + 1);
    std::vector<short> og2(nogene + 1);
    if (nooffspring > 0) {
        for (int j = 1; j<= nooffspring; ++j) {
            for (int k = 1; k<= nogene; ++k) {
                // inheritance genes from mother and father
                og1[k] = randombit() ? gene1[k] : gene2[k];
                og2[k] = randombit() ? male.gene1[k] : male.gene2[k];
                ///////// mutation ////////
                const double mrr = (k > 10) ? mr : nmr;
                if (longurnd() < mrr) {
                    og1[k] = (og1[k] == 1) ? 0 : 1;
                }
                if (longurnd() < mrr) {
                    og2[k] = (og2[k] == 1) ? 0 : 1;
                }
            }
            const int gg = randombit();// determin offspring sex
            const double sddispersal = (gg == 0) ? fdis : mdis;
            long nx = 0;
            long ny = 0;
            do {//// desersal
                // random number from normal distribtion with sddispersal standard deviationa and 0 mean
                const double didi = nrnd(sddispersal);
                const long d = static_cast<long>(std::abs(didi));
                randomove(xp, yp, d, &nx, &ny);
                if (ny > yrange) ny = ny - yrange;
                if (ny <= 0) ny = yrange + ny;
            } while ((nx <= minxrange) || (nx > xrange ) || (ny <= 0) || (ny > yrange));
            // Initialize offspring position and sex
            Cball child(static_cast<int>(nx), static_cast<int>(ny), gg);
            for (int i = 1; i <= nogene; ++i) {
                child.Igene(i, og1[i], og2[i]);// offspring genes
            }
            ablist->push_back(child);// add offspring the list
        }
    }
}


// calculate fitness
void Cball::measurefitness(double RR, double gradient, double Vs, double K, double Range, double MS) {
//Vs fitness function variance, K carring capacity,
    const int xg = static_cast<int>(ceil(xp / 200.0));
    const int yg = static_cast<int>(ceil(yp / 200.0));
    int jj[4];
    if (yg - 1 < 1) {
        jj[1] = 1;
        jj[2] = 2;
        jj[3] = 5;
    } else if (yg + 1 > yrange / 200) {
        jj[1] = 4;
        jj[2] = 5;
        jj[3] = 1;
    } else {
        jj[1] = yg - 1;
        jj[2] = yg;
        jj[3] = yg + 1;
    }
    int tot = 0;
    for (int i = std::max(xg - 1, 1); i <= std::min(xg + 1, xrange / 200); ++i) {
        for (int j = 1; j <= 3; ++j) {
            const std::vector<std::list<Cball>::iterator> grid_ij = gridindiv[i][jj[j]];
            for (size_t k = 0; k < grid_ij.size(); ++k) {
                const double dist = distance(*grid_ij[k]);
                if (dist <= Range) {
                    //count number of individuals within the range
                    ++tot;
                }
                // count number of mailes within the ara of radius 200
                if (sexi == 0 && dist <= MS && grid_ij[k]->sexi == 1) {
                    // if this individul is female and called individual is male and if the distance between them is less than 200, then the called individual is recored as candidate males.
                    candidatemate.push_back(grid_ij[k]);
                    //NOTE: should be refreshed every generation?
                }
            }
        }
    }
    //// calculate fitness
    double Sx = 0.0;
    if (xmi > 0 || xma > 0) {// BK
        const int AA = xrange / 2 - xmi;
        const int BB = xma - xrange / 2;
        if (xp < xrange / 2 - AA) {
            Sx = 16 + gradient * (xp - 4000 + AA);
        } else if (xp > xrange / 2 + BB) {
            Sx = 16 + gradient * (xp - 4000 - BB);
        } else {
            Sx = 16 + gradient * (xrange / 2 - 4000);
        }
    } else {// 2008BK
        Sx = 64 + gradient * (xp - 16000) * (xp - 16000) * (xp - 16000) / 1000000.;
    }
    const double Sx_resource = Sx - ResourceM();
    dfitness = 2 + RR * (1 - tot / K) - Sx_resource * Sx_resource / (2 * Vs);
    if (dfitness < 0) dfitness = 0;
    fitness = prnd(dfitness);
}

Cball::Cball(const std::vector<int>& row)
: xp(row[0]), yp(row[1]), sexi(row[2]),
  nomating(0), fitness(0), dfitness(0) {
    for (size_t col = 3; col < row.size(); ++col) {
        const int col_8 = static_cast<int>(col) + 8;
        if (row[col] == 0) Igene(col_8, 0, 0);
        if (row[col] == 2) Igene(col_8, 1, 1);
        if (row[col] == 1) {
            if (randombit() == 0) {
                Igene(col_8, 1, 0);
            } else {
                Igene(col_8, 0, 1);
            }
        }
    }
    //////////// set genes for resource use ///////
    for (int j = 1; j <= 10; ++j) {
        Igene(j, randombit(), randombit());
    }
}

// Create new indiviaul
void Newball(const char* infile, std::list<Cball>* list1) {
    //// creat n individuals and initialize position and sex
    const std::vector<std::vector<int> > matrix = read_int_array(infile);
    nloci = static_cast<int>(matrix[0].size() - 3);
    for (size_t row = 0; row < matrix.size(); ++row) {
        list1->push_back(Cball(matrix[row]));
    }
}

void Newball2008(int n, int male, std::list<Cball>* list1, double fr) {
    short *tur = new short[n+1];
    std::vector<short> ge1(n + 1);
    std::vector<short> ge2(n + 1);
    // creat n individuals and initialize position and sex
    for (int i = 1; i <= n; ++i) {
        const int xx = rndfrom1(500) + xrange / 2 - 250;
        const int yy = rndfrom1(static_cast<short>(yrange));
        list1->push_back(Cball(xx, yy, i <= male));
    }
    const int aa = rounds(fr * fr * n);
    const int ab = rounds(fr * (1 - fr) * 2 * n);
    ////// set genes for resource use ///////
    for (int j = 1; j <= 10; ++j) {
        GerateRandomperm(static_cast<short>(n), tur);
        for (int i = 1; i <= n; ++i) {
            ge1[tur[i]] = (i <= aa + ab);
            ge2[tur[i]] = (i <= aa);
        }
        size_t i = 1;
        for (std::list<Cball>::iterator individual = list1->begin(); individual != list1->end(); ++individual, ++i) {
           individual->Igene(j, ge1[i], ge2[i]);
       }
    }
    for (int j = 11; j <= 10 + (nloci - nopoly) / 2; ++j) {
        for (std::list<Cball>::iterator individual = list1->begin(); individual != list1->end(); ++individual) {
            individual->Igene(j, 1, 1);
        }
    }
    // genes 1 and 0 at 3 loci are randomly allocated for all the individuals
    for (int j = 11 + (nloci - nopoly) / 2; j <= 10 + (nloci - nopoly) / 2 + nopoly; ++j) {
        GerateRandomperm(static_cast<short>(n), tur);
        for (int i = 1; i <= n; ++i) {
            ge1[tur[i]] = (i <= aa + ab);
            ge2[tur[i]] = (i <= aa);
        }
        size_t i = 1;
        for (std::list<Cball>::iterator individual = list1->begin(); individual != list1->end(); ++individual, ++i) {
            individual->Igene(j, ge1[i], ge2[i]);
        }
    }
    for (int j = 11 + (nloci - nopoly) / 2 + nopoly; j <= 10 + nloci; ++j) {
        for (std::list<Cball>::iterator individual = list1->begin(); individual != list1->end(); ++individual) {
            individual->Igene(j, 0, 0);
        }
    }
    delete[] tur;
}

///// serach for candidate mates
size_t Cball::matingcount(int matingsize) const {
    double total_fitness = 0.0;
    std::vector<size_t> indices;
    indices.reserve(candidatemate.size());
    for (size_t i = 0; i < candidatemate.size(); ++i) {
        if (candidatemate[i]->fitness > 0) {
            const double dist = distance(*candidatemate[i]);
            if (dist <= matingsize) {
                /// count and save candidate male
                total_fitness += candidatemate[i]->dfitness;
                indices.push_back(i);
            }
        }
    }
    if (total_fitness > 0.0) {
        double sum = 0.0;
        unsigned i = 0;
        const double r = urnd() * total_fitness;
        do {
            sum += candidatemate[indices[i]]->dfitness;
            ++i;
        } while (sum < r);
        return indices[--i];
    }
    return candidatemate.size();
}

// save the results as afile
void SaveF(const std::list<Cball>& clist, int g, int gg, size_t n, int nogene) {
    std::ostringstream oss;
    oss << "File" << gg << "-" << g;
    FILE* fp = fopen(oss.str().c_str(), "w");
    std::list<Cball>::const_iterator it = clist.begin();
    for (size_t i = 0; i < n; ++it, ++i) {
        const int m1 = it->sexi;
        const double m2 = it->ResourceM();
        const double m3 = (m1 == 0) ? it->fitness : it->dfitness;
        fprintf(fp, "%d\t %d\t %d\t %f\t  %7.3f\t", it->xp, it->yp, m1, m2, m3);
        for (long ge = 1; ge <= nogene; ++ge) {
            fprintf(fp, "%d\t ", it->gene1[ge]+ it->gene2[ge]);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
}

void SaveE(int g, int gg) {
    std::ostringstream oss;
    oss << "File" << gg << "-" << g;
    FILE* fp = fopen(oss.str().c_str(), "w");
    fprintf(fp, "%d", 0);
    fclose(fp);
}

void SaveA(const std::list<Cball>& clist, int g, int gg, size_t n, int clas) {
    std::ostringstream oss;
    oss << "Resl" << gg << "-" << g;
    FILE* fp = fopen(oss.str().c_str(), "w");
    std::vector<double> avfit(clas);
    std::vector<int> nof(clas);
    const double w = 32000.0 / clas;
    for (int x = 0; x < clas; ++x) {
        std::list<Cball>::const_iterator it = clist.begin();
        for (size_t i = 0; i < n; ++it, ++i) {
            if(it->sexi == 0) {
                if(w * x < it->xp && it->xp <= w * (x + 1)) {
                    avfit[x] += it->fitness;
                    nof[x] += 1;
                }
            }
        }
    }
    for (int x = 0; x < clas; ++x) {
        const int m1 = static_cast<int>(w * x);
        const int m2 = static_cast<int>(w * (x + 1));
        const double m3 = static_cast<double>(nof[x]);
        const double m4 = m3 > 0.0 ? avfit[x] / m3 : 0.0;
        fprintf(fp, "%d\t %d\t %7.4f\t %7.4f\t", m1, m2, m3, m4);
        fprintf(fp, "\n");
    }
    fclose(fp);
}

void AssignBucket(std::list<Cball>* list1) {
    const int nogx = xrange / 200;//the max number of x dimention of the bucket girds y = 3200 x = 160
    const int nogy = yrange / 200;//the max number of y dimention of the bucket girds y = 1000 y = 5
    for (int xc = 0; xc <= nogx; ++xc) {
        for (int yc = 0; yc <= nogy; ++yc) {
            //Initialize the array of the number of individuals in a bucket grid xc yc
            gridindiv[xc][yc].clear();
        }
    }
    for (std::list<Cball>::iterator it = list1->begin(); it != list1->end(); ++it) {
        const int xc = static_cast<int>(ceil(it->xp / 200.));
        const int yc = static_cast<int>(ceil(it->yp / 200.));
        gridindiv[xc][yc].push_back(it);//store the individuals into the array
    }
}
