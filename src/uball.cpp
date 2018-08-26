#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <sstream>
#include <iostream>
#include "uball.h"
#include "read_array.hpp"
#include "random.hpp"

extern int xrange;
extern int minxrange;
extern int yrange;
extern int xmi;
extern int xma;
extern double Vp;
extern int nloci;
extern int nopoly;
extern double allele_effect;
extern std::mt19937 engine;

std::vector<std::list<Cball>::iterator> gridindiv[162][7];

// calculate resouce use phenotypes from the genotype
void Cball::set_resource() {
    resource_ = 0.0;
    for (int i = 11; i <= 10 + nloci; ++i) {
        resource_ += allele_effect * (gene1[i] + gene2[i]);
    }
    std::normal_distribution<double> normal(0.0, std::sqrt(Vp));
    resource_ += normal(engine);
}

double Cball::distance(const Cball& other) const {
    const int y2 = (yp >= other.yp) ? (other.yp + yrange) : (other.yp - yrange);
    const double dist1 = std::sqrt((xp - other.xp) * (xp - other.xp) + (yp - other.yp) * (yp - other.yp));
    const double dist2 = std::sqrt((xp - other.xp) * (xp - other.xp) + (yp - y2) * (yp - y2));
    return std::min(dist1, dist2);
}

std::vector<short> Cball::make_gamete(const double mr, const double nmr) const {
    std::vector<short> gamete = gene1;
    for (size_t i=0; i<gamete.size(); ++i) {
        if (wtl::randombit(engine)) {
            gamete[i] = gene2[i];
        }
        ///////// mutation ////////
        const double mrr = (i > (10 + 1)) ? mr : nmr;
        if (wtl::bernoulli(mrr, engine)) {
            gamete[i] = (gamete[i] == 1) ? 0 : 1;
        }
    }
    return gamete;
}

// reproduction
void Cball::nreproduction(const Cball& male, std::list<Cball>* ablist, double mdis, double fdis, double mr, double nmr) {
    ++nomating;
    for (int i = 0; i < nooffspring; ++i) {
        const int gg = wtl::randombit(engine);// determin offspring sex
        const double sddispersal = (gg == 0) ? fdis : mdis;
        std::normal_distribution<double> normal(0.0, sddispersal);
        int nx = 0;
        int ny = 0;
        do {//// desersal
            // random number from normal distribtion with sddispersal standard deviationa and 0 mean
            nx = xp + normal(engine);
            ny = yp + normal(engine);
            if (ny > yrange) ny = ny - yrange;
            if (ny <= 0) ny = yrange + ny;
        } while ((nx <= minxrange) || (nx > xrange ) || (ny <= 0) || (ny > yrange));
        // Initialize offspring position and sex
        // add offspring the list
        ablist->push_back(Cball(static_cast<int>(nx), static_cast<int>(ny), gg,
                                make_gamete(mr, nmr), male.make_gamete(mr, nmr)));
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
    if (xmi > 0 || xma > 0) {// flat BK
        const int AA = xrange / 2 - xmi;
        const int BB = xma - xrange / 2;
        if (xp < xrange / 2 - AA) {
            Sx = 16 + gradient * (xp - 4000 + AA);
        } else if (xp > xrange / 2 + BB) {
            Sx = 16 + gradient * (xp - 4000 - BB);
        } else {
            Sx = 16 + gradient * (xrange / 2 - 4000);
        }
    } else if (xmi < 0) {// linear
        Sx = 16 + gradient * (xp - 4000);
    } else {// steep 2008BK
        Sx = 64 + gradient * (xp - 16000) * (xp - 16000) * (xp - 16000) / 1000000.;
    }
    dfitness = 2 + RR * (1 - tot / K) - (Sx - resource()) * (Sx - resource()) / (2 * Vs);
    if (dfitness < 0) dfitness = 0;
    std::poisson_distribution<int> poisson(dfitness);
    nooffspring = poisson(engine);
}

Cball::Cball(const std::vector<int>& row)
: xp(row[0]), yp(row[1]), sexi(row[2]),
  gene1(row.size() - 3 + 10 + 1), gene2(gene1.size()),
  nomating(0), dfitness(0.0), nooffspring(0), resource_(0.0) {
    for (size_t col = 3; col < row.size(); ++col) {
        const int col_8 = static_cast<int>(col) + 8;
        if (row[col] == 2) Igene(col_8, 1, 1);
        else if (row[col] == 1) {
            if (wtl::randombit(engine) == 0) {
                Igene(col_8, 1, 0);
            } else {
                Igene(col_8, 0, 1);
            }
        }
    }
    //////////// set genes for resource use ///////
    for (int j = 1; j <= 10; ++j) {
        Igene(j, wtl::randombit(engine), wtl::randombit(engine));
    }
    set_resource();
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
    std::uniform_int_distribution<int> uniform_x(1, 500);
    std::uniform_int_distribution<int> uniform_y(1, yrange);
    // creat n individuals and initialize position and sex
    for (int i = 1; i <= n; ++i) {
        const int xx = uniform_x(engine) + xrange / 2 - 250;
        const int yy = uniform_y(engine);
        list1->push_back(Cball(xx, yy, i <= male, std::vector<short>(nloci + 10 + 1), std::vector<short>(nloci + 10 + 1)));
    }
    const int aa = std::round(fr * fr * n);
    const int ab = std::round(fr * (1 - fr) * 2 * n);
    std::vector<short> ge1(n);
    std::vector<short> ge2(n);
    for (int i = 0; i < n; ++i) {
        ge1[i] = (i < aa + ab);
        ge2[i] = (i < aa);
    }
    ////// set genes for resource use ///////
    for (int j = 1; j <= 10; ++j) {
        std::shuffle(ge1.begin(), ge1.end(), engine);
        std::shuffle(ge2.begin(), ge2.end(), engine);
        size_t i = 0;
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
        std::shuffle(ge1.begin(), ge1.end(), engine);
        std::shuffle(ge2.begin(), ge2.end(), engine);
        size_t i = 0;
        for (std::list<Cball>::iterator individual = list1->begin(); individual != list1->end(); ++individual, ++i) {
            individual->Igene(j, ge1[i], ge2[i]);
        }
    }
    for (int j = 11 + (nloci - nopoly) / 2 + nopoly; j <= 10 + nloci; ++j) {
        for (std::list<Cball>::iterator individual = list1->begin(); individual != list1->end(); ++individual) {
            individual->Igene(j, 0, 0);
        }
    }
    for (std::list<Cball>::iterator individual = list1->begin(); individual != list1->end(); ++individual) {
        individual->set_resource();
    }
}

///// serach for candidate mates
size_t Cball::matingcount(int matingsize) const {
    double total_fitness = 0.0;
    std::vector<size_t> indices;
    indices.reserve(candidatemate.size());
    for (size_t i = 0; i < candidatemate.size(); ++i) {
        if (candidatemate[i]->nooffspring > 0) {
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
        std::uniform_real_distribution<double> uniform(0.0, total_fitness);
        const double r = uniform(engine);
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
        const double m2 = it->resource();
        const double m3 = (m1 == 0) ? it->nooffspring : it->dfitness;
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
                    avfit[x] += it->nooffspring;
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
