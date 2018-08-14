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

// reproduction
void Cball::nreproduction (const Cball& male, std::list<Cball>* ablist, int nogene, double mdis, double fdis, double mr, double nmr) {
    ++nomating;
    const long xx = male.xp;
    const long yy = male.yp;
    const long y2 = (yp >= male.yp) ? (male.yp + yrange) : (male.yp - yrange);
    const double dist1 = std::sqrt((xp - xx) * (xp - xx) + (yp - yy) * (yp - yy));
    const double dist2 = std::sqrt((xp - xx) * (xp - xx) + (yp - y2) * (yp - y2));
    mdistance = std::min(dist1, dist2);
    const int nooffspring = static_cast<int>(fitness);
    short og1[200], og2[200];
    if (nooffspring > 0) {
        for (int j = 1; j<= nooffspring; ++j) {
            for (int k = 1; k<= nogene; ++k) {
                // inheritance genes from mother and father
                int gg = randombit();// random number 1 or 0
                if (gg == 1) {
                    og1[k] = gene1[k];
                } else {
                    og1[k] = gene2[k];
                }
                gg = randombit();
                if (gg == 1) {
                    og2[k] = male.gene1[k];
                } else {
                    og2[k] = male.gene2[k];
                }
                ///////// mutation ////////
                double mrr = 0.0;
                if (k > 10) {
                    mrr = mr;
                } else {
                    mrr = nmr;
                }
                double mu = longurnd();
                if (mrr >= mu) {
                    if (og1[k] == 1) {
                        og1[k] = 0;
                    } else {
                        og1[k] = 1;
                    }
                }
                mu = longurnd();
                if (mrr >= mu) {
                    if (og2[k] == 1) {
                        og2[k] = 0;
                    } else {
                        og2[k] = 1;
                    }
                }
            }
            int gg = randombit();// determin offspring sex
            const double sddispersal = (gg == 0) ? fdis : mdis;
            long nx = 0;
            long ny = 0;
            do {//// desersal
                // random number from normal distribtion with sddispersal standard deviationa and 0 mean
                const double didi = nrnd(sddispersal);
                const long d = static_cast<long>(std::abs(didi));
                randomove(ix, iy, d, &nx, &ny);
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
    const int ff = std::max(xg - 1, 1);
    const int t = yg - 1;
    const int b = yg + 1;
    int r = xg + 1;
    if (r > xrange / 200) {r = xrange / 200;}
    int jj[4];
    if (t < 1) {
        jj[1] = 1;
        jj[2] = 2;
        jj[3] = 5;
    } else if (b > yrange / 200) {
        jj[1] = 4;
        jj[2] = 5;
        jj[3] = 1;
    } else {
        jj[1] = yg - 1;
        jj[2] = yg;
        jj[3] = yg + 1;
    }
    int tot = 0;
    for (int i = ff; i <= r; ++i) {
        for (int j = 1; j <= 3; ++j) {
            const std::vector<std::list<Cball>::iterator> grid_ij = gridindiv[i][jj[j]];
            for (size_t k = 0; k < grid_ij.size(); ++k) {
                const int xx = grid_ij[k]->xp;
                const int yy = grid_ij[k]->yp;
                int y2 = 0;
                if (iy >= yy) {
                    y2 = yy + yrange;
                } else {
                    y2 = yy - yrange;
                }
                const double dist1 = std::sqrt((ix-xx) * (ix-xx) + (iy-yy) * (iy-yy));
                const double dist2 = std::sqrt((ix-xx) * (ix-xx) + (iy-y2) * (iy-y2));
                const double dist = std::min(dist1, dist2);
                if (dist <= Range) {
                    //count number of individuals within the range
                    ++tot;
                }
                // count number of mailes within the ara of radius 200
                if (sexi == 0 && dist <= MS && grid_ij[k]->sexi == 1) {
                    // if this individul is female and called individual is male and if the distance between them is less than 200, then the called individual is recored as candidate males.
                    ++nocandiate;
                    // store the canditate males
                    candidatemate[nocandiate] = grid_ij[k];
                }
            }
        }
    }
    //// calculate fitness
    double Sx = 0.0;
    if (xmi > 0 || xma > 0) {// BK
        int AA = xrange / 2 - xmi;
        int BB = xma - xrange / 2;
        if (ix >= xrange / 2 - AA && ix <= xrange / 2 + BB ) Sx = 16 + gradient * (xrange / 2 - 4000);
        if (ix < xrange / 2 - AA) Sx = 16 + gradient * (ix - 4000 + AA);
        if (ix > xrange / 2 + BB) Sx = 16 + gradient * (ix - 4000 - BB);
    } else {// 2008BK
        Sx = 64 + gradient * (ix - 16000) * (ix - 16000) * (ix - 16000) / 1000000.;
    }
    dfitness = 2 + RR * (1 - tot / K) - (Sx - ResourceM()) * (Sx - ResourceM()) / (2 * Vs);
    if (dfitness <= 0) dfitness = 0;
    fitness = prnd(dfitness);
}

// Create new indiviaul
void Newball(std::list<Cball>* list1) {
    //// creat n individuals and initialize position and sex
    const std::vector<std::vector<int> > vecvecint = read_int_array("TestInput.txt");
    nloci = static_cast<int>(vecvecint[0].size() - 3);
    for (size_t row = 0; row < vecvecint.size(); ++row) {
        Cball child(vecvecint[row][0], vecvecint[row][1], vecvecint[row][2]);
        for (size_t col = 3; col < vecvecint[row].size(); ++col) {
            const int col_8 = static_cast<int>(col) + 8;
            if (vecvecint[row][col] == 0) child.Igene(col_8, 0, 0);
            if (vecvecint[row][col] == 2) child.Igene(col_8, 1, 1);
            if (vecvecint[row][col] == 1) {
                if (randombit() == 0) {
                    child.Igene(col_8, 1, 0);
                } else {
                    child.Igene(col_8, 0, 1);
                }
            }
        }
        list1->push_back(child);
    }
    //////////// set genes for resource use ///////
    for (std::list<Cball>::iterator individual = list1->begin(); individual != list1->end(); ++individual) {
        for (int j = 1; j <= 10; ++j) {
            individual->Igene(j, randombit(), randombit());
        }
    }
}

void Newball2008(int n, int male, std::list<Cball>* list1, double fr) {
    short *tur = new short[n+1];
    short *ge1 = new short[n+1];
    short *ge2 = new short[n+1];
    // creat n individuals and initialize position and sex
    for (int i = 1; i <= n - male; i++) {
        const int xx = rndfrom1(500) + xrange / 2 - 250;
        const int yy = rndfrom1(static_cast<short>(yrange));
        list1->push_back(Cball(xx, yy, 0));
    }
    for (int i = 1; i <= male; i++) {
        const int xx = rndfrom1(500) + xrange / 2 - 250;
        const int yy = rndfrom1(static_cast<short>(yrange));
        list1->push_back(Cball(xx, yy, 1));
    }
    const int aa = rounds(fr * fr * n);
    const int ab = rounds(fr * (1 - fr) * 2 * n);
    ////// set genes for resource use ///////
    for (int j = 1; j <= 10; ++j) {
        GerateRandomperm(static_cast<short>(n), tur);
        for (int i = 1; i <= aa; ++i) {
            ge1[tur[i]] = 1;
            ge2[tur[i]] = 1;
        }
        for (int i = aa + 1; i <= aa + ab; ++i) {
            ge1[tur[i]] = 1;
            ge2[tur[i]] = 0;
        }
        for (int i = aa + ab + 1; i <= n; ++i) {
            ge1[tur[i]] = 0;
            ge2[tur[i]] = 0;
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
        for (int i = 1; i <= aa; ++i) {
            ge1[tur[i]] = 1;
            ge2[tur[i]] = 1;
        }
        for (int i = aa + 1; i <= aa + ab; ++i) {
            ge1[tur[i]] = 1;
            ge2[tur[i]] = 0;
        }
        for (int i = aa + ab + 1; i <= n; ++i){
            ge1[tur[i]] = 0;
            ge2[tur[i]] = 0;
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
    delete[] ge1;
    delete[] ge2;
}

///// serach for candidate mates
int Cball::matingcount(std::list<Cball>::iterator* matp, int matingsize) const {
    int dens = 0;
    long k = 0;
    double total_fitness = 0.0;
    for (long i = 1; i <= nocandiate; ++i) {
        if (candidatemate[i]->fitness > 0) {
            const long xx = candidatemate[i]->xp;
            const long yy = candidatemate[i]->yp;
            long y2 = 0;
            if (yp >= yy) {
                y2 = yy + yrange;
            } else {
                y2 = yy - yrange;
            }
            const double dist1 = std::sqrt((xp - xx) * (xp - xx) + (yp - yy) * (yp - yy));
            const double dist2 = std::sqrt((xp - xx) * (xp - xx) + (yp - y2) * (yp - y2));
            const double dist = std::min(dist1, dist2);
            if (dist <= matingsize) {
                /// count and save candidate male
                ++k;
                total_fitness += candidatemate[i]->dfitness;
                matp[k] = candidatemate[i];
                ++dens;
            }
        }
    }
    if (dens > 0) {
        double sum = 0.0;
        int i = 1;
        const double r = urnd() * total_fitness;
        do {
            sum += matp[i]->dfitness;
            ++i;
        } while (sum < r && i <= k);
        dens = i - 1;
    }
    return dens;
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
