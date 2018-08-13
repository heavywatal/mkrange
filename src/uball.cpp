#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <sstream>
#include <iostream>
#include "ucrandom.h"
#include "uball.h"
#include "read_array.hpp"

extern short xrange;
extern short minxrange;
extern short yrange;
extern short xmi;
extern short xma;
extern double Vp;
extern short nloci;
extern short nopoly;
extern double allele_effect;

std::vector<std::list<Cball>::iterator> gridindiv[162][7];

// calculate resouce use phenotypes from the genotype
double Cball::ResourceM() const {
    double sum = 0;
    for (int i = 11; i <= 10 + nloci; ++i) {
        sum += allele_effect * (gene1[i] + gene2[i]);
    }
    sum += nrnd(sqrt(Vp));
    return sum;
}

// reproduction
void Cball::nreproduction (const Cball& male, std::list<Cball>* ablist, short nogene, double mdis, double fdis, double mr, double nmr) {
    ++nomating;
    const long x = xp;// position x for focal female
    const long y = yp;// position y for focal female
    const long xx = male.xp;
    const long yy = male.yp;
    const long y2 = (y >= male.yp) ? (male.yp + yrange) : (male.yp - yrange);
    const double dist1 = sqrt((x - xx) * (x - xx) + (y - yy) * (y - yy));
    const double dist2 = sqrt((x - xx) * (x - xx) + (y - y2) * (y - y2));
    mdistance = std::min(dist1, dist2);
    const short nooffspring = fitness;
    short og1[200], og2[200];
    if (nooffspring > 0) {
        for (short j = 1; j<= nooffspring; ++j) {
            for (short k = 1; k<= nogene; ++k) {
                // inheritance genes from mother and father
                short gg = randombit();// random number 1 or 0
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
            short gg = randombit();// determin offspring sex
            const double sddispersal = (gg == 0) ? fdis : mdis;
            long nx = 0;
            long ny = 0;
            do {//// desersal
                // random number from normal distribtion with sddispersal standard deviationa and 0 mean
                const double didi = nrnd(sddispersal);
                const long d = abs(didi);
                randomove(ix, iy, d, &nx, &ny);
                if (ny > yrange) ny = ny - yrange;
                if (ny <= 0) ny = yrange + ny;
            } while ((nx <= minxrange) || (nx > xrange ) || (ny <= 0) || (ny > yrange));
            // Initialize offspring position and sex
            Cball child(nx, ny, gg);
            for (short i = 1; i <= nogene; ++i) {
                child.Igene(i, og1[i], og2[i]);// offspring genes
            }
            ablist->push_back(child);// add offspring the list
        }
    }
}


// calculate fitness
void Cball::measurefitness(double RR, double gradient, double Vs, double K, double Range, double MS) {
//Vs fitness function variance, K carring capacity,
    const short xg = (short)ceil((double)(xp) / 200.);
    const short yg = (short)ceil((double)(yp) / 200.);
    const short ff = std::max(xg - 1, 1);
    const short t = yg - 1;
    const short b = yg + 1;
    short r = xg + 1;
    if (r > xrange / 200) {r = xrange / 200;}
    short jj[4];
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
    short tot = 0;
    for (short i = ff; i <= r; ++i) {
        for (short j = 1; j <= 3; ++j) {
            const std::vector<std::list<Cball>::iterator> grid_ij = gridindiv[i][jj[j]];
            for (size_t k = 0; k < grid_ij.size(); ++k) {
                const short xx = grid_ij[k]->xp;
                const short yy = grid_ij[k]->yp;
                short y2 = 0;
                if (iy >= yy) {
                    y2 = yy + yrange;
                } else {
                    y2 = yy - yrange;
                }
                const double dist1 = sqrt((ix-xx) * (ix-xx) + (iy-yy) * (iy-yy));
                const double dist2 = sqrt((ix-xx) * (ix-xx) + (iy-y2) * (iy-y2));
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
    nloci = vecvecint[0].size()-3;
    for (size_t row = 0; row < vecvecint.size(); ++row) {
        Cball child(vecvecint[row][0], vecvecint[row][1], vecvecint[row][2]);
        for (size_t col = 3; col < vecvecint[row].size(); ++col) {
            if (vecvecint[row][col] == 0) child.Igene(col+8, 0, 0);
            if (vecvecint[row][col] == 2) child.Igene(col+8, 1, 1);
            if (vecvecint[row][col] == 1) {
                if (randombit() == 0) {
                    child.Igene(col + 8, 1, 0);
                } else {
                    child.Igene(col + 8, 0, 1);
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

void Newball2008(short n, short male, std::list<Cball>* list1, double fr) {
    short *tur;
    short *ge1,*ge2;
    ge1 = new short[n+1];
    ge2 = new short[n+1];
    tur = new short[n+1];
    // creat n individuals and initialize position and sex
    for (short i = 1; i <= n - male; i++) {
        const short xx = rndfrom1(500) + xrange / 2 - 250;
        const short yy = rndfrom1(yrange);
        list1->push_back(Cball(xx, yy, 0));
    }
    for (short i = 1; i <= male; i++) {
        const short xx = rndfrom1(500) + xrange / 2 - 250;
        const short yy = rndfrom1(yrange);
        list1->push_back(Cball(xx, yy, 1));
    }
    const short aa = rounds(fr * fr * n);
    const short ab = rounds(fr * (1 - fr) * 2 * n);
    ////// set genes for resource use ///////
    for (short j = 1; j <= 10; ++j) {
        GerateRandomperm(n, tur);
        for (short i = 1; i <= aa; ++i) {
            ge1[tur[i]] = 1;
            ge2[tur[i]] = 1;
        }
        for (short i = aa + 1; i <= aa + ab; ++i) {
            ge1[tur[i]] = 1;
            ge2[tur[i]] = 0;
        }
        for (short i = aa + ab + 1; i <= n; ++i) {
            ge1[tur[i]] = 0;
            ge2[tur[i]] = 0;
        }
        size_t i = 1;
        for (std::list<Cball>::iterator individual = list1->begin(); individual != list1->end(); ++individual, ++i) {
           individual->Igene(j, ge1[i], ge2[i]);
       }
    }
    for (short j = 11; j <= 10 + (nloci - nopoly) / 2; ++j) {
        for (std::list<Cball>::iterator individual = list1->begin(); individual != list1->end(); ++individual) {
            individual->Igene(j, 1, 1);
        }
    }
    // genes 1 and 0 at 3 loci are randomly allocated for all the individuals
    for (short j = 11 + (nloci - nopoly) / 2; j <= 10 + (nloci - nopoly) / 2 + nopoly; ++j) {
        GerateRandomperm(n, tur);
        for (short i = 1; i <= aa; ++i) {
            ge1[tur[i]] = 1;
            ge2[tur[i]] = 1;
        }
        for (short i = aa + 1; i <= aa + ab; ++i) {
            ge1[tur[i]] = 1;
            ge2[tur[i]] = 0;
        }
        for (short i = aa + ab + 1; i <= n; ++i){
            ge1[tur[i]] = 0;
            ge2[tur[i]] = 0;
        }
        size_t i = 1;
        for (std::list<Cball>::iterator individual = list1->begin(); individual != list1->end(); ++individual, ++i) {
            individual->Igene(j, ge1[i], ge2[i]);
        }
    }
    for (short j = 11 + (nloci - nopoly) / 2 + nopoly; j <= 10 + nloci; ++j) {
        for (std::list<Cball>::iterator individual = list1->begin(); individual != list1->end(); ++individual) {
            individual->Igene(j, 0, 0);
        }
    }
    delete tur;
    delete ge1;
    delete ge2;
}

///// serach for candidate mates
short matingcount (std::list<Cball>::iterator focalindiv, std::list<Cball>::iterator* matp, short matingsize) {
    const long x = focalindiv->xp;// position x for focal female
    const long y = focalindiv->yp;// position y for focal female
    const long NM = focalindiv->nocandiate;
    short dens = 0;
    long k = 0;
    double total_fitness = 0.0;
    for (long i = 1; i <= NM; ++i) {
        if (focalindiv->candidatemate[i]->fitness >0 ) {
            const long xx = focalindiv->candidatemate[i]->xp;
            const long yy = focalindiv->candidatemate[i]->yp;
            long y2 = 0;
            if (y >= yy) {
                y2 = yy + yrange;
            } else {
                y2 = yy - yrange;
            }
            const double dist1 = sqrt((x - xx) * (x - xx) + (y - yy) * (y - yy));
            const double dist2 = sqrt((x - xx) * (x - xx) + (y - y2) * (y - y2));
            const double dist = std::min(dist1, dist2);
            if (dist <= matingsize) {
                /// count and save candidate male
                ++k;
                total_fitness += focalindiv->candidatemate[i]->dfitness;
                matp[k] = focalindiv->candidatemate[i];
                ++dens;
            }
        }
    }
    if (dens > 0) {
        double sum = 0.0;
        long i = 1;
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
void SaveF(const std::list<Cball>& clist, short g, short gg, long n, short nogene) {
    std::ostringstream oss;
    oss << "File" << gg << "-" << g;
    FILE* fp = fopen(oss.str().c_str(), "w");
    int i = 0;
    for (std::list<Cball>::const_iterator it = clist.begin(); i < n; ++it, ++i) {
        const long m1 = it->sexi;
        const double m2 = it->ResourceM();
        const double m3 = (m1 == 0) ? it->fitness : it->dfitness;
        fprintf(fp, "%d\t %d\t %ld\t %f\t  %7.3f\t", it->xp, it->yp, m1, m2, m3);
        for (long ge = 1; ge <= nogene; ++ge) {
            fprintf(fp, "%d\t ", it->gene1[ge]+ it->gene2[ge]);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
}

void SaveE(short g, short gg) {
    std::ostringstream oss;
    oss << "File" << gg << "-" << g;
    FILE* fp = fopen(oss.str().c_str(), "w");
    short m1 = 0;
    fprintf(fp, "%d", m1);
    fclose(fp);
}

void SaveA(const std::list<Cball>& clist, short g, short gg, long n, short clas ) {
    std::ostringstream oss;
    oss << "Resl" << gg << "-" << g;
    FILE* fp = fopen(oss.str().c_str(), "w");
    double avfit[322];
    short nof[322];
    for (long x = 0; x<= 321; ++x) {
        avfit[x] = 0;
        nof[x] = 0;
    }
    const double w = 32000 / (double) clas;
    for (long x = 1; x <= clas; ++x) {
        int i = 0;
        for (std::list<Cball>::const_iterator it = clist.begin(); i < n; ++it, ++i) {
            if(it->sexi == 0) {
                if(it->xp > w * (x - 1) && it->xp <= w * x) {
                    avfit[x] += it->fitness;
                    nof[x] += 1;
                }
            }
        }
    }
    for (long x = 1; x <= clas; ++x) {
        const long m1 = w * (x-1);
        const long m2 = w * x;
        const double m3 = nof[x];
        double m4 = 0.0;
        if (nof[x] > 0) {
            m4 = avfit[x] / (double)nof[x];
        }
        fprintf(fp, "%ld\t %ld\t %7.4f\t %7.4f\t", m1, m2, m3, m4);
        fprintf(fp, "\n");
    }
    fclose(fp);
}

void AssignBucket(std::list<Cball>* list1) {
    const size_t nogx = xrange / 200;//the max number of x dimention of the bucket girds y = 3200 x = 160
    const size_t nogy = yrange / 200;//the max number of y dimention of the bucket girds y = 1000 y = 5
    for (size_t xc = 0; xc <= nogx; ++xc) {
        for (size_t yc = 0; yc <= nogy; ++yc) {
            //Initialize the array of the number of individuals in a bucket grid xc yc
            gridindiv[xc][yc].clear();
        }
    }
    for (std::list<Cball>::iterator it = list1->begin(); it != list1->end(); ++it) {
        const size_t xc = ceil(it->xp / 200.);
        const size_t yc = ceil(it->yp / 200.);
        gridindiv[xc][yc].push_back(it);//store the individuals into the array
    }
}
