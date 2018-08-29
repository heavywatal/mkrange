#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <sstream>
#include <iostream>
#include "uball.h"
#include "global.hpp"
#include "read_array.hpp"
#include "random.hpp"

std::vector<std::list<Cball>::iterator> gridindiv[162][7];

Cball::Cball(const std::vector<int>& row)
: xp(row[0]), yp(row[1]), sexi(row[2]),
  gene1(nloci_neutral + row.size() - 3), gene2(gene1.size()),
  dfitness(0.0), nooffspring(0), nomating(0), resource_(0.0) {
    //////////// set genes for resource use ///////
    unsigned i = 0;
    for (; i < nloci_neutral; ++i) {
        gene1[i] = wtl::randombit(engine);
        gene2[i] = wtl::randombit(engine);
    }
    for (unsigned col = 3; col < row.size(); ++col, ++i) {
        if (row[col] == 2) {
            gene1[i] = 1;
            gene2[i] = 1;
        }
        else if (row[col] == 1) {
            if (wtl::randombit(engine)) {
                gene1[i] = 1;
            } else {
                gene2[i] = 1;
            }
        }
    }
    resource_ = calc_resource();
}

// calculate resouce use phenotypes from the genotype
double Cball::calc_resource() const {
    double x = 0.0;
    for (unsigned i = nloci_neutral; i < nloci_neutral + nloci; ++i) {
        x += allele_effect * (gene1[i] + gene2[i]);
    }
    std::normal_distribution<double> normal(0.0, std::sqrt(Vp));
    return x += normal(engine);
}

double Cball::distance(const Cball& other) const {
    const int y2 = (yp >= other.yp) ? (other.yp + yrange) : (other.yp - yrange);
    const double dist1 = std::sqrt((xp - other.xp) * (xp - other.xp) + (yp - other.yp) * (yp - other.yp));
    const double dist2 = std::sqrt((xp - other.xp) * (xp - other.xp) + (yp - y2) * (yp - y2));
    return std::min(dist1, dist2);
}

std::vector<short> Cball::make_gamete(const double mr, const double nmr) const {
    std::vector<short> gamete = gene1;
    for (unsigned i=0; i<gamete.size(); ++i) {
        if (wtl::randombit(engine)) {
            gamete[i] = gene2[i];
        }
        ///////// mutation ////////
        const double mrr = (i < nloci_neutral) ? nmr : mr;
        if (wtl::bernoulli(mrr, engine)) {
            gamete[i] = (gamete[i] == 1) ? 0 : 1;
        }
    }
    return gamete;
}

// reproduction
void Cball::nreproduction(const Cball& male, std::list<Cball>* ablist, double mdis, double fdis, double mr, double nmr) {
    ++nomating;
    constexpr double PI = 3.14159265358979323846;
    std::uniform_real_distribution<double> angle_dist(0.0, 2 * PI);
    for (unsigned i = 0; i < nooffspring; ++i) {
        const int gg = wtl::randombit(engine);// determin offspring sex
        const double sddispersal = (gg == 0) ? fdis : mdis;
        std::normal_distribution<double> radius_dist(0.0, sddispersal);
        int nx = 0;
        int ny = 0;
        do {//// desersal
            // random number from normal distribtion with sddispersal standard deviationa and 0 mean
            const double radius = radius_dist(engine);
            const double angle = angle_dist(engine);
            nx = xp + static_cast<int>(radius * std::cos(angle));
            ny = yp + static_cast<int>(radius * std::sin(angle));
            if (ny <= 0) ny += yrange;
            if (ny > yrange) ny -= yrange;
        } while ((nx <= minxrange) || (nx > xrange) || (ny <= 0) || (ny > yrange));
        // Initialize offspring position and sex
        // add offspring the list
        ablist->push_back(Cball(nx, ny, gg, make_gamete(mr, nmr), male.make_gamete(mr, nmr)));
    }
}


// calculate fitness
void Cball::measurefitness(double RR, double gradient, double Vs, double K, double Range, double MS) {
//Vs fitness function variance, K carring capacity,
    const int xg = static_cast<int>(ceil(xp / 200.0));
    const int yg = static_cast<int>(ceil(yp / 200.0));
    std::vector<int> jj(3);
    if (yg - 1 < 1) {
        jj[0] = 1;
        jj[1] = 2;
        jj[2] = 5;
    } else if (yg + 1 > yrange / 200) {
        jj[0] = 4;
        jj[1] = 5;
        jj[2] = 1;
    } else {
        jj[0] = yg - 1;
        jj[1] = yg;
        jj[2] = yg + 1;
    }
    unsigned tot = 0;
    for (int i = std::max(xg - 1, 1); i <= std::min(xg + 1, xrange / 200); ++i) {
        for (unsigned j = 0; j < 3; ++j) {
            const std::vector<std::list<Cball>::iterator> grid_ij = gridindiv[i][jj[j]];
            for (unsigned k = 0; k < grid_ij.size(); ++k) {
                const double dist = distance(*grid_ij[k]);
                if (dist <= Range) {
                    //count number of individuals within the range
                    ++tot;
                }
                // count number of mailes within the ara of radius 200
                if (sexi == 0 && dist <= MS && grid_ij[k]->sexi == 1) {
                    // if this individul is female and called individual is male and if the distance between them is less than 200, then the called individual is recored as candidate males.
                    candidatemate.push_back(grid_ij[k]);
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
    std::poisson_distribution<unsigned> poisson(dfitness);
    nooffspring = poisson(engine);
}

///// serach for candidate mates
unsigned Cball::matingcount(int matingsize) const {
    double total_fitness = 0.0;
    std::vector<unsigned> indices;
    indices.reserve(candidatemate.size());
    for (unsigned i = 0; i < candidatemate.size(); ++i) {
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
    return static_cast<unsigned>(candidatemate.size());
}

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
namespace {

template <class Matrix> inline
Matrix transpose(const Matrix& A) {
    const size_t nrow = A.size();
    if (nrow == 0u) return A;
    const size_t ncol = A[0u].size();
    Matrix out(ncol, typename Matrix::value_type(nrow));
    for (unsigned row=0; row < nrow; ++row) {
        for (unsigned col=0; col < ncol; ++col) {
            out[col][row] = A[row][col];
        }
    }
    return out;
}

inline std::vector<std::vector<short>> make_haplotypes(unsigned popsize, unsigned freq) {
    std::vector<short> bits(popsize, 0);
    for (unsigned i = 0; i < freq; ++i) {
        bits[i] = 1;
    }
    std::vector<std::vector<short>> sites;
    sites.reserve(nloci_neutral + nloci);
    ////// set genes for resource use ///////
    for (unsigned i = 0; i < nloci_neutral; ++i) {
        std::shuffle(bits.begin(), bits.end(), engine);
        sites.push_back(bits);
    }
    const unsigned num_monomorphic = nloci - nopoly;
    const unsigned num_1_fixed = num_monomorphic / 2;
    const unsigned num_0_fixed = num_monomorphic - num_1_fixed;
    for (unsigned i = 0; i < num_1_fixed; ++i) {
        sites.push_back(std::vector<short>(popsize, 1));
    }
    // genes 1 and 0 at 3 loci are randomly allocated for all the individuals
    for (unsigned i = 0; i < nopoly; ++i) {
        std::shuffle(bits.begin(), bits.end(), engine);
        sites.push_back(bits);
    }
    for (unsigned i = 0; i < num_0_fixed; ++i) {
        sites.push_back(std::vector<short>(popsize, 0));
    }
    return transpose(sites);
}

}
/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

void Newball(unsigned n, double fr, std::list<Cball>* list1) {
    const unsigned num_males = n / 2;
    const unsigned aa = static_cast<unsigned>(std::round(fr * fr * n));
    const unsigned ab = static_cast<unsigned>(std::round(fr * (1 - fr) * 2 * n));
    const std::vector<std::vector<short>> haplotypes1 = make_haplotypes(n, aa + ab);
    const std::vector<std::vector<short>> haplotypes2 = make_haplotypes(n, aa);
    std::uniform_int_distribution<int> uniform_x(1, 500);
    std::uniform_int_distribution<int> uniform_y(1, yrange);
    for (unsigned i = 0; i < n; ++i) {
        const int xx = uniform_x(engine) + xrange / 2 - 250;
        const int yy = uniform_y(engine);
        list1->push_back(Cball(xx, yy, i < num_males, haplotypes1[i], haplotypes2[i]));
    }
}

// Create new indiviaul
void Newball(const char* infile, std::list<Cball>* list1) {
    //// creat n individuals and initialize position and sex
    const std::vector<std::vector<int> > matrix = read_int_array(infile);
    nloci = static_cast<unsigned>(matrix[0].size()) - 3;
    for (unsigned row = 0; row < matrix.size(); ++row) {
        list1->push_back(Cball(matrix[row]));
    }
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

// save the results as afile
void SaveF(const std::list<Cball>& clist, unsigned g, unsigned gg, unsigned n) {
    std::ostringstream oss;
    oss << "File" << gg << "-" << g;
    FILE* fp = fopen(oss.str().c_str(), "w");
    std::list<Cball>::const_iterator it = clist.begin();
    for (unsigned i = 0; i < n; ++it, ++i) {
        const int m1 = it->sexi;
        const double m2 = it->resource();
        const double m3 = (m1 == 0) ? it->nooffspring : it->dfitness;
        fprintf(fp, "%d\t %d\t %d\t %f\t  %7.3f\t", it->xp, it->yp, m1, m2, m3);
        for (unsigned ge = 0; ge < nloci_neutral + nloci; ++ge) {
            fprintf(fp, "%d\t ", it->gene1[ge]+ it->gene2[ge]);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
}

void SaveA(const std::list<Cball>& clist, unsigned g, unsigned gg, unsigned n, unsigned clas) {
    std::ostringstream oss;
    oss << "Resl" << gg << "-" << g;
    FILE* fp = fopen(oss.str().c_str(), "w");
    std::vector<double> avfit(clas);
    std::vector<int> nof(clas);
    const double w = 32000.0 / clas;
    for (unsigned x = 0; x < clas; ++x) {
        std::list<Cball>::const_iterator it = clist.begin();
        for (unsigned i = 0; i < n; ++it, ++i) {
            if(it->sexi == 0) {
                if(w * x < it->xp && it->xp <= w * (x + 1)) {
                    avfit[x] += it->nooffspring;
                    nof[x] += 1;
                }
            }
        }
    }
    for (unsigned x = 0; x < clas; ++x) {
        const int m1 = static_cast<int>(w * x);
        const int m2 = static_cast<int>(w * (x + 1));
        const double m3 = static_cast<double>(nof[x]);
        const double m4 = m3 > 0.0 ? avfit[x] / m3 : 0.0;
        fprintf(fp, "%d\t %d\t %7.4f\t %7.4f\t", m1, m2, m3, m4);
        fprintf(fp, "\n");
    }
    fclose(fp);
}

void SaveE(unsigned g, unsigned gg) {
    std::ostringstream oss;
    oss << "File" << gg << "-" << g;
    FILE* fp = fopen(oss.str().c_str(), "w");
    fprintf(fp, "%d", 0);
    fclose(fp);
}
