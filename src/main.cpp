#include <iostream>
#include <fstream>
#include <string>
#include <list>
#include <random>
#include "uball.h"
#include "global.hpp"

const unsigned nloci_neutral = 10;
int xrange;
int minxrange;
int yrange;
int xmi = 0;
int xma = 0;
double Vp;
unsigned nloci;
unsigned nopoly;
double allele_effect;
std::mt19937_64 engine;

int main(int argc, char* argv[]) {
    if (argc < 2) throw std::runtime_error("too few arguments");
    const std::string gradient(argv[1]);
    int UNUSED;
    unsigned no;
    int homeranges;
    int sizemating;
    unsigned genS;
    double G;
    double VS;
    double CC;
    double reprate;
    double mutationr;
    double nem;
    unsigned norepeat;
    unsigned nogeneration;
    unsigned noclas;
    double mdispersal;
    {
        std::string buffer;
        std::ifstream fin("Inputfile");
        fin >> buffer >> xrange       // range of habitat
            >> buffer >> minxrange
            >> buffer >> yrange
            >> buffer >> UNUSED
            >> buffer >> UNUSED
            >> buffer >> UNUSED
            >> buffer >> UNUSED;
        fin >> buffer;
        if (buffer == "flatmin") {
            fin           >> xmi
                >> buffer >> xma
                >> buffer;
        } else if (gradient == "flat") {
            throw std::runtime_error("Inputfile lacks flatmin/flatmax");
        }
        fin           >> no           // initial number of individuals
            >> buffer >> Vp           // phenotypic variance
            >> buffer >> homeranges   // area for counting densities
            >> buffer >> sizemating   // area for searching mates
            >> buffer >> genS         // Save files every genS generation
            >> buffer >> G            // slope of environment
            >> buffer >> VS           // variation of fitness function
            >> buffer >> CC           // Carring capacity
            >> buffer >> reprate      // reproductive rate
            >> buffer >> mutationr    // mutation rate
            >> buffer >> nem
            >> buffer >> nloci
            >> buffer >> nopoly
            >> buffer >> allele_effect
            >> buffer >> norepeat     // No of reprilcate
            >> buffer >> nogeneration // No of generations for one simulation
            >> buffer >> noclas
            >> buffer >> mdispersal;  // dispersal distance of males
    }
    ///////////////////////////////////////
    if (gradient == "linear") {
        xmi = -1;// global flag for measurefitness()
    } else if (gradient != "flat" && gradient != "steep") {
        throw std::runtime_error("invalid gradient choice: " + gradient);
    }
    const double fdispersal = mdispersal;
    std::random_device seeder;
    engine.seed(seeder());
    for (unsigned prep = 1; prep <= 1; ++prep) {
        std::cout << "no of loci x effect x 2 =" << nloci*allele_effect*2 << std::endl;
        std::cout << "gradient x xmax(=x range)=" << G*xrange << std::endl;
        if (xmi > 0 || xma > 0) {// not 2008BK
            if (nloci * allele_effect * 2 - G * xrange > 1e-6) norepeat = 0;
        }
        for (unsigned gggg = 1; gggg <= norepeat; ++gggg) {
            std::list<Cball> alist;// creat list for control individuals
            // Creat individuals
            if (argc > 2) {
                Newball(argv[2], &alist);
            } else {
                Newball(no, 0.5, &alist);
            }
            int minx = 0;
            int maxx = 0;
            for (unsigned y = 1; y <= nogeneration; ++y) {// Generation
                const size_t itemP = alist.size();// check number of individuals
                if (itemP == 0) {
                    SaveE(prep, gggg);
                    break;
                }
                if (y == 1 || y % 2 == 0) {
                    std::cout << "Gener = " << y << ":Repeat = " << gggg << "\tNo= " << itemP <<"\tMin(x)= " << minx << "\tMax(x)= " << maxx << std::endl;
                    // print no of generations and no of individuals
                }
                AssignBucket(alist);
                // Measure fitness for each individual
                for (std::list<Cball>::iterator it = alist.begin(); it != alist.end(); ++it) {
                    it->measurefitness(reprate, G, VS, CC, homeranges, sizemating);
                }
                // reproduced by females /////
                maxx = 0;
                minx = xrange;
                const unsigned num_parents = static_cast<unsigned>(alist.size());
                for (std::list<Cball>::iterator it = alist.begin(); it != alist.end(); ++it) {
                    if (it->xp > maxx) maxx = it->xp;
                    if (it->xp < minx) minx = it->xp;
                    if (it->sexi == 0 && it->nooffspring > 0) {
                        // choose female having nonzero fitness
                        // female search candidate mates
                        const unsigned nofm = it->matingcount(sizemating);
                        if (nofm < it->candidatemate.size()) {
                            const Cball& male = *it->candidatemate[nofm];
                            // give birth to young
                            it->nreproduction(male, &alist, mdispersal, fdispersal, mutationr, nem);
                        }
                    }
                }
                if (y % genS == 0 || y == 1) {
                    SaveF(alist, y, gggg, num_parents);
                    SaveA(alist, y, gggg, num_parents, noclas);
                }
                std::list<Cball>::iterator it = alist.begin();
                for (unsigned i = 0; i < num_parents; ++i) {
                    it = alist.erase(it);// clear parents
                }
            }// y no generation
        }
    }
    std::cout << "x=" << xrange
              << ";y=" << yrange
              << ";n=" << no
              << ";Vp=" << Vp
              << ";HR=" << homeranges
              << ";mdis=" << mdispersal
              << ";fdis=" << fdispersal
              << ";mating=" << sizemating
              << ";G=" << G
              << ";Vs=" << VS
              << ";K=" << CC
              << ";r="<< reprate
              << ";mr=" << mutationr
              << ";ng=" << nloci_neutral + nloci
              << "\n";
}
