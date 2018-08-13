#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <ctime>
#include <iostream>
#include <fstream>
#include <list>
#include "uglobal.h"
#include "ucrandom.h"
#include "uball.h"
#include <unistd.h>

/*** v1 to v20 are double variables from dialog box ************/
double v1,v2,v3,v4,v5,v6,v7,v8,v9,v10;
double v11,v12,v13,v14,v15,v16,v17,v18,v19,v20;
/**********************************************/
extern short xrange, yrange, minxrange, generation, homeranges, xmi, xma;
double mutationrate, Vp, CC, revar[51], revar2[51];
short sizemating;
gene genen;
double reprate;

std::list<Cball>::iterator mp;
std::list<Cball> *alist;
double mdispersal, fdispersal, mutationr;
extern std::list<Cball>::iterator gridindiv[162][7][501];
short nogene;
Cball *indiv;
extern short noindiv[165][7];
short nloci, nopoly;
double allele_effect;


int main() {
    std::list<Cball>::iterator individual;
    std::list<Cball>::iterator *matp;
    alist = new std::list<Cball>;// creat list for control individuals
    matp = new std::list<Cball>::iterator[30000+1];// creat array for candidate mates
    indiv = new Cball;// creat new individual
    double G, VS, mdis, mdis1, nem;
    short nogeneration, gggg, norepeat, prep;
    long item, itemP, nn;
    long x, y, no, i, maxx, minx, noclas;
    short nofm, nomale, genS;
    short gfff;
    short leftp, rightp, Lphenotype, Rphenotype, LN, RN;
    char name_dir[64];
    getcwd(name_dir, 64);// get current directly
    chdir(name_dir);// change to current directly
    short numberofgenes;
    randomizec();
    lrandomizec();
    // //////////////Input parameter values
    /*
    xrange=5000;
    yrange=1000; // range of habitat
    no=500;     // initial number of individuals
    Vp=0;// phenotypic variance
    homeranges=50;// area for counting densities
    mdispersal=50;// dispersal distance of males
    fdispersal=50;// dispersal CompetitionRange=50;
    sizemating=70;// area for searching mates
    genS=10;    // Save files every genS generation.
    G=0.004;// slope of environment
    VS=2;// variation of fitness function
    CC=7; // Carring capacity
    reprate=1.6; // reproductive rate
    mutationr=0.0001;//mutation rate
    nogene=10;       // No of loci
    norepeat=1;      // No of reprilcate
    nogeneration=100;// No of generations for one simulation 2000;*/
    /////////////////File Input ////////////
    {
        std::string buffer;
        std::ifstream fin("Inputfile");
        fin >> buffer >> xrange
            >> buffer >> minxrange
            >> buffer >> yrange
            >> buffer >> Lphenotype
            >> buffer >> Rphenotype
            >> buffer >> LN
            >> buffer >> RN;
        fin >> buffer;
        if (buffer == "flatmin") {
            fin           >> xmi
                >> buffer >> xma
                >> buffer;
        }
        fin           >> no
            >> buffer >> Vp
            >> buffer >> homeranges
            >> buffer >> sizemating
            >> buffer >> genS
            >> buffer >> G
            >> buffer >> VS
            >> buffer >> CC
            >> buffer >> reprate
            >> buffer >> mutationr
            >> buffer >> nem
            >> buffer >> nloci
            >> buffer >> nopoly
            >> buffer >> allele_effect
            >> buffer >> norepeat
            >> buffer >> nogeneration
            >> buffer >> noclas
            >> buffer >> mdispersal;
    }
    ///////////////////////////////////////
    nomale = static_cast<short>(no / 2); // initial number of males
    nogene = nloci + 10;
    numberofgenes = nogene;
    fdispersal = mdispersal;
    gfff = 2;
    for (prep = 1; prep <= 1; prep++) {
        std::cout << "no of loci x effect x 2 =" << nloci*allele_effect*2 << std::endl;
        std::cout << "gradient x xmax(=x range)=" << G*xrange << std::endl;
        if (xmi > 0 || xma > 0) {// not 2008BK
          if (nloci*allele_effect*2 - G*xrange > 1e-6) norepeat = 0;
        }
        for (gggg = 1; gggg <= norepeat; gggg++) {
            randomizec();// initialize random number (short type)
            lrandomizec();// initialize random number (long type)
            // Creat individuals
            if (xmi > 0 || xma > 0) {
                Newball(no, alist);
            } else {
                Newball2008(no, nomale, alist, 0.5);
            }
            for (y = 1; y <= nogeneration; y++) {// Generation
                AssignBucket(alist);
                itemP = alist->size();// check number of individuals
                if (itemP <= 0) {
                    SaveE(prep, gggg);
                    y = nogeneration;
                }
                if (y == 1) {
                    std::cout << "Gener = " << y << ":Repeat = " << gggg << "\tNo= " << itemP << std::endl;
                }
                if (y % 2 == 0) {
                    std::cout << "Gener = " << y << ":Repeat = " << gggg << "\tNo= " << itemP <<"\tMin(x)= " << minx << "\tMax(x)= " << maxx << std::endl;
                    // print no of generations and no of individuals
                    leftp = 0;
                    rightp = 0;
                }
                if (itemP > 0) {// itemP= number of individuals
                    if (y % 10){
                        randomizec();// initialize random number (short type)
                        lrandomizec();
                    }
                    // Measure fitness for each individual
                    for (individual = alist->begin(); individual != alist->end(); individual++) {
                        individual->measurefitness(reprate, G, VS, CC, homeranges, sizemating);
                    }
                    // reproduced by females /////
                    maxx = 0;
                    minx = 40000;
                    itemP = alist->size();
                    individual = alist->begin();
                    for (i = 1; i<= itemP; i++, individual++) {
                        if (individual->xp > maxx) maxx = individual->xp;
                        if (individual->xp < minx) minx = individual->xp;
                        if (individual->sexi == 0 && individual->fitness > 0) {
                            // choose female having nonzero fitness
                            // Search candidate mates :
                            // female search candidate mates
                            matingcount(individual, matp, sizemating, &nofm);
                            if (nofm != 0) {
                                // if candidate males were not zero
                                mp = matp[nofm];
                                // mp is a chosen male
                                individual->nreproduction (alist, nogene, mdispersal, fdispersal,mutationr,nem);// give birth to young
                            }
                        }
                    }
                    item = alist->size();
                    if (y % genS == 0 || y == 1) {
                        SaveF(alist, y, gggg, itemP, nogene);
                        SaveA(alist, y, gggg, itemP, noclas);
                    }
                    mdis = 0;
                    mdis1 = 0;
                    nn = 0;
                    individual = alist->begin();
                    for (x = 1; x <= itemP; x++) {
                        if (individual->sexi == 0 && individual->fitness > 0) {
                            mdis += individual->mdistance*individual->mdistance;
                            nn++;
                        }
                        individual = alist->erase(individual);// crear parents
                    }
                    mdis = sqrt(mdis / (double)nn);
                }
            }// y no generation
            alist->clear();// clear all individuals
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
              << ";ng=" << nogene ;
}
