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
using namespace std;  //introduces namespace std
#include "uball.h"
#include <unistd.h>

//extern short nogene;
/*** v1 to v20 are double variables from dialog box ************/
double v1,v2,v3,v4,v5,v6,v7,v8,v9,v10;
double v11,v12,v13,v14,v15,v16,v17,v18,v19,v20;
/**********************************************/
extern short xrange,yrange, minxrange, generation, homeranges;
extern long no;
double mutationrate,Vp,CC,revar[51],revar2[51];
short sizemating;
gene genen;
double reprate;

list<Cball>::iterator mp;
list<Cball> *alist;
double mdispersal, fdispersal,mutationr;
extern list<Cball>::iterator gridindiv[162][7][501];
short nogene;
Cball *indiv;


short nloci, nopoly;
 double allele_effect;


int  main()
{
char a[100];

list<Cball>::iterator individual;
list<Cball>::iterator *matp;
alist=new list<Cball>;// creat list for control individuals
matp=new list<Cball>::iterator[30000+1];// creat array for candidate mates
indiv=new Cball;// creat new individual




double G, VS,mdis,mdis1,nem;
short nogeneration,gggg,norepeat,prep;
long item,itemP,nn;
long x,y, no,i, maxx, minx,noclas;
short nofm,nomale,genS;
short gfff;
short leftp, rightp,Lphenotype,Rphenotype,LN, RN;
char name_dir[64];
unsigned long ti;
ti=time(NULL);
getcwd(name_dir,64);// get current directly
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
ifstream fin("Inputfile");

fin
 >> a >> xrange
 >> a >> minxrange
 >> a >> yrange
 >> a >> Lphenotype
 >> a >> Rphenotype
 >> a >> LN
 >> a >> RN
 >> a >> no
 >> a >> Vp
  >>  a >> homeranges
  >>   a >> sizemating
  >> a >> genS
  >>  a >> G
  >> a >> VS
  >> a >> CC
  >>   a >> reprate
  >> a >> mutationr
  >> a >> nem
   >> a >> nloci
   >> a >> nopoly
   >> a >> allele_effect
   >> a >> norepeat
   >> a   >>nogeneration
  >> a >> noclas
>> a >> mdispersal;


fin.close();
nomale=(short)no*0.5; // initial number of males
///////////////////////////////////////
nogene=nloci+10;
numberofgenes=nogene;
 fdispersal=mdispersal;
gfff=2;
for(prep=1;prep<=1;prep++)// repeat
{
//fdispersal=revar[prep];
//mdispersal=revar[prep];
cout<< "no of loci x effect x 2 =" << nloci*allele_effect*2 << endl;
cout << "gradient x xmax(=x range)=" << G*xrange << endl;
if(nloci*allele_effect*2 !=  G*xrange)norepeat=0;
for(gggg=1;gggg<=norepeat;gggg++)// repeat
{
randomizec();// initialize random number  (short type)
lrandomizec();// initialize random number (long type


Newball(no, nomale, alist,  0.5);// Creat individuals



for (y=1;y<=nogeneration;y++)// Generation
{


AssignBucket(alist);



itemP=alist->size();// check number of individuals
if(itemP<=0)
 {
 SaveE(prep,gggg);
 y=nogeneration;
 }

if (y==1 )
   cout << "Gener = " << y << ":Repeat = " << gggg << "No= " << itemP <<  endl;

 if ( y % 2 ==0)
	 {
  	 cout << "Gener = " << y << ":Repeat = " << gggg << "No= " << itemP << "Min(x)=" << minx << "Max(x)=" << maxx << endl;
	// print no of generations and no of individuals
	leftp=0;
	rightp=0;
//	for(individual=alist->begin();individual!=alist->end();individual++)// Measure fitness for each individual
//		{
 //   	 		if(individual->ResourceM() == Lphenotype)leftp++;
  //  	 		if(individual->ResourceM() == Rphenotype)rightp++;
// 		}
 //	if ( leftp > LN  && rightp >RN)
 //		{
 //
 //		y=nogeneration;
 //		}
	}


if ( itemP > 0 )// itemP= number of individuals
{





  if (y % 10){ randomizec();// initialize random number  (short type)
lrandomizec();}
 	for(individual=alist->begin();individual!=alist->end();individual++)// Measure fitness for each individual
	{
     individual->measurefitness( alist, reprate,  G, VS, CC,homeranges,sizemating);
 	}


 	// reproduced by females /////
 	maxx=0; minx=40000;
    itemP=alist->size();
    individual=alist->begin();
 	for(i=1;i<=itemP;i++,individual++)
	{
    if(individual->xp > maxx)maxx=individual->xp;
    if(individual->xp <minx)minx=individual->xp;

 	if(individual->sexi==0 && individual->fitness > 0)// choose female having nonzero fitness
 		{
 		// Search candidate mates :
 		matingcount(alist, individual,matp,itemP, sizemating, &nofm);// female search candidate mates

 		if(nofm != 0)// if candidate males were not zero
 		{

 			mp=matp[nofm];	// mp is a chosen male
 			individual->nreproduction (alist,  nogene,  mdispersal, fdispersal,mutationr,nem);// give birth to young

 			}
 		}

 	 }
 	 item=alist->size();
 		if(y % genS ==0 )
 		{
 		SaveF( alist, y,gggg,mdis,itemP);
 		SaveA( alist,y,gggg,mdis, itemP,noclas);
 		}

 	mdis=0;
 	mdis1=0;
 	nn=0;
	// time=total number of individual = paranet + offspring
    individual=alist->begin();
 	for(x=1;x<=itemP;x++)
	{
 	     if (individual->sexi==0  && individual->fitness >0)
    		  {
 //   	    mdis1=individual->mdistance;
  	      mdis +=individual->mdistance*individual->mdistance;
 	      nn++;
	       }

	individual=alist->erase(individual);// crear parents
 	}

// 	mdis =sqrt(mdis/(double)itemP);
 	mdis =sqrt(mdis/(double)nn);

 }



 }// y no generation



 //--------------------------------
alist->clear();// clear all individuals
}

}



cout << "x=" << xrange <<  ";y=" << yrange << ";n=" << no << ";Vp=" << Vp <<  ";HR=" << homeranges << ";mdis=" << mdispersal <<
";fdis=" << fdispersal <<  ";mating=" << sizemating <<   ";G=" << G << ";Vs=" << VS << ";K=" <<
CC <<  ";r="<< reprate << ";mr=" << mutationr << ";ng=" << nogene ;


}


