#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include "ucrandom.h"
#include "uball.h"
#include "read_array.hpp"

extern short *Fresource;
extern short nogene,rn,doend;
short xrange,yrange,minxrange, generation, kfood,homeranges;
short xmi = 0, xma = 0;
long no;
extern double Vp;
extern gene genen;
extern double mutationrate,benefit,cost;
extern std::list<Cball>::iterator mp;
extern Cball *indiv;
extern double mdispersal;
extern short nloci, nopoly;
extern double allele_effect;
extern double reso[21];
short noindiv[165][7];
std::list<Cball>::iterator gridindiv[161][6][301];


//extern short resoN[21],KresoN[21];
void Cball::Iball(short x, short y, short s)// initialize position, sex, and fitness;
{

	sexi=s;
	xp=x;
	ix=x;
	yp=y;
	iy=y;
	mx=0;
	my=0;
fitness=0;
	nomating=0;
	nocandiate=0;
}

void Cball::Igene(short nthgene, short gf, short gl)// initialize genes
{
     gene1[nthgene]=gf;
     gene2[nthgene]=gl;
     genotype[nthgene]=gf+gl;
}

double Cball::ResourceM()// calculate resouce use phenotypes from the genotype
{
  short i;
  double sum;
  sum=0;
  for(i=11;i<=10+nloci;i++)
  {
   sum += allele_effect*(gene1[i]+gene2[i]);
   }
   sum +=nrnd(sqrt(Vp));
   return sum;
 }

/*

short Cball::choice()// calculate female choice phenotypes from the genotype
{
  short i,sum;
  sum=0;
  for(i=11;i<=20;i++)
  {
   sum += gene1[i]+gene2[i];
   }
   return sum;
 }
short Cball::choiceM()// calculate male sexual phenotypes from the genotype
{

  short i,sum;
  sum=0;
  for(i=21;i<=30;i++)
  {
   sum += gene1[i]+gene2[i];
   }
   return sum;
 }
short Cball::select()// calculate selectivity phenotypes from the genotype
{
  short i,sum;
  sum=0;
  for(i=31;i<=40;i++)
  {
   sum += gene1[i]+gene2[i];
   }
   return sum;
 }
*/

// reproduction
void Cball::nreproduction (std::list<Cball> *ablist,  short nogene, double mdis, double fdis, double mr, double nmr)
		{
		double didi, sddispersal,dist,dist2;
		long  j,i, nx, ny,  d,k,x,y,xx,yy,y2;
		short nooffspring,gg, og1[200], og2[200];

		double mu,mrr;
	nomating++;
	mateP=mp->ResourceM();
mx=mp->xp;
my=mp->yp;

x=xp;// position x for focal female
y=yp;// position y for focal female
xx=mp->xp;
yy=mp->yp;
if (y >= mp->yp)y2=yrange + mp->yp;else y2=mp->yp - yrange;
 dist=sqrt((x - xx)*(x - xx)+(y - yy)*(y - yy));
 dist2=sqrt((x-xx)*(x-xx)+(y-y2)*(y-y2));
 if(dist > dist2) dist=dist2;
mdistance=dist;
	//mdistance=sqrt( (xp- mp->xp)*(xp-mp->xp)+(yp- mp->yp)*(yp- mp->yp));









//	mx=mp->xp;
//	my=mp->yp;
//	if(mdistance > 10000 || mdistance < 0)
//	{
//	k++;


//	}

	nooffspring = fitness;

		if (nooffspring > 0 )
			{
				for (j = 1; j<= nooffspring;j++)
					{
						for(k=1;k<=nogene;k++)
						{

						///////  inheritance genes from mother and father
						gg = randombit();// random number 1 or 0
						if ( gg ==1 )
							{ og1[k] = gene1[k];}
						else
							{og1[k] = gene2[k];}
						gg = randombit();
						if ( gg==1 )
							{ og2[k] = mp->gene1[k];}
						else
						{og2[k] = mp->gene2[k];}
						///////// mutation ////////
						if(k>10)mrr=mr;else mrr=nmr;
						mu=longurnd();

						if (mrr >= mu)
							{
							if( og1[k]==1)
								{
								og1[k]=0;
								}
								else
								{
								og1[k]=1;
								}
							}
						mu=longurnd();
							if (mrr >= mu)
							{
							if( og2[k]==1)
								{
								og2[k]=0;
								}
								else
								{
								og2[k]=1;
								}
							}

						}

					gg = randombit();// determin offspring sex
					if (gg==0)
					{sddispersal=fdis;}
					if(gg==1)
					{sddispersal=mdis;}


							do{//// desersal
							didi = nrnd(sddispersal); // random number from normal distribtion with sddispersal standard deviationa and 0 mean

							d=abs(didi);
							randomove(ix, iy, d, &nx, &ny);
							if(ny > yrange) ny=ny - yrange;
							if(ny <=0) ny=yrange+ny;
							}while((nx <= minxrange) || (nx > xrange )   ||  (ny <= 0) || (ny > yrange ) );


						//indiv=new Cball;
						indiv->Iball(nx, ny,gg);// Initialize offspring position and sex
						for (i=1;i<=nogene;i++)
							{
							indiv->Igene(i, og1[i], og2[i]);// offspring genes
							}
							ablist->push_back(*indiv);// add offspring the list

				}
			}

	}


// calculate fitness
void Cball::measurefitness(std::list<Cball>  *clist,   double RR, double gradient,double Vs,double K,double Range,double MS)
{
//Vs fitness function variance, K carring capacity,
short  AA, AB, BB,tot,y2,t1,b1,xg,yg,ff,r,t,b,i,j,n,xx,yy,jj[4],sfitness;
double dist,Sx, dfit,sfit,dist2;
tot=0;
AA=0;
AB=0;
BB=0;

std::list<Cball>::iterator indiv;
t1=0;
b1=0;
xg=(short)ceil((double)(xp)/200.);
yg=(short)ceil((double)(yp)/200.);
ff=xg-1;
r=xg+1;
t=yg-1;
b=yg+1;
if(ff<1){ff=1;}
if(r> xrange/200){r=xrange/200;}

if(t<1){
	jj[1]=1;
	jj[2]=2;
	jj[3]=5;
	}else
if(b>yrange/200)
        {
	jj[1]=4;
	jj[2]=5;
	jj[3]=1;
}else
{
jj[1]=yg-1;
jj[2]=yg;
jj[3]=yg+1;
}





tot=0;
for(i=ff;i<=r;i++)
		{
	for(j=1;j<=3;j++)
		{

			for(n=1;n<=noindiv[i][jj[j]];n++)
			{
			xx=gridindiv[i][jj[j]][n]->xp;
			yy=gridindiv[i][jj[j]][n]->yp;
			if (iy >= yy)y2=yrange +yy;else y2=yy - yrange;
				dist=sqrt((ix-xx)*(ix-xx)+(iy-yy)*(iy-yy));
				dist2=sqrt((ix-xx)*(ix-xx)+(iy-y2)*(iy-y2));
				if(dist > dist2) dist=dist2;
					if (dist <= Range )
					{
					tot++;//count number of individuals within the range
					}
					if (sexi==0 && dist <= MS && gridindiv[i][jj[j]][n]->sexi==1)// count number of mailes within the ara of radius 200
					{
			 		nocandiate++;// if this individul is female and called individual is male and if the distance between them is less than 200, then the called individual is recored as candidate males.
			 		candidatemate[nocandiate]=gridindiv[i][jj[j]][n];// store the canditate males
					}
			}



		 }
		}

// count number of individual within the range
/*for(indiv=clist->begin();indiv !=clist->end();indiv++)
			{
			if (iy >= indiv->yp)y2=yrange + indiv->yp;else y2=indiv->yp - yrange;
			 dist=sqrt((ix-indiv->xp)*(ix-indiv->xp)+(iy-indiv->yp)*(iy-indiv->yp));
			 dist2=sqrt((ix-indiv->xp)*(ix-indiv->xp)+(iy-y2)*(iy-y2));
			 if(dist > dist2) dist=dist2;
				if (dist <= Range )tot++;
			 }

*/

	AA=xrange/2-xmi;
	BB=xma-xrange/2;


//// calculate fitness

if(ix >= xrange/2 - AA && ix <= xrange/2 + BB )Sx=16+gradient*(xrange/2-4000);
if(ix < xrange/2 - AA)Sx=16+gradient*(ix-4000+AA);
if(ix > xrange/2 + BB)Sx=16+gradient*(ix-4000-BB);



dfitness=2+RR*(1-tot/K)-(Sx-ResourceM())*(Sx-ResourceM())/(2*Vs);
if(dfitness <=0) dfitness=0;


fitness=prnd(dfitness);
//if(urnd() <= dfit) fitness=sfit+1;else fitness=sfit;// change double fit to integer fitness
//if(dfitness <=0) dfitness=0;


}//endomethod;
//*******************************************************************

      // Functions


//*******************************************************************


// Create new indiviaul

void Newball(short n, short male, std::list<Cball>  *list1, double fr)
{
	long i,j, aa, bb,ab;

	std::list<Cball>::iterator individual;


	//// creat n individuals and initialize position and sex
	std::vector<std::vector<int> > vecvecint;
	vecvecint = read_int_array("TestInput.txt");
	n=vecvecint.size();
	nloci=	vecvecint[0].size()-3;
	for (size_t row=0; row<vecvecint.size(); ++row) {
		indiv->Iball(vecvecint[row][0],vecvecint[row][1],vecvecint[row][2]);
		for (size_t col=3; col<vecvecint[row].size(); ++col)
		{
			if(vecvecint[row][col]==0)indiv->Igene(col+8,0, 0);
			if(vecvecint[row][col]==2)indiv->Igene(col+8,1, 1);
			if(vecvecint[row][col]==1){
				if(randombit() ==0)
				{
					indiv->Igene(col+8,1, 0); 												}
				else
				{
					indiv->Igene(col+8,0, 1);
				}
			}
		}
		list1->push_back(*indiv);
	}

	//////////// set genes for  resource use ///////

	for(individual=list1->begin();individual!=list1->end();individual++)
	{
		for( j=1;j<=10;j++)
		{
			if(randombit()==1)
			{aa=1;}
			else{
				aa=0;
			}
			if(randombit()==1)
			{bb=1;}
			else{
				bb=0;
			}
			individual->Igene(j,aa,bb);
		}
	}

	/*
for( j=1;j<=10;j++)
		{
		 GerateRandomperm (n, tur);
		 for(i=1;i<=aa;i++){ ge1[tur[i]]=1;ge2[tur[i]]=1;}
		 for(i=aa+1;i<=aa+ab;i++){ ge1[tur[i]]=1;ge2[tur[i]]=0;}
		 for(i=aa+ab+1;i<=n;i++){ ge1[tur[i]]=0;ge2[tur[i]]=0;}
		 for(individual=list1->begin(),i=1;individual!=list1->end();individual++,i++)
		individual->Igene(j,ge1[i],ge2[i]);
		}

for( j=11;j<=10+(nloci-nopoly)/2;j++)
		{
		for(individual=list1->begin();individual!=list1->end();individual++)
			individual->Igene(j,1, 1);
		}
for( j=11+(nloci-nopoly)/2;j<=10+(nloci-nopoly)/2+nopoly;j++)/// genes 1 and 0 at 3  loci  are randomly allocated for all the  individuals
		{
		 GerateRandomperm (n, tur);
		 for(i=1;i<=aa;i++){ ge1[tur[i]]=1;ge2[tur[i]]=1;}
		 for(i=aa+1;i<=aa+ab;i++){ ge1[tur[i]]=1;ge2[tur[i]]=0;}
		 for(i=aa+ab+1;i<=n;i++){ ge1[tur[i]]=0;ge2[tur[i]]=0;}
		 for(individual=list1->begin(),i=1;individual!=list1->end();individual++,i++)
			individual->Igene(j,ge1[i],ge2[i]);
		}
for( j=11+(nloci-nopoly)/2+nopoly;j<=10+nloci;j++)
		{
		for(individual=list1->begin();individual!=list1->end();individual++)
			individual->Igene(j,0, 0);
		}

	delete tur;
	delete ge1;
	delete ge2;
	 */


}



///// serach for candidate mates

void matingcount (std::list<Cball>  *clist, std::list<Cball>::iterator focalindiv,std::list<Cball>::iterator *matp, long noi,  short matingsize, short *dens)
	{
			long  i, k,xx, yy, x, y,NM,y2;
			double  dist,sum,r,dist2;
			double total_fitness=0;
			std::list<Cball>::iterator indiv, mk;
		x=focalindiv->xp;// position x for focal female
		y=focalindiv->yp;// position y for focal female
		NM=focalindiv->nocandiate;
		*dens = 0;
		k = 0;
		for(i=1;i <=NM;i++)
					{

						if (focalindiv->candidatemate[i]->fitness >0 )
							{


								xx=focalindiv->candidatemate[i]->xp;
								yy=focalindiv->candidatemate[i]->yp;

							      if (y >= yy)y2=yrange + yy;else y2=yy - yrange;

								 dist=sqrt((x - xx)*(x - xx)+(y - yy)*(y - yy));
								 dist2=sqrt((x-xx)*(x-xx)+(y-y2)*(y-y2));
								 if(dist > dist2) dist=dist2;




								if(dist <=matingsize  )
									{

									/// count and save candidate male
									 k++;
									 total_fitness +=focalindiv->candidatemate[i]->dfitness;
									// total_fitness +=indiv->fitness;
									 matp[k]=focalindiv->candidatemate[i];
									 *dens +=1;
				   					 }
							}

					}

		if(*dens > 0)
		{
		sum=0;i=1;
		r=urnd()*total_fitness;
			do{
			mk=matp[i];
			//sum +=mk->fitness;
			sum +=mk->dfitness;
			i++;
			}while(sum < r && i <= k);
		*dens=i-1;
		}
	mk=matp[*dens];


	}


// save the results as afile


void SaveF(std::list<Cball>  *clist,short g,short gg,double md,long n )
{
long  m1,i,ge;
double m2;
char ab[12]="File%%%", rch[10]="0000", nots[10]="0000";
std::list<Cball>::iterator indiv;
double m3;

FILE *fp;
sprintf(rch, "%d", g);
sprintf(nots, "%d", gg);
if(gg >= 10)
{
ab[7]=rch[0];
ab[8]=rch[1];
ab[9]=rch[2];
ab[10]=rch[3];
ab[11]=rch[4];
ab[6]=0x2D;
ab[4]=nots[0];
ab[5]=nots[1];
}
if(gg <10)
{
ab[6]=rch[0];
ab[7]=rch[1];
ab[8]=rch[2];
ab[9]=rch[3];
ab[10]=rch[4];
ab[5]=0x2D;
ab[4]=nots[0];
//ab[5]=not[1];
}

fp=fopen(ab,"w");
// fprintf(fp, "%7.3f", md);
//  fprintf(fp, "\n");
for(indiv=clist->begin(), i=1;i<=n;indiv++,i++)
	{
   m1=indiv->sexi;
   m2=indiv->ResourceM();
  if(m1==0) m3=indiv->fitness; else m3=indiv->dfitness;
    fprintf(fp, "%d\t %d\t %d\t %f\t  %7.3f\t", indiv->xp, indiv->yp,m1, m2, m3);

   for(ge=1;ge<=nogene;ge++)
        fprintf(fp, "%d\t ", indiv->gene1[ge]+ indiv->gene2[ge]);

 	 fprintf(fp, "\n");
  }

  fclose(fp);
 }



void SaveE(short g,short gg)
{
short  m1;
char ab[12]="File%%%", rch[10]="0000", nots[10]="0000";


FILE *fp;
sprintf(rch, "%d", g);
sprintf(nots, "%d", gg);
if(gg >= 10)
{
ab[7]=rch[0];
ab[8]=rch[1];
ab[9]=rch[2];
ab[10]=rch[3];
ab[11]=rch[4];
ab[6]=0x2D;
ab[4]=nots[0];
ab[5]=nots[1];
}
if(gg <10)
{
ab[6]=rch[0];
ab[7]=rch[1];
ab[8]=rch[2];
ab[9]=rch[3];
ab[10]=rch[4];
ab[5]=0x2D;
ab[4]=nots[0];
//ab[5]=not[1];
}

fp=fopen(ab,"w");

	{
   m1=0;

    fprintf(fp, "%d",m1);
 	// fprintf(fp, "\n");
  }

  fclose(fp);
 }


void SaveA(std::list<Cball>  *clist,short g,short gg,double md, long n, short clas )
{
long  m1,m2,i,x;
char ab[12]="Resl%%%", rch[10]="0000", nots[10]="0000";
std::list<Cball>::iterator indiv;
double m3,m4,w;
double avfit[322];
short nof[322];
for(x=0;x<=321;x++)
  {
    avfit[x]=0;
    nof[x]=0;
    }

FILE *fp;
sprintf(rch, "%d", g);
sprintf(nots, "%d", gg);
if(gg >= 10)
{
ab[7]=rch[0];
ab[8]=rch[1];
ab[9]=rch[2];
ab[10]=rch[3];
ab[11]=rch[4];
ab[6]=0x2D;
ab[4]=nots[0];
ab[5]=nots[1];
}
if(gg <10)
{
ab[6]=rch[0];
ab[7]=rch[1];
ab[8]=rch[2];
ab[9]=rch[3];
ab[10]=rch[4];
ab[5]=0x2D;
ab[4]=nots[0];
//ab[5]=not[1];
}

w=32000/(double) clas;


fp=fopen(ab,"w");
//  fprintf(fp, "%7.3f", md);
// fprintf(fp, "\n");
for(x=1; x<=clas;x++)
{
	for(indiv=clist->begin(), i=1;i<=n;indiv++,i++)
	{
	if(indiv->sexi == 0)
	{
	 if( indiv->xp > w*(x-1)  && indiv->xp <= w*x)
	   {
	      avfit[x] += indiv-> fitness;
	      nof[x] +=1;
	    }
	}
      }
   }



for(x=1; x<=clas;x++)
	{
   	m1=w*(x-1);
  	 m2=w*x;
   	m3=nof[x];
  	 if (nof[x]==0) m4=0; else m4= avfit[x]/(double)nof[x];
	   fprintf(fp, "%d\t %d\t %7.4f\t %7.4f\t",m1, m2, m3, m4);
	  fprintf(fp, "\n");
  }

  fclose(fp);
 }

void AssignBucket(std::list<Cball>  *list1)
{
short nogx,nogy, xc,yc;
std::list<Cball>::iterator indiv;
//std::list<Cball>::iterator *mk;

nogx=xrange/200;//the max number of x dimention of the bucket girds y=3200  x=160
nogy=yrange/200;//the max number of y dimention of the bucket girds　ｙ=1000 y=5
for(xc=0;xc<=nogx;xc++)
 for(yc=0;yc<=nogy;yc++)
		noindiv[xc][yc]=0;//Initialize the array of the number of individuals in a bucket grid xc yc

	for(indiv=list1->begin();indiv !=list1->end();indiv++)
		{
  			xc=(short)ceil(indiv->xp/200.);yc=(short)ceil(indiv->yp/200.);
  			noindiv[xc][yc]+=1;// count the number of individual in a bucket
  			gridindiv[xc][yc][noindiv[xc][yc]]=indiv;//store the individuals into the array
			}





}
