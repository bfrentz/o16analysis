#include <iostream>
#include <fstream>
#include <cstdio>
#include <math.h>
using namespace::std;

#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TStyle.h"

#include "evt_config.h"
#include "hit1.1.h"
#include "signal1.1.h"

/**************************************
* Utility Functions for Data Analysis
* Updated 24 Feb. 2016
* Wanpeng Tan (wtan@nd.edu)
* Armen Gyurjinyan (agyurjin@nd.edu)
* Ethan Sauer (esauer2@nd.edu)
* Bryce Frentz (bfrentz@nd.edu)
***************************************/

double a = 3.727379*pow(10,9);//6.64424*pow(10,-27);
double c = 1.117793*pow(10,10);//1.211255*pow(10,10);<-That's 13C //2.1592577*pow(10,-26);
double B = 25*1000000;
//double v1i = sqrt( 2*( B*1000*( 1.60218*pow(10,-19) ) )/a  );
//double t = 0;
//double p = 0;
//double F = 0;
double pi = 3.1415926;

double maxE = B;
double minE = 0;


/*
double cellXLen = 2;
double cellYLen = 2;
double startHeight = -(32*cellYLen)/2;

int iXPos(int dNum){
        if(dNum == 1){
                return 400;
                }
        if(dNum == 2){
                return 500;
                }
        if(dNum == 3){
               return 500;
                }
        if(dNum == 4){
               return 400;
          }
}

int iYPos(int dNum){
        if(dNum == 1){
                return 400;
                }
        if(dNum == 2){
                return 200;
                }
        if(dNum == 3){
                return -200;
                }
        if(dNum == 4){
                return -400;
                }
}

float detecAng(int dNum){
        if(dNum == 1){
                return 45*pi/180;
                }
        if(dNum == 2){
                return 85*pi/180;
                }
        if(dNum == 3){
                return -85*pi/180;
                }
        if(dNum == 4){
                return -45*pi/180;
                }
            }

void getAngles(float outAng[], int dNum, int ix, int iy){
        float xt = iXPos(dNum)+cellXLen*ix*cos(detecAng(dNum));
        float yt = iYPos(dNum)-cellXLen*ix*sin(detecAng(dNum));
        float xp = xt;
        float yp = startHeight+cellYLen*iy;
        /*float theta outAng[0] = atan(yt/xt);
        /*float phi outAng[1] = atan(yp/xp);
        /*return theta,phi;
            }*/

float CheckGeom(float the) {
    // int iy = 0;
    // int ix = 0;
    // int dNum = 2;
    
    // float ang[2];
    // getAngles(ang,dNum,ix,iy);
    // Double_t v1f1 = cos(phi)*(sqrt(pow(a,2)+pow(c,2)+2*a*c*cos(the))/(a+c))*v1i;
    // Double_t v1f2 = sin(phi)*(sqrt(pow(a,2)+pow(c,2)-2*a*c*cos(phi))/(a+c))*v1i;
    // Double_t v1f = sqrt(pow(v1f1,2)+pow(v1f2,2));
    // Double_t ejoules = 0.5*a*pow(v1f,2);
    // Double_t eV = ejoules*6.24150934*pow(10,15);
    Double_t gamma = the;// asin(sqrt(pow(sin(the),2)+pow(tan(phi),2))*cos(phi));
    Double_t p1 = sqrt(2*a*B);
    Double_t gamma2 = pi/2+gamma/2-0.5*asin((a/c)*sin(gamma))-gamma;
    Double_t p1p = (sin(gamma2)/sin(gamma))*sqrt((pow(p1,2)/(2*a))*pow((pow(sin(gamma2),2)/(2*a*pow(sin(gamma),2))+1/(2*c)),-1));
    Double_t eV = pow(p1p,2)/(2*a);
    Double_t keV = eV/1000;
    
    //        cout << "Gamma: "<<gamma<<"rad "<<"Gamma2: "<<gamma2<<"rad "<<"Alpha Momentum: "<<p1p<<endl;
    //	cout << "Energy: "<<keV<<"keV"<<endl;
    if(eV>maxE){
        cout << "Error! E too large!" << endl;
    }
    return keV;
}

// sorting in descending order using insertion method
void sortins2(Double_t arr[][ASIC_MAX_HIT], Int_t brr[][ASIC_MAX_HIT], Int_t n[])
{
	int id;
	for (id = 0; id < 4; id++) {
		int i, j;
		Double_t a;
		Int_t b;
		if (n[id] > 1) {
			for (j = 1; j < n[id]; j++) {
        			a = arr[id][j];
       				b = brr[id][j];
        			i = j;
       				while (i > 0 && arr[id][i-1] < a) {
       					arr[id][i] = arr[id][i-1];
					brr[id][i] = brr[id][i-1];
            				i--;
       				}
        			arr[id][i] = a;
        			brr[id][i] = b;
    			}
		}
	}
}

//e0 is energy that we want to split between e1 and e2
Double_t chimatch(Double_t e0, Double_t e1, Double_t e2)
{
	Double_t temp, chi2;
	Double_t step  = 0;

	Double_t diff = e0 - e1 - e2;

	chi2 = pow((step),2)/pow(e1,2) + pow((e0 - (e1+step) - e2),2)/pow(e2,2);

	if(diff >= 0) step += 0.01;
	else step -= 0.01;

	temp = pow((step),2)/pow(e1,2) + pow((e0 - (e1+step) - e2),2)/pow(e2,2); 

	while(chi2 > temp){
		chi2 = temp;

		if(diff >= 0) step += 0.01;
		else step -= 0.01;

		temp = pow((step),2)/pow(e1,2) + pow((e0 - (e1+step) - e2),2)/pow(e2,2); 
	}

	return step;
}

Double_t twobytwo(Double_t a,Double_t b, Double_t c, Double_t d)
{
	Double_t arr[4] = {0};
	Int_t arg = 0;
	if(a>0) { arr[arg] = a; arg++; }
	if(b>0) { arr[arg] = b; arg++; }
	if(c>0) { arr[arg] = c; arg++; }
	if(d>0) { arr[arg] = d; arg++; }

	Double_t temp = 1000;
	for(int i = 0; i < arg; i++) {
        if(arr[arg] < temp){
			temp = arr[arg];
        }
	}
	return temp;
}

Int_t three(Double_t a, Double_t b, Double_t c)
{
	Double_t arr[3] = {0};
	arr[0] = a;
    arr[1] = b;
    arr[2] = c;
	Double_t temp = 1000;
	Int_t CASE = 0;
	for(int i = 0; i < 3; i++){
		if(arr[i] < temp && arr[i] > 0) {
            temp = arr[i];
            CASE = i+1;
        }
	}
	return CASE;
}

/*
Double_t defaultOffSets[4][5]={0.0};	// default offset values for Det0-3 and Errcode0-4
defaultOffSets[0][1]=288.3;
defaultOffSets[1][1]=220.0;
defaultOffSets[2][1]=205.3;
defaultOffSets[3][1]=220.0;
defaultOffSets[0][3]=113.4;
defaultOffSets[1][3]=210.1;
defaultOffSets[2][3]=198.2;
defaultOffSets[3][3]=228.9;
defaultOffSets[0][2]=288.3/2;
defaultOffSets[1][2]=220.0/2;
defaultOffSets[2][2]=205.3/2;
defaultOffSets[3][2]=220.0/2;
defaultOffSets[0][4]=113.4/2;
defaultOffSets[1][4]=210.1/2;
defaultOffSets[2][4]=198.2/2;
defaultOffSets[3][4]=228.9/2;
*/

Int_t match(Int_t ndet, Int_t nhitf, Int_t nhitb, signal sigF[], signal sigB[], hit hit[], Double_t err0=200.0)
{
	Double_t var0=err0*err0;	// condition for a possible match
	Int_t ngood=0;	// number of good matched hits

	Double_t defaultOffSets[4][5] = {0.0};	// default offset values for Det0-3 and Errcode0-4
	defaultOffSets[0][1]=288.3;
	defaultOffSets[1][1]=220.0;
	defaultOffSets[2][1]=205.3;
	defaultOffSets[3][1]=220.0;
	defaultOffSets[0][3]=113.4;
	defaultOffSets[1][3]=210.1;
	defaultOffSets[2][3]=198.2;
	defaultOffSets[3][3]=228.9;
	defaultOffSets[0][2]=288.3/2;
	defaultOffSets[1][2]=220.0/2;
	defaultOffSets[2][2]=205.3/2;
	defaultOffSets[3][2]=220.0/2;
	defaultOffSets[0][4]=113.4/2;
	defaultOffSets[1][4]=210.1/2;
	defaultOffSets[2][4]=198.2/2;
	defaultOffSets[3][4]=228.9/2;

	if (nhitf<=0 || nhitb<=0) return 0;

	Int_t nf_ini[ASIC_MAX_HIT];
	Int_t nb_ini[ASIC_MAX_HIT];
	Int_t nf_fin[ASIC_MAX_HIT];
	Int_t nb_fin[ASIC_MAX_HIT];
	Int_t nf_min[ASIC_MAX_HIT];
	Int_t nb_min[ASIC_MAX_HIT];

	Double_t varmin=0;
//	nb_ini[0]=0;
	// initial rough matching
	for(int jf=0; jf<nhitf; jf++) {
		int jb0=0;
		if (ngood>0) jb0=nb_ini[ngood-1]+1;
		for(int jb=jb0; jb<nhitb; jb++) {
			Double_t var = (sigF[jf].GetEnergy()-sigB[jb].GetEnergy())*(sigF[jf].GetEnergy()-sigB[jb].GetEnergy());
			if (var < var0) {
				varmin += var;
				nf_ini[ngood] = jf;
				nb_ini[ngood] = jb;
				nf_min[ngood] = jf;
				nb_min[ngood] = jb;
				ngood++;
				if (jb==nhitb-1) goto endfor2;
				break;
			}
		}
	}
endfor2:

	if (ngood==0) goto doublehit;

	// find the range of back strips for a given front strip
	for(int jj=ngood;  jj>0; jj--) {
		int jf = nf_ini[jj-1];
		for (int jb=nb_ini[jj-1]; jb<nhitb; jb++) {
			Double_t var = (sigF[jf].GetEnergy()-sigB[jb].GetEnergy())*(sigF[jf].GetEnergy()-sigB[jb].GetEnergy());
			if (var < var0) nb_fin[jj-1] = jb;
		}
		if ((jj < ngood) && nb_fin[jj-1] >= nb_fin[jj]) nb_fin[jj-1] = nb_fin[jj] - 1;
	}
	// find the range of front strips for a given back strip
	for(int jj=ngood;  jj>0; jj--) {
		int jb = nb_ini[jj-1];
		for (int jf=nf_ini[jj-1]; jf<nhitf; jf++) {
			Double_t var = (sigF[jf].GetEnergy()-sigB[jb].GetEnergy())*(sigF[jf].GetEnergy()-sigB[jb].GetEnergy());
			if (var < var0) nf_fin[jj-1] = jf;
		}
		if ((jj < ngood) && nf_fin[jj-1] >= nf_fin[jj]) nf_fin[jj-1] = nf_fin[jj] - 1;
	}

	// chi2 calculation for at most 4 hits in one detector
	Int_t nf[ASIC_MAX_HIT];
	Int_t nb[ASIC_MAX_HIT];
	Double_t var,v0,v1,v2,v3;
	if (ngood <= 4) {

		// vary Front strips for a given back strip
		for (int jf0=nf_ini[0]; jf0<nf_fin[0]; jf0++) {
			nf[0]=jf0;
			int jb0 = nb_ini[0];
			v0 = (sigF[jf0].GetEnergy()-sigB[jb0].GetEnergy())*(sigF[jf0].GetEnergy()-sigB[jb0].GetEnergy());
			if (ngood < 2) {
				var = v0;
				if (var < varmin) {
					for (int j=0; j<ngood; j++) {
						nf_min[j] = nf[j];
						nb_min[j] = nb_ini[j];
						varmin = var;
					}
				}
			} else {
				int jf10 = max(nf[0]+1,nf_ini[1]);
				for (int jf1=jf10; jf1<nf_fin[1]; jf1++) {
					nf[1] = jf1;
					int jb1 = nb_ini[1];
					v1 = (sigF[jf1].GetEnergy()-sigB[jb1].GetEnergy())*(sigF[jf1].GetEnergy()-sigB[jb1].GetEnergy());
					if (ngood < 3) {
						var = v0 + v1;
						if (var < varmin) {
							for (int j=0; j<ngood; j++) {
								nf_min[j] = nf[j];
								nb_min[j] = nb_ini[j];
								varmin = var;
							}
						}
					} else {
						int jf20 = max(nf[1]+1,nf_ini[2]);
						for (int jf2=jf20; jf2<nf_fin[2]; jf2++) {
							nf[2] = jf2;
							int jb2 = nb_ini[2];
							v2 = (sigF[jf2].GetEnergy()-sigB[jb2].GetEnergy())*(sigF[jf2].GetEnergy()-sigB[jb2].GetEnergy());
							if (ngood < 4) {
								var = v0 + v1 + v2;
								if (var < varmin) {
									for (int j=0; j<ngood; j++) {
										nf_min[j] = nf[j];
										nb_min[j] = nb_ini[j];
										varmin = var;
									}
								}
							} else {
								int jf30 = max(nf[2]+1,nf_ini[3]);
								for (int jf3=jf30; jf3<nf_fin[3]; jf3++) {
									nf[3] = jf3;
									int jb3 = nb_ini[3];
									v3 = (sigF[jf3].GetEnergy()-sigB[jb3].GetEnergy())*(sigF[jf3].GetEnergy()-sigB[jb3].GetEnergy());
									if (ngood < 5) {
										var = v0 + v1 + v2 + v3;
										if (var < varmin) {
											for (int j=0; j<ngood; j++) {
												nf_min[j] = nf[j];
												nb_min[j] = nb_ini[j];
												varmin = var;
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}

		// vary back strips for a given front strip
		for (int jb0=nb_ini[0]; jb0<nb_fin[0]; jb0++) {
			nb[0]=jb0;
			int jf0 = nf_ini[0];
			v0 = (sigF[jf0].GetEnergy()-sigB[jb0].GetEnergy())*(sigF[jf0].GetEnergy()-sigB[jb0].GetEnergy());
			if (ngood < 2) {
				var = v0;
				if (var < varmin) {
					for (int j=0; j<ngood; j++) {
						nb_min[j] = nb[j];
						nf_min[j] = nf_ini[j];
						varmin = var;
					}
				}
			} else {
				int jb10 = max(nb[0]+1,nb_ini[1]);
				for (int jb1=jb10; jb1<nb_fin[1]; jb1++) {
					nb[1] = jb1;
					int jf1 = nf_ini[1];
					v1 = (sigF[jf1].GetEnergy()-sigB[jb1].GetEnergy())*(sigF[jf1].GetEnergy()-sigB[jb1].GetEnergy());
					if (ngood < 3) {
						var = v0 + v1;
						if (var < varmin) {
							for (int j=0; j<ngood; j++) {
								nb_min[j] = nb[j];
								nf_min[j] = nf_ini[j];
								varmin = var;
							}
						}
					} else {
						int jb20 = max(nb[1]+1,nb_ini[2]);
						for (int jb2=jb20; jb2<nb_fin[2]; jb2++) {
							nb[2] = jb2;
							int jf2 = nf_ini[2];
							v2 = (sigF[jf2].GetEnergy()-sigB[jb2].GetEnergy())*(sigF[jf2].GetEnergy()-sigB[jb2].GetEnergy());
							if (ngood < 4) {
								var = v0 + v1 + v2;
								if (var < varmin) {
									for (int j=0; j<ngood; j++) {
										nb_min[j] = nb[j];
										nf_min[j] = nf_ini[j];
										varmin = var;
									}
								}
							} else {
								int jb30 = max(nb[2]+1,nb_ini[3]);
								for (int jb3=jb30; jb3<nb_fin[3]; jb3++) {
									nb[3] = jb3;
									int jf3 = nf_ini[3];
									v3 = (sigF[jf3].GetEnergy()-sigB[jb3].GetEnergy())*(sigF[jf3].GetEnergy()-sigB[jb3].GetEnergy());
									if (ngood < 5) {
										var = v0 + v1 + v2 + v3;
										if (var < varmin) {
											for (int j=0; j<ngood; j++) {
												nb_min[j] = nb[j];
												nf_min[j] = nf_ini[j];
												varmin = var;
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}

	}

	// pair up good single hits using the above chi2 approach
	for (int j = 0; j < ngood; j++) {
		int jf = nf_min[j];
		int jb = nb_min[j];
		sigF[jf].SetStatus(1);
		sigB[jb].SetStatus(1);
		Double_t ene = (sigF[jf].GetEnergy() + sigB[jb].GetEnergy())/2;
		Double_t espread = sigF[jf].GetEnergy()-sigB[jb].GetEnergy();
		int nstripf = sigF[jf].GetStrip();
		int nstripb = sigB[jb].GetStrip();
		hit[j].SetValues(ene,espread,ndet,nstripf,nstripb,false,0);	// errcode=0 means normal hit
	}

doublehit:

	// found out all the unpaired signals
	int nhitf2 = 0;
	int nhitb2 = 0;
	signal sigF2[nhitf];
	signal sigB2[nhitb];
	for (int j = 0; j<nhitf; j++) {
		if (sigF[j].GetStatus()==0) {
			sigF2[nhitf2].Apply(sigF[j]);
			nhitf2++;
		}
	}
	for (int j = 0; j<nhitb; j++) {
		if (sigB[j].GetStatus()==0) {
			sigB2[nhitb2].Apply(sigB[j]);
			nhitb2++;
		}
	}

	// search double-hit like signals (one Ef and two Eb)
	if (nhitf2 > 0 && nhitb2 > 1) {
		for (int jf=0; jf<nhitf2; jf++) {
			Double_t ef=sigF2[jf].GetEnergy();
			for (int jb1=0; jb1<nhitb2-1; jb1++) {
				Double_t eb1=sigB2[jb1].GetEnergy();
				for (int jb2=jb1+1; jb2<nhitb2; jb2++) {
					Double_t eb2=sigB2[jb2].GetEnergy();
					Double_t var = (ef - eb1 - eb2 + defaultOffSets[ndet][1])*(ef - eb1 - eb2 + defaultOffSets[ndet][1]);
					if (var < var0) {
						sigF2[jf].SetStatus(1);
						sigB2[jb1].SetStatus(1);
						sigB2[jb2].SetStatus(1);
						int nstripf = sigF2[jf].GetStrip();
						int nstripb1 = sigB2[jb1].GetStrip();
						int nstripb2 = sigB2[jb2].GetStrip();
						if (abs(nstripb1-nstripb2)==1) {	// split hit
							// Double_t ene = (ef + eb1 + eb2)/2;	// split energy is not reliable
							Double_t ene = ef;
							Double_t espread = ef-eb1-eb2 + defaultOffSets[ndet][1];
							int nstripb = nstripb1;
							if (nstripb2 < nstripb1) nstripb = nstripb2;
							hit[ngood].SetValues(ene,espread,ndet,nstripf,nstripb,false,1);	// later the angle should be calculated from nstripb+0.5
							ngood++;
						} else { // double hit
							Double_t ef1 = ef * eb1/(eb1+eb2);
							Double_t ef2 = ef * eb2/(eb1+eb2);
							// Double_t e1 = (ef1 + eb1)/2;
							// Double_t e2 = (ef2 + eb2)/2;
							Double_t e1 = ef1;
							Double_t e2 = ef2;
							Double_t espread1 = ef1-eb1 + defaultOffSets[ndet][2];
							Double_t espread2 = ef2-eb2 + defaultOffSets[ndet][2];
							hit[ngood].SetValues(e1,espread1,ndet,nstripf,nstripb1,true,2);
							ngood++;
							hit[ngood].SetValues(e2,espread2,ndet,nstripf,nstripb2,true,2);
							ngood++;
						}
						goto end2hit;
					}
				}
			}
		}
	}
	// search double-hit like signals (one Eb and two Ef)
	if (nhitb2 > 0 && nhitf2 > 1) {
		for (int jb=0; jb<nhitb2; jb++) {
			Double_t eb=sigB2[jb].GetEnergy();
			for (int jf1=0; jf1<nhitf2-1; jf1++) {
				Double_t ef1=sigF2[jf1].GetEnergy();
				for (int jf2=jf1+1; jf2<nhitf2; jf2++) {
					Double_t ef2=sigF2[jf2].GetEnergy();
					Double_t var = (eb - ef1 - ef2 + defaultOffSets[ndet][3])*(eb - ef1 - ef2 + defaultOffSets[ndet][3]);
					if (var < var0) {
						sigB2[jb].SetStatus(1);
						sigF2[jf1].SetStatus(1);
						sigF2[jf2].SetStatus(1);
						int nstripb = sigB2[jb].GetStrip();
						int nstripf1 = sigF2[jf1].GetStrip();
						int nstripf2 = sigF2[jf2].GetStrip();
						if (abs(nstripf1-nstripf2)==1) {	// split hit
							// Double_t ene = (eb + ef1 + ef2)/2; // split energy is not reliable
							Double_t ene = eb;
							Double_t espread = eb-ef1-ef2 + defaultOffSets[ndet][3];
							int nstripf = nstripf1;
							if (nstripf2 < nstripf1) nstripf = nstripf2;
							hit[ngood].SetValues(ene,espread,ndet,nstripf,nstripb,false,3);	// later the angle should be calculated from nstripb+0.5
							ngood++;
						} else { // double hit
							Double_t eb1 = eb * ef1/(ef1+ef2);
							Double_t eb2 = eb * ef2/(ef1+ef2);
							// Double_t e1 = (eb1 + ef1)/2;
							// Double_t e2 = (eb2 + ef2)/2;
							Double_t e1 = eb1;
							Double_t e2 = eb2;
							Double_t espread1 = eb1-ef1 + defaultOffSets[ndet][4];
							Double_t espread2 = eb2-ef2 + defaultOffSets[ndet][4];
							hit[ngood].SetValues(e1,espread1,ndet,nstripf1,nstripb,true,4);
							ngood++;
							hit[ngood].SetValues(e2,espread2,ndet,nstripf2,nstripb,true,4);
							ngood++;
						}
						goto end2hit;
					}
				}
			}
		}
	}

end2hit:

	return ngood;
}
