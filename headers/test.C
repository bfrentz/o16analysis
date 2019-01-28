/*******************************************
 * The main macro for 13C(a,n)16O* analysis.
 * Usage with root (see also run_all.C):
 *       .L libExpEvent.so
 *       .L analysis_xtalk.C
 *       Run(ENERGY)
 *       -> Submit as a batch job with qsub
 *       -> Using the file run_all.C
 * Updated February 2017
 * This introduces the new crosstalk 
 * calibration to the analysis. 
 * Bryce Frentz (bfrentz@nd.edu)
 * Armen Gyurjinyan (agyurjin@nd.edu)
 * Wanpeng Tan (wtan@nd.edu)
 * Ethan Sauer (esauer2@nd.edu)
 ********************************************/

#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>
#include <sstream>
#include "TCanvas.h"
#include "TRandom2.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TCut.h"

#include "eventData.h"
#include "functions1.6.h" // New matching algorithm in 1.6, correcting for shifts in Espread
#include "hit1.1.h"
#include "signal1.1.h"

// from the evt_config.h file (included by eventData.h)
//#define ASIC_NUM_CB 9
//#define ASIC_MAX_HIT 32

using namespace std;


void Run(Int_t eGroup){
    clock_t start = clock();
    
    // create a random # generator [0,1]
    TRandom2 myRandom(0);
    
    // Create an output file for the histograms
    ostringstream stm;
    stm << eGroup;
    TFile file("/afs/crc.nd.edu/group/nsl/nuc/users/bfrentz/analysis/output/data/xtalk/data" + TString(stm.str()) + "MeV_xtalk.root","RECREATE");

    TChain chain("evtTree");
    
    // Initialize values
    Double_t ebeamMeV = 0;
    Int_t Be9Gate1 = 0;
    Int_t Be9Gate2 = 0;
    
    Int_t He5Gate0 = 0;
    Int_t He5Gate1 = 0;
    Int_t He5Gate2 = 0;
    Int_t He5Gate3 = 0;
    Int_t He5Gate4 = 0;
    
    Int_t Q8BeMin = 0;
    Int_t Q8BeMax = 0;
    Int_t Q16OMin = 0;
    Int_t Q16OMax = 0;
    Double_t Q16Olow = 0;
    Double_t Q16Ohigh = 0;
    
    Double_t L1Min = 0;
    Double_t L1Max = 0;
    Double_t L2Min = 0;
    Double_t L2Max = 0;
    Double_t L3Min = 0;
    Double_t L3Max = 0;
    Double_t L4Min = 0;
    Double_t L4Max = 0;
    Double_t L5Min = 0;
    Double_t L5Max = 0;
    Double_t L6Min = 0;
    Double_t L6Max = 0;
    
    Double_t match_cond = 600;

    // For adding files in chain
    ifstream fRuns;
    char dataf[255];
    int nRuns;
    int run_num;
    
    //Optimize cache and disk reading because this isn't default
    //in this version of ROOT.
    chain.SetCacheSize(1E8);
    chain.AddBranchToCache("*");
    
    // Define gates for each energy group
    // Add data files to the chain
    switch(eGroup){
        
        case 24: {
            ebeamMeV = 24.0;
            
            Be9Gate1 = 2400;
            Be9Gate2 = 2000;
            
            Q8BeMin = 76;
            Q8BeMax = 104;
            Q16OMin = -13500.0;
            Q16OMax = -11600.0;
            Q16Olow = -13632.2;
            Q16Ohigh = -11621.4;
            
            He5Gate1 = 1500;
            He5Gate2 = 3000;
            He5Gate3 = 4750;
            
            L1Min = 16.5404;
            L1Max = 17.4709;
            L2Min = 17.2541;
            L2Max = 18.5409;
            L3Min = 18.5409;
            L3Max = 18.9182;
            L4Min = 19.0454;
            L4Max = 19.4317;
            
            // Add data to the chain
            fRuns.open("/afs/crc.nd.edu/group/nsl/nuc/users/bfrentz/analysis/input/24MeV.in", ios::in);
            fRuns.seekg(0, ios::beg);
            fRuns >> nRuns;
            fRuns.seekg(1, ios::cur);
            for(int i = 0; i < nRuns; i++){
                fRuns >> run_num;
                sprintf(dataf,"/afs/crc.nd.edu/group/nsl/nuc/users/wtan/2015Mar_16O/evt2root/data/root/run%d.root", run_num);
                cout << "Data File: " << dataf << endl;
                chain.Add(dataf);
            }
            fRuns.close();
            //Number of events: 75135594
            /*
            chain.Add("/afs/crc.nd.edu/group/nsl/nuc/user/wtan/2015Mar_16O/evt2root/data/root/run82.root");
            chain.Add("/afs/crc.nd.edu/group/nsl/nuc/user/wtan/2015Mar_16O/evt2root/data/root/run83.root");
            cout<<"Data file: /afs/crc.nd.edu/group/nsl/nuc/user/wtan/2015Mar_16O/evt2root/data/root/run82.root" <<endl;
            cout<<"Data file: /afs/crc.nd.edu/group/nsl/nuc/user/wtan/2015Mar_16O/evt2root/data/root/run83.root" <<endl;
            */
            break;
        }

        case 25: {
            ebeamMeV = 25.0;

            Be9Gate1 = 4000;
            Be9Gate2 = 3500;
            
            Q8BeMin = 74;
            Q8BeMax = 104;
            Q16OMin = -13600;
            Q16OMax = -11500;
            Q16Olow = -13730.5;
            Q16Ohigh = -11507.1;
            
            He5Gate1 = 1750;
            He5Gate2 = 2750;
            He5Gate3 = 4600;
            He5Gate4 = 11750;
            
            L1Min = 16.9417;
            L1Max = 17.274;
            L2Min = 17.274;
            L2Max = 18.4438;
            L3Min = 18.2719;
            L3Max = 18.998;
            L4Min = 18.8696;
            L4Max = 20.0382;
            
            // Add data to the chain
            //Number of events: 728489588
            fRuns.open("/afs/crc.nd.edu/group/nsl/nuc/users/bfrentz/analysis/input/25MeV.in", ios::in);
            fRuns.seekg(0, ios::beg);
            fRuns >> nRuns;
            fRuns.seekg(1, ios::cur);
            for(int i = 0; i < nRuns; i++){
                fRuns >> run_num;
                sprintf(dataf,"/afs/crc.nd.edu/group/nsl/nuc/users/wtan/2015Mar_16O/evt2root/data/root/run%d.root", run_num);
                cout << "Data File: " << dataf << endl;
                chain.Add(dataf);
            }
            fRuns.close();
            /*
            for (Int_t i=6; i<9; i++) {
                sprintf(dataf,"/afs/crc.nd.edu/group/nsl/nuc/user/wtan/2015Mar_16O/evt2root/data/root/run%d.root",i);
                cout<<"Data file: "<< dataf <<endl;
                chain.Add(dataf);
            }
            chain.Add("/afs/crc.nd.edu/group/nsl/nuc/user/wtan/2015Mar_16O/evt2root/data/root/run10.root");
            cout<<"Data file: /afs/crc.nd.edu/group/nsl/nuc/user/wtan/2015Mar_16O/evt2root/data/root/run10.root" <<endl;
            for (Int_t i=12; i<16; i++) {
                sprintf(dataf,"/afs/crc.nd.edu/group/nsl/nuc/user/wtan/2015Mar_16O/evt2root/data/root/run%d.root",i);
                cout<<"Data file: "<< dataf <<endl;
                chain.Add(dataf);
            }
            chain.Add("/afs/crc.nd.edu/group/nsl/nuc/user/wtan/2015Mar_16O/evt2root/data/root/run17.root");
            cout<<"Data file: /afs/crc.nd.edu/group/nsl/nuc/user/wtan/2015Mar_16O/evt2root/data/root/run17.root" <<endl;
            for (Int_t i=20; i<25; i++) {
                sprintf(dataf,"/afs/crc.nd.edu/group/nsl/nuc/user/wtan/2015Mar_16O/evt2root/data/root/run%d.root",i);
                cout<<"Data file: "<< dataf<<endl;
                chain.Add(dataf);
            }
            for (Int_t i=27; i<33; i++) {
                sprintf(dataf,"/afs/crc.nd.edu/group/nsl/nuc/user/wtan/2015Mar_16O/evt2root/data/root/run%d.root",i);
                cout<<"Data file: "<< dataf<<endl;
                chain.Add(dataf);
            }
            */
            break;
        }
            
        case 27: {
            ebeamMeV = 27.75;
            
            Be9Gate1 = 3800;
            Be9Gate2 = 3800;
            
            Q8BeMin = 73;
            Q8BeMax = 105;
            Q16OMin = -13900.0;
            Q16OMax = -11400.0;
            Q16Olow = -14079.6;
            Q16Ohigh = -11136.1;
            
            He5Gate1 = 1500;
            He5Gate2 = 3300;
            He5Gate3 = 6500;
            
            L1Min = 16.9194;
            L1Max = 17.3265;
            L2Min = 17.6289;
            L2Max = 18.1291;
            L3Min = 18.3966;
            L3Max = 18.8154;
            L4Min = 18.8735;
            L4Max = 19.5167;
            L5Min = 20.1646;
            L5Max = 21.2755;
            
            // Add data to the chain
            //Number of events: 199211739
            fRuns.open("/afs/crc.nd.edu/group/nsl/nuc/users/bfrentz/analysis/input/27MeV.in", ios::in);
            fRuns.seekg(0, ios::beg);
            fRuns >> nRuns;
            fRuns.seekg(1, ios::cur);
            for(int i = 0; i < nRuns; i++){
                fRuns >> run_num;
                sprintf(dataf,"/afs/crc.nd.edu/group/nsl/nuc/users/wtan/2015Mar_16O/evt2root/data/root/run%d.root", run_num);
                cout << "Data File: " << dataf << endl;
                chain.Add(dataf);
            }
            fRuns.close();
            /*
            for (Int_t i=46; i<56; i++) {
                sprintf(dataf,"/afs/crc.nd.edu/group/nsl/nuc/user/wtan/2015Mar_16O/evt2root/data/root/run%d.root",i);
                cout<<"Data file: "<< dataf<<endl;
                chain.Add(dataf);
            }
            */
            break;
        }
            
        case 28: {
            ebeamMeV = 28.1;

            Be9Gate1 = 3600;
            Be9Gate2 = 3800;
            
            Q8BeMin = 74;
            Q8BeMax = 104;
            Q16OMin = -13800.0;
            Q16OMax = -11400.0;
            Q16Olow = -13890.5;
            Q16Ohigh = -11495.2;
            
            He5Gate1 = 1400;
            He5Gate2 = 3300;
            
            L1Min = 16.7848;
            L1Max = 17.4269;
            L2Min = 17.2674;
            L2Max = 18.1315;
            L3Min = 18.4199;
            L3Max = 18.8148;
            L4Min = 18.8148;
            L4Max = 19.6296;
            L5Min = 20.4299;
            L5Max = 21.2249;
            L6Min = 21.4243;
            L6Max = 21.6237;
            
            // Add data to chain
            //Number of events: 12445922
            fRuns.open("/afs/crc.nd.edu/group/nsl/nuc/users/bfrentz/analysis/input/28MeV.in", ios::in);
            fRuns.seekg(0, ios::beg);
            fRuns >> nRuns;
            fRuns.seekg(1, ios::cur);
            for(int i = 0; i < nRuns; i++){
                fRuns >> run_num;
                sprintf(dataf,"/afs/crc.nd.edu/group/nsl/nuc/users/wtan/2015Mar_16O/evt2root/data/root/run%d.root", run_num);
                cout << "Data File: " << dataf << endl;
                chain.Add(dataf);
            }
            fRuns.close();
            /*
            chain.Add("/afs/crc.nd.edu/group/nsl/nuc/user/wtan/2015Mar_16O/evt2root/data/root/run43.root");
            chain.Add("/afs/crc.nd.edu/group/nsl/nuc/user/wtan/2015Mar_16O/evt2root/data/root/run44.root");
            chain.Add("/afs/crc.nd.edu/group/nsl/nuc/user/wtan/2015Mar_16O/evt2root/data/root/run45.root");
            cout<<"Data file: /afs/crc.nd.edu/group/nsl/nuc/user/wtan/2015Mar_16O/evt2root/data/root/run43.root" <<endl;
            cout<<"Data file: /afs/crc.nd.edu/group/nsl/nuc/user/wtan/2015Mar_16O/evt2root/data/root/run44.root" <<endl;
            cout<<"Data file: /afs/crc.nd.edu/group/nsl/nuc/user/wtan/2015Mar_16O/evt2root/data/root/run45.root" <<endl;
            */
            break;
        }

        case 29: {
            ebeamMeV = 29.19;

            Be9Gate1 = 4000;
            Be9Gate2 = 3500;
            
            Q8BeMin = 74;
            Q8BeMax = 104;
            Q16OMin = -14100.0;
            Q16OMax = -11200.0;
            Q16Olow = -14050.1;
            Q16Ohigh = -11348.7;
            
            He5Gate1 = 2000;
            He5Gate2 = 4250;
            He5Gate3 = 7250;
            
            L1Min = 16.9258;
            L1Max = 17.3565;
            L2Min = 17.4229;
            L2Max = 18.2072;
            L3Min = 18.2072;
            L3Max = 19.7493;
            L4Min = 18.9889;
            L4Max = 19.4608;
            L5Min = 20.6056;
            L5Max = 21.1956;
            L6Min = 20.9457;
            L6Max = 22.0903;
            
            // Add data to chain
            //Number of events: 210687427
            fRuns.open("/afs/crc.nd.edu/group/nsl/nuc/users/bfrentz/analysis/input/29MeV.in", ios::in);
            fRuns.seekg(0, ios::beg);
            fRuns >> nRuns;
            fRuns.seekg(1, ios::cur);
            for(int i = 0; i < nRuns; i++){
                fRuns >> run_num;
                sprintf(dataf,"/afs/crc.nd.edu/group/nsl/nuc/users/wtan/2015Mar_16O/evt2root/data/root/run%d.root", run_num);
                cout << "Data File: " << dataf << endl;
                chain.Add(dataf);
            }
            fRuns.close();
            /*
            for (Int_t i=34; i<43; i++) {
                sprintf(dataf,"/afs/crc.nd.edu/group/nsl/nuc/user/wtan/2015Mar_16O/evt2root/data/root/run%d.root",i);
                cout<<"Data file: "<< dataf<<endl;
                chain.Add(dataf);
            }
            for (Int_t i=57; i<77; i++) {
                sprintf(dataf,"/afs/crc.nd.edu/group/nsl/nuc/user/wtan/2015Mar_16O/evt2root/data/root/run%d.root",i);
                cout<<"Data file: "<< dataf<<endl;
                chain.Add(dataf);
            }
            chain.Add("/afs/crc.nd.edu/group/nsl/nuc/user/wtan/2015Mar_16O/evt2root/data/root/run81.root");
            cout<<"Data file: /afs/crc.nd.edu/group/nsl/nuc/user/wtan/2015Mar_16O/evt2root/data/root/run81.root" <<endl;
            */
            break;
        }
        default:
            cout << "Unknown group number!\n";
            break;
    }
    
    // define gates
    TCut gate_alpha_peaks = "(asic_hite[0] > 2000 || asic_hite[1] > 2000) && asic_nhit > 1";
    
    
    //Define the histograms with binning here
    TH1F *h1m = new TH1F("h1m","Matched Hit Multiplicity",ASIC_MAX_HIT+1,-0.5,0.5+ASIC_MAX_HIT);
    TH1F *h1mef = new TH1F("h1mef","Multiplicity in EF",ASIC_MAX_HIT+1,-0.5,0.5+ASIC_MAX_HIT);
    TH1F *h1meb = new TH1F("h1meb","Multiplicity in EB",ASIC_MAX_HIT+1,-0.5,0.5+ASIC_MAX_HIT);
    TH1F *h1en = new TH1F("h1en","Calibrated Energy Spectrum(Matched Hits)",3200,0.5,32000.5);
    TH1F *h1errSingle = new TH1F("h1errSingle","Energy Spread (matched single hits)",1000, -100.0,100.0);
    TH1F *h1errDouble = new TH1F("h1errDouble","Energy Spread (matched double hits)",1000, -100.0,100.0);
    TH1F *h1chi2 = new TH1F("h1chi2","Chi square distribution",1000,0,1000);
    TH2F *hit2d = new TH2F("hit2d","Matched Hit Pattern Theta vs Phi",1110,-185.0,185.0,315,0.0,90.0);
    TH2F *hitxy = new TH2F("hitxy","Matched Hit Pattern Y vs X",1300,-130.0,130.0,144,-36.0,36.0);
    TH2F *hitec = new TH2F("hitec","asic channel vs calibrated_E",290,0.0,290.0,8191,0.5,8191.5);
    TH2F *h2efb = new TH2F("h2efb","EF_strip# vs EB_strip# (matched hits)",130,-1.0,129.0,36,-1.0,35.0);
    TH2F *h2epix = new TH2F("h2epix","Calibrated Energy vs Pixel (matched)",4100,-1.0,4099.0,900,0.5,9000.5);

    TH2F *hitxy2 = new TH2F("hitxy2","Matched Hit Pattern Y vs X for checking",1300,-130.0,130.0,144,-36.0,36.0);
    
    TH1F *h1qval = new TH1F("h1qval","Q-value [MeV] of 13C(a,n+8Be+8Be)",1600, -16.0,0.0);
    TH1F *h1ex16o = new TH1F("h1ex16o","Ex [MeV] in 16O",240, 14.0,26.0);
    TH1F *h1ex16ogated = new TH1F("h1ex16ogated","Ex [MeV] in 16O",240, 14.0,26.0);
    TH1F *h1ex16ogatedN = new TH1F("h1ex16ogatedN","Ex [MeV] in 16O",240, 14.0,26.0);
    
    TH1F *h1angdistL1 = new TH1F("h1angdistL1","Angular Distribution 1st State",180,0.5,180.5);
    TH1F *h1angdistL2 = new TH1F("h1angdistL2","Angular Distribution 2nd State",180,0.5,180.5);
    TH1F *h1angdistL3 = new TH1F("h1angdistL3","Angular Distribution 3rd State",180,0.5,180.5);
    TH1F *h1angdistL4 = new TH1F("h1angdistL4","Angular Distribution 4th State",180,0.5,180.5);
    TH1F *h1angdistL5 = new TH1F("h1angdistL5","Angular Distribution 5th State",180,0.5,180.5);
    TH1F *h1angdistL6 = new TH1F("h1angdistL6","Angular Distribution 6th State",180,0.5,180.5);
    TH1F *h1angdistBG = new TH1F("h1angdistBG","Angular Distribution Background",180,0.5,180.5);
    
    TH2F* angdistL1 = new TH2F("angdistL1","angle vs Energy",400,0.,20.,180,0.,180.);
    TH2F* angdistL2 = new TH2F("angdistL2","angle vs Energy",400,0.,20.,180,0.,180.);
    TH2F* angdistL3 = new TH2F("angdistL3","angle vs Energy",400,0.,20.,180,0.,180.);
    TH2F* angdistL4 = new TH2F("angdistL4","angle vs Energy",400,0.,20.,180,0.,180.);
    TH2F* angdistL5 = new TH2F("angdistL5","angle vs Energy",400,0.,20.,180,0.,180.);
    TH2F* angdistL6 = new TH2F("angdistL6","angle vs Energy",400,0.,20.,180,0.,180.);
    TH2F* angdistBG = new TH2F("angdistBG","angle vs Energy",400,0.,20.,180,0.,180.);
    
    TH1F *CheckG = new TH1F("CheckG","Energy Error",60000,-30000,30000);
    TH2F *AngDep = new TH2F("AngDep","Error vs Angle",1440,-360,360,3000,-3000,3000);
    TH2F *AngDepEn = new TH2F("AngDepEn","Energy vs Angle",1440,-360,360,3000,0,30000);
    
    TH1F *calc = new TH1F("calc","calc",15000,0,30000);
    
    TH1F *strip10 = new TH1F("strip10","strip10",15000,0,30000);
    TH1F *strip9 = new TH1F("strip9","strip9",15000,0,30000);
    TH1F *strip8 = new TH1F("strip8","strip8",15000,0,30000);
    TH1F *strip7 = new TH1F("strip7","strip7",15000,0,30000);
    TH1F *strip6 = new TH1F("strip6","strip6",15000,0,30000);
    TH1F *strip5 = new TH1F("strip5","strip5",15000,0,30000);
    TH1F *strip4 = new TH1F("strip4","strip4",15000,0,30000);
    TH1F *strip3 = new TH1F("strip3","strip3",15000,0,30000);
    TH1F *strip2 = new TH1F("strip2","strip2",15000,0,30000);
    TH1F *strip1 = new TH1F("strip1","strip1",15000,0,30000);
    TH1F *strip0 = new TH1F("strip0","strip0",15000,0,30000);
    
    TH1F *h1q8be = new TH1F("h1q8be","Q-value [keV] of 8Be decay",500, -100.0,350.0);
    TH2F *h2q8be = new TH2F("h2q8be","Sum_Energy vs. Q-value [keV] of 8Be decay",500,-100.0,350.0,1600,0.5,32000.5);
    
    TH1F* h1q9be = new TH1F("h1q9be","Excitation Energy [keV] of 9Be decay to two alphas",1000,0.0,4000.0);
    TH2F* h2q9be = new TH2F("h2q9be","Sum_Energy vs Excitation Energy [keV] of 9Be decay", 1000, 0.0, 4000.0,1600,0.5, 32000.5);
    
    TH1F *h1q12c = new TH1F("h1q12c", "Q-value [MeV] of 12C decay to three alphas", 5000, 0.0, 30.0);
    TH1F *h1ex12c = new TH1F("h1ex12c", "Excitation Energy [MeV] of 12C decay to three alphas", 7000, 0, 35.0);
    TH2F *h2q12c = new TH2F("h2q12c", "Sum_Energy vs. Excitation Energy [MeV] of 12C decay (in 16O reconstruction)", 7000, 0, 35.0, 8000, 0, 40.0);
    
    TH1F *h1qval2 = new TH1F("h1qval2","Q-value [MeV] of 13C(a, 12C_hs + a + n)",1500, -15.0,0.0);
    TH1F *h1ex16o2 = new TH1F("h1ex16o2","Excitation Energy [MeV] in 16O Reconstructed from 12C + a",250, 15.0,25.0);
    
    TH1F *h1ex16o_cut1 = new TH1F("h1ex16o_cut1","Excitation Energy [MeV] in 16O, First Helium Cut",250, 15.0,25.0);
    TH1F *h1ex16o_cut2 = new TH1F("h1ex16o_cut2","Excitation Energy [MeV] in 16O, Second Helium Cut",250, 15.0,25.0);
    TH1F *h1ex16o_cut3 = new TH1F("h1ex16o_cut3","Excitation Energy [MeV] in 16O, Third Helium Cut",250, 15.0,25.0);
    TH1F *h1ex16o_cut4 = new TH1F("h1ex16o_cut4","Excitation Energy [MeV] in 16O, Fourth Helium Cut",250, 15.0,25.0);
    
    TH1F *h1ex13c = new TH1F("h1ex13c", "Excitation Energy [MeV] of 13C decay to 3a+n", 6000, 0, 30.0);
    TH2F *h2ex13c = new TH2F("h2ex13c", "Sum_Energy vs. Excitation Energy [MeV] of 13C decay", 6000, 0, 24.0, 6000, 0.0, 32.0);
    
    TH1F *h1q5he = new TH1F("h1q5he", "Q-Value [MeV] of 5He decay to a+n", 600, -2, 28.0);
    TH2F *h2q5he = new TH2F("h2q5he", "Sum_Energy vs. Q-Value [MeV] of 5He decay", 5000, 0, 24.0, 6000, -2, 28.0);
    
    TH2F* h216Od1 = new TH2F("h216Od1","16O vs 9Be DP1",240,14.,26.,600,0.,30.);
    TH2F* h216Od2 = new TH2F("h216Od2","16O vs 9Be DP2",240,14.,26.,600,0.,30.);
    
    TH2F* hSpread = new TH2F("hSpread", "E vs eSpread (keV) - Total", 500, -1000.0, 1000, 500, 0.0, 30000.0);
    TH2F* hSpread0 = new TH2F("hSpread0", "E vs eSpread (keV) - Detector 0", 500, -1000.0, 1000, 500, 0.0, 30000.0);
    TH2F* hSpread1 = new TH2F("hSpread1", "E vs eSpread (keV) - Detector 1", 500, -1000.0, 1000, 500, 0.0, 30000.0);
    TH2F* hSpread2 = new TH2F("hSpread2", "E vs eSpread (keV) - Detector 2", 500, -1000.0, 1000, 500, 0.0, 30000.0);
    TH2F* hSpread3 = new TH2F("hSpread3", "E vs eSpread (keV) - Detector 3", 500, -1000.0, 1000, 500, 0.0, 30000.0);
    TH2F* hSpreadpix1 = new TH2F("hSpreadpix1", "E vs eSpread (keV) - Pixel (Weird)", 500, -1000.0, 1000, 500, 0.0, 30000.0);
    TH2F* hSpreadpix2 = new TH2F("hSpreadpix2", "E vs eSpread (keV) - Pixel (Ethan)", 500, -1000.0, 1000, 500, 0.0, 30000.0);
    TH2F* hSpreadpix3 = new TH2F("hSpreadpix3", "E vs eSpread (keV) - Double Pixels", 500, -1000.0, 1000, 500, 0.0, 30000.0);

    TH1F* hdiag01 = new TH1F("hdiag01", "Counts vs Strip Identifier - Diagonal upper", 64, -0.5, 63.5);
    TH1F* hdiag02 = new TH1F("hdiag02", "Counts vs Strip Identifier - Diagonal lower", 64, -0.5, 63.5);
    TH1F* hdiag11 = new TH1F("hdiag11", "Counts vs Strip Identifier - Diagonal ", 64, -0.5, 63.5);
    TH1F* hlat11 = new TH1F("hlat11", "Counts vs Strip Identifier - lateral", 64, -0.5, 63.5);
    TH1F* hdiag21 = new TH1F("hdiag21", "Counts vs Strip Identifier - Diagonal far", 64, -0.5, 63.5);
    TH1F* hdiag22 = new TH1F("hdiag22", "Counts vs Strip Identifier - Diagonal close", 64, -0.5, 63.5);
    TH1F* hdiag23 = new TH1F("hdiag23", "Counts vs Strip Identifier - Diagonal strip", 64, -0.5, 63.5);
    TH1F* hdiag31 = new TH1F("hdiag31", "Counts vs Strip Identifier - Diagonal upper", 64, -0.5, 63.5);
    TH1F* hdiag32 = new TH1F("hdiag32", "Counts vs Strip Identifier - Diagonal lower", 64, -0.5, 63.5);
    TH1F* hlat31 = new TH1F("hlat31", "Counts vs Strip Identifier - lateral lower", 64, -0.5, 63.5);
    
    // neutron psd ratio, e.g., ch01/ch17
    Double_t psd[14]={1.0,3753.0/2648.0,3768.0/2659.0,3796.0/2515.0,3724.0/2792.0,3810.0/2864.0,
        3676.0/2838.0,3432.0/2471.0,3796.0/2912.0,3844.0/3064.0,3816.0/2868.0,3808.0/2920.0,1.0,3148.0/2411.0};
    TH1F *h1n[14];
    for (Int_t i=1; i<14; i++) {
        char hname1[255],hname2[255];
        sprintf(hname1,"h1n%d",i);
        sprintf(hname2,"Gated N%d spectrum",i);
        h1n[i] = new TH1F(hname1,hname2,4096,-0.5,4095.5);
    }
    
    
    //neutron histograms for gates
    TH2F* h2n[12];
    for (Int_t i=0; i<12; i++) {
        char hname1[255],hname2[255];
        sprintf(hname1,"h2n%d",i+1);
        sprintf(hname2,"N%d 2D spectrum",i+1);
        h2n[i] = new TH2F(hname1,hname2,4096,-0.5,4095.5,4096,-0.5,4095.5);
    }
    
    TH1F* h1nLong8Be[12];
    for (Int_t i=0; i<12; i++) {
        char hname1[255],hname2[255];
        sprintf(hname1,"h1nLong8Be%d",i+1);
        sprintf(hname2,"N%d spectrum",i+1);
        h1nLong8Be[i] = new TH1F(hname1,hname2,4096,-0.5,4095.5);
    }
    
    TH1F* h1nShort8Be[12];
    for (Int_t i=0; i<12; i++) {
        char hname1[255],hname2[255];
        sprintf(hname1,"h1nShort8Be%d",i+1);
        sprintf(hname2,"N%d spectrum",i+1);
        h1nShort8Be[i] = new TH1F(hname1,hname2,4096,-0.5,4095.5);
    }
    
    TH1F* h1nLong16O[12];
    for (Int_t i=0; i<12; i++) {
        char hname1[255],hname2[255];
        sprintf(hname1,"h1nLong16O%d",i+1);
        sprintf(hname2,"N%d spectrum",i+1);
        h1nLong16O[i] = new TH1F(hname1,hname2,4096,-0.5,4095.5);
    }
    
    TH1F* h1nShort16O[12];
    for (Int_t i=0; i<12; i++) {
        char hname1[255],hname2[255];
        sprintf(hname1,"h1nShort16O%d",i+1);
        sprintf(hname2,"N%d spectrum",i+1);
        h1nShort16O[i] = new TH1F(hname1,hname2,4096,-0.5,4095.5);
    }
    
    TH1F* h1nLongTotal[12];
    for (Int_t i=0; i<12; i++) {
        char hname1[255],hname2[255];
        sprintf(hname1,"h1nLongTotal%d",i+1);
        sprintf(hname2,"N%d spectrum",i+1);
        h1nLongTotal[i] = new TH1F(hname1,hname2,4096,-0.5,4095.5);
    }
    
    TH1F* h1nShortTotal[12];
    for (Int_t i=0; i<12; i++) {
        char hname1[255],hname2[255];
        sprintf(hname1,"h1nShortTotal%d",i+1);
        sprintf(hname2,"N%d spectrum",i+1);
        h1nShortTotal[i] = new TH1F(hname1,hname2,4096,-0.5,4095.5);
    }
    
    
    // Initialize reconstruction info
    Double_t p8bex[2],p8bey[2],p8bez[2],e8be[2],q8be;                   // 8Be
    Double_t p16ox,p16oy,p16oz,e16o,q16o,ex16o;                         // 16O (from 8Be)
    Double_t pnx,pny,pnz,eneu,etot,qval,pbeam,ebeam, vbeam;             // n , beam, and reaction info
    Double_t p12cx,p12cy,p12cz,e12c,q12c,ex12c;                         // 12C
    Double_t p16ox2, p16oy2, p16oz2, e16o2, q16o2, ex16o2;              // 16O (from 12C)
    Double_t pn2x, pn2y, pn2z, eneu2, etot2, qval2;                     // General Neutron
    Double_t p9bex[2], p9bey[2], p9bez[2], e9be[2], q9be[2], ex9be[2];  // 9Be
    Double_t p5hex, p5hey, p5hez, e5he, q5he, ex5he;                    // 5He
    Double_t p13cx, p13cy, p13cz, e13c, q13c, ex13c;                    // 13C
    
    
    //Beam info
    ebeam = ebeamMeV*1000.0; // beam energy in keV
    pbeam = sqrt(2.0*4.002603*ebeam); // beam momentum
    vbeam = sqrt(2.0*ebeam/4.002603); //beam velocity
    
    
    //Hits
    hit hit[ASIC_MAX_HIT];
    
    
    // Read Energy Calibration
    Double_t a[ASIC_NUM_CB][32];           // offset
    Double_t b[ASIC_NUM_CB][32];           // linear
    Double_t c[ASIC_NUM_CB][32];           // quadratic
    Bool_t lgood[ASIC_NUM_CB][32] = {0};    //true if it is a good strip
    
    ifstream in1("./input/CalibDataNew.dat");
    if (!in1.is_open()) {
        cout <<"Error: can't open input file: CalibDataNew.dat"<<endl;
    }
    while (!in1.eof()) {
        Int_t i,j;  // chip board#(1-9), channel#(0-31)
        Int_t t;
        Double_t aa,bb,cc;
        in1 >> i >> j >> aa >> bb >> cc;
        if (!in1.good()) break;
        a[i-1][j] = aa; 
        b[i-1][j] = bb; // chipboard# start from zero
        c[i-1][j] = cc;
        if (aa != 0.0 && bb > 0.0) lgood[i-1][j]=1;
    }
    in1.close();
    
    
    // Read Threshold
    Double_t ET[ASIC_NUM_CB][32] = {0};  //Threshold converted to Energy(keV)
    
    ifstream in5("./input/ThresholdEnergy.dat");
    if(!in5.is_open()) {
        cout << "Error can't open input file: ThresholdEnergy.dat" << endl;
    }
    while(!in5.eof()) {
        Int_t i, j;
        Double_t t;
        in5 >> i >> j >> t;
        ET[i-1][j] = t;
    }
    in5.close();
    
    
    // Read slope and offset for crosstalk calculation
    Double_t slopex[ASIC_NUM_CB][32] = {0};
    
    ifstream in6("./input/CrosstalkNew.dat");
    if(!in6.is_open()) {
        cout << "Error can't open input file: CrosstalkNew.dat" << endl;
    }
    while(!in6.eof()) {
        Int_t i, j;
        Double_t s;
        in6 >> i >> j >> s;
        slopex[i-1][j] = s;
    }
    in6.close();
    
    
    // Read Bad ASIC/DSSD Channels
    ifstream in4("./input/bad.dat");
    if (!in4.is_open()) {
        cout <<"Error: can't open input file: bad.dat"<<endl;
    }
    while (!in4.eof()) {
        Int_t i,j;  // chip board#(1-9), channel#(0-31)
        in4 >> i >> j;
        if (!in4.good()) break;
        lgood[i-1][j]=0;
        cout << "Bad Channels: " << i << " " << j << endl;
    }
    in4.close();
    
    
    //Read Angular Information (z-axis = beam direction)
    Double_t xpos[4][32][32];   // x/r=sin(the)cos(phi)
    Double_t ypos[4][32][32];   // y/r=sin(the)sin(phi)
    Double_t zpos[4][32][32];   // z/r=cos(the)
    Double_t rdis[4][32][32];   // distance from target to pixel
    Double_t the[4][32][32];    // theta of a pixel
    Double_t phi[4][32][32];    // phi
    
    ifstream in2("./input/dssdangles_new.dat");
    if (!in2.is_open()) {
        cout <<"Error: can't open input file: dssdangles_new.dat"<<endl;
    }
    for (Int_t i=0;i<4;i++) {               // detectror #
        for (Int_t j=0;j<32;j++){           // front strip#
            for (Int_t k=0;k<32;k++){       // back strip#
                Double_t xx,yy,zz;
                in2 >> xx >> yy >> zz >> rdis[i][j][k] >> the[i][j][k] >> phi[i][j][k];
                xpos[i][j][k] = xx / rdis[i][j][k];
                ypos[i][j][k] = yy / rdis[i][j][k];
                zpos[i][j][k] = zz / rdis[i][j][k];
                if (!in2.good()) break;
            }
        }
    }
    in2.close();
    
    
    // Read ASIC - DSSD strip arrangement
    Int_t nhitef[4],nhiteb[4];                              // hit numbers of front/back strips
    Int_t nsi[ASIC_NUM_CB]={0};                             // dssd detector # (1,2,3,4) as a function of chip board #  (0 means not used)
    Int_t chef[4][ASIC_MAX_HIT], cheb[4][ASIC_MAX_HIT];     // channel numbers of front/back strips
    Int_t cbef[4][ASIC_MAX_HIT], cbeb[4][ASIC_MAX_HIT];     // chip board numbers of front/back strips
    Int_t stef[4][ASIC_MAX_HIT], steb[4][ASIC_MAX_HIT];     // strip numbers of front/back strips
    Bool_t lfb[ASIC_NUM_CB];                                // if it is for front strips (true) or back strips (false)
    Double_t ef[4][ASIC_MAX_HIT],eb[4][ASIC_MAX_HIT];       // calibrated energies of front/back strips
    Double_t sit[4][ASIC_MAX_HIT];                          // timing
    
    ifstream in3("./input/asic.dat");
    if (!in3.is_open()) {
        cout <<"Error: can't open input file: asic.dat"<<endl;
    }
    for (Int_t i=0;i<8;i++) {
        Int_t ncb, n;
        Bool_t lfront;
        in3 >> ncb >> n >> lfront;
        if (!in3.good()) break;
        nsi[ncb-1] = n;     // change chipboard# start from zero
        lfb[ncb-1] = lfront;
    }
    in3.close();
    
    
    //convert ch# to strip#
    //      strip#   30, 28, ..., 0, 31, 29, ..., 1
    //      channel#  0,  1, ...,15, 16, 17, ...,31
    Int_t nstrip[32];   // strip# corresponding to asic channel# of each chip board
    for (Int_t i=0;i<16;i++) {
        nstrip[i] = 30-i*2;
    }
    for (Int_t i=16;i<32;i++) {
        nstrip[i] = 31-(i-16)*2;
    }
    // DSSD3 and DSSD4 back strips seem to be reversely connected
    Int_t nstrip34b[32];    // strip# corresponding to asic channel# of each chip board
    for (Int_t i=0;i<16;i++) {
        nstrip34b[31-i] = 30-i*2;
    }
    for (Int_t i=16;i<32;i++) {
        nstrip34b[31-i] = 31-(i-16)*2;
    }
    
    //cout<<"Debug: finished reading input"<<endl;
    
    
    // get event structure from the chain
    eventData *event = new eventData();
    chain.SetBranchAddress("event", &event);
    
    //Get total event number for the chain
    ULong64_t nevent = chain.GetEntries();
    cout << "Number of Entries: " << nevent <<"\n";

    //Text file for hit/event output
    //ofstream ftext;
    //ofstream fhits;
    //char stext[255];
    //sprintf(stext,"/afs/crc.nd.edu/group/nsl/nuc/users/bfrentz/analysis/output/runbyrun/small_events%d.txt", nrun);
    //ftext.open(stext);
    //fhits.open("/afs/crc.nd.edu/group/nsl/nuc/user/bfrentz/analysis/output/hits13_new.txt");
    
    
    // Start the event loop
    for (ULong64_t i = 0; i < nevent; i++){ //6
        //  cout << "------------\n";
        chain.GetEvent(i);  // get an event
        
        //cout<<"Debug: start calib event#="<<i<<endl;
        
        // store calibrated silicon energies
        for (Int_t j=0; j < 4; j++) {
            nhitef[j]=0; nhiteb[j]=0;
        }
        
        for (Int_t j=0; j < event->asic_nhit; j++) {    // loop over all hits
            Int_t k=event->asic_ncb[j];                 // chip board number (0-8)
            Int_t l=event->asic_nch[j];                 // channel # (0-31)
            Int_t n=nsi[k]-1;                           // detector number (0,1,2,3)
            if (n>=0) {                                 // used chip board (nsi=0 for unused)
                if (lgood[k][l]) {  // a good strip
                    
                    //Calibrate energy
                    Double_t ecalib;
                    bool f = true;

                    // Apply calibration
                    ecalib = a[k][l] + b[k][l] * (event->asic_hite[j]+myRandom.Rndm()-0.5) + c[k][l]*pow((event->asic_hite[j]+myRandom.Rndm()-0.5),2);   
                    hitec->Fill(k*32+l, ecalib);

                    // Fix Det0 shift on 29 MeV runs
                    if(eGroup == 29 && k == 7 && (56505370 < i && i < 172396367)){
                        ecalib -= 123;
                    }

                    if (ecalib > ET[k][l]) {                      // set a threshold in keV
                        if (lfb[k]) {                             // is a front strip
                            ef[n][nhitef[n]] = ecalib;
                            stef[n][nhitef[n]] = nstrip[l];       //strip number
                            cbef[n][nhitef[n]] = k;               //chipboard number
                            chef[n][nhitef[n]] = l;               //channel number
                            nhitef[n]++;                          // hit# in ef
                            f = true;
                        }
                        else {                                      // is a back strip
                            eb[n][nhiteb[n]] = ecalib;
                            steb[n][nhiteb[n]] = nstrip[l];         //strip number
                            cbeb[n][nhiteb[n]] = k;                 //chipboard number
                            cheb[n][nhiteb[n]] = l;                 //channel number
                            
                            // special treatment for DSSD3 and DSSD4 back strips which were reversed
                            if (n == 2 || n == 3) steb[n][nhiteb[n]] = nstrip34b[l];
                            nhiteb[n]++;                            // hit# in eb
                            f = false;
                        }
                    }
                }
            }
        }
        
        //cout<<"Debug: finished calib"<<endl;
        
        
        //Calculate crosstalk
        Double_t efX[4][ASIC_MAX_HIT] = {0}, ebX[4][ASIC_MAX_HIT] = {0};
        
        for(int j = 0; j<4;j++) {
            //pick a strip
            for(int k = 0;k<nhitef[j];k++) {
                for(int l=0;l<nhitef[j];l++){
                    bool check_strip = (stef[j][k] == stef[j][l]+1 || stef[j][k] == stef[j][l]-1);
                    if(check_strip && l!=k){
                        efX[j][k]+= slopex[cbef[j][k]][chef[j][k]] * ef[j][l];
                    }
                }
            }
            for(int k = 0;k<nhiteb[j];k++) {
                for(int l=0;l<nhiteb[j];l++){
                    bool check_strip = (steb[j][k] == steb[j][l]+1 || steb[j][k] == steb[j][l]-1);
                    if(check_strip && l!=k){
                        ebX[j][k] += slopex[cbeb[j][k]][cheb[j][k]] * eb[j][l];
                    }
                }
            }
        }
        
        //We need new array to do crosstalk subtraction
        Double_t ef_new[4][ASIC_MAX_HIT], eb_new[4][ASIC_MAX_HIT];
        Int_t stef_new[4][ASIC_MAX_HIT], steb_new[4][ASIC_MAX_HIT];
        Int_t nhitef_new[4], nhiteb_new[4];
        
        for(int j=0;j<4;j++) {nhitef_new[j] = 0; nhiteb_new[j] = 0;}
        
        for(int k=0;k<4;k++) {
            for(int j=0;j<nhitef[k];j++) {
                if((ef[k][j] - efX[k][j])>ET[cbef[k][j]][chef[k][j]]) {
                    ef_new[k][nhitef_new[k]] = ef[k][j] - efX[k][j];
                    stef_new[k][nhitef_new[k]] = stef[k][j];
                    nhitef_new[k]++;
                }
            }
            for(int j=0;j<nhiteb[k];j++) {
                if((eb[k][j] - ebX[k][j])>ET[cbeb[k][j]][cheb[k][j]]) {
                    eb_new[k][nhiteb_new[k]] = eb[k][j] - ebX[k][j];
                    steb_new[k][nhiteb_new[k]] = steb[k][j];
                    nhiteb_new[k]++;
                }
            }
        }
        
        // sort ef-stef, eb-steb arrays in energy descending order
        sortins2(ef_new, stef_new, nhitef_new);
        sortins2(eb_new, steb_new, nhiteb_new);
        
        //cout<<"Debug: finished sorting"<<endl;
        
        Int_t neftot=0; //total number of EF hits
        Int_t nebtot=0; //total number of EB hits
        for (Int_t j=0; j < 4; j++) {neftot += nhitef_new[j]; nebtot += nhiteb_new[j];}
        h1mef->Fill(neftot);
        h1meb->Fill(nebtot);
        
        /*
        // OLD MATCHING
        // pair up EF-EB hits in each of four silicons
        //Best case
        Int_t igood=0;  // number of good EB-EF matched hits
        Int_t igoodl=0; // number of good EB-EF matched hits for DSSD 1, 2 (left)
        Int_t igoodr=0; // number of good EB-EF matched hits for DSSD 3, 4 (right)
        for (Int_t j=0; j < 4; j++) {   // loop over detectors
            if(nhitef_new[j]>0 && nhiteb_new[j]>0){
                int f = nhitef_new[j];
                int b = nhiteb_new[j];
                signal sigF[f];
                signal sigB[b];
                
                for(int l=0;l<f;l++)
                    sigF[l].SetValues(ef_new[j][l],stef_new[j][l]);
                for(int l=0;l<b;l++)
                    sigB[l].SetValues(eb_new[j][l],steb_new[j][l]);
                
                bool stop = false;
                while(!stop) {
                    if(fabs(sigF[0].GetEnergy()-sigB[0].GetEnergy()) < match_cond) {
                        hit[igood].SetValues((sigF[0].GetEnergy()+sigB[0].GetEnergy())/2,(sigF[0].GetEnergy()-sigB[0].GetEnergy()),j,sigF[0].GetStrip(),sigB[0].GetStrip(),false);
                        igood++;
                        f--;
                        b--;
                        if(f==0 || b==0) stop = true;
                        for(int l=0;l<f;l++)
                            sigF[l].Apply(sigF[l+1]);
                        for(int l=0;l<b;l++)
                            sigB[l].Apply(sigB[l+1]);
                    }
                    else if(sigF[0].GetEnergy() > sigB[0].GetEnergy()) {
                        bool loopcheck = false;
                        bool matched = false;
                        for(int l=0;l<b;l++){
                            for(int k=l+1;k<b;k++){
                                if(fabs(sigF[0].GetEnergy()-sigB[l].GetEnergy()-sigB[k].GetEnergy()) < match_cond) {
                                    double e1=0, e2=0;
                                    e1 = (sigB[l].GetEnergy() * sigF[0].GetEnergy()) / (sigB[l].GetEnergy()+sigB[k].GetEnergy());
                                    e2 = (sigB[k].GetEnergy()*sigF[0].GetEnergy())/(sigB[l].GetEnergy()+sigB[k].GetEnergy());
                                    hit[igood].SetValues((e1+sigB[l].GetEnergy())/2,(e1-sigB[l].GetEnergy()),j,sigF[0].GetStrip(),sigB[l].GetStrip(),true);
                                    igood++;
                                    hit[igood].SetValues((e2+sigB[k].GetEnergy())/2,(e2-sigB[k].GetEnergy()),j,sigF[0].GetStrip(),sigB[k].GetStrip(),true);
                                    igood++;
                                    f--;
                                    b-=2;
                                    loopcheck = true;
                                    matched = true;
                                    if(f==0 || b==0) stop = true;
                                    
                                    for(int m=0;m<f;m++)
                                        sigF[m].Apply(sigF[m+1]);
                                    
                                    int n=0;
                                    for(int m=0;m<b;m++){
                                        while(n==l || n==k) n++;
                                        
                                        sigB[m].Apply(sigB[n]);
                                        n++;
                                    }
                                    break;
                                }
                            }
                            if(loopcheck) break;
                        }
                        if(!matched){
                            f--;
                            if(f==0 || b==0) stop = true;
                            for(int l=0;l<f;l++)
                                sigF[l].Apply(sigF[l+1]);
                        }
                    }
                    else if(sigF[0].GetEnergy() < sigB[0].GetEnergy()) {
                        bool loopcheck = false;
                        bool matched = false;
                        for(int l=0;l<f;l++){
                            for(int k=l+1;k<f;k++){
                                if(fabs(sigB[0].GetEnergy()-sigF[l].GetEnergy()-sigF[k].GetEnergy()) < match_cond) {
                                    double e1, e2;
                                    e1 = sigF[l].GetEnergy()*(sigF[l].GetEnergy()+sigF[k].GetEnergy())*sigB[0].GetEnergy();
                                    e2 = sigF[k].GetEnergy()*(sigF[l].GetEnergy()+sigF[k].GetEnergy())*sigB[0].GetEnergy();
                                    
                                    hit[igood].SetValues((e1+sigF[l].GetEnergy())/2,(e1-sigF[l].GetEnergy()),j,sigF[l].GetStrip(),sigB[0].GetStrip(),true);
                                    igood++;
                                    hit[igood].SetValues((e2+sigF[k].GetEnergy())/2,(e2-sigF[k].GetEnergy()),j,sigF[k].GetStrip(),sigB[0].GetStrip(),true);
                                    igood++;
                                    f-=2;
                                    b--;
                                    loopcheck = true;
                                    matched = true;
                                    if(f==0 || b==0) stop = true;
                                    
                                    for(int m=0;m<b;m++)
                                        sigB[m].Apply(sigB[m+1]);
                                    
                                    int n=0;
                                    for(int m=0;m<f;m++){
                                        while(n==l || n==k) n++;
                                        
                                        sigF[m].Apply(sigF[n]);
                                        n++;
                                    }
                                    
                                    break;
                                }
                            }
                            if(loopcheck) break;
                        }
                        if(!matched){
                            b--;
                            if(f==0 || b==0) stop = true;
                            for(int l=0;l<b;l++)
                                sigB[l].Apply(sigB[l+1]);
                        }
                    }
                    else {continue;}
                }
            }
        }
        */

        // NEW, CORRECT MATCHING
        // pair up EF-EB hits in each of four silicons
        // use the new sophiscated matching algorithm match() in functions1.5.h
        Int_t igood=0;  // number of good EB-EF matched hits
        Int_t igoodl=0; // number of good EB-EF matched hits for DSSD 1, 2 (left)
        Int_t igoodr=0; // number of good EB-EF matched hits for DSSD 3, 4 (right)
        for (Int_t j=0; j < 4; j++) {   // loop over detectors
            int f = nhitef_new[j];
            int b = nhiteb_new[j];
            if (f>0 && b>0) {
                signal sigF[f];
                signal sigB[b];
                for(int l=0;l<f;l++) sigF[l].SetValues(ef_new[j][l],stef_new[j][l]);
                for(int l=0;l<b;l++) sigB[l].SetValues(eb_new[j][l],steb_new[j][l]);
                igood += match(j, nhitef_new[j], nhiteb_new[j], sigF, sigB, &hit[igood], match_cond);
            }
        }

        //Filling histograms
        double AvgChi2 = 0;
        int DetN = 0;
        int hitnumber=0;
        for(int j=0;j<igood;j++) {

            //if(i < 1000 || ((nevent/2) < i && i < ((nevent/2)+1000)) || (nevent-1000) < i){
            //ftext << i << "\t" << hit[j].GetDetectorN() << "\t" << hit[j].GetStripF() << "\t" << hit[j].GetStripB() << "\t" << hit[j].GetEnergy() << "\t" << hit[j].GetESpread() << endl;
            //}

            Double_t pp = sqrt(2*4.002603*hit[j].GetEnergy()); // momemtum assuming it is alpha
            Double_t px = xpos[hit[j].GetDetectorN()][hit[j].GetStripF()][hit[j].GetStripB()]*pp;
            Double_t py = ypos[hit[j].GetDetectorN()][hit[j].GetStripF()][hit[j].GetStripB()]*pp;
            Double_t pz = zpos[hit[j].GetDetectorN()][hit[j].GetStripF()][hit[j].GetStripB()]*pp;
            Double_t th = the[hit[j].GetDetectorN()][hit[j].GetStripF()][hit[j].GetStripB()];   // theta
            Double_t ph = phi[hit[j].GetDetectorN()][hit[j].GetStripF()][hit[j].GetStripB()];   // phi
            
            hit[j].SetAngles(th,ph);
            hit[j].SetMomentum(px,py,pz);
            
            //if(hit[j].GetDetectorN()==1 && hit[j].GetStripF()==16 && hit[j].GetStripB()==4){
            //  angles->Fill(th);
            //}
            
            if(hit[j].GetDouble()){
                h1errDouble->Fill(hit[j].GetESpread());
            }
            else{
                h1errSingle->Fill(hit[j].GetESpread());
            }
            
            if(DetN == hit[j].GetDetectorN()) {
                AvgChi2 += pow(hit[j].GetESpread(),2);
                hitnumber++;
            }
            else {
                DetN = hit[j].GetDetectorN();
                h1chi2->Fill(AvgChi2/hitnumber);
                AvgChi2 = 0;
                hitnumber = 0;
            }
            
            hit2d->Fill(ph, th);
            hitxy->Fill(xpos[hit[j].GetDetectorN()][hit[j].GetStripF()][hit[j].GetStripB()]*130,ypos[hit[j].GetDetectorN()][hit[j].GetStripF()][hit[j].GetStripB()]*130);

            if(hit[j].GetDetectorN() == 3 && -250 < hit[j].GetESpread() && hit[j].GetESpread() < -60 && 3000 < hit[j].GetEnergy() && hit[j].GetEnergy() < 15000){
                hitxy2->Fill(xpos[hit[j].GetDetectorN()][hit[j].GetStripF()][hit[j].GetStripB()]*130,ypos[hit[j].GetDetectorN()][hit[j].GetStripF()][hit[j].GetStripB()]*130);
            }

            h1en->Fill(hit[j].GetEnergy());
            h2efb->Fill(hit[j].GetStripF()+hit[j].GetDetectorN()*32, hit[j].GetStripB());
            h2epix->Fill(hit[j].GetStripF()+hit[j].GetDetectorN()*1024+hit[j].GetStripB()*32, hit[j].GetEnergy());
        }
        
        
        // Energy spread
        for (int k = 0; k < igood; k++) {
            
            // Total spread
            hSpread->Fill(hit[k].GetESpread(), hit[k].GetEnergy());
            
            // Individual detectors
            if (hit[k].GetDetectorN() == 0) {
                hSpread0->Fill(hit[k].GetESpread(), hit[k].GetEnergy());
            }
            
            if (hit[k].GetDetectorN() == 1) {
                hSpread1->Fill(hit[k].GetESpread(), hit[k].GetEnergy());
            }
            
            if (hit[k].GetDetectorN() == 2) {
                hSpread2->Fill(hit[k].GetESpread(), hit[k].GetEnergy());
            }
            
            if (hit[k].GetDetectorN() == 3) {
                hSpread3->Fill(hit[k].GetESpread(), hit[k].GetEnergy());
            }
            
            if (hit[k].GetDetectorN() == 0 && hit[k].GetStripF() == 4 && hit[k].GetStripB() == 16) {
                hSpreadpix1->Fill(hit[k].GetESpread(), hit[k].GetEnergy());
            }
            
            if (hit[k].GetDetectorN() == 1 && hit[k].GetStripF() == 4 && hit[k].GetStripB() == 16) {
                hSpreadpix2->Fill(hit[k].GetESpread(), hit[k].GetEnergy());
            }

            // Checking the bad strips with the new calibration and looking at the common pixels on a back strip to see if anything weird is happening
            if (hit[k].GetDetectorN() == 0 && hit[k].GetStripF() == 11 && hit[k].GetStripB() == 16) {
                hSpreadpix3->Fill(hit[k].GetESpread(), hit[k].GetEnergy());
            }
            if (hit[k].GetDetectorN() == 0 && hit[k].GetStripF() == 21 && hit[k].GetStripB() == 16) {
                hSpreadpix3->Fill(hit[k].GetESpread(), hit[k].GetEnergy());
            }


            // Checking for bad strips
            // Diagonals and laterals
            
            if (hit[k].GetDetectorN() == 0 && 250 < hit[k].GetESpread() && hit[k].GetESpread() < 460 && 3130 < hit[k].GetEnergy() && hit[k].GetEnergy() < 3500) {
                hdiag01->Fill(hit[k].GetStripF());
                hdiag01->Fill(hit[k].GetStripB() + 32);
            }

            if (hit[k].GetDetectorN() == 0 && 250 < hit[k].GetESpread() && hit[k].GetESpread() < 460 && 3130 < hit[k].GetEnergy() && hit[k].GetEnergy() < 3500) {
                hdiag02->Fill(hit[k].GetStripF());
                hdiag02->Fill(hit[k].GetStripB() + 32);
            }

            if (hit[k].GetDetectorN() == 1 && 250 < hit[k].GetESpread() && hit[k].GetESpread() < 460 && 3130 < hit[k].GetEnergy() && hit[k].GetEnergy() < 3500) {
                hdiag11->Fill(hit[k].GetStripF());
                hdiag11->Fill(hit[k].GetStripB() + 32);
            }

            if (hit[k].GetDetectorN() == 1 && 250 < hit[k].GetESpread() && hit[k].GetESpread() < 460 && 3130 < hit[k].GetEnergy() && hit[k].GetEnergy() < 3500) {
                hlat11->Fill(hit[k].GetStripF());
                hlat11->Fill(hit[k].GetStripB() + 32);
            }

            if (hit[k].GetDetectorN() == 2 && 250 < hit[k].GetESpread() && hit[k].GetESpread() < 460 && 3130 < hit[k].GetEnergy() && hit[k].GetEnergy() < 3500) {
                hdiag21->Fill(hit[k].GetStripF());
                hdiag21->Fill(hit[k].GetStripB() + 32);
            }

            if (hit[k].GetDetectorN() == 2 && 250 < hit[k].GetESpread() && hit[k].GetESpread() < 460 && 3130 < hit[k].GetEnergy() && hit[k].GetEnergy() < 3500) {
                hdiag22->Fill(hit[k].GetStripF());
                hdiag22->Fill(hit[k].GetStripB() + 32);
            }

            if (hit[k].GetDetectorN() == 2 && 250 < hit[k].GetESpread() && hit[k].GetESpread() < 460 && 3130 < hit[k].GetEnergy() && hit[k].GetEnergy() < 3500) {
                hdiag23->Fill(hit[k].GetStripF());
                hdiag23->Fill(hit[k].GetStripB() + 32);
            }

            if (hit[k].GetDetectorN() == 3 && 250 < hit[k].GetESpread() && hit[k].GetESpread() < 460 && 3130 < hit[k].GetEnergy() && hit[k].GetEnergy() < 3500) {
                hdiag31->Fill(hit[k].GetStripF());
                hdiag31->Fill(hit[k].GetStripB() + 32);
            }

            if (hit[k].GetDetectorN() == 3 && 250 < hit[k].GetESpread() && hit[k].GetESpread() < 460 && 3130 < hit[k].GetEnergy() && hit[k].GetEnergy() < 3500) {
                hdiag32->Fill(hit[k].GetStripF());
                hdiag32->Fill(hit[k].GetStripB() + 32);
            }

            if (hit[k].GetDetectorN() == 3 && 250 < hit[k].GetESpread() && hit[k].GetESpread() < 460 && 3130 < hit[k].GetEnergy() && hit[k].GetEnergy() < 3500) {
                hlat31->Fill(hit[k].GetStripF());
                hlat31->Fill(hit[k].GetStripB() + 32);
            }


        }
        
        
        
        
        
        

        // Angles
        for(int k = 0; k < igood; k++){
            Int_t detecNum = hit[k].GetDetectorN();
            Int_t j = hit[k].GetStripF();
            Int_t i = hit[k].GetStripB();
            if(detecNum > 1){
                i = hit[k].GetStripF();         //these are switched!!
                j = hit[k].GetStripB();
            }
            //Double_t theta = hit[k].GetAngleTheta(); //(pi/180)*(the[detecNum][j][i]);
            Double_t phiVal = hit[k].GetAnglePhi(); //(pi/180)*(phi[detecNum][j][i]);
            //Double_t theoVal = CheckGeom(theta,phiVal,detecNum, j, i);
            Double_t e1 = hit[k].GetEnergy();
            
            //cout <<"D: "<<detecNum<<" Vert Strip: " <<j<<" Horiz Strip: "<<i<<" Theta: "<<theta<<" rad "<<"Phi: "<<phiVal<<" rad "<<end$
            
            /*
             if(detecNum == 2){
             cout <<"D: "<<detecNum<<" Vert Strip: " <<j<<" Horiz Strip: "<<i<<" Theta: "<<theta<<" rad "<<"Phi: "<<phiVal<<" ra$
             }
             */
            
            //Double_t valDiff = (theoVal - e1);
            
            //if(hit[k].GetEnergy()>23800 && hit[k].GetEnergy()<24100 && i==16 && j==15 && detecNum == 2){

            if(i==16 && j==10 && detecNum == 1){
                strip10->Fill(e1);
                Double_t theta = hit[k].GetAngleTheta();
                Double_t theoVal = CheckGeom(theta*pi/180);
                strip10->Fill(theoVal);
                Double_t valDiff = (theoVal - e1);
                //CheckG->Fill(valDiff);
                AngDep->Fill(theta,valDiff);
                AngDepEn->Fill(theta,e1);
            }
            
            if(i==16 && j==9 && detecNum == 1){
                strip9->Fill(e1);
                Double_t theta = hit[k].GetAngleTheta();
                Double_t theoVal = CheckGeom(theta*pi/180);
                strip9->Fill(theoVal);
                Double_t valDiff = (theoVal - e1);
                //CheckG->Fill(valDiff);
                AngDep->Fill(theta,valDiff);
                AngDepEn->Fill(theta,e1);
            }
            
            if(i==16 && j==8 && detecNum == 1){
                strip8->Fill(e1);
                Double_t theta = hit[k].GetAngleTheta();
                Double_t theoVal = CheckGeom(theta*pi/180);
                //cout << "Vert Strip Num: " << j<<"Theta: "<<theta<<"Theo E: "<<theoVal<<endl;
                strip8->Fill(theoVal);
                Double_t valDiff = (theoVal - e1);
                //CheckG->Fill(valDiff);
                AngDep->Fill(theta,valDiff);
                AngDepEn->Fill(theta,e1);
            }
            
            if(i==16 && j==7 && detecNum == 1){ //&& hit[k].GetEnergy()>22600 && hit[k].GetEnergy()<23100){
                strip7->Fill(e1);
                Double_t theta = hit[k].GetAngleTheta();
                Double_t theoVal = CheckGeom(theta*pi/180);
                //cout << "Vert Strip Num: " << j<<"Theta: "<<theta<<"Theo E: "<<theoVal<<endl;
                strip7->Fill(theoVal);
                Double_t valDiff = (theoVal - e1);
                //CheckG->Fill(valDiff);
                AngDep->Fill(theta,valDiff);
                AngDepEn->Fill(theta,e1);
            }
            
            if(i==16 && j==6 && detecNum == 1){
                strip6->Fill(e1);
                Double_t theta = hit[k].GetAngleTheta();
                Double_t theoVal = CheckGeom(theta*pi/180);
                //cout << "Vert Strip Num: " << j<<"Theta: "<<theta<<"Theo E: "<<theoVal<<endl;
                strip6->Fill(theoVal);
                Double_t valDiff = (theoVal - e1);
                //CheckG->Fill(valDiff);
                AngDep->Fill(theta,valDiff);
                AngDepEn->Fill(theta,e1);
            }
            
            if(i==16 && j==5 && detecNum == 1){
                strip5->Fill(e1);
                Double_t theta = hit[k].GetAngleTheta();
                Double_t theoVal = CheckGeom(theta*pi/180);
                //cout << "Vert Strip Num: " << j<<"Theta: "<<theta<<"Theo E: "<<theoVal<<endl;
                strip5->Fill(theoVal);
                Double_t valDiff = (theoVal - e1);
                //CheckG->Fill(valDiff);
                AngDep->Fill(theta,valDiff);
                AngDepEn->Fill(theta,e1);
            }
            
            if(i==16 && j==4 && detecNum == 1){
                strip4->Fill(e1);
                Double_t theta = hit[k].GetAngleTheta(); //(pi/180)*(the[detecNum][j][i]);
                Double_t theoVal = CheckGeom(theta*pi/180);
                //cout << "Vert Strip Num: " << j<<"Theta: "<<theta<<"Theo E: "<<theoVal<<endl;
                strip4->Fill(theoVal);
                Double_t valDiff = (theoVal - e1);
                CheckG->Fill(valDiff);
                AngDep->Fill(theta,valDiff);
                AngDepEn->Fill(theta,e1);
            }
            
            if(i==16 && j==3 && detecNum == 1){
                strip3->Fill(e1);
                Double_t theta = hit[k].GetAngleTheta();
                //cout <<   "Vert Strip Num: " << j<<"Theta: "<<theta<<endl;
                Double_t theoVal = CheckGeom(theta*pi/180);
                strip3->Fill(theoVal);
                Double_t valDiff = (theoVal - e1);
                //CheckG->Fill(valDiff);
                AngDep->Fill(theta,valDiff);
                AngDepEn->Fill(theta,e1);
            }
            
            if(i==16 && j==2 && detecNum == 1){
                //tree->Draw("asic_e[2][16]>>strip16","asic_e[2][16]>23800 && asic_e[2][16]<24100");
                strip2->Fill(e1);
                Double_t theta = hit[k].GetAngleTheta();
                //cout << "Vert Strip Num: " << j<<"Theta: "<<theta<<endl;
                Double_t theoVal = CheckGeom(theta*pi/180);
                strip2->Fill(theoVal);
                Double_t valDiff = (theoVal - e1);
                //CheckG->Fill(valDiff);
                AngDep->Fill(theta,valDiff);
                AngDepEn->Fill(theta,e1);
                /*
                 cout << "Theta: "<<theta<<"rad "<<"Phi: "<<phiVal<<"rad "<<endl;
                 cout << "Calculated: "<<theoVal<<"keV"<<endl;
                 cout << "Measured: "<<hit[k].GetEnergy()<<"keV"<<endl;
                 cout << "--------------------------------------------------------------------"<<endl;
                 cout << "Difference: "<<valDiff<<"keV"<<endl;
                 */
            }

            if(i==16 && j==1 && detecNum == 1){
                strip1->Fill(e1);
                Double_t theta = hit[k].GetAngleTheta();
                //cout <<   "Vert Strip Num: " << j<<"Theta: "<<theta<<endl;
                Double_t theoVal = CheckGeom(theta*pi/180);
                strip1->Fill(theoVal);
                Double_t valDiff = (theoVal - e1);
                CheckG->Fill(valDiff);
                AngDep->Fill(theta,valDiff);
                AngDepEn->Fill(theta,e1);
            }
            
            if(i==16 && j==0 && detecNum == 1){
                strip0->Fill(e1);
                Double_t theta = hit[k].GetAngleTheta();
                //cout <<   "Vert Strip Num: " << j<<"Theta: "<<theta<<endl;
                Double_t theoVal = CheckGeom(theta*pi/180);
                strip0->Fill(theoVal);
                Double_t valDiff = (theoVal - e1);
                CheckG->Fill(valDiff);
                AngDep->Fill(theta,valDiff);
                AngDepEn->Fill(theta,e1);
            }
            /*
             CheckG->Fill(valDiff);
             cout << "j: "<<j<<" i: "<<i<<" Detec: "<<detecNum<<endl;
             cout << "Theta: "<<theta<<"rad "<<"Phi: "<<phiVal<<"rad "<<endl;
             cout << "Calculated: "<<theoVal<<"keV"<<endl;
             cout << "Measured: "<<hit[k].GetEnergy()<<"keV"<<endl;
             cout << "Difference: "<<valDiff<<"keV"<<endl;
             */
        }
        
        
        //****************RECONSTRUCTION PART***********************
        // MASS EXCESS
        // n = 8071.31710;
        // 4He = 2424.91565;
        // 5He = 11386.233;
        // 8Be = 4941.672;
        // 9Be = 11347.648;
        // 12C = 0.0;
        // 13C = 3125.01129;
        
        //reconstruct 8Be
        Int_t igood8Be = 0;
        Int_t a[2][2];
        
        if (igood>1) {
            for (Int_t k=0; k<igood-1; k++) {
                if (igood8Be < 2) { // no more than two 8Be
                    for (Int_t l=k+1; l<igood; l++) {
                        
                        // Momenta
                        p8bex[igood8Be]=hit[k].GetMomentumX()+hit[l].GetMomentumX();
                        p8bey[igood8Be]=hit[k].GetMomentumY()+hit[l].GetMomentumY();
                        p8bez[igood8Be]=hit[k].GetMomentumZ()+hit[l].GetMomentumZ();
                        
                        // kinetic energy of 8Be
                        e8be[igood8Be] = ( pow(p8bex[igood8Be],2) + pow(p8bey[igood8Be],2) + pow(p8bez[igood8Be],2) )/2.0/8.005305;
                        q8be = hit[k].GetEnergy()+hit[l].GetEnergy() - e8be[igood8Be];  // Q-value of 8Be decay
                        
                        if ( fabs(q8be-91.84) < 160.0 ) { // q-value  matched in keV
                            h1q8be->Fill(q8be);
                            h2q8be->Fill(q8be,hit[k].GetEnergy()+hit[l].GetEnergy());
                            if ( q8be > Q8BeMin && q8be < Q8BeMax) { // q-value cut that we want to apply
                                a[igood8Be][0] = k;
                                a[igood8Be][1] = l;
                                igood8Be++;
                                break;
                            }
                        }
                    }
                }
            }
        }
        
        //Neutron gates
        double s[12] = {1.4615384615, 1.4615384615, 1.4615384615, 1.3068965517, 1.3107638889, 1.2925170068, 1.3996350365,
            1.2908474576, 1.2756666667, 1.3660714286, 11.351590106, 1.3327526132};
        double o[12] = {0, 0, 0, 10, 25, 0, -35, -8, -27, -25, -25, -25};
        
        //generate 1D Neutron histograms with and without gates
        Double_t ratio;
        if (igood8Be>0) {
            for (Int_t j=1; j<14; j++) {
                if (event->qdc[j] > 0 && event->qdc[j+16] > 0 && event->qdc[j] < 3820) {
                    ratio = (Double_t) event->qdc[j] / (Double_t) event->qdc[j+16];
                    if (j==7 && event->qdc[j] < 1687 && ratio > 1687.0/1290.0 || ratio > psd[j]) {
                        h1n[j]->Fill(event->qdc[j]);
                    }
                }
            }
            
            //Neutron detectors gated with 8Be
            for(int i=0;i<12;i++) {
                if(i==11) {
                    if(event->qdc[i+2] > (s[i]*event->qdc[i+18]+o[i])) {
                        if(event->qdc[i+2]<3800){
                            if(event->qdc[i+2]>0){
                                h1nShort8Be[i]->Fill(event->qdc[i+18]);
                                h1nLong8Be[i]->Fill(event->qdc[i+2]);
                            }
                        }
                    }
                }
                else{
                    if(event->qdc[i+1] > (s[i]*event->qdc[i+17]+o[i])) {
                        if(event->qdc[i+1]<3800){
                            if(event->qdc[i+1]>0){
                                h1nShort8Be[i]->Fill(event->qdc[i+17]);
                                h1nLong8Be[i]->Fill(event->qdc[i+1]);
                            }
                        }
                    }
                }
            }
        }
        
        //generate 2D neutron histograms
        for(int i=0;i<12;i++) {
            if(i==11)
                h2n[i]->Fill(event->qdc[i+18],event->qdc[i+2]);
            else
                h2n[i]->Fill(event->qdc[i+17],event->qdc[i+1]);
        }
        
        //reconstruct 16O from 8Be + 8Be
        if (igood8Be==2) {
            // Momenta
            p16ox = p8bex[0] + p8bex[1];
            p16oy = p8bey[0] + p8bey[1];
            p16oz = p8bez[0] + p8bez[1];
            
            // Energies
            e16o = (p16ox*p16ox+p16oy*p16oy+p16oz*p16oz)/2.0/15.9949146;    // kinetic energy of 16O
            q16o = e8be[0] + e8be[1] - e16o;                                // Q-value (g.s. Q-value = -14.6203MeV)
            ex16o = q16o + 14620.3;                                         // excitation energy
            
            /*
             //Angle of 16O in lab frame with radians
             Double_t angleLab = atan(sqrt(pow(p16ox,2) + pow(p16oy,2))/p16oz);
             //Velocity of 16O in lab frame
             Double_t v16OLab = sqrt(pow(p16ox,2)+pow(p16oy,2)+pow(p16oz,2))/15.9949146;
             //Velocity of center of momentum frame in final stage
             Double_t VCoM = (4.00260325415/(1.00866491574 + 15.99491461956))*vbeam;
             //Angle of 16O in center of mass frame in degrees
             Double_t angleCoM = acos((v16OLab*cos(angleLab)-VCoM)/v16OLab)*180/3.14159265359;
             */
            
            //Angle of 16O in lab frame with radians
            Double_t tgThetaLab = sqrt(pow(p16ox,2) + pow(p16oy,2))/p16oz;
            Double_t cosThetaLab;
            if(tgThetaLab>=0)
                cosThetaLab = 1/sqrt(1+pow(tgThetaLab,2));
            else
                cosThetaLab = -1/sqrt(1+pow(tgThetaLab,2));
            
            //Velocity of 16O in lab frame
            Double_t v16OLab = sqrt(pow(p16ox,2)+pow(p16oy,2)+pow(p16oz,2))/15.9949146;
            
            //Velocity of center of momentum frame in final stage
            Double_t VCoM = (4.00260325415/(1.00866491574 + 15.99491461956))*vbeam;
            
            //Velocity of 16O in center of mass frame
            Double_t V16OCM = sqrt(pow(v16OLab,2)+pow(VCoM,2)-2*v16OLab*VCoM*cosThetaLab);
            
            //Angle of 16O in CoM
            Double_t angleCoM = acos((v16OLab*cosThetaLab - VCoM)/V16OCM)*180/3.14159265359;
            
            bool isL1 = (ex16o>L1Min*1000 && ex16o<L1Max*1000 && L1Max != 0);
            bool isL2 = (ex16o>L2Min*1000 && ex16o<L2Max*1000 && L2Max != 0);
            bool isL3 = (ex16o>L3Min*1000 && ex16o<L3Max*1000 && L3Max != 0);
            bool isL4 = (ex16o>L4Min*1000 && ex16o<L4Max*1000 && L4Max != 0);
            bool isL5 = (ex16o>L5Min*1000 && ex16o<L5Max*1000 && L5Max != 0);
            bool isL6 = (ex16o>L6Min*1000 && ex16o<L6Max*1000 && L6Max != 0);
            bool isBG = !isL1 && !isL2 && !isL3 && !isL4 && !isL5 && !isL6;
            if(isL1){
                h1angdistL1->Fill(angleCoM);
                angdistL1->Fill(e16o/1000.0, angleCoM);
            }
            if(isL2){
                h1angdistL2->Fill(angleCoM);
                angdistL2->Fill(e16o/1000.0, angleCoM);
            }
            if(isL3){
                h1angdistL3->Fill(angleCoM);
                angdistL3->Fill(e16o/1000.0, angleCoM);
            }
            if(isL4){
                h1angdistL4->Fill(angleCoM);
                angdistL4->Fill(e16o/1000.0, angleCoM);
            }
            if(isL5){
                h1angdistL5->Fill(angleCoM);
                angdistL5->Fill(e16o/1000.0, angleCoM);
            }
            if(isL6){
                h1angdistL6->Fill(angleCoM);
                angdistL6->Fill(e16o/1000.0, angleCoM);
            }
            if(isBG){
                h1angdistBG->Fill(angleCoM);
                angdistBG->Fill(e16o/1000.0, angleCoM);
            }
            
            // reconstruct neutron
            pnx = -p16ox;
            pny = -p16oy;
            pnz = pbeam-p16oz;
            eneu = (pnx*pnx+pny*pny+pnz*pnz)/2.0/1.0086649;
            
            
            // reaction Q-value
            etot = e8be[0] + e8be[1] + eneu;
            qval = etot - ebeam;                        // 13C(4he,n+8Be+8Be) Q-value = -12404.734 keV
            h1qval->Fill(qval/1000.0);
            if (qval > Q16OMin && qval < Q16OMax) {     // q-value  matched in keV
                h1ex16o->Fill(ex16o/1000.0);
            }
            
            
            //Q-value(9Be->8Be+n) = −1665.3411
            for(int p=0;p<igood8Be;p++){
                p9bex[p] = p8bex[p]+pnx;
                p9bey[p] = p8bey[p]+pny;
                p9bez[p] = p8bez[p]+pnz;
                e9be[p] = (pow(p9bex[p],2)+pow(p9bey[p],2)+pow(p9bez[p],2))/2.0/9.012182201;
                q9be[p] = e8be[p] + eneu - e9be[p];
                ex9be[p] = q9be[p] + 1665.34;
                
                h1q9be->Fill(q9be[p]);
                h2q9be->Fill(q9be[p],e8be[p] + eneu);
            }
            if(qval > Q16OMin && qval < Q16OMax) {
                h216Od1->Fill(ex16o/1000.0,ex9be[0]/1000.0);
                h216Od2->Fill(ex16o/1000.0,ex9be[1]/1000.0);
            }
            //Best Gate: 4000 && 3500
            if(qval > Q16OMin && qval < Q16OMax && ex9be[0]>Be9Gate1 && ex9be[1]>Be9Gate2){
                h1ex16ogated->Fill(ex16o/1000.0);
                
                //Neutron detectors gated on 16O
                for(int i=0;i<12;i++) {
                    if(i==11) {
                        if(event->qdc[i+2] > (s[i]*event->qdc[i+18]+o[i])) {
                            if(event->qdc[i+2]<3800){
                                if(event->qdc[i+2]>0){
                                    h1nShort16O[i]->Fill(event->qdc[i+18]);
                                    h1nLong16O[i]->Fill(event->qdc[i+2]);
                                }
                            }
                        }
                    }
                    else{
                        if(event->qdc[i+1] > (s[i]*event->qdc[i+17]+o[i])) {
                            if(event->qdc[i+1]<3800){
                                if(event->qdc[i+1]>0){
                                    h1nShort16O[i]->Fill(event->qdc[i+17]);
                                    h1nLong16O[i]->Fill(event->qdc[i+1]);
                                }
                            }
                        }
                    }
                }
            }
        }
        
        //Neutron detectors gated
        for(int i=0;i<12;i++) {
            if(i==11) {
                if(event->qdc[i+2] > (s[i]*event->qdc[i+18]+o[i])) {
                    if(event->qdc[i+2]<3800){
                        if(event->qdc[i+2]>0){
                            h1nShortTotal[i]->Fill(event->qdc[i+18]);
                            h1nLongTotal[i]->Fill(event->qdc[i+2]);
                        }
                    }
                }
            }
            else {
                if(event->qdc[i+1] > (s[i]*event->qdc[i+17]+o[i])) {
                    if(event->qdc[i+1]<3800){
                        if(event->qdc[i+1]>0){
                            h1nShortTotal[i]->Fill(event->qdc[i+17]);
                            h1nLongTotal[i]->Fill(event->qdc[i+1]);
                        }
                    }
                }
            }
        }
        
        
        // Reconstruct 12C
        Int_t igood12C = 0;
        Int_t n12ca1 = 0;
        Int_t n12ca2 = 0;
        Int_t n12ca3 = 0;
        if (igood>2) {
            for (Int_t k=0; k<igood-2; k++) {
                if(igood12C == 1) break;
                else {
                    for (Int_t l=k+1; l<igood-1; l++) {
                        if(igood12C == 1) break;
                        else {
                            for (Int_t j = l+1; j < igood; j++) {
                                // Sum of Momenta
                                p12cx = hit[k].GetMomentumX() + hit[l].GetMomentumX() + hit[j].GetMomentumX();
                                p12cy = hit[k].GetMomentumY() + hit[l].GetMomentumY() + hit[j].GetMomentumY();
                                p12cz = hit[k].GetMomentumZ() + hit[l].GetMomentumZ() + hit[j].GetMomentumZ();
                                // Energies
                                e12c = (p12cx*p12cx + p12cy*p12cy + p12cz*p12cz)/2.0/12.0;                                  // Kinetic energy of 12C
                                q12c = (hit[k].GetEnergy() + hit[l].GetEnergy() + hit[j].GetEnergy()) - e12c;               // Measured Q value (g.s. Q-value = -7274.7 keV)
                                h1q12c->Fill(q12c/1000.0);
                                ex12c = (hit[k].GetEnergy() + hit[l].GetEnergy() + hit[j].GetEnergy()) - e12c + 7274.7;     // Excitation Energy
                                
                                // Fill histogram with excitation energy
                                h1ex12c->Fill(ex12c/1000.0);
                                
                                // Fill 2D histogram
                                h2q12c->Fill(ex12c/1000.0, (hit[k].GetEnergy() + hit[l].GetEnergy() + hit[j].GetEnergy())/1000.0);
                                
                                // Give a good 12C if the energy is a Hoyle state (Excitation Energy = 7.6542 MeV)
                                // Peak cut is 7.6 MeV and 7.68 MeV
                                if (ex12c > 7600.0 && ex12c < 7680.0) {
                                    n12ca1 = k;
                                    n12ca2 = l;
                                    n12ca3 = j;
                                    igood12C++;
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }
        
        // Reconstruct 16O from 12C
        if (igood12C == 1 && igood > 3) {
            for (Int_t i = 0; i < igood; i++) {
                
                // Only use particles that aren't used in 12C reconstruction
                if(i != n12ca1 && i != n12ca2 && i != n12ca3){
                    // Sum of momenta
                    p16ox2 = p12cx + hit[i].GetMomentumX();
                    p16oy2 = p12cy + hit[i].GetMomentumY();
                    p16oz2 = p12cz + hit[i].GetMomentumZ();
                    //Energies
                    e16o2 = (p16ox2*p16ox2 + p16oy2*p16oy2 + p16oz2*p16oz2)/2.0/15.9949146;     // KE of 16O
                    q16o2 = e12c + hit[i].GetEnergy() - e16o2;                                  // Measured Q-value of reaction (Should be -7161.91698 keV)
                    ex16o2 = q16o2 + 7161.91698 + 7654.20;                                      // Excited state energies (ground + hoyle state)
                    
                    // Reconstruct neutron
                    pn2x = -p16ox2;
                    pn2y = -p16oy2;
                    pn2z = pbeam - p16oz2;
                    eneu2 = (pn2x*pn2x+pn2y*pn2y+pn2z*pn2z)/2.0/1.0086649;
                    
                    // Reaction Q-Value
                    etot2 = e12c + hit[i].GetEnergy() + eneu2;
                    qval2 = etot2 - ebeam;                            // 13C(a, 12C_hs + a + n) Q-value = -12600.508
                    h1qval2->Fill(qval2/1000.0);
                    if (qval2 < Q16Ohigh && qval2 > Q16Olow) {
                        h1ex16o2->Fill(ex16o2/1000);
                    }
                    
                    
                    // Reconstruct 13C
                    // Q-value = -4946.30581
                    p13cx = p12cx + pn2x;
                    p13cy = p12cy + pn2y;
                    p13cz = p12cz + pn2z;
                    e13c = (pow(p13cx, 2) + pow(p13cy, 2) + pow(p13cz, 2))/2.0/13.00335483778;
                    q13c = e12c + eneu2 - e13c; // Q-Value = -4946.30581
                    ex13c = q13c + 4946.30581;
                    
                    // Fill histograms
                    h1ex13c->Fill(ex13c/1000.0);
                    h2ex13c->Fill(ex13c/1000.0, (e12c + eneu2)/1000.0);
                    if(qval2 < Q16Ohigh && qval2 > Q16Olow) {
                        //h216Od1->Fill(ex16o2/1000.0, ex13c/1000.0);
                    }
                    
                    
                    // Reconstruct 5He
                    p5hex = hit[i].GetMomentumX() + pn2x;
                    p5hey = hit[i].GetMomentumY() + pn2y;
                    p5hez = hit[i].GetMomentumZ() + pn2z;
                    e5he = ((p5hex * p5hex) + (p5hey * p5hey) + (p5hez * p5hez))/2.0/5.012057224;
                    q5he = hit[i].GetEnergy() + eneu2 - e5he; // Q-Value = 735.00025 kev
                    ex5he = q5he - 735.00025;
                    
                    // Fill histogram
                    h1q5he->Fill(q5he/1000.0);
                    h2q5he->Fill(q5he/1000.0, (hit[i].GetEnergy() + eneu2)/1000.0);
                    if(qval2 < Q16Ohigh && qval2 > Q16Olow) {
                        //h216Od2->Fill(ex16o2/1000.0, q5he/1000.0);
                        
                        if (q5he > He5Gate0 && q5he < He5Gate2) {}
                        else {
                            h1ex16o_cut1->Fill(ex16o2/1000.0);
                        }
                        
                        if (q5he > He5Gate0 && q5he < He5Gate2) {}
                        else {
                            h1ex16o_cut2->Fill(ex16o2/1000.0);
                        }
                        
                        if (q5he > He5Gate0 && q5he < He5Gate3) {}
                        else {
                            h1ex16o_cut3->Fill(ex16o2/1000.0);
                        }
                        
                        if (q5he > He5Gate0 && q5he < He5Gate4) {}
                        else {
                            h1ex16o_cut4->Fill(ex16o2/1000.0);
                        }
                    }
                }
            }
        }
        
        h1m->Fill(igood);
    }
    
    file.Write();
    file.Close();

    //ftext.close(); 
    
    clock_t stop=clock();
    printf("Analysis Run Time: %10.3f (seconds) \n",(double) (stop-start)/CLOCKS_PER_SEC);
}
