/*******************************************
 * The main macro for 13C(a,n)16O* analysis.
 * Usage with root (see also run_all.C):
 *       .L libExpEvent.so
 *       .L analysis4.5.C
 *       Run(ENERGY)
 *       -> Submit as a batch job with qsub
 *       -> Using the file run_all.C
 * Updated 28 July 2017
 * This file is to be used in determining
 * the angular distribution of the reaction
 * Specifically to implement WPT corrections
 * to reconstruction of 12C+a
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
#include "functions1.7.h" // Update matching algorithm to correct matching of low energy signals when detector is swamped, thereby reducing background
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
    TFile file("/afs/crc.nd.edu/group/nsl/nuc/users/bfrentz/analysis/output/data/finalAngDist/twoDegree/data" + TString(stm.str()) + "MeV_correlation.root","RECREATE");

    TChain chain("evtTree");
    
    // Initialize values
    Double_t ebeamMeV = 0.0;

    Double_t Be9Gate1 = 0.0;
    Double_t Be9Gate2 = 0.0;
    Double_t Be9Gate3 = 0.0;
    Double_t Be9Gate4 = 0.0;

    Double_t Be9Gate_best = 0.0;

    Double_t Be9Gateb11 = 0.0;
    Double_t Be9Gateb12 = 0.0;
    Double_t Be9Gateb13 = 0.0;
    Double_t Be9Gateb14 = 0.0;
    Double_t Be9Gateb15 = 0.0;
    Double_t Be9Gateb16 = 0.0;
    Double_t Be9Gateb17 = 0.0;
    Double_t Be9Gateb18 = 0.0;
    Double_t Be9Gateb19 = 0.0;
    Double_t Be9Gateb20 = 0.0;
    
    Double_t He5Gate0 = 0.0;
    Double_t He5Gate1 = 0.0;
    Double_t He5Gate2 = 0.0;
    Double_t He5Gate3 = 0.0;
    Double_t He5Gate4 = 0.0;
    Double_t He5Gate5 = 0.0;
    Double_t He5Gate6 = 0.0;
    Double_t He5Gate7 = 0.0;

    Double_t C13Gate1 = 0.0;
    Double_t C13Gate2 = 0.0;
    Double_t C13Gate3 = 0.0;
    Double_t C13Gate4 = 0.0;
    Double_t C13Gate5 = 0.0;

    Double_t cGateCombo0 = 0.0;
    Double_t cGateCombo1 = 0.0;
    Double_t cGateCombo2 = 0.0;
    Double_t hGateCombo00 = 0.0;
    Double_t hGateCombo01 = 0.0;
    Double_t hGateCombo02 = 0.0;
    Double_t hGateCombo03 = 0.0;
    Double_t hGateCombo10 = 0.0;
    Double_t hGateCombo11 = 0.0;
    Double_t hGateCombo12 = 0.0;
    Double_t hGateCombo13 = 0.0;
    Double_t hGateCombo20 = 0.0;
    Double_t hGateCombo21 = 0.0;
    Double_t hGateCombo22 = 0.0;
    Double_t hGateCombo23 = 0.0;

    Double_t He5Gate_best = 0.0;
    Double_t C13Gate_best = 0.0;
    
    Double_t Q8BeMin = 0.0;
    Double_t Q8BeMax = 0.0;
    Double_t Q16OMin = 0.0;
    Double_t Q16OMax = 0.0;
    Double_t C12Gate1 = 0.0;
    Double_t C12Gate2 = 0.0;
    Double_t Q16Olow = 0.0;
    Double_t Q16Ohigh = 0.0;
    
    Double_t LBe_00Min = 14.5;
    Double_t LBe_00Max = 14.9;
    Double_t LBe_0Min = 0.0;
    Double_t LBe_0Max = 0.0;
    Double_t LBe_1Min = 0.0;
    Double_t LBe_1Max = 0.0;
    Double_t LBe_2Min = 0.0;
    Double_t LBe_2Max = 0.0;
    Double_t LBe_3Min = 0.0;
    Double_t LBe_3Max = 0.0;
    Double_t LBe_4Min = 0.0;
    Double_t LBe_4Max = 0.0;
    Double_t LBe_5Min = 0.0;
    Double_t LBe_5Max = 0.0;
    Double_t LBe_6Min = 0.0;
    Double_t LBe_6Max = 0.0;
    Double_t LBe_7Min = 0.0;
    Double_t LBe_7Max = 0.0;
    Double_t LBe_8Min = 0.0;
    Double_t LBe_8Max = 0.0;
    Double_t LBe_9Min = 0.0;
    Double_t LBe_9Max = 0.0;
    Double_t LBe_10Min = 0.0;
    Double_t LBe_10Max = 0.0;

    Double_t LHe_0Min = 0.0;
    Double_t LHe_0Max = 0.0;
    Double_t LHe_1Min = 0.0;
    Double_t LHe_1Max = 0.0;
    Double_t LHe_2Min = 0.0;
    Double_t LHe_2Max = 0.0;
    Double_t LHe_3Min = 0.0;
    Double_t LHe_3Max = 0.0;
    Double_t LHe_4Min = 0.0;
    Double_t LHe_4Max = 0.0;
    Double_t LHe_5Min = 0.0;
    Double_t LHe_5Max = 0.0;
    Double_t LHe_6Min = 0.0;
    Double_t LHe_6Max = 0.0;
    Double_t LHe_7Min = 0.0;
    Double_t LHe_7Max = 0.0;
    Double_t LHe_8Min = 0.0;
    Double_t LHe_8Max = 0.0;
    Double_t LHe_9Min = 0.0;
    Double_t LHe_9Max = 0.0;
    Double_t LHe_10Min = 0.0;
    Double_t LHe_10Max = 0.0;

    Double_t LC_0Min = 0.0;
    Double_t LC_0Max = 0.0;
    Double_t LC_1Min = 0.0;
    Double_t LC_1Max = 0.0;
    Double_t LC_2Min = 0.0;
    Double_t LC_2Max = 0.0;
    Double_t LC_3Min = 0.0;
    Double_t LC_3Max = 0.0;
    Double_t LC_4Min = 0.0;
    Double_t LC_4Max = 0.0;
    Double_t LC_5Min = 0.0;
    Double_t LC_5Max = 0.0;
    Double_t LC_6Min = 0.0;
    Double_t LC_6Max = 0.0;
    Double_t LC_7Min = 0.0;
    Double_t LC_7Max = 0.0;
    Double_t LC_8Min = 0.0;
    Double_t LC_8Max = 0.0;
    Double_t LC_9Min = 0.0;
    Double_t LC_9Max = 0.0;
    Double_t LC_10Min = 0.0;
    Double_t LC_10Max = 0.0;
    
    Double_t match_cond = 210.0;

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
            
            Be9Gate1 = 2000.0;
            Be9Gate2 = 2300.0;
            Be9Gate3 = 2500.0;
            Be9Gate4 = 2750.0;

            Be9Gateb11 = 1800.0;
            Be9Gateb12 = 1850.0;
            Be9Gateb13 = 1900.0;
            Be9Gateb14 = 1950.0;
            Be9Gateb15 = 2050.0;
            Be9Gateb16 = 2100.0;
            Be9Gateb17 = 2150.0;
            Be9Gateb18 = 2200.0;
            Be9Gateb19 = 2250.0;
            Be9Gateb20 = 2300.0;
            
            Q8BeMin = 71.0;
            Q8BeMax = 116.0;
            Q16OMin = -13667.0;
            Q16OMax = -10912.0;

            C12Gate1 = 7571.0;
            C12Gate2 = 7737.0;
            Q16Olow = -13789.0;
            Q16Ohigh = -11102.0;
            
            He5Gate1 = 750.0;
            He5Gate2 = 850.0;
            He5Gate3 = 950.0;
            He5Gate4 = 1050.0;
            He5Gate5 = 1150.0;
            He5Gate6 = 1250.0;
            He5Gate7 = 1350.0;

            C13Gate1 = 5800.0;
            C13Gate2 = 5900.0;
            C13Gate3 = 6000.0;
            C13Gate4 = 6100.0;
            C13Gate5 = 6200.0;

            cGateCombo0 = 5500.0;
    		cGateCombo1 = 5700.0;
            cGateCombo2 = 5900.0;
    		hGateCombo00 = 650.0;
    		hGateCombo01 = 750.0;
    		hGateCombo02 = 850.0;
            hGateCombo03 = 950.0;
            hGateCombo10 = 650.0;
    		hGateCombo11 = 750.0;
    		hGateCombo12 = 850.0;
    		hGateCombo13 = 950.0;
            hGateCombo20 = 650.0;
            hGateCombo21 = 750.0;
            hGateCombo22 = 850.0;
            hGateCombo23 = 950.0;

            Be9Gate_best = 2100;
            He5Gate_best = 650;         // Previously 850
            C13Gate_best = 5500;
            
            LBe_2Min = 17.025;
            LBe_2Max = 17.325;
            LBe_4Min = 17.675;
            LBe_4Max = 18.625;
            LBe_5Min = 18.675;
            LBe_5Max = 19.025;
            LBe_6Min = 19.175;
            LBe_6Max = 19.625;

            LHe_1Min = 16.225;
            LHe_1Max = 16.625;
            LHe_2Min = 17.025;
            LHe_2Max = 17.325;
            LHe_3Min = 17.450;
            LHe_3Max = 17.825;
            LHe_4Min = 17.825;
            LHe_4Max = 18.255;
            LHe_5Min = 18.625;
            LHe_5Max = 18.875;
            LHe_6Min = 18.525;
            LHe_6Max = 20.525;

            LC_1Min = 16.225;
            LC_1Max = 16.625;
            LC_2Min = 17.025;
            LC_2Max = 17.325;
            LC_3Min = 17.450;
            LC_3Max = 17.825;
            LC_4Min = 17.825;
            LC_4Max = 18.375;
            LC_5Min = 18.625;
            LC_5Max = 18.875;
            LC_6Min = 19.125;
            LC_6Max = 19.525;
            
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

            Be9Gate1 = 2000.0;
            Be9Gate2 = 2500.0;
            Be9Gate3 = 3000.0;
            Be9Gate4 = 3800.0;

            Be9Gateb11 = 2800.0;
            Be9Gateb12 = 2900.0;
            Be9Gateb13 = 3000.0;
            Be9Gateb14 = 3100.0;
            Be9Gateb15 = 3200.0;
            Be9Gateb16 = 3300.0;
            Be9Gateb17 = 3400.0;
            Be9Gateb18 = 3500.0;
            Be9Gateb19 = 3600.0;
            Be9Gateb20 = 3700.0;
            
            Q8BeMin = 71.0;
            Q8BeMax = 116.0;
            Q16OMin = -13839.0;
            Q16OMax = -10786.0;

            C12Gate1 = 7546.0;
            C12Gate2 = 7763.0;
            Q16Olow = -13801.0;
            Q16Ohigh = -10934.0;
            
            He5Gate1 = 500.0;
            He5Gate2 = 750.0;
            He5Gate3 = 1000.0;
            He5Gate4 = 1250.0;
            He5Gate5 = 1500.0;
            He5Gate6 = 1750.0;
            He5Gate7 = 2000.0;

            C13Gate1 = 5800.0;
            C13Gate2 = 6000.0;
            C13Gate3 = 6200.0;
            C13Gate4 = 6400.0;
            C13Gate5 = 6600.0;

            cGateCombo0 = 5600.0;
    		cGateCombo1 = 5800.0;
            cGateCombo2 = 6000.0;
    		hGateCombo00 = 750.0;
    		hGateCombo01 = 1050.0;
    		hGateCombo02 = 1350.0;
            hGateCombo03 = 1750.0;
            hGateCombo10 = 750.0;
            hGateCombo11 = 1050.0;
            hGateCombo12 = 1350.0;
            hGateCombo13 = 1750.0;
            hGateCombo20 = 750.0;
            hGateCombo21 = 1050.0;
            hGateCombo22 = 1350.0;
            hGateCombo23 = 1750.0;

            Be9Gate_best = 2800;
            He5Gate_best = 750;         // Previously 1750
            C13Gate_best = 6000;        // Previously 6200
            
            LBe_2Min = 16.975;
            LBe_2Max = 17.375;
            LBe_4Min = 17.975;
            LBe_4Max = 18.225;
            LBe_5Min = 18.325;
            LBe_5Max = 19.125;
            LBe_6Min = 19.125;
            LBe_6Max = 19.875;

            LHe_0Min = 15.225;
            LHe_0Max = 15.475;
            LHe_1Min = 16.375;
            LHe_1Max = 16.625;
            LHe_2Min = 16.925;
            LHe_2Max = 17.325;
            LHe_3Min = 17.375;
            LHe_3Max = 17.775;
            LHe_4Min = 17.825;
            LHe_4Max = 18.275;
            LHe_5Min = 18.725;
            LHe_5Max = 18.875;
            LHe_6Min = 19.275;
            LHe_6Max = 19.875;

            LC_0Min = 15.225;
            LC_0Max = 15.475;
            LC_1Min = 16.375;
            LC_1Max = 16.625;
            LC_2Min = 16.925;
            LC_2Max = 17.325;
            LC_3Min = 17.525;
            LC_3Max = 17.775;
            LC_4Min = 17.925;
            LC_4Max = 18.275;
            LC_5Min = 18.725;
            LC_5Max = 18.875;
            LC_6Min = 19.275;
            LC_6Max = 19.675;
            
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
            
            Be9Gate1 = 1000.0;
            Be9Gate2 = 2500.0;
            Be9Gate3 = 3000.0;
            Be9Gate4 = 3500.0;

            Be9Gateb11 = 1200.0;
            Be9Gateb12 = 1400.0;
            Be9Gateb13 = 1600.0;
            Be9Gateb14 = 1800.0;
            Be9Gateb15 = 2000.0;
            Be9Gateb16 = 3100.0;
            Be9Gateb17 = 3200.0;
            Be9Gateb18 = 3300.0;
            Be9Gateb19 = 3400.0;
            Be9Gateb20 = 3500.0;
            
            Q8BeMin = 71.0;
            Q8BeMax = 116.0;
            Q16OMin = -13834.0;
            Q16OMax = -10681.0;

            C12Gate1 = 7592.0;
            C12Gate2 = 7712.0;
            Q16Olow = -14015.0;
            Q16Ohigh = -10832.0;
            
            He5Gate1 = 700.0;
            He5Gate2 = 900.0;
            He5Gate3 = 1200.0;
            He5Gate4 = 1400.0;
            He5Gate5 = 1600.0;
            He5Gate6 = 1800.0;
            He5Gate7 = 2000.0;

            C13Gate1 = 5200.0;
            C13Gate2 = 5400.0;
            C13Gate3 = 5600.0;
            C13Gate4 = 5800.0;
            C13Gate5 = 6000.0;

            cGateCombo0 = 5350.0;
            cGateCombo1 = 5575.0;
    		cGateCombo2 = 5800.0;
    		hGateCombo00 = 500.0;
    		hGateCombo01 = 600.0;
    		hGateCombo02 = 700.0;
            hGateCombo03 = 800.0;
            hGateCombo10 = 500.0;
            hGateCombo11 = 600.0;
            hGateCombo12 = 700.0;
            hGateCombo13 = 800.0;
            hGateCombo20 = 500.0;
            hGateCombo21 = 600.0;
            hGateCombo22 = 700.0;
            hGateCombo23 = 800.0;

            Be9Gate_best = 2100;
            He5Gate_best = 700;
            C13Gate_best = 5800;
            
            LBe_2Min = 17.025;
            LBe_2Max = 17.325;
            LBe_3Min = 17.375;
            LBe_3Max = 17.575;
            LBe_4Min = 17.825;
            LBe_4Max = 18.375;
            LBe_5Min = 18.675;
            LBe_5Max = 18.975;
            LBe_6Min = 19.125;
            LBe_6Max = 19.675;
            LBe_7Min = 20.625;
            LBe_7Max = 21.475;
            LBe_8Min = 21.525;
            LBe_8Max = 22.125;

            LHe_1Min = 16.325;
            LHe_1Max = 16.675;
            LHe_2Min = 16.975;
            LHe_2Max = 17.325;
            LHe_3Min = 17.325;
            LHe_3Max = 17.825;
            LHe_4Min = 17.975;
            LHe_4Max = 18.325;
            LHe_5Min = 18.575;
            LHe_5Max = 19.125;
            LHe_6Min = 19.325;
            LHe_6Max = 19.675;
            LHe_7Min = 20.875;
            LHe_7Max = 21.175;
            LHe_8Min = 21.575;
            LHe_8Max = 21.975;
            LHe_9Min = 21.975;
            LHe_9Max = 22.725;

            LC_1Min = 16.325;
            LC_1Max = 16.675;
            LC_2Min = 16.975;
            LC_2Max = 17.325;
            LC_3Min = 17.325;
            LC_3Max = 17.825;
            LC_4Min = 17.975;
            LC_4Max = 18.325;
            LC_5Min = 18.575;
            LC_5Max = 19.125;
            LC_6Min = 19.325;
            LC_6Max = 19.675;
            LC_7Min = 20.825;
            LC_7Max = 21.275;
            LC_8Min = 21.575;
            LC_8Max = 21.975;
            
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

            Be9Gate1 = 2000.0;
            Be9Gate2 = 2500.0;
            Be9Gate3 = 3000.0;
            Be9Gate4 = 3500.0;

            Be9Gateb11 = 2300.0;
            Be9Gateb12 = 2400.0;
            Be9Gateb13 = 2500.0;
            Be9Gateb14 = 2600.0;
            Be9Gateb15 = 2700.0;
            Be9Gateb16 = 2800.0;
            Be9Gateb17 = 2900.0;
            Be9Gateb18 = 3000.0;
            Be9Gateb19 = 3100.0;
            Be9Gateb20 = 3200.0;
            
            Q8BeMin = 71.0;
            Q8BeMax = 117.0;
            Q16OMin = -13890.0;
            Q16OMax = -10782.0;

            C12Gate1 = 7581.0;
            C12Gate2 = 7717.0;
            Q16Olow = -14045.0;
            Q16Ohigh = -10781.0;
            
            He5Gate1 = 1200.0;
            He5Gate2 = 1300.0;
            He5Gate3 = 1400.0;
            He5Gate4 = 1500.0;
            He5Gate5 = 1600.0;
            He5Gate6 = 1700.0;
            He5Gate7 = 1800.0;

            C13Gate1 = 5600.0;
            C13Gate2 = 5700.0;
            C13Gate3 = 5800.0;
            C13Gate4 = 5900.0;
            C13Gate5 = 6000.0;

            cGateCombo0 = 5400.0;
    		cGateCombo1 = 5500.0;
            cGateCombo2 = 5600.0;
    		hGateCombo00 = 700.0;
    		hGateCombo01 = 900.0;
    		hGateCombo02 = 1100.0;
            hGateCombo03 = 1300.0;
            hGateCombo10 = 700.0;
            hGateCombo11 = 900.0;
            hGateCombo12 = 1100.0;
            hGateCombo13 = 1300.0;
            hGateCombo20 = 700.0;
            hGateCombo21 = 900.0;
            hGateCombo22 = 1100.0;
            hGateCombo23 = 1300.0;

            Be9Gate_best = 3100;
            He5Gate_best = 1200;
            C13Gate_best = 5600;
            
            LBe_2Min = 16.875;
            LBe_2Max = 17.3249;
            LBe_3Min = 17.325;
            LBe_3Max = 17.5249;
            LBe_4Min = 17.525;
            LBe_4Max = 18.375;
            LBe_5Min = 18.625;
            LBe_5Max = 19.025;
            LBe_6Min = 19.125;
            LBe_6Max = 19.675;
            LBe_7Min = 20.775;
            LBe_7Max = 21.125;
            LBe_8Min = 21.325;
            LBe_8Max = 21.825;

            LHe_2Min = 17.025;
            LHe_2Max = 17.325;
            LHe_3Min = 17.375;
            LHe_3Max = 17.775;
            LHe_4Min = 17.925;
            LHe_4Max = 18.325;
            LHe_5Min = 18.425;
            LHe_5Max = 19.275;
            LHe_6Min = 19.325;
            LHe_6Max = 19.675;
            LHe_7Min = 20.775;
            LHe_7Max = 21.225;
            LHe_8Min = 21.525;
            LHe_8Max = 21.925;
            LHe_9Min = 21.925;
            LHe_9Max = 22.375;
            LHe_10Min = 22.775;
            LHe_10Max = 24.175;

            LC_2Min = 17.025;
            LC_2Max = 17.325;
            LC_3Min = 17.375;
            LC_3Max = 17.775;
            LC_4Min = 17.925;
            LC_4Max = 18.325;
            LC_5Min = 18.425;
            LC_5Max = 19.275;
            LC_6Min = 19.325;
            LC_6Max = 19.675;
            LC_7Min = 20.775;
            LC_7Max = 21.225;
            LC_8Min = 21.525;
            LC_8Max = 21.925;
            LC_9Min = 21.925;
            LC_9Max = 22.375;
            LC_10Min = 22.575;
            LC_10Max = 24.125;
            
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

            Be9Gate1 = 2100.0;
            Be9Gate2 = 3400.0;
            Be9Gate3 = 3700.0;
            Be9Gate4 = 4000.0;

            Be9Gateb11 = 1500.0;
            Be9Gateb12 = 1700.0;
            Be9Gateb13 = 1900.0;
            Be9Gateb14 = 2100.0;
            Be9Gateb15 = 2300.0;
            Be9Gateb16 = 2500.0;
            Be9Gateb17 = 2700.0;
            Be9Gateb18 = 2900.0;
            Be9Gateb19 = 3100.0;
            Be9Gateb20 = 3300.0;
            
            Q8BeMin = 70.0;
            Q8BeMax = 117.0;
            Q16OMin = -14230.0;
            Q16OMax = -10525.0;

            C12Gate1 = 7577.0;
            C12Gate2 = 7719.0;
            Q16Olow = -14147.0;
            Q16Ohigh = -10721.0;
            
            He5Gate1 = 300.0;
            He5Gate2 = 500.0;
            He5Gate3 = 700.0;
            He5Gate4 = 900.0;
            He5Gate5 = 1100.0;
            He5Gate6 = 1300.0;
            He5Gate7 = 1500.0;

            C13Gate1 = 5500.0;
            C13Gate2 = 5700.0;
            C13Gate3 = 5900.0;
            C13Gate4 = 6100.0;
            C13Gate5 = 6300.0;

            cGateCombo0 = 5300.0;
    		cGateCombo1 = 5500.0;
            cGateCombo2 = 5700.0;
    		hGateCombo00 = 150.0;
    		hGateCombo01 = 300.0;
    		hGateCombo02 = 450.0;
            hGateCombo03 = 600.0;
            hGateCombo10 = 150.0;
            hGateCombo11 = 300.0;
            hGateCombo12 = 450.0;
            hGateCombo13 = 600.0;
            hGateCombo20 = 150.0;
            hGateCombo21 = 300.0;
            hGateCombo22 = 450.0;
            hGateCombo23 = 600.0;

            Be9Gate_best = 2900;
            He5Gate_best = 300;
            C13Gate_best = 5500;
            
            LBe_2Min = 16.975;
            LBe_2Max = 17.325;
            LBe_4Min = 17.625;
            LBe_4Max = 18.375;
            LBe_5Min = 18.575;
            LBe_5Max = 18.925;
            LBe_6Min = 19.125;
            LBe_6Max = 19.675;
            LBe_7Min = 20.775;
            LBe_7Max = 21.3749;
            LBe_8Min = 21.375;
            LBe_8Max = 22.125;
            LBe_10Min = 22.775;
            LBe_10Max = 23.225;

            LHe_2Min = 17.025;
            LHe_2Max = 17.375;
            LHe_3Min = 17.375;
            LHe_3Max = 17.775;
            LHe_4Min = 17.925;
            LHe_4Max = 18.325;
            LHe_5Min = 18.425;
            LHe_5Max = 19.275;
            LHe_6Min = 19.175;
            LHe_6Max = 19.775;
            LHe_7Min = 20.775;
            LHe_7Max = 21.225;
            LHe_8Min = 21.375;
            LHe_8Max = 21.925;
            LHe_9Min = 21.925;
            LHe_9Max = 22.375;
            LHe_10Min = 22.775;
            LHe_10Max = 24.175;

            LC_2Min = 17.025;
            LC_2Max = 17.375;
            LC_3Min = 17.375;
            LC_3Max = 17.775;
            LC_4Min = 17.925;
            LC_4Max = 18.325;
            LC_5Min = 18.425;
            LC_5Max = 19.275;
            LC_6Min = 19.175;
            LC_6Max = 19.775;
            LC_7Min = 20.775;
            LC_7Max = 21.225;
            LC_8Min = 21.375;
            LC_8Max = 21.925;
            LC_9Min = 21.925;
            LC_9Max = 22.375;
            LC_10Min = 22.775;
            LC_10Max = 24.175;
            
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

        // For the alpha source runs
        case 1: {
            // These gates won't matter for the alpha source runs, so I copied the data from 24MeV
            // this is the case because they are only used in the reconstructions, which we don't 
            // have for the alpha source
            ebeamMeV = 1.0;
            
            Be9Gate1 = 2400.0;
            Be9Gate2 = 2000.0;
            
            Q8BeMin = 76.0;
            Q8BeMax = 104.0;
            Q16OMin = -13500.0;
            Q16OMax = -11600.0;

            C12Gate1 = 7600.0;
            C12Gate2 = 7680.0;
            Q16Olow = -13632.0;
            Q16Ohigh = -11621.0;
            
            He5Gate1 = 1500.0;
            He5Gate2 = 3000.0;
            He5Gate3 = 4750.0;
            
            
            // Add data to the chain
            fRuns.open("/afs/crc.nd.edu/group/nsl/nuc/users/bfrentz/analysis/input/source.in", ios::in);
            fRuns.seekg(0, ios::beg);
            fRuns >> nRuns;
            fRuns.seekg(1, ios::cur);
            for(int i = 0; i < nRuns; i++){
                fRuns >> run_num;
                sprintf(dataf,"/afs/crc.nd.edu/group/nsl/nuc/users/wtan/2015Mar_16O/evt2root/data/asic/root/run%d.root", run_num);
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
    TH1F *h1errSingle = new TH1F("h1errSingle","Energy Spread (matched single hits)",6000, -600.0,600.0);
    TH1F *h1errDouble = new TH1F("h1errDouble","Energy Spread (matched double hits)",6000, -600.0,600.0);
    TH1F *h1chi2 = new TH1F("h1chi2","Chi square distribution",1000,0,1000);
    TH2F *hit2d = new TH2F("hit2d","Matched Hit Pattern Theta vs Phi",1110,-185.0,185.0,315,0.0,90.0);
    TH2F *hitxy = new TH2F("hitxy","Matched Hit Pattern Y vs X",1300,-130.0,130.0,144,-36.0,36.0);
    TH2F *hitec = new TH2F("hitec","asic channel vs calibrated_E",290,0.0,290.0,8191,0.5,8191.5);
    TH2F *h2efb = new TH2F("h2efb","EF_strip# vs EB_strip# (matched hits)",130,-1.0,129.0,36,-1.0,35.0);
    TH2F *h2epix = new TH2F("h2epix","Calibrated Energy vs Pixel (matched)",4100,-1.0,4099.0,900,0.5,9000.5);

    TH2F *hitxy2 = new TH2F("hitxy2","Matched Hit Pattern Y vs X for checking",1300,-130.0,130.0,144,-36.0,36.0);
    
    TH1F *h1qval = new TH1F("h1qval","Q-value [MeV] of 13C(a,n+8Be+8Be)",1600, -16.0,0.0);
    TH1F *h1ex16o = new TH1F("h1ex16o","Ex [MeV] in 16O reconstructed from 2 8Be",240, 14.0,26.0);
    TH1F *h1ex16ogated = new TH1F("h1ex16ogated","Ex [MeV] in 16O with 9Be cut1",240, 14.0,26.0);
    TH1F *h1ex16ogated2 = new TH1F("h1ex16ogated2","Ex [MeV] in 16O with 9Be cut2",240, 14.0,26.0);
    TH1F *h1ex16ogated3 = new TH1F("h1ex16ogated3","Ex [MeV] in 16O with 9Be cut3",240, 14.0,26.0);
    TH1F *h1ex16ogated4 = new TH1F("h1ex16ogated4","Ex [MeV] in 16O with 9Be cut4",240, 14.0,26.0);
    TH1F *h1ex16ogatedN = new TH1F("h1ex16ogatedN","Ex [MeV] in 16O",240, 14.0,26.0);

    TH1F *h1ex16ogated_b11 = new TH1F("h1ex16ogated_b11","Ex [MeV] in 16O with 9Be cut11",240, 14.0,26.0);
    TH1F *h1ex16ogated_b12 = new TH1F("h1ex16ogated_b12","Ex [MeV] in 16O with 9Be cut12",240, 14.0,26.0);
    TH1F *h1ex16ogated_b13 = new TH1F("h1ex16ogated_b13","Ex [MeV] in 16O with 9Be cut13",240, 14.0,26.0);
    TH1F *h1ex16ogated_b14 = new TH1F("h1ex16ogated_b14","Ex [MeV] in 16O with 9Be cut14",240, 14.0,26.0);
    TH1F *h1ex16ogated_b15 = new TH1F("h1ex16ogated_b15","Ex [MeV] in 16O with 9Be cut15",240, 14.0,26.0);
    TH1F *h1ex16ogated_b16 = new TH1F("h1ex16ogated_b16","Ex [MeV] in 16O with 9Be cut16",240, 14.0,26.0);
    TH1F *h1ex16ogated_b17 = new TH1F("h1ex16ogated_b17","Ex [MeV] in 16O with 9Be cut17",240, 14.0,26.0);
    TH1F *h1ex16ogated_b18 = new TH1F("h1ex16ogated_b18","Ex [MeV] in 16O with 9Be cut18",240, 14.0,26.0);
    TH1F *h1ex16ogated_b19 = new TH1F("h1ex16ogated_b19","Ex [MeV] in 16O with 9Be cut19",240, 14.0,26.0);
    TH1F *h1ex16ogated_b20 = new TH1F("h1ex16ogated_b20","Ex [MeV] in 16O with 9Be cut20",240, 14.0,26.0);

    TH1F *h1q8be = new TH1F("h1q8be","Q-value [keV] of 8Be decay",500, -100.0,350.0);
    TH2F *h2q8be = new TH2F("h2q8be","Sum_Energy vs. Q-value [keV] of 8Be decay",500,-100.0,350.0,1600,0.5,32000.5);
    
    TH1F *h1q8be_restrict = new TH1F("h1q8be_restrict","Q-value [keV] of 8Be decay", 500, 50.0, 150.0);

    TH1F* h1q9be = new TH1F("h1q9be","Excitation Energy [keV] of 9Be decay to two alphas",1000,0.0,4000.0);
    TH2F* h2q9be = new TH2F("h2q9be","Sum_Energy vs Excitation Energy [keV] of 9Be decay", 1000, 0.0, 4000.0,1600,0.5, 32000.5);
    
    TH1F *h1q12c = new TH1F("h1q12c", "Q-value [MeV] of 12C decay to three alphas", 5000, 0.0, 30.0);
    TH1F *h1ex12c = new TH1F("h1ex12c", "Excitation Energy [MeV] of 12C decay to three alphas", 3500, 0, 35.0);
    TH1F *h1erel12c = new TH1F("h1erel12c", "12C decay to three alphas, Erel", 750, 0, 15.0);
    TH1F *h1ex12c4 = new TH1F("h1ex12c4", "Excitation Energy [MeV] of 12C decay to three alphas requiring 4 matched hits", 1750, 0, 35.0);
    TH1F *h1ex12c2 = new TH1F("h1ex12c2", "Excitation Energy [MeV] of 12C decay to three alphas, requiring 4 matched hits on reaction q-value", 3500, 0, 35.0);
    TH2F *h2q12c = new TH2F("h2q12c", "Sum_Energy vs. Excitation Energy [MeV] of 12C decay (in 16O reconstruction)", 7000, 0, 35.0, 8000, 0, 40.0);
    
    TH1F *h1qval2 = new TH1F("h1qval2","Q-value [MeV] of 13C(a, 12C_hs + a + n)",1500, -15.0,0.0);
    TH1F *h1qval4a = new TH1F("h1qval4a","Q-value [MeV] of 13C(a, 4a + n)",1500, -15.0,0.0);
    TH1F *h1ex16o2 = new TH1F("h1ex16o2","Excitation Energy [MeV] in 16O Reconstructed from 12C + a", 240, 14.0, 26.0);
    
    TH1F *h1ex16o_h1 = new TH1F("h1ex16o_h1","Excitation Energy [MeV] in 16O, First Helium5 Cut", 240, 14.0, 26.0);
    TH1F *h1ex16o_h2 = new TH1F("h1ex16o_h2","Excitation Energy [MeV] in 16O, Second Helium5 Cut", 240, 14.0, 26.0);
    TH1F *h1ex16o_h3 = new TH1F("h1ex16o_h3","Excitation Energy [MeV] in 16O, Third Helium5 Cut", 240, 14.0, 26.0);
    TH1F *h1ex16o_h4 = new TH1F("h1ex16o_h4","Excitation Energy [MeV] in 16O, Fourth Helium5 Cut", 240, 14.0, 26.0);
    TH1F *h1ex16o_h5 = new TH1F("h1ex16o_h5","Excitation Energy [MeV] in 16O, Fifth Helium5 Cut", 240, 14.0, 26.0);
    TH1F *h1ex16o_h6 = new TH1F("h1ex16o_h6","Excitation Energy [MeV] in 16O, Sixth Helium5 Cut", 240, 14.0, 26.0);
    TH1F *h1ex16o_h7 = new TH1F("h1ex16o_h7","Excitation Energy [MeV] in 16O, Seventh Helium5 Cut", 240, 14.0, 26.0);
    TH1F *h1ex16o_c1 = new TH1F("h1ex16o_c1","Excitation Energy [MeV] in 16O, First Carbon13 Cut", 240, 14.0, 26.0);
    TH1F *h1ex16o_c2 = new TH1F("h1ex16o_c2","Excitation Energy [MeV] in 16O, Second Carbon13 Cut", 240, 14.0, 26.0);
    TH1F *h1ex16o_c3 = new TH1F("h1ex16o_c3","Excitation Energy [MeV] in 16O, Third Carbon13 Cut", 240, 14.0, 26.0);
    TH1F *h1ex16o_c4 = new TH1F("h1ex16o_c4","Excitation Energy [MeV] in 16O, Fourth Carbon13 Cut", 240, 14.0, 26.0);
    TH1F *h1ex16o_c5 = new TH1F("h1ex16o_c5","Excitation Energy [MeV] in 16O, Fifth Carbon13 Cut", 240, 14.0, 26.0);

    TH1F *h1ex16o_combined00 = new TH1F("h1ex16o_combined00", "Excitation Energy [MeV] in 16O, combining 5He and 13C 00", 240, 14.0, 26.0);
    TH1F *h1ex16o_combined01 = new TH1F("h1ex16o_combined01", "Excitation Energy [MeV] in 16O, combining 5He and 13C 01", 240, 14.0, 26.0);
    TH1F *h1ex16o_combined02 = new TH1F("h1ex16o_combined02", "Excitation Energy [MeV] in 16O, combining 5He and 13C 02", 240, 14.0, 26.0);
    TH1F *h1ex16o_combined03 = new TH1F("h1ex16o_combined03", "Excitation Energy [MeV] in 16O, combining 5He and 13C 10", 240, 14.0, 26.0);
    TH1F *h1ex16o_combined10 = new TH1F("h1ex16o_combined10", "Excitation Energy [MeV] in 16O, combining 5He and 13C 10", 240, 14.0, 26.0);
    TH1F *h1ex16o_combined11 = new TH1F("h1ex16o_combined11", "Excitation Energy [MeV] in 16O, combining 5He and 13C 11", 240, 14.0, 26.0);
    TH1F *h1ex16o_combined12 = new TH1F("h1ex16o_combined12", "Excitation Energy [MeV] in 16O, combining 5He and 13C 12", 240, 14.0, 26.0);
    TH1F *h1ex16o_combined13 = new TH1F("h1ex16o_combined13", "Excitation Energy [MeV] in 16O, combining 5He and 13C 10", 240, 14.0, 26.0);
    TH1F *h1ex16o_combined20 = new TH1F("h1ex16o_combined20", "Excitation Energy [MeV] in 16O, combining 5He and 13C 10", 240, 14.0, 26.0);
    TH1F *h1ex16o_combined21 = new TH1F("h1ex16o_combined21", "Excitation Energy [MeV] in 16O, combining 5He and 13C 11", 240, 14.0, 26.0);
    TH1F *h1ex16o_combined22 = new TH1F("h1ex16o_combined22", "Excitation Energy [MeV] in 16O, combining 5He and 13C 12", 240, 14.0, 26.0);
    TH1F *h1ex16o_combined23 = new TH1F("h1ex16o_combined23", "Excitation Energy [MeV] in 16O, combining 5He and 13C 10", 240, 14.0, 26.0);

    TH1F *h1ex16o_bebest = new TH1F("h1ex16o_bebest","Excitation Energy [MeV] in 16O, Best Beryllium9 Cut", 240, 14.0, 26.0);
    TH1F *h1ex16o_hbest = new TH1F("h1ex16o_hbest","Excitation Energy [MeV] in 16O, Best Helium5 Cut", 240, 14.0, 26.0);
    TH1F *h1ex16o_cbest = new TH1F("h1ex16o_cbest","Excitation Energy [MeV] in 16O, Best Carbon13 Cut", 240, 14.0, 26.0);
    TH1F *h1ex16o_2best = new TH1F("h1ex16o_2best","Excitation Energy [MeV] in 16O, Best Combination of C13 and He5 Cut", 240, 14.0, 26.0);

    TH1F *h1ex16o_beAll = new TH1F("h1ex16o_beAll","Excitation Energy [MeV] in 16O", 240, 14.0, 26.0);
    TH1F *h1ex16o_combAll = new TH1F("h1ex16o_comAll","Excitation Energy [MeV] in 16O", 240, 14.0, 26.0);

    TH1F *h1_8be_corr = new TH1F("h1_8be_corr","Correlation function from Event Mixing w/ 8Be's", 240, 0.0, 12.0);
    TH1F *h1_12ca_corr = new TH1F("h1_12ca_corr","Correlation function from Event Mixing w/ 12C + a's", 240, 0.0, 12.0);
    TH1F *h1_12c_corr = new TH1F("h1_12c_corr","Correlation function of 12C from Event Mixing 3a's", 750, 0, 15.0);

    TH1F *h1ex16o_em = new TH1F("h1ex16o_em","Ex [MeV] in 16O from Event Mixing w/ 8Be's", 240, 0.0, 12.0);
    TH1F *h1ex16o_em1 = new TH1F("h1ex16o_em1","Ex [MeV] in 16O from Event Mixing w/ 12C + a's", 240, 0.0, 12.0);
    TH1F *h1ex12c_em = new TH1F("h1ex12c_em","Ex [MeV] in 12C from Event Mixing w/ 3a's", 750, 0, 15.0);
    TH1F *h1ex4a_c_em = new TH1F("h1ex4a_c_em", "Ex [MeV] in 12C from event mixing with 4a's", 750, 0, 15.0);
    TH1F *h1erel12c_em = new TH1F("h1erel12c_em", "12C decay to three alphas, Erel from event mixing", 750, 0, 15.0);

    TH1F *h1q_em = new TH1F("h1q_em","Q-value [MeV] of event mixing",1600, -16.0,0.0);
    TH1F *h1ex_em = new TH1F("h1ex_em","Excitation Energy [MeV] in 16O, Mixing tests", 240, 14.0, 26.0);
    TH1F *h1e1_em = new TH1F("h1e1_em","Excitation Energy E1 [MeV] in 16O, Mixing tests", 2000, -50.0, 50.0);
    TH1F *h1eneu_em = new TH1F("h1eneu_em","Excitation Energy Neutron [MeV] in 16O, Mixing tests", 2000, -50.0, 50.0);

    TH1F *h1qval4hit = new TH1F("h1qval4hit","Q-value [MeV] of 13C(a, 4a + n)",1500, -15.0,0.0);

    TH1F *h1ex13c = new TH1F("h1ex13c", "Excitation Energy [MeV] of 13C decay to 3a+n", 6000, 0, 30.0);
    TH2F *h2ex13c = new TH2F("h2ex13c", "Sum_Energy vs. Excitation Energy [MeV] of 13C decay", 6000, 0, 24.0, 6000, 0.0, 32.0);
    
    TH1F *h1q5he = new TH1F("h1q5he", "Q-Value [MeV] of 5He decay to a+n", 600, -2, 28.0);
    TH2F *h2q5he = new TH2F("h2q5he", "Sum_Energy vs. Q-Value [MeV] of 5He decay", 5000, 0, 24.0, 6000, -2, 28.0);
    
    TH2F* h216Od1 = new TH2F("h216Od1","16O vs 9Be DP1",240,14.,26.,600,0.,30.);
    TH2F* h216Od2 = new TH2F("h216Od2","16O vs 9Be DP2",240,14.,26.,600,0.,30.);
    TH2F* h216Od3 = new TH2F("h216Od3","16O vs 5He DP",240,14.,26.,600,0.,30.);
    TH2F* h216Od4 = new TH2F("h216Od4","16O vs 13C DP",240,14.,26.,600,0.,30.);

    TH2F* h28ben = new TH2F("h28ben","8Be1 + n (x-axis) vs 8Be2 + n", 180, 1.0, 10.0, 180, 1.0, 10.0);
    
    TH2F* hSpread0 = new TH2F("hSpread0", "E vs eSpread (keV) - Detector 0", 500, -1000.0, 1000, 500, 0.0, 30000.0);
    TH2F* hSpread1 = new TH2F("hSpread1", "E vs eSpread (keV) - Detector 1", 500, -1000.0, 1000, 500, 0.0, 30000.0);
    TH2F* hSpread2 = new TH2F("hSpread2", "E vs eSpread (keV) - Detector 2", 500, -1000.0, 1000, 500, 0.0, 30000.0);
    TH2F* hSpread3 = new TH2F("hSpread3", "E vs eSpread (keV) - Detector 3", 500, -1000.0, 1000, 500, 0.0, 30000.0);

    TH2F* hSpreaderr01 = new TH2F("hSpreaderr01", "Detector 0: Error Code 1 - Split back", 500, -1000.0, 1000, 500, 0.0, 30000.0);
    TH2F* hSpreaderr02 = new TH2F("hSpreaderr02", "Detector 0: Error Code 2 - Double back", 500, -1000.0, 1000, 500, 0.0, 30000.0);
    TH2F* hSpreaderr03 = new TH2F("hSpreaderr03", "Detector 0: Error Code 3 - Split front", 500, -1000.0, 1000, 500, 0.0, 30000.0);
    TH2F* hSpreaderr04 = new TH2F("hSpreaderr04", "Detector 0: Error Code 4 - Double front", 500, -1000.0, 1000, 500, 0.0, 30000.0);

    TH2F* hSpreaderr11 = new TH2F("hSpreaderr11", "Detector 1: Error Code 1 - Split back", 500, -1000.0, 1000, 500, 0.0, 30000.0);
    TH2F* hSpreaderr12 = new TH2F("hSpreaderr12", "Detector 1: Error Code 2 - Double back", 500, -1000.0, 1000, 500, 0.0, 30000.0);
    TH2F* hSpreaderr13 = new TH2F("hSpreaderr13", "Detector 1: Error Code 3 - Split front", 500, -1000.0, 1000, 500, 0.0, 30000.0);
    TH2F* hSpreaderr14 = new TH2F("hSpreaderr14", "Detector 1: Error Code 4 - Double front", 500, -1000.0, 1000, 500, 0.0, 30000.0);

    TH2F* hSpreaderr21 = new TH2F("hSpreaderr21", "Detector 2: Error Code 1 - Split back", 500, -1000.0, 1000, 500, 0.0, 30000.0);
    TH2F* hSpreaderr22 = new TH2F("hSpreaderr22", "Detector 2: Error Code 2 - Double back", 500, -1000.0, 1000, 500, 0.0, 30000.0);
    TH2F* hSpreaderr23 = new TH2F("hSpreaderr23", "Detector 2: Error Code 3 - Split front", 500, -1000.0, 1000, 500, 0.0, 30000.0);
    TH2F* hSpreaderr24 = new TH2F("hSpreaderr24", "Detector 2: Error Code 4 - Double front", 500, -1000.0, 1000, 500, 0.0, 30000.0);

    TH2F* hSpreaderr31 = new TH2F("hSpreaderr31", "Detector 3: Error Code 1 - Split back", 500, -1000.0, 1000, 500, 0.0, 30000.0);
    TH2F* hSpreaderr32 = new TH2F("hSpreaderr32", "Detector 3: Error Code 2 - Double back", 500, -1000.0, 1000, 500, 0.0, 30000.0);
    TH2F* hSpreaderr33 = new TH2F("hSpreaderr33", "Detector 3: Error Code 3 - Split front", 500, -1000.0, 1000, 500, 0.0, 30000.0);
    TH2F* hSpreaderr34 = new TH2F("hSpreaderr34", "Detector 3: Error Code 4 - Double front", 500, -1000.0, 1000, 500, 0.0, 30000.0);

    TH1F* hStrips = new TH1F("hStrips", "Counts vs Strip Identifier", 64, -0.5, 63.5);
    TH1F* hStrips2 = new TH1F("hStrips2", "Counts vs Strip Identifier", 64, -0.5, 63.5);
    
    // Angular Distributions
    //Peaks
    TH1F *h1angdistL00_Be = new TH1F("h1angdistL00_Be","Angular Distribution 14.6 MeV State (Be reconstruction)",90,0.5,180.5);
    TH1F *h1angdistL0_Be = new TH1F("h1angdistL0_Be","Angular Distribution 15.4 MeV State (Be reconstruction)",90,0.5,180.5);
    TH1F *h1angdistL1_Be = new TH1F("h1angdistL1_Be","Angular Distribution 16.5 MeV State (Be reconstruction)",90,0.5,180.5);
    TH1F *h1angdistL2_Be = new TH1F("h1angdistL2_Be","Angular Distribution 17.2 MeV State (Be reconstruction)",90,0.5,180.5);
    TH1F *h1angdistL3_Be = new TH1F("h1angdistL3_Be","Angular Distribution 17.5 MeV State (Be reconstruction)",90,0.5,180.5);
    TH1F *h1angdistL4_Be = new TH1F("h1angdistL4_Be","Angular Distribution 18.1 MeV State (Be reconstruction)",90,0.5,180.5);
    TH1F *h1angdistL5_Be = new TH1F("h1angdistL5_Be","Angular Distribution 18.8 MeV State (Be reconstruction)",90,0.5,180.5);
    TH1F *h1angdistL6_Be = new TH1F("h1angdistL6_Be","Angular Distribution 19.4 MeV State (Be reconstruction)",90,0.5,180.5);
    TH1F *h1angdistL7_Be = new TH1F("h1angdistL7_Be","Angular Distribution 21.0 MeV State (Be reconstruction)",90,0.5,180.5);
    TH1F *h1angdistL8_Be = new TH1F("h1angdistL8_Be","Angular Distribution 21.7 MeV State (Be reconstruction)",90,0.5,180.5);
    TH1F *h1angdistL9_Be = new TH1F("h1angdistL9_Be","Angular Distribution 22.3 MeV State (Be reconstruction)",90,0.5,180.5);
    TH1F *h1angdistL10_Be = new TH1F("h1angdistL10_Be","Angular Distribution 23.3 MeV State (Be reconstruction)",90,0.5,180.5);
    TH1F *h1angdistBG_Be = new TH1F("h1angdistBG_Be","Angular Distribution Background (Be reconstruction)",90,0.5,180.5);
    TH1F *h1angdist9Be = new TH1F("h1angdist9Be", "Angular Distribution 9Be Background", 90, 0.5, 180.5);

    //Backgrounds
    TH1F *h1angdistL0_BeB = new TH1F("h1angdistL0_BeB","Angular Distribution Below 15.4 MeV State (Be reconstruction)",90,0.5,180.5);
    TH1F *h1angdistL1_BeB = new TH1F("h1angdistL1_BeB","Angular Distribution Above 15.4 Below 16.5 MeV State (Be reconstruction)",90,0.5,180.5);
    TH1F *h1angdistL2_BeB = new TH1F("h1angdistL2_BeB","Angular Distribution Above 16.5 Below 17.2 MeV State (Be reconstruction)",90,0.5,180.5);
    TH1F *h1angdistL3_BeB = new TH1F("h1angdistL3_BeB","Angular Distribution Above 17.2 Below 17.5 MeV State (Be reconstruction)",90,0.5,180.5);
    TH1F *h1angdistL4_BeB = new TH1F("h1angdistL4_BeB","Angular Distribution Above 17.5 Below 18.1 MeV State (Be reconstruction)",90,0.5,180.5);
    TH1F *h1angdistL5_BeB = new TH1F("h1angdistL5_BeB","Angular Distribution Above 18.1 Below 18.8 MeV State (Be reconstruction)",90,0.5,180.5);
    TH1F *h1angdistL6_BeB = new TH1F("h1angdistL6_BeB","Angular Distribution Above 18.8 Below 19.4 MeV State (Be reconstruction)",90,0.5,180.5);
    TH1F *h1angdistL7_BeB = new TH1F("h1angdistL7_BeB","Angular Distribution Above 19.4 Below 21.0 MeV State (Be reconstruction)",90,0.5,180.5);
    TH1F *h1angdistL8_BeB = new TH1F("h1angdistL8_BeB","Angular Distribution Above 21.0 Below 21.7 MeV State (Be reconstruction)",90,0.5,180.5);
    TH1F *h1angdistL9_BeB = new TH1F("h1angdistL9_BeB","Angular Distribution Above 21.7 Below 22.3 MeV State (Be reconstruction)",90,0.5,180.5);
    TH1F *h1angdistL10_BeB = new TH1F("h1angdistL10_BeB","Angular Distribution Above 22.3 Below 23.3 MeV State (Be reconstruction)",90,0.5,180.5);
    TH1F *h1angdistL11_BeB = new TH1F("h1angdistL11_BeB","Angular Distribution Above 23.3 MeV State (Be reconstruction)",90,0.5,180.5);



    TH1F *h1EMangdistL0_Be = new TH1F("h1EMangdistL0_Be","Angular Distribution 15.4 MeV State (Event Mixing Be reconstruction)",90,0.5,180.5);
    TH1F *h1EMangdistL1_Be = new TH1F("h1EMangdistL1_Be","Angular Distribution 16.5 MeV State (Event Mixing Be reconstruction)",90,0.5,180.5);
    TH1F *h1EMangdistL2_Be = new TH1F("h1EMangdistL2_Be","Angular Distribution 17.2 MeV State (Event Mixing Be reconstruction)",90,0.5,180.5);
    TH1F *h1EMangdistL3_Be = new TH1F("h1EMangdistL3_Be","Angular Distribution 17.5 MeV State (Event Mixing Be reconstruction)",90,0.5,180.5);
    TH1F *h1EMangdistL4_Be = new TH1F("h1EMangdistL4_Be","Angular Distribution 18.1 MeV State (Event Mixing Be reconstruction)",90,0.5,180.5);
    TH1F *h1EMangdistL5_Be = new TH1F("h1EMangdistL5_Be","Angular Distribution 18.8 MeV State (Event Mixing Be reconstruction)",90,0.5,180.5);
    TH1F *h1EMangdistL6_Be = new TH1F("h1EMangdistL6_Be","Angular Distribution 19.4 MeV State (Event Mixing Be reconstruction)",90,0.5,180.5);
    TH1F *h1EMangdistL7_Be = new TH1F("h1EMangdistL7_Be","Angular Distribution 21.0 MeV State (Event Mixing Be reconstruction)",90,0.5,180.5);
    TH1F *h1EMangdistL8_Be = new TH1F("h1EMangdistL8_Be","Angular Distribution 21.7 MeV State (Event Mixing Be reconstruction)",90,0.5,180.5);
    TH1F *h1EMangdistL9_Be = new TH1F("h1EMangdistL9_Be","Angular Distribution 22.3 MeV State (Event Mixing Be reconstruction)",90,0.5,180.5);
    TH1F *h1EMangdistL10_Be = new TH1F("h1EMangdistL10_Be","Angular Distribution 23.3 MeV State (Event Mixing Be reconstruction)",90,0.5,180.5);

    // Backgrounds
    // Actually combined
    TH1F *h1angdistL0_HeB = new TH1F("h1angdistL0_HeB","Angular Distribution Below 15.4 MeV State (Combined cuts)",90,0.5,180.5);
    TH1F *h1angdistL1_HeB = new TH1F("h1angdistL1_HeB","Angular Distribution Above 15.4 Below 16.5 MeV State (Combined cuts)",90,0.5,180.5);
    TH1F *h1angdistL2_HeB = new TH1F("h1angdistL2_HeB","Angular Distribution Above 16.5 Below 17.2 MeV State (Combined cuts)",90,0.5,180.5);
    TH1F *h1angdistL3_HeB = new TH1F("h1angdistL3_HeB","Angular Distribution Above 17.2 Below 17.5 MeV State (Combined cuts)",90,0.5,180.5);
    TH1F *h1angdistL4_HeB = new TH1F("h1angdistL4_HeB","Angular Distribution Above 17.5 Below 18.1 MeV State (Combined cuts)",90,0.5,180.5);
    TH1F *h1angdistL5_HeB = new TH1F("h1angdistL5_HeB","Angular Distribution Above 18.1 Below 18.8 MeV State (Combined cuts)",90,0.5,180.5);
    TH1F *h1angdistL6_HeB = new TH1F("h1angdistL6_HeB","Angular Distribution Above 18.8 Below 19.4 MeV State (Combined cuts)",90,0.5,180.5);
    TH1F *h1angdistL7_HeB = new TH1F("h1angdistL7_HeB","Angular Distribution Above 19.4 Below 21.0 MeV State (Combined cuts)",90,0.5,180.5);
    TH1F *h1angdistL8_HeB = new TH1F("h1angdistL8_HeB","Angular Distribution Above 21.0 Below 21.7 MeV State (Combined cuts)",90,0.5,180.5);
    TH1F *h1angdistL9_HeB = new TH1F("h1angdistL9_HeB","Angular Distribution Above 21.7 Below 22.3 MeV State (Combined cuts)",90,0.5,180.5);
    TH1F *h1angdistL10_HeB = new TH1F("h1angdistL10_HeB","Angular Distribution Above 22.3 Below 23.3 MeV State (Combined cuts)",90,0.5,180.5);
    TH1F *h1angdistL11_HeB = new TH1F("h1angdistL11_HeB","Angular Distribution Above 23.3 MeV State (Combined cuts)",90,0.5,180.5);
    TH1F *h1angdistBG_He = new TH1F("h1angdistBG_He","Angular Distribution Background (Combined cuts)",90,0.5,180.5);
    TH1F *h1angdist5He = new TH1F("h1angdist5He", "Angular Distribution 5He Background", 90, 0.5, 180.5);

    //Peak angular distributions
    // Actually Combined
    TH1F *h1angdistL0_He = new TH1F("h1angdistL0_He","Angular Distribution 15.4 MeV State (C reconstruction and 5He cuts)",90,0.5,180.5);
    TH1F *h1angdistL1_He = new TH1F("h1angdistL1_He","Angular Distribution 16.5 MeV State (C reconstruction and 5He cuts)",90,0.5,180.5);
    TH1F *h1angdistL2_He = new TH1F("h1angdistL2_He","Angular Distribution 17.2 MeV State (C reconstruction and 5He cuts)",90,0.5,180.5);
    TH1F *h1angdistL3_He = new TH1F("h1angdistL3_He","Angular Distribution 17.5 MeV State (C reconstruction and 5He cuts)",90,0.5,180.5);
    TH1F *h1angdistL4_He = new TH1F("h1angdistL4_He","Angular Distribution 18.1 MeV State (C reconstruction and 5He cuts)",90,0.5,180.5);
    TH1F *h1angdistL5_He = new TH1F("h1angdistL5_He","Angular Distribution 18.8 MeV State (C reconstruction and 5He cuts)",90,0.5,180.5);
    TH1F *h1angdistL6_He = new TH1F("h1angdistL6_He","Angular Distribution 19.4 MeV State (C reconstruction and 5He cuts)",90,0.5,180.5);
    TH1F *h1angdistL7_He = new TH1F("h1angdistL7_He","Angular Distribution 21.0 MeV State (C reconstruction and 5He cuts)",90,0.5,180.5);
    TH1F *h1angdistL8_He = new TH1F("h1angdistL8_He","Angular Distribution 21.7 MeV State (C reconstruction and 5He cuts)",90,0.5,180.5);
    TH1F *h1angdistL9_He = new TH1F("h1angdistL9_He","Angular Distribution 22.3 MeV State (C reconstruction and 5He cuts)",90,0.5,180.5);
    TH1F *h1angdistL10_He = new TH1F("h1angdistL10_He","Angular Distribution 23.3 MeV State (C reconstruction and 5He cuts)",90,0.5,180.5);


    TH1F *h1angdistL0_C = new TH1F("h1angdistL0_C","Angular Distribution 15.4 MeV State (C reconstruction and 13C cuts)",90,0.5,180.5);
    TH1F *h1angdistL1_C = new TH1F("h1angdistL1_C","Angular Distribution 16.5 MeV State (C reconstruction and 13C cuts)",90,0.5,180.5);
    TH1F *h1angdistL2_C = new TH1F("h1angdistL2_C","Angular Distribution 17.2 MeV State (C reconstruction and 13C cuts)",90,0.5,180.5);
    TH1F *h1angdistL3_C = new TH1F("h1angdistL3_C","Angular Distribution 17.5 MeV State (C reconstruction and 13C cuts)",90,0.5,180.5);
    TH1F *h1angdistL4_C = new TH1F("h1angdistL4_C","Angular Distribution 18.1 MeV State (C reconstruction and 13C cuts)",90,0.5,180.5);
    TH1F *h1angdistL5_C = new TH1F("h1angdistL5_C","Angular Distribution 18.8 MeV State (C reconstruction and 13C cuts)",90,0.5,180.5);
    TH1F *h1angdistL6_C = new TH1F("h1angdistL6_C","Angular Distribution 19.4 MeV State (C reconstruction and 13C cuts)",90,0.5,180.5);
    TH1F *h1angdistL7_C = new TH1F("h1angdistL7_C","Angular Distribution 21.0 MeV State (C reconstruction and 13C cuts)",90,0.5,180.5);
    TH1F *h1angdistL8_C = new TH1F("h1angdistL8_C","Angular Distribution 21.7 MeV State (C reconstruction and 13C cuts)",90,0.5,180.5);
    TH1F *h1angdistL9_C = new TH1F("h1angdistL9_C","Angular Distribution 22.3 MeV State (C reconstruction and 13C cuts)",90,0.5,180.5);
    TH1F *h1angdistL10_C = new TH1F("h1angdistL10_C","Angular Distribution 23.3 MeV State (C reconstruction and 13C cuts)",90,0.5,180.5);
    TH1F *h1angdistBG_C = new TH1F("h1angdistBG_C","Angular Distribution Background (C reconstruction and 13C cuts)",90,0.5,180.5);
    TH1F *h1angdist13C = new TH1F("h1angdist13C", "Angular Distribution 13C Background", 90, 0.5, 180.5);


    TH1F *h1EMangdistL0_He = new TH1F("h1EMangdistL0_He","Angular Distribution 15.4 MeV State (Event Mixing C reconstruction and 5He cuts)",90,0.5,180.5);
    TH1F *h1EMangdistL1_He = new TH1F("h1EMangdistL1_He","Angular Distribution 16.5 MeV State (Event Mixing C reconstruction and 5He cuts)",90,0.5,180.5);
    TH1F *h1EMangdistL2_He = new TH1F("h1EMangdistL2_He","Angular Distribution 17.2 MeV State (Event Mixing C reconstruction and 5He cuts)",90,0.5,180.5);
    TH1F *h1EMangdistL3_He = new TH1F("h1EMangdistL3_He","Angular Distribution 17.5 MeV State (Event Mixing C reconstruction and 5He cuts)",90,0.5,180.5);
    TH1F *h1EMangdistL4_He = new TH1F("h1EMangdistL4_He","Angular Distribution 18.1 MeV State (Event Mixing C reconstruction and 5He cuts)",90,0.5,180.5);
    TH1F *h1EMangdistL5_He = new TH1F("h1EMangdistL5_He","Angular Distribution 18.8 MeV State (Event Mixing C reconstruction and 5He cuts)",90,0.5,180.5);
    TH1F *h1EMangdistL6_He = new TH1F("h1EMangdistL6_He","Angular Distribution 19.4 MeV State (Event Mixing C reconstruction and 5He cuts)",90,0.5,180.5);
    TH1F *h1EMangdistL7_He = new TH1F("h1EMangdistL7_He","Angular Distribution 21.0 MeV State (Event Mixing C reconstruction and 5He cuts)",90,0.5,180.5);
    TH1F *h1EMangdistL8_He = new TH1F("h1EMangdistL8_He","Angular Distribution 21.7 MeV State (Event Mixing C reconstruction and 5He cuts)",90,0.5,180.5);
    TH1F *h1EMangdistL9_He = new TH1F("h1EMangdistL9_He","Angular Distribution 22.3 MeV State (Event Mixing C reconstruction and 5He cuts)",90,0.5,180.5);
    TH1F *h1EMangdistL10_He = new TH1F("h1EMangdistL10_He","Angular Distribution 23.3 MeV State (Event Mixing C reconstruction and 5He cuts)",90, 0.5, 180.5);
    
    TH2F* angdistL0_Be = new TH2F("angdistL0_Be","angle vs Energy_Be0 15.4 MeV State",400,0.,20.,180,0.,180.);
    TH2F* angdistL1_Be = new TH2F("angdistL1_Be","angle vs Energy_Be1 16.5 MeV State",400,0.,20.,180,0.,180.);
    TH2F* angdistL2_Be = new TH2F("angdistL2_Be","angle vs Energy_Be2 17.2 MeV State",400,0.,20.,180,0.,180.);
    TH2F* angdistL3_Be = new TH2F("angdistL3_Be","angle vs Energy_Be3 17.5 MeV State",400,0.,20.,180,0.,180.);
    TH2F* angdistL4_Be = new TH2F("angdistL4_Be","angle vs Energy_Be4 18.1 MeV State",400,0.,20.,180,0.,180.);
    TH2F* angdistL5_Be = new TH2F("angdistL5_Be","angle vs Energy_Be5 18.8 MeV State",400,0.,20.,180,0.,180.);
    TH2F* angdistL6_Be = new TH2F("angdistL6_Be","angle vs Energy_Be6 19.4 MeV State",400,0.,20.,180,0.,180.);
    TH2F* angdistL7_Be = new TH2F("angdistL7_Be","angle vs Energy_Be7 21.0 MeV State",400,0.,20.,180,0.,180.);
    TH2F* angdistL8_Be = new TH2F("angdistL8_Be","angle vs Energy_Be8 21.7 MeV State",400,0.,20.,180,0.,180.);
    TH2F* angdistL9_Be = new TH2F("angdistL9_Be","angle vs Energy_Be9 22.3 MeV State",400,0.,20.,180,0.,180.);
    TH2F* angdistL10_Be = new TH2F("angdistL10_Be","angle vs Energy_Be10 23.3 MeV State",400,0.,20.,180,0.,180.);
    TH2F* angdistBG_Be = new TH2F("angdistBG_Be","angle vs Energy_BeBG",400,0.,20.,180,0.,180.);

    TH2F* angdistL0_He = new TH2F("angdistL0_He","angle vs Energy_He0 15.4 MeV State",400,0.,20.,180,0.,180.);
    TH2F* angdistL1_He = new TH2F("angdistL1_He","angle vs Energy_He1 16.5 MeV State",400,0.,20.,180,0.,180.);
    TH2F* angdistL2_He = new TH2F("angdistL2_He","angle vs Energy_He2 17.2 MeV State",400,0.,20.,180,0.,180.);
    TH2F* angdistL3_He = new TH2F("angdistL3_He","angle vs Energy_He3 17.5 MeV State",400,0.,20.,180,0.,180.);
    TH2F* angdistL4_He = new TH2F("angdistL4_He","angle vs Energy_He4 18.1 MeV State",400,0.,20.,180,0.,180.);
    TH2F* angdistL5_He = new TH2F("angdistL5_He","angle vs Energy_He5 18.8 MeV State",400,0.,20.,180,0.,180.);
    TH2F* angdistL6_He = new TH2F("angdistL6_He","angle vs Energy_He6 19.4 MeV State",400,0.,20.,180,0.,180.);
    TH2F* angdistL7_He = new TH2F("angdistL7_He","angle vs Energy_He7 21.0 MeV State",400,0.,20.,180,0.,180.);
    TH2F* angdistL8_He = new TH2F("angdistL8_He","angle vs Energy_He8 21.7 MeV State",400,0.,20.,180,0.,180.);
    TH2F* angdistL9_He = new TH2F("angdistL9_He","angle vs Energy_He9 22.3 MeV State",400,0.,20.,180,0.,180.);
    TH2F* angdistL10_He = new TH2F("angdistL10_He","angle vs Energy_He10 23.3 MeV State",400,0.,20.,180,0.,180.);
    TH2F* angdistBG_He = new TH2F("angdistBG_He","angle vs Energy_HeBG",400,0.,20.,180,0.,180.);

    TH2F* angdistL0_C = new TH2F("angdistL0_C","angle vs Energy_C0 15.4 MeV State",400,0.,20.,180,0.,180.);
    TH2F* angdistL1_C = new TH2F("angdistL1_C","angle vs Energy_C1 16.5 MeV State",400,0.,20.,180,0.,180.);
    TH2F* angdistL2_C = new TH2F("angdistL2_C","angle vs Energy_C2 17.2 MeV State",400,0.,20.,180,0.,180.);
    TH2F* angdistL3_C = new TH2F("angdistL3_C","angle vs Energy_C3 17.5 MeV State",400,0.,20.,180,0.,180.);
    TH2F* angdistL4_C = new TH2F("angdistL4_C","angle vs Energy_C4 18.1 MeV State",400,0.,20.,180,0.,180.);
    TH2F* angdistL5_C = new TH2F("angdistL5_C","angle vs Energy_C5 18.8 MeV State",400,0.,20.,180,0.,180.);
    TH2F* angdistL6_C = new TH2F("angdistL6_C","angle vs Energy_C6 19.4 MeV State",400,0.,20.,180,0.,180.);
    TH2F* angdistL7_C = new TH2F("angdistL7_C","angle vs Energy_C7 21.0 MeV State",400,0.,20.,180,0.,180.);
    TH2F* angdistL8_C = new TH2F("angdistL8_C","angle vs Energy_C8 21.7 MeV State",400,0.,20.,180,0.,180.);
    TH2F* angdistL9_C = new TH2F("angdistL9_C","angle vs Energy_C9 22.3 MeV State",400,0.,20.,180,0.,180.);
    TH2F* angdistL10_C = new TH2F("angdistL10_C","angle vs Energy_C10 23.3 MeV State",400,0.,20.,180,0.,180.);
    TH2F* angdistBG_C = new TH2F("angdistBG_C","angle vs Energy_CBG",400,0.,20.,180,0.,180.);
    
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
    Double_t p8bex[2] = {0.0}, p8bey[2] = {0.0}, p8bez[2] = {0.0}, e8be[2] = {0.0}, q8be = 0.0;                               // 8Be
    Double_t p16ox = 0.0, p16oy = 0.0, p16oz = 0.0, e16o = 0.0, q16o = 0.0, ex16o = 0.0;                                      // 16O (from 8Be)
    Double_t pnx = 0.0, pny = 0.0, pnz = 0.0, eneu = 0.0, etot = 0.0, qval = 0.0, pbeam = 0.0, ebeam = 0.0,  vbeam = 0.0;     // n, beam, & reaction info
    Double_t p12cx = 0.0, p12cy = 0.0, p12cz = 0.0, e12c = 0.0, q12c = 0.0, ex12c = 0.0;                                      // 12C
    Double_t p12cx1 = 0.0, p12cy1 = 0.0, p12cz1 = 0.0, e12c1 = 0.0, q12c1 = 0.0, ex12c1 = 0.0;                                // 12C for event mixing
    Double_t p16ox2 = 0.0, p16oy2 = 0.0, p16oz2 = 0.0, e16o2 = 0.0, q16o2 = 0.0, ex16o2 = 0.0;                                // 16O (from 12C)
    Double_t p16ox4 = 0.0, p16oy4 = 0.0, p16oz4 = 0.0, e16o4 = 0.0, q16o4 = 0.0, ex16o4 = 0.0;                                // 16O (from 12C 4 hit)
    Double_t pn2x = 0.0, pn2y = 0.0, pn2z = 0.0, eneu2 = 0.0, etot2 = 0.0, qval2 = 0.0;                                       // General Neutron
    Double_t pn2x4 = 0.0, pn2y4 = 0.0, pn2z4 = 0.0, eneu4 = 0.0, etot4 = 0.0, qval4 = 0.0;                                    // General Neutron 4 hit
    Double_t p9bex[2] = {0.0}, p9bey[2] = {0.0}, p9bez[2] = {0.0}, e9be[2] = {0.0}, q9be[2] = {0.0}, ex9be[2] = {0.0};        // 9Be
    Double_t p5hex = 0.0, p5hey = 0.0, p5hez = 0.0, e5he = 0.0, q5he = 0.0, ex5he = 0.0;                                      // 5He
    Double_t p13cx = 0.0, p13cy = 0.0, p13cz = 0.0, e13c = 0.0, q13c = 0.0, ex13c = 0.0;                                      // 13C
    
    
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
    
    
    // Read slope for crosstalk calculation
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
        lgood[i-1][j] = 0;
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
    ofstream ftext;
    ofstream fhits;
    char stext[255];
    char htext[255];
    sprintf(stext,"/afs/crc.nd.edu/group/nsl/nuc/users/bfrentz/analysis/output/data/check_run%d_event.txt", eGroup);
    sprintf(htext,"/afs/crc.nd.edu/group/nsl/nuc/users/bfrentz/analysis/output/data/check_run%d_hits.txt", eGroup);
    ftext.open(stext);
    fhits.open(htext);


    // ************
    // EVENT MIXING
    // The idea is to take hits/particles from different events and mix them toether in order to create uncorrelated data
    // For the different reconstruction methods i mix them differently.
    // For 8Be + 8Be, I take the 8Be's and use the ones from different events together.
    // For 12C + a, I combine the 12C and a's from different events together
    // for 12C, I take alpha's from different events, requiring the total reaction q-value to be that of 13C(a, n)4a
    
    // Create the arrays for event mixing
    Int_t em_countBe = 0;
    Int_t em_countO = 0;
    Double_t eMixingBe_px[6][2] = {0.0};
    Double_t eMixingBe_py[6][2] = {0.0};
    Double_t eMixingBe_pz[6][2] = {0.0};
    Double_t eMixingBe_e[6][2] = {0.0};
    Double_t eMixingO_px[480] = {0.0};
    Double_t eMixingO_py[480] = {0.0};
    Double_t eMixingO_pz[480] = {0.0};
    Double_t eMixingO_e[480] = {0.0};
    Double_t eMixingO_e1[480] = {0.0};
    Double_t eMixingO_q[480] = {0.0};
    Double_t eMixingO_ex[480] = {0.0};
    Double_t eMixingN_e[480] = {0.0};
    Double_t eMixingQ[480] = {0.0};

    Int_t em_countCa = 0;
    Int_t em_count1 = 0;
    Double_t eMixingCa_Cx[10] = {0.0};
    Double_t eMixingCa_Cy[10] = {0.0};
    Double_t eMixingCa_Cz[10] = {0.0};
    Double_t eMixingCa_Ce[10] = {0.0};
    Double_t eMixingCa_ax[10][7] = {0.0};
    Double_t eMixingCa_ay[10][7] = {0.0};
    Double_t eMixingCa_az[10][7] = {0.0};
    Double_t eMixingCa_ae[10][7] = {0.0};
    Double_t eMixingO_px1 = 0.0;
    Double_t eMixingO_py1 = 0.0;
    Double_t eMixingO_pz1 = 0.0;
    Double_t eMixingO_e_1 = 0.0;
    Double_t eMixingO_e11 = 0.0;
    Double_t eMixingO_q1 = 0.0;
    Double_t eMixingO_ex1 = 0.0;
    Double_t eMixingN_e1 = 0.0;
    Double_t eMixingQ1 = 0.0;

    Int_t em_countC = 0;
    Int_t em_count2 = 0;
    Double_t eMixing3a_px[3][16] = {0.0};
    Double_t eMixing3a_py[3][16] = {0.0};
    Double_t eMixing3a_pz[3][16] = {0.0};
    Double_t eMixing3a_e[3][16] = {0.0};
    Double_t eMixingC_px = 0.0;
    Double_t eMixingC_py = 0.0;
    Double_t eMixingC_pz = 0.0;
    Double_t eMixingC_e = 0.0;
    Double_t eMixingC_e1 = 0.0;
    Double_t eMixingC_q = 0.0;
    Double_t eMixingC_ex = 0.0;

    Int_t em_countC4 = 0;
    Int_t em_count4 = 0;
    Double_t eMixing4a_px[4][10] = {0.0};
    Double_t eMixing4a_py[4][10] = {0.0};
    Double_t eMixing4a_pz[4][10] = {0.0};
    Double_t eMixing4a_e[4][10] = {0.0};
    Double_t eMixing4a_cx = 0.0;
    Double_t eMixing4a_cy = 0.0;
    Double_t eMixing4a_cz = 0.0;
    Double_t eMixing4a_ce = 0.0;
    Double_t eMixing4a_q = 0.0;
    Double_t eMixing4a_cex = 0.0;
    Double_t eMixing4a_ox = 0.0;
    Double_t eMixing4a_oy = 0.0;
    Double_t eMixing4a_oz = 0.0;
    Double_t eMixing4a_oe = 0.0;
    Double_t eMixing4a_oe1 = 0.0;
    Double_t eMixing4a_oq = 0.0;
    Double_t eMixing4a_oex = 0.0;

    Int_t igoodMax = 0;
    Int_t igoodev[3] = {0};
    Int_t igoodev4[4] = {0};
    Int_t igoodevCa[10] = {0};
    
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
        
        

        // NEW, CORRECT MATCHING
        // pair up EF-EB hits in each of four silicons
        // use the new sophiscated matching algorithm match() in functions1.5.h
        Int_t igood=0;  // number of good EB-EF matched hits
        Int_t igoodl=0; // number of good EB-EF matched hits for DSSD 1, 2 (left)
        Int_t igoodr=0; // number of good EB-EF matched hits for DSSD 3, 4 (right)
        Int_t icheck=0; // counter for checking
        for (Int_t j=0; j < 4; j++) {   // loop over detectors
            int f = nhitef_new[j];
            int b = nhiteb_new[j];
            if (f>0 && b>0) {
                signal sigF[f];
                signal sigB[b];
                for(int l=0;l<f;l++){
                    sigF[l].SetValues(ef_new[j][l],stef_new[j][l]);
                }
                for(int l=0;l<b;l++){
                    sigB[l].SetValues(eb_new[j][l],steb_new[j][l]);
                }
                Int_t igood0 = match(j, nhitef_new[j], nhiteb_new[j], sigF, sigB, &hit[igood], match_cond);
                if (igood0 > 4)
                {
                    igood0 = 4;
                }
                igood += igood0;

                // if igood > 4 in any detector, set igood =
            }
        }

        //Find max igood
        if(igood > igoodMax) igoodMax = igood;


        for(int k = 0; k < igood; k++){
            if(hit[k].GetDetectorN() == 3 && -275 < hit[k].GetESpread() && hit[k].GetESpread() < -125 && hit[k].GetEnergy() > 2000 && hit[k].GetEnergy() < 12000){
                if(i < 1000000){
                        for (Int_t j=0; j < 4; j++) {   // loop over detectors
                            int f = nhitef_new[j];
                            int b = nhiteb_new[j];
                            if (f>0 && b>0) {
                                    signal sigF[f];
                                    signal sigB[b];
                                    for(int l=0;l<f;l++){
                                        sigF[l].SetValues(ef_new[j][l],stef_new[j][l]);
                                        if(j == 3){
                                                fhits << i << "\t" << "F" << "\t" << stef_new[j][l] << "\t" << ef_new[j][l] << endl;
                                        }
                                    }
                                    for(int l=0;l<b;l++){
                                        sigB[l].SetValues(eb_new[j][l],steb_new[j][l]);
                                        if(j == 3){
                                                fhits << i << "\t" << "B" << "\t" << steb_new[j][l] << "\t" << eb_new[j][l] << endl;
                                        }
                                    }
                            }
                        }
                }
            }
        }
        

        //Filling histograms
        double AvgChi2 = 0;
        int DetN = 0;
        int hitnumber=0;
        for(int j=0;j<igood;j++) {

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

            // ESpread for the error codes
            if (hit[k].GetDetectorN() == 0)
            {   
                if (hit[k].GetErrCode() == 1)
                {
                    hSpreaderr01->Fill(hit[k].GetESpread(), hit[k].GetEnergy());
                }

                if (hit[k].GetErrCode() == 2)
                {
                    hSpreaderr02->Fill(hit[k].GetESpread(), hit[k].GetEnergy());
                }

                if (hit[k].GetErrCode() == 3)
                {
                    hSpreaderr03->Fill(hit[k].GetESpread(), hit[k].GetEnergy());
                }

                if (hit[k].GetErrCode() == 4)
                {
                    hSpreaderr04->Fill(hit[k].GetESpread(), hit[k].GetEnergy());
                }
            }

            if (hit[k].GetDetectorN() == 1)
            {   
                if (hit[k].GetErrCode() == 1)
                {
                    hSpreaderr11->Fill(hit[k].GetESpread(), hit[k].GetEnergy());
                }

                if (hit[k].GetErrCode() == 2)
                {
                    hSpreaderr12->Fill(hit[k].GetESpread(), hit[k].GetEnergy());
                }

                if (hit[k].GetErrCode() == 3)
                {
                    hSpreaderr13->Fill(hit[k].GetESpread(), hit[k].GetEnergy());
                }

                if (hit[k].GetErrCode() == 4)
                {
                    hSpreaderr14->Fill(hit[k].GetESpread(), hit[k].GetEnergy());
                }
            }

            if (hit[k].GetDetectorN() == 2)
            {   
                if (hit[k].GetErrCode() == 1)
                {
                    hSpreaderr21->Fill(hit[k].GetESpread(), hit[k].GetEnergy());
                }

                if (hit[k].GetErrCode() == 2)
                {
                    hSpreaderr22->Fill(hit[k].GetESpread(), hit[k].GetEnergy());
                }

                if (hit[k].GetErrCode() == 3)
                {
                    hSpreaderr23->Fill(hit[k].GetESpread(), hit[k].GetEnergy());
                }

                if (hit[k].GetErrCode() == 4)
                {
                    hSpreaderr24->Fill(hit[k].GetESpread(), hit[k].GetEnergy());
                }
            }

            if (hit[k].GetDetectorN() == 3)
            {   
                if (hit[k].GetErrCode() == 1)
                {
                    hSpreaderr31->Fill(hit[k].GetESpread(), hit[k].GetEnergy());
                }

                if (hit[k].GetErrCode() == 2)
                {
                    hSpreaderr32->Fill(hit[k].GetESpread(), hit[k].GetEnergy());
                }

                if (hit[k].GetErrCode() == 3)
                {
                    hSpreaderr33->Fill(hit[k].GetESpread(), hit[k].GetEnergy());
                }

                if (hit[k].GetErrCode() == 4)
                {
                    hSpreaderr34->Fill(hit[k].GetESpread(), hit[k].GetEnergy());
                }
            }

            // Checking for the bad strip(s)
            if (hit[k].GetDetectorN() == 2 && -350 < hit[k].GetESpread() && hit[k].GetESpread() < -250 && 5000 < hit[k].GetEnergy() && hit[k].GetEnergy() < 7500) {
                hStrips->Fill(hit[k].GetStripF());
                hStrips->Fill(hit[k].GetStripB() + 32);
            }
            
            if (hit[k].GetDetectorN() == 2 && -500 < hit[k].GetESpread() && hit[k].GetESpread() < -340 && 15750 < hit[k].GetEnergy() && hit[k].GetEnergy() < 18250) {
                hStrips2->Fill(hit[k].GetStripF());
                hStrips2->Fill(hit[k].GetStripB() + 32);
            }



            // Check the individual events for the matched hits
            if (hit[k].GetDetectorN() == 3 && -275 < hit[k].GetESpread() && hit[k].GetESpread() < -125 && hit[k].GetEnergy() > 2000 && hit[k].GetEnergy() < 12000){
                //ftext << i << "\t" << hit[k].GetDouble() << "\t" << hit[k].GetStripF() << "\t" << hit[k].GetStripB() << "\t" << hit[k].GetEnergy() << "\t" << hit[k].GetESpread() << endl;
            }

        }
        
        
        
        
        
        

        // Angles
        for(int k = 0; k < igood; k++){
            Int_t detecNum = hit[k].GetDetectorN();
            Int_t j = hit[k].GetStripF();
            Int_t l = hit[k].GetStripB();
            if(detecNum > 1){
                l = hit[k].GetStripF();         //these are switched!!
                j = hit[k].GetStripB();
            }
            //Double_t theta = hit[k].GetAngleTheta(); //(pi/180)*(the[detecNum][j][l]);
            Double_t phiVal = hit[k].GetAnglePhi(); //(pi/180)*(phi[detecNum][j][l]);
            //Double_t theoVal = CheckGeom(theta,phiVal,detecNum, j, l);
            Double_t e1 = hit[k].GetEnergy();
            
            //cout <<"D: "<<detecNum<<" Vert Strip: " <<j<<" Horiz Strip: "<<l<<" Theta: "<<theta<<" rad "<<"Phi: "<<phiVal<<" rad "<<end$
            
            /*
             if(detecNum == 2){
             cout <<"D: "<<detecNum<<" Vert Strip: " <<j<<" Horiz Strip: "<<l<<" Theta: "<<theta<<" rad "<<"Phi: "<<phiVal<<" ra$
             }
             */
            
            //Double_t valDiff = (theoVal - e1);
            
            //if(hit[k].GetEnergy()>23800 && hit[k].GetEnergy()<24100 && l==16 && j==15 && detecNum == 2){

            if(l==16 && j==10 && detecNum == 1){
                strip10->Fill(e1);
                Double_t theta = hit[k].GetAngleTheta();
                Double_t theoVal = CheckGeom(theta*pi/180);
                strip10->Fill(theoVal);
                Double_t valDiff = (theoVal - e1);
                //CheckG->Fill(valDiff);
                AngDep->Fill(theta,valDiff);
                AngDepEn->Fill(theta,e1);
            }
            
            if(l==16 && j==9 && detecNum == 1){
                strip9->Fill(e1);
                Double_t theta = hit[k].GetAngleTheta();
                Double_t theoVal = CheckGeom(theta*pi/180);
                strip9->Fill(theoVal);
                Double_t valDiff = (theoVal - e1);
                //CheckG->Fill(valDiff);
                AngDep->Fill(theta,valDiff);
                AngDepEn->Fill(theta,e1);
            }
            
            if(l==16 && j==8 && detecNum == 1){
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
            
            if(l==16 && j==7 && detecNum == 1){ //&& hit[k].GetEnergy()>22600 && hit[k].GetEnergy()<23100){
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
            
            if(l==16 && j==6 && detecNum == 1){
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
            
            if(l==16 && j==5 && detecNum == 1){
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
            
            if(l==16 && j==4 && detecNum == 1){
                strip4->Fill(e1);
                Double_t theta = hit[k].GetAngleTheta(); //(pi/180)*(the[detecNum][j][l]);
                Double_t theoVal = CheckGeom(theta*pi/180);
                //cout << "Vert Strip Num: " << j<<"Theta: "<<theta<<"Theo E: "<<theoVal<<endl;
                strip4->Fill(theoVal);
                Double_t valDiff = (theoVal - e1);
                CheckG->Fill(valDiff);
                AngDep->Fill(theta,valDiff);
                AngDepEn->Fill(theta,e1);
            }
            
            if(l==16 && j==3 && detecNum == 1){
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
            
            if(l==16 && j==2 && detecNum == 1){
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

            if(l==16 && j==1 && detecNum == 1){
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
            
            if(l==16 && j==0 && detecNum == 1){
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
             cout << "j: "<<j<<" l: "<<l<<" Detec: "<<detecNum<<endl;
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
        //Int_t a[2][2];
        if (igood>1) {
            Int_t l1 = -1;  // index of 2nd alpha for 1st 8Be
            for (Int_t k=0; k<igood-1; k++) {
                if (l1 != k && igood8Be < 2) { // no more than two 8Be
                    for (Int_t l=k+1; l<igood; l++) {
                        if (l1 != l) {
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
                                    //a[igood8Be][0] = k;
                                    //a[igood8Be][1] = l;
                                    h1q8be_restrict->Fill(q8be);
                                    
                                    // Filling the Event Mixing arrays with 8Be's
                                    eMixingBe_px[em_countBe][igood8Be] = p8bex[igood8Be];
                                    eMixingBe_py[em_countBe][igood8Be] = p8bey[igood8Be];
                                    eMixingBe_pz[em_countBe][igood8Be] = p8bez[igood8Be];
                                    eMixingBe_e[em_countBe][igood8Be] = e8be[igood8Be];
                                    igood8Be++; l1=l;
                                    break;
                                }
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
            for(int j=0;j<12;j++) {
                if(j==11) {
                    if(event->qdc[j+2] > (s[j]*event->qdc[j+18]+o[j])) {
                        if(event->qdc[j+2]<3800){
                            if(event->qdc[j+2]>0){
                                h1nShort8Be[j]->Fill(event->qdc[j+18]);
                                h1nLong8Be[j]->Fill(event->qdc[j+2]);
                            }
                        }
                    }
                }
                else{
                    if(event->qdc[j+1] > (s[j]*event->qdc[j+17]+o[j])) {
                        if(event->qdc[j+1]<3800){
                            if(event->qdc[j+1]>0){
                                h1nShort8Be[j]->Fill(event->qdc[j+17]);
                                h1nLong8Be[j]->Fill(event->qdc[j+1]);
                            }
                        }
                    }
                }
            }
        }
        
        //generate 2D neutron histograms
        for(int j=0;j<12;j++) {
            if(j==11)
                h2n[j]->Fill(event->qdc[j+18],event->qdc[j+2]);
            else
                h2n[j]->Fill(event->qdc[j+17],event->qdc[j+1]);
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
            
            //Peaks
            bool isL00_Be = (ex16o>LBe_00Min*1000 && ex16o<LBe_00Max*1000 && LBe_00Max != 0);
            bool isL0_Be = (ex16o>LBe_0Min*1000 && ex16o<LBe_0Max*1000 && LBe_0Max != 0);
            bool isL1_Be = (ex16o>LBe_1Min*1000 && ex16o<LBe_1Max*1000 && LBe_1Max != 0);
            bool isL2_Be = (ex16o>LBe_2Min*1000 && ex16o<LBe_2Max*1000 && LBe_2Max != 0);
            bool isL3_Be = (ex16o>LBe_3Min*1000 && ex16o<LBe_3Max*1000 && LBe_3Max != 0);
            bool isL4_Be = (ex16o>LBe_4Min*1000 && ex16o<LBe_4Max*1000 && LBe_4Max != 0);
            bool isL5_Be = (ex16o>LBe_5Min*1000 && ex16o<LBe_5Max*1000 && LBe_5Max != 0);
            bool isL6_Be = (ex16o>LBe_6Min*1000 && ex16o<LBe_6Max*1000 && LBe_6Max != 0);
            bool isL7_Be = (ex16o>LBe_7Min*1000 && ex16o<LBe_7Max*1000 && LBe_7Max != 0);
            bool isL8_Be = (ex16o>LBe_8Min*1000 && ex16o<LBe_8Max*1000 && LBe_8Max != 0);
            bool isL9_Be = (ex16o>LBe_9Min*1000 && ex16o<LBe_9Max*1000 && LBe_9Max != 0);
            bool isL10_Be = (ex16o>LBe_10Min*1000 && ex16o<LBe_10Max*1000 && LBe_10Max != 0);
            bool isBG_Be = !isL00_Be && !isL0_Be && !isL1_Be && !isL2_Be && !isL3_Be && !isL4_Be && !isL5_Be && !isL6_Be && !isL7_Be && !isL8_Be && !isL9_Be && !isL10_Be;

            //Backgrounds?
            bool isL0_BeB = (ex16o<LBe_0Min*1000 && LBe_0Max != 0);
            bool isL1_BeB = (ex16o<LBe_1Min*1000 && ex16o>LBe_0Max*1000 && LBe_1Min != 0 && LBe_0Max != LBe_1Min);
            bool isL2_BeB = (ex16o<LBe_2Min*1000 && ex16o>LBe_1Max*1000 && LBe_2Min != 0 && LBe_1Max != LBe_2Min);
            bool isL3_BeB = (ex16o<LBe_3Min*1000 && ex16o>LBe_2Max*1000 && LBe_3Min != 0 && LBe_2Max != LBe_3Min);
            bool isL4_BeB = (ex16o<LBe_4Min*1000 && ex16o>LBe_3Max*1000 && LBe_4Min != 0 && LBe_3Max != LBe_4Min);
            bool isL5_BeB = (ex16o<LBe_5Min*1000 && ex16o>LBe_4Max*1000 && LBe_5Min != 0 && LBe_4Max != LBe_5Min);
            bool isL6_BeB = (ex16o<LBe_6Min*1000 && ex16o>LBe_5Max*1000 && LBe_6Min != 0 && LBe_5Max != LBe_6Min);
            bool isL7_BeB = (ex16o<LBe_7Min*1000 && ex16o>LBe_6Max*1000 && LBe_7Min != 0 && LBe_6Max != LBe_7Min);
            bool isL8_BeB = (ex16o<LBe_8Min*1000 && ex16o>LBe_7Max*1000 && LBe_8Min != 0 && LBe_7Max != LBe_8Min);
            bool isL9_BeB = (ex16o<LBe_9Min*1000 && ex16o>LBe_8Max*1000 && LBe_9Min != 0 && LBe_8Max != LBe_9Min);
            bool isL10_BeB = (ex16o<LBe_10Min*1000 && ex16o>LBe_9Max*1000 && LBe_10Min != 0 && LBe_9Max != LBe_10Min);
            bool isL11_BeB = (ex16o>LBe_10Max*1000 && LBe_10Max != 0);

           
            
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
                h1_8be_corr->Fill(q16o/1000.0);
                if(em_countBe < 6) em_countBe++;
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

                h28ben->Fill(ex9be[0]/1000.0, ex9be[1]/1000.0);
            }
            //Best Gate: 4000 && 3500
            if(qval > Q16OMin && qval < Q16OMax && ex9be[0]>Be9Gate1 && ex9be[1]>Be9Gate1){
                h1ex16ogated->Fill(ex16o/1000.0);

                if(ex9be[0] > Be9Gate_best && ex9be[1] > Be9Gate_best){
                    h1ex16o_bebest->Fill(ex16o/1000.0);

                    if(isL00_Be){
                        h1angdistL00_Be->Fill(angleCoM);
                    }
                    if(isL0_Be){
                        h1angdistL0_Be->Fill(angleCoM);
                        angdistL0_Be->Fill(e16o/1000.0, angleCoM);
                    }
                    if(isL1_Be){
                        h1angdistL1_Be->Fill(angleCoM);
                        angdistL1_Be->Fill(e16o/1000.0, angleCoM);
                    }
                    if(isL2_Be){
                        h1angdistL2_Be->Fill(angleCoM);
                        angdistL2_Be->Fill(e16o/1000.0, angleCoM);
                    }
                    if(isL3_Be){
                        h1angdistL3_Be->Fill(angleCoM);
                        angdistL3_Be->Fill(e16o/1000.0, angleCoM);
                    }
                    if(isL4_Be){
                        h1angdistL4_Be->Fill(angleCoM);
                        angdistL4_Be->Fill(e16o/1000.0, angleCoM);
                    }
                    if(isL5_Be){
                        h1angdistL5_Be->Fill(angleCoM);
                        angdistL5_Be->Fill(e16o/1000.0, angleCoM);
                    }
                    if(isL6_Be){
                        h1angdistL6_Be->Fill(angleCoM);
                        angdistL6_Be->Fill(e16o/1000.0, angleCoM);
                    }
                    if(isL7_Be){
                        h1angdistL7_Be->Fill(angleCoM);
                        angdistL7_Be->Fill(e16o/1000.0, angleCoM);
                    }
                    if(isL8_Be){
                        h1angdistL8_Be->Fill(angleCoM);
                        angdistL8_Be->Fill(e16o/1000.0, angleCoM);
                    }
                    if(isL9_Be){
                        h1angdistL9_Be->Fill(angleCoM);
                        angdistL9_Be->Fill(e16o/1000.0, angleCoM);
                    }
                    if(isL10_Be){
                        h1angdistL10_Be->Fill(angleCoM);
                        angdistL10_Be->Fill(e16o/1000.0, angleCoM);
                    }
        
                    //Backgrounds
                    if(isL0_BeB){
                        h1angdistL0_BeB->Fill(angleCoM);
                    }
                    if(isL1_BeB){
                        h1angdistL1_BeB->Fill(angleCoM);
                    }
                    if(isL2_BeB){
                        h1angdistL2_BeB->Fill(angleCoM);
                    }
                    if(isL3_BeB){
                        h1angdistL3_BeB->Fill(angleCoM);
                    }
                    if(isL4_BeB){
                        h1angdistL4_BeB->Fill(angleCoM);
                    }
                    if(isL5_BeB){
                        h1angdistL5_BeB->Fill(angleCoM);
                    }
                    if(isL6_BeB){
                        h1angdistL6_BeB->Fill(angleCoM);
                    }
                    if(isL7_BeB){
                        h1angdistL7_BeB->Fill(angleCoM);
                    }
                    if(isL8_BeB){
                        h1angdistL8_BeB->Fill(angleCoM);
                    }
                    if(isL9_BeB){
                        h1angdistL9_BeB->Fill(angleCoM);
                    }
                    if(isL10_BeB){
                        h1angdistL10_BeB->Fill(angleCoM);
                    }
                    if(isL11_BeB){
                        h1angdistL11_BeB->Fill(angleCoM);
                    }
                }
                else h1angdist9Be->Fill(angleCoM);

                if(ex9be[0] > Be9Gate2 && ex9be[1] > Be9Gate2){
                    h1ex16ogated2->Fill(ex16o/1000.0);
                }

                if(ex9be[0] > Be9Gate3 && ex9be[1] > Be9Gate3){
                    h1ex16ogated3->Fill(ex16o/1000.0);
                }

                if(ex9be[0] > Be9Gate4 && ex9be[1] > Be9Gate4){
                    h1ex16ogated4->Fill(ex16o/1000.0);
                }

                if(ex9be[0] > Be9Gateb11 && ex9be[1] > Be9Gateb11){
                    h1ex16ogated_b11->Fill(ex16o/1000.0);
                }

                if(ex9be[0] > Be9Gateb12 && ex9be[1] > Be9Gateb12){
                    h1ex16ogated_b12->Fill(ex16o/1000.0);
                }

                if(ex9be[0] > Be9Gateb13 && ex9be[1] > Be9Gateb13){
                    h1ex16ogated_b13->Fill(ex16o/1000.0);
                }

                if(ex9be[0] > Be9Gateb14 && ex9be[1] > Be9Gateb14){
                    h1ex16ogated_b14->Fill(ex16o/1000.0);
                }

                if(ex9be[0] > Be9Gateb15 && ex9be[1] > Be9Gateb15){
                    h1ex16ogated_b15->Fill(ex16o/1000.0);
                }

                if(ex9be[0] > Be9Gateb16 && ex9be[1] > Be9Gateb16){
                    h1ex16ogated_b16->Fill(ex16o/1000.0);
                }

                if(ex9be[0] > Be9Gateb17 && ex9be[1] > Be9Gateb17){
                    h1ex16ogated_b17->Fill(ex16o/1000.0);
                }

                if(ex9be[0] > Be9Gateb18 && ex9be[1] > Be9Gateb18){
                    h1ex16ogated_b18->Fill(ex16o/1000.0);
                }

                if(ex9be[0] > Be9Gateb19 && ex9be[1] > Be9Gateb19){
                    h1ex16ogated_b19->Fill(ex16o/1000.0);
                }

                if(ex9be[0] > Be9Gateb20 && ex9be[1] > Be9Gateb20){
                    h1ex16ogated_b20->Fill(ex16o/1000.0);
                }


                
                //Neutron detectors gated on 16O
                for(int j=0;j<12;j++) {
                    if(j==11) {
                        if(event->qdc[j+2] > (s[j]*event->qdc[j+18]+o[j])) {
                            if(event->qdc[j+2]<3800){
                                if(event->qdc[j+2]>0){
                                    h1nShort16O[j]->Fill(event->qdc[j+18]);
                                    h1nLong16O[j]->Fill(event->qdc[j+2]);
                                }
                            }
                        }
                    }
                    else{
                        if(event->qdc[j+1] > (s[j]*event->qdc[j+17]+o[j])) {
                            if(event->qdc[j+1]<3800){
                                if(event->qdc[j+1]>0){
                                    h1nShort16O[j]->Fill(event->qdc[j+17]);
                                    h1nLong16O[j]->Fill(event->qdc[j+1]);
                                }
                            }
                        }
                    }
                }
            }
        }

        //Event mixing 8Be
        //Create 16O
        if(em_countBe == 6){
        	em_countO = 0;
        	em_countBe = 0;
        	for(int j = 0; j < 5; j++){
        		for(int k = 0; k < 2; k++){
        			for(int l = j+1; l < 6; l++){
        				for(int m = 0; m < 2; m++){
        					eMixingO_px[em_countO] = eMixingBe_px[j][k] + eMixingBe_px[l][m];
        					eMixingO_py[em_countO] = eMixingBe_py[j][k] + eMixingBe_py[l][m];
        					eMixingO_pz[em_countO] = eMixingBe_pz[j][k] + eMixingBe_pz[l][m];
        					eMixingO_e[em_countO] = (eMixingO_px[em_countO]*eMixingO_px[em_countO] + eMixingO_py[em_countO]*eMixingO_py[em_countO] + eMixingO_pz[em_countO]*eMixingO_pz[em_countO])/2.0/15.9949146;
        					eMixingO_q[em_countO] = eMixingBe_e[j][k] + eMixingBe_e[l][m] - eMixingO_e[em_countO];
        					eMixingO_ex[em_countO] = eMixingO_q[em_countO] + 14620.3;
        					eMixingO_e1[em_countO] = eMixingBe_e[j][k] + eMixingBe_e[l][m];
        					em_countO++;
        				} 
        			}
        		}
        	}

        	// EM Angles
        	for(int j = 0; j < 480; j++){
        		//Angle of 16O in lab frame with radians
            	Double_t tgThetaLab_em = sqrt(pow(eMixingO_px[j],2) + pow(eMixingO_py[j],2))/eMixingO_pz[j];
            	Double_t cosThetaLab_em;
            	if(tgThetaLab_em >= 0)
            	    cosThetaLab_em = 1/sqrt(1+pow(tgThetaLab_em,2));
            	else
            	    cosThetaLab_em = -1/sqrt(1+pow(tgThetaLab_em,2));
            	
            	//Velocity of 16O in lab frame
            	Double_t v16OLab_em = sqrt(pow(eMixingO_px[j],2)+pow(eMixingO_py[j],2)+pow(eMixingO_pz[j],2))/15.9949146;
            	
            	//Velocity of center of momentum frame in final stage
            	Double_t VCoM_em = (4.00260325415/(1.00866491574 + 15.99491461956))*vbeam;
            	
            	//Velocity of 16O in center of mass frame
            	Double_t V16OCM_em = sqrt(pow(v16OLab_em,2)+pow(VCoM_em,2)-2*v16OLab_em*VCoM_em*cosThetaLab_em);
            	
            	//Angle of 16O in CoM
            	Double_t angleCoM_em = acos((v16OLab_em*cosThetaLab_em - VCoM_em)/V16OCM_em)*180/3.14159265359;

            	//Fill histograms
            	//h1angdistBG_Be->Fill(angleCoM_em);
            	//angdistBG_Be->Fill(eMixingO_ex[j]/1000.0, angleCoM_em);

            	// reconstruct neutron
            	eMixingN_e[j] = ((eMixingO_px[j]*eMixingO_px[j]) + (eMixingO_py[j]*eMixingO_py[j]) + ((pbeam - eMixingO_pz[j])*(pbeam - eMixingO_pz[j])))/2.0/1.0086649;
            	eMixingQ[j] = eMixingO_e1[j] + eMixingN_e[j] - ebeam;

            	h1ex_em->Fill(eMixingO_ex[j]/1000.0);
            	h1q_em->Fill(eMixingQ[j]/1000.0);
            	h1eneu_em->Fill(eMixingN_e[j]/1000.0);
            	h1e1_em->Fill(eMixingO_e1[j]/1000.0);

            	if(eMixingQ[j] > Q16OMin && eMixingQ[j] < Q16OMax){
            		h1ex16o_em->Fill(eMixingO_q[j]/1000.0);
            		h1angdistBG_Be->Fill(angleCoM_em);
            		angdistBG_Be->Fill(eMixingO_q[j]/1000.0, angleCoM_em);

                    bool isL0_Be = (eMixingO_ex[j]>LBe_0Min*1000 && eMixingO_ex[j]<LBe_0Max*1000 && LBe_0Max != 0);
                    bool isL1_Be = (eMixingO_ex[j]>LBe_1Min*1000 && eMixingO_ex[j]<LBe_1Max*1000 && LBe_1Max != 0);
                    bool isL2_Be = (eMixingO_ex[j]>LBe_2Min*1000 && eMixingO_ex[j]<LBe_2Max*1000 && LBe_2Max != 0);
                    bool isL3_Be = (eMixingO_ex[j]>LBe_3Min*1000 && eMixingO_ex[j]<LBe_3Max*1000 && LBe_3Max != 0);
                    bool isL4_Be = (eMixingO_ex[j]>LBe_4Min*1000 && eMixingO_ex[j]<LBe_4Max*1000 && LBe_4Max != 0);
                    bool isL5_Be = (eMixingO_ex[j]>LBe_5Min*1000 && eMixingO_ex[j]<LBe_5Max*1000 && LBe_5Max != 0);
                    bool isL6_Be = (eMixingO_ex[j]>LBe_6Min*1000 && eMixingO_ex[j]<LBe_6Max*1000 && LBe_6Max != 0);
                    bool isL7_Be = (eMixingO_ex[j]>LBe_7Min*1000 && eMixingO_ex[j]<LBe_7Max*1000 && LBe_7Max != 0);
                    bool isL8_Be = (eMixingO_ex[j]>LBe_8Min*1000 && eMixingO_ex[j]<LBe_8Max*1000 && LBe_8Max != 0);
                    bool isL9_Be = (eMixingO_ex[j]>LBe_9Min*1000 && eMixingO_ex[j]<LBe_9Max*1000 && LBe_9Max != 0);
                    bool isL10_Be = (eMixingO_ex[j]>LBe_10Min*1000 && eMixingO_ex[j]<LBe_10Max*1000 && LBe_10Max != 0);
                    if(isL0_Be){
                        h1EMangdistL0_Be->Fill(angleCoM_em);
                    }
                    if(isL1_Be){
                        h1EMangdistL1_Be->Fill(angleCoM_em);
                    }
                    if(isL2_Be){
                        h1EMangdistL2_Be->Fill(angleCoM_em);
                    }
                    if(isL3_Be){
                        h1EMangdistL3_Be->Fill(angleCoM_em);
                    }
                    if(isL4_Be){
                        h1EMangdistL4_Be->Fill(angleCoM_em);
                    }
                    if(isL5_Be){
                        h1EMangdistL5_Be->Fill(angleCoM_em);
                    }
                    if(isL6_Be){
                        h1EMangdistL6_Be->Fill(angleCoM_em);
                    }
                    if(isL7_Be){
                        h1EMangdistL7_Be->Fill(angleCoM_em);
                    }
                    if(isL8_Be){
                        h1EMangdistL8_Be->Fill(angleCoM_em);
                    }
                    if(isL9_Be){
                        h1EMangdistL9_Be->Fill(angleCoM_em);
                    }
                    if(isL10_Be){
                        h1EMangdistL10_Be->Fill(angleCoM_em);
                    }

                    
            	}
        	}

        	//Clearing arrays
			for(int j = 0; j < 6; j++){
        		for(int k = 0; k < 2; k++){
        			eMixingBe_px[j][k] = 0.0;
        			eMixingBe_py[j][k] = 0.0;
        			eMixingBe_pz[j][k] = 0.0;
        			eMixingBe_e[j][k] = 0.0;
        		}
        	}

        	for(int j = 0; j < 480; j++){
        		eMixingO_px[j] = 0.0;
				eMixingO_py[j] = 0.0;
				eMixingO_pz[j] = 0.0;
				eMixingO_e[j] = 0.0;
                eMixingO_e1[j] = 0.0;
				eMixingO_q[j] = 0.0;
				eMixingO_ex[j] = 0.0;
                eMixingN_e[j] = 0.0;
                eMixingQ[j] = 0.0;
        	}
        }

        // END THE EVENT MIXING OF 8Be + 8Be
        // *********************************

        
        //Neutron detectors gated
        for(int j=0;j<12;j++) {
            if(j==11) {
                if(event->qdc[j+2] > (s[j]*event->qdc[j+18]+o[j])) {
                    if(event->qdc[j+2]<3800){
                        if(event->qdc[j+2]>0){
                            h1nShortTotal[j]->Fill(event->qdc[j+18]);
                            h1nLongTotal[j]->Fill(event->qdc[j+2]);
                        }
                    }
                }
            }
            else {
                if(event->qdc[j+1] > (s[j]*event->qdc[j+17]+o[j])) {
                    if(event->qdc[j+1]<3800){
                        if(event->qdc[j+1]>0){
                            h1nShortTotal[j]->Fill(event->qdc[j+17]);
                            h1nLongTotal[j]->Fill(event->qdc[j+1]);
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
        Int_t n12ca4 = 0;
        bool em12c = false;

        //reconstruct 12C for 4a+n channel
        Int_t igood12C4a = 0;
        if (igood>3) {
            for (Int_t m=0; m < igood-3; m++) { if(igood12C4a == 1) break;
                for (Int_t k=m+1; k < igood-2; k++) { if(igood12C4a == 1) break; if (m == k) continue;
                    for (Int_t l=k+1; l < igood-1; l++) { if(igood12C4a == 1) break; if (m == l) continue;
                        for (Int_t j=l+1; j < igood; j++) { if (m == j) continue;

                            // Sum of Momenta
                            Double_t p12cx2 = hit[k].GetMomentumX()+hit[l].GetMomentumX()+hit[j].GetMomentumX();
                            Double_t p12cy2 = hit[k].GetMomentumY()+hit[l].GetMomentumY()+hit[j].GetMomentumY();
                            Double_t p12cz2 = hit[k].GetMomentumZ()+hit[l].GetMomentumZ()+hit[j].GetMomentumZ();
                            // Energy
                            Double_t e12c2 = (p12cx2*p12cx2+p12cy2*p12cy2+p12cz2*p12cz2)/2.0/12.0;
                            // Q-value and Kinetic Energy (g.s. Q-value = -7274.7 keV)
                            Double_t q12c2 = hit[k].GetEnergy() + hit[l].GetEnergy() + hit[j].GetEnergy() - e12c2; 
                            Double_t ex12c2 = (hit[k].GetEnergy() + hit[l].GetEnergy() + hit[j].GetEnergy()) - e12c2 + 7274.7;

                            // reconstruct neutron
                            Double_t pnx2 = -p12cx2-hit[m].GetMomentumX();
                            Double_t pny2 = -p12cy2-hit[m].GetMomentumY();
                            Double_t pnz2 = pbeam-p12cz2-hit[m].GetMomentumZ();
                            Double_t eneu2 = (pnx2*pnx2+pny2*pny2+pnz2*pnz2)/2.0/1.0086649;

                            // reaction Q-value
                            Double_t etot2 = hit[k].GetEnergy() + hit[l].GetEnergy() + hit[j].GetEnergy() + hit[m].GetEnergy() + eneu2;    //           4a+n
                            Double_t qval2 = etot2 - ebeam;    // 13C+a->4a+n Qvalue = -12221.05522 keV
                            h1qval4a->Fill(qval2/1000.0);

                            if (qval2 < Q16Ohigh && qval2 > Q16Olow) {
                                h1q12c->Fill(q12c2/1000.0);
                                igood12C4a++;
                                h1_12c_corr->Fill(q12c2/1000.0);
                                h1ex12c2->Fill(ex12c2/1000.0);

                                

                                break;
                            }
                        }
                    }
                }
            }
        }

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
                                //h1q12c->Fill(q12c/1000.0);
                                ex12c = (hit[k].GetEnergy() + hit[l].GetEnergy() + hit[j].GetEnergy()) - e12c + 7274.7;     // Excitation Energy
                                
                                // Fill histogram with excitation energy
                                h1ex12c->Fill(ex12c/1000.0);
                                
                                // Fill 2D histogram
                                h2q12c->Fill(ex12c/1000.0, (hit[k].GetEnergy() + hit[l].GetEnergy() + hit[j].GetEnergy())/1000.0);

                                // Give a good 12C if the energy is a Hoyle state (Excitation Energy = 7.6542 MeV)
                                // Peak cut is 7.6 MeV and 7.68 MeV
                                if (ex12c > C12Gate1 && ex12c < C12Gate2) {
                                    n12ca1 = k;
                                    n12ca2 = l;
                                    n12ca3 = j;
                                    igood12C++;
                                }

                                if(igood12C == 1) break;
                            }
                        }
                    }
                }
            }
        }

        // Reconstruct 16O from 12C
        if (igood12C == 1 && igood > 3) {

            for (Int_t j = 0; j < igood; j++) {
                // Only use particles that aren't used in 12C reconstruction
                if(j != n12ca1 && j != n12ca2 && j != n12ca3){
                    // Sum of momenta
                    p16ox2 = p12cx + hit[j].GetMomentumX();
                    p16oy2 = p12cy + hit[j].GetMomentumY();
                    p16oz2 = p12cz + hit[j].GetMomentumZ();
                    //Energies
                    e16o2 = (p16ox2*p16ox2 + p16oy2*p16oy2 + p16oz2*p16oz2)/2.0/15.9949146;     // KE of 16O
                    q16o2 = e12c + hit[j].GetEnergy() - e16o2;                                  // Measured Q-value of reaction (Should be -7161.91698 keV)
                    ex16o2 = q16o2 + 7161.91698 + 7654.20;                                      // Excited state energies (ground + hoyle state)
                    
                    
                    // Reconstruct neutron
                    pn2x = -p16ox2;
                    pn2y = -p16oy2;
                    pn2z = pbeam - p16oz2;
                    eneu2 = (pn2x*pn2x+pn2y*pn2y+pn2z*pn2z)/2.0/1.0086649;
                    
                    // Reaction Q-Value
                    etot2 = e12c + hit[j].GetEnergy() + eneu2;
                    qval2 = etot2 - ebeam;                            // 13C(a, 12C_hs + a + n) Q-value = -12600.508
                    h1qval2->Fill(qval2/1000.0);
                    if (qval2 < Q16Ohigh && qval2 > Q16Olow) {
                        h1ex16o2->Fill(ex16o2/1000);
                        //if(em_countCa < 10) em_countCa++;
                    }

                    // Angular distributions for this reconstruction method
                    //Angle of 16O in lab frame with radians
                    Double_t tgThetaLab2 = sqrt(pow(p16ox2,2) + pow(p16oy2,2))/p16oz2;
                    Double_t cosThetaLab2;
                    if(tgThetaLab2>=0)
                        cosThetaLab2 = 1/sqrt(1+pow(tgThetaLab2,2));
                    else
                        cosThetaLab2 = -1/sqrt(1+pow(tgThetaLab2,2));
                    
                    //Velocity of 16O in lab frame
                    Double_t v16OLab2 = sqrt(pow(p16ox2,2)+pow(p16oy2,2)+pow(p16oz2,2))/15.9949146;
                    
                    //Velocity of center of momentum frame in final stage
                    Double_t VCoM2 = (4.00260325415/(1.00866491574 + 15.99491461956))*vbeam;
                    
                    //Velocity of 16O in center of mass frame
                    Double_t V16OCM2 = sqrt(pow(v16OLab2,2)+pow(VCoM2,2)-2*v16OLab2*VCoM2*cosThetaLab2);
                    
                    //Angle of 16O in CoM
                    Double_t angleCoM2 = acos((v16OLab2*cosThetaLab2 - VCoM2)/V16OCM2)*180/3.14159265359;
                    
                    // Peaks
                    // USING HE FOR COMBINED
                    bool isL0_He = (ex16o2>LHe_0Min*1000 && ex16o2<LHe_0Max*1000 && LHe_0Max != 0);
                    bool isL1_He = (ex16o2>LHe_1Min*1000 && ex16o2<LHe_1Max*1000 && LHe_1Max != 0);
                    bool isL2_He = (ex16o2>LHe_2Min*1000 && ex16o2<LHe_2Max*1000 && LHe_2Max != 0);
                    bool isL3_He = (ex16o2>LHe_3Min*1000 && ex16o2<LHe_3Max*1000 && LHe_3Max != 0);
                    bool isL4_He = (ex16o2>LHe_4Min*1000 && ex16o2<LHe_4Max*1000 && LHe_4Max != 0);
                    bool isL5_He = (ex16o2>LHe_5Min*1000 && ex16o2<LHe_5Max*1000 && LHe_5Max != 0);
                    bool isL6_He = (ex16o2>LHe_6Min*1000 && ex16o2<LHe_6Max*1000 && LHe_6Max != 0);
                    bool isL7_He = (ex16o2>LHe_7Min*1000 && ex16o2<LHe_7Max*1000 && LHe_7Max != 0);
                    bool isL8_He = (ex16o2>LHe_8Min*1000 && ex16o2<LHe_8Max*1000 && LHe_8Max != 0);
                    bool isL9_He = (ex16o2>LHe_9Min*1000 && ex16o2<LHe_9Max*1000 && LHe_9Max != 0);
                    bool isL10_He = (ex16o2>LHe_10Min*1000 && ex16o2<LHe_10Max*1000 && LHe_10Max != 0);
                    bool isBG_He = !isL0_He && !isL1_He && !isL2_He && !isL3_He && !isL4_He && !isL5_He && !isL6_He && !isL7_He && !isL8_He && !isL9_He && !isL10_He;

                    //Backgrounds?
                    bool isL0_HeB = (ex16o2<LHe_0Min*1000 && LHe_0Min != 0);
                    bool isL1_HeB = (ex16o2<LHe_1Min*1000 && ex16o2>LHe_0Max*1000 && LHe_1Min != 0 && LHe_0Max != LHe_1Min);
                    bool isL2_HeB = (ex16o2<LHe_2Min*1000 && ex16o2>LHe_1Max*1000 && LHe_2Min != 0 && LHe_1Max != LHe_2Min);
                    bool isL3_HeB = (ex16o2<LHe_3Min*1000 && ex16o2>LHe_2Max*1000 && LHe_3Min != 0 && LHe_2Max != LHe_3Min);
                    bool isL4_HeB = (ex16o2<LHe_4Min*1000 && ex16o2>LHe_3Max*1000 && LHe_4Min != 0 && LHe_3Max != LHe_4Min);
                    bool isL5_HeB = (ex16o2<LHe_5Min*1000 && ex16o2>LHe_4Max*1000 && LHe_5Min != 0 && LHe_4Max != LHe_5Min);
                    bool isL6_HeB = (ex16o2<LHe_6Min*1000 && ex16o2>LHe_5Max*1000 && LHe_6Min != 0 && LHe_5Max != LHe_6Min);
                    bool isL7_HeB = (ex16o2<LHe_7Min*1000 && ex16o2>LHe_6Max*1000 && LHe_7Min != 0 && LHe_6Max != LHe_7Min);
                    bool isL8_HeB = (ex16o2<LHe_8Min*1000 && ex16o2>LHe_7Max*1000 && LHe_8Min != 0 && LHe_7Max != LHe_8Min);
                    bool isL9_HeB = (ex16o2<LHe_9Min*1000 && ex16o2>LHe_8Max*1000 && LHe_9Min != 0 && LHe_8Max != LHe_9Min);
                    bool isL10_HeB = (ex16o2<LHe_10Min*1000 && ex16o2>LHe_9Max*1000 && LHe_10Min != 0 && LHe_9Max != LHe_10Min);
                    bool isL11_HeB = (ex16o2>LHe_10Max*1000 && LHe_10Max != 0);

                    

                    bool isL0_C = (ex16o2>LC_0Min*1000 && ex16o2<LC_0Max*1000 && LC_0Max != 0);
                    bool isL1_C = (ex16o2>LC_1Min*1000 && ex16o2<LC_1Max*1000 && LC_1Max != 0);
                    bool isL2_C = (ex16o2>LC_2Min*1000 && ex16o2<LC_2Max*1000 && LC_2Max != 0);
                    bool isL3_C = (ex16o2>LC_3Min*1000 && ex16o2<LC_3Max*1000 && LC_3Max != 0);
                    bool isL4_C = (ex16o2>LC_4Min*1000 && ex16o2<LC_4Max*1000 && LC_4Max != 0);
                    bool isL5_C = (ex16o2>LC_5Min*1000 && ex16o2<LC_5Max*1000 && LC_5Max != 0);
                    bool isL6_C = (ex16o2>LC_6Min*1000 && ex16o2<LC_6Max*1000 && LC_6Max != 0);
                    bool isL7_C = (ex16o2>LC_7Min*1000 && ex16o2<LC_7Max*1000 && LC_7Max != 0);
                    bool isL8_C = (ex16o2>LC_8Min*1000 && ex16o2<LC_8Max*1000 && LC_8Max != 0);
                    bool isL9_C = (ex16o2>LC_9Min*1000 && ex16o2<LC_9Max*1000 && LC_9Max != 0);
                    bool isL10_C = (ex16o2>LC_10Min*1000 && ex16o2<LC_10Max*1000 && LC_10Max != 0);
                    bool isBG_C = !isL0_C && !isL1_C && !isL2_C && !isL3_C && !isL4_C && !isL5_C && !isL6_C && !isL7_C && !isL8_C && !isL9_C && !isL10_C;
                    if(isL0_C){
                        h1angdistL0_C->Fill(angleCoM2);
                        angdistL0_C->Fill(e16o2/1000.0, angleCoM2);
                    }
                    if(isL1_C){
                        h1angdistL1_C->Fill(angleCoM2);
                        angdistL1_C->Fill(e16o2/1000.0, angleCoM2);
                    }
                    if(isL2_C){
                        h1angdistL2_C->Fill(angleCoM2);
                        angdistL2_C->Fill(e16o2/1000.0, angleCoM2);
                    }
                    if(isL3_C){
                        h1angdistL3_C->Fill(angleCoM2);
                        angdistL3_C->Fill(e16o2/1000.0, angleCoM2);
                    }
                    if(isL4_C){
                        h1angdistL4_C->Fill(angleCoM2);
                        angdistL4_C->Fill(e16o2/1000.0, angleCoM2);
                    }
                    if(isL5_C){
                        h1angdistL5_C->Fill(angleCoM2);
                        angdistL5_C->Fill(e16o2/1000.0, angleCoM2);
                    }
                    if(isL6_C){
                        h1angdistL6_C->Fill(angleCoM2);
                        angdistL6_C->Fill(e16o2/1000.0, angleCoM2);
                    }
                    if(isL7_C){
                        h1angdistL7_C->Fill(angleCoM2);
                        angdistL7_C->Fill(e16o2/1000.0, angleCoM2);
                    }
                    if(isL8_C){
                        h1angdistL8_C->Fill(angleCoM2);
                        angdistL8_C->Fill(e16o2/1000.0, angleCoM2);
                    }
                    if(isL9_C){
                        h1angdistL9_C->Fill(angleCoM2);
                        angdistL9_C->Fill(e16o2/1000.0, angleCoM2);
                    }
                    if(isL10_C){
                        h1angdistL10_C->Fill(angleCoM2);
                        angdistL10_C->Fill(e16o2/1000.0, angleCoM2);
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
                        h216Od4->Fill(ex16o2/1000.0, ex13c/1000.0);

                        if (ex13c > C13Gate1) {
                            h1ex16o_c1->Fill(ex16o2/1000.0);
                        }

                        if (ex13c > C13Gate2) {
                            h1ex16o_c2->Fill(ex16o2/1000.0);
                        }

                        if (ex13c > C13Gate3) {
                            h1ex16o_c3->Fill(ex16o2/1000.0);
                        }

                        if (ex13c > C13Gate4) {
                            h1ex16o_c4->Fill(ex16o2/1000.0);
                        }

                        if (ex13c > C13Gate5) {
                            h1ex16o_c5->Fill(ex16o2/1000.0);
                        }

                        if (ex13c > C13Gate_best) {
                            h1ex16o_cbest->Fill(ex16o2/1000.0);
                        }
                        else h1angdist13C->Fill(angleCoM2);
                    }
                    
                    
                    // Reconstruct 5He
                    p5hex = hit[j].GetMomentumX() + pn2x;
                    p5hey = hit[j].GetMomentumY() + pn2y;
                    p5hez = hit[j].GetMomentumZ() + pn2z;
                    e5he = ((p5hex * p5hex) + (p5hey * p5hey) + (p5hez * p5hez))/2.0/5.012057224;
                    q5he = hit[j].GetEnergy() + eneu2 - e5he; // Q-Value = 735.00025 kev
                    ex5he = q5he - 735.00025;
                    
                    // Fill histogram
                    h1q5he->Fill(q5he/1000.0);
                    h2q5he->Fill(q5he/1000.0, (hit[j].GetEnergy() + eneu2)/1000.0);
                    if(qval2 < Q16Ohigh && qval2 > Q16Olow) {
                        h216Od3->Fill(ex16o2/1000.0, q5he/1000.0);
                        
                        if (q5he > He5Gate0 && q5he < He5Gate1) {}
                        else {
                            h1ex16o_h1->Fill(ex16o2/1000.0);
                        }
                        
                        if (q5he > He5Gate0 && q5he < He5Gate2) {}
                        else {
                            h1ex16o_h2->Fill(ex16o2/1000.0);
                        }
                        
                        if (q5he > He5Gate0 && q5he < He5Gate3) {}
                        else {
                            h1ex16o_h3->Fill(ex16o2/1000.0);
                        }

                        if (q5he > He5Gate0 && q5he < He5Gate4) {}
                        else {
                            h1ex16o_h4->Fill(ex16o2/1000.0);
                        }

                        if (q5he > He5Gate0 && q5he < He5Gate5) {}
                        else {
                            h1ex16o_h5->Fill(ex16o2/1000.0);
                        }

                        if (q5he > He5Gate0 && q5he < He5Gate6) {}
                        else {
                            h1ex16o_h6->Fill(ex16o2/1000.0);
                        }

                        if (q5he > He5Gate0 && q5he < He5Gate7) {}
                        else {
                            h1ex16o_h7->Fill(ex16o2/1000.0);
                        }

                        if (q5he > He5Gate0 && q5he < He5Gate_best) {
                            h1angdist5He->Fill(angleCoM2);
                        }
                        else {
                            h1ex16o_hbest->Fill(ex16o2/1000.0);
                        }
                    }

                    // Combined 5He and 13C
                    if(qval2 < Q16Ohigh && qval2 > Q16Olow) {
                        h1_12ca_corr->Fill(q16o2/1000.0);

                        if (ex13c >  C13Gate_best && ex5he > He5Gate_best) {

                            h1ex16o_2best->Fill(ex16o2/1000.0);

                            // Create the angular distributions from the data
                            if(isL0_He){
                                h1angdistL0_He->Fill(angleCoM2);
                                angdistL0_He->Fill(e16o2/1000.0, angleCoM2);
                            }
                            if(isL1_He){
                                h1angdistL1_He->Fill(angleCoM2);
                                angdistL1_He->Fill(e16o2/1000.0, angleCoM2);
                            }
                            if(isL2_He){
                                h1angdistL2_He->Fill(angleCoM2);
                                angdistL2_He->Fill(e16o2/1000.0, angleCoM2);
                            }
                            if(isL3_He){
                                h1angdistL3_He->Fill(angleCoM2);
                                angdistL3_He->Fill(e16o2/1000.0, angleCoM2);
                            }
                            if(isL4_He){
                                h1angdistL4_He->Fill(angleCoM2);
                                angdistL4_He->Fill(e16o2/1000.0, angleCoM2);
                            }
                            if(isL5_He){
                                h1angdistL5_He->Fill(angleCoM2);
                                angdistL5_He->Fill(e16o2/1000.0, angleCoM2);
                            }
                            if(isL6_He){
                                h1angdistL6_He->Fill(angleCoM2);
                                angdistL6_He->Fill(e16o2/1000.0, angleCoM2);
                            }
                            if(isL7_He){
                                h1angdistL7_He->Fill(angleCoM2);
                                angdistL7_He->Fill(e16o2/1000.0, angleCoM2);
                            }
                            if(isL8_He){
                                h1angdistL8_He->Fill(angleCoM2);
                                angdistL8_He->Fill(e16o2/1000.0, angleCoM2);
                            }
                            if(isL9_He){
                                h1angdistL9_He->Fill(angleCoM2);
                                angdistL9_He->Fill(e16o2/1000.0, angleCoM2);
                            }
                            if(isL10_He){
                                h1angdistL10_He->Fill(angleCoM2);
                                angdistL10_He->Fill(e16o2/1000.0, angleCoM2);
                            }

                            //Backgrounds
                            if(isL0_HeB){
                                h1angdistL0_HeB->Fill(angleCoM2);
                            }
                            if(isL1_HeB){
                                h1angdistL1_HeB->Fill(angleCoM2);
                            }
                            if(isL2_HeB){
                                h1angdistL2_HeB->Fill(angleCoM2);
                            }
                            if(isL3_HeB){
                                h1angdistL3_HeB->Fill(angleCoM2);
                            }
                            if(isL4_HeB){
                                h1angdistL4_HeB->Fill(angleCoM2);
                            }
                            if(isL5_HeB){
                                h1angdistL5_HeB->Fill(angleCoM2);
                            }
                            if(isL6_HeB){
                                h1angdistL6_HeB->Fill(angleCoM2);
                            }
                            if(isL7_HeB){
                                h1angdistL7_HeB->Fill(angleCoM2);
                            }
                            if(isL8_HeB){
                                h1angdistL8_HeB->Fill(angleCoM2);
                            }
                            if(isL9_HeB){
                                h1angdistL9_HeB->Fill(angleCoM2);
                            }
                            if(isL10_HeB){
                                h1angdistL10_HeB->Fill(angleCoM2);
                            }
                            if(isL11_HeB){
                                h1angdistL11_HeB->Fill(angleCoM2);
                            }
                        }

                        if (ex13c > cGateCombo0 && ex5he > hGateCombo00) {
                            h1ex16o_combined00->Fill(ex16o2/1000.0);
                        }
                        if (ex13c > cGateCombo0 && ex5he > hGateCombo01) {
                            h1ex16o_combined01->Fill(ex16o2/1000.0);
                        }
                        if (ex13c > cGateCombo0 && ex5he > hGateCombo02) {
                            h1ex16o_combined02->Fill(ex16o2/1000.0);
                        }
                        if (ex13c > cGateCombo0 && ex5he > hGateCombo03) {
                            h1ex16o_combined03->Fill(ex16o2/1000.0);
                        }

                        if (ex13c > cGateCombo1 && ex5he > hGateCombo10) {
                            h1ex16o_combined10->Fill(ex16o2/1000.0);
                        }
                        if (ex13c > cGateCombo1 && ex5he > hGateCombo11) {
                            h1ex16o_combined11->Fill(ex16o2/1000.0);
                        }
                        if (ex13c > cGateCombo1 && ex5he > hGateCombo12) {
                            h1ex16o_combined12->Fill(ex16o2/1000.0);
                        }
                        if (ex13c > cGateCombo1 && ex5he > hGateCombo12) {
                            h1ex16o_combined13->Fill(ex16o2/1000.0);
                        }

                        if (ex13c > cGateCombo2 && ex5he > hGateCombo20) {
                            h1ex16o_combined20->Fill(ex16o2/1000.0);
                        }
                        if (ex13c > cGateCombo2 && ex5he > hGateCombo21) {
                            h1ex16o_combined21->Fill(ex16o2/1000.0);
                        }
                        if (ex13c > cGateCombo2 && ex5he > hGateCombo22) {
                            h1ex16o_combined22->Fill(ex16o2/1000.0);
                        }
                        if (ex13c > cGateCombo2 && ex5he > hGateCombo23) {
                            h1ex16o_combined23->Fill(ex16o2/1000.0);
                        }
                    }
                }
            }
        }

        //cout << "Debug: Finished 16O Construction" << endl;

        // Fix event mixing 12C + a
        if(igood12C == 1 && igood > 3){
            //12C
            eMixingCa_Cx[em_countCa] = p12cx;
            eMixingCa_Cy[em_countCa] = p12cy;
            eMixingCa_Cz[em_countCa] = p12cz;
            eMixingCa_Ce[em_countCa] = e12c;

            //All the event alphas
            em_count1 = 0;

            for(int j = 0; j < igood; j++){
                if(j != n12ca1 && j != n12ca2 && j != n12ca3){
                    eMixingCa_ax[em_countCa][em_count1] = hit[j].GetMomentumX();
                    eMixingCa_ay[em_countCa][em_count1] = hit[j].GetMomentumY();
                    eMixingCa_az[em_countCa][em_count1] = hit[j].GetMomentumZ();
                    eMixingCa_ae[em_countCa][em_count1] = hit[j].GetEnergy();
                    em_count1++;
                }
            }
            if(em_countCa < 10) em_countCa++;
            igoodevCa[em_countCa] = igood;
        }

        //cout << "Debug: Finished 12C+a filling" << endl;


        // Event Mixing 3a -> 12C Possiblities
        // Filling arrays without bias
        if (igood>2) {
            for (Int_t k=0; k<igood; k++) {
                eMixing3a_px[em_countC][k] = hit[k].GetMomentumX();
                eMixing3a_py[em_countC][k] = hit[k].GetMomentumY();
                eMixing3a_pz[em_countC][k] = hit[k].GetMomentumZ();
                eMixing3a_e[em_countC][k] = hit[k].GetEnergy();
            }

            igoodev[em_countC] = igood;

            em_countC++;
        }

        if(igood > 3){
            for(Int_t k = 0; k < igood; k++){
                eMixing4a_px[em_countC4][k] = hit[k].GetMomentumX();
                eMixing4a_py[em_countC4][k] = hit[k].GetMomentumY();
                eMixing4a_pz[em_countC4][k] = hit[k].GetMomentumZ();
                eMixing4a_e[em_countC4][k] = hit[k].GetEnergy();
            }

            igoodev4[em_countC4] = igood;
            em_countC4++;
        }


        // Now create the mixed 12C
        if(em_countC == 3){
            em_countC = 0;
            em_count2 = 0;

            for(int j = 0; j < igoodev[0]; j++){
                for(int k = 0; k < igoodev[1]; k++){
                    for(int l = 0; l < igoodev[2]; l++){
                        eMixingC_px = eMixing3a_px[0][j] + eMixing3a_px[1][k] + eMixing3a_px[2][l];
                        eMixingC_py = eMixing3a_py[0][j] + eMixing3a_py[1][k] + eMixing3a_py[2][l];
                        eMixingC_pz = eMixing3a_pz[0][j] + eMixing3a_pz[1][k] + eMixing3a_pz[2][l];
                        eMixingC_e = ((eMixingC_px*eMixingC_px) + (eMixingC_py*eMixingC_py) + (eMixingC_pz*eMixingC_pz))/2.0/12.0; 
                        eMixingC_q = eMixing3a_e[0][j] + eMixing3a_e[1][k] + eMixing3a_e[2][l] - eMixingC_e;
                        eMixingC_ex = (eMixing3a_e[0][j] + eMixing3a_e[1][k] + eMixing3a_e[2][l]) - eMixingC_e + 7274.7;

                        //h1ex12c_em->Fill(eMixingC_ex/1000.0);
                        //if (eMixingC_q > 0) h1erel12c_em->Fill(eMixingC_q/1000.0);
                        
                        em_count2++;
                    }
                }
            }

            // Zero the arrays
            for(int j = 0; j < 3; j++){
                for(int k = 0; k < 16; k++){
                    eMixing3a_px[j][k] = 0.0;
                    eMixing3a_py[j][k] = 0.0;
                    eMixing3a_pz[j][k] = 0.0;
                    eMixing3a_e[j][k] = 0.0;
                }

                igoodev[j] = 0;
            }
            
            eMixingC_px = 0.0;
            eMixingC_py = 0.0;
            eMixingC_pz = 0.0;
            eMixingC_e = 0.0;
            eMixingC_e1 = 0.0;
            eMixingC_q = 0.0;
            eMixingC_ex = 0.0;
        }

        // Event Mixing 12C w/ 4 alpha reaction q-value
        if(em_countC4 == 4){

            em_countC4 = 0;

            for(int j = 0; j < igoodev4[0]; j++){
                for(int k = 0; k < igoodev4[1]; k++){
                    for(int l = 0; l < igoodev4[2]; l++){
                        for(int m = 0; m < igoodev4[3]; m++){
                            
                            for (int n = 0; n < 4; n++){
                                Int_t imod0 = n % 4;
                                Int_t imod1 = (n+1) % 4;
                                Int_t imod2 = (n+2) % 4;
                                Int_t imod3 = (n+3) % 4;

                                eMixing4a_cx = eMixing4a_px[imod0][j] + eMixing4a_px[imod1][k] + eMixing4a_px[imod2][l];
                                eMixing4a_cy = eMixing4a_py[imod0][j] + eMixing4a_py[imod1][k] + eMixing4a_py[imod2][l];
                                eMixing4a_cz = eMixing4a_pz[imod0][j] + eMixing4a_pz[imod1][k] + eMixing4a_pz[imod2][l];
                                eMixing4a_ce = ((eMixing4a_cx*eMixing4a_cx) + (eMixing4a_cy*eMixing4a_cy) + (eMixing4a_cz*eMixing4a_cz))/2.0/12.0;
                                eMixing4a_q = (eMixing4a_e[imod0][j] + eMixing4a_e[imod1][k] + eMixing4a_e[imod2][l]) - eMixing4a_ce;

                                pn2x4 = -(eMixing4a_cx + eMixing4a_px[imod3][m]);
                                pn2y4 = -(eMixing4a_cy + eMixing4a_py[imod3][m]);
                                pn2z4 = pbeam - (eMixing4a_cz + eMixing4a_pz[imod3][m]);
                                eneu4 = (pn2x4*pn2x4+pn2y4*pn2y4+pn2z4*pn2z4)/2.0/1.0086649;
                                        etot4 = eMixing4a_e[imod0][j] + eMixing4a_e[imod1][k] + eMixing4a_e[imod2][l] + eMixing4a_e[imod3][m] + eneu4;
                                qval4 = etot4 - ebeam;                        
                                        // 13C(a, 4a + n) Q-value = -12221.055
                                if(qval4 < Q16Ohigh && qval4 > Q16Olow){
                                    h1erel12c_em->Fill(eMixing4a_q/1000.0);
                                }
                            }
                        }
                    }
                }
            }
        }

        //cout << "Debug: Before EM Reconstruct" << endl;

        // Event Mixing 12C + a
        if(em_countCa == 10){
            em_countCa = 0;

            //cout << "Debug: Started EM 12C+a" << endl;

            // Constructing 16O
            for(int j = 0; j < 10; j++){
                for(int k = 0; k < 10; k++){
                    for(int l = 0; l < igoodevCa[k]; l++){
                        if(j != k){
                            // Mix from different events the 12C and a
                            eMixingO_px1 = eMixingCa_Cx[j] + eMixingCa_ax[k][l];
                            eMixingO_py1 = eMixingCa_Cy[j] + eMixingCa_ay[k][l];
                            eMixingO_pz1 = eMixingCa_Cz[j] + eMixingCa_az[k][l];
                            eMixingO_e_1 = ((eMixingO_px1*eMixingO_px1) + (eMixingO_py1*eMixingO_py1) + (eMixingO_pz1*eMixingO_pz1))/2.0/15.9949146;
                            eMixingO_q1 = eMixingCa_Ce[j] + eMixingCa_ae[k][l] - eMixingO_e_1;
                            eMixingO_ex1 = eMixingO_q1 + 7161.91698 + 7654.20;
                            eMixingO_e11 = eMixingCa_Ce[j] + eMixingCa_ae[k][l];
    
                            // reconstruct neutron
                            eMixingN_e1 = ((eMixingO_px1*eMixingO_px1) + (eMixingO_py1*eMixingO_py1) + ((pbeam - eMixingO_pz1)*(pbeam - eMixingO_pz1)))/2.0/1.0086649;
                            eMixingQ1 = eMixingO_e11 + eMixingN_e1 - ebeam;
    
                            //Angle of 16O in lab frame with radians
                            Double_t tgThetaLab_em1 = sqrt(pow(eMixingO_px1,2) + pow(eMixingO_py1,2))/eMixingO_pz1;
                            Double_t cosThetaLab_em1;
                            if(tgThetaLab_em1 >= 0)
                                cosThetaLab_em1 = 1/sqrt(1+pow(tgThetaLab_em1,2));
                            else
                                cosThetaLab_em1 = -1/sqrt(1+pow(tgThetaLab_em1,2));
                            
                            //Velocity of 16O in lab frame
                            Double_t v16OLab_em1 = sqrt(pow(eMixingO_px1,2)+pow(eMixingO_py1,2)+pow(eMixingO_pz1,2))/15.9949146;
                            
                            //Velocity of center of momentum frame in final stage
                            Double_t VCoM_em1 = (4.00260325415/(1.00866491574 + 15.99491461956))*vbeam;
                            
                            //Velocity of 16O in center of mass frame
                            Double_t V16OCM_em1 = sqrt(pow(v16OLab_em1,2)+pow(VCoM_em1,2)-2*v16OLab_em1*VCoM_em1*cosThetaLab_em1);
                            
                            //Angle of 16O in CoM
                            Double_t angleCoM_em1 = acos((v16OLab_em1*cosThetaLab_em1 - VCoM_em1)/V16OCM_em1)*180/3.14159265359;
            
                            // If mixed "reaction" q-value matches, fill histograms
                            if(eMixingQ1 > Q16Olow && eMixingQ1 < Q16Ohigh){
                                h1ex16o_em1->Fill(eMixingO_q1/1000.0);
                                h1angdistBG_C->Fill(angleCoM_em1);
                                h1angdistBG_He->Fill(angleCoM_em1);
                                angdistBG_C->Fill(eMixingO_q1/1000.0, angleCoM_em1);
                                angdistBG_He->Fill(eMixingO_q1/1000.0, angleCoM_em1);
            
                                bool isL0_He = (eMixingO_ex1>LHe_0Min*1000 && eMixingO_ex1<LHe_0Max*1000 && LHe_0Max != 0);
                                bool isL1_He = (eMixingO_ex1>LHe_1Min*1000 && eMixingO_ex1<LHe_1Max*1000 && LHe_1Max != 0);
                                bool isL2_He = (eMixingO_ex1>LHe_2Min*1000 && eMixingO_ex1<LHe_2Max*1000 && LHe_2Max != 0);
                                bool isL3_He = (eMixingO_ex1>LHe_3Min*1000 && eMixingO_ex1<LHe_3Max*1000 && LHe_3Max != 0);
                                bool isL4_He = (eMixingO_ex1>LHe_4Min*1000 && eMixingO_ex1<LHe_4Max*1000 && LHe_4Max != 0);
                                bool isL5_He = (eMixingO_ex1>LHe_5Min*1000 && eMixingO_ex1<LHe_5Max*1000 && LHe_5Max != 0);
                                bool isL6_He = (eMixingO_ex1>LHe_6Min*1000 && eMixingO_ex1<LHe_6Max*1000 && LHe_6Max != 0);
                                bool isL7_He = (eMixingO_ex1>LHe_7Min*1000 && eMixingO_ex1<LHe_7Max*1000 && LHe_7Max != 0);
                                bool isL8_He = (eMixingO_ex1>LHe_8Min*1000 && eMixingO_ex1<LHe_8Max*1000 && LHe_8Max != 0);
                                bool isL9_He = (eMixingO_ex1>LHe_9Min*1000 && eMixingO_ex1<LHe_9Max*1000 && LHe_9Max != 0);
                                bool isL10_He = (eMixingO_ex1>LHe_10Min*1000 && eMixingO_ex1<LHe_10Max*1000 && LHe_10Max != 0);
                                bool isBG_He = !isL0_He && !isL1_He && !isL2_He && !isL3_He && !isL4_He && !isL5_He && !isL6_He && !    isL7_He && !isL8_He && !isL9_He && !isL10_He;
                                if(isL0_He){
                                    h1EMangdistL0_He->Fill(angleCoM_em1);
                                }
                                if(isL1_He){
                                    h1EMangdistL1_He->Fill(angleCoM_em1);
                                }
                                if(isL2_He){
                                    h1EMangdistL2_He->Fill(angleCoM_em1);
                                }
                                if(isL3_He){
                                    h1EMangdistL3_He->Fill(angleCoM_em1);
                                }
                                if(isL4_He){
                                    h1EMangdistL4_He->Fill(angleCoM_em1);
                                }
                                if(isL5_He){
                                    h1EMangdistL5_He->Fill(angleCoM_em1);
                                }
                                if(isL6_He){
                                    h1EMangdistL6_He->Fill(angleCoM_em1);
                                }
                                if(isL7_He){
                                    h1EMangdistL7_He->Fill(angleCoM_em1);
                                }
                                if(isL8_He){
                                    h1EMangdistL8_He->Fill(angleCoM_em1);
                                }
                                if(isL9_He){
                                    h1EMangdistL9_He->Fill(angleCoM_em1);
                                }
                                if(isL10_He){
                                    h1EMangdistL10_He->Fill(angleCoM_em1);
                                }
            
                                h1q_em->Fill(eMixingQ1/1000.0);
                                h1eneu_em->Fill(eMixingN_e1/1000.0);
                                h1e1_em->Fill(eMixingO_e_1/1000.0);
                            }
                        }
                    }
                }
            }

            //cout << "Debug: Finished Event Mixed 16O Construction" << endl;

            // Clearing arrays
            for(int j=0; j < 10; j++){
                eMixingCa_Cx[j] = 0.0;
                eMixingCa_Cy[j] = 0.0;
                eMixingCa_Cz[j] = 0.0;
                eMixingCa_Ce[j] = 0.0;
                for(int k=0; k < igoodevCa[j]; k++){
                    eMixingCa_ax[j][k] = 0.0;
                    eMixingCa_ay[j][k] = 0.0;
                    eMixingCa_az[j][k] = 0.0;
                    eMixingCa_ae[j][k] = 0.0;
                }
            }

            //cout << "Debug: Cleared Arrays" << endl;
            

            eMixingO_px1 = 0.0;
            eMixingO_py1 = 0.0;
            eMixingO_pz1 = 0.0;
            eMixingO_e_1 = 0.0;
            eMixingO_e11 = 0.0;
            eMixingO_q1 = 0.0;
            eMixingO_ex1 = 0.0;
            eMixingN_e1 = 0.0;
            eMixingQ1 = 0.0;
            
        }

        
        h1m->Fill(igood);
    }

    h1_8be_corr->Divide(h1ex16o_em);
    h1_12ca_corr->Divide(h1ex16o_em1);
    h1_12c_corr->Divide(h1erel12c_em);

    cout << igoodMax << endl;
    
    file.Write();
    file.Close();

    ftext.close(); 
    fhits.close();

    clock_t stop=clock();
    printf("Analysis Run Time: %10.3f (seconds) \n",(double) (stop-start)/CLOCKS_PER_SEC);
}
    