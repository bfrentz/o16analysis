/*
	Program to do angular distribution backgound subtraction and normalization. 
	For each energy, the user needs to change IO file names, histogram write names, and user input the data about the peak/background statistics in order to properly normalize the histograms
	Created by Bryce Frentz
	Updated 10/27/17
*/

{
	// 24MeV: 2, 4, 5, 6
	// 25MeV: 2, 4, 5, 6
	// 27MeV: 2, 3, 4, 5, 6, 7, 8
	// 28MeV: not included, too little statistics
	// 29MeV: 2, 4, 5, 6, 7, 8, 10

	// Open the root file and get the histograms
	// USER CHANGE: Beam Energy
	TFile *f = new TFile("./twoDegree/data29MeV_correlation.root");
	//TH1F *h0 =  (TH1F*)f->Get("h1angdistL0_Be");						// L0 = 15.36
	//TH1F *h1 =  (TH1F*)f->Get("h1angdistL1_Be");						// L1 = 16.46
	TH1F *h2 =  (TH1F*)f->Get("h1angdistL2_Be");						// L2 = 17.20
	//TH1F *h3 =  (TH1F*)f->Get("h1angdistL3_Be");						// L3 = 17.50
	TH1F *h4 =  (TH1F*)f->Get("h1angdistL4_Be");						// L4 = 18.10
	TH1F *h5 =  (TH1F*)f->Get("h1angdistL5_Be");						// L5 = 18.80
	TH1F *h6 =  (TH1F*)f->Get("h1angdistL6_Be");						// L6 = 19.40
	TH1F *h7 =  (TH1F*)f->Get("h1angdistL7_Be");						// L7 = 21.00
	TH1F *h8 =  (TH1F*)f->Get("h1angdistL8_Be");						// L8 = 21.70
	//TH1F *h9 =  (TH1F*)f->Get("h1angdistL9_Be");						// L9 = 22.34
	TH1F *h10 = (TH1F*)f->Get("h1angdistL10_Be");				    	// L10 = 23.30
	//TH1F *hbg = (TH1F*)f->Get("h1angdist9Be");
	//TH1F *hL0 = (TH1F*) h0->Clone();
	//TH1F *hL1 = (TH1F*) h1->Clone();
	TH1F *hL2 = (TH1F*) h2->Clone();
	//TH1F *hL3 = (TH1F*) h3->Clone();
	TH1F *hL4 = (TH1F*) h4->Clone();
	TH1F *hL5 = (TH1F*) h5->Clone();
	TH1F *hL6 = (TH1F*) h6->Clone();
	TH1F *hL7 = (TH1F*) h7->Clone();
	TH1F *hL8 = (TH1F*) h8->Clone();
	//TH1F *hL9 = (TH1F*) h9->Clone();
	TH1F *hL10 = (TH1F*) h10->Clone();
	//TH1F *h9Be = (TH1F*) hbg->Clone();

	//Rebin
	//hL0->Rebin(5);
	//hL1->Rebin(5);
	hL2->Rebin(5);
	//hL3->Rebin(5);
	hL4->Rebin(5);
	hL5->Rebin(5);
	hL6->Rebin(5);
	hL7->Rebin(5);
	hL8->Rebin(5);
	//hL9->Rebin(5);
	hL10->Rebin(5);
	
	// Get bin data
	Int_t distribution[18][8];

	for(int i = 0; i < 18; i++){
		distribution[i][0] = (int)round(hL2->GetBinCenter(i+1)) - 1;
	}

	for(int i = 0; i < 18; i++){
		distribution[i][1] = hL2->GetBinContent(i+1);
	}
	for(int i = 0; i < 18; i++){
		distribution[i][2] = hL4->GetBinContent(i+1);
	}
	for(int i = 0; i < 18; i++){
		distribution[i][3] = hL5->GetBinContent(i+1);
	}
	for(int i = 0; i < 18; i++){
		distribution[i][4] = hL6->GetBinContent(i+1);
	}
	for(int i = 0; i < 18; i++){
		distribution[i][5] = hL7->GetBinContent(i+1);
	}
	for(int i = 0; i < 18; i++){
		distribution[i][6] = hL8->GetBinContent(i+1);
	}
	for(int i = 0; i < 18; i++){
		distribution[i][7] = hL10->GetBinContent(i+1);
	}

	// Output file for bin contents
	// USER CHANGE energy
    std::ofstream fOut("./distributionData29.csv", ios::app);

    // Write data
    for(int i = 0; i < 18; i++){
    	for(int j = 0; j < 8; j++){
    		fOut << distribution[i][j] << ",";
    	}

    	fOut << endl;
    }

    // Close the file
    fOut.close();


    // Create TGraphs for the data
    TGraphErrors *g2 = new TGraphErrors(18);
    g2->SetTitle("Angular Distribution (E_{beam}=29MeV, 17.2MeV state)");
	g2->GetXaxis()->SetTitle("Center of Mass Angle");
	g2->GetYaxis()->SetTitle("Counts per 10 degrees");
	g2->SetMarkerStyle(20);
	g2->SetMarkerSize(1);

	/*
	TGraphErrors *g3 = new TGraphErrors(18);
    g3->SetTitle("Angular Distribution (E_{beam}=27MeV, 17.5MeV state)");
	g3->GetXaxis()->SetTitle("Center of Mass Angle");
	g3->GetYaxis()->SetTitle("Counts per 10 degrees");
	g3->SetMarkerStyle(20);
	g3->SetMarkerSize(1);
	*/

    TGraphErrors *g4 = new TGraphErrors(18);
    g4->SetTitle("Angular Distribution (E_{beam}=29MeV, 18.1MeV state)");
	g4->GetXaxis()->SetTitle("Center of Mass Angle");
	g4->GetYaxis()->SetTitle("Counts per 10 degrees");
	g4->SetMarkerStyle(20);
	g4->SetMarkerSize(1);

    TGraphErrors *g5 = new TGraphErrors(18);
    g5->SetTitle("Angular Distribution (E_{beam}=29MeV, 18.8MeV state)");
	g5->GetXaxis()->SetTitle("Center of Mass Angle");
	g5->GetYaxis()->SetTitle("Counts per 10 degrees");
	g5->SetMarkerStyle(20);
	g5->SetMarkerSize(1);

    TGraphErrors *g6 = new TGraphErrors(18);
    g6->SetTitle("Angular Distribution (E_{beam}=29MeV, 19.4MeV state)");
	g6->GetXaxis()->SetTitle("Center of Mass Angle");
	g6->GetYaxis()->SetTitle("Counts per 10 degrees");
	g6->SetMarkerStyle(20);
	g6->SetMarkerSize(1);

	TGraphErrors *g7 = new TGraphErrors(18);
    g7->SetTitle("Angular Distribution (E_{beam}=29MeV, 21.0MeV state)");
	g7->GetXaxis()->SetTitle("Center of Mass Angle");
	g7->GetYaxis()->SetTitle("Counts per 10 degrees");
	g7->SetMarkerStyle(20);
	g7->SetMarkerSize(1);

    TGraphErrors *g8 = new TGraphErrors(18);
    g8->SetTitle("Angular Distribution (E_{beam}=29MeV, 21.7MeV state)");
	g8->GetXaxis()->SetTitle("Center of Mass Angle");
	g8->GetYaxis()->SetTitle("Counts per 10 degrees");
	g8->SetMarkerStyle(20);
	g8->SetMarkerSize(1);

	TGraphErrors *g10 = new TGraphErrors(18);
    g10->SetTitle("Angular Distribution (E_{beam}=29MeV, 23.3MeV state)");
	g10->GetXaxis()->SetTitle("Center of Mass Angle");
	g10->GetYaxis()->SetTitle("Counts per 10 degrees");
	g10->SetMarkerStyle(20);
	g10->SetMarkerSize(1);


    for(int i = 0; i < 18; i++){
    	g2->SetPoint(i, distribution[i][0], distribution[i][1]);
    	g2->SetPointError(i, 0.05, sqrt(distribution[i][1]));
    }
    /*
    for(int i = 0; i < 18; i++){
    	g3->SetPoint(i, distribution[i][0], distribution[i][2]);
    	g3->SetPointError(i, 0.05, sqrt(distribution[i][2]));
    }
    */
    for(int i = 0; i < 18; i++){
    	g4->SetPoint(i, distribution[i][0], distribution[i][2]);
    	g4->SetPointError(i, 0.05, sqrt(distribution[i][2]));
    }
    for(int i = 0; i < 18; i++){
    	g5->SetPoint(i, distribution[i][0], distribution[i][3]);
    	g5->SetPointError(i, 0.05, sqrt(distribution[i][3]));
    }
    for(int i = 0; i < 18; i++){
    	g6->SetPoint(i, distribution[i][0], distribution[i][4]);
    	g6->SetPointError(i, 0.05, sqrt(distribution[i][4]));
    }
    for(int i = 0; i < 18; i++){
    	g7->SetPoint(i, distribution[i][0], distribution[i][5]);
    	g7->SetPointError(i, 0.05, sqrt(distribution[i][5]));
    }
    for(int i = 0; i < 18; i++){
    	g8->SetPoint(i, distribution[i][0], distribution[i][6]);
    	g8->SetPointError(i, 0.05, sqrt(distribution[i][6]));
    }
    for(int i = 0; i < 18; i++){
    	g10->SetPoint(i, distribution[i][0], distribution[i][7]);
    	g10->SetPointError(i, 0.05, sqrt(distribution[i][7]));
    }

	// Create output file
	TFile *fOutRoot = TFile::Open("./finalAngularDistributions.root","UPDATE");


	// Set new names for write out
	// USER CHANGE: Beam Energy
	//hL0->SetName( "h24_L0");				// copied
	//hL1->SetName( "h24_L1");
	hL2->SetName( "h29_L2");
	//hL3->SetName( "h29_L3");
	hL4->SetName( "h29_L4");
	hL5->SetName( "h29_L5");
	hL6->SetName( "h29_L6");
	hL7->SetName( "h29_L7");
	hL8->SetName( "h29_L8");
	//hL9->SetName( "h29_L9");
	hL10->SetName("h29_L10");

	g2->SetName( "g29MeV_L2");
	//g3->SetName( "g29MeV_L3");
	g4->SetName( "g29MeV_L4");
	g5->SetName( "g29MeV_L5");
	g6->SetName( "g29MeV_L6");
	g7->SetName( "g29MeV_L7");
	g8->SetName( "g29MeV_L8");
	g10->SetName("g29MeV_L10");
	

	// Write out
	//hL0->Write( "h24_L0",  TObject::kWriteDelete);				
	//hL1->Write( "h24_L1",  TObject::kWriteDelete);				
	hL2->Write( "h29_L2",  TObject::kWriteDelete);
	//hL3->Write( "h29_L3",  TObject::kWriteDelete);
	hL4->Write( "h29_L4",  TObject::kWriteDelete);
	hL5->Write( "h29_L5",  TObject::kWriteDelete);
	hL6->Write( "h29_L6",  TObject::kWriteDelete);
	hL7->Write( "h29_L7",  TObject::kWriteDelete);
	hL8->Write( "h29_L8",  TObject::kWriteDelete);
	//hL9->Write( "h29_L9",  TObject::kWriteDelete);
	hL10->Write("h29_L10", TObject::kWriteDelete);

	g2->Write( "g29MeV_L2", TObject::kWriteDelete);
	//g3->Write( "g29MeV_L3", TObject::kWriteDelete);
	g4->Write( "g29MeV_L4", TObject::kWriteDelete);
	g5->Write( "g29MeV_L5", TObject::kWriteDelete);
	g6->Write( "g29MeV_L6", TObject::kWriteDelete);
	g7->Write( "g29MeV_L7", TObject::kWriteDelete);
	g8->Write( "g29MeV_L8", TObject::kWriteDelete);
	g10->Write("g29MeV_L10", TObject::kWriteDelete);

	fOutRoot->Close();
	
}