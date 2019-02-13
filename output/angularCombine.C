/*
	Program to do angular distribution backgound subtraction and normalization. 
	For each energy, the user needs to change IO file names, histogram write names, and user input the data about the peak/background statistics in order to properly normalize the histograms
	Created by Bryce Frentz
	Updated 10/27/17
*/

{
	// Open the root file and get the histograms
	// USER CHANGE: Beam Energy
	TFile *f = new TFile("./tenDegree/data24MeV_correlation.root");
	TH1F *h0 =  (TH1F*)f->Get("h1angdistL0_Be");						// L0 = 15.36
	TH1F *h1 =  (TH1F*)f->Get("h1angdistL1_Be");						// L1 = 16.46
	TH1F *h2 =  (TH1F*)f->Get("h1angdistL2_Be");						// L2 = 17.20
	TH1F *h3 =  (TH1F*)f->Get("h1angdistL3_Be");						// L3 = 17.50
	TH1F *h4 =  (TH1F*)f->Get("h1angdistL4_Be");						// L4 = 18.10
	TH1F *h5 =  (TH1F*)f->Get("h1angdistL5_Be");						// L5 = 18.80
	TH1F *h6 =  (TH1F*)f->Get("h1angdistL6_Be");						// L6 = 19.40
	TH1F *h7 =  (TH1F*)f->Get("h1angdistL7_Be");						// L7 = 21.00
	TH1F *h8 =  (TH1F*)f->Get("h1angdistL8_Be");						// L8 = 21.70
	TH1F *h9 =  (TH1F*)f->Get("h1angdistL9_Be");						// L9 = 22.34
	TH1F *h10 = (TH1F*)f->Get("h1angdistL10_Be");				    	// L10 = 23.30
	TH1F *hbg = (TH1F*)f->Get("h1angdist9Be");
	TH1F *hL0 = (TH1F*) h0->Clone();
	TH1F *hL1 = (TH1F*) h1->Clone();
	TH1F *hL2 = (TH1F*) h2->Clone();
	TH1F *hL3 = (TH1F*) h3->Clone();
	TH1F *hL4 = (TH1F*) h4->Clone();
	TH1F *hL5 = (TH1F*) h5->Clone();
	TH1F *hL6 = (TH1F*) h6->Clone();
	TH1F *hL7 = (TH1F*) h7->Clone();
	TH1F *hL8 = (TH1F*) h8->Clone();
	TH1F *hL9 = (TH1F*) h9->Clone();
	TH1F *hL10 = (TH1F*) h10->Clone();
	TH1F *h9Be = (TH1F*) hbg->Clone();
	

	// Create output file
	TFile *fout = TFile::Open("./compareDistributions24.root","UPDATE");

	// Set new histogram names for write out
	// USER CHANGE: Beam Energy
	hL0->SetName( "h24_tenDegree_L0");				// copied
	hL1->SetName( "h24_tenDegree_L1");
	hL2->SetName( "h24_tenDegree_L2");
	hL3->SetName( "h24_tenDegree_L3");
	hL4->SetName( "h24_tenDegree_L4");
	hL5->SetName( "h24_tenDegree_L5");
	hL6->SetName( "h24_tenDegree_L6");
	hL7->SetName( "h24_tenDegree_L7");
	hL8->SetName( "h24_tenDegree_L8");
	hL9->SetName( "h24_tenDegree_L9");
	hL10->SetName("h24_tenDegree_L10");
	

	// Write out
	hL0->Write( "h24_tenDegree_L0",  TObject::kWriteDelete);				
	hL1->Write( "h24_tenDegree_L1",  TObject::kWriteDelete);				
	hL2->Write( "h24_tenDegree_L2",  TObject::kWriteDelete);
	hL3->Write( "h24_tenDegree_L3",  TObject::kWriteDelete);
	hL4->Write( "h24_tenDegree_L4",  TObject::kWriteDelete);
	hL5->Write( "h24_tenDegree_L5",  TObject::kWriteDelete);
	hL6->Write( "h24_tenDegree_L6",  TObject::kWriteDelete);
	hL7->Write( "h24_tenDegree_L7",  TObject::kWriteDelete);
	hL8->Write( "h24_tenDegree_L8",  TObject::kWriteDelete);
	hL9->Write( "h24_tenDegree_L9",  TObject::kWriteDelete);
	hL10->Write("h24_tenDegree_L10", TObject::kWriteDelete);
	
}