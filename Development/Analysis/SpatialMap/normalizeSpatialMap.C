void normalizeSpatialMap(TString mapFile, Int_t nrSamples=100000){
	// Input (unnormalized maps)
	TFile * file = TFile::Open(mapFile, "READ");
	TH1D * prInnerDet = (TH1D*) file->Get("PrInnerDet");
	TH2D * innerMap = (TH2D*) file->Get("InnerMap");
	TH2D * outerMap = (TH2D*) file->Get("OuterMap");
	// Parameters
	double min_r=innerMap->GetXaxis()->GetXmin(), max_r=innerMap->GetXaxis()->GetXmax();
        double min_angle=innerMap->GetYaxis()->GetXmin(), max_angle=innerMap->GetYaxis()->GetXmax();
	int rbins=innerMap->GetXaxis()->GetNbins(), angle_bins=innerMap->GetYaxis()->GetNbins();
	// Output
	TFile * output = TFile::Open("Normalized_" + mapFile, "RECREATE");
	TH1D * normPrInnerD = new TH1D("PrInnerDet", "Pr ~ Fract. Inner/(Inner+Outer) Detections", rbins, min_r, max_r);
        TH2D * normInnerMap = new TH2D("InnerMap", "R-Angle Inner Shroud Map", rbins, min_r, max_r, angle_bins, min_angle, max_angle);
        TH2D * normOuterMap = new TH2D("OuterMap", "R-Angle Outer Shroud Map", rbins, min_r, max_r, angle_bins, min_angle, max_angle);
	// Create normalized maps
	for(int r=min_r; r<=max_r; r++){
		std::cout << "\rR: " << r << std::flush;
		// Copy value of PrInnerDet map
		double pr = prInnerDet->GetBinContent(prInnerDet->GetBin(r));
		normPrInnerD->Fill(r, pr);
		// Inner map normalization
		TH1D * innerProjR = innerMap->ProjectionY("py", r, r+1);
		for(int i=0; i<nrSamples; i++){
			double angle = innerProjR->GetRandom();
			normInnerMap->Fill(r, angle);
		}
		// Outer map normalization
		TH1D * outerProjR = outerMap->ProjectionY("py", r, r+1);
		for(int i=0; i<nrSamples; i++){
			double angle = outerProjR->GetRandom();
			normOuterMap->Fill(r, angle);
		}
	}
	std::cout << std::endl;
	// Write normalized maps
	normPrInnerD->Write();
	normInnerMap->Write();
	normOuterMap->Write();
	output->Close();
	file->Close();
}
