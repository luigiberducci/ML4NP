void plot_edep_gbyevent(const char * fileIn){
	TFile *file = new TFile(fileIn, "READ");
	TTree *fTree = (TTree*) file->Get("fTree");
	double energydeposition;
	int eventnumber;
	fTree->SetBranchAddress("energydeposition", &energydeposition);
	fTree->SetBranchAddress("eventnumber", &eventnumber);

	double sum_edep = 0;
	int last_event = 0;
	TFile *output = new TFile("sum_edep_gby_event.root", "recreate");
	TTree *oTree = new TTree("fTree", "");
	oTree->Branch("energydeposition", &sum_edep);
	for(int i=0; i<fTree->GetEntries(); i++){
		if(i%1000000==0)
			cout << "\r" << i << "/" << fTree->GetEntries() << flush;
		fTree->GetEntry(i);
		if(last_event>0 && eventnumber!=last_event){
			// write
			oTree->Fill();
			// reset
			sum_edep = 0;
		}
		if(energydeposition>0)
			sum_edep += energydeposition;
		last_event = eventnumber;
	}
	cout << "\r" << fTree->GetEntries() << "/" << fTree->GetEntries() << flush;
	cout << endl;
	// last fill and write
	oTree->Fill();
	oTree->Write();
	output->Close();
	file->Close();
}
