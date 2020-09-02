#define LIM_PRINTOUT 10000
void check_edep_vs_starting_energy(const char * fileIn, const char *fileOut="output_startingEnergy_vs_Edep"){
	int event, last_event=0, PID;
	Double_t kineticenergy=0, energydeposition=0, sum_edep=0, starting_e=0;

	TFile *input = new TFile(fileIn, "READ");
	TTree *fTree = (TTree*) input->Get("fTree");

	fTree->SetBranchAddress("PID", &PID);
	fTree->SetBranchAddress("energydeposition", &energydeposition);
	fTree->SetBranchAddress("kineticenergy", &kineticenergy);
	fTree->SetBranchAddress("eventnumber", &event);

	TFile *output = new TFile(fileOut, "RECREATE");
	TTree *oTree = new TTree("fTree", "");
	oTree->Branch("sum_edep", &sum_edep);
	oTree->Branch("starting_e", &starting_e);

	int nentries = fTree->GetEntries();
	for(int i=0; i<nentries; i++){
		fTree->GetEntry(i);
		if(i%LIM_PRINTOUT==0){
			cout << "\r" << i << "/" << nentries << flush;
		}
		if(last_event>0 && event!=last_event){
			//write
			oTree->Fill();
			//reset
			sum_edep = 0;
			starting_e = 0;
			//take starting energy
			if(PID!=2112)
				cout << "WARNING - Initial energy taken from not neutron\n";
			starting_e = kineticenergy;
		}
		sum_edep += energydeposition;
		last_event = event;
	}
	oTree->Fill();
	oTree->Write();
	output->Close();
}
