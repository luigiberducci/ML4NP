bool isRootFile(TString fileName, TString prefix){
    // Consider files `prefixXXXXXX.root`
    TString begin = fileName(0, prefix.Length());
    if(begin.EqualTo(prefix) & fileName.EndsWith("root"))
        return true;
    return false;
}

void plot_LArDeposition_Gecoincidence(const char * dirIn="./", const char * dirOut="./", int geThreshold=200){
	char* fullDirIn = gSystem->ExpandPathName(dirIn);
	char* fullDirOut = gSystem->ExpandPathName(dirOut);
	TString outputfilename(dirOut);
	outputfilename += "GeCoincidence_Spectra.root";
	TFile * output = TFile::Open(outputfilename, "RECREATE");
	// Create histograms
	double totalargonenergy = 0, totalgeenergy = 0;
	int totalpe = 0;
	TTree * oTree = new TTree("oTree", "Event-level information");
	oTree->Branch("energy", &totalargonenergy);
	oTree->Branch("NPE", &totalpe);
	oTree->Branch("GeThreshold", &geThreshold);

	void* dirp = gSystem->OpenDirectory(fullDirIn);
	const char* entry;
	int kFiles = 0, kProducedEntries = 0;
	bool debug = false;
	while((entry = (char*)gSystem->GetDirEntry(dirp))) {
		TString fileName = entry;
		if(!isRootFile(fileName, "SlicedDetections"))    continue;
		kFiles++;
		if(debug && kFiles>5)
			break;
		cout << "\t" << dirIn + fileName << "\n";

      		double energydeposition = 0;
      		int pedetected = 0;
      		int eventnumber = 0;
      		string *material = NULL;
		TFile *input = new TFile(fileName,"READ");
		TTree *fTree = (TTree*)input->Get("fTree");
     		fTree->SetBranchAddress("energydeposition",&energydeposition);
     		fTree->SetBranchAddress("pedetected",&pedetected);
     		fTree->SetBranchAddress("eventnumber",&eventnumber);
     		fTree->SetBranchAddress("material",&material);

		int entries = fTree->GetEntries();
      		cout <<"\tNumber of entries in current file: "<< entries << endl;
      		int storedeventnumber = -1;

     		fTree->GetEntry(0);
     		storedeventnumber = eventnumber;
     		for (int i=i;i<entries;i++){//For all entries
			fTree->GetEntry(i);
         		if(eventnumber!=storedeventnumber){//Event is complete, begin energy summing of another event
                		if(totalgeenergy>geThreshold){//Total Ge energy deposition has to exceed 200 keV
                        		oTree->Fill();
					kProducedEntries++;
                		}
                		totalpe = 0;
                		totalgeenergy = 0;
                		totalargonenergy = 0;
                		storedeventnumber = eventnumber;
         		}

		        if(strstr(material->c_str(), "ArgonLiquid"))
                		totalargonenergy+=energydeposition;
		        if(strstr(material->c_str(), "GermaniumEnriched"))
                		totalgeenergy+=energydeposition;
		        totalpe += pedetected;
		}
		cout << "\n";
		input->Close();
	}
	cout << "[Info] Writing the results in " << output->GetName() << "\n";
	cout << "[Info] Produced " << kProducedEntries << " entries\n";
	output->cd();
	oTree->Write();
	output->Close();
}

