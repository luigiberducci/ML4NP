//
// Created by luigi on 25/05/20.
//

bool isRootFile(TString fileName, TString prefix){
    // Consider files `prefixXXXXXX.root`
    TString begin = fileName(0, prefix.Length());
    if(begin.EqualTo(prefix) & fileName.EndsWith("root"))
        return true;
    return false;
}


void run_inout_dirs(const char * dirIn, const char * dirOut){
    cout << "[Info] Extract Muons Conditions...\n";
    // IO directories
    char* fullDirIn = gSystem->ExpandPathName(dirIn);
    char* fullDirOut = gSystem->ExpandPathName(dirOut);
    TString outputfilename(dirOut);
    outputfilename += "MuonsStartingConditions.root";
    TFile * output = TFile::Open(outputfilename, "RECREATE");
    TTree * fTree = new TTree("fTree", "initial condition of simulated muons");
    Int_t output_eventnumber = 0;
    Double_t x, y, z, px, py, pz, energy;
    fTree->Branch("eventnumber", &output_eventnumber, "eventnumber/I");
    fTree->Branch("x", &x, "x/D");
    fTree->Branch("y", &y, "y/D");
    fTree->Branch("z", &z, "z/D");
    fTree->Branch("energy", &energy, "energy/D");
    fTree->Branch("px", &px, "px/D");
    fTree->Branch("py", &py, "py/D");
    fTree->Branch("pz", &pz, "pz/D");

    void* dirp = gSystem->OpenDirectory(fullDirIn);
    const char* entry;
    while((entry = (char*)gSystem->GetDirEntry(dirp))) {
        TString fileName = entry;
        if(!isRootFile(fileName, "output"))    continue;
        cout << "\t" << dirIn + fileName << endl;

        // loop file entries
        TFile *f = TFile::Open(fullDirIn + fileName);
        TTree *originalTree = (TTree*) f->Get("fTree");
        Int_t previous_eventnumber = -1, eventnumber;
        originalTree->SetBranchAddress("muonx", &x);
        originalTree->SetBranchAddress("muony", &y);
        originalTree->SetBranchAddress("muonz", &z);
        originalTree->SetBranchAddress("muonenergy", &energy);
        originalTree->SetBranchAddress("muonpx", &px);
        originalTree->SetBranchAddress("muonpy", &py);
        originalTree->SetBranchAddress("muonpz", &pz);
        originalTree->SetBranchAddress("eventnumber", &eventnumber);
        for(Int_t i = 0; i < originalTree->GetEntries(); i++){
            originalTree->GetEntry(i);
            if(eventnumber != previous_eventnumber){
                output_eventnumber++;
                previous_eventnumber = eventnumber;
                fTree->Fill();  // Write the starting condition of the muon on each new event
            }
        }
	originalTree->Delete();
	f->Close();
	f->Delete();
    }
    output->cd();
    cout << "[Info] Complete computation. Writing in " << fullDirOut << output->GetName() << endl;
    fTree->Write();
    fTree->Delete();
    cout << "[Info] Tree written and deleted\n";
    output->Close();
    output->Delete();
    cout << "[Info] File closed and deleted\n";
    gSystem->FreeDirectory(dirp);   // Free pointer to dir
}

void extract_initial_muons(){
    run_inout_dirs("./", "./");
}
