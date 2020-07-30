/*
 * =====================================================================================
 *
 *       Filename:  check_space_muons.C
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  28/07/2020 09:55:14
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Luigi Berducci (), berducci.1647887@studenti.uniroma1.it
 *        Company:  University of Rome La Sapienza
 *
 * =====================================================================================
 */

bool isRootFile(TString fileName, TString prefix){
    // Consider files `prefixXXXXXX.root`
    TString begin = fileName(0, prefix.Length());
    if(begin.EqualTo(prefix) & fileName.EndsWith("root"))
        return true;
    return false;
}

void check_space_muons(){
    char* fullDirIn = gSystem->ExpandPathName("./");
    void* dirp = gSystem->OpenDirectory(fullDirIn);
    const char* entry;
    Int_t filePartID = 0, kFiles=0;
    bool debug = false;
    // Variables
    double r=0, z=0;
    string * material = 0;
    TFile * oFile = TFile::Open("RZExtraction.root", "RECREATE");
    TTree * oTree = new TTree("fTree", "");
    oTree->Branch("r", &r);
    oTree->Branch("z", &z);
    oTree->Branch("material", &material);
    while((entry = (char*)gSystem->GetDirEntry(dirp))) {
        TString fileName = entry;
        if(!isRootFile(entry, "SlicedDetections"))   continue;
        cout << "\t" << fullDirIn + fileName << "\n";
        kFiles++;
        if(debug && kFiles>5){
            cout << "[Info] Interrupt by DEBUG MODE.\n";
            break;
        }
        TFile *f = TFile::Open(fullDirIn + fileName);
        TTree *fTree = (TTree*) f->Get("fTree");
        fTree->SetBranchAddress("r", &r);
        fTree->SetBranchAddress("z", &z);
        fTree->SetBranchAddress("material", &material);

        for(int i=0; i<fTree->GetEntries(); i++){
            fTree->GetEntry(i);
            oTree->Fill();
        }
    }
    oFile->cd();
    oTree->Write();
    oFile->Close();
    cout << "[Info] Done.\n";
}
