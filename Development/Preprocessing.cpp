#include <iostream>
#include <TFile.h>
#include <TString.h>
#include <TTree.h>
#include <TBranch.h>
#include <TChain.h>
#include <TH3D.h>
#include <TROOT.h>
#include <TSystem.h>

using namespace std;

void data_cleaning(const char * dirIn, const char * dirOut, const char * mapDir){
    cout << "[Info] Data cleaning...\n";
    // IO directories
    char* fullDirIn = gSystem->ExpandPathName(dirIn);
    char* fullDirOut = gSystem->ExpandPathName(dirOut);
    char* fullMapDir = gSystem->ExpandPathName(mapDir);
    // Loop files in input dir
    void* dirp = gSystem->OpenDirectory(fullDirIn);
    const char* entry;
    // Get Neil's map
    TString mapFileName = "OpticalMapL200XeD.14String.5mm";
    TFile* mapFile = TFile::Open(fullMapDir + mapFileName + TString(".root"));
    TH3D *hMap;
    mapFile->GetObject("ProbMapInterior", hMap);

    TTree *mergedTree;
    TList *listTree = new TList;
    Long64_t k_entries = 0, entry_x_file = 1000000;
    Int_t id_part = 0;
    while((entry = (char*)gSystem->GetDirEntry(dirp))) {
        // Consider files `outputXXXXXX.root`
        TString fileName = entry;
        TString begin = fileName(0, 6);
        if(!(begin.EqualTo("output") & fileName.EndsWith("root")))
            continue;
        std::cout << "\t[Info] Processing " << dirIn + fileName << endl;
        // Open input file and take original tree
        TFile* file = TFile::Open(dirIn + fileName);
        TTree* originalTree = (TTree*)file->Get("fTree");
        TFile* output = TFile::Open(dirOut + fileName.Copy().Replace(0, 6, "cutROI"), "RECREATE");
        // cut: ROI & Edep>0
        // ROI def: x, y in [-500, +500], z in [-1000, +1000]
        TString cutX("(x >= -500.) & (x <= 500.)");
        TString cutY("(y >= -500.) & (y <= 500.)");
        TString cutZ("(z >= -1000) & (z <= 1000)");
        TString cutE("(energydeposition > 0)");
        // Copy tree
        TTree* reducedTree = originalTree->CopyTree(cutX + " & " + cutY + " & " + cutZ);

        // add branch detection efficiency from Neil's map
        Double_t x, y, z, detectionefficiency;
        TBranch *bDEff = reducedTree->Branch("detectionefficiency",
                                              &detectionefficiency,
                                              "detectionefficiency/D");
        reducedTree->SetBranchAddress("x",&x);
        reducedTree->SetBranchAddress("y",&y);
        reducedTree->SetBranchAddress("z",&z);
        for (Long64_t i = 0; i < reducedTree->GetEntries(); i++) {
            reducedTree->GetEntry(i);
            Int_t bin = hMap->FindBin(x, y, z);
            detectionefficiency = hMap->GetBinContent(bin);
            if(detectionefficiency>1){
                cout << "i: " << i << ", ";
                cout << "x: " << x << ", ";
                cout << "y: " << y << ", ";
                cout << "z: " << z << ", ";
                cout << "deff: " << detectionefficiency << endl;
            }
            bDEff->Fill();
        }
        // save cleaned files
        cout << "\t[Info] original entries: " << originalTree->GetEntries() << ", ";
        cout << "ROI entries: " << reducedTree->GetEntries() << ", ";
        cout << "perc cut: " << (1 - ((double)reducedTree->GetEntries() / originalTree->GetEntries())) << " %\n\n";

        // merge files based on number of entries
        listTree->Add(reducedTree);

    }
    gSystem->FreeDirectory(dirp);   // Free pointer to dir
    cout << "[Info] Data cleaning: completed.\n";
}

void data_preparation(){
    cout << "[Info] Data preparation...\n";
    // generate ntuple with SiPM readouts
    cout << "[Info] Data preparation: completed.\n";
}

int main(){
    cout << "[Info] Preprocessing...\n";
    // Data cleaning
    const char * dirIn = "/home/luigi/Development/ML4NP/Development/Data/";
    const char * dirOut = "/home/luigi/Development/ML4NP/Development/Out/";
    const char * mapDir = "/home/luigi/Development/ML4NP/Data/root/";
    data_cleaning(dirIn, dirOut, mapDir);
    // Data preparation
    data_preparation();
    cout << "[Info] End.\n";
}