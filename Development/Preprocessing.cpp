#include <iostream>
#include <TFile.h>
#include <TString.h>
#include <TTree.h>
#include <TBranch.h>
#include <TChain.h>
#include <TH3D.h>
#include <TROOT.h>
#include <TSystem.h>

#define TmpROIFilePrefix "tmproi"
#define OutFilePrefix "roi"
#define CJSimFilePrefix "output"

using namespace std;

TH3D* getOpticalMap(const char * mapDir){
    char* fullMapDir = gSystem->ExpandPathName(mapDir);
    // Get Neil's map
    TString mapFileName = "OpticalMapL200XeD.14String.5mm";
    TFile* mapFile = TFile::Open(fullMapDir + mapFileName + TString(".root"));
    TH3D *hMap;
    mapFile->GetObject("ProbMapInterior", hMap);
    return hMap;
}

bool isSimulationFile(TString fileName){
    // Consider files `outputXXXXXX.root`
    TString begin = fileName(0, 6);
    if(begin.EqualTo(CJSimFilePrefix) & fileName.EndsWith("root"))
        return true;
    return false;
}

bool isTmpROIFile(TString fileName){
    // Consider files `tmproiXXXXXX.root`
    TString begin = fileName(0, 6);
    if(begin.EqualTo(TmpROIFilePrefix) & fileName.EndsWith("root"))
        return true;
    return false;
}

TTree* getReducedTree(const char * dirIn, const char * dirOut, TString fileName){
    char* fullDirIn = gSystem->ExpandPathName(dirIn);
    // Open input file and take original tree
    TFile* file = TFile::Open(dirIn + fileName);
    TFile* output = TFile::Open(dirOut + fileName.Copy().Replace(0, 6, TmpROIFilePrefix), "RECREATE");
    TTree* originalTree = (TTree*)file->Get("fTree");
    // cut: ROI & Edep>0
    TString cutX("(x >= -500.) & (x <= 500.)");
    TString cutY("(y >= -500.) & (y <= 500.)");
    TString cutZ("(z >= -1000) & (z <= 1000)");
    TString cutE("(energydeposition > 0)");
    TTree* reducedTree = originalTree->CopyTree(cutX + " & " + cutY + " & " + cutZ);
    reducedTree->SetDirectory(output);      // otherwise remain related to infile
    delete originalTree;    // I cannot delete file, otwise seg fault
    return reducedTree;
}

void addDetEfficiencyBranch(TTree * reducedTree, TH3D * hMap){
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
        bDEff->Fill();
    }
}

void data_cleaning(const char * dirIn, const char * dirOut, const char * mapDir){
    cout << "[Info] Data cleaning: Only ROI and Edep>0. Add OP DetEff.\n";
    // IO directories
    char* fullDirIn = gSystem->ExpandPathName(dirIn);
    char* fullDirOut = gSystem->ExpandPathName(dirOut);
    TH3D* hMap = getOpticalMap(mapDir);

    void* dirp = gSystem->OpenDirectory(fullDirIn);
    const char* entry;
    while((entry = (char*)gSystem->GetDirEntry(dirp))) {
        TString fileName = entry;
        if(!isSimulationFile(fileName))    continue;
        cout << "\t" << dirIn + fileName << endl;

        TFile* output = TFile::Open(dirOut + fileName.Copy().Replace(0, 6, TmpROIFilePrefix), "RECREATE");
        TTree* reducedTree = getReducedTree(fullDirIn, fullDirOut, fileName);
        addDetEfficiencyBranch(reducedTree, hMap);
        reducedTree->Write();
        output->Close();
        delete reducedTree;
        delete output;
    }
    delete hMap;
    gSystem->FreeDirectory(dirp);   // Free pointer to dir
    cout << "[Info] Data cleaning: completed.\n";
}

void mergeChain(TChain* ch, const char * dirOut, int id_group){
    TString fileOut(dirOut);
    fileOut += OutFilePrefix;
    fileOut += "_part";
    fileOut += id_group;
    fileOut += ".root";
    ch->Merge(fileOut);
    // Debug
    cout << "\tMerged " << ch->GetEntries() << " entries in " << fileOut << endl << endl;
}

void compact_data(const char * dirIn, const char * dirOut){
    cout << "[Info] Merge data...\n";
    char* fullDirIn = gSystem->ExpandPathName(dirIn);
    char* fullDirOut = gSystem->ExpandPathName(dirOut);
    void* dirp = gSystem->OpenDirectory(fullDirIn);
    const char* entry;
    TChain* ch = new TChain("fTree");
    Long64_t k_entries = 0, entry_x_file = 2000000;
    Int_t id_group = 0;
    while((entry = (char*)gSystem->GetDirEntry(dirp))) {
        TString fileName = entry;
        if(!isTmpROIFile(fileName))    continue;
        cout << "\t" << dirIn + fileName << endl;

        if(k_entries >= entry_x_file){
            mergeChain(ch, fullDirOut, ++id_group);
            k_entries = 0;     // reset counter
            ch->Reset();
        }
        ch->Add(dirIn + fileName);
        k_entries = ch->GetEntries();
    }
    if(k_entries > 0){
        mergeChain(ch, fullDirOut, ++id_group);
        ch->Reset();
    }
    delete ch;
    gSystem->FreeDirectory(dirp);   // Free pointer to dir
    cout << "[Info] Merge data: completed.\n";
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
    compact_data(dirOut, dirOut);
    // Data preparation
    data_preparation();
    cout << "[Info] End.\n";
}