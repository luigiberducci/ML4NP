#include <iostream>
#include <TFile.h>
#include <TString.h>
#include <TTree.h>
#include <TBranch.h>
#include <TChain.h>
#include <TH3D.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TRandom.h>
#include <array>
#include <assert.h>

#define PI 3.14159265

#define CJSimFilePrefix "output"
#define TmpROIFilePrefix "tmproi"
#define OutROIFilePrefix "roi"
#define TmpSiPMFilePrefix "tmpsipm"
#define SiPMFilePrefix "sipm"

using namespace std;

pair<TFile*, TH3D*> getOpticalMap(const char * mapDir){
    TFile *mapFile;
    TH3D *hMap;
    char* fullMapDir = gSystem->ExpandPathName(mapDir);
    // Get Neil's map
    TString mapFileName = "OpticalMapL200XeD.14String.5mm";
    mapFile = TFile::Open(fullMapDir + mapFileName + TString(".root"));
    mapFile->GetObject("ProbMapInterior", hMap);
    return make_pair(mapFile, hMap);
}

bool isRootFile(TString fileName, TString prefix){
    // Consider files `prefixXXXXXX.root`
    TString begin = fileName(0, prefix.Length());
    if(begin.EqualTo(prefix) & fileName.EndsWith("root"))
        return true;
    return false;
}

bool isSimulationFile(TString fileName){
    // Consider files `outputXXXXXX.root`
    TString prefix = CJSimFilePrefix;
    return isRootFile(fileName, prefix);
}

bool isTmpROIFile(TString fileName){
    // Consider files `tmproiXXXXXX.root`
    TString prefix = TmpROIFilePrefix;
    return isRootFile(fileName, prefix);
}

bool isOutFile(TString fileName){
    // Consider files `roiXXXXXX.root`
    TString prefix = OutROIFilePrefix;
    return isRootFile(fileName, prefix);
}

pair<TFile*, TTree*> getReducedTree(const char * dirIn, const char * dirOut, TString fileName){
    TFile *file, *output;
    TTree *originalTree, *reducedTree;
    char* fullDirIn = gSystem->ExpandPathName(dirIn);
    // Open input file and take original tree
    file = TFile::Open(dirIn + fileName);
    output = TFile::Open(dirOut + fileName.Copy().Replace(0, 6, TmpROIFilePrefix), "RECREATE");
    originalTree = (TTree*)file->Get("fTree");
    // cut: ROI & Edep>0
    TString cutX("(x >= -500.) & (x <= 500.)");
    TString cutY("(y >= -500.) & (y <= 500.)");
    TString cutZ("(z >= -1000) & (z <= 1000)");
    TString cutE("(energydeposition > 0)");
    // cut original tree
    reducedTree = originalTree->CopyTree(cutX + " & " + cutY + " & " + cutZ);
    originalTree->Delete();    // I cannot delete file, otwise seg fault
    file->Close();
    file->Delete();
    return make_pair(output, reducedTree);
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
    pair<TFile*, TH3D*> pair_file_hmap = getOpticalMap(mapDir);
    TFile* mapFile = pair_file_hmap.first;
    TH3D* hMap = pair_file_hmap.second;

    void* dirp = gSystem->OpenDirectory(fullDirIn);
    const char* entry;
    while((entry = (char*)gSystem->GetDirEntry(dirp))) {
        TString fileName = entry;
        if(!isSimulationFile(fileName))    continue;
        cout << "\t" << dirIn + fileName;

        pair<TFile*, TTree*> file_tree_pair = getReducedTree(fullDirIn, fullDirOut, fileName);
        TFile * output = file_tree_pair.first;
        TTree * reducedTree = file_tree_pair.second;
        addDetEfficiencyBranch(reducedTree, hMap);
        
        cout << " -> " << fullDirOut << output->GetName() << endl;
        reducedTree->Write();
        reducedTree->Delete();
        output->Close();
        output->Delete();
    }
    hMap->Delete();
    mapFile->Close();
    mapFile->Delete();
    gSystem->FreeDirectory(dirp);   // Free pointer to dir
    cout << "[Info] Data cleaning: completed.\n";
}

void mergeChain(TChain* ch, TString dirOut, TString prefix, int id_group){
    TString fileOut(dirOut);
    fileOut += prefix;
    fileOut += "_part";
    fileOut += id_group;
    fileOut += ".root";
    ch->Merge(fileOut);
    // Debug
    cout << "\tMerged " << ch->GetEntries() << " entries in " << fileOut << endl << endl;
}

void compact_data(const char * dirIn, const char * dirOut, TString prefixIn, TString prefixOut){
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
        if(!isRootFile(fileName, prefixIn))    continue;
        cout << "\t" << fullDirIn + fileName << endl;

        if(k_entries >= entry_x_file){
            mergeChain(ch, fullDirOut, prefixOut, ++id_group);
            k_entries = 0;     // reset counter
            ch->Reset();
        }
        ch->Add(fullDirIn + fileName);
        k_entries = ch->GetEntries();
    }
    if(k_entries > 0){
        mergeChain(ch, fullDirOut, prefixOut, ++id_group);
        ch->Reset();
    }
    delete ch;
    gSystem->FreeDirectory(dirp);   // Free pointer to dir
    cout << "[Info] Merge data: completed.\n";
}

void data_preparation(const char * dirIn, const char * dirOut){
    cout << "[Info] Data preparation...\n";
    double m=0, s=5;        // tuned for 72 sipms
    const int opYield=40;   // 40 ops / KeV
    const int nSiPM = 72;
    const double thetaSegment = 2 * PI / nSiPM;
    char* fullDirIn = gSystem->ExpandPathName(dirIn);
    char* fullDirOut = gSystem->ExpandPathName(dirOut);
    void* dirp = gSystem->OpenDirectory(fullDirIn);
    const char* entry;
    while((entry = (char*)gSystem->GetDirEntry(dirp))) {
        TString fileName = entry;
        if(!isOutFile(entry))   continue;
        cout << "\t" << fullDirIn + fileName;

        TFile *f = TFile::Open(fullDirIn + fileName);
        TTree *fTree = (TTree*) f->Get("fTree");

        Double_t x, y, z, time, Edep, deteff;
        Int_t eventnumber;
        array<Long64_t, nSiPM> readouts;
        fTree->SetBranchAddress("x", &x);
        fTree->SetBranchAddress("y", &y);
        fTree->SetBranchAddress("z", &z);
        fTree->SetBranchAddress("time", &time);
        fTree->SetBranchAddress("eventnumber", &eventnumber);
        fTree->SetBranchAddress("energydeposition", &Edep);
        fTree->SetBranchAddress("detectionefficiency", &deteff);

        TString SiPMFilePrefixString = TmpSiPMFilePrefix;
        SiPMFilePrefixString += nSiPM;
        SiPMFilePrefixString += "SiPMs_";
        SiPMFilePrefixString += opYield;
        SiPMFilePrefixString += "Yield";
        TFile *output = TFile::Open(dirOut + fileName.Copy().Replace(0, 3, SiPMFilePrefixString), "RECREATE");
        TTree SiPMTree("fTree", "");
        TBranch *bN = SiPMTree.Branch("eventnumber", &eventnumber, "eventnumber/I");
        TBranch *bT = SiPMTree.Branch("time", &time, "time/D");
        TBranch *bX = SiPMTree.Branch("x", &x, "x/D");
        TBranch *bY = SiPMTree.Branch("y", &y, "y/D");
        TBranch *bZ = SiPMTree.Branch("z", &z, "z/D");
        TBranch *bE = SiPMTree.Branch("energydeposition", &Edep, "energydeposition/D");
        TBranch *bDE = SiPMTree.Branch("detectionefficiency", &deteff, "detectionefficiency/D");
        array<TBranch*, nSiPM> branchSiPM;
        for(int sipm = 0; sipm < nSiPM; sipm++){
            TString branchName = "SiPM";
            branchName += sipm;
            TString branchDesc = branchName + "/I";
            branchSiPM[sipm] = SiPMTree.Branch(branchName, &readouts[sipm], branchDesc);
        }

        for(Long64_t i = 0; i < fTree->GetEntries(); i++){
            fTree->GetEntry(i);
            double angle = atan2(y, x);
            if(angle < 0)
                angle += 2 * PI;
            int segment = angle / thetaSegment;
            assert((segment >= 0) & (segment <nSiPM));
            TRandom rnd = TRandom();
            for(int sipm = 0; sipm < nSiPM; sipm++){
                readouts[sipm] = 0;
            }
            Long64_t opDetected = ceil(Edep * opYield * deteff);
            for(Long64_t op = 0; op < opDetected; op++){
                int r = round(rnd.Gaus(m, s));
                if (segment + r < 0)
                    readouts[segment + r + nSiPM]++;
                else
                    readouts[segment + r]++;
            }
            if(opDetected > 0)
                SiPMTree.Fill();
        }

        cout << " -> " << fullDirOut << output->GetName() << endl;
        SiPMTree.Write();
        output->Close();
        output->Delete();
    }
    cout << "[Info] Data preparation: completed.\n";
}

int main(){
    cout << "[Info] Preprocessing...\n";
    // Data cleaning
    /* const char * dirIn = "/home/data/"; */
    /* const char * dirOut = "/home/data/ROI/"; */
    /* const char * mapDir = "/home/data/"; */
    // Local
    const char * dirIn = "Data/";
    const char * dirOut = "Out/";
    const char * mapDir = "../Data/root/";
    // Data cleaning
    data_cleaning(dirIn, dirOut, mapDir);
    compact_data(dirOut, dirOut, TmpROIFilePrefix, OutROIFilePrefix);
    // Data preparation
    data_preparation(dirOut, dirOut);
    compact_data(dirOut, dirOut, TmpSiPMFilePrefix, SiPMFilePrefix);
    cout << "[Info] End.\n";
}
