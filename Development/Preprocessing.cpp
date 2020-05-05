#include <iostream>
#include <TEntryList.h>
#include <TParameter.h>
#include <TFile.h>
#include <TString.h>
#include <TTree.h>
#include <TBranch.h>
#include <TChain.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TRandom.h>
#include <array>
#include <vector>
#include <set>
#include <assert.h>

#define PI 3.14159265

#define CJSimFilePrefix "output"
#define TmpROIFilePrefix "tmproi"
#define OutROIFilePrefix "roi"
#define TmpSiPMFilePrefix "tmpsipm"
#define SiPMFilePrefix "sipm"

using namespace std;

pair<TFile*, pair<TH3D*, TH2D*>> getOpticalMap(const char * mapDir){
    TFile *mapFile;
    TH3D *hMap;
    TH2D *oMap;
    char* fullMapDir = gSystem->ExpandPathName(mapDir);
    // Get Neil's map
    TString mapFileName = "OpticalMapL200XeD.14String.5mm";
    mapFile = TFile::Open(fullMapDir + mapFileName + TString(".root"));
    mapFile->GetObject("ProbMapInterior", hMap);
    mapFile->GetObject("ProbMapExterior", oMap);
    pair<TH3D*, TH2D*> map_pair = make_pair(hMap, oMap);
    pair<TFile*, pair<TH3D*, TH2D*>> filemap_pair = make_pair(mapFile, map_pair);
    return filemap_pair;
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
    // Debug
    // cout << "[Cut ROI] File: " << fileName << ", original entries: " << originalTree->GetEntries() << ", reduced entries: " << reducedTree->GetEntries() << endl;
    originalTree->Delete();    // I cannot delete file, otwise seg fault
    file->Close();
    file->Delete();
    return make_pair(output, reducedTree);
}

void addBranches(TTree * reducedTree, TH3D * hMap, TH2D * oMap, Int_t event_x_file, Int_t file_number){
    Double_t x, y, z, detectionefficiency;
    TBranch *bDEff = reducedTree->Branch("detectionefficiency",
                                         &detectionefficiency,
                                         "detectionefficiency/D");
    Int_t eventnumber, increment_eventnumber;
    TBranch *bIncENumber = reducedTree->Branch("inc_eventnumber",
                        	               &increment_eventnumber,
                                	       "inc_eventnumber/I");
    reducedTree->SetBranchAddress("x",&x);
    reducedTree->SetBranchAddress("y",&y);
    reducedTree->SetBranchAddress("z",&z);
    reducedTree->SetBranchAddress("eventnumber", &eventnumber);
    for (Long64_t i = 0; i < reducedTree->GetEntries(); i++) {
        reducedTree->GetEntry(i);
        Double_t r = sqrt(x*x + y*y);    // Euclidean distance ignoring Z
        /* if (r >= 300){      // WAIT: exterior map is bugged */
        /*     Int_t bin = oMap->FindBin(r, z); */
        /*     detectionefficiency = hMap->GetBinContent(bin); */
        /* }else{ */
            /* Int_t bin = hMap->FindBin(x, y, z); */
            /* detectionefficiency = hMap->GetBinContent(bin); */
        /* } */
        Int_t bin = hMap->FindBin(x, y, z);
        detectionefficiency = hMap->GetBinContent(bin);
        increment_eventnumber = file_number * event_x_file + eventnumber;
        bDEff->Fill();
        bIncENumber->Fill();	// for ar39 event number overlaps
    }
}

void data_cleaning(const char * dirIn, const char * dirOut, const char * mapDir){
    cout << "[Info] Data cleaning: Only ROI and Edep>0. Add OP DetEff.\n";
    // IO directories
    char* fullDirIn = gSystem->ExpandPathName(dirIn);
    char* fullDirOut = gSystem->ExpandPathName(dirOut);
    pair<TFile*, pair<TH3D*, TH2D*>> pair_file_map = getOpticalMap(mapDir);
    TFile* mapFile = pair_file_map.first;
    pair<TH3D*, TH2D*> map_pair = pair_file_map.second;
    TH3D* hMap = map_pair.first;
    TH2D* oMap = map_pair.second;

    void* dirp = gSystem->OpenDirectory(fullDirIn);
    const char* entry;
    Int_t event_x_file = 10000, ar39_file_counter = 0;	// to solve the repeatition of event number of ar39 sims
    while((entry = (char*)gSystem->GetDirEntry(dirp))) {
        TString fileName = entry;
        if(!isSimulationFile(fileName))    continue;
        cout << "\t" << dirIn + fileName;

        pair<TFile*, TTree*> file_tree_pair = getReducedTree(fullDirIn, fullDirOut, fileName);
        TFile * output = file_tree_pair.first;
        TTree * reducedTree = file_tree_pair.second;
        addBranches(reducedTree, hMap, oMap, event_x_file, ar39_file_counter);
        
        cout << " -> " << fullDirOut << output->GetName() << endl;
        reducedTree->Write();
        reducedTree->Delete();
        output->Close();
        output->Delete();
	ar39_file_counter++;	// increment file counter to avoid future events' overlaps
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

void compact_data(const char * dirIn, const char * dirOut, TString prefixIn, TString prefixOut, Long64_t entry_x_file=2000000){
    cout << "[Info] Merge data...\n";
    char* fullDirIn = gSystem->ExpandPathName(dirIn);
    char* fullDirOut = gSystem->ExpandPathName(dirOut);
    void* dirp = gSystem->OpenDirectory(fullDirIn);
    const char* entry;
    TChain* ch = new TChain("fTree");
    Long64_t k_entries = 0;
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
        TParameter<Int_t> *NSiPM = new TParameter<Int_t>("NSiPM", nSiPM);
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
            TString branchDesc = branchName + "/L";
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
        NSiPM->Write();
        output->Close();
        output->Delete();
    }
    cout << "[Info] Data preparation: completed.\n";
}

TString createDatasetOutFile(const char * dirOut, TString prefixOut, Double_t nDeltaT, Double_t DeltaTns){
    char* fullDirOut = gSystem->ExpandPathName(dirOut);
    TString outFileName = fullDirOut + prefixOut + "_T";
    outFileName += nDeltaT * DeltaTns;
    outFileName += "_DT";
    outFileName += DeltaTns;
    outFileName += ".csv";
    return outFileName;
}

set<Int_t> getSetOfEventNumbers(TTree *fTree){
    Int_t eventnumber;
    fTree->SetBranchAddress("eventnumber", &eventnumber);
    // Collect distinct event numbers
    set<Int_t> eventnumbers;
    for(Long64_t i = 0; i < fTree->GetEntries(); i++){
        fTree->GetEntry(i);
        eventnumbers.insert(eventnumber);
    }
    return eventnumbers;
}

vector<vector<Long64_t>> newDatasetEventInstance(Int_t nSiPM, Int_t nDeltaT){
    vector<vector<Long64_t>> TSiPMEvent;
    for(int dt = 0; dt < nDeltaT; dt++){
        vector<Long64_t> dTSiPMSnapshot(nSiPM);
        TSiPMEvent.push_back(dTSiPMSnapshot);
    }
    return TSiPMEvent;
}

TEntryList* getEntryListOfEvent(TTree *fTree, Int_t eventnumber){
    // Select entries of event
    TString selection = "eventnumber==";
    selection += eventnumber;
    fTree->Draw(">>entries", selection, "entrylist");
    TEntryList *elist = (TEntryList*)gDirectory->Get("entries");
    return elist;
}

void produce_time_dataset(const char * dirIn, const char * dirOut, TString prefixIn, TString prefixOut){
    Double_t DeltaTns = 2, shift = 0;    // T, dT in ns
    Int_t nDeltaT = 40;     // number of Dt for each event
    cout << "[Info] Create dataset wt T=" << nDeltaT * DeltaTns << ", dT=" << DeltaTns << "...\n";
    char* fullDirIn = gSystem->ExpandPathName(dirIn);
    // Create output file stream
    ofstream outCSV;
    TString outFile(createDatasetOutFile(dirOut, prefixOut, nDeltaT, DeltaTns));
    outCSV.open(outFile);
    // Write CSV header: number of Dt, Dt in ns, resulting time T
    outCSV << nDeltaT << "," << DeltaTns << "," << nDeltaT * DeltaTns << "\n";
    // Loop files in input directory
    void* dirp = gSystem->OpenDirectory(fullDirIn);
    const char* entry;
    while((entry = (char*)gSystem->GetDirEntry(dirp))) {
        TString fileName = entry;
        if(!isRootFile(fileName, prefixIn))   continue;
        cout << "\t" << fullDirIn + fileName << endl;
        // Open file, get tree and number of sipm
        TFile *f = TFile::Open(fullDirIn + fileName);
        TTree *fTree = (TTree*) f->Get("fTree");
        TParameter<Int_t> *NSiPM = (TParameter<Int_t>*) f->Get("NSiPM");
        const int nSiPM = (const int) NSiPM->GetVal();
        // Collect distinct event numbers
        set<Int_t> eventnumbers = getSetOfEventNumbers(fTree);
        // Connect branches
        Int_t eventnumber;
        Double_t time;
        fTree->SetBranchAddress("eventnumber", &eventnumber);
        fTree->SetBranchAddress("time", &time);
        vector<Long64_t> SiPM(nSiPM);
        for(int sipm = 0; sipm < nSiPM; sipm++){
            TString branchName = "SiPM";
            branchName += sipm;
            fTree->SetBranchAddress(branchName, &SiPM[sipm]);
        }
        // Loop over events
        for(auto event : eventnumbers){
            // Create dataset instance struct
            vector<vector<Long64_t>> TSiPMEvent = newDatasetEventInstance(nSiPM, nDeltaT);
            // Loop over event entries
            TEntryList *elist = getEntryListOfEvent(fTree, event);
            Long64_t eventEntry;
            while((eventEntry = elist->Next()) >= 0){   // when list ends -1
                fTree->GetEntry(eventEntry);
                int idDTSnapshot = floor((time + shift) / DeltaTns);
                if(idDTSnapshot < 0 || idDTSnapshot >= nDeltaT)
                    continue;   // all the others are bigger than time T (eventually overflow)
                // Integrate according to Dt
                for(int sipm = 0; sipm < nSiPM; sipm++){
                    TSiPMEvent[idDTSnapshot][sipm] += SiPM[sipm];
                }
            }
            // Write the event to the out csv file
            for(auto snapshot : TSiPMEvent){
                for(auto sipm : snapshot)
                    outCSV << sipm << ",";
                outCSV << endl;
            }
            outCSV << endl;
        }
        // Close file and tree
        fTree->Delete();
        f->Close();
        f->Delete();
    }
    cout << "\tWritten in file " << outFile << endl;
    gSystem->FreeDirectory(dirp);
}

int main(){
    cout << "[Info] Preprocessing...\n";
    // Data cleaning
    const char * dirIn = "/home/data/Ar39/";
    const char * dirOut = "/home/data/Ar39Preproc/";
    const char * mapDir = "/home/data/";
    // Local
    // const char * dirIn = "Data/";
    // const char * dirOut = "Out/";
    // const char * mapDir = "../Data/root/";
    // Data cleaning
    Long64_t entry_x_file = 3000000;	//Compact root files to have this number of entries
    data_cleaning(dirIn, dirOut, mapDir);
    compact_data(dirOut, dirOut, TmpROIFilePrefix, OutROIFilePrefix, entry_x_file);
    // Data preparation
    /* data_preparation(dirOut, dirOut); */
    //produce_time_dataset(dirOut, dirOut, TmpSiPMFilePrefix, "dataset");
    cout << "[Info] End.\n";
}
