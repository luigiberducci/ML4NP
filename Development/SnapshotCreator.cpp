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
    // Parameters
    Double_t DeltaT = 2;    // integration time
    Int_t nDeltaT = 40;     // number of successive integrations
    int num_shiftings = 100;	// Number of instances for each event
    int min_shifting = 0;	// Interval shifting (min)
    int max_shifting = 100;	// Interval shifting (max)
    // Let's go
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
    cout << "[Info] Snapshot Dataset Creator...\n";
    // Data cleaning
    const char * dirIn = "/home/data/Ar39/";
    const char * dirOut = "/home/data/Ar39Preproc/";
    // Local
    // const char * dirIn = "Data/";
    // const char * dirOut = "Out/";
    // Data cleaning
    Long64_t entry_x_file = 3000000;	//Compact root files to have this number of entries
    //produce_time_dataset(dirOut, dirOut, TmpSiPMFilePrefix, "dataset");
    cout << "[Info] End.\n";
}
