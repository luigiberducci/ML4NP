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
#define INF 999999999999999999999999999999999999999999999

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

TString createDatasetFilename(const char * dirOut, TString prefixOut, Double_t nDeltaT, Double_t DeltaTns){
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
    Double_t DeltaT = 10000;    // integration time (ns)
    Int_t nDeltaT = 5;     	// number of successive integrations
    int num_shiftings = 100;	// Number of instances for each event
    int min_shifting = 0;	// Interval shifting (min)
    int max_shifting = 100;	// Interval shifting (max)
    // Debug
    cout << "[Info] From files " << dirIn << "/" << prefixIn << "*" <<  endl;
    cout << "[Info] Create snapshot of T=" << nDeltaT * DeltaT << " ns, ";
    cout << "dT=" << DeltaT << "\n";
    // IO management
    ofstream outCSV;
    TString outFile(createDatasetFilename(dirOut, prefixOut, nDeltaT, DeltaT));
    outCSV.open(outFile);
    // Write CSV header: number of Dt, Dt in ns, resulting time T
    outCSV << nDeltaT << "," << DeltaT << "," << nDeltaT * DeltaT << "\n";

    // Loop files in input directory
    char* fullDirIn = gSystem->ExpandPathName(dirIn);
    void* dirp = gSystem->OpenDirectory(fullDirIn);
    const char* entry;
    while((entry = (char*)gSystem->GetDirEntry(dirp))) {
        TString fileName = entry;
        if(!isRootFile(fileName, prefixIn))   continue;
        cout << "\t" << fullDirIn + fileName << endl;	// Debug
        // Open file, get tree and number of sipm
        TFile *f = TFile::Open(fullDirIn + fileName);
        TTree *fTree = (TTree*) f->Get("fTree");
        TParameter<Int_t> *NSiPM = (TParameter<Int_t>*) f->Get("NSiPM");
        const int nSiPM = (const int) NSiPM->GetVal();

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
        set<Int_t> eventnumbers = getSetOfEventNumbers(fTree);
	Int_t nevents = eventnumbers.size();
	Int_t ecounter = 0;
        for(auto event : eventnumbers){
	    ecounter++;
	    cout << "\t\rat:" << ecounter << "/" << nevents << std::flush;	// Debug
	    // TODO introduce shifting here!
            // Create empty snapshot struct
            vector<vector<Long64_t>> TSiPMEvent = newDatasetEventInstance(nSiPM, nDeltaT);
            // Loop over event entries
            TEntryList *elist = getEntryListOfEvent(fTree, event);
            Long64_t eventEntry;
	    // Compute time of first deposit (time0)
	    Long64_t time0 = INF;
            while((eventEntry = elist->Next()) >= 0){   // when list ends -1
		fTree->GetEntry(eventEntry);
		if (time < time0){
		    time0 = time;
		}
	    }
	    // Loop entries
	    elist = getEntryListOfEvent(fTree, event);
            while((eventEntry = elist->Next()) >= 0){   // when list ends -1
                fTree->GetEntry(eventEntry);
                int id_time_bin = floor((time - time0 + min_shifting) / DeltaT);
                if(id_time_bin < 0 || id_time_bin >= nDeltaT)
                    continue;   // all the others are bigger than time T (eventually overflow)
                // Integrate in the time bin
                for(int sipm = 0; sipm < nSiPM; sipm++){
                    TSiPMEvent[id_time_bin][sipm] += SiPM[sipm];
                }
            }
            // Write the event to the out csv file
	    // TODO: add check empty snapshot (all zero non write it)
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
    // LGND Docker filesystem
    // const char * dirIn = "/home/data/Ar39/";
    // const char * dirOut = "/home/data/Ar39Preproc/";
    // Local
    const char * dirIn = "./";
    const char * dirOut = "Out/";
    produce_time_dataset(dirIn, dirOut, "SiPMTrace_Ar39_", "Ar39_Snapshots_");
    cout << "[Info] End.\n";
}
