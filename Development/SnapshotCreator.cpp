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

map<Int_t, Double_t> getFirstTimeOfEvents(TTree *fTree){
    Int_t eventnumber;
    Double_t time;
    fTree->SetBranchAddress("eventnumber", &eventnumber);
    fTree->SetBranchAddress("time", &time);
    // Collect distinct event numbers
    set<Int_t> eventnumbers;
    map<Int_t, Double_t> map_event_t0;
    for(Long64_t i = 0; i < fTree->GetEntries(); i++){
        fTree->GetEntry(i);
	    map<Int_t, Double_t>::iterator it = map_event_t0.find(eventnumber);
	    if ((it == map_event_t0.end()) || (it->second > time)) {
            map_event_t0.insert(make_pair(eventnumber, time));
        }
    }
    return map_event_t0;
}

map<Int_t, Double_t> getRndOffsetPerEvents(map<Int_t, Double_t> map_event_t0, Double_t max_offset){
    map<Int_t, Double_t> map_event_offset;
    TRandom rnd = TRandom();
    for(auto event_t0 : map_event_t0){
	Double_t offset = rnd.Uniform(max_offset);    // unif. rnd [0, maxoffset]
        map_event_offset[event_t0.first] = offset;
    }
    return map_event_offset;
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
    Double_t DeltaT = 4;    	// integration time (ns) - 4ns is Dt of FlashADC
    Int_t nDeltaT = 2500;     	// number of successive integrations - 2500 integrations of 4ns are 10us
    nDeltaT = 25;           	// number of successive integrations - 2500 integrations of 4ns are 10us
    Double_t margin = 5;	    // Margin at the end to avoid partial events -> approx. event length
    int num_shiftings = 100;	// Number of instances for each event
    num_shiftings = 2;	// Number of instances for each event
    int min_shifting = 0;	// Interval shifting (min)
    int max_shifting = nDeltaT * DeltaT - margin;	// Interval shifting (max)
    int group_events = 2;	// Number of events to be grouped in the same snapshot
    // Debug
    cout << "[Info] From files " << dirIn << "/" << prefixIn << "*" <<  endl;
    cout << "[Info] Create snapshot of T=" << nDeltaT * DeltaT << " ns, ";
    cout << "dT=" << DeltaT << "\n";
    // IO management
    ofstream outCSV;
    TString outFile(createDatasetFilename(dirOut, prefixOut, nDeltaT, DeltaT));
    outCSV.open(outFile);
    cout << "[Info] Writing in " << outFile << "...\n";
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
        cout << "[Debug] Loaded Tree: nentries " << fTree->GetEntries() << endl;

        // Loop over events
        map<Int_t, Double_t> map_event_t0 = getFirstTimeOfEvents(fTree);
        cout << "[Debug] Computed T0: nevents " << map_event_t0.size() << endl;
	    cout << "Firs time event 11 " << map_event_t0[11] << endl;
	    exit(0);
        map<Int_t, Double_t> map_event_offset = getRndOffsetPerEvents(map_event_t0, max_shifting);

        // Connect branches
        Int_t eventnumber;
        Double_t time;
        fTree->SetBranchAddress("eventnumber", &eventnumber);
        fTree->SetBranchAddress("time", &time);
        vector<Long64_t> SiPM(nSiPM);
        for(int sipm = 0; sipm < nSiPM; sipm++){
            TString branchName = "Slice";
            branchName += sipm;
            fTree->SetBranchAddress(branchName, &SiPM[sipm]);
        }

	    Int_t nevents = map_event_t0.size();
	    Int_t ecounter = 0, last_event = -1;
        vector<vector<Long64_t>> TSiPMEvent = newDatasetEventInstance(nSiPM, nDeltaT);
        for(Int_t i=0; i < fTree->GetEntries(); i++){
	    fTree->GetEntry(i);
            cout << "Event: " << eventnumber << ", Time: " << time << endl;
	    if(eventnumber != last_event){	// New Event
	    	last_event = eventnumber;
		ecounter++;
	        if(ecounter % group_events == 1){
		    if(ecounter > 1){
		        // Write output
    		        for(auto snapshot : TSiPMEvent){
			    for(auto sipm : snapshot){
			        outCSV << sipm << ",";
			        outCSV << endl;
			    }
		        }
                        outCSV << endl;
                    }
            	    TSiPMEvent = newDatasetEventInstance(nSiPM, nDeltaT);
		}
	    }
	    Double_t shifted_time = time - map_event_t0[eventnumber] + map_event_offset[eventnumber];
            int id_time_bin = floor((shifted_time) / DeltaT);
            if(id_time_bin < 0 || id_time_bin >= nDeltaT){
	        cout << "[Info] Skipped entry out-of-time. Entry: " << i << ", Event: " << eventnumber << endl; 
                continue;   // all the others are bigger than time T (eventually overflow)
	    }
            // Integrate in the time bin
            for(int sipm = 0; sipm < nSiPM; sipm++){
            	TSiPMEvent[id_time_bin][sipm] += SiPM[sipm];
            }
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
    // const char * dirIn = "/home/data/Ar39Preproc/";
    // const char * dirOut = "/home/data/Ar39Preproc/";
    // Local
    const char * dirIn = "./";
    const char * dirOut = "Out/";
    produce_time_dataset(dirIn, dirOut, "SiPMTrace_Ar39_", "Ar39_Snapshots");
    cout << "[Info] End.\n";
}
