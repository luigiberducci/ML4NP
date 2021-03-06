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

// Params
#define RNDSEED 123456789
#define PI 3.14159265
#define INF 999999999

#define ARGONMATERIAL "ArgonLiquid"

// Filename prefixes
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

TString createDatasetFilename(const char * dirOut, TString prefixOut, Double_t nDeltaT, Double_t DeltaTns,
                              Int_t nGroupedEvents, Int_t filePartID){
    char* fullDirOut = gSystem->ExpandPathName(dirOut);
    TString outFileName = fullDirOut + prefixOut + "_T";
    outFileName += nDeltaT * DeltaTns;
    outFileName += "_DT";
    outFileName += DeltaTns;
    outFileName += "_Grp";
    outFileName += nGroupedEvents;
    outFileName += "_Seed";
    outFileName += RNDSEED;
    outFileName += "_Part";
    outFileName += filePartID;
    outFileName += ".csv";
    return outFileName;
}

map<Int_t, Double_t> getFirstDetectionInLArPerEvents(TTree * fTree){
    // Assumption: entries with other material (!=LAr) have NPE==0
    Int_t eventnumber, pedetected;
    Double_t time;
    fTree->SetBranchAddress("eventnumber", &eventnumber);
    fTree->SetBranchAddress("time", &time);
    fTree->SetBranchAddress("pedetected", &pedetected);
    // Collect distinct event numbers
    set<Int_t> eventnumbers;
    map<Int_t, Double_t> map_event_t0;
    for(Long64_t i = 0; i < fTree->GetEntries(); i++){
        fTree->GetEntry(i);
        if(pedetected <= 0)
            continue;
	    map<Int_t, Double_t>::iterator it = map_event_t0.find(eventnumber);
	    if (it == map_event_t0.end()){
		    map_event_t0.insert(make_pair(eventnumber, time));
    	}else{
	    	if (it->second > time)
		    	it->second = time;	//eventually update time with min val
            }
        }
    return map_event_t0;
}

map<Int_t, Double_t> getRndOffsetPerEvents(map<Int_t, Double_t> map_event_t0, Double_t min_offset, Double_t max_offset){
    map<Int_t, Double_t> map_event_offset;
    TRandom rnd = TRandom(RNDSEED);
    for(auto event_t0 : map_event_t0){
	Double_t offset = rnd.Uniform(min_offset, max_offset);    // unif. rnd [minoffset, maxoffset]
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

void write_header(ofstream& outCSV, Int_t nDeltaT, Double_t DeltaT, Int_t filePartID, Int_t nInnerSlices, Int_t nOuterSlices){
    // Write CSV header: number of Dt, Dt in ns, resulting time T
    outCSV << "# File Part: " << filePartID << ", ";
    outCSV << "# N DeltaT: " << nDeltaT << ", ";
    outCSV << "DeltaT: " << DeltaT << ", ";
    outCSV << "T: " << nDeltaT * DeltaT << "\n";
    outCSV << "# Each event is a sequence of NDeltaT rows.";
    outCSV << "The events are separated by '#'.\n";
    outCSV << "# Row Format: ID, InnerSlice0, ..., InnerSlice(N-1), OuterSlice0, ..., OuterSlice(N-1)\n";
    outCSV << "#\n# Columns:\n";
    outCSV << "eventnumber,energydeposition,pedetected,";
    for(Int_t i=0; i < nInnerSlices; i++){
        outCSV << "InnerSlice" << i << ",";
    }
    for(Int_t i=0; i < nOuterSlices; i++){
        outCSV << "OuterSlice" << i << ",";
    }
    outCSV << "\n";
}

template<size_t NDt>
void write_produced_event_in_outfile(ofstream& outCSV, Int_t prodEventID, 
				     vector<vector<Long64_t>> TSiPMEvent_inner, vector<vector<Long64_t>> TSiPMEvent_outer,
                                     vector<Int_t> TSiPMEvent_events, vector<Double_t> TSiPMEvent_offsets,
                                     array<Int_t, NDt> PEDetectedInDt, array<Double_t, NDt> EDepositedInDt){
    // Write output
    outCSV << "# Grouped Events: ";
    for(auto event_id : TSiPMEvent_events)
        outCSV << event_id << ",";
    outCSV << " Event Offsets: ";
    for(auto offset : TSiPMEvent_offsets)
        outCSV << offset << ",";
    outCSV << endl;
    for(int iDt=0; iDt < TSiPMEvent_inner.size(); iDt++){
        auto inner_snapshot = TSiPMEvent_inner[iDt];
        auto outer_snapshot = TSiPMEvent_outer[iDt];
        outCSV << prodEventID << ",";
        outCSV << EDepositedInDt[iDt] << ",";
        outCSV << PEDetectedInDt[iDt] << ",";
        for(auto sipm : inner_snapshot){
            outCSV << sipm << ",";
        }
        for(auto sipm : outer_snapshot){
            outCSV << sipm << ",";
        }
        outCSV << endl;
    }
}

void produce_time_dataset(const char * dirIn, const char * dirOut, TString prefixIn, TString prefixOut, Int_t group_events, Int_t min_npe, Int_t max_npe){
    // Parameters
    const Double_t DeltaT = 10000;    	// integration time (ns) - 4ns is Dt of FlashADC, 10000=10us
    const Int_t nDeltaT = 1;     	// number of successive integrations - e.g. 25 integrations of 4ns are 100ns
    const Double_t margin = 20;         // Margin at the end to avoid partial events -> approx. event length
    const int min_shifting = 0;	        // Interval shifting (min)
    //const int max_shifting = nDeltaT * DeltaT - margin;	// Interval shifting (max)
    const int max_shifting = 0; // No shift -> Trigger on first non-zero detection
    const int max_event_x_file = 250000;  // Max number of output events (snapshots) per file
    // Debug
    cout << "[Info] From files " << dirIn << "/" << prefixIn << "*" <<  endl;
    cout << "[Info] Create snapshot of T=" << nDeltaT * DeltaT << " ns, ";
    cout << "dT=" << DeltaT << "ns, ";
    cout << "grouping " << group_events << " events\n";
    // IO management
    ofstream outCSV;
    Int_t file_part_id = 0;
    Int_t skipped_entries = 0;

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
        TParameter<Int_t> *NInnerSlicesParam= (TParameter<Int_t>*) f->Get("NInnerSlices");
        TParameter<Int_t> *NOuterSlicesParam= (TParameter<Int_t>*) f->Get("NOuterSlices");
        //const int nSiPM = (const int) NSiPM->GetVal();
        const int nInnerSlices = NInnerSlicesParam->GetVal();
        const int nOuterSlices = NOuterSlicesParam->GetVal();
        cout << "[Debug] Loaded Tree: nentries " << fTree->GetEntries() << endl;

        // Compute time of first deposit in lar, and random offset for event-shifting
        map<Int_t, Double_t> map_event_t0 = getFirstDetectionInLArPerEvents(fTree);
        cout << "[Debug] Computed T0 for nevents " << map_event_t0.size() << endl;
        map<Int_t, Double_t> map_event_offset = getRndOffsetPerEvents(map_event_t0, min_shifting, max_shifting);
        cout << "[Debug] Computed Rnd Offsets for nevents " << map_event_offset.size() << endl;

        // Connect branches
        Int_t eventnumber;
        Double_t time, energydeposition;
        fTree->SetBranchAddress("eventnumber", &eventnumber);
        fTree->SetBranchAddress("energydeposition", &energydeposition);
        fTree->SetBranchAddress("time", &time);
        vector<Int_t> innerSlices(nInnerSlices);
        vector<Int_t> outerSlices(nOuterSlices);
        for(int sipm = 0; sipm < nInnerSlices; sipm++){
            TString branchName = "InnerSlice";
            branchName += sipm;
            fTree->SetBranchAddress(branchName, &innerSlices[sipm]);
        }
        for(int sipm = 0; sipm < nOuterSlices; sipm++){
            TString branchName = "OuterSlice";
            branchName += sipm;
            fTree->SetBranchAddress(branchName, &outerSlices[sipm]);
        }

    	Int_t nevents = map_event_t0.size();
	    Int_t ecounter = 0, last_event = -1, prodevents_counter=0, written_events_counter=0;
    	// Create data struct
        vector<vector<Long64_t>> TSiPMEvent_inner = newDatasetEventInstance(nInnerSlices, nDeltaT);
        vector<vector<Long64_t>> TSiPMEvent_outer = newDatasetEventInstance(nOuterSlices, nDeltaT);
        array<Double_t, nDeltaT> EDepositedInDt = {0};     // cumulative Edep in Lar, over each Dt
        array<Int_t, nDeltaT> PEDetectedInDt = {0};        // cumulative Nr of PE detected, over each Dt
        vector<Int_t> TSiPMEvent_events;
        vector<Double_t> TSiPMEvent_offsets;
        Long64_t nEntries = fTree->GetEntries();
        for(Int_t i=0; i < nEntries; i++){
	    fTree->GetEntry(i);

	    if(eventnumber != last_event){	// New Event
	    	    last_event = eventnumber;
                    ecounter++;
	            // Note: we have to check group_event==1 every time
	            if((ecounter % group_events == 1) || (group_events == 1)){
                    if((ecounter > 1) || (group_events == 1)){
                            prodevents_counter++;
                            if (prodevents_counter % max_event_x_file == 1) {
                                outCSV.close();
                                file_part_id++; // Increment the id of file part (part1, part2, ..)
                                TString outFile(createDatasetFilename(dirOut, prefixOut, nDeltaT, DeltaT, group_events, file_part_id));
                                outCSV.open(outFile);
                                write_header(outCSV, nDeltaT, DeltaT, file_part_id, nInnerSlices, nOuterSlices);
                                cout << "[Info] Event " << prodevents_counter << " - Writing in " << outFile << "...\n";
                            }

                            // Check filter on Tot NPE
                            Int_t tot_npe = 0;
                            for(auto pe : PEDetectedInDt)
                                tot_npe += pe;
                            if(tot_npe >= min_npe && tot_npe <= max_npe){
                                written_events_counter++;
                                write_produced_event_in_outfile(outCSV, written_events_counter, TSiPMEvent_inner, TSiPMEvent_outer, 
								TSiPMEvent_events, TSiPMEvent_offsets, PEDetectedInDt, EDepositedInDt);
                            }
                    }
			        TSiPMEvent_inner = newDatasetEventInstance(nInnerSlices, nDeltaT);
			        TSiPMEvent_outer = newDatasetEventInstance(nOuterSlices, nDeltaT);
                    TSiPMEvent_events.clear();
                    TSiPMEvent_offsets.clear();
                    PEDetectedInDt.fill(0);
                    EDepositedInDt.fill(0);
		        }
                TSiPMEvent_events.push_back(eventnumber);
                TSiPMEvent_offsets.push_back(map_event_offset[eventnumber]);
	    }
	    // Debug
            if(i % 10000 == 0)
                cout << "\rentry: " << i << "/" << nEntries << std::flush;
	    Double_t shifted_time = time - map_event_t0[eventnumber] + map_event_offset[eventnumber];
            int id_time_bin = floor((shifted_time) / DeltaT);

            if(id_time_bin < 0 || id_time_bin >= nDeltaT){
	            //cout << "[Info] Skipped entry out-of-time. Entry: " << i << ", Event: " << eventnumber << endl;
                skipped_entries++;
                continue;   // all the others are bigger than time T (eventually overflow)
	    }
            // Integrate in the time bin
	    for(int sipm = 0; sipm < nInnerSlices; sipm++){
            	TSiPMEvent_inner[id_time_bin][sipm] += innerSlices[sipm];
                PEDetectedInDt[id_time_bin] += innerSlices[sipm];
            }
	    for(int sipm = 0; sipm < nOuterSlices; sipm++){
            	TSiPMEvent_outer[id_time_bin][sipm] += outerSlices[sipm];
                PEDetectedInDt[id_time_bin] += outerSlices[sipm];
            }
            EDepositedInDt[id_time_bin] += energydeposition;
        }
        // Close file and tree
        fTree->Delete();
        f->Close();
        f->Delete();
        cout << "[Info] Produced snapshots: " << written_events_counter << ", ";
	cout << "skipped entries for time: " << skipped_entries << "\n";
    }
    gSystem->FreeDirectory(dirp);
}

int main(){
    cout << "[Info] Snapshot Dataset Creator...\n";
    // LGND Docker filesystem
    // const char * dirIn = "/home/data/Ar39Preproc/";
    // const char * dirOut = "/home/data/Ar39Preproc/";
    // Local
    //const char * dirIn = "../Data/ar39/06-14-2020-10M/";
    //const char * dirOut = "Out/T10us/Ar39_1to7Pileups/";
    //const int group_events = 1;	        // Number of events to be grouped in the same snapshot
    // User input
    TString dirIn, dirOut, outputPrefix;
    Int_t groupEvents, min_npe, max_npe;
    cout << "[Input] What is the INPUT directory? (where are the files `SlicedDetections*root`)" << endl;
    cin >> dirIn;
    cout << "[Input] What is the OUTPUT directory?" << endl;
    cin >> dirOut;
    cout << "[Input] What is the output PREFIX?" << endl;
    cin >> outputPrefix;
    cout << "[Input] How many events to GROUP for each snapshot?" << endl;
    cin >> groupEvents;
    cout << "[Input] Filter event by NPE: Min NPE?" << endl;
    cin >> min_npe;
    cout << "[Input] Filter event by NPE: Max NPE? (-1 to unbounded)" << endl;
    cin >> max_npe;
    if(max_npe<0)
	    max_npe = INF;
    produce_time_dataset(dirIn, dirOut, "SlicedDetections", outputPrefix, groupEvents, min_npe, max_npe);
    cout << "[Info] End.\n";
}
