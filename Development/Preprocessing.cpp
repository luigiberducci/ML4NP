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

#define RNDSEED 123456789
#define PI 3.14159265
#define EVENTxFILE 100000

#define CJSimFilePrefix "output"
#define TmpROIFilePrefix "tmproi"
#define OutROIFilePrefix "roi"
#define TmpSiPMFilePrefix "SlicedDeposits"

using namespace std;

pair<TFile*, pair<TH3D*, TH2D*>> getOpticalMap(const char * mapDir){
    TFile *mapFile;
    TH3D *hMap;
    TH2D *oMap;
    char* fullMapDir = gSystem->ExpandPathName(mapDir);
    // Get Neil's map
    // TString mapFileName = "OpticalMapL200XeD.14String.5mm";
    // TString intMapName = "ProbMapInterior";
    // TString extMapName = "ProbMapExterior";	// this is bugged
    // Get Old name
    TString mapFileName = "LGND200_14_OpticalMapRealisticDetectorLayout";
    TString intMapName = "OpticalMap_Scaled";
    TString extMapName = "OpticalMap_Scaled";	// same map
    // Open and load maps
    mapFile = TFile::Open(fullMapDir + mapFileName + TString(".root"));
    mapFile->GetObject(intMapName, hMap);
    mapFile->GetObject(extMapName, oMap);
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
    TTree *originalTree, *reducedTree, *checkTree;
    char* fullDirIn = gSystem->ExpandPathName(dirIn);
    // Open input file and take original tree
    file = TFile::Open(dirIn + fileName);
    output = TFile::Open(dirOut + fileName.Copy().Replace(0, 6, TmpROIFilePrefix), "RECREATE");
    originalTree = (TTree*)file->Get("fTree");
    // Check Event consistency
    TString cutEvent("eventnumber > ");  // If this cut has some entries -> BREAK!
    cutEvent += EVENTxFILE;
    checkTree = originalTree->CopyTree(cutEvent);
    assert(checkTree->GetEntries() == 0); // Assert: risk event overlap
    // cut: ROI
    TString cutX("(x >= -700.) & (x <= 700.)");
    TString cutY("(y >= -700.) & (y <= 700.)");
    TString cutZ("(z >= -845.) & (z <= 845.)");
    // cut: edep>0 (deprecated)
    //TString cutE("(energydeposition > 0)");
    // cut original tree
    reducedTree = originalTree->CopyTree(cutX + " & " + cutY + " & " + cutZ);
    // Debug
    cout << " | Cut ROI: original entries: " << originalTree->GetEntries() << ", reduced entries: " << reducedTree->GetEntries();
    originalTree->Delete();    // I cannot delete output file, otwise seg fault
    file->Close();
    file->Delete();
    checkTree->Delete();
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
	// Use interior/exterior map, based on coordinate
        if (r >= 300){
            Int_t bin = hMap->FindBin(x, y, z);
            detectionefficiency = hMap->GetBinContent(bin);
        }else{
            Int_t bin = hMap->FindBin(x, y, z);
            detectionefficiency = hMap->GetBinContent(bin);
        }
        increment_eventnumber = file_number * event_x_file + eventnumber;
        bDEff->Fill();
        bIncENumber->Fill();
    }
}

void data_cleaning(const char * dirIn, const char * dirOut, const char * mapDir){
    cout << "[Info] Data cleaning: Only ROI. Add OP DetEff.\n";
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
    Int_t event_x_file = EVENTxFILE, input_file_counter = 0;	// to solve the repeatition of event number in input files
    while((entry = (char*)gSystem->GetDirEntry(dirp))) {
        TString fileName = entry;
        if(!isSimulationFile(fileName))    continue;
        cout << "\t" << dirIn + fileName;

        pair<TFile*, TTree*> file_tree_pair = getReducedTree(fullDirIn, fullDirOut, fileName);
        TFile * output = file_tree_pair.first;
        TTree * reducedTree = file_tree_pair.second;
        addBranches(reducedTree, hMap, oMap, event_x_file, input_file_counter);
        
        cout << " -> " << fullDirOut << output->GetName() << endl;
        reducedTree->Write();
        reducedTree->Delete();
        output->Close();
        output->Delete();
	input_file_counter++;	// increment file counter to avoid future events' overlaps
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
    fileOut += "_RndSeed";
    fileOut += RNDSEED;
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

TString createDatasetFilename(const char * dirOut, TString prefixOut, Int_t nSlices, Double_t opYield, Double_t quantumEff, Int_t filePartID){
    char* fullDirOut = gSystem->ExpandPathName(dirOut);
    TString outFileName = fullDirOut + prefixOut + "_Slices";
    outFileName += nSlices;
    outFileName += "_Yield";
    outFileName += opYield;
    outFileName += "_QuantumEff";
    outFileName += quantumEff;
    outFileName += "_Seed";
    outFileName += RNDSEED;
    outFileName += "_Part";
    outFileName += filePartID;
    outFileName += ".root";
    return outFileName;
}

void convert_to_sliced_detections(const char * dirIn, const char * dirOut){
    cout << "[Info] Data preparation...\n";
    Double_t m=0, s=5;        // tuned for 72 sipms
    const Int_t opYield=40;   // 40 ops / KeV
    const Double_t quantumEff = .40;   // 40% quantum efficiency
    const Int_t nSiPM = 72;
    const Double_t thetaSegment = 2 * PI / nSiPM;
    char* fullDirIn = gSystem->ExpandPathName(dirIn);
    char* fullDirOut = gSystem->ExpandPathName(dirOut);
    void* dirp = gSystem->OpenDirectory(fullDirIn);
    const char* entry;
    set<Int_t> original_events;
    set<Int_t> produced_events;
    Int_t filePartID = 0;
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
        fTree->SetBranchAddress("inc_eventnumber", &eventnumber);    // To avoid overlap of events
        fTree->SetBranchAddress("energydeposition", &Edep);
        fTree->SetBranchAddress("detectionefficiency", &deteff);
	
	// Create new output file
	filePartID++;
        TString outFilepath(createDatasetFilename(dirOut, TmpSiPMFilePrefix, nSiPM, opYield, quantumEff, filePartID));
        TFile *output = TFile::Open(outFilepath, "RECREATE");
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
            TString branchName = "Slice";
            branchName += sipm;
            TString branchDesc = branchName + "/L";
            branchSiPM[sipm] = SiPMTree.Branch(branchName, &readouts[sipm], branchDesc);
        }

        for(Long64_t i = 0; i < fTree->GetEntries(); i++){
            fTree->GetEntry(i);
            original_events.insert(eventnumber);
            double angle = atan2(y, x);
            if(angle < 0)
                angle += 2 * PI;
            int segment = angle / thetaSegment;
            assert((segment >= 0) & (segment <nSiPM));
            TRandom rnd = TRandom(RNDSEED);
            for(int sipm = 0; sipm < nSiPM; sipm++){
                readouts[sipm] = 0;
            }
            Long64_t opDetected = round(Edep * opYield * deteff * quantumEff);
	    if(eventnumber==107111){
	    	cout << "\t[Debug] Event: " << eventnumber << ",";
	    	cout << "\t[Debug] Edep: " << Edep << ",";
	    	cout << "\t[Debug] DetEff: " << deteff << ",";
	    	cout << "\t[Debug] Quantum: " << quantumEff<< ",";
	    	cout << "\t[Debug] opDetected: " << opDetected << ",";
	    }
            assert(opDetected>=0);
            for(Long64_t op = 0; op < opDetected; op++){
                int r = round(rnd.Gaus(m, s));
	        if(eventnumber==107111){
			cout << "[Debug] Segment: " << segment << ", offset r: " << r << endl;
                if (segment + r < 0)
                    readouts[segment + r + nSiPM]++;
                else
                    readouts[segment + r]++;
            }
            if(opDetected > 0){
                SiPMTree.Fill();
                produced_events.insert(eventnumber);
            }
        }
	// Debug
        cout << " -> " << fullDirOut << output->GetName() << endl;
        cout << "\tOriginal Events: " << original_events.size() << ",\n";
        cout << "\tProduced Events: " << produced_events.size() << "\n\n";
        // Write output
        SiPMTree.Write();
        NSiPM->Write();
        output->Close();
        output->Delete();
    }
    cout << "[Info] Data preparation: completed.\n";
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

int main(){
    cout << "[Info] Preprocessing...\n";
    // Data cleaning
    const char * dirIn = "/home/data/Muons/";
    const char * dirOut = "/home/data/MuonsPreproc/";
    const char * mapDir = "/home/data/OpticalMaps/";
    // Local
    // const char * dirIn = "../Data/muons/MuonsROI/";
    // const char * dirOut = "../Data/muons/MuonsROI/";
    // const char * mapDir = "../Data/OpticalMaps/";
    // Data cleaning
    Long64_t entry_x_file = 20000000;	//Compact root files to have this number of entries
    data_cleaning(dirIn, dirOut, mapDir);
    compact_data(dirOut, dirOut, TmpROIFilePrefix, OutROIFilePrefix, entry_x_file);
    // Data preparation
    convert_to_sliced_detections(dirOut, dirOut);
    cout << "[Info] End.\n";
}
