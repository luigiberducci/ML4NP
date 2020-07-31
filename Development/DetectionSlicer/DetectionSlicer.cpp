// File: DetectionSlicer.cpp
// Author: Luigi Berducci
// Purpose: generate a spatial distribution of detections,
//          starting from Energy deposited in LAr.
//
// Usage (compiled):
// $>g++ DetectionSlicer.cpp -o DetectionSlicer.exe `root-config --cflags --glibs`
// $>./DetectionSlicer.exe
// Usage (interactive):
// $>root
// root [0].L DetectionSlicer.cpp
// root [0]DetectionSlicer("Muons", "output", "Output")

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

#define ROI_MINX -700
#define ROI_MAXX +700
#define ROI_MINY -700
#define ROI_MAXY +700
#define ROI_MINZ -845
#define ROI_MAXZ +845
#define ARGONMATERIAL "ArgonLiquid"
#define MAP_MAX_RADIUS 705

using namespace std;

// Flat for input format
const Bool_t useMaterialBranch = true;              // If `true`, use `material` branch to compute NPE for Argon only
// Flag for output scheme
const Bool_t writeOnlyROIEntries = false;           // If `true`, all the entries outside the ROI are ignored
const Bool_t writeOnlyNonZeroDetections = false;     // If `true`, all the entries wt NPE=0 are ignored
const Bool_t writeOnlyArgonLiquidEntries = false;     // If `true`, all the entries wt material!='ArgonLiquid' are ignored
// Parameter for spatial distribution of detections
const Int_t N_INNERSLICES = 12;		// Number of Slices to segment the X-Y plane of Inner Shroud
const Int_t N_OUTERSLICES = 20;		// Number of Slices to segment the X-Y plane of Outer Shroud
const Int_t OPYIELD = 40;		// Number of Optical Photons per KeV
const Double_t QUANTUMEFF = 0.20;	// Quantum Efficiency of SiPMs
// Optical Map
TFile * mapFile;
TFile * spatMapFile;
TH3D * hMap;
TH1D * spatialInnerDet;
TH2D * spatialInnerMap;
TH2D * spatialOuterMap;
array<TH1D*, MAP_MAX_RADIUS> spatialInnerMapProjections;
array<TH1D*, MAP_MAX_RADIUS> spatialOuterMapProjections;

Bool_t isRootFile(TString fileName, TString prefix){
    // Return `true` if the filename `fileName` is a ROOT file starting wt `prefix`, otherwise `false`.
    // Note: path is not managed, this function expect only the filename.
    TString begin = fileName(0, prefix.Length());
    if(begin.EqualTo(prefix) & fileName.EndsWith("root"))
        return true;
    return false;
}

TString createOutputFilepath(TString inFileName, TString inFilePrefix, const char * fullDirOut, const char * outBasePrefix="SlicedDetections"){
    // Create the path of the output file, starting from the input filename.
    // This function replaces the input prefix `inFilePrefix` with the output prefix `outBasePrefix`,
    // and includes the parameters used to reproduce the detections.
    // Params:  `inFileName` is the filename of the input file (path in not managed)
    //          `inFilePrefix` is the prefix of the input file
    //          `fullDirOut` is the path of output directory
    //          `outBasePrefix` is the prefix of the output files
    // Returns: the path of the output file.
    TString outFilePrefix(TString(inFileName(0, inFileName.Length()-5)).Replace(0, inFilePrefix.Length(), outBasePrefix));
    TString outFilePath;
    const char * outTemplate = "%s/%s_Slices%d_%d_Yield%d_QuantumEff%f_Seed%d.root";
    outFilePath.Form(outTemplate, fullDirOut, outFilePrefix.Data(), N_INNERSLICES, N_OUTERSLICES, OPYIELD, QUANTUMEFF, RNDSEED);
    return outFilePath;
}

Bool_t isInROI(Double_t x, Double_t y, Double_t z){
    // Return `true` if the given coordinate is contained in the ROI, otherwise `false`.
    if((x >= ROI_MINX & x <= ROI_MAXX) & (y >= ROI_MINY & y <= ROI_MAXY) & (z >= ROI_MINZ & z <= ROI_MAXZ))
        return true;
    return false;
}

Int_t getShiftedSliceID(Double_t x, Double_t y, Double_t shift_angle, Int_t nSlices){
    // Return the ID in [0, nSlices) of the slice that contains the (x,y) coordinate, shifted with angle `shift_angle`.
    Double_t thetaSlice = 2 * PI / nSlices;
    Double_t angle = atan2(y, x);
    if(angle < 0)
        angle += 2 * PI;
    angle += shift_angle;
    if(angle < 0)
        angle += 2 * PI;
    return (int) round(angle / thetaSlice) % nSlices;
}

Int_t getCurrentSliceID(Double_t x, Double_t y, Int_t nSlices){
    // Return the ID in [0, nSlices) of the slice that contains the (x,y) coordinate.
    Double_t thetaSlice = 2 * PI / nSlices;
    Double_t angle = atan2(y, x);
    if(angle < 0)
        angle += 2 * PI;
    return angle / thetaSlice;
}

Double_t getDetectionEfficiency(Double_t x, Double_t y, Double_t z){
    // Return the Detection Efficiency in (x, y, z).
    if(!isInROI(x, y, z))
        return 0;
    Int_t bin = hMap->FindBin(x, y, z);
    return hMap->GetBinContent(bin);
}

void convertSingleFile(TString inFilePath, TString outFilePath, TString treeName="fTree"){
    // Open file and get tree
    TFile *f = TFile::Open(inFilePath);
    TTree *simTree = (TTree*) f->Get(treeName);
    // Connect branches
    Double_t x, y, z, r, time, Edep;
    Int_t eventnumber, PID;
    std::string * material = 0;
    simTree->SetBranchAddress("PID", &PID);
    simTree->SetBranchAddress("x", &x);
    simTree->SetBranchAddress("y", &y);
    simTree->SetBranchAddress("z", &z);
    simTree->SetBranchAddress("time", &time);
    simTree->SetBranchAddress("eventnumber", &eventnumber);    // To avoid overlap of events
    simTree->SetBranchAddress("energydeposition", &Edep);
    if(useMaterialBranch==true)
        simTree->SetBranchAddress("material", &material);
    // Create new tree
    cout << "\tProcessing " << inFilePath << "...\n";

    // Create out file and out tree
    TFile *output = TFile::Open(outFilePath, "RECREATE");
    TTree SlicedTree(treeName, "Sliced Detections Tree");
    // New variables
    Double_t deteff, quantumEff = QUANTUMEFF;
    Int_t pedetected;
    array<Int_t, N_INNERSLICES> inner_readouts;
    array<Int_t, N_OUTERSLICES> outer_readouts;
    TParameter<Int_t> *nInnerSlicesParam = new TParameter<Int_t>("NInnerSlices", N_INNERSLICES);
    TParameter<Int_t> *nOuterSlicesParam = new TParameter<Int_t>("NOuterSlices", N_OUTERSLICES);
    TParameter<Bool_t> *flagUseMaterialBranch = new TParameter<Bool_t>("FlagUseMaterialBranch", useMaterialBranch);
    TParameter<Bool_t> *flagWriteOnlyROIEntries = new TParameter<Bool_t>("FlagWriteOnlyROIEntries", writeOnlyROIEntries);
    TParameter<Bool_t> *flagWriteOnlyNonZeroDet = new TParameter<Bool_t>("FlagWriteOnlyNonZeroDet", writeOnlyNonZeroDetections);
    TParameter<Bool_t> *flagWriteOnlyLArEntries = new TParameter<Bool_t>("FlagWriteOnlyLArEntries", writeOnlyArgonLiquidEntries);
    TBranch *bN = SlicedTree.Branch("eventnumber", &eventnumber, "eventnumber/I");
    TBranch *bPID = SlicedTree.Branch("PID", &PID, "PID/I");
    TBranch *bT = SlicedTree.Branch("time", &time, "time/D");
    TBranch *bX = SlicedTree.Branch("x", &x, "x/D");
    TBranch *bY = SlicedTree.Branch("y", &y, "y/D");
    TBranch *bZ = SlicedTree.Branch("z", &z, "z/D");
    TBranch *bR = SlicedTree.Branch("r", &r, "r/D");
    TBranch *bM = SlicedTree.Branch("material", &material);
    TBranch *bE = SlicedTree.Branch("energydeposition", &Edep, "energydeposition/D");
    TBranch *bPE = SlicedTree.Branch("pedetected", &pedetected, "pedetected/I");
    TBranch *bDE = SlicedTree.Branch("detectionefficiency", &deteff, "detectionefficiency/D");
    TBranch *bQE = SlicedTree.Branch("quantumefficiency", &quantumEff, "quantumefficiency/D");
    array<TBranch*, N_INNERSLICES> branchInnerSlices;
    array<TBranch*, N_OUTERSLICES> branchOuterSlices;
    for(int slice = 0; slice < N_INNERSLICES; slice++){
        TString branchName = "InnerSlice";
        branchName += slice;
        TString branchDesc = branchName + "/I";
        branchInnerSlices[slice] = SlicedTree.Branch(branchName, &inner_readouts[slice], branchDesc);
    }
    for(int slice = 0; slice < N_OUTERSLICES; slice++){
        TString branchName = "OuterSlice";
        branchName += slice;
        TString branchDesc = branchName + "/I";
        branchOuterSlices[slice] = SlicedTree.Branch(branchName, &outer_readouts[slice], branchDesc);
    }
    // Loop over tree entries
    TRandom rnd = TRandom(RNDSEED);
    Int_t kAllEvents = 0, lastReadEventNr = -1;
    Int_t kROIEvents = 0, lastReadROIEventNr = -1;
    Int_t kWrittenEvents = 0, lastWrittenEventNr = -1;
    Int_t kArgonEvents = 0, lastArgonEvent = -1;
    Int_t nEntries = simTree->GetEntries();
    for(Long64_t i = 0; i < nEntries; i++){
	if(i % 10000 == 0)
		cout << "\r\tentry: " << i << "/" << nEntries << std::flush;
        simTree->GetEntry(i);
        Bool_t entryIsInROI = isInROI(x, y, z);
        // Debug (Just some event statistics)
        if(eventnumber > lastReadEventNr){
            kAllEvents++;
            lastReadEventNr = eventnumber;
        }else if(eventnumber < lastReadEventNr){
	    cout << "\n[Warning] Found a decreasing order of event number: risk of overlapping the entries.\n";
	}
        if((writeOnlyROIEntries) & (!entryIsInROI)) continue;
        if((useMaterialBranch) && (writeOnlyArgonLiquidEntries) && (*material!=ARGONMATERIAL)) continue;
	if((useMaterialBranch) && (*material==ARGONMATERIAL) && (eventnumber > lastArgonEvent)){
	    kArgonEvents++;
	    lastArgonEvent = eventnumber;
	}
        // Compute Detection Efficiency
        if(entryIsInROI) {
            if (eventnumber > lastReadROIEventNr) {
                kROIEvents++;
                lastReadROIEventNr = eventnumber;
            }
            deteff = getDetectionEfficiency(x, y, z);
        }else{
            deteff = 0.0;
        }
	if(useMaterialBranch && *material!=ARGONMATERIAL)
            deteff = 0.0;	// force 0 detection for Ge entries
        assert(deteff >= 0.0 & deteff <= 1.0);
        // Compute the radius
        r = sqrt(x*x + y*y);
        // Spread detections
        pedetected = round(Edep * OPYIELD * deteff * QUANTUMEFF);
        assert(pedetected>=0);
        for(int slice = 0; slice < N_INNERSLICES; slice++){
            inner_readouts[slice] = 0;     // Reset readouts
        }
        for(int slice = 0; slice < N_OUTERSLICES; slice++){
            outer_readouts[slice] = 0;     // Reset readouts
        }
	// Compute number of PE detected from inner and outer shrouds
	Double_t inner_fraction = spatialInnerDet->GetBinContent(spatialInnerDet->GetBin(r));
	//cout << inner_fraction<< endl;
	Int_t inner_pe = round(inner_fraction * pedetected);
	Int_t outer_pe = round((1 - inner_fraction) * pedetected);
	if(inner_pe + outer_pe > pedetected){	//because of rounding error
		if(outer_pe > inner_pe)
			outer_pe--;
		else
			inner_pe--;
	}
	// take the distribution based on distance from origin
	if(r <= MAP_MAX_RADIUS){
		TH1D * inner_hitspace_dist = spatialInnerMapProjections[(Int_t)round(r)];
		TH1D * outer_hitspace_dist = spatialOuterMapProjections[(Int_t)round(r)];
		assert(inner_pe==0 || inner_hitspace_dist->ComputeIntegral()>0);
		assert(outer_pe==0 || outer_hitspace_dist->ComputeIntegral()>0);
		// Spread on Inner Shroud
		for(Long64_t pe = 0; pe < inner_pe; pe++){
		    // Compute random detection angle accordin to spatial distribution
		    Double_t det_angle = inner_hitspace_dist->GetRandom();
		    Int_t activated_slice = getShiftedSliceID(x, y, det_angle, N_INNERSLICES);
		    assert((activated_slice >= 0) & (activated_slice < N_INNERSLICES));
		    inner_readouts[activated_slice]++;
		}
		// Spread on Outer Shroud
		for(Long64_t pe = 0; pe < outer_pe; pe++){
		    // Compute random detection angle accordin to spatial distribution
		    Double_t det_angle = outer_hitspace_dist->GetRandom();
		    Int_t activated_slice = getShiftedSliceID(x, y, det_angle, N_OUTERSLICES);
		    assert((activated_slice >= 0) & (activated_slice < N_OUTERSLICES));
		    outer_readouts[activated_slice]++;
		}
		// DEBUG COUNT
		Int_t innerPE = 0, outerPE = 0;
		for(int slice = 0; slice < N_INNERSLICES; slice++){
		    innerPE += inner_readouts[slice];
		}
		for(int slice = 0; slice < N_OUTERSLICES; slice++){
		    outerPE += outer_readouts[slice];
		}
		//cout << "PEDet: " << pedetected << ", ";
		//cout << "Written PE: " << innerPE<< " + " << outerPE << "= " << innerPE+outerPE << ",\n";
		assert(innerPE+outerPE== pedetected);
	}
	// To allow cut on Ge entries, we write the entries in other materials (!=Argon).
	// Since Det.Eff. is 0 for other materials, these entries will have 0 PE
	// If WriteOnlyNonZeroDetecions => (PE>0 or Germanium)
        if((!writeOnlyNonZeroDetections) || (pedetected > 0 || (useMaterialBranch && *material!=ARGONMATERIAL))){
            SlicedTree.Fill();
            if(eventnumber > lastWrittenEventNr){
                kWrittenEvents++;
                lastWrittenEventNr = eventnumber;
            }
        }
        // Safety check for writeOnlyNonZeroDetections==false and writeOnlyArgonEntries==false
        // 1) writeOnlyROIEntries==true => kWrittenEvents == kROIEvents
        // 2) writeOnlyROIEntries==false => kWrittenEvents == kAllEvents
        if(writeOnlyNonZeroDetections==false && writeOnlyArgonLiquidEntries==false){
            assert((!writeOnlyROIEntries) || (kWrittenEvents == kROIEvents));
            assert((writeOnlyROIEntries) || (kWrittenEvents == kAllEvents));
        }
    }
    // Write output
    nInnerSlicesParam->Write();
    nOuterSlicesParam->Write();
    flagUseMaterialBranch->Write();
    flagWriteOnlyROIEntries->Write();
    flagWriteOnlyNonZeroDet->Write();
    flagWriteOnlyLArEntries->Write();
    SlicedTree.Write();
    output->Close();
    output->Delete();
    // Debug
    cout << endl;
    cout << "\t[Info] Number of events in simulation file: " << kAllEvents << " ";
    cout << "(" << kROIEvents << " in ROI, ";
    cout << kArgonEvents << " wt Argon entries)" << endl;
    cout << "\t[Info] Number of events produced: " << kWrittenEvents << endl;
    cout << "\t[Info] Output written in " << output->GetName() << endl << endl;
}

void loadOpticalMap(const char * mapDir="OpticalMaps",
                    const char * mapFileName = "LGND200_14_OpticalMapRealisticDetectorLayout",
                    const char * mapObjectName = "OpticalMap_Scaled"){
    // Note: this method works with the OP map in Regular LAr
    // TODO: extend with interior/exterior map for Xenon map
    TString mapFilePath;
    mapFilePath.Form("%s/%s.root", mapDir, mapFileName);
    mapFile = TFile::Open(mapFilePath);
    mapFile->GetObject(mapObjectName, hMap);
    cout << "[Info] Loaded Optical Map " << mapFilePath << endl;
}

void loadSpatialMaps(const char * mapDir="OpticalMaps",
                    const char * mapFileName = "ToySpatialMap_711R_100AngleSlices_10000ops_AttLen500.000000_NewSampling_Last",
                    const char * PrInnerDetObjectName = "PrInnerDet",
                    const char * InnermapObjectName = "InnerMap",
                    const char * OutermapObjectName = "OuterMap"){
    TString mapFilePath;
    mapFilePath.Form("%s/%s.root", mapDir, mapFileName);
    spatMapFile = TFile::Open(mapFilePath);
    spatMapFile->GetObject(PrInnerDetObjectName, spatialInnerDet);
    spatMapFile->GetObject(InnermapObjectName, spatialInnerMap);
    spatMapFile->GetObject(OutermapObjectName, spatialOuterMap);
    cout << "[Info] Loaded Spatial Maps from " << mapFilePath << endl;
    // Create Projections to save time
    for(Int_t r=0; r<MAP_MAX_RADIUS; r++){
        spatialInnerMapProjections[r] = spatialInnerMap->ProjectionY("py", r, r+1);
        spatialOuterMapProjections[r] = spatialOuterMap->ProjectionY("py", r, r+1);
    }
    cout << "[Info] Loaded Radius-Projections from Spatial Maps\n";
}

void writeHeaderInfo(const char* fullDirIn, const char * inFilePrefix="output"){
    // Print info for the user about the output scheme
    TString header = "";
    header += "   ------------------------------------------------------------------\n";
    header += "  |                         Detection Slicer                          |\n";
    header += "  |                                                                   |\n";
    header += "  |        To convert Energy Depositions in LAr to Detections         |\n";
    header += "  |   using the Optical Map and sampling the hit-space distribution   |\n";
    header += "  |-------------------------------------------------------------------|\n";
    header += "  |   Author: Luigi Berducci, INFN RomaTre                            |\n";
    header += "  |   Date:   May 2020                                                |\n";
    header += "  |   For any questions: luigi.berducci@roma3.infn.it                 |\n";
    header += "   ------------------------------------------------------------------\n";

    // Templates
    TString inout_template("[Info] Simulation Files: %s/%s*.root");
    TString roi_template("[Info] CUT ROI %s: %s.");
    TString npe_template("[Info] CUT NPE %s: %s.");
    TString arg_template("[Info] CUT ARG %s: %s.");
    // Create printout strings
    TString inout_printout, roi_printout, arg_printout, npe_printout;
    inout_printout.Form(inout_template, fullDirIn, inFilePrefix);
    if(writeOnlyROIEntries)
        roi_printout.Form(roi_template, "ENABLED ", "Only the entries in ROI will be considered");
    else
        roi_printout.Form(roi_template, "DISABLED", "All the entries (even outside the ROI) will be considered");
    if(writeOnlyArgonLiquidEntries)
        arg_printout.Form(arg_template, "ENABLED ", "Only the entries in Argon material will be considered");
    else
        arg_printout.Form(arg_template, "DISABLED", "All the entries (even other materials) will be considered");
    if(writeOnlyNonZeroDetections)
        npe_printout.Form(npe_template, "ENABLED ", "Only the entries wt NPE>0 will be written (or other materials, if enabled)");
    else
        npe_printout.Form(npe_template, "DISABLED", "All the entries (even NPE==0) will be written");
    // Print them out
    cout << header << endl;
    cout << inout_printout << endl;
    cout << roi_printout << endl;
    cout << arg_printout << endl;
    cout << npe_printout << endl;
    cout << endl;
}

void DetectionSlicer(const char * dirIn="./Input", const char * inFilePrefix="output", const char * dirOut="./Output"){
    // Produce the "sliced" detections for each of the specified files.
    // Params:  `dirin` is the path of the input directory
    //          `inFilePrefix` is the prefix of the simulation files (e.g. "output" for output123456789.root)
    //          `dirOut` is the directory where you want to write the output files
    gSystem->SetAclicMode(TSystem::kDebug);
    auto start = std::chrono::system_clock::now();
    loadOpticalMap();
    loadSpatialMaps();
    char* fullDirIn = gSystem->ExpandPathName(dirIn);
    char* fullDirOut = gSystem->ExpandPathName(dirOut);
    void* dirp = gSystem->OpenDirectory(fullDirIn);
    const char* entry;
    writeHeaderInfo(fullDirIn, inFilePrefix);
    while((entry = (char*)gSystem->GetDirEntry(dirp))) {
        if(!isRootFile(entry, inFilePrefix))   continue;
        TString fileName = entry;
        TString inFilePath;
        inFilePath.Form("%s/%s", fullDirIn, entry);
        TString outFilePath(createOutputFilepath(fileName, inFilePrefix, fullDirOut));
        convertSingleFile(inFilePath, outFilePath);
    }
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    cout << "[Info] Convertion completed in " << elapsed_seconds.count() << " seconds.\n";
    cout << "[Info] Look for results in " << fullDirOut << "/\n";
}

int main(){
	DetectionSlicer();
}
