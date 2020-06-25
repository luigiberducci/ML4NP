// File: DetectionSlicer.cpp
// Author: Luigi Berducci
// Purpose: generate a spatial distribution of detections,
//          starting from Energy deposited in LAr.
//
// Usage:
// $>root
// root [0].L DetectionSlicer.cpp
// root [0]DetectionSlicer("Muons", "output", "Output")

#define RNDSEED 123456789
#define PI 3.14159265

#define ROI_MINX -700
#define ROI_MAXX +700
#define ROI_MINY -700
#define ROI_MAXY +700
#define ROI_MINZ -845
#define ROI_MAXZ +845

// Flag for output scheme
const Bool_t writeOnlyROIEntries = false;           // If `true`, all the entries outside the ROI are ignored
const Bool_t writeOnlyNonZeroDetections = true;     // If `true`, all the entries wt NPE=0 are ignored
// Parameter for spatial distribution of detections
const Int_t NSLICES = 72;           // Number of Slices to segment the X-Y plane
const Int_t OPYIELD = 40000;           // Number of Optical Photons per KeV
const Double_t QUANTUMEFF = 0.40;   // Quantum Efficiency of SiPMs
const Double_t M = 0, S = 5;        // Mean, StdDev of Gaussian used to sample where PE are detected
// Optical Map
TFile * mapFile;
TFile * spatMapFile;
TH3D * hMap;
TH2D * sMap;

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
    const char * outTemplate = "%s/%s_Slices%d_Yield%d_QuantumEff%f_Seed%d.root";
    outFilePath.Form(outTemplate, fullDirOut, outFilePrefix.Data(), NSLICES, OPYIELD, QUANTUMEFF, RNDSEED);
    return outFilePath;
}

Bool_t isInROI(Double_t x, Double_t y, Double_t z){
    // Return `true` if the given coordinate is contained in the ROI, otherwise `false`.
    if((x >= ROI_MINX & x <= ROI_MAXX) & (y >= ROI_MINY & y <= ROI_MAXY) & (z >= ROI_MINZ & z <= ROI_MAXZ))
        return true;
    return false;
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

// TODO doc
void convertSingleFile(TString inFilePath, TString outFilePath, TString treeName="fTree"){
    // Open file and get tree
    TFile *f = TFile::Open(inFilePath);
    TTree *simTree = (TTree*) f->Get(treeName);
    // Connect branches
    Double_t x, y, z, r, time, Edep;
    Int_t eventnumber;
    simTree->SetBranchAddress("x", &x);
    simTree->SetBranchAddress("y", &y);
    simTree->SetBranchAddress("z", &z);
    simTree->SetBranchAddress("time", &time);
    simTree->SetBranchAddress("inc_eventnumber", &eventnumber);    // To avoid overlap of events
    simTree->SetBranchAddress("energydeposition", &Edep);
    // Create new tree
    cout << "\tProcessing " << inFilePath << "..." << endl;

    // Create out file and out tree
    TFile *output = TFile::Open(outFilePath, "RECREATE");
    TTree SlicedTree(treeName, "Sliced Detections Tree");
    // New variables
    Double_t deteff, quantumEff = QUANTUMEFF;
    Int_t pedetected;
    array<Int_t, NSLICES> readouts;
    TParameter<Int_t> *nSlicesParam = new TParameter<Int_t>("NSlices", NSLICES);
    TBranch *bN = SlicedTree.Branch("eventnumber", &eventnumber, "eventnumber/I");
    TBranch *bT = SlicedTree.Branch("time", &time, "time/D");
    TBranch *bX = SlicedTree.Branch("x", &x, "x/D");
    TBranch *bY = SlicedTree.Branch("y", &y, "y/D");
    TBranch *bZ = SlicedTree.Branch("z", &z, "z/D");
    TBranch *bR = SlicedTree.Branch("r", &r, "r/D");
    TBranch *bE = SlicedTree.Branch("energydeposition", &Edep, "energydeposition/D");
    TBranch *bPE = SlicedTree.Branch("pedetected", &pedetected, "pedetected/I");
    TBranch *bDE = SlicedTree.Branch("detectionefficiency", &deteff, "detectionefficiency/D");
    TBranch *bQE = SlicedTree.Branch("quantumefficiency", &quantumEff, "quantumefficiency/D");
    array<TBranch*, NSLICES> branchSlices;
    for(int slice = 0; slice < NSLICES; slice++){
        TString branchName = "Slice";
        branchName += slice;
        TString branchDesc = branchName + "/I";
        branchSlices[slice] = SlicedTree.Branch(branchName, &readouts[slice], branchDesc);
    }
    // Loop over tree entries
    TRandom rnd = TRandom(RNDSEED);
    Int_t kAllEvents = 0, lastReadEventNr = -1;
    Int_t kROIEvents = 0, lastReadROIEventNr = -1;
    Int_t kWrittenEvents = 0, lastWrittenEventNr = -1;
    for(Long64_t i = 0; i < simTree->GetEntries(); i++){
        simTree->GetEntry(i);
        Bool_t entryIsInROI = isInROI(x, y, z);
        // Debug (Just some event statistics)
        if(eventnumber > lastReadEventNr){
            kAllEvents++;
            lastReadEventNr = eventnumber;
        }else if(eventnumber < lastReadEventNr){
	    cout << "[Warning] Found a decreasing order of event number: risk of overlapping the entries.\n";
	}
        if((writeOnlyROIEntries) & (!entryIsInROI)) continue;
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
        assert(deteff >= 0.0 & deteff <= 1.0);
        // Compute the radius
        r = sqrt(x*x + y*y);
        // Compute the current slice
        int currentSliceID = getCurrentSliceID(x, y, NSLICES);
        assert((currentSliceID >= 0) & (currentSliceID < NSLICES));
        // Spread detections
        pedetected = round(Edep * OPYIELD * deteff * QUANTUMEFF);
        assert(pedetected>=0);
        for(int slice = 0; slice < NSLICES; slice++){
            readouts[slice] = 0;     // Reset readouts
        }
        TH1D * hitspace_dist = sMap->ProjectionY("py", (Int_t) round(r));    // take the distribution based on distance from origin
        assert(pedetected==0 || hitspace_dist->ComputeIntegral()>0);
        for(Long64_t pe = 0; pe < pedetected; pe++){
            int activated_slice = (currentSliceID + (int)round(hitspace_dist->GetRandom())) % NSLICES;
            assert((activated_slice > -NSLICES) & (activated_slice < NSLICES));
            if (activated_slice < 0)
                readouts[activated_slice + NSLICES]++;
            else
                readouts[activated_slice]++;
        }
        if((!writeOnlyNonZeroDetections) || (pedetected > 0)){
            SlicedTree.Fill();
            if(eventnumber > lastWrittenEventNr){
                kWrittenEvents++;
                lastWrittenEventNr = eventnumber;
            }
        }
        // Safety check for writeOnlyNonZeroDetections==false
        // 1) writeOnlyROIEntries==true => kWrittenEvents == kROIEvents
        // 2) writeOnlyROIEntries==false => kWrittenEvents == kAllEvents
        if(writeOnlyNonZeroDetections==false){
            assert((!writeOnlyROIEntries) || (kWrittenEvents == kROIEvents));
            assert((writeOnlyROIEntries) || (kWrittenEvents == kAllEvents));
        }
    }
    // Debug
    cout << "\t[Info] Number of events in simulation file: " << kAllEvents << " ";
    cout << "(" << kROIEvents << " in ROI)" << endl;
    cout << "\t[Info] Number of events produced: " << kWrittenEvents << endl;
    cout << "\t[Info] Output written in " << output->GetName() << endl << endl;
    // Write output
    SlicedTree.Write();
    nSlicesParam->Write();
    output->Close();
    output->Delete();
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

void loadSpatialMap(const char * mapDir="OpticalMaps",
                    const char * mapFileName = "ToySpatialMap_710R_72Slices_1000ops",
                    const char * mapObjectName = "SpatialMap"){
    TString mapFilePath;
    mapFilePath.Form("%s/%s.root", mapDir, mapFileName);
    spatMapFile = TFile::Open(mapFilePath);
    spatMapFile->GetObject(mapObjectName, sMap);
    cout << "[Info] Loaded Spatial Map " << mapFilePath << endl;
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
    // Create printout strings
    TString inout_printout, roi_printout, npe_printout;
    inout_printout.Form(inout_template, fullDirIn, inFilePrefix);
    if(writeOnlyROIEntries)
        roi_printout.Form(roi_template, "ENABLED ", "Only the entries in ROI will be considered");
    else
        roi_printout.Form(roi_template, "DISABLED", "All the entries (even outside the ROI) will be considered");
    if(writeOnlyNonZeroDetections)
        npe_printout.Form(npe_template, "ENABLED ", "Only the entries wt NPE>0 will be written");
    else
        npe_printout.Form(npe_template, "DISABLED", "All the entries (even NPE==0) will be written");
    // Print them out
    cout << header << endl;
    cout << inout_printout << endl;
    cout << roi_printout << endl;
    cout << npe_printout << endl;
    cout << endl;
}

void DetectionSlicer(const char * dirIn="./Muons", const char * inFilePrefix="output", const char * dirOut="./Output"){
    // Produce the "sliced" detections for each of the specified files.
    // Params:  `dirin` is the path of the input directory
    //          `inFilePrefix` is the prefix of the simulation files (e.g. "output" for output123456789.root)
    //          `dirOut` is the directory where you want to write the output files
    auto start = std::chrono::system_clock::now();
    loadOpticalMap();
    loadSpatialMap();
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
    cout << "[Info] Results written in " << fullDirOut << "/\n";
}
