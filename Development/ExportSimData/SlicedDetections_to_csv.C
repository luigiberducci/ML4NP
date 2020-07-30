/*
 * =====================================================================================
 *
 *       Filename:  macro.C
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  12/02/2020 12:58:21
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Luigi Berducci (), berducci.1647887@studenti.uniroma1.it
 *        Company:  University of Rome La Sapienza
 *
 * =====================================================================================
 */
#include<assert.h>

bool isRootFile(TString fileName, TString prefix){
    // Consider files `prefixXXXXXX.root`
    TString begin = fileName(0, prefix.Length());
    if(begin.EqualTo(prefix) & fileName.EndsWith("root"))
        return true;
    return false;
}

void produceSingleExport(TString inFilepath){
    // IO params
    TFile fileIn(inFilepath);
    // Create out file name
    TString outFileName(TString(inFilepath(0, inFilepath.Length()-5)) + ".csv");
    // Tree
    TTree* fTree = (TTree*) fileIn.Get("fTree");
    // Params
    TParameter<Int_t> * NInnerSlicesParam = (TParameter<Int_t>*) fileIn.Get("NInnerSlices");
    TParameter<Int_t> * NOuterSlicesParam = (TParameter<Int_t>*) fileIn.Get("NOuterSlices");
    // Branches
    Double_t x, y, z, r, time, energydeposition, detectionefficiency, quantumefficiency;
    Int_t eventnumber, PID, pedetected;
    string * material = 0;
    vector<Int_t> innerSlices, outerSlices;
    for(int j=0; j<NInnerSlicesParam->GetVal(); j++)
	innerSlices.push_back(0);
    for(int j=0; j<NOuterSlicesParam->GetVal(); j++)
	outerSlices.push_back(0);
    // Connect branches
    fTree->SetBranchAddress("eventnumber", &eventnumber);
    fTree->SetBranchAddress("PID", &PID);
    fTree->SetBranchAddress("time", &time);
    fTree->SetBranchAddress("x", &x);
    fTree->SetBranchAddress("y", &y);
    fTree->SetBranchAddress("z", &z);
    fTree->SetBranchAddress("r", &r);
    fTree->SetBranchAddress("material", &material);
    fTree->SetBranchAddress("energydeposition", &energydeposition);
    fTree->SetBranchAddress("pedetected", &pedetected);
    fTree->SetBranchAddress("detectionefficiency", &detectionefficiency);
    fTree->SetBranchAddress("quantumefficiency", &quantumefficiency);
    for(int j=0; j<NInnerSlicesParam->GetVal(); j++){
        TString branchName = "InnerSlice";
        branchName += j;
        fTree->SetBranchAddress(branchName, &innerSlices[j]);
    }
    for(int j=0; j<NOuterSlicesParam->GetVal(); j++){
        TString branchName = "OuterSlice";
        branchName += j;
        fTree->SetBranchAddress(branchName, &outerSlices[j]);
    }
    // Loop over entries
    std::ofstream out;
    for(int i=0; i<fTree->GetEntries(); i++){
        fTree->GetEntry(i);
        // Check if we have a new event and reached the max num in curr file
        if(i == 0){
            cout << "\t[Info] Creating " << outFileName << endl;
            if(i > 0)    out.close();
	        out.open(outFileName);
            //Print header
            out << "eventnumber,PID,time,x,y,z,r,material,";
            out << "energydeposition,pedetected,detectionefficiency,quantumefficiency,";
            for(int j=0; j<NInnerSlicesParam->GetVal(); j++)
                out << "InnerSlice" << j << ",";
            for(int j=0; j<NOuterSlicesParam->GetVal(); j++)
                out << "OuterSlice" << j << ",";
            out << "\n";
	    }
        out << eventnumber << "," << PID << "," << time << "," << x << "," << y << "," << z << "," << r << ",";
        out << *material << "," << setprecision(15) << energydeposition << "," << pedetected << ",";
        out << detectionefficiency << "," << quantumefficiency << ",";
        for(int j=0; j<NInnerSlicesParam->GetVal(); j++)
            out << innerSlices[j] <<  ",";
        for(int j=0; j<NOuterSlicesParam->GetVal(); j++)
            out << outerSlices[j] <<  ",";
	out << "\n";
   }
}

void SlicedDetections_to_csv(const char * dirIn="./", const char * prefixIn="SlicedDetections"){
    // Loop files in input directory
    char* fullDirIn = gSystem->ExpandPathName(dirIn);
    void* dirp = gSystem->OpenDirectory(fullDirIn);
    const char* entry;
    while((entry = (char*)gSystem->GetDirEntry(dirp))) {
        TString fileName = entry;
        if(!isRootFile(fileName, prefixIn))   continue;
        cout << "\t[Info] Loading " << fullDirIn + fileName << endl;   // Debug
        produceSingleExport(fileName);
    }
}
