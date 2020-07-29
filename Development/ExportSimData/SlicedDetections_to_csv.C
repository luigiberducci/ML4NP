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

void produce_single_export(TString inFilepath){
    // IO params
    TFile fileIn(inFilepath);
    TObjArray *filepath = inFilepath.Tokenize("/");     // Assume unix-like path
    TString outFileBase = ((TObjString *)(filepath->At(filepath->GetEntries() - 1)))->String();
    TObjArray *filename= outFileBase.Tokenize(".");
    outFileBase = ((TObjString *)(filename->At(0)))->String();
    // Tree
    TTree* fTree = (TTree*) fileIn.Get("fTree");
    // Branches
    TParameter<Int_t> * NInnerSlicesParam = (TParameter<Int_t>*) fileIn.Get("NInnerSlices");
    TParameter<Int_t> * NOuterSlicesParam = (TParameter<Int_t>*) fileIn.Get("NOuterSlices");
    Int_t nInnerSlices = NInnerSlicesParam->GetVal();
    Int_t nOuterSlices = NOuterSlicesParam->GetVal();
    // Loop over entries
    std::ofstream out;
    for(int i=0; i<fTree->GetEntries(); i++){
        fTree->GetEntry(i);
        // Check if we have a new event and reached the max num in curr file
        if(i == 0){
            std::cout << "[Info] Creating " << outFileBase << ".csv" << endl;
            if(i > 0)    out.close();
	        out.open(outFileBase + ".csv");
            //Print header
            /* out << "PID,ParentTrackID,energydeposition,kineticenergy,time,x,y,z,"; */
            /* out << "px,py,pz,eventnumber,tracknumber,creatorprocess,process,parentnucleusPID,"; */
            /* out << "muonx,muony,muonz,muonpx,muonpy,muonpz,muonenergy,"; */
            /* out << "material,detectornumber"; */
            /* out << endl; */
	    }
        /* out << *pid << "," << *p_trace_id << ","; */
        /* out << setprecision(20) <<  *energydeposition << "," << *kineticenergy << "," << *time << ","; */
        /* out << *x << "," << *y << "," << *z << ","; */
        /* out << *px << "," << *py << "," << *pz << ","; */
        /* out << *eventnumber << "," << *tracknumber  << "," << *creatorprocess << "," << *process << ","; */
        /* out << *parentnucleusPID << ","; */
        /* out << *muonx << "," << *muony << "," << *muonz << ","; */
        /* out << *muonpx << "," << *muonpy << "," << *muonpz << ","; */
        /* out << *muone << ","; */
        /* out << *material << "," << *detectornumber; */
        /* out << endl; */
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
        cout << "\t" << fullDirIn + fileName << endl;   // Debug
        produce_single_export(fileName);
    }
}
