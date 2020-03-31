/*
 * =====================================================================================
 *
 *       Filename:  provaROI.C
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  31/03/2020 16:21:40
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Luigi Berducci (), berducci.1647887@studenti.uniroma1.it
 *        Company:  University of Rome La Sapienza
 *
 * =====================================================================================
 */

void cut_output_ROI() {
    // IO params
    const char * dirIn = "../Data/root/";
    const char * dirOut = "../Data/root/ROI/";
    char* fullDirIn = gSystem->ExpandPathName(dirIn);
    char* fullDirOut = gSystem->ExpandPathName(dirOut);
    // Loop files in input dir
    void* dirp = gSystem->OpenDirectory(fullDirIn);
    const char* entry;
    while((entry = (char*)gSystem->GetDirEntry(dirp))) {
        // Filter only output####.root files
        TString fileName = entry;
        TString begin = fileName(0, 6);
        if(!(begin.EqualTo("output") & fileName.EndsWith("root")))
            continue;
        // Debug
        std::cout << "[Info] Processing " << dirIn + fileName << endl;
        // Open input file and take original tree
        TFile* file = TFile::Open(dirIn + fileName);
        TTree* originalTree = (TTree*)file->Get("fTree");
        // Create (or overwrite) output file
        TFile* output = TFile::Open(dirOut + fileName.Copy().Replace(0, 6, "cutROI"), "RECREATE");
        // Define cut w.r.t. coordinates
        // ROI def: x, y in [-500, +500], z in [-1000, +1000]
        TString cutX("(x >= -500.) & (x <= 500.)");
        TString cutY("(y >= -500.) & (y <= 500.)");
        TString cutZ("(z >= -1000) & (z <= 1000)");
        // Copy tree
        TTree* selectedTree = originalTree->CopyTree(cutX + " & " + cutY + " & " + cutZ);
        std::cout << "[Info] original entries: " << originalTree->GetEntries() << ", ";
        std::cout << "ROI entries: " << selectedTree->GetEntries() << ", ";
        std::cout << "perc cut: " << (1 - ((double)selectedTree->GetEntries() / originalTree->GetEntries())) << " %\n\n";
        // Write out file and close
        output->Write();
        output->Close();
    }
    gSystem->FreeDirectory(dirp);   // Free pointer to dir
}
