/*
 * =====================================================================================
 *
 *       Filename:  merge_trees.C
 *
 *    Description:  merge trees from multiple root files
 *
 *        Version:  1.0
 *        Created:  31/03/2020 19:33:45
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Luigi Berducci (), berducci.1647887@studenti.uniroma1.it
 *        Company:  University of Rome La Sapienza
 *
 * =====================================================================================
 */

void merge_trees(){
    // IO params
    const char * dirIn = "../Data/cutROI/root/";
    char* fullDirIn = gSystem->ExpandPathName(dirIn);
    // Define chain
    TChain ch("fTree");
    // Loop files in input dir
    void* dirp = gSystem->OpenDirectory(fullDirIn);
    const char* entry;
    int file_k = 0, file_group = 20, n_group = 0;
    while((entry = (char*)gSystem->GetDirEntry(dirp))) {
        // Filter only output####.root files
        TString fileName = entry;
        TString begin = fileName(0, 6);
        if(!(begin.EqualTo("cutROI") & fileName.EndsWith("root")))
            continue;
        if(file_k >= file_group){
            n_group++;      // increment group id
            file_k = 0;     // reset counter
            TString fileOut(dirIn);
            fileOut += "output2eROI_part";
            fileOut += n_group;
            fileOut += ".root";
            ch.Merge(fileOut);
            ch.Reset();
            // Debug
            std::cout << "[Info] Merged in " << fileOut << endl << endl;
        }
        ch.Add(dirIn + fileName);
        file_k++;
        // Debug
        std::cout << "[Info] Add " << dirIn + fileName << endl;
    }
    n_group++;
    TString fileOut(dirIn);
    fileOut += "output2eROI_part";
    fileOut += n_group;
    fileOut += ".root";
    ch.Merge(fileOut);
    // Debug
    std::cout << "[Info] Merged in " << fileOut << endl << endl;
}
