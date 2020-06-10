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


void output_to_csv(TString inFilepath){
    // IO params
    TFile fileIn(inFilepath);
    TObjArray *filepath = inFilepath.Tokenize("/");     // Assume unix-like path
    TString outFileBase = ((TObjString *)(filepath->At(filepath->GetEntries() - 1)))->String();
    TObjArray *filename= outFileBase.Tokenize(".");
    outFileBase = ((TObjString *)(filename->At(0)))->String();
    // Tree
    TTree* theTree = nullptr;
    TTreeReader theReader("fTree", &fileIn);
    // Branches
    TTreeReaderValue<Int_t> pid(theReader, "PID");
    TTreeReaderValue<Int_t> p_trace_id(theReader, "ParentTrackID");
    TTreeReaderValue<Double_t> energydeposition(theReader, "energydeposition");
    TTreeReaderValue<Double_t> kineticenergy(theReader, "kineticenergy");
    TTreeReaderValue<Double_t> time(theReader, "time");
    TTreeReaderValue<Double_t> x(theReader, "x");
    TTreeReaderValue<Double_t> y(theReader, "y");
    TTreeReaderValue<Double_t> z(theReader, "z");
    TTreeReaderValue<Double_t> px(theReader, "px");
    TTreeReaderValue<Double_t> py(theReader, "py");
    TTreeReaderValue<Double_t> pz(theReader, "pz");
    TTreeReaderValue<Double_t> muonx(theReader, "muonx");
    TTreeReaderValue<Double_t> muony(theReader, "muony");
    TTreeReaderValue<Double_t> muonz(theReader, "muonz");
    TTreeReaderValue<Double_t> muonpx(theReader, "muonpx");
    TTreeReaderValue<Double_t> muonpy(theReader, "muonpy");
    TTreeReaderValue<Double_t> muonpz(theReader, "muonpz");
    TTreeReaderValue<Double_t> muone(theReader, "muonenergy");
    TTreeReaderValue<Double_t> detectionefficiency(theReader, "detectionefficiency");
    TTreeReaderValue<Int_t> eventnumber(theReader, "inc_eventnumber");
    TTreeReaderValue<Int_t> tracknumber(theReader, "tracknumber");
    TTreeReaderValue<string> creatorprocess(theReader, "creatorprocess");
    // TTreeReaderValue<string> creatorprocess(theReader, "process");
    TTreeReaderValue<Int_t> parentnucleusPID(theReader, "parentnucleusPID");

    // Loop over entries
    int i=0;
    Int_t maxNumberEntriesInFile = 2500000;
    Int_t outFileNumber = 0;          // part id (for out filename)
    std::ofstream out;
    while(theReader.Next()){
        Int_t cur_event = *eventnumber;
        // Check if we have a new event and reached the max num in curr file
	if(i % maxNumberEntriesInFile == 0){
            outFileNumber += 1;
            std::cout << "[Info] Creating " << outFileBase + "_part" + outFileNumber + ".csv" << endl;
            if(i > 0)    out.close();
	    out.open(outFileBase + "_part" + outFileNumber + ".csv");
            //Print header
            out << "PID,ParentTrackID,energydeposition,kineticenergy,time,x,y,z,";
            out << "px,py,pz,eventnumber,tracknumber,creatorprocess,parentnucleusPID,detectionefficiency,";
            out << "muonx,muony,muonz,muonpx,muonpy,muonpz,muonenergy,";
            out << endl;
	    }
        i++;
        out << *pid << "," << *p_trace_id << ",";
        out << setprecision(20) <<  *energydeposition << "," << *kineticenergy << "," << *time << ",";
        out << *x << "," << *y << "," << *z << ",";
        out << *px << "," << *py << "," << *pz << ",";
        out << *eventnumber << "," << *tracknumber  << "," << *creatorprocess << ",";
        out << *parentnucleusPID << "," << *detectionefficiency << ",";
        out << *muonx << "," << *muony << "," << *muonz << ",";
        out << *muonpx << "," << *muonpy << "," << *muonpz << ",";
        out << *muone << ",";
        out << endl;
   }
   std::cout << "\n[Info] Written " << i << " entries in " << outFileNumber << " files.\n";
}
