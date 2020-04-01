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


void output_to_csv(){
    // IO params
    TString inFileBase = "output2eROI_part5";
    TFile fileIn("../Data/cutROI/root/" + inFileBase + ".root");
    TString outFileBase = inFileBase;
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
    TTreeReaderValue<Int_t> eventnumber(theReader, "eventnumber");
    TTreeReaderValue<Int_t> tracknumber(theReader, "tracknumber");
    TTreeReaderValue<string> creatorprocess(theReader, "creatorprocess");
    TTreeReaderValue<Int_t> parentnucleusPID(theReader, "parentnucleusPID");

    // Loop over entries
    int i=0;
    Int_t maxNumberEventInFile = 20000;
    Int_t lastEventIDInFile = -1;     // last eventID in the current file
    Int_t counterEventInFile = 0;     // num of event in the current file
    Int_t outFileNumber = 0;          // part id (for out filename)
    std::ofstream out;
    while(theReader.Next()){
        i++;
        Int_t cur_event = *eventnumber;
        assert(cur_event >= lastEventID);	// we assume eventIDs are ordered
        // Check if we have a new event and reached the max num in curr file
	    if((cur_event > lastEventIDInFile) & (counterEventInFile % maxNumberEventInFile == 0)){
            outFileNumber += 1;
            std::cout << "[Info] Creating " << outFileBase + "_part" + outFileNumber + ".csv" << endl;
            if(i > 0)    out.close();
	        out.open(outFileBase + "_part" + outFileNumber + ".csv");
            //Print header
            out << "PID,ParentTrackID,energydeposition,kineticenergy,time,x,y,z,";
            out << "px,py,pz,eventnumber,tracknumber,creatorprocess,parentnucleusPID";
            out << endl;
            counterEventInFile = 0;
	    }
        if(cur_event > lastEventIDInFile){
            counterEventInFile += 1;
            lastEventIDInFile = cur_event;
        }
        out << *pid << "," << *p_trace_id << ",";
        out << *energydeposition << "," << *kineticenergy << "," << *time << ",";
        out << *x << "," << *y << "," << *z << ",";
        out << *px << "," << *py << "," << *pz << ",";
        out << *eventnumber << "," << *tracknumber  << "," << *creatorprocess << ",";
        out << *parentnucleusPID << ",";
        out << endl;
   }
   std::cout << "\n[Info] Written " << i << " entries in " << outFileNumber << " files.\n";
}
