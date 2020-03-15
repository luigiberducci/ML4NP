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

void nuclear_to_csv(){
    std::ofstream out("nuclear123456789.csv");
    //   Connect file generated in $ROOTSYS/test
    TFile fileIn("nuclearinfo123456789.root");
    TTree* theTree = nullptr;    
    TTreeReader theReader("nuclearTree", &fileIn);

    TTreeReaderValue<Int_t> pid(theReader, "PID");
    TTreeReaderValue<Double_t> x(theReader, "x");
    TTreeReaderValue<Double_t> y(theReader, "y");
    TTreeReaderValue<Double_t> z(theReader, "z");
    TTreeReaderValue<Double_t> muonpx(theReader, "muonpx");
    TTreeReaderValue<Double_t> muonpy(theReader, "muonpy");
    TTreeReaderValue<Double_t> muonpz(theReader, "muonpz");
    TTreeReaderValue<Int_t> tracknumber(theReader, "tracknumber");
    TTreeReaderValue<Int_t> randomseed(theReader, "randomseed");
    TTreeReaderValue<Int_t> primaryPID(theReader, "primaryPID");
    TTreeReaderValue<string> creatorprocess(theReader, "creatorprocess");
    TTreeReaderValue<Double_t> kineticenergy(theReader, "startingkineticenergy");
    TTreeReaderValue<Double_t> time(theReader, "time");
    
    //Print header
    out << "PID,x,y,z,muonpx,muonpy,muonpz,tracknumber,randomseed,primaryPID,creatorprocess,startingkineticenergy,time";
    out << endl;
    // Print each event by row
    int i=0;
    while(theReader.Next()){
        i++;
        out << *pid << ",";
        out << *x << "," << *y << "," << *z << ",";
        out << *muonpx << "," << *muonpy << "," << *muonpz << ",";
        out << *tracknumber  << "," << *randomseed << "," << *primaryPID << ",";
        out << *creatorprocess << "," << *kineticenergy << "," << *time << ",";
        out << endl;
   }
}
