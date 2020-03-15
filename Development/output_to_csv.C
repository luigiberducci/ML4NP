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

void output_to_csv(){
    std::ofstream out("output123456789.csv");
    //   Connect file generated in $ROOTSYS/test
    TFile fileIn("output123456789.root");
    TTree* theTree = nullptr;    
    TTreeReader theReader("fTree", &fileIn);

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
    TTreeReaderValue<Double_t> muonenergy(theReader, "muonenergy");
    TTreeReaderValue<Int_t> eventnumber(theReader, "eventnumber");
    TTreeReaderValue<Int_t> stepnumber(theReader, "stepnumber");
    TTreeReaderValue<Double_t> steplength(theReader, "steplength");
    TTreeReaderValue<Int_t> tracknumber(theReader, "tracknumber");
    TTreeReaderValue<Int_t> detectornumber(theReader, "detectornumber");
    TTreeReaderValue<string> creatorprocess(theReader, "creatorprocess");
    TTreeReaderValue<Double_t> startx(theReader, "startx");
    TTreeReaderValue<Double_t> starty(theReader, "starty");
    TTreeReaderValue<Double_t> startz(theReader, "startz");
    TTreeReaderValue<Int_t> randomseed(theReader, "randomseed");
    TTreeReaderValue<Int_t> parentnucleusPID(theReader, "parentnucleusPID");
    TTreeReaderValue<Double_t> nuclearx(theReader, "nuclearx");
    TTreeReaderValue<Double_t> nucleary(theReader, "nucleary");
    TTreeReaderValue<Double_t> nuclearz(theReader, "nuclearz");

    //Print header
    out << "PID,ParentTrackID,energydeposition,kineticenergy,time,x,y,z,";
    out << "px,py,pz,muonx,muony,muonz,muonpx,muonpy,muonpz,muonenergy,eventnumber,stepnumber,steplength,tracknumber,detectornumber,creatorprocess,";
    out << "startx,starty,startz,randomseed,parentnucleusPID,nuclearx,nucleary,nuclearz";
    out << endl;
    // Print each event by row
    int i=0;
    while(theReader.Next()){
        i++;
        out << *pid << "," << *p_trace_id << ",";
        out << *energydeposition << "," << *kineticenergy << "," << *time << ",";
        out << *x << "," << *y << "," << *z << ",";
        out << *px << "," << *py << "," << *pz << ",";
        out << *muonx << "," << *muony << "," << *muonz << ",";
        out << *muonpx << "," << *muonpy << "," << *muonpz << ",";
        out << *muonenergy << ",";
        out << *eventnumber << "," << *stepnumber << "," << *steplength << ",";
        out << *tracknumber  << "," << *detectornumber  << "," << *creatorprocess << ",";
        out << *startx << "," << *starty << "," << *startz << ",";
        out << *randomseed << "," << *parentnucleusPID << ",";
        out << *nuclearx << "," << *nucleary << "," << *nuclearz << ",";
        out << endl;
   }
}
