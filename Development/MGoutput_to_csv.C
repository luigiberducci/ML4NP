/*
 * =====================================================================================
 *
 *       Filename:  MGoutput_to_csv.C
 *
 *    Description:  csv export from MG output scheme
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

void MGoutput_to_csv(){
    // IO params
    TString dir = "../Data/SimsForGiuseppe/root/";
    TString fileName = "LGND20014String_bulk_A228_Z90_in_Cables_100000_-132362.root";
    std::ofstream out("output_AmbRad_A228_Z90_100000.csv");
    // create the tree from input file
    TChain *fTree = new TChain("fTree");
    fTree->Add(dir + fileName);
    // extract num entries
    Long64_t nentries = (Long64_t)fTree->GetEntries();
    // create step branches
    MGTMCEventSteps *eventSteps = 0;
    MGTMCEventSteps *eventPrimaries = 0;
    const MGTMCStepData *step,*primaries;
    if(fTree != NULL){
        fTree->SetBranchAddress("eventSteps", &eventSteps);
        fTree->SetBranchAddress("eventPrimaries", &eventPrimaries);
    }
    else{
        cout<<"NULL fTree"<<endl;
        return;
    }
    //Print header
    out << "PID,ParentTrackID,energydeposition,kineticenergy,time,x,y,z,";
    out << "px,py,pz,eventnumber,tracknumber,creatorprocess,parentnucleusPID";
    out << endl;
    for(Int_t i = 0; i < nentries ; i++){
        fTree->GetEntry(i);
        TString physName;
        primaries = eventPrimaries->GetStep(0);
        Int_t eventnumber = eventSteps->GetEventID();
        for (Int_t j = 0; j < eventSteps->GetNSteps();j++){
            step = eventSteps->GetStep(j);

            Int_t pid = step->GetParticleID();
            Int_t p_trace_id = step->GetParentTrackID();
            Double_t energydeposition = step->GetEdep();
            Double_t kineticenergy = step->GetKineticE();
            Double_t time = step->GetT();
            Double_t x = step->GetX();
            Double_t y = step->GetY();
            Double_t z = step->GetZ();
            Double_t px = step->GetPx();
            Double_t py = step->GetPy();
            Double_t pz = step->GetPz();
            Int_t tracknumber = step->GetTrackID();
            TString creatorprocess = step->GetProcessName();
            Int_t parentnucleusPID = 666;   // to mark that is not implemented

            // Print each event by row
            out << pid << "," << p_trace_id << ",";
            out << energydeposition << "," << kineticenergy << "," << time << ",";
            out << x << "," << y << "," << z << ",";
            out << px << "," << py << "," << pz << ",";
            out << eventnumber << "," << tracknumber  << "," << creatorprocess << ",";
            out << parentnucleusPID << ",";
            out << endl;
        }
    }
}
