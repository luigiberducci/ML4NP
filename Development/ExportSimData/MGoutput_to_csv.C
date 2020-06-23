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
#include<assert.h>

void MGoutput_to_csv(){
    // IO params
    TString dir = "./";
    TString fileName = "LGND_200_example_Danila.root";
    TString outFileBase = "LGND_200_example_Danila";
    // create the tree from input file
    TChain *fTree = new TChain("fTree");
    fTree->Add(dir + fileName);
    // extract num entries
    Long64_t nentries = (Long64_t)fTree->GetEntries();
    // create step branches
    MGMCRun *mcRun = 0;
    MGTMCEventHeader *eventHeader = 0;
    MGTMCEventSteps *eventSteps = 0;
    MGTMCEventSteps *eventPrimaries = 0;
    const MGTMCStepData *step,*primaries;
    if(fTree != NULL){
        fTree->SetBranchAddress("fMCRun", &mcRun);
        fTree->SetBranchAddress("eventHeader", &eventHeader);
        fTree->SetBranchAddress("eventSteps", &eventSteps);
        fTree->SetBranchAddress("eventPrimaries", &eventPrimaries);
    }
    else{
        cout<<"NULL fTree"<<endl;
        return;
    }
    std::cout << "[Info] Total number of entries: " << nentries << endl << endl;
    Int_t maxNumberEventInFile = 25000;   
    Int_t lastEventIDInFile = -1;     // last eventID in the current file
    Int_t counterEventInFile = 0;     // num of event in the current file
    Int_t outFileNumber = 0;          // part id (for out filename)
    std::ofstream out;
    for(Int_t i = 0; i < nentries ; i++){
	    fTree->GetEntry(i);
        Int_t eventnumber = eventSteps->GetEventID();
        // Check if we have a new event and reached the max num in curr file
	    if((eventnumber > lastEventIDInFile) & (counterEventInFile % maxNumberEventInFile == 0)){
            outFileNumber += 1;
            std::cout << "[Info] Creating " << outFileBase + "_part" + outFileNumber + ".csv" << endl;
            if(i > 0)    out.close();
	        out.open(outFileBase + "_part" + outFileNumber + ".csv");
            //Print header
            out << "PID,ParentTrackID,energydeposition,kineticenergy,time,x,y,z,";
            out << "px,py,pz,eventnumber,tracknumber,creatorprocess,volumeID,parentnucleusPID";
            out << endl;
            counterEventInFile = 0;
	    }
        if(eventnumber > lastEventIDInFile)
            counterEventInFile += 1;
        // Extract primary event
        primaries = eventPrimaries->GetStep(0);     // should be only 1 step with id=0
	    if(primaries == NULL)
	        continue;
        Int_t pri_pid = primaries->GetParticleID();
        Int_t pri_p_trace_id = primaries->GetParentTrackID();
        Double_t pri_energydeposition = primaries->GetEdep();
        Double_t pri_kineticenergy = primaries->GetKineticE();
        Double_t pri_time = primaries->GetT();
        Double_t pri_x = primaries->GetX();
        Double_t pri_y = primaries->GetY();
        Double_t pri_z = primaries->GetZ();
        Double_t pri_px = primaries->GetPx();
        Double_t pri_py = primaries->GetPy();
        Double_t pri_pz = primaries->GetPz();
        Int_t pri_tracknumber = primaries->GetTrackID();
        TString pri_creatorprocess = primaries->GetProcessName();
        Int_t pri_sensitivevolID= primaries->GetSensitiveVolumeID();
        Int_t pri_parentnucleusPID = 666;   // to mark that is not implemented
        // Print each event by row
        out << pri_pid << "," << pri_p_trace_id << ",";
        out << setprecision(20) << pri_energydeposition << "," << pri_kineticenergy << "," << pri_time << ",";
        out << pri_x << "," << pri_y << "," << pri_z << ",";
        out << pri_px << "," << pri_py << "," << pri_pz << ",";
        out << eventnumber << "," << pri_tracknumber  << "," << pri_creatorprocess << ",";
        out << pri_sensitivevolID << "," << pri_parentnucleusPID << ",";
        out << endl;
	    // Extract steps
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
            Int_t sensitivevolID= step->GetSensitiveVolumeID();
            Int_t parentnucleusPID = 666;   // to mark that is not implemented

            // Print each event by row
            out << pid << "," << p_trace_id << ",";
            out << setprecision(20) << energydeposition << "," << kineticenergy << ",";
	    out << time << ",";
            out << x << "," << y << "," << z << ",";
            out << px << "," << py << "," << pz << ",";
            out << eventnumber << "," << tracknumber  << "," << creatorprocess << ",";
            out << sensitivevolID << "," << parentnucleusPID << ",";
            out << endl;
        }
        // Update last event ID in file
        assert(eventnumber >= lastEventID);	// we assume eventIDs are ordered
        lastEventIDInFile = eventnumber;
    }
    std::cout << endl << "[Info] End."<< endl;
}
