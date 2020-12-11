#include <set>

bool isRootFile(TString fileName, TString prefix){
    // Consider files `prefixXXXXXX.root`
    TString begin = fileName(0, prefix.Length());
    if(begin.EqualTo(prefix) & fileName.EndsWith("root"))
        return true;
    return false;
}

int processFile(const char* dirIn, const char* startfile){
      TFile *input = new TFile(Form("%s/%s", dirIn, startfile),"READ");
      TTree *fTree = (TTree*)input->Get("fTree");

      int entries = fTree->GetEntries();
      double energydeposition = 0;
      int pedetected = 0;
      int eventnumber = 0;
      int storedeventnumber = -1;
      double totalpe = 0;
      double totalargonenergy = 0;
      double totalgeenergy = 0;
      string *material = NULL;

     fTree->SetBranchAddress("energydeposition",&energydeposition);
     fTree->SetBranchAddress("pedetected",&pedetected);
     fTree->SetBranchAddress("eventnumber",&eventnumber);
     fTree->SetBranchAddress("material",&material);

     fTree->GetEntry(0);
     storedeventnumber = eventnumber;
     set<int> eventnumbers;
     int k_critical_events = 0;
     for (int i=i;i<entries;i++)
       {//For all entries
         fTree->GetEntry(i);

         if(eventnumber!=storedeventnumber)
           {//Event is complete, begin energy summing of another event

             if(totalgeenergy>=1839 & totalgeenergy<=2239 & totalargonenergy<=20000){ //Total Ge energy deposition has to exceed 200 keV
                cout << "\tEvent: " << storedeventnumber << ": GeEnergy: " << totalgeenergy << " KeV, NPE: " << totalpe << endl;
                k_critical_events++;
                eventnumbers.insert(storedeventnumber);
                }

             totalgeenergy = 0;
             totalpe = 0;
             totalargonenergy = 0;
             storedeventnumber = eventnumber;
           }

          if(strstr(material->c_str(),"ArgonLiquid"))
            totalargonenergy+=energydeposition;
          if(strstr(material->c_str(),"GermaniumEnriched"))
            totalgeenergy+=energydeposition;
        totalpe+=pedetected;
       }//For all entries
     if(totalgeenergy>=1839 & totalgeenergy<=2239 & totalargonenergy<=20000){ //Total Ge energy deposition has to exceed 200 keV
        cout << "\tEvent: " << storedeventnumber << ": GeEnergy: " << totalgeenergy << " KeV, NPE: " << totalpe << endl;
        k_critical_events++;
        eventnumbers.insert(storedeventnumber);
    }
    int i_e = 0;
    for(auto critical_event: eventnumbers){
        i_e++;
        TFile *out = new TFile(Form("ExportCriticalEvent%d_%d.root", i_e, critical_event), "RECREATE");
        TTree* oTree = fTree->CopyTree(Form("eventnumber==%d", critical_event));
        oTree->Write();
        out->Close();
    }
     input->Close();
    return k_critical_events;
}//EOF

void listGeEvents(const char * dirIn = "PostProc/"){
    void *dirp = gSystem->OpenDirectory(dirIn);
    const char *direntry;
    int k_critical_events = 0;
    while((direntry = (char*) gSystem->GetDirEntry(dirp))){
        TString fileName = direntry;
        if(!isRootFile(fileName, "SlicedDetections"))    continue;
        cout << "File: " << direntry << "\n";
        k_critical_events += processFile(dirIn, direntry);
    }
    cout << endl;
    cout << "Critical Events: " << k_critical_events << endl;
    gSystem->FreeDirectory(dirp);
}
