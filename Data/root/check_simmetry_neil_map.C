/*
 * =====================================================================================
 *
 *       Filename:  check_simmetry_neil_map.C
 *
 *    Description:  create ttree with detection efficiency for x, y, z coordinates.
 *                  the goal is to check if there is symmetry w.r.t. Z axis
 *
 *        Version:  1.0
 *        Created:  15/04/2020 18:17:48
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Luigi Berducci (), berducci.1647887@studenti.uniroma1.it
 *        Company:  University of Rome La Sapienza
 *
 * =====================================================================================
 */
#include<fstream>

void check_simmetry_neil_map(){
    TH3D *hMap;
    TString mapFile = "OpticalMapL200XeD.14String.5mm.root";
    TFile *f = TFile::Open(mapFile);
    f->GetObject("ProbMapInterior", hMap);
    cout << "Map loaded\n";

    Double_t x, y, z, de;
    TFile *o = TFile::Open("reflectZAxis.root", "RECREATE");
    cout << "Outfile created\n";
    TTree *fTree = new TTree();
    /* fTree->Branch("x", &x, "x/D"); */
    /* fTree->Branch("y", &y, "y/D"); */
    fTree->Branch("z", &z, "z/D");
    fTree->Branch("detectionefficiency", &de, "detectionefficiency/D");
    cout << "TTree created.\n";

    ofstream output("reflectZAxis.csv");
    output << "x,y,z,detectionefficiency,\n";
    for(Double_t xx=-500; xx<+500; xx=xx+5){
        for(Double_t yy=-500; yy<+500; yy=yy+5){
            for(Double_t zz=-1000; zz<+1000; zz=zz+5){
                x = xx;
                y = yy;
                if(zz < 0)
                    z = -zz;
                else
                    z = zz;
                de = hMap->GetBinContent(hMap->FindBin(xx, yy, z));
                z = zz;
                fTree->Fill();
                output << x << ",";
                output << y << ",";
                output << z << ",";
                output << de << ",\n";
            }
        }
    }
    fTree->Write();
    o->Close();
}
