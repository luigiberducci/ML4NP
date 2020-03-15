/*
 * =====================================================================================
 *
 *       Filename:  dat_to_tree.C
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  10/03/2020 17:56:18
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Luigi Berducci (), berducci.1647887@studenti.uniroma1.it
 *        Company:  University of Rome La Sapienza
 *
 * =====================================================================================
 */
#include <fstream>
#include <sstream>
#include <string>

void dat_to_tree(){
   std::ifstream inFile("test.dat");     // input dat file
   // create root file with root tree
   TFile f("musun.root", "recreate");
   TTree t1("t1", "tree with muon starting position");
   Float_t x, y, z, energy, theta, phi;
   Int_t ev;
   t1.Branch("x",&x,"x/F");
   t1.Branch("y",&y,"y/F");
   t1.Branch("z",&z,"z/F");
   t1.Branch("energy",&energy,"energy/F");
   t1.Branch("theta",&theta,"theta/F");
   t1.Branch("phi",&phi,"phi/F");
   t1.Branch("eventnumber",&ev,"eventnumber/I");

   // fill the tree
   std::string line;
   ev = 0;
   while (std::getline(inFile, line)) {
        // process according to line format
        std::stringstream ss(line);
        Int_t id, unknown;
        ss >> id >> unknown >> energy >> x >> y >> z >> theta >> phi;
        /* std::cout << "ev:" << id << " x:" << x << " y:" << y << " z:" << z << " energy:" << energy << " theta:" << theta << " phi:" << phi << endl; */
        t1.Fill();
        ev = ev + 1;
   }
   // save the Tree heade; the file will be automatically closed
   // when going out of the function scope
   t1.Write();
}
