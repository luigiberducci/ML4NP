/*
 * =====================================================================================
 *
 *       Filename:  data_cleaning.C
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06/04/2020 15:39:00
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Luigi Berducci (), berducci.1647887@studenti.uniroma1.it
 *        Company:  University of Rome La Sapienza
 *
 * =====================================================================================
 */
#include <math.h>

void data_cleaning(){
    double m=0, s=5;        // tuned for 72 sipms
    const int opYield=40;   // 40 ops / KeV
    const int nSiPM = 72;
    // input
    double x = 0, y = 0, z = 0;
    double Edep = 3500;  // KeV
    double opDE = 0.001;    // computed based on x,y,z
    int segment = 0;
    // artificial segmentation
    TRandom rnd = TRandom();
    array<int, nSiPM> readout;
    TH1* h3 = new TH1I("h1", "normal distribution", 100, -40, +40); //keep margin
    for(int i=0; i<ceil(Edep * opYield * opDE); i++){
        int r = round(rnd.Gaus(m, s));
        cout << r << " ";
        if (segment + r < 0)
            readout[segment + r + nSiPM]++;
        else
            readout[segment + r]++;

        h3->Fill(segment + r);
    }
    cout << endl;

    h3->SetFillColor(kGreen);
    THStack *hs = new THStack("hs","");
    /* hs->Add(h1); */
    /* hs->Add(h2); */
    hs->Add(h3);
    hs->Draw();
}
