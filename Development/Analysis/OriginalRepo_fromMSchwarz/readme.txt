The root file has a few different histograms for each layout. The 3D map is named OpticalMap* . 
The tag _Raw is the raw counts per bin and has not been normalized. 
The _Scaled tag has had each bin corrected for the average SiPM QE and been normalized. 
The untagged file has normalized bins, but the average QE has not been included.

Neil thought, it is a good idea to start a variable name with a number...
So one has to use when e.g. drawing this thing:
 TFile* eberhard = TFile::Open("LGND200_14_OpticalMapRealisticDetectorLayout.root", "READ");
dynamic_cast<TH2D*>(eberhard->Get("2DOpticalMap_XYScaled"))->Draw("colz")
