
#include "TFile.h"
#include "TTree.h"

/*
Use e.g.
 .L toyMC.cpp
 TFile::Open("peSpectrum.root", "RECREATE")
 TH1D* dist = getPhotonNrDist(1000000, 700, 845, "peSpectrum39Ar")
 dist->Write()
*/

const char* mapFile = "../Data/root/OpticalMapL200XeD.14String.5mm.root";

TH3D* getMap(){
	TFile* file = TFile::Open(mapFile, "READ");
	return dynamic_cast<TH3D*>(file->Get("OpticalMap_Scaled"));
}



TH1D* getBetaPDF(const std::string& isotope = "Ar39");
int getBin(double, double, double, int);
bool pointValid(TH3D* map, double xPos, double yPos, double zPos);
double readPoint(TH3D* map, double xPos, double yPos, double zPos);
int getPhotonNumber(const std::string& isotope = "Ar39");	//from random, using 39Ar power spectrum & light yield
void printPoint(TH3D* map, double xPos, double yPos, double zPos);

TRandom3 modnar;
TH1D* ar39Histo = getBetaPDF();
TH3D* heatMap = getMap();

/*
* Gets the distribution of the number of detected photons in the detector.
* Each event ( = entry in the histo) is one decay of the isotope of interest (39Ar now)
* events: number of total decays to be sampled
* the sampled region is a cylinder, with radius radius and height = halfHeight * 2
* full heat map has xmin/xmax: -700/700, ymin/ymax: -700/700, zmin/zmax: -850/845
*    --> call getPhotonNrDist(..., 700, 845) should span full map
*/
TH1D* getPhotonNrDist(uint32_t events, double radius, double halfHeight, const std::string& histName = ""){
	
	const char* name = (histName == "") ? std::to_string((int)(modnar.Rndm()*10000)).c_str() : histName.c_str();

	TH1D* distr = new TH1D(name, name, 250, 0., 249.);
	
	uint32_t currentEvents = 0;
	while(true){

		if(currentEvents%1000==0) std::cout << "\rat: " << currentEvents << std::flush;
	
		double x = modnar.Rndm() * radius * 2 - radius;			//sample in a box first
		double y =  modnar.Rndm() * radius * 2 - radius;
		double z = modnar.Rndm() * halfHeight * 2 - halfHeight;

		if(x*x + y*y > radius*radius) continue;				//cut cylinder within box
		//printPoint(heatMap, x, y, z);
		
		double mapEff = readPoint(heatMap, x, y, z);		//TODO check if the position makes sense; i.e. if it is filled w/ LAr
															//for now: assume, that complete volume is LAr, as approximation
		
		if(mapEff < 0.) continue;			

		int phNr = getPhotonNumber();

		int detPhNr = modnar.Binomial(phNr, mapEff);
		//std::cout << "  "<< phNr << "  " << mapEff << "  "<<detPhNr <<  std::endl;

		distr->Fill(detPhNr);
	
		currentEvents++;	
		if(currentEvents >= events) break;
	}
	return distr;
}






int getPhotonNumber(const std::string& isotope){
	TH1D* histo;
	if(isotope == "Ar39") histo = ar39Histo;

	if(histo == NULL) throw 123;
	
	double yield = 40;	// in photons per keV

	double energyInkeV = histo->GetRandom();
	
	double expectPhotons = energyInkeV * yield;

	//TODO check if we have to use poisson statistics for determining nr of photons produced in the scintillation process

	return modnar.Poisson(expectPhotons);
}

TH1D* getBetaPDF(const std::string& isotope){
	ifstream dataFile;
	if(isotope == "Ar39"){
		dataFile.open("../Data/ar39/input_Ar39_ESpectrum.txt");
	}
	double e, f, df;

	uint32_t validLines = 0;
	while(true){
		std::cout << "\rat: " << validLines << std::flush;
		if(dataFile.peek() == '#'){ 
			char dummy[1024];
			dataFile.getline(dummy, 1024);	//throw comment away
			continue;	//jump over comments 
		}
		dataFile >> e >> f >> df;
		if(!dataFile.good()) break;
		validLines++;
	}
	std::cout << std::endl << "valid lines: " << validLines << ", last e: " << e << std::endl;

	TH1D* histo = new TH1D(isotope.c_str(), isotope.c_str(), validLines, 0., e + 1);

	dataFile.clear();	//clear EOF marker
	dataFile.seekg(0, ios::beg);

	validLines = 0;
	while(true){
		std::cout << "\rat: " << validLines << std::flush;
		if(dataFile.peek() == '#'){ 
			char dummy[1024];
			dataFile.getline(dummy, 1024);	//throw comment away
			continue;	//jump over comments 
		}
		dataFile >> e >> f >> df;
		if(!dataFile.good()) break;
		//std::cout << validLines + 1  << "; " << f << "; " << df << std::endl;
		histo->SetBinContent(validLines + 1, f); //bin0=underflow, 1...nBins: from xlow(inc) to xhigh(exc), nBins+1=overflow
		histo->SetBinError(validLines + 1, df);
		validLines++;
	}
	std::cout << std::endl;

	dataFile.close();

	return histo;
}

bool pointValid(){
	return true;	//TODO
}

double readPoint(TH3D* map, double xPos, double yPos, double zPos){
	int x_bin = getBin(xPos, map->GetXaxis()->GetXmin(), map->GetXaxis()->GetXmax(), map->GetNbinsX());
	int y_bin = getBin(yPos, map->GetYaxis()->GetXmin(), map->GetYaxis()->GetXmax(), map->GetNbinsY());
	int z_bin = getBin(zPos, map->GetZaxis()->GetXmin(), map->GetZaxis()->GetXmax(), map->GetNbinsZ());
	
	return map->GetBinContent(x_bin, y_bin, z_bin);
}

void printPoint(TH3D* map, double xPos, double yPos, double zPos){
	int x_bin = getBin(xPos, map->GetXaxis()->GetXmin(), map->GetXaxis()->GetXmax(), map->GetNbinsX());
	int y_bin = getBin(yPos, map->GetYaxis()->GetXmin(), map->GetYaxis()->GetXmax(), map->GetNbinsY());
	int z_bin = getBin(zPos, map->GetZaxis()->GetXmin(), map->GetZaxis()->GetXmax(), map->GetNbinsZ());
	std::cout << "Point "<<"("<<xPos<<","<<yPos<<","<<zPos<<") aka ("<<x_bin<<","<<y_bin<<","<<z_bin<<")"<< " --> " <<map->GetBinContent(x_bin, y_bin, z_bin) << std::endl;
}

int getBin(double x, double xmin, double xmax, int bin_nr){
	if(x < xmin || x > xmax) return -1;
	return (x - xmin) / (xmax - xmin) * bin_nr;
}






double getSum(){
	TFile* file = TFile::Open(mapFile, "READ");
	//std::cout << file->GetName() << std::endl;

	TH3D* map_scaled = dynamic_cast<TH3D*>(file->Get("OpticalMap_Scaled"));
	
	double sum_var = 0.;

	for(int i = 0; i < map_scaled->GetNbinsX(); i++){
		for(int j = 0; j < map_scaled->GetNbinsY(); j++){
			for(int k = 0; k < map_scaled->GetNbinsZ(); k++){
				sum_var += map_scaled->GetBinContent(i, j, k);
			}
		}
	}

	return sum_var;
}

/*main*/
void toyMC(){
	//nix
}


