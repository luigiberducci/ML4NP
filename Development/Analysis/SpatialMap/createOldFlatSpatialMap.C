#define INNER_REGION -1
#define MIDDLE_REGION 0
#define OUTER_REGION +1

Double_t EPS = 0.00001;
Double_t INF = 1000000;
Double_t PI = 3.141592653589793;

pair<Double_t, Double_t> compute_intersection_slice(Double_t x1, Double_t y1, Double_t x2, Double_t y2, Double_t r){
	// Compute line parameters
	Double_t a = y1 - y2;
	Double_t b = x2 - x1;
	Double_t c = x1 * y2 - x2 * y1;
	// Compute the point in line closest to the origin
	double xf = -a*c/(a*a+b*b), yf = -b*c/(a*a+b*b);
	// Based on dist from origin and radius, we know multiplicity
	if (c*c > r*r*(a*a+b*b)+EPS){
	    xf = 0;
	    yf = 0;
	} else if (abs (c*c - r*r*(a*a+b*b)) >= EPS) {
	    double d = r*r - c*c/(a*a+b*b);
	    double mult = sqrt (d / (a*a+b*b));
	    double ax, ay, bx, by;
	    ax = xf + b * mult;
	    bx = xf - b * mult;
	    ay = yf - a * mult;
	    by = yf + a * mult;
	    // assert: a valid intersection has the same direction(both x, y)
	    // if 2 intersections can be 1) only 1 in the same direction, 2) both of them (line from outside)
	    double distance = INF;
	    if((x1-x2)*(x1-ax)>=0 && (y1-y2)*(y1-ay)>=0){
	    	xf = ax;	    
	    	yf = ay;
	        distance = sqrt((x1-ax)*(x1-ax) + (y1-ay)*(y1-ay));	
	    }
	    double distance_b = sqrt((x1-bx)*(x1-bx) + (y1-by)*(y1-by));
	    if((x1-x2)*(x1-bx)>=0 && (y1-y2)*(y1-by)>=0 && (distance_b<distance)){
			xf = bx;	    
			yf = by;
			distance = distance_b;
	    }
	    if(distance == INF){
	    	xf = 0;		// 2 intersection but none of them has the right direction
	    	yf = 0;
	    }
	} // else only 1 (xf, yf)
	pair<Double_t, Double_t> intersection(make_pair(xf, yf));
	return intersection;
}

Int_t getPointPlacement(Double_t x, Double_t y, Double_t inner_r, Double_t outer_r){
	Double_t dist0 = sqrt(x*x + y*y);
	if(dist0 <= inner_r)	// point in [0,inner_r]
		return -1;
	if(dist0 <= outer_r)	// point in (inner_r, outer_r]
		return 0;
	return 1;		// point in (outer_r, *)
}

TRandom * rnd = new TRandom();
const Int_t n_slices = 72;
//array<Double_t, n_slices> runToyOpticsFromPoint(Double_t x1, Double_t y1, Int_t nOptics=360){
void runToyOpticsFromPoint(Double_t x1, Double_t y1, Int_t nOptics, array<Int_t, n_slices> &sliced_detections){
	// Parameters
    Int_t m = 0, s=5;
	for(int i=0; i<n_slices; i++)
		sliced_detections[i]=0;
	for(int i=0; i<nOptics; i++){
		// Compute the slice
        int slice_id = round(rnd->Gaus(m, s));
        if(slice_id < 0)
            slice_id += n_slices;
        sliced_detections[slice_id]++;
	}
}

void createOldFlatSpatialMap(){
	Double_t min_r = 0, max_r = 710, rbins = max_r - min_r + 1;
	Int_t slice_bins = n_slices;
	Int_t min_scaled_s  = -(n_slices / 2);		// in the map, slice id for 72 slices is in [-36, +35]
	Int_t max_scaled_s = (n_slices / 2) - 1;
	Int_t nOpticsPerPoint = 1000;
	// Create map
    TString outfile;
    outfile.Form("ToySpatialMap_%dR_%dSlices_%dops.root", (int)max_r, n_slices, nOpticsPerPoint);
    TFile * file = new TFile(outfile, "RECREATE");
	TH2D * map = new TH2D("SpatialMap", "R-Slice Map", rbins, min_r, max_r, slice_bins, min_scaled_s, max_scaled_s);
	for(Double_t x=min_r; x<=max_r; x += 1){
		array<Int_t, n_slices> sliced_detections;
		runToyOpticsFromPoint(x, 0., nOpticsPerPoint, sliced_detections);
		for(int s=0; s<n_slices; s++){
			cout << sliced_detections[s] << ", ";
			if(sliced_detections[s]>0){
				// Scale to center current slice in 0 and the rest in +-N/2
				Int_t scaled_s = s;
				if(s >= n_slices/2)
					scaled_s = - (n_slices - s);
                for(int k=0; k<sliced_detections[s]; k++)
    				map->Fill(x, scaled_s);
			}
		}
		cout << endl;
	}
    map->Write();
    file->Close();
}
