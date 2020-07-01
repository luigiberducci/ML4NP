#include<assert.h>

#define INNER_REGION -1
#define MIDDLE_REGION 0
#define OUTER_REGION +1
#define INNER_SHROUD 0
#define OUTER_SHROUD 1
#define UNKWN_SHROUD -1

Double_t EPS = 0.00001;
Double_t INF = 1000000;
Double_t PI = 3.141592653589793;

pair<Double_t, Double_t> compute_closest_intersection(Double_t x1, Double_t y1, Double_t x2, Double_t y2, Double_t r){
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

pair<Double_t, Double_t> compute_farthest_intersection(Double_t x1, Double_t y1, Double_t x2, Double_t y2, Double_t r){
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
	    double distance = -INF;
	    if((x1-x2)*(x1-ax)>=0 && (y1-y2)*(y1-ay)>=0){
	    	xf = ax;	    
	    	yf = ay;
	        distance = sqrt((x1-ax)*(x1-ax) + (y1-ay)*(y1-ay));	
	    }
	    double distance_b = sqrt((x1-bx)*(x1-bx) + (y1-by)*(y1-by));
	    if((x1-x2)*(x1-bx)>=0 && (y1-y2)*(y1-by)>=0 && (distance_b>distance)){
			xf = bx;	    
			yf = by;
			distance = distance_b;
	    }
	    if(distance == -INF){
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

// Parameters
const Int_t n_inner_slices = 100;
const Int_t n_outer_slices = 100;
const Double_t inner_theta_slice = 2 * PI / n_inner_slices;    // angle for each slice
const Double_t outer_theta_slice = 2 * PI / n_outer_slices;    // angle for each slice
void runToyOpticsFromPoint(Double_t x1, Double_t y1, Int_t nOptics, TH1D * &prInnerD, TH2D* &innerMap, TH2D* &outerMap){
	cout << "Running at R=" << x1 << endl;
	// Fiber radius
	Double_t inner_r = 175, outer_r = 295;	// mm
	Double_t delta_phi = PI / nOptics;
	// Loop on angles
	Double_t phi = 0;
	Bool_t debug = false;
	Bool_t printout = false;
	Double_t kInnerHits = 0, kOuterHits = 0;
	while(phi < 2*PI){
		// Derived
		Double_t x2=x1+cos(phi), y2=y1+sin(phi);
		if(printout){
			cout << "X1: " << x1 << ", ";
			cout << "Y1: " << y1 << ", ";
			cout << "Phi: " << phi << ": ";
			cout << "X2: " << x2 << ", ";
			cout << "Y2: " << y2 << ", ";
		}
		
		Int_t pointPlacement = getPointPlacement(x1, y1, inner_r, outer_r);
		Double_t close_xf = 0, close_yf = 0;
		Double_t far_xf = 0, far_yf = 0;
		Int_t close_shroud = UNKWN_SHROUD, far_shroud = UNKWN_SHROUD;	// Flag to select the shroud
		if(pointPlacement==INNER_REGION){
			// 2 cases:
			pair<Double_t, Double_t> inner_int, outer_int;
			inner_int = compute_closest_intersection(x1, y1, x2, y2, inner_r);
			outer_int = compute_closest_intersection(x1, y1, x2, y2, outer_r);
			// The closest from inner region is always the inner shroud
			close_xf = inner_int.first;	    
			close_yf = inner_int.second;
			close_shroud = INNER_SHROUD;
			// 1) Inner+Outer: intersect both the shrouds
			if(outer_int.first!=0 || outer_int.second!=0){
				far_xf = outer_int.first; 
				far_yf = outer_int.second;
				far_shroud = OUTER_SHROUD;
			}else{
			// 2) Inner: intersect once the inner shroud
				far_shroud = UNKWN_SHROUD;
			}
		} else if(pointPlacement==OUTER_REGION){
			// 2 cases:
			pair<Double_t, Double_t> inner_int, outer_int, outer_int2;
			inner_int = compute_closest_intersection(x1, y1, x2, y2, inner_r);
			outer_int = compute_closest_intersection(x1, y1, x2, y2, outer_r);
			// The closest from outer region is always the outer shroud
			close_xf = outer_int.first; 
			close_yf = outer_int.second;
			close_shroud = OUTER_SHROUD;
			// 1) Outer+Inner: intersect both the shrouds
			if(inner_int.first!=0 || inner_int.second!=0){
				far_xf = inner_int.first; 
				far_yf = inner_int.second;
				far_shroud = INNER_SHROUD;
			}else{
			// 2) Outer: intersect once the outer shroud
				far_shroud = UNKWN_SHROUD;
			}
		} else {
			// MIDDLE REGION
			pair<Double_t, Double_t> inner_int, outer_int;
			inner_int = compute_closest_intersection(x1, y1, x2, y2, inner_r);
			outer_int = compute_closest_intersection(x1, y1, x2, y2, outer_r);
			// assert: only 1 of the two intersections has the same direction and min distance
			Double_t distance = INF;
			if((inner_int.first!=0 || inner_int.second!=0) && 
			   (x1-x2)*(x1-inner_int.first)>=0 && (y1-y2)*(y1-inner_int.second)>=0){
			    	close_xf = inner_int.first;	    
    				close_yf = inner_int.second;
    				close_shroud = INNER_SHROUD;
				if(printout)
					cout << "HIT INNER ";
				distance = sqrt((x1-close_xf)*(x1-close_xf) + (y1-close_yf)*(y1-close_yf));	
	    		}
			Double_t distance_b = sqrt((x1-outer_int.first)*(x1-outer_int.first) + (y1-outer_int.second)*(y1-outer_int.second));	
			if((outer_int.first!=0 || outer_int.second!=0) && 
			   (x1-x2)*(x1-outer_int.first)>=0 && (y1-y2)*(y1-outer_int.second)>=0 && 
			   (distance_b < distance)){
			    	close_xf = outer_int.first; 
	    			close_yf = outer_int.second;
    				close_shroud = OUTER_SHROUD;
				if(printout)
					cout << "HIT OUTER ";
	    		}
		}
		// Compute the slice of the closest hit
		Double_t WEIGHT_CLOSE_HIT = 0.97462302;	//Norm of: P(capt close fib) + P(pass close fib)*P(hit Ge) = .54 + .46*.9
		Double_t WEIGHT_FAR_HIT= .02484;	//Norm of: P(pass close fib)*P(pass Ge)*P(capt far fib) = .46*.1*.54
		Bool_t close_hit = (close_xf!=0) || (close_yf!=0);
		Bool_t far_hit = (far_xf!=0) || (far_yf!=0);
		assert(!far_hit || close_hit);	// far_hit->close_hit
		if(close_hit==true){
			double angle = atan2(close_yf, close_xf);
			if(close_shroud == INNER_SHROUD){
				innerMap->Fill(x1, angle, WEIGHT_CLOSE_HIT);
				kInnerHits += WEIGHT_CLOSE_HIT;
			}else if(close_shroud == OUTER_SHROUD){
				outerMap->Fill(x1, angle, WEIGHT_CLOSE_HIT);
				kOuterHits += WEIGHT_CLOSE_HIT;
			}else{
				cout << "ERROR: UKNW SHROUD BUT HIT!!!!\n";
				exit(-1);
			}
			if(printout){
				cout << "Close Hit on the fibers at: " << close_xf << ", " << close_yf << ", angle: " << angle << "\n";
			}
		}
		if(far_hit==true){
			double angle = atan2(far_yf, far_xf);
			if(far_shroud == INNER_SHROUD){
				innerMap->Fill(x1, angle, WEIGHT_FAR_HIT);
				kInnerHits += WEIGHT_FAR_HIT;
			}else if(far_shroud == OUTER_SHROUD){
				outerMap->Fill(x1, angle, WEIGHT_FAR_HIT);
				kOuterHits += WEIGHT_FAR_HIT;
			}else{
				cout << "ERROR: UKNW SHROUD BUT HIT!!!!\n";
				exit(-1);
			}
			if(printout){
				cout << "Far Hit on the fibers at: " << far_xf << ", " << far_yf << ", angle: " << angle << "\n";
			}
		}
		if(printout)
			cout << endl;
		phi += delta_phi;
		if(debug)
			break;
	}
	if(kInnerHits + kOuterHits > 0)
		prInnerD->Fill(x1, kInnerHits/(kInnerHits+kOuterHits));
	// return sliced_detections;
}

void createSeparatedSpatialMap(){
	Double_t min_r = 0, max_r = 710, rbins = max_r - min_r + 1;
	Int_t angle_bins = 100;
	Int_t min_angle = -ceil(PI);
	Int_t max_angle = +ceil(PI);
	Int_t nOpticsPerPoint = 5000;
	// Create map
	TString outfile;
    	outfile.Form("ToySpatialMap_%dR_%dAngleSlices_%dops.root", (int)max_r, angle_bins, nOpticsPerPoint);
    	TFile * file = new TFile(outfile, "RECREATE");
	TH1D * prInnerD = new TH1D("PrInnerDet", "Pr ~ Fract. Inner/(Inner+Outer) Detections", rbins, min_r, max_r);
	TH2D * innerMap = new TH2D("InnerMap", "R-Angle Inner Shroud Map", rbins, min_r, max_r, angle_bins, min_angle, max_angle);
	TH2D * outerMap = new TH2D("OuterMap", "R-Angle Outer Shroud Map", rbins, min_r, max_r, angle_bins, min_angle, max_angle);
	for(Double_t x=min_r; x<=max_r; x += 1){
		array<Int_t, n_inner_slices> inner_sliced_detections;
		array<Int_t, n_outer_slices> outer_sliced_detections;
		Int_t innerDetections=0, totDetections=0;
		runToyOpticsFromPoint(x, 0., nOpticsPerPoint, prInnerD, innerMap, outerMap);
	}
    	prInnerD->Write();
    	innerMap->Write();
    	outerMap->Write();
    	file->Close();
}
