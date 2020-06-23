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

const Int_t n_slices = 72;
//array<Double_t, n_slices> runToyOpticsFromPoint(Double_t x1, Double_t y1, Int_t nOptics=360){
void runToyOpticsFromPoint(Double_t x1, Double_t y1, Int_t nOptics, array<Int_t, n_slices> &sliced_detections){
	// Parameters
	Double_t theta_slice = 2 * PI / n_slices;    // angle for each slice
	Int_t sum_detections = 0;     // For normalization
	for(int i=0; i<n_slices; i++)
		sliced_detections[i]=0;
	// Fiber radius
	Double_t inner_r = 175, outer_r = 295;	// mm
	Double_t delta_phi = PI / nOptics;
	// Loop on angles
	Double_t phi = 0;
	Bool_t debug = false;
	Bool_t printout = false;
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
		Double_t xf = 0, yf = 0;
		if(pointPlacement==INNER_REGION){
			pair<Double_t, Double_t> inner_int;
			inner_int = compute_intersection_slice(x1, y1, x2, y2, inner_r);
			xf = inner_int.first;	    
			yf = inner_int.second;
		} else if(pointPlacement==OUTER_REGION){
			pair<Double_t, Double_t> outer_int;
			outer_int = compute_intersection_slice(x1, y1, x2, y2, outer_r);
			xf = outer_int.first; 
			yf = outer_int.second;
		} else {
			pair<Double_t, Double_t> inner_int, outer_int;
			inner_int = compute_intersection_slice(x1, y1, x2, y2, inner_r);
			outer_int = compute_intersection_slice(x1, y1, x2, y2, outer_r);
			// assert: only 1 of the two intersections has the same direction
			Double_t distance = INF;
			if((inner_int.first!=0 || inner_int.second!=0) && 
			   (x1-x2)*(x1-inner_int.first)>=0 && (y1-y2)*(y1-inner_int.second)>=0){
			    	xf = inner_int.first;	    
	    			yf = inner_int.second;
				if(printout)
					cout << "HIT INNER ";
				distance = sqrt((x1-xf)*(x1-xf) + (y1-yf)*(y1-yf));	
	    		}
			Double_t distance_b = sqrt((x1-outer_int.first)*(x1-outer_int.first) + (y1-outer_int.second)*(y1-outer_int.second));	
			if((outer_int.first!=0 || outer_int.second!=0) && 
			   (x1-x2)*(x1-outer_int.first)>=0 && (y1-y2)*(y1-outer_int.second)>=0 && 
			   (distance_b < distance)){
			    	xf = outer_int.first; 
	    			yf = outer_int.second;
				if(printout)
					cout << "HIT OUTER ";
	    		}
		}
		// Compute the slice
		Bool_t hit = (xf!=0) || (yf!=0);
		if(hit==true){
			double angle = atan2(yf, xf);
			if(angle < 0)
				angle += 2 * PI;
			int slice_id = (angle / theta_slice);
			sliced_detections[slice_id]++;
			sum_detections++;
			if(printout){
				cout << "Hit on the fibers at: " << xf << ", " << yf << ", angle: " << angle << "("<< slice_id << ")";
			}
		}
		if(printout)
			cout << endl;
		phi += delta_phi;
		if(debug)
			break;
	}
	// return sliced_detections;
}

void createToySpatialMap(){
	Double_t min_r = 0, max_r = 700, rbins = max_r - min_r + 1;
	Int_t slice_bins = n_slices;
	Int_t min_scaled_s  = -(n_slices / 2);		// in the map, slice id for 72 slices is in [-36, +35]
	Int_t max_scaled_s = (n_slices / 2) - 1;
	Int_t nOpticsPerPoint = 1000;
	// Create map
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
				map->Fill(x, scaled_s);
			}
		}
		cout << endl;
	}
}
