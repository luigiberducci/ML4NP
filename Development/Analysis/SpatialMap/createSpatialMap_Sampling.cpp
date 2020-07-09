#include <iostream>
#include <TEntryList.h>
#include <TParameter.h>
#include <TFile.h>
#include <TString.h>
#include <TTree.h>
#include <TBranch.h>
#include <TChain.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TRandom.h>
#include <array>
#include <vector>
#include <set>
#include <assert.h>

#define INNER_REGION -1
#define MIDDLE_REGION 0
#define OUTER_REGION +1
#define INNER_SHROUD 0
#define OUTER_SHROUD 1
#define UNKWN_SHROUD -1

#define TOPZ 	   +845
#define BOTTOMZ    -845
#define TOPZ_GE    +425
#define BOTTOMZ_GE -425
#define LEFT_GE    -(235+40)
#define RIGHT_GE   235+40

#define NONCOORD -666666666

using namespace std;

Double_t EPS = 0.00001;
Double_t INF = 1000000;
Double_t PI = 3.141592653589793;
// Rnd Generator
TRandom * rnd = new TRandom();
// Coordinate Ge Crystals
Int_t nGeCrystals = 14;
Double_t radiusGeCrystals = 235.0;
Double_t maxGeRotationAngle = 0.22;	// the rot angle ranges 0,.22
Double_t geRadius = 40.0;
std::vector<Double_t> geCenters_x, geCenters_y;

pair<Double_t, Double_t> compute_closest_intersection(Double_t x1, Double_t y1, Double_t x2, Double_t y2, Double_t r, Double_t maxDistance=INF){
	if(sqrt(x1*x1 + y1*y1) == r)	// check if point stay on circle boundary
		return make_pair(x1, y1);
	// Compute line parameters
	Double_t a = y1 - y2;
	Double_t b = x2 - x1;
	Double_t c = x1 * y2 - x2 * y1;
	// Compute the point in line closest to the origin
	double xf = -a*c/(a*a+b*b), yf = -b*c/(a*a+b*b);
	// Based on dist from origin and radius, we know multiplicity
	if (c*c > r*r*(a*a+b*b)+EPS){
	    xf = NONCOORD;
	    yf = NONCOORD;
	} else if (abs (c*c - r*r*(a*a+b*b)) >= EPS) {
	    double d = r*r - c*c/(a*a+b*b);
	    double mult = sqrt (d / (a*a+b*b));
	    double ax, ay, bx, by;
	    ax = xf + b * mult;
	    bx = xf - b * mult;
	    ay = yf - a * mult;
	    by = yf + a * mult;
	    // Consider only interesection within maxDistance
	    double distance_a = sqrt((ax-x1)*(ax-x1) + (ay-y1)*(ay-y1));
	    double distance_b = sqrt((bx-x1)*(bx-x1) + (by-y1)*(by-y1));
	    // assert: a valid intersection has the same direction(both x, y)
	    // if 2 intersections can be 1) only 1 in the same direction, 2) both of them (line from outside)
	    double distance = INF;
	    if((x1-x2)*(x1-ax)>=0 && (y1-y2)*(y1-ay)>=0 && (distance_a<=maxDistance)){
	    	xf = ax;	    
	    	yf = ay;
	        distance = sqrt((x1-ax)*(x1-ax) + (y1-ay)*(y1-ay));	
	    }
	    if((x1-x2)*(x1-bx)>=0 && (y1-y2)*(y1-by)>=0 && (distance_b<=maxDistance) && (distance_b<distance)){
			xf = bx;	    
			yf = by;
			distance = distance_b;
	    }
	    if(distance == INF){
	    	xf = NONCOORD;		// 2 intersection but none of them has the right direction
	    	yf = NONCOORD;
	    }
	}else{ // else only 1 (xf, yf)
	    double distance = sqrt((xf-x1)*(xf-x1) + (yf-y1)*(yf-y1));
	    if(distance > maxDistance){	// the only intersection if farther than limit
	    	xf = NONCOORD;
	    	yf = NONCOORD;
	    }
	}
	pair<Double_t, Double_t> intersection(make_pair(xf, yf));
	return intersection;
}

pair<Double_t, Double_t> compute_farthest_intersection(Double_t x1, Double_t y1, Double_t x2, Double_t y2, Double_t r, Double_t maxDistance=INF){
	// Compute line parameters
	Double_t a = y1 - y2;
	Double_t b = x2 - x1;
	Double_t c = x1 * y2 - x2 * y1;
	// Compute the point in line closest to the origin
	double xf = -a*c/(a*a+b*b), yf = -b*c/(a*a+b*b);
	// Based on dist from origin and radius, we know multiplicity
	if (c*c > r*r*(a*a+b*b)+EPS){
	    xf = NONCOORD;
	    yf = NONCOORD;
	} else if (abs (c*c - r*r*(a*a+b*b)) >= EPS) {
	    double d = r*r - c*c/(a*a+b*b);
	    double mult = sqrt (d / (a*a+b*b));
	    double ax, ay, bx, by;
	    ax = xf + b * mult;
	    bx = xf - b * mult;
	    ay = yf - a * mult;
	    by = yf + a * mult;
	    // Maintain only point close enough
	    double distance_a = sqrt((x1-ax)*(x1-ax) + (y1-ay)*(y1-ay));
	    double distance_b = sqrt((x1-bx)*(x1-bx) + (y1-by)*(y1-by));
	    // assert: a valid intersection has the same direction(both x, y)
	    // if 2 intersections can be 1) only 1 in the same direction, 2) both of them (line from outside)
	    double distance = -INF;
	    if((x1-x2)*(x1-ax)>=0 && (y1-y2)*(y1-ay)>=0 && distance_a<=maxDistance){
	    	xf = ax;	    
	    	yf = ay;
	        distance = sqrt((x1-ax)*(x1-ax) + (y1-ay)*(y1-ay));	
	    }
	    if((x1-x2)*(x1-bx)>=0 && (y1-y2)*(y1-by)>=0 && distance_b<=maxDistance && distance_b>distance){
			xf = bx;	    
			yf = by;
			distance = distance_b;
	    }
	    if(distance == -INF){
	    	xf = NONCOORD;		// 2 intersection but none of them has the right direction
	    	yf = NONCOORD;
	    }
	}else{ // else only 1 (xf, yf)
	    double distance = sqrt((xf-x1)*(xf-x1) + (yf-y1)*(yf-y1));
	    if(distance > maxDistance){	// the only intersection if farther than limit
	    	xf = NONCOORD;
	    	yf = NONCOORD;
	    }
	}
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

Int_t initializePositionGeCrystals(){
	// Note: rotation_angle has been computed to disalign the Ge crystals from the x-axis
	// 0.22 is computed for 14 crystals of radius 40
	Double_t geAngle = 0;
	cout << "[Info] Initialization Ge Crystals...\n";
	for(Int_t ige=0; ige<nGeCrystals; ige++){
		// Compute center of Ge crystal
		Double_t x0 = radiusGeCrystals * cos(geAngle);
		Double_t y0 = radiusGeCrystals * sin(geAngle);
		// Rotate it
		cout << "\tCrystal " << ige+1 << "(Initially not shifted): ";
		cout << "Centered in " << x0 << ", " << y0 << " => R: " << sqrt(x0*x0+y0*y0) << endl;
		// Append coordinates
		geCenters_x.push_back(x0);
		geCenters_y.push_back(y0);
		geAngle += 2*PI / nGeCrystals;
	}
}

pair<Double_t, Double_t> getHitOnGermanium(Double_t x1, Double_t y1, Double_t x2, Double_t y2, Bool_t enableGe){
	// Idea: the method to find an intersection  assumes the circle to be centered in 0,0
	// We can translate the points to have each Ge crystal centered
	Double_t distance = INF;
	pair<Double_t, Double_t> closest_ge_int = make_pair(NONCOORD, NONCOORD);
	if(!enableGe)
		return closest_ge_int;	// If Ge disable, return always NoHit!
	// We aim to average the result over a slice, then sample the rotation angle
	for(int i=0; i<geCenters_x.size(); i++){
		// Shift the point coordinates w.r.t. the center of crystal
		Double_t xx1 = x1 - geCenters_x[i];
		Double_t yy1 = y1 - geCenters_y[i];
		Double_t xx2 = x2 - geCenters_x[i];
		Double_t yy2 = y2 - geCenters_y[i];
		// Once Ge has been centered, check hit
		pair<Double_t, Double_t> ge_int = compute_closest_intersection(xx1, yy1, xx2, yy2, geRadius);
		if(ge_int.first!=NONCOORD || ge_int.second!=NONCOORD){
			ge_int.first += geCenters_x[i];
			ge_int.second += geCenters_y[i];
			Double_t new_distance = sqrt((x1-ge_int.first)*(x1-ge_int.first) + (y1-ge_int.second)*(y1-ge_int.second));
			if(new_distance < distance){
				distance = new_distance;
				closest_ge_int.first = ge_int.first;
				closest_ge_int.second = ge_int.second;
			}	
		}
	}
	return closest_ge_int;
}

Double_t getAngleBtw3Points(Double_t cx, Double_t cy, Double_t ax, Double_t ay, Double_t bx, Double_t by){
	// Compute the angle between point A, B w.r.t. C (center)
	Double_t angle_a = atan2(ay - cy, ax - cx);
	if(angle_a < 0)
		angle_a += 2 * PI;
	Double_t angle_b = atan2(by - cy, bx - cx);
	if(angle_b < 0)
		angle_b += 2 * PI;
	return angle_b - angle_a;
}

Double_t computeProbTrajectoryToGermanium(Double_t radius, Double_t z){
	//Idea: given a position in r, z return the ratio of angle hitting the Ge over the angle hitting the shroud
	// Points: A highest intersection with ROI
	// 	   B highest intersection Ge region 
	// 	   C lowest intersection Ge region
	// 	   D lowest intersection with ROI
	Double_t bx, cx;
	Double_t ax = RIGHT_GE, ay = TOPZ, by = TOPZ_GE;
	Double_t cy = BOTTOMZ_GE, dx = RIGHT_GE, dy = BOTTOMZ;
	if(z > TOPZ_GE){
		bx = LEFT_GE;
		cx = RIGHT_GE;
	}else if(z > BOTTOMZ_GE){
		bx = RIGHT_GE;
		cx = RIGHT_GE;
	}else{
		bx = RIGHT_GE;
		cx = LEFT_GE;
	}
	Double_t angle_ad = getAngleBtw3Points(radius, z, ax, ay, dx, dy);
	Double_t angle_bc = getAngleBtw3Points(radius, z, bx, by, cx, cy);
	return angle_bc/angle_ad;
}

Bool_t checkIfEnableGe(Double_t x1, Double_t y1, Double_t z1, Double_t x2, Double_t y2, Double_t z2){
	// Compute intersection with the Ge volume, meaning the cylinder that includes the Ge Strings
	Double_t ge_r = 235 + 40;
	pair<Double_t, Double_t> ge_int = compute_closest_intersection(x1, y1, x2, y2, ge_r);
	if(ge_int.first==NONCOORD || ge_int.second==NONCOORD)
		return false;	// No interesection at all, no difference between the two maps
	// If there is an interesection, look at where it occurs in Z
	// Use parametrization of line: line = P + t D, where P is a vector (point), D direction vector
	Double_t t = (ge_int.first - x1) / (x2-x1);	// From param: t = (x-x1)/dx
	Double_t z_int = z1 + t * (z2-z1);   		// z = z1 + t * dz, where z1 point, dz direction vector
	if(BOTTOMZ_GE < z_int && z_int < TOPZ_GE)
		return true;	// Interesection with Z in [BOTTOMGE,TOPGE], enable crystals
	return false;		// Intersection with Z outside the Ge volume, disable crystals
}

// Parameters
const Int_t n_inner_slices = 100;
const Int_t n_outer_slices = 100;
const Double_t inner_theta_slice = 2 * PI / n_inner_slices;    // angle for each slice
const Double_t outer_theta_slice = 2 * PI / n_outer_slices;    // angle for each slice
void runToyOpticsFromPoint(Double_t radius, Int_t nOptics, TH1D * &prInnerD, TH2D* &innerMap, TH2D* &outerMap){
	cout << "Running at R=" << radius << endl;
	// Fiber radius
	Double_t inner_r = 175, outer_r = 295;	// mm
	// Flag sample Ge
	Bool_t enableGe = true;
	// Loop on angles
	Double_t phi, theta;
	Bool_t debug = false;
	Bool_t printout = false;
	Double_t kInnerHits = 0, kOuterHits = 0;
	// Get coordinate
	Double_t x1 = radius, y1 = 0, z1 = 0;
	while(nOptics>0){
		// Get random point (y1 according to angle shifting, z1 random) 
		Double_t phi1 = -maxGeRotationAngle + rnd->Rndm() * 2 * maxGeRotationAngle;
		y1 = x1 * sin(phi1);
		z1 = BOTTOMZ + rnd->Rndm() * TOPZ;
		// Compute Pr of Ge Region
		Double_t prTrajectoryToGe = computeProbTrajectoryToGermanium(x1, z1);
		// Get random angle
		theta = rnd->Rndm() * PI;
		phi = rnd->Rndm() * 2 * PI;
		// Derived
		Double_t x2 = x1 + 1 * sin(theta) * cos(phi);
		Double_t y2 = y1 + 1 * sin(theta) * sin(phi);
		Double_t z2 = z1 + 1 * cos(theta);
		Double_t lenTrajectory = 0;	// distance to hit the plane delimited by upper or lower boundary
		if(theta <= PI/2){
			lenTrajectory = (TOPZ-z1)/cos(theta);	    // trajectory towards top Z
		}else{
			lenTrajectory = (BOTTOMZ-z1)/cos(theta);    // trajectory towards bottom Z
		}
		enableGe = checkIfEnableGe(x1, y1, z1, x2, y2, z2);
		// Note: Ignore z
		if(printout){
			cout << "X1: " << x1 << ", ";
			cout << "Y1: " << y1 << ", ";
			cout << "Phi: " << phi << ": ";
			cout << "X2: " << x2 << ", ";
			cout << "Y2: " << y2 << ", ";
			cout << "Hitting Ge: " << enableGe << ", ";
			cout << "Max Len Traj: " << lenTrajectory << "\n";
		}
		
		Int_t pointPlacement = getPointPlacement(x1, y1, inner_r, outer_r);
		Double_t close_xf = 0, close_yf = 0;
		Double_t far_xf = 0, far_yf = 0;
		Double_t farfar_xf = 0, farfar_yf = 0;
		Int_t close_shroud = UNKWN_SHROUD, far_shroud = UNKWN_SHROUD, farfar_shroud = UNKWN_SHROUD;	// Flag to select the shroud
		if(pointPlacement==INNER_REGION){
			// 2 cases:
			pair<Double_t, Double_t> inner_int, outer_int, ge_int;
			inner_int = compute_closest_intersection(x1, y1, x2, y2, inner_r, lenTrajectory);
			ge_int = getHitOnGermanium(x1, y1, x2, y2, enableGe);
			outer_int = compute_closest_intersection(x1, y1, x2, y2, outer_r, lenTrajectory);
			// The closest from inner region is always the inner shroud
			// Not hit on inner shroud is not possible
			if(inner_int.first!=NONCOORD || inner_int.second!=NONCOORD){
				close_xf = inner_int.first;	    
				close_yf = inner_int.second;
				close_shroud = INNER_SHROUD;
				if(printout)
					cout << "HIT INNER ";
				// 1) Inner+Outer: intersect both the shrouds, no Ge hit in the middle
				if((outer_int.first!=NONCOORD || outer_int.second!=NONCOORD) ){
					if(ge_int.first==NONCOORD && ge_int.second==NONCOORD){
						far_xf = outer_int.first; 
						far_yf = outer_int.second;
						far_shroud = OUTER_SHROUD;
						if(printout)
							cout << "2 HIT OUTER ";
					}
				}else{
				// 2) Inner: intersect once the inner shroud (the outer is not reached or hit Ge before)
					far_shroud = UNKWN_SHROUD;
				}
			}
		} else if(pointPlacement==OUTER_REGION){
			// assert: if hit shroud, always the outer first
			// assert: if 2nd hit, check if no Ge before
			// 2 cases:
			pair<Double_t, Double_t> inner_int, inner_int2, outer_int, outer_int2, ge_int;
			inner_int = compute_closest_intersection(x1, y1, x2, y2, inner_r, lenTrajectory);
			ge_int = getHitOnGermanium(x1, y1, x2, y2, enableGe);
			outer_int = compute_closest_intersection(x1, y1, x2, y2, outer_r, lenTrajectory);
			outer_int2 = compute_farthest_intersection(x1, y1, x2, y2, outer_r, lenTrajectory);
			// The closest from outer region is always the outer shroud
			if(outer_int.first!=NONCOORD || outer_int.second!=NONCOORD){	// Discard no hit on outer shroud
				close_xf = outer_int.first; 
				close_yf = outer_int.second;
				close_shroud = OUTER_SHROUD;
				if(printout)
					cout << "1 HIT OUTER ";
				// The 2nd hit is the one at min dist from the point
				Double_t distance = INF;
				if((outer_int2.first!=NONCOORD || outer_int2.second!=NONCOORD) && 
				   (x1-x2)*(x1-outer_int2.first)>=0 && (y1-y2)*(y1-outer_int2.second)>=0){
					far_xf = outer_int2.first;	    
					far_yf = outer_int2.second;
					far_shroud = OUTER_SHROUD;
					if(printout)
						cout << "2 HIT OUTER ";
					distance = sqrt((x1-far_xf)*(x1-far_xf) + (y1-far_xf)*(y1-far_xf));
				}
				Double_t distance_in = sqrt((x1-inner_int.first)*(x1-inner_int.first) + (y1-inner_int.second)*(y1-inner_int.second));	
				if((inner_int.first!=NONCOORD || inner_int.second!=NONCOORD) && 
				   (x1-x2)*(x1-inner_int.first)>=0 && (y1-y2)*(y1-inner_int.second)>=0 && 
				   (distance_in < distance)){
					far_xf = inner_int.first; 
					far_yf = inner_int.second;
					far_shroud = INNER_SHROUD;
					if(printout)
						cout << "2 HIT INNER ";
					// Check 2nd hit on the inner shroud
					inner_int2 = compute_farthest_intersection(x1, y1, x2, y2, inner_r, lenTrajectory);
					if((inner_int2.first!=NONCOORD || inner_int2.second!=NONCOORD) &&
					   (x1-x2)*(x1-inner_int2.first)>=0 && (y1-y2)*(y1-inner_int2.second)>=0){
						farfar_xf = inner_int2.first; 
						farfar_yf = inner_int2.second;
						farfar_shroud = INNER_SHROUD;
					}
					distance = distance_in;
				}
				Double_t distance_ge = sqrt((x1-ge_int.first)*(x1-ge_int.first) + (y1-ge_int.second)*(y1-ge_int.second));
				if((ge_int.first!=NONCOORD || ge_int.second!=NONCOORD) && 
				   (x1-x2)*(x1-ge_int.first)>=0 && (y1-y2)*(y1-ge_int.second)>=0 && 
				   (distance_ge < distance)){
					// if hit ge, always after the outer shroud, then keep close_shroud
					far_shroud = UNKWN_SHROUD;
					farfar_shroud = UNKWN_SHROUD;
					if(printout)
						cout << "HIT GE ";
				}
			}
		} else {
			// MIDDLE REGION
			pair<Double_t, Double_t> inner_int, outer_int, ge_int;
			inner_int = compute_closest_intersection(x1, y1, x2, y2, inner_r, lenTrajectory);
			ge_int = getHitOnGermanium(x1, y1, x2, y2, enableGe);
			outer_int = compute_closest_intersection(x1, y1, x2, y2, outer_r, lenTrajectory);
			// assert: only 1 of the intersections (inner, outer, ge) has the same direction and min distance
			Double_t distance = INF;
			// assert: only 1 eventual hit in the outer shroud
			if((outer_int.first!=NONCOORD || outer_int.second!=NONCOORD) && 
			   (x1-x2)*(x1-outer_int.first)>=0 && (y1-y2)*(y1-outer_int.second)>=0){
			    	close_xf = outer_int.first;	    
    				close_yf = outer_int.second;
    				close_shroud = OUTER_SHROUD;
				if(printout)
					cout << "HIT OUTER ";
				distance = sqrt((x1-close_xf)*(x1-close_xf) + (y1-close_yf)*(y1-close_yf));	
	    		}
			Double_t distance_in = sqrt((x1-inner_int.first)*(x1-inner_int.first) + (y1-inner_int.second)*(y1-inner_int.second));	
			// assert: at most 2 hits in the inner shroud
			if((inner_int.first!=NONCOORD || inner_int.second!=NONCOORD) && 
			   (x1-x2)*(x1-inner_int.first)>=0 && (y1-y2)*(y1-inner_int.second)>=0 && 
			   (distance_in < distance)){
			    	close_xf = inner_int.first; 
	    			close_yf = inner_int.second;
    				close_shroud = INNER_SHROUD;
				if(printout)
					cout << "HIT INNER";
				// Check 2nd hit on the inner shroud
				inner_int = compute_farthest_intersection(x1, y1, x2, y2, inner_r, lenTrajectory);
				if((inner_int.first!=NONCOORD || inner_int.second!=NONCOORD) &&
				   (x1-x2)*(x1-inner_int.first)>=0 && (y1-y2)*(y1-inner_int.second)>=0){
					far_xf = inner_int.first; 
					far_yf = inner_int.second;
					far_shroud = INNER_SHROUD;
				}
				distance = distance_in;
	    		}
			Double_t distance_ge = sqrt((x1-ge_int.first)*(x1-ge_int.first) + (y1-ge_int.second)*(y1-ge_int.second));
			if((ge_int.first!=NONCOORD || ge_int.second!=NONCOORD) && 
			   (x1-x2)*(x1-ge_int.first)>=0 && (y1-y2)*(y1-ge_int.second)>=0 && 
			   (distance_ge < distance)){
    				close_shroud = UNKWN_SHROUD;	// reset eventual other hits recorded
    				far_shroud = UNKWN_SHROUD;
				if(printout)
					cout << "HIT GE";
	    		}
		}
		Bool_t close_hit = close_shroud != UNKWN_SHROUD;
		Bool_t far_hit = far_shroud != UNKWN_SHROUD;
		Bool_t farfar_hit = farfar_shroud != UNKWN_SHROUD;
		assert(!far_hit || close_hit);	// far_hit->close_hit
		assert(!farfar_hit || far_hit);	// farfar_hit->far_hit
		// Compute probabilities of multiple hits
		Double_t NORM_FACTOR;
		if(farfar_hit==true)
			NORM_FACTOR = .54 + .46*.54 + .46*.46*.54;
		else if(far_hit==true)
			NORM_FACTOR = .54 + .46*.54;
		else
			NORM_FACTOR = .54;
		Double_t WEIGHT_CLOSE_HIT = .54/NORM_FACTOR;	//P(captured closest shroud)
		Double_t WEIGHT_FAR_HIT= .46*.54/NORM_FACTOR;	//P(not capt closest shroud) P(captured second shroud)
		Double_t WEIGHT_FARFAR_HIT= .46*.46*.54/NORM_FACTOR;	//P(not capt closest) P(not capt second) P(captured third shroud)
		// Sampling the hit
		Double_t rndChoice = rnd->Rndm();
		Bool_t fillCloseHit = false, fillFarHit = false, fillFarFarHit = false;
		if(farfar_hit==true){
			if(rndChoice > 1 - WEIGHT_FARFAR_HIT)
				fillFarFarHit = true;
			else if(rndChoice > 1 - (WEIGHT_FAR_HIT+WEIGHT_FARFAR_HIT))
				fillFarHit = true;
			else
				fillCloseHit = true;
		}else if(far_hit==true){
			if(rndChoice > 1 - WEIGHT_FAR_HIT)
				fillFarHit = true;
			else
				fillCloseHit = true;
		}else if(close_hit==true){
			fillCloseHit = true;
		}
		// Only one must be true
		assert(!fillCloseHit || (!fillFarHit && !fillFarFarHit));
		assert(!fillFarHit || (!fillCloseHit && !fillFarFarHit));
		assert(!fillFarFarHit || (!fillFarHit && !fillCloseHit));
		if(fillCloseHit==true){
			assert(close_xf!=NONCOORD || close_yf!=NONCOORD);
			double angle = atan2(close_yf, close_xf);
			if(close_shroud == INNER_SHROUD){
				innerMap->Fill(x1, angle);
				kInnerHits += 1;
			}else if(close_shroud == OUTER_SHROUD){
				outerMap->Fill(x1, angle);
				kOuterHits += 1;
			}
			if(printout){
				cout << "Close Hit on the fibers at: " << close_xf << ", " << close_yf << ", angle: " << angle << "\n";
			}
		}
		if(fillFarHit==true){
			assert(far_xf!=NONCOORD || far_yf!=NONCOORD);
			double angle = atan2(far_yf, far_xf);
			if(far_shroud == INNER_SHROUD){
				innerMap->Fill(x1, angle);
				kInnerHits += 1;
			}else if(far_shroud == OUTER_SHROUD){
				outerMap->Fill(x1, angle);
				kOuterHits += 1;
			}
			if(printout){
				cout << "Far Hit on the fibers at: " << far_xf << ", " << far_yf << ", angle: " << angle << "\n";
			}
		}
		if(fillFarFarHit==true){
			assert(farfar_xf!=NONCOORD || farfar_yf!=NONCOORD);
			double angle = atan2(farfar_yf, farfar_xf);
			if(farfar_shroud == INNER_SHROUD){
				innerMap->Fill(x1, angle);
				kInnerHits += 1;
			}else if(farfar_shroud == OUTER_SHROUD){
				outerMap->Fill(x1, angle);
				kOuterHits += 1;
			}
			if(printout){
				cout << "Far Far Hit on the fibers at: " << farfar_xf << ", " << farfar_yf << ", angle: " << angle << "\n";
			}
		}
		if(printout)
			cout << endl;
		if(close_hit){
			nOptics--;	// decrement n optics only when hit a shroud
		}
		if(debug)
			break;
	}
	if(kInnerHits + kOuterHits > 0)
		prInnerD->Fill(x1, kInnerHits/(kInnerHits+kOuterHits));
	// return sliced_detections;
}

void createSpatialMap_Sampling(){
	Double_t min_r = 0, max_r = 710, rbins = max_r - min_r + 1;
	Int_t angle_bins = 100;
	Int_t min_angle = -ceil(PI);
	Int_t max_angle = +ceil(PI);
	Int_t nOpticsPerPoint = 1000000;
	// Initialize Ge Crystals
	initializePositionGeCrystals();
	// Create map
	TString outfile;
    	outfile.Form("ToySpatialMap_%dR_%dAngleSlices_%dops_ShortHitSampling.root", (int)rbins, angle_bins, nOpticsPerPoint);
    	TFile * file = new TFile(outfile, "RECREATE");
	TH1D * prInnerD = new TH1D("PrInnerDet", "Pr ~ Fract. Inner/(Inner+Outer) Detections", rbins, min_r, max_r);
	TH2D * innerMap = new TH2D("InnerMap", "R-Angle Inner Shroud Map", rbins, min_r, max_r, angle_bins, min_angle, max_angle);
	TH2D * outerMap = new TH2D("OuterMap", "R-Angle Outer Shroud Map", rbins, min_r, max_r, angle_bins, min_angle, max_angle);
	for(Double_t x=min_r; x<=max_r; x += 1){
		array<Int_t, n_inner_slices> inner_sliced_detections;
		array<Int_t, n_outer_slices> outer_sliced_detections;
		Int_t innerDetections=0, totDetections=0;
		runToyOpticsFromPoint(x, nOpticsPerPoint, prInnerD, innerMap, outerMap);
	}
    	prInnerD->Write();
    	innerMap->Write();
    	outerMap->Write();
    	file->Close();
}

int main(){
	createSpatialMap_Sampling();
}
