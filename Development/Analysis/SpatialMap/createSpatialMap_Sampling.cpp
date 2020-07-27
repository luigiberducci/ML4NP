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
#include <TMath.h>
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
#define GE_VOLUME -2

#define TOPZ 	   +845
#define BOTTOMZ    -845
#define TOPZ_GE    +425
#define BOTTOMZ_GE -425
#define LEFT_GE    -(235+40)
#define RIGHT_GE   235+40

#define INNER_SHROUD_RADIUS 175
#define OUTER_SHROUD_RADIUS 295
#define FIBER_THICKNESS 0.06

#define NONCOORD -666666666
#define NONDIST 666666666

using namespace std;

Double_t EPS = 0.00001;
Double_t INF = 1000000;
Double_t PI = 3.141592653589793;
// Rnd Generator
TRandom * rnd = new TRandom();
// Attenuation Lenght, Geometric Coverage
Double_t attenuationLen = 500;     //mm Attenuation Leng
Double_t shroudCapturePr = .54;    //fiber coverage over real cylinder surface
// Coordinate Ge Crystals
Int_t nGeCrystals = 14;
Double_t radiusGeCrystals = 235.0;
Double_t shiftAngleForXYSampling = 0.22;	// the rot angle ranges 0,.22
Double_t geRadius = 40.0;
std::vector<Double_t> geCenters_x, geCenters_y;

// Struct 3d point
struct Point{
	Double_t x;
	Double_t y;
	Double_t z;
    Int_t shroud;       // Optional
    Double_t distance;     // Optional
	Point(Double_t xx=NONCOORD, Double_t yy=NONCOORD, Double_t zz=NONCOORD,
	      Int_t shrd=UNKWN_SHROUD, Double_t dist=NONDIST){
		x = xx;
		y = yy;
		z = zz;
		shroud = shrd;
		distance = dist;
	}
	void reset(){
		x = NONCOORD;
		y = NONCOORD;
		z = NONCOORD;
		shroud = UNKWN_SHROUD;
		distance = NONDIST;
	}
	bool isDefined(){
		return (x!=NONCOORD || y!=NONCOORD || z!=NONCOORD);
	}
	bool checkSameDirection(Double_t x2, Double_t y2, struct Point p){
		return (x-x2)*(x-p.x)>=0 && (y-y2)*(y-p.y)>=0;
	}
    bool operator < (const Point& point){
        return distance < point.distance;
    }
    bool operator == (const Point& point){
        return abs(x - point.x)<=EPS || abs(y - point.y)<=EPS || abs(z - point.z)<=EPS;
    }
};

std::ostream& operator<<(std::ostream& os, const Point& point){
    os << "X: " << point.x << ", ";
    os << "Y: " << point.y << ", ";
    os << "Z: " << point.z << " | ";
    os << "Shroud: " << point.shroud << ", ";
    os << "Dist_prod: " << point.distance;
    return os;
}

Double_t getPointDistance(struct Point a, struct Point b){
    return sqrt((a.x-b.x)*(a.x-b.x) + (a.y-b.y)*(a.y-b.y) + (a.z-b.z)*(a.z-b.z));
}

vector<Point> compute_intersection(Point prodPoint, Point dirPoint, Double_t r, Double_t maxDistance=INF) {
    vector<Point> hits;
    if (abs(sqrt(prodPoint.x * prodPoint.x + prodPoint.y * prodPoint.y) - r) > FIBER_THICKNESS/2) {    // check if point stay on circle boundary
        // Compute line parameters
        Double_t a = prodPoint.y - dirPoint.y;
        Double_t b = dirPoint.x - prodPoint.x;
        Double_t c = prodPoint.x * dirPoint.y - dirPoint.x * prodPoint.y;
        // Compute the point in line closest to the origin
        Point hit = Point();
        hit.x = -a * c / (a * a + b * b), hit.y = -b * c / (a * a + b * b);
        // Based on dist from origin and radius, we know multiplicity
        if (c * c > r * r * (a * a + b * b) + EPS) {
            hit.reset();
        } else if (abs(c * c - r * r * (a * a + b * b)) >= EPS) {
            double d = r * r - c * c / (a * a + b * b);
            double mult = sqrt(d / (a * a + b * b));
            double ax, ay, bx, by;
            ax = hit.x + b * mult;
            bx = hit.x - b * mult;
            ay = hit.y - a * mult;
            by = hit.y + a * mult;
            // Consider only interesection within maxDistance
            double distance_a = sqrt((ax - prodPoint.x) * (ax - prodPoint.x) + (ay - prodPoint.y) * (ay - prodPoint.y));
            double distance_b = sqrt((bx - prodPoint.x) * (bx - prodPoint.x) + (by - prodPoint.y) * (by - prodPoint.y));
            // assert: a valid intersection has the same direction(both x, y)
            // if 2 intersections can be 1) only 1 in the same direction, 2) both of them (line from outside)
            Point hit = Point(ax, ay);
            if (prodPoint.checkSameDirection(dirPoint.x, dirPoint.y, hit) && (distance_a <= maxDistance))
                hits.push_back(hit);
            hit = Point(bx, by);
            if (prodPoint.checkSameDirection(dirPoint.x, dirPoint.y, hit) && (distance_b <= maxDistance))
                hits.push_back(hit);
        } else { // else only 1 (xf, yf)
            double distance = sqrt((hit.x - prodPoint.x) * (hit.x - prodPoint.x) + (hit.y - prodPoint.y) * (hit.y - prodPoint.y));
            if (prodPoint.checkSameDirection(dirPoint.x, dirPoint.y, hit) && (distance <= maxDistance))     // the only intersection if farther than limit
                hits.push_back(hit);
        }
    } else {      // point on the circle boundary
        hits.push_back(prodPoint);
    }
    for (auto &hit: hits) {
        if (hit.isDefined()) {
            // If there is an interesection, look at where it occurs in Z
            // Use parametrization of line: line = P + t D, where P is a vector (point), D direction vector
            Double_t t = (hit.x - prodPoint.x) / (dirPoint.x - prodPoint.x);     // From param: t = (x-x1)/dx
            hit.z = prodPoint.z + t * (dirPoint.z - prodPoint.z);              // z = z1 + t * dz, where z1 point, dz direction vector
            hit.distance = getPointDistance(prodPoint, hit);
        }
    }
    std::sort(hits.begin(), hits.end());    // sort ascending distance
    return hits;
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

struct Point getHitOnGermanium(Point point, Point dirPoint, Bool_t enableGe){
	// Idea: the method to find an intersection  assumes the circle to be centered in 0,0
	// We can translate the points to have each Ge crystal centered
	Double_t distance = INF;
	struct Point closestGeHit = Point();
	if(!enableGe)
		return closestGeHit;	// If Ge disable, return always NoHit!
	// We aim to average the result over a slice, then sample the rotation angle
    Double_t original_x1 = point.x, original_y1 = point.y;
    Double_t original_x2 = dirPoint.x, original_y2 = dirPoint.y;
	for(int i=0; i<geCenters_x.size(); i++){
		// Shift the point coordinates w.r.t. the center of crystal
		point.x = point.x - geCenters_x[i];
        point.y = point.y - geCenters_y[i];
		dirPoint.x = dirPoint.x - geCenters_x[i];
		dirPoint.y = dirPoint.y - geCenters_y[i];
		// Once Ge has been centered, check hit
		/*cout << "SHIFTING: \n";
        cout << "[ProdPoint] " << point << endl;
        cout << "[DirPoint] " << dirPoint << endl;
        cout << "[GeRadius] " << geRadius << endl;
        */
        vector<Point> geHits = compute_intersection(point, dirPoint, geRadius);
		if(geHits.size() > 0){
		    Point geHit = geHits[0];    // closest hit on Ge crystal
			geHit.x += geCenters_x[i];
            geHit.y += geCenters_y[i];
            // Note: no shift on Z
			if(geHit.distance < distance){
				distance = geHit.distance;
                closestGeHit = geHit;
			}
		}
		// Restore original coord
        point.x = original_x1;
        point.y = original_y1;
        dirPoint.x = original_x2;
        dirPoint.y = original_y2;
	}
	/*
    cout << "[ProdPoint] " << point << endl;
    cout << "[DirPoint] " << dirPoint << endl;
    cout << "[Hit] " << closestGeHit << endl;
    cout << "Hit Defined? " << closestGeHit.isDefined() << endl << endl;
    */
	return closestGeHit;
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

Bool_t checkIfEnableGe(Point prodPoint, Point dirPoint){
	// Compute intersection with the Ge volume, meaning the cylinder that includes the Ge Strings
	Double_t ge_r = RIGHT_GE;
	vector<Point> geHits = compute_intersection(prodPoint, dirPoint, ge_r);
	if(geHits.size()==0)
		return false;	// No interesection at all, no difference between the two maps
	// If there is an interesection, look at where it occurs in Z
	// Use parametrization of line: line = P + t D, where P is a vector (point), D direction vector
	Point geHit = geHits[0];
	Double_t t = (geHit.x - prodPoint.x) / (dirPoint.x-prodPoint.x);	// From param: t = (x-x1)/dx
	Double_t z_int = prodPoint.z + t * (dirPoint.z-prodPoint.z);   		// z = z1 + t * dz, where z1 point, dz direction vector
	if(BOTTOMZ_GE < z_int && z_int < TOPZ_GE)
		return true;	// Interesection with Z in [BOTTOMGE,TOPGE], enable crystals
	return false;		// Intersection with Z outside the Ge volume, disable crystals
}

struct Point generateRandomPointInAngleSlice(Double_t radius, Double_t angleShift=shiftAngleForXYSampling,
                                             Double_t minZ=BOTTOMZ, Double_t maxZ=TOPZ){
    Double_t phi = -angleShift + rnd->Rndm() * (2*angleShift);  // Random Phi in angle
    struct Point point = Point(radius * cos(phi),
                               radius * sin(phi),
                               minZ + rnd->Rndm() * maxZ); // Random Z position
    return point;
}

pair<struct Point, Double_t> generatePhotonEmission(struct Point point){
    Double_t theta = rnd->Rndm() * PI;
    Double_t phi = rnd->Rndm() * 2 * PI;
    // generate direction point
    Point dirPoint = Point(point.x + 1 * sin(theta) * cos(phi),
                           point.y + 1 * sin(theta) * sin(phi),
                           point.z + 1 * cos(theta));
    Double_t lenTrajectory = 0;	// distance to hit the plane delimited by upper or lower boundary
    if(theta <= PI/2){
        lenTrajectory = (TOPZ-point.z)/cos(theta);	    // trajectory towards top Z
    }else{
        lenTrajectory = (BOTTOMZ-point.z)/cos(theta);    // trajectory towards bottom Z
    }
    return make_pair(dirPoint, lenTrajectory);
}

vector<struct Point> computeAllHits(struct Point prodPoint, struct Point dirPoint, Double_t lenTrajectory, Bool_t enableGe){
    // Possible scenarios: outer hit (close), inner hit (close), inner hit (far), outer hit(far).
    // assert: if hit shroud, always the outer first
    // assert: if 2nd hit, check if no Ge before
    vector<struct Point> hits;
    // 2 cases:
    struct Point geHit;
    geHit = getHitOnGermanium(prodPoint, dirPoint, enableGe);
    // First, check for germanium hit
    Double_t distanceToGeHit = INF;
    if(geHit.isDefined()){
        distanceToGeHit = geHit.distance;
    }
    for(auto shroud : {INNER_SHROUD, OUTER_SHROUD}) {
        Int_t shroudRadius;
        if(shroud==INNER_SHROUD)
            shroudRadius = INNER_SHROUD_RADIUS;
        else
            shroudRadius = OUTER_SHROUD_RADIUS;
        vector<Point> shroudHits = compute_intersection(prodPoint, dirPoint, shroudRadius, lenTrajectory);
        for(auto &hit : shroudHits){
            if(hit.isDefined() && prodPoint.checkSameDirection(dirPoint.x, dirPoint.y, hit) && hit.distance<distanceToGeHit) {    // Discard no hit on outer shroud
                hit.shroud = shroud;
                hits.push_back(hit);
            }
        }
    }
    std::sort(hits.begin(), hits.end());    // sort ascending distance
    return hits;
}

// Parameters
const Int_t n_inner_slices = 100;
const Int_t n_outer_slices = 100;
void runToyOpticsFromPoint(Double_t radius, Int_t nOptics, TH1D * &prInnerD, TH2D* &innerMap, TH2D* &outerMap){
	cout << "Running at R=" << radius << endl;
	// Fiber radius
	Double_t inner_r = 175, outer_r = 295;	// mm
	// Flag sample Ge
	Bool_t enableGe;
	// Loop on angles
	Bool_t debug = false;
	Bool_t printout = false;
	Double_t kInnerHits = 0, kOuterHits = 0;
	// Get coordinate
	struct Point prodPoint = Point();
	Int_t kGe = 0;
	while(nOptics>0){
		// Get random point (y1 according to angle shifting, z1 random) 
		prodPoint = generateRandomPointInAngleSlice(radius);
		pair<struct Point, Double_t> photonEmission = generatePhotonEmission(prodPoint);
		struct Point dirPoint = photonEmission.first;
		Double_t lenTrajectory = photonEmission.second;
		enableGe = checkIfEnableGe(prodPoint, dirPoint);
		// Debug
        assert(lenTrajectory>=0);
		if(printout){
            cout << "[ProdPoint]" << prodPoint << endl;
            cout << "[DirPoint]" << dirPoint << endl;
			cout << "[Other] Enable Ge: " << enableGe << ", ";
			cout << "Max Len Traj: " << lenTrajectory << "\n";
		}
        // Compute the possible hits
		Int_t pointPlacement = getPointPlacement(prodPoint.x, prodPoint.y, inner_r, outer_r);
		vector<Point> hits = computeAllHits(prodPoint, dirPoint, lenTrajectory, enableGe);
        if(radius>=295 && printout){
            cout << "[ProdPoint]" << prodPoint << endl;
            cout << "[DirPoint]" << dirPoint << endl;
            for(auto hit : hits){
                cout << "[Hit]" << hit;
            }
            cout << endl;
        }
        assert(pointPlacement!=INNER_REGION || hits.size()<=2);     //INNER_REGION => nr hits <=2
        assert(pointPlacement!=MIDDLE_REGION || hits.size()<=3);    //MIDDLE_REGION => nr hits <=3
        assert(pointPlacement!=OUTER_REGION || hits.size()<=4);     //OUTER_REGION => nr hits <=4

        if(hits.size() > 0) {
            vector <Double_t> probabilities, norm_probabilities;
            Int_t i = 0;
            Double_t NORM_FACTOR = 0;
            for (auto hit : hits) {
                // prob: (1-shroudCaptPr)^(i-1) * shroudCaptPr * Exp(-dist/L)
                i++;    // avoid i=0
                Double_t pr = pow(1 - shroudCapturePr, i - 1) * shroudCapturePr * TMath::Exp(-hit.distance / attenuationLen);
                probabilities.push_back(pr);
                NORM_FACTOR += pr;
            }
            for (auto pr : probabilities) {
                norm_probabilities.push_back(pr / NORM_FACTOR);
            }
            // Sampling the hit
            Double_t rndChoice = rnd->Rndm();
            Int_t id_currentHit = 0;
            while (rndChoice > norm_probabilities[id_currentHit]) {
                rndChoice -= norm_probabilities[id_currentHit];
                id_currentHit++;
            }
            assert(id_currentHit < hits.size());
            Point sampledHit = hits[id_currentHit];
            // Fill map with the sampled hit
            assert(sampledHit.isDefined());
            assert(sampledHit.shroud == INNER_SHROUD || sampledHit.shroud == OUTER_SHROUD);
            double angle = atan2(sampledHit.y, sampledHit.x);
            if (sampledHit.shroud == INNER_SHROUD) {
                innerMap->Fill(radius, angle);
                kInnerHits += 1;
            } else if (sampledHit.shroud == OUTER_SHROUD) {
                outerMap->Fill(radius, angle);
                kOuterHits += 1;
            }
            // DEBUG
            if (printout) {
                cout << "Hit on the fibers at: " << sampledHit.x << ", " << sampledHit.y << ", " << sampledHit.z
                     << " | angle: " << angle << "\n\n";
            }
            // Update counters and debug    
    	    nOptics--;
	    if (enableGe)
                kGe++;
        }
        if(debug)
            break;
	}
	if(kInnerHits + kOuterHits > 0)
		prInnerD->Fill(radius, kInnerHits/(kInnerHits+kOuterHits));
	// cout << "DEBUG - Enabled Ge: " << kGe << endl;
	// return sliced_detections;
}

void createSpatialMap_Sampling(){
	Double_t min_r = 0, max_r = 710, rbins = max_r - min_r + 1;
	Int_t angle_bins = 100;
	Int_t min_angle = -ceil(PI);
	Int_t max_angle = +ceil(PI);
	Int_t nOpticsPerPoint = 10000;
	// Initialize Ge Crystals
	initializePositionGeCrystals();
	// Create map
	TString outfile;
    	outfile.Form("ToySpatialMap_%dR_%dAngleSlices_%dops_AttLen%f_NewSampling_Last.root", (int)rbins, angle_bins, nOpticsPerPoint, attenuationLen);
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
