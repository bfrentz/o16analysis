#ifndef HIT_H
#define HIT_H

/* hit1.1.h for functions1.5.h and analysis3.2.C or above */

class hit {

private:
	double energy;
	double espread;
	int stripF;
	int stripB;
	int detectorN;
	bool twohits;
// errorcode: 0 -- normal, 
// 1-- split hit on back strip, 2-- double hit on back strip,
// 3-- split hit on front strip, 4-- double hit on front strip, later more 5 - saturation ...etc
	int errcode;

	double theta;
	double phi;
	double px;
	double py;
	double pz;

public:
	hit();
	//Get Methodes
	double GetEnergy() { return energy; }
	double GetESpread() { return espread; }
	int GetStripF() { return stripF; }
	int GetStripB() { return stripB; }
	int GetDetectorN() { return detectorN; }
	double GetAnglePhi() { return phi; }
	double GetAngleTheta() { return theta; }
	double GetMomentumX() { return px; }
	double GetMomentumY() { return py; }
	double GetMomentumZ() { return pz; }
	bool GetDouble() { return twohits; }
	int GetErrCode() { return errcode; }
	//Set Methodes
	void SetEnergy(double a) { energy = a; }
	void SetESpread(double a) { espread = a; }
	void SetStripF(int a) { stripF = a; }
	void SetStripB(int a) { stripB = a; }
	void SetDetectorN(int a) { detectorN = a; }
	void SetValues(double, double, int, int, int, bool, int = 0);
	void SetAngles(double a,double b) { theta = a; phi = b; }
	void SetMomentum(double a, double b, double c) { px = a; py = b; pz = c; }
	void Apply(hit);
};

hit::hit() {
	energy = -1;
	espread = -1;
	detectorN = -1;
	stripF = -1;
	stripB = -1;
	theta = 0;
	phi = 0;
	px = 0;
	py = 0;
	pz = 0;
	twohits = false;
}
//SetValues([Energy],[ESpread],[DetectorNumber],[FrontStrip],[BackStrip],[Is it double hits])
void hit::SetValues(double a, double b, int c, int d, int e, bool f, int g) {
	energy = a;
	espread = b;
	detectorN = c;
	stripF = d;
	stripB = e;
	twohits = f;
	errcode = g;
}

void hit::Apply(hit a) {
	energy = a.GetEnergy();
	espread = a.GetESpread();
	detectorN = a.GetDetectorN();
	stripF = a.GetStripF();
	stripB = a.GetStripB();
}

#endif
