#ifndef SIGNAL_H
#define SIGNAL_H

/* signal1.1.h for functions1.5.h and analysis3.2.C or above */

class signal {

private:
	double energy;
	int strip;
	int status;
public:
	signal() {
		energy = -1;
		strip = -1;
		status = 0; //0 means not matched, 1 is matching candidate
	}
	void SetStatus(int c=0) {
		status = c;
	}
	void SetValues(double a, int b, int c=0) {
		energy = a;
		strip = b;
		status = c;
	}
	void Apply(signal a) {	
		energy = a.GetEnergy();
		strip = a.GetStrip();
		status = a.GetStatus();
	}
	double GetEnergy( ) { return energy; }
	int GetStrip() {return strip; }
	int GetStatus() {return status; }
};

#endif
