#ifndef MyAnalysis_OBJ_H
#define MyAnalysis_OBJ_H

//<Common Header>
#include <TLorentzVector.h>
#include <string>
#include <vector>
using namespace std;

//<Predefined Constants>
const double PI=3.1415926;
const double ZMass=91.1876e3; //MeV
//const double m_mass=105.658367; //MeV
//const double e_mass = 0.510998910; //MeV
const double m_mass=105.658; //MeV
const double e_mass = 0.511; //MeV
const double Unit_GeV = 1000.; //MeV

// Truth Object
// Truth_Electron, Truth_Muon

//<Object Class>
//<General, Muon, Electron, Jet, MET>
class OBJ_MUON {
    public:
        int charge;
        TLorentzVector L;
};

class OBJ_ELECTRON {
    public:
        double charge;
        TLorentzVector L;
};       

class OBJ_NEUTRINO {
    public:
	TLorentzVector L;
};

class OBJ_JET {
    public:
        TLorentzVector L;
};

class OBJ_Z {
    public:
	TLorentzVector L;
};

typedef vector<OBJ_MUON> VOmuon;
typedef vector<OBJ_ELECTRON> VOelectron;
typedef vector<OBJ_JET> VOjet;
typedef vector<OBJ_Z>  VOz;
typedef vector<OBJ_NEUTRINO> VOneutrino ;

#endif
