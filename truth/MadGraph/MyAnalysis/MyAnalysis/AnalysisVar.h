//<AnalysisVar.h>
//<Aim to provide variables and objects for analysis>
//<Yusheng WU, April 2011, Ann Arbor>

#ifndef AnalysisVar_h
#define AnalysisVar_h

#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <math.h> 
#include <string>
#include <map>
#include <vector>
#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>

//<Use STD map datatypes to manage counting and histograms>
using namespace std;

// <Map type STRING:STRING...?
typedef map<string, map<string, string> > MapType2_String;
//<Map type STRING:INT...>
typedef map<string, int> MapType_Int;
typedef map<string, map<string,int> > MapType2_Int;
typedef map<string, map<string, vector<int> > > MapType2_VInt;

//<Map type STRING: UNSIGNED LONG LONG INT...>
typedef map<string, map<string, unsigned long long> > MapType2_ULong64;

//<Map type STRING:FLOAT...>
typedef map<string, float> MapType_Float;
typedef map<string, map<string, float> > MapType2_Float;
typedef map<string, map<string, vector<float> > > MapType2_VFloat;

//<Map type STRING:DOUBLE...>
typedef map<string, double> MapType_Double;
typedef map<string, map<string, double> > MapType2_Double;
typedef map<string, map<string, vector<double> > > MapType2_VDouble;


//<Map type STRING:TTree...>


//<AnalysisVar Class>
class AnalysisVar {

    public:        
        MapType2_Int TreeIntVar;
        MapType2_ULong64 TreeUloVar;
        MapType2_String TreeStrVar;
        MapType2_Float TreeFltVar;
        MapType2_VFloat TreeFltVVar;
        MapType2_VInt TreeIntVVar;
        MapType2_VDouble TreeDblVVar;

	TTree * Tree;
};

#endif
