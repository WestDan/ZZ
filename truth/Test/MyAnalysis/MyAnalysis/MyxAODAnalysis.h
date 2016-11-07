#ifndef MyAnalysis_MyxAODAnalysis_H
#define MyAnalysis_MyxAODAnalysis_H

#include <EventLoop/Algorithm.h>

// Infrastructure include(s):
#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"
#include "AnalysisVar.h"


// #include "PATInterfaces/SystematicRegistry.h"
// #include "PATInterfaces/SystematicCode.h"
// #include "PATInterfaces/CorrectionCode.h"
// #include "PATInterfaces/SystematicsUtil.h"
// #include "PATInterfaces/SystematicVariation.h" 


#include <TTree.h>
#include <TH1.h>
#include <TLorentzVector.h>
#include <time.h>
#include <TMacro.h>
using namespace std;

class MyxAODAnalysis : public EL::Algorithm, public AnalysisVar  // no definition even in AnalysisVar.h ??? --Zhang
{
  // put your configuration variables here as public variables.
  // that way they can be set directly from CINT and python.
public:

  time_t start, end; //!

  xAOD::TEvent *m_event;  //!
  //xAOD::EventInfo* eventInfo; //!
  // count number of events
  bool isMC;
  int m_eventCounter; //!
  int m_fiducial; //!

  double m_sumOfWeights;  //! sumofWeights

  float sumOfWeights; //!
  float sumOfWeightsSquared; //!

  string set; //!
  string m_FileName; //!
  

  // variables that don't get filled at submission time should be
  // protected from being send from the submission node to the worker
  // node (done by the //!)
public:

  // this is a standard constructor
  MyxAODAnalysis ();
  MyxAODAnalysis (string );

  // these are the functions inherited from Algorithm
  virtual EL::StatusCode setupJob (EL::Job& job);
  virtual EL::StatusCode fileExecute ();
  virtual EL::StatusCode histInitialize ();
  virtual EL::StatusCode changeInput (bool firstFile);
  virtual EL::StatusCode initialize ();
  virtual EL::StatusCode execute ();
  virtual EL::StatusCode postExecute ();
  virtual EL::StatusCode finalize ();
  virtual EL::StatusCode histFinalize ();

  // this is needed to distribute the algorithm to the workers
  ClassDef(MyxAODAnalysis, 1);

  void InitStrVec(vector<string>& out, string in, string de=",");

  void InitTreeVar(string varlist, string type);
  void AddVarIntoTree(TTree *tree);
};

#endif
