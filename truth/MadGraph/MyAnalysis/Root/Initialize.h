#include <iostream>
#include <fstream>
#include "TSystem.h"
using namespace std;

#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>

#include "xAODMuon/Muon.h"  // only muon needed??? --Zhang

#include "MyAnalysis/MyxAODAnalysis.h"
#include "MyAnalysis/AnalysisVar.h"

#include "ElectronPhotonSelectorTools/AsgElectronLikelihoodTool.h"

#include "AssociationUtils/OverlapRemovalTool.h"
#include "PileupReweighting/PileupReweightingTool.h"

#include "xAODCutFlow/CutBookkeeper.h"
#include "xAODCutFlow/CutBookkeeperContainer.h"

//ort::inputAccessor_t selectAcc("selected");
//ort::inputDecorator_t selectDec("selected");
//ort::outputAccessor_t overlapAcc("overlaps");
//ort::objLinkAccessor_t objLinkAcc("overlapObject");

EL::StatusCode MyxAODAnalysis :: initialize ()
{
    // Here you do everything that you need to do after the first input
    // file has been connected and before the first event is precessed,
    // e.g. create additional histograms based on which variables are
    // available in the input files. You can also create all of your
    // histograms and trees in here, but be aware that this method
    // doesn't get called if no events are processed. So any method
    // you create here won't be available in the output if you have no 
    // input events.   --Zhang

  
  time(&start);

  m_event = wk()->xaodEvent();
  Info("initialize()", "Number of events = %lli", m_event->getEntries() );
  m_FileName = wk()->inputFile()->GetName();
  cout << endl << endl << endl
      << "File name: " << m_FileName << endl;

  const xAOD::EventInfo* eventInfo = 0;
  if( !m_event->retrieve( eventInfo, "EventInfo").isSuccess() )
  {
    Error("execute ()", "Failed to retrieve EventInfo. Exiting." );
    return EL::StatusCode::FAILURE;
  }
  isMC = eventInfo->eventType( xAOD::EventInfo::IS_SIMULATION );
  // eventType :
  // IS_SIMULATION | true: simulation, false: data
  // IS_TESTBEAM   | true: testbeam,   false: full detector
  // IS_CALIBRATION| true: calibration,false: physics

  uint64_t nEventsProcessed=0;
  sumOfWeights = 0.;
  sumOfWeightsSquared = 0.;

  TFile *file1 = wk()->getOutputFile ("tree_output");
  Tree = new TTree("Madgraph", "Madgraph_truth");
  Tree->SetDirectory(file1);
  AddVarIntoTree(Tree);

  m_eventCounter = 0;
  m_fiducial = 0;

  m_sumOfWeights = 0.0;

  return EL::StatusCode::SUCCESS;
}


void MyxAODAnalysis :: InitStrVec(vector<string>& out, string in, string de) {
    int pos=0, pos_pre=0;
    while(true) {
        pos=in.find(de,pos_pre);
        if(pos==-1) {out.push_back(in.substr(pos_pre,in.size()-pos_pre)); break;}
        else  out.push_back(in.substr(pos_pre,pos-pos_pre)); // this depend on that no space between literal word and delimiter.
        pos_pre=pos+1;
    }
}

void MyxAODAnalysis :: InitTreeVar(string varlist, string type) {
    vector<string> variables;
    InitStrVec(variables, varlist, ",");

    if(type=="I") {
      for(int i=0; i<(int)variables.size(); i++) {
        TreeIntVar[variables[i]]["Value"]=-9999;
      }
    }
    if(type=="U") {
      for(int i=0; i<(int)variables.size(); i++) {
        TreeUloVar[variables[i]]["Value"]=0;
      }
    }
    if(type=="F") {
      for(int i=0; i<(int)variables.size(); i++) {
        TreeFltVar[variables[i]]["Value"]=-9999.0;
      }
    }
    if(type=="S") {
      for(int i=0; i<(int)variables.size(); i++) {
        TreeStrVar[variables[i]]["Value"]=" ";
      }
    }
    if(type=="vector<F>") {
      for(int i=0; i<(int)variables.size(); i++) {
        TreeFltVVar[variables[i]]["Value"].clear();
      }
    }

    if(type=="vector<I>") {
      for(int i=0; i<(int)variables.size(); i++) {
        TreeIntVVar[variables[i]]["Value"].clear();
      }
    }
}

void MyxAODAnalysis :: AddVarIntoTree(TTree *tree) {
    
    char buf[1000];
    getcwd(buf, sizeof(buf));
    setenv("LOCAL", buf, 1);
    FILE* fp;
    char result_buf[1000];
    fp = popen("find $LOCAL -name MiniTree.txt", "r");
    fgets(result_buf, sizeof(result_buf), fp);

    string varfile(result_buf);
    size_t len = varfile.size();
    varfile.erase(len-1);
    ifstream file;
    file.open(varfile.c_str(), ios::out);

    if (file.is_open())  {
      char line[256];
      while (!file.eof() )  {
        string varname, type;

        file.getline (line,100);
        string sline(line);

        if(sline.find("Int_t")!=string::npos) {
          type = "I";
          varname = sline.substr(9);
          type = varname+"/"+type;
          InitTreeVar(varname,"I");
          tree->Branch(varname.c_str(),&TreeIntVar[varname]["Value"], type.c_str());
        }
        if(sline.find("Ulo_t")!=string::npos) {
          type = "U";
          varname = sline.substr(9);
          type = varname+"/"+type;
          InitTreeVar(varname, "U");
          tree->Branch(varname.c_str(),&TreeUloVar[varname]["Value"]);
        }
        if(sline.find("Float_t")!=string::npos) {
          type = "F";
          varname = sline.substr(9);
          type = varname+"/"+type;
          InitTreeVar(varname, "F");
          tree->Branch(varname.c_str(),&TreeFltVar[varname]["Value"], type.c_str());
        }
        if(sline.find("Str")!=string::npos) {
          type = "S";
          varname = sline.substr(9);
          type = varname+"/"+type;
          InitTreeVar(varname, "S");
          tree->Branch(varname.c_str(),&TreeStrVar[varname]["Value"]);
        }
        if(sline.find("Vector_I")!=string::npos) {
          varname = sline.substr(9);
          InitTreeVar(varname, "Vector_I");
          tree->Branch(varname.c_str(), &TreeIntVVar[varname]["Value"]);
        }
        if(sline.find("Vector_F")!=string::npos) {
          varname = sline.substr(9);
          InitTreeVar(varname, "Vector_F");
          tree->Branch(varname.c_str(), &TreeFltVVar[varname]["Value"]);
        }
      }
    }
}
