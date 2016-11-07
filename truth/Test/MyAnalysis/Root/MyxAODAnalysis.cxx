#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>
#include <EventLoop/OutputStream.h>
#include <EventLoopAlgs/NTupleSvc.h>
#include <EventLoopAlgs/AlgSelect.h>

#include "MyAnalysis/MyxAODAnalysis.h"
#include "MyAnalysis/OBJ_Base.h"

#include "xAODJet/JetContainer.h"

#include <xAODMuon/Muon.h>
#include <xAODMuon/MuonContainer.h>
#include <xAODMuon/MuonAuxContainer.h>

#include "xAODEventInfo/EventInfo.h"

#include "xAODTruth/TruthEvent.h"
#include "xAODTruth/TruthEventContainer.h"
#include "xAODTruth/TruthParticle.h"
#include "xAODTruth/TruthParticleContainer.h"
#include "xAODTruth/TruthVertex.h"
#include "xAODTruth/TruthVertexContainer.h"

#include "xAODRootAccess/TStore.h"
#include "xAODCore/ShallowCopy.h"

#include "Initialize.h"


using namespace std;


// this is needed to distribute the algorithm to the workers
ClassImp(MyxAODAnalysis)	


MyxAODAnalysis :: MyxAODAnalysis ()
{
  // Here you put any code for the base initialization of variables,
  // e.g. initialize all pointers to 0.  Note that you should only put
  // the most basic initialization here, since this method will be
  // called on both the submission and the worker node.  Most of your
  // initialization code will go into histInitialize() and
  // initialize().
  //m_setSysList=new std::vector<std::string>();
}

MyxAODAnalysis :: MyxAODAnalysis (string treename="physics")
{

  set = treename;
  //m_setSysList=new std::vector<std::string>();
  //

}



EL::StatusCode MyxAODAnalysis :: setupJob (EL::Job& job)
{
  // Here you put code that sets up the job on the submission object
  // so that it is ready to work with your algorithm, e.g. you can
  // request the D3PDReader service or add output files.  Any code you
  // put here could instead also go into the submission script.  The
  // sole advantage of putting it here is that it gets automatically
  // activated/deactivated when you add/remove the algorithm from your
  // job, which may or may not be of value to you.
  EL::OutputStream output1("tree_output");
  job.outputAdd (output1);

  job.useXAOD ();   // filetype is xAOD, tell EventLoop that we actually want to use the xAODRootAccess   --Zhang

  // let's initialize the algorithm to use the xAODRootAccess package
  xAOD::Init( "MyxAODAnalysis" ).ignore(); // call before opening first file   
  // xAOD::Init(const char * appname = "xAOD::Init") return a TReturnCode type var, which call the ignore method to ignore the return code, mark it as checked, Why so??? --Zhang

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode MyxAODAnalysis :: histInitialize ()
{
    // Here you do everything that needs to be done at the very
    // beginning on each worker node, e.g. create histograms and 
    // output trees. This method gets called before any input files
    // are connected.   --Zhang
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode MyxAODAnalysis :: fileExecute ()
{
    // Here you do everything you needs to be done exactly once for 
    // every single file, e.g. collect a list of all lumi-blocks
    // processed     -- Zhnag
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode MyxAODAnalysis :: changeInput (bool firstFile)
{
  // Here you do everything you need to do when we change input files,
  // e.g. resetting branch addresses on trees.  If you are using
  // D3PDReader or a similar service this method is not needed.
  return EL::StatusCode::SUCCESS;
}

// < move the "initialize ()" to a separate file "initialize.h"

EL::StatusCode MyxAODAnalysis :: execute ()
{
    m_eventCounter++;
    if(m_eventCounter==1 || m_eventCounter%1000==0)  
        cout << "event counter " << m_eventCounter << endl;
     
//    if(m_eventCounter == 100) return EL::StatusCode::FAILURE;

    // get eventinfo  
    const xAOD::EventInfo* eventInfo = 0;
    if( !m_event->retrieve( eventInfo, "EventInfo").isSuccess() )
    {
      Error("execute ()", "Failed to retrieve EventInfo. Exiting." );
      return EL::StatusCode::FAILURE;
    }

    uint32_t run = 0;
    unsigned long long event = 0;
//  /  float mu = eventInfo->actualInteractionsPerCrossing(); 
    vector<float> vw;
    float mcWeight = 1.0;
    if(isMC) {
      run = eventInfo->mcChannelNumber(); 
      event = eventInfo->mcEventNumber();
      vw = eventInfo->mcEventWeights(); 
      mcWeight = vw.size()>0?vw[0]:1.0;
      if(mcWeight>5) mcWeight = 1;
    }
    //
    // END OF EVENTINFO
    // 

    //
    // muon, electron, jet, met
    //
    VOmuon goodm;	
    VOelectron goode;
    VOjet goodj;
    VOneutrino goodn;

    // get TruthEvent
    const xAOD::TruthEventContainer* TruthEvtContainer = 0;
    if( !m_event->retrieve( TruthEvtContainer, "TruthEvents").isSuccess() ){
      Error("excute()", "Failed to retrieve TruthParticleContainer container. Exiting." );
      return EL::StatusCode::FAILURE;
    }

    //
    // Truth Particles Dealing 
    //

    // retrieve truth leptons
    const xAOD::TruthParticleContainer* mc_particle = 0;
    if ( ! m_event->retrieve(mc_particle,"TruthParticles").isSuccess() ) {
        Error("execute()", "Failed to retrieve TruthParticleContainer container. Exiting.");
        return EL::StatusCode::FAILURE;
    }

    cout << "event: " << m_eventCounter << endl;
    for (auto truth = mc_particle->begin(); truth != mc_particle->end(); ++truth ) {
      int pdg = (*truth)->pdgId();
      int status = (*truth)->status();
      if( pdg == 23 && status == 62 ) 
      { nz++; }
      if( abs(pdg) == 11 )
      {
	  if(status == 1)
	  { ne_1++; }
	  else if(status == 23)
	  { ne_23++; }
	  else 
	  { 
	      ne_other++;
	      cout << "electron's other status: " << pdg << "\t" << status << endl;
  	  }
      }
      else if( abs(pdg) == 13 )
      {
	  if(status == 1)
	  { nm_1++; }
	  else if(status == 23)
	  { nm_23++; }
	  else 
	  { 
	      nm_other++;	
	      cout << "muon's other status: " << pdg << "\t" << status << endl;
	  }
      }
      else if ( abs(pdg) == 12 )
      {
	  if(status == 1)
	  { nve++;  }
	  else 
	  {
	      nve_other++;
	      cout << "electron neutrino's other status: " << pdg << "\t" << status << endl;
	  }
      }
      else if ( abs(pdg) == 14 )
      {
	  if(status == 1)
	  { nvm++;  }
	  else 
	  {
	      nvm_other++;
	      cout << "muon neutrino's other status: " << pdg << "\t" << status << endl;
	  }
      }
      else if ( abs(pdg) == 16 )
      {
	  if(status == 1)
	  { nvt++;  }
	  else 
	  {
	      nvt_other++;
	      cout << "tau neutrino's other status: " << pdg << "\t" << status << endl;
	  }
      }
//      int barcode = (*truth)->barcode();
    }
    cout << "total leptons: " << ne_1+ne23+nm_1+nm_23 << "\t total neutrinos: " << nve+nvm+nvt << endl
	 << "electron status 1: " << ne_1 << "\t status 23: " << ne_23 << endl 
	 << "muon status 1: " << nm_1 << "\t status 23: " << nm_23 << endl 
	 << "electron status other: " << ne_other << "\t muon status other: " << nm_other << endl 
	 << "electron neutrino status 1: " << nve << "\t status other: " << nve_other << endl
	 << "muon neutrino status 1: " << nvm << "\t status other: " << nvm_other << endl
	 << "tau neutrino status 1: " << nvt << "\t status other: " << nvt_other << endl
	 << endl;

  // retrieve truth jets
  const xAOD::JetContainer* TruthJets = 0;
  if(! m_event->retrieve(TruthJets, "AntiKt4TruthJets").isSuccess())
  {
    Error("excute()", "Failed to retrieve Truth Jets info. Exiting.");
    return EL::StatusCode::FAILURE;
  }

  goodj.clear();
  xAOD::JetContainer::const_iterator trJets_itr = TruthJets->begin();
  xAOD::JetContainer::const_iterator trJets_end = TruthJets->end();
  for(; trJets_itr != trJets_end; ++trJets_itr)
  {
    if(((*trJets_itr)->pt() > 20.e3 && fabs((*trJets_itr)->eta()) < 4.5))
    {
      OBJ_JET jetInfo;

      double px = (*trJets_itr)->px();
      double py = (*trJets_itr)->py();
      double pz = (*trJets_itr)->pz();
      double e  = (*trJets_itr)->e();
      jetInfo.L.SetPxPyPzE(px, py, pz, e);
      goodj.push_back(jetInfo);
    }
  }
  cout << "jet number: " << goodj.size() << endl;

}


EL::StatusCode MyxAODAnalysis :: postExecute ()
{
  // Here you do everything that needs to be done after the main event
  // processing.  This is typically very rare, particularly in user
  // code.  It is mainly used in implementing the NTupleSvc.
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode MyxAODAnalysis :: finalize ()
{

  // This method is the mirror image of initialize(), meaning it gets
  // called after the last event has been processed on the worker node
  // and allows you to finish up any objects you created in
  // initialize() before they are written to disk.  This is actually
  // fairly rare, since this happens separately for each worker node.
  // Most of the time you want to do your post-processing on the
  // submission node after all your histogram outputs have been
  // merged.  This is different from histFinalize() in that it only
  // gets called on worker nodes that processed input events.

  cout << "####" << endl;
  printf("Finalize : %i events have been processed !\n", m_eventCounter);
  printf("Fiducial : %i events", m_fiducial);
  printf("Finalize MyxAODAnalysis !\n");
  cout << "####" << endl;
  cout << endl;

  time(&end);  // time(&start) in initial() file
  double timedif = difftime(end,start);
  if(timedif>3600) { printf("Finalize : Time Cost: %f hours\n", timedif/3600.); }
  else if(timedif>60) { printf("Finalize : Time Cost: %f minutes\n", timedif/60.); }
  else { printf("Finalize : Time Cost: %f second\n", timedif); }

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode MyxAODAnalysis :: histFinalize ()
{


  // This method is the mirror image of histInitialize(), meaning it
  // gets called after the last event has been processed on the worker
  // node and allows you to finish up any objects you created in
  // histInitialize() before they are written to disk.  This is
  // actually fairly rare, since this happens separately for each
  // worker node.  Most of the time you want to do your
  // post-processing on the submission node after all your histogram
  // outputs have been merged.  This is different from finalize() in
  // that it gets called on all worker nodes regardless of whether
  // they processed input events.
  return EL::StatusCode::SUCCESS;
}
