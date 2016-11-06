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
  if(m_eventCounter==1 || m_eventCounter%1000==0)  // for monitering the dealing precess
    cout << "event counter " << m_eventCounter << endl;
   
//  if(m_eventCounter == 100) return EL::StatusCode::FAILURE;

  // get eventinfo  
  const xAOD::EventInfo* eventInfo = 0;
  if( !m_event->retrieve( eventInfo, "EventInfo").isSuccess() )
  {
    Error("execute ()", "Failed to retrieve EventInfo. Exiting." );
    return EL::StatusCode::FAILURE;
  }

//  float aveIntPC = eventInfo->actualInteractionsPerCrossing() + eventInfo->averageInteractionsPerCrossing();
  uint32_t run = 0;
  unsigned long long event = 0;
///  float mu = eventInfo->actualInteractionsPerCrossing(); 
  vector<float> vw;
  float mcWeight = 1.0;
  if(isMC) {
    run = eventInfo->mcChannelNumber(); 
    event = eventInfo->mcEventNumber();
    vw = eventInfo->mcEventWeights(); 
    // SETNAME[0] = "physics"
    mcWeight = vw.size()>0?vw[0]:1.0;
    if(mcWeight>5) mcWeight = 1;
  }
  //
  // END OF EVENTINFO
  // 

  //
  // muon, electron, jet, met
  //
  VOmuon goodm;	// !!! truth muon 
  VOelectron goode;
  VOjet goodj;
  VOneutrino goodn;

  // get TruthEvent
  const xAOD::TruthEventContainer* TruthEvtContainer = 0;
  if( !m_event->retrieve( TruthEvtContainer, "TruthEvents").isSuccess() ){
    Error("excute()", "Failed to retrieve Truth info. Exiting." );
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

    goodm.clear();
    goode.clear();
    goodn.clear();
    OBJ_Z truthz1, truthz2;
    truthz1.L.SetPtEtaPhiE(9999.e3, -10, 0, 99999.e3);
    truthz2.L.SetPtEtaPhiE(9999.e3, -10, 0, 99999.e3);
    for (auto truth = mc_particle->begin(); truth != mc_particle->end(); ++truth ) {
      int pdg = (*truth)->pdgId();
      int status = (*truth)->status();
//      int barcode = (*truth)->barcode();

      /*
      if( abs(pdg) == 23 && status == 62 )  // !!! Z in Madgraph
      {
	  double pt  = (*truth)->pt();
	  double eta = (*truth)->eta();
	  double phi = (*truth)->phi();
	  double e   = (*truth)->e();
	  OBJ_Z zInfo;
	  zInfo.L.SetPtEtaPhiE(pt, eta, phi, e);
	  continue;
      }
       */

      // choose prompt(barcode) stable (status == 1) muon (pdg == 13) !!! status 
      if( ( abs(pdg) != 12 && abs(pdg) != 14 && abs(pdg) != 16 && abs(pdg)!=13 && abs(pdg)!=11 ) || status!=1 ) continue;

      // must have a parent and eventually pointing back to the Z boson
      int nparent = (*truth)->nParents();
      if(nparent!=1) continue;

      bool findZ = false;
      const xAOD::TruthParticle* theparent = (*truth)->parent(0);

      while(nparent == 1) {
        int pdgparent = theparent->pdgId();
// cout << pdgparent << "\t" << theparent->status() << endl;
	if(pdgparent == 23) 
	{
	    findZ = true;
	}
	else {
	  nparent = theparent->nParents();
	  theparent = theparent->parent(0);
	}
	if(findZ) break;
      }

      if (!findZ) continue;

      // find Z
      int charge = (*truth)->charge();
      double pt = (*truth)->pt();
      double eta = (*truth)->eta();
      double phi = (*truth)->phi();
      double e   = (*truth)->e();
      if( pt < 7.e3 || abs(eta) > 2.5) continue;    // !!! Is it valid for neutrino

      TLorentzVector temp;
      temp.SetPtEtaPhiE(pt, eta, phi, e);   // !!!! set E, not set M

      // dressing
      for(auto truth2 = mc_particle->begin(); truth2 != mc_particle->end(); ++truth2){
	int pdg2 = (*truth2)->pdgId();
	int status2 = (*truth2)->status();
//	int barcode2 = (*truth2)->barcode();
	if(pdg2!=22 || status2!=1 ) continue;

	double pt = (*truth2)->pt() ;
	double eta = (*truth2)->eta();
	double phi = (*truth2)->phi();
	TLorentzVector temp_photon;
	temp_photon.SetPtEtaPhiM(pt, eta, phi, 0);
	double dR = temp_photon.DeltaR(temp);
	if(dR<0.1) temp += temp_photon;
      }

      double z_pt  = theparent->pt();
      double z_eta = theparent->eta();
      double z_phi = theparent->phi();
      double z_e   = theparent->e();

      if(abs(pdg) == 13)
      { // muon
	  OBJ_MUON muonInfo;
	  muonInfo.charge = charge;
	  muonInfo.L = temp;
	  goodm.push_back(muonInfo);

	  truthz1.L.SetPtEtaPhiE(z_pt, z_eta, z_phi, z_e);
      }
      else if (abs(pdg) == 11)
      {	// electron
	  OBJ_ELECTRON elecInfo;
	  elecInfo.charge = charge;
	  elecInfo.L = temp;
	  goode.push_back(elecInfo);

	  truthz1.L.SetPtEtaPhiE(z_pt, z_eta, z_phi, z_e);
      }
      else 
      {	// neutrino
	  OBJ_NEUTRINO neuInfo;
	  neuInfo.L = temp;
	  goodn.push_back(neuInfo);

	  truthz2.L.SetPtEtaPhiE(z_pt, z_eta, z_phi, z_e);
      }
    }

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

    int nmuon = goodm.size();
    int nele  = goode.size();
    int nneu  = goodn.size();
    if( ! ( (nmuon == 2 && nele == 0 && nneu == 2) || ( nmuon == 0  && nele == 2 && nneu == 2 ) ) )    // llvv channel 
    {	
	return EL::StatusCode::SUCCESS;	    
    }
if(truthz1.L.Pt() > 1.e6)
{   
    cout << "event: " << m_eventCounter << "truthz_ll "
         << "\t\t" << nmuon << "\t" << nele << "\t\t"
	 << endl;
}
else if( truthz2.L.Pt() > 1.e6)
{   
    cout << "event: " << m_eventCounter << "truthz_vv "
         << "\t\t" << nneu << "\t\t"
	 << endl;
}

    // extract leptons
    
    float l1_pt = -9999.0, l1_eta = -9999.0, l1_phi = -9999.0, l1_e = -9999.0; 
    float l2_pt = -9999.0, l2_eta = -9999.0, l2_phi = -9999.0, l2_e = -9999.0; 
    float l3_pt = -9999.0, l3_eta = -9999.0, l3_phi = -9999.0, l3_e = -9999.0; 
    float l4_pt = -9999.0, l4_eta = -9999.0, l4_phi = -9999.0, l4_e = -9999.0; 
    float delta_R_ll_1 = -9999.0, delta_R_ll_2 = -9999.0;
    float delta_R_zz = -9999.0, delta_R_jj = -9999.0;
    float delta_R_zz_jj = -9999.0;
    float leadingj_rap = -9999., subleadingj_rap = -9999., zz_rap = -9999.;
    float delta_eta_jj = -9999.0;
    float centrality = -9999.0;

    // sort leptons
    OBJ_Z z[2];
    if(nmuon == 2)
    {
	if(goodm[0].charge == goodm[1].charge)
	{
	    cout << "Event: " << m_eventCounter << " is same charged muons" << endl;
	    return EL::StatusCode::SUCCESS;
	}
	if(goodm[0].L.Pt() < goodm[1].L.Pt() )
	{
	    OBJ_MUON temp = goodm[0];
	    goodm[0] = goodm[1];
	    goodm[1] = temp;
	}
	z[0].L = goodm[0].L + goodm[1].L;

	l1_pt  = goodm[0].L.Pt();
	l1_eta = goodm[0].L.Eta();
	l1_phi = goodm[0].L.Phi();
	l1_e   = goodm[0].L.E();
	l2_pt  = goodm[1].L.Pt();
	l2_eta = goodm[1].L.Eta();
	l2_phi = goodm[1].L.Phi();
	l2_e   = goodm[1].L.E();
    }
    else if (nele == 2)
    {
	if(goode[0].charge == goode[1].charge)
	{
	    cout << "Event: " << m_eventCounter << " is same charged electrons" << endl;
	    return EL::StatusCode::SUCCESS;
	}
	if( goode[0].L.Pt() < goode[1].L.Pt() )
	{
	    OBJ_ELECTRON temp = goode[0];
	    goode[0] = goode[1];
	    goode[1] = temp;
	}
	z[0].L = goode[0].L + goode[1].L;

	l1_pt  = goode[0].L.Pt();
	l1_eta = goode[0].L.Eta();
	l1_phi = goode[0].L.Phi();
	l1_e   = goode[0].L.E();
	l2_pt  = goode[1].L.Pt();
	l2_eta = goode[1].L.Eta();
	l2_phi = goode[1].L.Phi();
	l2_e   = goode[1].L.E();
    }

    delta_R_ll_1 = sqrt( pow( (l1_eta - l2_eta), 2) + pow( (l1_phi - l2_phi), 2) );

    {
	if (goodn[0].L.Pt() < goodn[1].L.Pt() )
	{
	    OBJ_NEUTRINO temp = goodn[0];
	    goodn[0] = goodn[1];
	    goodn[1] = temp;
	}
	z[1].L = goodn[0].L + goodn[1].L;

	l3_pt  = goodn[0].L.Pt();
	l3_eta = goodn[0].L.Eta();
	l3_phi = goodn[0].L.Phi();
	l3_e   = goodn[0].L.E();
	l4_pt  = goodn[1].L.Pt();
	l4_eta = goodn[1].L.Eta();
	l4_phi = goodn[1].L.Phi();
	l4_e   = goodn[1].L.E();
    }
    delta_R_ll_2 = sqrt( pow( (l3_eta - l4_eta), 2) + pow( (l3_phi - l4_phi), 2) );

    TLorentzVector zz = z[0].L + z[1].L;
    delta_R_zz = sqrt( pow( (z[0].L.Eta() - z[1].L.Eta()), 2) + pow( (z[0].L.Phi() - z[1].L.Phi()), 2));
    zz_rap = 0.5*log( (zz.E() + zz.Pz()) / (zz.E() - zz.Pz()));

    TLorentzVector truthzz = truthz1.L + truthz2.L;

    // truth jets
    int njet = goodj.size();
    OBJ_JET leadingj, subleadingj;
    if ( njet < 2)
    {	return EL::StatusCode::SUCCESS ;    }

    for ( int i=0; i<njet; i++)
    {
	if( goodj[i].L.Pt() > leadingj.L.Pt() )
	{
	    subleadingj = leadingj;
	    leadingj = goodj[i];
	}
	else if ( goodj[i].L.Pt() > subleadingj.L.Pt() )
	{
	    subleadingj = goodj[i];
	}
    }

    TLorentzVector jj = leadingj.L + subleadingj.L;
    delta_eta_jj = fabs( leadingj.L.Eta() - subleadingj.L.Eta() );
    delta_R_jj = sqrt( pow(delta_eta_jj, 2) + pow( (leadingj.L.Phi() - subleadingj.L.Phi()), 2));
    leadingj_rap = 0.5*log( (leadingj.L.E() + leadingj.L.Pz()) / ( leadingj.L.E() - leadingj.L.Pz()));
    subleadingj_rap = 0.5*log( (subleadingj.L.E() + subleadingj.L.Pz()) / ( subleadingj.L.E() - subleadingj.L.Pz()));
    centrality = (zz_rap - (leadingj_rap + subleadingj_rap) / 2) / fabs( leadingj_rap - subleadingj_rap );

    TreeStrVar["filename"]["Value"]     = m_FileName;
    TreeIntVar["run"]["Value"]          = run;
    TreeUloVar["event"]["Value"]        = event;
    TreeIntVar["njet"]["Value"]         = njet;
    TreeIntVar["nmuon"]["Value"]        = nmuon;
    TreeIntVar["nele"]["Value"]         = nele;
    TreeFltVar["mcweight"]["Value"]     = mcWeight;

    TreeFltVar["l1_pt"]["Value"]  = l1_pt;
    TreeFltVar["l1_eta"]["Value"] = l1_eta;
    TreeFltVar["l1_phi"]["Value"] = l1_phi;
    TreeFltVar["l1_e"]["Value"]   = l1_e;
    TreeFltVar["l2_pt"]["Value"]  = l2_pt;
    TreeFltVar["l2_eta"]["Value"] = l2_eta;
    TreeFltVar["l2_phi"]["Value"] = l2_phi;
    TreeFltVar["l2_e"]["Value"]   = l2_e;
    TreeFltVar["l3_pt"]["Value"]  = l3_pt;
    TreeFltVar["l3_eta"]["Value"] = l3_eta;
    TreeFltVar["l3_phi"]["Value"] = l3_phi;
    TreeFltVar["l3_e"]["Value"]   = l3_e;
    TreeFltVar["l4_pt"]["Value"]  = l4_pt;
    TreeFltVar["l4_eta"]["Value"] = l4_eta;
    TreeFltVar["l4_phi"]["Value"] = l4_phi;
    TreeFltVar["l4_e"]["Value"]   = l4_e;
    TreeFltVar["leadingj_pt"]["Value"]  = leadingj.L.Pt();
    TreeFltVar["leadingj_eta"]["Value"] = leadingj.L.Eta();
    TreeFltVar["leadingj_phi"]["Value"] = leadingj.L.Phi();
    TreeFltVar["leadingj_e"]["Value"]   = leadingj.L.E();
    TreeFltVar["subleadingj_pt"]["Value"]  = subleadingj.L.Pt();
    TreeFltVar["subleadingj_eta"]["Value"] = subleadingj.L.Eta();
    TreeFltVar["subleadingj_phi"]["Value"] = subleadingj.L.Phi();
    TreeFltVar["subleadingj_e"]["Value"]   = subleadingj.L.E();
    TreeFltVar["z1_pt"]["Value"]        = z[0].L.Pt();
    TreeFltVar["z1_eta"]["Value"]       = z[0].L.Eta();
    TreeFltVar["z1_phi"]["Value"]       = z[0].L.Phi();
    TreeFltVar["z1_e"]["Value"]         = z[0].L.E();
    TreeFltVar["z1_m"]["Value"]         = z[0].L.M();
    TreeFltVar["z2_pt"]["Value"]        = z[1].L.Pt();
    TreeFltVar["z2_eta"]["Value"]       = z[1].L.Eta();
    TreeFltVar["z2_phi"]["Value"]       = z[1].L.Phi();
    TreeFltVar["z2_e"]["Value"]         = z[1].L.E();
    TreeFltVar["z2_m"]["Value"]         = z[1].L.M();
    TreeFltVar["truthz1_pt"]["Value"]   = truthz1.L.Pt();
    TreeFltVar["truthz1_eta"]["Value"]  = truthz1.L.Eta();
    TreeFltVar["truthz1_phi"]["Value"]  = truthz1.L.Phi();
    TreeFltVar["truthz1_e"]["Value"]    = truthz1.L.E();
    TreeFltVar["truthz1_m"]["Value"]    = truthz1.L.M();
    TreeFltVar["truthz2_pt"]["Value"]   = truthz2.L.Pt();
    TreeFltVar["truthz2_eta"]["Value"]  = truthz2.L.Eta();
    TreeFltVar["truthz2_phi"]["Value"]  = truthz2.L.Phi();
    TreeFltVar["truthz2_e"]["Value"]    = truthz2.L.E();
    TreeFltVar["truthz2_m"]["Value"]    = truthz2.L.M();
    TreeFltVar["zz_pt"]["Value"]        = zz.Pt();
    TreeFltVar["zz_eta"]["Value"]       = zz.Eta();
    TreeFltVar["zz_phi"]["Value"]       = zz.Phi();
    TreeFltVar["zz_e"]["Value"]         = zz.E();
    TreeFltVar["zz_m"]["Value"]         = zz.M();
    TreeFltVar["truthzz_pt"]["Value"]   = truthzz.Pt();
    TreeFltVar["truthzz_eta"]["Value"]  = truthzz.Eta();
    TreeFltVar["truthzz_phi"]["Value"]  = truthzz.Phi();
    TreeFltVar["truthzz_e"]["Value"]    = truthzz.E();
    TreeFltVar["truthzz_m"]["Value"]    = truthzz.M();
    TreeFltVar["jj_pt"]["Value"]        = jj.Pt();
    TreeFltVar["jj_eta"]["Value"]       = jj.Eta();
    TreeFltVar["jj_phi"]["Value"]       = jj.Phi();
    TreeFltVar["jj_e"]["Value"]         = jj.E();
    TreeFltVar["jj_m"]["Value"]         = jj.M();
    TreeFltVar["delta_eta_jj"]["Value"] = delta_eta_jj;
    TreeFltVar["delta_R_ll_1"]["Value"] = delta_R_ll_1;
    TreeFltVar["delta_R_ll_2"]["Value"] = delta_R_ll_2;
    TreeFltVar["delta_R_zz"]["Value"]   = delta_R_zz;
    TreeFltVar["delta_R_jj"]["Value"]   = delta_R_jj;
    TreeFltVar["delta_R_zz_jj"]["Value"]= delta_R_zz_jj;
    TreeFltVar["centrality"]["Value"]   = centrality;
    Tree->Fill();

    m_fiducial++;
    return EL::StatusCode::SUCCESS;
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
