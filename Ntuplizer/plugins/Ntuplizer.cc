// -*- C++ -*-
//
// Package:    Analyzer/Ntuplizer
// Class:      Ntuplizer
// 
/**\class Ntuplizer Ntuplizer.cc Analyzer/Ntuplizer/plugins/Ntuplizer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Christophe ochando
//         Created:  Mon, 10 Mar 2014 14:51:20 GMT
//
//


// MY include
#include "Ntuplizer.h"

// C++ include files
#include <memory>

// CMSSW include
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

//#include <DataFormats/Common/interface/MergeableCounter.h>
//#include <DataFormats/Common/interface/View.h>
//#include <DataFormats/Candidate/interface/Candidate.h>
//#include <DataFormats/PatCandidates/interface/CompositeCandidate.h>
//#include <DataFormats/MuonReco/interface/Muon.h>
//#include <DataFormats/MuonReco/interface/MuonFwd.h>
//
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
//
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
//CMSSW/ RecoEcal/ EgammaCoreTools/ interface/ EcalClusterLazyTools.h
#include "FWCore/Utilities/interface/isFinite.h"
//#include "DataFormats/PatCandidates/interface/Muon.h"
//#include "DataFormats/PatCandidates/interface/MET.h"
//#include "DataFormats/PatCandidates/interface/Electron.h"
//#include "DataFormats/PatCandidates/interface/Jet.h"
//#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
//"
// MET
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/METFwd.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"
//
//#include "DataFormats/Math/interface/LorentzVector.h"
//#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "RecoVertex/PrimaryVertexProducer/interface/PrimaryVertexSorter.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"

// Trigger
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

//
// class declaration
//

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//using namespace std;
using namespace reco;
using namespace edm;

// =============================================================================================
// constructor
Ntuplizer::Ntuplizer(const edm::ParameterSet& iConfig) :
// ==============================================================================================
EleTag_ (iConfig.getParameter<edm::InputTag> ("EleTag")),
VerticesTag_(iConfig.getParameter<edm::InputTag> ("VerticesTag")),
HLTTag_(iConfig.getParameter<edm::InputTag> ("HLTTag")),
isMC_ (iConfig.getParameter<bool>("isMC")),
ispythia6_ (iConfig.getParameter<bool>("ispythia6")),
PileupSrc_ ("addPileupInfo"),
MVAidCollection_ (iConfig.getParameter<std::vector<edm::InputTag> >("MVAId")),
runGsfRefitter (iConfig.getParameter<bool>("runGsfRefitter")),
GSFTrajColl (iConfig.getParameter<std::string>("GSFTrajectoryInput"))

// std::vector<edm::InputTag> MVAidCollection_;
{


}

// =============================================================================================
// destructor
Ntuplizer::~Ntuplizer()
// =============================================================================================
{
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  delete m_electrons ;
//   if (fill_L1trigger) {
//     delete m_L1emIso;
//     delete m_L1emNonIso;
//   }
//   delete m_muons;
//   delete _m_jets_pf;
//   delete m_photons;
  
  if(isMC_ ) {
    delete _m_MC_gen_V;
    //delete _m_MC_gen_photons;
    delete _m_MC_gen_leptons;
    delete _m_MC_gen_Higgs;
    delete _m_MC_gen_leptons_status1;
    delete _m_MC_gen_leptons_status2;
  } // if MC
  
  
}

// =============================================================================================
// ------------ method called once each job just before starting event loop  ------------
void Ntuplizer::beginJob()
//=============================================================================================
{
  // Book histograms
  edm::Service<TFileService> fs ;
  _mytree  = fs->make <TTree>("simpleRootTree","simpleRootTree"); 
  
  //// Counters
  //_mytree->Branch("Nevt_Gen",&Nevt_Gen,"Nevt_Gen/I");
  //_mytree->Branch("Nevt_Skim",&Nevt_afterSkim,"Nevt_Skim/I");
  
  // Global
  _mytree->Branch("nEvent",&_nEvent,"nEvent/I");
  _mytree->Branch("nRun",&_nRun,"nRun/I");
  _mytree->Branch("nLumi",&_nLumi,"nLumi/I");
  
  // Pile UP
  _mytree->Branch("PU_N",&_PU_N,"PU_N/I");
  
  // Trigger
  _mytree->Branch("trig_fired_names",&trig_fired_names,"trig_fired_names[10000]/C");

  // Vertices
  _mytree->Branch("vtx_N",&_vtx_N,"vtx_N/I");

  // Electrons
  _mytree->Branch("ele_N",&ele_N,"ele_N/I");
  m_electrons = new TClonesArray ("TLorentzVector");
  _mytree->Branch ("electrons", "TClonesArray", &m_electrons, 256000,0);


 
  _mytree->Branch("ele_signedEstimateSumPred", &ele_signedEstimateSumPred);//;&ele_signedEstimateSumPred,"ele_signedEstimateSumPred[50]/D");
  _mytree->Branch("ele_signedEstimateSumPred_A", &ele_signedEstimateSumPred_A,"ele_signedEstimateSumPred_A[50]/F");
  _mytree->Branch("ele_propagatorSignedEstimateSumPred", &ele_propagatorSignedEstimateSumPred);



 
  _mytree->Branch("ele_echarge",&ele_echarge,"ele_echarge[50]/I");
  //
  _mytree->Branch("ele_he",&ele_he,"ele_he[50]/D");
  _mytree->Branch("ele_hebc",&ele_hebc,"ele_hebc[50]/D");
  _mytree->Branch("ele_oldhe",&ele_oldhe,"ele_oldhe[50]/D");
  _mytree->Branch("ele_oldhebc",&ele_oldhebc,"ele_oldhebc[50]/D");
  //
  _mytree->Branch("ele_eseedpout",&ele_eseedpout,"ele_eseedpout[50]/D");
  _mytree->Branch("ele_ep",&ele_ep,"ele_ep[50]/D");
  _mytree->Branch("ele_eseedp",&ele_eseedp,"ele_eseedp[50]/D");
  _mytree->Branch("ele_eelepout",&ele_eelepout,"ele_eelepout[50]/D");
  //
  _mytree->Branch("ele_pin_mode",&ele_pin_mode,"ele_pin_mode[50]/D");
  _mytree->Branch("ele_pout_mode",&ele_pout_mode,"ele_pout_mode[50]/D");
  _mytree->Branch("ele_pTin_mode",&ele_pTin_mode,"ele_pTin_mode[50]/D");
  _mytree->Branch("ele_pTout_mode",&ele_pTout_mode,"ele_pTout_mode[50]/D");
  //
  _mytree->Branch("ele_deltaetaseed",&ele_deltaetaseed,"ele_deltaetaseed[50]/D");
  _mytree->Branch("ele_deltaphiseed",&ele_deltaphiseed,"ele_deltaphiseed[50]/D");
  _mytree->Branch("ele_deltaetaele",&ele_deltaetaele,"ele_deltaetaele[50]/D");
  _mytree->Branch("ele_deltaphiele",&ele_deltaphiele,"ele_deltaphiele[50]/D");
  _mytree->Branch("ele_deltaetain",&ele_deltaetain,"ele_deltaetain[50]/D");
  _mytree->Branch("ele_deltaphiin",&ele_deltaphiin,"ele_deltaphiin[50]/D");
  //
  _mytree->Branch("ele_sigmaietaieta",&ele_sigmaietaieta,"ele_sigmaietaieta[50]/D");
  _mytree->Branch("ele_sigmaetaeta",&ele_sigmaetaeta,"ele_sigmaetaeta[50]/D");
  _mytree->Branch("ele_sigmaiphiiphi",&ele_sigmaiphiiphi,"ele_sigmaiphiiphi[50]/D");
  _mytree->Branch("ele_e15",&ele_e15,"ele_e15[50]/D");
  _mytree->Branch("ele_e25max",&ele_e25max,"ele_e25max[50]/D");
  _mytree->Branch("ele_e55",&ele_e55,"ele_e55[50]/D");
  //_mytree->Branch("ele_e1",&ele_e1,"ele_e1[50]/D");
  _mytree->Branch("ele_r9",&ele_r9,"ele_r9[50]/D");
  //
  _mytree->Branch("ele_oldsigmaietaieta",&ele_oldsigmaietaieta,"ele_oldsigmaietaieta[50]/D");
  _mytree->Branch("ele_oldsigmaetaeta",&ele_oldsigmaetaeta,"ele_oldsigmaetaeta[50]/D");
  _mytree->Branch("ele_oldsigmaiphiiphi",&ele_oldsigmaiphiiphi,"ele_oldsigmaiphiiphi[50]/D");
  _mytree->Branch("ele_oldsigmaietaiphi",&ele_oldsigmaietaiphi,"ele_oldsigmaietaiphi[50]/D");
  _mytree->Branch("ele_olde15",&ele_olde15,"ele_olde15[50]/D");
  _mytree->Branch("ele_olde25max",&ele_olde25max,"ele_olde25max[50]/D");
  _mytree->Branch("ele_olde55",&ele_olde55,"ele_olde55[50]/D");
  //_mytree->Branch("ele_olde1",&ele_olde1,"ele_olde1[50]/D");
  _mytree->Branch("ele_oldr9",&ele_oldr9,"ele_oldr9[50]/D");
  //
  //_mytree->Branch("ele_e33",&ele_e33,"ele_e33[50]/D");  ---> No ECAL Reduced Collection
  //_mytree->Branch("ele_e2overe9",&ele_e2overe9,"ele_e2overe9[50]/D");  ---> No ECAL Reduced Collection
  // 
  _mytree->Branch("ele_fbrem",&ele_fbrem,"ele_fbrem[50]/D");
  _mytree->Branch("ele_trackfbrem",&ele_trackfbrem,"ele_trackfbrem[50]/D");
  _mytree->Branch("ele_SCfbrem",&ele_SCfbrem,"ele_SCfbrem[50]/D");
  _mytree->Branch("ele_pfSCfbrem",&ele_pfSCfbrem,"ele_pfSCfbrem[50]/D");
  _mytree->Branch("ele_nbrem",&ele_nbrem,"ele_nbrem[50]/I");
  _mytree->Branch("ele_eClass",&ele_eClass,"ele_eClass[50]/I");
  //
  _mytree->Branch("ele_mva",&ele_mva,"ele_mva[50]/D");
  //
  _mytree->Branch("ele_isbarrel",&ele_isbarrel,"ele_isbarrel[50]/I");
  _mytree->Branch("ele_isendcap",&ele_isendcap,"ele_isendcap[50]/I");
  _mytree->Branch("ele_isEBetaGap",&ele_isEBetaGap,"ele_isEBetaGap[50]/I");
  _mytree->Branch("ele_isEBphiGap",&ele_isEBphiGap,"ele_isEBphiGap[50]/I");
  _mytree->Branch("ele_isEEdeeGap",&ele_isEEdeeGap,"ele_isEEdeeGap[50]/I");
  _mytree->Branch("ele_isEEringGap",&ele_isEEringGap,"ele_isEEringGap[50]/I");
  _mytree->Branch("ele_isecalDriven",&ele_isecalDriven,"ele_isecalDriven[50]/I");
  _mytree->Branch("ele_istrackerDriven",&ele_istrackerDriven,"ele_istrackerDriven[50]/I");
  //
  _mytree->Branch("ele_valid_hits",&ele_valid_hits,"ele_valid_hits[50]/I");
  _mytree->Branch("ele_lost_hits",&ele_lost_hits,"ele_lost_hits[50]/I");
  _mytree->Branch("ele_gsfchi2",&ele_gsfchi2,"ele_gsfchi2[50]/D");	
  //
  _mytree->Branch("ele_dxy",&ele_dxy,"ele_dxy[50]/D");
  _mytree->Branch("ele_dxyB",&ele_dxyB,"ele_dxyB[50]/D");
  _mytree->Branch("ele_dz",&ele_dz,"ele_dz[50]/D");
  _mytree->Branch("ele_dzB",&ele_dzB,"ele_dzB[50]/D");
  _mytree->Branch("ele_dsz",&ele_dsz,"ele_dsz[50]/D");
  _mytree->Branch("ele_dszB",&ele_dszB,"ele_dszB[50]/D");
  //
  _mytree->Branch("ele_dzPV", &ele_dzPV, "ele_dzPV[50]/D");
  _mytree->Branch("ele_d0", &ele_d0, "ele_d0[50]/D");
  _mytree->Branch("ele_d0err",&ele_d0err, "ele_d0err[50]/D");
  //
  _mytree->Branch("ele_IP",&ele_IP,"ele_IP[50]/D");
  _mytree->Branch("ele_IPError",&ele_IPError,"ele_IPError[50]/D");	
  _mytree->Branch("ele_SIP",&ele_SIP,"ele_SIP[50]/D");
  //  _mytree->Branch("ele_tkSumPt_dr03",&ele_tkSumPt_dr03,"ele_tkSumPt_dr03[50]/D"); 
  //   _mytree->Branch("ele_ecalRecHitSumEt_dr03",&ele_ecalRecHitSumEt_dr03,"ele_ecalRecHitSumEt_dr03[50]/D"); 
  //   _mytree->Branch("ele_hcalDepth1TowerSumEt_dr03",&ele_hcalDepth1TowerSumEt_dr03,"ele_hcalDepth1TowerSumEt_dr03[50]/D"); 
  //   _mytree->Branch("ele_hcalDepth2TowerSumEt_dr03",&ele_hcalDepth2TowerSumEt_dr03,"ele_hcalDepth2TowerSumEt_dr03[50]/D"); 
  //   _mytree->Branch("ele_tkSumPt_dr04",&ele_tkSumPt_dr04,"ele_tkSumPt_dr04[50]/D"); 
  //   _mytree->Branch("ele_ecalRecHitSumEt_dr04",&ele_ecalRecHitSumEt_dr04,"ele_ecalRecHitSumEt_dr04[50]/D"); 
  //   _mytree->Branch("ele_hcalDepth1TowerSumEt_dr04",&ele_hcalDepth1TowerSumEt_dr04,"ele_hcalDepth1TowerSumEt_dr04[50]/D"); 
  //   _mytree->Branch("ele_hcalDepth2TowerSumEt_dr04",&ele_hcalDepth2TowerSumEt_dr04,"ele_hcalDepth2TowerSumEt_dr04[50]/D"); 
  //
  _mytree->Branch("ele_conv_dist",&ele_conv_dist,"ele_conv_dist[50]/D");
  _mytree->Branch("ele_conv_dcot",&ele_conv_dcot,"ele_conv_dcot[50]/D");
  _mytree->Branch("ele_conv_radius",&ele_conv_radius,"ele_conv_radius[50]/D");
  _mytree->Branch("ele_expected_inner_hits",&ele_expected_inner_hits,"ele_expected_inner_hits[50]/I");
  _mytree->Branch("ele_vtxconv",&ele_vtxconv,"ele_vtxconv[50]/I");
  //
  _mytree->Branch("ele_pfChargedHadIso",&ele_pfChargedHadIso,"ele_pfChargedHadIso[50]/D");
  _mytree->Branch("ele_pfNeutralHadIso",&ele_pfNeutralHadIso,"ele_pfNeutralHadIso[50]/D");
  _mytree->Branch("ele_pfPhotonIso",&ele_pfPhotonIso,"ele_pfPhotonIso[50]/D");

  _mytree->Branch("ele_pfChargedIso", &ele_pfChargedIso, "ele_pfChargedIso[50]/D");
  _mytree->Branch("ele_pfSumPUIso", &ele_pfSumPUIso, "ele_pfSumPUIso[50]/D");
  //  _mytree->Branch("ele_pfChargedHadPUIso",&ele_pfChargedHadPUIso,"ele_pfChargedHadPUIso[50]/D");
  //   _mytree->Branch("ele_pfCombRelIso",&ele_pfCombRelIso,"ele_pfCombRelIso[50]/D");
  // 
 
  //
  _mytree->Branch("ele_sclRawE", &ele_sclRawE, "ele_sclRawE[50]/D");
  _mytree->Branch("ele_sclE",    &ele_sclE,    "ele_sclE[50]/D");
  _mytree->Branch("ele_sclEt",   &ele_sclEt,   "ele_sclEt[50]/D");
  _mytree->Branch("ele_sclEta",  &ele_sclEta,  "ele_sclEta[50]/D");
  _mytree->Branch("ele_sclPhi",  &ele_sclPhi,  "ele_sclPhi[50]/D");
  _mytree->Branch("ele_sclNclus",  &ele_sclNclus,  "ele_sclNclus[50]/I");
  _mytree->Branch("ele_sclphiwidth", &ele_sclphiwidth, "ele_sclphiwidth[50]/D");
  _mytree->Branch("ele_scletawidth", &ele_scletawidth, "ele_scletawidth[50]/D");
  //
  _mytree->Branch("ele_sclsubE", &ele_sclsubE, "ele_sclsubE[50][20]/D");
  _mytree->Branch("ele_sclsubEta", &ele_sclsubEta, "ele_sclsubEta[50][20]/D");
  _mytree->Branch("ele_sclsubPhi", &ele_sclsubPhi, "ele_sclsubPhi[50][20]/D");
  _mytree->Branch("ele_sclsubisseed", &ele_sclsubisseed, "ele_sclsubisseed[50][20]/I");
  //
  _mytree->Branch("ele_psE", &ele_psE, "ele_psE[50]/D");
  //
  _mytree->Branch("ele_ecalE",   &ele_ecalE,   "ele_ecalE[50]/D");
  _mytree->Branch("ele_ecalErr",   &ele_ecalErr,   "ele_ecalErr[50]/D");
  _mytree->Branch("ele_trackErr",  &ele_trackErr,  "ele_trackErr[50]/D");
  _mytree->Branch("ele_combErr",   &ele_combErr,   "ele_combErr[50]/D");
  _mytree->Branch("ele_PFcombErr", &ele_PFcombErr, "ele_PFcombErr[50]/D");
  
  // Electron ID
  _mytree->Branch("ele_mvaphys14",   &ele_mvaphys14,   "ele_mvaphys14[50]/D");
  _mytree->Branch("ele_mvaphys14fix",   &ele_mvaphys14fix,   "ele_mvaphys14fix[50]/D");
    

  //  _mytree->Branch("ele_ecalRegressionEnergy",   ele_ecalRegressionEnergy,   "ele_ecalRegressionEnergy[50]/D");
  //   _mytree->Branch("ele_ecalRegressionError", ele_ecalRegressionError, "ele_ecalRegressionError[50]/D");
  //   _mytree->Branch("ele_ecalTrackRegressionEnergy",&ele_ecalTrackRegressionEnergy,"ele_ecalTrackRegressionEnergy[50]/D");
  //   _mytree->Branch("ele_ecalTrackRegressionError",&ele_ecalTrackRegressionError,"ele_ecalTrackRegressionError[50]/D");
  //   _mytree->Branch("ele_ecalScale",ele_ecalScale,"ele_ecalScale[50]/D");
  //   _mytree->Branch("ele_ecalSmear",ele_ecalSmear,"ele_ecalSmear[50]/D");
  //   _mytree->Branch("ele_ecalRegressionScale",ele_ecalRegressionScale,"ele_ecalRegressionScale[50]/D");
  //   _mytree->Branch("ele_ecalRegressionSmear",ele_ecalRegressionSmear,"ele_ecalRegressionSmear[50]/D");
  //   _mytree->Branch("ele_ecalTrackRegressionScale",ele_ecalTrackRegressionScale,"ele_ecalTrackRegressionScale[50]/D");
  //   _mytree->Branch("ele_ecalTrackRegressionSmear",ele_ecalTrackRegressionSmear,"ele_ecalTrackRegressionSmear[50]/D");
 	
  // Variables for the mva. Most of them are duplicated, but since they are corrected at analysis level, it could be dangerous
  //  _mytree->Branch("ele_mvafbrem", ele_mvafbrem,"ele_mvafbrem[50]/D");
  //   _mytree->Branch("ele_mvadetain", ele_mvadetain,"ele_mvadetain[50]/D");
  //   _mytree->Branch("ele_mvadphiin", ele_mvadphiin,"ele_mvadphiin[50]/D");
  //   _mytree->Branch("ele_mvasieie", ele_mvasieie,"ele_mvasiesie[50]/D");
  //   _mytree->Branch("ele_mvahoe", ele_mvahoe,"ele_mvahoe[50]/D");
  //   _mytree->Branch("ele_mvaeop", ele_mvaeop,"ele_mvaeop[50]/D");
  //   _mytree->Branch("ele_mvae1x5e5x5", ele_mvae1x5e5x5,"ele_mvae1x5e5x5[50]/D");
  //   _mytree->Branch("ele_mvaeleopout", ele_mvaeleopout,"ele_mvaeleopout[50]/D");
  _mytree->Branch("ele_kfchi2",&ele_kfchi2,"ele_kfchi2[50]/D");
  _mytree->Branch("ele_kfhits",&ele_kfhits,"ele_kfhits[50]/I");
  _mytree->Branch("ele_gsfhits",&ele_gsfhits,"ele_gsfhits[50]/I");
   
  //   _mytree->Branch("ele_mvamishits", ele_mvamishits,"ele_mvamisthits[50]/I");
  //   _mytree->Branch("ele_mvadist", ele_mvadist,"ele_mvadist[50]/D");
  //   _mytree->Branch("ele_mvadcot", ele_mvadcot,"ele_mvadcot[50]/D");
  //   _mytree->Branch("ele_mvaeta", ele_mvaeta,"ele_mvaeta[50]/D");
  //   _mytree->Branch("ele_mvapt", ele_mvapt,"ele_mvapt[50]/D");
  //   _mytree->Branch("ele_mvaecalseed", ele_mvaecalseed,"ele_mvaecalseed[50]/I");
 
  
  // PFMET
  _mytree->Branch("met_pf_et",&_met_pf_et,"met_pf_et/D");
  _mytree->Branch("met_pf_px",&_met_pf_px,"met_pf_px/D");
  _mytree->Branch("met_pf_py",&_met_pf_py,"met_pf_py/D");
  _mytree->Branch("met_pf_phi",&_met_pf_phi,"met_pf_phi/D");
  _mytree->Branch("met_pf_set",&_met_pf_set,"met_pf_set/D");
  _mytree->Branch("met_pf_sig",&_met_pf_sig,"met_pf_sig/D");

  // Truth Leptons
  //cout << "truth leptons" << endl;
  _m_MC_gen_V = new TClonesArray ("TLorentzVector");
  _mytree->Branch ("MC_gen_V", "TClonesArray", &_m_MC_gen_V, 256000,0);
  _mytree->Branch ("MC_gen_V_pdgid",&_MC_gen_V_pdgid, "MC_gen_V_pdgid[10]/D");
  //
  _m_MC_gen_Higgs = new TClonesArray ("TLorentzVector");
  _mytree->Branch ("MC_gen_Higgs", "TClonesArray", &_m_MC_gen_Higgs, 256000,0);
  //_mytree->Branch ("MC_gen_Higgs_pdgid",&_MC_gen_Higgs_pdgid, "MC_gen_Higgs_pdgid[10]/D");
  //
  _m_MC_gen_leptons         = new TClonesArray ("TLorentzVector");
  _m_MC_gen_leptons_status1 = new TClonesArray ("TLorentzVector");
  _m_MC_gen_leptons_status2 = new TClonesArray ("TLorentzVector");
  _mytree->Branch ("MC_gen_leptons", "TClonesArray", &_m_MC_gen_leptons, 256000,0);
  _mytree->Branch ("MC_gen_leptons_status1", "TClonesArray", &_m_MC_gen_leptons_status1, 256000,0);
  _mytree->Branch ("MC_gen_leptons_status2", "TClonesArray", &_m_MC_gen_leptons_status2, 256000,0);
  _mytree->Branch ("MC_gen_leptons_pdgid",&_MC_gen_leptons_pdgid, "MC_gen_leptons_pdgid[30]/D");
  _mytree->Branch ("MC_gen_leptons_status1_pdgid",&_MC_gen_leptons_status1_pdgid, "MC_gen_leptons_status1_pdgid[30]/D");
  _mytree->Branch ("MC_gen_leptons_status1_FromWZ",&_MC_gen_leptons_status1_FromWZ, "MC_gen_leptons_status1_FromWZ[30]/I");
  _mytree->Branch ("MC_gen_leptons_status1_FromTaus",&_MC_gen_leptons_status1_FromTaus, "MC_gen_leptons_status1_FromTaus[30]/I");
  _mytree->Branch ("MC_gen_leptons_status1_FromNonPrompt",&_MC_gen_leptons_status1_FromNonPrompt, "MC_gen_leptons_status1_FromNonPrompt[30]/I");
  _mytree->Branch ("MC_gen_leptons_status1_FromBC",&_MC_gen_leptons_status1_FromBC, "MC_gen_leptons_status1_FromBC[30]/I");
  _mytree->Branch ("MC_gen_leptons_status2_pdgid",&_MC_gen_leptons_status2_pdgid, "MC_gen_leptons_status2_pdgid[30]/D");
  //_mytree->Branch ("MC_pthat",&_MC_pthat,"MC_pthat/D");
  _mytree->Branch ("MC_flavor",&_MC_flavor,"MC_flavor[2]/I");
  //_m_MC_gen_photons         = new TClonesArray ("TLorentzVector");
  //_mytree->Branch ("MC_gen_photons", "TClonesArray", &_m_MC_gen_photons, 256000,0);
  //_mytree->Branch ("MC_gen_photons_isFSR",&_MC_gen_photons_isFSR,"MC_gen_photons_isFSR[5000]/I");
  
 
}


//
// member functions
//
// =============================================================================================
// ------------ method called for each event  ------------
void Ntuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
// =============================================================================================
{
   using namespace edm;

   Init();

   //cout << "salut" << endl;

   FillEvent(iEvent, iSetup);
   
   FillVertices(iEvent, iSetup);

   m_electrons -> Clear();
   FillElectrons(iEvent, iSetup);
   
   FillMET (iEvent, iSetup);


   if(isMC_ ) {
     //	cout << "truth2" << endl;
     _m_MC_gen_V->Clear();
     _m_MC_gen_Higgs->Clear();
     //_m_MC_gen_photons->Clear();
     _m_MC_gen_leptons->Clear();
     _m_MC_gen_leptons_status1->Clear();
     _m_MC_gen_leptons_status2->Clear();
     //cout << "truth" << endl;
     FillTruth(iEvent, iSetup);
   }

   _mytree->Fill();

// #ifdef THIS_IS_AN_EVENT_EXAMPLE
//    Handle<ExampleData> pIn;
//    iEvent.getByLabel("example",pIn);
// #endif
   
// #ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
//    ESHandle<SetupData> pSetup;
//    iSetup.get<SetupRecord>().get(pSetup);
// #endif
}

// =============================================================================================
void Ntuplizer::FillEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup)
//=============================================================================================
{

  _nEvent = iEvent.id().event();
  _nRun   = iEvent.id().run();
  _nLumi  = iEvent.luminosityBlock();

  // ----------------
  // Fired Triggers
  // ----------------
  //cout << "fired" << endl;
  Handle<edm::TriggerResults> triggerResultsHandle;
  iEvent.getByLabel (HLTTag_,triggerResultsHandle);
  const edm::TriggerNames & triggerNames = iEvent.triggerNames(*triggerResultsHandle);
  
  //Get List of available Triggers
  //for (int in=0;in<(int)triggerNames.size();in++) {
  //cout << " Trigger Names " << in << " = " << triggerNames.triggerName(in) << endl;
  //} // for loop in triggernames
  
  // LOOP Over Trigger Results
  char trig_fired_names_local[10000];
  strcpy(trig_fired_names_local,"*");
  for (int iHLT = 0 ; 
       iHLT<static_cast<int>(triggerResultsHandle->size()); 
       ++iHLT) {	
    
    if (triggerResultsHandle->accept (iHLT)) {
      if ( strlen(trig_fired_names_local) <= 9950) {
	{
	  const char* c_str();
	  string hlt_string = triggerNames.triggerName(iHLT);
	  strcat(trig_fired_names_local,hlt_string.c_str());
	  strcat(trig_fired_names_local,"*");
	}
      }
    } // if HLT
  }
  strcpy(trig_fired_names,trig_fired_names_local);
  
  //cout << "trig = " << trig_fired_names << endl;
  //cout << "" << endl;

} // end of FillEvent


// =============================================================================================
void Ntuplizer::FillVertices(const edm::Event& iEvent, const edm::EventSetup& iSetup)
//=============================================================================================
{
  
  Handle<vector<reco::Vertex> >  recoPrimaryVertexCollection;
  //   //iEvent.getByLabel("goodPrimaryVertices",recoPrimaryVertexCollection);
  iEvent.getByLabel(VerticesTag_, recoPrimaryVertexCollection);

  //const reco::VertexCollection & vertices = *recoPrimaryVertexCollection.product();
  
//   // 	edm::Handle<reco::BeamSpot> recoBeamSpotHandle;
//   // 	iEvent.getByType(recoBeamSpotHandle);
//   // 	const reco::BeamSpot bs = *recoBeamSpotHandle;
  
  //int vtx_counter=0;
  _vtx_N = recoPrimaryVertexCollection->size();
  
} // end of FillVertices

// =============================================================================================
void Ntuplizer::FillElectrons(const edm::Event& iEvent, const edm::EventSetup& iSetup)
//=============================================================================================
{
  
  edm::Handle<reco::GsfElectronCollection> electronsCol;
  iEvent.getByLabel(EleTag_, electronsCol);
  
  InputTag  vertexLabel(string("offlinePrimaryVertices"));
  Handle<reco::VertexCollection> thePrimaryVertexColl;
  iEvent.getByLabel(VerticesTag_ ,thePrimaryVertexColl);

  Vertex dummy;
  const Vertex *pv = &dummy;
  if (thePrimaryVertexColl->size() != 0) {
    pv = &*thePrimaryVertexColl->begin();
  } else { // create a dummy PV                                                                                                                                                                                             
    Vertex::Error e;
    e(0, 0) = 0.0015 * 0.0015;
    e(1, 1) = 0.0015 * 0.0015;
    e(2, 2) = 15. * 15.;
    Vertex::Point p(0, 0, 0);
    dummy = Vertex(p, e, 0, 0, 0);
  }

  edm::ESHandle<TransientTrackBuilder> builder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", builder);
  TransientTrackBuilder thebuilder = *(builder.product());

  //load the conversion collection
  edm::Handle<reco::ConversionCollection> conversions_h;
  iEvent.getByLabel("allConversions", conversions_h);
  
  // get the beam spot
  edm::Handle<reco::BeamSpot> beamspot_h;
  iEvent.getByLabel("offlineBeamSpot", beamspot_h);
  const reco::BeamSpot &beamSpot = *(beamspot_h.product());

  // get the value map for eiD
  edm::Handle<edm::ValueMap<float> >  mapMVAcollection_phys14;
  edm::Handle<edm::ValueMap<float> >  mapMVAcollection_phys14fix;
  if (MVAidCollection_.size()>0 ) iEvent.getByLabel(MVAidCollection_[0] , mapMVAcollection_phys14);
  if (MVAidCollection_.size()>1 ) iEvent.getByLabel(MVAidCollection_[1], mapMVAcollection_phys14fix);
//  const edm::ValueMap<float> & mapMVA_phys14    = *mapMVAcollection_phys14;
//  const edm::ValueMap<float> & mapMVA_phys14fix = *mapMVAcollection_phys14fix;
  
  
  TClonesArray & electrons = *m_electrons;
  int counter = 0;
  ele_N = electronsCol->size();

  //cout << "ele N = " << ele_N << endl;
  ele_signedEstimateSumPred.clear();
  ele_propagatorSignedEstimateSumPred.clear();
  for (reco::GsfElectronCollection::const_iterator ielectrons=electronsCol->begin(); ielectrons!=electronsCol->end();++ielectrons) {
    if(counter>49) continue;



    float tmp_ele_signedEstimateSumPred = -200.;
    float tmp_ele_propagatorSignedEstimateSumPred = -200.;
    const Trajectory* theGSFTraj = nullptr;
    if(runGsfRefitter) {
        edm::LogInfo("Test") << "GSF refit was run";
        edm::Handle<std::vector<Trajectory> > GSFTrajCollectionHandle;
        iEvent.getByLabel(GSFTrajColl, GSFTrajCollectionHandle);
        theGSFTraj = GetTrajectoryFromTrack(GSFTrajCollectionHandle.product(), ielectrons->gsfTrack().get()); 

edm::LogInfo("test") << "test";
    }

    if(theGSFTraj == nullptr) {
        edm::LogWarning("No GSF Traj") << " Did not find the GSF Trajectory for track";
    } else {
        edm::LogInfo("Ntup") << "Found Trajectory with hits " << theGSFTraj->foundHits();

        tmp_ele_signedEstimateSumPred = ReducedChi2AsymmetryNormalizedNValidHits(theGSFTraj, ielectrons->gsfTrack().get());

        tmp_ele_propagatorSignedEstimateSumPred = GetSignedChiPropagator(theGSFTraj, *ielectrons, iEvent, iSetup);    
    }
    ele_signedEstimateSumPred.push_back(tmp_ele_signedEstimateSumPred);
    ele_signedEstimateSumPred_A[counter] = tmp_ele_signedEstimateSumPred;
    ele_propagatorSignedEstimateSumPred.push_back(tmp_ele_propagatorSignedEstimateSumPred);


    setMomentum(myvector, ielectrons->p4());
    new (electrons[counter]) TLorentzVector(myvector);

    ele_echarge[counter] = ielectrons->charge(); 
    //
    ele_he[counter]        = ielectrons->hcalOverEcal(); //hadronicOverEm() ;
    ele_hebc[counter]      = ielectrons-> hcalOverEcalBc();

    // TrackCluster Matching
    ele_eseedpout[counter] = ielectrons->eSeedClusterOverPout();
    ele_ep[counter]        = ielectrons->eSuperClusterOverP() ;        
    ele_eseedp[counter]    = ielectrons->eSeedClusterOverP() ;         
    ele_eelepout[counter]  = ielectrons->eEleClusterOverPout() ;       
    //
    ele_pin_mode[counter]    = ielectrons->trackMomentumAtVtx().R() ; 
    ele_pout_mode[counter]   = ielectrons->trackMomentumOut().R() ; 
    ele_pTin_mode[counter]   = ielectrons->trackMomentumAtVtx().Rho() ; 
    ele_pTout_mode[counter]  = ielectrons->trackMomentumOut().Rho() ; 
    //
    ele_deltaetaseed[counter] = ielectrons->deltaEtaSeedClusterTrackAtCalo() ; 
    ele_deltaphiseed[counter] = ielectrons->deltaPhiSeedClusterTrackAtCalo() ;  
    ele_deltaetaele[counter]  = ielectrons->deltaEtaEleClusterTrackAtCalo() ;  
    ele_deltaphiele[counter]  = ielectrons->deltaPhiEleClusterTrackAtCalo() ; 
    ele_deltaetain[counter]   = ielectrons->deltaEtaSuperClusterTrackAtVtx();
    ele_deltaphiin[counter]   = ielectrons->deltaPhiSuperClusterTrackAtVtx();   

    // Shower Shape
    ele_sigmaietaieta[counter] = (ielectrons->showerShape()).sigmaIetaIeta; //  ielectrons->
    ele_sigmaetaeta[counter]   = (ielectrons->showerShape()).sigmaEtaEta ; //ielectrons->sigmaEtaEta() ;
    ele_sigmaiphiiphi[counter] = (ielectrons->showerShape()).sigmaIphiIphi;
    ele_e15[counter]           = (ielectrons->showerShape()).e1x5; // ;ielectrons->e1x5() ;
    ele_e25max[counter]        = (ielectrons->showerShape()).e2x5Max; // ;ielectrons->e2x5Max() ;
    ele_e55[counter]           = (ielectrons->showerShape()).e5x5; // ;ielectrons->e5x5() ;
    ele_r9[counter]            = (ielectrons->showerShape()).r9; 
    //ele_e1[counter]            = FIXME ;
    //ele_e33[counter]           = FIXME ;
    //
    // Old-Style ShowerShape
    edm::InputTag  reducedBarrelRecHitCollection("reducedEcalRecHitsEB");
    edm::InputTag  reducedEndcapRecHitCollection("reducedEcalRecHitsEE");
    //noZS::
    // noZS::EcalClusterLazyTools lazyToolsNoZS(iEvent, iSetup, reducedBarrelRecHitCollection, reducedEndcapRecHitCollection);

    // const auto & seedCluster = ielectrons->superCluster()->seed();
    // std::vector<float> vCov = lazyToolsNoZS.localCovariances(*seedCluster);
    // std::vector<float> Cov  = lazyToolsNoZS.covariances(*seedCluster);
    
    ele_oldsigmaetaeta[counter]   =  ielectrons->full5x5_sigmaEtaEta();    //( !edm::isNotFinite(Cov[0]) ) ? sqrt(Cov[0]) : 0;
    ele_oldsigmaietaieta[counter] =  ielectrons->full5x5_sigmaIetaIeta();   //( !edm::isNotFinite(vCov[0]) ) ? sqrt(vCov[0]) : 0;
    ele_oldsigmaiphiiphi[counter] =  ielectrons->full5x5_sigmaIphiIphi();   //( !edm::isNotFinite(vCov[2]) ) ? sqrt(vCov[2]) : 0;
    //ele_oldsigmaietaiphi[counter] = vCov[1]; // this is missing in the struct
    ele_oldr9[counter]              =  ielectrons->full5x5_r9();  //lazyToolsNoZS.e3x3(*seedCluster) / ielectrons->superCluster()->rawEnergy() ;
    ele_olde15[counter]           =  ielectrons->full5x5_e1x5(); //lazyToolsNoZS.e1x5(*seedCluster);
    ele_olde25max[counter]   =  ielectrons->full5x5_e2x5Max(); //lazyToolsNoZS.e2x5Max(*seedCluster);
    ele_olde55[counter]           =  ielectrons->full5x5_e5x5();       // lazyToolsNoZS.e5x5(*seedCluster);
    ele_oldhe[counter]             =  ielectrons->full5x5_hcalOverEcal();
    ele_oldhebc[counter]        =  ielectrons->full5x5_hcalOverEcalBc();
            // hcal stuff is not filled
    //electron.full5x5_setShowerShape(ss);
    //electron.full5x5_setSigmaIetaIphi(sigmaIetaIphi);


    // E/P combination
    ele_ecalE[counter]     = ielectrons->ecalEnergy();
    ele_ecalErr[counter]   = ielectrons->ecalEnergyError();
    ele_trackErr[counter]  = ielectrons->trackMomentumError();
    ele_combErr[counter]   = ielectrons->p4Error(GsfElectron::P4_COMBINATION);
    ele_PFcombErr[counter] = ielectrons->p4Error(GsfElectron::P4_PFLOW_COMBINATION);
    //cout << "Errors (ecal/track/p4comb/PFprcomb) :" <<ele_ecalErr[counter] <<" " << ele_trackErr[counter]<<" "<< ele_combErr[counter]<<" "<< ele_PFcombErr[counter] <<endl;
    
    //regression
    //ele_ecalRegressionEnergy[counter]  = ielectrons->ecalRegressionEnergy();
    //ele_ecalRegressionError[counter] = ielectrons->ecalRegressionError();
    // ele_ecalTrackRegressionEnergy[counter] = ielectrons->ecalTrackRegressionEnergy();
    //     ele_ecalTrackRegressionError[counter] = ielectrons->ecalTrackRegressionError();
    //     ele_ecalScale[counter] = ielectrons->ecalScale();               
    //     ele_ecalSmear[counter] = ielectrons->ecalSmear();                
    //     ele_ecalRegressionScale[counter] = ielectrons->ecalRegressionScale();     
    //     ele_ecalRegressionSmear[counter] = ielectrons->ecalRegressionSmear();     
    //     ele_ecalTrackRegressionScale[counter] = ielectrons->ecalTrackRegressionScale();
    //     ele_ecalTrackRegressionSmear[counter] = ielectrons->ecalTrackRegressionSmear();
    
    
    // 		cout << "REGRESSION " << ele_ecalRegressionEnergy[counter]<< " " << ele_ecalRegressionError[counter] << endl;
    //
    ele_mva[counter]   = 0; //ielectrons->mva() ;
    //
    if (ielectrons->isEB()) ele_isbarrel[counter] = 1 ; 
    else  ele_isbarrel[counter] = 0 ;
    if (ielectrons->isEE()) ele_isendcap[counter] = 1 ; 
    else  ele_isendcap[counter] = 0 ;
    if (ielectrons->isEBEtaGap()) ele_isEBetaGap[counter] = 1 ;  
    if (ielectrons->isEBPhiGap()) ele_isEBphiGap[counter] = 1 ;  
    if (ielectrons->isEEDeeGap()) ele_isEEdeeGap[counter] = 1 ;  
    if (ielectrons->isEERingGap()) ele_isEEringGap[counter] = 1 ;
    if (ielectrons->ecalDrivenSeed()) ele_isecalDriven[counter] = 1 ;
    if (ielectrons->trackerDrivenSeed()) ele_istrackerDriven[counter] = 1 ;
    //

    // -----------------------------------------------------------------
    // Tracking Variables
    // -----------------------------------------------------------------
    ele_lost_hits[counter]   = ielectrons->gsfTrack()->lost(); //numberOfLostHits();
    ele_valid_hits[counter] = ielectrons->gsfTrack()->found(); //numberOfValidHits() ;
    ele_gsfchi2[counter]    = ielectrons->gsfTrack()->normalizedChi2() ;
    
    //double hits = ielectrons->gsfTrack()->hitPattern().trackerLayersWithMeasurement();
    //cout << "hits = " << hits << " valid = " << ele_valid_hits[counter] << endl;
    ele_gsfhits[counter] = ielectrons->gsfTrack()->hitPattern().trackerLayersWithMeasurement();

    bool validKF=false;
    reco::TrackRef myTrackRef = ielectrons->closestCtfTrackRef();
    validKF = myTrackRef.isNonnull();
    ele_kfchi2[counter] = validKF ? myTrackRef->normalizedChi2() : 0 ; //ielectrons->track()->normalizedChi2() : 0 ;
    ele_kfhits[counter] = validKF ? myTrackRef->hitPattern().trackerLayersWithMeasurement() : -1. ; //ielectrons->track()->hitPattern().trackerLayersWithMeasurement() : 0 ;
    //
    // 		ele_dxyB[counter] = ielectrons->gsfTrack()->dxy(bs.position()) ;
    ele_dxy[counter]  = ielectrons->gsfTrack()->dxy() ;
    // 		ele_dzB[counter]  = ielectrons->gsfTrack()->dz(bs.position()) ;
    ele_dz[counter]   = ielectrons->gsfTrack()->dz() ;
    // 		ele_dszB[counter] = ielectrons->gsfTrack()->dsz(bs.position()) ;
    ele_dsz[counter]  = ielectrons->gsfTrack()->dsz() ;
    ele_dzPV[counter] = ielectrons->gsfTrack()->dz(pv->position());
    //
    // Isolation variables
    // ele_tkSumPt_dr03[counter]              = ielectrons->dr03TkSumPt() ;
    //     ele_ecalRecHitSumEt_dr03[counter]      = ielectrons->dr03EcalRecHitSumEt() ;
    //     ele_hcalDepth1TowerSumEt_dr03[counter] = ielectrons->dr03HcalDepth1TowerSumEt() ;
    //     ele_hcalDepth2TowerSumEt_dr03[counter] = ielectrons->dr03HcalDepth2TowerSumEt() ;
    //     ele_tkSumPt_dr04[counter]              = ielectrons->dr04TkSumPt() ;
    //     ele_ecalRecHitSumEt_dr04[counter]      = ielectrons->dr04EcalRecHitSumEt() ;
    //     ele_hcalDepth1TowerSumEt_dr04[counter] = ielectrons->dr04HcalDepth1TowerSumEt() ;
    //     ele_hcalDepth2TowerSumEt_dr04[counter] = ielectrons->dr04HcalDepth2TowerSumEt() ;
    //
    // Custom HCAL
    //  m_calotowers = new edm::Handle<CaloTowerCollection>() ;
    //     if (!iEvent.getByLabel("towerMaker",*m_calotowers)) //hcalTowers_
    //       { edm::LogError("ElectronHcalHelper::readEvent")<<"failed to get the hcal towers of label "<<hcalTowers_ ; }
    
    //     edm::Ref<pat::ElectronCollection> electronEdmRef(electronsCol, counter);
    //     double egHcalIsoConeSizeOutSmall = 0.3;
    //     double egHcalIsoConeSizeIn       = 0.0, egHcalIsoPtMin=0.0;
    //     int egHcalDepth1 = 1; 
    //     int egHcalDepth2 = 2;
    //hadDepth1Isolation03_  = new EgammaTowerIsolation(egHcalIsoConeSizeOutSmall,egHcalIsoConeSizeIn,egHcalIsoPtMin,egHcalDepth1,m_calotowers->product()) ;
    //hadDepth2Isolation03_  = new EgammaTowerIsolation(egHcalIsoConeSizeOutSmall,egHcalIsoConeSizeIn,egHcalIsoPtMin,egHcalDepth2,m_calotowers->product()) ;
    //double hcalDepth1TowerSumEt03 = hadDepth1Isolation03_->getTowerEtSum(&(*electronEdmRef)); //ielectrons)); //electronRef));
    //double hcalDepth2TowerSumEt03 = hadDepth2Isolation03_->getTowerEtSum(&(*electronEdmRef)); // electronEdmRefielectrons)); //electronRef));
    //ele_HCALFullConeSum[counter]  = hcalDepth1TowerSumEt03+hcalDepth2TowerSumEt03;
    //&((*EleHandle)[i])
    //
    // Conversion Rejection
    ele_conv_dcot[counter]   = ielectrons->convDist(); //userFloat("dcot");
    ele_conv_dist[counter]   = ielectrons->convDcot(); //ielectrons->userFloat("dist");
    ele_conv_radius[counter] = ielectrons->convRadius();
    
    // FIXME: Always returns 0 :(
    //cout << " dcot = " << ielectrons->userFloat("dcot") << " dist = " << ielectrons->userFloat("dist") << endl;
    
    ele_expected_inner_hits[counter] = ielectrons->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS);
    //ielectrons->gsfTrack()->trackerExpectedHitsInner().numberOfHits();
    
    
    bool vtxFitConversion = ConversionTools::hasMatchedConversion(*ielectrons, conversions_h, beamSpot.position());
    ele_vtxconv[counter] = vtxFitConversion;
    
    // 
    
    // -----------------------------------------------------------------
    // PF Isolation for electrons in the cone 0.4
    // -----------------------------------------------------------------
    ele_pfChargedHadIso[counter]   = (ielectrons->pfIsolationVariables()).sumChargedHadronPt ; //chargedHadronIso();
    ele_pfNeutralHadIso[counter]   = (ielectrons->pfIsolationVariables()).sumNeutralHadronEt ; //neutralHadronIso();
    ele_pfPhotonIso[counter]       = (ielectrons->pfIsolationVariables()).sumPhotonEt; //photonIso();
    //
    ele_pfChargedIso[counter]      = (ielectrons->pfIsolationVariables()).sumChargedParticlePt;
    ele_pfSumPUIso[counter]        = (ielectrons->pfIsolationVariables()).sumPUPt;
    //ele_pfCombRelIso[counter]      = LeptonIsoHelper::combRelIsoPF(lepton_setup, lepton_setup, _PU_Elerho, *ielectrons);
    //FIXME these two for applying the beta corrections
    //ele_pfChargedHadPUIso[counter]  = ielectrons->chargedAllIso();
    //ele_pfChargedHadPUIso[counter]  = ielectrons->puChargedHadronIso();
    
    // -----------------------------------------------------------------
    // SIP3D
    // -----------------------------------------------------------------
    //default values for IP3D
    double ele_ip3D     = -999.0; 
    double ele_ip3D_err = 999.; //fMVAVar_ip3dSig = 0.0;
    double d0_corr = 999;
    double d0_err   = 999;

    if (ielectrons->gsfTrack().isNonnull()) {
      const double gsfsign   = ( (-ielectrons->gsfTrack()->dxy(pv->position()))   >=0 ) ? 1. : -1.;
      
      const reco::TransientTrack &tt = thebuilder.build(ielectrons->gsfTrack()); //transientTrackBuilder.build(ielectrons->gsfTrack()); 
      const std::pair <bool,Measurement1D> &ip3dpv =  IPTools::absoluteImpactParameter3D(tt, *pv);
      std::pair<bool,Measurement1D> result   = IPTools::absoluteTransverseImpactParameter(tt, * pv);
      
      if (ip3dpv.first) {
	double ip3d = gsfsign*ip3dpv.second.value();
	double ip3derr = ip3dpv.second.error();  
	
	ele_ip3D     = ip3d; 
	ele_ip3D_err = ip3derr;
	// fMVAVar_ip3dSig = ip3d/ip3derr;
      } // if ip3dpv.first

      if(result.first) {
	d0_corr = result.second.value();
	d0_err   = result.second.error();
      } // if d0
    } // if gsftrack.isNonnull
    
    ele_IP[counter]      = ele_ip3D; //fabs(ielectrons->dB(pat::Electron::PV3D));
    ele_IPError[counter] = ele_ip3D_err; //ielectrons->edB(pat::Electron::PV3D);	
    double ele_sip = -999; if(ele_ip3D_err!=0) ele_sip = ele_ip3D / ele_ip3D_err;
    ele_SIP[counter]     = ele_sip; //ele_IP[counter]/ele_IPError[counter];
    
    ele_d0[counter] = d0_corr;
    ele_d0err[counter] = d0_err;

    // -----------------------------------------------------------------
    // Get SuperCluster Informations
    // -----------------------------------------------------------------
    //cout << " SuperCluster "<< endl;
    //	if(ielectrons->ecalDrivenSeed()) {
    reco::SuperClusterRef sclRef = ielectrons->superCluster();
    //cout << " SuperClusterRef" << endl;
    ////math::XYZPoint sclPos        = ielectrons->superClusterPosition();
    //      cout << " pflow" << endl;
    
    //if (!ielectrons->ecalDrivenSeed() && ielectrons->trackerDrivenSeed()) 
    //sclRef = ielectrons->pflowSuperCluster();
    
    double R  = TMath::Sqrt(sclRef->x()*sclRef->x() + sclRef->y()*sclRef->y() +sclRef->z()*sclRef->z());
    double Rt = TMath::Sqrt(sclRef->x()*sclRef->x() + sclRef->y()*sclRef->y());
    ele_sclRawE[counter]  = sclRef->rawEnergy() ;
    
    ele_sclE[counter]     = sclRef->energy() ;
    // 		ele_sclE[counter]     = ielectrons->correctedEcalEnergy();  //for 5XY
    ele_sclEt[counter]    = sclRef->energy()*(Rt/R) ;
    ele_sclEta[counter]   = sclRef->eta() ;
    ele_sclPhi[counter]   = sclRef->phi() ;
    ele_sclNclus[counter] = sclRef->clustersSize();
    
    //cout << " etawidth = " << sclRef->etaWidth() << endl;

    ele_sclphiwidth[counter] = sclRef->phiWidth();
    ele_scletawidth[counter] = sclRef->etaWidth();

    // --------------
    // Sub-Clusters
    // --------------
    int countersub = 0;
    //eSubClusters_ = 0.;
    // Store subclusters
    reco::CaloCluster_iterator itscl  = sclRef->clustersBegin();
    reco::CaloCluster_iterator itsclE = sclRef->clustersEnd();
    
    for(; itscl < itsclE ; ++itscl) {
      bool isseed = false;
      if((*itscl)==ielectrons->superCluster()->seed()) isseed=true; // continue; // skip seed cluster
      //theBasicClusters_.push_back(&(**itscl));  
      //eSubClusters_ += (*itscl)->energy();
      ele_sclsubE[counter][countersub]      = (*itscl)->energy();
      ele_sclsubEta[counter][countersub]    = (*itscl)->eta();
      ele_sclsubPhi[counter][countersub]    = (*itscl)->phi();
      ele_sclsubisseed[counter][countersub] = isseed;
      //N subclusters?->sclNclus... to be checked
      countersub++;
    }
    // sort subclusters

    //	} // if ECAL driven
    
    //cout << " Etawidth = " <<  ele_scletawidth[counter] << endl;
    //cout << "" << endl;
    
    // -----------------------------------------------------------------
    // Get PreShower Informations
    // -----------------------------------------------------------------
    ele_psE[counter] = sclRef->preshowerEnergy();
    //T_Elec_PreShowerOverRaw->push_back((eleIt->superCluster()->rawEnergy()>0) ? eleIt->superCluster()->preshowerEnergy() / eleIt->superCluster()->rawEnergy() : -1)

    // -----------------------------------------------------------------
    //fbrem
    // -----------------------------------------------------------------
    ele_fbrem[counter]      = ielectrons->fbrem();
    ele_trackfbrem[counter] = ielectrons->trackFbrem();
    ele_SCfbrem[counter]    = ielectrons->superClusterFbrem();
    ele_pfSCfbrem[counter]  = ielectrons->pfSuperClusterFbrem(); // Should be identical to the previous one...
    ele_eClass[counter]     = ielectrons->classification() ;
    ele_nbrem[counter]      = ielectrons->numberOfBrems();
    
    // FOR MVA...
    //  ele_mvafbrem[counter] = ielectrons->fbrem();
    //     ele_mvadetain[counter] = ielectrons->deltaEtaSuperClusterTrackAtVtx();
    //     ele_mvadphiin[counter] = ielectrons->deltaPhiSuperClusterTrackAtVtx();
    //     ele_mvasieie[counter] = ielectrons->sigmaIetaIeta();
    //     ele_mvahoe[counter] = ielectrons->hcalOverEcal();
    //     ele_mvaeop[counter] = ielectrons->eSuperClusterOverP();
    //     ele_mvae1x5e5x5[counter] = (ielectrons->e5x5()) !=0. ? ielectrons->e1x5()/ielectrons->e5x5() : -1. ;
    //     ele_mvaeleopout[counter] = ielectrons->eEleClusterOverPout();
    //     bool validKF= (ielectrons->track().isNonnull());
    //     ele_mvakfchi2[counter] = validKF ? ielectrons->track()->normalizedChi2() : 0 ;
    //     ele_mvakfhits[counter] = validKF ? ielectrons->track()->hitPattern().trackerLayersWithMeasurement() : 0 ;
    //     ele_mvamishits[counter] = ielectrons->gsfTrack()->trackerExpectedHitsInner().numberOfLostHits();
    //     ele_mvadist[counter] = ielectrons->convDist();
    //     ele_mvadcot[counter] = ielectrons->convDcot();
    //     ele_mvaeta[counter] = ielectrons->eta();
    //     ele_mvapt[counter] = ielectrons->pt();
    //     ele_mvaecalseed[counter] = ielectrons->ecalDrivenSeed();

    // -----------------------------------------------------------------
    // Electron ID electronsCol
    // -----------------------------------------------------------------
    edm::Ref<reco::GsfElectronCollection> electronRef(electronsCol, counter); //electronsCollection,i); i++; //reference to the electron
    //cout << "MVA = " << mapMVA[electronRef] << endl;
//    ele_mvaphys14[counter] = mapMVA_phys14[electronRef] ;
//    ele_mvaphys14fix[counter] = mapMVA_phys14fix[electronRef] ;

    //T_Elec_MVAoutput->push_back(mapMVA[electronRef]);

    ++counter;
  } // for loop on gsfelectrons

  if(counter>49) { ele_N = 50; cout << "Number of electrons>49, electrons_N set to 50" << endl;}
  
} // end of FillElectrons

// ====================================================================================
void Ntuplizer::FillMET (const edm::Event& iEvent, const edm::EventSetup& iSetup)
// ====================================================================================
{
	
  // caloMET object (negative vector sum of calorimeter towers)
  //edm::Handle< edm::View<reco::CaloMET> > caloMEThandle;
  //iEvent.getByLabel("met", caloMEThandle);
  
  // MET object that corrects the basic calorimeter MET for muons
  // edm::Handle< edm::View<reco::CaloMET> > muCorrMEThandle;
  //   iEvent.getByLabel("corMetGlobalMuons", muCorrMEThandle);
  
  // MET object that corrects the basic calorimeter MET for muons and tracks
  //edm::Handle< edm::View<reco::MET> > tcMEThandle;
  //iEvent.getByLabel("tcMet", tcMEThandle);
  
  // MET object built as the (negative) vector sum of all particles (PFCandidates) reconstructed in the event
  // 	edm::Handle< edm::View<pat::MET> > pfMEThandle;
  // 	iEvent.getByLabel("patMETs", pfMEThandle);
  
  //edm::Handle< edm::View<cmg::BaseMET> > pfMEThandle;
  edm::Handle< edm::View<reco::PFMET> > pfMEThandle;
  iEvent.getByLabel("pfMet", pfMEThandle);
  
  
  // PFMET
  _met_pf_et  = (pfMEThandle->front() ).et();
  _met_pf_px  = (pfMEThandle->front() ).px();
  _met_pf_py  = (pfMEThandle->front() ).py();
  _met_pf_phi = (pfMEThandle->front() ).phi();
  _met_pf_set = (pfMEThandle->front() ).sumEt();
  _met_pf_sig = (pfMEThandle->front() ).mEtSig();
  
} // end of Fill MET



// ====================================================================================
void Ntuplizer::FillTruth(const edm::Event& iEvent, const edm::EventSetup& iSetup)
// ====================================================================================
{
  
  //edm::Handle< GenEventInfoProduct > HepMCEvt;
  //iEvent.getByLabel(MCTag_, HepMCEvt);
  //if(HepMCEvt->hasBinningValues()) _MC_pthat = (HepMCEvt->binningValues())[0];
  //else  _MC_pthat = 0.0;
  
  edm::Handle<View<Candidate> > genCandidatesCollection;
  //iEvent.getByLabel("prunedGen", genCandidatesCollection);
  iEvent.getByLabel("genParticles", genCandidatesCollection); //genParticlesPruned
  
  TClonesArray &MC_gen_V               = *_m_MC_gen_V;
  TClonesArray &MC_gen_Higgs           = *_m_MC_gen_Higgs;
  //cout << " photon" << endl;
  //TClonesArray &MC_gen_photons         = *_m_MC_gen_photons;
  TClonesArray &MC_gen_leptons         = *_m_MC_gen_leptons;
  TClonesArray &MC_gen_leptons_status2 = *_m_MC_gen_leptons_status2;
  TClonesArray &MC_gen_leptons_status1 = *_m_MC_gen_leptons_status1;
  
  int counter             = 0;
  int counter_higgs       = 0;
  int counter_daughters   = 0;
  int counter_lep_status2 = 0;
  int counter_lep_status1 = 0;
  //int counter_photon      = 0;
  
  // ----------------------------
  //      Loop on particles
  // ----------------------------
  for( View<Candidate>::const_iterator p = genCandidatesCollection->begin();p != genCandidatesCollection->end(); ++ p ) {
    
    // %%%%%%%%%%%%%%%%%%
    // If Higgs
    // %%%%%%%%%%%%%%%%%%
    if (p->pdgId() == 25 && p->status()==3) {
      setMomentum (myvector,p->p4());
      // 		  cout << "Higgs PdgId=" << p->pdgId() << " Higgs status=" << p->status() << " Mass=" << myvector.M() << endl;
      new (MC_gen_Higgs[counter_higgs]) TLorentzVector(myvector);
      counter_higgs++;
    } // if Higgs
    
    // %%%%%%%%%%%%%%%%%%
    // If Leptons from Z
    // %%%%%%%%%%%%%%%%%%
    if(fabs(p->pdgId())==11) { // || fabs(p->pdgId())==13 ||  fabs(p->pdgId())==15) {
      
      // cout << "Status pdgid = " << fabs(p->pdgId()) << " status = "<< p->status()  << " Pt = "<< p->pt() << " Eta = " << p->eta() <<  endl;
      
      //       if(p->status()==1) {
      // 	cout << "Status1 pdgid = " << fabs(p->pdgId()) << " status = "<< p->status()  << " Pt = "<< p->pt() << " Eta = " << p->eta() <<  endl;
      // 	cout << " Nmother = " << p->numberOfMothers() << endl;
      // 	for(unsigned int i=0;i<p->numberOfMothers();i++) {
      // 	  cout << " mother pdgid = " << fabs(p->mother(i)->pdgId()) << " status = " << p->mother(i)->status() << endl;
      
      // 	  for(unsigned int j=0;j<p->mother(i)->numberOfMothers();j++) {
      // 	    cout << " Grandmother pdgid = " << fabs(p->mother(i)->mother(j)->pdgId()) << " status = " << p->mother(i)->mother(j)->status() << endl;
      // 	  }
      // 	} // for loop on mothers
      // 	//cout << "" << endl;
      //       } // if status==1
      
      if(p->status()==1) {
	//if(p->numberOfMothers()>0) { // Need a mother...
	
	setMomentum(myvector, p->p4());
	new (MC_gen_leptons_status1[counter_lep_status1]) TLorentzVector(myvector);
	_MC_gen_leptons_status1_pdgid[counter_lep_status1] = p->pdgId();
	
	
	// Let's look at the mothers ... [works really only for electrons...=>eid studies]
	if(p->numberOfMothers()>0) {
	  //bool foundW = false;
	  //if (p.numberOfMothers()==0) return foundW;
	  //const reco::Candidate *mother= p->moter();
	  
	  const reco::Candidate  *Mom = p->mother();
	  while(Mom!=0) {
	    if(fabs(Mom->pdgId()) == 23 || fabs(Mom->pdgId()) == 24) { 
	      //cout << "W or Z !" << endl; 
	      _MC_gen_leptons_status1_FromWZ[counter_lep_status1] = 1;
	      break;}
	    else if(fabs(Mom->pdgId()) == 15) { 
	      //cout << " Taus  !" << endl; 
	      _MC_gen_leptons_status1_FromTaus[counter_lep_status1] = 1;
	      break;}
	    else if(fabs(Mom->pdgId()) > 50 ) { 
	      //cout << " Non Prompt  ! pdgID = " << fabs(Mom->pdgId()) << " Nmother = " <<  Mom->numberOfMothers() << endl; 
	      _MC_gen_leptons_status1_FromNonPrompt[counter_lep_status1] = 1;

	      const reco::Candidate  *GrandMom = Mom; //->mother();
	      while(GrandMom!=0) {
		for(int unsigned jj=0;jj<GrandMom->numberOfMothers();jj++) {
		  if(fabs(GrandMom->mother(jj)->pdgId()) == 4 || fabs(GrandMom->mother(jj)->pdgId()) == 5) { 
		    //cout << " B or C (loop)  !" << endl; 
		    _MC_gen_leptons_status1_FromBC[counter_lep_status1] = 1;
		    break;}
		} // for loop on grand-mothers, when they are more than 1...
		//cout << " Non Prompt GM = " << fabs(GrandMom->pdgId()) << endl;
		if(fabs(GrandMom->pdgId()) == 4 || fabs(GrandMom->pdgId()) == 5) { 
		  //cout << " B or C  !" << endl; 
		  _MC_gen_leptons_status1_FromBC[counter_lep_status1] = 1;
		  break;}
		GrandMom = GrandMom ->mother();
	      } // while loop on grand mom
	      
	      break;} // if pdgid>50
	    
	    Mom = Mom ->mother();
	  } // while loop on mothers
	  
	  
	  // 	  while (p->numberOfMothers()>0) {
	  // 	    const reco::Candidate  *Mom = ->mother();
	  // 	    if ((fabs(MomPart->pdgId())>=22)&&(fabs(MomPart->pdgId())<=24)){
	  // 	  }// while loop on mother
	  
	  // 	  const reco::Candidate  *part = (p.mother());
	  // 	  // loop on the mother particles to check if is has a W has mother
	  // 	  while ((part->numberOfMothers()>0)) {
	  // 	    const reco::Candidate  *MomPart =part->mother();
	  // 	    if ((fabs(MomPart->pdgId())>=22)&&(fabs(MomPart->pdgId())<=24)){
	  // 	      foundW = true;
	  //             break;
	  //         }
	  //         part = MomPart;
	  
	  
	  
	  //if(p->mother(0)->pdgId()== p->pdgId()) {
	  
	  //}// if pdgid
	  
	} // if mother
	
	counter_lep_status1++;
      } // if status 1
      //cout << "" << endl;

      if(p->status()==3) {
	if(p->numberOfMothers()>0) { // Need a mother...
	  if(p->mother(0)->pdgId()==23) {  // If Mother is a Z 
	    
	    //cout << "number of daughters = " << p->numberOfDaughters() << " mother id = " << p->pdgId() << endl;
	    
	    if(p->numberOfDaughters()>0) { // Need a daughter...
	      
	      //cout << " status of daughter = " << p->daughter(0)->status() << " pdgid = " << p->daughter(0)->pdgId() << endl;
	      
	      // Status 2 Leptons
	      if(p->daughter(0)->pdgId()==p->pdgId() && p->daughter(0)->status()==2) { // if daughter is lepton & status 2
		setMomentum(myvector, p->daughter(0)->p4());
		new (MC_gen_leptons_status2[counter_lep_status2]) TLorentzVector(myvector);
		_MC_gen_leptons_status2_pdgid[counter_lep_status2] = p->daughter(0)->pdgId();
		counter_lep_status2++;
		
		//cout << "dod = " << p->daughter(0)->daughter(0)->pdgId() << " status = " << p->daughter(0)->daughter(0)->status() << endl;
		
		// 		if(p->daughter(0)->daughter(0)->status()==2) {
		// 		  //cout << "Ndodod = " << p->daughter(0)->daughter(0)->numberOfDaughters() << endl;
		// 		  for(unsigned int i=0;i<p->numberOfDaughters();i++) {
		// 		    cout << " Dodod pdgid = " << 
		// 		  }
		// 		}
		
		// 	// Status 1 Leptons, from Status 2
		// 		if(p->daughter(0)->daughter(0)->pdgId()==p->pdgId() && p->daughter(0)->daughter(0)->status()==1) {
		// 		  setMomentum(myvector, p->daughter(0)->daughter(0)->p4());
		// 		  new (MC_gen_leptons_status1[counter_lep_status1]) TLorentzVector(myvector);
		// 		  _MC_gen_leptons_status1_pdgid[counter_lep_status1] = p->daughter(0)->daughter(0)->pdgId();
		// 		  counter_lep_status1++;
		// 		} // if status 1 from status 2 from status 3
		// 		else {
		// 		  if(fabs(p->pdgId())==11)
		// 		    //cout << "(status2) mother   id ? " << p->daughter(0)->pdgId() << " status ? " << p->daughter(0)->status() << endl;
		// 		    //cout << "(status2) daughter id ? " << p->daughter(0)->daughter(0)->pdgId() << " status ? " << p->daughter(0)->daughter(0)->status() << endl;
		// 		}
		
	      } // if Daughter Status = 2
	      
	      //  // Status 1 Leptons, from Status 3
	      // 	      if(p->daughter(0)->pdgId()==p->pdgId() && p->daughter(0)->status()==1) {
	      // 		setMomentum(myvector, p->daughter(0)->p4());
	      // 		new (MC_gen_leptons_status1[counter_lep_status1]) TLorentzVector(myvector);
	      // 		_MC_gen_leptons_status1_pdgid[counter_lep_status1] = p->daughter(0)->pdgId();
	      // 		counter_lep_status1++;
	      // 	      } // if Daughters Status =1 (from status 3)
	      // 	      else {
	      // 		if(fabs(p->pdgId())==11)
	      // 		  //cout << "(status3) daughter id ? " << p->daughter(0)->pdgId() << " status ? " << p->daughter(0)->status() << endl;
	      // 	      }
	      
	      
	    } // Need a daughter
	  } // if MOther = Z
	} // if Nmother (status3) >0
      } // if hard scatter electron (status3)
    } // if leptons
    
    //  if(p->status()==2 && p->mother(i)->pdgId() ==  p->pdgId() && p->mother(i)->status()==3 
    // 	 &&  p->mother(i)->mother(0)==23) { 
    
    //       } // status 2 && mother, same pdgid and status 3
    
    
    //       if(p->status()==1 && p->mother(i)->pdgId() == p->pdgId() 
    // 	 && (p->mother(i)->status()==2 && p->mother(i)->mother(0)->pdgId()== p->pdgId())
    // || p->mother(i)->status()==3) {
    
    //       } // if 
    
    //     } // if leptons
    
    
    //       if(p->status()==1) { 
    // 	cout << "Ele Status 1, Nmother =  " << p->numberOfMothers() << " Pt = " << p->p4().Pt() << endl;
    // 	for(unsigned int i=0;i<p->numberOfMothers();i++) {
    // 	  cout << " Mother1 pdgid = " << p->mother(i)->pdgId() << " status = " << p->mother(i)->status() << endl;
    // 	} //mother
    //        } // status
    
    //       if(p->status()==2) { 
    // 	cout << "Ele Status 2, Nmother =  " << p->numberOfMothers() << " Pt = " << p->p4().Pt() << endl;
    // 	for(unsigned int i=0;i<p->numberOfMothers();i++) {
    // 	  cout << " Mother2 pdgid = " << p->mother(i)->pdgId() << " status = " << p->mother(i)->status() << endl;
    // 	} //mother
    //       } // status
    
    //       if(p->status()==3) { 
    // 	cout << "Ele Status 3, Nmother =  " << p->numberOfMothers() << " Pt = " << p->p4().Pt() << endl;
    // 	for(unsigned int i=0;i<p->numberOfMothers();i++) {
    // 	  cout << " Mother3 pdgid = " << p->mother(i)->pdgId() << " status = " << p->mother(i)->status() << endl;
    // 	} //mother
    //       } // status
    
    //     } // if electron
    
    // %%%%%%%%%%%%%%%%%%
    //     If W or Z
    // %%%%%%%%%%%%%%%%%%
    if (p->pdgId() == 23 || fabs(p->pdgId())==24) {
      
      //cout << "Z status= " << p->status() << " mass = " << p->mass() << endl;
      
      //  if(p->status()==2) {
      // 	//cout << "status2,NB daugthers ? " << p->numberOfDaughters() << endl;
      // 	for(unsigned int i=0;i<p->numberOfDaughters();i++) {
      // 	  //cout << "Status2 id daughters = " << p->daughter(i)->pdgId() << " Status = " << p->daughter(i)->status()  << " Pt = " << p->daughter(i)->p4().Pt() << endl;
      // 	} // loop on daughters
      //       }
      
      bool pass_status=false;
      if(p->status()==3 && ispythia6_==true) pass_status=true; //{ // if PYTHIA6
      if(p->status() && p->numberOfDaughters()==2 && ispythia6_==false) pass_status=true; //{ // if PYTHIA 8
      // UGLY !!!! To uniformize at some point...
      
      if(pass_status) {
	// Fill truth W,Z
	setMomentum (myvector,p->p4());
	new (MC_gen_V[counter]) TLorentzVector(myvector);
	_MC_gen_V_pdgid[counter] = p->pdgId();
	
	//size_t nfirstdau = p->numberOfDaughters();
	
	//cout << "Z status= " << p->status() << " mass = " << p->mass() << " n filles = " <<  p->numberOfDaughters() << endl;

	// Loop on daughters
	for(unsigned int i=0;i<p->numberOfDaughters();i++) {
	  //cout << "fille " << i << " id = " << fabs(p->daughter(i)->pdgId()) << endl;
	  bool islep = false;
	  if(fabs(p->daughter(i)->pdgId())==11) { _MC_flavor[counter] = 0; islep=true;} // electron
	  if(fabs(p->daughter(i)->pdgId())==13) { _MC_flavor[counter] = 1; islep=true;} // muon
	  if(fabs(p->daughter(i)->pdgId())==15) { _MC_flavor[counter] = 2; islep=true;} // taus
	  
	  //cout << "Status3 id daughters = " << p->daughter(i)->pdgId() << " Status = " << p->daughter(i)->status()  << " Pt = " << p->daughter(i)->p4().Pt() << endl;
	  
	  
	  //  for(unsigned int j=0;j<p->daughter(i)->numberOfDaughters();j++) {
	  // 	    //cout << "DoD status = " <<  p->daughter(i)->daughter(j)->status() << " Pt = " << p->daughter(i)->daughter(j)->p4().Pt() << endl;
	  
	  // 	    for(unsigned int k=0;k<p->daughter(i)->daughter(j)->numberOfDaughters();k++) {
	  // 	      //cout << "DoDoD status = " <<  p->daughter(i)->daughter(j)->daughter(k)->status() << " Pt = " << p->daughter(i)->daughter(j)->daughter(k)->p4().Pt() << endl;
	  // 	    } // loop on k
	  
	  // 	  } //
	  
	  if(islep) { // p->daughter(i)->status()==1) { ?!
	    setMomentum(myvector, p->daughter(i)->p4());
	    new (MC_gen_leptons[counter_daughters]) TLorentzVector(myvector);
	    _MC_gen_leptons_pdgid[counter_daughters] = p->daughter(i)->pdgId();
	    
	    counter_daughters++;
	  } // if is lepton
	} // for loop on daughters
	counter++;
      } // if status stable
    } // if W or Z
    
    // // -------------
    // 		//   Photons
    // 		// -------------
    
    // 		if(fabs(p->pdgId())==22) { 
    // 			if(p->status()==1) {
    
    // 				//cout << " if photon" << endl;
    
    // 				setMomentum(myvector, p->p4());
    // 				new (MC_gen_photons[counter_photon]) TLorentzVector(myvector);
    
    // 				const Candidate * gen_photon = 0; 
    // 				//const GenParticle * gen_photon2 = 0; 
    // 				gen_photon  = &*p;
    // 				//gen_photon2 =  &*p;
    // 				//cout << "is fsr ?" << endl;
    // 				//const reco::GenParticle* photon = getParent(gen_photon);
    // 				bool fsr = isFSR(gen_photon); //const reco::GenParticle* genLep)
    
    // 				//cout << "fsr = " << fsr << endl;
    // 				_MC_gen_photons_isFSR[counter_photon] = fsr;
    // 				//" pT = " <<p->p4().Pt()  << endl;
    
    // 				//_MC_photon_isFSR[counter_photon] = isFSR(p, 
    
    // 				counter_photon++;
    // 			} // if status 1
    // 		} // if photon
    
    
  } // for loop on particles
  
  
  
  
  
  // 	// -------------------------
  // 	// Only 4e events
  // 	// -------------------------
  // 	if(_MC_flavor[0]==0 && _MC_flavor[1]==0) {
// 		// ----------------------------
// 		//      Loop on particles
// 		// ----------------------------
// 		for( View<Candidate>::const_iterator p = genCandidatesCollection->begin();p != genCandidatesCollection->end(); ++ p ) {
			
// 			if(fabs(p->pdgId())==11) { 
				
// 				if(p->status()==1) {
					
// 					//cout << "N mother = " << p->numberOfMothers() << endl;
// 					for(unsigned int i=0;i<p->numberOfMothers();i++) {
// 						//cout << " mother pdgid = " << p->mother(i)->pdgId() << " status = " << p->mother(i)->status() << endl;
// 					} // for loop on mothers
// 				} // if status 1
				
// 				if(p->status()==2) {
					
// 					//cout << "N daughters (status2) = " << p->numberOfDaughters() << endl;
// 					for(unsigned int i=0;i<p->numberOfDaughters();i++) {
// 						//cout << " daughter pdgid = " << p->daughter(i)->pdgId() << " status = " << p->daughter(i)->status() << endl;
// 					} // for loop on mothers
// 				} // if status 2
				
// 				if(p->status()==3) {
					
// 					//cout << "N daughters (status3) = " << p->numberOfDaughters() << endl;
// 					for(unsigned int i=0;i<p->numberOfDaughters();i++) {
// 						//cout << " daughter pdgid = " << p->daughter(i)->pdgId() << " status = " << p->daughter(i)->status() << endl;
// 					} // for loop on mothers
// 				} // if status 3
				
				
				
// 			} // if electron
// 		} // for loop on gen particles
		
// 	} // if 4e event
	
	
	
	
	
	
	
	
} // end of FillTruth




// ------------ method called once each job just after ending the event loop  ------------
void 
Ntuplizer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
void 
Ntuplizer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
Ntuplizer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
Ntuplizer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
Ntuplizer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/


 // ====================================================================================================
void Ntuplizer::setMomentum(TLorentzVector &myvector, const LorentzVector & mom)
// ====================================================================================================
{
  
  myvector.SetPx (mom.Px());
  myvector.SetPy (mom.Py());
  myvector.SetPz (mom.Pz());
  myvector.SetE (mom.E());
  
}


 // ====================================================================================================
void Ntuplizer::Init()
// ====================================================================================================
{

  _PU_N = 0;
  
  _vtx_N = 0;

  ele_N = 0;
  
  _met_pf_et  = 0.;
  _met_pf_px  = 0.; 
  _met_pf_py  = 0.; 
  _met_pf_phi = 0.; 
  _met_pf_set = 0.; 
  _met_pf_sig = 0.; 
  
  for (int i = 0 ; i < 50 ; ++i) {
    //electrons
    ele_echarge[i] = 0 ;
    ele_he[i] = 0 ;  
    ele_hebc[i] = 0 ;  
    ele_oldhe[i] = 0 ;  
    ele_oldhebc[i] = 0 ;  
    ele_eseedpout[i] = 0 ;  
    ele_ep[i] = 0 ;  
    ele_eseedp[i] = 0 ;  
    ele_eelepout[i] = 0 ;        
    ele_deltaetaseed[i] = 0 ;  
    ele_deltaetaele[i] = 0 ;  
    ele_deltaphiseed[i] = 0 ;  
    ele_deltaphiele[i] = 0 ;  
    ele_deltaetain[i] = 0 ;  
    ele_deltaphiin[i] = 0 ; 
    ele_sigmaietaieta[i] = 0 ;  
    ele_sigmaiphiiphi[i] = 0 ;  
    ele_sigmaetaeta[i] = 0 ;  
    ele_e15[i] = 0 ;  
    ele_e25max[i] = 0 ;  
    ele_e55[i] = 0 ;  
    ele_e1[i] = 0 ;  
    ele_r9[i] = 0;
    //
    ele_oldsigmaetaeta[i]   = 0;
    ele_oldsigmaietaieta[i] = 0;
    ele_oldsigmaiphiiphi[i] = 0;
    ele_oldsigmaietaiphi[i] = 0; 
    ele_oldr9[i]            = 0; 
    ele_olde15[i]           = 0; 
    ele_olde25max[i]        = 0; 
    ele_olde55[i]           = 0; 
    //
    //ele_e33[i] = 0 ;  
    //ele_e2overe9[i] = 0 ; 
    ele_pin_mode[i] = 0 ;  
    ele_pout_mode[i] = 0 ;  
    ele_pTin_mode[i] = 0 ;  
    ele_pTout_mode[i] = 0 ;  
    //
    ele_fbrem[i] = 0 ;  
    ele_SCfbrem[i] = 0 ;  
    ele_pfSCfbrem[i] = 0 ;  
    ele_trackfbrem[i] = 0 ;  
    ele_nbrem[i] = 0;
    //
    ele_mva[i] = 0 ; 
    ele_isbarrel[i] = 0 ;  
    ele_isendcap[i] = 0 ;  
    ele_isEBetaGap[i] = 0 ;  
    ele_isEBphiGap[i] = 0 ;  
    ele_isEEdeeGap[i] = 0 ;  
    ele_isEEringGap[i] = 0 ; 
    ele_isecalDriven[i] = 0 ;  
    ele_istrackerDriven[i] = 0 ; 
    ele_eClass[i] = 0 ;
    ele_valid_hits[i] = 0 ; 
    ele_lost_hits[i] = 0 ; 
    ele_gsfchi2[i] = 0 ; 
    ele_dxyB[i] = 0 ; 
    ele_dxy[i] = 0 ; 
    ele_dzB[i] = 0 ; 
    ele_dz[i] = 0 ; 
    ele_dszB[i] = 0 ; 
    ele_dsz[i] = 0 ;              
    //  ele_tkSumPt_dr03[i] = 0 ;  
    //     ele_ecalRecHitSumEt_dr03[i] = 0 ;  
    //     ele_hcalDepth1TowerSumEt_dr03[i] = 0 ;  
    //     ele_hcalDepth2TowerSumEt_dr03[i] = 0 ; 
    //     ele_tkSumPt_dr04[i] = 0 ;  
    //     ele_ecalRecHitSumEt_dr04[i] = 0 ;  
    //     ele_hcalDepth1TowerSumEt_dr04[i] = 0 ;  
    //     ele_hcalDepth2TowerSumEt_dr04[i] = 0 ; 
    ele_conv_dcot[i] = 0 ;
    ele_conv_dist[i] = 0 ;
    ele_conv_radius[i] = 0;
    ele_expected_inner_hits[i] = 0;
    ele_vtxconv[i] = 0;
    //
    ele_pfChargedHadIso[i]   = 0;
    ele_pfNeutralHadIso[i]   = 0;
    ele_pfPhotonIso[i]       = 0;

    ele_pfChargedIso[i] = 0;
    ele_pfSumPUIso[i]   = 0;
    //ele_pfChargedHadPUIso[i] = 0;
    //ele_pfCombRelIso[i] = 0;
    
    ele_IP[i] = 0 ; 
    ele_IPError[i] = 0 ; 
    ele_SIP[i] = 0 ; 

    ele_dzPV[i] = 0;
    ele_d0[i] = 0;
    ele_d0err[i] = 0;

    //
    ele_sclRawE[i]=0;
    ele_sclE[i]=0;
    ele_sclEt[i]=0;
    ele_sclEta[i]=0;
    ele_sclPhi[i]=0;
    ele_sclNclus[i]=0;
    ele_sclphiwidth[i]=0;
    ele_scletawidth[i]=0;
    
    for(int j=0;j<20;j++) {
      ele_sclsubE[i][j] = 0; //counter][countersub];
      ele_sclsubEta[i][j] = 0;
      ele_sclsubPhi[i][j] = 0;
      ele_sclsubisseed[i][j] = 0;
    } // for loop on subclusters
    
    ele_ecalE[i]=0;
    ele_ecalErr[i]=0;
    ele_trackErr[i]=0;
    ele_combErr[i]=0;
    ele_PFcombErr[i]=0;
    
    ele_ecalRegressionEnergy[i]  = 0;
    ele_ecalRegressionError[i] = 0;
    ele_ecalTrackRegressionEnergy[i]  = 0;
    ele_ecalTrackRegressionError[i]  = 0;
    ele_ecalScale[i]  = 0;
    ele_ecalSmear[i]  = 0;
    ele_ecalRegressionScale[i]  = 0;
    ele_ecalRegressionSmear[i]  = 0;
    ele_ecalTrackRegressionScale[i]  = 0;
    ele_ecalTrackRegressionSmear[i]  = 0;
    
    //  ele_mvafbrem[i]=0;
    //     ele_mvadetain[i]=0;
    //     ele_mvadphiin[i]=0;
    //     ele_mvasieie[i]=0;
    //     ele_mvahoe[i]=0;
    //     ele_mvaeop[i]=0;
    //     ele_mvae1x5e5x5[i]=0;
    //     ele_mvaeleopout[i]=0;
    ele_kfchi2[i]=0;
    ele_kfhits[i]=-1;
    ele_gsfhits[i] = -1;
    //     ele_mvamishits[i]=0;
    //     ele_mvadist[i]=0;
    //     ele_mvadcot[i]=0;
    //     ele_mvaeta[i]=0;
    //     ele_mvapt[i]=0;
    //     ele_mvaecalseed[i]=0;
    
    //ele_sclphiwidth[i]= 0;
    //ele_scletawidth[i]= 0;	
    ele_psE[i] = 0.;

    ele_mvaphys14[i] = 0.;
    ele_mvaphys14fix[i] = 0.;

  } // for loop on electrons


  // Gen Informations
  for(int ii=0;ii<10;ii++) {
    _MC_gen_V_pdgid[ii]               = 0;
  }

  for(int ii=0;ii<30;ii++) {
    _MC_gen_leptons_pdgid[ii]         = 0;
    _MC_gen_leptons_status1_pdgid[ii] = 0;
    _MC_gen_leptons_status2_pdgid[ii] = 0;
    
    _MC_gen_leptons_status1_FromWZ[ii]         = 0;
    _MC_gen_leptons_status1_FromTaus[ii]       = 0;
    _MC_gen_leptons_status1_FromNonPrompt[ii]  = 0;
    _MC_gen_leptons_status1_FromBC[ii]         = 0;
  } // for loop on gen leptons
  
  _MC_flavor[0] = 0;
  _MC_flavor[1] = 0;
    

}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Ntuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Ntuplizer);
