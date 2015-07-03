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
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

//
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
//
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "FWCore/Utilities/interface/isFinite.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
//"
// MET
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/METFwd.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"
//
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
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

conf(iConfig),
inFileType(inputFileTypes::UNDEF),
HLTTag_(iConfig.getParameter<edm::InputTag> ("HLTTag")),
isMC_ (iConfig.getParameter<bool>("isMC")),
ispythia6_ (iConfig.getParameter<bool>("ispythia6")),
PileupSrc_ ("addPileupInfo"),
MVAidCollection_ (iConfig.getParameter<std::vector<edm::InputTag> >("MVAId"))
//runGsfRefitter (iConfig.getParameter<bool>("runGsfRefitter")),
//GSFTrajColl (iConfig.getParameter<std::string>("GSFTrajectoryInput")),
//runKfWithGsfRefitter (iConfig.getParameter<bool>("runKfWithGsfRefitter")),
//CKFTrajColl (iConfig.getParameter<std::string>("CKFTrajectoryInput"))


// std::vector<edm::InputTag> MVAidCollection_;
{

  std::string inFileType_s = iConfig.getUntrackedParameter<std::string>("inputFileFormat");

  if(inFileType_s == "AOD")     inFileType = inputFileTypes::AOD; 
  if(inFileType_s == "MiniAOD") inFileType = inputFileTypes::MINIAOD; 

  if(inFileType == inputFileTypes::UNDEF) LogError("") << "Did not recognize input file format!";

  beamSpotToken_    = consumes<reco::BeamSpot> 
                        (iConfig.getParameter <edm::InputTag>
                        ("beamSpot"));

  

  if(inFileType == inputFileTypes::AOD) {
    electronsToken_ = mayConsume<edm::View<reco::GsfElectron> >
                        (iConfig.getParameter<edm::InputTag>
                        ("electronsAOD"));
    conversionsToken_ = mayConsume< reco::ConversionCollection >
                        (iConfig.getParameter<edm::InputTag>
                        ("conversionsAOD"));
    vertexToken_ = mayConsume<reco::VertexCollection>
                        (iConfig.getParameter<edm::InputTag>
                        ("verticesAOD"));

    pfMETToken_ = mayConsume<edm::View<reco::MET> >
                        (iConfig.getParameter<edm::InputTag>
                        ("PFMETAOD"));
    genParticleToken_ = mayConsume<vector<reco::GenParticle> >
                        (iConfig.getParameter<edm::InputTag>
                        ("genParticlesAOD"));

  }

  if(inFileType == inputFileTypes::MINIAOD) {
    electronsToken_ = mayConsume<edm::View<reco::GsfElectron> >
                        (iConfig.getParameter<edm::InputTag>
                        ("electronsMiniAOD"));
    conversionsToken_ = mayConsume< reco::ConversionCollection >
                        (iConfig.getParameter<edm::InputTag>
                        ("conversionsMiniAOD"));
    vertexToken_ = mayConsume<reco::VertexCollection>
                        (iConfig.getParameter<edm::InputTag>
                        ("verticesMiniAOD"));

    pfMETToken_ = mayConsume<edm::View<reco::MET> >
                        (iConfig.getParameter<edm::InputTag>
                        ("PFMETMiniAOD"));
    genParticleToken_ = mayConsume<vector<reco::GenParticle> >
                        (iConfig.getParameter<edm::InputTag>
                        ("genParticlesMiniAOD"));

  }



}

// =============================================================================================
// destructor
Ntuplizer::~Ntuplizer()
// =============================================================================================
{
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  delete m_electrons ;
 
  if(isMC_ ) {
    delete _m_MC_gen_V;
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


/*
  _mytree->Branch("ele_trackHitPDGID",&trackHitPDGID);
  _mytree->Branch("ele_lastHitPt",&lastHitPt);

  _mytree->Branch("ele_foundGSFTraj", &ele_foundGSFTraj);
  _mytree->Branch("ele_foundCKFTraj", &ele_foundCKFTraj);

  _mytree->Branch("ele_signedEstimateSumPred", &ele_signedEstimateSumPred);//;&ele_signedEstimateSumPred,"ele_signedEstimateSumPred[50]/D");
  //_mytree->Branch("ele_signedEstimateSumPred_A", &ele_signedEstimateSumPred_A,"ele_signedEstimateSumPred_A[50]/F");
  _mytree->Branch("ele_propagatorSignedEstimateSumPred", &ele_propagatorSignedEstimateSumPred);
  _mytree->Branch("ele_signSumPredNormVH", &ele_signSumPredNormVH);

  _mytree->Branch("ele_signedEstimateSumPredCKF", &ele_signedEstimateSumPredCKF);
  _mytree->Branch("ele_reducedChi2CKF", &ele_reducedChi2CKF);
*/
  _mytree->Branch("ele_conversionVertexFitProbability", &ele_conversionVertexFitProbability);

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
  _mytree->Branch("ele_r9",&ele_r9,"ele_r9[50]/D");
  //
  _mytree->Branch("ele_oldsigmaietaieta",&ele_oldsigmaietaieta,"ele_oldsigmaietaieta[50]/D");
  _mytree->Branch("ele_oldsigmaetaeta",&ele_oldsigmaetaeta,"ele_oldsigmaetaeta[50]/D");
  _mytree->Branch("ele_oldsigmaiphiiphi",&ele_oldsigmaiphiiphi,"ele_oldsigmaiphiiphi[50]/D");
  _mytree->Branch("ele_oldsigmaietaiphi",&ele_oldsigmaietaiphi,"ele_oldsigmaietaiphi[50]/D");
  _mytree->Branch("ele_olde15",&ele_olde15,"ele_olde15[50]/D");
  _mytree->Branch("ele_olde25max",&ele_olde25max,"ele_olde25max[50]/D");
  _mytree->Branch("ele_olde55",&ele_olde55,"ele_olde55[50]/D");
  _mytree->Branch("ele_oldr9",&ele_oldr9,"ele_oldr9[50]/D");
  //
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
    
  _mytree->Branch("ele_kfchi2",&ele_kfchi2,"ele_kfchi2[50]/D");
  _mytree->Branch("ele_kfhits",&ele_kfhits,"ele_kfhits[50]/I");
  _mytree->Branch("ele_gsfhits",&ele_gsfhits,"ele_gsfhits[50]/I");
   

  
  // PFMET
  _mytree->Branch("met_pf_et",&_met_pf_et,"met_pf_et/D");
  _mytree->Branch("met_pf_px",&_met_pf_px,"met_pf_px/D");
  _mytree->Branch("met_pf_py",&_met_pf_py,"met_pf_py/D");
  _mytree->Branch("met_pf_phi",&_met_pf_phi,"met_pf_phi/D");
  _mytree->Branch("met_pf_set",&_met_pf_set,"met_pf_set/D");
  _mytree->Branch("met_pf_sig",&_met_pf_sig,"met_pf_sig/D");

  // Truth Leptons
  _m_MC_gen_V = new TClonesArray ("TLorentzVector");
  _mytree->Branch ("MC_gen_V", "TClonesArray", &_m_MC_gen_V, 256000,0);
  _mytree->Branch ("MC_gen_V_pdgid",&_MC_gen_V_pdgid, "MC_gen_V_pdgid[10]/D");
  //
  _m_MC_gen_Higgs = new TClonesArray ("TLorentzVector");
  _mytree->Branch ("MC_gen_Higgs", "TClonesArray", &_m_MC_gen_Higgs, 256000,0);
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
  _mytree->Branch ("MC_flavor",&_MC_flavor,"MC_flavor[2]/I");
  
  _mytree->Branch ("MC_TrueNumInteractions",&_MC_TrueNumInteractions,"MC_TrueNumInteractions/I");

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
   LogDebug("") << "After FillEvent"; 
   FillVertices(iEvent, iSetup);
   LogDebug("") << "After FillVertices";
   m_electrons -> Clear();
   FillElectrons(iEvent, iSetup);
   LogDebug("") << "After FillElectrons";
   
   FillMET (iEvent, iSetup);
   LogDebug("") << "After FillMET";


   if(isMC_ ) {
     _MC_TrueNumInteractions = 0;
     _m_MC_gen_V->Clear();
     _m_MC_gen_Higgs->Clear();
     _m_MC_gen_leptons->Clear();
     _m_MC_gen_leptons_status1->Clear();
     _m_MC_gen_leptons_status2->Clear();
     FillTruth(iEvent, iSetup);
     LogDebug("") << "After FillTruth";
   }

   _mytree->Fill();

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
  

} // end of FillEvent


// =============================================================================================
void Ntuplizer::FillVertices(const edm::Event& iEvent, const edm::EventSetup& iSetup)
//=============================================================================================
{
  Handle<vector<reco::Vertex> >  recoPrimaryVertexCollection;
  iEvent.getByToken(vertexToken_, recoPrimaryVertexCollection);

  _vtx_N = recoPrimaryVertexCollection->size();
  
} // end of FillVertices

// =============================================================================================
void Ntuplizer::FillElectrons(const edm::Event& iEvent, const edm::EventSetup& iSetup)
//=============================================================================================
{
  LogVerbatim("") << "Running FillElectrons";
  LogDebug("") << "Ntuplizer::FillElectrons";

  edm::ESHandle<TransientTrackBuilder> builder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", builder);
  TransientTrackBuilder thebuilder = *(builder.product());

  //load the conversion collection

      edm::Handle<edm::View<reco::GsfElectron>> electronsColl_h;
      edm::Handle<reco::VertexCollection> primaryVertexColl_h;
      edm::Handle<reco::ConversionCollection> conversions_h;
      edm::Handle<reco::BeamSpot> beamspot_h;


  iEvent.getByToken(electronsToken_, electronsColl_h);
  iEvent.getByToken(vertexToken_, primaryVertexColl_h);
  iEvent.getByToken(conversionsToken_, conversions_h);
  iEvent.getByToken(beamSpotToken_, beamspot_h);

  Vertex dummy;
  const Vertex *pv = &dummy;
  if (primaryVertexColl_h->size() != 0) {
    pv = &*primaryVertexColl_h->begin();
  } else { // create a dummy PV                                                                                                                                                                                             
    Vertex::Error e;
    e(0, 0) = 0.0015 * 0.0015;
    e(1, 1) = 0.0015 * 0.0015;
    e(2, 2) = 15. * 15.;
    Vertex::Point p(0, 0, 0);
    dummy = Vertex(p, e, 0, 0, 0);
  }




  
  // get the beam spot
  const reco::BeamSpot &beamSpot = *(beamspot_h.product());

  // get the value map for eiD
  edm::Handle<edm::ValueMap<float> >  mapMVAcollection_phys14;
  edm::Handle<edm::ValueMap<float> >  mapMVAcollection_phys14fix;
  if (MVAidCollection_.size()>0 ) iEvent.getByLabel(MVAidCollection_[0] , mapMVAcollection_phys14);
  if (MVAidCollection_.size()>1 ) iEvent.getByLabel(MVAidCollection_[1], mapMVAcollection_phys14fix);
  
  
  TClonesArray & electrons = *m_electrons;
  int counter = 0;
  ele_N = electronsColl_h->size();

  //cout << "ele N = " << ele_N << endl;
  /*
  ele_signedEstimateSumPred.clear();
  ele_propagatorSignedEstimateSumPred.clear();
  ele_signSumPredNormVH.clear();

  ele_signedEstimateSumPredCKF.clear();
  ele_reducedChi2CKF.clear();
  ele_foundGSFTraj.clear();
  ele_foundCKFTraj.clear();
 
  trackHitPDGID.clear();
  lastHitPt.clear();
  */
  ele_conversionVertexFitProbability.clear();


  for (auto ielectrons=electronsColl_h->begin(); ielectrons != electronsColl_h->end(); ++ielectrons) {
    if(counter>49) { continue; } 
        edm::LogVerbatim("") << "Processing new electron";
/*
        const reco::GsfTrack* gsfTrack = ielectrons->gsfTrack().get();
        edm::LogVerbatim("") << "Processinf track with Pt=" << gsfTrack->pt() << " and eta=" << gsfTrack->eta();
        RecSimHitMatcher* theRecSimHitMatcher = nullptr;
        theRecSimHitMatcher = new RecSimHitMatcher(iEvent, iSetup, conf, gsfTrack);
        theRecSimHitMatcher->InitializeSimTrackRecHitAssociations(true);
        edm::LogInfo("test") << "After InitializeSimTrackRecHitAssociations";
        if(!theRecSimHitMatcher->AssociationSucceeded()) {
            edm::LogWarning("Association failed") << "The Association RecHit<>SimTrack failed";
            theRecSimHitMatcher = nullptr;
        }

        std::vector<int> pdgId(gsfTrack->found(), 0);
        float _lastHitPt = -1.;
        if(theRecSimHitMatcher != nullptr) {
            _lastHitPt = theRecSimHitMatcher->GetSimTrackToLastRecHit()->momentum().Pt();            
            pdgId = theRecSimHitMatcher->GetVectorOfHitPDGID();
            edm::LogVerbatim("") <<"PDGID list";
            for(const auto& id :pdgId) edm::LogVerbatim("") <<"PDG ID = " << id; 

        }
    trackHitPDGID.push_back(pdgId);
    lastHitPt.push_back(_lastHitPt);

    // variables that rely on GsfRefit
    float tmp_ele_signedEstimateSumPred = -200.;
    float tmp_ele_propagatorSignedEstimateSumPred = -200.;
    float tmp_ele_signSumPredNormVH = -200.;
    bool foundGSFTraj = false;
    const Trajectory* theGSFTraj = nullptr;
    if(runGsfRefitter) {
        edm::LogInfo("Test") << "GSF refit was run";
        edm::Handle<std::vector<Trajectory> > GSFTrajCollectionHandle;
        iEvent.getByLabel(GSFTrajColl, GSFTrajCollectionHandle);
        theGSFTraj = GetTrajectoryFromTrack(GSFTrajCollectionHandle.product(), ielectrons->gsfTrack().get()); 
    }

    if(theGSFTraj == nullptr) {
        edm::LogWarning("No GSF Traj") << " Did not find the GSF Trajectory for track";
    } else {
        foundGSFTraj = true;
        edm::LogInfo("Ntup") << "Found Trajectory with hits " << theGSFTraj->foundHits();

        tmp_ele_signedEstimateSumPred = ReducedChi2AsymmetryNormalizedNValidHits(theGSFTraj, ielectrons->gsfTrack().get());
        tmp_ele_signSumPredNormVH = SignAsymmetryNormalizedNValidHits(theGSFTraj, gsfTrack);
        tmp_ele_propagatorSignedEstimateSumPred = GetSignedChiPropagator(theGSFTraj, *ielectrons, iEvent, iSetup);    
    }
    ele_foundGSFTraj.push_back(foundGSFTraj);
    ele_signedEstimateSumPred.push_back(tmp_ele_signedEstimateSumPred);
    //ele_signedEstimateSumPred_A[counter] = tmp_ele_signedEstimateSumPred;
    ele_propagatorSignedEstimateSumPred.push_back(tmp_ele_propagatorSignedEstimateSumPred);
    ele_signSumPredNormVH.push_back(tmp_ele_signSumPredNormVH);

    // variables that rely on KfTrack from GsfTrack fit

    float tmp_ele_signedEstimateSumPredCKF = -200.;
    float tmp_ele_reducedChi2CKF = -200.;
    bool foundCKFTraj = false;
   
    if(runKfWithGsfRefitter) {
        edm::Handle<std::vector<Trajectory> > KfWithGsfTrajCollectionHandle;
        iEvent.getByLabel(CKFTrajColl, KfWithGsfTrajCollectionHandle);
        edm::LogInfo("Demo") << "Size of trajectory collection " << KfWithGsfTrajCollectionHandle.product()->size(); 
        //h_NumberKfWithGsfTracks->Fill(KfWithGsfTrajCollectionHandle->size());
         
        const Trajectory* theCKFTraj  = GetTrajectoryFromTrack(KfWithGsfTrajCollectionHandle.product(), gsfTrack);
        if(theCKFTraj == nullptr) {
            edm::LogWarning("No CKF Traj") << " Did not find the CKF Trajectory for track";
        } else {
            foundCKFTraj = true;
            tmp_ele_reducedChi2CKF = theCKFTraj->chiSquared() / theCKFTraj->ndof();
            tmp_ele_signedEstimateSumPredCKF = ReducedChi2AsymmetryNormalizedNValidHits(theCKFTraj, gsfTrack);
        }
    }

    ele_foundCKFTraj.push_back(foundCKFTraj);
    ele_signedEstimateSumPredCKF.push_back(tmp_ele_signedEstimateSumPredCKF);
    ele_reducedChi2CKF.push_back(tmp_ele_reducedChi2CKF);
*/
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
   
    ele_oldsigmaetaeta[counter]   =  ielectrons->full5x5_sigmaEtaEta();    //( !edm::isNotFinite(Cov[0]) ) ? sqrt(Cov[0]) : 0;
    ele_oldsigmaietaieta[counter] =  ielectrons->full5x5_sigmaIetaIeta();   //( !edm::isNotFinite(vCov[0]) ) ? sqrt(vCov[0]) : 0;
    ele_oldsigmaiphiiphi[counter] =  ielectrons->full5x5_sigmaIphiIphi();   //( !edm::isNotFinite(vCov[2]) ) ? sqrt(vCov[2]) : 0;
    ele_oldr9[counter]              =  ielectrons->full5x5_r9();  //lazyToolsNoZS.e3x3(*seedCluster) / ielectrons->superCluster()->rawEnergy() ;
    ele_olde15[counter]           =  ielectrons->full5x5_e1x5(); //lazyToolsNoZS.e1x5(*seedCluster);
    ele_olde25max[counter]   =  ielectrons->full5x5_e2x5Max(); //lazyToolsNoZS.e2x5Max(*seedCluster);
    ele_olde55[counter]           =  ielectrons->full5x5_e5x5();       // lazyToolsNoZS.e5x5(*seedCluster);
    ele_oldhe[counter]             =  ielectrons->full5x5_hcalOverEcal();
    ele_oldhebc[counter]        =  ielectrons->full5x5_hcalOverEcalBc();


    // E/P combination
    ele_ecalE[counter]     = ielectrons->ecalEnergy();
    ele_ecalErr[counter]   = ielectrons->ecalEnergyError();
    ele_trackErr[counter]  = ielectrons->trackMomentumError();
    ele_combErr[counter]   = ielectrons->p4Error(GsfElectron::P4_COMBINATION);
    ele_PFcombErr[counter] = ielectrons->p4Error(GsfElectron::P4_PFLOW_COMBINATION);
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
    LogDebug("") << "Start looking at track vars";
    ele_lost_hits[counter]   = ielectrons->gsfTrack()->lost(); //numberOfLostHits();
    ele_valid_hits[counter] = ielectrons->gsfTrack()->found(); //numberOfValidHits() ;
    ele_gsfchi2[counter]    = ielectrons->gsfTrack()->normalizedChi2() ;
    LogDebug("") << "After gsfTrack";
   
    ele_gsfhits[counter] = ielectrons->gsfTrack()->hitPattern().trackerLayersWithMeasurement();

    bool validKF=false;
    reco::TrackRef myTrackRef = ielectrons->closestCtfTrackRef();
    LogDebug("") << "After clostste ctdf";

    validKF = myTrackRef.isAvailable();
    LogDebug("") << "isAvailable : " << validKF;

    validKF &= myTrackRef.isNonnull();
    LogDebug("") << "isNonnull&isAvailable : " << validKF;
  
    ele_kfchi2[counter] = validKF ? myTrackRef->normalizedChi2() : 0 ; //ielectrons->track()->normalizedChi2() : 0 ;
    ele_kfhits[counter] = validKF ? myTrackRef->hitPattern().trackerLayersWithMeasurement() : -1. ; //ielectrons->track()->hitPattern().trackerLayersWithMeasurement() : 0 ;

    LogDebug("") << "After ctfTrack";
    //
    // 		ele_dxyB[counter] = ielectrons->gsfTrack()->dxy(bs.position()) ;
    ele_dxy[counter]  = ielectrons->gsfTrack()->dxy() ;
    // 		ele_dzB[counter]  = ielectrons->gsfTrack()->dz(bs.position()) ;
    ele_dz[counter]   = ielectrons->gsfTrack()->dz() ;
    // 		ele_dszB[counter] = ielectrons->gsfTrack()->dsz(bs.position()) ;
    ele_dsz[counter]  = ielectrons->gsfTrack()->dsz() ;
    ele_dzPV[counter] = ielectrons->gsfTrack()->dz(pv->position());
    //

    // Conversion Rejection
    ele_conv_dcot[counter]   = ielectrons->convDist(); //userFloat("dcot");
    ele_conv_dist[counter]   = ielectrons->convDcot(); //ielectrons->userFloat("dist");
    ele_conv_radius[counter] = ielectrons->convRadius();
    
    
    ele_expected_inner_hits[counter] = ielectrons->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS);
    
   

    reco::ConversionRef conv_ref = ConversionTools::matchedConversion(*ielectrons, conversions_h, beamSpot.position()); 
    bool vtxFitConversion = ConversionTools::hasMatchedConversion(*ielectrons, conversions_h, beamSpot.position());
    ele_vtxconv[counter] = vtxFitConversion;

    double vertexFitProbability = -1.;;
    if(!conv_ref.isNull()) {
        const reco::Vertex &vtx = conv_ref.get()->conversionVertex();
   
        //vertex validity
        if (vtx.isValid()) {
   
            //fit probability
            vertexFitProbability = TMath::Prob( vtx.chi2(),  vtx.ndof());
            
        }
    }
    ele_conversionVertexFitProbability.push_back(vertexFitProbability);
    LogVerbatim("") << "Vertex fit probability: " << vertexFitProbability;
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
   
    double R  = TMath::Sqrt(sclRef->x()*sclRef->x() + sclRef->y()*sclRef->y() +sclRef->z()*sclRef->z());
    double Rt = TMath::Sqrt(sclRef->x()*sclRef->x() + sclRef->y()*sclRef->y());
    ele_sclRawE[counter]  = sclRef->rawEnergy() ;
    
    ele_sclE[counter]     = sclRef->energy() ;
    ele_sclEt[counter]    = sclRef->energy()*(Rt/R) ;
    ele_sclEta[counter]   = sclRef->eta() ;
    ele_sclPhi[counter]   = sclRef->phi() ;
    ele_sclNclus[counter] = sclRef->clustersSize();
    

    ele_sclphiwidth[counter] = sclRef->phiWidth();
    ele_scletawidth[counter] = sclRef->etaWidth();

    LogDebug("") << "After access to SC ref";
    // --------------
    // Sub-Clusters
    // --------------
    int countersub = 0;
   
    if(inFileType == inputFileTypes::AOD) { 
      reco::CaloCluster_iterator itscl  = sclRef->clustersBegin();
      reco::CaloCluster_iterator itsclE = sclRef->clustersEnd();
      LogDebug("") << "After reco::CaloCluster_iterator itsclE = sclRef->clustersEnd();"; 
      for(; itscl < itsclE ; ++itscl) {
        bool isseed = false;
        if((*itscl)==ielectrons->superCluster()->seed()) isseed=true; // continue; // skip seed cluster
        LogDebug("") << "After if((*itscl)==ielectrons->superCluster()->seed())";
        ele_sclsubE[counter][countersub]      = (*itscl)->energy();
        ele_sclsubEta[counter][countersub]    = (*itscl)->eta();
        ele_sclsubPhi[counter][countersub]    = (*itscl)->phi();
        ele_sclsubisseed[counter][countersub] = isseed;
        countersub++;
      }
    }
    LogDebug("") << "After SubClusters";

    // -----------------------------------------------------------------
    // Get PreShower Informations
    // -----------------------------------------------------------------
    ele_psE[counter] = sclRef->preshowerEnergy();

    // -----------------------------------------------------------------
    //fbrem
    // -----------------------------------------------------------------
    ele_fbrem[counter]      = ielectrons->fbrem();
    ele_trackfbrem[counter] = ielectrons->trackFbrem();
    ele_SCfbrem[counter]    = ielectrons->superClusterFbrem();
    // GsfElectron definition changed in 7X, 
    ele_pfSCfbrem[counter]  = 0.;//ielectrons->pfSuperClusterFbrem(); // Should be identical to the previous one...
    ele_eClass[counter]     = ielectrons->classification() ;
    ele_nbrem[counter]      = ielectrons->numberOfBrems();
    
    // -----------------------------------------------------------------
    // Electron ID electronsCol
    // -----------------------------------------------------------------
    

    // Does this do anything - doesnt look like it
    //edm::Ref<reco::GsfElectronCollection> electronRef(electronsColl_h, counter); //electronsCollection,i); i++; //reference to the electron


    ++counter;
  } // for loop on gsfelectrons

  if(counter>49) { ele_N = 50; cout << "Number of electrons>49, electrons_N set to 50" << endl;}
  
} // end of FillElectrons

// ====================================================================================
void Ntuplizer::FillMET (const edm::Event& iEvent, const edm::EventSetup& iSetup)
// ====================================================================================
{
	

  edm::Handle< edm::View<reco::MET> > pfMEThandle;
  iEvent.getByToken(pfMETToken_, pfMEThandle);
  
  
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
  LogDebug("") << "Ntuplizer::FillTruth"; 

  Handle<std::vector< PileupSummaryInfo > >  PupInfo;
  iEvent.getByLabel(edm::InputTag("addPileupInfo"), PupInfo);
  
  std::vector<PileupSummaryInfo>::const_iterator PVI;
  
  int Tnpv = -1;
  for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
  
     int BX = PVI->getBunchCrossing();
  
     if(BX == 0) { 
       Tnpv = PVI->getTrueNumInteractions();
       continue;
     }
  
  }
  _MC_TrueNumInteractions = Tnpv; 
 
  edm::Handle<vector<reco::GenParticle> > genCandidatesCollection;
  iEvent.getByToken(genParticleToken_, genCandidatesCollection); //genParticlesPruned
  
  TClonesArray &MC_gen_V               = *_m_MC_gen_V;
  TClonesArray &MC_gen_Higgs           = *_m_MC_gen_Higgs;
  TClonesArray &MC_gen_leptons         = *_m_MC_gen_leptons;
  TClonesArray &MC_gen_leptons_status2 = *_m_MC_gen_leptons_status2;
  TClonesArray &MC_gen_leptons_status1 = *_m_MC_gen_leptons_status1;
  
  int counter             = 0;
  int counter_higgs       = 0;
  int counter_daughters   = 0;
  int counter_lep_status2 = 0;
  int counter_lep_status1 = 0;
  
  // ----------------------------
  //      Loop on particles
  // ----------------------------
  for( auto p = genCandidatesCollection->begin();p != genCandidatesCollection->end(); ++ p ) {
    
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
      
     
        if(p->status()==1) {
	
	    setMomentum(myvector, p->p4());
	    new (MC_gen_leptons_status1[counter_lep_status1]) TLorentzVector(myvector);
	    _MC_gen_leptons_status1_pdgid[counter_lep_status1] = p->pdgId();
	
	
	    // Let's look at the mothers ... [works really only for electrons...=>eid studies]
	    if(p->numberOfMothers()>0) {
	        const reco::Candidate  *Mom = p->mother();
	        while(Mom!=0) {
	            if(fabs(Mom->pdgId()) == 23 || fabs(Mom->pdgId()) == 24) { 
	                _MC_gen_leptons_status1_FromWZ[counter_lep_status1] = 1;
	                break;
                }
	            else if(fabs(Mom->pdgId()) == 15) { 
	                _MC_gen_leptons_status1_FromTaus[counter_lep_status1] = 1;
	                break;
                }
	            else if(fabs(Mom->pdgId()) > 50 ) { 
	                _MC_gen_leptons_status1_FromNonPrompt[counter_lep_status1] = 1;

	                const reco::Candidate  *GrandMom = Mom; //->mother();
	                while(GrandMom!=0) {
	                    for(int unsigned jj=0;jj<GrandMom->numberOfMothers();jj++) {
	                            if(fabs(GrandMom->mother(jj)->pdgId()) == 4 || fabs(GrandMom->mother(jj)->pdgId()) == 5) { 
	                                _MC_gen_leptons_status1_FromBC[counter_lep_status1] = 1;
	                                break;
                                }
	                    } // for loop on grand-mothers, when they are more than 1...
	                    if(fabs(GrandMom->pdgId()) == 4 || fabs(GrandMom->pdgId()) == 5) { 
	                        _MC_gen_leptons_status1_FromBC[counter_lep_status1] = 1;
	                        break;
                        }
	                    GrandMom = GrandMom ->mother();
	                } // while loop on grand mom
	                break;
                } // if pdgid>50
	        
	            Mom = Mom ->mother();
	        } // while loop on mothers
	    } // if mother
	
	    counter_lep_status1++;
    } // if status 1

    if(p->status()==3) {
	    if(p->numberOfMothers()>0) { // Need a mother...
	        if(p->mother(0)->pdgId()==23) {  // If Mother is a Z 
	            if(p->numberOfDaughters()>0) { // Need a daughter...
	      
	      // Status 2 Leptons
	                if(p->daughter(0)->pdgId()==p->pdgId() && p->daughter(0)->status()==2) { // if daughter is lepton & status 2
		                setMomentum(myvector, p->daughter(0)->p4());
		                new (MC_gen_leptons_status2[counter_lep_status2]) TLorentzVector(myvector);
		                _MC_gen_leptons_status2_pdgid[counter_lep_status2] = p->daughter(0)->pdgId();
		                counter_lep_status2++;
		
	                } // if Daughter Status = 2
	            } // Need a daughter
	        } // if MOther = Z
	    } // if Nmother (status3) >0
    } // if hard scatter electron (status3)
} // if leptons
   
    // %%%%%%%%%%%%%%%%%%
    //     If W or Z
    // %%%%%%%%%%%%%%%%%%
    if (p->pdgId() == 23 || fabs(p->pdgId())==24) {
     
        bool pass_status=false;
        if(p->status()==3 && ispythia6_==true) pass_status=true; //{ // if PYTHIA6
        if(p->status() && p->numberOfDaughters()==2 && ispythia6_==false) pass_status=true; //{ // if PYTHIA 8
      
        if(pass_status) {
	        // Fill truth W,Z
	        setMomentum (myvector,p->p4());
	        new (MC_gen_V[counter]) TLorentzVector(myvector);
	        _MC_gen_V_pdgid[counter] = p->pdgId();
	        
	        // Loop on daughters
	        for(unsigned int i=0;i<p->numberOfDaughters();i++) {
	            bool islep = false;
	            if(fabs(p->daughter(i)->pdgId())==11) { _MC_flavor[counter] = 0; islep=true;} // electron
	            if(fabs(p->daughter(i)->pdgId())==13) { _MC_flavor[counter] = 1; islep=true;} // muon
	            if(fabs(p->daughter(i)->pdgId())==15) { _MC_flavor[counter] = 2; islep=true;} // taus
  
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
} // for loop on particles
	
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
    

    ele_kfchi2[i]=0;
    ele_kfhits[i]=-1;
    ele_gsfhits[i] = -1;
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
