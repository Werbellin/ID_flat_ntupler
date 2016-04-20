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
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

//
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
//
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
//"
// MET
#include "DataFormats/METReco/interface/MET.h"
//
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

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

enum ElectronMatchType {UNMATCHED = 0, 
              TRUE_PROMPT_ELECTRON, 
              TRUE_ELECTRON_FROM_TAU,
              TRUE_NON_PROMPT_ELECTRON}; // The last does not include tau parents

//
// static data member definitions
//

//using namespace std;
using namespace reco;
using namespace edm;

void findFirstNonElectronMother(const reco::Candidate *particle,
                         int &ancestorPID, int &ancestorStatus){

  if( particle == 0 ){
    printf("SimpleElectronNtupler: ERROR! null candidate pointer, this should never happen\n");
    return;
  }

  // Is this the first non-electron parent? If yes, return, otherwise
  // go deeper into recursion
  if( abs(particle->pdgId()) == 11 ){
    findFirstNonElectronMother(particle->mother(0), ancestorPID, ancestorStatus);
  }else{
    ancestorPID = particle->pdgId();
    ancestorStatus = particle->status();
  }

  return;
}


const reco::Candidate* GetClosestGenParticle(const edm::Ptr<reco::GsfElectron> el, 
                  const edm::Handle<edm::View<reco::GenParticle>> &prunedGenParticles){

  // 
  // Explicit loop and geometric matching method (advised by Josh Bendavid)
  //

  // Find the closest status 1 gen electron to the reco electron
  double dR = 999;
  const reco::Candidate *closestParticle = nullptr;
  for(size_t i=0; i<prunedGenParticles->size();i++){
    const reco::Candidate *particle = &(*prunedGenParticles)[i];
    // Drop everything that is not electron or not status 1
    if( particle->status() != 1 )
      continue;
    //
    double dRtmp = ROOT::Math::VectorUtil::DeltaR( el->p4(), particle->p4() );
    if( dRtmp < dR ){
      dR = dRtmp;
      closestParticle = particle;
    }
  }
  if( dR > 0.1 ) {
    return nullptr;//UNMATCHED;
  } else {
    return closestParticle;
  }
}

float photon_E_over_electron_E(const edm::Ptr<reco::GsfElectron> el, 
                               const edm::Handle<edm::View<reco::GenParticle>> &prunedGenParticles) {
  // Find the closest status 1 gen electron to the reco electron

 float photon_pt = 0.;
 for(size_t i=0; i<prunedGenParticles->size();i++){
    const reco::Candidate *particle = &(*prunedGenParticles)[i];
    // Drop everything that is not electron or not status 1
    if( abs(particle->pdgId()) != 22 || particle->status() != 1 )
      continue;
    //
    //
    float dR_max = 0.1;
    double dRtmp = ROOT::Math::VectorUtil::DeltaR( el->p4(), particle->p4() );
    if( dRtmp < dR_max ){
      photon_pt += particle->pt();
    }
  }

  float photon_E_over_electron_E = 0.;
   if(el->pt() > 0.)
      photon_E_over_electron_E = photon_pt / el->pt();
   else photon_E_over_electron_E = -1;
  
  return photon_E_over_electron_E;
}


int matchToTruth(const edm::Ptr<reco::GsfElectron> el, 
                  const edm::Handle<edm::View<reco::GenParticle>> &prunedGenParticles){

  // 
  // Explicit loop and geometric matching method (advised by Josh Bendavid)
  //

  // Find the closest status 1 gen electron to the reco electron
  double dR = 999;
  const reco::Candidate *closestElectron = 0;
  for(size_t i=0; i<prunedGenParticles->size();i++){
    const reco::Candidate *particle = &(*prunedGenParticles)[i];
    // Drop everything that is not electron or not status 1
    if( abs(particle->pdgId()) != 11 || particle->status() != 1 )
      continue;
    //
    double dRtmp = ROOT::Math::VectorUtil::DeltaR( el->p4(), particle->p4() );
    if( dRtmp < dR ){
      dR = dRtmp;
      closestElectron = particle;
    }
  }
  // See if the closest electron (if it exists) is close enough.
  // If not, no match found.
  if( !(closestElectron != 0 && dR < 0.1) ) {
    return UNMATCHED;
  }

  // 
  int ancestorPID = -999; 
  int ancestorStatus = -999;
  findFirstNonElectronMother(closestElectron, ancestorPID, ancestorStatus);

  if( ancestorPID == -999 && ancestorStatus == -999 ){
    // No non-electron parent??? This should never happen.
    // Complain.
    printf("SimpleElectronNtupler: ERROR! Electron does not apper to have a non-electron parent\n");
    return UNMATCHED;
  }
  
  if( abs(ancestorPID) > 50 && ancestorStatus == 2 )
    return TRUE_NON_PROMPT_ELECTRON;

  if( abs(ancestorPID) == 15 && ancestorStatus == 2 )
    return TRUE_ELECTRON_FROM_TAU;

  // What remains is true prompt electrons
  return TRUE_PROMPT_ELECTRON;
}

string electronIDBaseName(string fullname, string producer_prefix, string producer_suffix) {

        LogError("") << fullname;
      string::size_type i = fullname.find(producer_prefix);  
      if(i != std::string::npos)
          fullname.erase(i, producer_prefix.length());
        LogError("") << fullname;

      i = fullname.find(producer_suffix);  
      if(i != std::string::npos)
          fullname.erase(i, producer_suffix.length());
         LogError("") << fullname;
   
      return fullname;
}

// =============================================================================================
// constructor
Ntuplizer::Ntuplizer(const edm::ParameterSet& iConfig) :
// ==============================================================================================

conf(iConfig),
inFileType(inputFileTypes::UNDEF),
isMC_ (iConfig.getParameter<bool>("isMC")),
ID1_use_userFloat_ (iConfig.getParameter<bool>("ID1_use_userFloat"))
{

  std::string inFileType_s = iConfig.getUntrackedParameter<std::string>("inputFileFormat");

  if(inFileType_s == "AOD")     inFileType = inputFileTypes::AOD; 
  else {
      LogError("") << "MINIAOD!";
      inFileType = inputFileTypes::MINIAOD; 
  }
      inFileType = inputFileTypes::MINIAOD; 

  if(inFileType == inputFileTypes::UNDEF) LogError("") << "Did not recognize input file format!";

  if(inFileType != inputFileTypes::MINIAOD && ID1_use_userFloat_) LogError("") << "Trying to get ID from userFloat but file input is not miniAOD!";

  beamSpotToken_    = consumes<reco::BeamSpot> 
                        (iConfig.getParameter <edm::InputTag>
                        ("beamSpot"));
  HLTToken = consumes<edm::TriggerResults>(iConfig.getParameter <edm::InputTag>
                        ("HLTTag")); 

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
    genParticlesToken_CB = mayConsume<edm::View<reco::GenParticle> > 
                        (iConfig.getParameter<edm::InputTag>
                        ("genParticlesAOD"));
    genEventInfoProductTagToken_ = consumes<GenEventInfoProduct>
                        (iConfig.getParameter<edm::InputTag>
                        ("genEventInfoProductAOD"));
  }

  if(inFileType == inputFileTypes::MINIAOD) {
    LogInfo("") << "Running on miniAOD";
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
    genParticlesToken_CB = mayConsume<edm::View<reco::GenParticle> > 
                        (iConfig.getParameter<edm::InputTag>
                        ("genParticlesMiniAOD"));
    genEventInfoProductTagToken_ = consumes<GenEventInfoProduct>
                        (iConfig.getParameter<edm::InputTag>
                        ("genEventInfoProductMiniAOD"));
 
  }

  electronEcalPFClusterIsolationProducerToken_ = mayConsume<ValueMap<float>>
                        (InputTag
                        ("electronEcalPFClusterIsolationProducer"));
  
  electronHcalPFClusterIsolationProducerToken_ = mayConsume<ValueMap<float>>
                        (InputTag
                        ("electronHcalPFClusterIsolationProducer"));



  electronID1_name = (iConfig.getParameter<string>("electronID1"));
  //electronID1_name = electronIDBaseName(electronID1_name, "electronMVAValueMapProducer:", "Values");

  electronID2_name = (iConfig.getParameter<string>("electronID2"));
  //electronID2_name = electronIDBaseName(electronID2_name, "electronMVAValueMapProducer:", "Values");



  if(!ID1_use_userFloat_) {

  electronID1Token_ = mayConsume<ValueMap<float>>
                        (InputTag
                        ("electronMVAValueMapProducer:" + electronID1_name + "Values"));

  electronID1CatToken_ = mayConsume<ValueMap<int>>
                        (InputTag
                        ("electronMVAValueMapProducer:" + electronID1_name + "Categories"));

  electronID1_pass_Token_ = mayConsume<ValueMap<bool>>
                        (iConfig.getParameter<edm::InputTag>
                        ("electronID1_pass"));
  }

  electronID2Token_ = mayConsume<ValueMap<float>>
                        (InputTag
                        ("electronMVAValueMapProducer:" + electronID2_name + "Values"));

  electronID2CatToken_ = mayConsume<ValueMap<int>>
                        (InputTag
                        ("electronMVAValueMapProducer:" + electronID2_name + "Categories"));

  electronID2_pass_Token_ = mayConsume<ValueMap<bool>>
                        (iConfig.getParameter<edm::InputTag>
                        ("electronID2_pass"));

  PUinfoToken = mayConsume<std::vector<PileupSummaryInfo>>
                        (InputTag
                        ("slimmedAddPileupInfo"));
}

// =============================================================================================
// destructor
Ntuplizer::~Ntuplizer()
// =============================================================================================
{
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  delete m_electrons ;
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
  _mytree->Branch("event_trig_fired", &event_trig_fired);
  _mytree->Branch("ele_trig_passed_filter", &ele_trig_passed_filter);
  _mytree->Branch("ele_pass_hltEle27WP75GsfTrackIsoFilter", &ele_pass_hltEle27WP75GsfTrackIsoFilter);



  // Vertices
  _mytree->Branch("vtx_N",&_vtx_N,"vtx_N/I");

  // Electrons
  _mytree->Branch("ele_N",&ele_N,"ele_N/I");
  _mytree->Branch("ele_N_saved",&ele_N_saved,"ele_N_saved/I");

  m_electrons = new TClonesArray ("TLorentzVector");
  _mytree->Branch ("electrons", "TClonesArray", &m_electrons, 256000,0);

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
  
  _mytree->Branch("ele_kfchi2",&ele_kfchi2,"ele_kfchi2[50]/D");
  _mytree->Branch("ele_kfhits",&ele_kfhits,"ele_kfhits[50]/I");
  _mytree->Branch("ele_gsfhits",&ele_gsfhits,"ele_gsfhits[50]/I");
   
  
  // Truth Leptons

  _mytree->Branch("mc_event_weight",&_mc_event_weight,"mc_event_weight/D");
  //
  _mytree->Branch ("MC_TrueNumInteractions",&_MC_TrueNumInteractions,"MC_TrueNumInteractions/I");

  _mytree->Branch("mc_ele_isPromptFinalState", &mc_ele_isPromptFinalState);
  _mytree->Branch("mc_ele_isDirectPromptTauDecayProductFinalState", &mc_ele_isDirectPromptTauDecayProductFinalState);
  _mytree->Branch("mc_ele_matchedFromCB", &mc_ele_matchedFromCB);
  _mytree->Branch("mc_ele_closestGenParticlePDGID", &mc_ele_matchMother_PDGID);
  _mytree->Branch("mc_ele_photon_over_ele_pt", &mc_ele_photon_over_ele_pt);

  _mytree->Branch("ele_dr03EcalRecHitSumEt", &ele_dr03EcalRecHitSumEt);
  _mytree->Branch("ele_dr03HcalTowerSumEt", &ele_dr03HcalTowerSumEt);
  _mytree->Branch("ele_dr03TkSumPt", &ele_dr03TkSumPt);
  _mytree->Branch("ele_pt", &ele_pt);


  _mytree->Branch("mc_gen_ele_p4", &_mc_gen_ele_p4);
  //_mytree->Branch("", &);
  _mytree->Branch("ele_electronEcalPFClusterIsolationProducer", &ele_electronEcalPFClusterIsolationProducer);
  _mytree->Branch("ele_electronHcalPFClusterIsolationProducer", &ele_electronHcalPFClusterIsolationProducer);
  _mytree->Branch("ele_full5x5_hcalOverEcal", &ele_full5x5_hcalOverEcal);

  _mytree->Branch("ele_ID1", &ele_ID1);
  _mytree->Branch("ele_ID2", &ele_ID2);

  _mytree->Branch("ele_ID1_pass", &ele_ID1_pass);
  _mytree->Branch("ele_ID2_pass", &ele_ID2_pass);

  _mytree->Branch("ele_ID1_cat", &ele_ID1_cat);
  _mytree->Branch("ele_ID2_cat", &ele_ID2_cat);


 _mytree->Branch("ele_index", &ele_index);

}


//
// member functions
//
// =============================================================================================
// ------------ method called for each event  ------------
void Ntuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
// =============================================================================================
{

   Init();

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
     FillTruth(iEvent, iSetup);
     LogDebug("") << "After FillTruth";
   }

   _mytree->Fill();

}

// =============================================================================================
void Ntuplizer::FillEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup)
//=============================================================================================
{
  _selectedObjects.clear();
  event_trig_fired.clear();

  _nEvent = iEvent.id().event();
  _nRun   = iEvent.id().run();
  _nLumi  = iEvent.luminosityBlock();

  // ----------------
  // Fired Triggers
  // ----------------
  Handle<edm::TriggerResults> triggerResultsHandle;
  iEvent.getByToken (HLTToken, triggerResultsHandle);
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
	  event_trig_fired.push_back(hlt_string);
      strcat(trig_fired_names_local,hlt_string.c_str());
	  strcat(trig_fired_names_local,"*");
	}
      }
    } // if HLT
  }
  strcpy(trig_fired_names,trig_fired_names_local);
  
if(inFileType == inputFileTypes::AOD) {
    //open the trigger summary
    edm::InputTag triggerSummaryLabel_ = edm::InputTag("hltTriggerSummaryAOD", "", "HLT");
    edm::Handle<trigger::TriggerEvent> triggerSummary;
    iEvent.getByLabel(triggerSummaryLabel_, triggerSummary);
    //trigger object we want to match
    edm::InputTag filterTag = edm::InputTag("hltEle12CaloIdLTrackIdLIsoVLTrackIsoFilter", "", "HLT"); //find the index corresponding to the event
    size_t filterIndex = (*triggerSummary).filterIndex(filterTag);
    trigger::TriggerObjectCollection allTriggerObjects = triggerSummary->getObjects(); //trigger::TriggerObjectCollection selectedObjects;
    if (filterIndex < (*triggerSummary).sizeFilters()) { //check if the trigger object is present ! 
    const trigger::Keys &keys = (*triggerSummary).filterKeys(filterIndex);
        for(size_t j = 0; j < keys.size(); j++) {
            trigger::TriggerObject foundObject = (allTriggerObjects)[keys[j]];
            _selectedObjects.push_back(foundObject);
        }
    }
    {
    edm::InputTag filterTag = edm::InputTag("hltEle27WP75GsfTrackIsoFilter", "", "HLT"); //find the index corresponding to the event
    size_t filterIndex = (*triggerSummary).filterIndex(filterTag);
    trigger::TriggerObjectCollection allTriggerObjects = triggerSummary->getObjects(); //trigger::TriggerObjectCollection hltEle27WP75GsfTrackIsoFilter;
    if (filterIndex < (*triggerSummary).sizeFilters()) { //check if the trigger object is present ! 
    const trigger::Keys &keys = (*triggerSummary).filterKeys(filterIndex);
        for(size_t j = 0; j < keys.size(); j++) {
            trigger::TriggerObject foundObject = (allTriggerObjects)[keys[j]];
            _hltEle27WP75GsfTrackIsoFilter.push_back(foundObject);
        }
    }
    }
}


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

  edm::Handle<edm::View<reco::GenParticle> > genCandidatesCollection;
  iEvent.getByToken(genParticlesToken_CB, genCandidatesCollection);

  
  // get the beam spot
  const reco::BeamSpot &beamSpot = *(beamspot_h.product());

   //electron PF isolation
  edm::Handle< edm::ValueMap<float> > electronECALIsoMapH;

  iEvent.getByToken(electronEcalPFClusterIsolationProducerToken_, electronECALIsoMapH);

  edm::Handle<edm::ValueMap<float> > ID1_map;
  edm::Handle<edm::ValueMap<int> > ID1_cat_map;

  if(!ID1_use_userFloat_) {
      iEvent.getByToken(electronID1Token_, ID1_map); 
      iEvent.getByToken(electronID1CatToken_, ID1_cat_map); 

  //iEvent.getByToken(electronID1_pass_Token_, ID1_pass_map); 
  }

  edm::Handle<edm::ValueMap<float> > ID2_map;
  iEvent.getByToken(electronID2Token_, ID2_map); 

  edm::Handle<edm::ValueMap<int> > ID2_cat_map;
  iEvent.getByToken(electronID2CatToken_, ID2_cat_map); 

  edm::Handle<edm::ValueMap<bool> > ID1_pass_map;
  edm::Handle<edm::ValueMap<bool> > ID2_pass_map;
  //iEvent.getByToken(electronID2_pass_Token_, ID2_pass_map); 


  //iEvent.getByLabel(ecalPFclusterIsolation_,electronECALIsoMapH);
  
  edm::Handle< edm::ValueMap<float> > electronHCALIsoMapH;
  

  if(inFileType == inputFileTypes::AOD) { 
    iEvent.getByToken(electronEcalPFClusterIsolationProducerToken_, electronECALIsoMapH);
    iEvent.getByToken(electronHcalPFClusterIsolationProducerToken_, electronHCALIsoMapH);
  }
 
  TClonesArray & electrons = *m_electrons;
  int counter = 0;
  ele_N = electronsColl_h->size();
  ele_N_saved = 0;

  ele_conversionVertexFitProbability.clear();
  ele_dr03EcalRecHitSumEt.clear();
  ele_dr03HcalTowerSumEt.clear();
  ele_dr03TkSumPt.clear();
  ele_pt.clear();
  ele_trackMomentumAtVtx_R.clear();

  ele_electronEcalPFClusterIsolationProducer.clear();
  ele_electronHcalPFClusterIsolationProducer.clear();

  ele_trig_passed_filter.clear();
  ele_pass_hltEle27WP75GsfTrackIsoFilter.clear();
  ele_full5x5_hcalOverEcal.clear();


  ele_ID1.clear();
  ele_ID2.clear();
  ele_index.clear();
  ele_ID1_pass.clear();
  ele_ID2_pass.clear();
  ele_ID1_cat.clear();
  ele_ID2_cat.clear();
 

  for(size_t i_ele = 0;  i_ele <  electronsColl_h->size(); ++i_ele) {
    if(counter>49) { continue; } 

    LogDebug("") << "Processing new electron";

    const auto ielectrons =  electronsColl_h->ptrAt(i_ele); 

    //const pat::electron pat_ele = nullptr;
    const edm::Ptr<pat::Electron> elePatPtr(ielectrons);

    if(inFileType == inputFileTypes::MINIAOD && elePatPtr.get() == NULL) {
        LogError("") << "Failed to get pointer to pat electron!";
    }
    if(ielectrons->pt() < 4.9) continue;

    ++ele_N_saved;


    float ID1_value = 999.;
    int ID1_cat = 999.; 

    float ID2_value = (*ID2_map)[ielectrons];
    int ID2_cat = (*ID2_cat_map)[ielectrons];


    if(ID1_use_userFloat_) {
        ID1_value = elePatPtr->userFloat(electronID1_name + "Values"); 
        ID1_cat = elePatPtr->userInt(electronID1_name + "Categories"); 
    } else {
        ID1_value = (*ID1_map)[ielectrons];
        ID1_cat = (*ID1_cat_map)[ielectrons];
    }

    ele_ID1.push_back(ID1_value);
    ele_ID2.push_back(ID2_value);

    ele_ID1_cat.push_back(ID1_cat);
    ele_ID2_cat.push_back(ID2_cat);

    //ele_ID1_pass.push_back((*ID1_pass_map)[ielectrons]);
    //ele_ID2_pass.push_back((*ID2_pass_map)[ielectrons]);


    bool _ele_trig_passed_filter = false;
    for(size_t t = 0; t < _selectedObjects.size(); ++t) {
        float deltaR = sqrt(pow(_selectedObjects[t].eta() - ielectrons->eta(), 2) + pow(acos(cos(_selectedObjects[t].phi() - ielectrons->phi())), 2));
        if(deltaR < 0.1) {//matching successfull
            _ele_trig_passed_filter = true;
            continue;
        }
    }
    ele_trig_passed_filter.push_back(_ele_trig_passed_filter);

    bool _ele_pass_hltEle27WP75GsfTrackIsoFilter = false;
    for(size_t t = 0; t < _hltEle27WP75GsfTrackIsoFilter.size(); ++t) {
        float deltaR = sqrt(pow(_hltEle27WP75GsfTrackIsoFilter[t].eta() - ielectrons->eta(), 2) + pow(acos(cos(_hltEle27WP75GsfTrackIsoFilter[t].phi() - ielectrons->phi())), 2));
        if(deltaR < 0.1) {//matching successfull
            _ele_pass_hltEle27WP75GsfTrackIsoFilter = true;
            continue;
        }
    }
    ele_pass_hltEle27WP75GsfTrackIsoFilter.push_back(_ele_pass_hltEle27WP75GsfTrackIsoFilter);

    ele_index.push_back(i_ele);

    ele_dr03EcalRecHitSumEt.push_back(ielectrons->dr03EcalRecHitSumEt());
    ele_dr03HcalTowerSumEt.push_back(ielectrons->dr03HcalTowerSumEt());
    ele_dr03TkSumPt.push_back(ielectrons->dr03TkSumPt());
    ele_pt.push_back(ielectrons->pt());
    ele_trackMomentumAtVtx_R.push_back(ielectrons->trackMomentumAtVtx().R());
    //edm::Ref<reco::GsfElectronCollection> electronRef(electronsCol, i_ele);
    float ECALIso = 999.;
    float HCALIso = 999.;

    if(inFileType == inputFileTypes::AOD) {
        const edm::ValueMap<float> electronHCALIsoMap = *(electronHCALIsoMapH);
        const edm::ValueMap<float> electronECALIsoMap = *(electronECALIsoMapH);
 
        ECALIso = electronECALIsoMap[ielectrons];
        HCALIso = electronHCALIsoMap[ielectrons];
    } else {
        ECALIso = elePatPtr->ecalPFClusterIso();
        HCALIso = elePatPtr->hcalPFClusterIso();
    }

    ele_electronEcalPFClusterIsolationProducer.push_back(ECALIso);
    ele_electronHcalPFClusterIsolationProducer.push_back(HCALIso);
    ele_full5x5_hcalOverEcal.push_back(ielectrons->full5x5_hcalOverEcal());

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
   
    ele_oldsigmaetaeta[counter]   =  ielectrons->full5x5_sigmaEtaEta();    
    ele_oldsigmaietaieta[counter] =  ielectrons->full5x5_sigmaIetaIeta();  
    ele_oldsigmaiphiiphi[counter] =  ielectrons->full5x5_sigmaIphiIphi();  
    ele_oldr9[counter]              =  ielectrons->full5x5_r9();  
    ele_olde15[counter]           =  ielectrons->full5x5_e1x5();
    ele_olde25max[counter]   =  ielectrons->full5x5_e2x5Max();
    ele_olde55[counter]           =  ielectrons->full5x5_e5x5();
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
    ele_dxy[counter]  = ielectrons->gsfTrack()->dxy() ;
    ele_dz[counter]   = ielectrons->gsfTrack()->dz() ;
    ele_dsz[counter]  = ielectrons->gsfTrack()->dsz() ;
    ele_dzPV[counter] = ielectrons->gsfTrack()->dz(pv->position());
    //

    // Conversion Rejection
    ele_conv_dcot[counter]   = ielectrons->convDist();
    ele_conv_dist[counter]   = ielectrons->convDcot();
    ele_conv_radius[counter] = ielectrons->convRadius();
    
    
    ele_expected_inner_hits[counter] = ielectrons->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS);
    
   

    reco::ConversionRef conv_ref = ConversionTools::matchedConversion(*ielectrons, conversions_h, beamSpot.position()); 
    bool vtxFitConversion = ConversionTools::hasMatchedConversion(*ielectrons, conversions_h, beamSpot.position());
    ele_vtxconv[counter] = vtxFitConversion;

    double vertexFitProbability = -1.;;
    if(!conv_ref.isNull()) {
        const reco::Vertex &vtx = conv_ref.get()->conversionVertex();
        if (vtx.isValid()) {
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

    ++counter;
  } // for loop on gsfelectrons

  if(counter>49) { ele_N = 50; cout << "Number of electrons>49, electrons_N set to 50" << endl;}
  
} // end of FillElectrons

// ====================================================================================
void Ntuplizer::FillMET (const edm::Event& iEvent, const edm::EventSetup& iSetup)
// ====================================================================================
{
	

  //edm::Handle< edm::View<reco::MET> > pfMEThandle;
  //iEvent.getByToken(pfMETToken_, pfMEThandle);
  
} // end of Fill MET



// ====================================================================================
void Ntuplizer::FillTruth(const edm::Event& iEvent, const edm::EventSetup& iSetup)
// ====================================================================================
{
  LogDebug("") << "Ntuplizer::FillTruth"; 

  Handle<std::vector< PileupSummaryInfo > >  PupInfo;
  iEvent.getByToken(PUinfoToken, PupInfo);
  
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


  Handle<GenEventInfoProduct> genEvtInfo;
  iEvent.getByToken(genEventInfoProductTagToken_, genEvtInfo );


  _mc_event_weight = genEvtInfo->weight();
  LogDebug("MC") << "MC generator event weight: " << _mc_event_weight;
 
  edm::Handle<vector<reco::GenParticle> > genCandidatesCollection;
  iEvent.getByToken(genParticleToken_, genCandidatesCollection); //genParticlesPruned
 
  mc_ele_isPromptFinalState.clear();
  mc_ele_isDirectPromptTauDecayProductFinalState.clear();
  _mc_gen_ele_p4.clear();

  mc_ele_matchedFromCB.clear();
  mc_ele_matchMother_PDGID.clear();
  mc_ele_photon_over_ele_pt.clear();

  Handle<View<GsfElectron>> electronsColl_h;
  iEvent.getByToken(electronsToken_, electronsColl_h);

  edm::Handle<edm::View<reco::GenParticle> > genCandidatesCollection_CB;
  iEvent.getByToken(genParticlesToken_CB, genCandidatesCollection_CB);


    for(size_t i_ele = 0;  i_ele <  electronsColl_h->size(); ++i_ele) {
        const auto ielectrons =  electronsColl_h->ptrAt(i_ele); 
        if(ielectrons->pt() < 4.9) continue;

        //ElectronMatchType
        int  matchType = matchToTruth(ielectrons, genCandidatesCollection_CB); 

        mc_ele_matchedFromCB.push_back(matchType);
        mc_ele_photon_over_ele_pt.push_back(photon_E_over_electron_E(ielectrons, genCandidatesCollection_CB));

        const reco::Candidate* genParticle = nullptr; 
        if(matchType == UNMATCHED) {
            genParticle = GetClosestGenParticle(ielectrons, genCandidatesCollection_CB);
        }
        int GEN_PDGID = 0;
        if(genParticle != nullptr) {
            GEN_PDGID = genParticle->pdgId();    
        }
        mc_ele_matchMother_PDGID.push_back(GEN_PDGID); 
    }

    // To be re-implemented
    //mc_ele_isPromptFinalState.push_back(p->isPromptFinalState());
    //mc_ele_isDirectPromptTauDecayProductFinalState.push_back(p->isDirectPromptTauDecayProductFinalState());	
	


    // ----------------------------
    //      Loop on particles
    // ----------------------------
    //for( auto p = genCandidatesCollection->begin();p != genCandidatesCollection->end(); ++ p ) {
    
    //} // for loop on particles
	
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

  } // for loop on electrons


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
