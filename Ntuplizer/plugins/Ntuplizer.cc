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
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
//"
// MET
#include "DataFormats/METReco/interface/MET.h"
//
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

// Trigger
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"


#include <LeptonIsolationHelper/Helper/interface/LeptonIsoHelper.h>

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
void findFirstNonPhotonMother(const reco::Candidate *particle,
                         int &ancestorPID, int &ancestorStatus){

  if( particle == 0 ){
    printf("SimpleElectronNtupler: ERROR! null candidate pointer, this should never happen\n");
    return;
  }

  // Is this the first non-electron parent? If yes, return, otherwise
  // go deeper into recursion
  if( abs(particle->pdgId()) == 22 ){
    findFirstNonPhotonMother(particle->mother(0), ancestorPID, ancestorStatus);
  }else{
    ancestorPID = particle->pdgId();
    ancestorStatus = particle->status();
  }

  return;
}
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


const reco::Candidate* GetClosestGenParticle(const edm::Ptr<reco::Candidate> el, 
                  const edm::Handle<edm::View<reco::GenParticle>> &prunedGenParticles){

  // 
  // Explicit loop and geometric matching method (advised by Josh Bendavid)
  //

  // Find the closest status 1 gen electron to the reco electron
  float dR = 999;
  const reco::Candidate *closestParticle = nullptr;
  for(size_t i=0; i<prunedGenParticles->size();i++){
    const reco::Candidate *particle = &(*prunedGenParticles)[i];
    // Drop everything that is not electron or not status 1
    if( particle->status() != 1 )
      continue;
    //
    float dRtmp = ROOT::Math::VectorUtil::DeltaR( el->p4(), particle->p4() );
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

float photon_E_over_electron_E(const edm::Ptr<reco::Candidate> el, 
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
    float dRtmp = ROOT::Math::VectorUtil::DeltaR( el->p4(), particle->p4() );
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


int matchToTruth(const edm::Ptr<reco::Candidate> el, 
                  const edm::Handle<edm::View<reco::GenParticle>> &prunedGenParticles, float max_dR = 0.1){

  // 
  // Explicit loop and geometric matching method (advised by Josh Bendavid)
  //

  // Find the closest status 1 gen electron to the reco electron
  float dR = 999;
  const reco::Candidate *closestElectron = 0;
  for(size_t i=0; i<prunedGenParticles->size();i++){
    const reco::Candidate *particle = &(*prunedGenParticles)[i];
    // Drop everything that is not electron or not status 1
    if( abs(particle->pdgId()) != 11 || particle->status() != 1 )
      continue;
    //
    float dRtmp = ROOT::Math::VectorUtil::DeltaR( el->p4(), particle->p4() );
    if( dRtmp < dR ){
      dR = dRtmp;
      closestElectron = particle;
    }
  }
  // See if the closest electron (if it exists) is close enough.
  // If not, no match found.
  if( !(closestElectron != 0 && dR < max_dR) ) {
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

  if(inFileType_s == "AOD") {
     inFileType = inputFileTypes::AOD; 
     isAOD = true;
  } else {
      LogError("") << "MINIAOD!";
      inFileType = inputFileTypes::MINIAOD; 
    isAOD = false;
  }

  if(inFileType == inputFileTypes::UNDEF) LogError("") << "Did not recognize input file format!";

  if(inFileType != inputFileTypes::MINIAOD && ID1_use_userFloat_) LogError("") << "Trying to get ID from userFloat but file input is not miniAOD!";

  do_TLE = iConfig.getParameter<bool>("do_TLE");
  do_signal = iConfig.getParameter<bool>("do_signal");
  do_min_dR_cut = iConfig.getParameter<bool>("do_min_dR_cut");

  beamSpotToken_    = consumes<reco::BeamSpot> 
                        (iConfig.getParameter <edm::InputTag>
                        ("beamSpot"));
  HLTToken = consumes<edm::TriggerResults>(iConfig.getParameter <edm::InputTag>
                        ("HLTTag")); 

  sampleType = 2015;
  setup = sampleType;
  rhoForEleToken = consumes<double>(LeptonIsoHelper::getEleRhoTag(sampleType, setup));

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
    PUinfoToken = mayConsume<std::vector<PileupSummaryInfo>>
                        (InputTag
                        ("addPileupInfo"));

    ebReducedRecHitCollection_        = mayConsume<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>
                                       ("ebReducedRecHitCollectionAOD"));

    eeReducedRecHitCollection_        = mayConsume<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>
                                       ("eeReducedRecHitCollectionAOD"));

    esReducedRecHitCollection_        = mayConsume<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>
                                       ("esReducedRecHitCollectionAOD"));

    phoChargedIsolation_CITKToken_ = mayConsume<ValueMap<float>>
                        (InputTag
                        ("egmPhotonIsolationAOD:h+-DR030-"));

    phoNeutralHadronIsolation_CITKToken_ = mayConsume<ValueMap<float>>
                        (InputTag
                        ("egmPhotonIsolationAOD:h0-DR030-"));

    phoPhotonIsolation_CITKToken_ = mayConsume<ValueMap<float>>
                        (InputTag
                        ("egmPhotonIsolationAOD:gamma-DR030-"));


  }

  if(inFileType == inputFileTypes::MINIAOD) {
    LogInfo("") << "Running on miniAOD";
    electronsToken_ = mayConsume<edm::View<reco::GsfElectron> >
                        (iConfig.getParameter<edm::InputTag>
                        ("electronsMiniAOD"));
    photonsToken_ = mayConsume<edm::View<reco::Photon> >
                        (iConfig.getParameter<edm::InputTag>
                        ("photonsMiniAOD"));

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
    PUinfoToken = mayConsume<std::vector<PileupSummaryInfo>>
                        (InputTag
                        ("slimmedAddPileupInfo"));
    ebReducedRecHitCollection_        = mayConsume<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>
                                       ("ebReducedRecHitCollectionMiniAOD"));

    eeReducedRecHitCollection_        = mayConsume<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>
                                       ("eeReducedRecHitCollectionMiniAOD"));

    esReducedRecHitCollection_        = mayConsume<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>
                                       ("esReducedRecHitCollectionMiniAOD"));

    phoChargedIsolation_CITKToken_ = mayConsume<ValueMap<float>>
                        (InputTag
                        ("egmPhotonIsolationMiniAOD:h+-DR030-"));

    phoNeutralHadronIsolation_CITKToken_ = mayConsume<ValueMap<float>>
                        (InputTag
                        ("egmPhotonIsolationMiniAOD:h0-DR030-"));

    phoPhotonIsolation_CITKToken_ = mayConsume<ValueMap<float>>
                        (InputTag
                        ("egmPhotonIsolationMiniAOD:gamma-DR030-"));
  }


  electronID1_name = (iConfig.getParameter<string>("electronID1"));
  //electronID1_name = electronIDBaseName(electronID1_name, "electronMVAValueMapProducer:", "Values");

  electronID2_name = (iConfig.getParameter<string>("electronID2"));
  //electronID2_name = electronIDBaseName(electronID2_name, "electronMVAValueMapProducer:", "Values");

  if(iConfig.exists("photonID1") && iConfig.exists("photonID2")) {
    photonID1_name = (iConfig.getParameter<string>("photonID1"));
    photonID2_name = (iConfig.getParameter<string>("photonID2"));
  }

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
  photonID1Token_ = mayConsume<ValueMap<float>>
                        (InputTag
                        ("photonMVAValueMapProducer:" + photonID1_name + "Values"));
  photonID2Token_ = mayConsume<ValueMap<float>>
                        (InputTag
                        ("photonMVAValueMapProducer:" + photonID2_name + "Values"));

  outputFile   = iConfig.getParameter<std::string>("outputFile");
  outputPath   = iConfig.getParameter<std::string>("outputPath");

  string postfix;
  if(do_signal) postfix = "_signal";
  else postfix = "_background";
  file = TFile::Open((outputPath + outputFile + postfix + ".root").c_str(),"RECREATE");

 }

// =============================================================================================
// destructor
Ntuplizer::~Ntuplizer()
// =============================================================================================
{
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  delete m_electrons ;
  delete file;
}

// =============================================================================================
// ------------ method called once each job just before starting event loop  ------------
void Ntuplizer::beginJob()
//=============================================================================================
{
  // Book histograms
  //edm::Service<TFileService> fs ;
  //_mytree  = fs->make <TTree>("tree","tree"); 

  file->cd(); 
  _mytree = new TTree("tree", "tree");

  //// Counters
  
  // Global
  _mytree->Branch("nEvent",&_nEvent,"nEvent/I");
  _mytree->Branch("nRun",&_nRun,"nRun/I");
  _mytree->Branch("nLumi",&_nLumi,"nLumi/I");
  
  // Pile UP
  _mytree->Branch("PU_N",&_PU_N,"PU_N/I");
  _mytree->Branch("ele_rho",&ele_rho,"ele_rho/F");
 


  // Vertices
  _mytree->Branch("vtx_N",&_vtx_N,"vtx_N/I");

  // Electrons
  _mytree->Branch("ele_N",&ele_N,"ele_N/I");
//  _mytree->Branch("ele_N_saved",&ele_N_saved,"ele_N_saved/I");

  //m_electrons = new TClonesArray ("TLorentzVector");
  //_mytree->Branch ("electrons", "TClonesArray", &m_electrons, 256000,0);

 //
  //_mytree->Branch("ele_he",&ele_he,"ele_he/F");
  //_mytree->Branch("ele_hebc",&ele_hebc,"ele_hebc/F");
 _mytree->Branch("ele_hadronicOverEm",&ele_hadronicOverEm,"ele_hadronicOverEm/F");

  //
 /*
  _mytree->Branch("ele_sigmaietaieta",&ele_sigmaietaieta,"ele_sigmaietaieta/F");
  _mytree->Branch("ele_sigmaetaeta",&ele_sigmaetaeta,"ele_sigmaetaeta/F");
  _mytree->Branch("ele_sigmaiphiiphi",&ele_sigmaiphiiphi,"ele_sigmaiphiiphi/F");
  _mytree->Branch("ele_e15",&ele_e15,"ele_e15/F");
  _mytree->Branch("ele_e25max",&ele_e25max,"ele_e25max/F");
  _mytree->Branch("ele_e55",&ele_e55,"ele_e55/F");
  _mytree->Branch("ele_r9",&ele_r9,"ele_r9/F");
  */
  _mytree->Branch("ele_oldsigmaietaieta",&ele_oldsigmaietaieta,"ele_oldsigmaietaieta/F");
//  _mytree->Branch("ele_oldsigmaetaeta",&ele_oldsigmaetaeta,"ele_oldsigmaetaeta/F");
  _mytree->Branch("ele_oldsigmaiphiiphi",&ele_oldsigmaiphiiphi,"ele_oldsigmaiphiiphi/F");
  //_mytree->Branch("ele_oldsigmaietaiphi",&ele_oldsigmaietaiphi,"ele_oldsigmaietaiphi/F");
  _mytree->Branch("ele_olde15",&ele_olde15,"ele_olde15/F");
  _mytree->Branch("ele_olde25max",&ele_olde25max,"ele_olde25max/F");
  _mytree->Branch("ele_olde55",&ele_olde55,"ele_olde55/F");
  _mytree->Branch("ele_oldr9",&ele_oldr9,"ele_oldr9/F");
  _mytree->Branch("ele_oldcircularity",&ele_oldcircularity,"ele_oldcircularity/F");

  //
  // 
 //
  _mytree->Branch("ele_isEB",&ele_isEB,"ele_isEB/O");
  _mytree->Branch("ele_isEE",&ele_isEE,"ele_isEE/O");
  _mytree->Branch("ele_isEBEtaGap",&ele_isEBEtaGap,"ele_isEBEtaGap/O");
  _mytree->Branch("ele_isEBPhiGap",&ele_isEBPhiGap,"ele_isEBPhiGap/O");
  _mytree->Branch("ele_isEBEEGap", &ele_isEBEEGap);
  _mytree->Branch("ele_isEEDeeGap",&ele_isEEDeeGap,"ele_isEEDeeGap/O");
  _mytree->Branch("ele_isEERingGap",&ele_isEERingGap,"ele_isEERingGap/O");
 //
 //
  //_mytree->Branch("ele_dxy",&ele_dxy,"ele_dxy/F");
  //_mytree->Branch("ele_dxyB",&ele_dxyB,"ele_dxyB/F");
  //_mytree->Branch("ele_dzB",&ele_dzB,"ele_dzB/F");
  //_mytree->Branch("ele_dsz",&ele_dsz,"ele_dsz/F");
  //_mytree->Branch("ele_dszB",&ele_dszB,"ele_dszB/F");
  //
  //_mytree->Branch("ele_dzPV", &ele_dzPV, "ele_dzPV/F");
  //_mytree->Branch("ele_d0", &ele_d0, "ele_d0/F");
  //_mytree->Branch("ele_d0err",&ele_d0err, "ele_d0err/F");
  //
  //_mytree->Branch("ele_IP",&ele_IP,"ele_IP/F");
  //_mytree->Branch("ele_IPError",&ele_IPError,"ele_IPError/F");	
  //
  /*
  _mytree->Branch("ele_pfChargedHadIso",&ele_pfChargedHadIso,"ele_pfChargedHadIso/F");
  _mytree->Branch("ele_pfNeutralHadIso",&ele_pfNeutralHadIso,"ele_pfNeutralHadIso/F");
  _mytree->Branch("ele_pfPhotonIso",&ele_pfPhotonIso,"ele_pfPhotonIso/F");

  _mytree->Branch("ele_pfChargedIso", &ele_pfChargedIso, "ele_pfChargedIso/F");
  _mytree->Branch("ele_pfSumPUIso", &ele_pfSumPUIso, "ele_pfSumPUIso/F");
  */
  // 
  //
  _mytree->Branch("ele_sclRawE", &ele_sclRawE, "ele_sclRawE/F");
  _mytree->Branch("scl_E",    &ele_sclE);//,    "ele_sclE/F");
  _mytree->Branch("scl_Et",   &ele_sclEt);//,   "ele_sclEt/F");
  _mytree->Branch("scl_eta",  &scl_eta);//,  "scl_eta/F");
  _mytree->Branch("scl_phi",  &ele_sclPhi);//,  "ele_sclPhi/F");
  _mytree->Branch("ele_sclNclus",  &ele_sclNclus,  "ele_sclNclus/I");
  _mytree->Branch("ele_sclphiwidth", &ele_sclphiwidth, "ele_sclphiwidth/F");
  _mytree->Branch("ele_scletawidth", &ele_scletawidth, "ele_scletawidth/F");
  //
  _mytree->Branch("ele_psEoverEraw", &ele_psEoverEraw, "ele_psEoverEraw/F");
  _mytree->Branch("ele_oldsirir", &ele_oldsirir);
  //
  //_mytree->Branch("ele_ecalE",   &ele_ecalE,   "ele_ecalE/F");
  //_mytree->Branch("ele_ecalErr",   &ele_ecalErr,   "ele_ecalErr/F");

  _mytree->Branch("ele_HZZ_iso", &ele_HZZ_iso);

  _mytree->Branch("ele_oldhe",&ele_oldhe,"ele_oldhe/F");
  if(do_TLE) {
//    _mytree->Branch("pho_min_dR", &pho_min_dR);
    //_mytree->Branch("ele_POG_iso", &ele_POG_iso);
    _mytree->Branch("pho_phoPhotonIsolation_CITK", &pho_phoPhotonIsolation_CITK);
    _mytree->Branch("pho_phoNeutralHadronIsolation_CITK", &pho_phoNeutralHadronIsolation_CITK);
    _mytree->Branch("pho_phoChargedIsolation_CITK", &pho_phoChargedIsolation_CITK);

    //_mytree->Branch("ele_etOutsideMustache", &ele_etOutsideMustache);
    //_mytree->Branch("ele_nClusterOutsideMustache", &ele_nClusterOutsideMustache);
    _mytree->Branch("pho_mipNhitCone", &pho_mipNhitCone);
//    _mytree->Branch("pho_conversionTrackProvenance", &pho_conversionTrackProvenance);
    _mytree->Branch("pho_hasConversionTracks", &pho_hasConversionTracks);
    _mytree->Branch("pho_conversion_size", &pho_conversion_size);
    _mytree->Branch("pho_one_leg_conversion_size", &pho_one_leg_conversion_size);
    _mytree->Branch("pho_ID1_value", &pho_ID1_value);
    _mytree->Branch("pho_ID2_value", &pho_ID2_value);
    _mytree->Branch("pho_hasPixelSeed", &pho_hasPixelSeed);
    _mytree->Branch("pho_passElectronVeto", &pho_passElectronVeto);

  }

 
  if(!do_TLE) {
    _mytree->Branch("ele_echarge",&ele_echarge,"ele_echarge/I");

    _mytree->Branch("ele_ecalDrivenSeed",&ele_ecalDrivenSeed,"ele_ecalDrivenSeed/O");
    _mytree->Branch("ele_trackerDrivenSeed",&ele_trackerDrivenSeed,"ele_trackerDrivenSeed/O");
 
    _mytree->Branch("ele_valid_hits",&ele_valid_hits,"ele_valid_hits/I");
    _mytree->Branch("ele_lost_hits",&ele_lost_hits,"ele_lost_hits/I");
    _mytree->Branch("ele_gsfchi2",&ele_gsfchi2,"ele_gsfchi2/F");	

    _mytree->Branch("ele_kfchi2",&ele_kfchi2,"ele_kfchi2/F");
    _mytree->Branch("ele_kfhits",&ele_kfhits,"ele_kfhits/I");
    _mytree->Branch("ele_gsfhits",&ele_gsfhits,"ele_gsfhits/I");

    _mytree->Branch("ele_dz",&ele_dz,"ele_dz/F");
    _mytree->Branch("ele_SIP",&ele_SIP,"ele_SIP/F");
    _mytree->Branch("ele_conv_dist",&ele_conv_dist,"ele_conv_dist/F");
    _mytree->Branch("ele_conv_dcot",&ele_conv_dcot,"ele_conv_dcot/F");
    _mytree->Branch("ele_conv_radius",&ele_conv_radius,"ele_conv_radius/F");
    _mytree->Branch("ele_expected_inner_hits",&ele_expected_inner_hits,"ele_expected_inner_hits/I");
    _mytree->Branch("ele_vtxconv",&ele_vtxconv,"ele_vtxconv/I");
    _mytree->Branch("ele_conversionVertexFitProbability", &ele_conversionVertexFitProbability);
 
    _mytree->Branch("ele_fbrem",&ele_fbrem,"ele_fbrem/F");
    _mytree->Branch("ele_SCfbrem",&ele_SCfbrem,"ele_SCfbrem/F");
    _mytree->Branch("ele_nbrem",&ele_nbrem,"ele_nbrem/I");
    _mytree->Branch("ele_eClass",&ele_eClass,"ele_eClass/I");

    //_mytree->Branch("ele_eseedpout",&ele_eseedpout,"ele_eseedpout/F");
    _mytree->Branch("ele_ep",&ele_ep,"ele_ep/F");
    _mytree->Branch("ele_IoEmIop",&ele_IoEmIop);

    //_mytree->Branch("ele_eseedp",&ele_eseedp,"ele_eseedp/F");
    _mytree->Branch("ele_eelepout",&ele_eelepout,"ele_eelepout/F");
    //
    /*
    _mytree->Branch("ele_pin_mode",&ele_pin_mode,"ele_pin_mode/F");
    _mytree->Branch("ele_pout_mode",&ele_pout_mode,"ele_pout_mode/F");
    _mytree->Branch("ele_pTin_mode",&ele_pTin_mode,"ele_pTin_mode/F");
    _mytree->Branch("ele_pTout_mode",&ele_pTout_mode,"ele_pTout_mode/F");
    */
    //
    _mytree->Branch("ele_deltaetaseed",&ele_deltaetaseed,"ele_deltaetaseed/F");
    //_mytree->Branch("ele_deltaphiseed",&ele_deltaphiseed,"ele_deltaphiseed/F");
    //_mytree->Branch("ele_deltaetaele",&ele_deltaetaele,"ele_deltaetaele/F");
    //_mytree->Branch("ele_deltaphiele",&ele_deltaphiele,"ele_deltaphiele/F");
    _mytree->Branch("ele_deltaetain",&ele_deltaetain,"ele_deltaetain/F");
    _mytree->Branch("ele_deltaphiin",&ele_deltaphiin,"ele_deltaphiin/F");

    _mytree->Branch("ele_oldhebc",&ele_oldhebc,"ele_oldhebc/F");
 
    _mytree->Branch("ele_trackErr",  &ele_trackErr,  "ele_trackErr/F");
    _mytree->Branch("ele_combErr",   &ele_combErr,   "ele_combErr/F");
    _mytree->Branch("ele_PFcombErr", &ele_PFcombErr, "ele_PFcombErr/F");
 
    // Trigger
    //_mytree->Branch("trig_fired_names",&trig_fired_names,"trig_fired_names[10000]/C");
    //_mytree->Branch("event_trig_fired", &event_trig_fired);
    //_mytree->Branch("ele_trig_passed_filter", &ele_trig_passed_filter);
    //_mytree->Branch("ele_pass_hltEle27WP75GsfTrackIsoFilter", &ele_pass_hltEle27WP75GsfTrackIsoFilter);

    _mytree->Branch("ele_ID1", &ele_ID1);
    _mytree->Branch("ele_ID2", &ele_ID2);

    _mytree->Branch("ele_ID1_pass", &ele_ID1_pass);
    _mytree->Branch("ele_ID2_pass", &ele_ID2_pass);

    _mytree->Branch("ele_ID1_cat", &ele_ID1_cat);
    _mytree->Branch("ele_ID2_cat", &ele_ID2_cat);


    _mytree->Branch("ele_index", &ele_index);


  }   
  // Truth Leptons

  _mytree->Branch("mc_event_weight",&_mc_event_weight,"mc_event_weight/F");
  //
  _mytree->Branch ("MC_TrueNumInteractions",&_MC_TrueNumInteractions,"MC_TrueNumInteractions/I");

  /*
  _mytree->Branch("mc_ele_isPromptFinalState", &mc_ele_isPromptFinalState);
  _mytree->Branch("mc_ele_isDirectPromptTauDecayProductFinalState", &mc_ele_isDirectPromptTauDecayProductFinalState);
*/
  _mytree->Branch("mc_ele_matchedFromCB", &mc_ele_matchedFromCB);
  _mytree->Branch("mc_ele_matchedFromCB2", &mc_ele_matchedFromCB2);
  _mytree->Branch("mc_ele_matchMother_status", &mc_ele_matchMother_status);
  _mytree->Branch("mc_ele_matchMother_PDGID", &mc_ele_matchMother_PDGID);
  _mytree->Branch("mc_pho_mother_status", &mc_pho_mother_status);
  _mytree->Branch("mc_pho_mother_id", &mc_pho_mother_id);

//  _mytree->Branch("mc_ele_photon_over_ele_pt", &mc_ele_photon_over_ele_pt);

  //_mytree->Branch("ele_dr03EcalRecHitSumEt", &ele_dr03EcalRecHitSumEt);
  //_mytree->Branch("ele_dr03HcalTowerSumEt", &ele_dr03HcalTowerSumEt);
  //_mytree->Branch("ele_dr03TkSumPt", &ele_dr03TkSumPt);
  _mytree->Branch("ele_pt", &ele_pT);


 // _mytree->Branch("mc_gen_ele_p4", &_mc_gen_ele_p4);
  _mytree->Branch("mc_gen_pt", &mc_gen_pT);
  _mytree->Branch("mc_gen_eta", &mc_gen_eta);
  _mytree->Branch("mc_gen_ID", &mc_gen_ID);


  //_mytree->Branch("", &);
  //_mytree->Branch("electronEcalPFClusterIsolationProducer", &electronEcalPFClusterIsolationProducer);
  //_mytree->Branch("electronHcalPFClusterIsolationProducer", &electronHcalPFClusterIsolationProducer);
  _mytree->Branch("ele_full5x5_hcalOverEcal", &ele_full5x5_hcalOverEcal);

/* _mytree->Branch("ele_isEBEtaGap", &ele_isEBEtaGap);
 _mytree->Branch("ele_isEBPhiGap", &ele_isEBPhiGap);
 _mytree->Branch("ele_isEEDeeGap", &ele_isEEDeeGap);
 _mytree->Branch("ele_isEERingGap", &ele_isEERingGap);

*/

}


//
// member functions
//
// =============================================================================================
// ------------ method called for each event  ------------
void Ntuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
// =============================================================================================
{
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", builder);
  thebuilder = &*(builder.product());

  //load the conversion collection

  iEvent.getByToken(electronsToken_, electronsColl_h);

  if(do_TLE&&!isAOD) {
      iEvent.getByToken(photonsToken_, photonsColl_h);
      iEvent.getByToken(phoNeutralHadronIsolation_CITKToken_, phoNeutralHadronIsolation_CITK_map);
      iEvent.getByToken(phoChargedIsolation_CITKToken_, phoChargedIsolation_CITK_map);
      iEvent.getByToken(phoPhotonIsolation_CITKToken_, phoPhotonIsolation_CITK_map);
  }

  iEvent.getByToken(vertexToken_, primaryVertexColl_h);
  iEvent.getByToken(conversionsToken_, conversions_h);
  iEvent.getByToken(beamSpotToken_, beamspot_h);

  Vertex dummy;
  pv = &dummy;
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
  edm::Handle<double> rhoHandle;
  iEvent.getByToken(rhoForEleToken, rhoHandle);
  ele_rho = *rhoHandle;

 // iEvent.getByToken(genParticlesToken_CB, genCandidatesCollection);


  // The object that will compute 5x5 quantities  
  lazyToolnoZS = std::unique_ptr<noZS::EcalClusterLazyTools>(new noZS::EcalClusterLazyTools(iEvent, iSetup, 
                          ebReducedRecHitCollection_,
                          eeReducedRecHitCollection_,
                          esReducedRecHitCollection_));

  
  // get the beam spot
  beamSpot = &*(beamspot_h.product());



  if(!ID1_use_userFloat_) {
      iEvent.getByToken(electronID1Token_, ID1_map); 
      iEvent.getByToken(electronID1CatToken_, ele_ID1_cat_map); 
  }

  iEvent.getByToken(electronID2Token_, ID2_map); 

  iEvent.getByToken(electronID2CatToken_, ele_ID2_cat_map); 

  if(!isAOD) {
  iEvent.getByToken(photonID1Token_, pho_ID1_map);
  iEvent.getByToken(photonID2Token_, pho_ID2_map);
  }

  if(isMC_) {
    iEvent.getByToken(PUinfoToken, PupInfo);
    iEvent.getByToken(genEventInfoProductTagToken_, genEvtInfo );
//    iEvent.getByToken(genParticleToken_, genCandidatesCollection); //genParticlesPruned
    iEvent.getByToken(genParticlesToken_CB, genCandidatesCollection_CB);
  }
  FillEvent(iEvent, iSetup);
  LogDebug("") << "After FillEvent"; 
  FillVertices(iEvent, iSetup);
  LogDebug("") << "After FillVertices";
  FillMET (iEvent, iSetup);
  LogDebug("") << "After FillMET";

  if(!do_TLE) {
    for(size_t i_ele = 0;  i_ele <  electronsColl_h->size(); ++i_ele) {
      Init();
      const auto ielectron =  electronsColl_h->ptrAt(i_ele); 
      if(ielectron->pt() < 4.9) continue;
      //FillCandidate(ielectron);
      FillElectron(ielectron);
      LogDebug("") << "After FillElectrons";
      
      if(isMC_ ) {
       FillTruth(ielectron);
       LogDebug("") << "After FillTruth";
      }

      if(abs(scl_eta) > 2.6) continue; 
      if(do_signal && mc_ele_matchedFromCB == TRUE_PROMPT_ELECTRON) {
       _mytree->Fill();
      }

      if(!do_signal && (mc_ele_matchedFromCB == UNMATCHED || mc_ele_matchedFromCB == TRUE_NON_PROMPT_ELECTRON)) {
       _mytree->Fill();
      }
    }
  } else {
    for(size_t i_pho = 0;  i_pho <  photonsColl_h->size(); ++i_pho) {
      const auto iphoton =  photonsColl_h->ptrAt(i_pho);
      double min_dR = 999.;
      double curr_dR = 999.;
      for(size_t i_ele = 0;  i_ele <  electronsColl_h->size(); ++i_ele) {
          const auto ielectron =  electronsColl_h->ptrAt(i_ele);
          curr_dR = ROOT::Math::VectorUtil::DeltaR( ielectron->p4(), iphoton->p4() );

          if(curr_dR < min_dR) {
              min_dR = curr_dR;
          }
      }
      if(do_min_dR_cut && min_dR < 0.3) continue;

      FillPhoton(iphoton);

      if(isMC_ ) {
        FillTruth(iphoton);
        LogDebug("") << "After FillTruth";
      }

      if(abs(scl_eta) > 2.6) continue; 
      if(do_signal && mc_ele_matchedFromCB == TRUE_PROMPT_ELECTRON) {
       _mytree->Fill();
      }

      if(!do_signal && (mc_ele_matchedFromCB == UNMATCHED || mc_ele_matchedFromCB == TRUE_NON_PROMPT_ELECTRON)) {
       _mytree->Fill();
      }
    }
  }
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

/*
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
	  event_trig_fired = hlt_string);
      strcat(trig_fired_names_local,hlt_string.c_str());
	  strcat(trig_fired_names_local,"*");
	}
      }
    } // if HLT
  }
  strcpy(trig_fired_names,trig_fired_names_local);
 
if(false) { //inFileType == inputFileTypes::AOD) {
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
            _selectedObjects = foundObject);
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
            _hltEle27WP75GsfTrackIsoFilter = foundObject);
        }
    }
    }
}
*/

} // end of FillEvent


// =============================================================================================
void Ntuplizer::FillVertices(const edm::Event& iEvent, const edm::EventSetup& iSetup)
//=============================================================================================
{
  Handle<vector<reco::Vertex> >  recoPrimaryVertexCollection;
  iEvent.getByToken(vertexToken_, recoPrimaryVertexCollection);

  _vtx_N = recoPrimaryVertexCollection->size();
  
} // end of FillVertices


/*
// =============================================================================================
void Ntuplizer::FillCandidate(const edm::Ptr<reco::Candidate> icand)
//=============================================================================================
{
  LogDebug("") << "Ntuplizer::FillPhoton";
const reco::Candidate* icandidate = nullptr;

       icandidate = icand.get();
    if(icand->isElectron()) {
       LogError("") << "isElectron";
    }

    ele_isEB = icandidate->isEB(); 
    ele_isEE = icandidate->isEE(); 
    ele_isEBEEGap = icandidate->isEBEEGap();  
    ele_isEBEtaGap = icandidate->isEBEtaGap();  
    ele_isEBPhiGap = icandidate->isEBPhiGap();  
    ele_isEEDeeGap = icandidate->isEEDeeGap();  
    ele_isEERingGap = icandidate->isEERingGap();  
}
*/


// =============================================================================================
void Ntuplizer::FillPhoton(const edm::Ptr<reco::Photon> iphoton)
//=============================================================================================
{
  LogDebug("") << "Ntuplizer::FillPhoton";
//  const auto iphoton = photonsColl_h->ptrAt(i_pho); 

  pho_hasPixelSeed = iphoton->hasPixelSeed();
  const edm::Ptr<pat::Photon> phPatPtr(iphoton);
  pho_passElectronVeto = phPatPtr->passElectronVeto();
  ele_pT = iphoton->pt();

    ele_isEB = iphoton->isEB();
    ele_isEE = iphoton->isEE();
    ele_isEBEEGap = iphoton->isEBEEGap();
    ele_isEBEtaGap = iphoton->isEBEtaGap();
    ele_isEBPhiGap = iphoton->isEBPhiGap();
    ele_isEEDeeGap = iphoton->isEEDeeGap();
    ele_isEERingGap = iphoton->isEERingGap();

    // Shower Shape
    const auto& full5x5_pss = iphoton->full5x5_showerShapeVariables();
    ele_oldsigmaietaieta = full5x5_pss.sigmaIetaIeta;
    //ele_oldsigmaetaeta   = full5x5_pss.sigmaEtaEta;
    ele_oldsigmaiphiiphi = full5x5_pss.sigmaIphiIphi;
    ele_olde15           = full5x5_pss.e1x5;
    ele_olde25max        = full5x5_pss.e2x5Max;
    ele_olde55           = full5x5_pss.e5x5;
    ele_oldr9            = (float) iphoton->r9();
    ele_oldhe               = full5x5_pss.hcalDepth1OverEcal + full5x5_pss.hcalDepth2OverEcal; 

    //ele_oldhe            =  iphoton->full5x5_hcalOverEcal();
    //ele_oldhebc          =  iphoton->full5x5_hcalOverEcalBc();
    ele_oldcircularity   = ele_olde55 != 0. ? 1. - (ele_olde15 / ele_olde55) : -1.;
    ele_hadronicOverEm   = iphoton->hadronicOverEm();
 
    // -----------------------------------------------------------------
    // Get SuperCluster Informations
    // -----------------------------------------------------------------
    reco::SuperClusterRef sclRef = iphoton->superCluster();
   
    float R  = TMath::Sqrt(sclRef->x()*sclRef->x() + sclRef->y()*sclRef->y() +sclRef->z()*sclRef->z());
    float Rt = TMath::Sqrt(sclRef->x()*sclRef->x() + sclRef->y()*sclRef->y());
    ele_sclRawE  = sclRef->rawEnergy() ;
    
    ele_sclE     = sclRef->energy() ;
    ele_sclEt    = sclRef->energy()*(Rt/R) ;
    scl_eta   = sclRef->eta() ;
    ele_sclPhi   = sclRef->phi() ;
    ele_sclNclus = sclRef->clustersSize();
    

    ele_sclphiwidth = sclRef->phiWidth();
    ele_scletawidth = sclRef->etaWidth();
    // -----------------------------------------------------------------
    // Get PreShower Informations
    // -----------------------------------------------------------------
    ele_psE = sclRef->preshowerEnergy();
    ele_psEoverEraw = ele_sclRawE != 0. ? ele_psE / ele_sclRawE : -1; 
    if (iphoton->isEE()) {
        ele_oldsirir = lazyToolnoZS->eseffsirir( *(iphoton->superCluster()));
    }
 

    float fsr = 0;
    
    ele_HZZ_iso = isAOD ? -1 : LeptonIsoHelper::combRelIsoPF(sampleType, setup, ele_rho, *iphoton, fsr);
  //  LogPrint("") << "-Werror=unused-but-set-variabl:" << rhoForEle;
    //float pog_iso = 0;//(*ISO_map)[iphotons];

    float phoChargedIsolation_CITK = isAOD ? -1 :(*phoChargedIsolation_CITK_map)[iphoton];
    float phoNeutralHadronIsolation_CITK = isAOD ? -1 : (*phoNeutralHadronIsolation_CITK_map)[iphoton];
    float phoPhotonIsolation_CITK = isAOD ? -1 : (*phoPhotonIsolation_CITK_map)[iphoton];

    pho_phoPhotonIsolation_CITK = phoPhotonIsolation_CITK;
    pho_phoNeutralHadronIsolation_CITK = phoNeutralHadronIsolation_CITK;
    pho_phoChargedIsolation_CITK = phoChargedIsolation_CITK;

    pho_mipNhitCone = iphoton->mipNhitCone();
//    pho_conversionTrackProvenance;
    pho_conversion_size = iphoton->conversions().size();
    pho_one_leg_conversion_size = iphoton->conversionsOneLeg().size();

    pho_hasConversionTracks = iphoton->hasConversionTracks();
   //ele_POG_iso = pog_iso;

    pho_ID1_value = isAOD ? -1 : (*pho_ID1_map)[iphoton];
    pho_ID2_value = isAOD ? -1 : (*pho_ID2_map)[iphoton];
   //    ele_ID1_cat = (*ele_ID1_cat_map)[ielectron];
    

//    ele_ID2 = ele_ID2_value;




}
// =============================================================================================
void Ntuplizer::FillElectron(const edm::Ptr<reco::GsfElectron> ielectron)
//=============================================================================================
{
  LogDebug("") << "Ntuplizer::FillElectrons";


  ele_N = electronsColl_h->size();
  //ele_N_saved = 0;


    LogDebug("") << "Processing new electron";


    const pat::Electron* pat_ele = nullptr;
    
    if(!isAOD) {
      const edm::Ptr<pat::Electron> elePatPtr(ielectron);
      if(elePatPtr.get() == NULL) {
          LogError("") << "Failed to get pointer to pat electron!";
      } else {
          pat_ele = elePatPtr.get();
      }
    }

 //   ++ele_N_saved;


    float ele_ID1_value = 999.;
    int ele_ID1_cat = 999.; 

    float ele_ID2_value = (*ID2_map)[ielectron];
    int ele_ID2_cat = (*ele_ID2_cat_map)[ielectron];


    if(!isAOD&&ID1_use_userFloat_) {
        ele_ID1_value = pat_ele->userFloat(electronID1_name + "Values"); 
        ele_ID1_cat = pat_ele->userInt(electronID1_name + "Categories"); 
    } else {
        ele_ID1_value = (*ID1_map)[ielectron];
        ele_ID1_cat = (*ele_ID1_cat_map)[ielectron];
    }

    ele_ID1 = ele_ID1_value;
    ele_ID2 = ele_ID2_value;

    ele_ID1_cat = ele_ID1_cat;
    ele_ID2_cat = ele_ID2_cat;

    //ele_ID1_pass = (*ID1_pass_map)[ielectron]);
    //ele_ID2_pass = (*ID2_pass_map)[ielectron]);


//    ele_index = i_ele;

    float fsr = 0;

    ele_HZZ_iso = isAOD ?-1 :LeptonIsoHelper::combRelIsoPF(sampleType, setup, ele_rho, *ielectron, fsr);
  //  LogPrint("") << "-Werror=unused-but-set-variabl:" << rhoForEle;
  //      ele_HZZ_iso = iso;
  //

    ele_pT = ielectron->pt();
    ele_trackMomentumAtVtx_R = ielectron->trackMomentumAtVtx().R();

    ele_full5x5_hcalOverEcal = ielectron->full5x5_hcalOverEcal();

//    setMomentum(myvector, ielectron->p4());
//    new (electrons[counter]) TLorentzVector(myvector);

    ele_echarge = ielectron->charge(); 
//    ele_rho = 

    //
    //ele_he        = ielectron->hcalOverEcal(); //hadronicOverEm() ;
    //ele_hebc      = ielectron-> hcalOverEcalBc();

    // TrackCluster Matching
    ele_eseedpout = ielectron->eSeedClusterOverPout();
    ele_ep        = ielectron->eSuperClusterOverP() ;        
    ele_eseedp    = ielectron->eSeedClusterOverP() ;         
    ele_eelepout  = ielectron->eEleClusterOverPout() ;       
    //
    ele_pin_mode    = ielectron->trackMomentumAtVtx().R() ; 
    ele_pout_mode   = ielectron->trackMomentumOut().R() ; 
    ele_pTin_mode   = ielectron->trackMomentumAtVtx().Rho() ; 
    ele_pTout_mode  = ielectron->trackMomentumOut().Rho() ; 
    //
    ele_deltaetaseed = ielectron->deltaEtaSeedClusterTrackAtCalo() ; 
    ele_deltaphiseed = ielectron->deltaPhiSeedClusterTrackAtCalo() ;  
    ele_deltaetaele  = ielectron->deltaEtaEleClusterTrackAtCalo() ;  
    ele_deltaphiele  = ielectron->deltaPhiEleClusterTrackAtCalo() ; 
    ele_deltaetain   = ielectron->deltaEtaSuperClusterTrackAtVtx();
    ele_deltaphiin   = ielectron->deltaPhiSuperClusterTrackAtVtx();   

    // Shower Shape
    /*
    ele_sigmaietaieta = (ielectron->showerShape()).sigmaIetaIeta; //  ielectron->
    ele_sigmaetaeta   = (ielectron->showerShape()).sigmaEtaEta ; //ielectron->sigmaEtaEta() ;
    ele_sigmaiphiiphi = (ielectron->showerShape()).sigmaIphiIphi;
    ele_e15           = (ielectron->showerShape()).e1x5; // ;ielectron->e1x5() ;
    ele_e25max        = (ielectron->showerShape()).e2x5Max; // ;ielectron->e2x5Max() ;
    ele_e55           = (ielectron->showerShape()).e5x5; // ;ielectron->e5x5() ;
    ele_r9            = (ielectron->showerShape()).r9; 
    */
    //
   
    //ele_oldsigmaetaeta   =  ielectron->full5x5_sigmaEtaEta();    
    ele_oldsigmaietaieta =  ielectron->full5x5_sigmaIetaIeta();  
    ele_oldsigmaiphiiphi =  ielectron->full5x5_sigmaIphiIphi();  
    ele_oldr9              =  ielectron->full5x5_r9();  
    ele_olde15           =  ielectron->full5x5_e1x5();
    ele_olde25max   =  ielectron->full5x5_e2x5Max();
    ele_olde55           =  ielectron->full5x5_e5x5();
    ele_oldhe             =  ielectron->full5x5_hcalOverEcal();
    ele_oldhebc        =  ielectron->full5x5_hcalOverEcalBc();
    ele_oldcircularity = ele_olde55 != 0. ? 1. - (ele_olde15 / ele_olde55) : -1.;
    ele_hadronicOverEm = ielectron->hadronicOverEm();
    // E/P combination
    ele_ecalE     = ielectron->ecalEnergy();
    ele_ecalErr   = ielectron->ecalEnergyError();
    ele_trackErr  = ielectron->trackMomentumError();
    ele_combErr   = ielectron->p4Error(GsfElectron::P4_COMBINATION);
    ele_PFcombErr = ielectron->p4Error(GsfElectron::P4_PFLOW_COMBINATION);
    ele_IoEmIop   = -1;
    if(ele_ecalE != 0 || ele_pin_mode != 0) {
      ele_IoEmIop = 1.0 / ele_ecalE - (1.0 / ele_pin_mode);
    }
   //
    //
    if (ielectron->isEE()) {
        ele_oldsirir = lazyToolnoZS->eseffsirir( *(ielectron->superCluster()));
    }
    ele_isEB = ielectron->isEB();
    ele_isEE = ielectron->isEE();
    ele_isEBEEGap = ielectron->isEBEEGap();
    ele_isEBEtaGap = ielectron->isEBEtaGap();
    ele_isEBPhiGap = ielectron->isEBPhiGap();
    ele_isEEDeeGap = ielectron->isEEDeeGap();
    ele_isEERingGap = ielectron->isEERingGap();
    ele_ecalDrivenSeed = ielectron->ecalDrivenSeed();
    ele_trackerDrivenSeed = ielectron->trackerDrivenSeed();
    //

    // -----------------------------------------------------------------
    // Tracking Variables
    // -----------------------------------------------------------------
    LogDebug("") << "Start looking at track vars";
    ele_lost_hits   = ielectron->gsfTrack()->lost(); //numberOfLostHits();
    ele_valid_hits = ielectron->gsfTrack()->found(); //numberOfValidHits() ;
    ele_gsfchi2    = ielectron->gsfTrack()->normalizedChi2() ;
    LogDebug("") << "After gsfTrack";
   
    ele_gsfhits = ielectron->gsfTrack()->hitPattern().trackerLayersWithMeasurement();

    bool validKF=false;
    reco::TrackRef myTrackRef = ielectron->closestCtfTrackRef();
    //const edm::Ptr<pat::Electron> elePatPtr(eleRecoPtr);
    // Check if this is really a pat::Electron, and if yes, get the track ref from this new
    // pointer instead
    if( pat_ele != nullptr )
        myTrackRef = pat_ele->closestCtfTrackRef();


    //LogVerbatim("") << "After clostste ctdf";

    validKF = myTrackRef.isAvailable();
    //LogVerbatim("") << "isAvailable : " << validKF;

    validKF &= myTrackRef.isNonnull();
    //LogVerbatim("") << "isNonnull&isAvailable : " << validKF;
  
    ele_kfchi2 = validKF ? myTrackRef->normalizedChi2() : 0 ; //ielectron->track()->normalizedChi2() : 0 ;
    ele_kfhits = validKF ? myTrackRef->hitPattern().trackerLayersWithMeasurement() : -1. ; //ielectron->track()->hitPattern().trackerLayersWithMeasurement() : 0 ;

    LogDebug("") << "After ctfTrack";
    //
    ele_dxy  = ielectron->gsfTrack()->dxy() ;
    ele_dz   = ielectron->gsfTrack()->dz() ;
    ele_dsz  = ielectron->gsfTrack()->dsz() ;
    ele_dzPV = ielectron->gsfTrack()->dz(pv->position());
    //

    // Conversion Rejection
    ele_conv_dcot   = ielectron->convDist();
    ele_conv_dist   = ielectron->convDcot();
    ele_conv_radius = ielectron->convRadius();
    
    ele_expected_inner_hits = ielectron->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS);
    
   

    reco::ConversionRef conv_ref = ConversionTools::matchedConversion(*ielectron, conversions_h, beamSpot->position()); 
    bool vtxFitConversion = ConversionTools::hasMatchedConversion(*ielectron, conversions_h, beamSpot->position());
    ele_vtxconv = vtxFitConversion;

    float vertexFitProbability = -1.;
    if(!conv_ref.isNull()) {
        const reco::Vertex &vtx = conv_ref.get()->conversionVertex();
        if (vtx.isValid()) {
            vertexFitProbability = TMath::Prob( vtx.chi2(),  vtx.ndof());
        }
    }
    ele_conversionVertexFitProbability = vertexFitProbability;
//    LogVerbatim("") << "Vertex fit probability: " << vertexFitProbability;
    // 
   /* 
    // -----------------------------------------------------------------
    // PF Isolation for electrons in the cone 0.4
    // -----------------------------------------------------------------
    ele_pfChargedHadIso[counter]   = (ielectron->pfIsolationVariables()).sumChargedHadronPt ; //chargedHadronIso();
    ele_pfNeutralHadIso[counter]   = (ielectron->pfIsolationVariables()).sumNeutralHadronEt ; //neutralHadronIso();
    ele_pfPhotonIso[counter]       = (ielectron->pfIsolationVariables()).sumPhotonEt; //photonIso();
    //
    ele_pfChargedIso[counter]      = (ielectron->pfIsolationVariables()).sumChargedParticlePt;
    ele_pfSumPUIso[counter]        = (ielectron->pfIsolationVariables()).sumPUPt;
   */
    // -----------------------------------------------------------------
    // SIP3D
    // -----------------------------------------------------------------
    //default values for IP3D
    float ele_ip3D     = -999.0; 
    float ele_ip3D_err = 999.; //fMVAVar_ip3dSig = 0.0;
    float d0_corr = 999;
    float d0_err   = 999;

    if (ielectron->gsfTrack().isNonnull()) {
      const float gsfsign   = ( (-ielectron->gsfTrack()->dxy(pv->position()))   >=0 ) ? 1. : -1.;
      
      const reco::TransientTrack &tt = thebuilder->build(ielectron->gsfTrack()); //transientTrackBuilder.build(ielectron->gsfTrack()); 
      const std::pair <bool,Measurement1D> &ip3dpv =  IPTools::absoluteImpactParameter3D(tt, *pv);
      std::pair<bool,Measurement1D> result   = IPTools::absoluteTransverseImpactParameter(tt, * pv);
      
      if (ip3dpv.first) {
	float ip3d = gsfsign*ip3dpv.second.value();
	float ip3derr = ip3dpv.second.error();  
	
	ele_ip3D     = ip3d; 
	ele_ip3D_err = ip3derr;
      } // if ip3dpv.first

      if(result.first) {
	d0_corr = result.second.value();
	d0_err   = result.second.error();
      } // if d0
    } // if gsftrack.isNonnull
    
    ele_IP      = ele_ip3D; //fabs(ielectron->dB(pat::Electron::PV3D));
    ele_IPError = ele_ip3D_err; //ielectron->edB(pat::Electron::PV3D);	
    float ele_sip = -999; if(ele_ip3D_err!=0) ele_sip = ele_ip3D / ele_ip3D_err;
    ele_SIP     = ele_sip; //ele_IP[counter]/ele_IPError[counter];
    
    ele_d0 = d0_corr;
    ele_d0err = d0_err;

    // -----------------------------------------------------------------
    // Get SuperCluster Informations
    // -----------------------------------------------------------------
    reco::SuperClusterRef sclRef = ielectron->superCluster();
   
    float R  = TMath::Sqrt(sclRef->x()*sclRef->x() + sclRef->y()*sclRef->y() +sclRef->z()*sclRef->z());
    float Rt = TMath::Sqrt(sclRef->x()*sclRef->x() + sclRef->y()*sclRef->y());
    ele_sclRawE  = sclRef->rawEnergy() ;
    
    ele_sclE     = sclRef->energy() ;
    ele_sclEt    = sclRef->energy()*(Rt/R) ;
    scl_eta   = sclRef->eta() ;
    ele_sclPhi   = sclRef->phi() ;
    ele_sclNclus = sclRef->clustersSize();
    

    ele_sclphiwidth = sclRef->phiWidth();
    ele_scletawidth = sclRef->etaWidth();

    // -----------------------------------------------------------------
    // Get PreShower Informations
    // -----------------------------------------------------------------
    ele_psE = sclRef->preshowerEnergy();
    ele_psEoverEraw = ele_sclRawE != 0. ? ele_psE / ele_sclRawE : -1; 
    // -----------------------------------------------------------------
    //fbrem
    // -----------------------------------------------------------------
    ele_fbrem      = ielectron->fbrem();
    ele_SCfbrem    = ielectron->superClusterFbrem();
    // GsfElectron definition changed in 7X, 
    ele_eClass     = ielectron->classification() ;
    ele_nbrem      = ielectron->numberOfBrems();
    
    // -----------------------------------------------------------------
    // Electron ID electronsCol
    // -----------------------------------------------------------------

} // end of FillElectrons

// ====================================================================================
void Ntuplizer::FillMET (const edm::Event& iEvent, const edm::EventSetup& iSetup)
// ====================================================================================
{
	

  //edm::Handle< edm::View<reco::MET> > pfMEThandle;
  //iEvent.getByToken(pfMETToken_, pfMEThandle);
  
} // end of Fill MET



// ====================================================================================
void Ntuplizer::FillTruth(const edm::Ptr<reco::Candidate> icandidate)
// ====================================================================================
{
  LogDebug("") << "Ntuplizer::FillTruth"; 

  
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


  _mc_event_weight = genEvtInfo->weight();
  LogDebug("MC") << "MC generator event weight: " << _mc_event_weight;

  //ElectronMatchType
  int  matchType = matchToTruth(icandidate, genCandidatesCollection_CB); 
  int  matchType2 = matchToTruth(icandidate, genCandidatesCollection_CB, 0.03); 

  mc_ele_matchedFromCB = matchType;
  mc_ele_matchedFromCB2 = matchType2;

  mc_ele_photon_over_ele_pt = photon_E_over_electron_E(icandidate, genCandidatesCollection_CB);

  const reco::Candidate* genParticle = nullptr; 
  int mother_ID = 0;
  int mother_status = -1;

  genParticle = GetClosestGenParticle(icandidate, genCandidatesCollection_CB);

  if(genParticle != nullptr) {
      findFirstNonElectronMother(genParticle, mother_ID, mother_status);
      findFirstNonPhotonMother(genParticle, mc_pho_mother_id, mc_pho_mother_status);
  }
  
  float gen_pT = -1;
  float gen_eta = -999;
  float gen_ID = 0;
  if(genParticle != nullptr) {
      gen_pT = genParticle->p4().Pt();
      gen_eta = genParticle->p4().Eta();
      gen_ID = genParticle->pdgId();
  }
  mc_ele_matchMother_PDGID = mother_ID; 
  mc_ele_matchMother_status = mother_status; 

  mc_gen_pT = gen_pT;
  mc_gen_eta = gen_eta;
  mc_gen_ID = gen_ID;
    // To be re-implemented
    //mc_ele_isPromptFinalState = p->isPromptFinalState());
    //mc_ele_isDirectPromptTauDecayProductFinalState = p->isDirectPromptTauDecayProductFinalState());	

} // end of FillTruth




// ------------ method called once each job just after ending the event loop  ------------
void Ntuplizer::endJob() {
  file->cd();
  file->Write();
  file->Close();

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
  ele_rho = 0; 
  //electrons
  ele_echarge = 0 ;
  //ele_he = 0 ;  
  //ele_hebc = 0 ;  
  ele_oldhe = 0 ;  
  ele_oldhebc = 0 ;  
  ele_hadronicOverEm = 0;
  ele_eseedpout = 0 ;  
  ele_ep = 0 ;  
  ele_eseedp = 0 ;  
  ele_eelepout = 0 ;        
  ele_deltaetaseed = 0 ;  
  ele_deltaetaele = 0 ;  
  ele_deltaphiseed = 0 ;  
  ele_deltaphiele = 0 ;  
  ele_deltaetain = 0 ;  
  ele_deltaphiin = 0 ; 
  /*
  ele_sigmaietaieta = 0 ;  
  ele_sigmaiphiiphi = 0 ;  
  ele_sigmaetaeta = 0 ;  
  ele_e15 = 0 ;  
  ele_e25max = 0 ;  
  ele_e55 = 0 ;  
  ele_e1 = 0 ;  
  ele_r9 = 0;
  */
  //
  //ele_oldsigmaetaeta   = 0;
  ele_oldsigmaietaieta = 0;
  ele_oldsigmaiphiiphi = 0;
  ele_oldsigmaietaiphi = 0; 
  ele_oldr9            = 0; 
  ele_olde15           = 0; 
  ele_olde25max        = 0; 
  ele_olde55           = 0; 
  ele_oldcircularity   = 0;

  //
  ele_pin_mode = 0 ;  
  ele_pout_mode = 0 ;  
  ele_pTin_mode = 0 ;  
  ele_pTout_mode = 0 ;  
  //
  ele_fbrem = 0 ;  
  ele_SCfbrem = 0 ;  
  ele_nbrem = 0;
  //
  ele_isEB = 0 ;  
  ele_isEE = 0 ;  
  ele_isEBEEGap = 0 ;  
  ele_isEBEtaGap = 0 ;  
  ele_isEBPhiGap = 0 ;  
  ele_isEEDeeGap = 0 ;  
  ele_isEERingGap = 0 ; 
  ele_ecalDrivenSeed = 0 ;  
  ele_trackerDrivenSeed = 0 ; 
  ele_eClass = 0 ;
  ele_valid_hits = 0 ; 
  ele_lost_hits = 0 ; 
  ele_gsfchi2 = 0 ; 
  ele_dxyB = 0 ; 
  ele_dxy = 0 ; 
  ele_dzB = 0 ; 
  ele_dz = 0 ; 
  ele_dszB = 0 ; 
  ele_dsz = 0 ;              

  ele_conv_dcot = 0 ;
  ele_conv_dist = 0 ;
  ele_conv_radius = 0;
  ele_expected_inner_hits = 0;
  ele_vtxconv = 0;
  //
  ele_pfChargedHadIso   = 0;
  ele_pfNeutralHadIso   = 0;
  ele_pfPhotonIso       = 0;

  ele_pfChargedIso = 0;
  ele_pfSumPUIso   = 0;
  
  ele_IP = 0 ; 
  ele_IPError = 0 ; 
  ele_SIP = 0 ; 

  ele_dzPV = 0;
  ele_d0 = 0;
  ele_d0err = 0;

  //
  ele_sclRawE=0;
  ele_sclE=0;
  ele_sclEt=0;
  scl_eta=0;
  ele_sclPhi=0;
  ele_sclNclus=0;
  ele_sclphiwidth=0;
  ele_scletawidth=0;
 
  ele_ecalE=0;
  ele_ecalErr=0;
  ele_trackErr=0;
  ele_combErr=0;
  ele_PFcombErr=0;
  ele_IoEmIop = 0; 
 
  ele_ecalRegressionEnergy  = 0;
  ele_ecalRegressionError = 0;
  ele_ecalTrackRegressionEnergy  = 0;
  ele_ecalTrackRegressionError  = 0;
  ele_ecalScale  = 0;
  ele_ecalSmear  = 0;
  ele_ecalRegressionScale  = 0;
  ele_ecalRegressionSmear  = 0;
  ele_ecalTrackRegressionScale  = 0;
  ele_ecalTrackRegressionSmear  = 0;
  

  ele_kfchi2=0;
  ele_kfhits=-1;
  ele_gsfhits = -1;
  ele_psE = 0.;
  ele_psEoverEraw = 0.;
  ele_conversionVertexFitProbability = 0;
  ele_pT = 0;
  ele_trackMomentumAtVtx_R = 0;


  //ele_trig_passed_filter = 0;
  //ele_pass_hltEle27WP75GsfTrackIsoFilter = 0;
  ele_full5x5_hcalOverEcal = 0;


  ele_ID1 = 0;
  ele_ID2 = 0;
  ele_index = 0;
  ele_ID1_pass = 0;
  ele_ID2_pass = 0;
  ele_ID1_cat = 0;
  ele_ID2_cat = 0;
 
  ele_oldsirir = 0;

  _MC_TrueNumInteractions = 0;
  mc_ele_isPromptFinalState = 0;
  mc_ele_isDirectPromptTauDecayProductFinalState = 0;
  //_mc_gen_ele_p4 = 0;
  mc_gen_pT = 0;
  mc_gen_eta = 0;
  mc_ele_matchedFromCB = 0;
  mc_ele_matchedFromCB2 = 0;

  mc_ele_matchMother_PDGID = 0;
  mc_ele_matchMother_status = 0;

  mc_pho_mother_id = 0;
  mc_pho_mother_status = 0;

  mc_ele_photon_over_ele_pt = 0;


  pho_hasPixelSeed = false;
  pho_passElectronVeto = false;


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
