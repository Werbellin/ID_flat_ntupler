#ifndef Ntuplizer_H
#define Ntuplizer_H

// CMSSW
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include <DataFormats/HepMCCandidate/interface/GenParticle.h>
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"

#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"

#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include <FWCore/Framework/interface/ESHandle.h>


#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"


#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"


// ROOT
#include "TTree.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TVector3.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "Math/VectorUtil.h"

// C++
#include<memory>
#include<vector>

using namespace std;
using namespace edm;
using namespace reco;

enum class inputFileTypes {AOD, MINIAOD, UNDEF};

class Ntuplizer : public edm::EDAnalyzer {
   public:
      explicit Ntuplizer(const edm::ParameterSet&);
      ~Ntuplizer();

      typedef math::XYZTLorentzVector LorentzVector ;
      
      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

      void Init();
      void FillEvent(const edm::Event&, const edm::EventSetup&);
      void FillCandidate(const edm::Ptr<reco::Candidate> icandidate);
      void FillElectron(const edm::Ptr<reco::GsfElectron>);
      void FillPhoton(const edm::Ptr<reco::Photon> iphoton);
      void FillVertices(const edm::Event&, const edm::EventSetup&);
      void FillTruth(const edm::Ptr<reco::Candidate>);
      void FillMET (const edm::Event& iEvent, const edm::EventSetup& iSetup);
      
      void setMomentum(TLorentzVector & myvector, const LorentzVector & mom) ;


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
      const edm::ParameterSet&  conf;
        
      inputFileTypes inFileType;   
      bool isAOD; 
      edm::EDGetToken electronsToken_;
      edm::EDGetToken photonsToken_;
      edm::Handle<edm::View<reco::Photon>> photonsColl_h;

      edm::EDGetTokenT<edm::View<reco::GenParticle> > genParticlesToken_;
      edm::EDGetTokenT<reco::VertexCollection> vertexToken_;
      edm::EDGetTokenT<reco::ConversionCollection> conversionsToken_;
      edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
      edm::EDGetToken pfMETToken_;
      edm::EDGetTokenT<vector<reco::GenParticle> > genParticleToken_;
      edm::EDGetTokenT<edm::View<reco::GenParticle> > genParticlesToken_CB;

      edm::EDGetTokenT<GenEventInfoProduct> genEventInfoProductTagToken_;

      edm::EDGetTokenT<ValueMap<float>> electronID1Token_;
      edm::EDGetTokenT<ValueMap<float>> electronID2Token_;

      edm::EDGetTokenT<ValueMap<float>> photonID1Token_;
      edm::Handle<edm::ValueMap<float> > pho_ID1_map;

      edm::EDGetTokenT<ValueMap<float>> photonID2Token_;
      edm::Handle<edm::ValueMap<float> > pho_ID2_map;


      edm::EDGetTokenT<ValueMap<int>> electronID1CatToken_;
      edm::EDGetTokenT<ValueMap<int>> electronID2CatToken_;



      edm::EDGetTokenT<ValueMap<bool>> electronID1_pass_Token_;
      edm::EDGetTokenT<ValueMap<bool>> electronID2_pass_Token_;

      string electronID1_name;
      string electronID2_name;

      string photonID1_name;
      string photonID2_name;

      std::unique_ptr<noZS::EcalClusterLazyTools> lazyToolnoZS;
      edm::EDGetTokenT<EcalRecHitCollection> ebReducedRecHitCollection_;
      edm::EDGetTokenT<EcalRecHitCollection> eeReducedRecHitCollection_;
      edm::EDGetTokenT<EcalRecHitCollection> esReducedRecHitCollection_;

      edm::ESHandle<TransientTrackBuilder> builder;
      const TransientTrackBuilder* thebuilder;

      edm::Handle<edm::View<reco::GsfElectron>> electronsColl_h;
      edm::Handle<reco::VertexCollection> primaryVertexColl_h;
      edm::Handle<reco::ConversionCollection> conversions_h;
      edm::Handle<reco::BeamSpot> beamspot_h;
      const reco::BeamSpot* beamSpot; // = *(beamspot_h.product());
      edm::Handle<edm::ValueMap<float> > ID1_map;
      edm::Handle<edm::ValueMap<int> > ele_ID1_cat_map;

      edm::EDGetTokenT<ValueMap<float>> phoChargedIsolation_CITKToken_;
      edm::EDGetTokenT<ValueMap<float>> phoNeutralHadronIsolation_CITKToken_;
      edm::EDGetTokenT<ValueMap<float>> phoPhotonIsolation_CITKToken_;

  edm::Handle<edm::ValueMap<float> > phoChargedIsolation_CITK_map;
  edm::Handle<edm::ValueMap<float> > phoPhotonIsolation_CITK_map;
  edm::Handle<edm::ValueMap<float> > phoNeutralHadronIsolation_CITK_map;
     edm::Handle<edm::ValueMap<float> > ID2_map;

      edm::Handle<edm::ValueMap<int> > ele_ID2_cat_map;

      edm::Handle<edm::ValueMap<bool> > ID1_pass_map;
      edm::Handle<edm::ValueMap<bool> > ID2_pass_map;
 
      Handle<GenEventInfoProduct> genEvtInfo;
      edm::Handle<vector<reco::GenParticle> > genCandidatesCollection;
      edm::Handle<edm::View<reco::GenParticle> > genCandidatesCollection_CB;
      Handle<std::vector< PileupSummaryInfo > >  PupInfo;
 

      edm::EDGetTokenT<double> rhoForEleToken;
      // Trigger Stuff
      edm::EDGetTokenT<edm::TriggerResults> HLTToken;
      bool isMC_;	
      bool ID1_use_userFloat_;

      bool do_TLE;
      bool do_signal;
      bool do_min_dR_cut;

      int sampleType, setup;

      vector<trigger::TriggerObject> _selectedObjects;
      vector<trigger::TriggerObject> _hltEle27WP75GsfTrackIsoFilter;
      
      edm::InputTag PileupSrc_;	
      edm::EDGetTokenT<std::vector<PileupSummaryInfo> >  PUinfoToken; 
      //tree
      TFile* file;
      TTree *_mytree;
      string outputFile;
      string outputPath;
      TLorentzVector myvector ;  

      //global variables
      int _nEvent, _nRun, _nLumi;
      //pile-up
      int _PU_N;

      float _mc_event_weight;

      typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float>> NewLorentzVector;
//      vector<NewLorentzVector> _mc_gen_ele_p4;

      //vertices
      int _vtx_N;
      const reco::Vertex *pv;
 
      //trigger fired names
      char trig_fired_names[10000];

      // MET
      float _met_pf_et,_met_pf_px, _met_pf_py, _met_pf_phi, _met_pf_set, _met_pf_sig; 

      float ele_rho;

      //electrons
      int ele_N;
      int ele_N_saved;
      TClonesArray * m_electrons;


      float ele_conversionVertexFitProbability;
      int  mc_ele_isPromptFinalState;
      int  mc_ele_isDirectPromptTauDecayProductFinalState;
      int mc_ele_matchedFromCB;
      int mc_ele_matchedFromCB2;
      int mc_ele_matchMother_status;

      int mc_ele_matchMother_PDGID;
      float mc_ele_photon_over_ele_pt;
        
      int mc_pho_mother_id;
      int mc_pho_mother_status;
      float mc_gen_pT;
      float mc_gen_eta;
      float mc_gen_ID;

      float ele_pT;
      float ele_trackMomentumAtVtx_R;
      float ele_ID1;
      float ele_ID2;

      int ele_ID1_pass;
      int ele_ID2_pass;

      int ele_ID1_cat;
      int ele_ID2_cat;

      int ele_index;
      
      bool ele_isEE;
      bool ele_isEB;
      bool ele_isEBEEGap;
      bool ele_isEBEtaGap;
      bool ele_isEBPhiGap;
      bool ele_isEEDeeGap;
      bool ele_isEERingGap;
      bool ele_ecalDrivenSeed;
      bool ele_trackerDrivenSeed;


      float ele_oldsirir;


      string event_trig_fired;
      bool ele_trig_passed_filter;
      bool ele_pass_hltEle27WP75GsfTrackIsoFilter;
      float ele_full5x5_hcalOverEcal;

      int ele_echarge;
      //float ele_he, ele_hebc
      float  ele_eseedpout , ele_ep , ele_eseedp , ele_eelepout ;       
      float ele_deltaetaseed , ele_deltaetaele , ele_deltaphiseed , ele_deltaphiele , ele_deltaetain , ele_deltaphiin ;
      //float ele_sigmaietaieta , ele_sigmaetaeta , ele_sigmaiphiiphi, ele_e15 , ele_e25max , ele_e55 , ele_e1, ele_r9 ;

      //float ele_oldsigmaetaeta
      float  ele_oldsigmaietaieta, ele_oldsigmaiphiiphi, ele_oldsigmaietaiphi, ele_oldr9, ele_olde15,           
	ele_olde25max, ele_olde55;
      float ele_oldcircularity;
      float ele_oldhe, ele_oldhebc, ele_hadronicOverEm;


      float ele_pin_mode , ele_pout_mode , ele_pTin_mode , ele_pTout_mode ; 
      //
      float ele_fbrem,ele_SCfbrem;
      int ele_nbrem;
      //
      //
	  int ele_eClass;
      int ele_valid_hits, ele_lost_hits; 
      float ele_gsfchi2; 
      float ele_dxyB, ele_dxy, ele_dzB, ele_dz, ele_dszB, ele_dsz;              
      float ele_tkSumPt_dr03 , ele_ecalRecHitSumEt_dr03 , ele_hcalDepth1TowerSumEt_dr03 , ele_hcalDepth2TowerSumEt_dr03 ,
	ele_tkSumPt_dr04 , ele_ecalRecHitSumEt_dr04 , ele_hcalDepth1TowerSumEt_dr04 , ele_hcalDepth2TowerSumEt_dr04 ;
      //
      float ele_conv_dcot;
      float ele_conv_dist;
      float ele_conv_radius;
      int ele_expected_inner_hits;
      int ele_vtxconv;
      //
      float ele_pfChargedHadIso, ele_pfNeutralHadIso, ele_pfPhotonIso, ele_pfChargedIso, ele_pfSumPUIso;
      float ele_HZZ_iso;
      //
      float ele_dzPV, ele_d0, ele_d0err;
      float ele_IP, ele_IPError, ele_SIP ;
      float ele_sclE, ele_sclEt, scl_eta, ele_sclPhi, ele_sclRawE;
      int ele_sclNclus;
      
    
      float ele_ecalE, ele_ecalErr, ele_trackErr, ele_combErr, ele_PFcombErr;
      float ele_ecalRegressionEnergy,ele_ecalRegressionError, ele_ecalTrackRegressionEnergy,ele_ecalTrackRegressionError,ele_ecalScale,
	ele_ecalSmear,ele_ecalRegressionScale,ele_ecalRegressionSmear,
	ele_ecalTrackRegressionScale,ele_ecalTrackRegressionSmear;
      
      float ele_sclphiwidth, ele_scletawidth;
      float ele_IoEmIop;


      float pho_ID1_value;
      float pho_ID2_value;

      //
      float ele_psE;
      float ele_psEoverEraw;  
      float ele_kfchi2;
      int ele_kfhits;
      int ele_gsfhits;

      bool pho_hasPixelSeed, pho_passElectronVeto;
      float pho_phoPhotonIsolation_CITK;
      float pho_phoNeutralHadronIsolation_CITK;
      float pho_phoChargedIsolation_CITK;

      //float pho_etOutsideMustache;
      //``int pho_nClusterOutsideMustache;
      int pho_mipNhitCone;
      int pho_conversionTrackProvenance;
      int pho_conversion_size;
      int pho_one_leg_conversion_size;

      bool pho_hasConversionTracks;
      int mc_pho_photon_parent_PDGID;
      int mc_pho_photon_parent_status;
      


      	//MC
    int _MC_TrueNumInteractions;
};
#endif
