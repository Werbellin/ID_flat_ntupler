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

enum class inputFileTypes {AOD, MINIAOD, UNDEF};

class Ntuplizer : public edm::EDAnalyzer {
   public:
      explicit Ntuplizer(const edm::ParameterSet&);
      ~Ntuplizer();

      typedef math::XYZTLorentzVector LorentzVector ;
      
      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

      void Init();
      void FillEvent(const edm::Event&, const edm::EventSetup&);
      void FillElectrons(const edm::Event&, const edm::EventSetup&);
      void FillVertices(const edm::Event&, const edm::EventSetup&);
      void FillTruth(const edm::Event&, const edm::EventSetup&);
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
      
      edm::EDGetToken electronsToken_;
      edm::EDGetTokenT<edm::View<reco::GenParticle> > genParticlesToken_;
      edm::EDGetTokenT<reco::VertexCollection> vertexToken_;
      edm::EDGetTokenT<reco::ConversionCollection> conversionsToken_;
      edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
      edm::EDGetToken pfMETToken_;
      edm::EDGetTokenT<vector<reco::GenParticle> > genParticleToken_;
      edm::EDGetTokenT<edm::View<reco::GenParticle> > genParticlesToken_CB;
      edm::EDGetTokenT<GenEventInfoProduct> genEventInfoProductTagToken_;
      edm::EDGetTokenT<ValueMap<float>> electronEcalPFClusterIsolationProducerToken_;
      edm::EDGetTokenT<ValueMap<float>> electronHcalPFClusterIsolationProducerToken_;

      edm::EDGetTokenT<ValueMap<float>> electronID1Token_;
      edm::EDGetTokenT<ValueMap<float>> electronID2Token_;

      // Trigger Stuff
      edm::InputTag HLTTag_; 
      bool isMC_;	

      vector<trigger::TriggerObject> _selectedObjects;
      vector<trigger::TriggerObject> _hltEle27WP75GsfTrackIsoFilter;
      
      edm::InputTag PileupSrc_;	
      
      //tree
      TTree *_mytree;
      TLorentzVector myvector ;  

      //global variables
      int _nEvent, _nRun, _nLumi;
      //pile-up
      int _PU_N;

      double _mc_event_weight;

      typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> NewLorentzVector;
      vector<NewLorentzVector> _mc_gen_ele_p4;

      //vertices
      int _vtx_N;
      
      //trigger fired names
      char trig_fired_names[10000];

      // MET
      double _met_pf_et,_met_pf_px, _met_pf_py, _met_pf_phi, _met_pf_set, _met_pf_sig; 

      //electrons
      int ele_N;
      int ele_N_saved;
      TClonesArray * m_electrons;


      vector<float> ele_conversionVertexFitProbability;
      vector<int>  mc_ele_isPromptFinalState;
      vector<int>  mc_ele_isDirectPromptTauDecayProductFinalState;
      vector<int> mc_ele_matchedFromCB;
      vector<int> mc_ele_matchMother_PDGID;

      vector<float> ele_dr03EcalRecHitSumEt;
      vector<float> ele_dr03HcalTowerSumEt;
      vector<float> ele_dr03TkSumPt;
      vector<float> ele_pt;
      vector<float> ele_trackMomentumAtVtx_R;
      vector<float> ele_electronEcalPFClusterIsolationProducer;
      vector<float> ele_electronHcalPFClusterIsolationProducer;
      vector<float> ele_ID1;
      vector<float> ele_ID2;

      vector<int> ele_index;

      vector<string> event_trig_fired;
      vector<bool> ele_trig_passed_filter;
      vector<bool> ele_pass_hltEle27WP75GsfTrackIsoFilter;
      vector<float> ele_full5x5_hcalOverEcal;

      int ele_echarge[50];
      double ele_he[50], ele_hebc[50], ele_eseedpout[50] , ele_ep[50] , ele_eseedp[50] , ele_eelepout[50] ;       
      double ele_deltaetaseed[50] , ele_deltaetaele[50] , ele_deltaphiseed[50] , ele_deltaphiele[50] , ele_deltaetain[50] , ele_deltaphiin[50] ;
      double ele_sigmaietaieta[50] , ele_sigmaetaeta[50] , ele_sigmaiphiiphi[50], ele_e15[50] , ele_e25max[50] , ele_e55[50] , ele_e1[50], ele_r9[50] ;

      double ele_oldsigmaetaeta[50], ele_oldsigmaietaieta[50], ele_oldsigmaiphiiphi[50], ele_oldsigmaietaiphi[50], ele_oldr9[50], ele_olde15[50],           
	ele_olde25max[50], ele_olde55[50];
      double ele_oldhe[50], ele_oldhebc[50];


      double ele_pin_mode[50] , ele_pout_mode[50] , ele_pTin_mode[50] , ele_pTout_mode[50] ; 
      //
      double ele_fbrem[50],ele_SCfbrem[50],ele_pfSCfbrem[50], ele_trackfbrem[50];
      int ele_nbrem[50];
      //
      double ele_mva[50] ;
      //
      int ele_isbarrel[50] , ele_isendcap[50] , 
	ele_isEBetaGap[50] , ele_isEBphiGap[50] , ele_isEEdeeGap[50] , ele_isEEringGap[50] ,
	ele_isecalDriven[50] , ele_istrackerDriven[50] ,
	ele_eClass[50];
      int ele_valid_hits[50], ele_lost_hits[50]; 
      double ele_gsfchi2[50]; 
      double ele_dxyB[50], ele_dxy[50], ele_dzB[50], ele_dz[50], ele_dszB[50], ele_dsz[50];              
      double ele_tkSumPt_dr03[50] , ele_ecalRecHitSumEt_dr03[50] , ele_hcalDepth1TowerSumEt_dr03[50] , ele_hcalDepth2TowerSumEt_dr03[50] ,
	ele_tkSumPt_dr04[50] , ele_ecalRecHitSumEt_dr04[50] , ele_hcalDepth1TowerSumEt_dr04[50] , ele_hcalDepth2TowerSumEt_dr04[50] ;
      //
      double ele_conv_dcot[50];
      double ele_conv_dist[50];
      double ele_conv_radius[50];
      int ele_expected_inner_hits[50];
      int ele_vtxconv[50];
      //
      double ele_eidVeryLoose[50], ele_eidLoose[50], ele_eidMedium[50], ele_eidTight[50]; 
      double ele_eidHZZVeryLoose[50], ele_eidHZZLoose[50], ele_eidHZZMedium[50], ele_eidHZZTight[50], ele_eidHZZSuperTight[50]; 
      double ele_eidMVATrig[50], ele_eidMVANoTrig[50];
      double ele_HZZisoTk[50],ele_HZZisoTk5[50],  ele_HZZisoEcal[50], ele_HZZisoHcal[50],ele_HZZisoComb[50] ;
      //
      double ele_pfChargedHadIso[50], ele_pfNeutralHadIso[50], ele_pfPhotonIso[50], ele_pfChargedIso[50], ele_pfSumPUIso[50];
      //
      double ele_dzPV[50], ele_d0[50], ele_d0err[50];
      double ele_IP[50], ele_IPError[50], ele_SIP[50] ;
      double ele_sclE[50], ele_sclEt[50], ele_sclEta[50], ele_sclPhi[50], ele_sclRawE[50];
      int ele_sclNclus[50];
      
      double ele_sclsubE[50][20], ele_sclsubEta[50][20], ele_sclsubPhi[50][20];
      int ele_sclsubisseed[50][20];
    
      double ele_ecalE[50], ele_ecalErr[50], ele_trackErr[50], ele_combErr[50], ele_PFcombErr[50];
      double ele_ecalRegressionEnergy[50],ele_ecalRegressionError[50], ele_ecalTrackRegressionEnergy[50],ele_ecalTrackRegressionError[50],ele_ecalScale[50],
	ele_ecalSmear[50],ele_ecalRegressionScale[50],ele_ecalRegressionSmear[50],
	ele_ecalTrackRegressionScale[50],ele_ecalTrackRegressionSmear[50];
      
      double ele_HCALFullConeSum[50];
      double ele_sclphiwidth[50], ele_scletawidth[50];
      //
      double ele_psE[50];
      
      double ele_kfchi2[50];
      int ele_kfhits[50];
      int ele_gsfhits[50];

      	//MC
    int _MC_TrueNumInteractions;
};
#endif
