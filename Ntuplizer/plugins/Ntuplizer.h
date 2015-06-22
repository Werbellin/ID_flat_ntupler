#ifndef Ntuplizer_H
#define Ntuplizer_H

// Own
//#include "fBremTest/DemoAnalyzer/interface/TrajectoryAnalysisHelper.h"
//#include "fBremTest/DemoAnalyzer/interface/RecSimHitMatcher.h"


// CMSSW
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include <DataFormats/HepMCCandidate/interface/GenParticleFwd.h>
#include "DataFormats/Math/interface/LorentzVector.h"
#include <FWCore/Framework/interface/ESHandle.h>

// ROOT
#include "TTree.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TParticle.h"
#include "TVector3.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"

// C++
#include<memory>
#include<vector>

using namespace std;
using std::vector;

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
      
      //void setMomentum(TLorentzVector & myvector, const LorentzVector & mom);
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
      //inputTag
      edm::InputTag EleTag_;
      //edm::InputTag MuonTag_;
      //edm::InputTag JetTag_;
      //edm::InputTag PhotonTag_;
      edm::InputTag VerticesTag_;
      // Trigger Stuff
      edm::InputTag HLTTag_; 
      //edm::InputTag triggerEventTag_;
      //edm::InputTag MCTag_ ;
      bool isMC_;	
      bool ispythia6_;
      //int lepton_setup;
      
      //edm::InputTag MuRhoCorrection_;
      //edm::InputTag EleRhoCorrection_;
      //edm::InputTag SigmaRhoCorrection_;
      edm::InputTag PileupSrc_;	
      
      std::vector<edm::InputTag> MVAidCollection_;
	
      //std::vector<edm::InputTag > HLT_Filters_;
      //edm::InputTag SCTag_;
	
      //tree
      TTree *_mytree;
      TLorentzVector myvector ;  

      //global variables
      int _nEvent, _nRun, _nLumi;
      //pile-up
      int _PU_N;

      //vertices
      int _vtx_N;
      
      //trigger fired names
      char trig_fired_names[10000];

      // MET
      double _met_pf_et,_met_pf_px, _met_pf_py, _met_pf_phi, _met_pf_set, _met_pf_sig; 

      //electrons
      int ele_N;
      TClonesArray * m_electrons;


      //vector < vector<int> > trackHitPDGID;
      //vector <float> lastHitPt;


      //vector<bool>  ele_foundGSFTraj;
      //vector<float> ele_signedEstimateSumPred;
      //float ele_signedEstimateSumPred_A[50];
      //vector<float> ele_propagatorSignedEstimateSumPred;
      //vector<float> ele_signSumPredNormVH;

      //vector<bool>  ele_foundCKFTraj;
      //vector<float> ele_signedEstimateSumPredCKF;
      //vector<float> ele_reducedChi2CKF;
      //vector<float> ele_conversionVertexFitProbability;
    



      int ele_echarge[50];
      double ele_he[50], ele_hebc[50], ele_eseedpout[50] , ele_ep[50] , ele_eseedp[50] , ele_eelepout[50] ;       
      double ele_deltaetaseed[50] , ele_deltaetaele[50] , ele_deltaphiseed[50] , ele_deltaphiele[50] , ele_deltaetain[50] , ele_deltaphiin[50] ;
      double ele_sigmaietaieta[50] , ele_sigmaetaeta[50] , ele_sigmaiphiiphi[50], ele_e15[50] , ele_e25max[50] , ele_e55[50] , ele_e1[50], ele_r9[50] ;

      double ele_oldsigmaetaeta[50], ele_oldsigmaietaieta[50], ele_oldsigmaiphiiphi[50], ele_oldsigmaietaiphi[50], ele_oldr9[50], ele_olde15[50],           
	ele_olde25max[50], ele_olde55[50];
      double ele_oldhe[50], ele_oldhebc[50];


      //, ele_e33[50] , ele_e2overe9[50] ;
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
      //ele_pfChargedHadPUIso[50],ele_pfCombRelIso[50] ;
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

      double ele_mvaphys14[50];
      double ele_mvaphys14fix[50];
      /* double ele_mvafbrem[50], ele_mvadetain[50], ele_mvadphiin[50], ele_mvasieie[50], ele_mvahoe[50], ele_mvaeop[50],  */
      /* 	ele_mvae1x5e5x5[50], ele_mvaeleopout[50], ele_mvakfchi2[50], ele_mvadist[50],  ele_mvadcot[50], ele_mvaeta[50], */
      /* 	ele_mvapt[50]; */
      /*       int ele_mvakfhits[50], ele_mvamishits[50], ele_mvaecalseed	[50]; */
      
      // MC
      	//MC
	TClonesArray * _m_MC_gen_V;
	TClonesArray * _m_MC_gen_Higgs;
	//TClonesArray * _m_MC_gen_photons;
	TClonesArray * _m_MC_gen_leptons;
	TClonesArray * _m_MC_gen_leptons_status1;
	TClonesArray * _m_MC_gen_leptons_status2;
	double _MC_gen_V_pdgid[10];
	double _MC_gen_leptons_pdgid[30];
	double _MC_gen_leptons_status1_pdgid[30];
	double _MC_gen_leptons_status2_pdgid[30];

	int _MC_gen_leptons_status1_FromWZ[30];
	int _MC_gen_leptons_status1_FromTaus[30];
	int _MC_gen_leptons_status1_FromNonPrompt[30];
	int _MC_gen_leptons_status1_FromBC[30];

	//double _MC_pthat;
	int _MC_flavor[2];

	//int _MC_gen_photons_isFSR[5000];
	
      bool runGsfRefitter;
      std::string GSFTrajColl;

      bool runKfWithGsfRefitter;
      std::string CKFTrajColl;

};
#endif
