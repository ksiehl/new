// -*- C++ -*-
//
// Package:    Analysis/B2GTTbarTreeMaker
// Class:      B2GTTbarTreeMaker
// 
/**\class B2GTTbarTreeMaker B2GTTbarTreeMaker.cc Analysis/B2GTTbarTreeMaker/plugins/B2GTTbarTreeMaker.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  James Dolen
//         Created:  Sat, 30 Apr 2016 17:40:42 GMT
//
//

//--------------------------
// To add:
// - NNPDF3weight  
// - electron quality cuts
// - trigger
//
// Note. The following items should be applied at the tree reader level:
// - MU ID (HIP)
// -- TFile* f_muID = TFile::Open("MuonID_Z_RunBCD_prompt80X_7p65.root", "read");
// -- TH1F* h_muID = (TH1F*) f_muID->Get("MC_NUM_MediumID_DEN_genTracks_PAR_eta/eta_ratio")->Clone();
// -- float SF_muID = h_muID->GetBinContent(h_muID->FindBin(eta););
// - BTagCalibrationReader
//--------------------------


// system include files
#include <memory>
#include <iostream>    
#include <algorithm>   
#include <bitset>   

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

// DataFormats
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/MET.h"

// TFileService
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// Gen particle
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

// JEC
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

// JER
#include "JetMETCorrections/Modules/interface/JetResolution.h"

// Electron
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

// Trigger
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

// Vertex
#include "DataFormats/VertexReco/interface/Vertex.h"

// Pileup
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

// LHE weights
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

// Utilities
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"

// root
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TLorentzVector.h"

//for files
#include <fstream>

#include "Analysis/B2GTTbar/metZcalc/METzCalculator.cc"
//std::ofstream outfile;

Int_t fakes = 0;

Float_t costheta1 = 0.0;
Float_t costheta2 = 0.0;
Float_t phi = 0.0;
Float_t costhetastar = 0.0;
Float_t phistar1 = 0.0;
Float_t phistar2 = 0.0;

const Float_t Wmass = 80.2;

Float_t FourBodyMass(TLorentzVector, TLorentzVector);

void calculateAngles(TLorentzVector thep4M11, TLorentzVector thep4M12, TLorentzVector thep4M21, TLorentzVector thep4M22, Float_t& costheta1, Float_t& costheta2, Float_t& phi, Float_t& costhetastar, Float_t& phistar1, Float_t& phistar2);

//RS gluon PDF weights
namespace LHAPDF
{
  void initPDFSet(int nset, const std::string& filename, int member = 0);
  int numberPDF(int nset);
  void usePDFMember(int nset, int member);
  double xfx(int nset, double x, double Q, int fl);
  double getXmin(int nset, int member);
  double getXmax(int nset, int member);
  double getQ2min(int nset, int member);
  double getQ2max(int nset, int member);
  void extrapolate(bool extrapolate = true);
}



//
// class declaration
//

class B2GTTbarTreeMaker : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
public:
  explicit B2GTTbarTreeMaker(const edm::ParameterSet&);
  ~B2GTTbarTreeMaker();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  std::ofstream outfile;
  

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  // ----------member data ---------------------------
  edm::EDGetTokenT<pat::JetCollection> ak4jetToken_;
  edm::EDGetTokenT<pat::JetCollection> ak8jetToken_;
  edm::EDGetTokenT<pat::JetCollection> puppijetToken_;
  edm::EDGetTokenT<pat::JetCollection> ak8CHSSoftDropSubjetsToken_;
  edm::EDGetTokenT<pat::JetCollection> ak8PuppiSoftDropSubjetsToken_;
  edm::EDGetTokenT<reco::GenJetCollection> ak4genjetToken_;
  edm::EDGetTokenT<reco::GenJetCollection> ak8genjetToken_;
  edm::EDGetTokenT<edm::View<reco::GenParticle> > prunedGenToken_;
  edm::EDGetTokenT<double> rhoToken_;  
  edm::EDGetTokenT<std::vector<reco::Vertex> > vtxToken_;
  edm::EDGetTokenT<edm::TriggerResults> triggerResultsMETFilterToken_;
  edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
  edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;
  edm::EDGetTokenT<bool> badMuonFilterToken_;
  edm::EDGetTokenT<bool> badChargedCandidateFilterToken_;
  edm::EDGetTokenT<pat::MuonCollection> muonToken_;
  edm::EDGetTokenT<pat::ElectronCollection> electronToken_;
  edm::EDGetTokenT<pat::METCollection> metToken_;
  edm::EDGetTokenT<std::vector<PileupSummaryInfo> > pileupInfoToken_;
  edm::EDGetTokenT<LHEEventProduct> theSrc_;
  edm::EDGetTokenT<GenEventInfoProduct> pdfToken_;


  
  bool useToolbox_;
  bool verbose_;
  bool verboseGen_;
  bool runGenLoop_;
  bool isZprime_;
  bool isttbar_;
  bool isRSG_;
  std::vector<std::string>  jecPayloadsAK4chs_;
  std::vector<std::string>  jecPayloadsAK8chs_;
  std::vector<std::string>  jecPayloadsAK4pup_;
  std::vector<std::string>  jecPayloadsAK8pup_;
  boost::shared_ptr<FactorizedJetCorrector> JetCorrectorAK4chs;
  boost::shared_ptr<FactorizedJetCorrector> JetCorrectorAK8chs;
  boost::shared_ptr<FactorizedJetCorrector> JetCorrectorAK4pup;
  boost::shared_ptr<FactorizedJetCorrector> JetCorrectorAK8pup;
  boost::shared_ptr<JetCorrectionUncertainty> JetCorrUncertAK4chs;
  boost::shared_ptr<JetCorrectionUncertainty> JetCorrUncertAK8chs;
  boost::shared_ptr<JetCorrectionUncertainty> JetCorrUncertAK4pup;
  boost::shared_ptr<JetCorrectionUncertainty> JetCorrUncertAK8pup;
      
  std::string jerSFtext_;

  //MY ADDITION
  /*edm::EDGetTokenT<double> mgreweight1_;
  edm::EDGetTokenT<double> mgreweight2_;
  edm::EDGetTokenT<double> mgreweight3_;
  edm::EDGetTokenT<double> mgreweight4_;
  edm::EDGetTokenT<double> mgreweight5_;
  edm::EDGetTokenT<double> mgreweight6_;
  edm::EDGetTokenT<double> mgreweight7_;
  edm::EDGetTokenT<double> mgreweight8_;
  edm::EDGetTokenT<double> mgreweight9_;
  edm::EDGetTokenT<double> mgreweight10_;
  edm::EDGetTokenT<double> mgreweight11_;
  edm::EDGetTokenT<double> mgreweight12_;
  edm::EDGetTokenT<double> mgreweight13_;
  edm::EDGetTokenT<double> mgreweight14_;
  edm::EDGetTokenT<double> mgreweight15_;
  edm::EDGetTokenT<double> mgreweight16_;
  edm::EDGetTokenT<double> mgreweight17_;
  edm::EDGetTokenT<double> mgreweight18_;
  edm::EDGetTokenT<double> mgreweight19_;
  edm::EDGetTokenT<double> mgreweight20_;
  edm::EDGetTokenT<double> mgreweight21_;
  edm::EDGetTokenT<double> mgreweight22_;
  edm::EDGetTokenT<double> mgreweight23_;
  edm::EDGetTokenT<double> mgreweight24_;
  edm::EDGetTokenT<double> mgreweight25_;
  edm::EDGetTokenT<double> mgreweight26_;
  edm::EDGetTokenT<double> mgreweight27_;
  edm::EDGetTokenT<double> mgreweight28_;
  edm::EDGetTokenT<double> mgreweight29_;
  edm::EDGetTokenT<double> mgreweight30_;
  edm::EDGetTokenT<double> mgreweight31_;
  edm::EDGetTokenT<double> mgreweight32_;
  edm::EDGetTokenT<double> mgreweight33_;
  edm::EDGetTokenT<double> mgreweight34_;
  edm::EDGetTokenT<double> mgreweight35_;
  edm::EDGetTokenT<double> mgreweight36_;
  edm::EDGetTokenT<double> mgreweight37_;
  edm::EDGetTokenT<double> mgreweight38_;
  edm::EDGetTokenT<double> mgreweight39_;
  edm::EDGetTokenT<double> mgreweight40_;
  edm::EDGetTokenT<double> mgreweight41_;
  edm::EDGetTokenT<double> mgreweight42_;
  edm::EDGetTokenT<double> mgreweight43_;
  edm::EDGetTokenT<double> mgreweight44_;
  edm::EDGetTokenT<double> mgreweight45_;
  edm::EDGetTokenT<double> mgreweight46_;
  edm::EDGetTokenT<double> mgreweight47_;
  edm::EDGetTokenT<double> mgreweight48_;
  edm::EDGetTokenT<double> mgreweight49_;
  edm::EDGetTokenT<double> mgreweight50_;*/

  //END MY ADDITION

  TFile* fPUweight;
  TH1D* hPUweight;
  TH1D* hPUweight_MBup;
  TH1D* hPUweight_MBdn;

  int count_GenTruth_semileptonic;
  int count_nMu_gt1; 
  int count_nEl_gt1; 
  int count_nMu_e1; 
  int count_nEl_e1; 
  int count_nLep_e1; 
  int count_JetPt300; 
  int count_JetPt300Eta; 
  int count_JetPt300Eta_AK4; 
  int count_JetPt300Eta_muPt40; 
  int count_JetPt300Eta_muPt40_MET40; 
  int count_JetPt300Eta_muPt40_MET40_AK4; 


  TH1D * h_ak8puppi_softDropMass;
  TH1D * h_ak8chs_softDropMass;
  TH1D * h_ak8chs_softDropMass_reweighted;
  TH1D * h_ak8chs_pt;
  TH1D * h_ak8chs_pt_reweighted;
  TH1D * h_NtrueIntPU;
  TH1D * h_NPV;               
  TH1D * h_NPVreweighted;                 

  //
  //       d8888 888 888        888    888               888     88888888888                           
  //      d88888 888 888        888    888               888         888                               
  //     d88P888 888 888        888    888               888         888                               
  //    d88P 888 888 888        8888888888  8888b.   .d88888         888     888d888  .d88b.   .d88b.  
  //   d88P  888 888 888        888    888     "88b d88" 888         888     888P"   d8P  Y8b d8P  Y8b 
  //  d88P   888 888 888 888888 888    888 .d888888 888  888         888     888     88888888 88888888 
  // d8888888888 888 888        888    888 888  888 Y88b 888         888     888     Y8b.     Y8b.     
  //d88P     888 888 888        888    888 "Y888888  "Y88888         888     888      "Y8888   "Y8888  
  //                                                                                                   
    
         

  TTree *TreeAllHad;   

  std::vector<std::string> *AllHadTrigNames     = new std::vector<std::string>;
  std::vector<int> *AllHadTrigPrescales = new std::vector<int>;
  std::vector<bool> *AllHadTrigPass    = new std::vector<bool>;

  std::string AllHadTrigAcceptBits;

  Int_t   PassMETFilters;
  Float_t Jet0PtRaw;
  Float_t Jet0EtaRaw;
  Float_t Jet0PhiRaw;
  Float_t Jet0MassRaw;
  Float_t Jet0P;
  Float_t Jet0Pt;
  Float_t Jet0Eta;
  Float_t Jet0Phi;
  Float_t Jet0Rap;
  Float_t Jet0Energy;
  Float_t Jet0Mass;
  Float_t Jet0Area;
  Float_t Jet0SDmass;
  Float_t Jet0SDmassRaw;
  Float_t Jet0SDmassCorrL23;    
  Float_t Jet0SDmassCorrL23Up;    
  Float_t Jet0SDmassCorrL23Dn;     
  Float_t Jet0SDmassCorrL23Smear;    
  Float_t Jet0SDmassCorrL23SmearUp;    
  Float_t Jet0SDmassCorrL23SmearDn;  
  Float_t Jet0SDmassCorrL123;    
  Float_t Jet0SDmassCorrL123Up;    
  Float_t Jet0SDmassCorrL123Dn;     
  Float_t Jet0SDptRaw;    
  Float_t Jet0SDptCorrL23;    
  Float_t Jet0SDptCorrL23Up;    
  Float_t Jet0SDptCorrL23Dn;    
  Float_t Jet0SDptCorrL123;    
  Float_t Jet0SDptCorrL123Up;    
  Float_t Jet0SDptCorrL123Dn;    
  Float_t Jet0SDptCorrL23Smear;    
  Float_t Jet0SDptCorrL23SmearUp;    
  Float_t Jet0SDptCorrL23SmearDn;    
  Float_t Jet0SDetaRaw;    
  Float_t Jet0SDphiRaw;    
  Float_t Jet0MassPruned;
  Float_t Jet0MassTrimmed;
  Float_t Jet0Tau1;
  Float_t Jet0Tau2;
  Float_t Jet0Tau3;
  Float_t Jet0Tau4;
  Float_t Jet0Tau32;
  Float_t Jet0Tau21;
  Float_t Jet0SDsubjet0bdisc;
  Float_t Jet0SDsubjet1bdisc;
  Float_t Jet0SDmaxbdisc;
  Float_t Jet0SDmaxbdiscflavHadron;
  Float_t Jet0SDmaxbdiscflavParton;
  Float_t Jet0SDsubjet0pt;
  Float_t Jet0SDsubjet0mass;
  Float_t Jet0SDsubjet0eta;
  Float_t Jet0SDsubjet0phi;
  Float_t Jet0SDsubjet0area;
  Float_t Jet0SDsubjet0flavHadron;
  Float_t Jet0SDsubjet0flavParton;
  Float_t Jet0SDsubjet0tau1;
  Float_t Jet0SDsubjet0tau2;
  Float_t Jet0SDsubjet0tau3;
  Float_t Jet0SDsubjet1pt;
  Float_t Jet0SDsubjet1mass;
  Float_t Jet0SDsubjet1eta;
  Float_t Jet0SDsubjet1phi;
  Float_t Jet0SDsubjet1area;
  Float_t Jet0SDsubjet1flavHadron;
  Float_t Jet0SDsubjet1flavParton;
  Float_t Jet0SDsubjet1tau1;
  Float_t Jet0SDsubjet1tau2;
  Float_t Jet0SDsubjet1tau3;
  Float_t Jet0PuppiP;
  Float_t Jet0PuppiPt;
  Float_t Jet0PuppiEta;
  Float_t Jet0PuppiPhi;
  Float_t Jet0PuppiMass;

  Float_t Jet0PuppiSDmass;
  Float_t Jet0PuppiSDmassCorr;
  Float_t Jet0PuppiSDmassCorrUp;
  Float_t Jet0PuppiSDmassCorrDn;
  Float_t Jet0PuppiSDmassCorrL23Smear;
  Float_t Jet0PuppiSDmassCorrL23SmearUp;
  Float_t Jet0PuppiSDmassCorrL23SmearDn;
  Float_t Jet0PuppiSDpt;
  Float_t Jet0PuppiSDptCorr;
  Float_t Jet0PuppiSDptCorrUp;
  Float_t Jet0PuppiSDptCorrDn;
  Float_t Jet0PuppiSDptCorrL23Smear;
  Float_t Jet0PuppiSDptCorrL23SmearUp;
  Float_t Jet0PuppiSDptCorrL23SmearDn;
  Float_t Jet0PuppiSDeta;
  Float_t Jet0PuppiSDphi;

  Float_t Jet0PuppiTau1;
  Float_t Jet0PuppiTau2;
  Float_t Jet0PuppiTau3;
  Float_t Jet0PuppiTau4;
  Float_t Jet0PuppiTau32;
  Float_t Jet0PuppiTau21;
  Float_t Jet0PuppiSDsubjet0bdisc;
  Float_t Jet0PuppiSDsubjet1bdisc;
  Float_t Jet0PuppiSDmaxbdisc;
  Float_t Jet0PuppiSDmaxbdiscflavHadron;
  Float_t Jet0PuppiSDmaxbdiscflavParton;
  Float_t Jet0PuppiSDsubjet0pt;
  Float_t Jet0PuppiSDsubjet0mass;
  Float_t Jet0PuppiSDsubjet0eta;
  Float_t Jet0PuppiSDsubjet0phi;
  Float_t Jet0PuppiSDsubjet0area;
  Float_t Jet0PuppiSDsubjet0flavHadron;
  Float_t Jet0PuppiSDsubjet0flavParton;
  Float_t Jet0PuppiSDsubjet0tau1;
  Float_t Jet0PuppiSDsubjet0tau2;
  Float_t Jet0PuppiSDsubjet0tau3;
  Float_t Jet0PuppiSDsubjet1pt;
  Float_t Jet0PuppiSDsubjet1mass;
  Float_t Jet0PuppiSDsubjet1eta;
  Float_t Jet0PuppiSDsubjet1phi;
  Float_t Jet0PuppiSDsubjet1area;
  Float_t Jet0PuppiSDsubjet1flavHadron;
  Float_t Jet0PuppiSDsubjet1flavParton;
  Float_t Jet0PuppiSDsubjet1tau1;
  Float_t Jet0PuppiSDsubjet1tau2;
  Float_t Jet0PuppiSDsubjet1tau3;
  Float_t Jet0CHF;
  Float_t Jet0NHF;
  Float_t Jet0CM;
  Float_t Jet0NM;
  Float_t Jet0NEF;
  Float_t Jet0CEF;
  Float_t Jet0MF;
  Float_t Jet0Mult;
  Float_t Jet0PuppiCHF;
  Float_t Jet0PuppiNHF;
  Float_t Jet0PuppiCM;
  Float_t Jet0PuppiNM;
  Float_t Jet0PuppiNEF;
  Float_t Jet0PuppiCEF;
  Float_t Jet0PuppiMF;
  Float_t Jet0PuppiMult;
  Float_t Jet0MassCorrFactor;
  Float_t Jet0MassCorrFactorUp;
  Float_t Jet0MassCorrFactorDn;
  Float_t Jet0CorrFactor;
  Float_t Jet0CorrFactorUp;
  Float_t Jet0CorrFactorDn;
  Float_t Jet0PtSmearFactor;
  Float_t Jet0PtSmearFactorUp;
  Float_t Jet0PtSmearFactorDn;
  Float_t Jet0PuppiMassCorrFactor;
  Float_t Jet0PuppiMassCorrFactorUp;
  Float_t Jet0PuppiMassCorrFactorDn;
  Float_t Jet0PuppiCorrFactor;
  Float_t Jet0PuppiCorrFactorUp;
  Float_t Jet0PuppiCorrFactorDn;
  Float_t Jet0PuppiPtSmearFactor;
  Float_t Jet0PuppiPtSmearFactorUp;
  Float_t Jet0PuppiPtSmearFactorDn;
  Float_t Jet0EtaScaleFactor;
  Float_t Jet0PhiScaleFactor;
  // Float_t Jet0MatchedGenJetDR;
  Float_t Jet0MatchedGenJetPt;
  Float_t Jet0MatchedGenJetMass;

  Int_t   Jet0GenMatched_TopHadronic;
  Float_t Jet0GenMatched_TopPt;
  Float_t Jet0GenMatched_TopEta;
  Float_t Jet0GenMatched_TopPhi;
  Float_t Jet0GenMatched_TopMass;
  Float_t Jet0GenMatched_bPt;
  Float_t Jet0GenMatched_WPt;
  Float_t Jet0GenMatched_Wd1Pt;
  Float_t Jet0GenMatched_Wd2Pt;
  Float_t Jet0GenMatched_Wd1ID;
  Float_t Jet0GenMatched_Wd2ID;
  Float_t Jet0GenMatched_MaxDeltaRPartonTop;
  Float_t Jet0GenMatched_MaxDeltaRWPartonTop;
  Float_t Jet0GenMatched_MaxDeltaRWPartonW;
  Float_t Jet0GenMatched_DeltaR_t_b;
  Float_t Jet0GenMatched_DeltaR_t_W;
  Float_t Jet0GenMatched_DeltaR_t_Wd1;
  Float_t Jet0GenMatched_DeltaR_t_Wd2;
  Float_t Jet0GenMatched_DeltaR_W_b1;
  Float_t Jet0GenMatched_DeltaR_W_Wd1;
  Float_t Jet0GenMatched_DeltaR_W_Wd2;
  Float_t Jet0GenMatched_DeltaR_Wd1_Wd2;
  Float_t Jet0GenMatched_DeltaR_Wd1_b;
  Float_t Jet0GenMatched_DeltaR_Wd2_b;
  Float_t Jet0GenMatched_DeltaR_jet_t;
  Float_t Jet0GenMatched_DeltaR_jet_W;
  Float_t Jet0GenMatched_DeltaR_jet_b;
  Float_t Jet0GenMatched_DeltaR_jet_Wd1;
  Float_t Jet0GenMatched_DeltaR_jet_Wd2;
  Float_t Jet0GenMatched_DeltaR_pup0_b;
  Float_t Jet0GenMatched_DeltaR_pup0_Wd1;
  Float_t Jet0GenMatched_DeltaR_pup0_Wd2;
  Float_t Jet0GenMatched_DeltaR_pup1_b;
  Float_t Jet0GenMatched_DeltaR_pup1_Wd1;
  Float_t Jet0GenMatched_DeltaR_pup1_Wd2;
  Float_t Jet0GenMatched_partonPt;
  Float_t Jet0GenMatched_partonEta;
  Float_t Jet0GenMatched_partonPhi;
  Float_t Jet0GenMatched_partonMass;
  Float_t Jet0GenMatched_partonID;
  Float_t Jet0GenMatched_DeltaRjetParton;

  Float_t Jet1PtRaw;
  Float_t Jet1EtaRaw;
  Float_t Jet1PhiRaw;
  Float_t Jet1MassRaw;
  Float_t Jet1P;
  Float_t Jet1Pt;
  Float_t Jet1Eta;
  Float_t Jet1Phi;
  Float_t Jet1Rap;
  Float_t Jet1Energy;
  Float_t Jet1Mass;
  Float_t Jet1Area;
  Float_t Jet1SDmass;
  Float_t Jet1SDmassRaw;
  Float_t Jet1SDmassCorrL23;
  Float_t Jet1SDmassCorrL23Up;
  Float_t Jet1SDmassCorrL23Dn;
  Float_t Jet1SDmassCorrL123;
  Float_t Jet1SDmassCorrL123Up;
  Float_t Jet1SDmassCorrL123Dn;
  Float_t Jet1SDmassCorrL23Smear;
  Float_t Jet1SDmassCorrL23SmearUp;
  Float_t Jet1SDmassCorrL23SmearDn;
  Float_t Jet1SDptRaw;
  Float_t Jet1SDptCorrL23;
  Float_t Jet1SDptCorrL23Up;
  Float_t Jet1SDptCorrL23Dn;
  Float_t Jet1SDptCorrL123;
  Float_t Jet1SDptCorrL123Up;
  Float_t Jet1SDptCorrL123Dn;
  Float_t Jet1SDptCorrL23Smear;
  Float_t Jet1SDptCorrL23SmearUp;
  Float_t Jet1SDptCorrL23SmearDn;
  Float_t Jet1SDetaRaw;
  Float_t Jet1SDphiRaw; 
  Float_t Jet1MassPruned;
  Float_t Jet1MassTrimmed;
  Float_t Jet1Tau1;
  Float_t Jet1Tau2;
  Float_t Jet1Tau3;
  Float_t Jet1Tau4;
  Float_t Jet1Tau32;
  Float_t Jet1Tau21;
  Float_t Jet1SDsubjet0bdisc;
  Float_t Jet1SDsubjet1bdisc;
  Float_t Jet1SDmaxbdisc;
  Float_t Jet1SDmaxbdiscflavHadron;
  Float_t Jet1SDmaxbdiscflavParton;
  Float_t Jet1SDsubjet0pt;
  Float_t Jet1SDsubjet0eta;
  Float_t Jet1SDsubjet0phi;
  Float_t Jet1SDsubjet0mass;
  Float_t Jet1SDsubjet0area;
  Float_t Jet1SDsubjet0flavHadron;
  Float_t Jet1SDsubjet0flavParton;
  Float_t Jet1SDsubjet0tau1;
  Float_t Jet1SDsubjet0tau2;
  Float_t Jet1SDsubjet0tau3;
  Float_t Jet1SDsubjet1pt;
  Float_t Jet1SDsubjet1eta;
  Float_t Jet1SDsubjet1phi;
  Float_t Jet1SDsubjet1mass;
  Float_t Jet1SDsubjet1area;
  Float_t Jet1SDsubjet1flavHadron;
  Float_t Jet1SDsubjet1flavParton;
  Float_t Jet1SDsubjet1tau1;
  Float_t Jet1SDsubjet1tau2;
  Float_t Jet1SDsubjet1tau3;
  Float_t Jet1PuppiP;
  Float_t Jet1PuppiPt;
  Float_t Jet1PuppiEta;
  Float_t Jet1PuppiPhi;
  Float_t Jet1PuppiMass;

  Float_t Jet1PuppiSDmass;
  Float_t Jet1PuppiSDmassCorr;
  Float_t Jet1PuppiSDmassCorrUp;
  Float_t Jet1PuppiSDmassCorrDn;
  Float_t Jet1PuppiSDmassCorrL23Smear;
  Float_t Jet1PuppiSDmassCorrL23SmearUp;
  Float_t Jet1PuppiSDmassCorrL23SmearDn;
  Float_t Jet1PuppiSDpt;
  Float_t Jet1PuppiSDptCorr;
  Float_t Jet1PuppiSDptCorrUp;
  Float_t Jet1PuppiSDptCorrDn;
  Float_t Jet1PuppiSDptCorrL23Smear;
  Float_t Jet1PuppiSDptCorrL23SmearUp;
  Float_t Jet1PuppiSDptCorrL23SmearDn;
  Float_t Jet1PuppiSDeta;
  Float_t Jet1PuppiSDphi;

  Float_t Jet1PuppiTau1;
  Float_t Jet1PuppiTau2;
  Float_t Jet1PuppiTau3;
  Float_t Jet1PuppiTau4;
  Float_t Jet1PuppiTau32;
  Float_t Jet1PuppiTau21;
  Float_t Jet1PuppiSDsubjet0bdisc;
  Float_t Jet1PuppiSDsubjet1bdisc;
  Float_t Jet1PuppiSDmaxbdisc;
  Float_t Jet1PuppiSDmaxbdiscflavHadron;
  Float_t Jet1PuppiSDmaxbdiscflavParton;
  Float_t Jet1PuppiSDsubjet0pt;
  Float_t Jet1PuppiSDsubjet0eta;
  Float_t Jet1PuppiSDsubjet0phi;
  Float_t Jet1PuppiSDsubjet0mass;
  Float_t Jet1PuppiSDsubjet0area;
  Float_t Jet1PuppiSDsubjet0flavHadron;
  Float_t Jet1PuppiSDsubjet0flavParton;
  Float_t Jet1PuppiSDsubjet0tau1;
  Float_t Jet1PuppiSDsubjet0tau2;
  Float_t Jet1PuppiSDsubjet0tau3;
  Float_t Jet1PuppiSDsubjet1pt;
  Float_t Jet1PuppiSDsubjet1eta;
  Float_t Jet1PuppiSDsubjet1phi;
  Float_t Jet1PuppiSDsubjet1mass;
  Float_t Jet1PuppiSDsubjet1area;
  Float_t Jet1PuppiSDsubjet1flavHadron;
  Float_t Jet1PuppiSDsubjet1flavParton;
  Float_t Jet1PuppiSDsubjet1tau1;
  Float_t Jet1PuppiSDsubjet1tau2;
  Float_t Jet1PuppiSDsubjet1tau3;
  Float_t Jet1CHF;
  Float_t Jet1NHF;
  Float_t Jet1CM;
  Float_t Jet1NM;
  Float_t Jet1NEF;
  Float_t Jet1CEF;
  Float_t Jet1MF;
  Float_t Jet1Mult;
  Float_t Jet1PuppiCHF;
  Float_t Jet1PuppiNHF;
  Float_t Jet1PuppiCM;
  Float_t Jet1PuppiNM;
  Float_t Jet1PuppiNEF;
  Float_t Jet1PuppiCEF;
  Float_t Jet1PuppiMF;
  Float_t Jet1PuppiMult;
  Float_t Jet1MassCorrFactor;
  Float_t Jet1MassCorrFactorUp;
  Float_t Jet1MassCorrFactorDn;
  Float_t Jet1CorrFactor;
  Float_t Jet1CorrFactorUp;
  Float_t Jet1CorrFactorDn;
  Float_t Jet1PtSmearFactor;
  Float_t Jet1PtSmearFactorUp;
  Float_t Jet1PtSmearFactorDn;
  Float_t Jet1PuppiMassCorrFactor;
  Float_t Jet1PuppiMassCorrFactorUp;
  Float_t Jet1PuppiMassCorrFactorDn;
  Float_t Jet1PuppiCorrFactor;
  Float_t Jet1PuppiCorrFactorUp;
  Float_t Jet1PuppiCorrFactorDn;
  Float_t Jet1PuppiPtSmearFactor;
  Float_t Jet1PuppiPtSmearFactorUp;
  Float_t Jet1PuppiPtSmearFactorDn;
  Float_t Jet1EtaScaleFactor;
  Float_t Jet1PhiScaleFactor;
  // Float_t Jet1MatchedGenJetDR;
  Float_t Jet1MatchedGenJetPt;
  Float_t Jet1MatchedGenJetMass;

  Int_t   Jet1GenMatched_TopHadronic;
  Float_t Jet1GenMatched_TopPt;
  Float_t Jet1GenMatched_TopEta;
  Float_t Jet1GenMatched_TopPhi;
  Float_t Jet1GenMatched_TopMass;
  Float_t Jet1GenMatched_bPt;
  Float_t Jet1GenMatched_WPt;
  Float_t Jet1GenMatched_Wd1Pt;
  Float_t Jet1GenMatched_Wd2Pt;
  Float_t Jet1GenMatched_Wd1ID;
  Float_t Jet1GenMatched_Wd2ID;
  Float_t Jet1GenMatched_MaxDeltaRPartonTop;
  Float_t Jet1GenMatched_MaxDeltaRWPartonTop;
  Float_t Jet1GenMatched_MaxDeltaRWPartonW;
  Float_t Jet1GenMatched_DeltaR_t_b;
  Float_t Jet1GenMatched_DeltaR_t_W;
  Float_t Jet1GenMatched_DeltaR_t_Wd1;
  Float_t Jet1GenMatched_DeltaR_t_Wd2;
  Float_t Jet1GenMatched_DeltaR_W_b1;
  Float_t Jet1GenMatched_DeltaR_W_Wd1;
  Float_t Jet1GenMatched_DeltaR_W_Wd2;
  Float_t Jet1GenMatched_DeltaR_Wd1_Wd2;
  Float_t Jet1GenMatched_DeltaR_Wd1_b;
  Float_t Jet1GenMatched_DeltaR_Wd2_b;
  Float_t Jet1GenMatched_DeltaR_jet_t;
  Float_t Jet1GenMatched_DeltaR_jet_W;
  Float_t Jet1GenMatched_DeltaR_jet_b;
  Float_t Jet1GenMatched_DeltaR_jet_Wd1;
  Float_t Jet1GenMatched_DeltaR_jet_Wd2;
  Float_t Jet1GenMatched_DeltaR_pup0_b;
  Float_t Jet1GenMatched_DeltaR_pup0_Wd1;
  Float_t Jet1GenMatched_DeltaR_pup0_Wd2;
  Float_t Jet1GenMatched_DeltaR_pup1_b;
  Float_t Jet1GenMatched_DeltaR_pup1_Wd1;
  Float_t Jet1GenMatched_DeltaR_pup1_Wd2;      
  Float_t Jet1GenMatched_partonPt;
  Float_t Jet1GenMatched_partonEta;
  Float_t Jet1GenMatched_partonPhi;
  Float_t Jet1GenMatched_partonMass;
  Float_t Jet1GenMatched_partonID;
  Float_t Jet1GenMatched_DeltaRjetParton;

  Float_t AllHadMETpx;           
  Float_t AllHadMETpy;           
  Float_t AllHadMETpt;           
  Float_t AllHadMETphi;           
  Float_t AllHadMETsumET;           
  Float_t AllHadNvtx;           
  Float_t AllHadNPUtrue;           
  Float_t AllHadRho;           
  Float_t AllHadEventWeight;    
  Float_t AllHadPUweight; 
  Float_t AllHadPUweight_MBup; 
  Float_t AllHadPUweight_MBdn;        
  Float_t DijetMass;           
  Float_t DijetMassPuppi;           
  Float_t DijetDeltaR;           
  Float_t DijetDeltaPhi;           
  Float_t DijetDeltaRap;           
  Float_t DiGenJetMass;           
  Float_t GenTTmass;           
  Float_t HT;           
  Float_t HT_CorrDn;           
  Float_t HT_CorrUp;           
  Float_t HT_PtSmearNom;           
  Float_t HT_PtSmearUp;           
  Float_t HT_PtSmearDn;           
  Float_t Q2weight_CorrDn;           
  Float_t Q2weight_CorrUp;           
  Float_t NNPDF3weight_CorrDn;           
  Float_t NNPDF3weight_CorrUp;           
  Float_t AllHadRunNum;           
  Float_t AllHadLumiBlock;           
  Float_t AllHadEventNum;    



  // 
  //  .d8888b.                         d8b        888                        888        88888888888                           
  // d88P  Y88b                        Y8P        888                        888            888                               
  // Y88b.                                        888                        888            888                               
  //  "Y888b.    .d88b.  88888b.d88b.  888        888       .d88b.  88888b.  888888         888     888d888  .d88b.   .d88b.  
  //     "Y88b. d8P  Y8b 888 "888 "88b 888        888      d8P  Y8b 888 "88b 888            888     888P"   d8P  Y8b d8P  Y8b 
  //       "888 88888888 888  888  888 888 888888 888      88888888 888  888 888            888     888     88888888 88888888 
  // Y88b  d88P Y8b.     888  888  888 888        888      Y8b.     888 d88P Y88b.          888     888     Y8b.     Y8b.     
  //  "Y8888P"   "Y8888  888  888  888 888        88888888  "Y8888  88888P"   "Y888         888     888      "Y8888   "Y8888  
  //                                                                888                                                       
  //                                                                888                                                       
  //                                                                888       
         
  TTree *TreeSemiLept;
  std::vector<std::string> *SemiLeptTrigNames     = new std::vector<std::string>;
  std::vector<int> *SemiLeptTrigPrescales = new std::vector<int>;
  std::vector<bool> *SemiLeptTrigPass    = new std::vector<bool>;



  std::string SemiLeptTrigAcceptBits;

  Float_t JetPtRaw;      
  Float_t JetEtaRaw;
  Float_t JetPhiRaw;
  Float_t JetMassRaw;
  Float_t JetP;
  Float_t JetPt;
  Float_t JetEta;
  Float_t JetPhi;
  Float_t JetRap;
  Float_t JetEnergy;
  Float_t JetMass;
  Float_t JetArea;
  Float_t JetSDmass;
  Float_t JetSDmassRaw;
  Float_t JetSDmassCorrL23;
  Float_t JetSDmassCorrL23Up;
  Float_t JetSDmassCorrL23Dn;
  Float_t JetSDmassCorrL123;
  Float_t JetSDmassCorrL123Up;
  Float_t JetSDmassCorrL123Dn;
  Float_t JetSDmassCorrL23Smear;
  Float_t JetSDmassCorrL23SmearUp;
  Float_t JetSDmassCorrL23SmearDn;
  Float_t JetSDptRaw;
  Float_t JetSDptCorrL23;
  Float_t JetSDptCorrL23Up;
  Float_t JetSDptCorrL23Dn;
  Float_t JetSDptCorrL123;
  Float_t JetSDptCorrL123Up;
  Float_t JetSDptCorrL123Dn;
  Float_t JetSDptCorrL23Smear;
  Float_t JetSDptCorrL23SmearUp;
  Float_t JetSDptCorrL23SmearDn;
  Float_t JetSDetaRaw;
  Float_t JetSDphiRaw;
  Float_t JetMassPruned;
  Float_t JetMassTrimmed;
  Float_t JetTau1;
  Float_t JetTau2;
  Float_t JetTau3;
  Float_t JetTau4;
  Float_t JetTau32;
  Float_t JetTau21;
  Float_t JetSDsubjet0bdisc;
  Float_t JetSDsubjet1bdisc;
  Float_t JetSDmaxbdisc;
  Float_t JetSDmaxbdiscflavHadron;
  Float_t JetSDmaxbdiscflavParton;
  Float_t JetSDsubjet0pt;
  Float_t JetSDsubjet0mass;
  Float_t JetSDsubjet0eta;
  Float_t JetSDsubjet0phi;
  Float_t JetSDsubjet0area;
  Float_t JetSDsubjet0flavHadron;
  Float_t JetSDsubjet0flavParton;
  Float_t JetSDsubjet0tau1;
  Float_t JetSDsubjet0tau2;
  Float_t JetSDsubjet0tau3;
  Float_t JetSDsubjet1pt;
  Float_t JetSDsubjet1mass;
  Float_t JetSDsubjet1eta;
  Float_t JetSDsubjet1phi;
  Float_t JetSDsubjet1area;
  Float_t JetSDsubjet1flavHadron;
  Float_t JetSDsubjet1flavParton;
  Float_t JetSDsubjet1tau1;
  Float_t JetSDsubjet1tau2;
  Float_t JetSDsubjet1tau3;
  Float_t JetPuppiP;
  Float_t JetPuppiPt;
  Float_t JetPuppiEta;
  Float_t JetPuppiPhi;
  Float_t JetPuppiMass;
  Float_t JetPuppiSDmass;
  Float_t JetPuppiSDmassCorr;
  Float_t JetPuppiSDmassCorrUp;
  Float_t JetPuppiSDmassCorrDn;
  Float_t JetPuppiSDmassCorrL23Smear;
  Float_t JetPuppiSDmassCorrL23SmearUp;
  Float_t JetPuppiSDmassCorrL23SmearDn;
  Float_t JetPuppiSDpt;
  Float_t JetPuppiSDptCorr;
  Float_t JetPuppiSDptCorrUp;
  Float_t JetPuppiSDptCorrDn;
  Float_t JetPuppiSDptCorrL23Smear;
  Float_t JetPuppiSDptCorrL23SmearUp;
  Float_t JetPuppiSDptCorrL23SmearDn;
  Float_t JetPuppiSDeta;
  Float_t JetPuppiSDphi;
  Float_t JetPuppiTau1;
  Float_t JetPuppiTau2;
  Float_t JetPuppiTau3;
  Float_t JetPuppiTau4;
  Float_t JetPuppiTau32;
  Float_t JetPuppiTau21;
  Float_t JetPuppiSDsubjet0bdisc;
  Float_t JetPuppiSDsubjet1bdisc;
  Float_t JetPuppiSDmaxbdisc;
  Float_t JetPuppiSDmaxbdiscflavHadron;
  Float_t JetPuppiSDmaxbdiscflavParton;
  Float_t JetPuppiSDsubjet0pt;
  Float_t JetPuppiSDsubjet0mass;
  Float_t JetPuppiSDsubjet0eta;
  Float_t JetPuppiSDsubjet0phi;
  Float_t JetPuppiSDsubjet0area;
  Float_t JetPuppiSDsubjet0flavHadron;
  Float_t JetPuppiSDsubjet0flavParton;
  Float_t JetPuppiSDsubjet0tau1;
  Float_t JetPuppiSDsubjet0tau2;
  Float_t JetPuppiSDsubjet0tau3;
  Float_t JetPuppiSDsubjet1pt;
  Float_t JetPuppiSDsubjet1mass;
  Float_t JetPuppiSDsubjet1eta;
  Float_t JetPuppiSDsubjet1phi;
  Float_t JetPuppiSDsubjet1area;
  Float_t JetPuppiSDsubjet1flavHadron;
  Float_t JetPuppiSDsubjet1flavParton;
  Float_t JetPuppiSDsubjet1tau1;
  Float_t JetPuppiSDsubjet1tau2;
  Float_t JetPuppiSDsubjet1tau3;
  Float_t JetCHF;
  Float_t JetNHF;
  Float_t JetCM;
  Float_t JetNM;
  Float_t JetNEF;
  Float_t JetCEF;
  Float_t JetMF;
  Float_t JetMult;
  Float_t JetPuppiCHF;
  Float_t JetPuppiNHF;
  Float_t JetPuppiCM;
  Float_t JetPuppiNM;
  Float_t JetPuppiNEF;
  Float_t JetPuppiCEF;
  Float_t JetPuppiMF;
  Float_t JetPuppiMult;
  Float_t JetMassCorrFactor;
  Float_t JetMassCorrFactorUp;
  Float_t JetMassCorrFactorDn;
  Float_t JetCorrFactor;
  Float_t JetCorrFactorUp;
  Float_t JetCorrFactorDn;
  Float_t JetPtSmearFactor;
  Float_t JetPtSmearFactorUp;
  Float_t JetPtSmearFactorDn;
  Float_t JetPuppiMassCorrFactor;
  Float_t JetPuppiMassCorrFactorUp;
  Float_t JetPuppiMassCorrFactorDn;
  Float_t JetPuppiCorrFactor;
  Float_t JetPuppiCorrFactorUp;
  Float_t JetPuppiCorrFactorDn;
  Float_t JetPuppiPtSmearFactor;
  Float_t JetPuppiPtSmearFactorUp;
  Float_t JetPuppiPtSmearFactorDn;
  Float_t JetEtaScaleFactor;
  Float_t JetPhiScaleFactor;
  // Float_t JetMatchedGenJetDR;
  Float_t JetMatchedGenJetPt;
  Float_t JetMatchedGenJetMass;
  Int_t   JetGenMatched_TopHadronic;
  Float_t JetGenMatched_TopPt;
  Float_t JetGenMatched_TopEta;
  Float_t JetGenMatched_TopPhi;
  Float_t JetGenMatched_TopMass;
  Float_t JetGenMatched_bPt;
  Float_t JetGenMatched_WPt;
  Float_t JetGenMatched_Wd1Pt;
  Float_t JetGenMatched_Wd2Pt;
  Float_t JetGenMatched_Wd1ID;
  Float_t JetGenMatched_Wd2ID;
  Float_t JetGenMatched_MaxDeltaRPartonTop;
  Float_t JetGenMatched_MaxDeltaRWPartonTop;
  Float_t JetGenMatched_MaxDeltaRWPartonW;
  Float_t JetGenMatched_DeltaR_t_b;
  Float_t JetGenMatched_DeltaR_t_W;
  Float_t JetGenMatched_DeltaR_t_Wd1;
  Float_t JetGenMatched_DeltaR_t_Wd2;
  Float_t JetGenMatched_DeltaR_W_b1;
  Float_t JetGenMatched_DeltaR_W_Wd1;
  Float_t JetGenMatched_DeltaR_W_Wd2;
  Float_t JetGenMatched_DeltaR_Wd1_Wd2;
  Float_t JetGenMatched_DeltaR_Wd1_b;
  Float_t JetGenMatched_DeltaR_Wd2_b;
  Float_t JetGenMatched_DeltaR_jet_t;
  Float_t JetGenMatched_DeltaR_jet_W;
  Float_t JetGenMatched_DeltaR_jet_b;
  Float_t JetGenMatched_DeltaR_jet_Wd1;
  Float_t JetGenMatched_DeltaR_jet_Wd2;
  Float_t JetGenMatched_DeltaR_pup0_b;
  Float_t JetGenMatched_DeltaR_pup0_Wd1;
  Float_t JetGenMatched_DeltaR_pup0_Wd2;
  Float_t JetGenMatched_DeltaR_pup1_b;
  Float_t JetGenMatched_DeltaR_pup1_Wd1;
  Float_t JetGenMatched_DeltaR_pup1_Wd2;
  Float_t JetGenMatched_partonPt;
  Float_t JetGenMatched_partonEta;
  Float_t JetGenMatched_partonPhi;
  Float_t JetGenMatched_partonMass;
  Float_t JetGenMatched_partonID;
  Float_t JetGenMatched_DeltaRjetParton;
  Float_t SemiLeptMETpx;
  Float_t SemiLeptMETpy;
  Float_t SemiLeptMETpt;
  Float_t SemiLeptMETphi;
  Float_t SemiLeptMETsumET;
  Float_t SemiLeptNvtx;
  Float_t SemiLeptNPUtrue;
  Float_t SemiLeptRho;
  Float_t SemiLeptEventWeight;
  Float_t SemiLeptPUweight;
  Float_t SemiLeptPUweight_MBup;
  Float_t SemiLeptPUweight_MBdn;


  Float_t SemiLeptGenTTmass;

  Float_t HTlep;
  Float_t ST;                
  Float_t ST_CorrDn;                
  Float_t ST_CorrUp;                
  Float_t ST_PtSmearNom;                
  Float_t ST_PtSmearUp;                
  Float_t ST_PtSmearDn;   

  Float_t SemiLeptQ2weight_CorrDn;
  Float_t SemiLeptQ2weight_CorrUp;
  Float_t SemiLeptNNPDF3weight_CorrDn;
  Float_t SemiLeptNNPDF3weight_CorrUp;
  Float_t SemiLeptRunNum;
  Float_t SemiLeptLumiBlock;
  Float_t SemiLeptEventNum;
  Int_t   SemiLeptPassMETFilters;       
        
  Float_t AK4dRminPt;
  Float_t AK4dRminEta;
  Float_t AK4dRminPhi;
  Float_t AK4dRminMass;
  Float_t AK4dRminBdisc;
  Float_t AK4dRminLep;
  Float_t AK4BtagdRminPt;
  Float_t AK4BtagdRminBdisc;
  Float_t AK4BtagdRminLep;
  Int_t   LepHemiContainsAK4BtagLoose;
  Int_t   LepHemiContainsAK4BtagMedium;
  Int_t   LepHemiContainsAK4BtagTight;


  Float_t LeptonPhi;
  Float_t LeptonPt;
  Float_t LeptonEta;
  Float_t LeptonMass;
  Float_t PtRel;
  Int_t   LeptonIsMu;
  Int_t   MuTight;
  Int_t   MuMedium;
  Float_t DeltaRJetLep; 
  Float_t DeltaPhiJetLep; 
  Float_t MuIso;
  Float_t Elecron_absiso;
  Float_t Elecron_relIsoWithDBeta;
  Float_t Elecron_absiso_EA;
  Float_t Elecron_relIsoWithEA;

  //MY ADDITION
  Int_t extra_muon = 0;
  Int_t extra_electron = 0;
  Int_t no_leptons = 0;
  Int_t loose_muon = 0;
  Int_t loose_electron = 0;
  Int_t no_jets = 0;
  Int_t lonesome_jet = 0;
  Int_t loose_jet = 0;
  Int_t extra_merged_jet = 0;
  Int_t met_small = 0;
  Int_t bad_dijet = 0;
  Int_t undefined_met = 0;
  Int_t no_merged_jet = 0;

  Int_t good_event = -1;
  
  
  double anom_weight1;
  double anom_weight2;
  double anom_weight3;
  double anom_weight4;
  double anom_weight5;
  double anom_weight6;
  double anom_weight7;
  double anom_weight8;
  double anom_weight9;
  double anom_weight10;
  double anom_weight11;
  double anom_weight12;
  double anom_weight13;
  double anom_weight14;
  double anom_weight15;
  double anom_weight16;
  double anom_weight17;
  double anom_weight18;
  double anom_weight19;
  double anom_weight20;
  double anom_weight21;
  double anom_weight22;
  double anom_weight23;
  double anom_weight24;
  double anom_weight25;
  double anom_weight26;
  double anom_weight27;
  double anom_weight28;
  double anom_weight29;
  double anom_weight30;
  double anom_weight31;
  double anom_weight32;
  double anom_weight33;
  double anom_weight34;
  double anom_weight35;
  double anom_weight36;
  double anom_weight37;
  double anom_weight38;
  double anom_weight39;
  double anom_weight40;
  double anom_weight41;
  double anom_weight42;
  double anom_weight43;
  double anom_weight44;
  double anom_weight45;
  double anom_weight46;
  double anom_weight47;
  double anom_weight48;
  double anom_weight49;
  double anom_weight50;
  
  Float_t Wplus_pt = -99.9;
  Float_t Wminus_pt = -99.9;
  Float_t Zneutral_pt = -99.9;
  Float_t quark_pt = -99.9;
  //Float_t q_pt = -99.9;
  Float_t antiquark_pt = -99.9;
  //Float_t aq_pt = -99.9;
  Float_t lepton_pt = -99.9;
  Float_t lep_pt = -99.9;
  Float_t neutrino_pt = -99.9;
  //Float_t neu_pt = -99.9;
  Float_t reco_jet0_pt = -99.9;
  Float_t reco_jet1_pt = -99.9;
  Float_t ak8jet0_pt = -99.9;
  Float_t ak8jet1_pt = -99.9;
  Float_t ak8jetONLY_pt = -99.9;
  //Float_t pupjet0_pt = -99.9;
  //Float_t pupjet1_pt = -99.9;
  //Float_t ak4genjet0_pt = -99.9;
  //Float_t ak4genjet1_pt = -99.9;
  Float_t ak8genjet0_pt = -99.9;
  Float_t ak8genjet1_pt = -99.9;
  Float_t reco_electron_pt = -99.9;
  Float_t reco_muon_pt = -99.9;
  Float_t reco_lepton_pt = -99.9;
  Float_t reco_met_pt = -99.9;
  //Float_t subjet0_pt = -99.9;
  //Float_t subjet1_pt = -99.9;
  
  //Float_t subjet2_pt = -99.9;
  //Float_t subjet3_pt = -99.9;
  //Float_t subjet4_pt = -99.9;
  //Float_t subjet5_pt = -99.9;
  //Float_t subjet6_pt = -99.9;
  
  Float_t ak8jet2_pt = -99.9;
  Float_t ak8jet3_pt = -99.9;

  Float_t Wplus_px = -99.9;
  Float_t Wminus_px = -99.9;
  Float_t Zneutral_px = -99.9;
  Float_t quark_px = -99.9;
  Float_t antiquark_px = -99.9;
  Float_t lepton_px = -99.9;
  Float_t neutrino_px = -99.9;
  Float_t reco_jet0_px = -99.9;
  Float_t reco_jet1_px = -99.9;
  Float_t ak8jet0_px = -99.9;
  Float_t ak8jet1_px = -99.9;
  Float_t ak8jetONLY_px = -99.9;
  //Float_t pupjet0_px = -99.9;
  //Float_t pupjet1_px = -99.9;
  //Float_t ak4genjet0_px = -99.9;
  //Float_t ak4genjet1_px = -99.9;
  Float_t ak8genjet0_px = -99.9;
  Float_t ak8genjet1_px = -99.9;
  Float_t reco_electron_px = -99.9;
  Float_t reco_muon_px = -99.9;
  Float_t reco_lepton_px = -99.9;
  Float_t reco_met_px = -99.9;
  //Float_t subjet0_px = -99.9;
  //Float_t subjet1_px = -99.9;

  //Float_t subjet2_px = -99.9;
  //Float_t subjet3_px = -99.9;
  //Float_t subjet4_px = -99.9;
  //Float_t subjet5_px = -99.9;
  //Float_t subjet6_px = -99.9;
  
  Float_t ak8jet2_px = -99.9;
  Float_t ak8jet3_px = -99.9;
  

  Float_t Wplus_py = -99.9;
  Float_t Wminus_py = -99.9;
  Float_t Zneutral_py = -99.9;
  Float_t quark_py = -99.9;
  Float_t antiquark_py = -99.9;
  Float_t lepton_py = -99.9;
  Float_t neutrino_py = -99.9;
  Float_t reco_jet0_py = -99.9;
  Float_t reco_jet1_py = -99.9;
  Float_t ak8jet0_py = -99.9;
  Float_t ak8jet1_py = -99.9;
  Float_t ak8jetONLY_py = -99.9;
  //Float_t pupjet0_py = -99.9;
  //Float_t pupjet1_py = -99.9;
  //Float_t ak4genjet0_py = -99.9;
  //Float_t ak4genjet1_py = -99.9;
  Float_t ak8genjet0_py = -99.9;
  Float_t ak8genjet1_py = -99.9;
  Float_t reco_electron_py = -99.9;
  Float_t reco_muon_py = -99.9;
  Float_t reco_lepton_py = -99.9;
  Float_t reco_met_py = -99.9;
  //Float_t subjet0_py = -99.9;
  //Float_t subjet1_py = -99.9;

  //Float_t subjet2_py = -99.9;
  //Float_t subjet3_py = -99.9;
  //Float_t subjet4_py = -99.9;
  //Float_t subjet5_py = -99.9;
  //Float_t subjet6_py = -99.9;
  
  Float_t ak8jet2_py = -99.9;
  Float_t ak8jet3_py = -99.9;

  Float_t Wplus_pz = -99.9;
  Float_t Wminus_pz = -99.9;
  Float_t Zneutral_pz = -99.9;
  Float_t quark_pz = -99.9;
  Float_t antiquark_pz = -99.9;
  Float_t lepton_pz = -99.9;
  Float_t neutrino_pz = -99.9;
  Float_t reco_jet0_pz = -99.9;
  Float_t reco_jet1_pz = -99.9;
  Float_t ak8jet0_pz = -99.9;
  Float_t ak8jet1_pz = -99.9;
  Float_t ak8jetONLY_pz = -99.9;
  //Float_t pupjet0_pz = -99.9;
  //Float_t pupjet1_pz = -99.9;
  //Float_t ak4genjet0_pz = -99.9;
  //Float_t ak4genjet1_pz = -99.9;
  Float_t ak8genjet0_pz = -99.9;
  Float_t ak8genjet1_pz = -99.9;
  Float_t reco_electron_pz = -99.9;
  Float_t reco_muon_pz = -99.9;
  Float_t reco_lepton_pz = -99.9;
  Float_t reco_met_pz = -99.9;
  //Float_t subjet0_pz = -99.9;
  //Float_t subjet1_pz = -99.9;

  //Float_t subjet2_pz = -99.9;
  //Float_t subjet3_pz = -99.9;
  //Float_t subjet4_pz = -99.9;
  //Float_t subjet5_pz = -99.9;
  //Float_t subjet6_pz = -99.9;
  
  Float_t ak8jet2_pz = -99.9;
  Float_t ak8jet3_pz = -99.9;

  Float_t Wplus_e = -99.9;
  Float_t Wminus_e = -99.9;
  Float_t Zneutral_e = -99.9;
  Float_t quark_e = -99.9;
  Float_t antiquark_e = -99.9;
  Float_t lepton_e = -99.9;
  Float_t neutrino_e = -99.9;
  Float_t reco_jet0_e = -99.9;
  Float_t reco_jet1_e = -99.9;
  Float_t ak8jet0_e = -99.9;
  Float_t ak8jet1_e = -99.9;
  Float_t ak8jetONLY_e = -99.9;
  //Float_t pupjet0_e = -99.9;
  //Float_t pupjet1_e = -99.9;
  //Float_t ak4genjet0_e = -99.9;
  //Float_t ak4genjet1_e = -99.9;
  Float_t ak8genjet0_e = -99.9;
  Float_t ak8genjet1_e = -99.9;
  Float_t reco_electron_e = -99.9;
  Float_t reco_muon_e = -99.9;
  Float_t reco_lepton_e = -99.9;
  Float_t reco_met_e = -99.9;
  //Float_t subjet0_e = -99.9;
  //Float_t subjet1_e = -99.9;

  //Float_t subjet2_e = -99.9;
  //Float_t subjet3_e = -99.9;
  //Float_t subjet4_e = -99.9;
  //Float_t subjet5_e = -99.9;
  //Float_t subjet6_e = -99.9;

  Float_t ak8jet2_e = -99.9;
  Float_t ak8jet3_e = -99.9;

  Float_t Wplus_eta = -99.9;
  Float_t Wminus_eta = -99.9;
  Float_t Zneutral_eta = -99.9;
  Float_t quark_eta = -99.9;
  Float_t antiquark_eta = -99.9;
  Float_t lepton_eta = -99.9;
  Float_t neutrino_eta = -99.9;
  Float_t reco_jet0_eta = -99.9;
  Float_t reco_jet1_eta = -99.9;
  Float_t ak8jet0_eta = -99.9;
  Float_t ak8jet1_eta = -99.9;
  Float_t ak8jetONLY_eta = -99.9;
  //Float_t pupjet0_eta = -99.9;
  //Float_t pupjet1_eta = -99.9;
  //Float_t ak4genjet0_eta = -99.9;
  //Float_t ak4genjet1_eta = -99.9;
  Float_t ak8genjet0_eta = -99.9;
  Float_t ak8genjet1_eta = -99.9;
  Float_t reco_electron_eta = -99.9;
  Float_t reco_muon_eta = -99.9;
  Float_t reco_lepton_eta = -99.9;
  Float_t reco_met_eta = -99.9;
  //Float_t subjet0_eta = -99.9;
  //Float_t subjet1_eta = -99.9;

  //Float_t subjet2_eta = -99.9;
  //Float_t subjet3_eta = -99.9;
  //Float_t subjet4_eta = -99.9;
  //Float_t subjet5_eta = -99.9;
  //Float_t subjet6_eta = -99.9;

  Float_t Wplus_et = -99.9;
  Float_t Wminus_et = -99.9;
  Float_t Zneutral_et = -99.9;
  Float_t quark_et = -99.9;
  Float_t antiquark_et = -99.9;
  Float_t lepton_et = -99.9;
  Float_t neutrino_et = -99.9;
  Float_t reco_jet0_et = -99.9;
  Float_t reco_jet1_et = -99.9;
  Float_t ak8jet0_et = -99.9;
  Float_t ak8jet1_et = -99.9;
  Float_t ak8jetONLY_et = -99.9;
  //Float_t pupjet0_et = -99.9;
  //Float_t pupjet1_et = -99.9;
  //Float_t ak4genjet0_et = -99.9;
  //Float_t ak4genjet1_et = -99.9;
  Float_t ak8genjet0_et = -99.9;
  Float_t ak8genjet1_et = -99.9;
  Float_t reco_electron_et = -99.9;
  Float_t reco_muon_et = -99.9;
  Float_t reco_lepton_et = -99.9;
  Float_t reco_met_et = -99.9;
  //Float_t subjet0_et = -99.9;
  //Float_t subjet1_et = -99.9;

  //Float_t subjet2_et = -99.9;
  //Float_t subjet3_et = -99.9;
  //Float_t subjet4_et = -99.9;
  //Float_t subjet5_et = -99.9;
  //Float_t subjet6_et = -99.9;

  Float_t ak8jet2_eta = -99.9;
  Float_t ak8jet3_eta = -99.9;

  Float_t ak8jet2_et = -99.9;
  Float_t ak8jet3_et = -99.9;

  Float_t Wplus_theta = -99.9;
  Float_t Wminus_theta = -99.9;
  Float_t Zneutral_theta = -99.9;
  Float_t quark_theta = -99.9;
  Float_t antiquark_theta = -99.9;
  Float_t lepton_theta = -99.9;
  Float_t neutrino_theta = -99.9;
  Float_t reco_jet0_theta = -99.9;
  Float_t reco_jet1_theta = -99.9;
  Float_t ak8jet0_theta = -99.9;
  Float_t ak8jet1_theta = -99.9;
  Float_t ak8jetONLY_theta = -99.9;
  //Float_t pupjet0_theta = -99.9;
  //Float_t pupjet1_theta = -99.9;
  //Float_t ak4genjet0_theta = -99.9;
  //Float_t ak4genjet1_theta = -99.9;
  Float_t ak8genjet0_theta = -99.9;
  Float_t ak8genjet1_theta = -99.9;
  Float_t reco_electron_theta = -99.9;
  Float_t reco_muon_theta = -99.9;
  Float_t reco_lepton_theta = -99.9;
  Float_t reco_met_theta = -99.9;
  //Float_t subjet0_theta = -99.9;
  //Float_t subjet1_theta = -99.9;

  //Float_t subjet2_theta = -99.9;
  //Float_t subjet3_theta = -99.9;
  //Float_t subjet4_theta = -99.9;
  //Float_t subjet5_theta = -99.9;
  //Float_t subjet6_theta = -99.9;

  Float_t ak8jet2_theta = -99.9;
  Float_t ak8jet3_theta = -99.9;

  Float_t Wplus_phi = -99.9;
  Float_t Wminus_phi = -99.9;
  Float_t Zneutral_phi = -99.9;
  Float_t quark_phi = -99.9;
  Float_t antiquark_phi = -99.9;
  Float_t lepton_phi = -99.9;
  Float_t neutrino_phi = -99.9;
  Float_t reco_jet0_phi = -99.9;
  Float_t reco_jet1_phi = -99.9;
  Float_t ak8jet0_phi = -99.9;
  Float_t ak8jet1_phi = -99.9;
  Float_t ak8jetONLY_phi = -99.9;
  //Float_t pupjet0_phi = -99.9;
  //Float_t pupjet1_phi = -99.9;
  //Float_t ak4genjet0_phi = -99.9;
  //Float_t ak4genjet1_phi = -99.9;
  Float_t ak8genjet0_phi = -99.9;
  Float_t ak8genjet1_phi = -99.9;
  Float_t reco_electron_phi = -99.9;
  Float_t reco_muon_phi = -99.9;
  Float_t reco_lepton_phi = -99.9;
  Float_t reco_met_phi = -99.9;
  //Float_t subjet0_phi = -99.9;
  //Float_t subjet1_phi = -99.9;

  //Float_t subjet2_phi = -99.9;
  //Float_t subjet3_phi = -99.9;
  //Float_t subjet4_phi = -99.9;
  //Float_t subjet5_phi = -99.9;
  //Float_t subjet6_phi = -99.9;

  Float_t ak8jet2_phi = -99.9;
  Float_t ak8jet3_phi = -99.9;

  Float_t Wplus_y = -99.9;
  Float_t Wminus_y = -99.9;
  Float_t Zneutral_y = -99.9;
  Float_t quark_y = -99.9;
  Float_t antiquark_y = -99.9;
  Float_t lepton_y = -99.9;
  Float_t neutrino_y = -99.9;
  Float_t reco_jet0_y = -99.9;
  Float_t reco_jet1_y = -99.9;
  Float_t ak8jet0_y = -99.9;
  Float_t ak8jet1_y = -99.9;
  Float_t ak8jetONLY_y = -99.9;
  //Float_t pupjet0_y = -99.9;
  //Float_t pupjet1_y = -99.9;
  //Float_t ak4genjet0_y = -99.9;
  //Float_t ak4genjet1_y = -99.9;
  Float_t ak8genjet0_y = -99.9;
  Float_t ak8genjet1_y = -99.9;
  Float_t reco_electron_y = -99.9;
  Float_t reco_muon_y = -99.9;
  Float_t reco_lepton_y = -99.9;
  Float_t reco_met_y = -99.9;
  //Float_t subjet0_y = -99.9;
  //Float_t subjet1_y = -99.9;

  //Float_t subjet2_y = -99.9;
  //Float_t subjet3_y = -99.9;
  //Float_t subjet4_y = -99.9;
  //Float_t subjet5_y = -99.9;
  //Float_t subjet6_y = -99.9;

  Float_t ak8jet2_y = -99.9;
  Float_t ak8jet3_y = -99.9;

  Float_t interquark_delta = -99.9;

  //Float_t reco_ak4_costheta1 = -99.9;
  //Float_t reco_ak4_costheta2 = -99.9;
  //Float_t reco_ak4_phi = -99.9;
  //Float_t reco_ak4_costhetastar = -99.9;
  //Float_t reco_ak4_phistar1 = -99.9;
  //Float_t reco_ak4_phistar2 = -99.9;

  //Float_t gen_ak4_costheta1 = -99.9;
  //Float_t gen_ak4_costheta2 = -99.9;
  //Float_t gen_ak4_phi = -99.9;
  //Float_t gen_ak4_costhetastar = -99.9;
  //Float_t gen_ak4_phistar1 = -99.9;
  //Float_t gen_ak4_phistar2 = -99.9;

  Float_t reco_costheta1 = -99.9;
  Float_t reco_costheta2 = -99.9;
  Float_t reco_phi = -99.9;
  Float_t reco_costhetastar = -99.9;
  Float_t reco_phistar1 = -99.9;
  Float_t reco_phistar2 = -99.9;

  Float_t reco_ak8_costheta1 = -99.9;
  Float_t reco_ak8_costheta2 = -99.9;
  Float_t reco_ak8_phi = -99.9;
  Float_t reco_ak8_costhetastar = -99.9;
  Float_t reco_ak8_phistar1 = -99.9;
  Float_t reco_ak8_phistar2 = -99.9;

  Float_t gen_ak8_costheta1 = -99.9;
  Float_t gen_ak8_costheta2 = -99.9;
  Float_t gen_ak8_phi = -99.9;
  Float_t gen_ak8_costhetastar = -99.9;
  Float_t gen_ak8_phistar1 = -99.9;
  Float_t gen_ak8_phistar2 = -99.9;

  Float_t reco_puppi_costheta1 = -99.9;
  Float_t reco_puppi_costheta2 = -99.9;
  Float_t reco_puppi_phi = -99.9;
  Float_t reco_puppi_costhetastar = -99.9;
  Float_t reco_puppi_phistar1 = -99.9;
  Float_t reco_puppi_phistar2 = -99.9;

  Float_t lhe_gen_costheta1 = -99.9;
  Float_t lhe_gen_costheta2 = -99.9;
  Float_t lhe_gen_phi = -99.9;
  Float_t lhe_gen_costhetastar = -99.9;
  Float_t lhe_gen_phistar1 = -99.9;
  Float_t lhe_gen_phistar2 = -99.9;

  //Float_t subjet_costheta1 = -99.9;
  //Float_t subjet_costheta2 = -99.9;
  //Float_t subjet_phi = -99.9;
  //Float_t subjet_costhetastar = -99.9;
  //Float_t subjet_phistar1 = -99.9;
  //Float_t subjet_phistar2 = -99.9;

  Float_t goodjet_costheta1 = -99.9;
  Float_t goodjet_costheta2 = -99.9;
  Float_t goodjet_phi = -99.9;
  Float_t goodjet_costhetastar = -99.9;
  Float_t goodjet_phistar1 = -99.9;
  Float_t goodjet_phistar2 = -99.9;

  //Float_t goodsubjet_costheta1 = -99.9;
  //Float_t goodsubjet_costheta2 = -99.9;
  //Float_t goodsubjet_phi = -99.9;
  //Float_t goodsubjet_costhetastar = -99.9;
  //Float_t goodsubjet_phistar1 = -99.9;
  //Float_t goodsubjet_phistar2 = -99.9;

  Int_t lep_gen = 0;
  Int_t reco_lep_gen = 0;

  Float_t lep_delta = -99.9;
  Float_t met_delta = -99.9;
  
  Float_t reco_jet_delta0 = -99.9;
  Float_t reco_jet_delta1 = -99.9;

  Float_t ak8jet_delta0 = -99.9;
  Float_t ak8jet_delta1 = -99.9;

  //Float_t pupjet_delta0 = -99.9;
  //Float_t pupjet_delta1 = -99.9;

  //Float_t ak4genjet_delta0 = -99.9;
  //Float_t ak4genjet_delta1 = -99.9;

  Float_t ak8genjet_delta0 = -99.9;
  Float_t ak8genjet_delta1 = -99.9;

  //Float_t subjet_delta0 = -99.9;
  //Float_t subjet_delta1 = -99.9;

  Int_t muon_count = -99;
  Int_t electron_count = -99;
  
  //Int_t reco_ak4_jet_count = -99;
  Int_t reco_ak8_jet_count = -99;
  //Int_t reco_puppi_jet_count = -99;
  
  //Int_t gen_ak4_jet_count = -99;
  Int_t gen_ak8_jet_count = -99;

  //Int_t subjet_count = -99;

  Bool_t good_jet_event;

  Float_t masswidth = -99.9;

  /*TEMPLATE;
  Float_t _pt = -99.9;
  Float_t _px = -99.9;
  Float_t _py = -99.9;
  Float_t _pz = -99.9;
  Float_t _e = -99.9;
  Float_t _eta = -99.9;
  Float_t _et = -99.9;
  Float_t _theta = -99.9;
  Float_t _phi = -99.9;
  Float_t _y = -99.9;
  END TEMPLATE;*/

  Float_t goodjet0_pt = -99.9;
  Float_t goodjet0_px = -99.9;
  Float_t goodjet0_py = -99.9;
  Float_t goodjet0_pz = -99.9;
  Float_t goodjet0_e = -99.9;
  Float_t goodjet0_eta = -99.9;
  Float_t goodjet0_et = -99.9;
  Float_t goodjet0_theta = -99.9;
  Float_t goodjet0_phi = -99.9;
  Float_t goodjet0_y = -99.9;

  Float_t goodjet1_pt = -99.9;
  Float_t goodjet1_px = -99.9;
  Float_t goodjet1_py = -99.9;
  Float_t goodjet1_pz = -99.9;
  Float_t goodjet1_e = -99.9;
  Float_t goodjet1_eta = -99.9;
  Float_t goodjet1_et = -99.9;
  Float_t goodjet1_theta = -99.9;
  Float_t goodjet1_phi = -99.9;
  Float_t goodjet1_y = -99.9;

  Float_t goodjetonly_pt = -99.9;
  Float_t goodjetonly_px = -99.9;
  Float_t goodjetonly_py = -99.9;
  Float_t goodjetonly_pz = -99.9;
  Float_t goodjetonly_e = -99.9;
  Float_t goodjetonly_eta = -99.9;
  Float_t goodjetonly_et = -99.9;
  Float_t goodjetonly_theta = -99.9;
  Float_t goodjetonly_phi = -99.9;
  Float_t goodjetonly_y = -99.9;

  //Float_t goodsubjet0_pt = -99.9;
  //Float_t goodsubjet0_px = -99.9;
  //Float_t goodsubjet0_py = -99.9;
  //Float_t goodsubjet0_pz = -99.9;
  //Float_t goodsubjet0_e = -99.9;
  //Float_t goodsubjet0_eta = -99.9;
  //Float_t goodsubjet0_theta = -99.9;
  //Float_t goodsubjet0_phi = -99.9;
  //Float_t goodsubjet0_y = -99.9;

  /*Float_t goodsubjet1_pt = -99.9;
  Float_t goodsubjet1_px = -99.9;
  Float_t goodsubjet1_py = -99.9;
  Float_t goodsubjet1_pz = -99.9;
  Float_t goodsubjet1_e = -99.9;
  Float_t goodsubjet1_eta = -99.9;
  Float_t goodsubjet1_theta = -99.9;
  Float_t goodsubjet1_phi = -99.9;
  Float_t goodsubjet1_y = -99.9;*/

  TLorentzVector GenJetMatched0;
  TLorentzVector GenJetMatched1;


 //DELTA TREES

 /* TTree *DeltaTree0;
  TTree *DeltaTree1;
  TTree *DeltaTree2;
  TTree *DeltaTree3;*/
  //TTree *DeltaTree4;
  //TTree *DeltaTree5;

  /* Float_t ak8_jet0lowq = -99.9;
  Float_t ak8_jet0highq = -99.9;

  Float_t ak8_jet1lowq = -99.9;
  Float_t ak8_jet1highq = -99.9;

  Float_t ak8_jet2lowq = -99.9;
  Float_t ak8_jet2highq = -99.9;

  Float_t ak8_jet3lowq = -99.9;
  Float_t ak8_jet3highq = -99.9;*/

  //Float_t ak8_jet4lowq = -99.9;
  //Float_t ak8_jet4highq = -99.9;
  
  //Float_t ak8_jet5lowq = -99.9;
  //Float_t ak8_jet5highq = -99.9;

  Float_t ak8jet_delta2 = -99.9;
  Float_t ak8jet_delta3 = -99.9;

  //Float_t subjet_delta2 = -99.9;
  //Float_t subjet_delta3 = -99.9;
  //Float_t subjet_delta4 = -99.9;
  //Float_t subjet_delta5 = -99.9;
  //Float_t subjet_delta6 = -99.9;
};
//std::cout << "we just finished class delcaration.\n";

//
// constructors and destructor
//
B2GTTbarTreeMaker::B2GTTbarTreeMaker(const edm::ParameterSet& iConfig):
  ak4jetToken_(consumes<pat::JetCollection>(edm::InputTag("slimmedJets"))),
  ak8jetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("ak8chsInput"))),//edm::InputTag("slimmedJetsAK8"))),
  puppijetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("ak8puppiInput"))),
  ak8CHSSoftDropSubjetsToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("ak8chsSubjetsInput"))),
  ak8PuppiSoftDropSubjetsToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("ak8puppiSubjetsInput"))),
  ak4genjetToken_(consumes<reco::GenJetCollection>(edm::InputTag("slimmedGenJets"))),
  ak8genjetToken_(consumes<reco::GenJetCollection>(edm::InputTag("slimmedGenJetsAK8"))),
  prunedGenToken_(consumes<edm::View<reco::GenParticle> >(edm::InputTag("prunedGenParticles"))),
  rhoToken_(consumes<double>(edm::InputTag("fixedGridRhoFastjetAll"))),
  vtxToken_(consumes<std::vector<reco::Vertex> >(edm::InputTag("offlineSlimmedPrimaryVertices"))),
  triggerResultsMETFilterToken_(consumes<edm::TriggerResults>(edm::InputTag("TriggerResults", "", "RECO"))),  //"PAT"
  triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerBits"))),//"TriggerResults", "", "HLT2"))),
  triggerObjects_(consumes<pat::TriggerObjectStandAloneCollection>(edm::InputTag("patTrigger"))),
  triggerPrescales_(consumes<pat::PackedTriggerPrescales>(edm::InputTag("patTrigger"))),  //   selectedPatTrigger))),
  badMuonFilterToken_(consumes<bool>(edm::InputTag("BadPFMuonFilter", ""))),
  badChargedCandidateFilterToken_(consumes<bool>(edm::InputTag("BadChargedCandidateFilter", ""))),
  muonToken_(consumes<pat::MuonCollection>(edm::InputTag("slimmedMuons"))),
  electronToken_(consumes<pat::ElectronCollection>(edm::InputTag("slimmedElectrons"))),
  metToken_(consumes<pat::METCollection>(edm::InputTag("slimmedMETs"))),
  pileupInfoToken_(consumes<std::vector<PileupSummaryInfo>>(edm::InputTag("slimmedAddPileupInfo"))),
  theSrc_(consumes<LHEEventProduct>(iConfig.getParameter<edm::InputTag>("theSrc"))),
  pdfToken_(consumes<GenEventInfoProduct>(edm::InputTag("generator"))),
  useToolbox_(iConfig.getParameter<bool>  ("useToolbox")),
  verbose_(iConfig.getParameter<bool>  ("verbose")),
  verboseGen_(iConfig.getParameter<bool>  ("verboseGen")),
  runGenLoop_(iConfig.getParameter<bool>  ("runGenLoop")),
  isZprime_(iConfig.getParameter<bool>  ("isZprime")),
  isttbar_(iConfig.getParameter<bool>  ("isttbar")),
  isRSG_(iConfig.getParameter<bool>  ("isRSG")),
  jecPayloadsAK4chs_ (iConfig.getParameter<std::vector<std::string> >  ("jecPayloadsAK4chs")),
  jecPayloadsAK8chs_ (iConfig.getParameter<std::vector<std::string> >  ("jecPayloadsAK8chs")),
  jecPayloadsAK4pup_ (iConfig.getParameter<std::vector<std::string> >  ("jecPayloadsAK4pup")),
  jecPayloadsAK8pup_ (iConfig.getParameter<std::vector<std::string> >  ("jecPayloadsAK8pup")),
  jerSFtext_ (iConfig.getParameter<std::string>  ("jerSFtext"))//,
  //MY ADDITION
  /*mgreweight1_(consumes<double>(edm::InputTag("allWeights", "mgreweight1", "HLT"))),
  mgreweight2_(consumes<double>(edm::InputTag("allWeights", "mgreweight2", "HLT"))),
  mgreweight3_(consumes<double>(edm::InputTag("allWeights", "mgreweight3", "HLT"))),
  mgreweight4_(consumes<double>(edm::InputTag("allWeights", "mgreweight4", "HLT"))),
  mgreweight5_(consumes<double>(edm::InputTag("allWeights", "mgreweight5", "HLT"))),
  mgreweight6_(consumes<double>(edm::InputTag("allWeights", "mgreweight6", "HLT"))),
  mgreweight7_(consumes<double>(edm::InputTag("allWeights", "mgreweight7", "HLT"))),
  mgreweight8_(consumes<double>(edm::InputTag("allWeights", "mgreweight8", "HLT"))),
  mgreweight9_(consumes<double>(edm::InputTag("allWeights", "mgreweight9", "HLT"))),
  mgreweight10_(consumes<double>(edm::InputTag("allWeights", "mgreweight10", "HLT"))),
  mgreweight11_(consumes<double>(edm::InputTag("allWeights", "mgreweight11", "HLT"))),
  mgreweight12_(consumes<double>(edm::InputTag("allWeights", "mgreweight12", "HLT"))),
  mgreweight13_(consumes<double>(edm::InputTag("allWeights", "mgreweight13", "HLT"))),
  mgreweight14_(consumes<double>(edm::InputTag("allWeights", "mgreweight14", "HLT"))),
  mgreweight15_(consumes<double>(edm::InputTag("allWeights", "mgreweight15", "HLT"))),
  mgreweight16_(consumes<double>(edm::InputTag("allWeights", "mgreweight16", "HLT"))),
  mgreweight17_(consumes<double>(edm::InputTag("allWeights", "mgreweight17", "HLT"))),
  mgreweight18_(consumes<double>(edm::InputTag("allWeights", "mgreweight18", "HLT"))),
  mgreweight19_(consumes<double>(edm::InputTag("allWeights", "mgreweight19", "HLT"))),
  mgreweight20_(consumes<double>(edm::InputTag("allWeights", "mgreweight20", "HLT"))),
  mgreweight21_(consumes<double>(edm::InputTag("allWeights", "mgreweight21", "HLT"))),
  mgreweight22_(consumes<double>(edm::InputTag("allWeights", "mgreweight22", "HLT"))),
  mgreweight23_(consumes<double>(edm::InputTag("allWeights", "mgreweight23", "HLT"))),
  mgreweight24_(consumes<double>(edm::InputTag("allWeights", "mgreweight24", "HLT"))),
  mgreweight25_(consumes<double>(edm::InputTag("allWeights", "mgreweight25", "HLT"))),
  mgreweight26_(consumes<double>(edm::InputTag("allWeights", "mgreweight26", "HLT"))),
  mgreweight27_(consumes<double>(edm::InputTag("allWeights", "mgreweight27", "HLT"))),
  mgreweight28_(consumes<double>(edm::InputTag("allWeights", "mgreweight28", "HLT"))),
  mgreweight29_(consumes<double>(edm::InputTag("allWeights", "mgreweight29", "HLT"))),
  mgreweight30_(consumes<double>(edm::InputTag("allWeights", "mgreweight30", "HLT"))),
  mgreweight31_(consumes<double>(edm::InputTag("allWeights", "mgreweight31", "HLT"))),
  mgreweight32_(consumes<double>(edm::InputTag("allWeights", "mgreweight32", "HLT"))),
  mgreweight33_(consumes<double>(edm::InputTag("allWeights", "mgreweight33", "HLT"))),
  mgreweight34_(consumes<double>(edm::InputTag("allWeights", "mgreweight34", "HLT"))),
  mgreweight35_(consumes<double>(edm::InputTag("allWeights", "mgreweight35", "HLT"))),
  mgreweight36_(consumes<double>(edm::InputTag("allWeights", "mgreweight36", "HLT"))),
  mgreweight37_(consumes<double>(edm::InputTag("allWeights", "mgreweight37", "HLT"))),
  mgreweight38_(consumes<double>(edm::InputTag("allWeights", "mgreweight38", "HLT"))),
  mgreweight39_(consumes<double>(edm::InputTag("allWeights", "mgreweight39", "HLT"))),
  mgreweight40_(consumes<double>(edm::InputTag("allWeights", "mgreweight40", "HLT"))),
  mgreweight41_(consumes<double>(edm::InputTag("allWeights", "mgreweight41", "HLT"))),
  mgreweight42_(consumes<double>(edm::InputTag("allWeights", "mgreweight42", "HLT"))),
  mgreweight43_(consumes<double>(edm::InputTag("allWeights", "mgreweight43", "HLT"))),
  mgreweight44_(consumes<double>(edm::InputTag("allWeights", "mgreweight44", "HLT"))),
  mgreweight45_(consumes<double>(edm::InputTag("allWeights", "mgreweight45", "HLT"))),
  mgreweight46_(consumes<double>(edm::InputTag("allWeights", "mgreweight46", "HLT"))),
  mgreweight47_(consumes<double>(edm::InputTag("allWeights", "mgreweight47", "HLT"))),
  mgreweight48_(consumes<double>(edm::InputTag("allWeights", "mgreweight48", "HLT"))),
  mgreweight49_(consumes<double>(edm::InputTag("allWeights", "mgreweight49", "HLT"))),
  mgreweight50_(consumes<double>(edm::InputTag("allWeights", "mgreweight50", "HLT")))*/
{
  std::cout << "we are starting the constructor.\n";
  std::cout << "B2GTTbarTreeMaker::B2GTTbarTreeMaker" << std::endl;

      
  std::cout << "use tool box is " << useToolbox_ << std::endl;

  outfile.open("special_output.txt");
  
  //RS gluon PDF weights
  LHAPDF::initPDFSet(1, "NNPDF30_lo_as_0130");

  usesResource("TFileService");

  edm::Service<TFileService> fs;

  h_ak8puppi_softDropMass            =  fs->make<TH1D>("h_ak8puppi_softDropMass"           , "", 200, 0, 400);
  h_ak8chs_softDropMass              =  fs->make<TH1D>("h_ak8chs_softDropMass"             , "", 200, 0, 400);
  h_ak8chs_softDropMass_reweighted   =  fs->make<TH1D>("h_ak8chs_softDropMass_reweighted"  , "", 200, 0, 400);
  h_ak8chs_pt                        =  fs->make<TH1D>("h_ak8chs_pt"                       , "", 200, 0, 4000);
  h_ak8chs_pt_reweighted             =  fs->make<TH1D>("h_ak8chs_pt_reweighted"            , "", 200, 0, 4000);
  h_NtrueIntPU                       =  fs->make<TH1D>("h_NtrueIntPU"                      , "", 200, 0, 200);
  h_NPV                              =  fs->make<TH1D>("h_NPV"                             , "", 200, 0, 200);
  h_NPVreweighted                    =  fs->make<TH1D>("h_NPVreweighted"                   , "", 200, 0, 200);



  //
  //       d8888 888 888        888    888               888     88888888888                           
  //      d88888 888 888        888    888               888         888                               
  //     d88P888 888 888        888    888               888         888                               
  //    d88P 888 888 888        8888888888  8888b.   .d88888         888     888d888  .d88b.   .d88b.  
  //   d88P  888 888 888        888    888     "88b d88" 888         888     888P"   d8P  Y8b d8P  Y8b 
  //  d88P   888 888 888 888888 888    888 .d888888 888  888         888     888     88888888 88888888 
  // d8888888888 888 888        888    888 888  888 Y88b 888         888     888     Y8b.     Y8b.     
  //d88P     888 888 888        888    888 "Y888888  "Y88888         888     888      "Y8888   "Y8888  
  //                                                                                                   
        
  TreeAllHad = new TTree("TreeAllHad", "TreeAllHad"); 


  TreeAllHad->Branch("AllHadTrigNames"    , "vector<std::string>", &AllHadTrigNames);
  /*TreeAllHad->Branch("AllHadTrigPrescales"   , "vector<int>", &AllHadTrigPrescales);
  TreeAllHad->Branch("AllHadTrigPass"        , "vector<bool>", &AllHadTrigPass);
  TreeAllHad->Branch("AllHadTrigAcceptBits"  , &AllHadTrigAcceptBits);


  TreeAllHad->Branch("PassMETFilters"                        , & PassMETFilters                     ,    "PassMETFilters/I");                                  
  TreeAllHad->Branch("Jet0PtRaw"                             , & Jet0PtRaw                          ,    "Jet0PtRaw/F");                                  
  TreeAllHad->Branch("Jet0EtaRaw"                            , & Jet0EtaRaw                         ,    "Jet0EtaRaw/F");                                   
  TreeAllHad->Branch("Jet0PhiRaw"                            , & Jet0PhiRaw                         ,    "Jet0PhiRaw/F");                                   
  TreeAllHad->Branch("Jet0MassRaw"                           , & Jet0MassRaw                        ,    "Jet0MassRaw/F");                                    
  TreeAllHad->Branch("Jet0P"                                 , & Jet0P                              ,    "Jet0P/F");                              
  TreeAllHad->Branch("Jet0Pt"                                , & Jet0Pt                             ,    "Jet0Pt/F");                               
  TreeAllHad->Branch("Jet0Eta"                               , & Jet0Eta                            ,    "Jet0Eta/F");                                
  TreeAllHad->Branch("Jet0Phi"                               , & Jet0Phi                            ,    "Jet0Phi/F");                                
  TreeAllHad->Branch("Jet0Rap"                               , & Jet0Rap                            ,    "Jet0Rap/F");                                
  TreeAllHad->Branch("Jet0Energy"                            , & Jet0Energy                         ,    "Jet0Energy/F");                                   
  TreeAllHad->Branch("Jet0Mass"                              , & Jet0Mass                           ,    "Jet0Mass/F");                                 
  TreeAllHad->Branch("Jet0Area"                              , & Jet0Area                           ,    "Jet0Area/F");                                 
  TreeAllHad->Branch("Jet0SDmass"                            , & Jet0SDmass                         ,    "Jet0SDmass/F");                                         
  TreeAllHad->Branch("Jet0SDmassRaw"                         , & Jet0SDmassRaw                      ,    "Jet0SDmassRaw/F");                                               
  TreeAllHad->Branch("Jet0SDmassCorrL23"                     , & Jet0SDmassCorrL23                  ,    "Jet0SDmassCorrL23/F");                                                    
  TreeAllHad->Branch("Jet0SDmassCorrL23Up"                   , & Jet0SDmassCorrL23Up                ,    "Jet0SDmassCorrL23Up/F");                                                      
  TreeAllHad->Branch("Jet0SDmassCorrL23Dn"                   , & Jet0SDmassCorrL23Dn                ,    "Jet0SDmassCorrL23Dn/F");                                                      
  TreeAllHad->Branch("Jet0SDmassCorrL123"                    , & Jet0SDmassCorrL123                 ,    "Jet0SDmassCorrL123/F");                                                      
  TreeAllHad->Branch("Jet0SDmassCorrL123Up"                  , & Jet0SDmassCorrL123Up               ,    "Jet0SDmassCorrL123Up/F");                                                        
  TreeAllHad->Branch("Jet0SDmassCorrL123Dn"                  , & Jet0SDmassCorrL123Dn               ,    "Jet0SDmassCorrL123Dn/F");                                                        
  TreeAllHad->Branch("Jet0SDmassCorrL23Smear"                , & Jet0SDmassCorrL23Smear             ,    "Jet0SDmassCorrL23Smear/F");                                                     
  TreeAllHad->Branch("Jet0SDmassCorrL23SmearUp"              , & Jet0SDmassCorrL23SmearUp           ,    "Jet0SDmassCorrL23SmearUp/F");                                                       
  TreeAllHad->Branch("Jet0SDmassCorrL23SmearDn"              , & Jet0SDmassCorrL23SmearDn           ,    "Jet0SDmassCorrL23SmearDn/F");   
  TreeAllHad->Branch("Jet0SDptRaw"                           , & Jet0SDptRaw                        ,    "Jet0SDptRaw/F");                                               
  TreeAllHad->Branch("Jet0SDptCorrL23"                       , & Jet0SDptCorrL23                    ,    "Jet0SDptCorrL23/F");                                                    
  TreeAllHad->Branch("Jet0SDptCorrL23Up"                     , & Jet0SDptCorrL23Up                  ,    "Jet0SDptCorrL23Up/F");                                                      
  TreeAllHad->Branch("Jet0SDptCorrL23Dn"                     , & Jet0SDptCorrL23Dn                  ,    "Jet0SDptCorrL23Dn/F");                                                      
  TreeAllHad->Branch("Jet0SDptCorrL123"                      , & Jet0SDptCorrL123                   ,    "Jet0SDptCorrL123/F");                                                      
  TreeAllHad->Branch("Jet0SDptCorrL123Up"                    , & Jet0SDptCorrL123Up                 ,    "Jet0SDptCorrL123Up/F");                                                        
  TreeAllHad->Branch("Jet0SDptCorrL123Dn"                    , & Jet0SDptCorrL123Dn                 ,    "Jet0SDptCorrL123Dn/F");                                                        
  TreeAllHad->Branch("Jet0SDptCorrL23Smear"                  , & Jet0SDptCorrL23Smear               ,    "Jet0SDptCorrL23Smear/F");                                                     
  TreeAllHad->Branch("Jet0SDptCorrL23SmearUp"                , & Jet0SDptCorrL23SmearUp             ,    "Jet0SDptCorrL23SmearUp/F");                                                       
  TreeAllHad->Branch("Jet0SDptCorrL23SmearDn"                , & Jet0SDptCorrL23SmearDn             ,    "Jet0SDptCorrL23SmearDn/F");                                                     
  TreeAllHad->Branch("Jet0SDetaRaw"                          , & Jet0SDetaRaw                       ,    "Jet0SDetaRaw/F");                                               
  TreeAllHad->Branch("Jet0SDphiRaw"                          , & Jet0SDphiRaw                       ,    "Jet0SDphiRaw/F");                                               
  TreeAllHad->Branch("Jet0MassPruned"                        , & Jet0MassPruned                     ,    "Jet0MassPruned/F");                                       
  TreeAllHad->Branch("Jet0MassTrimmed"                       , & Jet0MassTrimmed                    ,    "Jet0MassTrimmed/F");                                       
  TreeAllHad->Branch("Jet0Tau1"                              , & Jet0Tau1                           ,    "Jet0Tau1/F");                                 
  TreeAllHad->Branch("Jet0Tau2"                              , & Jet0Tau2                           ,    "Jet0Tau2/F");                                 
  TreeAllHad->Branch("Jet0Tau3"                              , & Jet0Tau3                           ,    "Jet0Tau3/F");                                 
  TreeAllHad->Branch("Jet0Tau4"                              , & Jet0Tau4                           ,    "Jet0Tau4/F");                                 
  TreeAllHad->Branch("Jet0Tau32"                             , & Jet0Tau32                          ,    "Jet0Tau32/F");                                  
  TreeAllHad->Branch("Jet0Tau21"                             , & Jet0Tau21                          ,    "Jet0Tau21/F");                                                                      
  TreeAllHad->Branch("Jet0SDmaxbdisc"                        , & Jet0SDmaxbdisc                     ,    "Jet0SDmaxbdisc/F");                                       
  TreeAllHad->Branch("Jet0SDmaxbdiscflavHadron"              , & Jet0SDmaxbdiscflavHadron           ,    "Jet0SDmaxbdiscflavHadron/F");                                           
  TreeAllHad->Branch("Jet0SDmaxbdiscflavParton"              , & Jet0SDmaxbdiscflavParton           ,    "Jet0SDmaxbdiscflavParton/F");                                           
  TreeAllHad->Branch("Jet0SDsubjet0pt"                       , & Jet0SDsubjet0pt                    ,    "Jet0SDsubjet0pt/F");                                        
  TreeAllHad->Branch("Jet0SDsubjet0mass"                     , & Jet0SDsubjet0mass                  ,    "Jet0SDsubjet0mass/F"); 
  TreeAllHad->Branch("Jet0SDsubjet0eta"                      , & Jet0SDsubjet0eta                   ,    "Jet0SDsubjet0eta/F");
  TreeAllHad->Branch("Jet0SDsubjet0phi"                      , & Jet0SDsubjet0phi                   ,    "Jet0SDsubjet0phi/F");                                         
  TreeAllHad->Branch("Jet0SDsubjet0area"                     , & Jet0SDsubjet0area                  ,    "Jet0SDsubjet0area/F");                                          
  TreeAllHad->Branch("Jet0SDsubjet0flavHadron"               , & Jet0SDsubjet0flavHadron            ,    "Jet0SDsubjet0flavHadron/F");                                          
  TreeAllHad->Branch("Jet0SDsubjet0flavParton"               , & Jet0SDsubjet0flavParton            ,    "Jet0SDsubjet0flavParton/F"); 
  TreeAllHad->Branch("Jet0SDsubjet0tau1"                     , & Jet0SDsubjet0tau1                  ,    "Jet0SDsubjet0tau1/F");
  TreeAllHad->Branch("Jet0SDsubjet0tau2"                     , & Jet0SDsubjet0tau2                  ,    "Jet0SDsubjet0tau2/F");
  TreeAllHad->Branch("Jet0SDsubjet0tau3"                     , & Jet0SDsubjet0tau3                  ,    "Jet0SDsubjet0tau3/F"); 
  TreeAllHad->Branch("Jet0SDsubjet0bdisc"                    , & Jet0SDsubjet0bdisc                 ,    "Jet0SDsubjet0bdisc/F");                                     
  TreeAllHad->Branch("Jet0SDsubjet1pt"                       , & Jet0SDsubjet1pt                    ,    "Jet0SDsubjet1pt/F");                                        
  TreeAllHad->Branch("Jet0SDsubjet1mass"                     , & Jet0SDsubjet1mass                  ,    "Jet0SDsubjet1mass/F");  
  TreeAllHad->Branch("Jet0SDsubjet1eta"                      , & Jet0SDsubjet1eta                   ,    "Jet0SDsubjet1eta/F");
  TreeAllHad->Branch("Jet0SDsubjet1phi"                      , & Jet0SDsubjet1phi                   ,    "Jet0SDsubjet1phi/F");                                        
  TreeAllHad->Branch("Jet0SDsubjet1area"                     , & Jet0SDsubjet1area                  ,    "Jet0SDsubjet1area/F");                                          
  TreeAllHad->Branch("Jet0SDsubjet1flavHadron"               , & Jet0SDsubjet1flavHadron            ,    "Jet0SDsubjet1flavHadron/F");                                          
  TreeAllHad->Branch("Jet0SDsubjet1flavParton"               , & Jet0SDsubjet1flavParton            ,    "Jet0SDsubjet1flavParton/F"); 
  TreeAllHad->Branch("Jet0SDsubjet1tau1"                     , & Jet0SDsubjet1tau1                  ,    "Jet0SDsubjet1tau1/F");
  TreeAllHad->Branch("Jet0SDsubjet1tau2"                     , & Jet0SDsubjet1tau2                  ,    "Jet0SDsubjet1tau2/F");
  TreeAllHad->Branch("Jet0SDsubjet1tau3"                     , & Jet0SDsubjet1tau3                  ,    "Jet0SDsubjet1tau3/F");  
  TreeAllHad->Branch("Jet0SDsubjet1bdisc"                    , & Jet0SDsubjet1bdisc                 ,    "Jet0SDsubjet1bdisc/F");                                                                                    
  TreeAllHad->Branch("Jet0PuppiP"                            , & Jet0PuppiP                         ,    "Jet0PuppiP/F");                                    
  TreeAllHad->Branch("Jet0PuppiPt"                           , & Jet0PuppiPt                        ,    "Jet0PuppiPt/F");                                    
  TreeAllHad->Branch("Jet0PuppiEta"                          , & Jet0PuppiEta                       ,    "Jet0PuppiEta/F");                                     
  TreeAllHad->Branch("Jet0PuppiPhi"                          , & Jet0PuppiPhi                       ,    "Jet0PuppiPhi/F");                                     
  TreeAllHad->Branch("Jet0PuppiMass"                         , & Jet0PuppiMass                      ,    "Jet0PuppiMass/F");                                      
  
  TreeAllHad->Branch("Jet0PuppiSDmass"                       , & Jet0PuppiSDmass                    ,   "Jet0PuppiSDmass/F");
  TreeAllHad->Branch("Jet0PuppiSDmassCorr"                   , & Jet0PuppiSDmassCorr                ,   "Jet0PuppiSDmassCorr/F");
  TreeAllHad->Branch("Jet0PuppiSDmassCorrUp"                 , & Jet0PuppiSDmassCorrUp              ,   "Jet0PuppiSDmassCorrUp/F");
  TreeAllHad->Branch("Jet0PuppiSDmassCorrDn"                 , & Jet0PuppiSDmassCorrDn              ,   "Jet0PuppiSDmassCorrDn/F");
  TreeAllHad->Branch("Jet0PuppiSDmassCorrL23Smear"           , & Jet0PuppiSDmassCorrL23Smear        ,   "Jet0PuppiSDmassCorrL23Smear/F");
  TreeAllHad->Branch("Jet0PuppiSDmassCorrL23SmearUp"         , & Jet0PuppiSDmassCorrL23SmearUp      ,   "Jet0PuppiSDmassCorrL23SmearUp/F");
  TreeAllHad->Branch("Jet0PuppiSDmassCorrL23SmearDn"         , & Jet0PuppiSDmassCorrL23SmearDn      ,   "Jet0PuppiSDmassCorrL23SmearDn/F");
  TreeAllHad->Branch("Jet0PuppiSDpt"                         , & Jet0PuppiSDpt                      ,   "Jet0PuppiSDpt/F");
  TreeAllHad->Branch("Jet0PuppiSDptCorr"                     , & Jet0PuppiSDptCorr                  ,   "Jet0PuppiSDptCorr/F");
  TreeAllHad->Branch("Jet0PuppiSDptCorrUp"                   , & Jet0PuppiSDptCorrUp                ,   "Jet0PuppiSDptCorrUp/F");
  TreeAllHad->Branch("Jet0PuppiSDptCorrDn"                   , & Jet0PuppiSDptCorrDn                ,   "Jet0PuppiSDptCorrDn/F");
  TreeAllHad->Branch("Jet0PuppiSDptCorrL23Smear"             , & Jet0PuppiSDptCorrL23Smear          ,   "Jet0PuppiSDptCorrL23Smear/F");
  TreeAllHad->Branch("Jet0PuppiSDptCorrL23SmearUp"           , & Jet0PuppiSDptCorrL23SmearUp        ,   "Jet0PuppiSDptCorrL23SmearUp/F");
  TreeAllHad->Branch("Jet0PuppiSDptCorrL23SmearDn"           , & Jet0PuppiSDptCorrL23SmearDn        ,   "Jet0PuppiSDptCorrL23SmearDn/F");
  TreeAllHad->Branch("Jet0PuppiSDeta"                        , & Jet0PuppiSDeta                     ,   "Jet0PuppiSDeta/F");
  TreeAllHad->Branch("Jet0PuppiSDphi"                        , & Jet0PuppiSDphi                     ,   "Jet0PuppiSDphi/F");
                                                           
  TreeAllHad->Branch("Jet0PuppiTau1"                         , & Jet0PuppiTau1                      ,    "Jet0PuppiTau1/F");                                      
  TreeAllHad->Branch("Jet0PuppiTau2"                         , & Jet0PuppiTau2                      ,    "Jet0PuppiTau2/F");                                      
  TreeAllHad->Branch("Jet0PuppiTau3"                         , & Jet0PuppiTau3                      ,    "Jet0PuppiTau3/F");                                      
  TreeAllHad->Branch("Jet0PuppiTau4"                         , & Jet0PuppiTau4                      ,    "Jet0PuppiTau4/F");                                      
  TreeAllHad->Branch("Jet0PuppiTau32"                        , & Jet0PuppiTau32                     ,    "Jet0PuppiTau32/F");                                       
  TreeAllHad->Branch("Jet0PuppiTau21"                        , & Jet0PuppiTau21                     ,    "Jet0PuppiTau21/F");                                                                                
  TreeAllHad->Branch("Jet0PuppiSDmaxbdisc"                   , & Jet0PuppiSDmaxbdisc                ,    "Jet0PuppiSDmaxbdisc/F");                                            
  TreeAllHad->Branch("Jet0PuppiSDmaxbdiscflavHadron"         , & Jet0PuppiSDmaxbdiscflavHadron      ,    "Jet0PuppiSDmaxbdiscflavHadron/F");                                                
  TreeAllHad->Branch("Jet0PuppiSDmaxbdiscflavParton"         , & Jet0PuppiSDmaxbdiscflavParton      ,    "Jet0PuppiSDmaxbdiscflavParton/F");                                                
  TreeAllHad->Branch("Jet0PuppiSDsubjet0pt"                  , & Jet0PuppiSDsubjet0pt               ,    "Jet0PuppiSDsubjet0pt/F");                                             
  TreeAllHad->Branch("Jet0PuppiSDsubjet0mass"                , & Jet0PuppiSDsubjet0mass             ,    "Jet0PuppiSDsubjet0mass/F");    
  TreeAllHad->Branch("Jet0PuppiSDsubjet0eta"                 , & Jet0PuppiSDsubjet0eta              ,    "Jet0PuppiSDsubjet0eta/F");
  TreeAllHad->Branch("Jet0PuppiSDsubjet0phi"                 , & Jet0PuppiSDsubjet0phi              ,    "Jet0PuppiSDsubjet0phi/F");                                           
  TreeAllHad->Branch("Jet0PuppiSDsubjet0area"                , & Jet0PuppiSDsubjet0area             ,    "Jet0PuppiSDsubjet0area/F");                                               
  TreeAllHad->Branch("Jet0PuppiSDsubjet0flavHadron"          , & Jet0PuppiSDsubjet0flavHadron       ,    "Jet0PuppiSDsubjet0flavHadron/F");                                               
  TreeAllHad->Branch("Jet0PuppiSDsubjet0flavParton"          , & Jet0PuppiSDsubjet0flavParton       ,    "Jet0PuppiSDsubjet0flavParton/F"); 
  TreeAllHad->Branch("Jet0PuppiSDsubjet0tau1"                , & Jet0PuppiSDsubjet0tau1             ,    "Jet0PuppiSDsubjet0tau1/F");
  TreeAllHad->Branch("Jet0PuppiSDsubjet0tau2"                , & Jet0PuppiSDsubjet0tau2             ,    "Jet0PuppiSDsubjet0tau2/F");
  TreeAllHad->Branch("Jet0PuppiSDsubjet0tau3"                , & Jet0PuppiSDsubjet0tau3             ,    "Jet0PuppiSDsubjet0tau3/F"); 
  TreeAllHad->Branch("Jet0PuppiSDsubjet0bdisc"               , & Jet0PuppiSDsubjet0bdisc            ,    "Jet0PuppiSDsubjet0bdisc/F");                                          
  TreeAllHad->Branch("Jet0PuppiSDsubjet1pt"                  , & Jet0PuppiSDsubjet1pt               ,    "Jet0PuppiSDsubjet1pt/F");                                             
  TreeAllHad->Branch("Jet0PuppiSDsubjet1mass"                , & Jet0PuppiSDsubjet1mass             ,    "Jet0PuppiSDsubjet1mass/F");  
  TreeAllHad->Branch("Jet0PuppiSDsubjet1eta"                 , & Jet0PuppiSDsubjet1eta              ,    "Jet0PuppiSDsubjet1eta/F");
  TreeAllHad->Branch("Jet0PuppiSDsubjet1phi"                 , & Jet0PuppiSDsubjet1phi              ,    "Jet0PuppiSDsubjet1phi/F");                                             
  TreeAllHad->Branch("Jet0PuppiSDsubjet1area"                , & Jet0PuppiSDsubjet1area             ,    "Jet0PuppiSDsubjet1area/F");                                               
  TreeAllHad->Branch("Jet0PuppiSDsubjet1flavHadron"          , & Jet0PuppiSDsubjet1flavHadron       ,    "Jet0PuppiSDsubjet1flavHadron/F");                                               
  TreeAllHad->Branch("Jet0PuppiSDsubjet1flavParton"          , & Jet0PuppiSDsubjet1flavParton       ,    "Jet0PuppiSDsubjet1flavParton/F");    
  TreeAllHad->Branch("Jet0PuppiSDsubjet1tau1"                , & Jet0PuppiSDsubjet1tau1             ,    "Jet0PuppiSDsubjet1tau1/F");
  TreeAllHad->Branch("Jet0PuppiSDsubjet1tau2"                , & Jet0PuppiSDsubjet1tau2             ,    "Jet0PuppiSDsubjet1tau2/F");
  TreeAllHad->Branch("Jet0PuppiSDsubjet1tau3"                , & Jet0PuppiSDsubjet1tau3             ,    "Jet0PuppiSDsubjet1tau3/F"); 
  TreeAllHad->Branch("Jet0PuppiSDsubjet1bdisc"               , & Jet0PuppiSDsubjet1bdisc            ,    "Jet0PuppiSDsubjet1bdisc/F");                                                                                          
  TreeAllHad->Branch("Jet0CHF"                               , & Jet0CHF                            ,    "Jet0CHF/F");                                
  TreeAllHad->Branch("Jet0NHF"                               , & Jet0NHF                            ,    "Jet0NHF/F");                                
  TreeAllHad->Branch("Jet0CM"                                , & Jet0CM                             ,    "Jet0CM/F");                               
  TreeAllHad->Branch("Jet0NM"                                , & Jet0NM                             ,    "Jet0NM/F");                               
  TreeAllHad->Branch("Jet0NEF"                               , & Jet0NEF                            ,    "Jet0NEF/F");                                
  TreeAllHad->Branch("Jet0CEF"                               , & Jet0CEF                            ,    "Jet0CEF/F");                                
  TreeAllHad->Branch("Jet0MF"                                , & Jet0MF                             ,    "Jet0MF/F");                               
  TreeAllHad->Branch("Jet0Mult"                              , & Jet0Mult                           ,    "Jet0Mult/F");                                 


  TreeAllHad->Branch("Jet0PuppiCHF"                         , & Jet0PuppiCHF                        ,    "Jet0PuppiCHF/F");                                
  TreeAllHad->Branch("Jet0PuppiNHF"                         , & Jet0PuppiNHF                        ,    "Jet0PuppiNHF/F");                                
  TreeAllHad->Branch("Jet0PuppiCM"                          , & Jet0PuppiCM                         ,    "Jet0PuppiCM/F");                               
  TreeAllHad->Branch("Jet0PuppiNM"                          , & Jet0PuppiNM                         ,    "Jet0PuppiNM/F");                               
  TreeAllHad->Branch("Jet0PuppiNEF"                         , & Jet0PuppiNEF                        ,    "Jet0PuppiNEF/F");                                
  TreeAllHad->Branch("Jet0PuppiCEF"                         , & Jet0PuppiCEF                        ,    "Jet0PuppiCEF/F");                                
  TreeAllHad->Branch("Jet0PuppiMF"                          , & Jet0PuppiMF                         ,    "Jet0PuppiMF/F");                               
  TreeAllHad->Branch("Jet0PuppiMult"                        , & Jet0PuppiMult                       ,    "Jet0PuppiMult/F");  

  TreeAllHad->Branch("Jet0MassCorrFactor"                    , & Jet0MassCorrFactor                 ,    "Jet0MassCorrFactor/F");                                           
  TreeAllHad->Branch("Jet0MassCorrFactorUp"                  , & Jet0MassCorrFactorUp               ,    "Jet0MassCorrFactorUp/F");                                             
  TreeAllHad->Branch("Jet0MassCorrFactorDn"                  , & Jet0MassCorrFactorDn               ,    "Jet0MassCorrFactorDn/F");                                             
  TreeAllHad->Branch("Jet0CorrFactor"                        , & Jet0CorrFactor                     ,    "Jet0CorrFactor/F");                                       
  TreeAllHad->Branch("Jet0CorrFactorUp"                      , & Jet0CorrFactorUp                   ,    "Jet0CorrFactorUp/F");                                         
  TreeAllHad->Branch("Jet0CorrFactorDn"                      , & Jet0CorrFactorDn                   ,    "Jet0CorrFactorDn/F");                                         
  TreeAllHad->Branch("Jet0PtSmearFactor"                     , & Jet0PtSmearFactor                  ,    "Jet0PtSmearFactor/F");                                          
  TreeAllHad->Branch("Jet0PtSmearFactorUp"                   , & Jet0PtSmearFactorUp                ,    "Jet0PtSmearFactorUp/F");                                            
  TreeAllHad->Branch("Jet0PtSmearFactorDn"                   , & Jet0PtSmearFactorDn                ,    "Jet0PtSmearFactorDn/F");                                            
  TreeAllHad->Branch("Jet0PuppiMassCorrFactor"               , & Jet0PuppiMassCorrFactor            ,    "Jet0PuppiMassCorrFactor/F");                                                
  TreeAllHad->Branch("Jet0PuppiMassCorrFactorUp"             , & Jet0PuppiMassCorrFactorUp          ,    "Jet0PuppiMassCorrFactorUp/F");                                                  
  TreeAllHad->Branch("Jet0PuppiMassCorrFactorDn"             , & Jet0PuppiMassCorrFactorDn          ,    "Jet0PuppiMassCorrFactorDn/F");                                                  
  TreeAllHad->Branch("Jet0PuppiCorrFactor"                   , & Jet0PuppiCorrFactor                ,    "Jet0PuppiCorrFactor/F");                                            
  TreeAllHad->Branch("Jet0PuppiCorrFactorUp"                 , & Jet0PuppiCorrFactorUp              ,    "Jet0PuppiCorrFactorUp/F");                                              
  TreeAllHad->Branch("Jet0PuppiCorrFactorDn"                 , & Jet0PuppiCorrFactorDn              ,    "Jet0PuppiCorrFactorDn/F");                                              
  TreeAllHad->Branch("Jet0PuppiPtSmearFactor"                , & Jet0PuppiPtSmearFactor             ,    "Jet0PuppiPtSmearFactor/F");                                               
  TreeAllHad->Branch("Jet0PuppiPtSmearFactorUp"              , & Jet0PuppiPtSmearFactorUp           ,    "Jet0PuppiPtSmearFactorUp/F");                                                 
  TreeAllHad->Branch("Jet0PuppiPtSmearFactorDn"              , & Jet0PuppiPtSmearFactorDn           ,    "Jet0PuppiPtSmearFactorDn/F");                                                 
  TreeAllHad->Branch("Jet0EtaScaleFactor"                    , & Jet0EtaScaleFactor                 ,    "Jet0EtaScaleFactor/F");                                           
  TreeAllHad->Branch("Jet0PhiScaleFactor"                    , & Jet0PhiScaleFactor                 ,    "Jet0PhiScaleFactor/F");                                           
  // TreeAllHad->Branch("Jet0MatchedGenJetDR"                   , & Jet0MatchedGenJetDR                ,    "Jet0MatchedGenJetDR/F");                                            
  TreeAllHad->Branch("Jet0MatchedGenJetPt"                   , & Jet0MatchedGenJetPt                ,    "Jet0MatchedGenJetPt/F");                                            
  TreeAllHad->Branch("Jet0MatchedGenJetMass"                 , & Jet0MatchedGenJetMass              ,    "Jet0MatchedGenJetMass/F");  

  TreeAllHad->Branch("Jet0GenMatched_TopHadronic"            , & Jet0GenMatched_TopHadronic         ,    "Jet0GenMatched_TopHadronic/I");      
  TreeAllHad->Branch("Jet0GenMatched_TopPt"                  , & Jet0GenMatched_TopPt               ,    "Jet0GenMatched_TopPt/F");      
  TreeAllHad->Branch("Jet0GenMatched_TopEta"                 , & Jet0GenMatched_TopEta              ,    "Jet0GenMatched_TopEta/F");      
  TreeAllHad->Branch("Jet0GenMatched_TopPhi"                 , & Jet0GenMatched_TopPhi              ,    "Jet0GenMatched_TopPhi/F");      
  TreeAllHad->Branch("Jet0GenMatched_TopMass"                , & Jet0GenMatched_TopMass             ,    "Jet0GenMatched_TopMass/F");      
  TreeAllHad->Branch("Jet0GenMatched_bPt"                    , & Jet0GenMatched_bPt                 ,    "Jet0GenMatched_bPt/F");      
  TreeAllHad->Branch("Jet0GenMatched_WPt"                    , & Jet0GenMatched_WPt                 ,    "Jet0GenMatched_WPt/F");      
  TreeAllHad->Branch("Jet0GenMatched_Wd1Pt"                  , & Jet0GenMatched_Wd1Pt               ,    "Jet0GenMatched_Wd1Pt/F");      
  TreeAllHad->Branch("Jet0GenMatched_Wd2Pt"                  , & Jet0GenMatched_Wd2Pt               ,    "Jet0GenMatched_Wd2Pt/F");      
  TreeAllHad->Branch("Jet0GenMatched_Wd1ID"                  , & Jet0GenMatched_Wd1ID               ,    "Jet0GenMatched_Wd1ID/F");      
  TreeAllHad->Branch("Jet0GenMatched_Wd2ID"                  , & Jet0GenMatched_Wd2ID               ,    "Jet0GenMatched_Wd2ID/F");      
  TreeAllHad->Branch("Jet0GenMatched_MaxDeltaRPartonTop"     , & Jet0GenMatched_MaxDeltaRPartonTop  ,    "Jet0GenMatched_MaxDeltaRPartonTop/F");      
  TreeAllHad->Branch("Jet0GenMatched_MaxDeltaRWPartonTop"    , & Jet0GenMatched_MaxDeltaRWPartonTop ,    "Jet0GenMatched_MaxDeltaRWPartonTop/F");      
  TreeAllHad->Branch("Jet0GenMatched_MaxDeltaRWPartonW"      , & Jet0GenMatched_MaxDeltaRWPartonW   ,    "Jet0GenMatched_MaxDeltaRWPartonW/F");      
  TreeAllHad->Branch("Jet0GenMatched_DeltaR_t_b"             , & Jet0GenMatched_DeltaR_t_b          ,    "Jet0GenMatched_DeltaR_t_b/F");      
  TreeAllHad->Branch("Jet0GenMatched_DeltaR_t_W"             , & Jet0GenMatched_DeltaR_t_W          ,    "Jet0GenMatched_DeltaR_t_W/F");      
  TreeAllHad->Branch("Jet0GenMatched_DeltaR_t_Wd1"           , & Jet0GenMatched_DeltaR_t_Wd1        ,    "Jet0GenMatched_DeltaR_t_Wd1/F");      
  TreeAllHad->Branch("Jet0GenMatched_DeltaR_t_Wd2"           , & Jet0GenMatched_DeltaR_t_Wd2        ,    "Jet0GenMatched_DeltaR_t_Wd2/F");      
  TreeAllHad->Branch("Jet0GenMatched_DeltaR_W_b1"            , & Jet0GenMatched_DeltaR_W_b1         ,    "Jet0GenMatched_DeltaR_W_b1/F");      
  TreeAllHad->Branch("Jet0GenMatched_DeltaR_W_Wd1"           , & Jet0GenMatched_DeltaR_W_Wd1        ,    "Jet0GenMatched_DeltaR_W_Wd1/F");      
  TreeAllHad->Branch("Jet0GenMatched_DeltaR_W_Wd2"           , & Jet0GenMatched_DeltaR_W_Wd2        ,    "Jet0GenMatched_DeltaR_W_Wd2/F");      
  TreeAllHad->Branch("Jet0GenMatched_DeltaR_Wd1_Wd2"         , & Jet0GenMatched_DeltaR_Wd1_Wd2      ,    "Jet0GenMatched_DeltaR_Wd1_Wd2/F");      
  TreeAllHad->Branch("Jet0GenMatched_DeltaR_Wd1_b"           , & Jet0GenMatched_DeltaR_Wd1_b        ,    "Jet0GenMatched_DeltaR_Wd1_b/F");      
  TreeAllHad->Branch("Jet0GenMatched_DeltaR_Wd2_b"           , & Jet0GenMatched_DeltaR_Wd2_b        ,    "Jet0GenMatched_DeltaR_Wd2_b/F");      
  TreeAllHad->Branch("Jet0GenMatched_DeltaR_jet_t"           , & Jet0GenMatched_DeltaR_jet_t        ,    "Jet0GenMatched_DeltaR_jet_t/F");      
  TreeAllHad->Branch("Jet0GenMatched_DeltaR_jet_W"           , & Jet0GenMatched_DeltaR_jet_W        ,    "Jet0GenMatched_DeltaR_jet_W/F");      
  TreeAllHad->Branch("Jet0GenMatched_DeltaR_jet_b"           , & Jet0GenMatched_DeltaR_jet_b        ,    "Jet0GenMatched_DeltaR_jet_b/F");      
  TreeAllHad->Branch("Jet0GenMatched_DeltaR_jet_Wd1"         , & Jet0GenMatched_DeltaR_jet_Wd1      ,    "Jet0GenMatched_DeltaR_jet_Wd1/F");      
  TreeAllHad->Branch("Jet0GenMatched_DeltaR_jet_Wd2"         , & Jet0GenMatched_DeltaR_jet_Wd2      ,    "Jet0GenMatched_DeltaR_jet_Wd2/F");      
  TreeAllHad->Branch("Jet0GenMatched_DeltaR_pup0_b"          , & Jet0GenMatched_DeltaR_pup0_b       ,    "Jet0GenMatched_DeltaR_pup0_b/F");      
  TreeAllHad->Branch("Jet0GenMatched_DeltaR_pup0_Wd1"        , & Jet0GenMatched_DeltaR_pup0_Wd1     ,    "Jet0GenMatched_DeltaR_pup0_Wd1/F");      
  TreeAllHad->Branch("Jet0GenMatched_DeltaR_pup0_Wd2"        , & Jet0GenMatched_DeltaR_pup0_Wd2     ,    "Jet0GenMatched_DeltaR_pup0_Wd2/F");      
  TreeAllHad->Branch("Jet0GenMatched_DeltaR_pup1_b"          , & Jet0GenMatched_DeltaR_pup1_b       ,    "Jet0GenMatched_DeltaR_pup1_b/F");      
  TreeAllHad->Branch("Jet0GenMatched_DeltaR_pup1_Wd1"        , & Jet0GenMatched_DeltaR_pup1_Wd1     ,    "Jet0GenMatched_DeltaR_pup1_Wd1/F");      
  TreeAllHad->Branch("Jet0GenMatched_DeltaR_pup1_Wd2"        , & Jet0GenMatched_DeltaR_pup1_Wd2     ,    "Jet0GenMatched_DeltaR_pup1_Wd2/F");               
  TreeAllHad->Branch("Jet0GenMatched_partonPt"               , & Jet0GenMatched_partonPt            ,    "Jet0GenMatched_partonPt/F");      
  TreeAllHad->Branch("Jet0GenMatched_partonEta"              , & Jet0GenMatched_partonEta           ,    "Jet0GenMatched_partonEta/F");      
  TreeAllHad->Branch("Jet0GenMatched_partonPhi"              , & Jet0GenMatched_partonPhi           ,    "Jet0GenMatched_partonPhi/F");      
  TreeAllHad->Branch("Jet0GenMatched_partonMass"             , & Jet0GenMatched_partonMass          ,    "Jet0GenMatched_partonMass/F");      
  TreeAllHad->Branch("Jet0GenMatched_partonID"               , & Jet0GenMatched_partonID            ,    "Jet0GenMatched_partonID/F");      
  TreeAllHad->Branch("Jet0GenMatched_DeltaRjetParton"        , & Jet0GenMatched_DeltaRjetParton     ,    "Jet0GenMatched_DeltaRjetParton/F");      
  std::cout << "Setup Jet0 in tree" << std::endl;
  
  TreeAllHad->Branch("Jet1PtRaw"                             , & Jet1PtRaw                          ,    "Jet1PtRaw/F");                                  
  TreeAllHad->Branch("Jet1EtaRaw"                            , & Jet1EtaRaw                         ,    "Jet1EtaRaw/F");                                   
  TreeAllHad->Branch("Jet1PhiRaw"                            , & Jet1PhiRaw                         ,    "Jet1PhiRaw/F");                                   
  TreeAllHad->Branch("Jet1MassRaw"                           , & Jet1MassRaw                        ,    "Jet1MassRaw/F");                                    
  TreeAllHad->Branch("Jet1P"                                 , & Jet1P                              ,    "Jet1P/F");                              
  TreeAllHad->Branch("Jet1Pt"                                , & Jet1Pt                             ,    "Jet1Pt/F");                               
  TreeAllHad->Branch("Jet1Eta"                               , & Jet1Eta                            ,    "Jet1Eta/F");                                
  TreeAllHad->Branch("Jet1Phi"                               , & Jet1Phi                            ,    "Jet1Phi/F");                                
  TreeAllHad->Branch("Jet1Rap"                               , & Jet1Rap                            ,    "Jet1Rap/F");                                
  TreeAllHad->Branch("Jet1Energy"                            , & Jet1Energy                         ,    "Jet1Energy/F");                                   
  TreeAllHad->Branch("Jet1Mass"                              , & Jet1Mass                           ,    "Jet1Mass/F");                                 
  TreeAllHad->Branch("Jet1Area"                              , & Jet1Area                           ,    "Jet1Area/F");                                 
  TreeAllHad->Branch("Jet1SDmass"                            , & Jet1SDmass                         ,    "Jet1SDmass/F");                                         
  TreeAllHad->Branch("Jet1SDmassRaw"                         , & Jet1SDmassRaw                      ,    "Jet1SDmassRaw/F");                                               
  TreeAllHad->Branch("Jet1SDmassCorrL23"                     , & Jet1SDmassCorrL23                  ,    "Jet1SDmassCorrL23/F");                                                    
  TreeAllHad->Branch("Jet1SDmassCorrL23Up"                   , & Jet1SDmassCorrL23Up                ,    "Jet1SDmassCorrL23Up/F");                                                      
  TreeAllHad->Branch("Jet1SDmassCorrL23Dn"                   , & Jet1SDmassCorrL23Dn                ,    "Jet1SDmassCorrL23Dn/F");                                                      
  TreeAllHad->Branch("Jet1SDmassCorrL123"                    , & Jet1SDmassCorrL123                 ,    "Jet1SDmassCorrL123/F");                                                      
  TreeAllHad->Branch("Jet1SDmassCorrL123Up"                  , & Jet1SDmassCorrL123Up               ,    "Jet1SDmassCorrL123Up/F");                                                        
  TreeAllHad->Branch("Jet1SDmassCorrL123Dn"                  , & Jet1SDmassCorrL123Dn               ,    "Jet1SDmassCorrL123Dn/F");                                                        
  TreeAllHad->Branch("Jet1SDmassCorrL23Smear"                   , & Jet1SDmassCorrL23Smear                ,    "Jet1SDmassCorrL23Smear/F");                                                     
  TreeAllHad->Branch("Jet1SDmassCorrL23SmearUp"                 , & Jet1SDmassCorrL23SmearUp              ,    "Jet1SDmassCorrL23SmearUp/F");                                                       
  TreeAllHad->Branch("Jet1SDmassCorrL23SmearDn"                 , & Jet1SDmassCorrL23SmearDn              ,    "Jet1SDmassCorrL23SmearDn/F");   
  TreeAllHad->Branch("Jet1SDptRaw"                           , & Jet1SDptRaw                        ,    "Jet1SDptRaw/F");                                               
  TreeAllHad->Branch("Jet1SDptCorrL23"                       , & Jet1SDptCorrL23                    ,    "Jet1SDptCorrL23/F");                                                    
  TreeAllHad->Branch("Jet1SDptCorrL23Up"                     , & Jet1SDptCorrL23Up                  ,    "Jet1SDptCorrL23Up/F");                                                      
  TreeAllHad->Branch("Jet1SDptCorrL23Dn"                     , & Jet1SDptCorrL23Dn                  ,    "Jet1SDptCorrL23Dn/F");                                                      
  TreeAllHad->Branch("Jet1SDptCorrL123"                      , & Jet1SDptCorrL123                   ,    "Jet1SDptCorrL123/F");                                                      
  TreeAllHad->Branch("Jet1SDptCorrL123Up"                    , & Jet1SDptCorrL123Up                 ,    "Jet1SDptCorrL123Up/F");                                                        
  TreeAllHad->Branch("Jet1SDptCorrL123Dn"                    , & Jet1SDptCorrL123Dn                 ,    "Jet1SDptCorrL123Dn/F");                                                        
  TreeAllHad->Branch("Jet1SDptCorrL23Smear"                     , & Jet1SDptCorrL23Smear                  ,    "Jet1SDptCorrL23Smear/F");                                                     
  TreeAllHad->Branch("Jet1SDptCorrL23SmearUp"                   , & Jet1SDptCorrL23SmearUp                ,    "Jet1SDptCorrL23SmearUp/F");                                                       
  TreeAllHad->Branch("Jet1SDptCorrL23SmearDn"                   , & Jet1SDptCorrL23SmearDn                ,    "Jet1SDptCorrL23SmearDn/F");                                                     
  TreeAllHad->Branch("Jet1SDetaRaw"                          , & Jet1SDetaRaw                       ,    "Jet1SDetaRaw/F");                                               
  TreeAllHad->Branch("Jet1SDphiRaw"                          , & Jet1SDphiRaw                       ,    "Jet1SDphiRaw/F");  

  TreeAllHad->Branch("Jet1MassPruned"                        , & Jet1MassPruned                     ,    "Jet1MassPruned/F");                                       
  TreeAllHad->Branch("Jet1MassTrimmed"                       , & Jet1MassTrimmed                    ,    "Jet1MassTrimmed/F");                                       
  TreeAllHad->Branch("Jet1Tau1"                              , & Jet1Tau1                           ,    "Jet1Tau1/F");                                 
  TreeAllHad->Branch("Jet1Tau2"                              , & Jet1Tau2                           ,    "Jet1Tau2/F");                                 
  TreeAllHad->Branch("Jet1Tau3"                              , & Jet1Tau3                           ,    "Jet1Tau3/F");                                 
  TreeAllHad->Branch("Jet1Tau4"                              , & Jet1Tau4                           ,    "Jet1Tau4/F");                                 
  TreeAllHad->Branch("Jet1Tau32"                             , & Jet1Tau32                          ,    "Jet1Tau32/F");                                  
  TreeAllHad->Branch("Jet1Tau21"                             , & Jet1Tau21                          ,    "Jet1Tau21/F");                                                                      
  TreeAllHad->Branch("Jet1SDmaxbdisc"                        , & Jet1SDmaxbdisc                     ,    "Jet1SDmaxbdisc/F");                                       
  TreeAllHad->Branch("Jet1SDmaxbdiscflavHadron"              , & Jet1SDmaxbdiscflavHadron           ,    "Jet1SDmaxbdiscflavHadron/F");                                           
  TreeAllHad->Branch("Jet1SDmaxbdiscflavParton"              , & Jet1SDmaxbdiscflavParton           ,    "Jet1SDmaxbdiscflavParton/F");                                           
  TreeAllHad->Branch("Jet1SDsubjet0pt"                       , & Jet1SDsubjet0pt                    ,    "Jet1SDsubjet0pt/F");                                        
  TreeAllHad->Branch("Jet1SDsubjet0mass"                     , & Jet1SDsubjet0mass                  ,    "Jet1SDsubjet0mass/F"); 
  TreeAllHad->Branch("Jet1SDsubjet0eta"                      , & Jet1SDsubjet0eta                   ,    "Jet1SDsubjet0eta/F");
  TreeAllHad->Branch("Jet1SDsubjet0phi"                      , & Jet1SDsubjet0phi                   ,    "Jet1SDsubjet0phi/F");                                         
  TreeAllHad->Branch("Jet1SDsubjet0area"                     , & Jet1SDsubjet0area                  ,    "Jet1SDsubjet0area/F");                                          
  TreeAllHad->Branch("Jet1SDsubjet0flavHadron"               , & Jet1SDsubjet0flavHadron            ,    "Jet1SDsubjet0flavHadron/F");                                          
  TreeAllHad->Branch("Jet1SDsubjet0flavParton"               , & Jet1SDsubjet0flavParton            ,    "Jet1SDsubjet0flavParton/F"); 
  TreeAllHad->Branch("Jet1SDsubjet0tau1"                     , & Jet1SDsubjet0tau1                  ,    "Jet1SDsubjet0tau1/F");
  TreeAllHad->Branch("Jet1SDsubjet0tau2"                     , & Jet1SDsubjet0tau2                  ,    "Jet1SDsubjet0tau2/F");
  TreeAllHad->Branch("Jet1SDsubjet0tau3"                     , & Jet1SDsubjet0tau3                  ,    "Jet1SDsubjet0tau3/F");
  TreeAllHad->Branch("Jet1SDsubjet0bdisc"                    , & Jet1SDsubjet0bdisc                 ,    "Jet1SDsubjet0bdisc/F");                                     
  TreeAllHad->Branch("Jet1SDsubjet1pt"                       , & Jet1SDsubjet1pt                    ,    "Jet1SDsubjet1pt/F");                                        
  TreeAllHad->Branch("Jet1SDsubjet1mass"                     , & Jet1SDsubjet1mass                  ,    "Jet1SDsubjet1mass/F");  
  TreeAllHad->Branch("Jet1SDsubjet1eta"                      , & Jet1SDsubjet1eta                   ,    "Jet1SDsubjet1eta/F");
  TreeAllHad->Branch("Jet1SDsubjet1phi"                      , & Jet1SDsubjet1phi                   ,    "Jet1SDsubjet1phi/F");                                        
  TreeAllHad->Branch("Jet1SDsubjet1area"                     , & Jet1SDsubjet1area                  ,    "Jet1SDsubjet1area/F");                                          
  TreeAllHad->Branch("Jet1SDsubjet1flavHadron"               , & Jet1SDsubjet1flavHadron            ,    "Jet1SDsubjet1flavHadron/F");                                          
  TreeAllHad->Branch("Jet1SDsubjet1flavParton"               , & Jet1SDsubjet1flavParton            ,    "Jet1SDsubjet1flavParton/F"); 
  TreeAllHad->Branch("Jet1SDsubjet1tau1"                     , & Jet1SDsubjet1tau1                  ,    "Jet1SDsubjet1tau1/F");
  TreeAllHad->Branch("Jet1SDsubjet1tau2"                     , & Jet1SDsubjet1tau2                  ,    "Jet1SDsubjet1tau2/F");
  TreeAllHad->Branch("Jet1SDsubjet1tau3"                     , & Jet1SDsubjet1tau3                  ,    "Jet1SDsubjet1tau3/F"); 
  TreeAllHad->Branch("Jet1SDsubjet1bdisc"                    , & Jet1SDsubjet1bdisc                 ,    "Jet1SDsubjet1bdisc/F");                                                                                    
  TreeAllHad->Branch("Jet1PuppiP"                            , & Jet1PuppiP                         ,    "Jet1PuppiP/F");                                    
  TreeAllHad->Branch("Jet1PuppiPt"                           , & Jet1PuppiPt                        ,    "Jet1PuppiPt/F");                                    
  TreeAllHad->Branch("Jet1PuppiEta"                          , & Jet1PuppiEta                       ,    "Jet1PuppiEta/F");                                     
  TreeAllHad->Branch("Jet1PuppiPhi"                          , & Jet1PuppiPhi                       ,    "Jet1PuppiPhi/F");                                     
  TreeAllHad->Branch("Jet1PuppiMass"                         , & Jet1PuppiMass                      ,    "Jet1PuppiMass/F");                                      
  
  TreeAllHad->Branch("Jet1PuppiSDmass"                       , & Jet1PuppiSDmass                    ,   "Jet1PuppiSDmass/F");
  TreeAllHad->Branch("Jet1PuppiSDmassCorr"                   , & Jet1PuppiSDmassCorr                ,   "Jet1PuppiSDmassCorr/F");
  TreeAllHad->Branch("Jet1PuppiSDmassCorrUp"                 , & Jet1PuppiSDmassCorrUp              ,   "Jet1PuppiSDmassCorrUp/F");
  TreeAllHad->Branch("Jet1PuppiSDmassCorrDn"                 , & Jet1PuppiSDmassCorrDn              ,   "Jet1PuppiSDmassCorrDn/F");
  TreeAllHad->Branch("Jet1PuppiSDmassCorrL23Smear"           , & Jet1PuppiSDmassCorrL23Smear        ,   "Jet1PuppiSDmassCorrL23Smear/F");
  TreeAllHad->Branch("Jet1PuppiSDmassCorrL23SmearUp"         , & Jet1PuppiSDmassCorrL23SmearUp      ,   "Jet1PuppiSDmassCorrL23SmearUp/F");
  TreeAllHad->Branch("Jet1PuppiSDmassCorrL23SmearDn"         , & Jet1PuppiSDmassCorrL23SmearDn      ,   "Jet1PuppiSDmassCorrL23SmearDn/F");
  TreeAllHad->Branch("Jet1PuppiSDpt"                         , & Jet1PuppiSDpt                      ,   "Jet1PuppiSDpt/F");
  TreeAllHad->Branch("Jet1PuppiSDptCorr"                     , & Jet1PuppiSDptCorr                  ,   "Jet1PuppiSDptCorr/F");
  TreeAllHad->Branch("Jet1PuppiSDptCorrUp"                   , & Jet1PuppiSDptCorrUp                ,   "Jet1PuppiSDptCorrUp/F");
  TreeAllHad->Branch("Jet1PuppiSDptCorrDn"                   , & Jet1PuppiSDptCorrDn                ,   "Jet1PuppiSDptCorrDn/F");
  TreeAllHad->Branch("Jet1PuppiSDptCorrL23Smear"             , & Jet1PuppiSDptCorrL23Smear          ,   "Jet1PuppiSDptCorrL23Smear/F");
  TreeAllHad->Branch("Jet1PuppiSDptCorrL23SmearUp"           , & Jet1PuppiSDptCorrL23SmearUp        ,   "Jet1PuppiSDptCorrL23SmearUp/F");
  TreeAllHad->Branch("Jet1PuppiSDptCorrL23SmearDn"           , & Jet1PuppiSDptCorrL23SmearDn        ,   "Jet1PuppiSDptCorrL23SmearDn/F");
  TreeAllHad->Branch("Jet1PuppiSDeta"                        , & Jet1PuppiSDeta                     ,   "Jet1PuppiSDeta/F");
  TreeAllHad->Branch("Jet1PuppiSDphi"                        , & Jet1PuppiSDphi                     ,   "Jet1PuppiSDphi/F");
                         
  TreeAllHad->Branch("Jet1PuppiTau1"                         , & Jet1PuppiTau1                      ,    "Jet1PuppiTau1/F");                                      
  TreeAllHad->Branch("Jet1PuppiTau2"                         , & Jet1PuppiTau2                      ,    "Jet1PuppiTau2/F");                                      
  TreeAllHad->Branch("Jet1PuppiTau3"                         , & Jet1PuppiTau3                      ,    "Jet1PuppiTau3/F");                                      
  TreeAllHad->Branch("Jet1PuppiTau4"                         , & Jet1PuppiTau4                      ,    "Jet1PuppiTau4/F");                                      
  TreeAllHad->Branch("Jet1PuppiTau32"                        , & Jet1PuppiTau32                     ,    "Jet1PuppiTau32/F");                                       
  TreeAllHad->Branch("Jet1PuppiTau21"                        , & Jet1PuppiTau21                     ,    "Jet1PuppiTau21/F");                                       
  TreeAllHad->Branch("Jet1PuppiSDmaxbdisc"                   , & Jet1PuppiSDmaxbdisc                ,    "Jet1PuppiSDmaxbdisc/F");                                            
  TreeAllHad->Branch("Jet1PuppiSDmaxbdiscflavHadron"         , & Jet1PuppiSDmaxbdiscflavHadron      ,    "Jet1PuppiSDmaxbdiscflavHadron/F");                                                
  TreeAllHad->Branch("Jet1PuppiSDmaxbdiscflavParton"         , & Jet1PuppiSDmaxbdiscflavParton      ,    "Jet1PuppiSDmaxbdiscflavParton/F");                                                
  TreeAllHad->Branch("Jet1PuppiSDsubjet0pt"                  , & Jet1PuppiSDsubjet0pt               ,    "Jet1PuppiSDsubjet0pt/F");                                             
  TreeAllHad->Branch("Jet1PuppiSDsubjet0mass"                , & Jet1PuppiSDsubjet0mass             ,    "Jet1PuppiSDsubjet0mass/F");    
  TreeAllHad->Branch("Jet1PuppiSDsubjet0eta"                 , & Jet1PuppiSDsubjet0eta              ,    "Jet1PuppiSDsubjet0eta/F");
  TreeAllHad->Branch("Jet1PuppiSDsubjet0phi"                 , & Jet1PuppiSDsubjet0phi              ,    "Jet1PuppiSDsubjet0phi/F");                                           
  TreeAllHad->Branch("Jet1PuppiSDsubjet0area"                , & Jet1PuppiSDsubjet0area             ,    "Jet1PuppiSDsubjet0area/F");                                               
  TreeAllHad->Branch("Jet1PuppiSDsubjet0flavHadron"          , & Jet1PuppiSDsubjet0flavHadron       ,    "Jet1PuppiSDsubjet0flavHadron/F");                                               
  TreeAllHad->Branch("Jet1PuppiSDsubjet0flavParton"          , & Jet1PuppiSDsubjet0flavParton       ,    "Jet1PuppiSDsubjet0flavParton/F"); 
  TreeAllHad->Branch("Jet1PuppiSDsubjet0tau1"                , & Jet1PuppiSDsubjet0tau1             ,    "Jet1PuppiSDsubjet0tau1/F");
  TreeAllHad->Branch("Jet1PuppiSDsubjet0tau2"                , & Jet1PuppiSDsubjet0tau2             ,    "Jet1PuppiSDsubjet0tau2/F");
  TreeAllHad->Branch("Jet1PuppiSDsubjet0tau3"                , & Jet1PuppiSDsubjet0tau3             ,    "Jet1PuppiSDsubjet0tau3/F");
  TreeAllHad->Branch("Jet1PuppiSDsubjet0bdisc"               , & Jet1PuppiSDsubjet0bdisc            ,    "Jet1PuppiSDsubjet0bdisc/F");                                                                                         
  TreeAllHad->Branch("Jet1PuppiSDsubjet1pt"                  , & Jet1PuppiSDsubjet1pt               ,    "Jet1PuppiSDsubjet1pt/F");                                             
  TreeAllHad->Branch("Jet1PuppiSDsubjet1mass"                , & Jet1PuppiSDsubjet1mass             ,    "Jet1PuppiSDsubjet1mass/F");  
  TreeAllHad->Branch("Jet1PuppiSDsubjet1eta"                 , & Jet1PuppiSDsubjet1eta              ,    "Jet1PuppiSDsubjet1eta/F");
  TreeAllHad->Branch("Jet1PuppiSDsubjet1phi"                 , & Jet1PuppiSDsubjet1phi              ,    "Jet1PuppiSDsubjet1phi/F");                                             
  TreeAllHad->Branch("Jet1PuppiSDsubjet1area"                , & Jet1PuppiSDsubjet1area             ,    "Jet1PuppiSDsubjet1area/F");                                               
  TreeAllHad->Branch("Jet1PuppiSDsubjet1flavHadron"          , & Jet1PuppiSDsubjet1flavHadron       ,    "Jet1PuppiSDsubjet1flavHadron/F");                                               
  TreeAllHad->Branch("Jet1PuppiSDsubjet1flavParton"          , & Jet1PuppiSDsubjet1flavParton       ,    "Jet1PuppiSDsubjet1flavParton/F");    
  TreeAllHad->Branch("Jet1PuppiSDsubjet1tau1"                , & Jet1PuppiSDsubjet1tau1             ,    "Jet1PuppiSDsubjet1tau1/F");
  TreeAllHad->Branch("Jet1PuppiSDsubjet1tau2"                , & Jet1PuppiSDsubjet1tau2             ,    "Jet1PuppiSDsubjet1tau2/F");
  TreeAllHad->Branch("Jet1PuppiSDsubjet1tau3"                , & Jet1PuppiSDsubjet1tau3             ,    "Jet1PuppiSDsubjet1tau3/F"); 
  TreeAllHad->Branch("Jet1PuppiSDsubjet1bdisc"               , & Jet1PuppiSDsubjet1bdisc            ,    "Jet1PuppiSDsubjet1bdisc/F");                                                                                       
  TreeAllHad->Branch("Jet1CHF"                               , & Jet1CHF                            ,    "Jet1CHF/F");                                
  TreeAllHad->Branch("Jet1NHF"                               , & Jet1NHF                            ,    "Jet1NHF/F");                                
  TreeAllHad->Branch("Jet1CM"                                , & Jet1CM                             ,    "Jet1CM/F");                               
  TreeAllHad->Branch("Jet1NM"                                , & Jet1NM                             ,    "Jet1NM/F");                               
  TreeAllHad->Branch("Jet1NEF"                               , & Jet1NEF                            ,    "Jet1NEF/F");                                
  TreeAllHad->Branch("Jet1CEF"                               , & Jet1CEF                            ,    "Jet1CEF/F");                                
  TreeAllHad->Branch("Jet1MF"                                , & Jet1MF                             ,    "Jet1MF/F");                               
  TreeAllHad->Branch("Jet1Mult"                              , & Jet1Mult                           ,    "Jet1Mult/F"); 
  TreeAllHad->Branch("Jet1PuppiCHF"                          , & Jet1PuppiCHF                       ,    "Jet1PuppiCHF/F");                                
  TreeAllHad->Branch("Jet1PuppiNHF"                          , & Jet1PuppiNHF                       ,    "Jet1PuppiNHF/F");                                
  TreeAllHad->Branch("Jet1PuppiCM"                           , & Jet1PuppiCM                        ,    "Jet1PuppiCM/F");                               
  TreeAllHad->Branch("Jet1PuppiNM"                           , & Jet1PuppiNM                        ,    "Jet1PuppiNM/F");                               
  TreeAllHad->Branch("Jet1PuppiNEF"                          , & Jet1PuppiNEF                       ,    "Jet1PuppiNEF/F");                                
  TreeAllHad->Branch("Jet1PuppiCEF"                          , & Jet1PuppiCEF                       ,    "Jet1PuppiCEF/F");                                
  TreeAllHad->Branch("Jet1PuppiMF"                           , & Jet1PuppiMF                        ,    "Jet1PuppiMF/F");                               
  TreeAllHad->Branch("Jet1PuppiMult"                         , & Jet1PuppiMult                      ,    "Jet1PuppiMult/F");                                     
  TreeAllHad->Branch("Jet1MassCorrFactor"                    , & Jet1MassCorrFactor                 ,    "Jet1MassCorrFactor/F");                                           
  TreeAllHad->Branch("Jet1MassCorrFactorUp"                  , & Jet1MassCorrFactorUp               ,    "Jet1MassCorrFactorUp/F");                                             
  TreeAllHad->Branch("Jet1MassCorrFactorDn"                  , & Jet1MassCorrFactorDn               ,    "Jet1MassCorrFactorDn/F");                                             
  TreeAllHad->Branch("Jet1CorrFactor"                        , & Jet1CorrFactor                     ,    "Jet1CorrFactor/F");                                       
  TreeAllHad->Branch("Jet1CorrFactorUp"                      , & Jet1CorrFactorUp                   ,    "Jet1CorrFactorUp/F");                                         
  TreeAllHad->Branch("Jet1CorrFactorDn"                      , & Jet1CorrFactorDn                   ,    "Jet1CorrFactorDn/F");                                         
  TreeAllHad->Branch("Jet1PtSmearFactor"                     , & Jet1PtSmearFactor                  ,    "Jet1PtSmearFactor/F");                                          
  TreeAllHad->Branch("Jet1PtSmearFactorUp"                   , & Jet1PtSmearFactorUp                ,    "Jet1PtSmearFactorUp/F");                                            
  TreeAllHad->Branch("Jet1PtSmearFactorDn"                   , & Jet1PtSmearFactorDn                ,    "Jet1PtSmearFactorDn/F");                                            
  TreeAllHad->Branch("Jet1PuppiMassCorrFactor"               , & Jet1PuppiMassCorrFactor            ,    "Jet1PuppiMassCorrFactor/F");                                                
  TreeAllHad->Branch("Jet1PuppiMassCorrFactorUp"             , & Jet1PuppiMassCorrFactorUp          ,    "Jet1PuppiMassCorrFactorUp/F");                                                  
  TreeAllHad->Branch("Jet1PuppiMassCorrFactorDn"             , & Jet1PuppiMassCorrFactorDn          ,    "Jet1PuppiMassCorrFactorDn/F");                                                  
  TreeAllHad->Branch("Jet1PuppiCorrFactor"                   , & Jet1PuppiCorrFactor                ,    "Jet1PuppiCorrFactor/F");                                            
  TreeAllHad->Branch("Jet1PuppiCorrFactorUp"                 , & Jet1PuppiCorrFactorUp              ,    "Jet1PuppiCorrFactorUp/F");                                              
  TreeAllHad->Branch("Jet1PuppiCorrFactorDn"                 , & Jet1PuppiCorrFactorDn              ,    "Jet1PuppiCorrFactorDn/F");                                              
  TreeAllHad->Branch("Jet1PuppiPtSmearFactor"                , & Jet1PuppiPtSmearFactor             ,    "Jet1PuppiPtSmearFactor/F");                                               
  TreeAllHad->Branch("Jet1PuppiPtSmearFactorUp"              , & Jet1PuppiPtSmearFactorUp           ,    "Jet1PuppiPtSmearFactorUp/F");                                                 
  TreeAllHad->Branch("Jet1PuppiPtSmearFactorDn"              , & Jet1PuppiPtSmearFactorDn           ,    "Jet1PuppiPtSmearFactorDn/F");                                                 
  TreeAllHad->Branch("Jet1EtaScaleFactor"                    , & Jet1EtaScaleFactor                 ,    "Jet1EtaScaleFactor/F");                                           
  TreeAllHad->Branch("Jet1PhiScaleFactor"                    , & Jet1PhiScaleFactor                 ,    "Jet1PhiScaleFactor/F");                                           
  // TreeAllHad->Branch("Jet1MatchedGenJetDR"                   , & Jet1MatchedGenJetDR                ,    "Jet1MatchedGenJetDR/F");                                            
  TreeAllHad->Branch("Jet1MatchedGenJetPt"                   , & Jet1MatchedGenJetPt                ,    "Jet1MatchedGenJetPt/F");                                            
  TreeAllHad->Branch("Jet1MatchedGenJetMass"                 , & Jet1MatchedGenJetMass              ,    "Jet1MatchedGenJetMass/F"); 
                        
  TreeAllHad->Branch("Jet1GenMatched_TopHadronic"            , & Jet1GenMatched_TopHadronic         ,    "Jet1GenMatched_TopHadronic/I");      
  TreeAllHad->Branch("Jet1GenMatched_TopPt"                  , & Jet1GenMatched_TopPt               ,    "Jet1GenMatched_TopPt/F");      
  TreeAllHad->Branch("Jet1GenMatched_TopEta"                 , & Jet1GenMatched_TopEta              ,    "Jet1GenMatched_TopEta/F");      
  TreeAllHad->Branch("Jet1GenMatched_TopPhi"                 , & Jet1GenMatched_TopPhi              ,    "Jet1GenMatched_TopPhi/F");      
  TreeAllHad->Branch("Jet1GenMatched_TopMass"                , & Jet1GenMatched_TopMass             ,    "Jet1GenMatched_TopMass/F");      
  TreeAllHad->Branch("Jet1GenMatched_bPt"                    , & Jet1GenMatched_bPt                 ,    "Jet1GenMatched_bPt/F");      
  TreeAllHad->Branch("Jet1GenMatched_WPt"                    , & Jet1GenMatched_WPt                 ,    "Jet1GenMatched_WPt/F");      
  TreeAllHad->Branch("Jet1GenMatched_Wd1Pt"                  , & Jet1GenMatched_Wd1Pt               ,    "Jet1GenMatched_Wd1Pt/F");      
  TreeAllHad->Branch("Jet1GenMatched_Wd2Pt"                  , & Jet1GenMatched_Wd2Pt               ,    "Jet1GenMatched_Wd2Pt/F");      
  TreeAllHad->Branch("Jet1GenMatched_Wd1ID"                  , & Jet1GenMatched_Wd1ID               ,    "Jet1GenMatched_Wd1ID/F");      
  TreeAllHad->Branch("Jet1GenMatched_Wd2ID"                  , & Jet1GenMatched_Wd2ID               ,    "Jet1GenMatched_Wd2ID/F");      
  TreeAllHad->Branch("Jet1GenMatched_MaxDeltaRPartonTop"     , & Jet1GenMatched_MaxDeltaRPartonTop  ,    "Jet1GenMatched_MaxDeltaRPartonTop/F");      
  TreeAllHad->Branch("Jet1GenMatched_MaxDeltaRWPartonTop"    , & Jet1GenMatched_MaxDeltaRWPartonTop ,    "Jet1GenMatched_MaxDeltaRWPartonTop/F");      
  TreeAllHad->Branch("Jet1GenMatched_MaxDeltaRWPartonW"      , & Jet1GenMatched_MaxDeltaRWPartonW   ,    "Jet1GenMatched_MaxDeltaRWPartonW/F");      
  TreeAllHad->Branch("Jet1GenMatched_DeltaR_t_b"             , & Jet1GenMatched_DeltaR_t_b          ,    "Jet1GenMatched_DeltaR_t_b/F");      
  TreeAllHad->Branch("Jet1GenMatched_DeltaR_t_W"             , & Jet1GenMatched_DeltaR_t_W          ,    "Jet1GenMatched_DeltaR_t_W/F");      
  TreeAllHad->Branch("Jet1GenMatched_DeltaR_t_Wd1"           , & Jet1GenMatched_DeltaR_t_Wd1        ,    "Jet1GenMatched_DeltaR_t_Wd1/F");      
  TreeAllHad->Branch("Jet1GenMatched_DeltaR_t_Wd2"           , & Jet1GenMatched_DeltaR_t_Wd2        ,    "Jet1GenMatched_DeltaR_t_Wd2/F");      
  TreeAllHad->Branch("Jet1GenMatched_DeltaR_W_b1"            , & Jet1GenMatched_DeltaR_W_b1         ,    "Jet1GenMatched_DeltaR_W_b1/F");      
  TreeAllHad->Branch("Jet1GenMatched_DeltaR_W_Wd1"           , & Jet1GenMatched_DeltaR_W_Wd1        ,    "Jet1GenMatched_DeltaR_W_Wd1/F");      
  TreeAllHad->Branch("Jet1GenMatched_DeltaR_W_Wd2"           , & Jet1GenMatched_DeltaR_W_Wd2        ,    "Jet1GenMatched_DeltaR_W_Wd2/F");      
  TreeAllHad->Branch("Jet1GenMatched_DeltaR_Wd1_Wd2"         , & Jet1GenMatched_DeltaR_Wd1_Wd2      ,    "Jet1GenMatched_DeltaR_Wd1_Wd2/F");      
  TreeAllHad->Branch("Jet1GenMatched_DeltaR_Wd1_b"           , & Jet1GenMatched_DeltaR_Wd1_b        ,    "Jet1GenMatched_DeltaR_Wd1_b/F");      
  TreeAllHad->Branch("Jet1GenMatched_DeltaR_Wd2_b"           , & Jet1GenMatched_DeltaR_Wd2_b        ,    "Jet1GenMatched_DeltaR_Wd2_b/F");      
  TreeAllHad->Branch("Jet1GenMatched_DeltaR_jet_t"           , & Jet1GenMatched_DeltaR_jet_t        ,    "Jet1GenMatched_DeltaR_jet_t/F");      
  TreeAllHad->Branch("Jet1GenMatched_DeltaR_jet_W"           , & Jet1GenMatched_DeltaR_jet_W        ,    "Jet1GenMatched_DeltaR_jet_W/F");      
  TreeAllHad->Branch("Jet1GenMatched_DeltaR_jet_b"           , & Jet1GenMatched_DeltaR_jet_b        ,    "Jet1GenMatched_DeltaR_jet_b/F");      
  TreeAllHad->Branch("Jet1GenMatched_DeltaR_jet_Wd1"         , & Jet1GenMatched_DeltaR_jet_Wd1      ,    "Jet1GenMatched_DeltaR_jet_Wd1/F");      
  TreeAllHad->Branch("Jet1GenMatched_DeltaR_jet_Wd2"         , & Jet1GenMatched_DeltaR_jet_Wd2      ,    "Jet1GenMatched_DeltaR_jet_Wd2/F");      
  TreeAllHad->Branch("Jet1GenMatched_DeltaR_pup0_b"          , & Jet1GenMatched_DeltaR_pup0_b       ,    "Jet1GenMatched_DeltaR_pup0_b/F");      
  TreeAllHad->Branch("Jet1GenMatched_DeltaR_pup0_Wd1"        , & Jet1GenMatched_DeltaR_pup0_Wd1     ,    "Jet1GenMatched_DeltaR_pup0_Wd1/F");      
  TreeAllHad->Branch("Jet1GenMatched_DeltaR_pup0_Wd2"        , & Jet1GenMatched_DeltaR_pup0_Wd2     ,    "Jet1GenMatched_DeltaR_pup0_Wd2/F");      
  TreeAllHad->Branch("Jet1GenMatched_DeltaR_pup1_b"          , & Jet1GenMatched_DeltaR_pup1_b       ,    "Jet1GenMatched_DeltaR_pup1_b/F");      
  TreeAllHad->Branch("Jet1GenMatched_DeltaR_pup1_Wd1"        , & Jet1GenMatched_DeltaR_pup1_Wd1     ,    "Jet1GenMatched_DeltaR_pup1_Wd1/F");      
  TreeAllHad->Branch("Jet1GenMatched_DeltaR_pup1_Wd2"        , & Jet1GenMatched_DeltaR_pup1_Wd2     ,    "Jet1GenMatched_DeltaR_pup1_Wd2/F");               
  TreeAllHad->Branch("Jet1GenMatched_partonPt"               , & Jet1GenMatched_partonPt            ,    "Jet1GenMatched_partonPt/F");      
  TreeAllHad->Branch("Jet1GenMatched_partonEta"              , & Jet1GenMatched_partonEta           ,    "Jet1GenMatched_partonEta/F");      
  TreeAllHad->Branch("Jet1GenMatched_partonPhi"              , & Jet1GenMatched_partonPhi           ,    "Jet1GenMatched_partonPhi/F");      
  TreeAllHad->Branch("Jet1GenMatched_partonMass"             , & Jet1GenMatched_partonMass          ,    "Jet1GenMatched_partonMass/F");      
  TreeAllHad->Branch("Jet1GenMatched_partonID"               , & Jet1GenMatched_partonID            ,    "Jet1GenMatched_partonID/F");      
  TreeAllHad->Branch("Jet1GenMatched_DeltaRjetParton"        , & Jet1GenMatched_DeltaRjetParton     ,    "Jet1GenMatched_DeltaRjetParton/F");      
  std::cout << "Setup Jet1 in tree" << std::endl;

  TreeAllHad->Branch("AllHadMETpx"                           , & AllHadMETpx                        ,    "AllHadMETpx/F");                                    
  TreeAllHad->Branch("AllHadMETpy"                           , & AllHadMETpy                        ,    "AllHadMETpy/F");                                    
  TreeAllHad->Branch("AllHadMETpt"                           , & AllHadMETpt                        ,    "AllHadMETpt/F");                                    
  TreeAllHad->Branch("AllHadMETphi"                          , & AllHadMETphi                       ,    "AllHadMETphi/F");                                     
  TreeAllHad->Branch("AllHadMETsumET"                        , & AllHadMETsumET                     ,    "AllHadMETsumET/F");                                     
  TreeAllHad->Branch("AllHadNvtx"                            , & AllHadNvtx                         ,    "AllHadNvtx/F");                                   
  TreeAllHad->Branch("AllHadRho"                             , & AllHadRho                          ,    "AllHadRho/F");                                  
  TreeAllHad->Branch("AllHadEventWeight"                     , & AllHadEventWeight                  ,    "AllHadEventWeight/F");                                          
  TreeAllHad->Branch("AllHadPUweight"                  , & AllHadPUweight               , "AllHadPUweight/F");
  TreeAllHad->Branch("AllHadPUweight_MBup"                  , & AllHadPUweight_MBup               , "AllHadPUweight_MBup/F");
  TreeAllHad->Branch("AllHadPUweight_MBdn"                  , & AllHadPUweight_MBdn               , "AllHadPUweight_MBdn/F");
       
  
  

  TreeAllHad->Branch("DijetMass"                             , & DijetMass                          ,    "DijetMass/F");                                  
  TreeAllHad->Branch("DijetMassPuppi"                        , & DijetMassPuppi                     ,    "DijetMassPuppi/F");                                  
  TreeAllHad->Branch("DijetDeltaR"                           , & DijetDeltaR                        ,    "DijetDeltaR/F");                                    
  TreeAllHad->Branch("DijetDeltaPhi"                         , & DijetDeltaPhi                      ,    "DijetDeltaPhi/F");                                      
  TreeAllHad->Branch("DijetDeltaRap"                         , & DijetDeltaRap                      ,    "DijetDeltaRap/F");                                      
  TreeAllHad->Branch("DiGenJetMass"                          , & DiGenJetMass                       ,    "DiGenJetMass/F");                                     
  TreeAllHad->Branch("GenTTmass"                             , & GenTTmass                          ,    "GenTTmass/F");                                  
  TreeAllHad->Branch("HT"                                    , & HT                                 ,    "HT/F");                           
  TreeAllHad->Branch("HT_CorrDn"                             , & HT_CorrDn                          ,    "HT_CorrDn/F");                                  
  TreeAllHad->Branch("HT_CorrUp"                             , & HT_CorrUp                          ,    "HT_CorrUp/F");                                  
  TreeAllHad->Branch("HT_PtSmearNom"                         , & HT_PtSmearNom                      ,    "HT_PtSmearNom/F");                                     
  TreeAllHad->Branch("HT_PtSmearUp"                          , & HT_PtSmearUp                       ,    "HT_PtSmearUp/F");                                     
  TreeAllHad->Branch("HT_PtSmearDn"                          , & HT_PtSmearDn                       ,    "HT_PtSmearDn/F");                                     
  TreeAllHad->Branch("Q2weight_CorrDn"                       , & Q2weight_CorrDn                    ,    "Q2weight_CorrDn/F");                                        
  TreeAllHad->Branch("Q2weight_CorrUp"                       , & Q2weight_CorrUp                    ,    "Q2weight_CorrUp/F");                                        
  TreeAllHad->Branch("NNPDF3weight_CorrDn"                   , & NNPDF3weight_CorrDn                ,    "NNPDF3weight_CorrDn/F");                                            
  TreeAllHad->Branch("NNPDF3weight_CorrUp"                   , & NNPDF3weight_CorrUp                ,    "NNPDF3weight_CorrUp/F");                                            
  TreeAllHad->Branch("AllHadRunNum"                          , & AllHadRunNum                       ,    "AllHadRunNum/F");                                     
  TreeAllHad->Branch("AllHadLumiBlock"                       , & AllHadLumiBlock                    ,    "AllHadLumiBlock/F");                                        
  TreeAllHad->Branch("AllHadEventNum"                        , & AllHadEventNum                     ,    "AllHadEventNum/F"); */                                      
  
  std::cout << "Setup all-had event tree" << std::endl;

  // 
  //  .d8888b.                         d8b        888                        888        88888888888                           
  // d88P  Y88b                        Y8P        888                        888            888                               
  // Y88b.                                        888                        888            888                               
  //  "Y888b.    .d88b.  88888b.d88b.  888        888       .d88b.  88888b.  888888         888     888d888  .d88b.   .d88b.  
  //     "Y88b. d8P  Y8b 888 "888 "88b 888        888      d8P  Y8b 888 "88b 888            888     888P"   d8P  Y8b d8P  Y8b 
  //       "888 88888888 888  888  888 888 888888 888      88888888 888  888 888            888     888     88888888 88888888 
  // Y88b  d88P Y8b.     888  888  888 888        888      Y8b.     888 d88P Y88b.          888     888     Y8b.     Y8b.     
  //  "Y8888P"   "Y8888  888  888  888 888        88888888  "Y8888  88888P"   "Y888         888     888      "Y8888   "Y8888  
  //                                                                888                                                       
  //                                                                888                                                       
  //                                                                888       
  
  TreeSemiLept = new TTree("TreeSemiLept", "TreeSemiLept"); 
       
  /*TreeSemiLept->Branch("SemiLeptTrigNames"    , "vector<std::string>", &SemiLeptTrigNames);
  TreeSemiLept->Branch("SemiLeptTrigPrescales" , "vector<int>",  &SemiLeptTrigPrescales);
  TreeSemiLept->Branch("SemiLeptTrigPass"      , "vector<bool>", &SemiLeptTrigPass);
  TreeSemiLept->Branch("SemiLeptTrigAcceptBits", &SemiLeptTrigAcceptBits);



  TreeSemiLept->Branch("JetPtRaw"                             , & JetPtRaw                          ,    "JetPtRaw/F");                                  
  TreeSemiLept->Branch("JetEtaRaw"                            , & JetEtaRaw                         ,    "JetEtaRaw/F");                                   
  TreeSemiLept->Branch("JetPhiRaw"                            , & JetPhiRaw                         ,    "JetPhiRaw/F");                                   
  TreeSemiLept->Branch("JetMassRaw"                           , & JetMassRaw                        ,    "JetMassRaw/F");                                    
  TreeSemiLept->Branch("JetP"                                 , & JetP                              ,    "JetP/F");                              
  TreeSemiLept->Branch("JetPt"                                , & JetPt                             ,    "JetPt/F");                               
  TreeSemiLept->Branch("JetEta"                               , & JetEta                            ,    "JetEta/F");                                
  TreeSemiLept->Branch("JetPhi"                               , & JetPhi                            ,    "JetPhi/F");                                
  TreeSemiLept->Branch("JetRap"                               , & JetRap                            ,    "JetRap/F");                                
  TreeSemiLept->Branch("JetEnergy"                            , & JetEnergy                         ,    "JetEnergy/F");                                   
  TreeSemiLept->Branch("JetMass"                              , & JetMass                           ,    "JetMass/F");                                 
  TreeSemiLept->Branch("JetArea"                              , & JetArea                           ,    "JetArea/F");                                 
  
  TreeSemiLept->Branch("JetSDmass"                            , & JetSDmass                         ,    "JetSDmass/F");                                         
  TreeSemiLept->Branch("JetSDmassRaw"                         , & JetSDmassRaw                      ,    "JetSDmassRaw/F");                                               
  TreeSemiLept->Branch("JetSDmassCorrL23"                     , & JetSDmassCorrL23                  ,    "JetSDmassCorrL23/F");                                                    
  TreeSemiLept->Branch("JetSDmassCorrL23Up"                   , & JetSDmassCorrL23Up                ,    "JetSDmassCorrL23Up/F");                                                      
  TreeSemiLept->Branch("JetSDmassCorrL23Dn"                   , & JetSDmassCorrL23Dn                ,    "JetSDmassCorrL23Dn/F");                                                      
  TreeSemiLept->Branch("JetSDmassCorrL123"                    , & JetSDmassCorrL123                 ,    "JetSDmassCorrL123/F");                                                      
  TreeSemiLept->Branch("JetSDmassCorrL123Up"                  , & JetSDmassCorrL123Up               ,    "JetSDmassCorrL123Up/F");                                                        
  TreeSemiLept->Branch("JetSDmassCorrL123Dn"                  , & JetSDmassCorrL123Dn               ,    "JetSDmassCorrL123Dn/F");                                                        
  TreeSemiLept->Branch("JetSDmassCorrL23Smear"                , & JetSDmassCorrL23Smear             ,    "JetSDmassCorrL23Smear/F");                                                     
  TreeSemiLept->Branch("JetSDmassCorrL23SmearUp"              , & JetSDmassCorrL23SmearUp           ,    "JetSDmassCorrL23SmearUp/F");                                                       
  TreeSemiLept->Branch("JetSDmassCorrL23SmearDn"              , & JetSDmassCorrL23SmearDn           ,    "JetSDmassCorrL23SmearDn/F");   
  TreeSemiLept->Branch("JetSDptRaw"                           , & JetSDptRaw                        ,    "JetSDptRaw/F");                                               
  TreeSemiLept->Branch("JetSDptCorrL23"                       , & JetSDptCorrL23                    ,    "JetSDptCorrL23/F");                                                    
  TreeSemiLept->Branch("JetSDptCorrL23Up"                     , & JetSDptCorrL23Up                  ,    "JetSDptCorrL23Up/F");                                                      
  TreeSemiLept->Branch("JetSDptCorrL23Dn"                     , & JetSDptCorrL23Dn                  ,    "JetSDptCorrL23Dn/F");                                                      
  TreeSemiLept->Branch("JetSDptCorrL123"                      , & JetSDptCorrL123                   ,    "JetSDptCorrL123/F");                                                      
  TreeSemiLept->Branch("JetSDptCorrL123Up"                    , & JetSDptCorrL123Up                 ,    "JetSDptCorrL123Up/F");                                                        
  TreeSemiLept->Branch("JetSDptCorrL123Dn"                    , & JetSDptCorrL123Dn                 ,    "JetSDptCorrL123Dn/F");                                                        
  TreeSemiLept->Branch("JetSDptCorrL23Smear"                  , & JetSDptCorrL23Smear               ,    "JetSDptCorrL23Smear/F");                                                     
  TreeSemiLept->Branch("JetSDptCorrL23SmearUp"                , & JetSDptCorrL23SmearUp             ,    "JetSDptCorrL23SmearUp/F");                                                       
  TreeSemiLept->Branch("JetSDptCorrL23SmearDn"                , & JetSDptCorrL23SmearDn             ,    "JetSDptCorrL23SmearDn/F");                                                     
  TreeSemiLept->Branch("JetSDetaRaw"                          , & JetSDetaRaw                       ,    "JetSDetaRaw/F");                                               
  TreeSemiLept->Branch("JetSDphiRaw"                          , & JetSDphiRaw                       ,    "JetSDphiRaw/F");  

  TreeSemiLept->Branch("JetMassPruned"                        , & JetMassPruned                     ,    "JetMassPruned/F");                                       
  TreeSemiLept->Branch("JetMassTrimmed"                       , & JetMassTrimmed                    ,    "JetMassTrimmed/F");                                       
  TreeSemiLept->Branch("JetTau1"                              , & JetTau1                           ,    "JetTau1/F");                                 
  TreeSemiLept->Branch("JetTau2"                              , & JetTau2                           ,    "JetTau2/F");                                 
  TreeSemiLept->Branch("JetTau3"                              , & JetTau3                           ,    "JetTau3/F");                                 
  TreeSemiLept->Branch("JetTau4"                              , & JetTau4                           ,    "JetTau4/F");                                 
  TreeSemiLept->Branch("JetTau32"                             , & JetTau32                          ,    "JetTau32/F");                                  
  TreeSemiLept->Branch("JetTau21"                             , & JetTau21                          ,    "JetTau21/F");                                  
  TreeSemiLept->Branch("JetSDmaxbdisc"                        , & JetSDmaxbdisc                     ,    "JetSDmaxbdisc/F");                                       
  TreeSemiLept->Branch("JetSDmaxbdiscflavHadron"              , & JetSDmaxbdiscflavHadron           ,    "JetSDmaxbdiscflavHadron/F");                                           
  TreeSemiLept->Branch("JetSDmaxbdiscflavParton"              , & JetSDmaxbdiscflavParton           ,    "JetSDmaxbdiscflavParton/F");  
                                         
  TreeSemiLept->Branch("JetSDsubjet0pt"                       , & JetSDsubjet0pt                    ,    "JetSDsubjet0pt/F");    
  TreeSemiLept->Branch("JetSDsubjet0mass"                     , & JetSDsubjet0mass                  ,    "JetSDsubjet0mass/F");
  TreeSemiLept->Branch("JetSDsubjet0eta"                      , & JetSDsubjet0eta                   ,    "JetSDsubjet0eta/F");
  TreeSemiLept->Branch("JetSDsubjet0phi"                      , & JetSDsubjet0phi                   ,    "JetSDsubjet0phi/F");
  TreeSemiLept->Branch("JetSDsubjet0area"                     , & JetSDsubjet0area                  ,    "JetSDsubjet0area/F");
  TreeSemiLept->Branch("JetSDsubjet0flavHadron"               , & JetSDsubjet0flavHadron            ,    "JetSDsubjet0flavHadron/F");
  TreeSemiLept->Branch("JetSDsubjet0flavParton"               , & JetSDsubjet0flavParton            ,    "JetSDsubjet0flavParton/F");
  TreeSemiLept->Branch("JetSDsubjet0tau1"                     , & JetSDsubjet0tau1                  ,    "JetSDsubjet0tau1/F");
  TreeSemiLept->Branch("JetSDsubjet0tau2"                     , & JetSDsubjet0tau2                  ,    "JetSDsubjet0tau2/F");
  TreeSemiLept->Branch("JetSDsubjet0tau3"                     , & JetSDsubjet0tau3                  ,    "JetSDsubjet0tau3/F"); 
  TreeSemiLept->Branch("JetSDsubjet0bdisc"                    , & JetSDsubjet0bdisc                 ,    "JetSDsubjet0bdisc/F");                                          
  TreeSemiLept->Branch("JetSDsubjet1pt"                       , & JetSDsubjet1pt                    ,    "JetSDsubjet1pt/F");    
  TreeSemiLept->Branch("JetSDsubjet1mass"                     , & JetSDsubjet1mass                  ,    "JetSDsubjet1mass/F");
  TreeSemiLept->Branch("JetSDsubjet1eta"                      , & JetSDsubjet1eta                   ,    "JetSDsubjet1eta/F");
  TreeSemiLept->Branch("JetSDsubjet1phi"                      , & JetSDsubjet1phi                   ,    "JetSDsubjet1phi/F");  
  TreeSemiLept->Branch("JetSDsubjet1area"                     , & JetSDsubjet1area                  ,    "JetSDsubjet1area/F");
  TreeSemiLept->Branch("JetSDsubjet1flavHadron"               , & JetSDsubjet1flavHadron            ,    "JetSDsubjet1flavHadron/F");
  TreeSemiLept->Branch("JetSDsubjet1flavParton"               , & JetSDsubjet1flavParton            ,    "JetSDsubjet1flavParton/F");
  TreeSemiLept->Branch("JetSDsubjet1tau1"                     , & JetSDsubjet1tau1                  ,    "JetSDsubjet1tau1/F");
  TreeSemiLept->Branch("JetSDsubjet1tau2"                     , & JetSDsubjet1tau2                  ,    "JetSDsubjet1tau2/F");
  TreeSemiLept->Branch("JetSDsubjet1tau3"                     , & JetSDsubjet1tau3                  ,    "JetSDsubjet1tau3/F");                                           
  TreeSemiLept->Branch("JetSDsubjet1bdisc"                    , & JetSDsubjet1bdisc                 ,    "JetSDsubjet1bdisc/F");                                     

  TreeSemiLept->Branch("JetPuppiP"                            , & JetPuppiP                         ,    "JetPuppiP/F");                                    
  TreeSemiLept->Branch("JetPuppiPt"                           , & JetPuppiPt                        ,    "JetPuppiPt/F");                                    
  TreeSemiLept->Branch("JetPuppiEta"                          , & JetPuppiEta                       ,    "JetPuppiEta/F");                                     
  TreeSemiLept->Branch("JetPuppiPhi"                          , & JetPuppiPhi                       ,    "JetPuppiPhi/F");                                     
  TreeSemiLept->Branch("JetPuppiMass"                         , & JetPuppiMass                      ,    "JetPuppiMass/F");                                      

  
  TreeSemiLept->Branch("JetPuppiSDmass"                         , & JetPuppiSDmass                    ,    "JetPuppiSDmass/F");
  TreeSemiLept->Branch("JetPuppiSDmassCorr"                     , & JetPuppiSDmassCorr                ,    "JetPuppiSDmassCorr/F");
  TreeSemiLept->Branch("JetPuppiSDmassCorrUp"                   , & JetPuppiSDmassCorrUp              ,    "JetPuppiSDmassCorrUp/F");
  TreeSemiLept->Branch("JetPuppiSDmassCorrDn"                   , & JetPuppiSDmassCorrDn              ,    "JetPuppiSDmassCorrDn/F");
  TreeSemiLept->Branch("JetPuppiSDmassCorrL23Smear"             , & JetPuppiSDmassCorrL23Smear        ,    "JetPuppiSDmassCorrL23Smear/F");
  TreeSemiLept->Branch("JetPuppiSDmassCorrL23SmearUp"           , & JetPuppiSDmassCorrL23SmearUp      ,    "JetPuppiSDmassCorrL23SmearUp/F");
  TreeSemiLept->Branch("JetPuppiSDmassCorrL23SmearDn"           , & JetPuppiSDmassCorrL23SmearDn      ,    "JetPuppiSDmassCorrL23SmearDn/F");
  TreeSemiLept->Branch("JetPuppiSDpt"                           , & JetPuppiSDpt                      ,    "JetPuppiSDpt/F");
  TreeSemiLept->Branch("JetPuppiSDptCorr"                       , & JetPuppiSDptCorr                  ,    "JetPuppiSDptCorr/F");
  TreeSemiLept->Branch("JetPuppiSDptCorrUp"                     , & JetPuppiSDptCorrUp                ,    "JetPuppiSDptCorrUp/F");
  TreeSemiLept->Branch("JetPuppiSDptCorrDn"                     , & JetPuppiSDptCorrDn                ,    "JetPuppiSDptCorrDn/F");
  TreeSemiLept->Branch("JetPuppiSDptCorrL23Smear"               , & JetPuppiSDptCorrL23Smear          ,    "JetPuppiSDptCorrL23Smear/F");
  TreeSemiLept->Branch("JetPuppiSDptCorrL23SmearUp"             , & JetPuppiSDptCorrL23SmearUp        ,    "JetPuppiSDptCorrL23SmearUp/F");
  TreeSemiLept->Branch("JetPuppiSDptCorrL23SmearDn"             , & JetPuppiSDptCorrL23SmearDn        ,    "JetPuppiSDptCorrL23SmearDn/F");
  TreeSemiLept->Branch("JetPuppiSDeta"                          , & JetPuppiSDeta                     ,    "JetPuppiSDeta/F");
  TreeSemiLept->Branch("JetPuppiSDphi"                          , & JetPuppiSDphi                     ,    "JetPuppiSDphi/F");
                         

  TreeSemiLept->Branch("JetPuppiTau1"                         , & JetPuppiTau1                      ,    "JetPuppiTau1/F");                                      
  TreeSemiLept->Branch("JetPuppiTau2"                         , & JetPuppiTau2                      ,    "JetPuppiTau2/F");                                      
  TreeSemiLept->Branch("JetPuppiTau3"                         , & JetPuppiTau3                      ,    "JetPuppiTau3/F");                                      
  TreeSemiLept->Branch("JetPuppiTau4"                         , & JetPuppiTau4                      ,    "JetPuppiTau4/F");                                      
  TreeSemiLept->Branch("JetPuppiTau32"                        , & JetPuppiTau32                     ,    "JetPuppiTau32/F");                                       
  TreeSemiLept->Branch("JetPuppiTau21"                        , & JetPuppiTau21                     ,    "JetPuppiTau21/F");                                       

  TreeSemiLept->Branch("JetPuppiSDmaxbdisc"                   , & JetPuppiSDmaxbdisc                ,    "JetPuppiSDmaxbdisc/F");                                            
  TreeSemiLept->Branch("JetPuppiSDmaxbdiscflavHadron"         , & JetPuppiSDmaxbdiscflavHadron      ,    "JetPuppiSDmaxbdiscflavHadron/F");                                                
  TreeSemiLept->Branch("JetPuppiSDmaxbdiscflavParton"         , & JetPuppiSDmaxbdiscflavParton      ,    "JetPuppiSDmaxbdiscflavParton/F");                                                
  TreeSemiLept->Branch("JetPuppiSDsubjet0pt"                  , & JetPuppiSDsubjet0pt               ,    "JetPuppiSDsubjet0pt/F");    
  TreeSemiLept->Branch("JetPuppiSDsubjet0mass"                , & JetPuppiSDsubjet0mass             ,    "JetPuppiSDsubjet0mass/F");
  TreeSemiLept->Branch("JetPuppiSDsubjet0eta"                 , & JetPuppiSDsubjet0eta              ,    "JetPuppiSDsubjet0eta/F");
  TreeSemiLept->Branch("JetPuppiSDsubjet0phi"                 , & JetPuppiSDsubjet0phi              ,    "JetPuppiSDsubjet0phi/F");
  TreeSemiLept->Branch("JetPuppiSDsubjet0area"                , & JetPuppiSDsubjet0area             ,    "JetPuppiSDsubjet0area/F");
  TreeSemiLept->Branch("JetPuppiSDsubjet0flavHadron"          , & JetPuppiSDsubjet0flavHadron       ,    "JetPuppiSDsubjet0flavHadron/F");
  TreeSemiLept->Branch("JetPuppiSDsubjet0flavParton"          , & JetPuppiSDsubjet0flavParton       ,    "JetPuppiSDsubjet0flavParton/F");
  TreeSemiLept->Branch("JetPuppiSDsubjet0tau1"                , & JetPuppiSDsubjet0tau1             ,    "JetPuppiSDsubjet0tau1/F");
  TreeSemiLept->Branch("JetPuppiSDsubjet0tau2"                , & JetPuppiSDsubjet0tau2             ,    "JetPuppiSDsubjet0tau2/F");
  TreeSemiLept->Branch("JetPuppiSDsubjet0tau3"                , & JetPuppiSDsubjet0tau3             ,    "JetPuppiSDsubjet0tau3/F");
  TreeSemiLept->Branch("JetPuppiSDsubjet0bdisc"               , & JetPuppiSDsubjet0bdisc            ,    "JetPuppiSDsubjet0bdisc/F");                                                
  TreeSemiLept->Branch("JetPuppiSDsubjet1pt"                  , & JetPuppiSDsubjet1pt               ,    "JetPuppiSDsubjet1pt/F");    
  TreeSemiLept->Branch("JetPuppiSDsubjet1mass"                , & JetPuppiSDsubjet1mass             ,    "JetPuppiSDsubjet1mass/F");
  TreeSemiLept->Branch("JetPuppiSDsubjet1eta"                 , & JetPuppiSDsubjet1eta              ,    "JetPuppiSDsubjet1eta/F");
  TreeSemiLept->Branch("JetPuppiSDsubjet1phi"                 , & JetPuppiSDsubjet1phi              ,    "JetPuppiSDsubjet1phi/F");  
  TreeSemiLept->Branch("JetPuppiSDsubjet1area"                , & JetPuppiSDsubjet1area             ,    "JetPuppiSDsubjet1area/F");
  TreeSemiLept->Branch("JetPuppiSDsubjet1flavHadron"          , & JetPuppiSDsubjet1flavHadron       ,    "JetPuppiSDsubjet1flavHadron/F");
  TreeSemiLept->Branch("JetPuppiSDsubjet1flavParton"          , & JetPuppiSDsubjet1flavParton       ,    "JetPuppiSDsubjet1flavParton/F");
  TreeSemiLept->Branch("JetPuppiSDsubjet1tau1"                , & JetPuppiSDsubjet1tau1             ,    "JetPuppiSDsubjet1tau1/F");
  TreeSemiLept->Branch("JetPuppiSDsubjet1tau2"                , & JetPuppiSDsubjet1tau2             ,    "JetPuppiSDsubjet1tau2/F");
  TreeSemiLept->Branch("JetPuppiSDsubjet1tau3"                , & JetPuppiSDsubjet1tau3             ,    "JetPuppiSDsubjet1tau3/F");  
  TreeSemiLept->Branch("JetPuppiSDsubjet1bdisc"               , & JetPuppiSDsubjet1bdisc            ,    "JetPuppiSDsubjet1bdisc/F");                                                                                                                        


  TreeSemiLept->Branch("JetCHF"                               , & JetCHF                            ,    "JetCHF/F");                                
  TreeSemiLept->Branch("JetNHF"                               , & JetNHF                            ,    "JetNHF/F");                                
  TreeSemiLept->Branch("JetCM"                                , & JetCM                             ,    "JetCM/F");                               
  TreeSemiLept->Branch("JetNM"                                , & JetNM                             ,    "JetNM/F");                               
  TreeSemiLept->Branch("JetNEF"                               , & JetNEF                            ,    "JetNEF/F");                                
  TreeSemiLept->Branch("JetCEF"                               , & JetCEF                            ,    "JetCEF/F");                                
  TreeSemiLept->Branch("JetMF"                                , & JetMF                             ,    "JetMF/F");                               
  TreeSemiLept->Branch("JetMult"                              , & JetMult                           ,    "JetMult/F");
  TreeSemiLept->Branch("JetPuppiCHF"                          , & JetPuppiCHF                       ,    "JetPuppiCHF/F");                                
  TreeSemiLept->Branch("JetPuppiNHF"                          , & JetPuppiNHF                       ,    "JetPuppiNHF/F");                                
  TreeSemiLept->Branch("JetPuppiCM"                           , & JetPuppiCM                        ,    "JetPuppiCM/F");                               
  TreeSemiLept->Branch("JetPuppiNM"                           , & JetPuppiNM                        ,    "JetPuppiNM/F");                               
  TreeSemiLept->Branch("JetPuppiNEF"                          , & JetPuppiNEF                       ,    "JetPuppiNEF/F");                                
  TreeSemiLept->Branch("JetPuppiCEF"                          , & JetPuppiCEF                       ,    "JetPuppiCEF/F");                                
  TreeSemiLept->Branch("JetPuppiMF"                           , & JetPuppiMF                        ,    "JetPuppiMF/F");                               
  TreeSemiLept->Branch("JetPuppiMult"                         , & JetPuppiMult                      ,    "JetPuppiMult/F");                                  
  TreeSemiLept->Branch("JetMassCorrFactor"                    , & JetMassCorrFactor                 ,    "JetMassCorrFactor/F");                                           
  TreeSemiLept->Branch("JetMassCorrFactorUp"                  , & JetMassCorrFactorUp               ,    "JetMassCorrFactorUp/F");                                             
  TreeSemiLept->Branch("JetMassCorrFactorDn"                  , & JetMassCorrFactorDn               ,    "JetMassCorrFactorDn/F");                                             
  TreeSemiLept->Branch("JetCorrFactor"                        , & JetCorrFactor                     ,    "JetCorrFactor/F");                                       
  TreeSemiLept->Branch("JetCorrFactorUp"                      , & JetCorrFactorUp                   ,    "JetCorrFactorUp/F");                                         
  TreeSemiLept->Branch("JetCorrFactorDn"                      , & JetCorrFactorDn                   ,    "JetCorrFactorDn/F");                                         
  TreeSemiLept->Branch("JetPtSmearFactor"                     , & JetPtSmearFactor                  ,    "JetPtSmearFactor/F");                                          
  TreeSemiLept->Branch("JetPtSmearFactorUp"                   , & JetPtSmearFactorUp                ,    "JetPtSmearFactorUp/F");                                            
  TreeSemiLept->Branch("JetPtSmearFactorDn"                   , & JetPtSmearFactorDn                ,    "JetPtSmearFactorDn/F");                                            
  TreeSemiLept->Branch("JetPuppiMassCorrFactor"               , & JetPuppiMassCorrFactor            ,    "JetPuppiMassCorrFactor/F");                                                
  TreeSemiLept->Branch("JetPuppiMassCorrFactorUp"             , & JetPuppiMassCorrFactorUp          ,    "JetPuppiMassCorrFactorUp/F");                                                  
  TreeSemiLept->Branch("JetPuppiMassCorrFactorDn"             , & JetPuppiMassCorrFactorDn          ,    "JetPuppiMassCorrFactorDn/F");                                                  
  TreeSemiLept->Branch("JetPuppiCorrFactor"                   , & JetPuppiCorrFactor                ,    "JetPuppiCorrFactor/F");                                            
  TreeSemiLept->Branch("JetPuppiCorrFactorUp"                 , & JetPuppiCorrFactorUp              ,    "JetPuppiCorrFactorUp/F");                                              
  TreeSemiLept->Branch("JetPuppiCorrFactorDn"                 , & JetPuppiCorrFactorDn              ,    "JetPuppiCorrFactorDn/F");                                              
  TreeSemiLept->Branch("JetPuppiPtSmearFactor"                , & JetPuppiPtSmearFactor             ,    "JetPuppiPtSmearFactor/F");                                               
  TreeSemiLept->Branch("JetPuppiPtSmearFactorUp"              , & JetPuppiPtSmearFactorUp           ,    "JetPuppiPtSmearFactorUp/F");                                                 
  TreeSemiLept->Branch("JetPuppiPtSmearFactorDn"              , & JetPuppiPtSmearFactorDn           ,    "JetPuppiPtSmearFactorDn/F");                                                 
  TreeSemiLept->Branch("JetEtaScaleFactor"                    , & JetEtaScaleFactor                 ,    "JetEtaScaleFactor/F");                                           
  TreeSemiLept->Branch("JetPhiScaleFactor"                    , & JetPhiScaleFactor                 ,    "JetPhiScaleFactor/F");                                           
  // TreeSemiLept->Branch("JetMatchedGenJetDR"                   , & JetMatchedGenJetDR                ,    "JetMatchedGenJetDR/F");                                            
  TreeSemiLept->Branch("JetMatchedGenJetPt"                   , & JetMatchedGenJetPt                ,    "JetMatchedGenJetPt/F");                                            
  TreeSemiLept->Branch("JetMatchedGenJetMass"                 , & JetMatchedGenJetMass              ,    "JetMatchedGenJetMass/F"); 
                        
  TreeSemiLept->Branch("JetGenMatched_TopHadronic"            , & JetGenMatched_TopHadronic         ,    "JetGenMatched_TopHadronic/I");      
  TreeSemiLept->Branch("JetGenMatched_TopPt"                  , & JetGenMatched_TopPt               ,    "JetGenMatched_TopPt/F");      
  TreeSemiLept->Branch("JetGenMatched_TopEta"                 , & JetGenMatched_TopEta              ,    "JetGenMatched_TopEta/F");      
  TreeSemiLept->Branch("JetGenMatched_TopPhi"                 , & JetGenMatched_TopPhi              ,    "JetGenMatched_TopPhi/F");      
  TreeSemiLept->Branch("JetGenMatched_TopMass"                , & JetGenMatched_TopMass             ,    "JetGenMatched_TopMass/F");      
  TreeSemiLept->Branch("JetGenMatched_bPt"                    , & JetGenMatched_bPt                 ,    "JetGenMatched_bPt/F");      
  TreeSemiLept->Branch("JetGenMatched_WPt"                    , & JetGenMatched_WPt                 ,    "JetGenMatched_WPt/F");      
  TreeSemiLept->Branch("JetGenMatched_Wd1Pt"                  , & JetGenMatched_Wd1Pt               ,    "JetGenMatched_Wd1Pt/F");      
  TreeSemiLept->Branch("JetGenMatched_Wd2Pt"                  , & JetGenMatched_Wd2Pt               ,    "JetGenMatched_Wd2Pt/F");      
  TreeSemiLept->Branch("JetGenMatched_Wd1ID"                  , & JetGenMatched_Wd1ID               ,    "JetGenMatched_Wd1ID/F");      
  TreeSemiLept->Branch("JetGenMatched_Wd2ID"                  , & JetGenMatched_Wd2ID               ,    "JetGenMatched_Wd2ID/F");      
  TreeSemiLept->Branch("JetGenMatched_MaxDeltaRPartonTop"     , & JetGenMatched_MaxDeltaRPartonTop  ,    "JetGenMatched_MaxDeltaRPartonTop/F");      
  TreeSemiLept->Branch("JetGenMatched_MaxDeltaRWPartonTop"    , & JetGenMatched_MaxDeltaRWPartonTop ,    "JetGenMatched_MaxDeltaRWPartonTop/F");      
  TreeSemiLept->Branch("JetGenMatched_MaxDeltaRWPartonW"      , & JetGenMatched_MaxDeltaRWPartonW   ,    "JetGenMatched_MaxDeltaRWPartonW/F");      
  TreeSemiLept->Branch("JetGenMatched_DeltaR_t_b"             , & JetGenMatched_DeltaR_t_b          ,    "JetGenMatched_DeltaR_t_b/F");      
  TreeSemiLept->Branch("JetGenMatched_DeltaR_t_W"             , & JetGenMatched_DeltaR_t_W          ,    "JetGenMatched_DeltaR_t_W/F");      
  TreeSemiLept->Branch("JetGenMatched_DeltaR_t_Wd1"           , & JetGenMatched_DeltaR_t_Wd1        ,    "JetGenMatched_DeltaR_t_Wd1/F");      
  TreeSemiLept->Branch("JetGenMatched_DeltaR_t_Wd2"           , & JetGenMatched_DeltaR_t_Wd2        ,    "JetGenMatched_DeltaR_t_Wd2/F");      
  TreeSemiLept->Branch("JetGenMatched_DeltaR_W_b1"            , & JetGenMatched_DeltaR_W_b1         ,    "JetGenMatched_DeltaR_W_b1/F");      
  TreeSemiLept->Branch("JetGenMatched_DeltaR_W_Wd1"           , & JetGenMatched_DeltaR_W_Wd1        ,    "JetGenMatched_DeltaR_W_Wd1/F");      
  TreeSemiLept->Branch("JetGenMatched_DeltaR_W_Wd2"           , & JetGenMatched_DeltaR_W_Wd2        ,    "JetGenMatched_DeltaR_W_Wd2/F");      
  TreeSemiLept->Branch("JetGenMatched_DeltaR_Wd1_Wd2"         , & JetGenMatched_DeltaR_Wd1_Wd2      ,    "JetGenMatched_DeltaR_Wd1_Wd2/F");      
  TreeSemiLept->Branch("JetGenMatched_DeltaR_Wd1_b"           , & JetGenMatched_DeltaR_Wd1_b        ,    "JetGenMatched_DeltaR_Wd1_b/F");      
  TreeSemiLept->Branch("JetGenMatched_DeltaR_Wd2_b"           , & JetGenMatched_DeltaR_Wd2_b        ,    "JetGenMatched_DeltaR_Wd2_b/F");      
  TreeSemiLept->Branch("JetGenMatched_DeltaR_jet_t"           , & JetGenMatched_DeltaR_jet_t        ,    "JetGenMatched_DeltaR_jet_t/F");      
  TreeSemiLept->Branch("JetGenMatched_DeltaR_jet_W"           , & JetGenMatched_DeltaR_jet_W        ,    "JetGenMatched_DeltaR_jet_W/F");      
  TreeSemiLept->Branch("JetGenMatched_DeltaR_jet_b"           , & JetGenMatched_DeltaR_jet_b        ,    "JetGenMatched_DeltaR_jet_b/F");      
  TreeSemiLept->Branch("JetGenMatched_DeltaR_jet_Wd1"         , & JetGenMatched_DeltaR_jet_Wd1      ,    "JetGenMatched_DeltaR_jet_Wd1/F");      
  TreeSemiLept->Branch("JetGenMatched_DeltaR_jet_Wd2"         , & JetGenMatched_DeltaR_jet_Wd2      ,    "JetGenMatched_DeltaR_jet_Wd2/F");      
  TreeSemiLept->Branch("JetGenMatched_DeltaR_pup0_b"          , & JetGenMatched_DeltaR_pup0_b       ,    "JetGenMatched_DeltaR_pup0_b/F");      
  TreeSemiLept->Branch("JetGenMatched_DeltaR_pup0_Wd1"        , & JetGenMatched_DeltaR_pup0_Wd1     ,    "JetGenMatched_DeltaR_pup0_Wd1/F");      
  TreeSemiLept->Branch("JetGenMatched_DeltaR_pup0_Wd2"        , & JetGenMatched_DeltaR_pup0_Wd2     ,    "JetGenMatched_DeltaR_pup0_Wd2/F");      
  TreeSemiLept->Branch("JetGenMatched_DeltaR_pup1_b"          , & JetGenMatched_DeltaR_pup1_b       ,    "JetGenMatched_DeltaR_pup1_b/F");      
  TreeSemiLept->Branch("JetGenMatched_DeltaR_pup1_Wd1"        , & JetGenMatched_DeltaR_pup1_Wd1     ,    "JetGenMatched_DeltaR_pup1_Wd1/F");      
  TreeSemiLept->Branch("JetGenMatched_DeltaR_pup1_Wd2"        , & JetGenMatched_DeltaR_pup1_Wd2     ,    "JetGenMatched_DeltaR_pup1_Wd2/F");               
  TreeSemiLept->Branch("JetGenMatched_partonPt"               , & JetGenMatched_partonPt            ,    "JetGenMatched_partonPt/F");      
  TreeSemiLept->Branch("JetGenMatched_partonEta"              , & JetGenMatched_partonEta           ,    "JetGenMatched_partonEta/F");      
  TreeSemiLept->Branch("JetGenMatched_partonPhi"              , & JetGenMatched_partonPhi           ,    "JetGenMatched_partonPhi/F");      
  TreeSemiLept->Branch("JetGenMatched_partonMass"             , & JetGenMatched_partonMass          ,    "JetGenMatched_partonMass/F");      
  TreeSemiLept->Branch("JetGenMatched_partonID"               , & JetGenMatched_partonID            ,    "JetGenMatched_partonID/F");      
  TreeSemiLept->Branch("JetGenMatched_DeltaRjetParton"        , & JetGenMatched_DeltaRjetParton     ,    "JetGenMatched_DeltaRjetParton/F");      
  std::cout << "Setup semi-lept jets in tree" << std::endl;


  TreeSemiLept->Branch("SemiLeptMETpx"                        , & SemiLeptMETpx                     , "SemiLeptMETpx/F");
  TreeSemiLept->Branch("SemiLeptMETpy"                        , & SemiLeptMETpy                     , "SemiLeptMETpy/F");
  TreeSemiLept->Branch("SemiLeptMETpt"                        , & SemiLeptMETpt                     , "SemiLeptMETpt/F");
  TreeSemiLept->Branch("SemiLeptMETphi"                       , & SemiLeptMETphi                    , "SemiLeptMETphi/F");
  TreeSemiLept->Branch("SemiLeptMETsumET"                     , & SemiLeptMETsumET                  , "SemiLeptMETsumET/F");
  TreeSemiLept->Branch("SemiLeptNvtx"                         , & SemiLeptNvtx                      , "SemiLeptNvtx/F");
  TreeSemiLept->Branch("SemiLeptRho"                          , & SemiLeptRho                       , "SemiLeptRho/F");
  TreeSemiLept->Branch("SemiLeptEventWeight"                  , & SemiLeptEventWeight               , "SemiLeptEventWeight/F");
  TreeSemiLept->Branch("SemiLeptPUweight"                  , & SemiLeptPUweight               , "SemiLeptPUweight/F");
  TreeSemiLept->Branch("SemiLeptPUweight_MBup"                  , & SemiLeptPUweight_MBup               , "SemiLeptPUweight_MBup/F");
  TreeSemiLept->Branch("SemiLeptPUweight_MBdn"                  , & SemiLeptPUweight_MBdn               , "SemiLeptPUweight_MBdn/F");
       
  
  
 
  TreeSemiLept->Branch("SemiLeptGenTTmass"                    , & SemiLeptGenTTmass                 , "SemiLeptGenTTmass/F");
  
  TreeSemiLept->Branch("HTlep"                                , & HTlep                             , "HTlep/F");
  TreeSemiLept->Branch("ST"                                   , & ST                                , "ST/F");
  TreeSemiLept->Branch("ST_CorrDn"                            , & ST_CorrDn                         , "ST_CorrDn/F");
  TreeSemiLept->Branch("ST_CorrUp"                            , & ST_CorrUp                         , "ST_CorrUp/F");
  TreeSemiLept->Branch("ST_PtSmearNom"                        , & ST_PtSmearNom                     , "ST_PtSmearNom/F");
  TreeSemiLept->Branch("ST_PtSmearUp"                         , & ST_PtSmearUp                      , "ST_PtSmearUp/F");
  TreeSemiLept->Branch("ST_PtSmearDn"                         , & ST_PtSmearDn                      , "ST_PtSmearDn/F");
  
  TreeSemiLept->Branch("SemiLeptQ2weight_CorrDn"              , & SemiLeptQ2weight_CorrDn           , "SemiLeptQ2weight_CorrDn/F");
  TreeSemiLept->Branch("SemiLeptQ2weight_CorrUp"              , & SemiLeptQ2weight_CorrUp           , "SemiLeptQ2weight_CorrUp/F");
  TreeSemiLept->Branch("SemiLeptNNPDF3weight_CorrDn"          , & SemiLeptNNPDF3weight_CorrDn       , "SemiLeptNNPDF3weight_CorrDn/F");
  TreeSemiLept->Branch("SemiLeptNNPDF3weight_CorrUp"          , & SemiLeptNNPDF3weight_CorrUp       , "SemiLeptNNPDF3weight_CorrUp/F");
  TreeSemiLept->Branch("SemiLeptRunNum"                       , & SemiLeptRunNum                    , "SemiLeptRunNum/F");
  TreeSemiLept->Branch("SemiLeptLumiBlock"                    , & SemiLeptLumiBlock                 , "SemiLeptLumiBlock/F");
  TreeSemiLept->Branch("SemiLeptEventNum"                     , & SemiLeptEventNum                  , "SemiLeptEventNum/F");
  TreeSemiLept->Branch("SemiLeptPassMETFilters"               , & SemiLeptPassMETFilters            , "SemiLeptPassMETFilters/I");
 
  TreeSemiLept->Branch("AK4dRminPt"                           , & AK4dRminPt                        , "AK4dRminPt/F");  
  TreeSemiLept->Branch("AK4dRminEta"                          , & AK4dRminEta                       , "AK4dRminEta/F");  
  TreeSemiLept->Branch("AK4dRminPhi"                          , & AK4dRminPhi                       , "AK4dRminPhi/F");  
  TreeSemiLept->Branch("AK4dRminMass"                         , & AK4dRminMass                     , "AK4dRminMass/F");  
  TreeSemiLept->Branch("AK4dRminBdisc"                        , & AK4dRminBdisc                     , "AK4dRminBdisc/F");  
  TreeSemiLept->Branch("AK4dRminLep"                          , & AK4dRminLep                       , "AK4dRminLep/F");  
  TreeSemiLept->Branch("AK4BtagdRminPt"                       , & AK4BtagdRminPt                    , "AK4BtagdRminPt/F");  
  TreeSemiLept->Branch("AK4BtagdRminBdisc"                    , & AK4BtagdRminBdisc                 , "AK4BtagdRminBdisc/F");  
  TreeSemiLept->Branch("AK4BtagdRminLep"                      , & AK4BtagdRminLep                   , "AK4BtagdRminLep/F");  
  TreeSemiLept->Branch("LepHemiContainsAK4BtagLoose"          , & LepHemiContainsAK4BtagLoose       , "LepHemiContainsAK4BtagLoose/I");  
  TreeSemiLept->Branch("LepHemiContainsAK4BtagMedium"         , & LepHemiContainsAK4BtagMedium      , "LepHemiContainsAK4BtagMedium/I");  
  TreeSemiLept->Branch("LepHemiContainsAK4BtagTight"          , & LepHemiContainsAK4BtagTight       , "LepHemiContainsAK4BtagTight/I");  

  TreeSemiLept->Branch("LeptonPhi"                            , &  LeptonPhi                        , "LeptonPhi/F"); 
  TreeSemiLept->Branch("LeptonPt"                             , &  LeptonPt                         , "LeptonPt/F"); 
  TreeSemiLept->Branch("LeptonEta"                            , &  LeptonEta                        , "LeptonEta/F"); 
  TreeSemiLept->Branch("LeptonMass"                           , &  LeptonMass                       , "LeptonMass/F"); 
  TreeSemiLept->Branch("PtRel"                                , &  PtRel                            , "PtRel/F"); 
  TreeSemiLept->Branch("LeptonIsMu"                           , &  LeptonIsMu                       , "LeptonIsMu/I"); 
  TreeSemiLept->Branch("MuMedium"                             , &  MuMedium                         , "MuMedium/I"); 
  TreeSemiLept->Branch("MuTight"                              , &  MuTight                          , "MuTight/I"); 
  TreeSemiLept->Branch("DeltaRJetLep"                         , &  DeltaRJetLep                     , "DeltaRJetLep/F"); 
  TreeSemiLept->Branch("DeltaPhiJetLep"                       , &  DeltaPhiJetLep                   , "DeltaPhiJetLep/F"); 
  

  TreeSemiLept->Branch("MuIso"                                , &  MuIso                            , "MuIso/F"); 
  TreeSemiLept->Branch("Elecron_absiso"                       , &  Elecron_absiso                   , "Elecron_absiso/F"); 
  TreeSemiLept->Branch("Elecron_relIsoWithDBeta"              , &  Elecron_relIsoWithDBeta          , "Elecron_relIsoWithDBeta/F"); 
  TreeSemiLept->Branch("Elecron_absiso_EA"                    , &  Elecron_absiso_EA                , "Elecron_absiso_EA/F"); 
  TreeSemiLept->Branch("Elecron_relIsoWithEA"                 , &  Elecron_relIsoWithEA             , "Elecron_relIsoWithEA/F"); */

  //MY ADDITION
  TreeSemiLept->Branch("extra_muon", &extra_muon, "extra_muon/I");
  TreeSemiLept->Branch("extra_electron", &extra_electron, "extra_electron/I");
  TreeSemiLept->Branch("no_leptons", &no_leptons, "no_leptons/I");
  TreeSemiLept->Branch("loose_muon", &loose_muon, "loose_muon/I");
  TreeSemiLept->Branch("loose_electron", &loose_electron, "loose_electron/I");
  TreeSemiLept->Branch("no_jets", &no_jets, "no_jets/I");
  TreeSemiLept->Branch("lonesome_jet", &lonesome_jet, "lonesome_jet/I");
  TreeSemiLept->Branch("loose_jet", &loose_jet, "loose_jet/I");
  TreeSemiLept->Branch("extra_merged_jet", &extra_merged_jet, "extra_merged_jet/I");
  TreeSemiLept->Branch("met_small", &met_small, "met_small/I");
  TreeSemiLept->Branch("bad_dijet", &bad_dijet, "bad_dijet/I");
  TreeSemiLept->Branch("undefined_met", &undefined_met, "undefined_met/I");
  TreeSemiLept->Branch("no_merged_jet", &no_merged_jet, "no_merged_jet/I");

  TreeSemiLept->Branch("good_event", &good_event, "good_event/I");

  
  /*TreeSemiLept->Branch("anom_weight1", &anom_weight1, "anom_weight1/D");
  TreeSemiLept->Branch("anom_weight2", &anom_weight2, "anom_weight2/D");
  TreeSemiLept->Branch("anom_weight3", &anom_weight3, "anom_weight3/D");
  TreeSemiLept->Branch("anom_weight4", &anom_weight4, "anom_weight4/D");
  TreeSemiLept->Branch("anom_weight5", &anom_weight5, "anom_weight5/D");
  TreeSemiLept->Branch("anom_weight6", &anom_weight6, "anom_weight6/D");
  TreeSemiLept->Branch("anom_weight7", &anom_weight7, "anom_weight7/D");
  TreeSemiLept->Branch("anom_weight8", &anom_weight8, "anom_weight8/D");
  TreeSemiLept->Branch("anom_weight9", &anom_weight9, "anom_weight9/D");
  TreeSemiLept->Branch("anom_weight10", &anom_weight10, "anom_weight10/D");
  TreeSemiLept->Branch("anom_weight11", &anom_weight11, "anom_weight11/D");
  TreeSemiLept->Branch("anom_weight12", &anom_weight12, "anom_weight12/D");
  TreeSemiLept->Branch("anom_weight13", &anom_weight13, "anom_weight13/D");
  TreeSemiLept->Branch("anom_weight14", &anom_weight14, "anom_weight14/D");
  TreeSemiLept->Branch("anom_weight15", &anom_weight15, "anom_weight15/D");
  TreeSemiLept->Branch("anom_weight16", &anom_weight16, "anom_weight16/D");
  TreeSemiLept->Branch("anom_weight17", &anom_weight17, "anom_weight17/D");
  TreeSemiLept->Branch("anom_weight18", &anom_weight18, "anom_weight18/D");
  TreeSemiLept->Branch("anom_weight19", &anom_weight19, "anom_weight19/D");
  TreeSemiLept->Branch("anom_weight20", &anom_weight20, "anom_weight20/D");
  TreeSemiLept->Branch("anom_weight21", &anom_weight21, "anom_weight21/D");
  TreeSemiLept->Branch("anom_weight22", &anom_weight22, "anom_weight22/D");
  TreeSemiLept->Branch("anom_weight23", &anom_weight23, "anom_weight23/D");
  TreeSemiLept->Branch("anom_weight24", &anom_weight24, "anom_weight24/D");
  TreeSemiLept->Branch("anom_weight25", &anom_weight25, "anom_weight25/D");
  TreeSemiLept->Branch("anom_weight26", &anom_weight26, "anom_weight26/D");
  TreeSemiLept->Branch("anom_weight27", &anom_weight27, "anom_weight27/D");
  TreeSemiLept->Branch("anom_weight28", &anom_weight28, "anom_weight28/D");
  TreeSemiLept->Branch("anom_weight29", &anom_weight29, "anom_weight29/D");
  TreeSemiLept->Branch("anom_weight30", &anom_weight30, "anom_weight30/D");
  TreeSemiLept->Branch("anom_weight31", &anom_weight31, "anom_weight31/D");
  TreeSemiLept->Branch("anom_weight32", &anom_weight32, "anom_weight32/D");
  TreeSemiLept->Branch("anom_weight33", &anom_weight33, "anom_weight33/D");
  TreeSemiLept->Branch("anom_weight34", &anom_weight34, "anom_weight34/D");
  TreeSemiLept->Branch("anom_weight35", &anom_weight35, "anom_weight35/D");
  TreeSemiLept->Branch("anom_weight36", &anom_weight36, "anom_weight36/D");
  TreeSemiLept->Branch("anom_weight37", &anom_weight37, "anom_weight37/D");
  TreeSemiLept->Branch("anom_weight38", &anom_weight38, "anom_weight38/D");
  TreeSemiLept->Branch("anom_weight39", &anom_weight39, "anom_weight39/D");
  TreeSemiLept->Branch("anom_weight40", &anom_weight40, "anom_weight40/D");
  TreeSemiLept->Branch("anom_weight41", &anom_weight41, "anom_weight41/D");
  TreeSemiLept->Branch("anom_weight42", &anom_weight42, "anom_weight42/D");
  TreeSemiLept->Branch("anom_weight43", &anom_weight43, "anom_weight43/D");
  TreeSemiLept->Branch("anom_weight44", &anom_weight44, "anom_weight44/D");
  TreeSemiLept->Branch("anom_weight45", &anom_weight45, "anom_weight45/D");
  TreeSemiLept->Branch("anom_weight46", &anom_weight46, "anom_weight46/D");
  TreeSemiLept->Branch("anom_weight47", &anom_weight47, "anom_weight47/D");
  TreeSemiLept->Branch("anom_weight48", &anom_weight48, "anom_weight48/D");
  TreeSemiLept->Branch("anom_weight49", &anom_weight49, "anom_weight49/D");
  TreeSemiLept->Branch("anom_weight50", &anom_weight50, "anom_weight50/D");*/
  
  TreeSemiLept->Branch("Wplus_pt", &Wplus_pt, "Wplus_pt/F");
  TreeSemiLept->Branch("Wminus_pt", &Wminus_pt, "Wminus_pt/F");
  TreeSemiLept->Branch("Zneutral_pt", &Zneutral_pt, "Zneutral_pt/F");
  TreeSemiLept->Branch("Wplus_px", &Wplus_px, "Wplus_px/F");
  TreeSemiLept->Branch("Wminus_px", &Wminus_px, "Wminus_px/F");
  TreeSemiLept->Branch("Zneutral_px", &Zneutral_px, "Zneutral_px/F");
  TreeSemiLept->Branch("Wplus_py", &Wplus_py, "Wplus_py/F");
  TreeSemiLept->Branch("Wminus_py", &Wminus_py, "Wminus_py/F");
  TreeSemiLept->Branch("Zneutral_py", &Zneutral_py, "Zneutral_py/F");
  TreeSemiLept->Branch("Wplus_pz", &Wplus_pz, "Wplus_pz/F");
  TreeSemiLept->Branch("Wminus_pz", &Wminus_pz, "Wminus_pz/F");
  TreeSemiLept->Branch("Zneutral_pz", &Zneutral_pz, "Zneutral_pz/F");
  TreeSemiLept->Branch("Wplus_e", &Wplus_e, "Wplus_e/F");
  TreeSemiLept->Branch("Wminus_e", &Wminus_e, "Wminus_e/F");
  TreeSemiLept->Branch("Zneutral_e", &Zneutral_e, "Zneutral_e/F");
  TreeSemiLept->Branch("Wplus_theta", &Wplus_theta, "Wplus_theta/F");
  TreeSemiLept->Branch("Wminus_theta", &Wminus_theta, "Wminus_theta/F");
  TreeSemiLept->Branch("Zneutral_theta", &Zneutral_theta, "Zneutral_theta/F");
  TreeSemiLept->Branch("Wplus_eta", &Wplus_eta, "Wplus_eta/F");
  TreeSemiLept->Branch("Wminus_eta", &Wminus_eta, "Wminus_eta/F");
  TreeSemiLept->Branch("Zneutral_eta", &Zneutral_eta, "Zneutral_eta/F");
   TreeSemiLept->Branch("Wplus_et", &Wplus_et, "Wplus_et/F");
  TreeSemiLept->Branch("Wminus_et", &Wminus_et, "Wminus_et/F");
  TreeSemiLept->Branch("Zneutral_et", &Zneutral_et, "Zneutral_et/F");
  TreeSemiLept->Branch("Wplus_phi", &Wplus_phi, "Wplus_phi/F");
  TreeSemiLept->Branch("Wminus_phi", &Wminus_phi, "Wminus_phi/F");
  TreeSemiLept->Branch("Zneutral_phi", &Zneutral_phi, "Zneutral_phi/F");
  TreeSemiLept->Branch("Wplus_y", &Wplus_y, "Wplus_y/F");
  TreeSemiLept->Branch("Wminus_y", &Wminus_y, "Wminus_y/F");
  TreeSemiLept->Branch("Zneutral_y", &Zneutral_y, "Zneutral_y/F");

  //////////////////////PARTICLE KINEMATICS
  
  TreeSemiLept->Branch("quark_pt", &quark_pt, "quark_pt/F"); 
  TreeSemiLept->Branch("antiquark_pt", &antiquark_pt, "antiquark_pt/F");
  TreeSemiLept->Branch("quark_px", &quark_px, "quark_px/F");
  TreeSemiLept->Branch("antiquark_px", &antiquark_px, "antiquark_px/F");
  TreeSemiLept->Branch("quark_py", &quark_py, "quark_py/F");
  TreeSemiLept->Branch("antiquark_py", &antiquark_py, "antiquark_py/F");
  TreeSemiLept->Branch("quark_pz", &quark_pz, "quark_pz/F");
  TreeSemiLept->Branch("antiquark_pz", &antiquark_pz, "antiquark_pz/F");
  TreeSemiLept->Branch("quark_e", &quark_e, "quark_e/F");
  TreeSemiLept->Branch("antiquark_e", &antiquark_e, "antiquark_e/F");
  TreeSemiLept->Branch("quark_eta", &quark_eta, "quark_eta/F");
  TreeSemiLept->Branch("antiquark_eta", &antiquark_eta, "antiquark_eta/F");
  TreeSemiLept->Branch("quark_et", &quark_et, "quark_et/F");
  TreeSemiLept->Branch("antiquark_et", &antiquark_et, "antiquark_et/F");
  TreeSemiLept->Branch("quark_theta", &quark_theta, "quark_theta/F");
  TreeSemiLept->Branch("antiquark_theta", &antiquark_theta, "antiquark_theta/F");
  TreeSemiLept->Branch("quark_phi", &quark_phi, "quark_phi/F");
  TreeSemiLept->Branch("antiquark_phi", &antiquark_phi, "antiquark_phi/F");
  TreeSemiLept->Branch("quark_y", &quark_y, "quark_y/F");
  TreeSemiLept->Branch("antiquark_y", &antiquark_y, "antiquark_y/F");

  
  TreeSemiLept->Branch("lepton_pt", &lepton_pt, "lepton_pt/F");
  TreeSemiLept->Branch("neutrino_pt", &neutrino_pt, "neutrino_pt/F");
  TreeSemiLept->Branch("lepton_px", &lepton_px, "lepton_px/F");
  TreeSemiLept->Branch("neutrino_px", &neutrino_px, "neutrino_px/F");
  TreeSemiLept->Branch("lepton_py", &lepton_py, "lepton_py/F");
  TreeSemiLept->Branch("neutrino_py", &neutrino_py, "neutrino_py/F");
  TreeSemiLept->Branch("lepton_pz", &lepton_pz, "lepton_pz/F");
  TreeSemiLept->Branch("neutrino_pz", &neutrino_pz, "neutrino_pz/F");
  TreeSemiLept->Branch("lepton_e", &lepton_e, "lepton_e/F");
  TreeSemiLept->Branch("neutrino_e", &neutrino_e, "neutrino_e/F");
  TreeSemiLept->Branch("lepton_eta", &lepton_eta, "lepton_eta/F");
  TreeSemiLept->Branch("neutrino_eta", &neutrino_eta, "neutrino_eta/F");
  TreeSemiLept->Branch("lepton_et", &lepton_et, "lepton_et/F");
  TreeSemiLept->Branch("neutrino_et", &neutrino_et, "neutrino_et/F");
  TreeSemiLept->Branch("lepton_theta", &lepton_theta, "lepton_theta/F");
  TreeSemiLept->Branch("neutrino_theta", &neutrino_theta, "neutrino_theta/F");
  TreeSemiLept->Branch("lepton_phi", &lepton_phi, "lepton_phi/F");
  TreeSemiLept->Branch("neutrino_phi", &neutrino_phi, "neutrino_phi/F");
  TreeSemiLept->Branch("lepton_y", &lepton_y, "lepton_y/F");
  TreeSemiLept->Branch("neutrino_y", &neutrino_y, "neutrino_y/F");


  //TreeSemiLept->Branch("ak4genjet0_pt", &ak4genjet0_pt, "ak4genjet0_pt/F");
  //TreeSemiLept->Branch("ak4genjet1_pt", &ak4genjet1_pt, "ak4genjet1_pt/F");
  TreeSemiLept->Branch("ak8genjet0_pt", &ak8genjet0_pt, "ak8genjet0_pt/F");
  TreeSemiLept->Branch("ak8genjet1_pt", &ak8genjet1_pt, "ak8genjet1_pt/F");
  //TreeSemiLept->Branch("ak4genjet0_px", &ak4genjet0_px, "ak4genjet0_px/F");
  //TreeSemiLept->Branch("ak4genjet1_px", &ak4genjet1_px, "ak4genjet1_px/F");
  TreeSemiLept->Branch("ak8genjet0_px", &ak8genjet0_px, "ak8genjet0_px/F");
  TreeSemiLept->Branch("ak8genjet1_px", &ak8genjet1_px, "ak8genjet1_px/F");
  //TreeSemiLept->Branch("ak4genjet0_py", &ak4genjet0_py, "ak4genjet0_py/F");
  //TreeSemiLept->Branch("ak4genjet1_py", &ak4genjet1_py, "ak4genjet1_py/F");
  TreeSemiLept->Branch("ak8genjet0_py", &ak8genjet0_py, "ak8genjet0_py/F");
  TreeSemiLept->Branch("ak8genjet1_py", &ak8genjet1_py, "ak8genjet1_py/F");
  //TreeSemiLept->Branch("ak4genjet0_pz", &ak4genjet0_pz, "ak4genjet0_pz/F");
  //TreeSemiLept->Branch("ak4genjet1_pz", &ak4genjet1_pz, "ak4genjet1_pz/F");
  TreeSemiLept->Branch("ak8genjet0_pz", &ak8genjet0_pz, "ak8genjet0_pz/F");
  TreeSemiLept->Branch("ak8genjet1_pz", &ak8genjet1_pz, "ak8genjet1_pz/F");
  //TreeSemiLept->Branch("ak4genjet0_e", &ak4genjet0_e, "ak4genjet0_e/F");
  //TreeSemiLept->Branch("ak4genjet1_e", &ak4genjet1_e, "ak4genjet1_e/F");
  TreeSemiLept->Branch("ak8genjet0_e", &ak8genjet0_e, "ak8genjet0_e/F");
  TreeSemiLept->Branch("ak8genjet1_e", &ak8genjet1_e, "ak8genjet1_e/F");
  //TreeSemiLept->Branch("ak4genjet0_eta", &ak4genjet0_eta, "ak4genjet0_eta/F");
  //TreeSemiLept->Branch("ak4genjet1_eta", &ak4genjet1_eta, "ak4genjet1_eta/F");
  TreeSemiLept->Branch("ak8genjet0_eta", &ak8genjet0_eta, "ak8genjet0_eta/F");
  TreeSemiLept->Branch("ak8genjet1_eta", &ak8genjet1_eta, "ak8genjet1_eta/F");
  //TreeSemiLept->Branch("ak4genjet0_et", &ak4genjet0_et, "ak4genjet0_et/F");
  //TreeSemiLept->Branch("ak4genjet1_et", &ak4genjet1_et, "ak4genjet1_et/F");
  TreeSemiLept->Branch("ak8genjet0_et", &ak8genjet0_et, "ak8genjet0_et/F");
  TreeSemiLept->Branch("ak8genjet1_et", &ak8genjet1_et, "ak8genjet1_et/F");
  //TreeSemiLept->Branch("ak4genjet0_theta", &ak4genjet0_theta, "ak4genjet0_theta/F");
  //TreeSemiLept->Branch("ak4genjet1_theta", &ak4genjet1_theta, "ak4genjet1_theta/F");
  TreeSemiLept->Branch("ak8genjet0_theta", &ak8genjet0_theta, "ak8genjet0_theta/F");
  TreeSemiLept->Branch("ak8genjet1_theta", &ak8genjet1_theta, "ak8genjet1_theta/F");
  //TreeSemiLept->Branch("ak4genjet0_phi", &ak4genjet0_phi, "ak4genjet0_phi/F");
  //TreeSemiLept->Branch("ak4genjet1_phi", &ak4genjet1_phi, "ak4genjet1_phi/F");
  TreeSemiLept->Branch("ak8genjet0_phi", &ak8genjet0_phi, "ak8genjet0_phi/F");
  TreeSemiLept->Branch("ak8genjet1_phi", &ak8genjet1_phi, "ak8genjet1_phi/F");
  //TreeSemiLept->Branch("ak4genjet0_y", &ak4genjet0_y, "ak4genjet0_y/F");
  //TreeSemiLept->Branch("ak4genjet1_y", &ak4genjet1_y, "ak4genjet1_y/F");
  TreeSemiLept->Branch("ak8genjet0_y", &ak8genjet0_y, "ak8genjet0_y/F");
  TreeSemiLept->Branch("ak8genjet1_y", &ak8genjet1_y, "ak8genjet1_y/F");

  
  ///RECO PARTICLE KINEMATICS

  TreeSemiLept->Branch("reco_muon_pt", &reco_muon_pt, "reco_muon_pt/F");
  TreeSemiLept->Branch("reco_electron_pt", &reco_electron_pt, "reco_electron_pt/F"); 
  TreeSemiLept->Branch("reco_muon_px", &reco_muon_px, "reco_muon_px/F");
  TreeSemiLept->Branch("reco_electron_px", &reco_electron_px, "reco_electron_px/F");
  TreeSemiLept->Branch("reco_muon_py", &reco_muon_py, "reco_muon_py/F");
  TreeSemiLept->Branch("reco_electron_py", &reco_electron_py, "reco_electron_py/F");
  TreeSemiLept->Branch("reco_muon_pz", &reco_muon_pz, "reco_muon_pz/F");
  TreeSemiLept->Branch("reco_electron_pz", &reco_electron_pz, "reco_electron_pz/F");
  TreeSemiLept->Branch("reco_muon_e", &reco_muon_e, "reco_muon_e/F");
  TreeSemiLept->Branch("reco_electron_e", &reco_electron_e, "reco_electron_e/F");
  TreeSemiLept->Branch("reco_muon_eta", &reco_muon_eta, "reco_muon_eta/F");
  TreeSemiLept->Branch("reco_electron_eta", &reco_electron_eta, "reco_electron_eta/F");
  TreeSemiLept->Branch("reco_muon_et", &reco_muon_et, "reco_muon_et/F");
  TreeSemiLept->Branch("reco_electron_et", &reco_electron_et, "reco_electron_et/F");
  TreeSemiLept->Branch("reco_muon_theta", &reco_muon_theta, "reco_muon_theta/F");
  TreeSemiLept->Branch("reco_electron_theta", &reco_electron_theta, "reco_electron_theta/F");
  TreeSemiLept->Branch("reco_muon_phi", &reco_muon_phi, "reco_muon_phi/F");
  TreeSemiLept->Branch("reco_electron_phi", &reco_electron_phi, "reco_electron_phi/F");
  TreeSemiLept->Branch("reco_muon_y", &reco_muon_y, "reco_muon_y/F");
  TreeSemiLept->Branch("reco_electron_y", &reco_electron_y, "reco_electron_y/F");
 
  
  TreeSemiLept->Branch("reco_lepton_pt", &reco_lepton_pt, "reco_lepton_pt/F");
  TreeSemiLept->Branch("reco_lepton_px", &reco_lepton_px, "reco_lepton_px/F");
  TreeSemiLept->Branch("reco_lepton_py", &reco_lepton_py, "reco_lepton_py/F");
  TreeSemiLept->Branch("reco_lepton_pz", &reco_lepton_pz, "reco_lepton_pz/F");
  TreeSemiLept->Branch("reco_lepton_e", &reco_lepton_e, "reco_lepton_e/F");
  TreeSemiLept->Branch("reco_lepton_eta", &reco_lepton_eta, "reco_lepton_eta/F");
  TreeSemiLept->Branch("reco_lepton_et", &reco_lepton_et, "reco_lepton_et/F");
  TreeSemiLept->Branch("reco_lepton_theta", &reco_lepton_theta, "reco_lepton_theta/F");
  TreeSemiLept->Branch("reco_lepton_phi", &reco_lepton_phi, "reco_lepton_phi/F");
  TreeSemiLept->Branch("reco_lepton_y", &reco_lepton_y, "reco_lepton_y/F");
 

  
  TreeSemiLept->Branch("reco_met_pt", &reco_met_pt, "reco_met_pt/F");
  TreeSemiLept->Branch("reco_met_px", &reco_met_px, "reco_met_px/F");
  TreeSemiLept->Branch("reco_met_py", &reco_met_py, "reco_met_py/F");
  TreeSemiLept->Branch("reco_met_pz", &reco_met_pz, "reco_met_pz/F");
  TreeSemiLept->Branch("reco_met_e", &reco_met_e, "reco_met_e/F");
  TreeSemiLept->Branch("reco_met_eta", &reco_met_eta, "reco_met_eta/F");
   TreeSemiLept->Branch("reco_met_et", &reco_met_et, "reco_met_et/F");
  TreeSemiLept->Branch("reco_met_theta", &reco_met_theta, "reco_met_theta/F");
  TreeSemiLept->Branch("reco_met_phi", &reco_met_phi, "reco_met_phi/F");
  TreeSemiLept->Branch("reco_met_y", &reco_met_y, "reco_met_y/F");

  
  TreeSemiLept->Branch("reco_jet0_pt", &reco_jet0_pt, "reco_jet0_pt/F");
  TreeSemiLept->Branch("reco_jet1_pt", &reco_jet1_pt, "reco_jet1_pt/F");
  TreeSemiLept->Branch("ak8jet0_pt", &ak8jet0_pt, "ak8jet0_pt/F");
  TreeSemiLept->Branch("ak8jet1_pt", &ak8jet1_pt, "ak8jet1_pt/F");
  TreeSemiLept->Branch("ak8jetONLY_pt", &ak8jetONLY_pt, "ak8jetONLY_pt/F");
  //TreeSemiLept->Branch("pupjet0_pt", &pupjet0_pt, "pupjet0_pt/F");
  //TreeSemiLept->Branch("pupjet1_pt", &pupjet1_pt, "pupjet1_pt/F"); 
  //TreeSemiLept->Branch("subjet0_pt", &subjet0_pt, "subjet0_pt/F");
  //TreeSemiLept->Branch("subjet1_pt", &subjet1_pt, "subjet1_pt/F");

  
 
  TreeSemiLept->Branch("reco_jet0_px", &reco_jet0_px, "reco_jet0_px/F");
  TreeSemiLept->Branch("reco_jet1_px", &reco_jet1_px, "reco_jet1_px/F");
  TreeSemiLept->Branch("ak8jet0_px", &ak8jet0_px, "ak8jet0_px/F");
  TreeSemiLept->Branch("ak8jet1_px", &ak8jet1_px, "ak8jet1_px/F");
  TreeSemiLept->Branch("ak8jetONLY_px", &ak8jetONLY_px, "ak8jetONLY_px/F");
  //TreeSemiLept->Branch("pupjet0_px", &pupjet0_px, "pupjet0_px/F");
  //TreeSemiLept->Branch("pupjet1_px", &pupjet1_px, "pupjet1_px/F"); 
  //TreeSemiLept->Branch("subjet0_px", &subjet0_px, "subjet0_px/F");
  //TreeSemiLept->Branch("subjet1_px", &subjet1_px, "subjet1_px/F");

 
 
  TreeSemiLept->Branch("reco_jet0_py", &reco_jet0_py, "reco_jet0_py/F");
  TreeSemiLept->Branch("reco_jet1_py", &reco_jet1_py, "reco_jet1_py/F");
  TreeSemiLept->Branch("ak8jet0_py", &ak8jet0_py, "ak8jet0_py/F");
  TreeSemiLept->Branch("ak8jet1_py", &ak8jet1_py, "ak8jet1_py/F");
  TreeSemiLept->Branch("ak8jetONLY_py", &ak8jetONLY_py, "ak8jetONLY_py/F");
  //TreeSemiLept->Branch("pupjet0_py", &pupjet0_py, "pupjet0_py/F");
  //TreeSemiLept->Branch("pupjet1_py", &pupjet1_py, "pupjet1_py/F"); 
  //TreeSemiLept->Branch("subjet0_py", &subjet0_py, "subjet0_py/F");
  //TreeSemiLept->Branch("subjet1_py", &subjet1_py, "subjet1_py/F");

  

 
  TreeSemiLept->Branch("reco_jet0_pz", &reco_jet0_pz, "reco_jet0_pz/F");
  TreeSemiLept->Branch("reco_jet1_pz", &reco_jet1_pz, "reco_jet1_pz/F");
  TreeSemiLept->Branch("ak8jet0_pz", &ak8jet0_pz, "ak8jet0_pz/F");
  TreeSemiLept->Branch("ak8jet1_pz", &ak8jet1_pz, "ak8jet1_pz/F");
  TreeSemiLept->Branch("ak8jetONLY_pz", &ak8jetONLY_pz, "ak8jetONLY_pz/F");
  //TreeSemiLept->Branch("pupjet0_pz", &pupjet0_pz, "pupjet0_pz/F");
  //TreeSemiLept->Branch("pupjet1_pz", &pupjet1_pz, "pupjet1_pz/F");
  //TreeSemiLept->Branch("subjet0_pz", &subjet0_pz, "subjet0_pz/F");
  //TreeSemiLept->Branch("subjet1_pz", &subjet1_pz, "subjet1_pz/F");

  


  TreeSemiLept->Branch("reco_jet0_e", &reco_jet0_e, "reco_jet0_e/F");
  TreeSemiLept->Branch("reco_jet1_e", &reco_jet1_e, "reco_jet1_e/F");
  TreeSemiLept->Branch("ak8jet0_e", &ak8jet0_e, "ak8jet0_e/F");
  TreeSemiLept->Branch("ak8jet1_e", &ak8jet1_e, "ak8jet1_e/F");
  TreeSemiLept->Branch("ak8jetONLY_e", &ak8jetONLY_e, "ak8jetONLY_e/F");
  //TreeSemiLept->Branch("pupjet0_e", &pupjet0_e, "pupjet0_e/F");
  //TreeSemiLept->Branch("pupjet1_e", &pupjet1_e, "pupjet1_e/F");
  //TreeSemiLept->Branch("subjet0_e", &subjet0_e, "subjet0_e/F");
  //TreeSemiLept->Branch("subjet1_e", &subjet1_e, "subjet1_e/F");



  TreeSemiLept->Branch("reco_jet0_eta", &reco_jet0_eta, "reco_jet0_eta/F");
  TreeSemiLept->Branch("reco_jet1_eta", &reco_jet1_eta, "reco_jet1_eta/F");
  TreeSemiLept->Branch("ak8jet0_eta", &ak8jet0_eta, "ak8jet0_eta/F");
  TreeSemiLept->Branch("ak8jet1_eta", &ak8jet1_eta, "ak8jet1_eta/F");
  TreeSemiLept->Branch("ak8jetONLY_eta", &ak8jetONLY_eta, "ak8jetONLY_eta/F");
  //TreeSemiLept->Branch("pupjet0_eta", &pupjet0_eta, "pupjet0_eta/F");
  //TreeSemiLept->Branch("pupjet1_eta", &pupjet1_eta, "pupjet1_eta/F");  
  //TreeSemiLept->Branch("subjet0_eta", &subjet0_eta, "subjet0_eta/F");
  //TreeSemiLept->Branch("subjet1_eta", &subjet1_eta, "subjet1_eta/F");

   TreeSemiLept->Branch("reco_jet0_et", &reco_jet0_et, "reco_jet0_et/F");
  TreeSemiLept->Branch("reco_jet1_et", &reco_jet1_et, "reco_jet1_et/F");
  TreeSemiLept->Branch("ak8jet0_et", &ak8jet0_et, "ak8jet0_et/F");
  TreeSemiLept->Branch("ak8jet1_et", &ak8jet1_et, "ak8jet1_et/F");
  TreeSemiLept->Branch("ak8jetONLY_et", &ak8jetONLY_et, "ak8jetONLY_et/F");
  //TreeSemiLept->Branch("pupjet0_et", &pupjet0_et, "pupjet0_et/F");
  //TreeSemiLept->Branch("pupjet1_et", &pupjet1_et, "pupjet1_et/F");  
  //TreeSemiLept->Branch("subjet0_et", &subjet0_et, "subjet0_et/F");
  //TreeSemiLept->Branch("subjet1_et", &subjet1_et, "subjet1_et/F");


  TreeSemiLept->Branch("reco_jet0_theta", &reco_jet0_theta, "reco_jet0_theta/F");
  TreeSemiLept->Branch("reco_jet1_theta", &reco_jet1_theta, "reco_jet1_theta/F");
  TreeSemiLept->Branch("ak8jet0_theta", &ak8jet0_theta, "ak8jet0_theta/F");
  TreeSemiLept->Branch("ak8jet1_theta", &ak8jet1_theta, "ak8jet1_theta/F");
  TreeSemiLept->Branch("ak8jetONLY_theta", &ak8jetONLY_theta, "ak8jetONLY_theta/F");
  //TreeSemiLept->Branch("pupjet0_theta", &pupjet0_theta, "pupjet0_theta/F");
  //TreeSemiLept->Branch("pupjet1_theta", &pupjet1_theta, "pupjet1_theta/F");
  //TreeSemiLept->Branch("subjet0_theta", &subjet0_theta, "subjet0_theta/F");
  //TreeSemiLept->Branch("subjet1_theta", &subjet1_theta, "subjet1_theta/F");



  TreeSemiLept->Branch("reco_jet0_phi", &reco_jet0_phi, "reco_jet0_phi/F");
  TreeSemiLept->Branch("reco_jet1_phi", &reco_jet1_phi, "reco_jet1_phi/F");
  TreeSemiLept->Branch("ak8jet0_phi", &ak8jet0_phi, "ak8jet0_phi/F");
  TreeSemiLept->Branch("ak8jet1_phi", &ak8jet1_phi, "ak8jet1_phi/F");
  TreeSemiLept->Branch("ak8jetONLY_phi", &ak8jetONLY_phi, "ak8jetONLY_phi/F");
  //TreeSemiLept->Branch("pupjet0_phi", &pupjet0_phi, "pupjet0_phi/F");
  //TreeSemiLept->Branch("pupjet1_phi", &pupjet1_phi, "pupjet1_phi/F");
  //TreeSemiLept->Branch("subjet0_phi", &subjet0_phi, "subjet0_phi/F");
  //TreeSemiLept->Branch("subjet1_phi", &subjet1_phi, "subjet1_phi/F");


  
  TreeSemiLept->Branch("reco_jet0_y", &reco_jet0_y, "reco_jet0_y/F");
  TreeSemiLept->Branch("reco_jet1_y", &reco_jet1_y, "reco_jet1_y/F");
  TreeSemiLept->Branch("ak8jet0_y", &ak8jet0_y, "ak8jet0_y/F");
  TreeSemiLept->Branch("ak8jet1_y", &ak8jet1_y, "ak8jet1_y/F");
  TreeSemiLept->Branch("ak8jetONLY_y", &ak8jetONLY_y, "ak8jetONLY_y/F");
  //TreeSemiLept->Branch("pupjet0_y", &pupjet0_y, "pupjet0_y/F");
  //TreeSemiLept->Branch("pupjet1_y", &pupjet1_y, "pupjet1_y/F");
  //TreeSemiLept->Branch("subjet0_y", &subjet0_y, "subjet0_y/F");
  //TreeSemiLept->Branch("subjet1_y", &subjet1_y, "subjet1_y/F");

  /////////////////////EXTRANIOUS JETS AND SUBJETS

  //TreeSemiLept->Branch("subjet2_pt", &subjet2_pt, "subjet2_pt/F");
  //TreeSemiLept->Branch("subjet3_pt", &subjet3_pt, "subjet3_pt/F");
  //TreeSemiLept->Branch("subjet4_pt", &subjet4_pt, "subjet4_pt/F");
  //TreeSemiLept->Branch("subjet5_pt", &subjet5_pt, "subjet5_pt/F");
  //TreeSemiLept->Branch("subjet6_pt", &subjet6_pt, "subjet6_pt/F");

  TreeSemiLept->Branch("ak8jet2_pt", &ak8jet2_pt, "ak8jet2_pt/F");
  TreeSemiLept->Branch("ak8jet3_pt", &ak8jet3_pt, "ak8jet3_pt/F");

  //////////

  //TreeSemiLept->Branch("subjet2_px", &subjet2_px, "subjet2_px/F");
  //TreeSemiLept->Branch("subjet3_px", &subjet3_px, "subjet3_px/F");
  //TreeSemiLept->Branch("subjet4_px", &subjet4_px, "subjet4_px/F");
  //TreeSemiLept->Branch("subjet5_px", &subjet5_px, "subjet5_px/F");
  //TreeSemiLept->Branch("subjet6_px", &subjet6_px, "subjet6_px/F");

  TreeSemiLept->Branch("ak8jet2_px", &ak8jet2_px, "ak8jet2_px/F");
  TreeSemiLept->Branch("ak8jet3_px", &ak8jet3_px, "ak8jet3_px/F");

  //////////

  //TreeSemiLept->Branch("subjet2_py", &subjet2_py, "subjet2_py/F");
  //TreeSemiLept->Branch("subjet3_py", &subjet3_py, "subjet3_py/F");
  //TreeSemiLept->Branch("subjet4_py", &subjet4_py, "subjet4_py/F");
  //TreeSemiLept->Branch("subjet5_py", &subjet5_py, "subjet5_py/F");
  //TreeSemiLept->Branch("subjet6_py", &subjet6_py, "subjet6_py/F");

  TreeSemiLept->Branch("ak8jet2_py", &ak8jet2_py, "ak8jet2_py/F");
  TreeSemiLept->Branch("ak8jet3_py", &ak8jet3_py, "ak8jet3_py/F");

  //////////

  //TreeSemiLept->Branch("subjet2_pz", &subjet2_pz, "subjet2_pz/F");
  //TreeSemiLept->Branch("subjet3_pz", &subjet3_pz, "subjet3_pz/F");
  //TreeSemiLept->Branch("subjet4_pz", &subjet4_pz, "subjet4_pz/F");
  //TreeSemiLept->Branch("subjet5_pz", &subjet5_pz, "subjet5_pz/F");
  //TreeSemiLept->Branch("subjet6_pz", &subjet6_pz, "subjet6_pz/F");

  TreeSemiLept->Branch("ak8jet2_pz", &ak8jet2_pz, "ak8jet2_pz/F");
  TreeSemiLept->Branch("ak8jet3_pz", &ak8jet3_pz, "ak8jet3_pz/F");

  //////////

  //TreeSemiLept->Branch("subjet2_e", &subjet2_e, "subjet2_e/F");
  //TreeSemiLept->Branch("subjet3_e", &subjet3_e, "subjet3_e/F");
  //TreeSemiLept->Branch("subjet4_e", &subjet4_e, "subjet4_e/F");
  //TreeSemiLept->Branch("subjet5_e", &subjet5_e, "subjet5_e/F");
  //TreeSemiLept->Branch("subjet6_e", &subjet6_e, "subjet6_e/F");

  TreeSemiLept->Branch("ak8jet2_e", &ak8jet2_e, "ak8jet2_e/F");
  TreeSemiLept->Branch("ak8jet3_e", &ak8jet3_e, "ak8jet3_e/F");

  //////////

  //TreeSemiLept->Branch("subjet2_eta", &subjet2_eta, "subjet2_eta/F");
  //TreeSemiLept->Branch("subjet3_eta", &subjet3_eta, "subjet3_eta/F");
  //TreeSemiLept->Branch("subjet4_eta", &subjet4_eta, "subjet4_eta/F");
  //TreeSemiLept->Branch("subjet5_eta", &subjet5_eta, "subjet5_eta/F");
  //TreeSemiLept->Branch("subjet6_eta", &subjet6_eta, "subjet6_eta/F");

  TreeSemiLept->Branch("ak8jet2_eta", &ak8jet2_eta, "ak8jet2_eta/F");
  TreeSemiLept->Branch("ak8jet3_eta", &ak8jet3_eta, "ak8jet3_eta/F");
  TreeSemiLept->Branch("ak8jet2_et", &ak8jet2_et, "ak8jet2_et/F");
  TreeSemiLept->Branch("ak8jet3_et", &ak8jet3_et, "ak8jet3_et/F");

  //////////

  //TreeSemiLept->Branch("subjet2_theta", &subjet2_theta, "subjet2_theta/F");
  //TreeSemiLept->Branch("subjet3_theta", &subjet3_theta, "subjet3_theta/F");
  //TreeSemiLept->Branch("subjet4_theta", &subjet4_theta, "subjet4_theta/F");
  //TreeSemiLept->Branch("subjet5_theta", &subjet5_theta, "subjet5_theta/F");
  //TreeSemiLept->Branch("subjet6_theta", &subjet6_theta, "subjet6_theta/F");

  TreeSemiLept->Branch("ak8jet2_theta", &ak8jet2_theta, "ak8jet2_theta/F");
  TreeSemiLept->Branch("ak8jet3_theta", &ak8jet3_theta, "ak8jet3_theta/F");

  //////////

  //TreeSemiLept->Branch("subjet2_phi", &subjet2_phi, "subjet2_phi/F");
  //TreeSemiLept->Branch("subjet3_phi", &subjet3_phi, "subjet3_phi/F");
  //TreeSemiLept->Branch("subjet4_phi", &subjet4_phi, "subjet4_phi/F");
  //TreeSemiLept->Branch("subjet5_phi", &subjet5_phi, "subjet5_phi/F");
  //TreeSemiLept->Branch("subjet6_phi", &subjet6_phi, "subjet6_phi/F");

  TreeSemiLept->Branch("ak8jet2_phi", &ak8jet2_phi, "ak8jet2_phi/F");
  TreeSemiLept->Branch("ak8jet3_phi", &ak8jet3_phi, "ak8jet3_phi/F");

  //////////

  //TreeSemiLept->Branch("subjet2_y", &subjet2_y, "subjet2_y/F");
  //TreeSemiLept->Branch("subjet3_y", &subjet3_y, "subjet3_y/F");
  //TreeSemiLept->Branch("subjet4_y", &subjet4_y, "subjet4_y/F");
  //TreeSemiLept->Branch("subjet5_y", &subjet5_y, "subjet5_y/F");
  //TreeSemiLept->Branch("subjet6_y", &subjet6_y, "subjet6_y/F");

  TreeSemiLept->Branch("ak8jet2_y", &ak8jet2_y, "ak8jet2_y/F");
  TreeSemiLept->Branch("ak8jet3_y", &ak8jet3_y, "ak8jet3_y/F");

  ///////////////////////////////ANGULAR VARIABLES
 
  
  //TreeSemiLept->Branch("reco_ak4_costheta1", &reco_ak4_costheta1, "reco_ak4_costheta1/F");
  //TreeSemiLept->Branch("reco_ak4_costheta2", &reco_ak4_costheta2, "reco_ak4_costheta2/F");
  //TreeSemiLept->Branch("reco_ak4_phi", &reco_ak4_phi, "reco_ak4_phi/F");
  //TreeSemiLept->Branch("reco_ak4_costhetastar", &reco_ak4_costhetastar, "reco_ak4_costhetastar/F");
  //TreeSemiLept->Branch("reco_ak4_phistar1", &reco_ak4_phistar1, "reco_ak4_phistar1/F");
  //TreeSemiLept->Branch("reco_ak4_phistar2", &reco_ak4_phistar2, "reco_ak4_phistar2/F");

  //TreeSemiLept->Branch("gen_ak4_costheta1", &gen_ak4_costheta1, "gen_ak4_costheta1/F");
  //TreeSemiLept->Branch("gen_ak4_costheta2", &gen_ak4_costheta2, "gen_ak4_costheta2/F");
  //TreeSemiLept->Branch("gen_ak4_phi", &gen_ak4_phi, "gen_ak4_phi/F");
  //TreeSemiLept->Branch("gen_ak4_costhetastar", &gen_ak4_costhetastar, "gen_ak4_costhetastar/F");
  //TreeSemiLept->Branch("gen_ak4_phistar1", &gen_ak4_phistar1, "gen_ak4_phistar1/F");
  //TreeSemiLept->Branch("gen_ak4_phistar2", &gen_ak4_phistar2, "gen_ak4_phistar2/F");

  TreeSemiLept->Branch("reco_costheta1", &reco_costheta1, "reco_costheta1/F");
  TreeSemiLept->Branch("reco_costheta2", &reco_costheta2, "reco_costheta2/F");
  TreeSemiLept->Branch("reco_phi", &reco_phi, "reco_phi/F");
  TreeSemiLept->Branch("reco_costhetastar", &reco_costhetastar, "reco_costhetastar/F");
  TreeSemiLept->Branch("reco_phistar1", &reco_phistar1, "reco_phistar1/F");
  TreeSemiLept->Branch("reco_phistar2", &reco_phistar2, "reco_phistar2/F");

  TreeSemiLept->Branch("gen_ak8_costheta1", &gen_ak8_costheta1, "gen_ak8_costheta1/F");
  TreeSemiLept->Branch("gen_ak8_costheta2", &gen_ak8_costheta2, "gen_ak8_costheta2/F");
  TreeSemiLept->Branch("gen_ak8_phi", &gen_ak8_phi, "gen_ak8_phi/F");
  TreeSemiLept->Branch("gen_ak8_costhetastar", &gen_ak8_costhetastar, "gen_ak8_costhetastar/F");
  TreeSemiLept->Branch("gen_ak8_phistar1", &gen_ak8_phistar1, "gen_ak8_phistar1/F");
  TreeSemiLept->Branch("gen_ak8_phistar2", &gen_ak8_phistar2, "gen_ak8_phistar2/F");

  //TreeSemiLept->Branch("reco_puppi_costheta1", &reco_puppi_costheta1, "reco_puppi_costheta1/F");
  //TreeSemiLept->Branch("reco_puppi_costheta2", &reco_puppi_costheta2, "reco_puppi_costheta2/F");
  //TreeSemiLept->Branch("reco_puppi_phi", &reco_puppi_phi, "reco_puppi_phi/F");
  //TreeSemiLept->Branch("reco_puppi_costhetastar", &reco_puppi_costhetastar, "reco_puppi_costhetastar/F");
  //TreeSemiLept->Branch("reco_puppi_phistar1", &reco_puppi_phistar1, "reco_puppi_phistar1/F");
  //TreeSemiLept->Branch("reco_puppi_phistar2", &reco_puppi_phistar2, "reco_puppi_phistar2/F");

  TreeSemiLept->Branch("lhe_gen_costheta1", &lhe_gen_costheta1, "lhe_gen_costheta1/F");
  TreeSemiLept->Branch("lhe_gen_costheta2", &lhe_gen_costheta2, "lhe_gen_costheta2/F");
  TreeSemiLept->Branch("lhe_gen_phi", &lhe_gen_phi, "lhe_gen_phi/F");
  TreeSemiLept->Branch("lhe_gen_costhetastar", &lhe_gen_costhetastar, "lhe_gen_costhetastar/F");
  TreeSemiLept->Branch("lhe_gen_phistar1", &lhe_gen_phistar1, "lhe_gen_phistar1/F");
  TreeSemiLept->Branch("lhe_gen_phistar2", &lhe_gen_phistar2, "lhe_gen_phistar2/F");

  //TreeSemiLept->Branch("subjet_costheta1", &subjet_costheta1, "subjet_costheta1/F");
  //TreeSemiLept->Branch("subjet_costheta2", &subjet_costheta2, "subjet_costheta2/F");
  //TreeSemiLept->Branch("subjet_phi", &subjet_phi, "subjet_phi/F");
  //TreeSemiLept->Branch("subjet_costhetastar", &subjet_costhetastar, "subjet_costhetastar/F");
  //TreeSemiLept->Branch("subjet_phistar1", &subjet_phistar1, "subjet_phistar1/F");
  //TreeSemiLept->Branch("subjet_phistar2", &subjet_phistar2, "subjet_phistar2/F");

  TreeSemiLept->Branch("goodjet_costheta1", &goodjet_costheta1, "goodjet_costheta1/F");
  TreeSemiLept->Branch("goodjet_costheta2", &goodjet_costheta2, "goodjet_costheta2/F");
  TreeSemiLept->Branch("goodjet_phi", &goodjet_phi, "goodjet_phi/F");
  TreeSemiLept->Branch("goodjet_costhetastar", &goodjet_costhetastar, "goodjet_costhetastar/F");
  TreeSemiLept->Branch("goodjet_phistar1", &goodjet_phistar1, "goodjet_phistar1/F");
  TreeSemiLept->Branch("goodjet_phistar2", &goodjet_phistar2, "goodjet_phistar2/F");

  //TreeSemiLept->Branch("goodsubjet_costheta1", &goodsubjet_costheta1, "goodsubjet_costheta1/F");
  //TreeSemiLept->Branch("goodsubjet_costheta2", &goodsubjet_costheta2, "goodsubjet_costheta2/F");
  //TreeSemiLept->Branch("goodsubjet_phi", &goodsubjet_phi, "goodsubjet_phi/F");
  //TreeSemiLept->Branch("goodsubjet_costhetastar", &goodsubjet_costhetastar, "goodsubjet_costhetastar/F");
  //TreeSemiLept->Branch("goodsubjet_phistar1", &goodsubjet_phistar1, "goodsubjet_phistar1/F");
  //TreeSemiLept->Branch("goodsubjet_phistar2", &goodsubjet_phistar2, "goodsubjet_phistar2/F");

  ///////////////////////////GENERATION
 
  TreeSemiLept->Branch("lep_gen", &lep_gen, "lep_gen/I");
  TreeSemiLept->Branch("reco_lep_gen", &reco_lep_gen, "reco_lep_gen/I");
  

  ///////////////////////////////DELTAS

  TreeSemiLept->Branch("lep_delta", &lep_delta, "lep_delta/F");
  TreeSemiLept->Branch("met_delta", &met_delta, "met_delta/F");

  TreeSemiLept->Branch("reco_jet_delta0", &reco_jet_delta0, "reco_jet_delta0/F");
  TreeSemiLept->Branch("reco_jet_delta1", &reco_jet_delta1, "reco_jet_delta1/F");

  TreeSemiLept->Branch("ak8jet_delta0", &ak8jet_delta0, "ak8jet_delta0/F");
  TreeSemiLept->Branch("ak8jet_delta1", &ak8jet_delta1, "ak8jet_delta1/F");
  TreeSemiLept->Branch("ak8jet_delta2", &ak8jet_delta2, "ak8jet_delta2/F");
  TreeSemiLept->Branch("ak8jet_delta3", &ak8jet_delta3, "ak8jet_delta3/F");

  //TreeSemiLept->Branch("pupjet_delta0", &pupjet_delta0, "pupjet_delta0/F");
  //TreeSemiLept->Branch("pupjet_delta1", &pupjet_delta1, "pupjet_delta1/F");


  //TreeSemiLept->Branch("ak4genjet_delta0", &ak4genjet_delta0, "ak4genjet_delta0/F");
  //TreeSemiLept->Branch("ak4genjet_delta1", &ak4genjet_delta1, "ak4genjet_delta1/F");

  TreeSemiLept->Branch("ak8genjet_delta0", &ak8genjet_delta0, "ak8genjet_delta0/F");
  TreeSemiLept->Branch("ak8genjet_delta1", &ak8genjet_delta1, "ak8genjet_delta1/F");

  TreeSemiLept->Branch("interquark_delta", &interquark_delta, "interquark_delta/F");

  //TreeSemiLept->Branch("subjet_delta0", &subjet_delta0, "subjet_delta0/F");
  //TreeSemiLept->Branch("subjet_delta1", &subjet_delta1, "subjet_delta1/F");

  //TreeSemiLept->Branch("subjet_delta2", &subjet_delta2, "subjet_delta2/F");
  //TreeSemiLept->Branch("subjet_delta3", &subjet_delta3, "subjet_delta3/F");
  //TreeSemiLept->Branch("subjet_delta4", &subjet_delta4, "subjet_delta4/F");
  //TreeSemiLept->Branch("subjet_delta5", &subjet_delta5, "subjet_delta5/F");
  //TreeSemiLept->Branch("subjet_delta6", &subjet_delta6, "subjet_delta6/F");

  ////////////////////////COUNTS

  TreeSemiLept->Branch("muon_count", &muon_count, "muon_count/I");
  TreeSemiLept->Branch("electron_count", &electron_count, "electron_count/I");
  
  //TreeSemiLept->Branch("reco_ak4_jet_count", &reco_ak4_jet_count, "reco_ak4_jet_count/I");
  TreeSemiLept->Branch("reco_ak8_jet_count", &reco_ak8_jet_count, "reco_ak8_jet_count/I");
  //TreeSemiLept->Branch("reco_puppi_jet_count", &reco_puppi_jet_count, "reco_puppi_jet_count/I");
  
  //TreeSemiLept->Branch("gen_ak4_jet_count", &gen_ak4_jet_count, "gen_ak4_jet_count/I");
  TreeSemiLept->Branch("gen_ak8_jet_count", &gen_ak8_jet_count, "gen_ak8_jet_count/I");

  //TreeSemiLept->Branch("subjet_count", &subjet_count, "subjet_count/I");

  /*TEMPLATE
    TreeSemiLept->Branch("value", &value, "value/F");
    TreeSemiLept->Branch("value", &value, "value/I");
    TreeSemiLept->Branch("value", &value, "value/O");
    TreeSemiLept->Branch("_pt", &_pt, "_pt/F");
    TreeSemiLept->Branch("_px", &_px, "_px/F");
    TreeSemiLept->Branch("_py", &_py, "_py/F");
    TreeSemiLept->Branch("_pz", &_pz, "_pz/F");
    TreeSemiLept->Branch("_e", &_e, "_e/F");
    TreeSemiLept->Branch("_eta", &_eta, "_eta/F");
    TreeSemiLept->Branch("_theta", &_theta, "_theta/F");
    TreeSemiLept->Branch("_phi", &_phi, "_phi/F");
    TreeSemiLept->Branch("_y", &_y, "_y/F");
    END TEMPLATE*/

  //4-BODY MASS SCHEME

  TreeSemiLept->Branch("masswidth", &masswidth, "masswidth/F");
  TreeSemiLept->Branch("good_jet_event", &good_jet_event, "good_jet_event/O");
   

  TreeSemiLept->Branch("goodjet0_pt", &goodjet0_pt, "goodjet0_pt/F");
  TreeSemiLept->Branch("goodjet0_px", &goodjet0_px, "goodjet0_px/F");
  TreeSemiLept->Branch("goodjet0_py", &goodjet0_py, "goodjet0_py/F");
  TreeSemiLept->Branch("goodjet0_pz", &goodjet0_pz, "goodjet0_pz/F");
  TreeSemiLept->Branch("goodjet0_e", &goodjet0_e, "goodjet0_e/F");
  TreeSemiLept->Branch("goodjet0_eta", &goodjet0_eta, "goodjet0_eta/F");
  TreeSemiLept->Branch("goodjet0_theta", &goodjet0_theta, "goodjet0_theta/F");
  TreeSemiLept->Branch("goodjet0_phi", &goodjet0_phi, "goodjet0_phi/F");
  TreeSemiLept->Branch("goodjet0_y", &goodjet0_y, "goodjet0_y/F");

  
  TreeSemiLept->Branch("goodjet1_pt", &goodjet1_pt, "goodjet1_pt/F");
  TreeSemiLept->Branch("goodjet1_px", &goodjet1_px, "goodjet1_px/F");
  TreeSemiLept->Branch("goodjet1_py", &goodjet1_py, "goodjet1_py/F");
  TreeSemiLept->Branch("goodjet1_pz", &goodjet1_pz, "goodjet1_pz/F");
  TreeSemiLept->Branch("goodjet1_e", &goodjet1_e, "goodjet1_e/F");
  TreeSemiLept->Branch("goodjet1_eta", &goodjet1_eta, "goodjet1_eta/F");
  TreeSemiLept->Branch("goodjet1_theta", &goodjet1_theta, "goodjet1_theta/F");
  TreeSemiLept->Branch("goodjet1_phi", &goodjet1_phi, "goodjet1_phi/F");
  TreeSemiLept->Branch("goodjet1_y", &goodjet1_y, "goodjet1_y/F");


  TreeSemiLept->Branch("goodjetonly_pt", &goodjetonly_pt, "goodjetonly_pt/F");
  TreeSemiLept->Branch("goodjetonly_px", &goodjetonly_px, "goodjetonly_px/F");
  TreeSemiLept->Branch("goodjetonly_py", &goodjetonly_py, "goodjetonly_py/F");
  TreeSemiLept->Branch("goodjetonly_pz", &goodjetonly_pz, "goodjetonly_pz/F");
  TreeSemiLept->Branch("goodjetonly_e", &goodjetonly_e, "goodjetonly_e/F");
  TreeSemiLept->Branch("goodjetonly_eta", &goodjetonly_eta, "goodjetonly_eta/F");
  TreeSemiLept->Branch("goodjetonly_theta", &goodjetonly_theta, "goodjetonly_theta/F");
  TreeSemiLept->Branch("goodjetonly_phi", &goodjetonly_phi, "goodjetonly_phi/F");
  TreeSemiLept->Branch("goodjetonly_y", &goodjetonly_y, "goodjetonly_y/F");

  ////////////////////

  //TreeSemiLept->Branch("goodsubjet0_pt", &goodsubjet0_pt, "goodsubjet0_pt/F");
  //TreeSemiLept->Branch("goodsubjet0_px", &goodsubjet0_px, "goodsubjet0_px/F");
  //TreeSemiLept->Branch("goodsubjet0_py", &goodsubjet0_py, "goodsubjet0_py/F");
  //TreeSemiLept->Branch("goodsubjet0_pz", &goodsubjet0_pz, "goodsubjet0_pz/F");
  //TreeSemiLept->Branch("goodsubjet0_e", &goodsubjet0_e, "goodsubjet0_e/F");
  //TreeSemiLept->Branch("goodsubjet0_eta", &goodsubjet0_eta, "goodsubjet0_eta/F");
  //TreeSemiLept->Branch("goodsubjet0_theta", &goodsubjet0_theta, "goodsubjet0_theta/F");
  //TreeSemiLept->Branch("goodsubjet0_phi", &goodsubjet0_phi, "goodsubjet0_phi/F");
  //TreeSemiLept->Branch("goodsubjet0_y", &goodsubjet0_y, "goodsubjet0_y/F");

  
  //TreeSemiLept->Branch("goodsubjet1_pt", &goodsubjet1_pt, "goodsubjet1_pt/F");
  //TreeSemiLept->Branch("goodsubjet1_px", &goodsubjet1_px, "goodsubjet1_px/F");
  //TreeSemiLept->Branch("goodsubjet1_py", &goodsubjet1_py, "goodsubjet1_py/F");
  //TreeSemiLept->Branch("goodsubjet1_pz", &goodsubjet1_pz, "goodsubjet1_pz/F");
  //TreeSemiLept->Branch("goodsubjet1_e", &goodsubjet1_e, "goodsubjet1_e/F");
  //TreeSemiLept->Branch("goodsubjet1_eta", &goodsubjet1_eta, "goodsubjet1_eta/F");
  //TreeSemiLept->Branch("goodsubjet1_theta", &goodsubjet1_theta, "goodsubjet1_theta/F");
  //TreeSemiLept->Branch("goodsubjet1_phi", &goodsubjet1_phi, "goodsubjet1_phi/F");
  //TreeSemiLept->Branch("goodsubjet1_y", &goodsubjet1_y, "goodsubjet1_y/F");


  
  std::cout << "Setup semi-lept tree" << std::endl;

  //DELTA TREE FILLING

  /* DeltaTree0 = new TTree("DeltaTree0", "DeltaTree0");
  DeltaTree1 = new TTree("DeltaTree1", "DeltaTree1");
  DeltaTree2 = new TTree("DeltaTree2", "DeltaTree2");
  DeltaTree3 = new TTree("DeltaTree3", "DeltaTree3");
  DeltaTree4 = new TTree("DeltaTree4", "DeltaTree4");
  DeltaTree5 = new TTree("DeltaTree5", "DeltaTree5");

  DeltaTree0->Branch("ak8_jet0lowq", &ak8_jet0lowq, "ak8_jet0lowq/F");
  DeltaTree0->Branch("ak8_jet0highq", &ak8_jet0highq, "ak8_jet0highq/F");

  DeltaTree0->Branch("ak8_jet1lowq", &ak8_jet1lowq, "ak8_jet1lowq/F");
  DeltaTree0->Branch("ak8_jet1highq", &ak8_jet1highq, "ak8_jet1highq/F");
      
  DeltaTree0->Branch("ak8_jet2lowq", &ak8_jet2lowq, "ak8_jet2lowq/F");
  DeltaTree0->Branch("ak8_jet2highq", &ak8_jet2highq, "ak8_jet2highq/F");

  DeltaTree0->Branch("ak8_jet3lowq", &ak8_jet3lowq, "ak8_jet3lowq/F");
  DeltaTree0->Branch("ak8_jet3highq", &ak8_jet3highq, "ak8_jet3highq/F");

  DeltaTree0->Branch("ak8_jet4lowq", &ak8_jet4lowq, "ak8_jet4lowq/F");
  DeltaTree0->Branch("ak8_jet4highq", &ak8_jet4highq, "ak8_jet4highq/F");

  DeltaTree0->Branch("ak8_jet5lowq", &ak8_jet5lowq, "ak8_jet5lowq/F");
  DeltaTree0->Branch("ak8_jet5highq", &ak8_jet5highq, "ak8_jet5highq/F");


  DeltaTree1->Branch("ak8_jet1lowq", &ak8_jet1lowq, "ak8_jet1lowq/F");
  DeltaTree1->Branch("ak8_jet1highq", &ak8_jet1highq, "ak8_jet1highq/F");
      
  DeltaTree1->Branch("ak8_jet2lowq", &ak8_jet2lowq, "ak8_jet2lowq/F");
  DeltaTree1->Branch("ak8_jet2highq", &ak8_jet2highq, "ak8_jet2highq/F");

  DeltaTree1->Branch("ak8_jet3lowq", &ak8_jet3lowq, "ak8_jet3lowq/F");
  DeltaTree1->Branch("ak8_jet3highq", &ak8_jet3highq, "ak8_jet3highq/F");

  DeltaTree1->Branch("ak8_jet4lowq", &ak8_jet4lowq, "ak8_jet4lowq/F");
  DeltaTree1->Branch("ak8_jet4highq", &ak8_jet4highq, "ak8_jet4highq/F");

  DeltaTree1->Branch("ak8_jet5lowq", &ak8_jet5lowq, "ak8_jet5lowq/F");
  DeltaTree1->Branch("ak8_jet5highq", &ak8_jet5highq, "ak8_jet5highq/F");


  DeltaTree2->Branch("ak8_jet2lowq", &ak8_jet2lowq, "ak8_jet2lowq/F");
  DeltaTree2->Branch("ak8_jet2highq", &ak8_jet2highq, "ak8_jet2highq/F");

  DeltaTree2->Branch("ak8_jet3lowq", &ak8_jet3lowq, "ak8_jet3lowq/F");
  DeltaTree2->Branch("ak8_jet3highq", &ak8_jet3highq, "ak8_jet3highq/F");

  DeltaTree2->Branch("ak8_jet4lowq", &ak8_jet4lowq, "ak8_jet4lowq/F");
  DeltaTree2->Branch("ak8_jet4highq", &ak8_jet4highq, "ak8_jet4highq/F");

  DeltaTree2->Branch("ak8_jet5lowq", &ak8_jet5lowq, "ak8_jet5lowq/F");
  DeltaTree2->Branch("ak8_jet5highq", &ak8_jet5highq, "ak8_jet5highq/F");


  DeltaTree3->Branch("ak8_jet3lowq", &ak8_jet3lowq, "ak8_jet3lowq/F");
  DeltaTree3->Branch("ak8_jet3highq", &ak8_jet3highq, "ak8_jet3highq/F");

  DeltaTree3->Branch("ak8_jet4lowq", &ak8_jet4lowq, "ak8_jet4lowq/F");
  DeltaTree3->Branch("ak8_jet4highq", &ak8_jet4highq, "ak8_jet4highq/F");

  DeltaTree3->Branch("ak8_jet5lowq", &ak8_jet5lowq, "ak8_jet5lowq/F");
  DeltaTree3->Branch("ak8_jet5highq", &ak8_jet5highq, "ak8_jet5highq/F");


  DeltaTree4->Branch("ak8_jet4lowq", &ak8_jet4lowq, "ak8_jet4lowq/F");
  DeltaTree4->Branch("ak8_jet4highq", &ak8_jet4highq, "ak8_jet4highq/F");

  DeltaTree4->Branch("ak8_jet5lowq", &ak8_jet5lowq, "ak8_jet5lowq/F");
  DeltaTree4->Branch("ak8_jet5highq", &ak8_jet5highq, "ak8_jet5highq/F");


  DeltaTree5->Branch("ak8_jet5lowq", &ak8_jet5lowq, "ak8_jet5lowq/F");
  DeltaTree5->Branch("ak8_jet5highq", &ak8_jet5highq, "ak8_jet5highq/F");*/

  

  std::cout << "Finished constructor" << std::endl;

}
//std::cout << "we defined the constructor.\n";

B2GTTbarTreeMaker::~B2GTTbarTreeMaker() {std::cout << "we are starting the destructor.\n";outfile.close();}

//std::cout "we defined the destructor.\n";
//
// member functions
//

// ------------ method called for each event  ------------

void B2GTTbarTreeMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  std::cout << "we are starting the analyze function.\n";
//#if 0
  using namespace edm;
  using namespace std;
  using namespace reco;
  using namespace pat;
  using namespace LHAPDF;

  std::cout << "verbose_ was " << verbose_ << endl;
  std::cout << "verboseGen_ was " << verboseGen_ << endl;

  verbose_ = true;
  verboseGen_ = true;
  
  if (verbose_)
    {
      cout << "----------------------------------------------------------------------------------------------------------" << endl;
      cout << "Analyze event " << iEvent.id().event() << " run " << iEvent.id().run() << " lumiblock " << iEvent.id().luminosityBlock() << endl;
    }
  // MARKER GEN PARTICLES
  //    .d8888b.  8888888888 888b    888     8888888b.                   888    d8b          888                   
  //   d88P  Y88b 888        8888b   888     888   Y88b                  888    Y8P          888                   
  //   888    888 888        88888b  888     888    888                  888                 888                   
  //   888        8888888    888Y88b 888     888   d88P  8888b.  888d888 888888 888  .d8888b 888  .d88b.  .d8888b  
  //   888  88888 888        888 Y88b888     8888888P"      "88b 888P"   888    888 d88P"    888 d8P  Y8b 88K      
  //   888    888 888        888  Y88888     888        .d888888 888     888    888 888      888 88888888 "Y8888b. 
  //   Y88b  d88P 888        888   Y8888     888        888  888 888     Y88b.  888 Y88b.    888 Y8b.          X88 
  //    "Y8888P88 8888888888 888    Y888     888        "Y888888 888      "Y888 888  "Y8888P 888  "Y8888   88888P' 
  //                                                                                                               
                                                                                                                                                                                                                        
  bool top1hadronic = false;
  bool top2hadronic = false;
  bool top1leptonic = false;
  bool top2leptonic = false;

  TLorentzVector t1_p4;
  TLorentzVector t2_p4;
  TLorentzVector finalt1_p4;
  TLorentzVector finalt2_p4;
  TLorentzVector b1_p4;
  TLorentzVector b2_p4;
  TLorentzVector W1_p4;
  TLorentzVector W2_p4;
  TLorentzVector W1d1_p4;
  TLorentzVector W1d2_p4;
  TLorentzVector W2d1_p4;
  TLorentzVector W2d2_p4;
  TLorentzVector resonantW1_p4;
  TLorentzVector resonantW2_p4;
  TLorentzVector Resonance_p4;
  TLorentzVector hardest_parton_hardScatterOutgoing_p4;
  TLorentzVector second_hardest_parton_hardScatterOutgoing_p4;

  double hardest_parton_hardScatterOutgoing_pt        = 0;
  double second_hardest_parton_hardScatterOutgoing_pt = 0;

  int parton1id = 0;
  int parton2id = 0;

  int W1d1_id = 0;
  int W1d2_id = 0;
  int W2d1_id = 0;
  int W2d2_id = 0;

  bool GenTruth_allhadronic = false;
  bool GenTruth_semileptonic = false;

  double deltaR_t1_t2       = 99; deltaR_t1_t2 = deltaR_t1_t2; //stop the ASININE compiler from throwing a damn fit when I comment out usage to debug things
  double deltaR_t1_b1       = 99; deltaR_t1_b1 = deltaR_t1_b1;
  double deltaR_t1_W1       = 99; deltaR_t1_W1 = deltaR_t1_W1;
  double deltaR_t1_W1d1     = 99; deltaR_t1_W1d1 = deltaR_t1_W1d1;
  double deltaR_t1_W1d2     = 99; deltaR_t1_W1d2 = deltaR_t1_W1d2;
  double deltaR_W1_b1       = 99; deltaR_W1_b1 = deltaR_W1_b1;
  double deltaR_W1_W1d1     = 99; deltaR_W1_W1d1 = deltaR_W1_W1d1;
  double deltaR_W1_W1d2     = 99; deltaR_W1_W1d2 = deltaR_W1_W1d2;
  double deltaR_W1d1_W1d2   = 99; deltaR_W1d1_W1d2 = deltaR_W1d1_W1d2;
  double deltaR_W1d1_b1     = 99; deltaR_W1d1_b1 = deltaR_W1d1_b1;
  double deltaR_W1d2_b1     = 99; deltaR_W1d2_b1 = deltaR_W1d2_b1;
  double deltaR_t2_b2       = 99; deltaR_t2_b2 =  deltaR_t2_b2;
  double deltaR_t2_W2       = 99; deltaR_t2_W2 = deltaR_t2_W2;
  double deltaR_t2_W2d1     = 99; deltaR_t2_W2d1 = deltaR_t2_W2d1;
  double deltaR_t2_W2d2     = 99; deltaR_t2_W2d2 = deltaR_t2_W2d2;
  double deltaR_W2_b2       = 99; deltaR_W2_b2 = deltaR_W2_b2;
  double deltaR_W2_W2d1     = 99; deltaR_W2_W2d1 = deltaR_W2_W2d1;
  double deltaR_W2_W2d2     = 99; deltaR_W2_W2d2 = deltaR_W2_W2d2;
  double deltaR_W2d1_W2d2   = 99; deltaR_W2d1_W2d2 = deltaR_W2d1_W2d2;
  double deltaR_W2d1_b2     = 99; deltaR_W2d1_b2 = deltaR_W2d1_b2;
  double deltaR_W2d2_b2     = 99; deltaR_W2d2_b2 = deltaR_W2d2_b2;

  double max_deltaR_parton_t1 = -1;
  double max_deltaR_parton_t2 = -1;
  double max_deltaR_Wparton_t1 = -1;
  double max_deltaR_Wparton_t2 = -1;
  double max_deltaR_Wparton_W1 = -1;
  double max_deltaR_Wparton_W2 = -1;

  double counttop = 0;
  if (!iEvent.isRealData() and runGenLoop_)
    {
      Handle<edm::View<reco::GenParticle> > genpart;
      iEvent.getByToken(prunedGenToken_, genpart);  

      // Classify the event based on the number of top quarks
      for (size_t i = 0; i < genpart->size(); i++)
	{
	  if (fabs((*genpart)[i].pdgId()) == 6 && (*genpart)[i].status() < 30 && (*genpart)[i].status() >= 20) counttop++;  // Z' events: status 22 = top from Z', status 52 with 2 daughters = the top that decays (after radiating a bunch of times)
	}
      if (verboseGen_) cout << "counttop " << counttop << endl;
   
      // Loop over all pruned gen particles and find the 4-vectors of the top, W, B, W duaghters
      double countW = 0;
      double countb = 0;
      for (size_t i = 0; i < genpart->size(); i++)
	{
	  int id        = (*genpart)[i].pdgId();
	  int status    = (*genpart)[i].status();
	  int ndau      = (*genpart)[i].numberOfDaughters();
	  double px     = (*genpart)[i].px();
	  double py     = (*genpart)[i].py();
	  double pz     = (*genpart)[i].pz();
	  double e      = (*genpart)[i].energy();
	  double m      = (*genpart)[i].mass();
	  double pt     = (*genpart)[i].pt();
	  double eta    = (*genpart)[i].eta();
	  double phi    = (*genpart)[i].phi();
	  //double nmothers = (*genpart)[i].numberOfMothers();

	  // Find the particles from the hard scatter (for QCD samples)
	  if (status == 23 && counttop == 0)
	    {
	      if (pt > hardest_parton_hardScatterOutgoing_pt)
		{
		  second_hardest_parton_hardScatterOutgoing_pt = hardest_parton_hardScatterOutgoing_pt;
		  second_hardest_parton_hardScatterOutgoing_p4 = hardest_parton_hardScatterOutgoing_p4;
		  hardest_parton_hardScatterOutgoing_pt = pt;
		  hardest_parton_hardScatterOutgoing_p4.SetPxPyPzE(px, py, pz, e);
		  parton1id = id;
		  if (verboseGen_) cout << "---------- pt>hardest_parton_hardScatterOutgoing_pt - parton1id = " << parton1id << endl;
		}
	      else if (pt > second_hardest_parton_hardScatterOutgoing_pt)
		{
		  second_hardest_parton_hardScatterOutgoing_pt = pt;
		  second_hardest_parton_hardScatterOutgoing_p4.SetPxPyPzE(px, py, pz, e); 
		  parton2id = id;
		  if (verboseGen_) cout << "---------- pt>second_hardest_parton_hardScatterOutgoing_pt - parton2id = " << parton2id << endl;
		}
	    }
	  // Find the the resonance mass for Z'
	  if (id > 1000000 && status == 22) 
	    {
	      Resonance_p4.SetPxPyPzE(px, py, pz, e); 
	      if (verboseGen_) cout << ".Resonant particle with mass " << m << endl; // RSGWW id = 5100039, Z' id = 6000047
	    }
	  // Get tops from hard subprocess (for ttbar samples)
	  if (id == 6 && status < 30 && status >= 20)
	    {
	      t1_p4.SetPxPyPzE(px, py, pz, e); 
	      parton1id = id;
	      if (verboseGen_) cout << "..top (hard)" << endl;//" with pt " << pt << " status " << status << " ndau " << ndau << " pt " << pt << " eta " << eta << " phi " << phi << " parton1id = " << parton1id << endl;
	    }
	  if (id == -6 && status < 30 && status >= 20)
	    {
	      t2_p4.SetPxPyPzE(px, py, pz, e); 
	      parton2id = id;
	      if (verboseGen_) cout << "..atop (hard)" << endl;//" with pt " << pt << " status " << status << " ndau " << ndau << " pt " << pt << " eta " << eta << " phi " << phi << " parton2id = " << parton2id << endl;
	    }

	  // Get the tops which decay - record b and W information
	  if (ndau == 2 && id == 6)
	    {
	      finalt1_p4.SetPxPyPzE(px, py, pz, e); 
	      if (verboseGen_) cout << "....two daughters top pt " << pt << " status " << status << " ndau " << ndau << " pt " << pt << " eta " << eta << " phi " << phi << endl;
	      for (int daught = 0; daught < 2; daught++)
		{
		  if (fabs((*genpart)[i].daughter(daught)->pdgId()) == 5)  b1_p4.SetPxPyPzE((*genpart)[i].daughter(daught)->px(), (*genpart)[i].daughter(daught)->py(), (*genpart)[i].daughter(daught)->pz(), (*genpart)[i].daughter(daught)->energy());
		  if (fabs((*genpart)[i].daughter(daught)->pdgId()) == 24) W1_p4.SetPxPyPzE((*genpart)[i].daughter(daught)->px(), (*genpart)[i].daughter(daught)->py(), (*genpart)[i].daughter(daught)->pz(), (*genpart)[i].daughter(daught)->energy());
		  if (verboseGen_) cout << "......top daughter ID " << (*genpart)[i].daughter(daught)->pdgId() << " pt " << (*genpart)[i].daughter(daught)->pt()  << endl;
		}
	    }
	  if (ndau == 2 && id == -6)
	    {
	      finalt2_p4.SetPxPyPzE(px, py, pz, e); 
	      if (verboseGen_) cout << "....two daughters atop pt " << pt << " status " << status << " ndau " << ndau << " pt " << pt << " eta " << eta << " phi " << phi << endl;
	      for (int daught = 0; daught < 2; daught++)
		{
		  if (fabs((*genpart)[i].daughter(daught)->pdgId()) == 5)  b2_p4.SetPxPyPzE((*genpart)[i].daughter(daught)->px(), (*genpart)[i].daughter(daught)->py(), (*genpart)[i].daughter(daught)->pz(), (*genpart)[i].daughter(daught)->energy());
		  if (fabs((*genpart)[i].daughter(daught)->pdgId()) == 24) W2_p4.SetPxPyPzE((*genpart)[i].daughter(daught)->px(), (*genpart)[i].daughter(daught)->py(), (*genpart)[i].daughter(daught)->pz(), (*genpart)[i].daughter(daught)->energy());
		  if (verboseGen_) cout << "......atop daughter ID " << (*genpart)[i].daughter(daught)->pdgId() << " pt " << (*genpart)[i].daughter(daught)->pt()  << endl;
		}
	    }

	  // Get the Ws which decay - record their daughter information
	  if (ndau == 2 && id == 24)
	    {
	      if (verboseGen_) cout << "....W+ with 2 daughters  id " << id << " status " << status << " ndau " << ndau << " pt " << pt << " eta " << eta << " phi " << phi << endl;
	      if (verboseGen_) cout << "......dd0 " << (*genpart)[i].daughter(0)->pdgId() << " ndau " << (*genpart)[i].daughter(0)->numberOfDaughters() << endl;
	      if (verboseGen_) cout << "......dd1 " << (*genpart)[i].daughter(1)->pdgId() << " ndau " << (*genpart)[i].daughter(1)->numberOfDaughters() << endl;
	      W1d1_p4.SetPxPyPzE((*genpart)[i].daughter(0)->px(), (*genpart)[i].daughter(0)->py(), (*genpart)[i].daughter(0)->pz(), (*genpart)[i].daughter(0)->energy());
	      W1d2_p4.SetPxPyPzE((*genpart)[i].daughter(1)->px(), (*genpart)[i].daughter(1)->py(), (*genpart)[i].daughter(1)->pz(), (*genpart)[i].daughter(1)->energy());
	      if (fabs((*genpart)[i].daughter(0)->pdgId()) < 6 && fabs((*genpart)[i].daughter(1)->pdgId()) < 6) top1hadronic = true;  
	      if (fabs((*genpart)[i].daughter(0)->pdgId()) <= 18 && fabs((*genpart)[i].daughter(0)->pdgId()) >= 11) top1leptonic = true;  
	      W1d1_id = (*genpart)[i].daughter(0)->pdgId();
	      W1d2_id = (*genpart)[i].daughter(1)->pdgId();
	    }
	  if (ndau == 2 && id == -24)
	    {
	      if (verboseGen_) cout << "....W- with 2 daughters  id " << id << " status " << status << " ndau " << ndau << " pt " << pt << " eta " << eta << " phi " << phi << endl;
	      if (verboseGen_) cout << "......dd0 " << (*genpart)[i].daughter(0)->pdgId() << " ndau " << (*genpart)[i].daughter(0)->numberOfDaughters() << endl;
	      if (verboseGen_) cout << "......dd1 " << (*genpart)[i].daughter(1)->pdgId() << " ndau " << (*genpart)[i].daughter(1)->numberOfDaughters() << endl;
	      W2d1_p4.SetPxPyPzE((*genpart)[i].daughter(0)->px(), (*genpart)[i].daughter(0)->py(), (*genpart)[i].daughter(0)->pz(), (*genpart)[i].daughter(0)->energy());
	      W2d2_p4.SetPxPyPzE((*genpart)[i].daughter(1)->px(), (*genpart)[i].daughter(1)->py(), (*genpart)[i].daughter(1)->pz(), (*genpart)[i].daughter(1)->energy());
	      if (fabs((*genpart)[i].daughter(0)->pdgId()) < 6 && fabs((*genpart)[i].daughter(1)->pdgId()) < 6) top2hadronic = true;  
	      if (fabs((*genpart)[i].daughter(0)->pdgId()) <= 18 && fabs((*genpart)[i].daughter(0)->pdgId()) >= 11) top2leptonic = true;  
	      W2d1_id = (*genpart)[i].daughter(0)->pdgId();
	      W2d2_id = (*genpart)[i].daughter(1)->pdgId();
	    }
	} // end genParticle loop

      // Generator truth all-hadronic or semileptonic
      if (top1hadronic && top2hadronic)  GenTruth_allhadronic = true;
      if (top1hadronic && !top2hadronic) GenTruth_semileptonic = true;
      if (!top1hadronic && top2hadronic) GenTruth_semileptonic = true;

      // Angles between particles
      deltaR_t1_t2       = t1_p4  .DeltaR(t2_p4);
      deltaR_t1_b1       = t1_p4  .DeltaR(b1_p4);
      deltaR_t1_W1       = t1_p4  .DeltaR(W1_p4);
      deltaR_t1_W1d1     = t1_p4  .DeltaR(W1d1_p4);
      deltaR_t1_W1d2     = t1_p4  .DeltaR(W1d2_p4);
      deltaR_W1_b1       = W1_p4  .DeltaR(b1_p4);
      deltaR_W1_W1d1     = W1_p4  .DeltaR(W1d1_p4);
      deltaR_W1_W1d2     = W1_p4  .DeltaR(W1d2_p4);
      deltaR_W1d1_W1d2   = W1d1_p4.DeltaR(W1d2_p4);
      deltaR_W1d1_b1     = W1d1_p4.DeltaR(b1_p4);
      deltaR_W1d2_b1     = W1d2_p4.DeltaR(b1_p4);
      deltaR_t2_b2       = t2_p4  .DeltaR(b2_p4);
      deltaR_t2_W2       = t2_p4  .DeltaR(W2_p4);
      deltaR_t2_W2d1     = t2_p4  .DeltaR(W2d1_p4);
      deltaR_t2_W2d2     = t2_p4  .DeltaR(W2d2_p4);
      deltaR_W2_b2       = W2_p4  .DeltaR(b2_p4);
      deltaR_W2_W2d1     = W2_p4  .DeltaR(W2d1_p4);
      deltaR_W2_W2d2     = W2_p4  .DeltaR(W2d2_p4);
      deltaR_W2d1_W2d2   = W2d1_p4.DeltaR(W2d2_p4);
      deltaR_W2d1_b2     = W2d1_p4.DeltaR(b2_p4);
      deltaR_W2d2_b2     = W2d2_p4.DeltaR(b2_p4);

      // Find top decay products which have the farthest angle from the top and from the W
      // max parton deltaR from t1
      if (deltaR_t1_b1   > max_deltaR_parton_t1) max_deltaR_parton_t1 = deltaR_t1_b1;
      if (deltaR_t1_W1d1 > max_deltaR_parton_t1) max_deltaR_parton_t1 = deltaR_t1_W1d1;
      if (deltaR_t1_W1d2 > max_deltaR_parton_t1) max_deltaR_parton_t1 = deltaR_t1_W1d2;

      // max parton deltaR from t2
      if (deltaR_t2_b2   > max_deltaR_parton_t2) max_deltaR_parton_t2 = deltaR_t2_b2;
      if (deltaR_t2_W2d1 > max_deltaR_parton_t2) max_deltaR_parton_t2 = deltaR_t2_W2d1;
      if (deltaR_t2_W2d2 > max_deltaR_parton_t2) max_deltaR_parton_t2 = deltaR_t2_W2d2;

      // max W parton deltaR from t1
      if (deltaR_t1_W1d1 > max_deltaR_Wparton_t1) max_deltaR_Wparton_t1 = deltaR_t1_W1d1;
      if (deltaR_t1_W1d2 > max_deltaR_Wparton_t1) max_deltaR_Wparton_t1 = deltaR_t1_W1d2;

      // max W parton deltaR from t2
      if (deltaR_t2_W2d1 > max_deltaR_Wparton_t2) max_deltaR_Wparton_t2 = deltaR_t2_W2d1;
      if (deltaR_t2_W2d2 > max_deltaR_Wparton_t2) max_deltaR_Wparton_t2 = deltaR_t2_W2d2;

      // max W parton deltaR from W1
      if (deltaR_W1_W1d1 > max_deltaR_Wparton_W1) max_deltaR_Wparton_W1 = deltaR_W1_W1d1;
      if (deltaR_W1_W1d2 > max_deltaR_Wparton_W1) max_deltaR_Wparton_W1 = deltaR_W1_W1d2;

      // max W parton deltaR from W2
      if (deltaR_W2_W2d1 > max_deltaR_Wparton_W2) max_deltaR_Wparton_W2 = deltaR_W2_W2d1;
      if (deltaR_W2_W2d2 > max_deltaR_Wparton_W2) max_deltaR_Wparton_W2 = deltaR_W2_W2d2;

      if (verboseGen_)
	{
	  cout << "second_hardest_parton_hardScatterOutgoing_pt " << second_hardest_parton_hardScatterOutgoing_pt      << endl;                
	  cout << "second_hardest_parton_hardScatterOutgoing_p4pt " << second_hardest_parton_hardScatterOutgoing_p4.Pt() << endl;                      
	  cout << "second_hardest_parton_hardScatterOutgoing_eta " << second_hardest_parton_hardScatterOutgoing_p4.Eta() << endl;                      
	  cout << "hardest_parton_hardScatterOutgoing_pt        " << hardest_parton_hardScatterOutgoing_pt             << endl;  
	  cout << "hardest_parton_hardScatterOutgoing_p4pt        " << hardest_parton_hardScatterOutgoing_p4.Pt()        << endl;       
	  cout << "hardest_parton_hardScatterOutgoing_eta        " << hardest_parton_hardScatterOutgoing_p4.Eta()        << endl;       
	  cout << "parton1id = " << parton1id << endl;
	  cout << "parton2id = " << parton1id << endl;

	  cout << "top1hadronic " << top1hadronic << endl;
	  cout << "top2hadronic " << top2hadronic << endl;
	  cout << "top1leptonic " << top1leptonic << endl;
	  cout << "top2leptonic " << top2leptonic << endl;
	  cout << "W1d1_id " << W1d1_id << endl;
	  cout << "W1d2_id " << W1d2_id << endl;
	  cout << "W2d1_id " << W2d1_id << endl;
	  cout << "W2d2_id " << W2d2_id << endl;

	  cout << "top1hadronic " << top1hadronic << endl;
	  cout << "top2hadronic " << top2hadronic << endl;
	  cout << "top1leptonic " << top1leptonic << endl;
	  cout << "top2leptonic " << top2leptonic << endl;
	  if (GenTruth_allhadronic)  cout << "allhadronic" << endl;
	  if (GenTruth_semileptonic) cout << "semileptonic" << endl;
    
	  cout << "t1_p4   Pt " << t1_p4  .Pt() << " Eta " << t1_p4  .Eta() << " Phi " << t1_p4  .Phi() << " M " << t1_p4  .M() << endl;
	  cout << "t2_p4   Pt " << t2_p4  .Pt() << " Eta " << t2_p4  .Eta() << " Phi " << t2_p4  .Phi() << " M " << t2_p4  .M() << endl;
	  cout << "b1_p4   Pt " << b1_p4  .Pt() << " Eta " << b1_p4  .Eta() << " Phi " << b1_p4  .Phi() << " M " << b1_p4  .M() << endl;
	  cout << "b2_p4   Pt " << b2_p4  .Pt() << " Eta " << b2_p4  .Eta() << " Phi " << b2_p4  .Phi() << " M " << b2_p4  .M() << endl;
	  cout << "W1_p4   Pt " << W1_p4  .Pt() << " Eta " << W1_p4  .Eta() << " Phi " << W1_p4  .Phi() << " M " << W1_p4  .M() << endl;
	  cout << "W2_p4   Pt " << W2_p4  .Pt() << " Eta " << W2_p4  .Eta() << " Phi " << W2_p4  .Phi() << " M " << W2_p4  .M() << endl;
	  cout << "W1d1_p4 Pt " << W1d1_p4.Pt() << " Eta " << W1d1_p4.Eta() << " Phi " << W1d1_p4.Phi() << " M " << W1d1_p4.M() << endl;
	  cout << "W1d2_p4 Pt " << W1d2_p4.Pt() << " Eta " << W1d2_p4.Eta() << " Phi " << W1d2_p4.Phi() << " M " << W1d2_p4.M() << endl;
	  cout << "W2d1_p4 Pt " << W2d1_p4.Pt() << " Eta " << W2d1_p4.Eta() << " Phi " << W2d1_p4.Phi() << " M " << W2d1_p4.M() << endl;
	  cout << "W2d2_p4 Pt " << W2d2_p4.Pt() << " Eta " << W2d2_p4.Eta() << " Phi " << W2d2_p4.Phi() << " M " << W2d2_p4.M() << endl;
	  cout << "resonantW1_p4   Pt " << resonantW1_p4.Pt() << " Eta " << resonantW1_p4.Eta() << " Phi " << resonantW1_p4.Phi() << " M " << resonantW1_p4.M() << endl;
	  cout << "resonantW2_p4   Pt " << resonantW2_p4.Pt() << " Eta " << resonantW2_p4.Eta() << " Phi " << resonantW2_p4.Phi() << " M " << resonantW2_p4.M() << endl;
    
	  cout << "deltaR_t1_t2       " << deltaR_t1_t2      << endl;
	  cout << "deltaR_t1_b1       " << deltaR_t1_b1      << endl; 
	  cout << "deltaR_t1_W1d1     " << deltaR_t1_W1d1    << endl; 
	  cout << "deltaR_t1_W1d2     " << deltaR_t1_W1d2    << endl; 
	  cout << "deltaR_t2_b2       " << deltaR_t2_b2      << endl; 
	  cout << "deltaR_t2_W2d1     " << deltaR_t2_W2d1    << endl; 
	  cout << "deltaR_t2_W2d2     " << deltaR_t2_W2d2    << endl; 
	  cout << "max_deltaR_parton_t1 " << max_deltaR_parton_t1 << endl;
	  cout << "max_deltaR_parton_t2 " << max_deltaR_parton_t2 << endl;
	  cout << "counttop " <<  counttop << " countW " << countW << " countb " << countb << endl;
	}
    }

  //   MARKER MY PARTICLES
  //   d88b.    .d88b   8b       d8          8888888b.                   888    d8b          888                   
  //   d88b.    .d888   88b     d88          888   Y88b                  888    Y8P          888                   
  //   d88b.9999.d888    88b   d88           888    888                  888                 888                   
  //   d88b. 99 .d888     88b d88            888   d88P  8888b.  888d888 888888 888  .d8888b 888  .d88b.  .d8888b  
  //   d88b.    .d888       8b8              8888888P"      "88b 888P"   888    888 d88P"    888 d8P  Y8b 88K      
  //   d88b.    .d888       888              888        .d888888 888     888    888 888      888 88888888 "Y8888b. 
  //   d88b.    .d888       Y8Y              888        888  888 888     Y88b.  888 Y88b.    888 Y8b.          X88 
  //   d88b.    .d888       888              888        "Y888888 888      "Y888 888  "Y8888P 888  "Y8888   88888P'         
  //                                                                                                                
 
  TLorentzVector Wplus;
  TLorentzVector Wminus;
  TLorentzVector Zneutral;
  TLorentzVector quark;
  TLorentzVector antiquark;
  TLorentzVector lepton;
  TLorentzVector neutrino;

  Wplus.SetPxPyPzE(-99.9,-99.9,-99.9,-99.9);
  Wminus.SetPxPyPzE(-99.9,-99.9,-99.9,-99.9);
  Zneutral.SetPxPyPzE(-99.9,-99.9,-99.9,-99.9);
  quark.SetPxPyPzE(-99.9,-99.9,-99.9,-99.9);
  antiquark.SetPxPyPzE(-99.9,-99.9,-99.9,-99.9);
  lepton.SetPxPyPzE(-99.9,-99.9,-99.9,-99.9);
  neutrino.SetPxPyPzE(-99.9,-99.9,-99.9,-99.9);

  //lepton.SetPxPyPzE(99.0, 99.0, 99.0, 99.0);

  //Int_t charge = 0;
  //charge = charge;//keep asinine compiler from throwing a fit if I later set it but don't use it.
  
  if (!iEvent.isRealData() and runGenLoop_) //if event is mc, not data
    {
      Handle<edm::View<reco::GenParticle>> genpart; //create pointer to .... particle array
      iEvent.getByToken(prunedGenToken_, genpart);  //connect said pointer to event

  //std::cout << "SEARCHING: About to enter for loop, looking for status of particles.\n";

  bool lepton_found = false;
  bool neutrino_found = false;
   
      for (size_t i = 0; i < genpart->size(); i++) //i is particle number, in a single event, sets, and resets, a bunch of quantities for the event (entire analyze function applied event-by-event)
	{
	  int id        = (*genpart)[i].pdgId(); //genpart[i] a particle
	  int status    = (*genpart)[i].status();
	  Float_t px     = (*genpart)[i].px();
	  Float_t py     = (*genpart)[i].py();
	  Float_t pz     = (*genpart)[i].pz();
	  Float_t e      = (*genpart)[i].energy();
          int ndau      = (*genpart)[i].numberOfDaughters();

	  std::string statsp = " ";
	  std::string daughtsp = " ";
	  std::string idsp = " ";
	  std::string xspace = " ";
	  std::string yspace = " ";
	  std::string zspace = " ";
	 
	  //EWK Gauge bosons

	  if (id ==  24)
	    {
	      if (status < 10) statsp = "  ";
	      if (status >= 10) statsp = " ";
	      if (ndau < 10) daughtsp = "  ";
	      if (ndau >= 10) daughtsp = " ";
	      if (id <= -10) idsp = " ";
	      if (id > -10 && id < 0) idsp = "  ";
	      if (id < 10 && id >= 0) idsp = "  +";
	      if (id >= 10) idsp = " +";
	      if (px < 0) xspace = " ";
	      if (px >= 0) xspace = " +";
	      if (py < 0) yspace = " ";
	      if (py >= 0) yspace = " +";
	      if (pz < 0) zspace = " ";
	      if (pz >= 0) zspace = " +";
	      Wplus.SetPxPyPzE(px, py, pz, e); 
	      //std::cout << "SEARCHING: status of Wplus  is" << statsp << status << " number of Wplus  daughters is" << daughtsp << ndau << " px is" << xspace << px << " py is" << yspace << py << " pz is" << zspace << pz << " e is " << e << endl;
	    }
	  
	  if (id == -24)
	    {
	      if (status < 10) statsp = "  ";
	      if (status >= 10) statsp = " ";
	      if (ndau < 10) daughtsp = "  ";
	      if (ndau >= 10) daughtsp = " ";
	      if (id <= -10) idsp = " ";
	      if (id > -10 && id < 0) idsp = "  ";
	      if (id < 10 && id >= 0) idsp = "  +";
	      if (id >= 10) idsp = " +";
	      if (px < 0) xspace = " ";
	      if (px >= 0) xspace = " +";
	      if (py < 0) yspace = " ";
	      if (py >= 0) yspace = " +";
	      if (pz < 0) zspace = " ";
	      if (pz >= 0) zspace = " +";
	      Wminus.SetPxPyPzE(px, py, pz, e);
	      //std::cout << "SEARCHING: status of Wminus is" << statsp << status << " number of Wminus daughters is" << daughtsp << ndau << " px is" << xspace << px << " py is" << yspace << py << " pz is" << zspace << pz << " e is " << e << endl;
	    }
	    
	  if (id ==  23)
	    {
	      if (status < 10) statsp = "  ";
	      if (status >= 10) statsp = " ";
	      if (ndau < 10) daughtsp = "  ";
	      if (ndau >= 10) daughtsp = " ";
	      if (id <= -10) idsp = " ";
	      if (id > -10 && id < 0) idsp = "  ";
	      if (id < 10 && id >= 0) idsp = "  +";
	      if (id >= 10) idsp = " +";
	      if (px < 0) xspace = " ";
	      if (px >= 0) xspace = " +";
	      if (py < 0) yspace = " ";
	      if (py >= 0) yspace = " +";
	      if (pz < 0) zspace = " ";
	      if (pz >= 0) zspace = " +";
	      Zneutral.SetPxPyPzE(px, py, pz, e);
	      //std::cout << "SEARCHING: status of Zneu   is" << statsp << status << " number of Zneu   daughters is" << daughtsp << ndau << " px is" << xspace << px << " py is" << yspace << py << " pz is" << zspace << pz << " e is " << e << endl;
	    }

	  //decay product fermions
	  
	  if (id == 11 && !lepton_found)
	    {
	      if (status < 10) statsp = "  ";
	      if (status >= 10) statsp = " ";
	      if (ndau < 10) daughtsp = "  ";
	      if (ndau >= 10) daughtsp = " ";
	      if (id <= -10) idsp = " ";
	      if (id > -10 && id < 0) idsp = "  ";
	      if (id < 10 && id >= 0) idsp = "  +";
	      if (id >= 10) idsp = " +";
	      if (px < 0) xspace = " ";
	      if (px >= 0) xspace = " +";
	      if (py < 0) yspace = " ";
	      if (py >= 0) yspace = " +";
	      if (pz < 0) zspace = " ";
	      if (pz >= 0) zspace = " +";
	      lepton.SetPxPyPzE(px, py, pz, e);
	      //if (id < 0) charge =  1;
	      //if (id > 0) charge = -1;
	      lepton_found = true;
	      lep_gen = 1;
	      //std::cout << "SEARCHING: status of  lepton is" << statsp << status << " number of  lepton daughters is" << daughtsp << ndau << " id is" << idsp << id << " px is" << xspace << px << " py is" << yspace << py << " pz is" << zspace << pz << " e is " << e << endl;
	    }


	  if (id == -11 && !lepton_found)
	    {
	      if (status < 10) statsp = "  ";
	      if (status >= 10) statsp = " ";
	      if (ndau < 10) daughtsp = "  ";
	      if (ndau >= 10) daughtsp = " ";
	      if (id <= -10) idsp = " ";
	      if (id > -10 && id < 0) idsp = "  ";
	      if (id < 10 && id >= 0) idsp = "  +";
	      if (id >= 10) idsp = " +";
	      if (px < 0) xspace = " ";
	      if (px >= 0) xspace = " +";
	      if (py < 0) yspace = " ";
	      if (py >= 0) yspace = " +";
	      if (pz < 0) zspace = " ";
	      if (pz >= 0) zspace = " +";
	      lepton.SetPxPyPzE(px, py, pz, e);
	      //if (id < 0) charge =  1;
	      //if (id > 0) charge = -1;
	      lepton_found = true;
	      lep_gen = -1;
	      //std::cout << "SEARCHING: status of  lepton is" << statsp << status << " number of  lepton daughters is" << daughtsp << ndau << " id is" << idsp << id << " px is" << xspace << px << " py is" << yspace << py << " pz is" << zspace << pz << " e is " << e << endl;
	    }

	  if (id == 13 && !lepton_found)
	    {
	      if (status < 10) statsp = "  ";
	      if (status >= 10) statsp = " ";
	      if (ndau < 10) daughtsp = "  ";
	      if (ndau >= 10) daughtsp = " ";
	      if (id <= -10) idsp = " ";
	      if (id > -10 && id < 0) idsp = "  ";
	      if (id < 10 && id >= 0) idsp = "  +";
	      if (id >= 10) idsp = " +";
	      if (px < 0) xspace = " ";
	      if (px >= 0) xspace = " +";
	      if (py < 0) yspace = " ";
	      if (py >= 0) yspace = " +";
	      if (pz < 0) zspace = " ";
	      if (pz >= 0) zspace = " +";
	      lepton.SetPxPyPzE(px, py, pz, e);
	      //if (id < 0) charge =  1;
	      //if (id > 0) charge = -1;
	      lepton_found = true;
	      lep_gen = 2;
	      //std::cout << "SEARCHING: status of  lepton is" << statsp << status << " number of  lepton daughters is" << daughtsp << ndau << " id is" << idsp << id << " px is" << xspace << px << " py is" << yspace << py << " pz is" << zspace << pz << " e is " << e << endl;
	    }

	  if (id == -13 && !lepton_found)
	    {
	      if (status < 10) statsp = "  ";
	      if (status >= 10) statsp = " ";
	      if (ndau < 10) daughtsp = "  ";
	      if (ndau >= 10) daughtsp = " ";
	      if (id <= -10) idsp = " ";
	      if (id > -10 && id < 0) idsp = "  ";
	      if (id < 10 && id >= 0) idsp = "  +";
	      if (id >= 10) idsp = " +";
	      if (px < 0) xspace = " ";
	      if (px >= 0) xspace = " +";
	      if (py < 0) yspace = " ";
	      if (py >= 0) yspace = " +";
	      if (pz < 0) zspace = " ";
	      if (pz >= 0) zspace = " +";
	      lepton.SetPxPyPzE(px, py, pz, e);
	      //if (id < 0) charge =  1;
	      //if (id > 0) charge = -1;
	      lepton_found = true;
	      lep_gen = -2;
	      //std::cout << "SEARCHING: status of  lepton is" << statsp << status << " number of  lepton daughters is" << daughtsp << ndau << " id is" << idsp << id << " px is" << xspace << px << " py is" << yspace << py << " pz is" << zspace << pz << " e is " << e << endl;
	    }

    	  if (id == 15 && !lepton_found)
	    {
	      if (status < 10) statsp = "  ";
	      if (status >= 10) statsp = " ";
	      if (ndau < 10) daughtsp = "  ";
	      if (ndau >= 10) daughtsp = " ";
	      if (id <= -10) idsp = " ";
	      if (id > -10 && id < 0) idsp = "  ";
	      if (id < 10 && id >= 0) idsp = "  +";
	      if (id >= 10) idsp = " +";
	      if (px < 0) xspace = " ";
	      if (px >= 0) xspace = " +";
	      if (py < 0) yspace = " ";
	      if (py >= 0) yspace = " +";
	      if (pz < 0) zspace = " ";
	      if (pz >= 0) zspace = " +";
	      lepton.SetPxPyPzE(px, py, pz, e);
	      //if (id < 0) charge =  1;
	      //if (id > 0) charge = -1;
	      lepton_found = true;
	      lep_gen = 3;
	      //std::cout << "SEARCHING: status of  lepton is" << statsp << status << " number of  lepton daughters is" << daughtsp << ndau << " id is" << idsp << id << " px is" << xspace << px << " py is" << yspace << py << " pz is" << zspace << pz << " e is " << e << endl;
	    }

	  if (id == -15 && !lepton_found)
	    {
	      if (status < 10) statsp = "  ";
	      if (status >= 10) statsp = " ";
	      if (ndau < 10) daughtsp = "  ";
	      if (ndau >= 10) daughtsp = " ";
	      if (id <= -10) idsp = " ";
	      if (id > -10 && id < 0) idsp = "  ";
	      if (id < 10 && id >= 0) idsp = "  +";
	      if (id >= 10) idsp = " +";
	      if (px < 0) xspace = " ";
	      if (px >= 0) xspace = " +";
	      if (py < 0) yspace = " ";
	      if (py >= 0) yspace = " +";
	      if (pz < 0) zspace = " ";
	      if (pz >= 0) zspace = " +";
	      lepton.SetPxPyPzE(px, py, pz, e);
	      //if (id < 0) charge =  1;
	      //if (id > 0) charge = -1;
	      lepton_found = true;
	      lep_gen = -3;
	      //std::cout << "SEARCHING: status of  lepton is" << statsp << status << " number of  lepton daughters is" << daughtsp << ndau << " id is" << idsp << id << " px is" << xspace << px << " py is" << yspace << py << " pz is" << zspace << pz << " e is " << e << endl;
	    }

	  
	  if (id == 12 && !neutrino_found)
	    {
	      if (status < 10) statsp = "  ";
	      if (status >= 10) statsp = " ";
	      if (ndau < 10) daughtsp = "  ";
	      if (ndau >= 10) daughtsp = " ";
	      if (id <= -10) idsp = " ";
	      if (id > -10 && id < 0) idsp = "  ";
	      if (id < 10 && id >= 0) idsp = "  +";
	      if (id >= 10) idsp = " +";
	      if (px < 0) xspace = " ";
	      if (px >= 0) xspace = " +";
	      if (py < 0) yspace = " ";
	      if (py >= 0) yspace = " +";
	      if (pz < 0) zspace= " ";
	      if (pz >= 0) zspace = " +";
	      neutrino.SetPxPyPzE(px, py, pz, e);
	      //if (id > 0) charge =  1;
	      //if (id < 0) charge = -1;
	      neutrino_found = true;
	      //lep_gen = -1;
	      //std::cout << "SEARCHING: status of neutrno is" << statsp << status << " number of neutrno daughters is" << daughtsp << ndau << " id is" << idsp << id << " px is" << xspace << px << " py is" << yspace << py << " pz is" << zspace << pz << " e is " << e << endl;
	    }

	  
	  if (id == -12 && !neutrino_found)
	    {
	      if (status < 10) statsp = "  ";
	      if (status >= 10) statsp = " ";
	      if (ndau < 10) daughtsp = "  ";
	      if (ndau >= 10) daughtsp = " ";
	      if (id <= -10) idsp = " ";
	      if (id > -10 && id < 0) idsp = "  ";
	      if (id < 10 && id >= 0) idsp = "  +";
	      if (id >= 10) idsp = " +";
	      if (px < 0) xspace = " ";
	      if (px >= 0) xspace = " +";
	      if (py < 0) yspace = " ";
	      if (py >= 0) yspace = " +";
	      if (pz < 0) zspace= " ";
	      if (pz >= 0) zspace = " +";
	      neutrino.SetPxPyPzE(px, py, pz, e);
	      //if (id > 0) charge =  1;
	      //if (id < 0) charge = -1;
	      neutrino_found = true;
	      //lep_gen = 1;
	      //std::cout << "SEARCHING: status of neutrno is" << statsp << status << " number of neutrno daughters is" << daughtsp << ndau << " id is" << idsp << id << " px is" << xspace << px << " py is" << yspace << py << " pz is" << zspace << pz << " e is " << e << endl;
	    }

	  
	    if (id == 14 && !neutrino_found)
	    {
	      if (status < 10) statsp = "  ";
	      if (status >= 10) statsp = " ";
	      if (ndau < 10) daughtsp = "  ";
	      if (ndau >= 10) daughtsp = " ";
	      if (id <= -10) idsp = " ";
	      if (id > -10 && id < 0) idsp = "  ";
	      if (id < 10 && id >= 0) idsp = "  +";
	      if (id >= 10) idsp = " +";
	      if (px < 0) xspace = " ";
	      if (px >= 0) xspace = " +";
	      if (py < 0) yspace = " ";
	      if (py >= 0) yspace = " +";
	      if (pz < 0) zspace= " ";
	      if (pz >= 0) zspace = " +";
	      neutrino.SetPxPyPzE(px, py, pz, e);
	      //if (id > 0) charge =  1;
	      //if (id < 0) charge = -1;
	      neutrino_found = true;
	      //lep_gen = -2;
	      //std::cout << "SEARCHING: status of neutrno is" << statsp << status << " number of neutrno daughters is" << daughtsp << ndau << " id is" << idsp << id << " px is" << xspace << px << " py is" << yspace << py << " pz is" << zspace << pz << " e is " << e << endl;
            }

	  
	    if (id == -14 && !neutrino_found)
	    {
	      if (status < 10) statsp = "  ";
	      if (status >= 10) statsp = " ";
	      if (ndau < 10) daughtsp = "  ";
	      if (ndau >= 10) daughtsp = " ";
	      if (id <= -10) idsp = " ";
	      if (id > -10 && id < 0) idsp = "  ";
	      if (id < 10 && id >= 0) idsp = "  +";
	      if (id >= 10) idsp = " +";
	      if (px < 0) xspace = " ";
	      if (px >= 0) xspace = " +";
	      if (py < 0) yspace = " ";
	      if (py >= 0) yspace = " +";
	      if (pz < 0) zspace= " ";
	      if (pz >= 0) zspace = " +";
	      neutrino.SetPxPyPzE(px, py, pz, e);
	      //if (id > 0) charge =  1;
	      //if (id < 0) charge = -1;
	      neutrino_found = true;
	      //lep_gen = 2;
	      //std::cout << "SEARCHING: status of neutrno is" << statsp << status << " number of neutrno daughters is" << daughtsp << ndau << " id is" << idsp << id << " px is" << xspace << px << " py is" << yspace << py << " pz is" << zspace << pz << " e is " << e << endl;
            }

	    
	    if (id == 16 && !neutrino_found)
	      {
		if (status < 10) statsp = "  ";
		if (status >= 10) statsp = " ";
		if (ndau < 10) daughtsp = "  ";
		if (ndau >= 10) daughtsp = " ";
		if (id <= -10) idsp = " ";
		if (id > -10 && id < 0) idsp = "  ";
		if (id < 10 && id >= 0) idsp = "  +";
		if (id >= 10) idsp = " +";
		if (px < 0) xspace = " ";
		if (px >= 0) xspace = " +";
		if (py < 0) yspace = " ";
		if (py >= 0) yspace = " +";
		if (pz < 0) zspace= " ";
		if (pz >= 0) zspace = " +";
		neutrino.SetPxPyPzE(px, py, pz, e);
		//if (id > 0) charge =  1;
		//if (id < 0) charge = -1;
		neutrino_found = true;
		//lep_gen = -3;
		//std::cout << "SEARCHING: status of neutrno is" << statsp << status << " number of neutrno daughters is" << daughtsp << ndau << " id is" << idsp << id << " px is" << xspace << px << " py is" << yspace << py << " pz is" << zspace << pz << " e is " << e << endl;
	      }

	    
	    if (id == -16 && !neutrino_found)
	      {
		if (status < 10) statsp = "  ";
		if (status >= 10) statsp = " ";
		if (ndau < 10) daughtsp = "  ";
		if (ndau >= 10) daughtsp = " ";
		if (id <= -10) idsp = " ";
		if (id > -10 && id < 0) idsp = "  ";
		if (id < 10 && id >= 0) idsp = "  +";
		if (id >= 10) idsp = " +";
		if (px < 0) xspace = " ";
		if (px >= 0) xspace = " +";
		if (py < 0) yspace = " ";
		if (py >= 0) yspace = " +";
		if (pz < 0) zspace= " ";
		if (pz >= 0) zspace = " +";
		neutrino.SetPxPyPzE(px, py, pz, e);
		//if (id > 0) charge =  1;
		//if (id < 0) charge = -1;
		neutrino_found = true;
		//lep_gen = 3;
		//std::cout << "SEARCHING: status of neutrno is" << statsp << status << " number of neutrno daughters is" << daughtsp << ndau << " id is" << idsp << id << " px is" << xspace << px << " py is" << yspace << py << " pz is" << zspace << pz << " e is " << e << endl;
	      }
	    

	    if (id >=  1 && id <=  6 && status == 23)
	      {
		if (status < 10) statsp = "  ";
		if (status >= 10) statsp = " ";
		if (ndau < 10) daughtsp = "  ";
		if (ndau >= 10) daughtsp = " ";
		if (id <= -10) idsp = " ";
		if (id > -10 && id < 0) idsp = "  ";
		if (id < 10 && id >= 0) idsp = "  +";
		if (id >= 10) idsp = " +";
		if (px < 0) xspace = " ";
		if (px >= 0) xspace = " +";
		if (py < 0) yspace = " ";
		if (py >= 0) yspace = " +";
		if (pz < 0) zspace = " ";
		if (pz >= 0) zspace = " +";
		quark.SetPxPyPzE(px, py, pz, e);
		//std::cout << "SEARCHING: status of   quark is" << statsp << status << " number of   quark daughters is" << daughtsp << ndau << " id is" << idsp << id << " px is" << xspace << px << " py is" << yspace << py << " pz is" << zspace << pz << " e is " << e << endl;
	    }
	      
	  if (id <= -1 && id >= -6 && status == 23)
	    {
	      if (status < 10) statsp = "  ";
	      if (status >= 10) statsp = " ";
	      if (ndau < 10) daughtsp = "  ";
	      if (ndau >= 10) daughtsp = " ";
	      if (id <= -10) idsp = " ";
	      if (id > -10 && id < 0) idsp = "  ";
	      if (id < 10 && id >= 0) idsp = "  +";
	      if (id >= 10) idsp = " +";
	      if (px < 0) xspace = " ";
	      if (px >= 0) xspace = " +";
	      if (py < 0) yspace = " ";
	      if (py >= 0) yspace = " +";
	      if (pz < 0) zspace = " ";
	      if (pz >= 0) zspace = " +";
	      antiquark.SetPxPyPzE(px, py, pz, e);
	      //std::cout << "SEARCHING: status of  aquark is" << statsp << status << " number of  aquark daughters is" << daughtsp << ndau << " id is" << idsp << id << " px is" << xspace << px << " py is" << yspace << py << " pz is" << zspace << pz << " e is " << e << endl;
	    }
	
	} // end genParticle loop
      
      //std::cout << "SEARCHING: Just exited for loop, looking for status of particles.\n";
      //marker
      
      Wplus_pt = Wplus.Pt();
      Wminus_pt = Wminus.Pt();
      Zneutral_pt = Zneutral.Pt();
      quark_pt = quark.Pt();
      antiquark_pt = antiquark.Pt();
      lepton_pt = lepton.Pt();
      neutrino_pt = neutrino.Pt();

      Wplus_px = Wplus.Px();
      Wminus_px = Wminus.Px();
      Zneutral_px = Zneutral.Px();
      quark_px = quark.Px();
      antiquark_px = antiquark.Px();
      lepton_px = lepton.Px();
      neutrino_px = neutrino.Px();

      Wplus_py = Wplus.Py();
      Wminus_py = Wminus.Py();
      Zneutral_py = Zneutral.Py();
      quark_py = quark.Py();
      antiquark_py = antiquark.Py();
      lepton_py = lepton.Py();
      neutrino_py = neutrino.Py();

      Wplus_pz = Wplus.Pz();
      Wminus_pz = Wminus.Pz();
      Zneutral_pz = Zneutral.Pz();
      quark_pz = quark.Pz();
      antiquark_pz = antiquark.Pz();
      lepton_pz = lepton.Pz();
      neutrino_pz = neutrino.Pz();

      Wplus_e = Wplus.E();
      Wminus_e = Wminus.E();
      Zneutral_e = Zneutral.E();
      quark_e = quark.E();
      antiquark_e = antiquark.E();
      lepton_e = lepton.E();
      neutrino_e = neutrino.E();

      Wplus_eta = Wplus.Eta();
      Wminus_eta = Wminus.Eta();
      Zneutral_eta = Zneutral.Eta();
      quark_eta = quark.Eta();
      antiquark_eta = antiquark.Eta();
      lepton_eta = lepton.Eta();
      neutrino_eta = neutrino.Eta();

      Wplus_et = Wplus.Et();
      Wminus_et = Wminus.Et();
      Zneutral_et = Zneutral.Et();
      quark_et = quark.Et();
      antiquark_et = antiquark.Et();
      lepton_et = lepton.Et();
      neutrino_et = neutrino.Et();

      Wplus_theta = Wplus.Theta();
      Wminus_theta = Wminus.Theta();
      Zneutral_theta = Zneutral.Theta();
      quark_theta = quark.Theta();
      antiquark_theta = antiquark.Theta();
      lepton_theta = lepton.Theta();
      neutrino_theta = neutrino.Theta();

      Wplus_phi = Wplus.Phi();
      Wminus_phi = Wminus.Phi();
      Zneutral_phi = Zneutral.Phi();
      quark_phi = quark.Phi();
      antiquark_phi = antiquark.Phi();
      lepton_phi = lepton.Phi();
      neutrino_phi = neutrino.Phi();

      Wplus_y = Wplus.Y();
      Wminus_y = Wminus.Y();
      Zneutral_y = Zneutral.Y();
      quark_y = quark.Y();
      antiquark_y = antiquark.Y();
      lepton_y = lepton.Y();
      neutrino_y = neutrino.Y();

      interquark_delta = quark.DeltaR(antiquark);

      if (lep_gen  > 0) calculateAngles(lepton, neutrino, quark, antiquark, costheta1, costheta2, phi, costhetastar, phistar1, phistar2);
      if (lep_gen  < 0) calculateAngles(neutrino, lepton, quark, antiquark, costheta1, costheta2, phi, costhetastar, phistar1, phistar2);

      if (lep_gen != 0) {lhe_gen_costheta1 = costheta1; lhe_gen_costheta2 = costheta2; lhe_gen_phi = phi; lhe_gen_costhetastar = costhetastar; lhe_gen_phistar1 = phistar1; lhe_gen_phistar2 = phistar2;}
    }
  
  // MARKER HLT
  // 888    888  888      88888888888 
  // 888    888  888          888     
  // 888    888  888          888     
  // 8888888888  888          888     
  // 888    888  888          888     
  // 888    888  888          888     
  // 888    888  888          888     
  // 888    888  88888888     888     
  // 
  
    edm::Handle<edm::TriggerResults> triggerBits;
    edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
    edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;

    iEvent.getByToken(triggerBits_, triggerBits);
    iEvent.getByToken(triggerObjects_, triggerObjects);
    iEvent.getByToken(triggerPrescales_, triggerPrescales);

    vector<string> trigsToRun;
    trigsToRun.push_back("HLT_PFHT300_v");
    trigsToRun.push_back("HLT_PFHT350_v");
    trigsToRun.push_back("HLT_PFHT400_v");
    trigsToRun.push_back("HLT_PFHT475_v");
    trigsToRun.push_back("HLT_PFHT600_v");
    trigsToRun.push_back("HLT_PFHT650_v");
    trigsToRun.push_back("HLT_PFHT800_v");
    trigsToRun.push_back("HLT_PFHT900_v");
    trigsToRun.push_back("HLT_PFJet320_v");
    trigsToRun.push_back("HLT_PFJet400_v");
    trigsToRun.push_back("HLT_PFJet450_v");
    trigsToRun.push_back("HLT_PFJet500_v");
    trigsToRun.push_back("HLT_AK8PFJet360_TrimMass30_v");
    trigsToRun.push_back("HLT_AK8PFHT700_TrimR0p1PT0p03Mass50_v");
    trigsToRun.push_back("HLT_AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV_p20_v");
    trigsToRun.push_back("HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20_v");
    trigsToRun.push_back("HLT_Mu45_eta2p1_v");
    trigsToRun.push_back("HLT_Mu50_v");
    trigsToRun.push_back("HLT_Mu55_v");
    trigsToRun.push_back("HLT_IsoMu22_eta2p1_v");
    trigsToRun.push_back("HLT_IsoMu24_v");
    trigsToRun.push_back("HLT_IsoMu27_v");
    trigsToRun.push_back("HLT_Mu30_eta2p1_PFJet150_PFJet50_v");
    trigsToRun.push_back("HLT_Mu40_eta2p1_PFJet200_PFJet50_v");
    trigsToRun.push_back("HLT_Ele32_eta2p1_WPTight_Gsf_v");
    trigsToRun.push_back("HLT_Ele35_WPLoose_Gsf_v");
    trigsToRun.push_back("HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50_v");
    trigsToRun.push_back("HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet140_v");
    trigsToRun.push_back("HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165_v");
    trigsToRun.push_back("HLT_Ele105_CaloIdVT_GsfTrkIdT_v");

    const int ntrigs = trigsToRun.size();
    if (verbose_) cout << "trigsToRun size " << ntrigs << endl;

    // do the same thing two different ways (to test)
    std::bitset<30> hltbit;
    vector<bool> trigAccept;

    AllHadTrigNames       ->clear();
    SemiLeptTrigNames     ->clear();
    AllHadTrigPrescales   ->clear();
    SemiLeptTrigPrescales ->clear();
    AllHadTrigPass        ->clear();
    SemiLeptTrigPass      ->clear();

    const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
    if (verbose_) std::cout << "\n === TRIGGER PATHS === " << std::endl;
    int counttrigs =0;
    // for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
    for (unsigned int j = 0; j < trigsToRun.size(); j++)
      {  // Running this loop first even though it is slow in order to preserve order (temporary solution)
	if (verbose_) cout << "try to find " << trigsToRun[j] << endl;
	for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i)
	  {
	    string name = names.triggerName(i);
	    if (verbose_) cout << " " << name << endl;
	    std::size_t found = name.find(trigsToRun[j]);
	    // cout << " Check: " << trigsToRun[j]  << " = " << name << " ?" << found << endl;
	    if (found != std::string::npos)
	      {
		int accept = triggerBits->accept(i);
		bool pass = false;
		if (accept == 1) pass = true;
		int prescale = triggerPrescales->getPrescaleForIndex(i);
        
		if (verbose_) std::cout << "  Found Trigger " << names.triggerName(i) << ", prescale " << triggerPrescales->getPrescaleForIndex(i) << ": " << (triggerBits->accept(i) ? "PASS" : "fail (or not run)") << std::endl; 
		trigAccept.push_back(pass);
		// trigStrings.push_back(name);
		// trigPrescale.push_back(prescale);
		AllHadTrigNames       ->push_back(name);
		SemiLeptTrigNames     ->push_back(name);
		AllHadTrigPrescales   ->push_back(prescale);
		SemiLeptTrigPrescales ->push_back(prescale);
		AllHadTrigPass        ->push_back(pass);
		SemiLeptTrigPass      ->push_back(pass);
		if (pass) hltbit[counttrigs]=1;  
		counttrigs++;
		break;
	      }
	  }
      }

    if (verbose_)
      {
	cout << "trig accept size " << trigAccept.size() << endl;
	for (unsigned int i = 0; i < trigAccept.size(); i++) {cout << trigAccept[trigAccept.size() - 1 - i];}
	cout << endl;
	cout << "hlt bit" << endl;
	cout << hltbit.to_string() << endl;
      }

    AllHadTrigAcceptBits   = hltbit.to_string();
    SemiLeptTrigAcceptBits = hltbit.to_string();
    
    // if (verbose_)
    // {
    //   cout << "trigStrings " << trigStrings.size() << endl;
    //   cout << "trigAccept " << trigAccept.size() << endl;
    //   cout << "trigPrescale " << trigPrescale.size() << endl;
    // }

    // MARKER MET NOISE FILTERS
    // 888b     d888 8888888888 88888888888     888b    888          d8b                       8888888888 d8b 888 888                              
    // 8888b   d8888 888            888         8888b   888          Y8P                       888        Y8P 888 888                              
    // 88888b.d88888 888            888         88888b  888                                    888            888 888                              
    // 888Y88888P888 8888888        888         888Y88b 888  .d88b.  888 .d8888b   .d88b.      8888888    888 888 888888  .d88b.  888d888 .d8888b  
    // 888 Y888P 888 888            888         888 Y88b888 d88""88b 888 88K      d8P  Y8b     888        888 888 888    d8P  Y8b 888P"   88K      
    // 888  Y8P  888 888            888         888  Y88888 888  888 888 "Y8888b. 88888888     888        888 888 888    88888888 888     "Y8888b. 
    // 888   "   888 888            888         888   Y8888 Y88..88P 888      X88 Y8b.         888        888 888 Y88b.  Y8b.     888          X88 
    // 888       888 8888888888     888         888    Y888  "Y88P"  888  88888P'  "Y8888      888        888 888  "Y888  "Y8888  888      88888P' 

    bool passMETfilters = false;
    if (iEvent.isRealData())
      {
	edm::Handle < edm::TriggerResults > metFilters;
	iEvent.getByToken(triggerResultsMETFilterToken_, metFilters);
	edm::TriggerNames const& filterTriggerNames = iEvent.triggerNames(*metFilters);

	int nMETfilters = metFilters->size();
	if (verbose_) cout << "nMETfilters " << nMETfilters << endl;

	//https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2#MiniAOD_8011_ICHEP_dataset
	// Flag_HBHENoiseFilter TO BE USED
	// Flag_HBHENoiseIsoFilter TO BE USED
	// Flag_EcalDeadCellTriggerPrimitiveFilter TO BE USED
	// Flag_goodVertices TO BE USED
	// Flag_eeBadScFilter TO BE USED
	// Flag_globalTightHalo2016Filter NEW TO BE USED
	// badMuon (NEW not available in 80X miniAOD v2, see snippet below) TO BE USED
	// badCharged hadron (NEW not available in 80X miniAOD v2, see snippet below) TO BE USED 

	// Recomendation:
	// primary vertex filter (available in miniAOD v2)
	// beam halo filter (NEW available in miniAOD v2)
	// HBHE noise filter (available in miniAOD v2)
	// HBHEiso noise filter (available in miniAOD v2)
	// ee badSC noise filter (available in miniAOD v2)
	// ECAL TP filter (available in miniAOD v2)
	// badMuon (NEW not available in miniAOD v2, see snippet below)
	// badCharged hadron (NEW not available in miniAOD v2, see snippet below) 

	vector<string> filterFlags;
	filterFlags.push_back("Flag_HBHENoiseFilter");
	filterFlags.push_back("Flag_HBHENoiseIsoFilter");
	filterFlags.push_back("Flag_EcalDeadCellTriggerPrimitiveFilter");
	filterFlags.push_back("Flag_goodVertices");
	filterFlags.push_back("Flag_eeBadScFilter");
	filterFlags.push_back("Flag_globalTightHalo2016Filter");

	unsigned int count_matched_accept = 0;
	for (int itrig = 0; itrig != nMETfilters; ++itrig)
	  {
	    std::string trigName = filterTriggerNames.triggerName(itrig);
	    bool accept = metFilters->accept(itrig);
	    if (verbose_) cout << trigName << "  " << accept;
	    if (std::find(filterFlags.begin(), filterFlags.end(), trigName) != filterFlags.end())
	      {
		if (verbose_) cout << "  -> matches filterFlags list (" << trigName << ")" << endl;
		if (accept) count_matched_accept++;
	      }
	    else if (verbose_) cout << endl;
	  }
	if (verbose_) cout << "filterFlags.size() " << filterFlags.size() << " count_matched_accept " << count_matched_accept << endl;
	if (count_matched_accept == filterFlags.size()) passMETfilters = true;
	if (verbose_) cout << "miniAOD Flags pass? " << passMETfilters << endl;
      }
    else passMETfilters = true;

    // RECO problemes -> apply to both data and MC
    Handle<bool> ifilterbadChCand;
    iEvent.getByToken(badChargedCandidateFilterToken_, ifilterbadChCand);
    bool  filterbadChCandidate = *ifilterbadChCand;
    if (verbose_) cout << "filterbadChCandidate " << filterbadChCandidate << endl;

    Handle<bool> ifilterbadPFMuon;
    iEvent.getByToken(badMuonFilterToken_, ifilterbadPFMuon);
    bool filterbadPFMuon = *ifilterbadPFMuon;
    if (verbose_) cout << "filterbadPFMuon " << filterbadPFMuon << endl;

    passMETfilters = passMETfilters && filterbadChCandidate && filterbadPFMuon;
    if (verbose_) cout << "passMETfilters = " << passMETfilters << endl;

    // MARKER JEC PAYLOADS
    //  888888 8888888888  .d8888b.      8888888b.                    888                        888          
    //    "88b 888        d88P  Y88b     888   Y88b                   888                        888          
    //     888 888        888    888     888    888                   888                        888          
    //     888 8888888    888            888   d88P  8888b.  888  888 888  .d88b.   8888b.   .d88888 .d8888b  
    //     888 888        888            8888888P"      "88b 888  888 888 d88""88b     "88b d88" 888 88K      
    //     888 888        888    888     888        .d888888 888  888 888 888  888 .d888888 888  888 "Y8888b. 
    //     88P 888        Y88b  d88P     888        888  888 Y88b 888 888 Y88..88P 888  888 Y88b 888      X88 
    //     888 8888888888  "Y8888P"      888        "Y888888  "Y88888 888  "Y88P"  "Y888888  "Y88888  88888P' 
    //   .d88P                                                    888                                         
    // .d88P"                                                Y8b d88P                                         
    //888P"                                                   "Y88P"                                          
    //

    // AK4chs JEC 
    std::vector<JetCorrectorParameters> vParAK4chs;
    for (std::vector<std::string>::const_iterator ipayload = jecPayloadsAK4chs_.begin(),
	   ipayloadEnd = jecPayloadsAK4chs_.end(); ipayload != ipayloadEnd - 1; ++ipayload)
      {
	JetCorrectorParameters pars(*ipayload);
	vParAK4chs.push_back(pars);
      }
    JetCorrectorAK4chs   = boost::shared_ptr<FactorizedJetCorrector>  (new FactorizedJetCorrector(vParAK4chs));
    JetCorrUncertAK4chs  = boost::shared_ptr<JetCorrectionUncertainty>(new JetCorrectionUncertainty(jecPayloadsAK4chs_.back()));
  
    // AK8chs JEC 
    std::vector<JetCorrectorParameters> vParAK8chs;
    for (std::vector<std::string>::const_iterator ipayload = jecPayloadsAK8chs_.begin(),
	   ipayloadEnd = jecPayloadsAK8chs_.end(); ipayload != ipayloadEnd - 1; ++ipayload)
      {
	JetCorrectorParameters pars(*ipayload);
	vParAK8chs.push_back(pars);
      }
    JetCorrectorAK8chs   = boost::shared_ptr<FactorizedJetCorrector>  (new FactorizedJetCorrector(vParAK8chs));
    JetCorrUncertAK8chs  = boost::shared_ptr<JetCorrectionUncertainty>(new JetCorrectionUncertainty(jecPayloadsAK8chs_.back()));

    // AK4pup JEC 
    std::vector<JetCorrectorParameters> vParAK4pup;
    for (std::vector<std::string>::const_iterator ipayload = jecPayloadsAK4pup_.begin(),
	   ipayloadEnd = jecPayloadsAK4pup_.end(); ipayload != ipayloadEnd - 1; ++ipayload)
      {
	JetCorrectorParameters pars(*ipayload);
	vParAK4pup.push_back(pars);
      }
    JetCorrectorAK4pup   = boost::shared_ptr<FactorizedJetCorrector>  (new FactorizedJetCorrector(vParAK4pup));
    JetCorrUncertAK4pup  = boost::shared_ptr<JetCorrectionUncertainty>(new JetCorrectionUncertainty(jecPayloadsAK4pup_.back()));
  
    // AK8pup JEC 
    std::vector<JetCorrectorParameters> vParAK8pup;
    for (std::vector<std::string>::const_iterator ipayload = jecPayloadsAK8pup_.begin(),
	   ipayloadEnd = jecPayloadsAK8pup_.end(); ipayload != ipayloadEnd - 1; ++ipayload)
      {
	JetCorrectorParameters pars(*ipayload);
	vParAK8pup.push_back(pars);
      }
    JetCorrectorAK8pup   = boost::shared_ptr<FactorizedJetCorrector>  (new FactorizedJetCorrector(vParAK8pup));
    JetCorrUncertAK8pup  = boost::shared_ptr<JetCorrectionUncertainty>(new JetCorrectionUncertainty(jecPayloadsAK8pup_.back()));
  
    // jet resolution scale factor from text files
    JME::JetResolutionScaleFactor jer_scaler;
    jer_scaler = JME::JetResolutionScaleFactor(jerSFtext_);


    // MARKER VERTICIES
    //  888     888                  888    d8b                            
    //  888     888                  888    Y8P                            
    //  888     888                  888                                   
    //  Y88b   d88P  .d88b.  888d888 888888 888  .d8888b  .d88b.  .d8888b  
    //   Y88b d88P  d8P  Y8b 888P"   888    888 d88P"    d8P  Y8b 88K      
    //    Y88o88P   88888888 888     888    888 888      88888888 "Y8888b. 
    //     Y888P    Y8b.     888     Y88b.  888 Y88b.    Y8b.          X88 
    //      Y8P      "Y8888  888      "Y888 888  "Y8888P  "Y8888   88888P' 
    //                                                                     
    //                                                                     
                                                                       
    edm::Handle<std::vector<reco::Vertex>> vertices;
    iEvent.getByToken(vtxToken_, vertices);
    int nvtx = vertices->size();
    if (vertices->empty()) return; // skip the event if no PV found
    const reco::Vertex &PV = vertices->front();  // save PV for tight muon ID

    // int nvtxgood =0;
    // for (VertexCollection::const_iterator vtx = vertices->begin(); vtx != vertices->end(); ++vtx, ++firstGoodVertexIdx)
    // {
    //   // Replace isFake() for miniAOD because it requires tracks and miniAOD vertices don't have tracks:
    //   // Vertex.h: bool isFake() const {return (chi2_==0 && ndof_ == 0 && tracks_.empty());}
    //   bool isFake = (vtx->chi2() == 0 && vtx->ndof() == 0);
    //   if (!isFake &&  vtx->ndof() >= 4. && vtx->position().Rho()<=2.0 && fabs(vtx->position().Z())<=24.0) nvtxgood++;
    // }
    // cout << "nvtx " << nvtx << " nvtxgood " << nvtxgood << endl;
    
    // MARKER PU WEIGHT
    //  8888888b.  888     888     888       888          d8b          888      888    
    //  888   Y88b 888     888     888   o   888          Y8P          888      888    
    //  888    888 888     888     888  d8b  888                       888      888    
    //  888   d88P 888     888     888 d888b 888  .d88b.  888  .d88b.  88888b.  888888 
    //  8888888P"  888     888     888d88888b888 d8P  Y8b 888 d88P"88b 888 "88b 888    
    //  888        888     888     88888P Y88888 88888888 888 888  888 888  888 888    
    //  888        Y88b. .d88P     8888P   Y8888 Y8b.     888 Y88b 888 888  888 Y88b.  
    //  888         "Y88888P"      888P     Y888  "Y8888  888  "Y88888 888  888  "Y888 
    //                                                             888                 
    //                                                        Y8b d88P                 
    //                                                         "Y88P"                  

    edm::Handle<std::vector<PileupSummaryInfo>> pileup;
    iEvent.getByToken(pileupInfoToken_, pileup);
    int nPU = 0;
    if (pileup.isValid())
      { // protection for data
	for (std::vector<PileupSummaryInfo>::const_iterator iPV = pileup->begin(); iPV != pileup->end(); ++iPV)
	  {
	    if (iPV->getBunchCrossing() == 0)
	      {
		nPU = iPV->getTrueNumInteractions();  
		//  numGenPV = iPV->getPU_NumInteractions();
		break;
	      }
	  }
      }

    double puweight   = hPUweight     ->GetBinContent(hPUweight     ->GetXaxis()->FindBin(nPU));
    double puweightUp = hPUweight_MBup->GetBinContent(hPUweight_MBup->GetXaxis()->FindBin(nPU));
    double puweightDn = hPUweight_MBdn->GetBinContent(hPUweight_MBdn->GetXaxis()->FindBin(nPU));

    if (verbose_) std::cout << "PU weight: " << puweight << std::endl;

    h_NtrueIntPU->Fill(nPU);
    h_NPV->Fill(nvtx);
    h_NPVreweighted->Fill(nvtx, puweight);

    // MARKER LHE WEIGHTS
    //  888      888    888 8888888888     888       888          d8b          888      888             
    //  888      888    888 888            888   o   888          Y8P          888      888             
    //  888      888    888 888            888  d8b  888                       888      888             
    //  888      8888888888 8888888        888 d888b 888  .d88b.  888  .d88b.  88888b.  888888 .d8888b  
    //  888      888    888 888            888d88888b888 d8P  Y8b 888 d88P"88b 888 "88b 888    88K      
    //  888      888    888 888            88888P Y88888 88888888 888 888  888 888  888 888    "Y8888b. 
    //  888      888    888 888            8888P   Y8888 Y8b.     888 Y88b 888 888  888 Y88b.       X88 
    //  88888888 888    888 8888888888     888P     Y888  "Y8888  888  "Y88888 888  888  "Y888  88888P' 
    //                                                                     888                          
    //                                                                Y8b d88P                          
    //                                                                 "Y88P"                           

    double Q2wgt_up = -999;
    double Q2wgt_down = -999;

    double NNPDF3wgt_up = -999;
    double NNPDF3wgt_down = -999;

    if (isZprime_ || isttbar_)
      {
	edm::Handle<LHEEventProduct> EvtHandle;
	iEvent.getByToken(theSrc_, EvtHandle);

	if (EvtHandle.isValid())
	  {
	    double centralWgt = EvtHandle->weights()[0].wgt;

	    //Q^2 uncertainties
	    if (verbose_) cout << "Calculating Q^2 uncertainties." << endl;
	    double maxQ2wgt_frac = 1;
	    double minQ2wgt_frac = 1;
      
	    if (verbose_) cout << "Q^2 loop" << endl;
	    for (unsigned int iLHE = 0; iLHE < 9; ++iLHE)
	      {
		if (iLHE != 5 && iLHE != 7)
		  {
		    double Q2wgt = EvtHandle->weights()[iLHE].wgt;
		    if (verbose_) cout << "Q^2 Weight: " << Q2wgt << endl;
		    double Q2wgt_frac = Q2wgt/(centralWgt);
		    if (verbose_) cout << "Fractional Q^2 Weight: " << Q2wgt_frac << endl;
		    maxQ2wgt_frac = max(maxQ2wgt_frac, Q2wgt_frac);
		    minQ2wgt_frac = min(minQ2wgt_frac, Q2wgt_frac);
		  }
	      }
      
	    Q2wgt_up = maxQ2wgt_frac;
	    Q2wgt_down = minQ2wgt_frac;

	    //NNPDF3 uncertainties
	    if (verbose_) cout << "Calculating NNPDF3 uncertainties." << endl;
	    double NNPDF3wgtAvg = 0.0;
	    double NNPDF3wgtRMS = 0.0;
	    double NNPDF3wgt = 0.0;
	    double NNPDF3wgt_frac = 0.0;

	    //ttbar
	    unsigned int PDFstart = 9;
	    unsigned int PDFend = 109;

	    //Zprime
	    if (isZprime_)
	      {
		PDFstart = 10;
		PDFend = 110;
	      }

	    //Making sure central PDF isn't zero                                                                                              
	    if (centralWgt == 0)
	      {
		NNPDF3wgt_up = 0.0;
		NNPDF3wgt_down = 0.0;
		if (verbose_) cout << "Unphysical: central PDF weight is zero!" << endl;
	      }
	    else
	      {
		for (unsigned int i_lhePDF = PDFstart; i_lhePDF < PDFend; ++i_lhePDF)
		  {
		    NNPDF3wgt = EvtHandle->weights()[i_lhePDF].wgt;
		    NNPDF3wgt_frac = NNPDF3wgt/(centralWgt);
		    NNPDF3wgtAvg += NNPDF3wgt_frac;
		    if (verbose_)
		      {
			cout << "-----" << endl;
			cout << i_lhePDF - PDFstart << endl;
			cout << "Fractional PDF weight: " << NNPDF3wgt_frac << endl;
			cout << "-----" << endl;
			cout << "" << endl;
		      }
		  }

		NNPDF3wgtAvg = NNPDF3wgtAvg/(PDFend - PDFstart);
		if (verbose_) cout << NNPDF3wgtAvg;
      
		for (unsigned int i_lhePDF = PDFstart; i_lhePDF < PDFend; ++i_lhePDF)
		  {
		    NNPDF3wgt = EvtHandle->weights()[i_lhePDF].wgt;
		    NNPDF3wgt_frac = NNPDF3wgt/(centralWgt);
		    NNPDF3wgtRMS += (NNPDF3wgt_frac - NNPDF3wgtAvg) * (NNPDF3wgt_frac - NNPDF3wgtAvg);
		  }
	
		NNPDF3wgtRMS = sqrt(NNPDF3wgtRMS/(PDFend - PDFstart - 1));
		NNPDF3wgt_up = 1.0 + NNPDF3wgtRMS;
		NNPDF3wgt_down = 1.0 - NNPDF3wgtRMS;
	      }
	  }
      }

    else if (isRSG_)
      {
	edm::Handle<GenEventInfoProduct> pdfstuff;
	iEvent.getByToken(pdfToken_, pdfstuff);

	LHAPDF::usePDFMember(1, 0);

	float q = pdfstuff->pdf()->scalePDF;

	int id1 = pdfstuff->pdf()->id.first;
	double x1 = pdfstuff->pdf()->x.first;
	//double pdf1 = pdfstuff->pdf()->xPDF.first;                                                                                           

	int id2 = pdfstuff->pdf()->id.second;
	double x2 = pdfstuff->pdf()->x.second;
	//double pdf2 = pdfstuff->pdf()->xPDF.second;                                                                                          

	double xpdf1 = LHAPDF::xfx(1, x1, q, id1);
	double xpdf2 = LHAPDF::xfx(1, x2, q, id2);
	double w0 = xpdf1 * xpdf2;
	double sumsq = 0.0;
	for (int i = 1; i <= 100; ++i)
	  {
	    LHAPDF::usePDFMember(1, i);
	    double xpdf1_new = LHAPDF::xfx(1, x1, q, id1);
	    double xpdf2_new = LHAPDF::xfx(1, x2, q, id2);
	    double weight = xpdf1_new * xpdf2_new / w0;
	    sumsq += (weight - w0) * (weight - w0);
	  }


	double rmsWt = sqrt((1./99.) * sumsq);

	if (rmsWt > 1.0)
	  {

	    NNPDF3wgt_up = rmsWt;
	    NNPDF3wgt_down = 2 - rmsWt;
	  }

	if (rmsWt < 1.0)
	  {

	    NNPDF3wgt_down = rmsWt;
	    NNPDF3wgt_up = 2 - rmsWt;
	  }
    
      }

    //MY ADDITION
    /*edm::Handle<double> weight1st;
    iEvent.getByToken(mgreweight1_, weight1st);
    anom_weight1 = *weight1st;

    edm::Handle<double> weight2nd;
    iEvent.getByToken(mgreweight2_, weight2nd);
    anom_weight2 = *weight2nd;

    edm::Handle<double> weight3rd;
    iEvent.getByToken(mgreweight3_, weight3rd);
    anom_weight3 = *weight3rd;

    edm::Handle<double> weight4th;
    iEvent.getByToken(mgreweight4_, weight4th);
    anom_weight4 = *weight4th;

    edm::Handle<double> weight5th;
    iEvent.getByToken(mgreweight5_, weight5th);
    anom_weight5 = *weight5th;

    edm::Handle<double> weight6th;
    iEvent.getByToken(mgreweight6_, weight6th);
    anom_weight6 = *weight6th;

    edm::Handle<double> weight7th;
    iEvent.getByToken(mgreweight7_, weight7th);
    anom_weight7 = *weight7th;

    edm::Handle<double> weight8th;
    iEvent.getByToken(mgreweight8_, weight8th);
    anom_weight8 = *weight8th;

    edm::Handle<double> weight9th;
    iEvent.getByToken(mgreweight9_, weight9th);
    anom_weight9 = *weight9th;

    edm::Handle<double> weight10th;
    iEvent.getByToken(mgreweight10_, weight10th);
    anom_weight10 = *weight10th;

    edm::Handle<double> weight11th;
    iEvent.getByToken(mgreweight11_, weight11th);
    anom_weight11 = *weight11th;

    edm::Handle<double> weight12th;
    iEvent.getByToken(mgreweight12_, weight12th);
    anom_weight12 = *weight12th;

    edm::Handle<double> weight13th;
    iEvent.getByToken(mgreweight13_, weight13th);
    anom_weight13 = *weight13th;

    edm::Handle<double> weight14th;
    iEvent.getByToken(mgreweight14_, weight14th);
    anom_weight14 = *weight14th;

    edm::Handle<double> weight15th;
    iEvent.getByToken(mgreweight15_, weight15th);
    anom_weight15 = *weight15th;

    edm::Handle<double> weight16th;
    iEvent.getByToken(mgreweight16_, weight16th);
    anom_weight16 = *weight16th;

    edm::Handle<double> weight17th;
    iEvent.getByToken(mgreweight17_, weight17th);
    anom_weight17 = *weight17th;

    edm::Handle<double> weight18th;
    iEvent.getByToken(mgreweight18_, weight18th);
    anom_weight18 = *weight18th;

    edm::Handle<double> weight19th;
    iEvent.getByToken(mgreweight19_, weight19th);
    anom_weight19 = *weight19th;

    edm::Handle<double> weight20th;
    iEvent.getByToken(mgreweight20_, weight20th);
    anom_weight20 = *weight20th;

    edm::Handle<double> weight21st;
    iEvent.getByToken(mgreweight21_, weight21st);
    anom_weight21 = *weight21st;

    edm::Handle<double> weight22nd;
    iEvent.getByToken(mgreweight22_, weight22nd);
    anom_weight22 = *weight22nd;

    edm::Handle<double> weight23rd;
    iEvent.getByToken(mgreweight23_, weight23rd);
    anom_weight23 = *weight23rd;

    edm::Handle<double> weight24th;
    iEvent.getByToken(mgreweight24_, weight24th);
    anom_weight24 = *weight24th;

    edm::Handle<double> weight25th;
    iEvent.getByToken(mgreweight25_, weight25th);
    anom_weight25 = *weight25th;

    edm::Handle<double> weight26th;
    iEvent.getByToken(mgreweight26_, weight26th);
    anom_weight26 = *weight26th;

    edm::Handle<double> weight27th;
    iEvent.getByToken(mgreweight27_, weight27th);
    anom_weight27 = *weight27th;

    edm::Handle<double> weight28th;
    iEvent.getByToken(mgreweight28_, weight28th);
    anom_weight28 = *weight28th;

    edm::Handle<double> weight29th;
    iEvent.getByToken(mgreweight29_, weight29th);
    anom_weight29 = *weight29th;

    edm::Handle<double> weight30th;
    iEvent.getByToken(mgreweight30_, weight30th);
    anom_weight30 = *weight30th;

    edm::Handle<double> weight31st;
    iEvent.getByToken(mgreweight31_, weight31st);
    anom_weight31 = *weight31st;

    edm::Handle<double> weight32nd;
    iEvent.getByToken(mgreweight32_, weight32nd);
    anom_weight32 = *weight32nd;

    edm::Handle<double> weight33rd;
    iEvent.getByToken(mgreweight33_, weight33rd);
    anom_weight33 = *weight33rd;

    edm::Handle<double> weight34th;
    iEvent.getByToken(mgreweight34_, weight34th);
    anom_weight34 = *weight34th;

    edm::Handle<double> weight35th;
    iEvent.getByToken(mgreweight35_, weight35th);
    anom_weight35 = *weight35th;

    edm::Handle<double> weight36th;
    iEvent.getByToken(mgreweight36_, weight36th);
    anom_weight36 = *weight36th;

    edm::Handle<double> weight37th;
    iEvent.getByToken(mgreweight37_, weight37th);
    anom_weight37 = *weight37th;

    edm::Handle<double> weight38th;
    iEvent.getByToken(mgreweight38_, weight38th);
    anom_weight38 = *weight38th;

    edm::Handle<double> weight39th;
    iEvent.getByToken(mgreweight39_, weight39th);
    anom_weight39 = *weight39th;

    edm::Handle<double> weight40th;
    iEvent.getByToken(mgreweight40_, weight40th);
    anom_weight40 = *weight40th;

    edm::Handle<double> weight41st;
    iEvent.getByToken(mgreweight41_, weight41st);
    anom_weight41 = *weight41st;

    edm::Handle<double> weight42nd;
    iEvent.getByToken(mgreweight42_, weight42nd);
    anom_weight42 = *weight42nd;

    edm::Handle<double> weight43rd;
    iEvent.getByToken(mgreweight43_, weight43rd);
    anom_weight43 = *weight43rd;

    edm::Handle<double> weight44th;
    iEvent.getByToken(mgreweight44_, weight44th);
    anom_weight44 = *weight44th;

    edm::Handle<double> weight45th;
    iEvent.getByToken(mgreweight45_, weight45th);
    anom_weight45 = *weight45th;

    edm::Handle<double> weight46th;
    iEvent.getByToken(mgreweight46_, weight46th);
    anom_weight46 = *weight46th;

    edm::Handle<double> weight47th;
    iEvent.getByToken(mgreweight47_, weight47th);
    anom_weight47 = *weight47th;

    edm::Handle<double> weight48th;
    iEvent.getByToken(mgreweight48_, weight48th);
    anom_weight48 = *weight48th;

    edm::Handle<double> weight49th;
    iEvent.getByToken(mgreweight49_, weight49th);
    anom_weight49 = *weight49th;

    edm::Handle<double> weight50th;
    iEvent.getByToken(mgreweight50_, weight50th);
    anom_weight50 = *weight50th;*/


    // MARKER RHO
    // 8888888b.  888               
    // 888   Y88b 888               
    // 888    888 888               
    // 888   d88P 88888b.   .d88b.  
    // 8888888P"  888 "88b d88""88b 
    // 888 T88b   888  888 888  888 
    // 888  T88b  888  888 Y88..88P 
    // 888   T88b 888  888  "Y88P"  
    //                

    Handle<double> rhoH;
    iEvent.getByToken(rhoToken_, rhoH);
    double rho = *rhoH;

    // MARKER MUON
    // 888b     d888                            
    // 8888b   d8888                            
    // 88888b.d88888                            
    // 888Y88888P888 888  888  .d88b.  88888b.  
    // 888 Y888P 888 888  888 d88""88b 888 "88b 
    // 888  Y8P  888 888  888 888  888 888  888 
    // 888   "   888 Y88b 888 Y88..88P 888  888 
    // 888       888  "Y88888  "Y88P"  888  888 
    //                                          
                                             
    edm::Handle<pat::MuonCollection> muons;
    iEvent.getByToken(muonToken_, muons);

    TLorentzVector mu0_p4;
    bool mu0_isTight = false;
    bool mu0_isMedium = false;
    double mu0_iso04 = 0;
    int count_mu = 0;
    for (const pat::Muon &mu : *muons)
      {
	if (mu.pt() < 200 || !mu.isLooseMuon() || fabs(mu.eta()) > 2.1) continue; //changed min pt from 30 to 200


	if (count_mu == 0)
	  {
	    mu0_p4.SetPtEtaPhiM(mu.pt(), mu.eta(), mu.phi(), mu.mass());
	    if (mu.isTightMuon(PV)) mu0_isTight = true;
	    // ICHEP HIP MEDIUM MUON ID https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideMuonIdRun2#Short_Term_Instructions_for_ICHE
	    bool goodGlob = mu.isGlobalMuon() && 
	      mu.globalTrack()->normalizedChi2() < 3 && 
						   mu.combinedQuality().chi2LocalPosition < 12 && 
											    mu.combinedQuality().trkKink < 20; 
	    bool isMedium = muon::isLooseMuon(mu) && 
															   mu.innerTrack()->validFraction() > 0.49 && 
						   muon::segmentCompatibility(mu) > (goodGlob ? 0.303 : 0.451); 
	    if (isMedium) mu0_isMedium = true;


	    // https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId#Accessing_PF_Isolation_from_reco
	    double sumChargedHadronPt = mu.pfIsolationR04().sumChargedHadronPt;
	    double sumNeutralHadronPt = mu.pfIsolationR04().sumNeutralHadronEt;
	    double sumPhotonPt        = mu.pfIsolationR04().sumPhotonEt;
	    double sumPUPt            = mu.pfIsolationR04().sumPUPt;
	    double pt                 = mu.pt();
	    double iso04 = (sumChargedHadronPt+TMath::Max(0., sumNeutralHadronPt+sumPhotonPt - 0.5 * sumPUPt))/pt;
	    mu0_iso04 = iso04;

	    if (verbose_) cout << "Muon pT " << mu.pt() << " iso04 " << iso04 << endl;
	  } 
	// printf("muon with pt %4.1f, dz(PV) %+5.3f, POG loose id %d, tight id %d\n",
	// mu.pt(), mu.muonBestTrack()->dz(PV.position()), mu.isLooseMuon(), mu.isTightMuon(PV));
	count_mu++;
      }


    // MARKER ELECTRON
    // 8888888888 888                   888                              
    // 888        888                   888                              
    // 888        888                   888                              
    // 8888888    888  .d88b.   .d8888b 888888 888d888  .d88b.  88888b.  
    // 888        888 d8P  Y8b d88P"    888    888P"   d88""88b 888 "88b 
    // 888        888 88888888 888      888    888     888  888 888  888 
    // 888        888 Y8b.     Y88b.    Y88b.  888     Y88..88P 888  888 
    // 8888888888 888  "Y8888   "Y8888P  "Y888 888      "Y88P"  888  888 
    //                                                                   
    edm::Handle<pat::ElectronCollection> electrons;
    iEvent.getByToken(electronToken_, electrons);

    TLorentzVector el0_p4;
    Float_t el0_absiso           = 0.0;
    Float_t el0_relIsoWithDBeta  = 0.0;
    Float_t el0_absiso_EA        = 0.0;
    Float_t el0_relIsoWithEA     = 0.0;
    int count_el = 0;
    for (const pat::Electron &el : *electrons)
      {
	if (el.pt() < 200.0 || fabs(el.eta()) > 2.4) continue; //changing min pt from 50 to 200

	float eta = el.eta();
	// Implemation of loose Quality cuts (need to study this and improve)
	// https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2
	// Recomended isolation variables are stored separately and not included in the loose quality cut


	bool passLoose = false;
	float ooEmooP_; 
	if (el.ecalEnergy() == 0.0)
	  {
	    if (verbose_) printf("Electron energy is zero!\n");
	    ooEmooP_ = 999.0;
	  }
	else if (!std::isfinite(el.ecalEnergy()))
	  {
	    if (verbose_) printf("Electron energy is not finite!\n");
	    ooEmooP_ = 998.0;
	  }
	else ooEmooP_ = fabs(1.0/el.ecalEnergy() - el.eSuperClusterOverP()/el.ecalEnergy());
	float missHits = el.gsfTrack()->hitPattern().numberOfLostTrackerHits(HitPattern::MISSING_INNER_HITS);
      

	GsfElectron::PflowIsolationVariables pfIso = el.pfIsolationVariables();
	float absiso = pfIso.sumChargedHadronPt + max(0.0, pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - 0.5 * pfIso.sumPUPt);
	float relIsoWithDBeta = absiso/el.pt();

	float effArea = 0.0;
	if (abs(eta) > 0.0 && abs(eta) <= 1.0) effArea = 0.1752;
	if (abs(eta) > 1.0 && abs(eta) <= 1.479) effArea = 0.1862;
	if (abs(eta) > 1.479 && abs(eta) <= 2.0) effArea = 0.1411;
	if (abs(eta) > 2.0 && abs(eta) <= 2.2) effArea = 0.1534;
	if (abs(eta) > 2.2 && abs(eta) <= 2.3) effArea = 0.1903;
	if (abs(eta) > 2.3 && abs(eta) <= 2.4) effArea = 0.2243;
	if (abs(eta) > 2.4 && abs(eta) <= 2.5) effArea = 0.2687;

	float absiso_EA = pfIso.sumChargedHadronPt + max(0.0, pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - rho * effArea);
	float relIsoWithEA = absiso_EA/el.pt();



	if (fabs(el.eta()) < 1.479)
	  {
	    if (el.full5x5_sigmaIetaIeta()                <      0.011      &&
		fabs(el.deltaEtaSuperClusterTrackAtVtx()) <      0.00477    &&
		fabs(el.deltaPhiSuperClusterTrackAtVtx()) <      0.222      &&
		el.hadronicOverEm()                       <      0.298      &&
		ooEmooP_                                  <      0.241      &&
		missHits                                  <=     1) passLoose = true;
	  }
	else
	  { 
	    if (el.full5x5_sigmaIetaIeta()                <      0.0314     &&
		fabs(el.deltaEtaSuperClusterTrackAtVtx()) <      0.00868    &&
		fabs(el.deltaPhiSuperClusterTrackAtVtx()) <      0.213      &&
		el.hadronicOverEm()                       <      0.101      &&
		ooEmooP_                                  <      0.14       &&
		missHits                                  <=     1) passLoose = true;
	  }
	if (!passLoose) continue;

	if (count_el == 0)
	  {
	    el0_p4.SetPtEtaPhiM(el.pt(), el.eta(), el.phi(), el.mass());

	    el0_absiso           = absiso;
	    el0_relIsoWithDBeta  = relIsoWithDBeta;
	    el0_absiso_EA        = absiso_EA;
	    el0_relIsoWithEA     = relIsoWithEA;

	    if (verbose_) cout << "Electron pT " << el.pt() << endl;
	  } 
	count_el++;
	//printf("elec with pt %4.1f, supercluster eta %+5.3f, sigmaIetaIeta %.3f (%.3f with full5x5 shower shapes), lost hits %d, pass conv veto %d\n",
	//              el.pt(), el.superCluster()->eta(), el.sigmaIetaIeta(), el.full5x5_sigmaIetaIeta(), el.gsfTrack()->trackerExpectedHitsInner().numberOfLostHits(), el.passConversionVeto());
      }

    TLorentzVector lep0_p4;
    if (count_mu == 1 &&  count_el == 0)      lep0_p4 = mu0_p4;
    else if (count_el == 1 && count_mu == 0) lep0_p4 = el0_p4;
    int count_lep = count_mu + count_el;

    if (verbose_)
      {
	cout << "count_mu  " << count_mu << endl;
	cout << "count_el  " << count_el << endl;
	cout << "count_lep " << count_lep << endl;
      }

 
    // MARKER MET
    // 888b     d888 8888888888 88888888888 
    // 8888b   d8888 888            888     
    // 88888b.d88888 888            888     
    // 888Y88888P888 8888888        888     
    // 888 Y888P 888 888            888     
    // 888  Y8P  888 888            888     
    // 888   "   888 888            888     
    // 888       888 8888888888     888     
    //                                      

    if (verbose_) cout << "debug: about to grab met" << endl;
    edm::Handle<pat::METCollection> mets;
    iEvent.getByToken(metToken_, mets);
    const pat::MET &met = mets->front();
    if (verbose_)
      {
	cout << "MET pt " << met.pt() << endl;
	cout << "MET phi " << met.phi() << endl;
	cout << "MET sumEt " << met.sumEt() << endl;
	if (!iEvent.isRealData())  cout << "genMET " << met.genMET()->pt() << endl;
      }
    // MARKER AK4 CHS JETS
    //        d8888 888    d8P      d8888       .d8888b.  888    888  .d8888b.         d8b          888             
    //       d88888 888   d8P      d8P888      d88P  Y88b 888    888 d88P  Y88b        Y8P          888             
    //      d88P888 888  d8P      d8P 888      888    888 888    888 Y88b.                          888             
    //     d88P 888 888d88K      d8P  888      888        8888888888  "Y888b.         8888  .d88b.  888888 .d8888b  
    //    d88P  888 8888888b    d88   888      888        888    888     "Y88b.       "888 d8P  Y8b 888    88K      
    //   d88P   888 888  Y88b   8888888888     888    888 888    888       "888        888 88888888 888    "Y8888b. 
    //  d8888888888 888   Y88b        888      Y88b  d88P 888    888 Y88b  d88P        888 Y8b.     Y88b.       X88 
    // d88P     888 888    Y88b       888       "Y8888P"  888    888  "Y8888P"         888  "Y8888   "Y888  88888P' 
    //                                                                                 888                          
    //                                                                                d88P                          

    //int count_AK4MINI = 0;
    TLorentzVector AK4_dRMinLep_p4;
    TLorentzVector AK4_btagged_dRMinLep_p4;
    double AK4_dRMinLep_bdisc = -99;
    double AK4_btagged_dRMinLep_bdisc = -99;
    double AK4_dRMinLep  = 99;
    double AK4_btagged_dRMinLep = 99;

    bool ak4_btag_loose  = false;
    bool ak4_btag_medium = false;
    bool ak4_btag_tight  = false;

    double HT_AK4_pt30           = 0;
    double HT_AK4_pt30_corrUp    = 0;
    double HT_AK4_pt30_corrDn    = 0;
    double HT_AK4_pt30_smearNom  = 0;
    double HT_AK4_pt30_smearUp   = 0;
    double HT_AK4_pt30_smearDn   = 0;

    edm::Handle<pat::JetCollection> AK4MINI;
    iEvent.getByToken(ak4jetToken_, AK4MINI);

    edm::Handle<reco::GenJetCollection> AK4GENJET;  
    iEvent.getByToken(ak4genjetToken_, AK4GENJET);

    if (verbose_) cout << "debug: about to grab ak4 jets" << endl;

    for (const pat::Jet &ijet : *AK4MINI)
      { //NOTE: ijet is an iterator over AK4MINI leading to AK4_dRMinLep_p4 
    
	if (ijet.pt() < 30 || fabs(ijet.eta()) > 10) {std::cout << "failed pt less than 30 or eta in range of 2.4\n"; continue;} //changing rapidity max from 2.4 to 10

	//------------------------------------
	// Noise jet ID
	//------------------------------------    
	double NHF       = ijet.neutralHadronEnergyFraction();
	double NEMF      = ijet.neutralEmEnergyFraction();
	double CHF       = ijet.chargedHadronEnergyFraction();
	// double MUF       = ijet.muonEnergyFraction();
	double CEMF      = ijet.chargedEmEnergyFraction();
	double NumConst  = ijet.chargedMultiplicity()+ijet.neutralMultiplicity();
	double NM        = ijet.neutralMultiplicity();
	double CM        = ijet.chargedMultiplicity(); 
	double eta       = ijet.eta(); 

	bool goodJet_looseJetID =  
	  (fabs(eta) <= 2.4 && NHF < 0.99 && NEMF < 0.99 && NumConst > 1 && CHF > 0.0  && CM > 0 && CEMF < 0.99) 
	  || (fabs(eta) <= 2.7 && fabs(eta) > 2.4 && NHF < 0.99 && NEMF < 0.99 && NumConst > 1) 
	  || (fabs(eta) <= 3.0 && fabs(eta) > 2.7 && NEMF < 0.9 && NM > 2) 
	  || (fabs(eta)  > 3.0 && NEMF < 0.9 && NM > 10);
	if (verbose_) cout << "  goodJet = " << goodJet_looseJetID << endl;

	if (!goodJet_looseJetID)
	  {
	    if (verbose_) cout << "  bad AK4 jet. skip.  (pt " << ijet.pt() << " eta " << ijet.eta() << " NumConst " << NumConst << ")" << endl;
	    continue;
	  }

	//------------------------------------
	// AK4CHS JEC correction 
	//------------------------------------
	reco::Candidate::LorentzVector uncorrJet = ijet.correctedP4(0);
	JetCorrectorAK4chs->setJetEta(uncorrJet.eta());
	JetCorrectorAK4chs->setJetPt (uncorrJet.pt());
	JetCorrectorAK4chs->setJetE  (uncorrJet.energy());
	JetCorrectorAK4chs->setJetA  (ijet.jetArea());
	JetCorrectorAK4chs->setRho   (rho);
	JetCorrectorAK4chs->setNPV   (nvtx);
	double corr = JetCorrectorAK4chs->getCorrection();

	reco::Candidate::LorentzVector corrJet = corr * uncorrJet;
	if (verbose_) cout << "uncorrected AK4 jet pt " << uncorrJet.pt() << " corrected jet pt " << corrJet.pt() << endl;
    
	//------------------------------------
	// AK4CHS JEC uncertainty
	//------------------------------------
	double corrDn = 1.0;
	JetCorrUncertAK4chs->setJetPhi(corrJet.phi());
	JetCorrUncertAK4chs->setJetEta(corrJet.eta());
	JetCorrUncertAK4chs->setJetPt(corrJet.pt());
	corrDn = corr - JetCorrUncertAK4chs->getUncertainty(0);
	double corrUp = 1.0;
	JetCorrUncertAK4chs->setJetPhi(corrJet.phi());
	JetCorrUncertAK4chs->setJetEta(corrJet.eta());
	JetCorrUncertAK4chs->setJetPt(corrJet.pt());
	corrUp = corr + JetCorrUncertAK4chs->getUncertainty(1);

	if (verbose_) cout << "  corrDn " << corrDn << " corrUp " << corrUp << endl;

	//------------------------------------
	// GenJet  matched
	//------------------------------------ 
	double genpt = 0;
	if (!iEvent.isRealData())
	  {
	    const reco::GenJet* genJet = ijet.genJet();
	    if (genJet)
	      {
		genpt = genJet->pt();
		if (verbose_) cout << "  ak4 genJet pt " << genJet->pt() << " mass " << genJet->mass() << endl;
	      }
	  }
	//------------------------------------
	// JER SF
	//------------------------------------
	double ptsmear   = 1;
	double ptsmearUp = 1;
	double ptsmearDn = 1;
	if (!iEvent.isRealData())
	  {
	    double jer_sf    = jer_scaler.getScaleFactor({{JME::Binning::JetEta, corrJet.eta()}});
	    double jer_sf_up = jer_scaler.getScaleFactor({{JME::Binning::JetEta, corrJet.eta()}}, Variation::UP);
	    double jer_sf_dn = jer_scaler.getScaleFactor({{JME::Binning::JetEta, corrJet.eta()}}, Variation::DOWN);
	    if (verbose_) std::cout << "  JER Scale factors (Nominal / Up / Down) : " << jer_sf << " / " << jer_sf_up << " / " << jer_sf_dn << std::endl;
	    double recopt    = corrJet.pt();
	    // double genpt     = GenJetMatched.Perp();
	    double deltapt   = (recopt-genpt) * (jer_sf - 1.0);
	    double deltaptUp = (recopt-genpt) * (jer_sf_up - 1.0);
	    double deltaptDn = (recopt-genpt) * (jer_sf_dn - 1.0);
	    ptsmear   = std::max((double)0.0, (recopt+deltapt)/recopt);
	    ptsmearUp = std::max((double)0.0, (recopt+deltaptUp)/recopt);
	    ptsmearDn = std::max((double)0.0, (recopt+deltaptDn)/recopt);
	  }

	if (verbose_) cout << "  ptsmear " << ptsmear << " ptsmearUp " << ptsmearDn << " ptsmearDn " << ptsmearUp << endl;


	//------------------------------------
	// AK4 variables 
	//------------------------------------
	double pt           = corrJet.pt(); std::cout << "pt for jet4 is " << pt << endl;
	double mass         = corrJet.mass(); std::cout << "mass for jet4 is " << mass << endl;
	eta                 = corrJet.eta(); std::cout << "eta for jet4 is " << eta << endl;
	double phi          = corrJet.phi(); std::cout << "phi for jet4 is " << phi << endl;
	//double rapidity     = ijet.rapidity();
	//double ndau         = ijet.numberOfDaughters();
	double bdisc        = ijet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"); 
	//double puid         = ijet.userFloat("pileupJetId:fullDiscriminant");
 
	//------------------------------------
	// HT calculation
	//------------------------------------

	if (corrJet.pt() > 30.0)             HT_AK4_pt30           +=   pt;
	if (corrUp * corrJet.pt() > 30.0)    HT_AK4_pt30_corrUp    +=   corrUp * corrJet.pt();
	if (corrDn * corrJet.pt() > 30.0)    HT_AK4_pt30_corrDn    +=   corrDn * corrJet.pt();
	if (ptsmear * corrJet.pt() > 30.0)   HT_AK4_pt30_smearNom  +=   ptsmear * corrJet.pt();
	if (ptsmearUp * corrJet.pt() > 30.0) HT_AK4_pt30_smearUp   +=   ptsmearUp * corrJet.pt();
	if (ptsmearDn * corrJet.pt() > 30.0) HT_AK4_pt30_smearDn   +=   ptsmearDn * corrJet.pt();

	//------------------------------------
	// Find AK4 jet closest to lepton
	//------------------------------------ 
	double deltaRlep = deltaR(corrJet.eta(), corrJet.phi(), lep0_p4.Eta(), lep0_p4.Phi());

	if (pt > 40.0 && fabs(eta) < 2.4 && goodJet_looseJetID)
	  {
	    if (deltaRlep<AK4_dRMinLep)
	      {
		AK4_dRMinLep = deltaRlep;
		AK4_dRMinLep_p4.SetPtEtaPhiM(pt, eta, phi, mass);
		AK4_dRMinLep_bdisc = bdisc;
	      }
	  }



	//------------------------------------
	// Find Loose b-tagged AK4 jet closest to lepton
	//------------------------------------ 
	if (pt > 40.0 && fabs(eta) < 2.4 && goodJet_looseJetID && bdisc > 0.460)
	  {
	    if (deltaRlep<AK4_btagged_dRMinLep)
	      {
		AK4_btagged_dRMinLep = deltaRlep;
		AK4_btagged_dRMinLep_p4.SetPtEtaPhiM(pt, eta, phi, mass); 
		AK4_btagged_dRMinLep_bdisc = bdisc;
	      }
	  }

	//------------------------------------
	// Check if there is a b-tagged AK4 jet in the lepton hemisphere
	//------------------------------------ 
	double deltaPhiLep = fabs(deltaPhi(phi, lep0_p4.Phi()));  
	if (pt > 40.0 && fabs(eta) < 2.4 && goodJet_looseJetID)
	  {              
	    if (deltaPhiLep <  3.14/2.0)
	      {
		if (bdisc > 0.460) ak4_btag_loose  = true;
		if (bdisc > 0.800) ak4_btag_medium = true;
		if (bdisc > 0.935) ak4_btag_tight  = true;
	      }
	  }
      } //end AK4 loop

    if (verbose_)
      {
	cout << "AK4 summary:" << endl;
	cout << "  closest ak4 jet to lepton:" << endl;
	cout << "    pt =  " << AK4_dRMinLep_p4.Perp() << endl;
	cout << "    bdisc =  " << AK4_dRMinLep_bdisc << endl;
	cout << "    dR  = " << AK4_dRMinLep << endl;
	cout << "  closest loose b-tagged ak4 jet to lepton:" << endl;
	cout << "    pt =  " << AK4_btagged_dRMinLep_p4.Perp() << endl;
	cout << "    bdisc =  " << AK4_btagged_dRMinLep_bdisc << endl;
	cout << "    dR  = " << AK4_btagged_dRMinLep << endl;
	cout << "  b-tagged jet in hemisphere around lepton?" << endl;
	cout << "    ak4_btag_loose  " << ak4_btag_loose  << endl; 
	cout << "    ak4_btag_medium " << ak4_btag_medium << endl; 
	cout << "    ak4_btag_tight  " << ak4_btag_tight  << endl; 
      }


    //      MARKER AK8 CHS JETS
    //         d8888 888    d8P   .d8888b.       .d8888b.  888    888  .d8888b.         d8b          888             
    //        d88888 888   d8P   d88P  Y88b     d88P  Y88b 888    888 d88P  Y88b        Y8P          888             
    //       d88P888 888  d8P    Y88b. d88P     888    888 888    888 Y88b.                          888             
    //      d88P 888 888d88K      "Y88888"      888        8888888888  "Y888b.         8888  .d88b.  888888 .d8888b  
    //     d88P  888 8888888b    .d8P""Y8b.     888        888    888     "Y88b.       "888 d8P  Y8b 888    88K      
    //    d88P   888 888  Y88b   888    888     888    888 888    888       "888        888 88888888 888    "Y8888b. 
    //   d8888888888 888   Y88b  Y88b  d88P     Y88b  d88P 888    888 Y88b  d88P        888 Y8b.     Y88b.       X88 
    //  d88P     888 888    Y88b  "Y8888P"       "Y8888P"  888    888  "Y8888P"         888  "Y8888   "Y888  88888P' 
    //                                                                                  888                          
    //                                                                               d88P                          
    //                                                                             888P"                                                                                              

    edm::Handle<pat::JetCollection> AK8CHS;
    iEvent.getByToken(ak8jetToken_, AK8CHS);

    edm::Handle<reco::GenJetCollection> AK8GENJET;  
    iEvent.getByToken(ak8genjetToken_, AK8GENJET);

    edm::Handle<pat::JetCollection> AK8CHSsub;
    edm::Handle<pat::JetCollection> AK8PUPPI;
    edm::Handle<pat::JetCollection> AK8PUPPIsub;
    if (useToolbox_)
      {
	iEvent.getByToken(ak8CHSSoftDropSubjetsToken_, AK8CHSsub);
	iEvent.getByToken(puppijetToken_, AK8PUPPI);
	iEvent.getByToken(ak8PuppiSoftDropSubjetsToken_, AK8PUPPIsub);
      }

    int count_AK8CHS = 0;
    int count_fill_leptTree = 0;

    TLorentzVector AK8jet_SemiLept_P4corr; //delta8s
    TLorentzVector AK8jet0_P4corr; //delta80
    TLorentzVector AK8jet1_P4corr; //delta81
    TLorentzVector PUPPIjet0_P4; //unused
    TLorentzVector PUPPIjet1_P4; //unused
    TLorentzVector PUPPIjet0_P4corr; //deltapup0
    TLorentzVector PUPPIjet1_P4corr; //deltapup1

    TLorentzVector GenJetMatched0;
    TLorentzVector GenJetMatched1;

    if (verbose_) cout << "debug: about to grab ak8 jets" << endl;
#if 0
    for (const pat::Jet &ijet : *AK8CHS)
      {  
	std::cout << "we've opened the loop.\n";
	if (count_AK8CHS > 1) break;
	std::cout << "we've gotten past one break statement.";
	if (count_AK8CHS == 0 && ijet.pt() < 250) break;
	std::cout << "we've gotten past two break statements.";
	if (verbose_) cout << "\nJet " << count_AK8CHS << " with pT " << ijet.pt() << " sdMass " << ijet.userFloat("ak8PFJetsCHSSoftDropMass") << endl;

	//------------------------------------
	// Noise jet ID
	//------------------------------------    
	double NHF       = ijet.neutralHadronEnergyFraction();
	double NEMF      = ijet.neutralEmEnergyFraction();
	double CHF       = ijet.chargedHadronEnergyFraction();
	// double MUF       = ijet.muonEnergyFraction();
	double CEMF      = ijet.chargedEmEnergyFraction();
	double NumConst  = ijet.chargedMultiplicity()+ijet.neutralMultiplicity();
	double NM        = ijet.neutralMultiplicity();
	double CM        = ijet.chargedMultiplicity(); 
	double eta       = ijet.eta();

	bool goodJet_looseJetID =  
	  (fabs(eta) <= 2.4 && NHF < 0.99 && NEMF < 0.99 && NumConst > 1 && CHF > 0.0  && CM > 0 && CEMF < 0.99) 
	  || (fabs(eta) <= 2.7 && fabs(eta) > 2.4 && NHF < 0.99 && NEMF < 0.99 && NumConst > 1) 
	  || (fabs(eta) <= 3.0 && fabs(eta) > 2.7 && NEMF < 0.9 && NM > 2) 
	  || (fabs(eta)  > 3.0 && NEMF < 0.9 && NM > 10);
	if (verbose_) cout << "  goodJet = " << goodJet_looseJetID << endl;

	if (!goodJet_looseJetID)
	  {
	    if (verbose_) cout << "  bad AK8 jet. skip.  (pt " << ijet.pt() << " eta " << ijet.eta() << " NumConst " << NumConst << ")" << endl;
	    continue;
	  }

	//------------------------------------
	// AK8CHS JEC correction 
	//------------------------------------
	reco::Candidate::LorentzVector uncorrJet = ijet.correctedP4(0);
	JetCorrectorAK8chs->setJetEta(uncorrJet.eta());
	JetCorrectorAK8chs->setJetPt (uncorrJet.pt());
	JetCorrectorAK8chs->setJetE  (uncorrJet.energy());
	JetCorrectorAK8chs->setJetA  (ijet.jetArea());
	JetCorrectorAK8chs->setRho   (rho);
	JetCorrectorAK8chs->setNPV   (nvtx);
	double corr = JetCorrectorAK8chs->getCorrection();

	reco::Candidate::LorentzVector corrJet = corr * uncorrJet;
	if (verbose_) cout << " uncorrected AK8 jet pt " << uncorrJet.pt() << " corrected jet pt " << corrJet.pt() << endl;
    
	if (count_AK8CHS == 0) AK8jet0_P4corr.SetPtEtaPhiM(corrJet.pt(), corrJet.eta(), corrJet.phi(), corrJet.mass());
	if (count_AK8CHS == 1) AK8jet1_P4corr.SetPtEtaPhiM(corrJet.pt(), corrJet.eta(), corrJet.phi(), corrJet.mass());

	//------------------------------------
	// AK8CHS JEC L23 correction
	//------------------------------------
	JetCorrectorAK8chs->setJetEta(uncorrJet.eta());
	JetCorrectorAK8chs->setJetPt (uncorrJet.pt());
	JetCorrectorAK8chs->setJetE  (uncorrJet.energy());
	JetCorrectorAK8chs->setJetA  (ijet.jetArea());
	JetCorrectorAK8chs->setRho   (rho);
	JetCorrectorAK8chs->setNPV   (nvtx);
	// getSubCorrections member function returns the vector of the subcorrections UP to the given level. For example in the example above, factors[0] is the L1 correction and factors[3] is the L1+L2+L3+Residual correction. 
	vector<float> factors = JetCorrectorAK8chs->getSubCorrections();
	float corr_factor_L1      = 1.0;
	float corr_factor_L12     = 1.0;
	float corr_factor_L123    = 1.0;
	float corr_factor_L123res = 1.0;
	if (factors.size() > 0) corr_factor_L1       = factors[0];
	if (factors.size() > 1) corr_factor_L12      = factors[1];
	if (factors.size() > 2) corr_factor_L123     = factors[2];
	if (factors.size() > 3) corr_factor_L123res  = factors[3];
	double corr_factor_L2 = corr_factor_L12/corr_factor_L1;
	double corr_factor_L3 = corr_factor_L123/corr_factor_L12;
	double corr_factor_res = corr_factor_L123res/corr_factor_L123;
	//double corr_factor_L23 = corr_factor_L2*corr_factor_L3;
	double corr_factor_L23res = corr_factor_L2*corr_factor_L3*corr_factor_res;

	//------------------------------------
	// AK8CHS JEC uncertainty
	//------------------------------------
	double corrDn_L23  = 1.0;
	double corrDn_L123 = 1.0;
	JetCorrUncertAK8chs->setJetPhi(corrJet.phi());
	JetCorrUncertAK8chs->setJetEta(corrJet.eta());
	JetCorrUncertAK8chs->setJetPt(corrJet.pt());
	double corrDn_temp1 = JetCorrUncertAK8chs->getUncertainty(0);
	corrDn_L23   = corr_factor_L23res - corrDn_temp1;
	corrDn_L123 = corr - corrDn_temp1;
	double corrUp_L23  = 1.0;
	double corrUp_L123 = 1.0;
	JetCorrUncertAK8chs->setJetPhi(corrJet.phi());
	JetCorrUncertAK8chs->setJetEta(corrJet.eta());
	JetCorrUncertAK8chs->setJetPt(corrJet.pt());
	double corrUp_temp1 = JetCorrUncertAK8chs->getUncertainty(1);
	corrUp_L23   = corr_factor_L23res + corrUp_temp1;
	corrUp_L123 = corr + corrUp_temp1;

	//------------------------------------
	// GenJet  matched
	//------------------------------------ 
	TLorentzVector GenJetMatched;
	if (!iEvent.isRealData())
	  {
	    const reco::GenJet* genJet = ijet.genJet();
	    if (genJet)
	      {
		GenJetMatched.SetPtEtaPhiM(genJet->pt(), genJet->eta(), genJet->phi(), genJet->mass());
		if (verbose_) cout << "  ak8 genJet pt " << genJet->pt() << " mass " << genJet->mass() << endl;
	      }
	  }

	if (count_AK8CHS == 0) GenJetMatched0 = GenJetMatched;
	if (count_AK8CHS == 1) GenJetMatched1 = GenJetMatched;


	//------------------------------------
	// JER SF
	//------------------------------------
	double ptsmear   = 1;
	double ptsmearUp = 1;
	double ptsmearDn = 1;
	if (!iEvent.isRealData())
	  {
	    double jer_sf    = jer_scaler.getScaleFactor({{JME::Binning::JetEta, corrJet.eta()}});
	    double jer_sf_up = jer_scaler.getScaleFactor({{JME::Binning::JetEta, corrJet.eta()}}, Variation::UP);
	    double jer_sf_dn = jer_scaler.getScaleFactor({{JME::Binning::JetEta, corrJet.eta()}}, Variation::DOWN);
	    if (verbose_) std::cout << " JER Scale factors (Nominal / Up / Down) : " << jer_sf << " / " << jer_sf_up << " / " << jer_sf_dn << std::endl;
	    double recopt    = corrJet.pt();
	    double genpt     = GenJetMatched.Perp();
	    double deltapt   = (recopt-genpt)*(jer_sf-1.0);
	    double deltaptUp = (recopt-genpt)*(jer_sf_up-1.0);
	    double deltaptDn = (recopt-genpt)*(jer_sf_dn-1.0);
	    ptsmear   = std::max((double)0.0, (recopt+deltapt)/recopt);
	    ptsmearUp = std::max((double)0.0, (recopt+deltaptUp)/recopt);
	    ptsmearDn = std::max((double)0.0, (recopt+deltaptDn)/recopt);
	  }

	//------------------------------------
	// AK8 variables from miniAOD
	//------------------------------------
	// double pt           = corrJet.pt();
	//double mass         = corrJet.mass();
	// double eta          = corrJet.eta();
	// double phi          = corrJet.phi();
	//double rapidity     = ijet.rapidity();
	//double ndau         = ijet.numberOfDaughters();

	double tau1         = 99;
	double tau2         = 99;
	double tau3         = 99;
	double tau4         = 99;
	double prunedMass   = ijet.userFloat("ak8PFJetsCHSPrunedMass");
	double softDropMass = ijet.userFloat("ak8PFJetsCHSSoftDropMass");
	double trimmedMass  = -1.0;

	if (useToolbox_)
	  {
	    tau1         = ijet.userFloat("NjettinessAK8CHS:tau1");
	    tau2         = ijet.userFloat("NjettinessAK8CHS:tau2");
	    tau3         = ijet.userFloat("NjettinessAK8CHS:tau3");
	    tau4         = ijet.userFloat("NjettinessAK8CHS:tau4");
	    trimmedMass  = ijet.userFloat("ak8PFJetsCHSTrimmedMass"); 
	  }
	else
	  {
	    //tau1         = ijet.userFloat("NjettinessAK8:tau1");//COMMENTED OUT FOR REMOVING JET TOOL
	    //tau2         = ijet.userFloat("NjettinessAK8:tau2");
	    //tau3         = ijet.userFloat("NjettinessAK8:tau3");
	  }
	double tau21        = 99;
	double tau32        = 99;

	double puppi_p              = -1;     
	double puppi_pt             = -1;     
	double puppi_mass           = -1;     
	double puppi_eta            = -1;     
	double puppi_phi            = -1;     
	double puppi_tau1           = -1;     
	double puppi_tau2           = -1;     
	double puppi_tau3           = -1;     
	double puppi_tau4           = -1;     
	// double puppi_prunedMass     = -1;     
	// double puppi_trimmedMass    = -1;     
	// double puppi_softDropMass   = -1;     

	double puppi_CHF    = -1;
	double puppi_NHF    = -1;
	double puppi_CM     = -1;
	double puppi_NM     = -1;
	double puppi_NEF    = -1;
	double puppi_CEF    = -1;
	double puppi_MF     = -1;
	double puppi_Mult   = -1;
 
	if (!useToolbox_)
	  {
	    /*puppi_pt           = ijet.userFloat("ak8PFJetsPuppiValueMap:pt");//COMMENTED OUT FOR REMOVING TOOLBOX
	    puppi_mass         = ijet.userFloat("ak8PFJetsPuppiValueMap:mass");
	    puppi_eta          = ijet.userFloat("ak8PFJetsPuppiValueMap:eta");
	    puppi_phi          = ijet.userFloat("ak8PFJetsPuppiValueMap:phi");
	    puppi_tau1         = ijet.userFloat("ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau1");
	    puppi_tau2         = ijet.userFloat("ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau2");
	    puppi_tau3         = ijet.userFloat("ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau3");*/
	  }
	if (useToolbox_)
	  {
	    for (const pat::Jet &ipup : *AK8PUPPI)
	      {  
		if (ipup.pt() < 180) break;    
		double deltaRpup = deltaR(ijet.eta(), ijet.phi(), ipup.eta(), ipup.phi());
		if (deltaRpup<0.8)
		  {
		    puppi_p            = ipup.p();
		    puppi_pt           = ipup.pt();
		    puppi_mass         = ipup.mass();
		    puppi_eta          = ipup.eta();
		    puppi_phi          = ipup.phi();
		    // puppi_prunedMass   = ipup.userFloat("ak8PFJetsPuppiPrunedMass");
		    // puppi_trimmedMass  = ipup.userFloat("ak8PFJetsPuppiTrimmedMass");
		    // puppi_softDropMass = ipup.userFloat("ak8PFJetsPuppiSoftDropMass");
		    puppi_tau1         = ipup.userFloat("NjettinessAK8Puppi:tau1");
		    puppi_tau2         = ipup.userFloat("NjettinessAK8Puppi:tau2");
		    puppi_tau3         = ipup.userFloat("NjettinessAK8Puppi:tau3");
		    puppi_tau4         = ipup.userFloat("NjettinessAK8Puppi:tau4");

		    puppi_CHF          = ipup.chargedHadronEnergy() / ipup.correctedP4(0).E();  
		    puppi_NHF          = ipup.neutralHadronEnergy() / ipup.correctedP4(0).E();  
		    puppi_CM           = ipup.chargedMultiplicity();                  
		    puppi_NM           = ipup.neutralMultiplicity();                  
		    puppi_NEF          = ipup.neutralEmEnergy() / ipup.correctedP4(0).E();      
		    puppi_CEF          = ipup.chargedEmEnergy() / ipup.correctedP4(0).E();      
		    puppi_MF           = ipup.muonEnergy() / ipup.correctedP4(0).E();           
		    puppi_Mult         = ipup.numberOfDaughters();   

		  }
	      }
	  }
	double puppi_tau21        = 99.0;
	double puppi_tau32        = 99.0;


	if (tau1!=0) tau21 = tau2/tau1;
	if (tau2!=0) tau32 = tau3/tau2;

	if (puppi_tau1 != 0) puppi_tau21 = puppi_tau2/puppi_tau1;
	if (puppi_tau2 != 0) puppi_tau32 = puppi_tau3/puppi_tau2;

	TLorentzVector jet_p4;
	jet_p4.SetPtEtaPhiM(corrJet.pt(), corrJet.eta(), corrJet.phi(), corrJet.mass());
  

	// //------------------------------------
	// // AK8PUPPI JEC L23 correction
	// //------------------------------------

	// TLorentzVector jet_p4;
	// pupppijet_p4.SetPtEtaPhiM(puppi_pt, puppi_eta, puppi_phi, puppi_mass);
  

	// JetCorrectorAK8chs->setJetEta(puppi_eta);
	// JetCorrectorAK8chs->setJetPt (puppi_pt);
	// JetCorrectorAK8chs->setJetE  (puppi_e);
	// JetCorrectorAK8chs->setJetA  (ijet.jetArea());
	// JetCorrectorAK8chs->setRho   (rho);
	// JetCorrectorAK8chs->setNPV   (nvtx);
	// // getSubCorrections member function returns the vector of the subcorrections UP to the given level. For example in the example above, factors[0] is the L1 correction and factors[3] is the L1+L2+L3+Residual correction. 
	// vector<float> factors = JetCorrectorAK8chs->getSubCorrections();
	// float corr_factor_L1      = 1.0;
	// float corr_factor_L12     = 1.0;
	// float corr_factor_L123    = 1.0;
	// float corr_factor_L123res = 1.0;
	// if (factors.size() > 0) corr_factor_L1       = factors[0];
	// if (factors.size() > 1) corr_factor_L12      = factors[1];
	// if (factors.size() > 2) corr_factor_L123     = factors[2];
	// if (factors.size() > 3) corr_factor_L123res  = factors[3];
	// double corr_factor_L2 = corr_factor_L12/corr_factor_L1;
	// double corr_factor_L3 = corr_factor_L123/corr_factor_L12;
	// double corr_factor_res = corr_factor_L123res/corr_factor_L123;
	// //double corr_factor_L23 = corr_factor_L2*corr_factor_L3;
	// double corr_factor_L23res = corr_factor_L2*corr_factor_L3*corr_factor_res;


	if (count_AK8CHS == 0) PUPPIjet0_P4corr.SetPtEtaPhiM(puppi_pt, puppi_eta, puppi_phi, puppi_mass);
	if (count_AK8CHS == 1) PUPPIjet1_P4corr.SetPtEtaPhiM(puppi_pt, puppi_eta, puppi_phi, puppi_mass);



	//------------------------------------
	// SoftDrop subjets
	//------------------------------------
	TLorentzVector sub0_P4_uncorr;
	TLorentzVector sub0_P4_L23res;
	TLorentzVector sub0_P4_L23resCorrUp;
	TLorentzVector sub0_P4_L23resCorrDn;
	TLorentzVector sub0_P4_L23resPtSmear;
	TLorentzVector sub0_P4_L23resPtSmearUp;
	TLorentzVector sub0_P4_L23resPtSmearDn;
	TLorentzVector sub0_P4_L123res;
	TLorentzVector sub0_P4_L123resCorrUp;
	TLorentzVector sub0_P4_L123resCorrDn;

	TLorentzVector sub1_P4_uncorr;
	TLorentzVector sub1_P4_L23res;
	TLorentzVector sub1_P4_L23resCorrUp;
	TLorentzVector sub1_P4_L23resCorrDn;
	TLorentzVector sub1_P4_L23resPtSmear;
	TLorentzVector sub1_P4_L23resPtSmearUp;
	TLorentzVector sub1_P4_L23resPtSmearDn;
	TLorentzVector sub1_P4_L123res;
	TLorentzVector sub1_P4_L123resCorrUp;
	TLorentzVector sub1_P4_L123resCorrDn;

	double sub0_area  = 0.0;
	double sub0_tau1  = 0.0;
	double sub0_tau2  = 0.0;
	double sub0_tau3  = 0.0;
	double sub0_flav_hadron  = 0.0;
	double sub0_flav_parton  = 0.0;
	double sub0_bdisc = 0.0;
	double sub1_area  = 0.0;
	double sub1_tau1  = 0.0;
	double sub1_tau2  = 0.0;
	double sub1_tau3  = 0.0;
	double sub1_flav_hadron  = 0.0;
	double sub1_flav_parton  = 0.0;
	double sub1_bdisc = 0.0;
	double mostMassiveSDsubjetMass = 0.0;
	int count_SD = 0;

	if (false) //COMMENTED OUT FOR NOT USING JET TOOLBOX (!useToolbox_)
	  {
	    auto const & sdSubjets = ijet.subjets("SoftDrop");
	    for (auto const & it : sdSubjets)
	      {
		double subjetPt       = it->correctedP4(0).pt();
		double subjetEta      = it->correctedP4(0).eta();
		double subjetPhi      = it->correctedP4(0).phi();
		double subjetMass     = it->correctedP4(0).mass();
		double subjetBdisc    = it->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"); 
		double deltaRsubjetJet = deltaR(ijet.eta(), ijet.phi(), subjetEta, subjetPhi);
      
		//------------------------------------
		// subjet JEC 
		//------------------------------------
		reco::Candidate::LorentzVector uncorrSubjet = it->correctedP4(0);
		JetCorrectorAK4chs->setJetEta(uncorrSubjet.eta());
		JetCorrectorAK4chs->setJetPt (uncorrSubjet.pt());
		JetCorrectorAK4chs->setJetE  (uncorrSubjet.energy());
		JetCorrectorAK4chs->setJetA  (it->jetArea());
		JetCorrectorAK4chs->setRho   (rho);
		JetCorrectorAK4chs->setNPV   (nvtx);
		double subjet_corr_factor_L123res_full = JetCorrectorAK4chs->getCorrection();
		reco::Candidate::LorentzVector corrSubjetL123res = subjet_corr_factor_L123res_full * uncorrSubjet;

		//------------------------------------
		// subjet L23 JEC 
		//------------------------------------
		JetCorrectorAK4chs->setJetEta(uncorrSubjet.eta());
		JetCorrectorAK4chs->setJetPt (uncorrSubjet.pt());
		JetCorrectorAK4chs->setJetE  (uncorrSubjet.energy());
		JetCorrectorAK4chs->setJetA  (it->jetArea());
		JetCorrectorAK4chs->setRho   (rho);
		JetCorrectorAK4chs->setNPV   (nvtx);
		// getSubCorrections member function returns the vector of the subcorrections UP to the given level. For example in the example above, factors[0] is the L1 correction and factors[3] is the L1+L2+L3+Residual correction. 
		vector<float> subjet_factors = JetCorrectorAK4chs->getSubCorrections();
		float subjet_corr_factor_L1      = 1.0;
		float subjet_corr_factor_L12     = 1.0;
		float subjet_corr_factor_L123    = 1.0;
		float subjet_corr_factor_L123res = 1.0;
		if (factors.size() > 0) subjet_corr_factor_L1      = subjet_factors[0];
		if (factors.size() > 1) subjet_corr_factor_L12     = subjet_factors[1];
		if (factors.size() > 2) subjet_corr_factor_L123    = subjet_factors[2];
		if (factors.size() > 3) subjet_corr_factor_L123res = subjet_factors[3];
		double subjet_corr_factor_L2     = subjet_corr_factor_L12     / subjet_corr_factor_L1;
		double subjet_corr_factor_L3     = subjet_corr_factor_L123    / subjet_corr_factor_L12;
		double subjet_corr_factor_res    = subjet_corr_factor_L123res / subjet_corr_factor_L123;
		double subjet_corr_factor_L23    = subjet_corr_factor_L2 * subjet_corr_factor_L3;
		double subjet_corr_factor_L23res = subjet_corr_factor_L2 * subjet_corr_factor_L3 * subjet_corr_factor_res;
		if (verbose_) cout << "subjet corr: L1 " << subjet_corr_factor_L1 << " L23 " << subjet_corr_factor_L23<< " L23res " << subjet_corr_factor_L23res << " L123res" << subjet_corr_factor_L123res << endl;
		reco::Candidate::LorentzVector corrSubjetL23res = subjet_corr_factor_L23res * uncorrSubjet;
        
		//------------------------------------
		// subjet JEC uncertainty
		//------------------------------------
		double subjet_corrDn_L23 =  1.0;
		double subjet_corrDn_L123 = 1.0;
		JetCorrUncertAK4chs->setJetPhi(corrSubjetL123res.phi());
		JetCorrUncertAK4chs->setJetEta(corrSubjetL123res.eta());
		JetCorrUncertAK4chs->setJetPt(corrSubjetL123res.pt());
		double corrDn_temp2 = JetCorrUncertAK4chs->getUncertainty(0);
		subjet_corrDn_L23   = subjet_corr_factor_L23res - corrDn_temp2;
		subjet_corrDn_L123  = subjet_corr_factor_L123res_full - corrDn_temp2;

		double subjet_corrUp_L23  = 1.0;
		double subjet_corrUp_L123 = 1.0;
		JetCorrUncertAK4chs->setJetPhi(corrSubjetL123res.phi());
		JetCorrUncertAK4chs->setJetEta(corrSubjetL123res.eta());
		JetCorrUncertAK4chs->setJetPt(corrSubjetL123res.pt());
		double corrUp_temp2 = JetCorrUncertAK4chs->getUncertainty(1);
		subjet_corrUp_L23   = subjet_corr_factor_L23res + corrUp_temp2;
		subjet_corrUp_L123  = subjet_corr_factor_L123res_full + corrUp_temp2;

		reco::Candidate::LorentzVector corrSubjetL123resCorrDn  = subjet_corrDn_L123  * uncorrSubjet;
		reco::Candidate::LorentzVector corrSubjetL123resCorrUp  = subjet_corrUp_L123  * uncorrSubjet;
		reco::Candidate::LorentzVector corrSubjetL23resCorrDn   = subjet_corrDn_L23   * uncorrSubjet;
		reco::Candidate::LorentzVector corrSubjetL23resCorrUp   = subjet_corrUp_L23   * uncorrSubjet;
     

		//------------------------------------
		// subjet values for Tree
		//------------------------------------
		if (count_SD == 0)
		  {
		    sub0_P4_uncorr            .SetPtEtaPhiM(subjetPt, subjetEta, subjetPhi, subjetMass);
		    sub0_P4_L123res           .SetPtEtaPhiM(corrSubjetL123res.pt()   , corrSubjetL123res.eta()   , corrSubjetL123res.phi()   , corrSubjetL123res.mass());
		    sub0_P4_L23res            .SetPtEtaPhiM(corrSubjetL23res.pt()    , corrSubjetL23res.eta()    , corrSubjetL23res.phi()    , corrSubjetL23res.mass());
		    sub0_P4_L123resCorrUp    .SetPtEtaPhiM(corrSubjetL123resCorrUp.pt() , corrSubjetL123resCorrUp.eta() , corrSubjetL123resCorrUp.phi() , corrSubjetL123resCorrUp.mass());
		    sub0_P4_L23resCorrUp     .SetPtEtaPhiM(corrSubjetL23resCorrUp.pt()  , corrSubjetL23resCorrUp.eta()  , corrSubjetL23resCorrUp.phi()  , corrSubjetL23resCorrUp.mass());
		    sub0_P4_L123resCorrDn    .SetPtEtaPhiM(corrSubjetL123resCorrDn.pt() , corrSubjetL123resCorrDn.eta() , corrSubjetL123resCorrDn.phi() , corrSubjetL123resCorrUp.mass());
		    sub0_P4_L23resCorrDn     .SetPtEtaPhiM(corrSubjetL23resCorrDn.pt()  , corrSubjetL23resCorrDn.eta()  , corrSubjetL23resCorrDn.phi()  , corrSubjetL23res.mass());
		    sub0_area   = it->jetArea();
		    sub0_flav_parton   = it->partonFlavour();
		    sub0_flav_hadron   = it->hadronFlavour();
		    sub0_bdisc  = subjetBdisc;
		  }
		if (count_SD == 1)
		  {
		    sub1_P4_uncorr          .SetPtEtaPhiM(subjetPt, subjetEta, subjetPhi, subjetMass);
		    sub1_P4_L123res         .SetPtEtaPhiM(corrSubjetL123res.pt()   , corrSubjetL123res.eta()   , corrSubjetL123res.phi()   , corrSubjetL123res.mass());
		    sub1_P4_L23res          .SetPtEtaPhiM(corrSubjetL23res.pt()    , corrSubjetL23res.eta()    , corrSubjetL23res.phi()    , corrSubjetL23res.mass());
		    sub1_P4_L123resCorrUp  .SetPtEtaPhiM(corrSubjetL123resCorrUp.pt() , corrSubjetL123resCorrUp.eta() , corrSubjetL123resCorrUp.phi() , corrSubjetL123resCorrUp.mass());
		    sub1_P4_L23resCorrUp   .SetPtEtaPhiM(corrSubjetL23resCorrUp.pt()  , corrSubjetL23resCorrUp.eta()  , corrSubjetL23resCorrUp.phi()  , corrSubjetL23resCorrUp.mass());
		    sub1_P4_L123resCorrDn  .SetPtEtaPhiM(corrSubjetL123resCorrDn.pt() , corrSubjetL123resCorrDn.eta() , corrSubjetL123resCorrDn.phi() , corrSubjetL123resCorrUp.mass());
		    sub1_P4_L23resCorrDn   .SetPtEtaPhiM(corrSubjetL23resCorrDn.pt()  , corrSubjetL23resCorrDn.eta()  , corrSubjetL23resCorrDn.phi()  , corrSubjetL23res.mass());
		    sub1_area   = it->jetArea();
		    sub1_flav_parton   = it->partonFlavour();
		    sub1_flav_hadron   = it->hadronFlavour();
		    sub1_bdisc  = subjetBdisc;
		  }
		if (subjetMass > mostMassiveSDsubjetMass) mostMassiveSDsubjetMass = subjetMass;

		if (verbose_) cout << " SD Subjet pt " << subjetPt << " Eta " << subjetEta << " deltaRsubjetJet " << deltaRsubjetJet << " Mass " << subjetMass << " Bdisc " << subjetBdisc << endl;
		count_SD++;
	      }
	  }
	if (useToolbox_)
	  {
	    for (const pat::Jet &isub : *AK8CHSsub)
	      {  
  
		double subjetPt       = isub.correctedP4(0).pt();
		double subjetEta      = isub.correctedP4(0).eta();
		double subjetPhi      = isub.correctedP4(0).phi();
		double subjetMass     = isub.correctedP4(0).mass();
		double subjetBdisc    = isub.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"); 

		double deltaRsubjetJet = deltaR(ijet.eta(), ijet.phi(), subjetEta, subjetPhi);
		if (deltaRsubjetJet < 0.8)
		  {
		    if (verbose_) cout << " matched subjet with mass " << subjetMass << endl;

		    //------------------------------------
		    // subjet JEC 
		    //------------------------------------
		    reco::Candidate::LorentzVector uncorrSubjet = isub.correctedP4(0);
		    JetCorrectorAK4chs->setJetEta(uncorrSubjet.eta());
		    JetCorrectorAK4chs->setJetPt (uncorrSubjet.pt());
		    JetCorrectorAK4chs->setJetE  (uncorrSubjet.energy());
		    JetCorrectorAK4chs->setJetA  (isub.jetArea());
		    JetCorrectorAK4chs->setRho   (rho);
		    JetCorrectorAK4chs->setNPV   (nvtx);
		    double subjet_corr_factor_L123res_full = JetCorrectorAK4chs->getCorrection();
		    reco::Candidate::LorentzVector corrSubjetL123res = subjet_corr_factor_L123res_full * uncorrSubjet;

		    //------------------------------------
		    // subjet L23 JEC 
		    //------------------------------------
		    JetCorrectorAK4chs->setJetEta(uncorrSubjet.eta());
		    JetCorrectorAK4chs->setJetPt (uncorrSubjet.pt());
		    JetCorrectorAK4chs->setJetE  (uncorrSubjet.energy());
		    JetCorrectorAK4chs->setJetA  (isub.jetArea());
		    JetCorrectorAK4chs->setRho   (rho);
		    JetCorrectorAK4chs->setNPV   (nvtx);
		    // getSubCorrections member function returns the vector of the subcorrections UP to the given level. For example in the example above, factors[0] is the L1 correction and factors[3] is the L1+L2+L3+Residual correction. 
		    vector<float> subjet_factors = JetCorrectorAK4chs->getSubCorrections();
		    float subjet_corr_factor_L1      = 1.0;
		    float subjet_corr_factor_L12     = 1.0;
		    float subjet_corr_factor_L123    = 1.0;
		    float subjet_corr_factor_L123res = 1.0;
		    if (factors.size() > 0) subjet_corr_factor_L1      = subjet_factors[0];
		    if (factors.size() > 1) subjet_corr_factor_L12     = subjet_factors[1];
		    if (factors.size() > 2) subjet_corr_factor_L123    = subjet_factors[2];
		    if (factors.size() > 3) subjet_corr_factor_L123res = subjet_factors[3];
		    double subjet_corr_factor_L2     = subjet_corr_factor_L12     / subjet_corr_factor_L1;
		    double subjet_corr_factor_L3     = subjet_corr_factor_L123    / subjet_corr_factor_L12;
		    double subjet_corr_factor_res    = subjet_corr_factor_L123res / subjet_corr_factor_L123;
		    double subjet_corr_factor_L23    = subjet_corr_factor_L2 * subjet_corr_factor_L3;
		    double subjet_corr_factor_L23res = subjet_corr_factor_L2 * subjet_corr_factor_L3 * subjet_corr_factor_res;
		    if (verbose_) cout << "  subjet corr: L1 " << subjet_corr_factor_L1 << " L23 " << subjet_corr_factor_L23 << " L23res " << subjet_corr_factor_L23res << " L123res " << subjet_corr_factor_L123res << endl;
		    reco::Candidate::LorentzVector corrSubjetL23res   = subjet_corr_factor_L23res * uncorrSubjet;
          
		    //------------------------------------
		    // subjet JEC uncertainty
		    //------------------------------------
		    double subjet_corrDn_L23 =  1.0;
		    double subjet_corrDn_L123 = 1.0;
		    JetCorrUncertAK4chs->setJetPhi(corrSubjetL123res.phi());
		    JetCorrUncertAK4chs->setJetEta(corrSubjetL123res.eta());
		    JetCorrUncertAK4chs->setJetPt(corrSubjetL123res.pt());
		    double corrDn_temp2 = JetCorrUncertAK4chs->getUncertainty(0);
		    subjet_corrDn_L23   = subjet_corr_factor_L23res - corrDn_temp2;
		    subjet_corrDn_L123  = subjet_corr_factor_L123res_full - corrDn_temp2;

		    double subjet_corrUp_L23  = 1.0;
		    double subjet_corrUp_L123 = 1.0;
		    JetCorrUncertAK4chs->setJetPhi(corrSubjetL123res.phi());
		    JetCorrUncertAK4chs->setJetEta(corrSubjetL123res.eta());
		    JetCorrUncertAK4chs->setJetPt(corrSubjetL123res.pt());
		    double corrUp_temp2 = JetCorrUncertAK4chs->getUncertainty(1);
		    subjet_corrUp_L23   = subjet_corr_factor_L23res + corrUp_temp2;
		    subjet_corrUp_L123  = subjet_corr_factor_L123res_full + corrUp_temp2;

		    reco::Candidate::LorentzVector corrSubjetL123resCorrDn  = subjet_corrDn_L123  * uncorrSubjet;
		    reco::Candidate::LorentzVector corrSubjetL123resCorrUp  = subjet_corrUp_L123  * uncorrSubjet;
		    reco::Candidate::LorentzVector corrSubjetL23resCorrDn   = subjet_corrDn_L23   * uncorrSubjet;
		    reco::Candidate::LorentzVector corrSubjetL23resCorrUp   = subjet_corrUp_L23   * uncorrSubjet;
         
		    //------------------------------------
		    // subjet JER SF
		    //------------------------------------
		    TLorentzVector GenSubJet;
		    double ptsmear   = 1;
		    double ptsmearUp = 1;
		    double ptsmearDn = 1;
		    if (!iEvent.isRealData())
		      {
			const reco::GenJet* genJet = isub.genJet();
			if (genJet)
			  {
			    GenSubJet.SetPtEtaPhiM(genJet->pt(), genJet->eta(), genJet->phi(), genJet->mass());
			    if (verbose_) cout << "  SD subjet genJet pt " << genJet->pt() << " mass " << genJet->mass() << " reco pt " << subjetPt << " reco mass " << subjetMass << endl;
			  }
			double jer_sf    = jer_scaler.getScaleFactor({{JME::Binning::JetEta, corrSubjetL23res.eta()}});
			double jer_sf_up = jer_scaler.getScaleFactor({{JME::Binning::JetEta, corrSubjetL23res.eta()}}, Variation::UP);
			double jer_sf_dn = jer_scaler.getScaleFactor({{JME::Binning::JetEta, corrSubjetL23res.eta()}}, Variation::DOWN);
			if (verbose_) std::cout << " SD subjet JER Scale factors (Nominal / Up / Down) : " << jer_sf << " / " << jer_sf_up << " / " << jer_sf_dn << std::endl;
			double recopt    = corrSubjetL23res.pt();
			double genpt     = GenJetMatched.Perp();
			double deltapt   = (recopt-genpt)*(jer_sf-1.0);
			double deltaptUp = (recopt-genpt)*(jer_sf_up-1.0);
			double deltaptDn = (recopt-genpt)*(jer_sf_dn-1.0);
			ptsmear   = std::max((double)0.0, (recopt+deltapt)/recopt);
			ptsmearUp = std::max((double)0.0, (recopt+deltaptUp)/recopt);
			ptsmearDn = std::max((double)0.0, (recopt+deltaptDn)/recopt);
			if (verbose_) std::cout << " SD subjet ptsmear " << ptsmear << " ptsmearUp " << ptsmearUp << " ptsmearDn " << ptsmearDn << endl;
		      }
		    reco::Candidate::LorentzVector corrSubjetL23resPtSmear   = ptsmear * corrSubjetL23res;
		    reco::Candidate::LorentzVector corrSubjetL23resPtSmearUp = ptsmearUp * corrSubjetL23res;
		    reco::Candidate::LorentzVector corrSubjetL23resPtSmearDn = ptsmearDn * corrSubjetL23res;

		    //------------------------------------
		    // subjet values for Tree
		    //------------------------------------
		    if (count_SD == 0)
		      {
			sub0_P4_uncorr            .SetPtEtaPhiM(subjetPt, subjetEta, subjetPhi, subjetMass);
			sub0_P4_L23res            .SetPtEtaPhiM(corrSubjetL23res          .pt() , corrSubjetL23res          .eta()  , corrSubjetL23res          .phi()  , corrSubjetL23res          .mass());
			sub0_P4_L23resCorrUp      .SetPtEtaPhiM(corrSubjetL23resCorrUp    .pt() , corrSubjetL23resCorrUp    .eta()  , corrSubjetL23resCorrUp    .phi()  , corrSubjetL23resCorrUp    .mass());
			sub0_P4_L23resCorrDn      .SetPtEtaPhiM(corrSubjetL23resCorrDn    .pt() , corrSubjetL23resCorrDn    .eta()  , corrSubjetL23resCorrDn    .phi()  , corrSubjetL23resCorrDn    .mass());
			sub0_P4_L23resPtSmear     .SetPtEtaPhiM(corrSubjetL23resPtSmear   .pt() , corrSubjetL23resPtSmear   .eta()  , corrSubjetL23resPtSmear   .phi()  , corrSubjetL23resPtSmear   .mass());
			sub0_P4_L23resPtSmearUp   .SetPtEtaPhiM(corrSubjetL23resPtSmearUp .pt() , corrSubjetL23resPtSmearUp .eta()  , corrSubjetL23resPtSmearUp .phi()  , corrSubjetL23resPtSmearUp .mass());
			sub0_P4_L23resPtSmearDn   .SetPtEtaPhiM(corrSubjetL23resPtSmearDn .pt() , corrSubjetL23resPtSmearDn .eta()  , corrSubjetL23resPtSmearDn .phi()  , corrSubjetL23resPtSmearDn .mass());
			sub0_P4_L123res           .SetPtEtaPhiM(corrSubjetL123res         .pt() , corrSubjetL123res         .eta()  , corrSubjetL123res         .phi()  , corrSubjetL123res         .mass());
			sub0_P4_L123resCorrUp     .SetPtEtaPhiM(corrSubjetL123resCorrUp   .pt() , corrSubjetL123resCorrUp   .eta()  , corrSubjetL123resCorrUp   .phi()  , corrSubjetL123resCorrUp   .mass());
			sub0_P4_L123resCorrDn     .SetPtEtaPhiM(corrSubjetL123resCorrDn   .pt() , corrSubjetL123resCorrDn   .eta()  , corrSubjetL123resCorrDn   .phi()  , corrSubjetL123resCorrDn   .mass());
 
			sub0_area          = isub.jetArea();
			sub0_flav_parton   = isub.partonFlavour();
			sub0_flav_hadron   = isub.hadronFlavour();
			sub0_bdisc         = subjetBdisc;
			// available from toolbox only (80X) THESE MUST BE COMMENTED OUT FOR NOT DEFINING TOOL BOX
			//sub0_tau1          = isub.userFloat("NsubjettinessAK8PFCHSSoftDropSubjets:tau1");
			//sub0_tau2          = isub.userFloat("NsubjettinessAK8PFCHSSoftDropSubjets:tau2");
			//sub0_tau3          = isub.userFloat("NsubjettinessAK8PFCHSSoftDropSubjets:tau3");
		      }
		    if (count_SD == 1)
		      {
			sub1_P4_uncorr            .SetPtEtaPhiM(subjetPt, subjetEta, subjetPhi, subjetMass);
			sub1_P4_L23res            .SetPtEtaPhiM(corrSubjetL23res          .pt() , corrSubjetL23res          .eta()  , corrSubjetL23res          .phi()  , corrSubjetL23res          .mass());
			sub1_P4_L23resCorrUp      .SetPtEtaPhiM(corrSubjetL23resCorrUp    .pt() , corrSubjetL23resCorrUp    .eta()  , corrSubjetL23resCorrUp    .phi()  , corrSubjetL23resCorrUp    .mass());
			sub1_P4_L23resCorrDn      .SetPtEtaPhiM(corrSubjetL23resCorrDn    .pt() , corrSubjetL23resCorrDn    .eta()  , corrSubjetL23resCorrDn    .phi()  , corrSubjetL23resCorrDn    .mass());
			sub1_P4_L23resPtSmear     .SetPtEtaPhiM(corrSubjetL23resPtSmear   .pt() , corrSubjetL23resPtSmear   .eta()  , corrSubjetL23resPtSmear   .phi()  , corrSubjetL23resPtSmear   .mass());
			sub1_P4_L23resPtSmearUp   .SetPtEtaPhiM(corrSubjetL23resPtSmearUp .pt() , corrSubjetL23resPtSmearUp .eta()  , corrSubjetL23resPtSmearUp .phi()  , corrSubjetL23resPtSmearUp .mass());
			sub1_P4_L23resPtSmearDn   .SetPtEtaPhiM(corrSubjetL23resPtSmearDn .pt() , corrSubjetL23resPtSmearDn .eta()  , corrSubjetL23resPtSmearDn .phi()  , corrSubjetL23resPtSmearDn .mass());
			sub1_P4_L123res           .SetPtEtaPhiM(corrSubjetL123res         .pt() , corrSubjetL123res         .eta()  , corrSubjetL123res         .phi()  , corrSubjetL123res         .mass());
			sub1_P4_L123resCorrUp     .SetPtEtaPhiM(corrSubjetL123resCorrUp   .pt() , corrSubjetL123resCorrUp   .eta()  , corrSubjetL123resCorrUp   .phi()  , corrSubjetL123resCorrUp   .mass());
			sub1_P4_L123resCorrDn     .SetPtEtaPhiM(corrSubjetL123resCorrDn   .pt() , corrSubjetL123resCorrDn   .eta()  , corrSubjetL123resCorrDn   .phi()  , corrSubjetL123resCorrDn   .mass());
 
			sub1_area          = isub.jetArea();
			sub1_flav_parton   = isub.partonFlavour();
			sub1_flav_hadron   = isub.hadronFlavour();
			sub1_bdisc         = subjetBdisc;
			// available from toolbox only (80X) THESE MUST BE COMMENTED OUT FOR NOT DEFINING TOOL BOX
			//sub1_tau1          = isub.userFloat("NsubjettinessAK8PFCHSSoftDropSubjets:tau1");
			//sub1_tau2          = isub.userFloat("NsubjettinessAK8PFCHSSoftDropSubjets:tau2");
			//sub1_tau3          = isub.userFloat("NsubjettinessAK8PFCHSSoftDropSubjets:tau3");
		      }
		    if (subjetMass > mostMassiveSDsubjetMass) mostMassiveSDsubjetMass = subjetMass;

		    if (verbose_) cout << "  SD Subjet pt " << subjetPt << " Eta " << subjetEta << " deltaRsubjetJet " << deltaRsubjetJet << " Mass " << subjetMass << " corrMass " << corrSubjetL23res.mass() << " Bdisc " << subjetBdisc << endl;
		    if (verbose_) cout << "    sub0_tau1 " << sub0_tau1 << " sub0_tau2 " << sub0_tau2<< " sub0_tau3 " << sub0_tau3 << endl;
		    count_SD++;

		  }
	      }
	  }

	TLorentzVector sumSDsubjets_P4_uncorr;
	TLorentzVector sumSDsubjets_P4_L23res;
	TLorentzVector sumSDsubjets_P4_L23resCorrUp;
	TLorentzVector sumSDsubjets_P4_L23resCorrDn;
	TLorentzVector sumSDsubjets_P4_L23resPtSmear;
	TLorentzVector sumSDsubjets_P4_L23resPtSmearUp;
	TLorentzVector sumSDsubjets_P4_L23resPtSmearDn;
	TLorentzVector sumSDsubjets_P4_L123res;
	TLorentzVector sumSDsubjets_P4_L123resCorrDn;
	TLorentzVector sumSDsubjets_P4_L123resCorrUp;

	if (count_SD > 1)
	  { 
	    sumSDsubjets_P4_uncorr             = sub0_P4_uncorr              + sub1_P4_uncorr; 
	    sumSDsubjets_P4_L23res             = sub0_P4_L23res              + sub1_P4_L23res; 
	    sumSDsubjets_P4_L23resCorrUp       = sub0_P4_L23resCorrUp        + sub1_P4_L23resCorrUp; 
	    sumSDsubjets_P4_L23resCorrDn       = sub0_P4_L23resCorrDn        + sub1_P4_L23resCorrDn; 
	    sumSDsubjets_P4_L23resPtSmear      = sub0_P4_L23resPtSmear       + sub1_P4_L23resPtSmear;
	    sumSDsubjets_P4_L23resPtSmearUp    = sub0_P4_L23resPtSmearUp     + sub1_P4_L23resPtSmearUp;
	    sumSDsubjets_P4_L23resPtSmearDn    = sub0_P4_L23resPtSmearDn     + sub1_P4_L23resPtSmearDn;
	    sumSDsubjets_P4_L123res            = sub0_P4_L123res             + sub1_P4_L123res; 
	    sumSDsubjets_P4_L123resCorrUp      = sub0_P4_L123resCorrUp       + sub1_P4_L123resCorrUp; 
	    sumSDsubjets_P4_L123resCorrDn      = sub0_P4_L123resCorrDn       + sub1_P4_L123resCorrDn; 
	  }  

	double maxbdisc = 0;
	double maxbdiscflav_hadron = 0;
	double maxbdiscflav_parton = 0;
	if (sub0_bdisc>=sub1_bdisc)
	  {
	    maxbdisc = sub0_bdisc;
	    maxbdiscflav_hadron = sub0_flav_hadron;
	    maxbdiscflav_parton = sub0_flav_parton;
	  } 
	else if (sub1_bdisc>sub0_bdisc)
	  {
	    maxbdisc = sub1_bdisc;
	    maxbdiscflav_hadron = sub1_flav_hadron;
	    maxbdiscflav_parton = sub1_flav_parton;
	  }  

	//------------------------------------
	// PUPPI SoftDrop subjets
	//------------------------------------ 
	TLorentzVector pup0_P4_uncorr;
	TLorentzVector pup0_P4_L23res;
	TLorentzVector pup0_P4_L23resCorrUp;
	TLorentzVector pup0_P4_L23resCorrDn;
	TLorentzVector pup0_P4_L23resPtSmear;
	TLorentzVector pup0_P4_L23resPtSmearUp;
	TLorentzVector pup0_P4_L23resPtSmearDn;

	TLorentzVector pup1_P4_uncorr;
	TLorentzVector pup1_P4_L23res;
	TLorentzVector pup1_P4_L23resCorrUp;
	TLorentzVector pup1_P4_L23resCorrDn;
	TLorentzVector pup1_P4_L23resPtSmear;
	TLorentzVector pup1_P4_L23resPtSmearUp;
	TLorentzVector pup1_P4_L23resPtSmearDn;


	double pup0_area  = 0.0;
	double pup0_tau1  = 0.0;
	double pup0_tau2  = 0.0;
	double pup0_tau3  = 0.0;
	double pup0_flav_hadron = 0.0;
	double pup0_flav_parton = 0.0;
	double pup0_bdisc = 0.0;
	double pup1_area  = 0.0;
	double pup1_tau1  = 0.0;
	double pup1_tau2  = 0.0;
	double pup1_tau3  = 0.0;
	double pup1_flav_hadron = 0.0;
	double pup1_flav_parton = 0.0;
	double pup1_bdisc = 0.0;
	double mostMassiveSDPUPPIsubjetMass = 0.0;
	int count_pup = 0;

	if (!true)//COMMENTED OUT FOR NOT USING TOOL BOX(!useToolbox_)
	  {
	    auto const &sdSubjetsPuppi = ijet.subjets("SoftDropPuppi");
	    for (auto const &it : sdSubjetsPuppi)
	      {
		double subjetPt       = it->correctedP4(0).pt();
		double subjetEta      = it->correctedP4(0).eta();
		double subjetPhi      = it->correctedP4(0).phi();
		double subjetMass     = it->correctedP4(0).mass();
		double subjetBdisc    = it->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"); 
		double deltaRsubjetJet = deltaR(ijet.eta(), ijet.phi(), subjetEta, subjetPhi);
		if (verbose_) cout << " SD Subjet pt " << subjetPt << " Eta " << subjetEta << " deltaRsubjetJet " << deltaRsubjetJet << " Mass " << subjetMass << " Bdisc " << subjetBdisc << endl; 
        
		//------------------------------------
		// PUPPI subjet JEC 
		//------------------------------------
		reco::Candidate::LorentzVector uncorrSubjet = it->correctedP4(0);
		JetCorrectorAK4pup->setJetEta(uncorrSubjet.eta());
		JetCorrectorAK4pup->setJetPt (uncorrSubjet.pt());
		JetCorrectorAK4pup->setJetE  (uncorrSubjet.energy());
		JetCorrectorAK4pup->setJetA  (it->jetArea());
		JetCorrectorAK4pup->setRho   (rho);
		JetCorrectorAK4pup->setNPV   (nvtx);
		double subjet_corr_factor_L23res_full = JetCorrectorAK4pup->getCorrection();
		reco::Candidate::LorentzVector corrSubjetL23res = subjet_corr_factor_L23res_full * uncorrSubjet;

		//------------------------------------
		// PUPPI subjet JEC uncertainty
		//------------------------------------
		double subjet_corrDn_L23 = 1.0;
		JetCorrUncertAK4pup->setJetPhi(corrSubjetL23res.phi());
		JetCorrUncertAK4pup->setJetEta(corrSubjetL23res.eta());
		JetCorrUncertAK4pup->setJetPt(corrSubjetL23res.pt());
		subjet_corrDn_L23 = subjet_corr_factor_L23res_full - JetCorrUncertAK4pup->getUncertainty(0);
		double subjet_corrUp_L23 = 1.0;
		JetCorrUncertAK4pup->setJetPhi(corrSubjetL23res.phi());
		JetCorrUncertAK4pup->setJetEta(corrSubjetL23res.eta());
		JetCorrUncertAK4pup->setJetPt(corrSubjetL23res.pt());
		subjet_corrUp_L23 = subjet_corr_factor_L23res_full + JetCorrUncertAK4pup->getUncertainty(1);

		reco::Candidate::LorentzVector corrSubjetL23resCorrDn   = subjet_corrDn_L23   * uncorrSubjet;
		reco::Candidate::LorentzVector corrSubjetL23resCorrUp   = subjet_corrUp_L23   * uncorrSubjet;
     
		//------------------------------------
		// subjet values for Tree
		//------------------------------------

		if (count_pup == 0)
		  {
		    pup0_P4_uncorr          .SetPtEtaPhiM(subjetPt, subjetEta, subjetPhi, subjetMass);
		    pup0_P4_L23res          .SetPtEtaPhiM(corrSubjetL23res.pt()    , corrSubjetL23res.eta()    , corrSubjetL23res.phi()    , corrSubjetL23res.mass());
		    pup0_P4_L23resCorrUp    .SetPtEtaPhiM(corrSubjetL23resCorrUp.pt()  , corrSubjetL23resCorrUp.eta()  , corrSubjetL23resCorrUp.phi()  , corrSubjetL23resCorrUp.mass());
		    pup0_P4_L23resCorrDn    .SetPtEtaPhiM(corrSubjetL23resCorrDn.pt()  , corrSubjetL23resCorrDn.eta()  , corrSubjetL23resCorrDn.phi()  , corrSubjetL23res.mass());
		    pup0_area = it->jetArea();
		    if (useToolbox_)
		      {
			pup0_tau1 = it->userFloat("ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau1");
			pup0_tau2 = it->userFloat("ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau2");
			pup0_tau3 = it->userFloat("ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau3");
		      }
		    pup0_flav_parton = it->partonFlavour();
		    pup0_flav_hadron = it->hadronFlavour();
		    pup0_bdisc  = subjetBdisc;
		  }
		if (count_pup == 1)
		  {
		    pup1_P4_uncorr          .SetPtEtaPhiM(subjetPt, subjetEta, subjetPhi, subjetMass);
		    pup1_P4_L23res          .SetPtEtaPhiM(corrSubjetL23res.pt()    , corrSubjetL23res.eta()    , corrSubjetL23res.phi()    , corrSubjetL23res.mass());
		    pup1_P4_L23resCorrUp    .SetPtEtaPhiM(corrSubjetL23resCorrUp.pt()  , corrSubjetL23resCorrUp.eta()  , corrSubjetL23resCorrUp.phi()  , corrSubjetL23resCorrUp.mass());
		    pup1_P4_L23resCorrDn    .SetPtEtaPhiM(corrSubjetL23resCorrDn.pt()  , corrSubjetL23resCorrDn.eta()  , corrSubjetL23resCorrDn.phi()  , corrSubjetL23res.mass());
		    pup1_area = it->jetArea();
		    if (useToolbox_)
		      {
			pup1_tau1 = it->userFloat("ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau1");
			pup1_tau2 = it->userFloat("ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau2");
			pup1_tau3 = it->userFloat("ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau3");
		      }
		    pup1_flav_parton = it->partonFlavour();
		    pup1_flav_hadron = it->hadronFlavour();
		    pup1_bdisc = subjetBdisc;
		  }

		if (subjetMass > mostMassiveSDPUPPIsubjetMass) mostMassiveSDPUPPIsubjetMass = subjetMass;
		count_pup++;
	      }
	  }
	if (useToolbox_)
	  {
	    for (const pat::Jet &isub : *AK8PUPPIsub)
	      {  
  
		double subjetPt        = isub.correctedP4(0).pt();
		double subjetEta       = isub.correctedP4(0).eta();
		double subjetPhi       = isub.correctedP4(0).phi();
		double subjetMass      = isub.correctedP4(0).mass();
		double subjetBdisc     = isub.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"); 
		double deltaRsubjetJet = deltaR(ijet.eta(), ijet.phi(), subjetEta, subjetPhi);

		if (deltaRsubjetJet < 0.8)
		  {
		    if (verbose_) cout << " matched puppi subjet with mass " << subjetMass << endl;

		    //------------------------------------
		    // PUPPI subjet JEC 
		    //------------------------------------
		    reco::Candidate::LorentzVector uncorrSubjet = isub.correctedP4(0);
		    JetCorrectorAK4pup->setJetEta(uncorrSubjet.eta());
		    JetCorrectorAK4pup->setJetPt (uncorrSubjet.pt());
		    JetCorrectorAK4pup->setJetE  (uncorrSubjet.energy());
		    JetCorrectorAK4pup->setJetA  (isub.jetArea());
		    JetCorrectorAK4pup->setRho   (rho);
		    JetCorrectorAK4pup->setNPV   (nvtx);
		    double subjet_corr_factor_L23res_full = JetCorrectorAK4pup->getCorrection();
		    reco::Candidate::LorentzVector corrSubjetL23res = subjet_corr_factor_L23res_full * uncorrSubjet;

		    //------------------------------------
		    // PUPPI subjet JEC uncertainty
		    //------------------------------------
		    double subjet_corrDn_L23 = 1.0;
		    JetCorrUncertAK4pup->setJetPhi(corrSubjetL23res.phi());
		    JetCorrUncertAK4pup->setJetEta(corrSubjetL23res.eta());
		    JetCorrUncertAK4pup->setJetPt(corrSubjetL23res.pt());
		    subjet_corrDn_L23   = subjet_corr_factor_L23res_full - JetCorrUncertAK4pup->getUncertainty(0);
		    double subjet_corrUp_L23 = 1.0;
		    JetCorrUncertAK4pup->setJetPhi(corrSubjetL23res.phi());
		    JetCorrUncertAK4pup->setJetEta(corrSubjetL23res.eta());
		    JetCorrUncertAK4pup->setJetPt(corrSubjetL23res.pt());
		    subjet_corrUp_L23 = subjet_corr_factor_L23res_full + JetCorrUncertAK4pup->getUncertainty(1);

		    reco::Candidate::LorentzVector corrSubjetL23resCorrDn   = subjet_corrDn_L23   * uncorrSubjet;
		    reco::Candidate::LorentzVector corrSubjetL23resCorrUp   = subjet_corrUp_L23   * uncorrSubjet;
       
		    //------------------------------------
		    // subjet JER SF
		    //------------------------------------
		    TLorentzVector GenSubJet;
		    double ptsmear   = 1.0;
		    double ptsmearUp = 1.0;
		    double ptsmearDn = 1.0;
		    if (!iEvent.isRealData())
		      {
			const reco::GenJet* genJet = isub.genJet();
			if (genJet)
			  {
			    GenSubJet.SetPtEtaPhiM(genJet->pt(), genJet->eta(), genJet->phi(), genJet->mass());
			    if (verbose_) cout << "  SD subjet genJet pt " << genJet->pt() << " mass " << genJet->mass() << " reco pt " << subjetPt << " reco mass " << subjetMass << endl;
			  }
			double jer_sf    = jer_scaler.getScaleFactor({{JME::Binning::JetEta, corrSubjetL23res.eta()}});
			double jer_sf_up = jer_scaler.getScaleFactor({{JME::Binning::JetEta, corrSubjetL23res.eta()}}, Variation::UP);
			double jer_sf_dn = jer_scaler.getScaleFactor({{JME::Binning::JetEta, corrSubjetL23res.eta()}}, Variation::DOWN);
			if (verbose_) std::cout << " SD subjet JER Scale factors (Nominal / Up / Down) : " << jer_sf << " / " << jer_sf_up << " / " << jer_sf_dn << std::endl;
			double recopt    = corrSubjetL23res.pt();
			double genpt     = GenJetMatched.Perp();
			double deltapt   = (recopt-genpt)*(jer_sf - 1.0);
			double deltaptUp = (recopt-genpt)*(jer_sf_up - 1.0);
			double deltaptDn = (recopt-genpt)*(jer_sf_dn - 1.0);
			ptsmear   = std::max((double)0.0, (recopt+deltapt)/recopt);
			ptsmearUp = std::max((double)0.0, (recopt+deltaptUp)/recopt);
			ptsmearDn = std::max((double)0.0, (recopt+deltaptDn)/recopt);
			if (verbose_) std::cout << " SD subjet ptsmear " << ptsmear << " ptsmearUp " << ptsmearUp << " ptsmearDn " << ptsmearDn << endl;
		      }
		    reco::Candidate::LorentzVector corrSubjetL23resPtSmear   = ptsmear   * corrSubjetL23res;
		    reco::Candidate::LorentzVector corrSubjetL23resPtSmearUp = ptsmearUp * corrSubjetL23res;
		    reco::Candidate::LorentzVector corrSubjetL23resPtSmearDn = ptsmearDn * corrSubjetL23res;

		    //------------------------------------
		    // subjet values for Tree
		    //------------------------------------

		    if (count_pup == 0)
		      {
			pup0_P4_uncorr            .SetPtEtaPhiM(subjetPt, subjetEta, subjetPhi, subjetMass);
			pup0_P4_L23res            .SetPtEtaPhiM(corrSubjetL23res          .pt() , corrSubjetL23res          .eta() , corrSubjetL23res          .phi() , corrSubjetL23res          .mass());
			pup0_P4_L23resCorrUp      .SetPtEtaPhiM(corrSubjetL23resCorrUp    .pt() , corrSubjetL23resCorrUp    .eta() , corrSubjetL23resCorrUp    .phi() , corrSubjetL23resCorrUp    .mass());
			pup0_P4_L23resCorrDn      .SetPtEtaPhiM(corrSubjetL23resCorrDn    .pt() , corrSubjetL23resCorrDn    .eta() , corrSubjetL23resCorrDn    .phi() , corrSubjetL23resCorrDn    .mass());
			pup0_P4_L23resPtSmear     .SetPtEtaPhiM(corrSubjetL23resPtSmear   .pt() , corrSubjetL23resPtSmear   .eta() , corrSubjetL23resPtSmear   .phi() , corrSubjetL23resPtSmear   .mass());
			pup0_P4_L23resPtSmearUp   .SetPtEtaPhiM(corrSubjetL23resPtSmearUp .pt() , corrSubjetL23resPtSmearUp .eta() , corrSubjetL23resPtSmearUp .phi() , corrSubjetL23resPtSmearUp .mass());
			pup0_P4_L23resPtSmearDn   .SetPtEtaPhiM(corrSubjetL23resPtSmearDn .pt() , corrSubjetL23resPtSmearDn .eta() , corrSubjetL23resPtSmearDn .phi() , corrSubjetL23resPtSmearDn .mass());

			pup0_tau1   = isub.userFloat("NsubjettinessAK8PFPuppiSoftDropSubjets:tau1");
			pup0_tau2   = isub.userFloat("NsubjettinessAK8PFPuppiSoftDropSubjets:tau2");
			pup0_tau3   = isub.userFloat("NsubjettinessAK8PFPuppiSoftDropSubjets:tau3");
          
			pup0_flav_parton   = isub.partonFlavour();
			pup0_flav_hadron   = isub.hadronFlavour();
			pup0_area          = isub.jetArea();
			pup0_bdisc         = subjetBdisc;
		      }
		    if (count_pup == 1)
		      {
			pup1_P4_uncorr            .SetPtEtaPhiM(subjetPt, subjetEta, subjetPhi, subjetMass);
			pup1_P4_L23res            .SetPtEtaPhiM(corrSubjetL23res          .pt() , corrSubjetL23res          .eta() , corrSubjetL23res          .phi() , corrSubjetL23res          .mass());
			pup1_P4_L23resCorrUp      .SetPtEtaPhiM(corrSubjetL23resCorrUp    .pt() , corrSubjetL23resCorrUp    .eta() , corrSubjetL23resCorrUp    .phi() , corrSubjetL23resCorrUp    .mass());
			pup1_P4_L23resCorrDn      .SetPtEtaPhiM(corrSubjetL23resCorrDn    .pt() , corrSubjetL23resCorrDn    .eta() , corrSubjetL23resCorrDn    .phi() , corrSubjetL23resCorrDn    .mass());
			pup1_P4_L23resPtSmear     .SetPtEtaPhiM(corrSubjetL23resPtSmear   .pt() , corrSubjetL23resPtSmear   .eta() , corrSubjetL23resPtSmear   .phi() , corrSubjetL23resPtSmear   .mass());
			pup1_P4_L23resPtSmearUp   .SetPtEtaPhiM(corrSubjetL23resPtSmearUp .pt() , corrSubjetL23resPtSmearUp .eta() , corrSubjetL23resPtSmearUp .phi() , corrSubjetL23resPtSmearUp .mass());
			pup1_P4_L23resPtSmearDn   .SetPtEtaPhiM(corrSubjetL23resPtSmearDn .pt() , corrSubjetL23resPtSmearDn .eta() , corrSubjetL23resPtSmearDn .phi() , corrSubjetL23resPtSmearDn .mass());

			pup1_tau1   = isub.userFloat("NsubjettinessAK8PFPuppiSoftDropSubjets:tau1");
			pup1_tau2   = isub.userFloat("NsubjettinessAK8PFPuppiSoftDropSubjets:tau2");
			pup1_tau3   = isub.userFloat("NsubjettinessAK8PFPuppiSoftDropSubjets:tau3");
          
			pup1_flav_parton   = isub.partonFlavour();
			pup1_flav_hadron   = isub.hadronFlavour();
			pup1_area          = isub.jetArea();
			pup1_bdisc         = subjetBdisc;
		      }

		    if (subjetMass > mostMassiveSDPUPPIsubjetMass) mostMassiveSDPUPPIsubjetMass = subjetMass;
		    count_pup++;
		    if (verbose_) cout << "  SD Subjet pt " << subjetPt << " Eta " << subjetEta << " deltaRsubjetJet " << deltaRsubjetJet << " Mass " << subjetMass << " Bdisc " << subjetBdisc << endl; 


		  }
	      }
	  }

	TLorentzVector sumPUPsubjets_P4_uncorr;
	TLorentzVector sumPUPsubjets_P4_L23res;
	TLorentzVector sumPUPsubjets_P4_L23resCorrUp;
	TLorentzVector sumPUPsubjets_P4_L23resCorrDn;
	TLorentzVector sumPUPsubjets_P4_L23resPtSmear;
	TLorentzVector sumPUPsubjets_P4_L23resPtSmearUp;
	TLorentzVector sumPUPsubjets_P4_L23resPtSmearDn;
	if (count_SD > 1)
	  { 
	    sumPUPsubjets_P4_uncorr            = pup0_P4_uncorr            + pup1_P4_uncorr; 
	    sumPUPsubjets_P4_L23res            = pup0_P4_L23res            + pup1_P4_L23res; 
	    sumPUPsubjets_P4_L23resCorrUp      = pup0_P4_L23resCorrUp      + pup1_P4_L23resCorrUp; 
	    sumPUPsubjets_P4_L23resCorrDn      = pup0_P4_L23resCorrDn      + pup1_P4_L23resCorrDn; 
	    sumPUPsubjets_P4_L23resPtSmear     = pup0_P4_L23resPtSmear     + pup1_P4_L23resPtSmear; 
	    sumPUPsubjets_P4_L23resPtSmearUp   = pup0_P4_L23resPtSmearUp   + pup1_P4_L23resPtSmearUp; 
	    sumPUPsubjets_P4_L23resPtSmearDn   = pup0_P4_L23resPtSmearDn   + pup1_P4_L23resPtSmearDn; 
	  } 


	double pup_maxbdisc = 0.0;
	double pup_maxbdiscflav_hadron = 0.0;
	double pup_maxbdiscflav_parton = 0.0;
	if (pup0_bdisc>=pup1_bdisc)
	  {
	    pup_maxbdisc = pup0_bdisc;
	    pup_maxbdiscflav_hadron = pup0_flav_hadron;
	    pup_maxbdiscflav_parton = pup0_flav_parton;
	  } 
	else if (pup1_bdisc>pup0_bdisc)
	  {
	    pup_maxbdisc = pup1_bdisc;
	    pup_maxbdiscflav_hadron = pup1_flav_hadron;
	    pup_maxbdiscflav_parton = pup1_flav_parton;
	  }  

	h_ak8chs_softDropMass  ->Fill(sumSDsubjets_P4_uncorr.M());
	h_ak8puppi_softDropMass->Fill(sumPUPsubjets_P4_uncorr.M());

	//------------------------------------
	// Gen particle info
	//------------------------------------ 
    
	double deltaR_jet_t1 = jet_p4.DeltaR(t1_p4);
	double deltaR_jet_t2 = jet_p4.DeltaR(t2_p4);
	bool jet_matched_t1 = false;
	bool jet_matched_t2 = false;
	if (deltaR_jet_t1<deltaR_jet_t2) jet_matched_t1 = true;
	if (deltaR_jet_t2<deltaR_jet_t1) jet_matched_t2 = true;

	double deltaR_jet_p1 =  jet_p4.DeltaR(hardest_parton_hardScatterOutgoing_p4);
	double deltaR_jet_p2 =  jet_p4.DeltaR(second_hardest_parton_hardScatterOutgoing_p4);
	bool jet_matched_p1 = false;
	bool jet_matched_p2 = false;
	if (deltaR_jet_p1<deltaR_jet_p2) jet_matched_p1 = true;
	if (deltaR_jet_p2<deltaR_jet_p1) jet_matched_p2 = true;
      
	//------------------------------------
	// Fill AllHadTree
	//------------------------------------ 
	if (verbose_) cout << "Fill AllHadTree " << endl;
	if (count_AK8CHS == 0)
	  {
	    Jet0PtRaw                              = uncorrJet.pt();                 
	    Jet0EtaRaw                             = uncorrJet.eta();                  
	    Jet0PhiRaw                             = uncorrJet.phi();   
	    Jet0MassRaw                            = uncorrJet.mass();                                           
	    Jet0P                                  = corrJet.P();        
	    Jet0Pt                                 = corrJet.pt();                  
	    Jet0Eta                                = corrJet.eta();                  
	    Jet0Phi                                = corrJet.phi();                  
	    Jet0Rap                                = corrJet.Rapidity();                  
	    Jet0Energy                             = corrJet.energy();                  
	    Jet0Mass                               = corrJet.mass();                    
	    Jet0Area                               = ijet.jetArea();                  
	    Jet0SDmass                             = softDropMass;                            
	    Jet0SDmassRaw                          = sumSDsubjets_P4_uncorr          .M(); // Full softDrop uncorrected 4-vector     
	    Jet0SDmassCorrL23                      = sumSDsubjets_P4_L23res          .M();   
	    Jet0SDmassCorrL23Up                    = sumSDsubjets_P4_L23resCorrUp    .M();   
	    Jet0SDmassCorrL23Dn                    = sumSDsubjets_P4_L23resCorrDn    .M(); 
	    Jet0SDmassCorrL123                     = sumSDsubjets_P4_L123res         .M();  
	    Jet0SDmassCorrL123Up                   = sumSDsubjets_P4_L123resCorrUp   .M();   
	    Jet0SDmassCorrL123Dn                   = sumSDsubjets_P4_L123resCorrDn   .M();  
	    Jet0SDmassCorrL23Smear                 = sumSDsubjets_P4_L23resPtSmear   .M();
	    Jet0SDmassCorrL23SmearUp               = sumSDsubjets_P4_L23resPtSmearUp .M();
	    Jet0SDmassCorrL23SmearDn               = sumSDsubjets_P4_L23resPtSmearDn .M();
	    Jet0SDptRaw                            = sumSDsubjets_P4_uncorr          .Perp();  // Full softDrop uncorrected 4-vector 
	    Jet0SDptCorrL23                        = sumSDsubjets_P4_L23res          .Perp();  
	    Jet0SDptCorrL23Up                      = sumSDsubjets_P4_L23resCorrUp    .Perp();  
	    Jet0SDptCorrL23Dn                      = sumSDsubjets_P4_L23resCorrDn    .Perp();  
	    Jet0SDptCorrL123                       = sumSDsubjets_P4_L123res         .Perp();  
	    Jet0SDptCorrL123Up                     = sumSDsubjets_P4_L123resCorrUp   .Perp();  
	    Jet0SDptCorrL123Dn                     = sumSDsubjets_P4_L123resCorrDn   .Perp();  
	    Jet0SDptCorrL23Smear                   = sumSDsubjets_P4_L23resPtSmear   .Perp();
	    Jet0SDptCorrL23SmearUp                 = sumSDsubjets_P4_L23resPtSmearUp .Perp();
	    Jet0SDptCorrL23SmearDn                 = sumSDsubjets_P4_L23resPtSmearDn .Perp();
	    Jet0SDetaRaw                           = sumSDsubjets_P4_uncorr.Eta();          // Full softDrop uncorrected 4-vector           
	    Jet0SDphiRaw                           = sumSDsubjets_P4_uncorr.Phi();          // Full softDrop uncorrected 4-vector           
	    Jet0MassPruned                         = prunedMass;     
	    Jet0MassTrimmed                        = trimmedMass;     
	    Jet0Tau1                               = tau1;  
	    Jet0Tau2                               = tau2;  
	    Jet0Tau3                               = tau3;  
	    Jet0Tau4                               = tau4;  
	    Jet0Tau32                              = tau32;  
	    Jet0Tau21                              = tau21;  
	    Jet0SDsubjet0bdisc                     = sub0_bdisc;  
	    Jet0SDsubjet1bdisc                     = sub1_bdisc;   
	    Jet0SDmaxbdisc                         = maxbdisc;
	    Jet0SDmaxbdiscflavHadron               = maxbdiscflav_hadron;  
	    Jet0SDmaxbdiscflavParton               = maxbdiscflav_parton;  
	    Jet0SDsubjet0pt                        = sub0_P4_uncorr.Pt();               
	    Jet0SDsubjet0mass                      = sub0_P4_uncorr.M();  
	    Jet0SDsubjet0eta                       = sub0_P4_uncorr.Eta();  
	    Jet0SDsubjet0phi                       = sub0_P4_uncorr.Phi();  
	    Jet0SDsubjet0area                      = sub0_area;  
	    Jet0SDsubjet0flavHadron                = sub0_flav_hadron;  
	    Jet0SDsubjet0flavParton                = sub0_flav_parton;  
	    Jet0SDsubjet0tau1                      = sub0_tau1;  
	    Jet0SDsubjet0tau2                      = sub0_tau2;  
	    Jet0SDsubjet0tau3                      = sub0_tau3;  
	    Jet0SDsubjet1pt                        = sub1_P4_uncorr.Pt();                    
	    Jet0SDsubjet1mass                      = sub1_P4_uncorr.M(); 
	    Jet0SDsubjet1eta                       = sub1_P4_uncorr.Eta();  
	    Jet0SDsubjet1phi                       = sub1_P4_uncorr.Phi();                     
	    Jet0SDsubjet1area                      = sub1_area;                    
	    Jet0SDsubjet1flavHadron                = sub1_flav_hadron;     
	    Jet0SDsubjet1flavParton                = sub1_flav_parton;     
	    Jet0SDsubjet1tau1                      = sub1_tau1;  
	    Jet0SDsubjet1tau2                      = sub1_tau2;  
	    Jet0SDsubjet1tau3                      = sub1_tau3; 

	    Jet0PuppiP                             = puppi_p;                  
	    Jet0PuppiPt                            = puppi_pt;                  
	    Jet0PuppiEta                           = puppi_eta;                   
	    Jet0PuppiPhi                           = puppi_phi;                  
	    Jet0PuppiMass                          = puppi_mass;                  
	    Jet0PuppiSDmass                        = sumPUPsubjets_P4_uncorr           .M();
	    Jet0PuppiSDmassCorr                    = sumPUPsubjets_P4_L23res           .M();
	    Jet0PuppiSDmassCorrUp                  = sumPUPsubjets_P4_L23resCorrUp     .M();
	    Jet0PuppiSDmassCorrDn                  = sumPUPsubjets_P4_L23resCorrDn     .M();
	    Jet0PuppiSDmassCorrL23Smear            = sumPUPsubjets_P4_L23resPtSmear    .M();
	    Jet0PuppiSDmassCorrL23SmearUp          = sumPUPsubjets_P4_L23resPtSmearUp  .M();
	    Jet0PuppiSDmassCorrL23SmearDn          = sumPUPsubjets_P4_L23resPtSmearDn  .M();
	    Jet0PuppiSDpt                          = sumPUPsubjets_P4_uncorr           .Perp();
	    Jet0PuppiSDptCorr                      = sumPUPsubjets_P4_L23res           .Perp();
	    Jet0PuppiSDptCorrUp                    = sumPUPsubjets_P4_L23resCorrUp     .Perp();
	    Jet0PuppiSDptCorrDn                    = sumPUPsubjets_P4_L23resCorrDn     .Perp();
	    Jet0PuppiSDptCorrL23Smear              = sumPUPsubjets_P4_L23resPtSmear    .Perp();
	    Jet0PuppiSDptCorrL23SmearUp            = sumPUPsubjets_P4_L23resPtSmearUp  .Perp();
	    Jet0PuppiSDptCorrL23SmearDn            = sumPUPsubjets_P4_L23resPtSmearDn  .Perp();
	    Jet0PuppiSDeta                         = sumPUPsubjets_P4_uncorr           .Eta();
	    Jet0PuppiSDphi                         = sumPUPsubjets_P4_uncorr           .Phi();
	    Jet0PuppiTau1                          = puppi_tau1;                  
	    Jet0PuppiTau2                          = puppi_tau2;                  
	    Jet0PuppiTau3                          = puppi_tau3;                  
	    Jet0PuppiTau4                          = puppi_tau4;                  
	    Jet0PuppiTau32                         = puppi_tau32;                  
	    Jet0PuppiTau21                         = puppi_tau21;                  
	    Jet0PuppiSDsubjet0bdisc                = pup0_bdisc;
	    Jet0PuppiSDsubjet1bdisc                = pup1_bdisc;
	    Jet0PuppiSDmaxbdisc                    = pup_maxbdisc;
	    Jet0PuppiSDmaxbdiscflavHadron          = pup_maxbdiscflav_hadron;
	    Jet0PuppiSDmaxbdiscflavParton          = pup_maxbdiscflav_parton;
	    Jet0PuppiSDsubjet0pt                   = pup0_P4_uncorr.Pt();
	    Jet0PuppiSDsubjet0mass                 = pup0_P4_uncorr.M();
	    Jet0PuppiSDsubjet0eta                  = pup0_P4_uncorr.Eta();
	    Jet0PuppiSDsubjet0phi                  = pup0_P4_uncorr.Phi();
	    Jet0PuppiSDsubjet0area                 = pup0_area;
	    Jet0PuppiSDsubjet0flavHadron           = pup0_flav_hadron;
	    Jet0PuppiSDsubjet0flavParton           = pup0_flav_parton;
	    Jet0PuppiSDsubjet0tau1                 = pup0_tau1;  
	    Jet0PuppiSDsubjet0tau2                 = pup0_tau2;  
	    Jet0PuppiSDsubjet0tau3                 = pup0_tau3; 
	    Jet0PuppiSDsubjet1pt                   = pup1_P4_uncorr.Pt();                 
	    Jet0PuppiSDsubjet1mass                 = pup1_P4_uncorr.M(); 
	    Jet0PuppiSDsubjet1eta                  = pup1_P4_uncorr.Eta();
	    Jet0PuppiSDsubjet1phi                  = pup1_P4_uncorr.Phi();             
	    Jet0PuppiSDsubjet1area                 = pup1_area;              
	    Jet0PuppiSDsubjet1flavHadron           = pup1_flav_hadron;   
	    Jet0PuppiSDsubjet1flavParton           = pup1_flav_parton;   
	    Jet0PuppiSDsubjet1tau1                 = pup1_tau1;  
	    Jet0PuppiSDsubjet1tau2                 = pup1_tau2;  
	    Jet0PuppiSDsubjet1tau3                 = pup1_tau3; 

	    Jet0CHF                                = ijet.chargedHadronEnergy() / uncorrJet.E();                        
	    Jet0NHF                                = ijet.neutralHadronEnergy() / uncorrJet.E();                         
	    Jet0CM                                 = ijet.chargedMultiplicity();                         
	    Jet0NM                                 = ijet.neutralMultiplicity();                          
	    Jet0NEF                                = ijet.neutralEmEnergy() / uncorrJet.E();                            
	    Jet0CEF                                = ijet.chargedEmEnergy() / uncorrJet.E();                          
	    Jet0MF                                 = ijet.muonEnergy() / uncorrJet.E();                         
	    Jet0Mult                               = ijet.numberOfDaughters();   

	    Jet0PuppiCHF                                = puppi_CHF; 
	    Jet0PuppiNHF                                = puppi_NHF; 
	    Jet0PuppiCM                                 = puppi_CM; 
	    Jet0PuppiNM                                 = puppi_NM; 
	    Jet0PuppiNEF                                = puppi_NEF; 
	    Jet0PuppiCEF                                = puppi_CEF; 
	    Jet0PuppiMF                                 = puppi_MF; 
	    Jet0PuppiMult                               = puppi_Mult; 

	    Jet0MassCorrFactor                     = corr_factor_L23res;        
	    Jet0MassCorrFactorUp                   = corrUp_L23;
	    Jet0MassCorrFactorDn                   = corrDn_L23;
	    Jet0CorrFactor                         = corr;        
	    Jet0CorrFactorUp                       = corrUp_L123;
	    Jet0CorrFactorDn                       = corrDn_L123;
	    Jet0PtSmearFactor                      = ptsmear;
	    Jet0PtSmearFactorUp                    = ptsmearUp;
	    Jet0PtSmearFactorDn                    = ptsmearDn;
	    Jet0PuppiMassCorrFactor                = 1;          
	    Jet0PuppiMassCorrFactorUp              = 1;          
	    Jet0PuppiMassCorrFactorDn              = 1;          
	    Jet0PuppiCorrFactor                    = 1;          
	    Jet0PuppiCorrFactorUp                  = 1;          
	    Jet0PuppiCorrFactorDn                  = 1;          
	    Jet0PuppiPtSmearFactor                 = 1;          
	    Jet0PuppiPtSmearFactorUp               = 1;          
	    Jet0PuppiPtSmearFactorDn               = 1;          
	    Jet0EtaScaleFactor                     = 1;          
	    Jet0PhiScaleFactor                     = 1;          
	    // Jet0MatchedGenJetDR                    = GenJetMatched_dRmin;             
	    Jet0MatchedGenJetPt                    = GenJetMatched.Perp();       
	    Jet0MatchedGenJetMass                  = GenJetMatched.M();   

	    if (!iEvent.isRealData() and runGenLoop_)
	      {
		if (counttop == 2 && jet_matched_t1)
		  {
		    if (top1hadronic) Jet0GenMatched_TopHadronic         = 1;
		    else Jet0GenMatched_TopHadronic                      = 0;
		    Jet0GenMatched_TopPt               = t1_p4.Perp();
		    Jet0GenMatched_TopEta              = t1_p4.Eta();
		    Jet0GenMatched_TopPhi              = t1_p4.Phi();
		    Jet0GenMatched_TopMass             = t1_p4.M();
		    Jet0GenMatched_bPt                 = b1_p4.Perp();
		    Jet0GenMatched_WPt                 = W1_p4.Perp();
		    Jet0GenMatched_Wd1Pt               = W1d1_p4.Perp();
		    Jet0GenMatched_Wd2Pt               = W1d2_p4.Perp();
		    Jet0GenMatched_Wd1ID               = W1d1_id;
		    Jet0GenMatched_Wd2ID               = W1d2_id;
		    Jet0GenMatched_MaxDeltaRPartonTop  = max_deltaR_parton_t1;
		    Jet0GenMatched_MaxDeltaRWPartonTop = max_deltaR_Wparton_t1;
		    Jet0GenMatched_MaxDeltaRWPartonW   = max_deltaR_Wparton_W1;
		    Jet0GenMatched_DeltaR_t_b          = deltaR_t1_b1;
		    Jet0GenMatched_DeltaR_t_W          = deltaR_t1_W1;
		    Jet0GenMatched_DeltaR_t_Wd1        = deltaR_t1_W1d1;
		    Jet0GenMatched_DeltaR_t_Wd2        = deltaR_t1_W1d2;
		    Jet0GenMatched_DeltaR_W_b1         = deltaR_W1_b1;
		    Jet0GenMatched_DeltaR_W_Wd1        = deltaR_W1_W1d1;
		    Jet0GenMatched_DeltaR_W_Wd2        = deltaR_W1_W1d2;
		    Jet0GenMatched_DeltaR_Wd1_Wd2      = deltaR_W1d1_W1d2;
		    Jet0GenMatched_DeltaR_Wd1_b        = deltaR_W1d1_b1;
		    Jet0GenMatched_DeltaR_Wd2_b        = deltaR_W1d2_b1;
		    Jet0GenMatched_DeltaR_jet_t        = deltaR_jet_t1;
		    Jet0GenMatched_DeltaR_jet_W        = jet_p4.DeltaR(W1_p4);
		    Jet0GenMatched_DeltaR_jet_b        = jet_p4.DeltaR(b1_p4);
		    Jet0GenMatched_DeltaR_jet_Wd1      = jet_p4.DeltaR(W1d1_p4);
		    Jet0GenMatched_DeltaR_jet_Wd2      = jet_p4.DeltaR(W1d2_p4);
		    Jet0GenMatched_DeltaR_pup0_b       = pup0_P4_L23res.DeltaR(b1_p4);
		    Jet0GenMatched_DeltaR_pup0_Wd1     = pup0_P4_L23res.DeltaR(W1d1_p4);
		    Jet0GenMatched_DeltaR_pup0_Wd2     = pup0_P4_L23res.DeltaR(W1d2_p4);
		    Jet0GenMatched_DeltaR_pup1_b       = pup1_P4_L23res.DeltaR(b1_p4);
		    Jet0GenMatched_DeltaR_pup1_Wd1     = pup1_P4_L23res.DeltaR(W1d1_p4);
		    Jet0GenMatched_DeltaR_pup1_Wd2     = pup1_P4_L23res.DeltaR(W1d2_p4);
		  }   
		if (counttop == 2 && jet_matched_t2)
		  {
		    if (top2hadronic) Jet0GenMatched_TopHadronic         = 1;
		    else Jet0GenMatched_TopHadronic                      = 0;
		    Jet0GenMatched_TopPt               = t2_p4.Perp();
		    Jet0GenMatched_TopEta              = t2_p4.Eta();
		    Jet0GenMatched_TopPhi              = t2_p4.Phi();
		    Jet0GenMatched_TopMass             = t2_p4.M();
		    Jet0GenMatched_bPt                 = b2_p4.Perp();
		    Jet0GenMatched_WPt                 = W2_p4.Perp();
		    Jet0GenMatched_Wd1Pt               = W2d1_p4.Perp();
		    Jet0GenMatched_Wd2Pt               = W2d2_p4.Perp();
		    Jet0GenMatched_Wd1ID               = W2d1_id;
		    Jet0GenMatched_Wd2ID               = W2d2_id;
		    Jet0GenMatched_MaxDeltaRPartonTop  = max_deltaR_parton_t2;
		    Jet0GenMatched_MaxDeltaRWPartonTop = max_deltaR_Wparton_t2;
		    Jet0GenMatched_MaxDeltaRWPartonW   = max_deltaR_Wparton_W2;
		    Jet0GenMatched_DeltaR_t_b          = deltaR_t2_b2;
		    Jet0GenMatched_DeltaR_t_W          = deltaR_t2_W2;
		    Jet0GenMatched_DeltaR_t_Wd1        = deltaR_t2_W2d1;
		    Jet0GenMatched_DeltaR_t_Wd2        = deltaR_t2_W2d2;
		    Jet0GenMatched_DeltaR_W_b1         = deltaR_W2_b2;
		    Jet0GenMatched_DeltaR_W_Wd1        = deltaR_W2_W2d1;
		    Jet0GenMatched_DeltaR_W_Wd2        = deltaR_W2_W2d2;
		    Jet0GenMatched_DeltaR_Wd1_Wd2      = deltaR_W2d1_W2d2;
		    Jet0GenMatched_DeltaR_Wd1_b        = deltaR_W2d1_b2;
		    Jet0GenMatched_DeltaR_Wd2_b        = deltaR_W2d2_b2;
		    Jet0GenMatched_DeltaR_jet_t        = deltaR_jet_t2;
		    Jet0GenMatched_DeltaR_jet_W        = jet_p4.DeltaR(W2_p4);
		    Jet0GenMatched_DeltaR_jet_b        = jet_p4.DeltaR(b2_p4);
		    Jet0GenMatched_DeltaR_jet_Wd1      = jet_p4.DeltaR(W2d1_p4);
		    Jet0GenMatched_DeltaR_jet_Wd2      = jet_p4.DeltaR(W2d2_p4);
		    Jet0GenMatched_DeltaR_pup0_b       = pup0_P4_L23res.DeltaR(b2_p4);
		    Jet0GenMatched_DeltaR_pup0_Wd1     = pup0_P4_L23res.DeltaR(W2d1_p4);
		    Jet0GenMatched_DeltaR_pup0_Wd2     = pup0_P4_L23res.DeltaR(W2d2_p4);
		    Jet0GenMatched_DeltaR_pup1_b       = pup1_P4_L23res.DeltaR(b2_p4);
		    Jet0GenMatched_DeltaR_pup1_Wd1     = pup1_P4_L23res.DeltaR(W2d1_p4);
		    Jet0GenMatched_DeltaR_pup1_Wd2     = pup1_P4_L23res.DeltaR(W2d2_p4);
		  }
		if (counttop == 0 && jet_matched_p1)
		  {
		    Jet0GenMatched_partonPt               = hardest_parton_hardScatterOutgoing_p4.Perp();
		    Jet0GenMatched_partonEta              = hardest_parton_hardScatterOutgoing_p4.Eta();
		    Jet0GenMatched_partonPhi              = hardest_parton_hardScatterOutgoing_p4.Phi();
		    Jet0GenMatched_partonMass             = hardest_parton_hardScatterOutgoing_p4.M();
		    Jet0GenMatched_partonID               = parton1id;
		    Jet0GenMatched_DeltaRjetParton        = deltaR_jet_p1;
		  }
		if (counttop == 0 && jet_matched_p2)
		  {
		    Jet0GenMatched_partonPt               = second_hardest_parton_hardScatterOutgoing_p4.Perp();
		    Jet0GenMatched_partonEta              = second_hardest_parton_hardScatterOutgoing_p4.Eta();
		    Jet0GenMatched_partonPhi              = second_hardest_parton_hardScatterOutgoing_p4.Phi();
		    Jet0GenMatched_partonMass             = second_hardest_parton_hardScatterOutgoing_p4.M();
		    Jet0GenMatched_partonID               = parton2id;
		    Jet0GenMatched_DeltaRjetParton        = deltaR_jet_p2;
		  }
	      }
	  }   
	if (count_AK8CHS == 1)
	  {
	    Jet1PtRaw                              = uncorrJet.pt();                 
	    Jet1EtaRaw                             = uncorrJet.eta();                  
	    Jet1PhiRaw                             = uncorrJet.phi();   
	    Jet1MassRaw                            = uncorrJet.mass();                                           
	    Jet1P                                  = corrJet.P();        
	    Jet1Pt                                 = corrJet.pt();                  
	    Jet1Eta                                = corrJet.eta();                  
	    Jet1Phi                                = corrJet.phi();                  
	    Jet1Rap                                = corrJet.Rapidity();                  
	    Jet1Energy                             = corrJet.energy();                  
	    Jet1Mass                               = corrJet.mass();                    
	    Jet1Area                               = ijet.jetArea();                  
	    Jet1SDmass                             = softDropMass;                           // Full softDrop uncorrected 4-vector      
	    Jet1SDmassRaw                          = sumSDsubjets_P4_uncorr          .M(); // Full softDrop uncorrected 4-vector     
	    Jet1SDmassCorrL23                      = sumSDsubjets_P4_L23res          .M();   
	    Jet1SDmassCorrL23Up                    = sumSDsubjets_P4_L23resCorrUp   .M();   
	    Jet1SDmassCorrL23Dn                    = sumSDsubjets_P4_L23resCorrDn   .M(); 
	    Jet1SDmassCorrL123                     = sumSDsubjets_P4_L123res         .M();  
	    Jet1SDmassCorrL123Up                   = sumSDsubjets_P4_L123resCorrUp  .M();   
	    Jet1SDmassCorrL123Dn                   = sumSDsubjets_P4_L123resCorrDn  .M();  
	    Jet1SDmassCorrL23Smear                    = sumSDsubjets_P4_L23resPtSmear   .M();
	    Jet1SDmassCorrL23SmearUp                  = sumSDsubjets_P4_L23resPtSmearUp .M();
	    Jet1SDmassCorrL23SmearDn                  = sumSDsubjets_P4_L23resPtSmearDn .M();
	    Jet1SDptRaw                            = sumSDsubjets_P4_uncorr          .Perp();  // Full softDrop uncorrected 4-vector 
	    Jet1SDptCorrL23                        = sumSDsubjets_P4_L23res          .Perp();  
	    Jet1SDptCorrL23Up                      = sumSDsubjets_P4_L23resCorrUp   .Perp();  
	    Jet1SDptCorrL23Dn                      = sumSDsubjets_P4_L23resCorrDn   .Perp();  
	    Jet1SDptCorrL123                       = sumSDsubjets_P4_L123res         .Perp();  
	    Jet1SDptCorrL123Up                     = sumSDsubjets_P4_L123resCorrUp  .Perp();  
	    Jet1SDptCorrL123Dn                     = sumSDsubjets_P4_L123resCorrDn  .Perp();  
	    Jet1SDptCorrL23Smear                      = sumSDsubjets_P4_L23resPtSmear   .Perp();
	    Jet1SDptCorrL23SmearUp                    = sumSDsubjets_P4_L23resPtSmearUp .Perp();
	    Jet1SDptCorrL23SmearDn                    = sumSDsubjets_P4_L23resPtSmearDn .Perp();
	    Jet1SDetaRaw                           = sumSDsubjets_P4_uncorr.Eta();  // Full softDrop uncorrected 4-vector           
	    Jet1SDphiRaw                           = sumSDsubjets_P4_uncorr.Phi();  // Full softDrop uncorrected 4-vector   
	    Jet1MassPruned                         = prunedMass;     
	    Jet1MassTrimmed                        = trimmedMass;     
	    Jet1Tau1                               = tau1;  
	    Jet1Tau2                               = tau2;  
	    Jet1Tau3                               = tau3;  
	    Jet1Tau4                               = tau4;  
	    Jet1Tau32                              = tau32;  
	    Jet1Tau21                              = tau21;  
	    Jet1SDsubjet0bdisc                     = sub0_bdisc;  
	    Jet1SDsubjet1bdisc                     = sub1_bdisc;   
	    Jet1SDmaxbdisc                         = maxbdisc;
	    Jet1SDmaxbdiscflavHadron               = maxbdiscflav_hadron;  
	    Jet1SDmaxbdiscflavParton               = maxbdiscflav_parton;  
	    Jet1SDsubjet0pt                        = sub0_P4_uncorr.Pt();               
	    Jet1SDsubjet0mass                      = sub0_P4_uncorr.M();  
	    Jet1SDsubjet0eta                       = sub0_P4_uncorr.Eta();  
	    Jet1SDsubjet0phi                       = sub0_P4_uncorr.Phi();  
	    Jet1SDsubjet0area                      = sub0_area;  
	    Jet1SDsubjet0flavHadron                = sub0_flav_hadron;  
	    Jet1SDsubjet0flavParton                = sub0_flav_parton;  
	    Jet1SDsubjet0tau1                      = sub0_tau1;  
	    Jet1SDsubjet0tau2                      = sub0_tau2;  
	    Jet1SDsubjet0tau3                      = sub0_tau3;  
	    Jet1SDsubjet1pt                        = sub1_P4_uncorr.Pt();                    
	    Jet1SDsubjet1mass                      = sub1_P4_uncorr.M();                    
	    Jet1SDsubjet1eta                       = sub1_P4_uncorr.Eta();                    
	    Jet1SDsubjet1phi                       = sub1_P4_uncorr.Phi();                    
	    Jet1SDsubjet1area                      = sub1_area;                    
	    Jet1SDsubjet1flavHadron                = sub1_flav_hadron;     
	    Jet1SDsubjet1flavParton                = sub1_flav_parton;     
	    Jet1SDsubjet1tau1                      = sub1_tau1;  
	    Jet1SDsubjet1tau2                      = sub1_tau2;  
	    Jet1SDsubjet1tau3                      = sub1_tau3;

	    Jet1PuppiP                             = puppi_p;                  
	    Jet1PuppiPt                            = puppi_pt;                  
	    Jet1PuppiEta                           = puppi_eta;                   
	    Jet1PuppiPhi                           = puppi_phi;                  
	    Jet1PuppiMass                          = puppi_mass;   

	    Jet1PuppiSDmass                        = sumPUPsubjets_P4_uncorr           .M();
	    Jet1PuppiSDmassCorr                    = sumPUPsubjets_P4_L23res           .M();
	    Jet1PuppiSDmassCorrUp                  = sumPUPsubjets_P4_L23resCorrUp     .M();
	    Jet1PuppiSDmassCorrDn                  = sumPUPsubjets_P4_L23resCorrDn     .M();
	    Jet1PuppiSDmassCorrL23Smear            = sumPUPsubjets_P4_L23resPtSmear    .M();
	    Jet1PuppiSDmassCorrL23SmearUp          = sumPUPsubjets_P4_L23resPtSmearUp  .M();
	    Jet1PuppiSDmassCorrL23SmearDn          = sumPUPsubjets_P4_L23resPtSmearDn  .M();
	    Jet1PuppiSDpt                          = sumPUPsubjets_P4_uncorr           .Perp();
	    Jet1PuppiSDptCorr                      = sumPUPsubjets_P4_L23res           .Perp();
	    Jet1PuppiSDptCorrUp                    = sumPUPsubjets_P4_L23resCorrUp     .Perp();
	    Jet1PuppiSDptCorrDn                    = sumPUPsubjets_P4_L23resCorrDn     .Perp();
	    Jet1PuppiSDptCorrL23Smear              = sumPUPsubjets_P4_L23resPtSmear    .Perp();
	    Jet1PuppiSDptCorrL23SmearUp            = sumPUPsubjets_P4_L23resPtSmearUp  .Perp();
	    Jet1PuppiSDptCorrL23SmearDn            = sumPUPsubjets_P4_L23resPtSmearDn  .Perp();
	    Jet1PuppiSDeta                         = sumPUPsubjets_P4_uncorr           .Eta();
	    Jet1PuppiSDphi                         = sumPUPsubjets_P4_uncorr           .Phi();

	    Jet1PuppiTau1                          = puppi_tau1;                  
	    Jet1PuppiTau2                          = puppi_tau2;                  
	    Jet1PuppiTau3                          = puppi_tau3;                  
	    Jet1PuppiTau4                          = puppi_tau4;                  
	    Jet1PuppiTau32                         = puppi_tau32;                  
	    Jet1PuppiTau21                         = puppi_tau21;                  
	    Jet1PuppiSDsubjet0bdisc                = pup0_bdisc;
	    Jet1PuppiSDsubjet1bdisc                = pup1_bdisc;
	    Jet1PuppiSDmaxbdisc                    = pup_maxbdisc;
	    Jet1PuppiSDmaxbdiscflavHadron          = pup_maxbdiscflav_hadron;
	    Jet1PuppiSDmaxbdiscflavParton          = pup_maxbdiscflav_parton;
	    Jet1PuppiSDsubjet0pt                   = pup0_P4_uncorr.Pt();
	    Jet1PuppiSDsubjet0mass                 = pup0_P4_uncorr.M();
	    Jet1PuppiSDsubjet0eta                  = pup0_P4_uncorr.Eta();
	    Jet1PuppiSDsubjet0phi                  = pup0_P4_uncorr.Phi();
	    Jet1PuppiSDsubjet0area                 = pup0_area;
	    Jet1PuppiSDsubjet0flavHadron           = pup0_flav_hadron;
	    Jet1PuppiSDsubjet0flavParton           = pup0_flav_parton;
	    Jet1PuppiSDsubjet0tau1                 = pup0_tau1;
	    Jet1PuppiSDsubjet0tau2                 = pup0_tau2;
	    Jet1PuppiSDsubjet0tau3                 = pup0_tau3;
	    Jet1PuppiSDsubjet1pt                   = pup1_P4_uncorr.Pt();                 
	    Jet1PuppiSDsubjet1mass                 = pup1_P4_uncorr.M();              
	    Jet1PuppiSDsubjet1eta                  = pup1_P4_uncorr.Eta();              
	    Jet1PuppiSDsubjet1phi                  = pup1_P4_uncorr.Phi();              
	    Jet1PuppiSDsubjet1area                 = pup1_area;              
	    Jet1PuppiSDsubjet1flavHadron           = pup1_flav_hadron;   
	    Jet1PuppiSDsubjet1flavParton           = pup1_flav_parton;   
	    Jet1PuppiSDsubjet1tau1                 = pup1_tau1;
	    Jet1PuppiSDsubjet1tau2                 = pup1_tau2;
	    Jet1PuppiSDsubjet1tau3                 = pup1_tau3;

	    Jet1CHF                                = ijet.chargedHadronEnergy() / uncorrJet.E();                        
	    Jet1NHF                                = ijet.neutralHadronEnergy() / uncorrJet.E();                         
	    Jet1CM                                 = ijet.chargedMultiplicity();                         
	    Jet1NM                                 = ijet.neutralMultiplicity();                          
	    Jet1NEF                                = ijet.neutralEmEnergy() / uncorrJet.E();                            
	    Jet1CEF                                = ijet.chargedEmEnergy() / uncorrJet.E();                          
	    Jet1MF                                 = ijet.muonEnergy() / uncorrJet.E();                         
	    Jet1Mult                               = ijet.numberOfDaughters();   


	    Jet1PuppiCHF                                = puppi_CHF; 
	    Jet1PuppiNHF                                = puppi_NHF; 
	    Jet1PuppiCM                                 = puppi_CM; 
	    Jet1PuppiNM                                 = puppi_NM; 
	    Jet1PuppiNEF                                = puppi_NEF; 
	    Jet1PuppiCEF                                = puppi_CEF; 
	    Jet1PuppiMF                                 = puppi_MF; 
	    Jet1PuppiMult                               = puppi_Mult; 


	    Jet1MassCorrFactor                     = corr_factor_L23res;        
	    Jet1MassCorrFactorUp                   = corrUp_L23;
	    Jet1MassCorrFactorDn                   = corrDn_L23;
	    Jet1CorrFactor                         = corr;        
	    Jet1CorrFactorUp                       = corrUp_L123;
	    Jet1CorrFactorDn                       = corrDn_L123;
	    Jet1PtSmearFactor                      = ptsmear;
	    Jet1PtSmearFactorUp                    = ptsmearUp;
	    Jet1PtSmearFactorDn                    = ptsmearDn;
	    Jet1PuppiMassCorrFactor                = 1;          
	    Jet1PuppiMassCorrFactorUp              = 1;          
	    Jet1PuppiMassCorrFactorDn              = 1;          
	    Jet1PuppiCorrFactor                    = 1;          
	    Jet1PuppiCorrFactorUp                  = 1;          
	    Jet1PuppiCorrFactorDn                  = 1;          
	    Jet1PuppiPtSmearFactor                 = 1;          
	    Jet1PuppiPtSmearFactorUp               = 1;          
	    Jet1PuppiPtSmearFactorDn               = 1;          
	    Jet1EtaScaleFactor                     = 1;          
	    Jet1PhiScaleFactor                     = 1;          
	    // Jet1MatchedGenJetDR                 = GenJetMatched_dRmin;             
	    Jet1MatchedGenJetPt                    = GenJetMatched.Perp();       
	    Jet1MatchedGenJetMass                  = GenJetMatched.M();   

	    if (!iEvent.isRealData() and runGenLoop_)
	      {
		if (counttop == 2 && jet_matched_t1)
		  {
		    if (top1hadronic) Jet1GenMatched_TopHadronic         = 1;
		    else Jet1GenMatched_TopHadronic                      = 0;
		    Jet1GenMatched_TopPt               = t1_p4.Perp();
		    Jet1GenMatched_TopEta              = t1_p4.Eta();
		    Jet1GenMatched_TopPhi              = t1_p4.Phi();
		    Jet1GenMatched_TopMass             = t1_p4.M();
		    Jet1GenMatched_bPt                 = b1_p4.Perp();
		    Jet1GenMatched_WPt                 = W1_p4.Perp();
		    Jet1GenMatched_Wd1Pt               = W1d1_p4.Perp();
		    Jet1GenMatched_Wd2Pt               = W1d2_p4.Perp();
		    Jet1GenMatched_Wd1ID               = W1d1_id;
		    Jet1GenMatched_Wd2ID               = W1d2_id;
		    Jet1GenMatched_MaxDeltaRPartonTop  = max_deltaR_parton_t1;
		    Jet1GenMatched_MaxDeltaRWPartonTop = max_deltaR_Wparton_t1;
		    Jet1GenMatched_MaxDeltaRWPartonW   = max_deltaR_Wparton_W1;
		    Jet1GenMatched_DeltaR_t_b          = deltaR_t1_b1;
		    Jet1GenMatched_DeltaR_t_W          = deltaR_t1_W1;
		    Jet1GenMatched_DeltaR_t_Wd1        = deltaR_t1_W1d1;
		    Jet1GenMatched_DeltaR_t_Wd2        = deltaR_t1_W1d2;
		    Jet1GenMatched_DeltaR_W_b1         = deltaR_W1_b1;
		    Jet1GenMatched_DeltaR_W_Wd1        = deltaR_W1_W1d1;
		    Jet1GenMatched_DeltaR_W_Wd2        = deltaR_W1_W1d2;
		    Jet1GenMatched_DeltaR_Wd1_Wd2      = deltaR_W1d1_W1d2;
		    Jet1GenMatched_DeltaR_Wd1_b        = deltaR_W1d1_b1;
		    Jet1GenMatched_DeltaR_Wd2_b        = deltaR_W1d2_b1;
		    Jet1GenMatched_DeltaR_jet_t        = deltaR_jet_t1;
		    Jet1GenMatched_DeltaR_jet_W        = jet_p4.DeltaR(W1_p4);
		    Jet1GenMatched_DeltaR_jet_b        = jet_p4.DeltaR(b1_p4);
		    Jet1GenMatched_DeltaR_jet_Wd1      = jet_p4.DeltaR(W1d1_p4);
		    Jet1GenMatched_DeltaR_jet_Wd2      = jet_p4.DeltaR(W1d2_p4);
		    Jet1GenMatched_DeltaR_pup0_b       = pup0_P4_L23res.DeltaR(b1_p4);
		    Jet1GenMatched_DeltaR_pup0_Wd1     = pup0_P4_L23res.DeltaR(W1d1_p4);
		    Jet1GenMatched_DeltaR_pup0_Wd2     = pup0_P4_L23res.DeltaR(W1d2_p4);
		    Jet1GenMatched_DeltaR_pup1_b       = pup1_P4_L23res.DeltaR(b1_p4);
		    Jet1GenMatched_DeltaR_pup1_Wd1     = pup1_P4_L23res.DeltaR(W1d1_p4);
		    Jet1GenMatched_DeltaR_pup1_Wd2     = pup1_P4_L23res.DeltaR(W1d2_p4);
		  }   
		if (counttop == 2 && jet_matched_t2)
		  {
		    if (top2hadronic) Jet1GenMatched_TopHadronic         = 1;
		    else              Jet1GenMatched_TopHadronic         = 0;
		    Jet1GenMatched_TopPt               = t2_p4.Perp();
		    Jet1GenMatched_TopEta              = t2_p4.Eta();
		    Jet1GenMatched_TopPhi              = t2_p4.Phi();
		    Jet1GenMatched_TopMass             = t2_p4.M();
		    Jet1GenMatched_bPt                 = b2_p4.Perp();
		    Jet1GenMatched_WPt                 = W2_p4.Perp();
		    Jet1GenMatched_Wd1Pt               = W2d1_p4.Perp();
		    Jet1GenMatched_Wd2Pt               = W2d2_p4.Perp();
		    Jet1GenMatched_Wd1ID               = W2d1_id;
		    Jet1GenMatched_Wd2ID               = W2d2_id;
		    Jet1GenMatched_MaxDeltaRPartonTop  = max_deltaR_parton_t2;
		    Jet1GenMatched_MaxDeltaRWPartonTop = max_deltaR_Wparton_t2;
		    Jet1GenMatched_MaxDeltaRWPartonW   = max_deltaR_Wparton_W2;
		    Jet1GenMatched_DeltaR_t_b          = deltaR_t2_b2;
		    Jet1GenMatched_DeltaR_t_W          = deltaR_t2_W2;
		    Jet1GenMatched_DeltaR_t_Wd1        = deltaR_t2_W2d1;
		    Jet1GenMatched_DeltaR_t_Wd2        = deltaR_t2_W2d2;
		    Jet1GenMatched_DeltaR_W_b1         = deltaR_W2_b2;
		    Jet1GenMatched_DeltaR_W_Wd1        = deltaR_W2_W2d1;
		    Jet1GenMatched_DeltaR_W_Wd2        = deltaR_W2_W2d2;
		    Jet1GenMatched_DeltaR_Wd1_Wd2      = deltaR_W2d1_W2d2;
		    Jet1GenMatched_DeltaR_Wd1_b        = deltaR_W2d1_b2;
		    Jet1GenMatched_DeltaR_Wd2_b        = deltaR_W2d2_b2;
		    Jet1GenMatched_DeltaR_jet_t        = deltaR_jet_t2;
		    Jet1GenMatched_DeltaR_jet_W        = jet_p4.DeltaR(W2_p4);
		    Jet1GenMatched_DeltaR_jet_b        = jet_p4.DeltaR(b2_p4);
		    Jet1GenMatched_DeltaR_jet_Wd1      = jet_p4.DeltaR(W2d1_p4);
		    Jet1GenMatched_DeltaR_jet_Wd2      = jet_p4.DeltaR(W2d2_p4);
		    Jet1GenMatched_DeltaR_pup0_b       = pup0_P4_L23res.DeltaR(b2_p4);
		    Jet1GenMatched_DeltaR_pup0_Wd1     = pup0_P4_L23res.DeltaR(W2d1_p4);
		    Jet1GenMatched_DeltaR_pup0_Wd2     = pup0_P4_L23res.DeltaR(W2d2_p4);
		    Jet1GenMatched_DeltaR_pup1_b       = pup1_P4_L23res.DeltaR(b2_p4);
		    Jet1GenMatched_DeltaR_pup1_Wd1     = pup1_P4_L23res.DeltaR(W2d1_p4);
		    Jet1GenMatched_DeltaR_pup1_Wd2     = pup1_P4_L23res.DeltaR(W2d2_p4);
		  }
		if (counttop == 0 && jet_matched_p1)
		  {
		    Jet1GenMatched_partonPt               = hardest_parton_hardScatterOutgoing_p4.Perp();
		    Jet1GenMatched_partonEta              = hardest_parton_hardScatterOutgoing_p4.Eta();
		    Jet1GenMatched_partonPhi              = hardest_parton_hardScatterOutgoing_p4.Phi();
		    Jet1GenMatched_partonMass             = hardest_parton_hardScatterOutgoing_p4.M();
		    Jet1GenMatched_partonID               = parton1id;
		    Jet1GenMatched_DeltaRjetParton        = deltaR_jet_p1;
		  }
		if (counttop == 0 && jet_matched_p2) {
		  Jet1GenMatched_partonPt               = second_hardest_parton_hardScatterOutgoing_p4.Perp();
		  Jet1GenMatched_partonEta              = second_hardest_parton_hardScatterOutgoing_p4.Eta();
		  Jet1GenMatched_partonPhi              = second_hardest_parton_hardScatterOutgoing_p4.Phi();
		  Jet1GenMatched_partonMass             = second_hardest_parton_hardScatterOutgoing_p4.M();
		  Jet1GenMatched_partonID               = parton2id;
		  Jet1GenMatched_DeltaRjetParton        = deltaR_jet_p2;
		}
	      }   
	  } 


	//------------------------------------
	// Fill SemiLeptTree
	//------------------------------------ 
	if (verbose_) cout << "Fill SemiLeptTree " << endl;
	double deltaPhi_lep_jet = fabs(deltaPhi(corrJet.phi(), lep0_p4.Phi()));
	// AK8 jet should be in opposite hemisphere from lepton. If leading jet matches then use it. If it doensn't then check the second leading jet.
	if (((count_AK8CHS == 0&& deltaPhi_lep_jet >=3.14/2) || (count_AK8CHS == 1&&deltaPhi_lep_jet >=3.14/2)) && count_lep >=1 && count_fill_leptTree == 0)
	  {
	    count_fill_leptTree++;
	    DeltaRJetLep                          = deltaR(corrJet.eta(), corrJet.phi(), lep0_p4.Eta(), lep0_p4.Phi());
	    DeltaPhiJetLep                        = deltaPhi_lep_jet;
	    JetPtRaw                              = uncorrJet.pt();                 
	    JetEtaRaw                             = uncorrJet.eta();                  
	    JetPhiRaw                             = uncorrJet.phi();   
	    JetMassRaw                            = uncorrJet.mass();                                           
	    JetP                                  = corrJet.P();        
	    JetPt                                 = corrJet.pt();                  
	    JetEta                                = corrJet.eta();                  
	    JetPhi                                = corrJet.phi();                  
	    JetRap                                = corrJet.Rapidity();                  
	    JetEnergy                             = corrJet.energy();                  
	    JetMass                               = corrJet.mass();                    
	    JetArea                               = ijet.jetArea();                  
	    JetSDmass                             = softDropMass;                           // Full softDrop uncorrected 4-vector      
	    JetSDmassRaw                          = sumSDsubjets_P4_uncorr          .M(); // Full softDrop uncorrected 4-vector     
	    JetSDmassCorrL23                      = sumSDsubjets_P4_L23res          .M();   
	    JetSDmassCorrL23Up                    = sumSDsubjets_P4_L23resCorrUp   .M();   
	    JetSDmassCorrL23Dn                    = sumSDsubjets_P4_L23resCorrDn   .M(); 
	    JetSDmassCorrL123                     = sumSDsubjets_P4_L123res         .M();  
	    JetSDmassCorrL123Up                   = sumSDsubjets_P4_L123resCorrUp  .M();   
	    JetSDmassCorrL123Dn                   = sumSDsubjets_P4_L123resCorrDn  .M();  
	    JetSDmassCorrL23Smear                    = sumSDsubjets_P4_L23resPtSmear   .M();
	    JetSDmassCorrL23SmearUp                  = sumSDsubjets_P4_L23resPtSmearUp .M();
	    JetSDmassCorrL23SmearDn                  = sumSDsubjets_P4_L23resPtSmearDn .M();
	    JetSDptRaw                            = sumSDsubjets_P4_uncorr          .Perp();  // Full softDrop uncorrected 4-vector 
	    JetSDptCorrL23                        = sumSDsubjets_P4_L23res          .Perp();  
	    JetSDptCorrL23Up                      = sumSDsubjets_P4_L23resCorrUp   .Perp();  
	    JetSDptCorrL23Dn                      = sumSDsubjets_P4_L23resCorrDn   .Perp();  
	    JetSDptCorrL123                       = sumSDsubjets_P4_L123res         .Perp();  
	    JetSDptCorrL123Up                     = sumSDsubjets_P4_L123resCorrUp  .Perp();  
	    JetSDptCorrL123Dn                     = sumSDsubjets_P4_L123resCorrDn  .Perp();  
	    JetSDptCorrL23Smear                      = sumSDsubjets_P4_L23resPtSmear   .Perp();
	    JetSDptCorrL23SmearUp                    = sumSDsubjets_P4_L23resPtSmearUp .Perp();
	    JetSDptCorrL23SmearDn                    = sumSDsubjets_P4_L23resPtSmearDn .Perp();
	    JetSDetaRaw                           = sumSDsubjets_P4_uncorr.Eta();      // Full softDrop uncorrected 4-vector           
	    JetSDphiRaw                           = sumSDsubjets_P4_uncorr.Phi();      // Full softDrop uncorrected 4-vector   
	    JetMassPruned                         = prunedMass;     
	    JetMassTrimmed                        = trimmedMass;     
	    JetTau1                               = tau1;  
	    JetTau2                               = tau2;  
	    JetTau3                               = tau3;  
	    JetTau4                               = tau4;  
	    JetTau32                              = tau32;  
	    JetTau21                              = tau21;  
	    JetSDsubjet0bdisc                     = sub0_bdisc;  
	    JetSDsubjet1bdisc                     = sub1_bdisc;   
	    JetSDmaxbdisc                         = maxbdisc;
	    JetSDmaxbdiscflavHadron               = maxbdiscflav_hadron;  
	    JetSDmaxbdiscflavParton               = maxbdiscflav_parton;  
	    JetSDsubjet0pt                        = sub0_P4_uncorr.Pt();               
	    JetSDsubjet0mass                      = sub0_P4_uncorr.M();  
	    JetSDsubjet0eta                       = sub0_P4_uncorr.Eta();  
	    JetSDsubjet0phi                       = sub0_P4_uncorr.Phi();  
	    JetSDsubjet0area                      = sub0_area;  
	    JetSDsubjet0flavHadron                = sub0_flav_hadron;  
	    JetSDsubjet0flavParton                = sub0_flav_parton;  
	    JetSDsubjet0tau1                      = sub0_tau1;
	    JetSDsubjet0tau2                      = sub0_tau2;
	    JetSDsubjet0tau3                      = sub0_tau3;
	    JetSDsubjet1pt                        = sub1_P4_uncorr.Pt();                    
	    JetSDsubjet1mass                      = sub1_P4_uncorr.M();                    
	    JetSDsubjet1eta                       = sub1_P4_uncorr.Eta();                    
	    JetSDsubjet1phi                       = sub1_P4_uncorr.Phi();                    
	    JetSDsubjet1area                      = sub1_area;                    
	    JetSDsubjet1flavHadron                = sub1_flav_hadron;     
	    JetSDsubjet1flavParton                = sub1_flav_parton;     
	    JetSDsubjet1tau1                      = sub1_tau1;
	    JetSDsubjet1tau2                      = sub1_tau2;
	    JetSDsubjet1tau3                      = sub1_tau3;

	    AK8jet_SemiLept_P4corr.SetPtEtaPhiM(corrJet.pt(), corrJet.eta(), corrJet.phi(), corrJet.mass());

	    JetPuppiP                             = puppi_p;                  
	    JetPuppiPt                            = puppi_pt;                  
	    JetPuppiEta                           = puppi_eta;                   
	    JetPuppiPhi                           = puppi_phi;                  
	    JetPuppiMass                          = puppi_mass;                  

	    Jet1PuppiSDmass                       = sumPUPsubjets_P4_uncorr           .M();
	    Jet1PuppiSDmassCorr                   = sumPUPsubjets_P4_L23res           .M();
	    Jet1PuppiSDmassCorrUp                 = sumPUPsubjets_P4_L23resCorrUp     .M();
	    Jet1PuppiSDmassCorrDn                 = sumPUPsubjets_P4_L23resCorrDn     .M();
	    Jet1PuppiSDmassCorrL23Smear           = sumPUPsubjets_P4_L23resPtSmear    .M();
	    Jet1PuppiSDmassCorrL23SmearUp         = sumPUPsubjets_P4_L23resPtSmearUp  .M();
	    Jet1PuppiSDmassCorrL23SmearDn         = sumPUPsubjets_P4_L23resPtSmearDn  .M();
	    Jet1PuppiSDpt                         = sumPUPsubjets_P4_uncorr           .Perp();
	    Jet1PuppiSDptCorr                     = sumPUPsubjets_P4_L23res           .Perp();
	    Jet1PuppiSDptCorrUp                   = sumPUPsubjets_P4_L23resCorrUp     .Perp();
	    Jet1PuppiSDptCorrDn                   = sumPUPsubjets_P4_L23resCorrDn     .Perp();
	    Jet1PuppiSDptCorrL23Smear             = sumPUPsubjets_P4_L23resPtSmear    .Perp();
	    Jet1PuppiSDptCorrL23SmearUp           = sumPUPsubjets_P4_L23resPtSmearUp  .Perp();
	    Jet1PuppiSDptCorrL23SmearDn           = sumPUPsubjets_P4_L23resPtSmearDn  .Perp();
	    Jet1PuppiSDeta                        = sumPUPsubjets_P4_uncorr           .Eta();
	    Jet1PuppiSDphi                        = sumPUPsubjets_P4_uncorr           .Phi();

	    JetPuppiTau1                          = puppi_tau1;                  
	    JetPuppiTau2                          = puppi_tau2;                  
	    JetPuppiTau3                          = puppi_tau3;                  
	    JetPuppiTau4                          = puppi_tau4;                  
	    JetPuppiTau32                         = puppi_tau32;                  
	    JetPuppiTau21                         = puppi_tau21;                  
	    JetPuppiSDsubjet0bdisc                = pup0_bdisc;
	    JetPuppiSDsubjet1bdisc                = pup1_bdisc;
	    JetPuppiSDmaxbdisc                    = pup_maxbdisc;   
	    JetPuppiSDmaxbdiscflavHadron          = pup_maxbdiscflav_hadron;
	    JetPuppiSDmaxbdiscflavParton          = pup_maxbdiscflav_parton;
	    JetPuppiSDsubjet0pt                   = pup0_P4_uncorr.Pt();
	    JetPuppiSDsubjet0mass                 = pup0_P4_uncorr.M();
	    JetPuppiSDsubjet0eta                  = pup0_P4_uncorr.Eta();
	    JetPuppiSDsubjet0phi                  = pup0_P4_uncorr.Phi();
	    JetPuppiSDsubjet0area                 = pup0_area;
	    JetPuppiSDsubjet0flavHadron           = pup0_flav_hadron;
	    JetPuppiSDsubjet0flavParton           = pup0_flav_parton;
	    JetPuppiSDsubjet0tau1                 = pup0_tau1;
	    JetPuppiSDsubjet0tau2                 = pup0_tau2;
	    JetPuppiSDsubjet0tau3                 = pup0_tau3;
	    JetPuppiSDsubjet1pt                   = pup1_P4_uncorr.Pt();                 
	    JetPuppiSDsubjet1mass                 = pup1_P4_uncorr.M();              
	    JetPuppiSDsubjet1eta                  = pup1_P4_uncorr.Eta();              
	    JetPuppiSDsubjet1phi                  = pup1_P4_uncorr.Phi();              
	    JetPuppiSDsubjet1area                 = pup1_area;              
	    JetPuppiSDsubjet1flavHadron           = pup1_flav_hadron;   
	    JetPuppiSDsubjet1flavParton           = pup1_flav_parton;   
	    JetPuppiSDsubjet1tau1                 = pup1_tau1;
	    JetPuppiSDsubjet1tau2                 = pup1_tau2;
	    JetPuppiSDsubjet1tau3                 = pup1_tau3;

	    JetCHF                                = ijet.chargedHadronEnergy() / uncorrJet.E();                        
	    JetNHF                                = ijet.neutralHadronEnergy() / uncorrJet.E();                         
	    JetCM                                 = ijet.chargedMultiplicity();                         
	    JetNM                                 = ijet.neutralMultiplicity();                          
	    JetNEF                                = ijet.neutralEmEnergy() / uncorrJet.E();                            
	    JetCEF                                = ijet.chargedEmEnergy() / uncorrJet.E();                          
	    JetMF                                 = ijet.muonEnergy() / uncorrJet.E();                         
	    JetMult                               = ijet.numberOfDaughters();   


	    JetPuppiCHF                                = puppi_CHF; 
	    JetPuppiNHF                                = puppi_NHF; 
	    JetPuppiCM                                 = puppi_CM; 
	    JetPuppiNM                                 = puppi_NM; 
	    JetPuppiNEF                                = puppi_NEF; 
	    JetPuppiCEF                                = puppi_CEF; 
	    JetPuppiMF                                 = puppi_MF; 
	    JetPuppiMult                               = puppi_Mult; 

	    JetMassCorrFactor                     = corr_factor_L23res;        
	    JetMassCorrFactorUp                   = corrUp_L23;
	    JetMassCorrFactorDn                   = corrDn_L23;
	    JetCorrFactor                         = corr;        
	    JetCorrFactorUp                       = corrUp_L123;
	    JetCorrFactorDn                       = corrDn_L123;
	    JetPtSmearFactor                      = ptsmear;
	    JetPtSmearFactorUp                    = ptsmearUp;
	    JetPtSmearFactorDn                    = ptsmearDn;
	    JetPuppiMassCorrFactor                = 1;          
	    JetPuppiMassCorrFactorUp              = 1;          
	    JetPuppiMassCorrFactorDn              = 1;          
	    JetPuppiCorrFactor                    = 1;          
	    JetPuppiCorrFactorUp                  = 1;          
	    JetPuppiCorrFactorDn                  = 1;          
	    JetPuppiPtSmearFactor                 = 1;          
	    JetPuppiPtSmearFactorUp               = 1;          
	    JetPuppiPtSmearFactorDn               = 1;          
	    JetEtaScaleFactor                     = 1;          
	    JetPhiScaleFactor                     = 1;          
	    // JetMatchedGenJetDR                    = GenJetMatched_dRmin;             
	    JetMatchedGenJetPt                    = GenJetMatched.Perp();       
	    JetMatchedGenJetMass                  = GenJetMatched.M();   

	    if (!iEvent.isRealData() and runGenLoop_)
	      {
		if (counttop == 2 && jet_matched_t1)
		  {
		    if (top1hadronic) JetGenMatched_TopHadronic         = 1;
		    else              JetGenMatched_TopHadronic         = 0;
		    JetGenMatched_TopPt               = t1_p4.Perp();
		    JetGenMatched_TopEta              = t1_p4.Eta();
		    JetGenMatched_TopPhi              = t1_p4.Phi();
		    JetGenMatched_TopMass             = t1_p4.M();
		    JetGenMatched_bPt                 = b1_p4.Perp();
		    JetGenMatched_WPt                 = W1_p4.Perp();
		    JetGenMatched_Wd1Pt               = W1d1_p4.Perp();
		    JetGenMatched_Wd2Pt               = W1d2_p4.Perp();
		    JetGenMatched_Wd1ID               = W1d1_id;
		    JetGenMatched_Wd2ID               = W1d2_id;
		    JetGenMatched_MaxDeltaRPartonTop  = max_deltaR_parton_t1;
		    JetGenMatched_MaxDeltaRWPartonTop = max_deltaR_Wparton_t1;
		    JetGenMatched_MaxDeltaRWPartonW   = max_deltaR_Wparton_W1;
		    JetGenMatched_DeltaR_t_b          = deltaR_t1_b1;
		    JetGenMatched_DeltaR_t_W          = deltaR_t1_W1;
		    JetGenMatched_DeltaR_t_Wd1        = deltaR_t1_W1d1;
		    JetGenMatched_DeltaR_t_Wd2        = deltaR_t1_W1d2;
		    JetGenMatched_DeltaR_W_b1         = deltaR_W1_b1;
		    JetGenMatched_DeltaR_W_Wd1        = deltaR_W1_W1d1;
		    JetGenMatched_DeltaR_W_Wd2        = deltaR_W1_W1d2;
		    JetGenMatched_DeltaR_Wd1_Wd2      = deltaR_W1d1_W1d2;
		    JetGenMatched_DeltaR_Wd1_b        = deltaR_W1d1_b1;
		    JetGenMatched_DeltaR_Wd2_b        = deltaR_W1d2_b1;
		    JetGenMatched_DeltaR_jet_t        = deltaR_jet_t1;
		    JetGenMatched_DeltaR_jet_W        = jet_p4.DeltaR(W1_p4);
		    JetGenMatched_DeltaR_jet_b        = jet_p4.DeltaR(b1_p4);
		    JetGenMatched_DeltaR_jet_Wd1      = jet_p4.DeltaR(W1d1_p4);
		    JetGenMatched_DeltaR_jet_Wd2      = jet_p4.DeltaR(W1d2_p4);
		    JetGenMatched_DeltaR_pup0_b       = pup0_P4_L23res.DeltaR(b1_p4);
		    JetGenMatched_DeltaR_pup0_Wd1     = pup0_P4_L23res.DeltaR(W1d1_p4);
		    JetGenMatched_DeltaR_pup0_Wd2     = pup0_P4_L23res.DeltaR(W1d2_p4);
		    JetGenMatched_DeltaR_pup1_b       = pup1_P4_L23res.DeltaR(b1_p4);
		    JetGenMatched_DeltaR_pup1_Wd1     = pup1_P4_L23res.DeltaR(W1d1_p4);
		    JetGenMatched_DeltaR_pup1_Wd2     = pup1_P4_L23res.DeltaR(W1d2_p4);
		  }   
		if (counttop == 2 && jet_matched_t2)
		  {
		    if (top2hadronic) JetGenMatched_TopHadronic         = 1;
		    else             JetGenMatched_TopHadronic         = 0;
		    JetGenMatched_TopPt               = t2_p4.Perp();
		    JetGenMatched_TopEta              = t2_p4.Eta();
		    JetGenMatched_TopPhi              = t2_p4.Phi();
		    JetGenMatched_TopMass             = t2_p4.M();
		    JetGenMatched_bPt                 = b2_p4.Perp();
		    JetGenMatched_WPt                 = W2_p4.Perp();
		    JetGenMatched_Wd1Pt               = W2d1_p4.Perp();
		    JetGenMatched_Wd2Pt               = W2d2_p4.Perp();
		    JetGenMatched_Wd1ID               = W2d1_id;
		    JetGenMatched_Wd2ID               = W2d2_id;
		    JetGenMatched_MaxDeltaRPartonTop  = max_deltaR_parton_t2;
		    JetGenMatched_MaxDeltaRWPartonTop = max_deltaR_Wparton_t2;
		    JetGenMatched_MaxDeltaRWPartonW   = max_deltaR_Wparton_W2;
		    JetGenMatched_DeltaR_t_b          = deltaR_t2_b2;
		    JetGenMatched_DeltaR_t_W          = deltaR_t2_W2;
		    JetGenMatched_DeltaR_t_Wd1        = deltaR_t2_W2d1;
		    JetGenMatched_DeltaR_t_Wd2        = deltaR_t2_W2d2;
		    JetGenMatched_DeltaR_W_b1         = deltaR_W2_b2;
		    JetGenMatched_DeltaR_W_Wd1        = deltaR_W2_W2d1;
		    JetGenMatched_DeltaR_W_Wd2        = deltaR_W2_W2d2;
		    JetGenMatched_DeltaR_Wd1_Wd2      = deltaR_W2d1_W2d2;
		    JetGenMatched_DeltaR_Wd1_b        = deltaR_W2d1_b2;
		    JetGenMatched_DeltaR_Wd2_b        = deltaR_W2d2_b2;
		    JetGenMatched_DeltaR_jet_t        = deltaR_jet_t2;
		    JetGenMatched_DeltaR_jet_W        = jet_p4.DeltaR(W2_p4);
		    JetGenMatched_DeltaR_jet_b        = jet_p4.DeltaR(b2_p4);
		    JetGenMatched_DeltaR_jet_Wd1      = jet_p4.DeltaR(W2d1_p4);
		    JetGenMatched_DeltaR_jet_Wd2      = jet_p4.DeltaR(W2d2_p4);
		    JetGenMatched_DeltaR_pup0_b       = pup0_P4_L23res.DeltaR(b2_p4);
		    JetGenMatched_DeltaR_pup0_Wd1     = pup0_P4_L23res.DeltaR(W2d1_p4);
		    JetGenMatched_DeltaR_pup0_Wd2     = pup0_P4_L23res.DeltaR(W2d2_p4);
		    JetGenMatched_DeltaR_pup1_b       = pup1_P4_L23res.DeltaR(b2_p4);
		    JetGenMatched_DeltaR_pup1_Wd1     = pup1_P4_L23res.DeltaR(W2d1_p4);
		    JetGenMatched_DeltaR_pup1_Wd2     = pup1_P4_L23res.DeltaR(W2d2_p4);
		  }
		if (counttop == 0 && jet_matched_p1)
		  {
		    JetGenMatched_partonPt               = hardest_parton_hardScatterOutgoing_p4.Perp();
		    JetGenMatched_partonEta              = hardest_parton_hardScatterOutgoing_p4.Eta();
		    JetGenMatched_partonPhi              = hardest_parton_hardScatterOutgoing_p4.Phi();
		    JetGenMatched_partonMass             = hardest_parton_hardScatterOutgoing_p4.M();
		    JetGenMatched_partonID               = parton1id;
		    JetGenMatched_DeltaRjetParton        = deltaR_jet_p1;
		  }
		if (counttop == 0 && jet_matched_p2)
		  {
		    JetGenMatched_partonPt               = second_hardest_parton_hardScatterOutgoing_p4.Perp();
		    JetGenMatched_partonEta              = second_hardest_parton_hardScatterOutgoing_p4.Eta();
		    JetGenMatched_partonPhi              = second_hardest_parton_hardScatterOutgoing_p4.Phi();
		    JetGenMatched_partonMass             = second_hardest_parton_hardScatterOutgoing_p4.M();
		    JetGenMatched_partonID               = parton2id;
		    JetGenMatched_DeltaRjetParton        = deltaR_jet_p2;
		  }
	      }
	  } 

	count_AK8CHS++;
      } //matches beginning of for loop where crash occurs
#endif
    // MARKER RECO PARTICLES
    //  ____  _____ ____ ___     ____   _    ____ _____ ___ ____ _     _____ ____  
    // |  _ \| ____/ ___/ _ \   |  _ \ / \  |  _ \_   _|_ _/ ___| |   | ____/ ___| 
    // | |_) |  _|| |  | | | |  | |_) / _ \ | |_) || |  | | |   | |   |  _| \___ \ |
    // |  _ <| |__| |__| |_| |  |  __/ ___ \|  _ < | |  | | |___| |___| |___ ___) |
    // |_| \_\_____\____\___/   |_| /_/   \_\_| \_\|_| |___\____|_____|_____|____/ 
    //                                                                          


    edm::Handle<pat::ElectronCollection> reco_electrons;
    iEvent.getByToken(electronToken_, reco_electrons);

    edm::Handle<pat::MuonCollection> reco_muons;
    iEvent.getByToken(muonToken_, reco_muons);

    edm::Handle<pat::JetCollection> reco_ak4_jets;
    iEvent.getByToken(ak8jetToken_, reco_ak4_jets);

    edm::Handle<pat::JetCollection> reco_ak8_jets;// max = 4; mode = 2; mean = 1.7;
    iEvent.getByToken(ak8jetToken_, reco_ak8_jets);

    //edm::Handle<pat::JetCollection> reco_puppi_jets;
    //iEvent.getByToken(puppijetToken_, reco_puppi_jets);

    //edm::Handle<reco::GenJetCollection> gen_ak4_jets;  
    //iEvent.getByToken(ak8genjetToken_, gen_ak4_jets);

    edm::Handle<reco::GenJetCollection> gen_ak8_jets;// max = 5; mode = 2; mean = 2.0;
    iEvent.getByToken(ak8genjetToken_, gen_ak8_jets);

    edm::Handle<pat::JetCollection> sub_jets;// max = 7; mode = 3; mean = 2.5;
    iEvent.getByToken(ak8CHSSoftDropSubjetsToken_, sub_jets);

   
    //Float_t reco_ak4_jet_pt_max = 0.0;
    //Float_t reco_ak4_jet_pt_next_to_max = 0.0;

    //Float_t reco_ak8_jet_pt_max = 0.0;
    //Float_t reco_ak8_jet_pt_next_to_max = 0.0;

    //Float_t reco_puppi_jet_pt_max = 0.0;
    //Float_t reco_puppi_jet_pt_next_to_max = 0.0;

    //Float_t gen_ak4_jet_pt_max = 0.0;
    //Float_t gen_ak4_jet_pt_next_to_max = 0.0;

    //Float_t gen_ak8_jet_pt_max = 0.0;
    //Float_t gen_ak8_jet_pt_next_to_max = 0.0;
    //Float_t num2aftermax = 0.0;
    //Float_t num3aftermax = 0.0;
    //Float_t num4aftermax = 0.0;
    //Float_t num5aftermax = 0.0;
  

    //Float_t muon_pt_max = 0.0;
    //Float_t electron_pt_max = 0.0;

    //Float_t subjet_pt_max = 0.0;
    //Float_t subjet_pt_next_to_max = 0.0;

    //TLorentzVector zeroth_reco_ak4_jet;
    //TLorentzVector first_reco_ak4_jet;

    TLorentzVector zeroth_reco_ak8_jet;
    TLorentzVector first_reco_ak8_jet;
    TLorentzVector second_reco_ak8_jet;
    TLorentzVector third_reco_ak8_jet;
    TLorentzVector fourth_reco_ak8_jet;
    TLorentzVector fifth_reco_ak8_jet;
    
    TLorentzVector only_reco_ak8_jet;
    


    //TLorentzVector zeroth_reco_puppi_jet;
    //TLorentzVector first_reco_puppi_jet;

    //TLorentzVector zeroth_gen_ak4_jet;
    //TLorentzVector first_gen_ak4_jet;

    TLorentzVector zeroth_gen_ak8_jet;
    TLorentzVector first_gen_ak8_jet;


    TLorentzVector reco_jet0, reco_jet1;
    TLorentzVector merged_jet;

    TLorentzVector reco_muon;
    TLorentzVector reco_electron;
    TLorentzVector reco_lepton;

    //TLorentzVector subjet0;
    //TLorentzVector subjet1;
    //TLorentzVector subjet2;
    //TLorentzVector subjet3;
    //TLorentzVector subjet4;
    //TLorentzVector subjet5;
    //TLorentzVector subjet6;

    //zeroth_reco_ak4_jet.SetPxPyPzE(-99.9, -99.9, -99.9, -99.9);
    //first_reco_ak4_jet.SetPxPyPzE(-99.9, -99.9, -99.9, -99.9);

    /*zeroth_reco_ak8_jet.SetPxPyPzE(-99.9, -99.9, -99.9, -99.9);
    first_reco_ak8_jet.SetPxPyPzE(-99.9, -99.9, -99.9, -99.9);
    second_reco_ak8_jet.SetPxPyPzE(-99.9, -99.9, -99.9, -99.9);
    third_reco_ak8_jet.SetPxPyPzE(-99.9, -99.9, -99.9, -99.9);
    fourth_reco_ak8_jet.SetPxPyPzE(-99.9, -99.9, -99.9, -99.9);
    fifth_reco_ak8_jet.SetPxPyPzE(-99.9, -99.9, -99.9, -99.9);
    
    only_reco_ak8_jet.SetPxPyPzE(-99.9, -99.9, -99.9, -99.9);*/


    //zeroth_reco_puppi_jet.SetPxPyPzE(-99.9, -99.9, -99.9, -99.9);
    //first_reco_puppi_jet.SetPxPyPzE(-99.9, -99.9, -99.9, -99.9);

    //zeroth_gen_ak4_jet.SetPxPyPzE(-99.9, -99.9, -99.9, -99.9);
    //first_gen_ak4_jet.SetPxPyPzE(-99.9, -99.9, -99.9, -99.9);

    //zeroth_gen_ak8_jet.SetPxPyPzE(-99.9, -99.9, -99.9, -99.9);
    //first_gen_ak8_jet.SetPxPyPzE(-99.9, -99.9, -99.9, -99.9);

    reco_muon.SetPxPyPzE(-99.9, -99.9, -99.9, -99.9);
    reco_electron.SetPxPyPzE(-99.9, -99.9, -99.9, -99.9);
    reco_lepton.SetPxPyPzE(-99.9, -99.9, -99.9, -99.9);
    reco_jet0.SetPxPyPzE(-99.9, -99.9, -99.9, -99.9);
    reco_jet1.SetPxPyPzE(-99.9, -99.9, -99.9, -99.9);
    merged_jet.SetPxPyPzE(-99.9, -99.9, -99.9, -99.9);

    /*subjet0.SetPxPyPzE(-99.9, -99.9, -99.9, -99.9);
    subjet1.SetPxPyPzE(-99.9, -99.9, -99.9, -99.9);
    subjet2.SetPxPyPzE(-99.9, -99.9, -99.9, -99.9);
    subjet3.SetPxPyPzE(-99.9, -99.9, -99.9, -99.9);
    subjet4.SetPxPyPzE(-99.9, -99.9, -99.9, -99.9);
    subjet5.SetPxPyPzE(-99.9, -99.9, -99.9, -99.9);
    subjet6.SetPxPyPzE(-99.9, -99.9, -99.9, -99.9);*/

    Int_t reco_el_charge = 0;
    Int_t reco_mu_charge = 0;
    Int_t counter = -1;
    Int_t subcounter = -1;
    Int_t mu_mark = -1;
    Int_t el_mark = -1;
    Int_t jet0_mark = -1;
    Int_t jet1_mark = -1;
    Int_t merged_mark = -1;

    Bool_t lepton_found = false;
    Bool_t bad_event = false;

    extra_muon = 0;
    extra_electron = 0;
    loose_muon = 0;
    loose_electron = 0;
    loose_jet = 0;
    no_leptons = 0;
    no_jets = 0;
    lonesome_jet = 0;
    extra_merged_jet = 0;
    met_small = 0;
    bad_dijet = 0;
    undefined_met = 0;

    const Float_t tight_muon_cutoff_pt = 53.0; //was 25.0
    const Float_t tight_electron_cutoff_pt = 50.0; //was 30.0
    const Float_t loose_muon_cutoff_pt = 20.0; //was ?? channel
    const Float_t loose_electron_cutoff_pt = 20.0; //35.0; //was 20.0 channel
    const Float_t muon_pseudorapidity = 2.4; //was 2.1 (2.5 loose)
    const Float_t electron_pseudorapidity = 2.5; //was 2.5
    const Float_t jet_pseudorapidity = 2.4; //was infinity
    const Float_t mu_channel_missing_et_cutoff = 40.0; //was 25.0
    const Float_t el_channel_missing_et_cutoff = 80.0; //was 25.0
    const Float_t goodjet_pt_cutoff = 200.0; //was 40.0 (35.0 second)
    const Float_t badjet_pt_cutoff = 30.0; //was 30.0
    
    Float_t loose_lepton_cutoff_pt = -1.0;
    Float_t missing_et_cutoff = -1.0;

    //electron eta cannot be between 1.44 and 1.57 (is this absolute value?)
    //std::cout << "COUNTER\n";

    for (const pat::Muon &mu_index : *reco_muons) //tight muons
      {
	//if (bad_event) break;
	counter++;
	if (mu_index.pt() > tight_muon_cutoff_pt && fabs(mu_index.eta()) < muon_pseudorapidity)
	  {
	    if (lepton_found) {bad_event = true; extra_muon = 1;}// std::cout << "BAD EVENT: more than one muon.\n";}
	    else
	      {
		lepton_found = true;
		reco_muon.SetPxPyPzE(mu_index.px(), mu_index.py(), mu_index.pz(), mu_index.energy());
		mu_mark = counter;
		reco_mu_charge = mu_index.charge();
	      }
	  }
      }
   //std::cout << "MUON COUNTER IS: " << counter << std::endl; 
    counter = -1;

    for (const pat::Electron &el_index : *reco_electrons) //tight electrons
      {
	//if (bad_event) break;
	counter++;
	if (el_index.et() > tight_electron_cutoff_pt && fabs(el_index.eta()) < electron_pseudorapidity)
	  {
	    if (lepton_found) {bad_event = true; extra_electron = 1;} //std::cout << "BAD EVENT: more than one electron, or electron and muon.\n";}
	    else
	      {
		lepton_found = true;
		reco_electron.SetPxPyPzE(el_index.px(), el_index.py(), el_index.pz(), el_index.energy());
		el_mark = counter;
		reco_el_charge = el_index.charge();
	      }
	  }
      }
   //std::cout << "ELECTRON COUNTER IS: " << counter << std::endl;
    if (!lepton_found) {bad_event = true;} //std::cout << "BAD EVENT: no leptons found.\n"; no_leptons = 1;}

    if (el_mark != -1 && reco_el_charge < 0)
      {
	reco_lep_gen = 1;
	loose_lepton_cutoff_pt = loose_electron_cutoff_pt;
	reco_lepton = reco_electron;
	missing_et_cutoff = el_channel_missing_et_cutoff;
      }
    if (el_mark != -1 && reco_el_charge > 0)
      {
	reco_lep_gen = -1;
	loose_lepton_cutoff_pt = loose_electron_cutoff_pt;
	reco_lepton = reco_electron;
	missing_et_cutoff = el_channel_missing_et_cutoff;
      }
    if (mu_mark != -1 && reco_mu_charge < 0)
      {
	reco_lep_gen = 2;
	loose_lepton_cutoff_pt = loose_muon_cutoff_pt;
	reco_lepton = reco_muon;
	missing_et_cutoff = mu_channel_missing_et_cutoff;
      }
    if (mu_mark != -1 && reco_mu_charge > 0)
      {
	reco_lep_gen = -2;
	loose_lepton_cutoff_pt = loose_muon_cutoff_pt;
	reco_lepton = reco_muon;
	missing_et_cutoff = mu_channel_missing_et_cutoff;
      }

    
    counter = -1;

    for (const pat::Muon &mu_index : *reco_muons) //loose muons
      {
	//if (bad_event) break;
	counter++;
	if (counter == mu_mark) continue; //avoid counting the tight muon as a loose muon
	if (mu_index.pt() > loose_lepton_cutoff_pt) {bad_event = true; loose_muon = 1;}// std::cout << "BAD EVENT: loose muon.\n";}	  
      }
    //std::cout << "LOOSE MUON COUNTER IS: " << counter << std::endl;
    counter = -1;

    for (const pat::Electron &el_index : *reco_electrons) //loose electrons
      {
	//if (bad_event) break;
	counter++;
	if (counter == el_mark) continue; //avoid counting the tight electron as a loose electron
	if (el_index.pt() > loose_lepton_cutoff_pt) {bad_event = true; loose_electron = 1;} // std::cout << "BAD EVENT: loose electron.\n";}
      }
    //std::cout << "LOOSE ELECTRON COUNTER IS: " << counter << std::endl;
    ///////////////////////////////////GOOD JETS
  
	counter = -1;
	Float_t leading_jet_pt = goodjet_pt_cutoff;

	for (const pat::Jet &jet_index : *reco_ak8_jets) //serach for leading jet
	  {
	    //if (bad_event) break;
	    counter++;
	    if (jet_index.pt() > leading_jet_pt && fabs(jet_index.eta()) < jet_pseudorapidity)
	      {
		leading_jet_pt = jet_index.pt();
		reco_jet0.SetPxPyPzE(jet_index.px(), jet_index.py(), jet_index.pz(), jet_index.energy());
		jet0_mark = counter;
	      }
	  }
        //std::cout << "JET COUNTER IS: " << counter << std::endl;
	counter = -1;
	Float_t sub_leading_jet_pt = goodjet_pt_cutoff;
    
	if (jet0_mark == -1) {bad_event = true;} //std::cout << "BAD EVENT: UNMERGED: no jets found.\n"; no_jets = 1;}

	for (const pat::Jet &jet_index : *reco_ak8_jets) //search for subleading jet
	  {
	    //if (bad_event) break;
	    counter++;
	    if (counter == jet0_mark) continue; //avoid treating zeroth jet as first jet
	    if (jet_index.pt() > sub_leading_jet_pt && fabs(jet_index.eta()) < jet_pseudorapidity)
	      {	   
		sub_leading_jet_pt = jet_index.pt();
		reco_jet1.SetPxPyPzE(jet_index.px(), jet_index.py(), jet_index.pz(), jet_index.energy());
		jet1_mark = counter;
	      }
	  }
        //std::cout << "SECOND JET COUNTER IS: " << counter << std::endl;
	//////////////////////////////BAD JETS

	counter = -1;
	Float_t background_jet_pt_max = badjet_pt_cutoff;

	if (jet1_mark == -1) {bad_event = true;}// std::cout << "BAD EVENT: UNMERGED: only one jet.\n"; lonesome_jet = 1;}

	for (const pat::Jet &jet_index : *reco_ak4_jets) //make sure no loose jets
	  {
	    //if (bad_event) break;
	    counter++;
	    if (counter == jet0_mark || counter == jet1_mark) continue; //avoid treating tight jets as loose jet
	    if (jet_index.pt() > background_jet_pt_max && fabs(jet_index.eta()) < jet_pseudorapidity) {bad_event = true;}// std::cout << "BAD EVENT: UNMERGED: loose jets.\n"; loose_jet = 1;}
	  }
         //std::cout << "BAD JET COUNTER IS: " << counter << std::endl;
         //std::cout << "COUNTER\n";
	//TLorentzVector dijet_system;
	//dijet_system.SetPxPyPzE(reco_jet0.Px() + reco_jet1.Px(), reco_jet0.Py() + reco_jet1.Py(), reco_jet0.Pz() + reco_jet1.Pz(), reco_jet0.E() + reco_jet1.E());

	//if (dijet_system.Pt() <= 70.0) {bad_event = true; std::cout << "BAD EVENT: UNMERGED: dijet system too small.\n"; bad_dijet = 1;}
    
	/*//////////////////////////////////////////////////////////////MERGED STUFF
	//if boosted (merged?), reconstructed w_pt must be over 200.0
	Float_t merged_jet_pt_max = 200.0;
	counter = -1;
	const Float_t merged_otherjet_pt_max = 150;//80.0;//seems to be a killer
	
	for (const pat::Jet &jet_index : *reco_ak8_jets)
	  {
	    //if (bad_event) break;
	    counter++;
	    if (jet_index.pt() > merged_jet_pt_max)
	      {
		merged_jet_pt_max = jet_index.pt();
		merged_jet.SetPxPyPzE(jet_index.px(), jet_index.py(), jet_index.pz(), jet_index.energy());
		merged_mark = counter;
	      }
	  }

	if (merged_mark == -1) {bad_event = true; std::cout << "BAD EVENT: MERGED: no merged jets found.\n"; no_merged_jet = 1;}

	counter = -1;

	for (const pat::Jet &jet_index : *reco_ak8_jets)
	  {
	    //if (bad_event) break;
	    counter++;
	    if (counter == merged_mark) continue; //skip over already found jet
	    if (jet_index.pt() > merged_otherjet_pt_max) {bad_event = true; std::cout << "BAD EVENT: MERGED: too many merged jets.\n"; extra_merged_jet = 1;}
	  }

	counter = -1;
	
	//now for the subjets
	Float_t subjet_pt_max = 0.0;
	Float_t subjet_pt_next_to_max = 0.0;
	for (const pat::Jet &jet_index : *reco_ak8_jets)
	  {
	    if (jet_index.pt() == 0) { } //keep compiler from warning about unused variable
	    //if (bad_event) break;
	    counter++;
	    if (counter != merged_mark) continue; //this time we want only the merged jet we found
	    
	    for (const pat::Jet &subjet_index : *sub_jets)
	      {
		subcounter++;
		if (subjet_index.pt() > subjet_pt_max)
		  {
		    subjet_pt_max = subjet_index.pt();
		    reco_jet0.SetPxPyPzE(subjet_index.px(), subjet_index.py(), subjet_index.pz(), subjet_index.energy());
		    jet0_mark = subcounter;
		  }
	      }

	    subcounter = -1;
	    
	    for (const pat::Jet &subjet_index : *sub_jets)
	      {
		subcounter++;
		if (subcounter == jet0_mark) continue; //don't want to count the same subjet twice
		if (subjet_index.pt() > subjet_pt_next_to_max)
		  {
		    subjet_pt_next_to_max = subjet_index.pt();
		    reco_jet1.SetPxPyPzE(subjet_index.px(), subjet_index.py(), subjet_index.pz(), subjet_index.energy());
		    jet1_mark = subcounter;
		  }
	      }
	  }
	*///////////////////////////////////////////////////////////////END MERGED STUFF	    	
       
    //MISSING ET
    
    TLorentzVector met4vector;
    met4vector.SetPxPyPzE(met.px(), met.py(), 0.0, met.energy());

    METzCalculator calc;
    calc.SetMET(met4vector);
    calc.SetLepton(reco_lepton);

    if (abs(reco_lep_gen) == 1) calc.SetLeptonType("electron");// std::cout << "DEBUGGING: lepton is an electron.\n";}
    if (abs(reco_lep_gen) == 2) calc.SetLeptonType("muon");// std::cout << "DEBUGGING: lepton is a muon.\n";}

    Float_t met_pz = 99.9;
  
    if (!bad_event) met_pz = calc.Calculate(0);

    if (met_pz != met_pz) {bad_event = true;}// std::cout << "BAD EVENT: met pz not well defined.\n"; undefined_met = 1;}

    if (!bad_event) met4vector.SetPxPyPzE(met.px(), met.py(), met_pz, met.energy());

    if (met4vector.Et() < missing_et_cutoff) {bad_event = true;}// std::cout << "BAD EVENT: met et too small.\n"; met_small = 1;}

    if (bad_event) 
      {
        //std::cout << "BAD EVENT: reco lep generation is: " << reco_lep_gen << endl;
        //std::cout << "BAD EVENT: electron mark is: " << el_mark << endl;
        //std::cout << "BAD EVENT: muon mark is: " << mu_mark << endl;
	//std::cout << "BAD EVENT: THAT'S ALL FOLKS!\n\n";
      }

    if (no_leptons == 1 || extra_muon == 1 || extra_electron == 1)
      {
	loose_electron = 0;
	loose_muon = 0;
	met_small = 0;
      }
    if (no_jets == 1)
      {
	lonesome_jet = 0;
	loose_jet = 0;
	met_small = 0;
      }
    if (lonesome_jet == 1)
      {
	loose_jet = 0;
	met_small = 0;
      }

    if (bad_event) good_event = 0;
    else good_event = 1;


    ////////SANITY CHECK
    if (no_leptons == 1)
      {
	//counter = -1;
	Float_t pt_mu_max = 0.0;
	Float_t pt_el_max = 0.0;
	Float_t eta_of_mu_max = 0.0;
	Float_t eta_of_el_max = 0.0;
	
	for (const pat::Muon &mu_index : *reco_muons)
	  {
	    if (mu_index.pt() > pt_mu_max)// && fabs(mu_index.eta()) < muon_pseudorapidity)
	      {
		pt_mu_max = mu_index.pt();
		eta_of_mu_max = mu_index.eta(); eta_of_mu_max = eta_of_mu_max;
	      }
	  }

	for (const pat::Electron &el_index : *reco_electrons)
	  {
	    if (el_index.et() > pt_el_max)// && fabs(el_index.eta()) < electron_pseudorapidity)
	      {
		pt_el_max = el_index.pt();
		eta_of_el_max = el_index.eta(); eta_of_el_max = eta_of_el_max;
	      }
	  }

	//std::cout << "\nSANITY CHECK\n";
	//std::cout << "SANITY CHECK: START\n";
	//std::cout << "SANITY CHECK: muon max pt is " << pt_mu_max << std::endl;
	//std::cout << "SANITY CHECK: muon eta is " << eta_of_mu_max << std::endl;
	//std::cout << "SANITY CHECK: electron max pt is " << pt_el_max << std::endl;
	//std::cout << "SANITY CHECK: electron eta is " << eta_of_el_max << std::endl;
	//std::cout << "SANITY CHECK: lepton pt is " << lepton_pt << std::endl;
	//std::cout << "SANITY CHECK: lepton eta is " << lepton_eta << std::endl;
	//std::cout << "SANITY CHECK: lep gen is " << lep_gen << std::endl;
	//std::cout << "SANITY CHECK: END\n";
	//std::cout << "SANITY CHECK\n";
    
      }

    //if (no_leptons == 1 && lepton_found) std::cout << "SANITY CHECK: HELP! SOMETHING'S NOT RIGHT! we found a lepton, but have none?\n";
    //if (no_leptons == 0 && !lepton_found) std::cout << "SANITY CHECK: HELP! SOMETHING'S NOT RIGHT! we didn't find a lepton, but have one?\n";
    /*
    muon_count = 0;
    for (const pat::Muon &mu_index : *reco_muons)
      {
	muon_count++;
	if (mu_index.pt() > muon_pt_max)
	  {
	    muon_pt_max = mu_index.pt();
	    reco_muon.SetPxPyPzE(mu_index.px(), mu_index.py(), mu_index.pz(), mu_index.energy());
	    reco_mu_charge = mu_index.charge();
	  }
      }
    electron_count = 0;
    for (const pat::Electron &el_index : *reco_electrons)
      {
	electron_count++;
	if (el_index.pt() > electron_pt_max)
	  {
	    electron_pt_max = el_index.pt();
	    reco_electron.SetPxPyPzE(el_index.px(), el_index.py(), el_index.pz(), el_index.energy());
	    reco_el_charge = el_index.charge();
	  }
      }
    */
    

    //std::cout << "electron charge is " << reco_el_charge << std::endl;
    //std::cout << "muon charge is " << reco_mu_charge << std::endl;

    /*reco_ak4_jet_count = 0;
      for (const pat::Jet &jet_index : *reco_ak4_jets)
      {
      reco_ak4_jet_count++;
      if (jet_index.pt() > reco_ak4_jet_pt_max)
      {
      reco_ak4_jet_pt_max = jet_index.pt();
      zeroth_reco_ak4_jet.SetPxPyPzE(jet_index.px(), jet_index.py(), jet_index.pz(), jet_index.energy());
      }
      else if (jet_index.pt() > reco_ak4_jet_pt_next_to_max)
      {
      reco_ak4_jet_pt_next_to_max = jet_index.pt();
      first_reco_ak4_jet.SetPxPyPzE(jet_index.px(), jet_index.py(), jet_index.pz(), jet_index.energy());
      }
      }*/
  
   
    /* for (const pat::Jet &jet_index : *reco_ak8_jets) // line 4633
       {
	
       if (jet_index.pt() > reco_ak8_jet_pt_max)
       {
       reco_ak8_jet_pt_max = jet_index.pt();
	  
       //auto const &subjetautoconst = jet_index.subjets("SoftDrop");
       //for (auto const &subjet_index : subjetautoconst) // line 4951
       // for (const pat::Jet &isub : *AK8CHSsub) -- line 5117
       for (const pat::Jet &subjet_index : *sub_jets)
       {
       subjet_count++;
       //if (subjet_index.correctedP4(0).pt() > subjet_pt_max)
       {
       subjet_pt_max = subjet_index.correctedP4(0).pt();
       subjet0.SetPxPyPzE(subjet_index.correctedP4(0).px(), subjet_index.correctedP4(0).py(), subjet_index.correctedP4(0).pz(), subjet_index.correctedP4(0).energy());
       }
       else if (subjet_index.correctedP4(0).pt() > subjet_pt_next_to_max)
       {
       subjet_pt_next_to_max = subjet_index.correctedP4(0).pt();
       subjet1.SetPxPyPzE(subjet_index.correctedP4(0).px(), subjet_index.correctedP4(0).py(), subjet_index.correctedP4(0).pz(), subjet_index.correctedP4(0).energy());
       }
       }
       }
       }*/

    //jets
    /*  reco_ak8_jet_count = 0;
    Float_t deltamin = 1.0;
    
    //subjets
    subjet_count = 0;
    Int_t subjet_count1 = 0; //counters to keep track of subsequent subjet searches
    Int_t subjet_count2 = 0;
    Int_t subjet_count3 = 0;
    Int_t subjet_count4 = 0;
    Int_t subjet_count5 = 0;
    Int_t subjet_count6 = 0;
    
    Float_t deltasubmin0 = 1.0;
    Float_t deltasubmin1 = 1.0;
    Float_t deltasubmin2 = 10.0;
    Float_t deltasubmin3 = 10.0;
    Float_t deltasubmin4 = 10.0;
    Float_t deltasubmin5 = 10.0;
    Float_t deltasubmin6 = 10.0;
    
    TString smallest = "fixme";
    
    Int_t mark0 = -1; //marks for the 0th, 1st, etc subjets, to avoid double counting
    Int_t mark1 = -1;
    Int_t mark2 = -1;
    Int_t mark3 = -1;
    Int_t mark4 = -1;
    Int_t mark5 = -1;
    Int_t mark6 = -1;

    for (const pat::Jet &jet_index : *reco_ak8_jets) // big for loop
      {
	reco_ak8_jet_count++;
	TLorentzVector dummy4vector;
	dummy4vector.SetPxPyPzE(jet_index.px(),jet_index.py(),jet_index.pz(),jet_index.energy());
	
	if (dummy4vector.DeltaR(quark) < deltamin || dummy4vector.DeltaR(antiquark) < deltamin) //this entire if statement, including all for loops, redone if better jet found within for loop
	  {
	    only_reco_ak8_jet.SetPxPyPzE(jet_index.px(), jet_index.py(), jet_index.pz(), jet_index.energy());
	    if (dummy4vector.DeltaR(quark) <  dummy4vector.DeltaR(antiquark)) deltamin = dummy4vector.DeltaR(quark);
	    else deltamin = dummy4vector.DeltaR(antiquark);
	    subjet_count = 0;
	    for (const pat::Jet &subjet_index : *sub_jets) //search within a single jet
	      {
		subjet_count++;
		TLorentzVector dummysub4vector;
		dummysub4vector.SetPxPyPzE(subjet_index.correctedP4(0).px(), subjet_index.correctedP4(0).py(), subjet_index.correctedP4(0).pz(), subjet_index.correctedP4(0).energy());
		if (dummysub4vector.DeltaR(quark) < deltasubmin0 || dummysub4vector.DeltaR(antiquark) < deltasubmin0) //if statement redone if better subjet is found
		  {
		    mark0 = subjet_count;
		    subjet0.SetPxPyPzE(subjet_index.correctedP4(0).px(), subjet_index.correctedP4(0).py(), subjet_index.correctedP4(0).pz(), subjet_index.correctedP4(0).energy());
		    if (dummysub4vector.DeltaR(quark) < dummysub4vector.DeltaR(antiquark)) {deltasubmin0 = dummysub4vector.DeltaR(quark); smallest = "quark";}
		    else {deltasubmin0 = dummysub4vector.DeltaR(antiquark); smallest = "antiquark";} //keep track of which parton it was closer to
		  }
	      } //end for loop for zeroth subjet
	    subjet_count1 = 0;
	    for (const pat::Jet &subjet_index : *sub_jets) //re-search within same jet, still inside same if statement, only look at parton not yet matched
	      {
		subjet_count1++;
		TLorentzVector dummysub4vector;
		dummysub4vector.SetPxPyPzE(subjet_index.correctedP4(0).px(), subjet_index.correctedP4(0).py(), subjet_index.correctedP4(0).pz(), subjet_index.correctedP4(0).energy());
		if (smallest == "quark" && dummysub4vector.DeltaR(antiquark) < deltasubmin1 && subjet_count1 != mark0) //quark found first, case
		  {
		    mark1 = subjet_count1;
		    subjet1.SetPxPyPzE(subjet_index.correctedP4(0).px(), subjet_index.correctedP4(0).py(), subjet_index.correctedP4(0).pz(), subjet_index.correctedP4(0).energy());
		    deltasubmin1 = dummysub4vector.DeltaR(antiquark);
		  }
		if (smallest == "antiquark" && dummysub4vector.DeltaR(quark) < deltasubmin1 && subjet_count1 != mark0) //other case, when antiquark found first
		  {
		    mark1 = subjet_count1;
		    subjet1.SetPxPyPzE(subjet_index.correctedP4(0).px(), subjet_index.correctedP4(0).py(), subjet_index.correctedP4(0).pz(), subjet_index.correctedP4(0).energy());
		    deltasubmin1 = dummysub4vector.DeltaR(quark);
		  }
	      } //end for loop for first subjet

	    ///////begin searching for extranious subjets
	    
	    subjet_count2 = 0;
	    for (const pat::Jet &subjet_index : *sub_jets) //further group by distance from either parton, whichever closest
	      {
		subjet_count2++;
		TLorentzVector dummysub4vector;
		dummysub4vector.SetPxPyPzE(subjet_index.correctedP4(0).px(), subjet_index.correctedP4(0).py(), subjet_index.correctedP4(0).pz(), subjet_index.correctedP4(0).energy());
		if (subjet_count2 != mark0 && subjet_count2 != mark1 &&
		    (dummysub4vector.DeltaR(quark) < deltasubmin2 || dummysub4vector.DeltaR(antiquark) < deltasubmin2))
		  {
		    mark2 = subjet_count2;
		    subjet2.SetPxPyPzE(subjet_index.correctedP4(0).px(), subjet_index.correctedP4(0).py(), subjet_index.correctedP4(0).pz(), subjet_index.correctedP4(0).energy());
		    if (dummysub4vector.DeltaR(quark) < dummysub4vector.DeltaR(antiquark)) deltasubmin2 = dummysub4vector.DeltaR(quark);
		    else deltasubmin2 = dummysub4vector.DeltaR(antiquark);
		  }
	      }
	    
	    subjet_count3 = 0;
	    for (const pat::Jet &subjet_index : *sub_jets) //further group by distance from either parton, whichever closest
	      {
		subjet_count3++;
		TLorentzVector dummysub4vector;
		dummysub4vector.SetPxPyPzE(subjet_index.correctedP4(0).px(), subjet_index.correctedP4(0).py(), subjet_index.correctedP4(0).pz(), subjet_index.correctedP4(0).energy());
		if (subjet_count3 != mark0 && subjet_count3 != mark1 && subjet_count3 != mark2 &&
		    (dummysub4vector.DeltaR(quark) < deltasubmin3 || dummysub4vector.DeltaR(antiquark) < deltasubmin3))
		  {
		    mark3 = subjet_count3;
		    subjet3.SetPxPyPzE(subjet_index.correctedP4(0).px(), subjet_index.correctedP4(0).py(), subjet_index.correctedP4(0).pz(), subjet_index.correctedP4(0).energy());
		    if (dummysub4vector.DeltaR(quark) < dummysub4vector.DeltaR(antiquark)) deltasubmin3 = dummysub4vector.DeltaR(quark);
		    else deltasubmin3 = dummysub4vector.DeltaR(antiquark);
		  }
	      }

	    subjet_count4 = 0;
	    for (const pat::Jet &subjet_index : *sub_jets) //further group by distance from either parton, whichever closest
	      {
		subjet_count4++;
		TLorentzVector dummysub4vector;
		dummysub4vector.SetPxPyPzE(subjet_index.correctedP4(0).px(), subjet_index.correctedP4(0).py(), subjet_index.correctedP4(0).pz(), subjet_index.correctedP4(0).energy());
		if (subjet_count4 != mark0 && subjet_count4 != mark1 && subjet_count4 != mark2 && subjet_count4 != mark3 &&
		    (dummysub4vector.DeltaR(quark) < deltasubmin4 || dummysub4vector.DeltaR(antiquark) < deltasubmin4))
		  {
		    mark4 = subjet_count4;
		    subjet4.SetPxPyPzE(subjet_index.correctedP4(0).px(), subjet_index.correctedP4(0).py(), subjet_index.correctedP4(0).pz(), subjet_index.correctedP4(0).energy());
		    if (dummysub4vector.DeltaR(quark) < dummysub4vector.DeltaR(antiquark)) deltasubmin4 = dummysub4vector.DeltaR(quark);
		    else deltasubmin4 = dummysub4vector.DeltaR(antiquark);
		  }
	      }

	    subjet_count5 = 0;
	    for (const pat::Jet &subjet_index : *sub_jets) //further group by distance from either parton, whichever closest
	      {
		subjet_count5++;
		TLorentzVector dummysub4vector;
		dummysub4vector.SetPxPyPzE(subjet_index.correctedP4(0).px(), subjet_index.correctedP4(0).py(), subjet_index.correctedP4(0).pz(), subjet_index.correctedP4(0).energy());
		if (subjet_count5 != mark0 && subjet_count5 != mark1 && subjet_count5 != mark2 && subjet_count5 != mark3 && subjet_count5 != mark4 &&
		    (dummysub4vector.DeltaR(quark) < deltasubmin5 || dummysub4vector.DeltaR(antiquark) < deltasubmin5))
		  {
		    mark5 = subjet_count5;
		    subjet5.SetPxPyPzE(subjet_index.correctedP4(0).px(), subjet_index.correctedP4(0).py(), subjet_index.correctedP4(0).pz(), subjet_index.correctedP4(0).energy());
		    if (dummysub4vector.DeltaR(quark) < dummysub4vector.DeltaR(antiquark)) deltasubmin5 = dummysub4vector.DeltaR(quark);
		    else deltasubmin5 = dummysub4vector.DeltaR(antiquark);
		  }
	      }

	    subjet_count6 = 0;
	    for (const pat::Jet &subjet_index : *sub_jets) //further group by distance from either parton, whichever closest
	      {
		subjet_count6++;
		TLorentzVector dummysub4vector;
		dummysub4vector.SetPxPyPzE(subjet_index.correctedP4(0).px(), subjet_index.correctedP4(0).py(), subjet_index.correctedP4(0).pz(), subjet_index.correctedP4(0).energy());
		if (subjet_count6 != mark0 && subjet_count6 != mark1 && subjet_count6 != mark2 && subjet_count6 != mark3 && subjet_count6 != mark4 && subjet_count6 != mark5 &&
		    (dummysub4vector.DeltaR(quark) < deltasubmin6 || dummysub4vector.DeltaR(antiquark) < deltasubmin6))
		  {
		    mark6 = subjet_count6;
		    subjet6.SetPxPyPzE(subjet_index.correctedP4(0).px(), subjet_index.correctedP4(0).py(), subjet_index.correctedP4(0).pz(), subjet_index.correctedP4(0).energy());
		    if (dummysub4vector.DeltaR(quark) < dummysub4vector.DeltaR(antiquark)) deltasubmin6 = dummysub4vector.DeltaR(quark);
		    else deltasubmin6 = dummysub4vector.DeltaR(antiquark);
		  }
	      }
	    
	  }//end big if statement
      }
    mark6 = mark6;	
    std::cout << "SUBJETS: smallest is " << smallest << std::endl;
    std::cout << "SUBJETS: deltamin is " << deltamin << std::endl;
    std::cout << "SUBJETS: deltasubmin0 is " << deltasubmin0 << std::endl;
    std::cout << "SUBJETS: deltasubmin1 is " << deltasubmin1 << std::endl;
    std::cout << "SUBJETS: mark0 is " << mark0 << std::endl;
    std::cout << "SUBJETS: mark1 is " << mark1 << std::endl;

    Int_t reco_ak8_jet_count0 = 0;
    Int_t reco_ak8_jet_count1 = 0;
    Int_t reco_ak8_jet_count2 = 0;
    Int_t reco_ak8_jet_count3 = 0;
    
    Float_t deltamin0 = 1.0;
    Float_t deltamin1 = 10.0;
    Float_t deltamin2 = 10.0;
    Float_t deltamin3 = 10.0;

    Int_t Mark0 = -1; //marks for the 0th, 1st, etc jets, to avoid double counting
    Int_t Mark1 = -1;
    Int_t Mark2 = -1;
    Int_t Mark3 = -1;
 

    for (const pat::Jet &jet_index : *reco_ak8_jets)
      {
	reco_ak8_jet_count0++;
	TLorentzVector dummy4vector;
	dummy4vector.SetPxPyPzE(jet_index.px(),jet_index.py(),jet_index.pz(),jet_index.energy());

	if (dummy4vector.DeltaR(quark) < deltamin0 || dummy4vector.DeltaR(antiquark) < deltamin0)
	  {
	    zeroth_reco_ak8_jet.SetPxPyPzE(jet_index.px(), jet_index.py(), jet_index.pz(), jet_index.energy());
	    if (dummy4vector.DeltaR(quark) < dummy4vector.DeltaR(antiquark)) deltamin0 = dummy4vector.DeltaR(quark);
	    else deltamin0 = dummy4vector.DeltaR(antiquark);
	    Mark0 = reco_ak8_jet_count0;
	  }				      	    
      }

    for (const pat::Jet &jet_index : *reco_ak8_jets)
      {
	reco_ak8_jet_count1++;
	TLorentzVector dummy4vector;
	dummy4vector.SetPxPyPzE(jet_index.px(),jet_index.py(),jet_index.pz(),jet_index.energy());

	if ((dummy4vector.DeltaR(quark) < deltamin1 || dummy4vector.DeltaR(antiquark) < deltamin1) &&
	    reco_ak8_jet_count1 != Mark0)
	  {
	    first_reco_ak8_jet.SetPxPyPzE(jet_index.px(), jet_index.py(), jet_index.pz(), jet_index.energy());
	    if (dummy4vector.DeltaR(quark) < dummy4vector.DeltaR(antiquark)) deltamin1 = dummy4vector.DeltaR(quark);
	    else deltamin1 = dummy4vector.DeltaR(antiquark);
	    Mark1 = reco_ak8_jet_count1;
	  }				      	    
      }

    
    for (const pat::Jet &jet_index : *reco_ak8_jets)
      {
	reco_ak8_jet_count2++;
	TLorentzVector dummy4vector;
	dummy4vector.SetPxPyPzE(jet_index.px(),jet_index.py(),jet_index.pz(),jet_index.energy());

	if ((dummy4vector.DeltaR(quark) < deltamin2 || dummy4vector.DeltaR(antiquark) < deltamin2) &&
	    reco_ak8_jet_count2 != Mark0 && reco_ak8_jet_count2 != Mark1)
	  {
	    second_reco_ak8_jet.SetPxPyPzE(jet_index.px(), jet_index.py(), jet_index.pz(), jet_index.energy());
	    if (dummy4vector.DeltaR(quark) < dummy4vector.DeltaR(antiquark)) deltamin2 = dummy4vector.DeltaR(quark);
	    else deltamin2 = dummy4vector.DeltaR(antiquark);
	    Mark2 = reco_ak8_jet_count2;
	  }				      	    
      }


    for (const pat::Jet &jet_index : *reco_ak8_jets)
      {
	reco_ak8_jet_count3++;
	TLorentzVector dummy4vector;
	dummy4vector.SetPxPyPzE(jet_index.px(),jet_index.py(),jet_index.pz(),jet_index.energy());

	if ((dummy4vector.DeltaR(quark) < deltamin3 || dummy4vector.DeltaR(antiquark) < deltamin3) &&
	    reco_ak8_jet_count3 != Mark0 && reco_ak8_jet_count3 != Mark1 && reco_ak8_jet_count3 != Mark2)
	  {
	    third_reco_ak8_jet.SetPxPyPzE(jet_index.px(), jet_index.py(), jet_index.pz(), jet_index.energy());
	    if (dummy4vector.DeltaR(quark) < dummy4vector.DeltaR(antiquark)) deltamin3 = dummy4vector.DeltaR(quark);
	    else deltamin3 = dummy4vector.DeltaR(antiquark);
	    Mark3 = reco_ak8_jet_count3;
	  }				      	    
      }
    Mark3 = Mark3;*/
/*
    //4-BODY MASS SCHEME
    
    Float_t mindiff = 1000.0;
    TLorentzVector goodjet0;
    TLorentzVector goodjet1;
    TLorentzVector goodjetonly;

    TLorentzVector goodsubjet0;
    TLorentzVector goodsubjet1;
    
    Int_t index0 = -1;
    Int_t index1 = -1;
    Int_t subindex0 = -1;
    Int_t subindex1 = -1;

    Bool_t found = false;
    Bool_t subfound = false;

    for (const pat::Jet &jet_index0 : *reco_ak8_jets)
      {
	index0++;
	TLorentzVector dummy0;
	dummy0.SetPxPyPzE(jet_index0.px(), jet_index0.py(), jet_index0.pz(), jet_index0.energy());
	for (const pat::Jet &jet_index1 : *reco_ak8_jets)
	  {
	    index1++;
	     if (index0 >= index1) continue;
	    TLorentzVector dummy1;
	    dummy1.SetPxPyPzE(jet_index1.px(), jet_index1.py(), jet_index1.pz(), jet_index1.energy());
	    if (fabs((FourBodyMass(dummy0, dummy1) - Wmass)) < mindiff)
	      {
		found = true;
		mindiff = fabs((FourBodyMass(dummy0, dummy1) - Wmass));
		goodjet0.SetPxPyPzE(jet_index0.px(), jet_index0.py(), jet_index0.pz(), jet_index0.energy());
		goodjet1.SetPxPyPzE(jet_index1.px(), jet_index1.py(), jet_index1.pz(), jet_index1.energy());
	      }
	  }
      }

    if (!found)
      {
	for (const pat::Jet jet_index : *reco_ak8_jets)
	  {
	    //possibly some condition here?
	    for (const pat::Jet subjet_index0 : *sub_jets)
	      {
		subindex0++;
		TLorentzVector dummy0;
		dummy0.SetPxPyPzE(subjet_index0.px(), subjet_index0.py(), subjet_index0.pz(), subjet_index0.energy());
		for (const pat::Jet subjet_index1 : *sub_jets)
		  {
		    subindex1++;
		    if (subindex0 >= subindex1) continue;
		    TLorentzVector dummy1;
		    dummy1.SetPxPyPzE(subjet_index1.px(), subjet_index1.py(), subjet_index1.pz(), subjet_index1.energy());
		    if (fabs((FourBodyMass(dummy0, dummy1) - Wmass)) < mindiff)
		      {
			subfound = true;
			mindiff = fabs((FourBodyMass(dummy0, dummy1) - Wmass));
			goodsubjet0.SetPxPyPzE(subjet_index0.px(), subjet_index0.py(), subjet_index0.pz(), subjet_index0.energy());
			goodsubjet1.SetPxPyPzE(subjet_index1.px(), subjet_index1.py(), subjet_index1.pz(), subjet_index1.energy());
			goodjetonly.SetPxPyPzE(jet_index.px(), jet_index.py(), jet_index.pz(), jet_index.energy());
		      }
		  }
	      }
	  }
	  
      }

    //for (const pat::Jet jet_index : *sub_jet)
    // {
    //	std::cout << "If this works, something's wrong.\n";
    //  }

    if (!(found || subfound)) good_jet_event = false;
    else good_jet_event = true;*/
	  
    /*reco_ak8_jet_pt_max = 0.0;
   
      for (const pat::Jet &jet_index : *reco_ak8_jets)
       {
	 //reco_ak8_jet_count++;
	 if (jet_index.pt() > reco_ak8_jet_pt_max)
	   {
	     reco_ak8_jet_pt_max = jet_index.pt();
	     zeroth_reco_ak8_jet.SetPxPyPzE(jet_index.px(), jet_index.py(), jet_index.pz(), jet_index.energy());
	   }
	 else if (jet_index.pt() > reco_ak8_jet_pt_next_to_max)
	   {
	     reco_ak8_jet_pt_next_to_max = jet_index.pt();
	     first_reco_ak8_jet.SetPxPyPzE(jet_index.px(), jet_index.py(), jet_index.pz(), jet_index.energy());
	   }
	 else if (jet_index.pt() > num2aftermax)
	   {
	     num2aftermax = jet_index.pt();
	     second_reco_ak8_jet.SetPxPyPzE(jet_index.px(), jet_index.py(), jet_index.pz(), jet_index.energy());
	   }
	 else if (jet_index.pt() > num3aftermax)
	   {
	     num3aftermax = jet_index.pt();
	     third_reco_ak8_jet.SetPxPyPzE(jet_index.px(), jet_index.py(), jet_index.pz(), jet_index.energy());
	   }
	 else if (jet_index.pt() > num4aftermax)
	   {
	     num4aftermax = jet_index.pt();
	     fourth_reco_ak8_jet.SetPxPyPzE(jet_index.px(), jet_index.py(), jet_index.pz(), jet_index.energy());
	   }
	 else if (jet_index.pt() > num5aftermax)
	   {
	     num5aftermax = jet_index.pt();
	     fifth_reco_ak8_jet.SetPxPyPzE(jet_index.px(), jet_index.py(), jet_index.pz(), jet_index.energy());
	   }
	   }*/

    /*reco_puppi_jet_count = 0;
      for (const pat::Jet &jet_index : *reco_puppi_jets)
      {
      reco_puppi_jet_count++;
      if (jet_index.pt() > reco_puppi_jet_pt_max)
      {
      reco_puppi_jet_pt_max = jet_index.pt();
      zeroth_reco_puppi_jet.SetPxPyPzE(jet_index.px(), jet_index.py(), jet_index.pz(), jet_index.energy());
      }
      else if (jet_index.pt() > reco_puppi_jet_pt_next_to_max)
      {
      reco_puppi_jet_pt_next_to_max = jet_index.pt();
      first_reco_puppi_jet.SetPxPyPzE(jet_index.px(), jet_index.py(), jet_index.pz(), jet_index.energy());
      }
      }*/
    
    /*gen_ak4_jet_count = 0;
      for (const pat::Jet &jet_index : *gen_ak4_jets)
      {
      gen_ak4_jet_count++;
      if (jet_index.pt() > gen_ak4_jet_pt_max)
      {
      gen_ak4_jet_pt_max = jet_index.pt();
      zeroth_gen_ak4_jet.SetPxPyPzE(jet_index.px(), jet_index.py(), jet_index.pz(), jet_index.energy());
      }
      else if (jet_index.pt() > gen_ak4_jet_pt_next_to_max)
      {
      gen_ak4_jet_pt_next_to_max = jet_index.pt();
      first_gen_ak4_jet.SetPxPyPzE(jet_index.px(), jet_index.py(), jet_index.pz(), jet_index.energy());
      }
      }*/
/*
    gen_ak8_jet_count = 0;
    for (const pat::Jet &jet_index : *gen_ak8_jets)
      {
	gen_ak8_jet_count++;
	if (jet_index.pt() > gen_ak8_jet_pt_max)
	  {
	    gen_ak8_jet_pt_max = jet_index.pt();
	    zeroth_gen_ak8_jet.SetPxPyPzE(jet_index.px(), jet_index.py(), jet_index.pz(), jet_index.energy());
	  }
	else if (jet_index.pt() > gen_ak8_jet_pt_next_to_max)
	  {
	    gen_ak8_jet_pt_next_to_max = jet_index.pt();
	    first_gen_ak8_jet.SetPxPyPzE(jet_index.px(), jet_index.py(), jet_index.pz(), jet_index.energy());
	  }
      }

*/
    //std::cout << "DEBUGGING:\n";
    //std::cout << "DEBUGGING:\n";
    //std::cout << "DEBUGGING:\n";
    //std::cout << "DEBUGGING: START.\n";
  
    //if (abs(lep_gen) == 1) std::cout << "DEBUGGING: gen lepton is an electron.\n";
    //if (abs(lep_gen) == 2) std::cout << "DEBUGGING: gen lepton is a muon.\n";
    //if (abs(lep_gen) == 3) std::cout << "DEBUGGING: gen lepton is a tau.\n";


    
    /*if (electron_count != 0 && muon_count == 0 && reco_el_charge < 0) {reco_lep_gen =  1; std::cout << "DEBUGGING: only found an electron.\n";}
    if (electron_count != 0 && muon_count == 0 && reco_el_charge > 0) {reco_lep_gen = -1; std::cout << "DEBUGGING: only found a  positron.\n";}
    if (electron_count == 0 && muon_count != 0 && reco_mu_charge < 0) {reco_lep_gen =  2; std::cout << "DEBUGGING: only found a  muon.\n";}
    if (electron_count == 0 && muon_count != 0 && reco_mu_charge > 0) {reco_lep_gen = -2; std::cout << "DEBUGGING: only found an antimuon.\n";}
    if (electron_count != 0 && muon_count != 0)
      {
	if (reco_muon.Pt() >  reco_electron.Pt() && reco_mu_charge < 0) {reco_lep_gen =  2; std::cout << "DEBUGGING: muon had higher pt.\n";}
	if (reco_muon.Pt() >  reco_electron.Pt() && reco_mu_charge > 0) {reco_lep_gen = -2; std::cout << "DEBUGGING: antimuon had higher pt.\n";}
	if (reco_muon.Pt() <  reco_electron.Pt() && reco_el_charge < 0) {reco_lep_gen =  1; std::cout << "DEBUGGING: electron had higher pt.\n";}
	if (reco_muon.Pt() <  reco_electron.Pt() && reco_el_charge > 0) {reco_lep_gen = -1; std::cout << "DEBUGGING: positron had higher pt.\n";}
	if (reco_muon.Pt() == reco_electron.Pt() && reco_mu_charge < 0) {reco_lep_gen =  2; std::cout << "DEBUGGING: they had the same pt.\n";}
	if (reco_muon.Pt() == reco_electron.Pt() && reco_mu_charge > 0) {reco_lep_gen = -2; std::cout << "DEBUGGING: they had the same pt.\n";}
      }
    if (electron_count == 0 && muon_count == 0) std::cout << "DEBUGGING: no leptons found.\n";*/
  


    //if (lep_gen == 0) {std::cout << "You screwed up. You didn't set the lepton generation.\n"; exit(2);}

    
    if (reco_lep_gen == lep_gen)
      {
	Float_t recoPT = reco_lepton.Pt(); 
	//std::cout << "DEBUGGING: good lepton pt is " << recoPT << std::endl;
	//std::cout << "DEBUGGING: gen generation is " << lep_gen << std::endl;
	//std::cout << "DEBUGGING: reco generation is " << reco_lep_gen << std::endl;
	//if (recoPT <= 200.0) std::cout << "DEBUGGING: we missed a lepton.\n";
      }
  
    if (reco_lep_gen != lep_gen && reco_lep_gen != 0)
      {
	Float_t recoPT = reco_lepton.Pt(); 
	//std::cout << "DEBUGGING: fake lepton pt is " << recoPT << std::endl;
	//std::cout << "DEBUGGING: reco generation is " << reco_lep_gen << std::endl;
	//std::cout << "DEBUGGING: gen generation is " << lep_gen << std::endl;
	if (recoPT > 200.0)
	  {
	    //if (abs(lep_gen) == 3) {std::cout << "DEBUGGING: a fake got past!\n"; fakes++;}
	    //if (abs(lep_gen) != 3) {std::cout << "DEBUGGING: a fake got past and we missed a real one!\n"; fakes++;}
	  }
	//else if (abs(lep_gen) != 3) std::cout << "DEBUGGING: we missed a lepton.\n";
      }

    if (reco_lep_gen == 0 && abs(lep_gen) != 3)
      {
	//std::cout << "DEBUGGING: we missed a lepton.\n";
	//std::cout << "DEBUGGING: gen generation is " << lep_gen << std::endl;
	//std::cout << "DEBUGGING: reco generation is " << reco_lep_gen << std::endl;
      }
  
    //std::cout << "DEBUGGING: END.\n";
    //std::cout << "DEBUGGING:\n";
    //std::cout << "DEBUGGING:\n";
    //std::cout << "DEBUGGING:\n";
 
    //std::cout << "number of fakes: " << fakes << std::endl;
  
    /*
  
      double debugx = reco_lepton.Px();
      double debugy = reco_lepton.Py();
      double debugz = reco_lepton.Pz();
      double debugt = reco_lepton.E();

      double debugxe = reco_electron.Px();
      double debugye = reco_electron.Py();
      double debugze = reco_electron.Pz();
      double debugte = reco_electron.E();

      double debugxm = reco_muon.Px();
      double debugym = reco_muon.Py();
      double debugzm = reco_muon.Pz();
      double debugtm = reco_muon.E();

      if (abs(lep_gen) != 3)
      {
      std::cout << "DEBUGGING: reco lepton is: px: " << debugx << std::endl;
      std::cout << "DEBUGGING: reco lepton is: py: " << debugy << std::endl;
      std::cout << "DEBUGGING: reco lepton is: pz: " << debugz << std::endl;
      std::cout << "DEBUGGING: reco lepton is: e: "  << debugt << std::endl;
      }

      if (abs(lep_gen) == 1)
      {
      std::cout << "DEBUGGING: reco electron is: px: " << debugxe << std::endl;
      std::cout << "DEBUGGING: reco electron is: py: " << debugye << std::endl;
      std::cout << "DEBUGGING: reco electron is: pz: " << debugze << std::endl;
      std::cout << "DEBUGGING: reco electron is: e: "  << debugte << std::endl;
      }
      if (abs(lep_gen) == 2)
      {
      std::cout << "DEBUGGING: reco muon is: px: " << debugxm << std::endl;
      std::cout << "DEBUGGING: reco muon is: py: " << debugym << std::endl;
      std::cout << "DEBUGGING: reco muon is: pz: " << debugzm << std::endl;
      std::cout << "DEBUGGING: reco muon is: e: "  << debugtm << std::endl;
      }
    */
  

  

    // MARKER ANGLES
    //    ##    #    #   ####   #       ######   ####
    //   #  #   ##   #  #    #  #       #       #
    //  #    #  # #  #  #       #       #####    ####
    //  ######  #  # #  #  ###  #       #            #
    //  #    #  #   ##  #    #  #       #       #    #
    //  #    #  #    #   ####   ######  ######   ####

    if (reco_lep_gen > 0) calculateAngles(reco_lepton, met4vector, reco_jet0, reco_jet1, costheta1, costheta2, phi, costhetastar, phistar1, phistar2);
    if (reco_lep_gen < 0) calculateAngles(met4vector, reco_lepton, reco_jet0, reco_jet1, costheta1, costheta2, phi, costhetastar, phistar1, phistar2);

    if (!bad_event) {reco_costheta1 = costheta1; reco_costheta2 = costheta2; reco_phi = phi; reco_costhetastar = costhetastar; reco_phistar1 = phistar1; reco_phistar2 = phistar2;}
 
    //if (reco_lep_gen  > 0) calculateAngles(reco_lepton, met4vector, zeroth_reco_ak4_jet, first_reco_ak4_jet, costheta1, costheta2, phi, costhetastar, phistar1, phistar2);
    //if (reco_lep_gen  < 0) calculateAngles(met4vector, reco_lepton, zeroth_reco_ak4_jet, first_reco_ak4_jet, costheta1, costheta2, phi, costhetastar, phistar1, phistar2);

    //if (reco_lep_gen != 0) {reco_ak4_costheta1 = costheta1; reco_ak4_costheta2 = costheta2; reco_ak4_phi = phi; reco_ak4_costhetastar = costhetastar; reco_ak4_phistar1 = phistar1; reco_ak4_phistar2 = phistar2;}

    ///////////////////////

    //if (reco_lep_gen  > 0) calculateAngles(reco_lepton, met4vector, zeroth_reco_ak8_jet, first_reco_ak8_jet, costheta1, costheta2, phi, costhetastar, phistar1, phistar2);
    //if (reco_lep_gen  < 0) calculateAngles(met4vector, reco_lepton, zeroth_reco_ak8_jet, first_reco_ak8_jet, costheta1, costheta2, phi, costhetastar, phistar1, phistar2);

    //if (reco_lep_gen != 0) {reco_ak8_costheta1 = costheta1; reco_ak8_costheta2 = costheta2; reco_ak8_phi = phi; reco_ak8_costhetastar = costhetastar; reco_ak8_phistar1 = phistar1; reco_ak8_phistar2 = phistar2;}

    ///////////////////////

    //if (reco_lep_gen  > 0) calculateAngles(reco_lepton, met4vector, zeroth_reco_puppi_jet, first_reco_puppi_jet, costheta1, costheta2, phi, costhetastar, phistar1, phistar2);
    //if (reco_lep_gen  < 0) calculateAngles(met4vector, reco_lepton, zeroth_reco_puppi_jet, first_reco_puppi_jet, costheta1, costheta2, phi, costhetastar, phistar1, phistar2);

    //if (reco_lep_gen != 0) {reco_puppi_costheta1 = costheta1; reco_puppi_costheta2 = costheta2; reco_puppi_phi = phi; reco_puppi_costhetastar = costhetastar; reco_puppi_phistar1 = phistar1; reco_puppi_phistar2 = phistar2;}

    ///////////////////////

    //if (reco_lep_gen  > 0) calculateAngles(reco_lepton, met4vector, zeroth_gen_ak4_jet, first_gen_ak4_jet, costheta1, costheta2, phi, costhetastar, phistar1, phistar2);
    //if (reco_lep_gen  < 0) calculateAngles(met4vector, reco_lepton, zeroth_gen_ak4_jet, first_gen_ak4_jet, costheta1, costheta2, phi, costhetastar, phistar1, phistar2);

    //if (reco_lep_gen != 0) {gen_ak4_costheta1 = costheta1; gen_ak4_costheta2 = costheta2; gen_ak4_phi = phi; gen_ak4_costhetastar = costhetastar; gen_ak4_phistar1 = phistar1; gen_ak4_phistar2 = phistar2;}

    ///////////////////////

    if (reco_lep_gen  > 0) calculateAngles(reco_lepton, met4vector, zeroth_gen_ak8_jet, first_gen_ak8_jet, costheta1, costheta2, phi, costhetastar, phistar1, phistar2);
    if (reco_lep_gen  < 0) calculateAngles(met4vector, reco_lepton, zeroth_gen_ak8_jet, first_gen_ak8_jet, costheta1, costheta2, phi, costhetastar, phistar1, phistar2);

    if (reco_lep_gen != 0) {gen_ak8_costheta1 = costheta1; gen_ak8_costheta2 = costheta2; gen_ak8_phi = phi; gen_ak8_costhetastar = costhetastar; gen_ak8_phistar1 = phistar1; gen_ak8_phistar2 = phistar2;}

    ///////////////////////

    //if (reco_lep_gen  > 0) calculateAngles(reco_lepton, met4vector, subjet0, subjet1, costheta1, costheta2, phi, costhetastar, phistar1, phistar2);
    //if (reco_lep_gen  < 0) calculateAngles(met4vector, reco_lepton, subjet0, subjet1, costheta1, costheta2, phi, costhetastar, phistar1, phistar2);

    //if (reco_lep_gen != 0) {subjet_costheta1 = costheta1; subjet_costheta2 = costheta2; subjet_phi = phi; subjet_costhetastar = costhetastar; subjet_phistar1 = phistar1; subjet_phistar2 = phistar2;}

    //////////////////////////MAKING ANGLES WITH 4-BODY MASS SCHEME
    /* if (found)
      {
	if (reco_lep_gen  > 0) calculateAngles(reco_lepton, met4vector, goodjet0, goodjet1, costheta1, costheta2, phi, costhetastar, phistar1, phistar2);
	if (reco_lep_gen  < 0) calculateAngles(met4vector, reco_lepton, goodjet0, goodjet1, costheta1, costheta2, phi, costhetastar, phistar1, phistar2);

	if (reco_lep_gen != 0) {goodjet_costheta1 = costheta1; goodjet_costheta2 = costheta2; goodjet_phi = phi; goodjet_costhetastar = costhetastar; goodjet_phistar1 = phistar1; goodjet_phistar2 = phistar2;}
      }

    else if (subfound)
      {
	if (reco_lep_gen  > 0) calculateAngles(reco_lepton, met4vector, goodsubjet0, goodsubjet1, costheta1, costheta2, phi, costhetastar, phistar1, phistar2);
	if (reco_lep_gen  < 0) calculateAngles(met4vector, reco_lepton, goodsubjet0, goodsubjet1, costheta1, costheta2, phi, costhetastar, phistar1, phistar2);

	if (reco_lep_gen != 0) {goodsubjet_costheta1 = costheta1; goodsubjet_costheta2 = costheta2; goodsubjet_phi = phi; goodsubjet_costhetastar = costhetastar; goodsubjet_phistar1 = phistar1; goodsubjet_phistar2 = phistar2;}
	}*/
    
    outfile << "============EVENT=====================\n";
  
    outfile << "electrons " << electron_count << std::endl;
    outfile << "muons " << muon_count << std::endl;
    outfile << "reco_ak8_jets " << reco_ak8_jet_count << std::endl;
    outfile << "gen_ak8 jets " << gen_ak8_jet_count << std::endl;
    //outfile << "reco_puppi_jets " << reco_puppi_jet_count << std::endl;
    if (reco_lep_gen == 0 || abs(lep_gen) == 3 || abs(reco_lep_gen) != abs(lep_gen) || reco_lepton.Pt() <= 200.0) outfile << "BAD EVENT.\n";
    else outfile << "GOOD EVENT.\n";

    outfile << "\n===================ENDEVENT===================\n" << endl;

    reco_jet0_pt = reco_jet0.Pt();
    reco_jet1_pt = reco_jet1.Pt(); 
    ak8jet0_pt = zeroth_reco_ak8_jet.Pt(); 
    ak8jet1_pt = first_reco_ak8_jet.Pt();
    ak8jet2_pt = second_reco_ak8_jet.Pt();
    ak8jet3_pt = third_reco_ak8_jet.Pt(); 
    ak8jetONLY_pt = only_reco_ak8_jet.Pt();
    //pupjet0_pt = zeroth_reco_puppi_jet.Pt(); 
    //pupjet1_pt = first_reco_puppi_jet.Pt(); 
    //ak4genjet0_pt = zeroth_gen_ak4_jet.Pt(); 
    //ak4genjet1_pt = first_gen_ak4_jet.Pt(); 
    ak8genjet0_pt = zeroth_gen_ak8_jet.Pt(); 
    ak8genjet1_pt = first_gen_ak8_jet.Pt(); 
    reco_electron_pt = reco_electron.Pt(); 
    reco_muon_pt = reco_muon.Pt(); 
    reco_lepton_pt = reco_lepton.Pt(); 
    reco_met_pt = met4vector.Pt();
    //subjet0_pt = subjet0.Pt();
    //subjet1_pt = subjet1.Pt();
    //subjet2_pt = subjet2.Pt();
    //subjet3_pt = subjet3.Pt();
    //subjet4_pt = subjet4.Pt();
    //subjet5_pt = subjet5.Pt();
    //subjet6_pt = subjet6.Pt();
 
    reco_jet0_px = reco_jet0.Px();
    reco_jet1_px = reco_jet1.Px(); 
    ak8jet0_px = zeroth_reco_ak8_jet.Px(); 
    ak8jet1_px = first_reco_ak8_jet.Px();
    ak8jet2_px = second_reco_ak8_jet.Px();
    ak8jet3_px = third_reco_ak8_jet.Px(); 
    ak8jetONLY_px = only_reco_ak8_jet.Px();
    //pupjet0_px = zeroth_reco_puppi_jet.Px(); 
    //pupjet1_px = first_reco_puppi_jet.Px(); 
    //ak4genjet0_px = zeroth_gen_ak4_jet.Px(); 
    //ak4genjet1_px = first_gen_ak4_jet.Px(); 
    ak8genjet0_px = zeroth_gen_ak8_jet.Px(); 
    ak8genjet1_px = first_gen_ak8_jet.Px(); 
    reco_electron_px = reco_electron.Px(); 
    reco_muon_px = reco_muon.Px(); 
    reco_lepton_px = reco_lepton.Px(); 
    reco_met_px = met4vector.Px();
    //subjet0_px = subjet0.Px();
    //subjet1_px = subjet1.Px();
    //subjet2_px = subjet2.Px();
    //subjet3_px = subjet3.Px();
    //subjet4_px = subjet4.Px();
    //subjet5_px = subjet5.Px();
    //subjet6_px = subjet6.Px();

    reco_jet0_py = reco_jet0.Py();
    reco_jet1_py = reco_jet1.Py(); 
    ak8jet0_py = zeroth_reco_ak8_jet.Py(); 
    ak8jet1_py = first_reco_ak8_jet.Py();
    ak8jet2_py = second_reco_ak8_jet.Py();
    ak8jet3_py = third_reco_ak8_jet.Py(); 
    ak8jetONLY_py = only_reco_ak8_jet.Py();
    //pupjet0_py = zeroth_reco_puppi_jet.Py(); 
    //pupjet1_py = first_reco_puppi_jet.Py(); 
    //ak4genjet0_py = zeroth_gen_ak4_jet.Py(); 
    //ak4genjet1_py = first_gen_ak4_jet.Py(); 
    ak8genjet0_py = zeroth_gen_ak8_jet.Py(); 
    ak8genjet1_py = first_gen_ak8_jet.Py(); 
    reco_electron_py = reco_electron.Py(); 
    reco_muon_py = reco_muon.Py(); 
    reco_lepton_py = reco_lepton.Py(); 
    reco_met_py = met4vector.Py();
    //subjet0_py = subjet0.Py();
    //subjet1_py = subjet1.Py();
    //subjet2_py = subjet2.Py();
    //subjet3_py = subjet3.Py();
    //subjet4_py = subjet4.Py();
    //subjet5_py = subjet5.Py();
    //subjet6_py = subjet6.Py();
   
    reco_jet0_pz = reco_jet0.Pz();
    reco_jet1_pz = reco_jet1.Pz(); 
    ak8jet0_pz = zeroth_reco_ak8_jet.Pz(); 
    ak8jet1_pz = first_reco_ak8_jet.Pz();
    ak8jet2_pz = second_reco_ak8_jet.Pz();
    ak8jet3_pz = third_reco_ak8_jet.Pz(); 
    ak8jetONLY_pz = only_reco_ak8_jet.Pz();
    //pupjet0_pz = zeroth_reco_puppi_jet.Pz(); 
    //pupjet1_pz = first_reco_puppi_jet.Pz(); 
    //ak4genjet0_pz = zeroth_gen_ak4_jet.Pz(); 
    //ak4genjet1_pz = first_gen_ak4_jet.Pz(); 
    ak8genjet0_pz = zeroth_gen_ak8_jet.Pz(); 
    ak8genjet1_pz = first_gen_ak8_jet.Pz(); 
    reco_electron_pz = reco_electron.Pz(); 
    reco_muon_pz = reco_muon.Pz(); 
    reco_lepton_pz = reco_lepton.Pz(); 
    reco_met_pz = met4vector.Pz();
    //subjet0_pz = subjet0.Pz();
    //subjet1_pz = subjet1.Pz();
    //subjet2_pz = subjet2.Pz();
    //subjet3_pz = subjet3.Pz();
    //subjet4_pz = subjet4.Pz();
    //subjet5_pz = subjet5.Pz();
    //subjet6_pz = subjet6.Pz();

    reco_jet0_e = reco_jet0.E();
    reco_jet1_e = reco_jet1.E(); 
    ak8jet0_e = zeroth_reco_ak8_jet.E(); 
    ak8jet1_e = first_reco_ak8_jet.E();
    ak8jet2_e = second_reco_ak8_jet.E();
    ak8jet3_e = third_reco_ak8_jet.E(); 
    ak8jetONLY_e = only_reco_ak8_jet.E();
    //pupjet0_e = zeroth_reco_puppi_jet.E(); 
    //pupjet1_e = first_reco_puppi_jet.E(); 
    //ak4genjet0_e = zeroth_gen_ak4_jet.E(); 
    //ak4genjet1_e = first_gen_ak4_jet.E(); 
    ak8genjet0_e = zeroth_gen_ak8_jet.E(); 
    ak8genjet1_e = first_gen_ak8_jet.E(); 
    reco_electron_e = reco_electron.E(); 
    reco_muon_e = reco_muon.E(); 
    reco_lepton_e = reco_lepton.E(); 
    reco_met_e = met4vector.E();
    //subjet0_e = subjet0.E();
    //subjet1_e = subjet1.E();
    //subjet2_e = subjet2.E();
    //subjet3_e = subjet3.E();
    //subjet4_e = subjet4.E();
    //subjet5_e = subjet5.E();
    //subjet6_e = subjet6.E();

    reco_jet0_eta = reco_jet0.Eta();
    reco_jet1_eta = reco_jet1.Eta(); 
    ak8jet0_eta = zeroth_reco_ak8_jet.Eta(); 
    ak8jet1_eta = first_reco_ak8_jet.Eta();
    ak8jet2_eta = second_reco_ak8_jet.Eta();
    ak8jet3_eta = third_reco_ak8_jet.Eta(); 
    ak8jetONLY_eta = only_reco_ak8_jet.Eta();
    reco_jet0_et = reco_jet0.Et();
    reco_jet1_et = reco_jet1.Et(); 
    ak8jet0_et = zeroth_reco_ak8_jet.Et(); 
    ak8jet1_et = first_reco_ak8_jet.Et();
    ak8jet2_et = second_reco_ak8_jet.Et();
    ak8jet3_et = third_reco_ak8_jet.Et(); 
    ak8jetONLY_et = only_reco_ak8_jet.Et();
    //pupjet0_eta = zeroth_reco_puppi_jet.Eta(); 
    //pupjet1_eta = first_reco_puppi_jet.Eta(); 
    //ak4genjet0_eta = zeroth_gen_ak4_jet.Eta(); 
    //ak4genjet1_eta = first_gen_ak4_jet.Eta(); 
    ak8genjet0_eta = zeroth_gen_ak8_jet.Eta(); 
    ak8genjet1_eta = first_gen_ak8_jet.Eta(); 
    reco_electron_eta = reco_electron.Eta(); 
    reco_muon_eta = reco_muon.Eta(); 
    reco_lepton_eta = reco_lepton.Eta(); 
    reco_met_eta = met4vector.Eta();
    ak8genjet0_et = zeroth_gen_ak8_jet.Et(); 
    ak8genjet1_et = first_gen_ak8_jet.Et(); 
    reco_electron_et = reco_electron.Et(); 
    reco_muon_et = reco_muon.Et(); 
    reco_lepton_et = reco_lepton.Et(); 
    reco_met_et = met4vector.Et();
    //subjet0_eta = subjet0.Eta();
    //subjet1_eta = subjet1.Eta();
    //subjet2_eta = subjet2.Eta();
    //subjet3_eta = subjet3.Eta();
    //subjet4_eta = subjet4.Eta();
    //subjet5_eta = subjet5.Eta();
    //subjet6_eta = subjet6.Eta();

    reco_jet0_theta = reco_jet0.Theta();
    reco_jet1_theta = reco_jet1.Theta(); 
    ak8jet0_theta = zeroth_reco_ak8_jet.Theta(); 
    ak8jet1_theta = first_reco_ak8_jet.Theta();
    ak8jet2_theta = second_reco_ak8_jet.Theta();
    ak8jet3_theta = third_reco_ak8_jet.Theta(); 
    ak8jetONLY_theta = only_reco_ak8_jet.Theta();
    //pupjet0_theta = zeroth_reco_puppi_jet.Theta(); 
    //pupjet1_theta = first_reco_puppi_jet.Theta(); 
    //ak4genjet0_theta = zeroth_gen_ak4_jet.Theta(); 
    //ak4genjet1_theta = first_gen_ak4_jet.Theta(); 
    ak8genjet0_theta = zeroth_gen_ak8_jet.Theta(); 
    ak8genjet1_theta = first_gen_ak8_jet.Theta(); 
    reco_electron_theta = reco_electron.Theta(); 
    reco_muon_theta = reco_muon.Theta(); 
    reco_lepton_theta = reco_lepton.Theta(); 
    reco_met_theta = met4vector.Theta();
    //subjet0_theta = subjet0.Theta();
    //subjet1_theta = subjet1.Theta();
    //subjet2_theta = subjet2.Theta();
    //subjet3_theta = subjet3.Theta();
    //subjet4_theta = subjet4.Theta();
    //subjet5_theta = subjet5.Theta();
    //subjet6_theta = subjet6.Theta();

    reco_jet0_phi = reco_jet0.Phi();
    reco_jet1_phi = reco_jet1.Phi(); 
    ak8jet0_phi = zeroth_reco_ak8_jet.Phi(); 
    ak8jet1_phi = first_reco_ak8_jet.Phi();
    ak8jet2_phi = second_reco_ak8_jet.Phi();
    ak8jet3_phi = third_reco_ak8_jet.Phi(); 
    ak8jetONLY_phi = only_reco_ak8_jet.Phi();
    //pupjet0_phi = zeroth_reco_puppi_jet.Phi(); 
    //pupjet1_phi = first_reco_puppi_jet.Phi(); 
    //ak4genjet0_phi = zeroth_gen_ak4_jet.Phi(); 
    //ak4genjet1_phi = first_gen_ak4_jet.Phi(); 
    ak8genjet0_phi = zeroth_gen_ak8_jet.Phi(); 
    ak8genjet1_phi = first_gen_ak8_jet.Phi(); 
    reco_electron_phi = reco_electron.Phi(); 
    reco_muon_phi = reco_muon.Phi(); 
    reco_lepton_phi = reco_lepton.Phi(); 
    reco_met_phi = met4vector.Phi();
    //subjet0_phi = subjet0.Phi();
    //subjet1_phi = subjet1.Phi();
    //subjet2_phi = subjet2.Phi();
    //subjet3_phi = subjet3.Phi();
    //subjet4_phi = subjet4.Phi();
    //subjet5_phi = subjet5.Phi();
    //subjet6_phi = subjet6.Phi();

    reco_jet0_y = reco_jet0.Y();
    reco_jet1_y = reco_jet1.Y(); 
    ak8jet0_y = zeroth_reco_ak8_jet.Y(); 
    ak8jet1_y = first_reco_ak8_jet.Y();
    ak8jet2_y = second_reco_ak8_jet.Y();
    ak8jet3_y = third_reco_ak8_jet.Y(); 
    ak8jetONLY_y = only_reco_ak8_jet.Y();
    //pupjet0_y = zeroth_reco_puppi_jet.Y(); 
    //pupjet1_y = first_reco_puppi_jet.Y(); 
    //ak4genjet0_y = zeroth_gen_ak4_jet.Y(); 
    //ak4genjet1_y = first_gen_ak4_jet.Y(); 
    ak8genjet0_y = zeroth_gen_ak8_jet.Y(); 
    ak8genjet1_y = first_gen_ak8_jet.Y(); 
    reco_electron_y = reco_electron.Y(); 
    reco_muon_y = reco_muon.Y(); 
    reco_lepton_y = reco_lepton.Y(); 
    reco_met_y = met4vector.Y();
    //subjet0_y = subjet0.Y();
    //subjet1_y = subjet1.Y();
    //subjet2_y = subjet2.Y();
    //subjet3_y = subjet3.Y();
    //subjet4_y = subjet4.Y();
    //subjet5_y = subjet5.Y();
    //subjet6_y = subjet6.Y();

    /*//////
      TEMPLATE
      _px = .Px();
      _py = .Py();
      _pz = .Pz();
      _e = .E();
      _eta = .Eta();
      _pt = .Pt();
      _theta = .Theta();
      _phi = .Phi();
      _y = .Y();
      END TEMPLATE */

  
    /*
    goodjet0_px = goodjet0.Px();
    goodjet0_py = goodjet0.Py();
    goodjet0_pz = goodjet0.Pz();
    goodjet0_e = goodjet0.E();
    goodjet0_eta = goodjet0.Eta();
    goodjet0_pt = goodjet0.Pt();
    goodjet0_theta = goodjet0.Theta();
    goodjet0_phi = goodjet0.Phi();
    goodjet0_y = goodjet0.Y();

    goodjet1_px = goodjet1.Px();
    goodjet1_py = goodjet1.Py();
    goodjet1_pz = goodjet1.Pz();
    goodjet1_e = goodjet1.E();
    goodjet1_eta = goodjet1.Eta();
    goodjet1_pt = goodjet1.Pt();
    goodjet1_theta = goodjet1.Theta();
    goodjet1_phi = goodjet1.Phi();
    goodjet1_y = goodjet1.Y();

    goodjetonly_px = goodjetonly.Px();
    goodjetonly_py = goodjetonly.Py();
    goodjetonly_pz = goodjetonly.Pz();
    goodjetonly_e = goodjetonly.E();
    goodjetonly_eta = goodjetonly.Eta();
    goodjetonly_pt = goodjetonly.Pt();
    goodjetonly_theta = goodjetonly.Theta();
    goodjetonly_phi = goodjetonly.Phi();
    goodjetonly_y = goodjetonly.Y();

    goodsubjet0_px = goodsubjet0.Px();
    goodsubjet0_py = goodsubjet0.Py();
    goodsubjet0_pz = goodsubjet0.Pz();
    goodsubjet0_e = goodsubjet0.E();
    goodsubjet0_eta = goodsubjet0.Eta();
    goodsubjet0_pt = goodsubjet0.Pt();
    goodsubjet0_theta = goodsubjet0.Theta();
    goodsubjet0_phi = goodsubjet0.Phi();
    goodsubjet0_y = goodsubjet0.Y();

    goodsubjet1_px = goodsubjet1.Px();
    goodsubjet1_py = goodsubjet1.Py();
    goodsubjet1_pz = goodsubjet1.Pz();
    goodsubjet1_e = goodsubjet1.E();
    goodsubjet1_eta = goodsubjet1.Eta();
    goodsubjet1_pt = goodsubjet1.Pt();
    goodsubjet1_theta = goodsubjet1.Theta();
    goodsubjet1_phi = goodsubjet1.Phi();
    goodsubjet1_y = goodsubjet1.Y();

    masswidth = mindiff;
    */
    // MARKER DELTAS
    //  mmmm          ""#      m                 
    //  #   "m  mmm     #    mm#mm   mmm    mmm  
    //  #    # #"  #    #      #    "   #  #   " 
    //  #    # #""""    #      #    m"""#   """m 
    //  #mmm"  "#mm"    "mm    "mm  "mm"#  "mmm" 
    //

 
    //TLorentzVector biggerquark;
    //TLorentzVector smallerquark;

    //if (quark_pt >= antiquark_pt) {biggerquark  = quark; smallerquark = antiquark;}
    //if (quark_pt <= antiquark_pt) {smallerquark = quark; biggerquark  = antiquark;}

    if (reco_lep_gen != 0) lep_delta = lepton.DeltaR(reco_lepton);
    else lep_delta = -99.9;

    if (reco_lep_gen != 0) met_delta = neutrino.DeltaR(met4vector);
    else met_delta = -99.9;

    /////////jets
   
    if (reco_ak8_jet_count > 0)
      {
	if (quark.DeltaR(zeroth_reco_ak8_jet) < antiquark.DeltaR(zeroth_reco_ak8_jet)) ak8jet_delta0 = quark.DeltaR(zeroth_reco_ak8_jet);
	else ak8jet_delta0 = antiquark.DeltaR(zeroth_reco_ak8_jet);
      }
    else ak8jet_delta0 = -99.9;

    ///

    if (reco_ak8_jet_count > 1)
      {
	if (quark.DeltaR(first_reco_ak8_jet) < antiquark.DeltaR(first_reco_ak8_jet)) ak8jet_delta1 = quark.DeltaR(first_reco_ak8_jet);
	else ak8jet_delta1 = antiquark.DeltaR(first_reco_ak8_jet);
      }
    else ak8jet_delta1 = -99.9;

    ///

    if (reco_ak8_jet_count > 2)
      {
	if (quark.DeltaR(zeroth_reco_ak8_jet) < antiquark.DeltaR(zeroth_reco_ak8_jet)) ak8jet_delta2 = quark.DeltaR(zeroth_reco_ak8_jet);
	else ak8jet_delta2 = antiquark.DeltaR(zeroth_reco_ak8_jet);
      }
    else ak8jet_delta2 = -99.9;

    ///
    
    if (reco_ak8_jet_count > 3)
      {
	if (quark.DeltaR(first_reco_ak8_jet) < antiquark.DeltaR(first_reco_ak8_jet)) ak8jet_delta3 = quark.DeltaR(first_reco_ak8_jet);
	else ak8jet_delta3 = antiquark.DeltaR(first_reco_ak8_jet);
      }
    else ak8jet_delta3 = -99.9;

    ///GEN JETS
  
    if (gen_ak8_jet_count > 0)
      {
	if (quark.DeltaR(zeroth_gen_ak8_jet) < antiquark.DeltaR(zeroth_gen_ak8_jet)) ak8genjet_delta0 = quark.DeltaR(zeroth_gen_ak8_jet);
	else ak8genjet_delta0 = antiquark.DeltaR(zeroth_gen_ak8_jet);
      }
    else ak8genjet_delta0 = -99.9;

    ///

    if (gen_ak8_jet_count > 1)
      {
	if (quark.DeltaR(first_gen_ak8_jet) < antiquark.DeltaR(first_gen_ak8_jet)) ak8genjet_delta1 = quark.DeltaR(first_gen_ak8_jet);
	else ak8genjet_delta1 = antiquark.DeltaR(first_gen_ak8_jet);
      }
    else ak8genjet_delta1 = -99.9;
   
    
    /* if (subjet_count > 0) subjet_delta0 = smallerquark.DeltaR(subjet0);
    else subjet_delta0 = -99.9;

    if (subjet_count > 1) subjet_delta1 = smallerquark.DeltaR(subjet1);
    else subjet_delta1 = -99.9;*/

    /*if (subjet0.DeltaR(quark) < subjet0.DeltaR(antiquark))
      {
	subjet_delta0 = quark.DeltaR(subjet0);
	subjet_delta1 = antiquark.DeltaR(subjet1);
      }
    else
      {
	subjet_delta0 = antiquark.DeltaR(subjet0);
	subjet_delta1 = quark.DeltaR(subjet1);
      }

    ////////extranious subjets 

    if (subjet_count > 2)
      {
	if (subjet2.DeltaR(quark) < subjet2.DeltaR(antiquark)) subjet_delta2 = quark.DeltaR(subjet2);
     
	else subjet_delta2 = antiquark.DeltaR(subjet2);
      }*/
    
      	
    //////////////////////
/*
    if (subjet_count > 3)
      {
	if (subjet3.DeltaR(quark) < subjet3.DeltaR(antiquark)) subjet_delta3 = quark.DeltaR(subjet3);
     
	else subjet_delta3 = antiquark.DeltaR(subjet3);
      }
    

    //////////////////////

    if (subjet_count > 4)
      {
	if (subjet4.DeltaR(quark) < subjet4.DeltaR(antiquark)) subjet_delta4 = quark.DeltaR(subjet4);
     
	else subjet_delta4 = antiquark.DeltaR(subjet4);
      }
    

    //////////////////////

    if (subjet_count > 5)
      {
	if (subjet5.DeltaR(quark) < subjet5.DeltaR(antiquark)) subjet_delta5 = quark.DeltaR(subjet5);
     
	else subjet_delta5 = antiquark.DeltaR(subjet5);
      }
    

    //////////////////////

    if (subjet_count > 6)
      {
	if (subjet6.DeltaR(quark) < subjet6.DeltaR(antiquark)) subjet_delta6 = quark.DeltaR(subjet6);
     
	else subjet_delta6 = antiquark.DeltaR(subjet6);
      }
 */
  
    /*if (reco_ak8_jet_count > 0) ak8_jet0lowq = smallerquark.DeltaR(zeroth_reco_ak8_jet);
    else ak8_jet0lowq  = -99.9;

    if (reco_ak8_jet_count > 0) ak8_jet0highq = biggerquark.DeltaR(zeroth_reco_ak8_jet);
    else ak8_jet0highq = -99.9;

  
    if (reco_ak8_jet_count > 1) ak8_jet1lowq = smallerquark.DeltaR(first_reco_ak8_jet);
    else ak8_jet1lowq  = -99.9;

    if (reco_ak8_jet_count > 1) ak8_jet1highq = biggerquark.DeltaR(first_reco_ak8_jet);
    else ak8_jet1highq = -99.9;

  
    if (reco_ak8_jet_count > 2) ak8_jet2lowq = smallerquark.DeltaR(second_reco_ak8_jet);
    else ak8_jet2lowq  = -99.9;

    if (reco_ak8_jet_count > 2) ak8_jet2highq = biggerquark.DeltaR(second_reco_ak8_jet);
    else ak8_jet2highq = -99.9;


    if (reco_ak8_jet_count > 3) ak8_jet3lowq = smallerquark.DeltaR(third_reco_ak8_jet);
    else ak8_jet3lowq  = -99.9;

    if (reco_ak8_jet_count > 3) ak8_jet3highq = biggerquark.DeltaR(third_reco_ak8_jet);
    else ak8_jet3highq = -99.9;*/


    /*if (reco_ak8_jet_count > 4) ak8_jet4lowq = smallerquark.DeltaR(fourth_reco_ak8_jet);
    else ak8_jet4lowq  = -99.9;

    if (reco_ak8_jet_count > 4) ak8_jet4highq = biggerquark.DeltaR(fourth_reco_ak8_jet);
    else ak8_jet4highq = -99.9;


    if (reco_ak8_jet_count > 5) ak8_jet5lowq = smallerquark.DeltaR(fifth_reco_ak8_jet);
    else ak8_jet5lowq  = -99.9;

    if (reco_ak8_jet_count > 5) ak8_jet5highq = biggerquark.DeltaR(fifth_reco_ak8_jet);
    else ak8_jet5highq = -99.9;*/
 

    /*
    if (reco_ak8_jet_count > 0) DeltaTree0->Fill();
    if (reco_ak8_jet_count > 1) DeltaTree1->Fill();
    if (reco_ak8_jet_count > 2) DeltaTree2->Fill();
    if (reco_ak8_jet_count > 3) DeltaTree3->Fill();
    if (reco_ak8_jet_count > 4) DeltaTree4->Fill();
    if (reco_ak8_jet_count > 5) DeltaTree5->Fill();
    */
   
    /*
      int index1 = -1;
      int index2 = -1;
      int match1 = -1;
      int match2 = -1;
      Float_t DeltaMinSqrd = 10000.0;
      TLorentzVector first_dummy;
      TLorentzVector second_dummy;
  
      for (const pat::Jet &first_index : *reco_ak8_jets2)
      {
      index1++; std::cout << "LOOPY: index1 " << index1 << std::endl;
      first_dummy.SetPxPyPzE(first_index.px(), first_index.py(), first_index.pz(), first_index.energy());
      for (const pat::Jet &second_index : *reco_ak8_jets2)
      {
      index2++; std::cout << "LOOPY: index2 " << index2 << std::endl;
      if (index1 == index2) continue;
      second_dummy.SetPxPyPzE(second_index.px(), second_index.py(), second_index.pz(), second_index.energy());
      if (quark.DeltaR(first_dummy) * quark.DeltaR(first_dummy) + antiquark.DeltaR(second_dummy) * antiquark.DeltaR(second_dummy) < DeltaMinSqrd)
      {
      Float_t newvalue = (quark.DeltaR(first_dummy) * quark.DeltaR(first_dummy) + antiquark.DeltaR(second_dummy) * antiquark.DeltaR(second_dummy));
      std::cout << "LOOPY: newvalue is smaller than Deltaminsquared " << newvalue << std::endl;
      std::cout << "LOOPY: DeltaminSqrd was " << DeltaMinSqrd << std::endl;
      DeltaMinSqrd = quark.DeltaR(first_dummy) * quark.DeltaR(first_dummy) + antiquark.DeltaR(second_dummy) * antiquark.DeltaR(second_dummy);
      match1 = index1; std::cout << "LOOPY: match1 is now: " << match1 << std::endl;
      match2 = index2; std::cout << "LOOPY: match2 is now: " << match2 << std::endl;
      }
      else
      {
      Float_t newvalue = (quark.DeltaR(first_dummy) * quark.DeltaR(first_dummy) + antiquark.DeltaR(second_dummy) * antiquark.DeltaR(second_dummy));
      std::cout << "LOOPY: newvalue is larger than or equal to Deltaminsquared " << newvalue << std::endl;
      std::cout << "LOOPY: DeltaminSqrd still is " << DeltaMinSqrd << std::endl;
      DeltaMinSqrd = quark.DeltaR(first_dummy) * quark.DeltaR(first_dummy) + antiquark.DeltaR(second_dummy) * antiquark.DeltaR(second_dummy);
      match1 = index1; std::cout << "LOOPY: match1 is still: " << match1 << std::endl;
      match2 = index2; std::cout << "LOOPY: match2 is still: " << match2 << std::endl;
      }
      index2 = -1;
      }
      }

      //std::cout << "LOOPY: index1 final: " << index1 << std::endl;
      //std::cout << "LOOPY: index2 final: " << index2 << std::endl;
      std::cout << "LOOPY: match1 final: " << match1 << std::endl;
      std::cout << "LOOPY: match2 final: " << match2 << std::endl;
      std::cout << "LOOPY: DeltaMinSqrdfinal: " << DeltaMinSqrd << std::endl;
    

    int index = 0;
    TLorentzVector dummy;
  
    for (const pat::Jet &only_index : *reco_ak8_jets)
      {
	std::cout << "LOOPY: index " << index << std::endl;
	index++;
	dummy.SetPxPyPzE(only_index.px(), only_index.py(), only_index.pz(), only_index.energy());
	Float_t DeltaqValue = quark.DeltaR(dummy);
	Float_t DeltaantiqValue = antiquark.DeltaR(dummy);
	std::cout << "LOOPY: Delta for quark: " << DeltaqValue << std::endl;
	std::cout << "LOOPY: Delta for antiquark: " << DeltaantiqValue << std::endl;
      }
    std::cout << "LOOPY: That's all for now.\n";*/
  
      
      
    // MARKER ALL-HAD TREE
    //        d8888 888 888        888    888               888     88888888888                           
    //       d88888 888 888        888    888               888         888                               
    //      d88P888 888 888        888    888               888         888                               
    //     d88P 888 888 888        8888888888  8888b.   .d88888         888     888d888  .d88b.   .d88b.  
    //    d88P  888 888 888        888    888     "88b d88" 888         888     888P"   d8P  Y8b d8P  Y8b 
    //   d88P   888 888 888 888888 888    888 .d888888 888  888         888     888     88888888 88888888 
    //  d8888888888 888 888        888    888 888  888 Y88b 888         888     888     Y8b.     Y8b.     
    // d88P     888 888 888        888    888 "Y888888  "Y88888         888     888      "Y8888   "Y8888  
    //                                                                                                    
         
    AllHadMETpx          = met.px();                   
    AllHadMETpy          = met.py();                   
    AllHadMETpt          = met.pt();                   
    AllHadMETphi         = met.phi();                   
    AllHadMETsumET       = met.sumEt();                   
    AllHadNvtx           = nvtx;    
    AllHadNPUtrue        = nPU;           
    AllHadRho            = rho;               
    AllHadEventWeight    = 1;   
    AllHadPUweight       = puweight; 
    AllHadPUweight_MBup  = puweightUp;
    AllHadPUweight_MBdn  = puweightDn;          
    DijetMass            = (AK8jet0_P4corr + AK8jet1_P4corr).M();                                                   
    DijetMassPuppi       = (PUPPIjet0_P4corr + PUPPIjet1_P4corr).M();                                                   
    DijetDeltaR          = deltaR(AK8jet0_P4corr.Eta(), AK8jet0_P4corr.Phi(), AK8jet1_P4corr.Eta(), AK8jet1_P4corr.Phi());               
    DijetDeltaPhi        = fabs(deltaPhi(AK8jet0_P4corr.Phi(), AK8jet1_P4corr.Phi()));                 
    DijetDeltaRap        = fabs(AK8jet0_P4corr.Rapidity() -  AK8jet1_P4corr.Rapidity());

    DiGenJetMass         = (GenJetMatched0 + GenJetMatched1).M();                   
    GenTTmass            = (t1_p4+t2_p4).M();               
    HT                   = HT_AK4_pt30;                
    HT_CorrDn            = HT_AK4_pt30_corrDn;                
    HT_CorrUp            = HT_AK4_pt30_corrUp;  
    HT_PtSmearNom        = HT_AK4_pt30_smearNom;          
    HT_PtSmearUp         = HT_AK4_pt30_smearUp;               
    HT_PtSmearDn         = HT_AK4_pt30_smearDn;               
    Q2weight_CorrDn      = Q2wgt_down;              
    Q2weight_CorrUp      = Q2wgt_up;              
    NNPDF3weight_CorrDn  = NNPDF3wgt_down;              
    NNPDF3weight_CorrUp  = NNPDF3wgt_up;              
    AllHadRunNum         = iEvent.id().run();              
    AllHadLumiBlock      = iEvent.id().luminosityBlock();              
    AllHadEventNum       = iEvent.id().event();  
    if (passMETfilters) PassMETFilters       = 1;
    else PassMETFilters                      = 0;

                   


    //------------------------------------
    // WRITE TREE WITH BASELINE PT CUT AND ETA CUT
    //------------------------------------
    if (AK8jet0_P4corr.Perp() > 300 && AK8jet1_P4corr.Perp() > 300 && fabs(AK8jet0_P4corr.Rapidity()) < 2.4 && fabs(AK8jet1_P4corr.Rapidity()) < 2.4) TreeAllHad->Fill(); 



    // MARKER SEMI-LEPT TREE
    //  .d8888b.                         d8b        888                        888        88888888888                           
    // d88P  Y88b                        Y8P        888                        888            888                               
    // Y88b.                                        888                        888            888                               
    //  "Y888b.    .d88b.  88888b.d88b.  888        888       .d88b.  88888b.  888888         888     888d888  .d88b.   .d88b.  
    //     "Y88b. d8P  Y8b 888 "888 "88b 888        888      d8P  Y8b 888 "88b 888            888     888P"   d8P  Y8b d8P  Y8b 
    //       "888 88888888 888  888  888 888 888888 888      88888888 888  888 888            888     888     88888888 88888888 
    // Y88b  d88P Y8b.     888  888  888 888        888      Y8b.     888 d88P Y88b.          888     888     Y8b.     Y8b.     
    //  "Y8888P"   "Y8888  888  888  888 888        88888888  "Y8888  88888P"   "Y888         888     888      "Y8888   "Y8888  
    //                                                                888                                                       
    //                                                                888                                                       
    //                                                                888       
     
    SemiLeptMETpx                = met.px();                   
    SemiLeptMETpy                = met.py();                   
    SemiLeptMETpt                = met.pt();                   
    SemiLeptMETphi               = met.phi();                   
    SemiLeptMETsumET             = met.sumEt();                   
    SemiLeptNvtx                 = nvtx;     
    SemiLeptNPUtrue              = nPU;     
    SemiLeptRho                  = rho;               
    SemiLeptEventWeight          = 1;              
    SemiLeptPUweight       = puweight; 
    SemiLeptPUweight_MBup  = puweightUp;
    SemiLeptPUweight_MBdn  = puweightDn;

    SemiLeptGenTTmass            = (t1_p4+t2_p4).M(); 
    

    double htlep = lep0_p4.Perp() + met.pt();
    HTlep                = htlep;
    ST                   = htlep + HT_AK4_pt30;                
    ST_CorrDn            = htlep + HT_AK4_pt30_corrDn;                
    ST_CorrUp            = htlep + HT_AK4_pt30_corrUp;                
    ST_PtSmearNom        = htlep + HT_AK4_pt30_smearNom;                
    ST_PtSmearUp         = htlep + HT_AK4_pt30_smearUp;                
    ST_PtSmearDn         = htlep + HT_AK4_pt30_smearDn;                
    SemiLeptQ2weight_CorrDn      = Q2wgt_down;              
    SemiLeptQ2weight_CorrUp      = Q2wgt_up;              
    SemiLeptNNPDF3weight_CorrDn  = NNPDF3wgt_down;              
    SemiLeptNNPDF3weight_CorrUp  = NNPDF3wgt_up;              
    SemiLeptRunNum               = iEvent.id().run();              
    SemiLeptLumiBlock            = iEvent.id().luminosityBlock();              
    SemiLeptEventNum             = iEvent.id().event();              
    if (passMETfilters) SemiLeptPassMETFilters  = 1;
    else SemiLeptPassMETFilters  = 0;

    AK4dRminPt        = AK4_dRMinLep_p4.Perp();
    AK4dRminEta       = AK4_dRMinLep_p4.Eta();
    AK4dRminPhi       = AK4_dRMinLep_p4.Phi();
    AK4dRminMass      = AK4_dRMinLep_p4.M();
    AK4dRminBdisc     = AK4_dRMinLep_bdisc;
    AK4dRminLep       = AK4_dRMinLep;
  
    AK4BtagdRminPt    = AK4_btagged_dRMinLep_p4.Perp();
    AK4BtagdRminBdisc = AK4_btagged_dRMinLep_bdisc;
    AK4BtagdRminLep   = AK4_btagged_dRMinLep;

 
    if (ak4_btag_loose)  LepHemiContainsAK4BtagLoose   = 1;
    else                 LepHemiContainsAK4BtagLoose   = 0;
    if (ak4_btag_medium) LepHemiContainsAK4BtagMedium  = 1;
    else                 LepHemiContainsAK4BtagMedium  = 0;
    if (ak4_btag_tight)  LepHemiContainsAK4BtagTight   = 1;
    else                 LepHemiContainsAK4BtagTight   = 0;

    LeptonPhi   = lep0_p4.Phi(); 
    LeptonPt    = lep0_p4.Perp();  
    LeptonEta   = lep0_p4.Eta(); 
    LeptonMass  = lep0_p4.M(); 


    if (count_mu == 1 && count_el == 0)      LeptonIsMu  = 1; 
    else if (count_mu == 0 && count_el == 1) LeptonIsMu  = 0; 
    else LeptonIsMu =0;

    PtRel  = AK4_dRMinLep_p4.Perp(lep0_p4.Vect());
    MuIso  = mu0_iso04;

    Elecron_absiso            = el0_absiso;  
    Elecron_relIsoWithDBeta   = el0_relIsoWithDBeta;  
    Elecron_absiso_EA         = el0_absiso_EA;  
    Elecron_relIsoWithEA      = el0_relIsoWithEA;  

    if (mu0_isMedium) MuMedium = 1;
    else             MuMedium = 0;

    if (mu0_isTight) MuTight = 1;
    else            MuTight = 0;

    //------------------------------------
    // WRITE TREE WITH BASELINE PT CUT AND ETA CUT
    //------------------------------------


    if (GenTruth_semileptonic)  count_GenTruth_semileptonic++;
    if (count_mu  >= 1)  count_nMu_gt1++; 
    if (count_el  >= 1)  count_nEl_gt1++; 
    if (count_mu  == 1)  count_nMu_e1++; 
    if (count_el  == 1)  count_nEl_e1++; 
    if (count_lep == 1)  count_nLep_e1++; 
    if (count_lep == 1  && AK8jet_SemiLept_P4corr.Perp() > 300)  count_JetPt300++; 
    if (count_lep == 1  && AK8jet_SemiLept_P4corr.Perp() > 300 && fabs(AK8jet_SemiLept_P4corr.Rapidity()) < 2.4)  count_JetPt300Eta++; 
    if (count_lep == 1  && AK8jet_SemiLept_P4corr.Perp() > 300 && fabs(AK8jet_SemiLept_P4corr.Rapidity()) < 2.4 && AK4_dRMinLep_p4.Perp() > 20)  count_JetPt300Eta_AK4++; 
    if (count_lep == 1  && AK8jet_SemiLept_P4corr.Perp() > 300 && fabs(AK8jet_SemiLept_P4corr.Rapidity()) < 2.4 &&  mu0_p4.Perp() > 40)  count_JetPt300Eta_muPt40++; 
    if (count_lep == 1  && AK8jet_SemiLept_P4corr.Perp() > 300 && fabs(AK8jet_SemiLept_P4corr.Rapidity()) < 2.4 &&  mu0_p4.Perp() > 40 && met.pt() > 40)  count_JetPt300Eta_muPt40_MET40++; 
    if (count_lep == 1  && AK8jet_SemiLept_P4corr.Perp() > 300 && fabs(AK8jet_SemiLept_P4corr.Rapidity()) < 2.4 &&  mu0_p4.Perp() > 40 && met.pt() > 40 &&  AK4_dRMinLep_p4.Perp() > 20)  count_JetPt300Eta_muPt40_MET40_AK4++; 

    if (count_lep == 1  && verbose_)
      {
	cout << " ak8pt " << AK8jet_SemiLept_P4corr.Perp() << endl;
	cout << " mu pt " << mu0_p4.Perp() << endl;
	cout << " el pt " << el0_p4.Perp() << endl;
	cout << " met " << met.pt() << endl;
	cout << " ak4 pt " << AK4_dRMinLep_p4.Perp() << endl;
      } 
    // if (count_lep == 1 && AK8jet_SemiLept_P4corr.Perp()>200 && fabs(AK8jet_SemiLept_P4corr.Rapidity()) < 2.4 
    //  &&  lep0_p4.Perp() > 30 && met.pt() > 30) {
    //if (reco_lep_gen != 0 && (reco_ak8_jet_count > 0 || gen_ak8_jet_count > 0) && good_value && fabs(reco_muon.Eta()) < 2.1) TreeSemiLept->Fill();
    //if (!bad_event)
    if (abs(lep_gen) != 3) TreeSemiLept->Fill();
    //  } 

//#endif
}

//std::cout << "we just defined the analyze function.\n";

// ------------ method called once each job just before starting event loop  ------------
void B2GTTbarTreeMaker::beginJob()
{
  std::cout << "we are starting the beginJob function.\n";
  fPUweight = new TFile("PUweight20160908.root");
  hPUweight      = (TH1D*) fPUweight->Get("PUweight_true");
  hPUweight_MBup = (TH1D*) fPUweight->Get("PUweight_true_MBup");
  hPUweight_MBdn = (TH1D*) fPUweight->Get("PUweight_true_MBdn");
       
  std::cout << "Test PU reweight file: " << hPUweight->GetBinContent(hPUweight->GetXaxis()->FindBin(30)) << std::endl;
    

  count_GenTruth_semileptonic = 0;
  count_nMu_gt1 = 0; 
  count_nEl_gt1 = 0; 
  count_nMu_e1 = 0; 
  count_nEl_e1 = 0; 
  count_nLep_e1 = 0; 
  count_JetPt300 = 0; 
  count_JetPt300Eta = 0; 
  count_JetPt300Eta_AK4 = 0; 
  count_JetPt300Eta_muPt40 = 0; 
  count_JetPt300Eta_muPt40_MET40 = 0; 
  count_JetPt300Eta_muPt40_MET40_AK4 = 0; 

}
//std::cout << "we just defined the beginJob function.\n";
// ------------ method called once each job just after ending the event loop  ------------
void B2GTTbarTreeMaker::endJob() 
{
  std::cout << "we are starting the endJob function.\n";
  std::cout << " nEvents GenTruth semileptonic  :" << count_GenTruth_semileptonic << std::endl;
  std::cout << " nEvents nMu =1   :" << count_nMu_e1 << std::endl;
  std::cout << " nEvents nEl =1   :" << count_nEl_e1 << std::endl;
  std::cout << " nEvents nLepton =1   :" << count_nLep_e1 << std::endl;
  std::cout << " nEvents nLepton =1 && JetPt300   :" << count_JetPt300 << std::endl;
  std::cout << " nEvents nLepton =1 && JetPt300Eta   :" << count_JetPt300Eta << std::endl;
  std::cout << " nEvents nLepton =1 && JetPt300Eta && AK4pt>20   :" << count_JetPt300Eta_AK4 << std::endl;
  std::cout << " nEvents nLepton =1 && JetPt300Eta && muPt40   :" << count_JetPt300Eta_muPt40 << std::endl;
  std::cout << " nEvents nLepton =1 && JetPt300Eta && muPt40 && MET40   :" << count_JetPt300Eta_muPt40_MET40 << std::endl;
  std::cout << " nEvents nLepton =1 && JetPt300Eta && muPt40 && MET40 && AK4pt>20  :" << count_JetPt300Eta_muPt40_MET40_AK4 << std::endl;

}
//std::cout << "we just defined the endJob function.\n";
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void B2GTTbarTreeMaker::fillDescriptions(edm::ConfigurationDescriptions& descriptions) 
{
  std::cout << "we are starting the fillDescriptions function.\n";
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//std::cout << "we just defined the fillDescriptions function.\n";
Float_t FourBodyMass(TLorentzVector jet0, TLorentzVector jet1)
{
  Float_t px = jet0.Px() + jet1.Px();
  Float_t py = jet0.Py() + jet1.Py();
  Float_t pz = jet0.Pz() + jet1.Pz();
  Float_t e  = jet0.E()  + jet1.E();

  Float_t fourbodymass = std::sqrt(e * e - px * px - py * py - pz * pz);
  return fourbodymass;
}

void calculateAngles(TLorentzVector thep4M11, TLorentzVector thep4M12, TLorentzVector thep4M21, TLorentzVector thep4M22, Float_t& costheta1, Float_t& costheta2, Float_t& phi, Float_t& costhetastar, Float_t& phistar1, Float_t& phistar2)
{


  TLorentzVector thep4H = thep4M11 + thep4M12 + thep4M21 + thep4M22;
  TLorentzVector thep4Z1 = thep4M11 + thep4M12;
  TLorentzVector thep4Z2 = thep4M21 + thep4M22;

  Float_t norm;

  TVector3 boostX = -(thep4H.BoostVector());
  TLorentzVector thep4Z1inXFrame(thep4Z1);
  TLorentzVector thep4Z2inXFrame(thep4Z2);      
  thep4Z1inXFrame.Boost(boostX);
  thep4Z2inXFrame.Boost(boostX);
  TVector3 theZ1X_p3 = TVector3(thep4Z1inXFrame.X(), thep4Z1inXFrame.Y(), thep4Z1inXFrame.Z());
  TVector3 theZ2X_p3 = TVector3(thep4Z2inXFrame.X(), thep4Z2inXFrame.Y(), thep4Z2inXFrame.Z());

  // calculate phi1, phi2, costhetastar
  ///phi1 = theZ1X_p3.Phi();
  ///phi2 = theZ2X_p3.Phi();

  ///////////////////////////////////////////////
  // check for z1/z2 convention, redefine all 4 vectors with convention
  /////////////////////////////////////////////// 
  TLorentzVector p4H, p4Z1, p4M11, p4M12, p4Z2, p4M21, p4M22;
  p4Z1 = thep4Z1; p4M11 = thep4M11; p4M12 = thep4M12;
  p4Z2 = thep4Z2; p4M21 = thep4M21; p4M22 = thep4M22;
  costhetastar = theZ1X_p3.CosTheta();

  // now helicity angles................................
  // ...................................................
  TVector3 boostZ1 = -(p4Z1.BoostVector());
  TLorentzVector p4Z2Z1(p4Z2);
  p4Z2Z1.Boost(boostZ1);
  //find the decay axis
  /////TVector3 unitx_1 = -Hep3Vector(p4Z2Z1);
  TVector3 unitx_1(-p4Z2Z1.X(), -p4Z2Z1.Y(), -p4Z2Z1.Z());
  norm = 1/(unitx_1.Mag());
  unitx_1*=norm;
  //boost daughters of z2
  TLorentzVector p4M21Z1(p4M21);
  TLorentzVector p4M22Z1(p4M22);
  p4M21Z1.Boost(boostZ1);
  p4M22Z1.Boost(boostZ1);
  //create z and y axes
  /////TVector3 unitz_1 = Hep3Vector(p4M21Z1).cross(Hep3Vector(p4M22Z1));
  TVector3 p4M21Z1_p3(p4M21Z1.X(), p4M21Z1.Y(), p4M21Z1.Z());
  TVector3 p4M22Z1_p3(p4M22Z1.X(), p4M22Z1.Y(), p4M22Z1.Z());
  TVector3 unitz_1 = p4M21Z1_p3.Cross(p4M22Z1_p3);
  norm = 1/(unitz_1.Mag());
  unitz_1 *= norm;
  TVector3 unity_1 = unitz_1.Cross(unitx_1);

  //caculate theta1
  TLorentzVector p4M11Z1(p4M11);
  p4M11Z1.Boost(boostZ1);
  TVector3 p3M11(p4M11Z1.X(), p4M11Z1.Y(), p4M11Z1.Z());
  TVector3 unitM11 = p3M11.Unit();
  Float_t x_m11 = unitM11.Dot(unitx_1); Float_t y_m11 = unitM11.Dot(unity_1); Float_t z_m11 = unitM11.Dot(unitz_1);
  TVector3 M11_Z1frame(y_m11, z_m11, x_m11);
  costheta1 = M11_Z1frame.CosTheta();
  //std::cout << "theta1: " << M11_Z1frame.Theta() << std::endl;
  //////-----------------------old way of calculating phi---------------/////////
  phi = M11_Z1frame.Phi();

  //set axes for other system
  TVector3 boostZ2 = -(p4Z2.BoostVector());
  TLorentzVector p4Z1Z2(p4Z1);
  p4Z1Z2.Boost(boostZ2);
  TVector3 unitx_2(-p4Z1Z2.X(), -p4Z1Z2.Y(), -p4Z1Z2.Z());
  norm = 1/(unitx_2.Mag());
  unitx_2 *= norm;
  //boost daughters of z2
  TLorentzVector p4M11Z2(p4M11);
  TLorentzVector p4M12Z2(p4M12);
  p4M11Z2.Boost(boostZ2);
  p4M12Z2.Boost(boostZ2);
  TVector3 p4M11Z2_p3(p4M11Z2.X(), p4M11Z2.Y(), p4M11Z2.Z());
  TVector3 p4M12Z2_p3(p4M12Z2.X(), p4M12Z2.Y(), p4M12Z2.Z());
  TVector3 unitz_2 = p4M11Z2_p3.Cross(p4M12Z2_p3);
  norm = 1/(unitz_2.Mag());
  unitz_2*=norm;
  TVector3 unity_2 = unitz_2.Cross(unitx_2);
  //calcuate theta2
  TLorentzVector p4M21Z2(p4M21);
  p4M21Z2.Boost(boostZ2);
  TVector3 p3M21(p4M21Z2.X(), p4M21Z2.Y(), p4M21Z2.Z());
  TVector3 unitM21 = p3M21.Unit();
  Float_t x_m21 = unitM21.Dot(unitx_2); Float_t y_m21 = unitM21.Dot(unity_2); Float_t z_m21 = unitM21.Dot(unitz_2);
  TVector3 M21_Z2frame(y_m21, z_m21, x_m21);
  costheta2 = M21_Z2frame.CosTheta();

  // calculate phi
  //calculating phi_n
  TLorentzVector n_p4Z1inXFrame(p4Z1);
  TLorentzVector n_p4M11inXFrame(p4M11);
  n_p4Z1inXFrame.Boost(boostX);
  n_p4M11inXFrame.Boost(boostX);        
  TVector3 n_p4Z1inXFrame_unit = n_p4Z1inXFrame.Vect().Unit();
  TVector3 n_p4M11inXFrame_unit = n_p4M11inXFrame.Vect().Unit();  
  TVector3 n_unitz_1(n_p4Z1inXFrame_unit);
  //// y-axis is defined by neg lepton cross z-axis
  //// the subtle part is here...
  //////////TVector3 n_unity_1 = n_p4M11inXFrame_unit.Cross(n_unitz_1);
  TVector3 n_unity_1 = n_unitz_1.Cross(n_p4M11inXFrame_unit);
  TVector3 n_unitx_1 = n_unity_1.Cross(n_unitz_1);

  TLorentzVector n_p4M21inXFrame(p4M21);
  n_p4M21inXFrame.Boost(boostX);
  TVector3 n_p4M21inXFrame_unit = n_p4M21inXFrame.Vect().Unit();
  //rotate into other plane
  TVector3 n_p4M21inXFrame_unitprime(n_p4M21inXFrame_unit.Dot(n_unitx_1), n_p4M21inXFrame_unit.Dot(n_unity_1), n_p4M21inXFrame_unit.Dot(n_unitz_1));

  ///////-----------------new way of calculating phi-----------------///////
  //Float_t phi_n =  n_p4M21inXFrame_unitprime.Phi();
  //
    //std::cout << "---------------------------" << std::endl;
    //std::cout << "phi: " << phi << std::endl;
    //std::cout << "phi_n: " << phi_n << std::endl;
    //std::cout << "phi + phi_n: " << (phi+phi_n) << std::endl;
  //
  /// and then calculate phistar1
  TVector3 n_p4PartoninXFrame_unit(0.0, 0.0, 1.0);
  TVector3 n_p4PartoninXFrame_unitprime(n_p4PartoninXFrame_unit.Dot(n_unitx_1), n_p4PartoninXFrame_unit.Dot(n_unity_1), n_p4PartoninXFrame_unit.Dot(n_unitz_1));
  // negative sign is for arrow convention in paper
  phistar1 = (n_p4PartoninXFrame_unitprime.Phi());

  // and the calculate phistar2
  TLorentzVector n_p4Z2inXFrame(p4Z2);
  n_p4Z2inXFrame.Boost(boostX);
  TVector3 n_p4Z2inXFrame_unit = n_p4Z2inXFrame.Vect().Unit();
  ///////TLorentzVector n_p4M21inXFrame(p4M21);
  //////n_p4M21inXFrame.Boost(boostX);        
  ////TVector3 n_p4M21inXFrame_unit = n_p4M21inXFrame.Vect().Unit();  
  TVector3 n_unitz_2(n_p4Z2inXFrame_unit);
  //// y-axis is defined by neg lepton cross z-axis
  //// the subtle part is here...
  //////TVector3 n_unity_2 = n_p4M21inXFrame_unit.Cross(n_unitz_2);
  TVector3 n_unity_2 = n_unitz_2.Cross(n_p4M21inXFrame_unit);
  TVector3 n_unitx_2 = n_unity_2.Cross(n_unitz_2);
  TVector3 n_p4PartoninZ2PlaneFrame_unitprime(n_p4PartoninXFrame_unit.Dot(n_unitx_2), n_p4PartoninXFrame_unit.Dot(n_unity_2), n_p4PartoninXFrame_unit.Dot(n_unitz_2));
  phistar2 = (n_p4PartoninZ2PlaneFrame_unitprime.Phi());

  //
    //Float_t phistar12_0 = phistar1 + phistar2;
    //if (phistar12_0 > TMath::Pi()) phistar12 = phistar12_0 - 2*TMath::Pi();
    //else if (phistar12_0 < (-1.)*TMath::Pi()) phistar12 = phistar12_0 + 2*TMath::Pi();
    //else phistar12 = phistar12_0;
  //

}

//define this as a plug-in
DEFINE_FWK_MODULE(B2GTTbarTreeMaker);

//  LocalWords:  NNPDF
