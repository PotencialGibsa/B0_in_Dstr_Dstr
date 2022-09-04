//код для D* D0 с катами
// интервал по массе для второго D расширяем .. например 1.77 ..... 1.97 GeV
// каты на отлёт и косинус для второго D
// -*- C++ -*-
//
// Package:    XbFrame/Xb_frame
// Class:      Xb_frame
//
/**\class Xb_frame Xb_frame.cc XbFrame/Xb_frame/plugins/Xb_frame.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Sergey Polikarpov
//         Created:  Tue, 15 Aug 2017 01:04:12 GMT
//
//

// system include files
#include <memory>

/// framework
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/FWLite/interface/EventBase.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/Handle.h"
/// triggers
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
// #include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
/// tracks
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "DataFormats/Candidate/interface/ShallowCloneCandidate.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h" //including

/// muons
#include "DataFormats/MuonReco/interface/MuonChamberMatch.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

/// vertex fits
#include "RecoVertex/VertexTools/interface/VertexDistance.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/V0Candidate/interface/V0Candidate.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "RecoVertex/VertexPrimitives/interface/BasicSingleVertexState.h"
#include "RecoVertex/VertexPrimitives/interface/VertexState.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "DataFormats/PatCandidates/interface/PackedCandidate.h" // for miniAOD
// there are three track collections in miniAOD (in addition to the muon
// and electron collections), embedded into particleflow objects:
//   *** update this ***
// packedPFCandidates allows to rebuild tracks (for pt>0.95/0.4 GeV) using
//                    pseudotrack() (object) or besttrack() (pointer)
// PackedPFCandidatesDiscarded presumably contains only discarded duplicate
//                    muon candidates -> do not use
// lostTracks (high purity only) might contain many of the non-vertex tracks
//                    needed for the slow pion measurement
#include "DataFormats/TrackReco/interface/Track.h" // for miniAOD

//// gen ??
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/CLHEP/interface/Migration.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Math/interface/Error.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"

#include <vector>
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "TLorentzVector.h"
#include <utility>
#include <string>
#include <map>


using namespace edm;
using namespace std;
using namespace reco;
//
///////
///// UPDATED TO 2018 PDG PARTICLES !!
/// {{{
double PDG_MUON_MASS    =   0.1056583745;
double PDG_PION_MASS    =   0.13957061;
double PDG_PIOZ_MASS    =   0.1349770;
double PDG_KAON_MASS    =   0.493677;
double PDG_PROTON_MASS  =   0.9382720813;
double PDG_KSHORT_MASS  =   0.497611;
double PDG_KSHORT_DM    =   0.000013;
double PDG_KSHORT_TIME  =   0.8954 * 0.0000000001;
double PDG_KS_MASS      =   PDG_KSHORT_MASS;
double PDG_LAMBDA_MASS  =   1.115683 ;
double PDG_LAMBDA_DM    =   0.000006 ;
double PDG_LAMBDA_TIME  =   2.632 * 0.0000000001;
double PDG_SIGMA0_MASS  =   1.192642 ;
double PDG_XImunus_MASS =   1.32171 ;
double PDG_XImunus_DM   =   0.00007 ;
double PDG_XImunus_TIME =   1.639 * 0.0000000001;
double PDG_OMmunus_MASS =   1.67245 ;
double PDG_OMmunus_DM   =   0.00029 ;
double PDG_OMmunus_TIME =   0.821 * 0.0000000001;
double PDG_DPM_MASS     =   1.86965 ;
double PDG_DPM_DM       =   0.00005 ;
double PDG_DPM_TIME     =   1.040 * 0.000000000001 ;
double PDG_DZ_MASS      =   1.86483 ;
double PDG_DZ_DM        =   0.00005 ;
double PDG_DZ_TIME      =   0.4101 * 0.000000000001 ;
double PDG_DS_MASS      =   1.96834 ;
double PDG_DS_DM        =   0.00007 ;
double PDG_DS_TIME      =   0.504 * 0.000000000001 ;
double PDG_LAMCZ_MASS   =   2.28646 ;
double PDG_LAMCZ_DM     =   0.00031 ;
double PDG_LAMCZ_TIME   =   2.00 * 0.0000000000001;
double PDG_XICZ_MASS    =   2.47087 ;
double PDG_XICZ_DM      =   0.00031 ;
double PDG_XICZ_TIME    =   1.12 * 0.0000000000001;
double PDG_XICP_MASS    =   2.46787 ;
double PDG_XICP_DM      =   0.00030 ;
double PDG_XICP_TIME    =   4.42 * 0.0000000000001;
double PDG_KSTARZ_MASS  =   0.89555;
double PDG_KSTARZ_GAMMA =   0.0473 ;
double PDG_KSTARP_MASS  =   0.89176;
double PDG_KSTARP_GAMMA =   0.0503 ;
double PDG_PHI_MASS     =   1.019461;
double PDG_PHI_GAMMA    =   0.004249;
double PDG_JPSI_MASS    =   3.096900;
double PDG_PSI2S_MASS   =   3.686097;
double PDG_X3872_MASS   =   3.87169;
double PDG_BU_MASS      =   5.27932;
double PDG_BU_TIME      =   1.638 * 0.000000000001;
double PDG_B0_MASS      =   5.27963;
double PDG_B0_TIME      =   1.520 * 0.000000000001;
double PDG_BS_MASS      =   5.36689;
double PDG_BS_TIME      =   1.509 * 0.000000000001;
double PDG_BC_MASS      =   6.2749;
double PDG_BC_TIME      =   0.507 * 0.000000000001;
double PDG_LB_MASS      =   5.61960;
double PDG_LB_TIME      =   1.470 * 0.000000000001;
double PDG_XIBZ_MASS    =   5.7919;
double PDG_XIBZ_TIME    =   1.479 * 0.000000000001;
double PDG_XIBM_MASS    =   5.7970;
double PDG_XIBM_TIME    =   1.571 * 0.000000000001;
double PDG_OMBM_MASS    =   6.0461;
double PDG_OMBM_TIME    =   1.64 * 0.000000000001;
double PDG_C            =   29979245800.; // in cm/c
//
/// }}}

ParticleMass PM_PDG_MUON_MASS = PDG_MUON_MASS;
ParticleMass PM_PDG_JPSI_MASS = PDG_JPSI_MASS;
ParticleMass PM_PDG_KAON_MASS = PDG_KAON_MASS;
ParticleMass PM_PDG_PION_MASS = PDG_PION_MASS;
ParticleMass PM_PDG_PROTON_MASS = PDG_PROTON_MASS;
ParticleMass PM_PDG_PSI2S_MASS = PDG_PSI2S_MASS;

//
// class declaration
//
class Xb_frame : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit Xb_frame(const edm::ParameterSet&);
      ~Xb_frame();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
// HLTConfigProvider hltConfig_;
// edm::EDGetTokenT<edm::TriggerResults> hlTriggerResults_;
edm::EDGetTokenT<reco::BeamSpot> thebeamspot_;
edm::EDGetTokenT<reco::VertexCollection> vtxSample;

// edm::EDGetTokenT<vector < pat::GenericParticle > > tracks_; /// AOD OLNY
edm::EDGetTokenT<pat::PackedCandidateCollection> trkTkn; /// miniAOD
edm::EDGetTokenT<pat::PackedCandidateCollection> trkTkndisc; /// miniAOD
edm::EDGetTokenT<pat::PackedCandidateCollection> trkTknlost; /// miniAOD

edm::EDGetTokenT<reco::VertexCompositePtrCandidateCollection> tok_v0_; /// miniAOD

edm::EDGetTokenT<edm::View<reco::GenParticle>> prunedGenToken_; /// miniAOD
edm::EDGetTokenT<edm::View<pat::PackedGenParticle> > packedGenToken_; /// miniAOD
edm::EDGetTokenT<edm::TriggerResults> triggerToken_;

std::vector<float> *B_mass          , *B_masd         ;
std::vector<float> *B_px            , *B_py           , *B_pz;
std::vector<float> *B_DecayVtxX     , *B_DecayVtxY    , *B_DecayVtxZ;
std::vector<float> *B_DecayVtxXE    , *B_DecayVtxYE   , *B_DecayVtxZE;
std::vector<float> *B_Prob          ;


std::vector<float> *D_mass,           *D_mass_CV;
std::vector<float> *D_px,             *D_py,            *D_pz ;
std::vector<float> *D_px_CV,          *D_py_CV,         *D_pz_CV ;
std::vector<float> *D_DecayVtxX,      *D_DecayVtxY,     *D_DecayVtxZ ;
std::vector<float> *D_DecayVtxXE,     *D_DecayVtxYE,    *D_DecayVtxZE ;
std::vector<float> *D_Prob ;


std::vector<float> *KA1_px    , *KA1_py     , *KA1_pz;
std::vector<float> *KA1_px_CV , *KA1_py_CV  , *KA1_pz_CV;
std::vector<float> *KA1_ips;
std::vector<int> *KA1_charge;

std::vector<float> *PI2_px    , *PI2_py     , *PI2_pz;
std::vector<float> *PI2_dr    , *PI2_dz;
std::vector<float> *PI2_px_CV , *PI2_py_CV  , *PI2_pz_CV;
std::vector<float> *PI2_ips;


std::vector<int>   *PIS1_charg , *PIS1_purit;
std::vector<float> *PIS1_px    , *PIS1_py     , *PIS1_pz;
std::vector<float> *PIS1_dr    , *PIS1_dz;
std::vector<float> *PIS1_px_CV , *PIS1_py_CV  , *PIS1_pz_CV;
std::vector<float> *PIS1_ips;

std::vector<int>   *PIS2_charg , *PIS2_purit;
std::vector<float> *PIS2_px    , *PIS2_py     , *PIS2_pz;
std::vector<float> *PIS2_dr    , *PIS2_dz;
std::vector<float> *PIS2_px_CV , *PIS2_py_CV  , *PIS2_pz_CV;
std::vector<float> *PIS2_ips;


std::vector<float> *E_mass,           *E_mass_CV;
std::vector<float> *E_px,             *E_py,            *E_pz ;
std::vector<float> *E_px_CV,          *E_py_CV,         *E_pz_CV ;
std::vector<float> *E_DecayVtxX,      *E_DecayVtxY,     *E_DecayVtxZ ;
std::vector<float> *E_DecayVtxXE,     *E_DecayVtxYE,    *E_DecayVtxZE ;
std::vector<float> *E_Prob ;


std::vector<float> *KA3_px    , *KA3_py     , *KA3_pz;
std::vector<float> *KA3_px_CV , *KA3_py_CV  , *KA3_pz_CV;
std::vector<float> *KA3_ips;
std::vector<int> *KA3_charge;

std::vector<float> *PI4_px    , *PI4_py     , *PI4_pz;
std::vector<float> *PI4_dr    , *PI4_dz;
std::vector<float> *PI4_px_CV , *PI4_py_CV  , *PI4_pz_CV;
std::vector<float> *PI4_ips;

std::vector<float> *PV_becos_XX , *PV_becos_YY  , *PV_becos_ZZ;
std::vector<float> *PV_becos_EX , *PV_becos_EY  , *PV_becos_EZ;
std::vector<float> *PV_becos_CL;
std::vector<int>   *PV_becos_dN;

std::vector<int>   *GEN_dectype;
std::vector<int>   *GEN_momID;
std::vector<float> *GEN_PIS_px  , *GEN_PIS_py   , *GEN_PIS_pz;
std::vector<float> *GEN_KS_px   , *GEN_KS_py    , *GEN_KS_pz;
std::vector<float> *GEN_KR_px   , *GEN_KR_py    , *GEN_KR_pz;
std::vector<float> *GEN_PI1_px  , *GEN_PI1_py   , *GEN_PI1_pz;
std::vector<float> *GEN_PI2_px  , *GEN_PI2_py   , *GEN_PI2_pz;

std::vector<bool>  *trig0_fi;
std::vector<bool>  *trig1_fi;
std::vector<bool>  *trig2_fi;
std::vector<bool>  *trig3_fi;
std::vector<bool>  *trig4_fi;
std::vector<bool>  *trig5_fi;
std::vector<bool>  *trig6_fi;
std::vector<bool>  *trig7_fi;
std::vector<bool>  *trig8_fi;
std::vector<bool>  *trig9_fi;

Int_t       nCand;

Int_t       run;
Int_t       event;
Float_t     lumi;

Float_t     BESP_x ;
Float_t     BESP_y ;
Float_t     BESP_z ;
Float_t     BESP_ex;
Float_t     BESP_ey;
Float_t     BESP_ez;

Int_t       numPV;
Int_t       numTrack;
Int_t       numV0;

TTree *wwtree;
TFile *f;
std::string fileName;

};

//
Xb_frame::Xb_frame(const edm::ParameterSet& iConfig) :

B_mass(0),          B_masd(0),
B_px(0),            B_py(0),            B_pz(0),
B_DecayVtxX(0),     B_DecayVtxY(0),     B_DecayVtxZ(0),
B_DecayVtxXE(0),    B_DecayVtxYE(0),    B_DecayVtxZE(0),
B_Prob(0),

D_mass(0),           D_mass_CV(0),
D_px(0),             D_py(0),            D_pz(0),
D_px_CV(0),          D_py_CV(0),         D_pz_CV(0),
D_DecayVtxX(0),      D_DecayVtxY(0),     D_DecayVtxZ(0),
D_DecayVtxXE(0),     D_DecayVtxYE(0),    D_DecayVtxZE(0),
D_Prob(0),

KA1_px(0)     , KA1_py(0)   , KA1_pz(0),
KA1_px_CV(0)  , KA1_py_CV(0), KA1_pz_CV(0),
KA1_ips(0),KA1_charge(0),

PI2_px(0)     , PI2_py(0)   , PI2_pz(0),
PI2_dr(0)     , PI2_dz(0) ,
PI2_px_CV(0)  , PI2_py_CV(0), PI2_pz_CV(0),
PI2_ips(0),


PIS1_charg(0)  , PIS1_purit(0),
PIS1_px(0)     , PIS1_py(0)   , PIS1_pz(0),
PIS1_dr(0)     , PIS1_dz(0) ,
PIS1_px_CV(0)  , PIS1_py_CV(0), PIS1_pz_CV(0),
PIS1_ips(0),

PIS2_charg(0)  , PIS2_purit(0),
PIS2_px(0)     , PIS2_py(0)   , PIS2_pz(0),
PIS2_dr(0)     , PIS2_dz(0) ,
PIS2_px_CV(0)  , PIS2_py_CV(0), PIS2_pz_CV(0),
PIS2_ips(0),

E_mass(0),           E_mass_CV(0),
E_px(0),             E_py(0),            E_pz(0),
E_px_CV(0),          E_py_CV(0),         E_pz_CV(0),
E_DecayVtxX(0),      E_DecayVtxY(0),     E_DecayVtxZ(0),
E_DecayVtxXE(0),     E_DecayVtxYE(0),    E_DecayVtxZE(0),
E_Prob(0),

KA3_px(0)     , KA3_py(0)   , KA3_pz(0),
KA3_px_CV(0)  , KA3_py_CV(0), KA3_pz_CV(0),
KA3_ips(0), KA3_charge(0),

PI4_px(0)     , PI4_py(0)   , PI4_pz(0),
PI4_dr(0)     , PI4_dz(0) ,
PI4_px_CV(0)  , PI4_py_CV(0), PI4_pz_CV(0),
PI4_ips(0),

PV_becos_XX(0)  , PV_becos_YY(0), PV_becos_ZZ(0),
PV_becos_EX(0)  , PV_becos_EY(0), PV_becos_EZ(0),
PV_becos_CL(0)  , PV_becos_dN(0),
//

GEN_dectype(0),
GEN_momID(0),
GEN_PIS_px(0)   , GEN_PIS_py(0) , GEN_PIS_pz(0),
GEN_KS_px(0)    , GEN_KS_py(0)  , GEN_KS_pz(0) ,
GEN_KR_px(0)    , GEN_KR_py(0)  , GEN_KR_pz(0) ,
GEN_PI1_px(0)   , GEN_PI1_py(0) , GEN_PI1_pz(0),
GEN_PI2_px(0)   , GEN_PI2_py(0) , GEN_PI2_pz(0),

trig0_fi(0),    trig1_fi(0),    trig2_fi(0),
trig3_fi(0),    trig4_fi(0),    trig5_fi(0),
trig6_fi(0),    trig7_fi(0),    trig8_fi(0),
trig9_fi(0),
//
nCand(0),

run(0),
event(0),
lumi(0),

BESP_x(0),
BESP_y(0),
BESP_z(0),
BESP_ex(0),
BESP_ey(0),
BESP_ez(0),

numPV(0),
numTrack(0),
numV0(0)

{
// fileName = iConfig.getUntrackedParameter<std::string>("fileName","BFinder.root");
fileName = iConfig.getUntrackedParameter<std::string>("fileName");
   //now do what ever initialization is needed
   usesResource("TFileService");
// tok_v0_     = consumes<reco::VertexCompositeCandidateCollection>(   edm::InputTag("generalV0Candidates:Lambda")); /// AOD

vtxSample       =   consumes<reco::VertexCollection>(edm::InputTag("offlineSlimmedPrimaryVertices")); /// miniAOD
thebeamspot_    =   consumes<reco::BeamSpot>(edm::InputTag("offlineBeamSpot"));  /// miniAOD
trkTkn          =   consumes<pat::PackedCandidateCollection>(edm::InputTag("packedPFCandidates")); // miniAOD
trkTkndisc      =   consumes<pat::PackedCandidateCollection>(edm::InputTag("packedPFCandidatesDiscarded")); /// miniAOD
trkTknlost      =   consumes<pat::PackedCandidateCollection>(edm::InputTag("lostTracks")); /// miniAOD
tok_v0_         =   consumes<reco::VertexCompositePtrCandidateCollection>(edm::InputTag("slimmedKshortVertices")); /// miniAOD

prunedGenToken_ =   consumes<edm::View<reco::GenParticle>> (edm::InputTag("prunedGenParticles"));
packedGenToken_ =   consumes<edm::View<pat::PackedGenParticle>>(edm::InputTag("packedGenParticles"));
triggerToken_   =   consumes<edm::TriggerResults>(edm::InputTag("TriggerResults", "", "HLT"));


}


Xb_frame::~Xb_frame()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
Xb_frame::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
using namespace edm;
using namespace reco;
using namespace std;
using reco::MuonCollection;

run   = iEvent.id().run();
event = iEvent.id().event();

lumi = 1;
lumi = iEvent.luminosityBlock();

//// {{{ GEN INFO// GEN INFO// GEN INFO// GEN INFO// GEN INFO  // GEN INFO
Handle<edm::View<reco::GenParticle>> pruned;
iEvent.getByToken(prunedGenToken_, pruned);
// std::cout << run << ":" << event << " val p  "<< pruned.isValid()  << std::endl;
Handle<edm::View<pat::PackedGenParticle>> packed;
iEvent.getByToken(packedGenToken_, packed);
// std::cout << run << ":" << event << " val pa "<< packed.isValid()  << std::endl;
Int_t genNumDzDau = 0;
if ( pruned.isValid() )
{
    int _numgenpr = pruned->size() ;
    // int _numgenpa = packed->size() ;

    // std::cout << run << ":" << event << " prun "<< _numgenpr << std::endl;
    // D*+ = 413, D0=421, KS=310, pi+=211 // B0 = 511, B+ = 521

        for(Int_t ii=0; ii< _numgenpr; ii++)
        {
            Int_t pdgId = (*pruned)[ii].pdgId();
            Int_t mysign = pdgId > 0 ? 1 : -1;
            if (abs(pdgId) == 413) // look for D*
            {
                const Candidate *genDST = &(*pruned)[ii];   //  D*
                const Candidate *genMom = genDST->mother(); //  its mother
                //
                // cout << 1 << std::endl;
                if (genDST->numberOfDaughters() < 2) continue;
                const Candidate *genDz  = genDST->daughter(0);  // Dst decays to D0(this) pi
                const Candidate *genPis = genDST->daughter(1);  // Dst decays to D0 pi(this)
                if (abs(genDz->pdgId()) != 421) continue;       // checking their IDS
                if (abs(genPis->pdgId()) != 211) continue;
                //
                // cout << 2 << std::endl;
                Int_t genGrMomId = 0;
                if (genMom->numberOfMothers() > 0)
                {
                    genGrMomId = genMom->mother()->pdgId();
                }
                if (genGrMomId == 211 ) {}; // just a silly check to not get "unused variable" error
                // cout << 3 << std::endl;
                ///
                genNumDzDau = genDz->numberOfDaughters() ;
                //std::cout << " dst m "<< genDST->mass()  << " pt "<< genDST->pt() << " Mom " << genMom->pdgId()  << " grMom " << genGrMomId << " dau0 " << genDz->pt() << " dau1 " << genPis->pt() << "nD0dau " << genNumDzDau << std::endl;
                if (genNumDzDau < 2) continue;
                //
                // removing extra QED photons here
                if (genNumDzDau == 5)
                {
                    if (genDz->daughter(4)->pdgId() == 22) genNumDzDau = 4;
                }
                if (genNumDzDau == 4)
                {
                    if (genDz->daughter(3)->pdgId() == 22) genNumDzDau = 3;
                }
                if (genNumDzDau == 3)
                {
                    if (genDz->daughter(2)->pdgId() == 22) genNumDzDau = 2;
                }
                // finished removing extra photons here
                // now understand which decay is there
                //
                if (genNumDzDau == 2)
                {
                    const Candidate *genH1 = genDz->daughter(0);    // D0 daughters
                    const Candidate *genH2 = genDz->daughter(1);    // D0 daughters
                    if ((abs(genH1->pdgId()) == 211) && (abs(genH2->pdgId()) == 211))
                    {
                        // HERE WE ARE D0 -> PI PI
                        // for every event fill all 'branches' and then figure out offline which decay and which branches to read (others zre 0)
                        GEN_dectype ->push_back(3*mysign);
                        GEN_momID   ->push_back(genMom->pdgId() );
                        GEN_PIS_px  ->push_back(genPis->px()    );
                        GEN_PIS_py  ->push_back(genPis->py()    );
                        GEN_PIS_pz  ->push_back(genPis->pz()    );
                        GEN_KS_px   ->push_back(0   );
                        GEN_KS_py   ->push_back(0   );
                        GEN_KS_pz   ->push_back(0   );
                        GEN_KR_px   ->push_back(0   );
                        GEN_KR_py   ->push_back(0   );
                        GEN_KR_pz   ->push_back(0   );
                        GEN_PI1_px  ->push_back(genH1->px()     );
                        GEN_PI1_py  ->push_back(genH1->py()     );
                        GEN_PI1_pz  ->push_back(genH1->pz()     );
                        GEN_PI2_px  ->push_back(genH1->px()     );
                        GEN_PI2_py  ->push_back(genH1->py()     );
                        GEN_PI2_pz  ->push_back(genH1->pz()     );
                    }
                    //
                    if ((abs(genH1->pdgId()) == 310) && (abs(genH2->pdgId()) == 310))
                    {
                        // HERE WE ARE D0 -> KS KS
                        GEN_dectype ->push_back(2*mysign);
                        GEN_momID   ->push_back(genMom->pdgId() );
                        GEN_PIS_px  ->push_back(genPis->px()    );
                        GEN_PIS_py  ->push_back(genPis->py()    );
                        GEN_PIS_pz  ->push_back(genPis->pz()    );
                        GEN_KS_px   ->push_back(genH1->px()     );
                        GEN_KS_py   ->push_back(genH1->py()     );
                        GEN_KS_pz   ->push_back(genH1->pz()     );
                        GEN_KR_px   ->push_back(genH2->px()     );
                        GEN_KR_py   ->push_back(genH2->py()     );
                        GEN_KR_pz   ->push_back(genH2->pz()     );
                        GEN_PI1_px  ->push_back(0   );
                        GEN_PI1_py  ->push_back(0   );
                        GEN_PI1_pz  ->push_back(0   );
                        GEN_PI2_px  ->push_back(0   );
                        GEN_PI2_py  ->push_back(0   );
                        GEN_PI2_pz  ->push_back(0   );
                    }
                    //
                    if ((abs(genH1->pdgId()) == 333) && (abs(genH2->pdgId()) == 310))
                    {
                        // HERE WE ARE D0 -> PHI KS
                        if (genH1->numberOfDaughters() < 2) continue;
                        const Candidate *genK1 = genH1->daughter(0);    // D0 daughters
                        const Candidate *genK2 = genH1->daughter(1);    // D0 daughters
                        if (abs(genK1->pdgId()) != 321) continue;
                        if (abs(genK2->pdgId()) != 321) continue;
                        //
                        GEN_dectype ->push_back(4*mysign);
                        GEN_momID   ->push_back(genMom->pdgId() );
                        GEN_PIS_px  ->push_back(genPis->px()    );
                        GEN_PIS_py  ->push_back(genPis->py()    );
                        GEN_PIS_pz  ->push_back(genPis->pz()    );
                        GEN_KS_px   ->push_back(genH2->px()     );
                        GEN_KS_py   ->push_back(genH2->py()     );
                        GEN_KS_pz   ->push_back(genH2->pz()     );
                        GEN_KR_px   ->push_back(0   );
                        GEN_KR_py   ->push_back(0   );
                        GEN_KR_pz   ->push_back(0   );
                        GEN_PI1_px  ->push_back(genK1->px()     );
                        GEN_PI1_py  ->push_back(genK1->py()     ); // HERE saving into pions but they are actually KAONS !! (we know it since GEN_dectype = +-4)
                        GEN_PI1_pz  ->push_back(genK1->pz()     );
                        GEN_PI2_px  ->push_back(genK2->px()     );
                        GEN_PI2_py  ->push_back(genK2->py()     );
                        GEN_PI2_pz  ->push_back(genK2->pz()     );
                    }
                    // if ((abs(genH1->pdgId()) == 310) && (abs(genH2->pdgId()) == 333))
                    // {// HERE WE ARE D0 -> KS PHI ---- THIS NEVER HAPPENS, ALWAYS FIRST PHI THEN KS
                    //    std::cout<<" 2 dz daus!!! " << genDz->daughter(0)->pdgId() << " "<< genDz->daughter(1)->pdgId() << std::endl<< std::endl<< std::endl;
                    // }
                }
                if (genNumDzDau == 3)
                {
                    const Candidate *genKS = genDz->daughter(0);    // D0 daughters in case of Ks PI PI the 1st one is always KS !!!
                    const Candidate *genP1 = genDz->daughter(1);    // D0 daughters
                    const Candidate *genP2 = genDz->daughter(2);    // D0 daughters
                    //
                    if (abs(genKS->pdgId()) != 310) continue;
                    if (abs(genP1->pdgId()) != 211) continue;
                    if (abs(genP2->pdgId()) != 211) continue;
                    //
                    // HERE WE ARE D0 -> KS PI PI
                    //
                    GEN_dectype ->push_back(1*mysign);
                    GEN_momID   ->push_back(genMom->pdgId() );
                    GEN_PIS_px  ->push_back(genPis->px()    );
                    GEN_PIS_py  ->push_back(genPis->py()    );
                    GEN_PIS_pz  ->push_back(genPis->pz()    );
                    GEN_KS_px   ->push_back(genKS->px()     );
                    GEN_KS_py   ->push_back(genKS->py()     );
                    GEN_KS_pz   ->push_back(genKS->pz()     );
                    GEN_KR_px   ->push_back(0   );
                    GEN_KR_py   ->push_back(0   );
                    GEN_KR_pz   ->push_back(0   );
                    GEN_PI1_px  ->push_back(genP1->px()     );
                    GEN_PI1_py  ->push_back(genP1->py()     );
                    GEN_PI1_pz  ->push_back(genP1->pz()     );
                    GEN_PI2_px  ->push_back(genP2->px()     );
                    GEN_PI2_py  ->push_back(genP2->py()     );
                    GEN_PI2_pz  ->push_back(genP2->pz()     );
                }
            }
        }
}   else
{
    // std::cout << run << ":" << event << " THIS IS DATA " << std::endl;
}
// if (GEN_dectype->size() ==0)
// {
//     std::cout << run << ":" << event << " num of gen-saved " << GEN_dectype->size() << " aa " << genNumDzDau << std::endl;
// }
// Decay MyD0
// 0.60  K_S0  pi+  pi-    PHSP;    310 211 211             DECAYTYPE = 1
// 0.30  K_S0  K_S0        PHSP;    310 310                 DECAYTYPE = 2
// 0.05  pi+   pi-         PHSP;    211 211                 DECAYTYPE = 3
// 0.05  MyPhi K_S0        SVS;     333 310 [333 -> 321 321]DECAYTYPE = 4

//// }}} GEN INFO // GEN INFO // GEN INFO // GEN INFO // GEN INFO // GEN INFO

//// {{{ TRIGGER

Handle<edm::TriggerResults> triggerResults;
iEvent.getByToken(triggerToken_, triggerResults);

const edm::TriggerNames &TrigNames = iEvent.triggerNames(*triggerResults);
// std::cout << run << ":" << event << " == TRIGGER PATHS= " << std::endl;

unsigned int NTRIGGERS = 10;
std::string TriggersToTest[NTRIGGERS] = {
"HLT_Mu12_IP6_p",       // 0
"HLT_Mu10p5_IP3p5_p",   // 1
"HLT_Mu9_IP6_p",        // 2
"HLT_Mu9_IP5_p",        // 3
"HLT_Mu9_IP4_p",        // 4
"HLT_Mu8p5_IP3p5_p",    // 5
"HLT_Mu8_IP6_p",        // 6
"HLT_Mu8_IP5_p",        // 7
"HLT_Mu8_IP3_p",        // 8
"HLT_Mu7_IP4_p",        // 9
};
unsigned short int TriggersFired[NTRIGGERS] = {0,0,0,0,0,0,0,0,0,0};

/// THIS WORKS !!!
bool FoundParkedTrigger = false;
for (unsigned int i = 0, n = triggerResults->size(); i < n; ++i)
{
    if (!triggerResults->accept(i)) continue;
    string hltName = TrigNames.triggerName(i);
    // std::cout << " Fired Path : " << hltName << std::endl;
    //
    // if( hltName.find("HLT_Mu7_IP") != string::npos || (hltName.find("HLT_Mu8_IP") != string::npos) || (hltName.find("HLT_Mu8p5_IP") != string::npos) || (hltName.find("HLT_Mu9_IP") != string::npos) || (hltName.find("HLT_Mu10p5_IP") != string::npos) || (hltName.find("HLT_Mu12_IP") != string::npos)  )
    // {
        // std::cout << " Fired Parking Path : " << hltName << std::endl;
        // FoundParkedTrigger = true;
    // }
    for (unsigned int k = 0; k < NTRIGGERS; k++)
    {
        if( hltName.find(TriggersToTest[k]) != string::npos )
        {
            // std::cout << " Fired Parking trigger " << k << std::endl;
            TriggersFired[k] = 1;
            FoundParkedTrigger = true;
        }
    }
}
if (!FoundParkedTrigger)
{
    /// this happens very often in MC but never in DATA
    // std::cout << run << ":" << event << "NO PARKING TRIGGER FIRED AAAAAAAAAAAAA" << std::endl << std::endl << std::endl << std::endl;
}
//// }}} TRIGGER

// declare new track builder for new Transient track collection for miniAOD ???
ESHandle<TransientTrackBuilder> theB;
iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);

//// BEAMSPOT
Handle<reco::BeamSpot> beamSpotHandle;
iEvent.getByToken(thebeamspot_, beamSpotHandle);
reco::BeamSpot vertexBeamSpot= *beamSpotHandle;
BESP_x = vertexBeamSpot.x0();
BESP_y = vertexBeamSpot.y0();
BESP_z = vertexBeamSpot.z0();
BESP_ex = vertexBeamSpot.BeamWidthX();
BESP_ey = vertexBeamSpot.BeamWidthY();
BESP_ez = vertexBeamSpot.sigmaZ();

//// PV PV  PV  PV  PV  PV  PV  PV  PV  PV  PV  PV  PV
Handle < VertexCollection > recVtxs;
iEvent.getByToken(vtxSample, recVtxs);
///
Vertex thePrimaryV;
thePrimaryV = Vertex(*(recVtxs->begin()));
const reco::VertexCollection & vertices = *recVtxs.product();

// Handle < vector < pat::GenericParticle > >thePATTrackHandle; /// AOD ONLY
// iEvent.getByToken(tracks_, thePATTrackHandle); /// AOD ONLY
Handle<pat::PackedCandidateCollection> tracks; /// miniAOD
Handle<pat::PackedCandidateCollection> discTracks; /// miniAOD
Handle<pat::PackedCandidateCollection> lostTracks; /// miniAOD
iEvent.getByToken(trkTkn, tracks); /// miniAOD
iEvent.getByToken(trkTkndisc, discTracks); /// miniAOD
iEvent.getByToken(trkTknlost, lostTracks); /// miniAOD

////  V0  V0  V0  V0  V0  V0  V0  V0  V0  V0  V0  V0
Handle<reco::VertexCompositePtrCandidateCollection> v0Coll; // miniAOD
iEvent.getByToken(tok_v0_, v0Coll);

numPV       = vertices.size();
numTrack    = tracks->size();
numV0       = v0Coll->size();

// if (run   < 316114 ) return;
// if (run   > 316114 ) return;  // for testing specific event
// if (event < 755293018 ) return;
// if (event > 755293018 ) return;

// if (numV0   < 1 ) return;
if (numTrack< 4 ) return;

// cout << "numV0 " << numV0 << " numTrack " << numTrack << std::endl;

nCand = 0;

float PM_sigma = 1.e-7;
KinematicParticleFactoryFromTransientTrack pFactory;
//
//
float chi = 0.;
float ndf = 0.;
int fitgood = 1;
//
/// Find first Kaon for first D0
//
for (pat::PackedCandidateCollection::const_iterator iTra1 = tracks->begin(); iTra1 != tracks->end(); iTra1++)
{
	if (!(iTra1->hasTrackDetails())) continue;
        auto iTrack1 = iTra1->bestTrack();
        if (iTrack1 == nullptr) continue;
        if (iTrack1->pt()<0.5) continue;
        if (!(iTrack1->quality(reco::TrackBase::highPurity))) continue;

        TLorentzVector p4ka1;
        p4ka1.SetPtEtaPhiM(iTrack1->pt(),iTrack1->eta(),iTrack1->phi(), PDG_KAON_MASS);

        if( fabs(iTrack1->eta()) > 2.6 ) continue;
        TransientTrack kaon1TT = theB->build( iTra1->pseudoTrack() );
        if (!kaon1TT.isValid()) continue;

        // Find first Pion for first D0
        for (pat::PackedCandidateCollection::const_iterator iTra2 = tracks->begin(); iTra2 != tracks->end(); iTra2++)
        {
            if (!(iTra2->hasTrackDetails())) continue;
            auto iTrack2 = iTra2->bestTrack();
            if (iTrack2 == nullptr) continue;

            if (iTrack2->pt()<0.50) continue;
            if (iTrack1->charge() * iTrack2->charge() > -0.5) continue;
            if (!(iTrack2->quality(reco::TrackBase::highPurity))) continue;

            if((iTrack2->pt()<0.6)&&(iTrack1->pt()<0.6)) continue;

            TLorentzVector p4pi2;
            p4pi2.SetPtEtaPhiM(iTrack2->pt(),iTrack2->eta(),iTrack2->phi(), PDG_PION_MASS);
            if ((p4ka1 + p4pi2).M() > 1.99) continue;
            if ((p4ka1 + p4pi2).M() < 1.76) continue;
            if ((p4ka1 + p4pi2).Pt() < 3.9) continue;

            if(iTrack1 == iTrack2) continue;

            if( fabs(iTrack2->eta()) > 2.6 ) continue;

            TransientTrack pion2TT = theB->build( iTra2->pseudoTrack() );
            if (!pion2TT.isValid()) continue;

            double pi2dz   = iTrack2->vz() - iTrack1->vz();
            double pi2dr   = sqrt( (iTrack2->vx() - iTrack1->vx())*(iTrack2->vx() - iTrack1->vx()) + (iTrack2->vy() - iTrack1->vy())*(iTrack2->vy() - iTrack1->vy()) );

            if (pi2dz > 3.0) continue;
            if (pi2dz < -3.0) continue;


            ///////////////////////////////  FIRST D0 VTX TMP /////////////////////////////////////////

            std::vector<RefCountedKinematicParticle> Q_candidate_init; //

            Q_candidate_init.push_back(pFactory.particle(kaon1TT, PM_PDG_KAON_MASS, chi,ndf, PM_sigma));
            Q_candidate_init.push_back(pFactory.particle(pion2TT, PM_PDG_PION_MASS, chi,ndf, PM_sigma));
            RefCountedKinematicTree qertexFitTree;


            std::vector<RefCountedKinematicParticle> Q_candidate = Q_candidate_init;
            KinematicParticleVertexFitter qFitter;
            fitgood = 1;
            try
            {
                qertexFitTree = qFitter.fit(Q_candidate);
            }
            catch (const VertexException &)
            {
                fitgood = 0;
            }
            if (fitgood == 0) continue;

            if (!qertexFitTree->isValid()) continue;
            qertexFitTree->movePointerToTheTop();

            RefCountedKinematicParticle qCand         = qertexFitTree->currentParticle();
            RefCountedKinematicVertex   qDecayVertex  = qertexFitTree->currentDecayVertex();
            double Q_Prob_tmp   = 0;
            fitgood = 1;
            try
            {
                if(qDecayVertex->chiSquared() < 0) fitgood = 0;
                Q_Prob_tmp   = TMath::Prob(qDecayVertex->chiSquared(), (int) qDecayVertex->degreesOfFreedom());
            }
            catch (const Exception &)
            {
                fitgood = 0;
            }
            if (fitgood == 0) continue;

            if(Q_Prob_tmp < 0.05) continue;
            if (!qDecayVertex->vertexIsValid())  continue;

            double Q_mass_tmp = qCand->currentState().mass();
            if(Q_mass_tmp < 1.82) continue;
            if(Q_mass_tmp > 1.91) continue;
            //
            ////
            reco::Vertex bestVtxBSIP;  /// using D0 fit and D0 vtx HERE !!
             // {{{ GET THE BEST PV BY CHOSING THE SMALLEST ANGLE
              // ********************* ****

            reco::Vertex vtxBSrf ;

            Double_t pVtxBSIPX_temp = -10000.0;
            Double_t pVtxBSIPY_temp = -10000.0;
            Double_t pVtxBSIPZ_temp = -10000.0;
            Double_t pVtxBSIPXE_temp = -10000.0;
            Double_t pVtxBSIPYE_temp = -10000.0;
            Double_t pVtxBSIPZE_temp = -10000.0;
            Double_t pVtxBSIPCL_temp = -10000.0;
            Double_t pVtxBSIPdN_temp = 0;
            Double_t lip = -100000.0;
            // Double_t ddd = -100000.0;
            //

            for(size_t i = 0; i < recVtxs->size(); ++i)
            {
                Double_t ptsum_ = 0;
                const Vertex &vtxBS = (*recVtxs)[i];
                vector<reco::TransientTrack> vertexTracks;
                for ( std::vector<TrackBaseRef >::const_iterator iTrack = vtxBS.tracks_begin(); iTrack != vtxBS.tracks_end(); ++iTrack)
                {
                    TrackRef trackRef = iTrack->castTo<TrackRef>();

                    if (  !( (iTrack1  ==trackRef.get() )
                           ||(iTrack2  ==trackRef.get() )) )
                    {

                        TransientTrack tt = theB->build( trackRef );
                        vertexTracks.push_back(tt);
                        ptsum_ += trackRef->pt();
                    }
                }

                // if no tracks in primary or no reco track included in primary then don't do anything
                vtxBSrf = vtxBS;
                GlobalPoint PVRfP = GlobalPoint( vtxBS.x(), vtxBS.y(), vtxBS.z() );
                if ( vertexTracks.size()>0 && (vtxBS.tracksSize()!=vertexTracks.size()) )
                {
                    AdaptiveVertexFitter theFitter;
                    TransientVertex v = theFitter.vertex(vertexTracks,PVRfP);
                    if ( v.isValid() )
                    {
                        vtxBSrf = reco::Vertex(v);
                    }
                }
                Double_t dx = (*qDecayVertex).position().x() - vtxBSrf.x();
                Double_t dy = (*qDecayVertex).position().y() - vtxBSrf.y();
                Double_t dz = (*qDecayVertex).position().z() - vtxBSrf.z();
                Double_t cosAlphaXYb = ( qCand->currentState().globalMomentum().x() * dx + qCand->currentState().globalMomentum().y()*dy + qCand->currentState().globalMomentum().z()*dz  )/( sqrt(dx*dx+dy*dy+dz*dz)* qCand->currentState().globalMomentum().mag() );


                if(cosAlphaXYb>lip)
                {

                    lip = cosAlphaXYb ;
                    pVtxBSIPX_temp     = vtxBSrf.x();
                    pVtxBSIPY_temp     = vtxBSrf.y();
                    pVtxBSIPZ_temp     = vtxBSrf.z();
                    pVtxBSIPXE_temp    = vtxBSrf.covariance(0, 0);
                    pVtxBSIPYE_temp    = vtxBSrf.covariance(1, 1);
                    pVtxBSIPZE_temp    = vtxBSrf.covariance(2, 2);
                    pVtxBSIPCL_temp    = (TMath::Prob(vtxBSrf.chi2(), (int)vtxBSrf.ndof()) );
                    pVtxBSIPdN_temp    = vtxBS.tracksSize() - vertexTracks.size();
                    bestVtxBSIP = vtxBSrf;
                }
            }

            // D0 vtx from PV
            double tmp_dx   = qDecayVertex->position().x() - pVtxBSIPX_temp ;
            double tmp_dy   = qDecayVertex->position().y() - pVtxBSIPY_temp ;
//             double tmp_dz   = qDecayVertex->position().z() - pVtxBSIPZ_temp ;
            double tmp_edx  = fabs(qDecayVertex->error().cxx()) + fabs(pVtxBSIPXE_temp) ;
            double tmp_edy  = fabs(qDecayVertex->error().cyy()) + fabs(pVtxBSIPYE_temp) ;
//             double tmp_edz  = qDecayVertex->error().czz() + pVtxBSIPZE_temp ;
            double tmp_px   = qCand->currentState().globalMomentum().x();
            double tmp_py   = qCand->currentState().globalMomentum().y();
//             double tmp_pz   = qCand->currentState().globalMomentum().z();
            //
            // double tmp_DS3  = sqrt(tmp_dx*tmp_dx/tmp_edx + tmp_dy*tmp_dy/tmp_edy + tmp_dz*tmp_dz/tmp_edz );
            double tmp_DS2  = sqrt((tmp_dx*tmp_dx)/tmp_edx + (tmp_dy*tmp_dy)/tmp_edy );
            if( tmp_DS2 < 3) continue;
            double tmp_Cos2 = (tmp_dx * tmp_px + tmp_dy * tmp_py) / (sqrt(tmp_dx * tmp_dx + tmp_dy * tmp_dy )*sqrt(tmp_px * tmp_px + tmp_py * tmp_py ));
            if( tmp_Cos2 < 0.5) continue;


            ///slow pion for first D* from first D0
            for (pat::PackedCandidateCollection::const_iterator iTraS1 = tracks->begin(); iTraS1 != tracks->end(); iTraS1++)
            {
                if (!(iTraS1->hasTrackDetails())) continue;
                auto iTrackS1 = iTraS1->bestTrack();
                if (iTrackS1 == nullptr) continue;
                //
                if (iTrackS1->charge() < 0.5) continue;
                if ( (iTrackS1->charge() * iTrack2->charge()) < 0.5) continue;
                //
                TLorentzVector p4piS1;
                p4piS1.SetPtEtaPhiM(iTrackS1->pt(),iTrackS1->eta(),iTrackS1->phi(), PDG_PION_MASS);
                if ((p4piS1 + p4ka1 + p4pi2).M() - (p4ka1 + p4pi2).M() + PDG_DZ_MASS > 2.080) continue; // D*(2010) selection, rough! better after vertex fit. However this is not too tight and not too loose
                //
                int piS1quality = 0; // 1 = loose, 2 = tight, 4 = highPurity
                if(iTrackS1->quality(reco::TrackBase::loose))
                {
                    piS1quality += 1;
                }
                if(iTrackS1->quality(reco::TrackBase::tight))
                {
                    piS1quality += 2;
                }
                if(iTrackS1->quality(reco::TrackBase::highPurity))
                {
                    piS1quality += 4;
                }
                //
                if(iTrackS1 == iTrack1) continue;
                if(iTrackS1 == iTrack2) continue;
                if (iTrackS1->charge() == 0) continue;
                //
                if( fabs(iTrackS1->eta()) > 2.6 ) continue;
                //
                TransientTrack pionSTT1 = theB->build( iTraS1->pseudoTrack() );
                if (!pionSTT1.isValid()) continue;

                ///////////////////////////////  FIRST D0 VTX TMP /////////////////////////////////////////

                std::vector<RefCountedKinematicParticle> W_candidate_init; //

                W_candidate_init.push_back(pFactory.particle(kaon1TT, PM_PDG_KAON_MASS, chi,ndf, PM_sigma));
                W_candidate_init.push_back(pFactory.particle(pion2TT, PM_PDG_PION_MASS, chi,ndf, PM_sigma));
                RefCountedKinematicTree wertexFitTree;

                std::vector<RefCountedKinematicParticle> W_candidate = W_candidate_init;
                KinematicParticleVertexFitter wFitter; //KinematicParticleVertexFitter
                fitgood = 1;
                try
                {
                    wertexFitTree = wFitter.fit(W_candidate);
                }
                catch (const VertexException &)
                {
                    fitgood = 0;
                }
                if (fitgood == 0) continue;
                if (!wertexFitTree->isValid()) continue;

                wertexFitTree->movePointerToTheTop();

                RefCountedKinematicParticle wCand         = wertexFitTree->currentParticle();
                RefCountedKinematicVertex   wDecayVertex  = wertexFitTree->currentDecayVertex();

                if(wDecayVertex->chiSquared() < 0) continue;
                if (!wDecayVertex->vertexIsValid())  continue;

                double piS1dz   = wDecayVertex->position().z() - iTrackS1->vz();
                double piS1dr   = sqrt( (wDecayVertex->position().x() - iTrackS1->vx())*(wDecayVertex->position().x() - iTrackS1->vx()) + (wDecayVertex->position().y() - iTrackS1->vy())*(wDecayVertex->position().y() - iTrackS1->vy()) );
                  ///////////////////////////////  FIRST D* VTX FIT /////////////////////////////////////////

                std::vector<RefCountedKinematicParticle> C_candidate_init; //

                C_candidate_init.push_back(wCand);
                C_candidate_init.push_back(pFactory.particle(pionSTT1, PM_PDG_PION_MASS, chi,ndf, PM_sigma));
                RefCountedKinematicTree certexFitTree;

                std::vector<RefCountedKinematicParticle> C_candidate = C_candidate_init;
                KinematicParticleVertexFitter cFitter; //KinematicParticleVertexFitter
                fitgood = 1;
                try
                {
                    certexFitTree = cFitter.fit(C_candidate);
                }
                catch (const VertexException &)
                {
                    fitgood = 0;
                }
                if (fitgood == 0) continue;
                //
                if (!certexFitTree->isValid()) continue;
                certexFitTree->movePointerToTheTop();

                RefCountedKinematicParticle cCand         = certexFitTree->currentParticle();
                RefCountedKinematicVertex   cDecayVertex  = certexFitTree->currentDecayVertex();
                if(cDecayVertex->chiSquared() < 0) continue;
                double C_Prob_tmp   = TMath::Prob(cDecayVertex->chiSquared(), (int) cDecayVertex->degreesOfFreedom());
                if(C_Prob_tmp < 0.01) continue;
                // get children from final first D* fit
                if (!cDecayVertex->vertexIsValid())  continue;

                double C_mass_tmp = cCand->currentState().mass();

                certexFitTree->movePointerToTheFirstChild();
                RefCountedKinematicParticle W0Cand    = certexFitTree->currentParticle();
                certexFitTree->movePointerToTheTop();
                //
		//
                if(C_mass_tmp - W0Cand->currentState().mass() + PDG_DZ_MASS > 2.025) continue;
                //
                //
    // Find second kaon for second D0
		for (pat::PackedCandidateCollection::const_iterator iTra3 = tracks->begin(); iTra3 != tracks->end(); iTra3++)
		{
			if (!(iTra3->hasTrackDetails())) continue;
			auto iTrack3 = iTra3->bestTrack();
			if (iTrack3 == nullptr) continue;
      if (iTrack3->charge() < 0.5) continue; // кат на заряд для второго отрицательного D*
			if(iTrack2 == iTrack3) continue;
			if(iTrack1 == iTrack3) continue;
			if(iTrackS1 == iTrack3) continue;
			//
			if (iTrack3->pt()<0.5) continue;
			if (!(iTrack3->quality(reco::TrackBase::highPurity))) continue;
			//
			TLorentzVector p4ka3;
			p4ka3.SetPtEtaPhiM(iTrack3->pt(),iTrack3->eta(),iTrack3->phi(), PDG_KAON_MASS);
			//
			if( fabs(iTrack3->eta()) > 2.6 ) continue;
			//
			TransientTrack kaon3TT = theB->build( iTra3->pseudoTrack() );
			if (!kaon3TT.isValid()) continue;
			// find second pion for second D0
			for (pat::PackedCandidateCollection::const_iterator iTra4 = tracks->begin(); iTra4 != tracks->end(); iTra4++)
			{
				if (!(iTra4->hasTrackDetails())) continue;
				auto iTrack4 = iTra4->bestTrack();
				if (iTrack4 == nullptr) continue;
				//
				if (iTrack4->pt()<0.50) continue;
				if (iTrack3->charge() * iTrack4->charge() > -0.5) continue; //OS tracks  negative
				if (!(iTrack4->quality(reco::TrackBase::highPurity))) continue;

				if((iTrack4->pt()<0.6)&&(iTrack3->pt()<0.6)) continue;
				//
				TLorentzVector p4pi4;
				p4pi4.SetPtEtaPhiM(iTrack4->pt(),iTrack4->eta(),iTrack4->phi(), PDG_PION_MASS);
				if ((p4ka3 + p4pi4).M() > 2.02) continue; // D0 selection, < 2.0 to get upto 1.91
				if ((p4ka3 + p4pi4).M() < 1.72) continue; // D0 selection, > 1.75 to get upto 1.82
				if ((p4ka3 + p4pi4).Pt() < 2.9) continue; // D0 pt > 3
				//
				if(iTrack3 == iTrack4) continue;
				if(iTrack2 == iTrack4) continue;
				if(iTrack1 == iTrack4) continue;
				if(iTrackS1 == iTrack4) continue;
				//
				if( fabs(iTrack4->eta()) > 2.6 ) continue;
				//
				TransientTrack pion4TT = theB->build( iTra4->pseudoTrack() );
				if (!pion4TT.isValid()) continue;
				//
				double pi4dz   = iTrack4->vz() - iTrack3->vz();
				double pi4dr   = sqrt( (iTrack4->vx() - iTrack3->vx())*(iTrack4->vx() - iTrack3->vx()) + (iTrack4->vy() - iTrack3->vy())*(iTrack4->vy() - iTrack3->vy()) );
				//
				if (pi4dz > 3.0) continue;
				if (pi4dz < -3.0) continue;


				  ///////////////////////////////  SECOND D0 FIT  /////////////////////////////////////////

				std::vector<RefCountedKinematicParticle> S_candidate_init; //

				S_candidate_init.push_back(pFactory.particle(kaon3TT, PM_PDG_KAON_MASS, chi,ndf, PM_sigma));
				S_candidate_init.push_back(pFactory.particle(pion4TT, PM_PDG_PION_MASS, chi,ndf, PM_sigma));
				RefCountedKinematicTree sertexFitTree;

				std::vector<RefCountedKinematicParticle> S_candidate = S_candidate_init;
				KinematicParticleVertexFitter sFitter; //KinematicParticleVertexFitter
				fitgood = 1;
				try
				{
				    sertexFitTree = sFitter.fit(S_candidate);
				}
				catch (const VertexException &)
				{
				    fitgood = 0;
				}
				if (fitgood == 0) continue;
				if (!sertexFitTree->isValid()) continue;
				sertexFitTree->movePointerToTheTop();

				RefCountedKinematicParticle sCand         = sertexFitTree->currentParticle();
				RefCountedKinematicVertex   sDecayVertex  = sertexFitTree->currentDecayVertex();
				double S_Prob_tmp   = 0;
				fitgood = 1;
				try
				{
				    if(sDecayVertex->chiSquared() < 0) fitgood = 0;
				    S_Prob_tmp   = TMath::Prob(sDecayVertex->chiSquared(), (int) sDecayVertex->degreesOfFreedom());
				}
				catch (const Exception &)
				{
				    fitgood = 0;
				}
				if (fitgood == 0) continue;
				//
				if(S_Prob_tmp < 0.05) continue;
				if (!sDecayVertex->vertexIsValid())  continue;

				double S_mass_tmp = sCand->currentState().mass();
				if(S_mass_tmp < 1.77) continue;
				if(S_mass_tmp > 1.97) continue;

        // D0 vtx from PV
                tmp_dx   = sDecayVertex->position().x() - pVtxBSIPX_temp ;
                tmp_dy   = sDecayVertex->position().y() - pVtxBSIPY_temp ;
//              double tmp_dz   = qDecayVertex->position().z() - pVtxBSIPZ_temp ;
                tmp_edx  = fabs(sDecayVertex->error().cxx()) + fabs(pVtxBSIPXE_temp) ;
                tmp_edy  = fabs(sDecayVertex->error().cyy()) + fabs(pVtxBSIPYE_temp) ;
//              double tmp_edz  = qDecayVertex->error().czz() + pVtxBSIPZE_temp ;
                tmp_px   = sCand->currentState().globalMomentum().x();
                tmp_py   = sCand->currentState().globalMomentum().y();
//             double tmp_pz   = qCand->currentState().globalMomentum().z();
        //
                // double tmp_DS3  = sqrt(tmp_dx*tmp_dx/tmp_edx + tmp_dy*tmp_dy/tmp_edy + tmp_dz*tmp_dz/tmp_edz );
                double tmp_ES2  = sqrt((tmp_dx*tmp_dx)/tmp_edx + (tmp_dy*tmp_dy)/tmp_edy );
                if( tmp_ES2 < 3) continue;
                double tmp_E_Cos2 = (tmp_dx * tmp_px + tmp_dy * tmp_py) / (sqrt(tmp_dx * tmp_dx + tmp_dy * tmp_dy )*sqrt(tmp_px * tmp_px + tmp_py * tmp_py ));
                if( tmp_E_Cos2 < 0.99) continue;


        //////////////////////////////ВОТ ТУТ ДОЛЖЕН СТОЯТЬ КОД ДЛЯ ВТОРОГО МЯГКОГО ПИОНА////////////////////////////
        ///slow pion for first D* from first D0
        for (pat::PackedCandidateCollection::const_iterator iTraS2 = tracks->begin(); iTraS2 != tracks->end(); iTraS2++)
        {
            if (!(iTraS2->hasTrackDetails())) continue;
            auto iTrackS2 = iTraS2->bestTrack();
            if (iTrackS2 == nullptr) continue;
            //
            if ( (iTrackS2->charge() * iTrack4->charge()) < 0.5) continue;
            //
            TLorentzVector p4piS2;
            p4piS2.SetPtEtaPhiM(iTrackS2->pt(),iTrackS2->eta(),iTrackS2->phi(), PDG_PION_MASS);
            if ((p4piS2 + p4ka3 + p4pi4).M() - (p4ka3 + p4pi4).M() + PDG_DZ_MASS > 2.080) continue; // D*(2010) selection, rough! better after vertex fit. However this is not too tight and not too loose
            //
            int piS2quality = 0; // 1 = loose, 2 = tight, 4 = highPurity
            if(iTrackS2->quality(reco::TrackBase::loose))
            {
                piS2quality += 1;
            }
            if(iTrackS2->quality(reco::TrackBase::tight))
            {
                piS2quality += 2;
            }
            if(iTrackS2->quality(reco::TrackBase::highPurity))
            {
                piS2quality += 4;
            }
            //
            if(iTrackS2 == iTrack3) continue;
            if(iTrackS2 == iTrack4) continue;
            if(iTrackS2 == iTrack1) continue;
            if(iTrackS2 == iTrack2) continue;
            if(iTrackS2 == iTrackS1) continue;
            if (iTrackS2->charge() == 0) continue;
            //
            if( fabs(iTrackS2->eta()) > 2.6 ) continue;
            //
            TransientTrack pionSTT2 = theB->build( iTraS2->pseudoTrack() );
            if (!pionSTT2.isValid()) continue;

            ///////////////////////////////  SECOND D0 VTX TMP /////////////////////////////////////////
            std::vector<RefCountedKinematicParticle> E_candidate_init; //
            E_candidate_init.push_back(pFactory.particle(kaon3TT, PM_PDG_KAON_MASS, chi,ndf, PM_sigma));
            E_candidate_init.push_back(pFactory.particle(pion4TT, PM_PDG_PION_MASS, chi,ndf, PM_sigma));
            RefCountedKinematicTree eertexFitTree;
            std::vector<RefCountedKinematicParticle> E_candidate = E_candidate_init;
            KinematicParticleVertexFitter eFitter; //KinematicParticleVertexFitter
            fitgood = 1;
            try
            {
                eertexFitTree = eFitter.fit(E_candidate);
            }
            catch (const VertexException &)
            {
                fitgood = 0;
            }
            if (fitgood == 0) continue;
            if (!eertexFitTree->isValid()) continue;
            eertexFitTree->movePointerToTheTop();
            RefCountedKinematicParticle eCand         = eertexFitTree->currentParticle();
            RefCountedKinematicVertex   eDecayVertex  = eertexFitTree->currentDecayVertex();
            if(eDecayVertex->chiSquared() < 0) continue;
            double E_Prob_tmp   = TMath::Prob(eDecayVertex->chiSquared(), (int) eDecayVertex->degreesOfFreedom());
            if (!eDecayVertex->vertexIsValid())  continue;
            double E_mass_tmp = eCand->currentState().mass();

            double piS2dz   = eDecayVertex->position().z() - iTrackS2->vz();
            double piS2dr   = sqrt( (eDecayVertex->position().x() - iTrackS2->vx())*(eDecayVertex->position().x() - iTrackS2->vx()) + (eDecayVertex->position().y() - iTrackS2->vy())*(eDecayVertex->position().y() - iTrackS2->vy()) );

            eertexFitTree->movePointerToTheFirstChild();
            RefCountedKinematicParticle KA3Cand    = eertexFitTree->currentParticle();
            eertexFitTree->movePointerToTheNextChild();
            RefCountedKinematicParticle PI4Cand    = eertexFitTree->currentParticle();
            eertexFitTree->movePointerToTheTop();
            GlobalVector KA3Cand_p = KA3Cand->currentState().kinematicParameters().momentum();
            GlobalVector PI4Cand_p = PI4Cand->currentState().kinematicParameters().momentum();

              ///////////////////////////////  SECOND D* VTX FIT /////////////////////////////////////////

            std::vector<RefCountedKinematicParticle> F_candidate_init; //

            F_candidate_init.push_back(eCand);
            F_candidate_init.push_back(pFactory.particle(pionSTT2, PM_PDG_PION_MASS, chi,ndf, PM_sigma));
            RefCountedKinematicTree fertexFitTree;

            std::vector<RefCountedKinematicParticle> F_candidate = F_candidate_init;
            KinematicParticleVertexFitter fFitter; //KinematicParticleVertexFitter
            fitgood = 1;
            try
            {
                fertexFitTree = fFitter.fit(F_candidate);
            }
            catch (const VertexException &)
            {
                fitgood = 0;
            }
            if (fitgood == 0) continue;
            //
            if (!fertexFitTree->isValid()) continue;
            fertexFitTree->movePointerToTheTop();

            RefCountedKinematicParticle fCand         = fertexFitTree->currentParticle();
            RefCountedKinematicVertex   fDecayVertex  = fertexFitTree->currentDecayVertex();
            if(fDecayVertex->chiSquared() < 0) continue;
            double F_Prob_tmp   = TMath::Prob(fDecayVertex->chiSquared(), (int) fDecayVertex->degreesOfFreedom());
            if(F_Prob_tmp < 0.01) continue;
            if (!fDecayVertex->vertexIsValid())  continue;

            double F_mass_tmp = fCand->currentState().mass();

            fertexFitTree->movePointerToTheFirstChild();
            RefCountedKinematicParticle F0Cand    = fertexFitTree->currentParticle();
            fertexFitTree->movePointerToTheTop();
            //
//
            if(F_mass_tmp - F0Cand->currentState().mass() + PDG_DZ_MASS > 2.025) continue;
            //
            //

        /////////////////////////////ВОТ ТУТ ДОЛЖЕН КОНЧАТЬСЯ КОД ДЛЯ ВТОРОГО D* ////////////////////////////////////

				//// refitting 1st d0 before final fit
				//////////////////////////////////////////////
				std::vector<RefCountedKinematicParticle> D_candidate_init; //
				D_candidate_init.push_back(pFactory.particle(kaon1TT, PM_PDG_KAON_MASS, chi,ndf, PM_sigma));
				D_candidate_init.push_back(pFactory.particle(pion2TT, PM_PDG_PION_MASS, chi,ndf, PM_sigma));
				RefCountedKinematicTree dertexFitTree;
				std::vector<RefCountedKinematicParticle> D_candidate = D_candidate_init;
				KinematicParticleVertexFitter dFitter; //KinematicParticleVertexFitter
				fitgood = 1;
				try
				{
				    dertexFitTree = dFitter.fit(D_candidate);
				}
				catch (const VertexException &)
				{
				    fitgood = 0;
				}
				if (fitgood == 0) continue;
				if (!dertexFitTree->isValid()) continue;
				dertexFitTree->movePointerToTheTop();
				RefCountedKinematicParticle dCand         = dertexFitTree->currentParticle();
				RefCountedKinematicVertex   dDecayVertex  = dertexFitTree->currentDecayVertex();
				if(dDecayVertex->chiSquared() < 0) continue;
				double D_Prob_tmp   = TMath::Prob(dDecayVertex->chiSquared(), (int) dDecayVertex->degreesOfFreedom());
				if (!dDecayVertex->vertexIsValid())  continue;
				double D_mass_tmp = dCand->currentState().mass();
				dertexFitTree->movePointerToTheFirstChild();
				RefCountedKinematicParticle KA1Cand    = dertexFitTree->currentParticle();
				dertexFitTree->movePointerToTheNextChild();
				RefCountedKinematicParticle PI2Cand    = dertexFitTree->currentParticle();
				dertexFitTree->movePointerToTheTop();
				GlobalVector KA1Cand_p = KA1Cand->currentState().kinematicParameters().momentum();
				GlobalVector PI2Cand_p = PI2Cand->currentState().kinematicParameters().momentum();
				//


				//// FINAL D* D0 fit
				std::vector<RefCountedKinematicParticle> B_candidate_init; //

				B_candidate_init.push_back(dCand);
				B_candidate_init.push_back(pFactory.particle(pionSTT1, PM_PDG_PION_MASS, chi,ndf, PM_sigma));
				B_candidate_init.push_back(eCand);
        B_candidate_init.push_back(pFactory.particle(pionSTT2, PM_PDG_PION_MASS, chi,ndf, PM_sigma));
				RefCountedKinematicTree vertexFitTree;

				std::vector<RefCountedKinematicParticle> B_candidate = B_candidate_init;
				KinematicParticleVertexFitter pFitter; //KinematicParticleVertexFitter
				fitgood = 1;
				try
				{
				    vertexFitTree = pFitter.fit(B_candidate);
				}
				catch (const VertexException &)
				{
				    fitgood = 0;
				}
				if (fitgood == 0) continue;
				//
				if (!vertexFitTree->isValid()) continue;
				vertexFitTree->movePointerToTheTop();

				RefCountedKinematicParticle bCand         = vertexFitTree->currentParticle();
				RefCountedKinematicVertex   bDecayVertex  = vertexFitTree->currentDecayVertex();
				if(bDecayVertex->chiSquared() < 0) continue;
				double B_Prob_tmp   = TMath::Prob(bDecayVertex->chiSquared(), (int) bDecayVertex->degreesOfFreedom());
				if(B_Prob_tmp < 0.05) continue;
				// get children from final B fit
				if (!bDecayVertex->vertexIsValid())  continue;

				double B_mass_tmp = bCand->currentState().mass();

				vertexFitTree->movePointerToTheFirstChild();
				RefCountedKinematicParticle D0Cand    = vertexFitTree->currentParticle();
				vertexFitTree->movePointerToTheNextChild();
				RefCountedKinematicParticle PIS1Cand   = vertexFitTree->currentParticle();
				vertexFitTree->movePointerToTheNextChild();
				RefCountedKinematicParticle E0Cand    = vertexFitTree->currentParticle();
				vertexFitTree->movePointerToTheNextChild();
        RefCountedKinematicParticle PIS2Cand   = vertexFitTree->currentParticle();
				vertexFitTree->movePointerToTheTop();
				//
				GlobalVector D0Cand_p  = D0Cand->currentState().kinematicParameters().momentum();
				GlobalVector PIS1Cand_p = PIS1Cand->currentState().kinematicParameters().momentum();
				GlobalVector E0Cand_p  = E0Cand->currentState().kinematicParameters().momentum();
        GlobalVector PIS2Cand_p = PIS2Cand->currentState().kinematicParameters().momentum();

				//////////////////   SAVE
				//
				B_mass          ->push_back(B_mass_tmp);
				B_masd          ->push_back(B_mass_tmp - D0Cand->currentState().mass() - E0Cand->currentState().mass() + 2*PDG_DZ_MASS);
				B_px            ->push_back(bCand->currentState().globalMomentum().x());
				B_py            ->push_back(bCand->currentState().globalMomentum().y());
				B_pz            ->push_back(bCand->currentState().globalMomentum().z());
				B_DecayVtxX     ->push_back(bDecayVertex->position().x());
				B_DecayVtxY     ->push_back(bDecayVertex->position().y());
				B_DecayVtxZ     ->push_back(bDecayVertex->position().z());
				B_DecayVtxXE    ->push_back(bDecayVertex->error().cxx());
				B_DecayVtxYE    ->push_back(bDecayVertex->error().cyy());
				B_DecayVtxZE    ->push_back(bDecayVertex->error().czz());
				B_Prob          ->push_back(B_Prob_tmp);

				D_mass          ->push_back(D_mass_tmp);
				D_mass_CV       ->push_back(D0Cand->currentState().mass() );
				D_px            ->push_back(dCand->currentState().globalMomentum().x());
				D_py            ->push_back(dCand->currentState().globalMomentum().y());
				D_pz            ->push_back(dCand->currentState().globalMomentum().z());
				D_px_CV         ->push_back(D0Cand_p.x());
				D_py_CV         ->push_back(D0Cand_p.y());
				D_pz_CV         ->push_back(D0Cand_p.z());
				D_DecayVtxX     ->push_back(dDecayVertex->position().x());
				D_DecayVtxY     ->push_back(dDecayVertex->position().y());
				D_DecayVtxZ     ->push_back(dDecayVertex->position().z());
				D_DecayVtxXE    ->push_back(dDecayVertex->error().cxx());
				D_DecayVtxYE    ->push_back(dDecayVertex->error().cyy());
				D_DecayVtxZE    ->push_back(dDecayVertex->error().czz());
				D_Prob          ->push_back(D_Prob_tmp);
				//
				//
				KA1_px          ->push_back( p4ka1.Px()                                             );
				KA1_py          ->push_back( p4ka1.Py()                                             );
				KA1_pz          ->push_back( p4ka1.Pz()                                             );
				KA1_px_CV       ->push_back( KA1Cand_p.x());
				KA1_py_CV       ->push_back( KA1Cand_p.y());
				KA1_pz_CV       ->push_back( KA1Cand_p.z());
				KA1_ips         ->push_back( fabs(iTrack1->dxy(bestVtxBSIP.position())) / (0.000001 + fabs(iTrack1->dxyError())) );
        KA1_charge      ->push_back(iTrack1->charge() );
				//
				PI2_px          ->push_back( p4pi2.Px()                                             );
				PI2_py          ->push_back( p4pi2.Py()                                             );
				PI2_pz          ->push_back( p4pi2.Pz()                                             );
				PI2_dr          ->push_back( pi2dr                                                  );
				PI2_dz          ->push_back( pi2dz                                                  );
				PI2_px_CV       ->push_back( PI2Cand_p.x());
				PI2_py_CV       ->push_back( PI2Cand_p.y());
				PI2_pz_CV       ->push_back( PI2Cand_p.z());
				PI2_ips         ->push_back( fabs(iTrack2->dxy(bestVtxBSIP.position())) / (0.000001 + fabs(iTrack2->dxyError())) );
				//
				//
				PIS1_charg       ->push_back( iTrackS1->charge()                                      );
				PIS1_purit       ->push_back( piS1quality                                             );
				PIS1_px          ->push_back( p4piS1.Px()                                             );
				PIS1_py          ->push_back( p4piS1.Py()                                             );
				PIS1_pz          ->push_back( p4piS1.Pz()                                             );
				PIS1_dr          ->push_back( piS1dr                                                  );
				PIS1_dz          ->push_back( piS1dz                                                  );
				PIS1_px_CV       ->push_back( PIS1Cand_p.x());
				PIS1_py_CV       ->push_back( PIS1Cand_p.y());
				PIS1_pz_CV       ->push_back( PIS1Cand_p.z());
				PIS1_ips         ->push_back( fabs(iTrackS1->dxy(bestVtxBSIP.position())) / (0.000001 + fabs(iTrackS1->dxyError())) );

        PIS2_charg       ->push_back( iTrackS2->charge()                                      );
        PIS2_purit       ->push_back( piS2quality                                             );
        PIS2_px          ->push_back( p4piS2.Px()                                             );
        PIS2_py          ->push_back( p4piS2.Py()                                             );
        PIS2_pz          ->push_back( p4piS2.Pz()                                             );
        PIS2_dr          ->push_back( piS2dr                                                  );
        PIS2_dz          ->push_back( piS2dz                                                  );
        PIS2_px_CV       ->push_back( PIS2Cand_p.x());
        PIS2_py_CV       ->push_back( PIS2Cand_p.y());
        PIS2_pz_CV       ->push_back( PIS2Cand_p.z());
        PIS2_ips         ->push_back( fabs(iTrackS2->dxy(bestVtxBSIP.position())) / (0.000001 + fabs(iTrackS2->dxyError())) );
				//
				E_mass          ->push_back(E_mass_tmp);
				E_mass_CV       ->push_back(E0Cand->currentState().mass() );
				E_px            ->push_back(eCand->currentState().globalMomentum().x());
				E_py            ->push_back(eCand->currentState().globalMomentum().y());
				E_pz            ->push_back(eCand->currentState().globalMomentum().z());
				E_px_CV         ->push_back(E0Cand_p.x());
				E_py_CV         ->push_back(E0Cand_p.y());
				E_pz_CV         ->push_back(E0Cand_p.z());
				E_DecayVtxX     ->push_back(eDecayVertex->position().x());
				E_DecayVtxY     ->push_back(eDecayVertex->position().y());
				E_DecayVtxZ     ->push_back(eDecayVertex->position().z());
				E_DecayVtxXE    ->push_back(eDecayVertex->error().cxx());
				E_DecayVtxYE    ->push_back(eDecayVertex->error().cyy());
				E_DecayVtxZE    ->push_back(eDecayVertex->error().czz());
				E_Prob          ->push_back(E_Prob_tmp);
				//
				//
				KA3_px          ->push_back( p4ka3.Px()                                             );
				KA3_py          ->push_back( p4ka3.Py()                                             );
				KA3_pz          ->push_back( p4ka3.Pz()                                             );
				KA3_px_CV       ->push_back( KA3Cand_p.x());
				KA3_py_CV       ->push_back( KA3Cand_p.y());
				KA3_pz_CV       ->push_back( KA3Cand_p.z());
				KA3_ips         ->push_back( fabs(iTrack3->dxy(bestVtxBSIP.position())) / (0.000001 + fabs(iTrack3->dxyError())) );
                KA3_charge      ->push_back(iTrack3->charge() );
				//
				PI4_px          ->push_back( p4pi4.Px()                                             );
				PI4_py          ->push_back( p4pi4.Py()                                             );
				PI4_pz          ->push_back( p4pi4.Pz()                                             );
				PI4_dr          ->push_back( pi4dr                                                  );
				PI4_dz          ->push_back( pi4dz                                                  );
				PI4_px_CV       ->push_back( PI4Cand_p.x());
				PI4_py_CV       ->push_back( PI4Cand_p.y());
				PI4_pz_CV       ->push_back( PI4Cand_p.z());
				PI4_ips         ->push_back( fabs(iTrack4->dxy(bestVtxBSIP.position())) / (0.000001 + fabs(iTrack4->dxyError())) );
				//
				PV_becos_XX     ->push_back( pVtxBSIPX_temp                                         );
				PV_becos_YY     ->push_back( pVtxBSIPY_temp                                         );
				PV_becos_ZZ     ->push_back( pVtxBSIPZ_temp                                         );
				PV_becos_EX     ->push_back( pVtxBSIPXE_temp                                        );
				PV_becos_EY     ->push_back( pVtxBSIPYE_temp                                        );
				PV_becos_EZ     ->push_back( pVtxBSIPZE_temp                                        );
				PV_becos_CL     ->push_back( pVtxBSIPCL_temp                                        );
				PV_becos_dN     ->push_back( pVtxBSIPdN_temp                                        );

				trig0_fi        ->push_back( TriggersFired[0]                                       );
				trig1_fi        ->push_back( TriggersFired[1]                                       );
				trig2_fi        ->push_back( TriggersFired[2]                                       );
				trig3_fi        ->push_back( TriggersFired[3]                                       );
				trig4_fi        ->push_back( TriggersFired[4]                                       );
				trig5_fi        ->push_back( TriggersFired[5]                                       );
				trig6_fi        ->push_back( TriggersFired[6]                                       );
				trig7_fi        ->push_back( TriggersFired[7]                                       );
				trig8_fi        ->push_back( TriggersFired[8]                                       );
				trig9_fi        ->push_back( TriggersFired[9]                                       );

				nCand ++;
				E_candidate_init.clear();
				E_candidate.clear();
				D_candidate_init.clear();
				D_candidate.clear();
				B_candidate_init.clear();
				B_candidate.clear();

    F_candidate_init.clear();
    F_candidate.clear();
      }// piS2
			} // pi4
		} //ka3
		C_candidate_init.clear();
		C_candidate.clear();
		W_candidate_init.clear();
		W_candidate.clear();
  } // piS1
            Q_candidate_init.clear();
            Q_candidate.clear();
        } // pi2
} // ka1

// for (reco::MuonCollection::const_iterator muon=muonColl->begin(), muonCollEnd=muonColl->end(); muon!=muonCollEnd;  ++muon)
// {
//     if(!(
//     muon->isGlobalMuon() &&
//     muon->globalTrack()->normalizedChi2() < 10 &&
//     muon->globalTrack()->hitPattern().numberOfValidMuonHits() > 0 &&
//     muon->numberOfMatchedStations() > 1 &&
//     muon->innerTrack()->hitPattern().numberOfValidPixelHits() > 0 &&
//     muon::isGoodMuon(*muon, muon::TMLastStationTight) &&
//     muon->innerTrack()->hitPattern().trackerLayersWithMeasurement() > 4 &&
//     muon->pt() > 3.0 && fabs(muon->eta()) < 2.4
//     )) continue;
//
//     mu_pt->push_back( muon->pt() ) ;
//     nCand ++;
// }

// ===================== END OF EVENT : WRITE ETC ++++++++++++++++++++++

if (nCand > 0)
{               wwtree->Fill();
} // nCand > 0

//
B_mass->clear();            B_masd->clear();
B_px->clear();              B_py->clear();              B_pz->clear();
B_DecayVtxX->clear();       B_DecayVtxY->clear();       B_DecayVtxZ->clear();
B_DecayVtxXE->clear();      B_DecayVtxYE->clear();      B_DecayVtxZE->clear();
B_Prob->clear();

D_mass->clear();            D_mass_CV->clear();
D_px->clear();              D_py->clear();              D_pz->clear();
D_px_CV->clear();           D_py_CV->clear();           D_pz_CV->clear();
D_DecayVtxX->clear();       D_DecayVtxY->clear();       D_DecayVtxZ->clear();
D_DecayVtxXE->clear();      D_DecayVtxYE->clear();      D_DecayVtxZE->clear();
D_Prob->clear();


KA1_px->clear();        KA1_py->clear();      KA1_pz->clear();
KA1_px_CV->clear();     KA1_py_CV->clear();   KA1_pz_CV->clear();
KA1_ips->clear();
KA1_charge->clear();

PI2_px->clear();        PI2_py->clear();      PI2_pz->clear();
PI2_dr->clear();        PI2_dz->clear();
PI2_px_CV->clear();     PI2_py_CV->clear();   PI2_pz_CV->clear();
PI2_ips->clear();


PIS1_charg->clear();     PIS1_purit->clear();
PIS1_px->clear();        PIS1_py->clear();      PIS1_pz->clear();
PIS1_dr->clear();        PIS1_dz->clear();
PIS1_px_CV->clear();     PIS1_py_CV->clear();   PIS1_pz_CV->clear();
PIS1_ips->clear();

PIS2_charg->clear();     PIS2_purit->clear();
PIS2_px->clear();        PIS2_py->clear();      PIS2_pz->clear();
PIS2_dr->clear();        PIS2_dz->clear();
PIS2_px_CV->clear();     PIS2_py_CV->clear();   PIS2_pz_CV->clear();
PIS2_ips->clear();

E_mass->clear();            E_mass_CV->clear();
E_px->clear();              E_py->clear();              E_pz->clear();
E_px_CV->clear();           E_py_CV->clear();           E_pz_CV->clear();
E_DecayVtxX->clear();       E_DecayVtxY->clear();       E_DecayVtxZ->clear();
E_DecayVtxXE->clear();      E_DecayVtxYE->clear();      E_DecayVtxZE->clear();
E_Prob->clear();


KA3_px->clear();        KA3_py->clear();      KA3_pz->clear();
KA3_px_CV->clear();     KA3_py_CV->clear();   KA3_pz_CV->clear();
KA3_ips->clear();
KA3_charge->clear();

PI4_px->clear();        PI4_py->clear();      PI4_pz->clear();
PI4_dr->clear();        PI4_dz->clear();
PI4_px_CV->clear();     PI4_py_CV->clear();   PI4_pz_CV->clear();
PI4_ips->clear();


PV_becos_XX->clear();   PV_becos_YY->clear();   PV_becos_ZZ->clear();
PV_becos_EX->clear();   PV_becos_EY->clear();   PV_becos_EZ->clear();
PV_becos_CL->clear();   PV_becos_dN->clear();

GEN_dectype->clear();
GEN_momID->clear();
GEN_PIS_px->clear();    GEN_PIS_py->clear();    GEN_PIS_pz->clear();
GEN_KS_px->clear();     GEN_KS_py->clear();     GEN_KS_pz->clear();
GEN_KR_px->clear();     GEN_KR_py->clear();     GEN_KR_pz->clear();
GEN_PI1_px->clear();    GEN_PI1_py->clear();    GEN_PI1_pz->clear();
GEN_PI2_px->clear();    GEN_PI2_py->clear();    GEN_PI2_pz->clear();

trig0_fi->clear();      trig1_fi->clear();      trig2_fi->clear();
trig3_fi->clear();      trig4_fi->clear();      trig5_fi->clear();
trig6_fi->clear();      trig7_fi->clear();      trig8_fi->clear();
trig9_fi->clear();

}


// ------------ method called once each job just before starting event loop  ------------
void Xb_frame::beginJob()
{
using namespace std;
using namespace reco;
//
cout << "------------------------------->>>>> Begin Job" << endl;

f = new TFile(fileName.c_str(), "RECREATE");
wwtree  = new TTree("wztree", "muons tree");

wwtree->Branch("nCand"              , &nCand            , "nCand/I"     );

wwtree->Branch("run"                , &run              , "run/I"       );
wwtree->Branch("event"              , &event            , "event/I"     );
wwtree->Branch("lumi"               , &lumi             , "lumi/F"      );

wwtree->Branch("BESP_x"             , &BESP_x           , "BESP_x/F"    );
wwtree->Branch("BESP_y"             , &BESP_y           , "BESP_y/F"    );
wwtree->Branch("BESP_z"             , &BESP_z           , "BESP_z/F"    );
wwtree->Branch("BESP_ex"            , &BESP_ex          , "BESP_ex/F"   );
wwtree->Branch("BESP_ey"            , &BESP_ey          , "BESP_ey/F"   );
wwtree->Branch("BESP_ez"            , &BESP_ez          , "BESP_ez/F"   );

wwtree->Branch("numPV"              , &numPV            , "numPV/I"     );
wwtree->Branch("numTrack"           , &numTrack         , "numTrack/I"  );
wwtree->Branch("numV0"              , &numV0            , "numV0/I"     );

wwtree->Branch("B_mass"             , &B_mass           );
wwtree->Branch("B_masd"             , &B_masd           );
wwtree->Branch("B_px"               , &B_px             );
wwtree->Branch("B_py"               , &B_py             );
wwtree->Branch("B_pz"               , &B_pz             );
wwtree->Branch("B_DecayVtxX"        , &B_DecayVtxX      );
wwtree->Branch("B_DecayVtxY"        , &B_DecayVtxY      );
wwtree->Branch("B_DecayVtxZ"        , &B_DecayVtxZ      );
wwtree->Branch("B_DecayVtxXE"       , &B_DecayVtxXE     );
wwtree->Branch("B_DecayVtxYE"       , &B_DecayVtxYE     );
wwtree->Branch("B_DecayVtxZE"       , &B_DecayVtxZE     );
wwtree->Branch("B_Prob"             , &B_Prob           );

wwtree->Branch("D_mass"             , &D_mass           );
wwtree->Branch("D_mass_CV"          , &D_mass_CV        );
wwtree->Branch("D_px"               , &D_px             );
wwtree->Branch("D_py"               , &D_py             );
wwtree->Branch("D_pz"               , &D_pz             );
wwtree->Branch("D_px_CV"            , &D_px_CV          );
wwtree->Branch("D_py_CV"            , &D_py_CV          );
wwtree->Branch("D_pz_CV"            , &D_pz_CV          );
wwtree->Branch("D_DecayVtxX"        , &D_DecayVtxX      );
wwtree->Branch("D_DecayVtxY"        , &D_DecayVtxY      );
wwtree->Branch("D_DecayVtxZ"        , &D_DecayVtxZ      );
wwtree->Branch("D_DecayVtxXE"       , &D_DecayVtxXE     );
wwtree->Branch("D_DecayVtxYE"       , &D_DecayVtxYE     );
wwtree->Branch("D_DecayVtxZE"       , &D_DecayVtxZE     );
wwtree->Branch("D_Prob"             , &D_Prob           );
//
//

wwtree->Branch("KA1_px"             , &KA1_px           );
wwtree->Branch("KA1_py"             , &KA1_py           );
wwtree->Branch("KA1_pz"             , &KA1_pz           );
wwtree->Branch("KA1_px_CV"          , &KA1_px_CV        );
wwtree->Branch("KA1_py_CV"          , &KA1_py_CV        );
wwtree->Branch("KA1_pz_CV"          , &KA1_pz_CV        );
wwtree->Branch("KA1_ips"            , &KA1_ips          );
wwtree->Branch("KA1_charge"            , &KA1_charge          );
//
wwtree->Branch("PI2_px"             , &PI2_px           );
wwtree->Branch("PI2_py"             , &PI2_py           );
wwtree->Branch("PI2_pz"             , &PI2_pz           );
wwtree->Branch("PI2_dr"             , &PI2_dr           );
wwtree->Branch("PI2_dz"             , &PI2_dz           );
wwtree->Branch("PI2_px_CV"          , &PI2_px_CV        );
wwtree->Branch("PI2_py_CV"          , &PI2_py_CV        );
wwtree->Branch("PI2_pz_CV"          , &PI2_pz_CV        );
wwtree->Branch("PI2_ips"            , &PI2_ips          );

//
wwtree->Branch("PIS1_charg"          , &PIS1_charg        );
wwtree->Branch("PIS1_purit"          , &PIS1_purit        );
wwtree->Branch("PIS1_px"             , &PIS1_px           );
wwtree->Branch("PIS1_py"             , &PIS1_py           );
wwtree->Branch("PIS1_pz"             , &PIS1_pz           );
wwtree->Branch("PIS1_dr"             , &PIS1_dr           );
wwtree->Branch("PIS1_dz"             , &PIS1_dz           );
wwtree->Branch("PIS1_px_CV"          , &PIS1_px_CV        );
wwtree->Branch("PIS1_py_CV"          , &PIS1_py_CV        );
wwtree->Branch("PIS1_pz_CV"          , &PIS1_pz_CV        );
wwtree->Branch("PIS1_ips"            , &PIS1_ips          );
//
wwtree->Branch("PIS2_charg"          , &PIS2_charg        );
wwtree->Branch("PIS2_purit"          , &PIS2_purit        );
wwtree->Branch("PIS2_px"             , &PIS2_px           );
wwtree->Branch("PIS2_py"             , &PIS2_py           );
wwtree->Branch("PIS2_pz"             , &PIS2_pz           );
wwtree->Branch("PIS2_dr"             , &PIS2_dr           );
wwtree->Branch("PIS2_dz"             , &PIS2_dz           );
wwtree->Branch("PIS2_px_CV"          , &PIS2_px_CV        );
wwtree->Branch("PIS2_py_CV"          , &PIS2_py_CV        );
wwtree->Branch("PIS2_pz_CV"          , &PIS2_pz_CV        );
wwtree->Branch("PIS2_ips"            , &PIS2_ips          );
//
wwtree->Branch("E_mass"             , &E_mass           );
wwtree->Branch("E_mass_CV"          , &E_mass_CV        );
wwtree->Branch("E_px"               , &E_px             );
wwtree->Branch("E_py"               , &E_py             );
wwtree->Branch("E_pz"               , &E_pz             );
wwtree->Branch("E_px_CV"            , &E_px_CV          );
wwtree->Branch("E_py_CV"            , &E_py_CV          );
wwtree->Branch("E_pz_CV"            , &E_pz_CV          );
wwtree->Branch("E_DecayVtxX"        , &E_DecayVtxX      );
wwtree->Branch("E_DecayVtxY"        , &E_DecayVtxY      );
wwtree->Branch("E_DecayVtxZ"        , &E_DecayVtxZ      );
wwtree->Branch("E_DecayVtxXE"       , &E_DecayVtxXE     );
wwtree->Branch("E_DecayVtxYE"       , &E_DecayVtxYE     );
wwtree->Branch("E_DecayVtxZE"       , &E_DecayVtxZE     );
wwtree->Branch("E_Prob"             , &E_Prob           );
//
//

wwtree->Branch("KA3_px"             , &KA3_px           );
wwtree->Branch("KA3_py"             , &KA3_py           );
wwtree->Branch("KA3_pz"             , &KA3_pz           );
wwtree->Branch("KA3_px_CV"          , &KA3_px_CV        );
wwtree->Branch("KA3_py_CV"          , &KA3_py_CV        );
wwtree->Branch("KA3_pz_CV"          , &KA3_pz_CV        );
wwtree->Branch("KA3_ips"            , &KA3_ips          );
wwtree->Branch("KA3_charge"            , &KA3_charge          );

//

wwtree->Branch("PI4_px"             , &PI4_px           );
wwtree->Branch("PI4_py"             , &PI4_py           );
wwtree->Branch("PI4_pz"             , &PI4_pz           );
wwtree->Branch("PI4_dr"             , &PI4_dr           );
wwtree->Branch("PI4_dz"             , &PI4_dz           );
wwtree->Branch("PI4_px_CV"          , &PI4_px_CV        );
wwtree->Branch("PI4_py_CV"          , &PI4_py_CV        );
wwtree->Branch("PI4_pz_CV"          , &PI4_pz_CV        );
wwtree->Branch("PI4_ips"            , &PI4_ips          );

wwtree->Branch("PV_becos_XX"        , &PV_becos_XX      );
wwtree->Branch("PV_becos_YY"        , &PV_becos_YY      );
wwtree->Branch("PV_becos_ZZ"        , &PV_becos_ZZ      );
wwtree->Branch("PV_becos_EX"        , &PV_becos_EX      );
wwtree->Branch("PV_becos_EY"        , &PV_becos_EY      );
wwtree->Branch("PV_becos_EZ"        , &PV_becos_EZ      );
wwtree->Branch("PV_becos_CL"        , &PV_becos_CL      );
wwtree->Branch("PV_becos_dN"        , &PV_becos_dN      );


wwtree->Branch("GEN_dectype"        , &GEN_dectype      );
wwtree->Branch("GEN_momID"          , &GEN_momID        );
wwtree->Branch("GEN_PIS_px"         , &GEN_PIS_px       );
wwtree->Branch("GEN_PIS_py"         , &GEN_PIS_py       );
wwtree->Branch("GEN_PIS_pz"         , &GEN_PIS_pz       );
wwtree->Branch("GEN_KS_px"          , &GEN_KS_px        );
wwtree->Branch("GEN_KS_py"          , &GEN_KS_py        );
wwtree->Branch("GEN_KS_pz"          , &GEN_KS_pz        );
wwtree->Branch("GEN_KR_px"          , &GEN_KR_px        );
wwtree->Branch("GEN_KR_py"          , &GEN_KR_py        );
wwtree->Branch("GEN_KR_pz"          , &GEN_KR_pz        );
wwtree->Branch("GEN_PI1_px"         , &GEN_PI1_px       );
wwtree->Branch("GEN_PI1_py"         , &GEN_PI1_py       );
wwtree->Branch("GEN_PI1_pz"         , &GEN_PI1_pz       );
wwtree->Branch("GEN_PI2_px"         , &GEN_PI2_px       );
wwtree->Branch("GEN_PI2_py"         , &GEN_PI2_py       );
wwtree->Branch("GEN_PI2_pz"         , &GEN_PI2_pz       );

wwtree->Branch("trig0_fi"           , &trig0_fi         );
wwtree->Branch("trig1_fi"           , &trig1_fi         );
wwtree->Branch("trig2_fi"           , &trig2_fi         );
wwtree->Branch("trig3_fi"           , &trig3_fi         );
wwtree->Branch("trig4_fi"           , &trig4_fi         );
wwtree->Branch("trig5_fi"           , &trig5_fi         );
wwtree->Branch("trig6_fi"           , &trig6_fi         );
wwtree->Branch("trig7_fi"           , &trig7_fi         );
wwtree->Branch("trig8_fi"           , &trig8_fi         );
wwtree->Branch("trig9_fi"           , &trig9_fi         );
}

// ------------ method called once each job just after ending the event loop  ------------
void
Xb_frame::endJob()
{

using namespace std;
cout << "------------------------------->>>>> End Job" << endl;
f->WriteTObject(wwtree);
delete wwtree;
f->Close();

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Xb_frame::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Xb_frame);
///
//  Path? : HLT_Mu12_IP6_part0_v2
//  Path? : HLT_Mu9_IP5_part0_v2
//  Path? : HLT_Mu9_IP6_part0_v3
//  Path? : HLT_Mu9_IP4_part0_v2
//  Path? : HLT_Mu8_IP5_part0_v2
//  Path? : HLT_Mu8_IP6_part0_v2
//  Path? : HLT_Mu8_IP3_part0_v3
//  Path? : HLT_Mu7_IP4_part0_v2
