//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Dec 18 09:44:35 2014 by ROOT version 5.34/23
// from TTree GenLevelTree/GenLevelTree
// found on file: GenLevelAnalyzer_v01/DiPhotonJetsBox_Pt_32_20_dR_04_7TeV_sherpa_AODSIM.root
//////////////////////////////////////////////////////////

#ifndef GenLevelAcceptance_h
#define GenLevelAcceptance_h

//#define _DEBUG

#include <TROOT.h>
#include <TChain.h>
#include <TTree.h>
#include <TMath.h>
#include <TFile.h>
#include <TString.h>
#include <assert.h>
#include <utility>
#include <iostream>

// Header file for the classes stored in the TTree if any.
#include <vector>

using namespace std;

typedef unsigned int uint;
static const float pi = 3.14159265358979312;

// Fixed size dimensions of array or collections stored in the TTree if any.

class GenLevelAcceptance {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   ULong64_t       run;
   ULong64_t       ls;
   ULong64_t       evt;
   vector<float>   *photon_pt;
   vector<float>   *photon_eta;
   vector<float>   *photon_phi;
   vector<int>     *photon_mother_id;
   vector<int>     *photon_mother_status;
   vector<float>   *photon_geniso04;
   vector<float>   *jet_pt;
   vector<float>   *jet_eta;
   vector<float>   *jet_phi;
   vector<float>   *jet_energy;
   vector<float>   *parton_pt;
   vector<float>   *parton_eta;
   vector<float>   *parton_phi;
   vector<int>     *parton_id;
   vector<int>     *parton_status;
   vector<float>   *weights_geninfo;
   vector<float>   *weights_lhe;
   Float_t         weight_geninfo;
   Float_t         weight_lhe;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_ls;   //!
   TBranch        *b_evt;   //!
   TBranch        *b_photon_pt;   //!
   TBranch        *b_photon_eta;   //!
   TBranch        *b_photon_phi;   //!
   TBranch        *b_photon_mother_id;   //!
   TBranch        *b_photon_mother_status;   //!
   TBranch        *b_photon_geniso04;   //!
   TBranch        *b_jet_pt;   //!
   TBranch        *b_jet_eta;   //!
   TBranch        *b_jet_phi;   //!
   TBranch        *b_jet_energy;   //!
   TBranch        *b_parton_pt;   //!
   TBranch        *b_parton_eta;   //!
   TBranch        *b_parton_phi;   //!
   TBranch        *b_parton_id;   //!
   TBranch        *b_parton_status;   //!
   TBranch        *b_weights_geninfo;   //!
   TBranch        *b_weights_lhe;   //!
   TBranch        *b_weight_geninfo;   //!
   TBranch        *b_weight_lhe;   //!

   GenLevelAcceptance(TString filename, float xsection, float sumw, TString generator_);
   virtual ~GenLevelAcceptance();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   float xsec;
   float sumweights;

   typedef pair<uint,float> OrderPair;
   struct IndexByPt {
     bool operator()(const OrderPair& j1, const OrderPair& j2 ) const {
       return j1.second > j2.second;
     }
   };
   IndexByPt ptComparator;

   typedef enum MyGenEngine {
     kUndefined,
     kSHERPA,
     kaMCatNLO
   } GenEngine;
   GenEngine generator;

   // functions

   void ApplyPhotonEtaCuts(vector<uint> &passing);
   void ApplyPhotonGenIso(vector<uint> &passing);
   void ApplyDiPhotonCuts(vector<uint> &passing);
   void ApplyJetPtCuts(vector<uint> &passing);
   void ApplyJetEtaCuts(vector<uint> &passing);
   void CleanPhotonJetOverlap(vector<uint> &passingpho, vector<uint> &passingjet);
   vector<uint> SortObjects(vector<float> *ptvector);
   float deltaR(float eta1, float phi1, float eta2, float phi2);
   void FillOutput(vector<uint> &passingpho, vector<uint> &passingjet);
   
   // output tree
   TFile *fOutput;
   TTree *lighttree;
   void InitOutputTree();

   Float_t event_luminormfactor;
   vector<float> *event_luminormfactor_scalevar;
   vector<float> *event_luminormfactor_pdfvar;
   Float_t event_base_luminormfactor;
   Float_t event_Kfactor;
   Float_t event_weight;

   Float_t pholead_GEN_pt, photrail_GEN_pt;
   Float_t pholead_GEN_eta, photrail_GEN_eta;
   Float_t pholead_GEN_phi, photrail_GEN_phi;

   static const uint global_maxN_jets=100;
   Int_t n_GEN_jets;
   Float_t jet_GEN_pt[global_maxN_jets];
   Float_t jet_GEN_eta[global_maxN_jets];
   Float_t jet_GEN_phi[global_maxN_jets];
   Float_t jet_GEN_energy[global_maxN_jets];

   Bool_t tree_gen_in_acc;
   Bool_t tree_reco_in_acc;
   Bool_t tree_matched;

};

#endif

#ifdef GenLevelAcceptance_cxx
GenLevelAcceptance::GenLevelAcceptance(TString filename, float xsection, float sumw, TString generator_) : fChain(0), xsec(xsection), sumweights(sumw)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
  TFile *f = new TFile(filename.Data(),"read");
  TTree *tree = 0;
  f->GetObject("GenLevelTree",tree);
  assert(tree);
  generator = kUndefined;
  if (generator_=="SHERPA") generator = kSHERPA;
  else if (generator_=="aMCatNLO") generator = kaMCatNLO;
  assert (generator != kUndefined);
  if (generator==kaMCatNLO) assert (xsec==1); // set sumweights = lhe_nevents for aMCatNLO
  Init(tree);
  InitOutputTree();
  Loop();
}

GenLevelAcceptance::~GenLevelAcceptance()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

Int_t GenLevelAcceptance::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t GenLevelAcceptance::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void GenLevelAcceptance::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   photon_pt = 0;
   photon_eta = 0;
   photon_phi = 0;
   photon_mother_id = 0;
   photon_mother_status = 0;
   photon_geniso04 = 0;
   jet_pt = 0;
   jet_eta = 0;
   jet_phi = 0;
   jet_energy = 0;
   parton_pt = 0;
   parton_eta = 0;
   parton_phi = 0;
   parton_id = 0;
   parton_status = 0;
   weights_geninfo = 0;
   weights_lhe = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("ls", &ls, &b_ls);
   fChain->SetBranchAddress("evt", &evt, &b_evt);
   fChain->SetBranchAddress("photon_pt", &photon_pt, &b_photon_pt);
   fChain->SetBranchAddress("photon_eta", &photon_eta, &b_photon_eta);
   fChain->SetBranchAddress("photon_phi", &photon_phi, &b_photon_phi);
   fChain->SetBranchAddress("photon_mother_id", &photon_mother_id, &b_photon_mother_id);
   fChain->SetBranchAddress("photon_mother_status", &photon_mother_status, &b_photon_mother_status);
   fChain->SetBranchAddress("photon_geniso04", &photon_geniso04, &b_photon_geniso04);
   fChain->SetBranchAddress("jet_pt", &jet_pt, &b_jet_pt);
   fChain->SetBranchAddress("jet_eta", &jet_eta, &b_jet_eta);
   fChain->SetBranchAddress("jet_phi", &jet_phi, &b_jet_phi);
   fChain->SetBranchAddress("jet_energy", &jet_energy, &b_jet_energy);
   fChain->SetBranchAddress("parton_pt", &parton_pt, &b_parton_pt);
   fChain->SetBranchAddress("parton_eta", &parton_eta, &b_parton_eta);
   fChain->SetBranchAddress("parton_phi", &parton_phi, &b_parton_phi);
   fChain->SetBranchAddress("parton_id", &parton_id, &b_parton_id);
   fChain->SetBranchAddress("parton_status", &parton_status, &b_parton_status);
   fChain->SetBranchAddress("weights_geninfo", &weights_geninfo, &b_weights_geninfo);
   fChain->SetBranchAddress("weights_lhe", &weights_lhe, &b_weights_lhe);
   fChain->SetBranchAddress("weight_geninfo", &weight_geninfo, &b_weight_geninfo);
   fChain->SetBranchAddress("weight_lhe", &weight_lhe, &b_weight_lhe);
   Notify();
}

Bool_t GenLevelAcceptance::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void GenLevelAcceptance::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t GenLevelAcceptance::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
  return (entry>-1);
}

void GenLevelAcceptance::InitOutputTree(){

  fOutput = new TFile("outfile.root","recreate");
  fOutput->cd();
  lighttree = new TTree("LightTreeGenReco_Default","LightTreeGenReco_Default");

  event_luminormfactor_scalevar = new vector<float>();
  event_luminormfactor_pdfvar   = new vector<float>();

  lighttree->Branch("event_luminormfactor",&event_luminormfactor,"event_luminormfactor/F");
  lighttree->Branch("event_luminormfactor_scalevar",&event_luminormfactor_scalevar);
  lighttree->Branch("event_luminormfactor_pdfvar",&event_luminormfactor_pdfvar);
  lighttree->Branch("event_Kfactor",&event_Kfactor,"event_Kfactor/F");
  lighttree->Branch("event_weight",&event_weight,"event_weight/F");
  lighttree->Branch("pholead_GEN_pt",&pholead_GEN_pt,"pholead_GEN_pt/F");
  lighttree->Branch("photrail_GEN_pt",&photrail_GEN_pt,"photrail_GEN_pt/F");
  lighttree->Branch("pholead_GEN_eta",&pholead_GEN_eta,"pholead_GEN_eta/F");
  lighttree->Branch("photrail_GEN_eta",&photrail_GEN_eta,"photrail_GEN_eta/F");
  lighttree->Branch("pholead_GEN_phi",&pholead_GEN_phi,"pholead_GEN_phi/F");
  lighttree->Branch("photrail_GEN_phi",&photrail_GEN_phi,"photrail_GEN_phi/F");
  lighttree->Branch("n_GEN_jets",&n_GEN_jets,"n_GEN_jets/I");
  lighttree->Branch("jet_GEN_pt",&jet_GEN_pt,"jet_GEN_pt[n_GEN_jets]/F");
  lighttree->Branch("jet_GEN_eta",&jet_GEN_eta,"jet_GEN_eta[n_GEN_jets]/F");
  lighttree->Branch("jet_GEN_phi",&jet_GEN_phi,"jet_GEN_phi[n_GEN_jets]/F");
  lighttree->Branch("jet_GEN_energy",&jet_GEN_energy,"jet_GEN_energy[n_GEN_jets]/F");
  lighttree->Branch("gen_in_acc",&tree_gen_in_acc,"gen_in_acc/O");
  lighttree->Branch("reco_in_acc",&tree_reco_in_acc,"reco_in_acc/O");
  lighttree->Branch("matched",&tree_matched,"matched/O");

  tree_reco_in_acc=false;
  tree_matched=false;

  event_base_luminormfactor = xsec*1e3/sumweights;
  event_Kfactor = 1;
  event_weight = 1;

}

void GenLevelAcceptance::FillOutput(vector<uint> &passingpho, vector<uint> &passingjet){

  event_luminormfactor_scalevar->clear();
  event_luminormfactor_pdfvar->clear();

  if (generator==kSHERPA){
    event_luminormfactor = event_base_luminormfactor*weight_geninfo;
  }
  else if (generator==kaMCatNLO){
    event_luminormfactor = event_base_luminormfactor*weight_lhe;
    for (uint i=0; i<8; i++) event_luminormfactor_scalevar->push_back(event_base_luminormfactor*weights_lhe->at(i+1));
    for (uint i=0; i<52; i++) event_luminormfactor_pdfvar->push_back(event_base_luminormfactor*weights_lhe->at(i+9));
  }
  else assert(false);

  tree_gen_in_acc = (passingpho.size()>=2);

  pholead_GEN_pt =   tree_gen_in_acc ? photon_pt->at(passingpho.at(0)) : -999;
  pholead_GEN_eta =  tree_gen_in_acc ? photon_eta->at(passingpho.at(0)) : -999;
  pholead_GEN_phi =  tree_gen_in_acc ? photon_phi->at(passingpho.at(0)) : -999;
  photrail_GEN_pt =  tree_gen_in_acc ? photon_pt->at(passingpho.at(1)) : -999;
  photrail_GEN_eta = tree_gen_in_acc ? photon_eta->at(passingpho.at(1)) : -999;
  photrail_GEN_phi = tree_gen_in_acc ? photon_phi->at(passingpho.at(1)) : -999;

  n_GEN_jets = 0;
  uint index=0;
  for (vector<uint>::const_iterator it = passingjet.begin(); it != passingjet.end() && index<global_maxN_jets; it++){
    jet_GEN_pt[index]=jet_pt->at(*it);
    jet_GEN_eta[index]=jet_eta->at(*it);
    jet_GEN_phi[index]=jet_phi->at(*it);
    jet_GEN_energy[index]=jet_energy->at(*it);
    n_GEN_jets++;
    index++;
  }
  while (index<global_maxN_jets) {
    jet_GEN_pt[index] =-999;
    jet_GEN_eta[index]=-999;
    jet_GEN_phi[index]=-999;
    jet_GEN_energy[index]=-999;
    index++;
  }

  lighttree->Fill();

}

float GenLevelAcceptance::deltaR(float eta1, float phi1, float eta2, float phi2){
  float deta = eta1-eta2;
  float dphi = phi1-phi2;
  while (dphi>pi) dphi-=2*pi;
  while (dphi<=-pi) dphi+=2*pi;
  return sqrt(deta*deta+dphi*dphi);
}

vector<uint> GenLevelAcceptance::SortObjects(vector<float> *ptvector){
  vector<pair<uint,float> > ordering;
  vector<uint> ordered;
  for (uint i=0; i<ptvector->size(); i++) ordering.push_back(pair<int,float>(i,ptvector->at(i)));
  sort(ordering.begin(),ordering.end(),ptComparator);
  for (uint i=0; i<ordering.size(); i++) ordered.push_back(ordering.at(i).first);
  return ordered;
}




#endif // #ifdef GenLevelAcceptance_cxx
