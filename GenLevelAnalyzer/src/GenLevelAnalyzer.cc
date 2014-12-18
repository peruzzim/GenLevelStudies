// -*- C++ -*-
//
// Package:    GenLevelAnalyzer
// Class:      GenLevelAnalyzer
// 
/**\class GenLevelAnalyzer GenLevelAnalyzer.cc GenLevelStudies/GenLevelAnalyzer/src/GenLevelAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Marco Peruzzi (ETHZ) [peruzzi]
//         Created:  Wed Dec 17 10:44:49 CET 2014
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "TLorentzVector.h"
#include "TVector.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include <vector>
#include <algorithm>

using namespace edm;
using namespace std;
using namespace reco;

//
// class declaration
//

class GenLevelAnalyzer : public edm::EDAnalyzer {
   public:
      explicit GenLevelAnalyzer(const edm::ParameterSet&);
      ~GenLevelAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

  void FillPhoton(GenParticleCollection::const_iterator p);
  void FillJet(GenJetCollection::const_iterator j);
  void FillParton(GenParticleCollection::const_iterator p);
  void FillWeight(const GenEventInfoProduct *info);
  void FillTree();
  void ClearEvent();
  void Error(TString err);

      // ----------member data ---------------------------

  std::string OutputFile_;
  TFile *fOutput;
  TTree *fTree;

  edm::InputTag genParticlesTag_;
  edm::InputTag genJetsTag_;
  double minPtPhotons_;
  double maxEtaPhotons_;
  double minPtJets_;
  double maxEtaJets_;
  vector<int> partonStatusList_;
  bool writeAllGenParticles;

  ULong64_t run;
  ULong64_t ls;
  ULong64_t evt;

  vector<float> *photon_pt = new vector<float>();
  vector<float> *photon_eta = new vector<float>();
  vector<float> *photon_phi = new vector<float>();
  vector<int>   *photon_mother_id = new vector<int>();
  vector<int>   *photon_mother_status = new vector<int>();
  vector<float> *photon_geniso04 = new vector<float>();
  vector<float> *jet_pt = new vector<float>();
  vector<float> *jet_eta = new vector<float>();
  vector<float> *jet_phi = new vector<float>();
  vector<float> *parton_pt = new vector<float>();
  vector<float> *parton_eta = new vector<float>();
  vector<float> *parton_phi = new vector<float>();
  vector<int>   *parton_id = new vector<int>();
  vector<int>   *parton_status = new vector<int>();
  vector<float> *weights = new vector<float>();

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
GenLevelAnalyzer::GenLevelAnalyzer(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed

  genParticlesTag_ = iConfig.getParameter<edm::InputTag>("genParticlesTag");
  genJetsTag_ = iConfig.getParameter<edm::InputTag>("genJetsTag");
  OutputFile_ = iConfig.getParameter<std::string>("OutputFile");

  minPtPhotons_ = iConfig.getParameter<double>("minPtPhotons");
  maxEtaPhotons_ = iConfig.getParameter<double>("maxEtaPhotons");
  minPtJets_ = iConfig.getParameter<double>("minPtJets");
  maxEtaJets_ = iConfig.getParameter<double>("maxEtaJets");
  partonStatusList_ = iConfig.getParameter<vector<int> >("partonStatusList");
  writeAllGenParticles = iConfig.getParameter<bool>("writeAllGenParticles");

  for (vector<int>::const_iterator it = partonStatusList_.begin(); it!=partonStatusList_.end(); it++){
    cout << "Particles with status " << *it << " will be inserted in the parton collection" << endl;
  }

  fOutput = new TFile(OutputFile_.c_str(),"recreate");
  fOutput->cd();
  fTree = new TTree("GenLevelTree","GenLevelTree");

  fTree->Branch("run",&run,"run/l");
  fTree->Branch("ls",&ls,"ls/l");
  fTree->Branch("evt",&evt,"evt/l");

  fTree->Branch("photon_pt",&photon_pt);
  fTree->Branch("photon_eta",&photon_eta);
  fTree->Branch("photon_phi",&photon_phi);
  fTree->Branch("photon_mother_id",&photon_mother_id);
  fTree->Branch("photon_mother_status",&photon_mother_status);
  fTree->Branch("photon_geniso04",&photon_geniso04);
  fTree->Branch("jet_pt",&jet_pt);
  fTree->Branch("jet_eta",&jet_eta);
  fTree->Branch("jet_phi",&jet_phi);
  fTree->Branch("parton_pt",&parton_pt);
  fTree->Branch("parton_eta",&parton_eta);
  fTree->Branch("parton_phi",&parton_phi);
  fTree->Branch("parton_id",&parton_id);
  fTree->Branch("parton_status",&parton_status);
  fTree->Branch("weights",&weights);

}


GenLevelAnalyzer::~GenLevelAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

  fOutput->Write();
  fOutput->Close();

}


//
// member functions
//

// ------------ method called for each event  ------------
void
GenLevelAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  run = iEvent.id().run();
  ls = iEvent.luminosityBlock();
  evt = iEvent.id().event();
  
  ClearEvent();

   edm::Handle<GenParticleCollection> genParticles;
   iEvent.getByLabel( genParticlesTag_, genParticles );
   const GenParticleCollection *parts = genParticles.product();

   edm::Handle<GenJetCollection> genJets;
   iEvent.getByLabel( genJetsTag_, genJets);

   edm::Handle<GenEventInfoProduct> genEvtInfo;
   iEvent.getByLabel("generator", genEvtInfo);

   int index=-1;
   for (GenParticleCollection::const_iterator p = genParticles->begin(); p!=genParticles->end(); p++){
     index++;
     if (writeAllGenParticles || find(partonStatusList_.begin(),partonStatusList_.end(),p->status())!=partonStatusList_.end()){
       FillParton(p);
     }
     if (p->status()==1){
       if (p->pdgId()!=22) continue;
       if (p->pt()<minPtPhotons_) continue;
       if (fabs(p->eta())>maxEtaPhotons_) continue;
       FillPhoton(p);
       float sum = 0;
       for (uint i=0; i<parts->size(); i++){
	 const GenParticle &p2 = parts->at(i);
	 if (p2.status()!=1) continue;
	 if (i==(uint)index) continue;
	 if (reco::deltaR(p->eta(),p->phi(),p2.eta(),p2.phi())<0.4){
	   sum+=p2.et(); // should probably be pt(), keep et() for consistency with NTupleProducer
	 }
       }
       photon_geniso04->push_back(sum);
     }
   }

   for (GenJetCollection::const_iterator j = genJets->begin(); j!=genJets->end(); j++){
     if (j->pt()<minPtJets_) continue;
     if (fabs(j->eta())>maxEtaJets_) continue;
     FillJet(j);
   }

   FillWeight(genEvtInfo.product());
   FillTree();

}

void GenLevelAnalyzer::FillPhoton(GenParticleCollection::const_iterator p){

  if (p->pdgId()!=22) {
    Error(Form("ERROR: Passed particle with pdgid %d to photon filler!",p->pdgId()));
    return;
  }
  if (p->status()!=1) {
    Error(Form("ERROR: Passed particle with status %d to photon filler!",p->status()));
    return;
  }

  photon_pt->push_back(p->pt());
  photon_eta->push_back(p->eta());
  photon_phi->push_back(p->phi());
  photon_mother_id->push_back(-999);
  photon_mother_status->push_back(-999);
  if (p->mother()){
    photon_mother_id->back()=p->mother()->pdgId();
    photon_mother_status->back()=p->mother()->status();
  }

}

void GenLevelAnalyzer::FillJet(GenJetCollection::const_iterator j){
  jet_pt->push_back(j->pt());
  jet_eta->push_back(j->eta());
  jet_phi->push_back(j->phi());
}

void GenLevelAnalyzer::FillParton(GenParticleCollection::const_iterator p){
  if (fabs(p->eta())>10) return;
  parton_pt->push_back(p->pt());
  parton_eta->push_back(p->eta());
  parton_phi->push_back(p->phi());
  parton_id->push_back(p->pdgId());
  parton_status->push_back(p->status());
}

void GenLevelAnalyzer::FillWeight(const GenEventInfoProduct *info){
  weights->push_back(info->weight());
};

void GenLevelAnalyzer::FillTree(){
  fTree->Fill();
}

void GenLevelAnalyzer::ClearEvent(){
  photon_pt->clear();
  photon_eta->clear();
  photon_phi->clear();
  photon_mother_id->clear();
  photon_mother_status->clear();
  photon_geniso04->clear();
  jet_pt->clear();
  jet_eta->clear();
  jet_phi->clear();
  parton_pt->clear();
  parton_eta->clear();
  parton_phi->clear();
  parton_id->clear();
  parton_status->clear();
  weights->clear();
}

void GenLevelAnalyzer::Error(TString err){
  cout << err.Data() << endl;
}

// ------------ method called once each job just before starting event loop  ------------
void 
GenLevelAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
GenLevelAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
GenLevelAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
GenLevelAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
GenLevelAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
GenLevelAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
GenLevelAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(GenLevelAnalyzer);
