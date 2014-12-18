#define GenLevelAcceptance_cxx
#include "GenLevelAcceptance.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void GenLevelAcceptance::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L GenLevelAcceptance.C
//      Root > GenLevelAcceptance t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      if (jentry%100000==0) cout << "Processing entry " << jentry << endl;

      vector<uint> genphotons = SortObjects(photon_pt);
      vector<uint> genjets = SortObjects(jet_pt);

      ApplyPhotonEtaCuts(genphotons);
      ApplyPhotonGenIso(genphotons);

      ApplyJetPtCuts(genjets);
      ApplyJetEtaCuts(genjets);
      
      ApplyDiPhotonCuts(genphotons);
      CleanPhotonJetOverlap(genphotons,genjets);

      FillOutput(genphotons,genjets);

   }

   lighttree->Write();
   fOutput->Close();


}





void GenLevelAcceptance::ApplyPhotonEtaCuts(vector<uint> &passing){
  for (vector<uint>::iterator it = passing.begin(); it != passing.end(); ){
#ifdef _DEBUG
    cout << "ApplyPhotonEtaCuts " << *it << " " << photon_eta->at(*it) << endl;
#endif
    float eta = fabs(photon_eta->at(*it));   
    bool pass = ((eta<1.4442) || (eta>1.566 && eta<2.5));
    if (!pass) it = passing.erase(it); else it++;
  }
}
void GenLevelAcceptance::ApplyPhotonGenIso(vector<uint> &passing){
  for (vector<uint>::iterator it = passing.begin(); it != passing.end(); ){
#ifdef _DEBUG
  cout << "ApplyPhotonGenIso " << *it << " " << photon_geniso04->at(*it) << endl;
#endif
    bool pass = (photon_geniso04->at(*it)<5);
    if (!pass) it = passing.erase(it); else it++;
  }
}
void GenLevelAcceptance::ApplyDiPhotonCuts(vector<uint> &passing){
  for (int i=0; i<int(passing.size())-1; i++) assert(photon_pt->at(passing.at(i))>=photon_pt->at(passing.at(i+1)));
#ifdef _DEBUG
  cout << "ApplyDiPhotonCuts" << endl;
#endif
  if (passing.size()>=2 && photon_pt->at(passing.at(0))>40 && photon_pt->at(passing.at(1))>25) passing.resize(2);
  else passing.clear();
  if (passing.size()>=2 && deltaR(photon_eta->at(passing.at(0)),photon_phi->at(passing.at(0)),photon_eta->at(passing.at(1)),photon_phi->at(passing.at(1)))<0.45) passing.clear();
}


void GenLevelAcceptance::ApplyJetPtCuts(vector<uint> &passing){
  for (vector<uint>::iterator it = passing.begin(); it != passing.end(); ){
#ifdef _DEBUG
    cout << "ApplyJetPtCuts " << *it << " " << jet_pt->at(*it) << endl;
#endif
    if (jet_pt->at(*it)<25) it = passing.erase(it); else it++;
  }
}
void GenLevelAcceptance::ApplyJetEtaCuts(vector<uint> &passing){
  for (vector<uint>::iterator it = passing.begin(); it != passing.end(); ){
#ifdef _DEBUG
    cout << "ApplyJetEtaCuts " << *it << " " << jet_eta->at(*it) << endl;
#endif
    float eta = fabs(jet_eta->at(*it));   
    bool pass = (eta<4.7);
    if (!pass) it = passing.erase(it); else it++;
  }
}


void GenLevelAcceptance::CleanPhotonJetOverlap(vector<uint> &passingpho, vector<uint> &passingjet){
  for (vector<uint>::iterator it = passingjet.begin(); it != passingjet.end(); ){
    bool bad = false;
    for (vector<uint>::iterator itpho = passingpho.begin(); itpho != passingpho.end() && !bad; itpho++){
      if (deltaR(jet_eta->at(*it),jet_phi->at(*it),photon_eta->at(*itpho),photon_phi->at(*itpho))<0.5) bad=true;
    }
    if (bad) it = passingjet.erase(it); else it++;
  }
#ifdef _DEBUG
  cout << "done CleanPhotonJetOverlap" << endl;
#endif

}


