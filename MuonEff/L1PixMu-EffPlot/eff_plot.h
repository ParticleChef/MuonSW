//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Oct  8 15:09:43 2018 by ROOT version 6.10/09
// from TTree t/t
// found on file: ../test.root
//////////////////////////////////////////////////////////

#ifndef eff_plot_h
#define eff_plot_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"

class eff_plot {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           totalEvent;
   Float_t         totalme0N;
   Int_t           ntnME02;
   vector<float>   *ntME0Et;
   vector<float>   *ntME0Eta;
   vector<float>   *ntME0Phi;
   vector<int>     *PiXTRKbit;
   vector<int>     *trigger_bit_width;
   vector<int>     *pix_comb;
   Int_t           nPix123_segments;
   Int_t           nPix124_segments;
   Int_t           nPix134_segments;
   Int_t           nPix234_segments;
   vector<bool>    *ntCl_match;
   vector<bool>    *isTrack_match;
   vector<float>   *chi2;
   vector<float>   *track_dr;
   vector<bool>    *withoutEM_match;
   vector<bool>    *withEM_match;
   vector<int>     *ntfirstPix;
   vector<int>     *ntsecondPix;
   vector<int>     *ntthirdPix;
   vector<int>     *ntfourthPix;
   Float_t         nt_genPhi;
   Float_t         nt_genEta;
   Float_t         nt_genPt;
   Float_t         matchedME0Et;
   Float_t         matchedME0Eta;
   Float_t         matchedME0Phi;
   Float_t         nt_lastSimtkpt;
   Float_t         nt_initialSimtkpt;
   Int_t           fired;

   // List of branches
   TBranch        *b_count_Entry;   //!
   TBranch        *b_me0N;   //!
   TBranch        *b_ntnME02;   //!
   TBranch        *b_ntME0Et;   //!
   TBranch        *b_ntME0Eta;   //!
   TBranch        *b_ntME0Phi;   //!
   TBranch        *b_PiXTRKbit;   //!
   TBranch        *b_trigger_bit_width;   //!
   TBranch        *b_pix_comb;   //!
   TBranch        *b_nPix123_segments;   //!
   TBranch        *b_nPix124_segments;   //!
   TBranch        *b_nPix134_segments;   //!
   TBranch        *b_nPix234_segments;   //!
   TBranch        *b_ntCl_match;   //!
   TBranch        *b_isTrack_match;   //!
   TBranch        *b_chi2;   //!
   TBranch        *b_track_dr;   //!
   TBranch        *b_withoutEM_match;   //!
   TBranch        *b_withEM_match;   //!
   TBranch        *b_ntfirstPix;   //!
   TBranch        *b_ntsecondPix;   //!
   TBranch        *b_ntthirdPix;   //!
   TBranch        *b_ntfourthPix;   //!
   TBranch        *b_nt_genPhi;   //!
   TBranch        *b_nt_genEta;   //!
   TBranch        *b_nt_genPt;   //!
   TBranch        *b_matchedME0Et;   //!
   TBranch        *b_matchedME0Eta;   //!
   TBranch        *b_matchedME0Phi;   //!
   TBranch        *b_nt_lastSimtkpt;   //!
   TBranch        *b_nt_initialSimtkpt;   //!
   TBranch        *b_fired;   //!

   eff_plot(TTree *tree=0);
   virtual ~eff_plot();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef eff_plot_cxx
eff_plot::eff_plot(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../test_200PU.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("../test_200PU.root");
      }
      f->GetObject("t",tree);

   }
   Init(tree);
}

eff_plot::~eff_plot()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t eff_plot::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t eff_plot::LoadTree(Long64_t entry)
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

void eff_plot::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   ntME0Et = 0;
   ntME0Eta = 0;
   ntME0Phi = 0;
   PiXTRKbit = 0;
   trigger_bit_width = 0;
   pix_comb = 0;
   ntCl_match = 0;
   isTrack_match = 0;
   chi2 = 0;
   track_dr = 0;
   withoutEM_match = 0;
   withEM_match = 0;
   ntfirstPix = 0;
   ntsecondPix = 0;
   ntthirdPix = 0;
   ntfourthPix = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("totalEvent", &totalEvent, &b_count_Entry);
   fChain->SetBranchAddress("totalme0N", &totalme0N, &b_me0N);
   fChain->SetBranchAddress("ntnME02", &ntnME02, &b_ntnME02);
   fChain->SetBranchAddress("ntME0Et", &ntME0Et, &b_ntME0Et);
   fChain->SetBranchAddress("ntME0Eta", &ntME0Eta, &b_ntME0Eta);
   fChain->SetBranchAddress("ntME0Phi", &ntME0Phi, &b_ntME0Phi);
   fChain->SetBranchAddress("PiXTRKbit", &PiXTRKbit, &b_PiXTRKbit);
   fChain->SetBranchAddress("trigger_bit_width", &trigger_bit_width, &b_trigger_bit_width);
   fChain->SetBranchAddress("pix_comb", &pix_comb, &b_pix_comb);
   fChain->SetBranchAddress("nPix123_segments", &nPix123_segments, &b_nPix123_segments);
   fChain->SetBranchAddress("nPix124_segments", &nPix124_segments, &b_nPix124_segments);
   fChain->SetBranchAddress("nPix134_segments", &nPix134_segments, &b_nPix134_segments);
   fChain->SetBranchAddress("nPix234_segments", &nPix234_segments, &b_nPix234_segments);
   fChain->SetBranchAddress("ntCl_match", &ntCl_match, &b_ntCl_match);
   fChain->SetBranchAddress("isTrack_match", &isTrack_match, &b_isTrack_match);
   fChain->SetBranchAddress("chi2", &chi2, &b_chi2);
   fChain->SetBranchAddress("track_dr", &track_dr, &b_track_dr);
   fChain->SetBranchAddress("withoutEM_match", &withoutEM_match, &b_withoutEM_match);
   fChain->SetBranchAddress("withEM_match", &withEM_match, &b_withEM_match);
   fChain->SetBranchAddress("ntfirstPix", &ntfirstPix, &b_ntfirstPix);
   fChain->SetBranchAddress("ntsecondPix", &ntsecondPix, &b_ntsecondPix);
   fChain->SetBranchAddress("ntthirdPix", &ntthirdPix, &b_ntthirdPix);
   fChain->SetBranchAddress("ntfourthPix", &ntfourthPix, &b_ntfourthPix);
   fChain->SetBranchAddress("nt_genPhi", &nt_genPhi, &b_nt_genPhi);
   fChain->SetBranchAddress("nt_genEta", &nt_genEta, &b_nt_genEta);
   fChain->SetBranchAddress("nt_genPt", &nt_genPt, &b_nt_genPt);
   fChain->SetBranchAddress("matchedME0Et", &matchedME0Et, &b_matchedME0Et);
   fChain->SetBranchAddress("matchedME0Eta", &matchedME0Eta, &b_matchedME0Eta);
   fChain->SetBranchAddress("matchedME0Phi", &matchedME0Phi, &b_matchedME0Phi);
   fChain->SetBranchAddress("nt_lastSimtkpt", &nt_lastSimtkpt, &b_nt_lastSimtkpt);
   fChain->SetBranchAddress("nt_initialSimtkpt", &nt_initialSimtkpt, &b_nt_initialSimtkpt);
   fChain->SetBranchAddress("fired", &fired, &b_fired);
   Notify();
}

Bool_t eff_plot::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void eff_plot::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t eff_plot::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef eff_plot_cxx
