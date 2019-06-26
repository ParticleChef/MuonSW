//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Jan  3 21:43:06 2019 by ROOT version 6.10/09
// from TTree L1PiXTRKTree/L1PiXTRKTree
// found on file: /xrootd/store/user/jhong/SingleMu_Pt2to200_Eta3p0_CMSSW_9_3_7_NoPU_D17test/crab_Muon0802/180802_123200/0000/SingleMuonNoPU_1.root
//////////////////////////////////////////////////////////

#ifndef plot_h
#define plot_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"
#include "vector"
#include "vector"

class plot {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           nVtx;
   Int_t           nMeanPU;
   Int_t           genPartN;
   vector<float>   *genPartE;
   vector<float>   *genPartPt;
   vector<float>   *genPartEta;
   vector<float>   *genPartPhi;
   vector<int>     *genPartCharge;
   vector<int>     *genPartId;
   vector<float>   *propgenElPartE;
   vector<float>   *propgenElPartPt;
   vector<float>   *propgenElPartEta;
   vector<float>   *propgenElPartPhi;
   vector<int>     *propgenElPartCharge;
   vector<float>   *propgenElPartx;
   vector<float>   *propgenElParty;
   vector<float>   *propgenElPartz;
   Int_t           simTrkN;
   vector<float>   *simTrkPt;
   vector<float>   *simTrkEta;
   vector<float>   *simTrkPhi;
   vector<int>     *simTrkId;
   vector<int>     *simTrkType;
   vector<float>   *simTrkVx;
   vector<float>   *simTrkVy;
   vector<float>   *simTrkVz;
   vector<float>   *simVx;
   vector<float>   *simVy;
   vector<float>   *simVz;
   Float_t         lastSimtkpt;
   Float_t         initialSimtkpt;
   Int_t           bremflag;
   vector<float>   *Brempos_radius;
   vector<float>   *Brem_eLoss;
   vector<float>   *Brem_ptLoss;
   vector<float>   *Brempos_x;
   vector<float>   *Brempos_y;
   vector<float>   *Brempos_z;
   vector<float>   *propgenPoPartE;
   vector<float>   *propgenPoPartPt;
   vector<float>   *propgenPoPartEta;
   vector<float>   *propgenPoPartPhi;
   vector<int>     *propgenPoPartCharge;
   vector<float>   *propgenPoPartx;
   vector<float>   *propgenPoParty;
   vector<float>   *propgenPoPartz;
   vector<int>     *bRecHitLayer;
   vector<int>     *bRecHitLadder;
   vector<int>     *bRecHitModule;
   vector<int>     *fRecHitDisk;
   vector<int>     *fRecHitBlade;
   vector<int>     *fRecHitSide;
   vector<int>     *fRecHitPanel;
   vector<int>     *fRecHitModule;
   Int_t           bRecHitN;
   Int_t           fRecHitN;
   vector<float>   *fRecHitGx;
   vector<float>   *fRecHitGy;
   vector<float>   *fRecHitGz;
   vector<float>   *fRhSize;
   vector<float>   *fRhSizeX;
   vector<float>   *fRhSizeY;
   vector<float>   *bRecHitGx;
   vector<float>   *bRecHitGy;
   vector<float>   *bRecHitGz;
   vector<float>   *bRhSize;
   vector<float>   *bRhSizeX;
   vector<float>   *bRhSizeY;
   Int_t           bfastsimHitN;
   Int_t           ffastsimHitN;
   vector<int>     *bfastsimHitLayer;
   vector<float>   *bfastsimHitGx;
   vector<float>   *bfastsimHitGy;
   vector<float>   *bfastsimHitGz;
   vector<int>     *ffastsimHitLayer;
   vector<float>   *ffastsimHitGx;
   vector<float>   *ffastsimHitGy;
   vector<float>   *ffastsimHitGz;
   Int_t           egCrysN;
   vector<float>   *egCrysE;
   vector<float>   *egCrysEt;
   vector<float>   *egCrysEta;
   vector<float>   *egCrysPhi;
   vector<float>   *egCrysGx;
   vector<float>   *egCrysGy;
   vector<float>   *egCrysGz;
   Int_t           egCrysClusterN;
   vector<float>   *egCrysClusterE;
   vector<float>   *egCrysClusterEt;
   vector<float>   *egCrysClusterEta;
   vector<float>   *egCrysClusterPhi;
   vector<float>   *egCrysClusterGx;
   vector<float>   *egCrysClusterGy;
   vector<float>   *egCrysClusterGz;
   vector<float>   *egCrysClusterPGx;
   vector<float>   *egCrysClusterPGy;
   vector<float>   *egCrysClusterPGz;
   vector<bool>    *isTrackMatched;
   vector<float>   *isoConeNTrack;
   vector<float>   *isoConePtTrack;
   vector<float>   *trackHighestPt;
   vector<float>   *trackHighestPtEta;
   vector<float>   *trackHighestPtPhi;
   vector<float>   *trackHighestPtChi2;
   vector<float>   *trackHighestPtCutChi2;
   vector<float>   *trackHighestPtCutChi2Eta;
   vector<float>   *trackHighestPtCutChi2Phi;
   vector<float>   *trackHighestPtCutChi2Chi2;
   vector<float>   *trackmatchingdR;
   vector<bool>    *hgcal_isTrackMatched;
   vector<float>   *hgcal_isoConeNTrack;
   vector<float>   *hgcal_isoConePtTrack;
   vector<float>   *hgcal_trackHighestPt;
   vector<float>   *hgcal_trackHighestPtEta;
   vector<float>   *hgcal_trackHighestPtPhi;
   vector<float>   *hgcal_trackHighestPtChi2;
   vector<float>   *hgcal_trackHighestPtCutChi2;
   vector<float>   *hgcal_trackHighestPtCutChi2Eta;
   vector<float>   *hgcal_trackHighestPtCutChi2Phi;
   vector<float>   *hgcal_trackHighestPtCutChi2Chi2;
   vector<float>   *hgcal_trackmatchingdR;
   Int_t           cl3d_n;
   vector<float>   *cl3d_pt;
   vector<float>   *cl3d_energy;
   vector<float>   *cl3d_eta;
   vector<float>   *cl3d_phi;
   vector<int>     *cl3d_nclu;
   vector<float>   *cl3d_x;
   vector<float>   *cl3d_y;
   vector<int>     *cl3d_z;
   vector<float>   *cl3d_hovere;
   vector<int>     *cl3d_showerlength;
   vector<int>     *cl3d_coreshowerlength;
   vector<int>     *cl3d_firstlayer;
   vector<int>     *cl3d_maxlayer;
   vector<float>   *cl3d_seetot;
   vector<float>   *cl3d_seemax;
   vector<float>   *cl3d_spptot;
   vector<float>   *cl3d_sppmax;
   vector<float>   *cl3d_szz;
   vector<float>   *cl3d_srrtot;
   vector<float>   *cl3d_srrmax;
   vector<float>   *cl3d_srrmean;
   vector<float>   *cl3d_emaxe;
   UShort_t        egN;
   vector<float>   *egEt;
   vector<float>   *egEta;
   vector<float>   *egPhi;
   vector<float>   *egGx;
   vector<float>   *egGy;
   vector<float>   *egGz;
   vector<short>   *egIEt;
   vector<short>   *egIEta;
   vector<short>   *egIPhi;
   vector<short>   *egIso;
   vector<short>   *egBx;
   vector<short>   *egTowerIPhi;
   vector<short>   *egTowerIEta;
   vector<short>   *egRawEt;
   vector<short>   *egIsoEt;
   vector<short>   *egFootprintEt;
   vector<short>   *egNTT;
   vector<short>   *egShape;
   vector<short>   *egTowerHoE;
   UInt_t          me0SegNum;
   vector<unsigned int> *me0SegDetId;
   vector<float>   *me0SegPosX;
   vector<float>   *me0SegPosY;
   vector<float>   *me0SegPosZ;
   vector<float>   *me0SegDirX;
   vector<float>   *me0SegDirY;
   vector<float>   *me0SegDirZ;
   vector<int>     *me0SegNumRecHit;
   vector<float>   *me0SegDeltaPhi;

   // List of branches
   TBranch        *b_nVtx;   //!
   TBranch        *b_nMeanPU;   //!
   TBranch        *b_genPartN;   //!
   TBranch        *b_genPartE;   //!
   TBranch        *b_genPartPt;   //!
   TBranch        *b_genPartEta;   //!
   TBranch        *b_genPartPhi;   //!
   TBranch        *b_genPartCharge;   //!
   TBranch        *b_genPartId;   //!
   TBranch        *b_propgenElPartE;   //!
   TBranch        *b_propgenElPartPt;   //!
   TBranch        *b_propgenElPartEta;   //!
   TBranch        *b_propgenElPartPhi;   //!
   TBranch        *b_propgenElPartCharge;   //!
   TBranch        *b_propgenElPartx;   //!
   TBranch        *b_propgenElParty;   //!
   TBranch        *b_propgenElPartz;   //!
   TBranch        *b_simTrkN;   //!
   TBranch        *b_simTrkPt;   //!
   TBranch        *b_simTrkEta;   //!
   TBranch        *b_simTrkPhi;   //!
   TBranch        *b_simTrkId;   //!
   TBranch        *b_simTrkType;   //!
   TBranch        *b_simTrkVx;   //!
   TBranch        *b_simTrkVy;   //!
   TBranch        *b_simTrkVz;   //!
   TBranch        *b_simVx;   //!
   TBranch        *b_simVy;   //!
   TBranch        *b_simVz;   //!
   TBranch        *b_lastSimtkpt;   //!
   TBranch        *b_initialSimtkpt;   //!
   TBranch        *b_bremflag;   //!
   TBranch        *b_Brempos_radius;   //!
   TBranch        *b_Brem_eLoss;   //!
   TBranch        *b_Brem_ptLoss;   //!
   TBranch        *b_Brempos_x;   //!
   TBranch        *b_Brempos_y;   //!
   TBranch        *b_Brempos_z;   //!
   TBranch        *b_propgenPoPartE;   //!
   TBranch        *b_propgenPoPartPt;   //!
   TBranch        *b_propgenPoPartEta;   //!
   TBranch        *b_propgenPoPartPhi;   //!
   TBranch        *b_propgenPoPartCharge;   //!
   TBranch        *b_propgenPoPartx;   //!
   TBranch        *b_propgenPoParty;   //!
   TBranch        *b_propgenPoPartz;   //!
   TBranch        *b_bRecHitLayer;   //!
   TBranch        *b_bRecHitLadder;   //!
   TBranch        *b_bRecHitModule;   //!
   TBranch        *b_fRecHitDisk;   //!
   TBranch        *b_fRecHitBlade;   //!
   TBranch        *b_fRecHitSide;   //!
   TBranch        *b_fRecHitPanel;   //!
   TBranch        *b_fRecHitModule;   //!
   TBranch        *b_bRecHitN;   //!
   TBranch        *b_fRecHitN;   //!
   TBranch        *b_fRecHitGx;   //!
   TBranch        *b_fRecHitGy;   //!
   TBranch        *b_fRecHitGz;   //!
   TBranch        *b_fRhSize;   //!
   TBranch        *b_fRhSizeX;   //!
   TBranch        *b_fRhSizeY;   //!
   TBranch        *b_bRecHitGx;   //!
   TBranch        *b_bRecHitGy;   //!
   TBranch        *b_bRecHitGz;   //!
   TBranch        *b_bRhSize;   //!
   TBranch        *b_bRhSizeX;   //!
   TBranch        *b_bRhSizeY;   //!
   TBranch        *b_bfastsimHitN;   //!
   TBranch        *b_ffastsimHitN;   //!
   TBranch        *b_bfastsimHitLayer;   //!
   TBranch        *b_bfastsimHitGx;   //!
   TBranch        *b_bfastsimHitGy;   //!
   TBranch        *b_bfastsimHitGz;   //!
   TBranch        *b_ffastsimHitLayer;   //!
   TBranch        *b_ffastsimHitGx;   //!
   TBranch        *b_ffastsimHitGy;   //!
   TBranch        *b_ffastsimHitGz;   //!
   TBranch        *b_egCrysN;   //!
   TBranch        *b_egCrysE;   //!
   TBranch        *b_egCrysEt;   //!
   TBranch        *b_egCrysEta;   //!
   TBranch        *b_egCrysPhi;   //!
   TBranch        *b_egCrysGx;   //!
   TBranch        *b_egCrysGy;   //!
   TBranch        *b_egCrysGz;   //!
   TBranch        *b_egCrysClusterN;   //!
   TBranch        *b_egCrysClusterE;   //!
   TBranch        *b_egCrysClusterEt;   //!
   TBranch        *b_egCrysClusterEta;   //!
   TBranch        *b_egCrysClusterPhi;   //!
   TBranch        *b_egCrysClusterGx;   //!
   TBranch        *b_egCrysClusterGy;   //!
   TBranch        *b_egCrysClusterGz;   //!
   TBranch        *b_egCrysClusterPGx;   //!
   TBranch        *b_egCrysClusterPGy;   //!
   TBranch        *b_egCrysClusterPGz;   //!
   TBranch        *b_isTrackMatched;   //!
   TBranch        *b_isoConeNTrack;   //!
   TBranch        *b_isoConePtTrack;   //!
   TBranch        *b_trackHighestPt;   //!
   TBranch        *b_trackHighestPtEta;   //!
   TBranch        *b_trackHighestPtPhi;   //!
   TBranch        *b_trackHighestPtChi2;   //!
   TBranch        *b_trackHighestPtCutChi2;   //!
   TBranch        *b_trackHighestPtCutChi2Eta;   //!
   TBranch        *b_trackHighestPtCutChi2Phi;   //!
   TBranch        *b_trackHighestPtCutChi2Chi2;   //!
   TBranch        *b_trackmatchingdR;   //!
   TBranch        *b_hgcal_isTrackMatched;   //!
   TBranch        *b_hgcal_isoConeNTrack;   //!
   TBranch        *b_hgcal_isoConePtTrack;   //!
   TBranch        *b_hgcal_trackHighestPt;   //!
   TBranch        *b_hgcal_trackHighestPtEta;   //!
   TBranch        *b_hgcal_trackHighestPtPhi;   //!
   TBranch        *b_hgcal_trackHighestPtChi2;   //!
   TBranch        *b_hgcal_trackHighestPtCutChi2;   //!
   TBranch        *b_hgcal_trackHighestPtCutChi2Eta;   //!
   TBranch        *b_hgcal_trackHighestPtCutChi2Phi;   //!
   TBranch        *b_hgcal_trackHighestPtCutChi2Chi2;   //!
   TBranch        *b_hgcal_trackmatchingdR;   //!
   TBranch        *b_cl3d_n;   //!
   TBranch        *b_cl3d_pt;   //!
   TBranch        *b_cl3d_energy;   //!
   TBranch        *b_cl3d_eta;   //!
   TBranch        *b_cl3d_phi;   //!
   TBranch        *b_cl3d_nclu;   //!
   TBranch        *b_cl3d_x;   //!
   TBranch        *b_cl3d_y;   //!
   TBranch        *b_cl3d_z;   //!
   TBranch        *b_cl3d_hovere;   //!
   TBranch        *b_cl3d_showerlength;   //!
   TBranch        *b_cl3d_coreshowerlength;   //!
   TBranch        *b_cl3d_firstlayer;   //!
   TBranch        *b_cl3d_maxlayer;   //!
   TBranch        *b_cl3d_seetot;   //!
   TBranch        *b_cl3d_seemax;   //!
   TBranch        *b_cl3d_spptot;   //!
   TBranch        *b_cl3d_sppmax;   //!
   TBranch        *b_cl3d_szz;   //!
   TBranch        *b_cl3d_srrtot;   //!
   TBranch        *b_cl3d_srrmax;   //!
   TBranch        *b_cl3d_srrmean;   //!
   TBranch        *b_cl3d_emaxe;   //!
   TBranch        *b_egN;   //!
   TBranch        *b_egEt;   //!
   TBranch        *b_egEta;   //!
   TBranch        *b_egPhi;   //!
   TBranch        *b_egGx;   //!
   TBranch        *b_egGy;   //!
   TBranch        *b_egGz;   //!
   TBranch        *b_egIEt;   //!
   TBranch        *b_egIEta;   //!
   TBranch        *b_egIPhi;   //!
   TBranch        *b_egIso;   //!
   TBranch        *b_egBx;   //!
   TBranch        *b_egTowerIPhi;   //!
   TBranch        *b_egTowerIEta;   //!
   TBranch        *b_egRawEt;   //!
   TBranch        *b_egIsoEt;   //!
   TBranch        *b_egFootprintEt;   //!
   TBranch        *b_egNTT;   //!
   TBranch        *b_egShape;   //!
   TBranch        *b_egTowerHoE;   //!
   TBranch        *b_me0SegNum;   //!
   TBranch        *b_me0SegDetId;   //!
   TBranch        *b_me0SegPosX;   //!
   TBranch        *b_me0SegPosY;   //!
   TBranch        *b_me0SegPosZ;   //!
   TBranch        *b_me0SegDirX;   //!
   TBranch        *b_me0SegDirY;   //!
   TBranch        *b_me0SegDirZ;   //!
   TBranch        *b_me0SegNumRecHit;   //!
   TBranch        *b_me0SegDeltaPhi;   //!

   plot(TTree *tree=0);
   virtual ~plot();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   inline float deltaPhi(float phi1, float phi2){
	   float result = phi1 - phi2;
	   while(result > float(M_PI)) result -= float(2*M_PI);
	   while(result <= -float(M_PI)) result += float(2*M_PI);
	   return result;
   }
};

#endif

#ifdef plot_cxx
plot::plot(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
/*
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/xrootd/store/user/jhong/SingleMu_Pt2to200_Eta3p0_CMSSW_9_3_7_NoPU_D17test/crab_Muon0802/180802_123200/0000/SingleMuonNoPU_1.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/xrootd/store/user/jhong/SingleMu_Pt2to200_Eta3p0_CMSSW_9_3_7_NoPU_D17test/crab_Muon0802/180802_123200/0000/SingleMuonNoPU_1.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("/xrootd/store/user/jhong/SingleMu_Pt2to200_Eta3p0_CMSSW_9_3_7_NoPU_D17test/crab_Muon0802/180802_123200/0000/SingleMuonNoPU_1.root:/l1PiXTRKTree");
      dir->GetObject("L1PiXTRKTree",tree);

   }
*/
   if (tree == 0) {
#ifdef SINGLE_TREE
	     TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/xrootd/store/user/jhong/SingleMu_Pt2to200_Eta3p0_CMSSW_9_3_7_NoPU_D17test/crab_Muon0802/180802_123200/*/");
	     if(!f || !f->IsOpen()){
		     f = new TFile("/xrootd/store/user/jhong/SingleMu_Pt2to200_Eta3p0_CMSSW_9_3_7_NoPU_D17test/crab_Muon0802/180802_123200/*/");
	     }
	     f->GetObject("l1PiXTRKTree/L1PiXTRKTree","");
#else
	     TChain *chain = new TChain("l1PiXTRKTree/L1PiXTRKTree","");

//	     chain->Add("/xrootd/store/user/jhong/SingleMu_Pt2to200_Eta3p0_CMSSW_9_3_7_NoPU_D17test/crab_Muon0802/180802_123200/*/*.root/L1PiXTRKTree");
	     chain->Add("root://cms-xrdr.private.lo:2094//xrd/store/user/jhong/SingleMu_Pt2to200_Eta3p0_CMSSW_9_3_7_NoPU_D17test/crab_Muon0802/180802_123200/0000/*.root/l1PiXTRKTree/L1PiXTRKTree");
//	     chain->Add("root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/jhong/SingleMu_Pt2to200_Eta3p0_CMSSW_9_3_7_NoPU_D17test/crab_Muon0802/180802_123200/0000/*.root/l1PiXTRKTree/L1PiXTRKTree");
//	     chain->Add("root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/jhong/SingleMu_Pt2to200_Eta3p0_CMSSW_9_3_7_NoPU_D17test/crab_Muon0802/180802_123200/0001/*.root/l1PiXTRKTree/L1PiXTRKTree");
//	     chain->Add("root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/jhong/SingleMu_Pt2to200_Eta3p0_CMSSW_9_3_7_NoPU_D17test/crab_Muon0802/180802_123200/0002/*.root/l1PiXTRKTree/L1PiXTRKTree");
//	     chain->Add("root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/jhong/SingleMu_Pt2to200_Eta3p0_CMSSW_9_3_7_NoPU_D17test/crab_Muon0802/180802_123200/0003/*.root/l1PiXTRKTree/L1PiXTRKTree");
//	     chain->Add("root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/jhong/SingleMu_Pt2to200_Eta3p0_CMSSW_9_3_7_NoPU_D17test/crab_Muon0802/180802_123200/0004/*.root/l1PiXTRKTree/L1PiXTRKTree");
       
	     tree = chain;
#endif
     }	
   Init(tree);
/*   
     if (tree == 0) {
#ifdef SINGLE_TREE
	     TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/xrootd/store/user/jhong/SingleMu_FlatPt-2to100/crab_SingleMu200PU/190122_053556/0000/");
	     if(!f || !f->IsOpen()){
			 f = new TFile("/xrootd/store/user/jhong/SingleMu_FlatPt-2to100/crab_SingleMu200PU/190122_053556/0000/");
	     }
	     f->GetObject("l1PiXTRKTree/L1PiXTRKTree","");
#else
	     TChain *chain = new TChain("l1PiXTRKTree/L1PiXTRKTree","");

//	     chain->Add("root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/jhong/SingleMu_FlatPt-2to100/crab_SingleMu200PU/190122_053556/0000/.root/l1PiXTRKTree/L1PiXTRKTree");
	     
	     tree = chain;
#endif
     }	
   Init(tree);
*/   
}

plot::~plot()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t plot::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t plot::LoadTree(Long64_t entry)
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

void plot::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   genPartE = 0;
   genPartPt = 0;
   genPartEta = 0;
   genPartPhi = 0;
   genPartCharge = 0;
   genPartId = 0;
   propgenElPartE = 0;
   propgenElPartPt = 0;
   propgenElPartEta = 0;
   propgenElPartPhi = 0;
   propgenElPartCharge = 0;
   propgenElPartx = 0;
   propgenElParty = 0;
   propgenElPartz = 0;
   simTrkPt = 0;
   simTrkEta = 0;
   simTrkPhi = 0;
   simTrkId = 0;
   simTrkType = 0;
   simTrkVx = 0;
   simTrkVy = 0;
   simTrkVz = 0;
   simVx = 0;
   simVy = 0;
   simVz = 0;
   Brempos_radius = 0;
   Brem_eLoss = 0;
   Brem_ptLoss = 0;
   Brempos_x = 0;
   Brempos_y = 0;
   Brempos_z = 0;
   propgenPoPartE = 0;
   propgenPoPartPt = 0;
   propgenPoPartEta = 0;
   propgenPoPartPhi = 0;
   propgenPoPartCharge = 0;
   propgenPoPartx = 0;
   propgenPoParty = 0;
   propgenPoPartz = 0;
   bRecHitLayer = 0;
   bRecHitLadder = 0;
   bRecHitModule = 0;
   fRecHitDisk = 0;
   fRecHitBlade = 0;
   fRecHitSide = 0;
   fRecHitPanel = 0;
   fRecHitModule = 0;
   fRecHitGx = 0;
   fRecHitGy = 0;
   fRecHitGz = 0;
   fRhSize = 0;
   fRhSizeX = 0;
   fRhSizeY = 0;
   bRecHitGx = 0;
   bRecHitGy = 0;
   bRecHitGz = 0;
   bRhSize = 0;
   bRhSizeX = 0;
   bRhSizeY = 0;
   bfastsimHitLayer = 0;
   bfastsimHitGx = 0;
   bfastsimHitGy = 0;
   bfastsimHitGz = 0;
   ffastsimHitLayer = 0;
   ffastsimHitGx = 0;
   ffastsimHitGy = 0;
   ffastsimHitGz = 0;
   egCrysE = 0;
   egCrysEt = 0;
   egCrysEta = 0;
   egCrysPhi = 0;
   egCrysGx = 0;
   egCrysGy = 0;
   egCrysGz = 0;
   egCrysClusterE = 0;
   egCrysClusterEt = 0;
   egCrysClusterEta = 0;
   egCrysClusterPhi = 0;
   egCrysClusterGx = 0;
   egCrysClusterGy = 0;
   egCrysClusterGz = 0;
   egCrysClusterPGx = 0;
   egCrysClusterPGy = 0;
   egCrysClusterPGz = 0;
   isTrackMatched = 0;
   isoConeNTrack = 0;
   isoConePtTrack = 0;
   trackHighestPt = 0;
   trackHighestPtEta = 0;
   trackHighestPtPhi = 0;
   trackHighestPtChi2 = 0;
   trackHighestPtCutChi2 = 0;
   trackHighestPtCutChi2Eta = 0;
   trackHighestPtCutChi2Phi = 0;
   trackHighestPtCutChi2Chi2 = 0;
   trackmatchingdR = 0;
   hgcal_isTrackMatched = 0;
   hgcal_isoConeNTrack = 0;
   hgcal_isoConePtTrack = 0;
   hgcal_trackHighestPt = 0;
   hgcal_trackHighestPtEta = 0;
   hgcal_trackHighestPtPhi = 0;
   hgcal_trackHighestPtChi2 = 0;
   hgcal_trackHighestPtCutChi2 = 0;
   hgcal_trackHighestPtCutChi2Eta = 0;
   hgcal_trackHighestPtCutChi2Phi = 0;
   hgcal_trackHighestPtCutChi2Chi2 = 0;
   hgcal_trackmatchingdR = 0;
   cl3d_pt = 0;
   cl3d_energy = 0;
   cl3d_eta = 0;
   cl3d_phi = 0;
   cl3d_nclu = 0;
   cl3d_x = 0;
   cl3d_y = 0;
   cl3d_z = 0;
   cl3d_hovere = 0;
   cl3d_showerlength = 0;
   cl3d_coreshowerlength = 0;
   cl3d_firstlayer = 0;
   cl3d_maxlayer = 0;
   cl3d_seetot = 0;
   cl3d_seemax = 0;
   cl3d_spptot = 0;
   cl3d_sppmax = 0;
   cl3d_szz = 0;
   cl3d_srrtot = 0;
   cl3d_srrmax = 0;
   cl3d_srrmean = 0;
   cl3d_emaxe = 0;
   egEt = 0;
   egEta = 0;
   egPhi = 0;
   egGx = 0;
   egGy = 0;
   egGz = 0;
   egIEt = 0;
   egIEta = 0;
   egIPhi = 0;
   egIso = 0;
   egBx = 0;
   egTowerIPhi = 0;
   egTowerIEta = 0;
   egRawEt = 0;
   egIsoEt = 0;
   egFootprintEt = 0;
   egNTT = 0;
   egShape = 0;
   egTowerHoE = 0;
   me0SegDetId = 0;
   me0SegPosX = 0;
   me0SegPosY = 0;
   me0SegPosZ = 0;
   me0SegDirX = 0;
   me0SegDirY = 0;
   me0SegDirZ = 0;
   me0SegNumRecHit = 0;
   me0SegDeltaPhi = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("nVtx", &nVtx, &b_nVtx);
   fChain->SetBranchAddress("nMeanPU", &nMeanPU, &b_nMeanPU);
   fChain->SetBranchAddress("genPartN", &genPartN, &b_genPartN);
   fChain->SetBranchAddress("genPartE", &genPartE, &b_genPartE);
   fChain->SetBranchAddress("genPartPt", &genPartPt, &b_genPartPt);
   fChain->SetBranchAddress("genPartEta", &genPartEta, &b_genPartEta);
   fChain->SetBranchAddress("genPartPhi", &genPartPhi, &b_genPartPhi);
   fChain->SetBranchAddress("genPartCharge", &genPartCharge, &b_genPartCharge);
   fChain->SetBranchAddress("genPartId", &genPartId, &b_genPartId);
   fChain->SetBranchAddress("propgenElPartE", &propgenElPartE, &b_propgenElPartE);
   fChain->SetBranchAddress("propgenElPartPt", &propgenElPartPt, &b_propgenElPartPt);
   fChain->SetBranchAddress("propgenElPartEta", &propgenElPartEta, &b_propgenElPartEta);
   fChain->SetBranchAddress("propgenElPartPhi", &propgenElPartPhi, &b_propgenElPartPhi);
   fChain->SetBranchAddress("propgenElPartCharge", &propgenElPartCharge, &b_propgenElPartCharge);
   fChain->SetBranchAddress("propgenElPartx", &propgenElPartx, &b_propgenElPartx);
   fChain->SetBranchAddress("propgenElParty", &propgenElParty, &b_propgenElParty);
   fChain->SetBranchAddress("propgenElPartz", &propgenElPartz, &b_propgenElPartz);
   fChain->SetBranchAddress("simTrkN", &simTrkN, &b_simTrkN);
   fChain->SetBranchAddress("simTrkPt", &simTrkPt, &b_simTrkPt);
   fChain->SetBranchAddress("simTrkEta", &simTrkEta, &b_simTrkEta);
   fChain->SetBranchAddress("simTrkPhi", &simTrkPhi, &b_simTrkPhi);
   fChain->SetBranchAddress("simTrkId", &simTrkId, &b_simTrkId);
   fChain->SetBranchAddress("simTrkType", &simTrkType, &b_simTrkType);
   fChain->SetBranchAddress("simTrkVx", &simTrkVx, &b_simTrkVx);
   fChain->SetBranchAddress("simTrkVy", &simTrkVy, &b_simTrkVy);
   fChain->SetBranchAddress("simTrkVz", &simTrkVz, &b_simTrkVz);
   fChain->SetBranchAddress("simVx", &simVx, &b_simVx);
   fChain->SetBranchAddress("simVy", &simVy, &b_simVy);
   fChain->SetBranchAddress("simVz", &simVz, &b_simVz);
   fChain->SetBranchAddress("lastSimtkpt", &lastSimtkpt, &b_lastSimtkpt);
   fChain->SetBranchAddress("initialSimtkpt", &initialSimtkpt, &b_initialSimtkpt);
   fChain->SetBranchAddress("bremflag", &bremflag, &b_bremflag);
   fChain->SetBranchAddress("Brempos_radius", &Brempos_radius, &b_Brempos_radius);
   fChain->SetBranchAddress("Brem_eLoss", &Brem_eLoss, &b_Brem_eLoss);
   fChain->SetBranchAddress("Brem_ptLoss", &Brem_ptLoss, &b_Brem_ptLoss);
   fChain->SetBranchAddress("Brempos_x", &Brempos_x, &b_Brempos_x);
   fChain->SetBranchAddress("Brempos_y", &Brempos_y, &b_Brempos_y);
   fChain->SetBranchAddress("Brempos_z", &Brempos_z, &b_Brempos_z);
   fChain->SetBranchAddress("propgenPoPartE", &propgenPoPartE, &b_propgenPoPartE);
   fChain->SetBranchAddress("propgenPoPartPt", &propgenPoPartPt, &b_propgenPoPartPt);
   fChain->SetBranchAddress("propgenPoPartEta", &propgenPoPartEta, &b_propgenPoPartEta);
   fChain->SetBranchAddress("propgenPoPartPhi", &propgenPoPartPhi, &b_propgenPoPartPhi);
   fChain->SetBranchAddress("propgenPoPartCharge", &propgenPoPartCharge, &b_propgenPoPartCharge);
   fChain->SetBranchAddress("propgenPoPartx", &propgenPoPartx, &b_propgenPoPartx);
   fChain->SetBranchAddress("propgenPoParty", &propgenPoParty, &b_propgenPoParty);
   fChain->SetBranchAddress("propgenPoPartz", &propgenPoPartz, &b_propgenPoPartz);
   fChain->SetBranchAddress("bRecHitLayer", &bRecHitLayer, &b_bRecHitLayer);
   fChain->SetBranchAddress("bRecHitLadder", &bRecHitLadder, &b_bRecHitLadder);
   fChain->SetBranchAddress("bRecHitModule", &bRecHitModule, &b_bRecHitModule);
   fChain->SetBranchAddress("fRecHitDisk", &fRecHitDisk, &b_fRecHitDisk);
   fChain->SetBranchAddress("fRecHitBlade", &fRecHitBlade, &b_fRecHitBlade);
   fChain->SetBranchAddress("fRecHitSide", &fRecHitSide, &b_fRecHitSide);
   fChain->SetBranchAddress("fRecHitPanel", &fRecHitPanel, &b_fRecHitPanel);
   fChain->SetBranchAddress("fRecHitModule", &fRecHitModule, &b_fRecHitModule);
   fChain->SetBranchAddress("bRecHitN", &bRecHitN, &b_bRecHitN);
   fChain->SetBranchAddress("fRecHitN", &fRecHitN, &b_fRecHitN);
   fChain->SetBranchAddress("fRecHitGx", &fRecHitGx, &b_fRecHitGx);
   fChain->SetBranchAddress("fRecHitGy", &fRecHitGy, &b_fRecHitGy);
   fChain->SetBranchAddress("fRecHitGz", &fRecHitGz, &b_fRecHitGz);
   fChain->SetBranchAddress("fRhSize", &fRhSize, &b_fRhSize);
   fChain->SetBranchAddress("fRhSizeX", &fRhSizeX, &b_fRhSizeX);
   fChain->SetBranchAddress("fRhSizeY", &fRhSizeY, &b_fRhSizeY);
   fChain->SetBranchAddress("bRecHitGx", &bRecHitGx, &b_bRecHitGx);
   fChain->SetBranchAddress("bRecHitGy", &bRecHitGy, &b_bRecHitGy);
   fChain->SetBranchAddress("bRecHitGz", &bRecHitGz, &b_bRecHitGz);
   fChain->SetBranchAddress("bRhSize", &bRhSize, &b_bRhSize);
   fChain->SetBranchAddress("bRhSizeX", &bRhSizeX, &b_bRhSizeX);
   fChain->SetBranchAddress("bRhSizeY", &bRhSizeY, &b_bRhSizeY);
   fChain->SetBranchAddress("bfastsimHitN", &bfastsimHitN, &b_bfastsimHitN);
   fChain->SetBranchAddress("ffastsimHitN", &ffastsimHitN, &b_ffastsimHitN);
   fChain->SetBranchAddress("bfastsimHitLayer", &bfastsimHitLayer, &b_bfastsimHitLayer);
   fChain->SetBranchAddress("bfastsimHitGx", &bfastsimHitGx, &b_bfastsimHitGx);
   fChain->SetBranchAddress("bfastsimHitGy", &bfastsimHitGy, &b_bfastsimHitGy);
   fChain->SetBranchAddress("bfastsimHitGz", &bfastsimHitGz, &b_bfastsimHitGz);
   fChain->SetBranchAddress("ffastsimHitLayer", &ffastsimHitLayer, &b_ffastsimHitLayer);
   fChain->SetBranchAddress("ffastsimHitGx", &ffastsimHitGx, &b_ffastsimHitGx);
   fChain->SetBranchAddress("ffastsimHitGy", &ffastsimHitGy, &b_ffastsimHitGy);
   fChain->SetBranchAddress("ffastsimHitGz", &ffastsimHitGz, &b_ffastsimHitGz);
   fChain->SetBranchAddress("egCrysN", &egCrysN, &b_egCrysN);
   fChain->SetBranchAddress("egCrysE", &egCrysE, &b_egCrysE);
   fChain->SetBranchAddress("egCrysEt", &egCrysEt, &b_egCrysEt);
   fChain->SetBranchAddress("egCrysEta", &egCrysEta, &b_egCrysEta);
   fChain->SetBranchAddress("egCrysPhi", &egCrysPhi, &b_egCrysPhi);
   fChain->SetBranchAddress("egCrysGx", &egCrysGx, &b_egCrysGx);
   fChain->SetBranchAddress("egCrysGy", &egCrysGy, &b_egCrysGy);
   fChain->SetBranchAddress("egCrysGz", &egCrysGz, &b_egCrysGz);
   fChain->SetBranchAddress("egCrysClusterN", &egCrysClusterN, &b_egCrysClusterN);
   fChain->SetBranchAddress("egCrysClusterE", &egCrysClusterE, &b_egCrysClusterE);
   fChain->SetBranchAddress("egCrysClusterEt", &egCrysClusterEt, &b_egCrysClusterEt);
   fChain->SetBranchAddress("egCrysClusterEta", &egCrysClusterEta, &b_egCrysClusterEta);
   fChain->SetBranchAddress("egCrysClusterPhi", &egCrysClusterPhi, &b_egCrysClusterPhi);
   fChain->SetBranchAddress("egCrysClusterGx", &egCrysClusterGx, &b_egCrysClusterGx);
   fChain->SetBranchAddress("egCrysClusterGy", &egCrysClusterGy, &b_egCrysClusterGy);
   fChain->SetBranchAddress("egCrysClusterGz", &egCrysClusterGz, &b_egCrysClusterGz);
   fChain->SetBranchAddress("egCrysClusterPGx", &egCrysClusterPGx, &b_egCrysClusterPGx);
   fChain->SetBranchAddress("egCrysClusterPGy", &egCrysClusterPGy, &b_egCrysClusterPGy);
   fChain->SetBranchAddress("egCrysClusterPGz", &egCrysClusterPGz, &b_egCrysClusterPGz);
   fChain->SetBranchAddress("isTrackMatched", &isTrackMatched, &b_isTrackMatched);
   fChain->SetBranchAddress("isoConeNTrack", &isoConeNTrack, &b_isoConeNTrack);
   fChain->SetBranchAddress("isoConePtTrack", &isoConePtTrack, &b_isoConePtTrack);
   fChain->SetBranchAddress("trackHighestPt", &trackHighestPt, &b_trackHighestPt);
   fChain->SetBranchAddress("trackHighestPtEta", &trackHighestPtEta, &b_trackHighestPtEta);
   fChain->SetBranchAddress("trackHighestPtPhi", &trackHighestPtPhi, &b_trackHighestPtPhi);
   fChain->SetBranchAddress("trackHighestPtChi2", &trackHighestPtChi2, &b_trackHighestPtChi2);
   fChain->SetBranchAddress("trackHighestPtCutChi2", &trackHighestPtCutChi2, &b_trackHighestPtCutChi2);
   fChain->SetBranchAddress("trackHighestPtCutChi2Eta", &trackHighestPtCutChi2Eta, &b_trackHighestPtCutChi2Eta);
   fChain->SetBranchAddress("trackHighestPtCutChi2Phi", &trackHighestPtCutChi2Phi, &b_trackHighestPtCutChi2Phi);
   fChain->SetBranchAddress("trackHighestPtCutChi2Chi2", &trackHighestPtCutChi2Chi2, &b_trackHighestPtCutChi2Chi2);
   fChain->SetBranchAddress("trackmatchingdR", &trackmatchingdR, &b_trackmatchingdR);
   fChain->SetBranchAddress("hgcal_isTrackMatched", &hgcal_isTrackMatched, &b_hgcal_isTrackMatched);
   fChain->SetBranchAddress("hgcal_isoConeNTrack", &hgcal_isoConeNTrack, &b_hgcal_isoConeNTrack);
   fChain->SetBranchAddress("hgcal_isoConePtTrack", &hgcal_isoConePtTrack, &b_hgcal_isoConePtTrack);
   fChain->SetBranchAddress("hgcal_trackHighestPt", &hgcal_trackHighestPt, &b_hgcal_trackHighestPt);
   fChain->SetBranchAddress("hgcal_trackHighestPtEta", &hgcal_trackHighestPtEta, &b_hgcal_trackHighestPtEta);
   fChain->SetBranchAddress("hgcal_trackHighestPtPhi", &hgcal_trackHighestPtPhi, &b_hgcal_trackHighestPtPhi);
   fChain->SetBranchAddress("hgcal_trackHighestPtChi2", &hgcal_trackHighestPtChi2, &b_hgcal_trackHighestPtChi2);
   fChain->SetBranchAddress("hgcal_trackHighestPtCutChi2", &hgcal_trackHighestPtCutChi2, &b_hgcal_trackHighestPtCutChi2);
   fChain->SetBranchAddress("hgcal_trackHighestPtCutChi2Eta", &hgcal_trackHighestPtCutChi2Eta, &b_hgcal_trackHighestPtCutChi2Eta);
   fChain->SetBranchAddress("hgcal_trackHighestPtCutChi2Phi", &hgcal_trackHighestPtCutChi2Phi, &b_hgcal_trackHighestPtCutChi2Phi);
   fChain->SetBranchAddress("hgcal_trackHighestPtCutChi2Chi2", &hgcal_trackHighestPtCutChi2Chi2, &b_hgcal_trackHighestPtCutChi2Chi2);
   fChain->SetBranchAddress("hgcal_trackmatchingdR", &hgcal_trackmatchingdR, &b_hgcal_trackmatchingdR);
   fChain->SetBranchAddress("cl3d_n", &cl3d_n, &b_cl3d_n);
   fChain->SetBranchAddress("cl3d_pt", &cl3d_pt, &b_cl3d_pt);
   fChain->SetBranchAddress("cl3d_energy", &cl3d_energy, &b_cl3d_energy);
   fChain->SetBranchAddress("cl3d_eta", &cl3d_eta, &b_cl3d_eta);
   fChain->SetBranchAddress("cl3d_phi", &cl3d_phi, &b_cl3d_phi);
   fChain->SetBranchAddress("cl3d_nclu", &cl3d_nclu, &b_cl3d_nclu);
   fChain->SetBranchAddress("cl3d_x", &cl3d_x, &b_cl3d_x);
   fChain->SetBranchAddress("cl3d_y", &cl3d_y, &b_cl3d_y);
   fChain->SetBranchAddress("cl3d_z", &cl3d_z, &b_cl3d_z);
   fChain->SetBranchAddress("cl3d_hovere", &cl3d_hovere, &b_cl3d_hovere);
   fChain->SetBranchAddress("cl3d_showerlength", &cl3d_showerlength, &b_cl3d_showerlength);
   fChain->SetBranchAddress("cl3d_coreshowerlength", &cl3d_coreshowerlength, &b_cl3d_coreshowerlength);
   fChain->SetBranchAddress("cl3d_firstlayer", &cl3d_firstlayer, &b_cl3d_firstlayer);
   fChain->SetBranchAddress("cl3d_maxlayer", &cl3d_maxlayer, &b_cl3d_maxlayer);
   fChain->SetBranchAddress("cl3d_seetot", &cl3d_seetot, &b_cl3d_seetot);
   fChain->SetBranchAddress("cl3d_seemax", &cl3d_seemax, &b_cl3d_seemax);
   fChain->SetBranchAddress("cl3d_spptot", &cl3d_spptot, &b_cl3d_spptot);
   fChain->SetBranchAddress("cl3d_sppmax", &cl3d_sppmax, &b_cl3d_sppmax);
   fChain->SetBranchAddress("cl3d_szz", &cl3d_szz, &b_cl3d_szz);
   fChain->SetBranchAddress("cl3d_srrtot", &cl3d_srrtot, &b_cl3d_srrtot);
   fChain->SetBranchAddress("cl3d_srrmax", &cl3d_srrmax, &b_cl3d_srrmax);
   fChain->SetBranchAddress("cl3d_srrmean", &cl3d_srrmean, &b_cl3d_srrmean);
   fChain->SetBranchAddress("cl3d_emaxe", &cl3d_emaxe, &b_cl3d_emaxe);
   fChain->SetBranchAddress("egN", &egN, &b_egN);
   fChain->SetBranchAddress("egEt", &egEt, &b_egEt);
   fChain->SetBranchAddress("egEta", &egEta, &b_egEta);
   fChain->SetBranchAddress("egPhi", &egPhi, &b_egPhi);
   fChain->SetBranchAddress("egGx", &egGx, &b_egGx);
   fChain->SetBranchAddress("egGy", &egGy, &b_egGy);
   fChain->SetBranchAddress("egGz", &egGz, &b_egGz);
   fChain->SetBranchAddress("egIEt", &egIEt, &b_egIEt);
   fChain->SetBranchAddress("egIEta", &egIEta, &b_egIEta);
   fChain->SetBranchAddress("egIPhi", &egIPhi, &b_egIPhi);
   fChain->SetBranchAddress("egIso", &egIso, &b_egIso);
   fChain->SetBranchAddress("egBx", &egBx, &b_egBx);
   fChain->SetBranchAddress("egTowerIPhi", &egTowerIPhi, &b_egTowerIPhi);
   fChain->SetBranchAddress("egTowerIEta", &egTowerIEta, &b_egTowerIEta);
   fChain->SetBranchAddress("egRawEt", &egRawEt, &b_egRawEt);
   fChain->SetBranchAddress("egIsoEt", &egIsoEt, &b_egIsoEt);
   fChain->SetBranchAddress("egFootprintEt", &egFootprintEt, &b_egFootprintEt);
   fChain->SetBranchAddress("egNTT", &egNTT, &b_egNTT);
   fChain->SetBranchAddress("egShape", &egShape, &b_egShape);
   fChain->SetBranchAddress("egTowerHoE", &egTowerHoE, &b_egTowerHoE);
   fChain->SetBranchAddress("me0SegNum", &me0SegNum, &b_me0SegNum);
   fChain->SetBranchAddress("me0SegDetId", &me0SegDetId, &b_me0SegDetId);
   fChain->SetBranchAddress("me0SegPosX", &me0SegPosX, &b_me0SegPosX);
   fChain->SetBranchAddress("me0SegPosY", &me0SegPosY, &b_me0SegPosY);
   fChain->SetBranchAddress("me0SegPosZ", &me0SegPosZ, &b_me0SegPosZ);
   fChain->SetBranchAddress("me0SegDirX", &me0SegDirX, &b_me0SegDirX);
   fChain->SetBranchAddress("me0SegDirY", &me0SegDirY, &b_me0SegDirY);
   fChain->SetBranchAddress("me0SegDirZ", &me0SegDirZ, &b_me0SegDirZ);
   fChain->SetBranchAddress("me0SegNumRecHit", &me0SegNumRecHit, &b_me0SegNumRecHit);
   fChain->SetBranchAddress("me0SegDeltaPhi", &me0SegDeltaPhi, &b_me0SegDeltaPhi);
   Notify();
}

Bool_t plot::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void plot::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t plot::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef plot_cxx
