//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Sep 28 21:00:10 2018 by ROOT version 6.10/09
// from TTree L1PiXTRKTree/L1PiXTRKTree
// found on file: /xrootd/store/user/jhong/SingleMu_Pt2to200_Eta3p0_CMSSW_9_3_7_NoPU_D17test/crab_Muon0802/180802_123200/0000/SingleMuonNoPU_1.root
//////////////////////////////////////////////////////////

#ifndef MuonEff_h
#define MuonEff_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"
#include "vector"
#include "vector"

#include <iostream>
#include <fstream>
#include <string>

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>

#include "withEM_SingleCrys_900pre6_v2.h"
#include "withoutEM_SingleCrys_900pre6.h"
#include "RegionOfInterest.h"

using namespace std;

class MuonEff {
private:

	map<TString, TH1*> maphist;

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

   int Ele, Pos;
   int skip;
   float me0N;
   int eta_region;

   int withoutEM_count_Ele, withEM_count_Ele;
   bool PixTrkPassed;
   int pass_count_wo4thPix, pass_count_wo3thPix, pass_count_wo2thPix, pass_count_wo1thPix;
   int woEM_pass_Ele_count_wo4thPix, woEM_pass_Ele_count_wo3thPix, woEM_pass_Ele_count_wo2thPix, woEM_pass_Ele_count_wo1thPix;
   int wEM_pass_Ele_count_wo4thPix, wEM_pass_Ele_count_wo3thPix, wEM_pass_Ele_count_wo2thPix, wEM_pass_Ele_count_wo1thPix;

   double all_cut_pass_eg;
   int all_cut_pass_Ele, withoutEM_pass_Ele, withEM_pass_Ele;
   int all_cut_pass_Pos;
   int fourth_layer_missing;
   int third_layer_missing;
   int second_layer_missing;
   int first_layer_missing;

   int bit1;
   int trigger_bit_width_;
   int pix_comb_;

   bool debug;

   double  L1_Dphi_cut1, L1_Dphi_cut2;
   double  L2_Dphi_cut1, L2_Dphi_cut2;
   double  L3_Dphi_cut1, L3_Dphi_cut2;
   double  L4_Dphi_cut1, L4_Dphi_cut2;
   double  D1_Dphi_cut1, D1_Dphi_cut2;
   double  D2_Dphi_cut1, D2_Dphi_cut2;
   double  D3_Dphi_cut1, D3_Dphi_cut2;

   double dPhi012;
   double dPhi013;
   double dPhi014;
   double dPhi023;
   double dPhi024;
   double dPhi034;

   double  L012_DPhi_cut1, L012_DPhi_cut2;
   
   double  L013_DPhi_cut1, L013_DPhi_cut2;
   
   double  L014_DPhi_cut1, L014_DPhi_cut2;
   
   double  L023_DPhi_cut1, L023_DPhi_cut2;
   
   double  L024_DPhi_cut1, L024_DPhi_cut2;
   
   double  L034_DPhi_cut1, L034_DPhi_cut2;
   
   double  L123_DPhi_cut1, L123_DPhi_cut2;
   double  L123_DEta_cut1, L123_DEta_cut2;
 
   double  L124_DPhi_cut1, L124_DPhi_cut2;
   double  L124_DEta_cut1, L124_DEta_cut2;
   
   double  L134_DPhi_cut1, L134_DPhi_cut2;
   double  L134_DEta_cut1, L134_DEta_cut2;
   
   double  L234_DPhi_cut1, L234_DPhi_cut2;
   double  L234_DEta_cut1, L234_DEta_cut2;

   double  L12_eta_upper, L13_eta_upper, L14_eta_upper, L23_eta_upper, L24_eta_upper, L34_eta_upper;
   double  L12_phi_upper, L13_phi_upper, L14_phi_upper, L23_phi_upper, L24_phi_upper, L34_phi_upper;
   double  L12_eta_bellow, L13_eta_bellow, L14_eta_bellow, L23_eta_bellow, L24_eta_bellow, L34_eta_bellow;
   double  L12_phi_bellow, L13_phi_bellow, L14_phi_bellow, L23_phi_bellow, L24_phi_bellow, L34_phi_bellow;
   double  L12_R_bellow, L13_R_bellow, L14_R_bellow, L23_R_bellow, L24_R_bellow, L34_R_bellow;

   double dPhi;
   double dEta;
   double dPhi_1, dPhi_2, dPhi_3;
   double dEta_1, dEta_2, dEta_3;
   TVector3 first_temp, second_temp;
   int _pass_Ele, _pass_Pos;

   int L012_pass_Ele, L012_pass_Pos; 
   int L013_pass_Ele, L013_pass_Pos;
   int L014_pass_Ele, L014_pass_Pos;
   int L023_pass_Ele, L023_pass_Pos;
   int L024_pass_Ele, L024_pass_Pos;
   int L034_pass_Ele, L034_pass_Pos;
   int L123_pass_Ele, L123_pass_Pos;
   int L124_pass_Ele, L124_pass_Pos;
   int L134_pass_Ele, L134_pass_Pos;
   int L234_pass_Ele, L234_pass_Pos;

   int L12_EM_Ele, L12_EM_Pos; 
   int L13_EM_Ele, L13_EM_Pos;
   int L14_EM_Ele, L14_EM_Pos;
   int L23_EM_Ele, L23_EM_Pos;
   int L24_EM_Ele, L24_EM_Pos;
   int L34_EM_Ele, L34_EM_Pos;

   std::vector<TVector3> first_disk_hits;
   std::vector<TVector3> second_disk_hits;
   std::vector<TVector3> third_disk_hits;
   std::vector<TVector3> fourth_disk_hits;

   double r; // r for radius of pixel tracker layer
   int layers[5];  // initialize as 0, layers contain # of hits on each pixel layer

   TVector3 emvector;
   float ME0Et;
   float ME0Eta;
   float ME0Phi;

   std::vector<int> first_disk_hits_Ele_or_Pos;
   std::vector<int> second_disk_hits_Ele_or_Pos;
   std::vector<int> third_disk_hits_Ele_or_Pos;
   std::vector<int> fourth_disk_hits_Ele_or_Pos;
   std::vector<int> hitted_layers;

   void MakeHistograms(TString hname, int nbins, float xmin, float xmax);
   TH1* GetHist(TString hname);
   void FillHist(TString histname, float value, float w, float xmin, float xmax, int nbins);

   void StorePixelHit( int region);
   double StandaloneDPhi( int first_hit, int second_hit, int third_hit, int which_first_hit, int which_second_hit, int which_third_hit );
   double StandaloneDEta( int first_hit, int second_hit, int third_hit, int which_first_hit, int which_second_hit, int which_third_hit );
   double EMmatchingDEta(TVector3& first_layer, TVector3& second_layer, TVector3& egvector);
   double EMmatchingDPhi(TVector3& first_layer, TVector3& second_layer, TVector3& egvector);
   int Signal_window_check( double upper, double value, double lower, int Ele_Pos);
   void FillCutFlow(TString cut, float weight);
   void SetROI(int region);
   void SetSignalBoundary(int region, double eg_dphi, double eg_deta, double sa_dphi, double sa_deta);
   void TriggeringWith_1st2ndPixel(int nthFirstHit, int nthSecondHit);
   void TriggeringWith_1st3rdPixel(int nthFirstHit, int nthSecondHit);
   void TriggeringWith_2nd3rdPixel(int nthFirstHit, int nthSecondHit);
   void TriggeringWithout_4thPixel(int nthFirstHit, int nthSecondHit, int nthThirdHit);
   void TriggeringWithout_3rdPixel(int nthFirstHit, int nthSecondHit, int nthThirdHit);
   void TriggeringWithout_2ndPixel(int nthFirstHit, int nthSecondHit, int nthThirdHit);
   void TriggeringWithout_1stPixel(int nthFirstHit, int nthSecondHit, int nthThirdHit);

   void TriggeringWith_1st2ndPixel_v2(int nthFirstHit, int nthSecondHit);
   void TriggeringWith_1st3rdPixel_v2(int nthFirstHit, int nthSecondHit);
   void TriggeringWith_1st4thPixel_v2(int nthFirstHit, int nthSecondHit);
   void TriggeringWith_2nd3rdPixel_v2(int nthFirstHit, int nthSecondHit);
   void TriggeringWith_2nd4thPixel_v2(int nthFirstHit, int nthSecondHit);
   void TriggeringWith_3rd4thPixel_v2(int nthFirstHit, int nthSecondHit);

   inline float deltaPhi(float phi1, float phi2) { 
     float result = phi1 - phi2;
     while (result > float(M_PI)) result -= float(2*M_PI);
     while (result <= -float(M_PI)) result += float(2*M_PI);
     return result;
   }


   MuonEff(TTree *tree=0);
   virtual ~MuonEff();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   TFile *outfile;
   TTree* pixtrk_tree;

   int count_Entry; 
   int pass_egobjects_check;
   int ntnME02; 
   int event_denominator; 
   int event_nominator; 

   int nPix123_segments;
   int nPix124_segments;
   int nPix134_segments;
   int nPix234_segments;

   vector<float> ntME0Et; 
   vector<float> ntME0Eta; 
   vector<float> ntME0Phi; 

   float matchedME0Et; 
   float matchedME0Eta; 
   float matchedME0Phi; 
   int   fired;

   vector<int> PiXTRKbit;
   vector<int> pix_comb;
   vector<int> trigger_bit_width;

   vector<bool> ntCl_match; 
   vector<bool> isTrack_match; 
   vector<float> chi2; 
   vector<float> track_dr;
   vector<bool> withoutEM_match; 
   vector<bool> withEM_match; 

   vector<int> ntfirstPix;  
   vector<int> ntsecondPix; 
   vector<int> ntthirdPix;  
   vector<int> ntfourthPix; 

   float nt_lastSimtkpt;
   float nt_initialSimtkpt;

   float nt_genPhi;
   float nt_genEta;
   float nt_genPt;

};

#endif

#ifdef MuonEff_cxx
MuonEff::MuonEff(TTree *tree) : fChain(0) 
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
   Init(tree);
*/

     if (tree == 0) {
#ifdef SINGLE_TREE
	    // TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/xrootd/store/user/jhong/SingleMu_Pt2to200_Eta3p0_CMSSW_9_3_7_NoPU_D17test/crab_Muon0802/180802_123200/*/");
	     TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/xrootd/store/user/jhong/SingleMu_FlatPt-2to100/crab_SingleMu200PU/190122_053556/0000/");
	     if(!f || !f->IsOpen()){
		    // f = new TFile("/xrootd/store/user/jhong/SingleMu_Pt2to200_Eta3p0_CMSSW_9_3_7_NoPU_D17test/crab_Muon0802/180802_123200/*/");
		     f = new TFile("xrootd/store/user/jhong/SingleMu_FlatPt-2to100/crab_SingleMu200PU/190122_053556/0000/");
	     }
	     f->GetObject("l1PiXTRKTree/L1PiXTRKTree","");
#else
	    // TChain *chain = new TChain("l1PiXTRKTree/L1PiXTRKTree","");
	     TChain *chain = new TChain("l1PiXTRKTree/L1PiXTRKTree","");

	     chain->Add("root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/jhong/SingleMu_FlatPt-2to100/crab_SingleMu200PU/190122_053556/0000/*.root/l1PiXTRKTree/L1PiXTRKTree");
//	     chain->Add("/xrootd/store/user/jhong/SingleMu_Pt2to200_Eta3p0_CMSSW_9_3_7_NoPU_D17test/crab_Muon0802/180802_123200/*/*.root/L1PiXTRKTree");
	   //  chain->Add("root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/jhong/SingleMu_Pt2to200_Eta3p0_CMSSW_9_3_7_NoPU_D17test/crab_Muon0802/180802_123200/0000/*.root/l1PiXTRKTree/L1PiXTRKTree");
	    // chain->Add("root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/jhong/SingleMu_Pt2to200_Eta3p0_CMSSW_9_3_7_NoPU_D17test/crab_Muon0802/180802_123200/0001/*.root/l1PiXTRKTree/L1PiXTRKTree");
	    // chain->Add("root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/jhong/SingleMu_Pt2to200_Eta3p0_CMSSW_9_3_7_NoPU_D17test/crab_Muon0802/180802_123200/0002/*.root/l1PiXTRKTree/L1PiXTRKTree");
	    // chain->Add("root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/jhong/SingleMu_Pt2to200_Eta3p0_CMSSW_9_3_7_NoPU_D17test/crab_Muon0802/180802_123200/0003/*.root/l1PiXTRKTree/L1PiXTRKTree");
	    // chain->Add("root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/jhong/SingleMu_Pt2to200_Eta3p0_CMSSW_9_3_7_NoPU_D17test/crab_Muon0802/180802_123200/0004/*.root/l1PiXTRKTree/L1PiXTRKTree");
	     
	     tree = chain;
#endif
     }	
   Init(tree);

   Ele = 1, Pos = 2;
   skip = 1;

//   outfile = new TFile("../output_tmp/Tree_output","recreate");
   outfile = new TFile("Eff_200PU.root","recreate");
   pixtrk_tree = new TTree("t","t");

   count_Entry = 1;
   pixtrk_tree->Branch("totalEvent", &count_Entry, "count_Entry/I");   pixtrk_tree->Branch("totalme0N", &me0N, "me0N/F");
   pixtrk_tree->Branch("ntnME02", &ntnME02, "ntnME02/I");

   pixtrk_tree->Branch("ntME0Et",&ntME0Et);
   pixtrk_tree->Branch("ntME0Eta",&ntME0Eta);
   pixtrk_tree->Branch("ntME0Phi",&ntME0Phi);

   pixtrk_tree->Branch("PiXTRKbit",&PiXTRKbit);
   pixtrk_tree->Branch("trigger_bit_width",&trigger_bit_width);
   pixtrk_tree->Branch("pix_comb",&pix_comb);

   pixtrk_tree->Branch("nPix123_segments",&nPix123_segments,"nPix123_segments/I");
   pixtrk_tree->Branch("nPix124_segments",&nPix124_segments,"nPix124_segments/I");
   pixtrk_tree->Branch("nPix134_segments",&nPix134_segments,"nPix134_segments/I");
   pixtrk_tree->Branch("nPix234_segments",&nPix234_segments,"nPix234_segments/I");

   pixtrk_tree->Branch("ntCl_match",&ntCl_match);
   pixtrk_tree->Branch("isTrack_match",&isTrack_match);
   pixtrk_tree->Branch("chi2",&chi2);
   pixtrk_tree->Branch("track_dr",&track_dr);
   pixtrk_tree->Branch("withoutEM_match",&withoutEM_match);
   pixtrk_tree->Branch("withEM_match",&withEM_match);

   pixtrk_tree->Branch("ntfirstPix",&ntfirstPix);
   pixtrk_tree->Branch("ntsecondPix",&ntsecondPix);
   pixtrk_tree->Branch("ntthirdPix",&ntthirdPix);
   pixtrk_tree->Branch("ntfourthPix",&ntfourthPix);

   pixtrk_tree->Branch("nt_genPhi",&nt_genPhi,"nt_genPhi/F");
   pixtrk_tree->Branch("nt_genEta",&nt_genEta,"nt_genEta/F");
   pixtrk_tree->Branch("nt_genPt",&nt_genPt,"nt_genPt/F");

   pixtrk_tree->Branch("matchedME0Et",&matchedME0Et,"matchedME0Et/F");
   pixtrk_tree->Branch("matchedME0Eta",&matchedME0Eta,"matchedME0Eta/F");
   pixtrk_tree->Branch("matchedME0Phi",&matchedME0Phi,"matchedME0Phi/F");

   pixtrk_tree->Branch("nt_lastSimtkpt",&nt_lastSimtkpt,"nt_lastSimtkpt/F");
   pixtrk_tree->Branch("nt_initialSimtkpt",&nt_initialSimtkpt,"nt_initialSimtkpt/F");

   pixtrk_tree->Branch("fired",&fired,"fired/I");


}

MuonEff::~MuonEff()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t MuonEff::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t MuonEff::LoadTree(Long64_t entry)
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

void MuonEff::Init(TTree *tree)
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

Bool_t MuonEff::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void MuonEff::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t MuonEff::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

void MuonEff::MakeHistograms(TString hname, int nbins, float xmin, float xmax){

 maphist[hname] =  new TH1F(hname.Data(),hname.Data(),nbins,xmin,xmax);

}

TH1* MuonEff::GetHist(TString hname){

 TH1* h = NULL;
 std::map<TString, TH1*>::iterator mapit = maphist.find(hname);
 if(mapit != maphist.end()) return mapit->second;

 return h;

}

void MuonEff::FillHist(TString histname, float value, float w, float xmin, float xmax, int nbins){

 if(GetHist(histname)) GetHist(histname)->Fill(value, w);
 else{
//     cout << "Making histogram..." << endl;
     MakeHistograms(histname, nbins, xmin, xmax);
     if(GetHist(histname)) GetHist(histname)->Fill(value, w);
 }
}

double MuonEff::StandaloneDPhi( int first_hit, int second_hit, int third_hit, int which_first_hit, int which_second_hit, int which_third_hit ){
  if( first_hit == 0 ){

     TVector3 temp;

     if( second_hit == 1 ) temp.SetXYZ( first_disk_hits[which_second_hit].X(), first_disk_hits[which_second_hit].Y(), first_disk_hits[which_second_hit].Z() );
     if( second_hit == 2 ) temp.SetXYZ( second_disk_hits[which_second_hit].X(), second_disk_hits[which_second_hit].Y(), second_disk_hits[which_second_hit].Z() );
     if( second_hit == 3 ) temp.SetXYZ( third_disk_hits[which_second_hit].X(), third_disk_hits[which_second_hit].Y(), third_disk_hits[which_second_hit].Z() );


     if( third_hit == 2 ) return deltaPhi( (second_disk_hits[which_third_hit] - temp).Phi(), temp.Phi());
     if( third_hit == 3 ) return deltaPhi( (third_disk_hits[which_third_hit] - temp).Phi(),  temp.Phi());
     if( third_hit == 4 ) return deltaPhi( (fourth_disk_hits[which_third_hit] - temp).Phi(), temp.Phi());

  }
  if( first_hit != 0 ){
    TVector3 temp_first_layer;
    TVector3 temp_second_layer;
    TVector3 temp_third_layer;

    if( first_hit == 1 ) temp_first_layer = first_disk_hits[which_first_hit];
    if( first_hit == 2 ) temp_first_layer = second_disk_hits[which_first_hit];

    if( second_hit == 2 ) temp_second_layer = second_disk_hits[which_second_hit];
    if( second_hit == 3 ) temp_second_layer = third_disk_hits[which_second_hit];


    if( third_hit == 3 ) temp_third_layer = third_disk_hits[which_third_hit];
    if( third_hit == 4 ) temp_third_layer = fourth_disk_hits[which_third_hit];

    return deltaPhi( (temp_third_layer - temp_second_layer).Phi(), (temp_second_layer - temp_first_layer).Phi());
  }
  return 0.;
}

double MuonEff::StandaloneDEta( int first_hit, int second_hit, int third_hit, int which_first_hit, int which_second_hit, int which_third_hit ){

    TVector3 temp_first_layer;
    TVector3 temp_second_layer;
    TVector3 temp_third_layer;

    if( first_hit == 1 ) temp_first_layer = first_disk_hits[which_first_hit];
    if( first_hit == 2 ) temp_first_layer = second_disk_hits[which_first_hit];

    if( second_hit == 2 ) temp_second_layer = second_disk_hits[which_second_hit];
    if( second_hit == 3 ) temp_second_layer = third_disk_hits[which_second_hit];


    if( third_hit == 3 ) temp_third_layer = third_disk_hits[which_third_hit];
    if( third_hit == 4 ) temp_third_layer = fourth_disk_hits[which_third_hit];

    return (temp_third_layer - temp_second_layer).PseudoRapidity() - (temp_second_layer - temp_first_layer).PseudoRapidity();
}

double MuonEff::EMmatchingDEta(TVector3& first_layer, TVector3& second_layer, TVector3& egvector){

    TVector3 pixelVector = second_layer - first_layer;
    TVector3 EM_pixelVector = egvector - second_layer;
    return EM_pixelVector.Eta() - pixelVector.Eta();

}

double MuonEff::EMmatchingDPhi(TVector3& first_layer, TVector3& second_layer, TVector3& egvector){


    TVector3 pixelVector = second_layer - first_layer;
    TVector3 EM_pixelVector = egvector - second_layer;

    return deltaPhi( EM_pixelVector.Phi(), pixelVector.Phi() );
}

int MuonEff::Signal_window_check( double upper, double value, double lower, int Ele_Pos){

  if( Ele_Pos == 1 ){ // 1 is Electron
    if( value <= upper && value >= lower){
      return true;
    }
    else
     return false;
  }
  if( Ele_Pos == 2 ){ // 2 is Positron
    if( value >= -upper && value <= -lower){
      return true;
    }
    else
     return false;
  }
 return 0;
}

void MuonEff::FillCutFlow(TString cut, float weight){


  if(GetHist("cutflow")) {
    GetHist("cutflow")->Fill(cut,weight);

  }
  else{
    MuonEff::MakeHistograms("cutflow", 6,0.,6.);

    GetHist("cutflow")->GetXaxis()->SetBinLabel(1,"NoCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(2,"MinEtCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(3,"EtaCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(4,"DRCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(5,"PtErrCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(6,"EvtCut");


  }
}

void MuonEff::StorePixelHit(int region){

        for(int a=0; a<fRecHitN; a++){
           int Dphi_Ele_pass = 0;
           int Dphi_Pos_pass = 0;
           double Dphi = 0.;
           int el_or_po = 0; // electron = 1, positron = 2, both = 3
           TVector3 current_hit;
           current_hit.SetXYZ( fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a) );
           Dphi = deltaPhi(current_hit.Phi(), ME0Phi);

           if( region == 2 ){

             if( fRecHitDisk->at(a) == 2 ){ // first disk
               if( Dphi < L1_Dphi_cut1 && Dphi > L1_Dphi_cut2){
                  Dphi_Ele_pass = 1; el_or_po = 1;
               }
               if( Dphi > -L1_Dphi_cut1 && Dphi < -L1_Dphi_cut2){
                  Dphi_Pos_pass = 1; el_or_po = el_or_po + 2;
               }
               if( Dphi_Ele_pass || Dphi_Pos_pass ){
                 layers[1]++;
                 first_disk_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a)));
                 first_disk_hits_Ele_or_Pos.push_back(el_or_po);
               }
             }

             if( fRecHitDisk->at(a) == 3 ){ // second disk
               if( Dphi < L2_Dphi_cut1 && Dphi > L2_Dphi_cut2){
                  Dphi_Ele_pass = 1; el_or_po = 1;
               }
               if( Dphi > -L2_Dphi_cut1 && Dphi < -L2_Dphi_cut2){
                  Dphi_Pos_pass = 1; el_or_po = el_or_po + 2;
               }
               if( Dphi_Ele_pass || Dphi_Pos_pass ){
                 layers[2]++;
                 second_disk_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a)));
                 second_disk_hits_Ele_or_Pos.push_back(el_or_po);
               }
             }
             if( fRecHitDisk->at(a) == 4 ){ // third disk
               if( Dphi < L3_Dphi_cut1 && Dphi > L3_Dphi_cut2){
                  Dphi_Ele_pass = 1; el_or_po = 1;
               }
               if( Dphi > -L3_Dphi_cut1 && Dphi < -L3_Dphi_cut2){
                  Dphi_Pos_pass = 1; el_or_po = el_or_po + 2;
               }
               if( Dphi_Ele_pass || Dphi_Pos_pass ){
                 layers[3]++;
                 third_disk_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a)));
                 third_disk_hits_Ele_or_Pos.push_back(el_or_po);
               }
             }

             if( fRecHitDisk->at(a) == 5 ){ // fourth disk
               if( Dphi < L4_Dphi_cut1 && Dphi > L4_Dphi_cut2){
                  Dphi_Ele_pass = 1; el_or_po = 1;
               }
               if( Dphi > -L4_Dphi_cut1 && Dphi < -L4_Dphi_cut2){
                  Dphi_Pos_pass = 1; el_or_po = el_or_po + 2;
               }
               if( Dphi_Ele_pass || Dphi_Pos_pass ){
                 layers[4]++;
                 fourth_disk_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a)));
                 fourth_disk_hits_Ele_or_Pos.push_back(el_or_po);
               }
             }

           }// region = 2

           if( region == 1 ){

             if( fRecHitDisk->at(a) == 1 ){ // first disk
               if( Dphi < L1_Dphi_cut1 && Dphi > L1_Dphi_cut2){
                  Dphi_Ele_pass = 1; el_or_po = 1;
               }
               if( Dphi > -L1_Dphi_cut1 && Dphi < -L1_Dphi_cut2){
                  Dphi_Pos_pass = 1; el_or_po = el_or_po + 2;
               }
               if( Dphi_Ele_pass || Dphi_Pos_pass ){
                 layers[1]++;
                 first_disk_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a)));
                 first_disk_hits_Ele_or_Pos.push_back(el_or_po);
               }
             }

             if( fRecHitDisk->at(a) == 2 ){ // second disk
               if( Dphi < L2_Dphi_cut1 && Dphi > L2_Dphi_cut2){
                  Dphi_Ele_pass = 1; el_or_po = 1;
               }
               if( Dphi > -L2_Dphi_cut1 && Dphi < -L2_Dphi_cut2){
                  Dphi_Pos_pass = 1; el_or_po = el_or_po + 2;
               }
               if( Dphi_Ele_pass || Dphi_Pos_pass ){
                 layers[2]++;
                 second_disk_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a)));
                 second_disk_hits_Ele_or_Pos.push_back(el_or_po);
               }
             }
             if( fRecHitDisk->at(a) == 3 ){ // third disk
               if( Dphi < L3_Dphi_cut1 && Dphi > L3_Dphi_cut2){
                  Dphi_Ele_pass = 1; el_or_po = 1;
               }
               if( Dphi > -L3_Dphi_cut1 && Dphi < -L3_Dphi_cut2){
                  Dphi_Pos_pass = 1; el_or_po = el_or_po + 2;
               }
               if( Dphi_Ele_pass || Dphi_Pos_pass ){
                 layers[3]++;
                 third_disk_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a)));
                 third_disk_hits_Ele_or_Pos.push_back(el_or_po);
               }
             }

             if( fRecHitDisk->at(a) == 4 ){ // fourth disk
               if( Dphi < L4_Dphi_cut1 && Dphi > L4_Dphi_cut2){
                  Dphi_Ele_pass = 1; el_or_po = 1;
               }
               if( Dphi > -L4_Dphi_cut1 && Dphi < -L4_Dphi_cut2){
                  Dphi_Pos_pass = 1; el_or_po = el_or_po + 2;
               }
               if( Dphi_Ele_pass || Dphi_Pos_pass ){
                 layers[4]++;
                 fourth_disk_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a)));
                 fourth_disk_hits_Ele_or_Pos.push_back(el_or_po);
               }
             }

           }// region = 1
        }

      ntfirstPix.push_back(layers[1]);
      ntsecondPix.push_back(layers[2]);
      ntthirdPix.push_back(layers[3]);
      ntfourthPix.push_back(layers[4]);
}

void MuonEff::SetROI(int region){


  float upper_width = 0.055;
  float lower_width = 0.055;

  L1_Dphi_cut1 = ROI_func(region, ME0Et);
  L1_Dphi_cut2 = ROI_func(region, ME0Et);

  L1_Dphi_cut1 = L1_Dphi_cut1 + upper_width;
  L1_Dphi_cut2 = L1_Dphi_cut2 - lower_width;

  L2_Dphi_cut1 = ROI_func(region, ME0Et);
  L2_Dphi_cut2 = ROI_func(region, ME0Et);

  L2_Dphi_cut1 = L2_Dphi_cut1 + upper_width;
  L2_Dphi_cut2 = L2_Dphi_cut2 - lower_width;

  L3_Dphi_cut1 = ROI_func(region, ME0Et);
  L3_Dphi_cut2 = ROI_func(region, ME0Et);

  L3_Dphi_cut1 = L3_Dphi_cut1 + upper_width;
  L3_Dphi_cut2 = L3_Dphi_cut2 - lower_width;

  L4_Dphi_cut1 = ROI_func(region, ME0Et);
  L4_Dphi_cut2 = ROI_func(region, ME0Et);

  L4_Dphi_cut1 = L4_Dphi_cut1 + upper_width;
  L4_Dphi_cut2 = L4_Dphi_cut2 - lower_width;

}

void MuonEff::SetSignalBoundary(int region, double eg_dphi, double eg_deta, double sa_dphi, double sa_deta){

      float EG_pixel_dphi_upper_width = eg_dphi;
      float EG_pixel_dphi_lower_width = eg_dphi;

      // pixel-ME0 
      L12_phi_upper =  SW_func2_dphi_v2(region, ME0Et);
      L12_phi_bellow = SW_func2_dphi_v2(region, ME0Et);
 
      L12_phi_upper = L12_phi_upper   + EG_pixel_dphi_upper_width;
      L12_phi_bellow = L12_phi_bellow - EG_pixel_dphi_lower_width;

      L13_phi_upper =  SW_func2_dphi_v2(region, ME0Et);
      L13_phi_bellow = SW_func2_dphi_v2(region, ME0Et);

      L13_phi_upper =  L13_phi_upper  + EG_pixel_dphi_upper_width;
      L13_phi_bellow = L13_phi_bellow - EG_pixel_dphi_lower_width;

      L14_phi_upper =  SW_func2_dphi_v2(region, ME0Et);
      L14_phi_bellow = SW_func2_dphi_v2(region, ME0Et);

      L14_phi_upper =  L14_phi_upper  + EG_pixel_dphi_upper_width;
      L14_phi_bellow = L14_phi_bellow - EG_pixel_dphi_lower_width;

      L23_phi_upper =  SW_func2_dphi_v2(region, ME0Et);
      L23_phi_bellow = SW_func2_dphi_v2(region, ME0Et);

      L23_phi_upper =  L23_phi_upper  + EG_pixel_dphi_upper_width;
      L23_phi_bellow = L23_phi_bellow - EG_pixel_dphi_lower_width;

      L24_phi_upper =  SW_func2_dphi_v2(region, ME0Et);
      L24_phi_bellow = SW_func2_dphi_v2(region, ME0Et);

      L24_phi_upper =  L24_phi_upper  + EG_pixel_dphi_upper_width;
      L24_phi_bellow = L24_phi_bellow - EG_pixel_dphi_lower_width;

      L34_phi_upper =  SW_func2_dphi_v2(region, ME0Et);
      L34_phi_bellow = SW_func2_dphi_v2(region, ME0Et);

      L34_phi_upper = L34_phi_upper   + EG_pixel_dphi_upper_width;
      L34_phi_bellow = L34_phi_bellow - EG_pixel_dphi_lower_width;

      float EG_pixel_deta_upper_width = eg_deta;
      float EG_pixel_deta_lower_width = eg_deta;

      L12_eta_upper =  SW_func2_deta_v2(ME0Et);
      L12_eta_bellow = SW_func2_deta_v2(ME0Et);

      L12_eta_upper =  L12_eta_upper  + EG_pixel_deta_upper_width;
      L12_eta_bellow = L12_eta_bellow - EG_pixel_deta_lower_width;

      L13_eta_upper =  SW_func2_deta_v2(ME0Et);
      L13_eta_bellow = SW_func2_deta_v2(ME0Et);

      L13_eta_upper =  L13_eta_upper  + EG_pixel_deta_upper_width;
      L13_eta_bellow = L13_eta_bellow - EG_pixel_deta_lower_width;

      L14_eta_upper =  SW_func2_deta_v2(ME0Et);
      L14_eta_bellow = SW_func2_deta_v2(ME0Et);

      L14_eta_upper =  L14_eta_upper  + EG_pixel_deta_upper_width;
      L14_eta_bellow = L14_eta_bellow - EG_pixel_deta_lower_width;

      L23_eta_upper =  SW_func2_deta_v2(ME0Et);
      L23_eta_bellow = SW_func2_deta_v2(ME0Et);

      L23_eta_upper =  L23_eta_upper  + EG_pixel_deta_upper_width;
      L23_eta_bellow = L23_eta_bellow - EG_pixel_deta_lower_width;

      L24_eta_upper =  SW_func2_deta_v2(ME0Et);
      L24_eta_bellow = SW_func2_deta_v2(ME0Et);

      L24_eta_upper =  L24_eta_upper  + EG_pixel_deta_upper_width;
      L24_eta_bellow = L24_eta_bellow - EG_pixel_deta_lower_width;

      L34_eta_upper =  SW_func2_deta_v2(ME0Et);
      L34_eta_bellow = SW_func2_deta_v2(ME0Et);

      L34_eta_upper =  L34_eta_upper  + EG_pixel_deta_upper_width;
      L34_eta_bellow = L34_eta_bellow - EG_pixel_deta_lower_width;

      // pixel-pixel 
      float pixel_pixel_dphi_upper_width = sa_dphi;
      float pixel_pixel_dphi_lower_width = sa_dphi;

      L012_DPhi_cut1 = SW_func1_dphi_v2(region, 0, ME0Et);
      L012_DPhi_cut2 = SW_func1_dphi_v2(region, 0, ME0Et);

      L012_DPhi_cut1 = L012_DPhi_cut1 + pixel_pixel_dphi_upper_width;
      L012_DPhi_cut2 = L012_DPhi_cut2 - pixel_pixel_dphi_lower_width;

      L013_DPhi_cut1 = SW_func1_dphi_v2(region, 1, ME0Et);
      L013_DPhi_cut2 = SW_func1_dphi_v2(region, 1, ME0Et);

      L013_DPhi_cut1 = L013_DPhi_cut1 + pixel_pixel_dphi_upper_width;
      L013_DPhi_cut2 = L013_DPhi_cut2 - pixel_pixel_dphi_lower_width;

      L014_DPhi_cut1 = SW_func1_dphi_v2(region, 2, ME0Et);
      L014_DPhi_cut2 = SW_func1_dphi_v2(region, 2, ME0Et);

      L014_DPhi_cut1 = L014_DPhi_cut1 + pixel_pixel_dphi_upper_width;
      L014_DPhi_cut2 = L014_DPhi_cut2 - pixel_pixel_dphi_lower_width;

      L023_DPhi_cut1 = SW_func1_dphi_v2(region, 3, ME0Et);
      L023_DPhi_cut2 = SW_func1_dphi_v2(region, 3, ME0Et);

      L023_DPhi_cut1 = L023_DPhi_cut1 + pixel_pixel_dphi_upper_width;
      L023_DPhi_cut2 = L023_DPhi_cut2 - pixel_pixel_dphi_lower_width;

      L024_DPhi_cut1 = SW_func1_dphi_v2(region, 4, ME0Et);
      L024_DPhi_cut2 = SW_func1_dphi_v2(region, 4, ME0Et);

      L024_DPhi_cut1 = L024_DPhi_cut1 + pixel_pixel_dphi_upper_width;
      L024_DPhi_cut2 = L024_DPhi_cut2 - pixel_pixel_dphi_lower_width;

      L034_DPhi_cut1 = SW_func1_dphi_v2(region, 5, ME0Et);
      L034_DPhi_cut2 = SW_func1_dphi_v2(region, 5, ME0Et);

      L034_DPhi_cut1 = L034_DPhi_cut1 + pixel_pixel_dphi_upper_width;
      L034_DPhi_cut2 = L034_DPhi_cut2 - pixel_pixel_dphi_lower_width;

      L123_DPhi_cut1 = SW_func1_dphi_v2(region, 6, ME0Et);
      L123_DPhi_cut2 = SW_func1_dphi_v2(region, 6, ME0Et);

      L123_DPhi_cut1 = L123_DPhi_cut1 + pixel_pixel_dphi_upper_width;
      L123_DPhi_cut2 = L123_DPhi_cut2 - pixel_pixel_dphi_lower_width;

      L124_DPhi_cut1 = SW_func1_dphi_v2(region, 7, ME0Et);
      L124_DPhi_cut2 = SW_func1_dphi_v2(region, 7 ,ME0Et);

      L124_DPhi_cut1 = L124_DPhi_cut1 + pixel_pixel_dphi_upper_width;
      L124_DPhi_cut2 = L124_DPhi_cut2 - pixel_pixel_dphi_lower_width;

      L134_DPhi_cut1 = SW_func1_dphi_v2(region, 8, ME0Et);
      L134_DPhi_cut2 = SW_func1_dphi_v2(region, 8, ME0Et);

      L134_DPhi_cut1 = L134_DPhi_cut1 + pixel_pixel_dphi_upper_width;
      L134_DPhi_cut2 = L134_DPhi_cut2 - pixel_pixel_dphi_lower_width;

      L234_DPhi_cut1 = SW_func1_dphi_v2(region, 9, ME0Et);
      L234_DPhi_cut2 = SW_func1_dphi_v2(region, 9, ME0Et);

      L234_DPhi_cut1 = L234_DPhi_cut1 + pixel_pixel_dphi_upper_width;
      L234_DPhi_cut2 = L234_DPhi_cut2 - pixel_pixel_dphi_lower_width;

      float pixel_pixel_deta_upper_width = sa_deta;
      float pixel_pixel_deta_lower_width = sa_deta;

      L123_DEta_cut1 = SW_func1_deta_v2(0, ME0Et);
      L123_DEta_cut2 = SW_func1_deta_v2(0, ME0Et);

      L123_DEta_cut1 = L123_DEta_cut1 + pixel_pixel_deta_upper_width;
      L123_DEta_cut2 = L123_DEta_cut2 - pixel_pixel_deta_lower_width;

      L124_DEta_cut1 = SW_func1_deta_v2(1, ME0Et);
      L124_DEta_cut2 = SW_func1_deta_v2(1 ,ME0Et);

      L124_DEta_cut1 = L124_DEta_cut1 + pixel_pixel_deta_upper_width;
      L124_DEta_cut2 = L124_DEta_cut2 - pixel_pixel_deta_lower_width;

      L134_DEta_cut1 = SW_func1_deta_v2(2, ME0Et);
      L134_DEta_cut2 = SW_func1_deta_v2(2, ME0Et);

      L134_DEta_cut1 = L134_DEta_cut1 + pixel_pixel_deta_upper_width;
      L134_DEta_cut2 = L134_DEta_cut2 - pixel_pixel_deta_lower_width;

      L234_DEta_cut1 = SW_func1_deta_v2(3, ME0Et);
      L234_DEta_cut2 = SW_func1_deta_v2(3, ME0Et);

      L234_DEta_cut1 = L234_DEta_cut1 + pixel_pixel_deta_upper_width;
      L234_DEta_cut2 = L234_DEta_cut2 - pixel_pixel_deta_lower_width;

}
void MuonEff::TriggeringWithout_4thPixel( int nthFirstHit, int nthSecondHit, int nthThirdHit){

dPhi012 = StandaloneDPhi( 0, 1, 2, 0, nthFirstHit, nthSecondHit);
dPhi013 = StandaloneDPhi( 0, 1, 3, 0, nthFirstHit, nthThirdHit );
dPhi023 = StandaloneDPhi( 0, 2, 3, 0, nthSecondHit, nthThirdHit);

dPhi_1 = EMmatchingDPhi(first_disk_hits[nthFirstHit], second_disk_hits[nthSecondHit], emvector);
dEta_1 = EMmatchingDEta(first_disk_hits[nthFirstHit], second_disk_hits[nthSecondHit], emvector);

dPhi_2 = EMmatchingDPhi(first_disk_hits[nthFirstHit], third_disk_hits[nthThirdHit], emvector);
dEta_2 = EMmatchingDEta(first_disk_hits[nthFirstHit], third_disk_hits[nthThirdHit], emvector);

dPhi_3 = EMmatchingDPhi(second_disk_hits[nthSecondHit], third_disk_hits[nthThirdHit], emvector);
dEta_3 = EMmatchingDEta(second_disk_hits[nthSecondHit], third_disk_hits[nthThirdHit], emvector);

L012_pass_Ele = Signal_window_check( L012_DPhi_cut1, dPhi012, L012_DPhi_cut2, Ele );
L012_pass_Pos = Signal_window_check( L012_DPhi_cut1, dPhi012, L012_DPhi_cut2, Pos );
L013_pass_Ele = Signal_window_check( L013_DPhi_cut1, dPhi013, L013_DPhi_cut2, Ele );
L013_pass_Pos = Signal_window_check( L013_DPhi_cut1, dPhi013, L013_DPhi_cut2, Pos );
L023_pass_Ele = Signal_window_check( L023_DPhi_cut1, dPhi023, L023_DPhi_cut2, Ele );
L023_pass_Pos = Signal_window_check( L023_DPhi_cut1, dPhi023, L023_DPhi_cut2, Pos );

if( Signal_window_check(L12_phi_upper, dPhi_1, L12_phi_bellow, Ele) && Signal_window_check(L12_eta_upper, dEta_1, L12_eta_bellow, Ele)  )
   L12_EM_Ele = 1;
if( Signal_window_check(L12_phi_upper, dPhi_1, L12_phi_bellow, Pos) && Signal_window_check(L12_eta_upper, dEta_1, L12_eta_bellow, Ele)  )
   L12_EM_Pos = 1;

if( Signal_window_check(L13_phi_upper, dPhi_2, L13_phi_bellow, Ele) && Signal_window_check(L13_eta_upper, dEta_2, L13_eta_bellow, Ele)  )
   L13_EM_Ele = 1;
if( Signal_window_check(L13_phi_upper, dPhi_2, L13_phi_bellow, Pos) && Signal_window_check(L13_eta_upper, dEta_2, L13_eta_bellow, Ele)  )
   L13_EM_Pos = 1;

if( Signal_window_check(L23_phi_upper, dPhi_3, L23_phi_bellow, Ele) && Signal_window_check(L23_eta_upper, dEta_3, L23_eta_bellow, Ele)  )
   L23_EM_Ele = 1;
if( Signal_window_check(L23_phi_upper, dPhi_3, L23_phi_bellow, Pos) && Signal_window_check(L23_eta_upper, dEta_3, L23_eta_bellow, Ele)  )
   L23_EM_Pos = 1;

if( Signal_window_check(L123_DPhi_cut1, dPhi, L123_DPhi_cut2, Ele ) && Signal_window_check(L123_DEta_cut1, dEta, L123_DEta_cut2, Ele )  )
   L123_pass_Ele = 1;
if( Signal_window_check(L123_DPhi_cut1, dPhi, L123_DPhi_cut2, Pos ) && Signal_window_check(L123_DEta_cut1, dEta, L123_DEta_cut2, Ele )  ) 
   L123_pass_Pos = 1;

}

void MuonEff::TriggeringWithout_3rdPixel( int nthFirstHit, int nthSecondHit, int nthThirdHit){

dPhi012 = StandaloneDPhi( 0, 1, 2, 0, nthFirstHit, nthSecondHit);
dPhi014 = StandaloneDPhi( 0, 1, 4, 0, nthFirstHit, nthThirdHit );
dPhi024 = StandaloneDPhi( 0, 2, 4, 0, nthSecondHit, nthThirdHit);

dPhi_1 = EMmatchingDPhi(first_disk_hits[nthFirstHit], second_disk_hits[nthSecondHit], emvector);
dEta_1 = EMmatchingDEta(first_disk_hits[nthFirstHit], second_disk_hits[nthSecondHit], emvector);

dPhi_2 = EMmatchingDPhi(first_disk_hits[nthFirstHit], fourth_disk_hits[nthThirdHit], emvector);
dEta_2 = EMmatchingDEta(first_disk_hits[nthFirstHit], fourth_disk_hits[nthThirdHit], emvector);

dPhi_3 = EMmatchingDPhi(second_disk_hits[nthSecondHit], fourth_disk_hits[nthThirdHit], emvector);
dEta_3 = EMmatchingDEta(second_disk_hits[nthSecondHit], fourth_disk_hits[nthThirdHit], emvector);

L012_pass_Ele = Signal_window_check( L012_DPhi_cut1, dPhi012, L012_DPhi_cut2, Ele );
L012_pass_Pos = Signal_window_check( L012_DPhi_cut1, dPhi012, L012_DPhi_cut2, Pos );
L014_pass_Ele = Signal_window_check( L014_DPhi_cut1, dPhi014, L014_DPhi_cut2, Ele );
L014_pass_Pos = Signal_window_check( L014_DPhi_cut1, dPhi014, L014_DPhi_cut2, Pos );
L024_pass_Ele = Signal_window_check( L024_DPhi_cut1, dPhi024, L024_DPhi_cut2, Ele );
L024_pass_Pos = Signal_window_check( L024_DPhi_cut1, dPhi024, L024_DPhi_cut2, Pos );

if( Signal_window_check(L12_phi_upper, dPhi_1, L12_phi_bellow, Ele) && Signal_window_check(L12_eta_upper, dEta_1, L12_eta_bellow, Ele)  ) 
   L12_EM_Ele = 1;
if( Signal_window_check(L12_phi_upper, dPhi_1, L12_phi_bellow, Pos) && Signal_window_check(L12_eta_upper, dEta_1, L12_eta_bellow, Ele)  )
   L12_EM_Pos = 1;

if( Signal_window_check(L14_phi_upper, dPhi_2, L14_phi_bellow, Ele) && Signal_window_check(L14_eta_upper, dEta_2, L14_eta_bellow, Ele)  )
   L14_EM_Ele = 1;
if( Signal_window_check(L14_phi_upper, dPhi_2, L14_phi_bellow, Pos) && Signal_window_check(L14_eta_upper, dEta_2, L14_eta_bellow, Ele)  )
   L14_EM_Pos = 1;

if( Signal_window_check(L24_phi_upper, dPhi_3, L24_phi_bellow, Ele) && Signal_window_check(L24_eta_upper, dEta_3, L24_eta_bellow, Ele)  ) 
   L24_EM_Ele = 1;
if( Signal_window_check(L24_phi_upper, dPhi_3, L24_phi_bellow, Pos) && Signal_window_check(L24_eta_upper, dEta_3, L24_eta_bellow, Ele)  )
   L24_EM_Pos = 1;

if( Signal_window_check(L124_DPhi_cut1, dPhi, L124_DPhi_cut2, Ele ) && Signal_window_check(L124_DEta_cut1, dEta, L124_DEta_cut2, Ele )  )
   L124_pass_Ele = 1;
if( Signal_window_check(L124_DPhi_cut1, dPhi, L124_DPhi_cut2, Pos ) && Signal_window_check(L124_DEta_cut1, dEta, L124_DEta_cut2, Ele )  )
   L124_pass_Pos = 1;
}

void MuonEff::TriggeringWithout_2ndPixel( int nthFirstHit, int nthSecondHit, int nthThirdHit){

dPhi013 = StandaloneDPhi( 0, 1, 3, 0, nthFirstHit, nthSecondHit);
dPhi014 = StandaloneDPhi( 0, 1, 4, 0, nthFirstHit, nthThirdHit );
dPhi034 = StandaloneDPhi( 0, 3, 4, 0, nthSecondHit, nthThirdHit);


dPhi_1 = EMmatchingDPhi(first_disk_hits[nthFirstHit], third_disk_hits[nthSecondHit], emvector);
dEta_1 = EMmatchingDEta(first_disk_hits[nthFirstHit], third_disk_hits[nthSecondHit], emvector);

dPhi_2 = EMmatchingDPhi(first_disk_hits[nthFirstHit], fourth_disk_hits[nthThirdHit], emvector);
dEta_2 = EMmatchingDEta(first_disk_hits[nthFirstHit], fourth_disk_hits[nthThirdHit], emvector);

dPhi_3 = EMmatchingDPhi(third_disk_hits[nthSecondHit], fourth_disk_hits[nthThirdHit], emvector);
dEta_3 = EMmatchingDEta(third_disk_hits[nthSecondHit], fourth_disk_hits[nthThirdHit], emvector);

L013_pass_Ele = Signal_window_check( L013_DPhi_cut1, dPhi013, L013_DPhi_cut2, Ele );
L013_pass_Pos = Signal_window_check( L013_DPhi_cut1, dPhi013, L013_DPhi_cut2, Pos );
L014_pass_Ele = Signal_window_check( L014_DPhi_cut1, dPhi014, L014_DPhi_cut2, Ele );
L014_pass_Pos = Signal_window_check( L014_DPhi_cut1, dPhi014, L014_DPhi_cut2, Pos );
L034_pass_Ele = Signal_window_check( L034_DPhi_cut1, dPhi034, L034_DPhi_cut2, Ele );
L034_pass_Pos = Signal_window_check( L034_DPhi_cut1, dPhi034, L034_DPhi_cut2, Pos );

if( Signal_window_check(L13_phi_upper, dPhi_1, L13_phi_bellow, Ele) && Signal_window_check(L13_eta_upper, dEta_1, L13_eta_bellow, Ele)  ) 
   L13_EM_Ele = 1;
if( Signal_window_check(L13_phi_upper, dPhi_1, L13_phi_bellow, Pos) && Signal_window_check(L13_eta_upper, dEta_1, L13_eta_bellow, Ele)  )
   L13_EM_Pos = 1;

if( Signal_window_check(L14_phi_upper, dPhi_2, L14_phi_bellow, Ele) && Signal_window_check(L14_eta_upper, dEta_2, L14_eta_bellow, Ele)  )
   L14_EM_Ele = 1;
if( Signal_window_check(L14_phi_upper, dPhi_2, L14_phi_bellow, Pos) && Signal_window_check(L14_eta_upper, dEta_2, L14_eta_bellow, Ele)  )
   L14_EM_Pos = 1;

if( Signal_window_check(L34_phi_upper, dPhi_3, L34_phi_bellow, Ele) && Signal_window_check(L34_eta_upper, dEta_3, L34_eta_bellow, Ele)  )
   L34_EM_Ele = 1;
if( Signal_window_check(L34_phi_upper, dPhi_3, L34_phi_bellow, Pos) && Signal_window_check(L34_eta_upper, dEta_3, L34_eta_bellow, Ele)  ) 
   L34_EM_Pos = 1;

if( Signal_window_check(L134_DPhi_cut1, dPhi, L134_DPhi_cut2, Ele ) && Signal_window_check(L134_DEta_cut1, dEta, L134_DEta_cut2, Ele )  )
  L134_pass_Ele = 1;
if( Signal_window_check(L134_DPhi_cut1, dPhi, L134_DPhi_cut2, Pos ) && Signal_window_check(L134_DEta_cut1, dEta, L134_DEta_cut2, Ele )  ) 
  L134_pass_Pos = 1;
}

void MuonEff::TriggeringWithout_1stPixel( int nthFirstHit, int nthSecondHit, int nthThirdHit){

dPhi023 = StandaloneDPhi( 0, 2, 3, 0, nthFirstHit, nthSecondHit);
dPhi024 = StandaloneDPhi( 0, 2, 4, 0, nthFirstHit, nthThirdHit );
dPhi034 = StandaloneDPhi( 0, 3, 4, 0, nthSecondHit, nthThirdHit);

dPhi_1 = EMmatchingDPhi(second_disk_hits[nthFirstHit], third_disk_hits[nthSecondHit], emvector);
dEta_1 = EMmatchingDEta(second_disk_hits[nthFirstHit], third_disk_hits[nthSecondHit], emvector);

dPhi_2 = EMmatchingDPhi(second_disk_hits[nthFirstHit], fourth_disk_hits[nthThirdHit], emvector);
dEta_2 = EMmatchingDEta(second_disk_hits[nthFirstHit], fourth_disk_hits[nthThirdHit], emvector);

dPhi_3 = EMmatchingDPhi(third_disk_hits[nthSecondHit], fourth_disk_hits[nthThirdHit], emvector);
dEta_3 = EMmatchingDEta(third_disk_hits[nthSecondHit], fourth_disk_hits[nthThirdHit], emvector);

L023_pass_Ele = Signal_window_check( L023_DPhi_cut1, dPhi023, L023_DPhi_cut2, Ele );
L023_pass_Pos = Signal_window_check( L023_DPhi_cut1, dPhi023, L023_DPhi_cut2, Pos );
L024_pass_Ele = Signal_window_check( L024_DPhi_cut1, dPhi024, L024_DPhi_cut2, Ele );
L024_pass_Pos = Signal_window_check( L024_DPhi_cut1, dPhi024, L024_DPhi_cut2, Pos );
L034_pass_Ele = Signal_window_check( L034_DPhi_cut1, dPhi034, L034_DPhi_cut2, Ele );
L034_pass_Pos = Signal_window_check( L034_DPhi_cut1, dPhi034, L034_DPhi_cut2, Pos );

if( Signal_window_check(L23_phi_upper, dPhi_1, L23_phi_bellow, Ele) && Signal_window_check(L23_eta_upper, dEta_1, L23_eta_bellow, Ele)  )
   L23_EM_Ele = 1;
if( Signal_window_check(L23_phi_upper, dPhi_1, L23_phi_bellow, Pos) && Signal_window_check(L23_eta_upper, dEta_1, L23_eta_bellow, Ele)  )
   L23_EM_Pos = 1;

if( Signal_window_check(L24_phi_upper, dPhi_2, L24_phi_bellow, Ele) && Signal_window_check(L24_eta_upper, dEta_2, L24_eta_bellow, Ele)  )
   L24_EM_Ele = 1;
if( Signal_window_check(L24_phi_upper, dPhi_2, L24_phi_bellow, Pos) && Signal_window_check(L24_eta_upper, dEta_2, L24_eta_bellow, Ele)  )
   L24_EM_Pos = 1;

if( Signal_window_check(L34_phi_upper, dPhi_3, L34_phi_bellow, Ele) && Signal_window_check(L34_eta_upper, dEta_3, L34_eta_bellow, Ele)  )
   L34_EM_Ele = 1;
if( Signal_window_check(L34_phi_upper, dPhi_3, L34_phi_bellow, Pos) && Signal_window_check(L34_eta_upper, dEta_3, L34_eta_bellow, Ele)  )
   L34_EM_Pos = 1;

if( Signal_window_check(L234_DPhi_cut1, dPhi, L234_DPhi_cut2, Ele ) && Signal_window_check(L234_DEta_cut1, dEta, L234_DEta_cut2, Ele )  )
   L234_pass_Ele = 1;
if( Signal_window_check(L234_DPhi_cut1, dPhi, L234_DPhi_cut2, Pos ) && Signal_window_check(L234_DEta_cut1, dEta, L234_DEta_cut2, Ele )  )
   L234_pass_Pos = 1;
}
void MuonEff::TriggeringWith_1st2ndPixel(int nthFirstHit, int nthSecondHit){
dPhi012 = StandaloneDPhi( 0, 1, 2, 0, nthFirstHit, nthSecondHit );
L012_pass_Ele = Signal_window_check( L012_DPhi_cut1, dPhi012, L012_DPhi_cut2, Ele );
L012_pass_Pos = Signal_window_check( L012_DPhi_cut1, dPhi012, L012_DPhi_cut2, Pos );

if( L012_pass_Ele &&
    (first_disk_hits_Ele_or_Pos[nthFirstHit] == 1 || first_disk_hits_Ele_or_Pos[nthFirstHit] ==3) &&
    (second_disk_hits_Ele_or_Pos[nthSecondHit] == 1 || second_disk_hits_Ele_or_Pos[nthSecondHit] ==3)) _pass_Ele = 1;
if( L012_pass_Pos &&
    (first_disk_hits_Ele_or_Pos[nthFirstHit] == 2 || first_disk_hits_Ele_or_Pos[nthFirstHit] ==3) &&
    (second_disk_hits_Ele_or_Pos[nthSecondHit] == 2 || second_disk_hits_Ele_or_Pos[nthSecondHit] ==3)) _pass_Pos = 1;
}

void MuonEff::TriggeringWith_1st3rdPixel(int nthFirstHit, int nthSecondHit){
dPhi013 = StandaloneDPhi( 0, 1, 3, 0, nthFirstHit, nthSecondHit );
L013_pass_Ele = Signal_window_check( L013_DPhi_cut1, dPhi013, L013_DPhi_cut2, Ele );
L013_pass_Pos = Signal_window_check( L013_DPhi_cut1, dPhi013, L013_DPhi_cut2, Pos );

 if( L013_pass_Ele &&
     (first_disk_hits_Ele_or_Pos[nthFirstHit] == 1 || first_disk_hits_Ele_or_Pos[nthFirstHit] ==3) &&
     (third_disk_hits_Ele_or_Pos[nthSecondHit] == 1 || third_disk_hits_Ele_or_Pos[nthSecondHit] ==3)) _pass_Ele = 1;
 if( L013_pass_Pos &&
     (first_disk_hits_Ele_or_Pos[nthFirstHit] == 2 || first_disk_hits_Ele_or_Pos[nthFirstHit] ==3) &&
     (third_disk_hits_Ele_or_Pos[nthSecondHit] == 2 || third_disk_hits_Ele_or_Pos[nthSecondHit] ==3)) _pass_Pos = 1;

}
void MuonEff::TriggeringWith_2nd3rdPixel(int nthFirstHit, int nthSecondHit){
dPhi023 = StandaloneDPhi( 0, 2, 3, 0, nthFirstHit, nthSecondHit );
L023_pass_Ele = Signal_window_check( L023_DPhi_cut1, dPhi023, L023_DPhi_cut2, Ele );
L023_pass_Pos = Signal_window_check( L023_DPhi_cut1, dPhi023, L023_DPhi_cut2, Pos );

 if( L023_pass_Ele &&
     (second_disk_hits_Ele_or_Pos[nthFirstHit] == 1 || second_disk_hits_Ele_or_Pos[nthFirstHit] ==3) &&
     (third_disk_hits_Ele_or_Pos[nthSecondHit] == 1 || third_disk_hits_Ele_or_Pos[nthSecondHit] ==3)) _pass_Ele = 1;
 if( L023_pass_Pos &&
     (second_disk_hits_Ele_or_Pos[nthFirstHit] == 2 || second_disk_hits_Ele_or_Pos[nthFirstHit] ==3) &&
     (third_disk_hits_Ele_or_Pos[nthSecondHit] == 2 || third_disk_hits_Ele_or_Pos[nthSecondHit] ==3)) _pass_Pos = 1;
}

void MuonEff::TriggeringWith_1st2ndPixel_v2(int nthFirstHit, int nthSecondHit){
dPhi012 = StandaloneDPhi( 0, 1, 2, 0, nthFirstHit, nthSecondHit );
L012_pass_Ele = Signal_window_check( L012_DPhi_cut1, dPhi012, L012_DPhi_cut2, Ele );
L012_pass_Pos = Signal_window_check( L012_DPhi_cut1, dPhi012, L012_DPhi_cut2, Pos );

dPhi_1 = EMmatchingDPhi(first_disk_hits[nthFirstHit], second_disk_hits[nthSecondHit], emvector);
dEta_1 = EMmatchingDEta(first_disk_hits[nthFirstHit], second_disk_hits[nthSecondHit], emvector);

if( Signal_window_check(L12_phi_upper, dPhi_1, L12_phi_bellow, Ele) && Signal_window_check(L12_eta_upper, dEta_1, L12_eta_bellow, Ele)  )
   L12_EM_Ele = 1;
if( Signal_window_check(L12_phi_upper, dPhi_1, L12_phi_bellow, Pos) && Signal_window_check(L12_eta_upper, dEta_1, L12_eta_bellow, Ele)  )
   L12_EM_Pos = 1;

if( L012_pass_Ele && L12_EM_Ele &&
    (first_disk_hits_Ele_or_Pos[nthFirstHit] == 1 || first_disk_hits_Ele_or_Pos[nthFirstHit] ==3) &&
    (second_disk_hits_Ele_or_Pos[nthSecondHit] == 1 || second_disk_hits_Ele_or_Pos[nthSecondHit] ==3)) _pass_Ele = 1;
if( L012_pass_Pos && L12_EM_Pos &&
    (first_disk_hits_Ele_or_Pos[nthFirstHit] == 2 || first_disk_hits_Ele_or_Pos[nthFirstHit] ==3) &&
    (second_disk_hits_Ele_or_Pos[nthSecondHit] == 2 || second_disk_hits_Ele_or_Pos[nthSecondHit] ==3)) _pass_Pos = 1;
}

void MuonEff::TriggeringWith_1st3rdPixel_v2(int nthFirstHit, int nthSecondHit){
dPhi013 = StandaloneDPhi( 0, 1, 2, 0, nthFirstHit, nthSecondHit );
L013_pass_Ele = Signal_window_check( L013_DPhi_cut1, dPhi013, L013_DPhi_cut2, Ele );
L013_pass_Pos = Signal_window_check( L013_DPhi_cut1, dPhi013, L013_DPhi_cut2, Pos );

dPhi_1 = EMmatchingDPhi(first_disk_hits[nthFirstHit], third_disk_hits[nthSecondHit], emvector);
dEta_1 = EMmatchingDEta(first_disk_hits[nthFirstHit], third_disk_hits[nthSecondHit], emvector);

if( Signal_window_check(L13_phi_upper, dPhi_1, L13_phi_bellow, Ele) && Signal_window_check(L13_eta_upper, dEta_1, L13_eta_bellow, Ele)  )
   L13_EM_Ele = 1;
if( Signal_window_check(L13_phi_upper, dPhi_1, L13_phi_bellow, Pos) && Signal_window_check(L13_eta_upper, dEta_1, L13_eta_bellow, Ele)  )
   L13_EM_Pos = 1;

if( L013_pass_Ele && L13_EM_Ele &&
    (first_disk_hits_Ele_or_Pos[nthFirstHit] == 1 || first_disk_hits_Ele_or_Pos[nthFirstHit] ==3) &&
    (third_disk_hits_Ele_or_Pos[nthSecondHit] == 1 || third_disk_hits_Ele_or_Pos[nthSecondHit] ==3)) _pass_Ele = 1;
if( L013_pass_Pos && L13_EM_Pos &&
    (first_disk_hits_Ele_or_Pos[nthFirstHit] == 2 || first_disk_hits_Ele_or_Pos[nthFirstHit] ==3) &&
    (third_disk_hits_Ele_or_Pos[nthSecondHit] == 2 || third_disk_hits_Ele_or_Pos[nthSecondHit] ==3)) _pass_Pos = 1;
}

void MuonEff::TriggeringWith_1st4thPixel_v2(int nthFirstHit, int nthSecondHit){
dPhi014 = StandaloneDPhi( 0, 1, 2, 0, nthFirstHit, nthSecondHit );
L014_pass_Ele = Signal_window_check( L014_DPhi_cut1, dPhi014, L014_DPhi_cut2, Ele );
L014_pass_Pos = Signal_window_check( L014_DPhi_cut1, dPhi014, L014_DPhi_cut2, Pos );

dPhi_1 = EMmatchingDPhi(first_disk_hits[nthFirstHit], fourth_disk_hits[nthSecondHit], emvector);
dEta_1 = EMmatchingDEta(first_disk_hits[nthFirstHit], fourth_disk_hits[nthSecondHit], emvector);

if( Signal_window_check(L14_phi_upper, dPhi_1, L14_phi_bellow, Ele) && Signal_window_check(L14_eta_upper, dEta_1, L14_eta_bellow, Ele)  )
   L14_EM_Ele = 1;
if( Signal_window_check(L14_phi_upper, dPhi_1, L14_phi_bellow, Pos) && Signal_window_check(L14_eta_upper, dEta_1, L14_eta_bellow, Ele)  )
   L14_EM_Pos = 1;

if( L014_pass_Ele && L14_EM_Ele &&
    (first_disk_hits_Ele_or_Pos[nthFirstHit] == 1 || first_disk_hits_Ele_or_Pos[nthFirstHit] ==3) &&
    (fourth_disk_hits_Ele_or_Pos[nthSecondHit] == 1 || fourth_disk_hits_Ele_or_Pos[nthSecondHit] ==3)) _pass_Ele = 1;
if( L014_pass_Pos && L14_EM_Pos &&
    (first_disk_hits_Ele_or_Pos[nthFirstHit] == 2 || first_disk_hits_Ele_or_Pos[nthFirstHit] ==3) &&
    (fourth_disk_hits_Ele_or_Pos[nthSecondHit] == 2 || fourth_disk_hits_Ele_or_Pos[nthSecondHit] ==3)) _pass_Pos = 1;
}

void MuonEff::TriggeringWith_2nd3rdPixel_v2(int nthFirstHit, int nthSecondHit){
dPhi023 = StandaloneDPhi( 0, 1, 2, 0, nthFirstHit, nthSecondHit );
L023_pass_Ele = Signal_window_check( L023_DPhi_cut1, dPhi023, L023_DPhi_cut2, Ele );
L023_pass_Pos = Signal_window_check( L023_DPhi_cut1, dPhi023, L023_DPhi_cut2, Pos );

dPhi_1 = EMmatchingDPhi(second_disk_hits[nthFirstHit], third_disk_hits[nthSecondHit], emvector);
dEta_1 = EMmatchingDEta(second_disk_hits[nthFirstHit], third_disk_hits[nthSecondHit], emvector);

if( Signal_window_check(L23_phi_upper, dPhi_1, L23_phi_bellow, Ele) && Signal_window_check(L23_eta_upper, dEta_1, L23_eta_bellow, Ele)  )
   L23_EM_Ele = 1;
if( Signal_window_check(L23_phi_upper, dPhi_1, L23_phi_bellow, Pos) && Signal_window_check(L23_eta_upper, dEta_1, L23_eta_bellow, Ele)  )
   L23_EM_Pos = 1;

if( L023_pass_Ele && L23_EM_Ele &&
    (first_disk_hits_Ele_or_Pos[nthFirstHit] == 1 || first_disk_hits_Ele_or_Pos[nthFirstHit] ==3) &&
    (third_disk_hits_Ele_or_Pos[nthSecondHit] == 1 || third_disk_hits_Ele_or_Pos[nthSecondHit] ==3)) _pass_Ele = 1;
if( L023_pass_Pos && L23_EM_Pos &&
    (first_disk_hits_Ele_or_Pos[nthFirstHit] == 2 || first_disk_hits_Ele_or_Pos[nthFirstHit] ==3) &&
    (third_disk_hits_Ele_or_Pos[nthSecondHit] == 2 || third_disk_hits_Ele_or_Pos[nthSecondHit] ==3)) _pass_Pos = 1;
}

void MuonEff::TriggeringWith_2nd4thPixel_v2(int nthFirstHit, int nthSecondHit){
dPhi024 = StandaloneDPhi( 0, 1, 2, 0, nthFirstHit, nthSecondHit );
L024_pass_Ele = Signal_window_check( L024_DPhi_cut1, dPhi024, L024_DPhi_cut2, Ele );
L024_pass_Pos = Signal_window_check( L024_DPhi_cut1, dPhi024, L024_DPhi_cut2, Pos );

dPhi_1 = EMmatchingDPhi(second_disk_hits[nthFirstHit], fourth_disk_hits[nthSecondHit], emvector);
dEta_1 = EMmatchingDEta(second_disk_hits[nthFirstHit], fourth_disk_hits[nthSecondHit], emvector);

if( Signal_window_check(L24_phi_upper, dPhi_1, L24_phi_bellow, Ele) && Signal_window_check(L24_eta_upper, dEta_1, L24_eta_bellow, Ele)  )
   L24_EM_Ele = 1;
if( Signal_window_check(L24_phi_upper, dPhi_1, L24_phi_bellow, Pos) && Signal_window_check(L24_eta_upper, dEta_1, L24_eta_bellow, Ele)  )
   L24_EM_Pos = 1;

if( L024_pass_Ele && L24_EM_Ele &&
    (first_disk_hits_Ele_or_Pos[nthFirstHit] == 1 || first_disk_hits_Ele_or_Pos[nthFirstHit] ==3) &&
    (third_disk_hits_Ele_or_Pos[nthSecondHit] == 1 || third_disk_hits_Ele_or_Pos[nthSecondHit] ==3)) _pass_Ele = 1;
if( L024_pass_Pos && L24_EM_Pos &&
    (first_disk_hits_Ele_or_Pos[nthFirstHit] == 2 || first_disk_hits_Ele_or_Pos[nthFirstHit] ==3) &&
    (third_disk_hits_Ele_or_Pos[nthSecondHit] == 2 || third_disk_hits_Ele_or_Pos[nthSecondHit] ==3)) _pass_Pos = 1;
}

void MuonEff::TriggeringWith_3rd4thPixel_v2(int nthFirstHit, int nthSecondHit){
dPhi034 = StandaloneDPhi( 0, 1, 2, 0, nthFirstHit, nthSecondHit );
L034_pass_Ele = Signal_window_check( L034_DPhi_cut1, dPhi034, L034_DPhi_cut2, Ele );
L034_pass_Pos = Signal_window_check( L034_DPhi_cut1, dPhi034, L034_DPhi_cut2, Pos );

dPhi_1 = EMmatchingDPhi(third_disk_hits[nthFirstHit], fourth_disk_hits[nthSecondHit], emvector);
dEta_1 = EMmatchingDEta(third_disk_hits[nthFirstHit], fourth_disk_hits[nthSecondHit], emvector);

if( Signal_window_check(L34_phi_upper, dPhi_1, L34_phi_bellow, Ele) && Signal_window_check(L34_eta_upper, dEta_1, L34_eta_bellow, Ele)  )
   L34_EM_Ele = 1;
if( Signal_window_check(L34_phi_upper, dPhi_1, L34_phi_bellow, Pos) && Signal_window_check(L34_eta_upper, dEta_1, L34_eta_bellow, Ele)  )
   L34_EM_Pos = 1;

if( L034_pass_Ele && L34_EM_Ele &&
    (first_disk_hits_Ele_or_Pos[nthFirstHit] == 1 || first_disk_hits_Ele_or_Pos[nthFirstHit] ==3) &&
    (third_disk_hits_Ele_or_Pos[nthSecondHit] == 1 || third_disk_hits_Ele_or_Pos[nthSecondHit] ==3)) _pass_Ele = 1;
if( L034_pass_Pos && L34_EM_Pos &&
    (first_disk_hits_Ele_or_Pos[nthFirstHit] == 2 || first_disk_hits_Ele_or_Pos[nthFirstHit] ==3) &&
    (third_disk_hits_Ele_or_Pos[nthSecondHit] == 2 || third_disk_hits_Ele_or_Pos[nthSecondHit] ==3)) _pass_Pos = 1;
}

#endif // #ifdef MuonEff_cxx
