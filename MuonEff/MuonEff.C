#define MuonEff_cxx
#include "MuonEff.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <string>
#include <TLorentzVector.h>
#include <TGraph.h>

using namespace std;

void MuonEff::Loop()
{
   if (fChain == 0) return;

//   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nentries = fChain->GetEntries();

   Long64_t nbytes = 0, nb = 0;

   float dr_cut = 0.1;

   int flag = 0;

   bit1 = 0x1;

   debug = false;

   const double ME0_PiX_dphi_width_[9] = {0.02, 0.04, 0.031, 0.033, 0.035, 0.037, 0.039, 0.04, 0.05}; //last: 0.02
   const double ME0_PiX_deta_width_[9] = {0.05, 0.03, 0.031, 0.033, 0.035, 0.037, 0.039, 0.04, 0.05}; //last: 0.05

   const double PiX_PiX_dphi_width_[9] = {0.002, 0.0035, 0.0031, 0.0033, 0.0035, 0.0037, 0.0039, 0.004, 0.005}; // last: 0.002
   const double PiX_PiX_deta_width_[9] = {0.005, 0.008, 0.0031, 0.0033, 0.0035, 0.0037, 0.0039, 0.004, 0.005}; // last: 0.005, 0.008

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

	  if(!(jentry%1)) cout<<"Processing entry " << jentry << "/" << nentries <<endl;
	  FillCutFlow("NoCut",1.);

      matchedME0Et  = -999.;
      matchedME0Eta = -999.;
      matchedME0Phi = -999.;
      fired        = 0;

      ntnME02 = 0;
      ntME0Et.clear();
      ntME0Eta.clear();
      ntME0Phi.clear();

      PiXTRKbit.clear();
      trigger_bit_width.clear();
      pix_comb.clear();

      ntCl_match.clear();
      isTrack_match.clear();
      chi2.clear();
      track_dr.clear();
      withoutEM_match.clear();
      withEM_match.clear();

      ntfirstPix.clear();
      ntsecondPix.clear();
      ntthirdPix.clear();
      ntfourthPix.clear();

      nt_genPhi = genPartPhi->at(0);
      nt_genEta = genPartEta->at(0);
      nt_genPt = genPartPt->at(0);

      nt_lastSimtkpt = lastSimtkpt;
      nt_initialSimtkpt = initialSimtkpt;
 
      float tempDR = 999.;
//      int   indx = -1;
//      int   indx_closestME0 = -1;

	  //find closest me0 object to the gen electron
	  me0N=me0SegPosX->size(); 
	  
	  float closest_dr = 9999.;
	  int closest_me0 = -1;
	  int me0Count = 0;

      for( int i = 0; i < me0N; i++){
	      TVector3 me0SegPos;
	      me0SegPos.SetXYZ(me0SegPosX->at(i), me0SegPosY->at(i), me0SegPosZ->at(i));
	      float me0Phi = me0SegPos.Phi();
	      float me0Eta = me0SegPos.Eta();

	      float dPhi = deltaPhi(genPartPhi->at(0), me0Phi);
	      float current_dr = sqrt( pow(dPhi,2) + pow(genPartEta->at(0) - me0Eta,2) );
	      if( genPartPt->at(0) < 10 ) continue; // !!!!! Use genPartPt instead me0 Pt !!!!!
	      me0Count++;
	      if(current_dr < closest_dr){
		      closest_dr = current_dr;
		      closest_me0 = i;
	      }
      }// end of loop to find the closest me0 the gen

	  pix_comb_ = 0x0;

	  nPix123_segments = 0;
	  nPix124_segments = 0;
	  nPix134_segments = 0;
	  nPix234_segments = 0;

	  if( closest_dr == 9999. || closest_me0 == -1 ) continue;

	  // find me0 passing pixtrk signal windows
	  if( closest_dr < dr_cut && closest_dr != 9999. ){
		  debug = false;

//		  indx++; //remove this variable

		  float me0X = me0SegPosX->at(closest_me0);
		  float me0Y = me0SegPosY->at(closest_me0);
		  float me0Z = me0SegPosZ->at(closest_me0);

		  TVector3 me0SegPos;
		  me0SegPos.SetXYZ(me0X, me0Y, me0Z);
		  ME0Et = genPartPt->at(0); // !!!!!!!! genPartPt instead me0 !!!!!!!
		  ME0Eta = me0SegPos.Eta();
		  ME0Phi = me0SegPos.Phi();

//cout<<"isTrackMatched->at "<<isTrackMatched->at(1)<<endl;
//		  isTrack_match.push_back(isTrackMatched->at(closest_me0));
//		  chi2.push_back(trackHighestPtCutChi2Chi2->at(closest_me0));
//		  track_dr.push_back(trackmatchingdR->at(closest_me0));

		  emvector.SetXYZ(me0X, me0Y, me0Z);

		  tempDR = closest_dr; 
		  matchedME0Eta = ME0Eta;
		  matchedME0Phi = ME0Phi;
		  matchedME0Et  = ME0Et;
//		  indx_closestME0 = indx;
cout<<"region "<<eta_region<<endl;
		  if( ME0Eta > 2.0 && ME0Eta < 2.4 ) eta_region = 1; 
		  if( ME0Eta > 2.4 && ME0Eta < 2.8 ) eta_region = 2;
		  if( ME0Eta < 2.0 || ME0Eta > 2.8 ) continue;

		  ntnME02++;
		  ntME0Et.push_back(ME0Et);
		  ntME0Eta.push_back(ME0Eta);
		  ntME0Phi.push_back(ME0Phi);

		  // set region of interest
		  SetROI(eta_region);

		  //initialize pixel hit variables
		  first_disk_hits.clear();
		  second_disk_hits.clear();
		  third_disk_hits.clear();
		  fourth_disk_hits.clear();

		  first_disk_hits_Ele_or_Pos.clear();
		  second_disk_hits_Ele_or_Pos.clear();
		  third_disk_hits_Ele_or_Pos.clear();
		  fourth_disk_hits_Ele_or_Pos.clear();
		  

		  hitted_layers.clear();

		  layers[0] = 1; // beam spot
		  layers[1] = 0; layers[2] = 0; layers[3] = 0; layers[4] = 0;
		  r = 0;

		  StorePixelHit(eta_region); // save pixel hits in Region of Interest for the given eta region

		  //check which pixel has hits
		  for( int i=1; i < 5; i++){
			  if( layers[i] != 0 ) hitted_layers.push_back(i);
		  }

		  int global_index_width = 0;
		  trigger_bit_width_ = 0x0;
		  // set pixtrk signal boundary
		  for( int nth_me0_pix_deta = 0; nth_me0_pix_deta < 9; nth_me0_pix_deta++){
			  if( nth_me0_pix_deta != 0) continue;
		  
		  if(eta_region == 1) SetSignalBoundary(eta_region, ME0_PiX_dphi_width_[0], ME0_PiX_deta_width_[0], PiX_PiX_dphi_width_[0], PiX_PiX_deta_width_[0]);
		  if(eta_region == 2) SetSignalBoundary(eta_region, ME0_PiX_dphi_width_[0], ME0_PiX_deta_width_[0], PiX_PiX_dphi_width_[0], PiX_PiX_deta_width_[1]);

		  // PixTRK algorithm 
		  PixTrkPassed = false;
		  withoutEM_count_Ele = 0, withEM_count_Ele = 0;

		  fourth_layer_missing = 0;
		  third_layer_missing = 0;
		  second_layer_missing = 0;
		  first_layer_missing = 0;

		  // loop over every 3 out of 4 pixel combination 
		  for( std::vector<int>::iterator first_hit = hitted_layers.begin(); first_hit != hitted_layers.end(); first_hit++){
			  for ( std::vector<int>::iterator second_hit = first_hit+1; second_hit != hitted_layers.end(); second_hit++){
				  for ( std::vector<int>::iterator third_hit = second_hit+1; third_hit != hitted_layers.end(); third_hit++){

					  // loop over every pixel hits in the given pixel combination
					  for( int k=0; k < layers[*first_hit]; k++){
						  for( int i=0; i < layers[*second_hit]; i++){
							  _pass_Ele = 0, _pass_Pos = 0;
							  L012_pass_Ele = 0, L012_pass_Pos = 0;
							  L013_pass_Ele = 0, L013_pass_Pos = 0;
							  L023_pass_Ele = 0, L023_pass_Pos = 0;

							  if( *first_hit == 1 && *second_hit == 2 )
								  TriggeringWith_1st2ndPixel(k,i);

							  if( *first_hit == 1 && *second_hit == 3 )
								  TriggeringWith_1st3rdPixel(k,i);

							  if( *first_hit == 2 && *second_hit == 3 )
								  TriggeringWith_2nd3rdPixel(k,i);

							  // skip only if both _pass_Ele and _pass_Pos are 0 i.e., both electron and positron signal window are not satisfied
							  if( !_pass_Ele && !_pass_Pos ) continue;

							  for( int j=0; j < layers[*third_hit]; j++){
								  all_cut_pass_Ele = 0, all_cut_pass_Pos = 0;
								  withoutEM_pass_Ele = 0, withEM_pass_Ele = 0;

								  L012_pass_Ele = 0, L012_pass_Pos = 0;
								  L013_pass_Ele = 0, L013_pass_Pos = 0;
								  L014_pass_Ele = 0, L014_pass_Pos = 0;
								  L023_pass_Ele = 0, L023_pass_Pos = 0;
								  L024_pass_Ele = 0, L024_pass_Pos = 0;
								  L034_pass_Ele = 0, L034_pass_Pos = 0;
								  L123_pass_Ele = 0, L123_pass_Pos = 0;
								  L124_pass_Ele = 0, L124_pass_Pos = 0;
								  L134_pass_Ele = 0, L134_pass_Pos = 0;
								  L234_pass_Ele = 0, L234_pass_Pos = 0;

								  L12_EM_Ele = 0, L12_EM_Pos = 0;
								  L13_EM_Ele = 0, L13_EM_Pos = 0;
								  L14_EM_Ele = 0, L14_EM_Pos = 0;
								  L23_EM_Ele = 0, L23_EM_Pos = 0;
								  L24_EM_Ele = 0, L24_EM_Pos = 0;
								  L34_EM_Ele = 0, L34_EM_Pos = 0;

								  dPhi = StandaloneDPhi( *first_hit, *second_hit, *third_hit, k, i, j );
								  dEta = StandaloneDEta( *first_hit, *second_hit, *third_hit, k, i, j );

								  if( *first_hit == 1 && *second_hit == 2 && *third_hit == 3 ){ // for efficiency counting  !!caution of the position of this codition
									  // This is for the case that the first hit is in the first pixel layer and the second hit is in the second pixel layer and the third hit is in the third layer. 
									  TriggeringWithout_4thPixel(k, i, j);

									  if( (first_disk_hits_Ele_or_Pos[k] == 1 || first_disk_hits_Ele_or_Pos[k] ==3) &&
											  (second_disk_hits_Ele_or_Pos[i] == 1 || second_disk_hits_Ele_or_Pos[i] ==3) &&
											  (third_disk_hits_Ele_or_Pos[j] == 1 || third_disk_hits_Ele_or_Pos[j] ==3)){
										  if(L012_pass_Ele && L013_pass_Ele && L023_pass_Ele && L123_pass_Ele && L12_EM_Ele && L13_EM_Ele && L23_EM_Ele)
											  all_cut_pass_Ele = 1; 
										  if(L012_pass_Ele && L013_pass_Ele && L023_pass_Ele && L123_pass_Ele)
											  withoutEM_pass_Ele = 1;
										  if(L12_EM_Ele && L13_EM_Ele && L23_EM_Ele)
											  withEM_pass_Ele = 1;
									  }

									  if( L012_pass_Pos && L013_pass_Pos && L023_pass_Pos && L123_pass_Pos && L12_EM_Pos && L13_EM_Pos && L23_EM_Pos &&
											  (first_disk_hits_Ele_or_Pos[k] == 2 || first_disk_hits_Ele_or_Pos[k] ==3) &&
											  (second_disk_hits_Ele_or_Pos[i] == 2 || second_disk_hits_Ele_or_Pos[i] ==3) &&
											  (third_disk_hits_Ele_or_Pos[j] == 2 || third_disk_hits_Ele_or_Pos[j] ==3)) all_cut_pass_Pos = 1; 

									  if( all_cut_pass_Ele == 1){
										  pix_comb_ = pix_comb_ | (bit1 << 1);
										  nPix123_segments++;
										  flag = 1;
									  }

									  if(skip){
										  if( all_cut_pass_Ele == 1 || all_cut_pass_Pos == 1 ){ // if pass exit loop
											  k = layers[*first_hit];
											  i = layers[*second_hit];
											  j = layers[*third_hit]; 
											  fourth_layer_missing = 1;
										  }
									  }
								  }// 1st,2nd,3rd case
								  if( *first_hit == 1 && *second_hit == 2 && *third_hit == 4 ){ // for efficiency counting  !!caution of the position of this codition
									  // This is for the case that the first hit is in the first pixel layer and the second is in the second layer and the third hit is in the fourth layer.
									  TriggeringWithout_3rdPixel(k, i, j);

									  if( (first_disk_hits_Ele_or_Pos[k] == 1 || first_disk_hits_Ele_or_Pos[k] ==3) &&
											  (second_disk_hits_Ele_or_Pos[i] == 1 || second_disk_hits_Ele_or_Pos[i] ==3) &&
											  (fourth_disk_hits_Ele_or_Pos[j] == 1 || fourth_disk_hits_Ele_or_Pos[j] ==3)){
										  if(L012_pass_Ele && L014_pass_Ele && L024_pass_Ele && L124_pass_Ele && L12_EM_Ele && L14_EM_Ele && L24_EM_Ele) all_cut_pass_Ele = 1;
										  if(L012_pass_Ele && L014_pass_Ele && L024_pass_Ele && L124_pass_Ele) withoutEM_pass_Ele = 1;
										  if(L12_EM_Ele && L14_EM_Ele && L24_EM_Ele) withEM_pass_Ele = 1;
									  } 

									  if( L012_pass_Pos && L014_pass_Pos && L024_pass_Pos && L124_pass_Pos && L12_EM_Pos && L14_EM_Pos && L24_EM_Pos &&
											  (first_disk_hits_Ele_or_Pos[k] == 2 || first_disk_hits_Ele_or_Pos[k] ==3) &&
											  (second_disk_hits_Ele_or_Pos[i] == 2 || second_disk_hits_Ele_or_Pos[i] ==3) &&
											  (fourth_disk_hits_Ele_or_Pos[j] == 2 || fourth_disk_hits_Ele_or_Pos[j] ==3)) all_cut_pass_Pos = 1; 

									  if( all_cut_pass_Ele == 1){
										  pix_comb_ = pix_comb_ | (bit1 << 2);
										  nPix124_segments++;
										  flag = 1;
									  }

									  if(skip){
										  if( all_cut_pass_Ele == 1 || all_cut_pass_Pos == 1 ){ // if pass exit the for loops
											  k = layers[*first_hit];
											  i = layers[*second_hit];
											  j = layers[*third_hit]; 
											  third_layer_missing = 1;
										  }
									  }
								  }// 1st,2nd,4th case

								  if( *first_hit == 1 && *second_hit == 3 && *third_hit == 4 ){ // for efficiency counting  !!caution of the position of this codition
									  TriggeringWithout_2ndPixel(k, i, j);

									  if( (first_disk_hits_Ele_or_Pos[k] == 1 || first_disk_hits_Ele_or_Pos[k] ==3) &&
											  (third_disk_hits_Ele_or_Pos[i] == 1 || third_disk_hits_Ele_or_Pos[i] ==3) &&
											  (fourth_disk_hits_Ele_or_Pos[j] == 1 || fourth_disk_hits_Ele_or_Pos[j] ==3)){
										  if(L013_pass_Ele && L014_pass_Ele && L034_pass_Ele && L134_pass_Ele && L13_EM_Ele && L14_EM_Ele && L34_EM_Ele)all_cut_pass_Ele = 1; 
										  if(L013_pass_Ele && L014_pass_Ele && L034_pass_Ele && L134_pass_Ele) withoutEM_pass_Ele = 1;
										  if(L13_EM_Ele && L14_EM_Ele && L34_EM_Ele) withEM_pass_Ele = 1;
									  }

									  if( L013_pass_Pos && L014_pass_Pos && L034_pass_Pos && L134_pass_Pos && L13_EM_Pos && L14_EM_Pos && L34_EM_Pos &&
											  (first_disk_hits_Ele_or_Pos[k] == 2 || first_disk_hits_Ele_or_Pos[k] ==3) &&
											  (third_disk_hits_Ele_or_Pos[i] == 2 || third_disk_hits_Ele_or_Pos[i] ==3) &&
											  (fourth_disk_hits_Ele_or_Pos[j] == 2 || fourth_disk_hits_Ele_or_Pos[j] ==3)) all_cut_pass_Pos = 1; 

									  if( all_cut_pass_Ele == 1){
										  pix_comb_ = pix_comb_ | (bit1 << 3);
										  nPix134_segments++;
										  flag = 1;
									  }
									  if(skip){
										  if( all_cut_pass_Ele == 1 || all_cut_pass_Pos == 1 ){ // if pass exit the for loops
											  k = layers[*first_hit];
											  i = layers[*second_hit];
											  j = layers[*third_hit]; 
											  second_layer_missing = 1; 
										  }
									  }
								  }// 1st,3rd,4th case

								  if( *first_hit == 2 && *second_hit == 3 && *third_hit == 4 ){ // for efficiency counting  !!caution of the position of this codition
									  TriggeringWithout_1stPixel(k, i, j);

									  if( (second_disk_hits_Ele_or_Pos[k] == 1 || second_disk_hits_Ele_or_Pos[k] ==3) &&
											  (third_disk_hits_Ele_or_Pos[i] == 1 || third_disk_hits_Ele_or_Pos[i] ==3) &&
											  (fourth_disk_hits_Ele_or_Pos[j] == 1 || fourth_disk_hits_Ele_or_Pos[j] ==3)){
										  if(L023_pass_Ele && L024_pass_Ele && L034_pass_Ele && L234_pass_Ele && L23_EM_Ele && L24_EM_Ele && L34_EM_Ele) all_cut_pass_Ele = 1;
										  if(L023_pass_Ele && L024_pass_Ele && L034_pass_Ele && L234_pass_Ele) withoutEM_pass_Ele = 1;
										  if(L23_EM_Ele && L24_EM_Ele && L34_EM_Ele) withEM_pass_Ele = 1;
									  } 

									  if( L023_pass_Pos && L024_pass_Pos && L034_pass_Pos && L234_pass_Pos && L23_EM_Pos && L24_EM_Pos && L34_EM_Pos &&
											  (second_disk_hits_Ele_or_Pos[k] == 2 || second_disk_hits_Ele_or_Pos[k] ==3) &&
											  (third_disk_hits_Ele_or_Pos[i] == 2 || third_disk_hits_Ele_or_Pos[i] ==3) &&
											  (fourth_disk_hits_Ele_or_Pos[j] == 2 || fourth_disk_hits_Ele_or_Pos[j] ==3)) all_cut_pass_Pos = 1; 

									  if( all_cut_pass_Ele == 1){
										  pix_comb_ = pix_comb_ | (bit1 << 4);
										  nPix234_segments++;
										  flag = 1;
									  } 
									  if(skip){
										  if( all_cut_pass_Ele == 1 || all_cut_pass_Pos == 1 ){ // if pass exit the for loops
											  k = layers[*first_hit];
											  i = layers[*second_hit];
											  j = layers[*third_hit]; 
											  first_layer_missing = 1;
										  } 
									  }
								  }// 2nd,3rd,4th case

								  if( all_cut_pass_Ele == 1 ) { PixTrkPassed = true;}
								  if( withoutEM_pass_Ele == 1 ) withoutEM_count_Ele = 1;
								  if( withEM_pass_Ele == 1 ) withEM_count_Ele = 1;

							  }// loop for third layer hit
						  }// loopt for second layer hit
					  }// loop for first layer hit

				  }// 3 out of 4
			  }// 3 out of 4
		  }// 3 out of 4

		  if( PixTrkPassed ){
			  trigger_bit_width_ = trigger_bit_width_| (bit1 << nth_me0_pix_deta);
		  }

		  global_index_width++;


	  }//closest dr cut
	  }// me0 loop
	  trigger_bit_width.push_back(trigger_bit_width_);
//	  trigger_bit_width_iso.push_back(trigger_bit_width_iso_);
	  pix_comb.push_back(pix_comb_);

	  pixtrk_tree->Fill();

	  
   }//evnet loop

   outfile->Write();

   outfile->Close();
}
