#define MakeSWMuon_cxx
#include "MakeSWMuon.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include <TLorentzVector.h>
#include <iostream>
#include <TGraph.h>
#include <TGraphErrors.h>

void MakeSWMuon::Loop()
{

	TFile *file = new TFile("Pix_Pix_median.root","recreate");

    if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   float low_et = 0;
   float high_et = 0;

   int binSize = 90;

   vector<float> x[10][2], median[10][2], x_err[10][2], median_err[10][2];
   vector<float> deta_x[4][2], deta_median[4][2], deta_x_err[4][2], deta_median_err[4][2];
/*
   TCanvas *cc = new TCanvas("cc","",800,700);
   TH1F *dPhidist[140]; 
*/ 
   // for 2D distribution plots
   TH2F* SADphi_dist[10][2];
   TH2F* SADeta_dist[4][2];

   for(int j = 0; j < 10; j++){
	   for(int i = 0; i < 2; i++){
		   TString nth_eg_sw_;
		   nth_eg_sw_.Form("%d", j + 1);
		   TString eta_region;
		   eta_region.Form("%d", i + 1);

		   SADphi_dist[j][i] = new TH2F("","",90,10,100,100,-0.02,0.02);
		   SADphi_dist[j][i]->SetName("SADphi_"+nth_eg_sw_+"_eta_region_"+eta_region);
	   }
   }

   for(int j = 0; j < 4; j++){
	   for(int i = 0; i < 2; i++){
		   TString nth_eg_sw_;
		   nth_eg_sw_.Form("%d", j + 1);
		   TString eta_region;
		   eta_region.Form("%d", i + 1);

		   SADeta_dist[j][i] = new TH2F("","",90,10,100,100,-0.02,0.02);
		   SADeta_dist[j][i]->SetName("SADeta_"+nth_eg_sw_+"_eta_region_"+eta_region);
		   
	   }
   }
   for( int nth = 0; nth < 189 ; nth++){
	   low_et = 10. + nth;
	   high_et = low_et + 1.;
	   //cout << "et: " << low_et << endl;

//	   Bool_t flag = false;
//	   if(low_et == 80 || low_et == 140 || low_et == 180 || low_et == 76 || low_et == 136 || low_et == 176 ) flag = true;

//	   if( !flag ) continue; 

	   cout << "et: " << low_et << endl;

	   vector<float> dPhi_pixV_pixV[10][2];
	   vector<float> dEta_pixV_pixV[4][2];
//	   float abs_pterr;
	   float closest_me0Et;
	   float closest_me0Eta;
	   float closest_me0Phi;
	   float closest_me0_dr;
   
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
//   for (Long64_t jentry=0; jentry<100;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      //
      //find closest me0 object to the gen
      int me0N_ = me0SegPosX->size(); // !!!!!! Use me0SegPosX instead me0 Pt !!!!!!
//      cout<<"me0N = "<<me0N_<<endl;
      if(me0N_ == 0 ) continue; //skip if there is no l1 me objects at all


      // closest me0
      float closest_dr = 9999;
      int closest_me0 = 0;
      int me0Count = 0;

      for( int i = 0; i < me0N_; i++){
	      TVector3 me0SegPos;
	      me0SegPos.SetXYZ(me0SegPosX->at(i), me0SegPosY->at(i), me0SegPosZ->at(i));
	      float me0Phi = me0SegPos.Phi();
	      float me0Eta = me0SegPos.Eta();

	      float dPhi = deltaPhi(genPartPhi->at(0), me0Phi);
	      float current_dr = sqrt( pow(dPhi,2) + pow(genPartEta->at(0) - me0Eta,2) );
	      if( genPartPt->at(0) > high_et || genPartPt->at(0) < low_et ) continue; // !!!!! Use genPartPt instead me0 Pt !!!!!
	      me0Count++;
	      if(current_dr < closest_dr){
		      closest_dr = current_dr;
		      closest_me0 = i;
	      }
      }// end of loop to find the closest me0 the gen

//      cout<<"me0Count = "<<me0Count<<endl;


      if(closest_dr == 9999) continue;
      if(me0Count == 0) continue;

      //cout<<"me0N_ : "<< me0N_ <<endl;
      //cout<<"genPartPt = "<<genPartPt->at(0)<<endl;
      //cout<<"closest_me0 = "<<closest_me0<<endl;
      //cout<<"closest_dr = "<<closest_dr<<endl;
      TVector3 emvector;

//      abs_pterr = fabs(genPartPt->at(0)-propgenElPartPt->at(0))/genPartPt->at(0); // !!!!! Use propgenElPartPt instead closest_me0 Pt !!!!!
      closest_me0Et = genPartPt->at(0); // !!!!! Use propgenElPartPt instead closest_me0 Pt !!!!!
      TVector3 me0SegPos;
      me0SegPos.SetXYZ(me0SegPosX->at(closest_me0), me0SegPosY->at(closest_me0), me0SegPosZ->at(closest_me0));
      closest_me0Eta = me0SegPos.Eta();
      closest_me0Phi = me0SegPos.Phi();
      closest_me0_dr = closest_dr;
      emvector.SetXYZ(me0SegPosX->at(closest_me0),me0SegPosY->at(closest_me0), me0SegPosZ->at(closest_me0));

      // store pixel hits into vectors
      std::vector<TVector3> first_disk_hits;
      std::vector<TVector3> second_disk_hits;
      std::vector<TVector3> third_disk_hits;
      std::vector<TVector3> fourth_disk_hits;
      std::vector<TVector3> fifth_disk_hits;
      std::vector<int> hitted_layers;

      int layers[6] = {}; // initialize as 0 for each event, save number of hits for each pixel layer
      layers[0] = 1; // beam spot

      std::vector<TVector3> PiXTRK_first_hits;
      std::vector<TVector3> PiXTRK_second_hits;
      std::vector<TVector3> PiXTRK_third_hits;
      std::vector<TVector3> PiXTRK_fourth_hits;

      int PiXTRK_layers[5] = {}; // initialize as 0 for each event, save number of hits for each pixel layer
      PiXTRK_layers[0] = 1; // beam spot

      int eta_region = 0;
      if( fabs(closest_me0Eta) < 2.4 && fabs(closest_me0Eta) > 2.0 ) eta_region = 1;
      if( fabs(closest_me0Eta) < 2.8 && fabs(closest_me0Eta) > 2.4 ) eta_region = 2;

      if( fabs(closest_me0Eta) < 2.0 || fabs(closest_me0Eta) > 2.8 ) continue;

      //set roi dphi cut
      float upper_roi = ROI_func(eta_region, closest_me0Eta) + 0.005;
      float lower_roi = ROI_func(eta_region, closest_me0Eta) - 0.005;

      int fpix_size = fRecHitGx->size();

      for( int a = 0; a < fpix_size; a++){
		
	      TVector3 current_hit;
	      current_hit.SetXYZ( fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a) );
	      double Dphi = deltaPhi(current_hit.Phi(), closest_me0Phi);

	      if( Dphi > upper_roi || Dphi < lower_roi ) continue;

	      if( fRecHitDisk->at(a) == 1 ){
		      layers[1]++;
		      first_disk_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a) ));
	      }
	      if( fRecHitDisk->at(a) == 2 ){
		      layers[2]++;
		      second_disk_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a) ));
	      }
	      if( fRecHitDisk->at(a) == 3 ){
		      layers[3]++;
		      third_disk_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a) ));
	      }
	      if( fRecHitDisk->at(a) == 4 ){
		      layers[4]++;
		      fourth_disk_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a) ));
	      }
	      if( fRecHitDisk->at(a) == 5 ){
		      layers[5]++;
		      fifth_disk_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a) ));
	      }
      }

      if( fabs(closest_me0Eta) < 2.4 && fabs(closest_me0Eta) > 2.0){

	      PiXTRK_first_hits = first_disk_hits;
	      PiXTRK_layers[1] = PiXTRK_first_hits.size();

	      PiXTRK_second_hits = second_disk_hits;
	      PiXTRK_layers[2] = PiXTRK_second_hits.size();

	      PiXTRK_third_hits = third_disk_hits;
	      PiXTRK_layers[3] = PiXTRK_third_hits.size();

	      PiXTRK_fourth_hits = fourth_disk_hits;
	      PiXTRK_layers[4] = PiXTRK_fourth_hits.size();
      }

      if( fabs(closest_me0Eta) < 2.8 && fabs(closest_me0Eta) > 2.4){

	      PiXTRK_first_hits = second_disk_hits;
	      PiXTRK_layers[1] = PiXTRK_first_hits.size();

	      PiXTRK_second_hits = third_disk_hits;
	      PiXTRK_layers[2] = PiXTRK_second_hits.size();

	      PiXTRK_third_hits = fourth_disk_hits;
	      PiXTRK_layers[3] = PiXTRK_third_hits.size();

	      PiXTRK_fourth_hits = fifth_disk_hits;
	      PiXTRK_layers[4] = PiXTRK_fourth_hits.size();
      }

      int n_pixels = 0; // initialize as 0 for each event
      hitted_layers.push_back(0); // 0 mean beam spot i.e., (0,0,0)
      for( int i = 1; i < 5; i++){
      	if( PiXTRK_layers[i] != 0 ){
      		hitted_layers.push_back(i); // check if hits on each barrel or disk exists
      		n_pixels++;
      	}
      }

	  if( n_pixels >= 3 ){
		  for( std::vector<int>::iterator first_hit = hitted_layers.begin(); first_hit != hitted_layers.end(); first_hit++){
			  for( std::vector<int>::iterator second_hit = first_hit+1; second_hit != hitted_layers.end(); second_hit++){
				  for( std::vector<int>::iterator third_hit = second_hit+1; third_hit != hitted_layers.end(); third_hit++){

					  for ( int k = 0; k < PiXTRK_layers[*first_hit]; k++){
						  for( int i = 0; i < PiXTRK_layers[*second_hit]; i++){
							  for( int j = 0; j < PiXTRK_layers[*third_hit]; j++){
								  double dPhi = 0;
								  double dEta = 0;

								  if( *first_hit == 0 && *second_hit == 1 && *third_hit == 2 ){
									  dPhi = StandaloneDPhi_BS( PiXTRK_first_hits[i], PiXTRK_second_hits[j]);
									  if( fabs(closest_me0Eta) > 2.0 && fabs(closest_me0Eta) < 2.4  ){ dPhi_pixV_pixV[0][0].push_back(dPhi); SADphi_dist[0][0]->Fill(closest_me0Et,dPhi);}
									  if( fabs(closest_me0Eta) > 2.4 && fabs(closest_me0Eta) < 2.8  ){ dPhi_pixV_pixV[0][1].push_back(dPhi); SADphi_dist[0][1]->Fill(closest_me0Et,dPhi);}
								  }
								  if( *first_hit == 0 && *second_hit == 1 && *third_hit == 3 ){
									  dPhi = StandaloneDPhi_BS( PiXTRK_first_hits[i], PiXTRK_third_hits[j]);
									  if( fabs(closest_me0Eta) > 2.0 && fabs(closest_me0Eta) < 2.4  ){ dPhi_pixV_pixV[1][0].push_back(dPhi); SADphi_dist[1][0]->Fill(closest_me0Et,dPhi);}
									  if( fabs(closest_me0Eta) > 2.4 && fabs(closest_me0Eta) < 2.8  ){ dPhi_pixV_pixV[1][1].push_back(dPhi); SADphi_dist[1][1]->Fill(closest_me0Et,dPhi);}
								  }
								  if( *first_hit == 0 && *second_hit == 1 && *third_hit == 4 ){
									  dPhi = StandaloneDPhi_BS( PiXTRK_first_hits[i], PiXTRK_fourth_hits[j]);
									  if( fabs(closest_me0Eta) > 2.0 && fabs(closest_me0Eta) < 2.4  ){ dPhi_pixV_pixV[2][0].push_back(dPhi); SADphi_dist[2][0]->Fill(closest_me0Et,dPhi);}
									  if( fabs(closest_me0Eta) > 2.4 && fabs(closest_me0Eta) < 2.8  ){ dPhi_pixV_pixV[2][1].push_back(dPhi); SADphi_dist[2][1]->Fill(closest_me0Et,dPhi);}
								  }
								  if( *first_hit == 0 && *second_hit == 2 && *third_hit == 3 ){
									  dPhi = StandaloneDPhi_BS( PiXTRK_second_hits[i], PiXTRK_third_hits[j]);
									  if( fabs(closest_me0Eta) > 2.0 && fabs(closest_me0Eta) < 2.4  ){ dPhi_pixV_pixV[3][0].push_back(dPhi); SADphi_dist[3][0]->Fill(closest_me0Et,dPhi);}
									  if( fabs(closest_me0Eta) > 2.4 && fabs(closest_me0Eta) < 2.8  ){ dPhi_pixV_pixV[3][1].push_back(dPhi); SADphi_dist[3][1]->Fill(closest_me0Et,dPhi);}
								  }
								  if( *first_hit == 0 && *second_hit == 2 && *third_hit == 4 ){
									  dPhi = StandaloneDPhi_BS( PiXTRK_second_hits[i], PiXTRK_fourth_hits[j]);
									  if( fabs(closest_me0Eta) > 2.0 && fabs(closest_me0Eta) < 2.4  ){ dPhi_pixV_pixV[4][0].push_back(dPhi); SADphi_dist[4][0]->Fill(closest_me0Et,dPhi);}
									  if( fabs(closest_me0Eta) > 2.4 && fabs(closest_me0Eta) < 2.8  ){ dPhi_pixV_pixV[4][1].push_back(dPhi); SADphi_dist[4][1]->Fill(closest_me0Et,dPhi);}
								  }
								  if( *first_hit == 0 && *second_hit == 3 && *third_hit == 4 ){
									  dPhi = StandaloneDPhi_BS( PiXTRK_third_hits[i], PiXTRK_fourth_hits[j]);
									  if( fabs(closest_me0Eta) > 2.0 && fabs(closest_me0Eta) < 2.4  ){ dPhi_pixV_pixV[5][0].push_back(dPhi); SADphi_dist[5][0]->Fill(closest_me0Et,dPhi);}
									  if( fabs(closest_me0Eta) > 2.4 && fabs(closest_me0Eta) < 2.8  ){ dPhi_pixV_pixV[5][1].push_back(dPhi); SADphi_dist[5][1]->Fill(closest_me0Et,dPhi);}
								  }

								  if( *first_hit == 1 && *second_hit == 2 && *third_hit == 3 ){
									  dPhi = StandaloneDPhi( PiXTRK_first_hits[k], PiXTRK_second_hits[i], PiXTRK_third_hits[j]);
									  dEta = StandaloneDEta( PiXTRK_first_hits[k], PiXTRK_second_hits[i], PiXTRK_third_hits[j]);
									  if( fabs(closest_me0Eta) > 2.0 && fabs(closest_me0Eta) < 2.4 ){
										  dPhi_pixV_pixV[6][0].push_back(dPhi); dEta_pixV_pixV[0][0].push_back(dEta);
										  SADphi_dist[6][0]->Fill(closest_me0Et,dPhi); SADeta_dist[0][0]->Fill(closest_me0Et,dEta);
									  }
									  if( fabs(closest_me0Eta) > 2.4 && fabs(closest_me0Eta) < 2.8 ){
										  dPhi_pixV_pixV[6][1].push_back(dPhi); dEta_pixV_pixV[0][1].push_back(dEta);
										  SADphi_dist[6][1]->Fill(closest_me0Et,dPhi); SADeta_dist[0][1]->Fill(closest_me0Et,dEta);
									  }
								  }
								  if( *first_hit == 1 && *second_hit == 2 && *third_hit == 4 ){
									  dPhi = StandaloneDPhi( PiXTRK_first_hits[k], PiXTRK_second_hits[i], PiXTRK_fourth_hits[j]);
									  dEta = StandaloneDEta( PiXTRK_first_hits[k], PiXTRK_second_hits[i], PiXTRK_fourth_hits[j]);
									  if( fabs(closest_me0Eta) > 2.0 && fabs(closest_me0Eta) < 2.4 ){
										  dPhi_pixV_pixV[7][0].push_back(dPhi); dEta_pixV_pixV[1][0].push_back(dEta);
										  SADphi_dist[7][0]->Fill(closest_me0Et,dPhi); SADeta_dist[1][0]->Fill(closest_me0Et,dEta);
									  }
									  if( fabs(closest_me0Eta) > 2.4 && fabs(closest_me0Eta) < 2.8 ){
										  dPhi_pixV_pixV[7][1].push_back(dPhi); dEta_pixV_pixV[1][1].push_back(dEta);
										  SADphi_dist[7][1]->Fill(closest_me0Et,dPhi); SADeta_dist[1][1]->Fill(closest_me0Et,dEta);
									  }
								  }
								  if( *first_hit == 1 && *second_hit == 3 && *third_hit == 4 ){
									  dPhi = StandaloneDPhi( PiXTRK_first_hits[k], PiXTRK_third_hits[i], PiXTRK_fourth_hits[j]);
									  dEta = StandaloneDEta( PiXTRK_first_hits[k], PiXTRK_third_hits[i], PiXTRK_fourth_hits[j]);
									  if( fabs(closest_me0Eta) > 2.0 && fabs(closest_me0Eta) < 2.4 ){
										  dPhi_pixV_pixV[8][0].push_back(dPhi); dEta_pixV_pixV[2][0].push_back(dEta);
										  SADphi_dist[8][0]->Fill(closest_me0Et,dPhi); SADeta_dist[2][0]->Fill(closest_me0Et,dEta);
									  }
									  if( fabs(closest_me0Eta) > 2.4 && fabs(closest_me0Eta) < 2.8 ){
										  dPhi_pixV_pixV[8][1].push_back(dPhi); dEta_pixV_pixV[2][1].push_back(dEta);
										  SADphi_dist[8][1]->Fill(closest_me0Et,dPhi); SADeta_dist[2][1]->Fill(closest_me0Et,dEta);
									  }
								  }
								  if( *first_hit == 2 && *second_hit == 3 && *third_hit == 4 ){
									  dPhi = StandaloneDPhi( PiXTRK_second_hits[k], PiXTRK_third_hits[i], PiXTRK_fourth_hits[j]);
									  dEta = StandaloneDEta( PiXTRK_second_hits[k], PiXTRK_third_hits[i], PiXTRK_fourth_hits[j]);
									  if( fabs(closest_me0Eta) > 2.0 && fabs(closest_me0Eta) < 2.4 ){
										  dPhi_pixV_pixV[9][0].push_back(dPhi); dEta_pixV_pixV[3][0].push_back(dEta);
										  SADphi_dist[9][0]->Fill(closest_me0Et,dPhi); SADeta_dist[3][0]->Fill(closest_me0Et,dEta);
									  }
									  if( fabs(closest_me0Eta) > 2.4 && fabs(closest_me0Eta) < 2.8 ){
										  dPhi_pixV_pixV[9][1].push_back(dPhi); dEta_pixV_pixV[3][1].push_back(dEta);
										  SADphi_dist[9][1]->Fill(closest_me0Et,dPhi); SADeta_dist[3][1]->Fill(closest_me0Et,dEta);
									  }
								  }

							  }// jth hit in third pixel disk
						  }// ith hit in second pixel disk
					  }// kth hit in first pixel disk
				  }
			  }
		  }
	  }

   }// event loop

   for(int nth_me0_sw = 0; nth_me0_sw < 10; nth_me0_sw++){
	   for(int i = 0; i < 2; i++){
		   std::sort (dPhi_pixV_pixV[nth_me0_sw][i].begin(), dPhi_pixV_pixV[nth_me0_sw][i].end());

		   if( dPhi_pixV_pixV[nth_me0_sw][i].size() != 0 ){
			   x[nth_me0_sw][i].push_back(low_et + 0.5);
			   median[nth_me0_sw][i].push_back(getMedian(dPhi_pixV_pixV[nth_me0_sw][i]));
			   x_err[nth_me0_sw][i].push_back(0.);
	//		   x_err[nth_me0_sw][i].push_back(0.);
	//		   median_err[nth_me0_sw][i].push_back(getMedianErr(dPhi_pixV_pixV[nth_me0_sw][i]));
			   median_err[nth_me0_sw][i].push_back(0.);
		   }
//// For drawing each pt range ////
/*******************************************8
		   cout<<"count = "<<count<<endl;

		   Char_t histname[20];
		   sprintf(histname, "hist_%dth", count);
		   dPhidist[count] = new TH1F(histname,"dPhidist",1000,-0.5,0.5);    
		   for( int a = 0; a < dPhi_l1me0_pixV[nth_me0_sw][i].size(); a++ ){
		   dPhidist[count]->Fill(dPhi_l1me0_pixV[nth_me0_sw][i].at(a));
		   }
		   TString nth_me0_sw_;
		   nth_me0_sw_.Form("%d", nth_me0_sw + 1 );
		   TString eta_region_;
		   eta_region_.Form("%d", i + 1);
		   TString pt_range_;		   
		   pt_range_.Form("%d", (int)low_et );
		   
		   // Draw andp save histograms
		   dPhidist[count]->SetTitle("test");
		   dPhidist[count]->Draw();
		   cc->SaveAs("dPhi_result/dPhi_"+nth_me0_sw_+"th_Pixel_"+eta_region_+"_region_"+pt_range_+"GeV.png");
		   cc->Clear();

		   count++;
*******************************************/
	   }
   }
   for(int nth_me0_sw = 0; nth_me0_sw < 4; nth_me0_sw++){
	   for(int i = 0; i < 2; i++){
		   std::sort (dEta_pixV_pixV[nth_me0_sw][i].begin(), dEta_pixV_pixV[nth_me0_sw][i].end());

		   if( dEta_pixV_pixV[nth_me0_sw][i].size() != 0 ){
			   deta_x[nth_me0_sw][i].push_back(low_et + 0.5);
			   deta_median[nth_me0_sw][i].push_back(getMedian(dEta_pixV_pixV[nth_me0_sw][i]));
		//	   deta_x_err[nth_me0_sw][i].push_back(0.);
			   deta_x_err[nth_me0_sw][i].push_back(0.);
		//	   deta_median_err[nth_me0_sw][i].push_back(getMedianErr(dEta_pixV_pixV[nth_me0_sw][i]));
			   deta_median_err[nth_me0_sw][i].push_back(0.);
		   }
	   }
   }
   
   }// et loop

   TGraphErrors* dphi_median_gr[10][2];
   for( int j = 0; j < 10; j++){
	   for( int i = 0; i < 2; i++){
		   int dphi_point_size = x[j][i].size();
		   dphi_median_gr[j][i] = new TGraphErrors(dphi_point_size, &x[j][i][0], &median[j][i][0], &x_err[j][i][0], &median_err[j][i][0]);

//		   TCanvas *c1 = new TCanvas();
//		   dphi_median_gr[j][i]->SetFillStyle(3002);
//		   dphi_median_gr[j][i]->SetFillColorAlpha(kBlue,0.2);

		   dphi_median_gr[j][i]->SetMarkerStyle(24);
		   dphi_median_gr[j][i]->SetMarkerSize(0.5);
		   TString nth_me0_sw_;
		   nth_me0_sw_.Form("%d", j + 1);
		   TString eta_region;
		   eta_region.Form("%d", i + 1);

		   dphi_median_gr[j][i]->SetName( "pix_pix_dphi_SW_"+nth_me0_sw_+"_eta_region_"+eta_region+"_median"); 

//		   dphi_median_gr[j][i]->Draw("apE3");
//		   c1->SaveAs("Pix_Pix_dPhi_SW_"+nth_me0_sw_+"_eta_region_"+eta_region+".png");

		   dphi_median_gr[j][i]->Write();
	   }
   }
   TGraphErrors* deta_median_gr[4][2];
   for( int j = 0; j < 4; j++){
	   for( int i = 0; i < 2; i++){
		   int deta_point_size = deta_x[j][i].size();
		   deta_median_gr[j][i] = new TGraphErrors(deta_point_size, &deta_x[j][i][0], &deta_median[j][i][0], &deta_x_err[j][i][0], &deta_median_err[j][i][0]);

//		   TCanvas *c2 = new TCanvas();
//		   c2->GetGridx();
//		   c2->GetGridy();

//		   deta_median_gr[j][i]->SetMarkerColor(kOrange);
//		   deta_median_gr[j][i]->SetFillStyle(3002);
//		   deta_median_gr[j][i]->SetFillColorAlpha(kBlue,0.2);

		   deta_median_gr[j][i]->SetMarkerStyle(24);
		   deta_median_gr[j][i]->SetMarkerSize(0.5);
		   TString nth_me0_sw_;
		   nth_me0_sw_.Form("%d", j + 1);
		   TString eta_region;
		   eta_region.Form("%d", i + 1);

		   deta_median_gr[j][i]->SetName( "pix_pix_deta_SW_"+nth_me0_sw_+"_eta_region_"+eta_region+"_median");

		   deta_median_gr[j][i]->GetYaxis()->SetRangeUser(-0.015,0.015);
//		   deta_median_gr[j][i]->Draw("apE3");
		   c2->SaveAs("Pix_Pix_dEta_sw_"+nth_me0_sw_+"_eta_region_"+eta_region+".png");
		   deta_median_gr[j][i]->Write();
	   }
   }

   // save all 2D plots
   for(int j = 0; j < 10; j++){
	   for(int i = 0; i < 2; i++){
		   SADphi_dist[j][i]->Write();
	   }
   }
   for(int j = 0; j < 4; j++){
	   for(int i = 0; i < 2; i++){
		   SADeta_dist[j][i]->Write();
	   
   }

   file->Write();
}
}
