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
//   In a ROOT session, you can do:
//      root> .L MakeSWMuon.C
//      root> MakeSWMuon t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
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

	TFile *file = new TFile("roi_median.root","recreate");


    if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   float low_et = 0;
   float high_et = 0;

   int binSize = 90;
   int count = 0;

   vector<float> x[6][2], median[6][2], x_err[6][2], median_err[6][2];
   vector<float> deta_x[6][2], deta_median[6][2], deta_x_err[6][2], deta_median_err[6][2];

   // For save each dPhi //
//   TCanvas *cc = new TCanvas("cc","",800,700);
//   TH1F *dPhidist[140]; 
   
   for( int nth = 0; nth < 189; nth++){
	   low_et = 10. + nth;
	   high_et = low_et + 1.;
	   //cout << "et: " << low_et << endl;

	  //Bool_t flag = false;
	  //if(low_et == 80 || low_et == 140 || low_et == 180 || low_et == 76 || low_et == 136 || low_et == 176 ) flag = true;

	  //if( !flag ) continue; 

	   cout << "et: " << low_et << endl;
	   vector<float> dPhi_l1me0_pixV[6][2];
	   vector<float> dEta_l1me0_pixV[6][2];
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

	      float dPhi = deltaPhi(propgenElPartPhi->at(0), me0Phi);
	      float current_dr = sqrt( pow(dPhi,2) + pow(propgenElPartEta->at(0) - me0Eta,2) );
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
      closest_me0Et = propgenElPartPt->at(0); // !!!!! Use propgenElPartPt instead closest_me0 Pt !!!!!
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
      float upper_roi = ROI_func(eta_region, closest_me0Eta) + 0.055;
      float lower_roi = ROI_func(eta_region, closest_me0Eta) - 0.055;

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
//      if( n_pixels < 3 ) continue;

      if( n_pixels >= 3 ){
      	for( std::vector<int>::iterator first_hit = hitted_layers.begin()+1; first_hit != hitted_layers.end(); first_hit++){
		for( std::vector<int>::iterator second_hit = first_hit+1; second_hit != hitted_layers.end(); second_hit++){
		
		for ( int k = 0; k < PiXTRK_layers[*first_hit]; k++){
			for( int i = 0; i < PiXTRK_layers[*second_hit]; i++){
			double dPhi = 0;
			double dEta = 0;

			TVector3 first_layer_;
			TVector3 second_layer_;

			if( *first_hit == 1 ) first_layer_ = PiXTRK_first_hits[k];
			if( *first_hit == 2 ) first_layer_ = PiXTRK_second_hits[k];
			if( *first_hit == 3 ) first_layer_ = PiXTRK_third_hits[k];

			if( *second_hit == 2 ) second_layer_ = PiXTRK_second_hits[i];
			if( *second_hit == 3 ) second_layer_ = PiXTRK_third_hits[i];
			if( *second_hit == 4 ) second_layer_ = PiXTRK_fourth_hits[i];

			TVector3 pixelVector = second_layer_ - first_layer_;
			TVector3 EM_pixelVector = emvector - second_layer_;

			dPhi = deltaPhi(EM_pixelVector.Phi(), pixelVector.Phi());
			dEta = EM_pixelVector.Eta() - pixelVector.Eta();

//			cout<<"closest_me0Eta = "<<closest_me0Eta<<endl;

			if( *first_hit == 1 && *second_hit == 2 ){
				if( fabs(closest_me0Eta) > 2.0 && fabs(closest_me0Eta) < 2.4  ){ dPhi_l1me0_pixV[0][0].push_back(dPhi); dEta_l1me0_pixV[0][0].push_back(dEta);}
				if( fabs(closest_me0Eta) > 2.4 && fabs(closest_me0Eta) < 2.8  ){ dPhi_l1me0_pixV[0][1].push_back(dPhi); dEta_l1me0_pixV[0][1].push_back(dEta);}
			}
			if( *first_hit == 1 && *second_hit == 3 ){
				if( fabs(closest_me0Eta) > 2.0 && fabs(closest_me0Eta) < 2.4  ){ dPhi_l1me0_pixV[1][0].push_back(dPhi); dEta_l1me0_pixV[1][0].push_back(dEta);}
				if( fabs(closest_me0Eta) > 2.4 && fabs(closest_me0Eta) < 2.8  ){ dPhi_l1me0_pixV[1][1].push_back(dPhi); dEta_l1me0_pixV[1][1].push_back(dEta);}
			}
			if( *first_hit == 1 && *second_hit == 4 ){
				if( fabs(closest_me0Eta) > 2.0 && fabs(closest_me0Eta) < 2.4  ){ dPhi_l1me0_pixV[2][0].push_back(dPhi); dEta_l1me0_pixV[2][0].push_back(dEta);}
				if( fabs(closest_me0Eta) > 2.4 && fabs(closest_me0Eta) < 2.8  ){ dPhi_l1me0_pixV[2][1].push_back(dPhi); dEta_l1me0_pixV[2][1].push_back(dEta);}
			}
			if( *first_hit == 2 && *second_hit == 3 ){
				if( fabs(closest_me0Eta) > 2.0 && fabs(closest_me0Eta) < 2.4  ){ dPhi_l1me0_pixV[3][0].push_back(dPhi); dEta_l1me0_pixV[3][0].push_back(dEta);}
				if( fabs(closest_me0Eta) > 2.4 && fabs(closest_me0Eta) < 2.8  ){ dPhi_l1me0_pixV[3][1].push_back(dPhi); dEta_l1me0_pixV[3][1].push_back(dEta);}
			}
			if( *first_hit == 2 && *second_hit == 4 ){
				if( fabs(closest_me0Eta) > 2.0 && fabs(closest_me0Eta) < 2.4  ){ dPhi_l1me0_pixV[4][0].push_back(dPhi); dEta_l1me0_pixV[4][0].push_back(dEta);}
				if( fabs(closest_me0Eta) > 2.4 && fabs(closest_me0Eta) < 2.8  ){ dPhi_l1me0_pixV[4][1].push_back(dPhi); dEta_l1me0_pixV[4][1].push_back(dEta);}
			}
			if( *first_hit == 3 && *second_hit == 4 ){
				if( fabs(closest_me0Eta) > 2.0 && fabs(closest_me0Eta) < 2.4  ){ dPhi_l1me0_pixV[5][0].push_back(dPhi); dEta_l1me0_pixV[5][0].push_back(dEta);}
				if( fabs(closest_me0Eta) > 2.4 && fabs(closest_me0Eta) < 2.8  ){ dPhi_l1me0_pixV[5][1].push_back(dPhi); dEta_l1me0_pixV[5][1].push_back(dEta);}
			}
		}
	}
      }
	}
      }

   }// event loop

   for(int nth_me0_sw = 0; nth_me0_sw < 6; nth_me0_sw++){
	   for(int i = 0; i < 2; i++){
		   std::sort (dPhi_l1me0_pixV[nth_me0_sw][i].begin(), dPhi_l1me0_pixV[nth_me0_sw][i].end());
		   std::sort (dEta_l1me0_pixV[nth_me0_sw][i].begin(), dEta_l1me0_pixV[nth_me0_sw][i].end());

//		   cout<<"dPhi_l1me0_pixV["<<nth_me0_sw<<"]["<<i<<"] = "<<dPhi_l1me0_pixV[nth_me0_sw][i].size()<<endl;

		   if( dPhi_l1me0_pixV[nth_me0_sw][i].size() != 0 ){
			   x[nth_me0_sw][i].push_back(low_et);
			   median[nth_me0_sw][i].push_back(getMedian(dPhi_l1me0_pixV[nth_me0_sw][i]));
			   x_err[nth_me0_sw][i].push_back(0.);
			   median_err[nth_me0_sw][i].push_back(getMedianErr(dPhi_l1me0_pixV[nth_me0_sw][i]));

			   deta_x[nth_me0_sw][i].push_back(low_et);
			   deta_median[nth_me0_sw][i].push_back(getMedian(dEta_l1me0_pixV[nth_me0_sw][i]));
			   deta_x_err[nth_me0_sw][i].push_back(0.);
			   deta_median_err[nth_me0_sw][i].push_back(getMedianErr(dEta_l1me0_pixV[nth_me0_sw][i]));
		   }


//		   cout<<"count = "<<count<<endl;
/************ Save each dPhi distribution ************/
/*****************************************************
		   Char_t histname[20];
		   sprintf(histname, "hist_%dth", count);
		   dPhidist[count] = new TH1F(histname,"dPhidist",1000,-0.1,0.1);    
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
		   dPhidist[count]->SetTitle("dPhi_"+nth_me0_sw_+"th_pixel_"+eta_region_+"_region_"+pt_range_+"GeV");
		   dPhidist[count]->Draw();
		   cc->SaveAs("dPhi_result/dPhi_"+nth_me0_sw_+"th_Pixel_"+eta_region_+"_region_"+pt_range_+"GeV.png");
		   cc->Clear();

		   count++;
*****************************************************/
	   }
   }

   
   
   }// et loop

   TGraphErrors* dphi_median_gr[6][2];
   TGraphErrors* deta_median_gr[6][2];
   for( int j = 0; j < 6; j++){
	   for( int i = 0; i < 2; i++){
		   int dphi_point_size = x[j][i].size();
		   int deta_point_size = deta_x[j][i].size();
		   dphi_median_gr[j][i] = new TGraphErrors(dphi_point_size, &x[j][i][0], &median[j][i][0], &x_err[j][i][0], &median_err[j][i][0]);
		   deta_median_gr[j][i] = new TGraphErrors(deta_point_size, &deta_x[j][i][0], &deta_median[j][i][0], &deta_x_err[j][i][0], &deta_median_err[j][i][0]);

		   dphi_median_gr[j][i]->SetMarkerStyle(24);
		   deta_median_gr[j][i]->SetMarkerSize(0.3);
		   TString nth_me0_sw_;
		   nth_me0_sw_.Form("%d", j + 1);
		   TString eta_region;
		   eta_region.Form("%d", i + 1);

		   dphi_median_gr[j][i]->SetTitle(nth_me0_sw_+"th_Pixel_"+eta_region+"_median");
		   dphi_median_gr[j][i]->SetName( "ME0_Pix_dphi_SW_"+nth_me0_sw_+"_eta_region_"+eta_region+"_median"); 

		   deta_median_gr[j][i]->SetMarkerStyle(24);
		   deta_median_gr[j][i]->SetMarkerSize(0.3);
		   deta_median_gr[j][i]->SetTitle(nth_me0_sw_+"th_Pixel_"+eta_region+"_median");
		   deta_median_gr[j][i]->SetName( "ME0_Pix_deta_SW_"+nth_me0_sw_+"_eta_region_"+eta_region+"_median");

		   dphi_median_gr[j][i]->Write();
		   deta_median_gr[j][i]->Write();
	   }
   }

   file->Write();
}
