#define PID_cxx
#include "PID.h"
#include <TH2.h>
#include <TH1.h>
#include "TF1.h"
#include <TStyle.h>
#include <TCanvas.h>
#include "TClonesArray.h"


#include<fstream>
#include<stdlib.h>
#include "TString.h"
#include <iostream>
#include "Riostream.h"
#include <string>
#include <time.h>
#include <sstream>
#include "TLine.h"

/*
This function/class is supposed to improve the PID before and after the target

We start with the empty run, as here should be almost no difference in after and before the target


*/

void PID::Loop(short iteration, short var){
  //create all histos for AOQ: 0,1,2,3
  //in dependending on
  /*
    beta
    angle a and b (x,y)(two times: first and second detector)
    x- and y- position (two times: first and second detector)
  */


//   In a ROOT session, you can do:
//      Root > .L PID.C
//      Root > PID t
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
   Histo(var);
   Read();
   Long64_t nbytes = 0, nb = 0;


   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      //global cut
      if(globalcut())continue;
      
      for(short irips = 0; irips < kMaxBigRIPSRIPS; irips++){//getting the variables

	if(BigRIPSRIPS_upstream_fpl[irips] == 3){
	  fin_var[0] = BigRIPSRIPS_angle[irips];
	  //	  angle_b_F3 = BigRIPSRIPS_angle[irips];
	}else if(BigRIPSRIPS_upstream_fpl[irips] == 5){
	  fin_var[1] = BigRIPSRIPS_angle[irips];
	}else if(BigRIPSRIPS_upstream_fpl[irips] == 8){
	  fout_var[0][0] = BigRIPSRIPS_angle[irips];
	  fvar[0][0] = BigRIPSRIPS_angle[irips];
	}

	/*if(BigRIPSRIPS_downstream_fpl[irips] == 9){
	  fout_var[0][0][0] = BigRIPSRIPS_angle[irips];
	  }else*/
	if(BigRIPSRIPS_downstream_fpl[irips] == 10){
	  fout_var[0][1] = BigRIPSRIPS_angle[irips];
	  fout_var[1][0] = BigRIPSRIPS_angle[irips];
	  fout_var[2][0] = BigRIPSRIPS_angle[irips];
	  fvar[1][0] = BigRIPSRIPS_angle[irips];
	}else if(BigRIPSRIPS_downstream_fpl[irips] == 11){
	  fout_var[1][1] = BigRIPSRIPS_angle[irips];
	  fout_var[2][1] = BigRIPSRIPS_angle[irips];
	  fvar[2][0] = BigRIPSRIPS_angle[irips];
	}else continue;

      }//for(RIPS)

      for(short ippac = 0; ippac < kMaxBigRIPSPPAC; ippac++){//getting the variables

	if(BigRIPSPPAC_fpl[ippac] == 3){
	  fin_var[2] = BigRIPSPPAC_fX[ippac];
	  fin_var[4] = BigRIPSPPAC_fY[ippac];
	}else if(BigRIPSPPAC_fpl[ippac] == 5){
	  fin_var[3] = BigRIPSPPAC_fX[ippac];
	  fin_var[5] = BigRIPSPPAC_fY[ippac];
	}else if(BigRIPSPPAC_fpl[ippac] == 8){
	  fout_var[0][2] = BigRIPSPPAC_fX[ippac];
	  fout_var[0][4] = BigRIPSPPAC_fY[ippac];
	  fvar[0][1] = BigRIPSPPAC_fX[ippac];
	  fvar[0][2] = BigRIPSPPAC_fY[ippac];
	}
	
	/*if(BigPPACPPAC_downstream_fpl[ippac] == 9){
	  pos_x_F9 = BigRIPSPPAC_fX[ippac];
	  pos_y_F9 = BigRIPSPPAC_fY[ippac];
	  }else*/
	if(BigRIPSPPAC_fpl[ippac] == 10){
	  fout_var[0][3] = BigRIPSPPAC_fX[ippac];
	  fout_var[0][5] = BigRIPSPPAC_fY[ippac];
	  fout_var[1][2] = BigRIPSPPAC_fX[ippac];
	  fout_var[1][4] = BigRIPSPPAC_fY[ippac];

	  fvar[1][1] = BigRIPSPPAC_fX[ippac];
	  fvar[1][2] = BigRIPSPPAC_fY[ippac];
	}else if(BigRIPSPPAC_fpl[ippac] == 11){
	  fout_var[1][3] = BigRIPSPPAC_fX[ippac];
	  fout_var[1][5] = BigRIPSPPAC_fY[ippac];

	  fvar[2][1] = BigRIPSPPAC_fX[ippac];
	  fvar[2][2] = BigRIPSPPAC_fY[ippac];
	}else continue;

      }//for(PPAC)
      if(var<6){
	cout<< "Hi5" <<endl;
	//fh2_aoq_0->Fill(fin_var[var],Correction(0,var));
	cout<< "Hi.55" <<endl;
	if( ! outcut()){
	  /*
	  fh2_aoq_1->Fill(fout_var[0][var],Correction(1,var));
	  fh2_aoq_2->Fill(fout_var[1][var],Correction(2,var));
	  */
	}
      }else if(var == 6){
	fh2_aoq_0->Fill(BETA[0],Correction(0,var));
	if( ! outcut()){
	  fh2_aoq_1->Fill(BETA[1],Correction(1,var));
	  fh2_aoq_2->Fill(BETA[2],Correction(2,var));
	}
      }
      if(var < 9)
	fh2_aoq_3->Fill(fvar[var/3][var%3],Correction_3(3,var));




      /*
      //_____incoming____________________________
     
      fh2_aoq_0_angle_F3->Fill(fin_var[0],Correction(0,100));
      fh2_aoq_0_angle_F5->Fill(fin_var[1],Correction(0,100));
      fh2_aoq_0_xpos_F3->Fill(fin_var[2],Correction(0,200));
      fh2_aoq_0_xpos_F5->Fill(fin_var[3],Correction(0,300));
      fh2_aoq_0_ypos_F3->Fill(fin_var[4],Correction(0,400));
      fh2_aoq_0_ypos_F5->Fill(fin_var[5],Correction(0,500));
      fh2_aoq_0_beta->Fill(BETA[0],Correction(0,600));

      //_____outgoing___________________________
      if( outcut() )continue;
      //AOQ[1]
      fh2_aoq_1_angle_F8->Fill(fout_var[0][0],Correction(1,100));
      fh2_aoq_1_angle_F10->Fill(fout_var[0][1],Correction(1,100));
      fh2_aoq_1_xpos_F8->Fill(fout_var[0][2],Correction(1,200));
      fh2_aoq_1_xpos_F10->Fill(fout_var[0][3],Correction(1,300));
      fh2_aoq_1_ypos_F8->Fill(fout_var[0][4],Correction(1,400));
      fh2_aoq_1_ypos_F10->Fill(fout_var[0][5],Correction(1,500));
      fh2_aoq_1_beta->Fill(BETA[1],Correction(1,600));

      //aoq[2]
      fh2_aoq_2_angle_F10->Fill(fout_var[1][0],Correction(2,100));
      fh2_aoq_2_angle_F11->Fill(fout_var[1][1],Correction(2,100));
      fh2_aoq_2_xpos_F10->Fill(fout_var[1][2],Correction(2,100));
      fh2_aoq_2_xpos_F11->Fill(fout_var[1][3],Correction(2,100));
      fh2_aoq_2_ypos_F10->Fill(fout_var[1][4],Correction(2,100));
      fh2_aoq_2_ypos_F11->Fill(fout_var[1][5],Correction(2,100));
      fh2_aoq_2_beta->Fill(BETA[2],Correction(2,100));
      //aoq[3]
      
      fh2_aoq_3_angle_F8 ->Fill(fvar[0][0],Correction_3(3,100));
      fh2_aoq_3_angle_F10->Fill(fvar[1][0],Correction_3(3,100));
      fh2_aoq_3_angle_F11->Fill(fvar[2][0],Correction_3(3,100));

      fh2_aoq_3_xpos_F8 ->Fill(fvar[0][1],Correction_3(3,100));
      fh2_aoq_3_xpos_F10->Fill(fvar[1][1],Correction_3(3,100));
      fh2_aoq_3_xpos_F11->Fill(fvar[1][2],Correction_3(3,100));

      fh2_aoq_3_ypos_F8 ->Fill(fvar[0][2],Correction_3(3,100));
      fh2_aoq_3_ypos_F10->Fill(fvar[1][2],Correction_3(3,100));
      fh2_aoq_3_ypos_F11->Fill(fvar[2][2],Correction_3(3,100));

      fh2_aoq_3_beta->Fill(BETA[3],Correction_3(3,100));
      
      // if (Cut(ientry) < 0) continue;
      */
   }//for(events)
   
   TCanvas* can = new TCanvas("can","canvas");

   can->Divide(2,2);
   can->cd(1);
   //   fh2_aoq_0->Draw("colz");

   can->cd(2);
   //fh2_aoq_1->Draw("colz");

   can->cd(3);
   //fh2_aoq_2->Draw("colz");

   can->cd(4);
   fh2_aoq_3->Draw("colz");
      
}//Loop

/** _________________________________________________________ **/
void PID::Fitting(TH2F* h2, short aoq, short var){

  TObjArray* dummy = new TObjArray();//needed for FitSlices to delete the histograms
  dummy->SetOwner(kTRUE);

  TF1* gausfunc = new TF1("gausfunc","gaus",2.5,2.7);//(0)+pol1(3)
  h2->FitSlicesY(gausfunc,0,-1,0,"QNRG4",dummy);
  TH1F* h1 = (TH1F*) dummy->FindObject(Form("fh2_aoq_%i_1",aoq));
  TF1* Linfunc = new TF1("Linfunc","[0]*(x+[1])+2.64",-10,10);
  Linfunc->SetParameters(0,0);
  Linfunc->SetParLimits(0,-1e-3,1e-3);
  Linfunc->SetParLimits(1,-1e-3,1e-3);
  h1->Fit(Linfunc);

  if(aoq!=3)
    fcorr_Lin[aoq][var] = Linfunc->GetParameter(0);
  else{
    short plane = var/3;
    short var_ = var%3;
    fcorr_Lin_3[plane][var_] =  Linfunc->GetParameter(0);
  }
  //  fcorr_Lin[aoq][var] = Linfunc->GetParameter(2);


  //dummy->Clear();

}//incoming

/** _________________________________________________________ **/
void PID::Outgoing(TH2F* h2_angle_a, TH2F* h2_angle_b,// TH2F* h2_angle_a_, TH2F* h2_angle_b_,
		   TH2F* h2_pos_x, TH2F* h2_pos_y, TH2F* h2_pos_x_, TH2F* h2_pos_y_,
		   TH2F* h2_beta){




}//Outgoing_1

/** _________________________________________________________ **/
Bool_t PID::globalcut(){

  if( TMath::Abs(ZET[0]-50)>0.3 )return kTRUE;
  if( TMath::Abs(AOQ[0]-2.64)>0.005)return kTRUE;

  return kFALSE;
}

Bool_t PID::outcut(){  

  if( TMath::Abs(delta[2] - delta[4] )> 0.2 ) return kTRUE;
  if( TMath::Abs(ZET[2]-45.3)>0.5 )return kTRUE;
  if( TMath::Abs(AOQ[2]-2.64)>0.005)return kTRUE;

  return kFALSE;
}

void PID::Histo(short var){

  if(var != 6){

    fh2_aoq_0 = new TH2F(Form("fh2_aoq_0_%i",var),"aoq vs ",200,-10,10,200,2.6,2.7);
    fh2_aoq_1 = new TH2F("fh2_aoq_111","aoq vs ",200,-10,10,200,2.6,2.7);
    fh2_aoq_2 = new TH2F("fh2_aoq_222","aoq vs ",200,-10,10,200,2.6,2.7);
    
    //fh2_aoq_4 = new TH2F("fh2_aoq_4","aoq vs ",200,-10,10,200,2.6,2.7);
  }else {
    fh2_aoq_0 = new TH2F("fh2_aoq_000","aoq vs ",200,0.5,0.6,200,2.6,2.7);
    fh2_aoq_1 = new TH2F("fh2_aoq_111","aoq vs ",200,0.5,0.6,200,2.6,2.7);
    fh2_aoq_2 = new TH2F("fh2_aoq_222","aoq vs ",200,0.5,0.6,200,2.6,2.7);

    //fh2_aoq_4 = new TH2F("fh2_aoq_4","aoq vs ",200,0.5,0.6,200,2.6,2.7);
  }
  if(var != 9)
    fh2_aoq_3 = new TH2F("fh2_aoq_333","aoq vs ",200,-10,10,200,2.6,2.7);
  else 
    fh2_aoq_3 = new TH2F("fh2_aoq_333","aoq vs ",200,0.5,0.6,200,2.6,2.7);

  fh2_aoq_0_angle_F3 = new TH2F("fh2_aoq_0_angle_F3","aoq vs angle F3",200,-10,10,200,2.6,2.7);
  fh2_aoq_0_angle_F5 = new TH2F("fh2_aoq_0_angle_F5","aoq vs angle F5",200,-10,10,200,2.6,2.7);
  fh2_aoq_0_xpos_F3 = new TH2F("fh2_aoq_0_xpos_F3","aoq vs xpos F3",200,-10,10,200,2.6,2.7); 
  fh2_aoq_0_xpos_F5 = new TH2F("fh2_aoq_0_xpos_F5","aoq vs xpos F5",200,-10,10,200,2.6,2.7); 
  fh2_aoq_0_ypos_F3 = new TH2F("fh2_aoq_0_ypos_F3","aoq vs ypos F3",200,-10,10,200,2.6,2.7); 
  fh2_aoq_0_ypos_F5 = new TH2F("fh2_aoq_0_ypos_F5","aoq vs ypos F5",200,-10,10,200,2.6,2.7); 
  fh2_aoq_0_beta = new TH2F("fh2_aoq_0_beta","aoq vs beta",200,0.5,0.6,200,2.6,2.7);	   
                   
  fh2_aoq_1_angle_F8 = new TH2F("fh2_aoq_1_angle_F8","aoq vs angle F8",200,-10,10,200,2.6,2.7);
  fh2_aoq_1_angle_F10 = new TH2F("fh2_aoq_1_angle_F10","aoq vs angle F10",200,-10,10,200,2.6,2.7);
  fh2_aoq_1_xpos_F8 = new TH2F("fh2_aoq_1_xpos_F8","aoq vs xpos F8",200,-10,10,200,2.6,2.7); 
  fh2_aoq_1_xpos_F10 = new TH2F("fh2_aoq_1_xpos_F10","aoq vs xpos F10",200,-10,10,200,2.6,2.7);
  fh2_aoq_1_ypos_F8 = new TH2F("fh2_aoq_1_ypos_F8","aoq vs ypos F8",200,-10,10,200,2.6,2.7); 
  fh2_aoq_1_ypos_F10 = new TH2F("fh2_aoq_1_ypos_F10","aoq vs ypos F10 ",200,-10,10,200,2.6,2.7);
  fh2_aoq_1_beta = new TH2F("fh2_aoq_1_beta","aoq vs beta",200,0.5,0.6,200,2.6,2.7);	   
                   
  fh2_aoq_2_angle_F10 = new TH2F("fh2_aoq_2_angle_F10","aoq vs angle F10",200,-10,10,200,2.6,2.7);
  fh2_aoq_2_angle_F11 = new TH2F("fh2_aoq_2_angle_F11","aoq vs angle F11",200,-10,10,200,2.6,2.7);
  fh2_aoq_2_xpos_F10 = new TH2F("fh2_aoq_2_xpos_F10","aoq vs xpos F10",200,-10,10,200,2.6,2.7);
  fh2_aoq_2_xpos_F11 = new TH2F("fh2_aoq_2_xpos_F11","aoq vs xpos F11",200,-10,10,200,2.6,2.7);
  fh2_aoq_2_ypos_F10 = new TH2F("fh2_aoq_2_ypos_F10","aoq vs ypos F10",200,-10,10,200,2.6,2.7);
  fh2_aoq_2_ypos_F11 = new TH2F("fh2_aoq_2_ypos_F11","aoq vs ypos F11",200,-10,10,200,2.6,2.7);
  fh2_aoq_2_beta = new TH2F("fh2_aoq_2_beta","aoq vs beta",200,0.5,0.6,200,2.6,2.7);	   

  fh2_aoq_3_angle_F8  = new TH2F("fh2_aoq_3_angle_F8","aoq vs angle F8",200,-10,10,200,2.6,2.7);
  fh2_aoq_3_angle_F10 = new TH2F("fh2_aoq_3_angle_F10","aoq vs angle F10",200,-10,10,200,2.6,2.7);
  fh2_aoq_3_angle_F11 = new TH2F("fh2_aoq_3_angle_F11","aoq vs angle F11",200,-10,10,200,2.6,2.7);
  fh2_aoq_3_xpos_F8   = new TH2F("fh2_aoq_3_xpos_F8","aoq vs xpos F8",200,-10,10,200,2.6,2.7); 
  fh2_aoq_3_xpos_F10  = new TH2F("fh2_aoq_3_xpos_F10","aoq vs xpos F10",200,-10,10,200,2.6,2.7);
  fh2_aoq_3_xpos_F11  = new TH2F("fh2_aoq_3_xpos_F11","aoq vs xpos F11",200,-10,10,200,2.6,2.7);
  fh2_aoq_3_ypos_F8   = new TH2F("fh2_aoq_3_ypos_F8","aoq vs ypos F8",200,-10,10,200,2.6,2.7); 
  fh2_aoq_3_ypos_F10  = new TH2F("fh2_aoq_3_ypos_F10","aoq vs ypos F10",200,-10,10,200,2.6,2.7);
  fh2_aoq_3_ypos_F11  = new TH2F("fh2_aoq_3_ypos_F11","aoq vs ypos F11",200,-10,10,200,2.6,2.7);
  fh2_aoq_3_beta      = new TH2F("fh2_aoq_3_beta","aoq vs beta",200,0.5,0.6,200,2.6,2.7);	   



}


Double_t PID::Correction_3(short aoq, short var){


  Double_t naoq = AOQ[aoq];
  
  for(short iplane = 0; iplane < 3; iplane++){
    for(short ivar = 0; ivar < 3; ivar++){
      if( iplane == var/3 && ivar ==  var%3) continue;
      naoq = naoq - fvar[iplane][ivar]*fcorr_Lin_3[iplane][ivar];
    }//for(iplane)
  }//for(ivar)
  
  if( var != 10 )
    naoq = naoq - BETA[3]*fcorr_Lin_3[4][0];


  return naoq;

}//Correction_3

Double_t PID::Correction(short aoq, short var){

  Double_t naoq = AOQ[aoq];
  if(aoq>0){
    for(short jj = 0; jj < 6; jj++){
      if( jj == var ) continue;
      naoq = naoq - fout_var[aoq-1][jj]*fcorr_Lin[aoq][jj];
    }
  }else{
    for(short jj = 0; jj < 6; jj++){
      if( jj == var ) continue;
      naoq = naoq - fin_var[jj]*fcorr_Lin[aoq][jj];
    }
  }
  if( var != 6 )
    naoq = naoq - BETA[aoq]*fcorr_Lin[aoq][6];


  return naoq;

}//correction()

void PID::Write(Char_t* name){

  TFile* out = new TFile(Form("%s",name),"RECREATE");

  fh2_aoq_0_angle_F3->Write();
  fh2_aoq_0_angle_F5->Write();
  fh2_aoq_0_xpos_F3 ->Write(); 
  fh2_aoq_0_xpos_F5 ->Write(); 
  fh2_aoq_0_ypos_F3 ->Write(); 
  fh2_aoq_0_ypos_F5 ->Write(); 
  fh2_aoq_0_beta    ->Write();	   
                   
  fh2_aoq_1_angle_F8 ->Write();
  fh2_aoq_1_angle_F10->Write();
  fh2_aoq_1_xpos_F8  ->Write(); 
  fh2_aoq_1_xpos_F10 ->Write();
  fh2_aoq_1_ypos_F8  ->Write(); 
  fh2_aoq_1_ypos_F10 ->Write();
  fh2_aoq_1_beta     ->Write();	   
                   
  fh2_aoq_2_angle_F10->Write();
  fh2_aoq_2_angle_F11->Write();
  fh2_aoq_2_xpos_F10 ->Write();
  fh2_aoq_2_xpos_F11 ->Write();
  fh2_aoq_2_ypos_F10 ->Write();
  fh2_aoq_2_ypos_F11 ->Write();
  fh2_aoq_2_beta     ->Write();	   
                 
  fh2_aoq_3_angle_F8 ->Write();
  fh2_aoq_3_angle_F10->Write();
  fh2_aoq_3_angle_F11->Write();
  fh2_aoq_3_xpos_F8  ->Write();
  fh2_aoq_3_xpos_F10 ->Write();
  fh2_aoq_3_xpos_F11 ->Write();
  fh2_aoq_3_ypos_F8  ->Write();
  fh2_aoq_3_ypos_F10 ->Write();
  fh2_aoq_3_ypos_F11 ->Write();
  fh2_aoq_3_beta     ->Write();	   

  for(short jj = 0; jj < 4; jj++)
    fh1_CorrAoQ[jj]->Write();

  out->Write();
                   
}

void PID::AoQ(){

  for(short ii = 0; ii < 4; ii++)
    fh1_CorrAoQ[ii] = new TH1F(Form("fh1_CorrAoQ_%i",ii),Form("Corr AOQ %i",ii),400,2.5,2.7);

  if (fChain == 0) return;
  
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    if( TMath::Abs( ZET[0] -50 ) > 0.4)continue;
    fh1_CorrAoQ[0]->Fill(Correction(0,100));
    if( globalcut() &&  TMath::Abs(ZET[2]-45.2)<0.7  ){
      for(short ii = 1; ii < 3; ii++)
	fh1_CorrAoQ[ii]->Fill(Correction(ii,100));
      
      fh1_CorrAoQ[3]->Fill(Correction_3(3,100));
    }//if
  
  }//for()

}//AoQ


void PID::Reset(){

  fh2_aoq_0->Delete();
  fh2_aoq_1->Delete();
  fh2_aoq_2->Delete();
  fh2_aoq_3->Delete();

  fh2_aoq_0_angle_F3->Delete();
  fh2_aoq_0_angle_F5->Delete();
  fh2_aoq_0_xpos_F3->Delete(); 
  fh2_aoq_0_xpos_F5->Delete(); 
  fh2_aoq_0_ypos_F3->Delete(); 
  fh2_aoq_0_ypos_F5->Delete(); 
  fh2_aoq_0_beta->Delete();	   
                   
  fh2_aoq_1_angle_F8->Delete();
  fh2_aoq_1_angle_F10->Delete();
  fh2_aoq_1_xpos_F8->Delete(); 
  fh2_aoq_1_xpos_F10->Delete();
  fh2_aoq_1_ypos_F8->Delete(); 
  fh2_aoq_1_ypos_F10->Delete();
  fh2_aoq_1_beta->Delete();	   
                   
  fh2_aoq_2_angle_F10->Delete();
  fh2_aoq_2_angle_F11->Delete();
  fh2_aoq_2_xpos_F10->Delete();
  fh2_aoq_2_xpos_F11->Delete();
  fh2_aoq_2_ypos_F10->Delete();
  fh2_aoq_2_ypos_F11->Delete();
  fh2_aoq_2_beta->Delete();	   


  
  fh2_aoq_3_angle_F8->Delete();
  fh2_aoq_3_angle_F10->Delete();
  fh2_aoq_3_angle_F11->Delete();

  fh2_aoq_3_xpos_F8->Delete();
  fh2_aoq_3_xpos_F10->Delete();
  fh2_aoq_3_xpos_F11->Delete();

  fh2_aoq_3_ypos_F8->Delete();
  fh2_aoq_3_ypos_F10->Delete();
  fh2_aoq_3_ypos_F11->Delete();

  fh2_aoq_3_beta->Delete();	   
  
}

void PID::Read(){

  ReadIncoming();
  ReadOutgoing_1();
  ReadOutgoing_2();
  ReadOutgoing_3();

}

void PID::ReadIncoming(){


  Char_t* path;
  path = "cutParameters/PID/Sn132";
  
  ifstream aoqfile;
  std::string aoq, line;
  aoqfile.open(Form("%s/AOQ_correction_in.txt",path));
  short ii =0;
  while(aoqfile.good()){
    getline(aoqfile,line);
    std::istringstream iss(line);
    iss >> aoq;
    fcorr_Lin[0][ii]=atof(aoq.c_str());
    ii++;
    if(ii == 6)continue;
  }
}

void PID::ReadOutgoing_1(){

  Char_t* path;
  path = "cutParameters/PID/Sn132";
  
  ifstream aoqfile;
  std::string aoq, line;
  aoqfile.open(Form("%s/AOQ_Correction_out1.txt",path));
  short ii =0;
  while(aoqfile.good()){
    getline(aoqfile,line);
    std::istringstream iss(line);
    iss >> aoq;
    fcorr_Lin[1][ii]=atof(aoq.c_str());
    ii++;
    if(ii == 6)continue;
  }


}

void PID::ReadOutgoing_2(){

  Char_t* path;
  path = "cutParameters/PID/Sn132";
  std::string aoq, line;
  ifstream aoqfile;
  aoqfile.open(Form("%s/AOQ_Correction_out2.txt",path));
  short ii =0;
  while(aoqfile.good()){

    getline(aoqfile,line);
    std::istringstream iss(line);
    iss >> aoq;
    fcorr_Lin[2][ii]=atof(aoq.c_str());
    ii++;
    if(ii == 6)continue;
  }


}

void PID::ReadOutgoing_3(){

  Char_t* path;
  path = "cutParameters/PID/Sn132";
  
  std::string var1,var2,var3,var4, line;
  ifstream aoqfile;
  aoqfile.open(Form("%s/AOQ_Correction_out3.txt",path));
  short ii =0;
  while(aoqfile.good()){

    getline(aoqfile,line);
    std::istringstream iss(line);
    iss >> var1 >> var2 >> var3 >> var4;
    fcorr_Lin_3[0][ii]=atof(var1.c_str());
    fcorr_Lin_3[1][ii]=atof(var2.c_str());
    fcorr_Lin_3[2][ii]=atof(var3.c_str());
    fcorr_Lin_3[3][ii]=atof(var4.c_str());

    ii++;
    if(ii == 3)continue;
  }



}



