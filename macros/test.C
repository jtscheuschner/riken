#define test_cxx
#include "test.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include "Riostream.h"
#include <string>
#include <time.h>
#include <sstream>

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TObject.h>
#include <TRandom3.h>

//#include "TArtStoreManager.hh"
//#include "TArtEventStore.hh"
#include "TClonesArray.h"
#include "TSystem.h"
#include "TChain.h"
#include "TH1F.h"
#include "TF1.h"

#include "TVector3.h"
#include "TMath.h"

#include "signal.h"
#include "TSpectrum.h"

//#include "hist.h"
#include "analyse.h"

#include<fstream>
#include<stdlib.h>
#include "TString.h"
#include "TSpectrum.h"
#include "TPaveStats.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TMath.h"
#include "TLine.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TLatex.h"
#include <vector>

void test::Loop(short aoq, short max){

//   In a ROOT session, you can do:
//      Root > .L test.C
//      Root > test t
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


  vector<Char_t*>names;

  names = vector<Char_t*>(variables);
  short plane_[2]={0x0}, delta_[2]={0x0}, beta_=0;

  if(aoq == 0){
    names[0]="F3A",  names[1]="F3X",  names[2]="F3Y",  names[3]="delta[0]";//First detector
    names[4]="F5A",  names[5]="F5X",  names[6]="F5Y",  names[7]="delta[1]";//Second detector
    names[8]="BETA[0]";                                                       //Beta before ([0]) or after ([1]) the target
    plane_[0] = 0;
    plane_[1] = 1;
    delta_[0] = 0;
    delta_[1] = 1;
    beta_= 0;
  }else if(aoq == 2){
    names[0]="F10A",  names[1]="F10X",  names[2]="F10Y",  names[3]="delta[2]";//First detector
    names[4]="F11A",  names[5]="F11X",  names[6]="F11Y",  names[7]="delta[4]";//Second detector
    names[8]="BETA[1]";                                                       //Beta before ([0]) or after ([1]) the target
    plane_[0] = 4;
    plane_[1] = 5;
    delta_[0] = 2;
    delta_[1] = 4;
    beta_= 1;
  }else if(aoq == 5){
    names[0]="F8A",  names[1]="F8X",  names[2]="F8Y",  names[3]="delta[2]";//First detector
    names[4]="F11A",  names[5]="F11X",  names[6]="F11Y",  names[7]="delta[4]";//Second detector
    names[8]="BETA[1]";                                                      //Beta before ([0]) or after ([1]) the target
    plane_[0] = 3;
    plane_[1] = 5;
    delta_[0] = 2;
    delta_[1] = 4;
    beta_= 1;
  }

  //h1 for testing
    std::vector<TH1F*> h1;
    h1 =std::vector<TH1F*> (4);
    for(short ii=0; ii< 4; ii++){
      h1[ii] = new TH1F(Form("h1_%i",ii), Form(" [%i] : %i ",aoq, ii),100,-5,5);
    }

  //TH2 to fill aoq vs variable
    std::vector<TH2F*> h2_AOQ;
    h2_AOQ =std::vector<TH2F*> (variables);

    for(short ii=0; ii< variables; ii++){
      if( ii ==8) h2_AOQ[ii] = new TH2F(Form("h2_AOQ_%i",ii), Form(" AOQ[%i] : %s ",aoq, names[ii]),100,0.52,0.6,400,2.5,2.7);
      else if( ii == 3 || ii == 7) h2_AOQ[ii] = new TH2F(Form("h2_AOQ_%i",ii), Form(" AOQ[%i] : %s ",aoq, names[ii]),100,-2,2,400,2.5,2.7);
      else h2_AOQ[ii] = new TH2F(Form("h2_AOQ_%i",ii), Form(" AOQ[%i] : %s ",aoq, names[ii]),100,-10,10,400,2.5,2.7);
    }



   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      //set the variables of interest into a dummy array to loop over those
	Double_t dummy_var[variables] = {0x0};
	//First Plane
	dummy_var[0] = FA[plane_[0]];
	dummy_var[1] = FX[plane_[0]];
	dummy_var[2] = FY[plane_[0]];

	//Second Plane
	dummy_var[4] = FA[plane_[1]];
	dummy_var[5] = FX[plane_[1]];
	dummy_var[6] = FY[plane_[1]];

	//delta
	dummy_var[3] = delta[0];
	dummy_var[7] = delta[1];
	//BETA
	dummy_var[8] = BETA[beta_];

	Float_t corr = 0;
	for(short ii = 0; ii<variables;ii++){
	  corr += +dummy_var[ii]*correction[0][ii]+dummy_var[ii]*dummy_var[ii]*correction[1][ii];
	}

	for(short jj = 0; jj< variables; jj++){

	  corr = corr - dummy_var[jj]*correction[0][jj]-dummy_var[jj]*dummy_var[jj]*correction[1][jj];//subtract current variable
	  h2_AOQ[jj]->Fill(dummy_var[jj],AOQ[aoq]);                                              //DRAW the histo
	   corr = corr + dummy_var[jj]*correction[0][jj]+dummy_var[jj]*dummy_var[jj]*correction[1][jj];//add again the correction

	}


   }//event

      TCanvas* can = new TCanvas("can","AOQ:var");
      can->Divide(3,3);

      TCanvas* can1 = new TCanvas("can1","var_Slices");
      can1->Divide(3,3);
      //Fit
      for(short ii=0; ii < variables; ii++){
	can->cd(ii+1);
	h2_AOQ[ii]->Draw("colz");
      }//filling canvas


    TCanvas* can_aoq = new TCanvas("can_aoq","aoq_corrected");
    can_aoq->Divide(2,2);
    for(short ii=0;ii<4;ii++){
      can_aoq->cd(ii+1);
      h1[ii]->Draw();
    }


}
