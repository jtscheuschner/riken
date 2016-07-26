/***********************************************************

This class optimizes the AOQ values in dependence of 

- BETA
                    and two times (First and Second PPAC)
- DELTA
- X
- Y
- Angle

The timeoffset has to be adapted to match the designevalue....


 **********************************************************/

#define Imp_PID_cxx
#include "Imp_PID.h"
#include <iostream>
#include "Riostream.h"
#include <string>
#include <time.h>
#include <sstream>

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TObject.h>
#include <TRandom3.h>

//#include "TArtStoreManager.hh"
//#include "TArtEventStore.hh"
#include "TClonesArray.h"
#include "TSystem.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TCanvas.h"

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

void Imp_PID::Loop(short aoq, short max){
  //Setting the variables according to the choice of AOQ

  LoadData(0,0);

  vector<Char_t*>names;

  names = vector<Char_t*>(variables);
  short plane_[2]={0x0}, delta_[2]={0x0}, beta_=0;

  if(aoq == 0){
    names[0]="FA[0]",  names[1]="FX[0]",  names[2]="FY[0]",  names[3]="delta[0]";//First detector
    names[4]="FA[1]",  names[5]="FX[1]",  names[6]="FY[1]",  names[7]="delta[1]";//Second detector
    names[8]="BETA[0]";                                                       //Beta before ([0]) or after ([1]) the target
    plane_[0] = 0;
    plane_[1] = 1;
    delta_[0] = 0;
    delta_[1] = 1;
    beta_= 0;
  }else if(aoq == 2){
    names[0]="FA[4]",  names[1]="FX[4]",  names[2]="FY[4]",  names[3]="delta[2]";//First detector
    names[4]="FA[5]",  names[5]="FX[5]",  names[6]="FY[5]",  names[7]="delta[4]";//Second detector
    names[8]="BETA[1]";                                                       //Beta before ([0]) or after ([1]) the target
    plane_[0] = 4;
    plane_[1] = 5;
    delta_[0] = 2;
    delta_[1] = 4;
    beta_= 1;
  }else if(aoq == 5){
    names[0]="FA[3]",  names[1]="FX[3]",  names[2]="FY[3]",  names[3]="delta[2]";//First detector
    names[4]="FA[5]",  names[5]="FX[5]",  names[6]="FY[5]",  names[7]="delta[4]";//Second detector
    names[8]="BETA[1]";                                                      //Beta before ([0]) or after ([1]) the target
    plane_[0] = 3;
    plane_[1] = 5;
    delta_[0] = 2;
    delta_[1] = 4;
    beta_= 1;
  }

    std::vector<TH1F*> h1;
    h1 =std::vector<TH1F*> (4);
    for(short ii=0; ii< 4; ii++){
      h1[ii] = new TH1F(Form("h1_%i",ii), Form(" [%i] : %i ",aoq, ii),100,-5,5);
    }
    std::vector<TH2F*> h2_AOQ;
    h2_AOQ =std::vector<TH2F*> (variables);

    for(short ii=0; ii< variables; ii++){
      if( ii ==8) h2_AOQ[ii] = new TH2F(Form("h2_AOQ_%i",ii), Form(" AOQ[%i] : %s ",aoq, names[ii]),100,0.52,0.6,400,2.5,2.7);
      else if( ii == 3 || ii == 7) h2_AOQ[ii] = new TH2F(Form("h2_AOQ_%i",ii), Form(" AOQ[%i] : %s ",aoq, names[ii]),100,-2,2,400,2.5,2.7);
      else h2_AOQ[ii] = new TH2F(Form("h2_AOQ_%i",ii), Form(" AOQ[%i] : %s ",aoq, names[ii]),100,-10,10,400,2.5,2.7);
    }

    //fChain->Draw("delta[1]>>h(100,-5,5)");

    for(short ii = 0; ii < variables; ii++)
      fChain->Draw(Form("AOQ[%i]:%s>>h2_AOQ_%i",aoq,names[ii],ii),"","colz");


    for(short iteration = 0; iteration<max;iteration++){

      if (fChain == 0) return;

      Long64_t nentries = fChain->GetEntriesFast();
      Long64_t nbytes = 0, nb = 0;

      for (Long64_t jentry=0; jentry<nentries;jentry++) {
	Long64_t ientry = LoadTree(jentry);

	//	std::cout<< delta[2] << "  "<< delta[1] << "  "<< delta[2] << "  "<< delta[3] << "  "<< delta[4] << "  " <<std::endl;
	if (ientry < 0) break;
	nb = fChain->GetEntry(jentry);   nbytes += nb;

	// if (Cut(ientry) < 0) continue;

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
	dummy_var[3] = delta[delta_[0]];
	dummy_var[7] = delta[delta_[1]];

	h1[1]->Fill(dummy_var[3]);
	h1[3]->Fill(dummy_var[7]);
	h1[0]->Fill(delta[0]);
	h1[2]->Fill(delta[1]);

	//BETA
	dummy_var[8] = BETA[beta_];

	//First Plane
	Float_t corr = 0;
	for(short ii = 0; ii<variables;ii++){
	  corr += +dummy_var[ii]*correction[0][ii]+dummy_var[ii]*dummy_var[ii]*correction[1][ii];
	}

	for(short jj = 0; jj< variables; jj++){

	  corr = corr - dummy_var[jj]*correction[0][jj]-dummy_var[jj]*dummy_var[jj]*correction[1][jj];//subtract current variable
	  /* if( jj != 7)*/ h2_AOQ[jj]->Fill(dummy_var[jj],AOQ[aoq]);                      //DRAW the histo
	  //else h2_AOQ[7]->Fill(delta[2],AOQ[aoq]);
	  corr = corr + dummy_var[jj]*correction[0][jj]+dummy_var[jj]*dummy_var[jj]*correction[1][jj];//add again the correction

	}
      
      }//for(event)
    
      TCanvas* can = new TCanvas("can","AOQ:var");
      can->Divide(3,3);

      TCanvas* can1 = new TCanvas("can1","var_Slices");
      can1->Divide(3,3);
      //Fit
      for(short ii=0; ii < variables; ii++){
	can->cd(ii+1);
	h2_AOQ[ii]->Draw("colz");
	can1->cd(ii+1);
	if(iteration+1 < max )Fitter(h2_AOQ[ii],ii,kTRUE);//Histo, variable number, clear

	else Fitter(h2_AOQ[ii],ii,kFALSE);//Histo, variable number, clear
	
      }
    }//vaiables


    TCanvas* can_aoq = new TCanvas("can_aoq","aoq_corrected");
    
    Char_t* corrLin = "", *corrSqu = "";
    for(short ii = 0; ii < variables; ii++){
      corrLin  = Form("%s+%s*%f",corrLin,names[ii],correction[0][ii]);
      corrSqu  = Form("%s+%s*%s*%f",corrSqu,names[ii],names[ii],correction[0][ii]);
    }
    TH1F* h1_AOQ = new TH1F("h1_AOQ","AOQ",400,2.5,2.7);
    std::cout<< corrLin <<std::endl;
    std::cout<< corrSqu <<std::endl;
    fChain->Draw(Form("AOQ[%i]-%s-%s>>h1_AOQ",aoq,corrLin,/*"0"*/corrSqu),"","colz");
    //    can_aoq->Divide(2,2);
    //for(short ii=0;ii<4;ii++){
    //can_aoq->cd(ii+1);
    //h1[ii]->Draw();
    //}



}


    void Imp_PID::Fitter(TH2F* h2, Int_t var, Bool_t clear){
  
      /************

This function Slices the h2 and then Fits the Slices with pol2

  ***********/
      if(var!=7)return;
  std::cout<< " var:  " << var <<std::endl;
  h2->GetYaxis()->SetRangeUser(2.55,2.57);
  TObjArray* dummy = new TObjArray();//this is needed to delete properly the fitslices histograms
  dummy->SetOwner(kTRUE);
  TF1* f1 = new TF1("f1","gaus",2.54,2.58);
  h2->FitSlicesY(f1,0,-1,5,"QNRG4",dummy);  //QNR default, Q: quiet, N: do not draw do not store drawing, R: use range Gx: Merges x(x number between 2 and 5) bins before Projections
  h1[var] =  (TH1F*) dummy->FindObject(Form("h2_AOQ_%i_1",var));

  TF1* func = new TF1("func","pol2",-9,9);
  //func->SetParameter(0,AOQ_value);
  //func->SetParLimits(0,2.556,2.564);
  func->SetParameters(2.56,correction[0][var],correction[1][var]);
  if(var!=8)h1[var]->Fit(func,"","",-9,9);
  else return;
  correction[0][var] = func->GetParameter(1);
  correction[1][var] = func->GetParameter(2);    

  if(clear)
    dummy->Clear();
}


