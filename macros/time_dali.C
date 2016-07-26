#include <iostream>
#include <string>
#include <time.h>
#include <sstream>
#include "Riostream.h"

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
#include "TDirectory.h"
#include "TList.h"


void time_dali(){

  FILE* file_dali = fopen("cutParameters/time_dali2.txt","w");
  TFile* _file = new TFile("data/Histo/Sn132/15.0_3.0_3.0_2.5_2.0_0.8_1.0_2.0_-1.0.root");

  TH2F* h2[150];
  TH1F* h1[150];
  TF1 * f1[150];
  TCanvas* can[10];
  for(Int_t jj = 0; jj<10;jj++){
    can[jj] = new TCanvas(Form("can_%i",jj),"");
    can[jj]->Divide(4,4);
  }
  Int_t ipad = 1, ican=0;
  for(Int_t ii = 0; ii<150;ii++){
    cout<< ii <<endl;
    f1[ii] = new TF1(Form("f%i",ii),"gaus(0)+pol0(3)",-20,20);
    //    f1[ii] = new TF1(Form("f%i",ii),"gaus",400,510);
    f1[ii]->SetParameters(50,50,10,100);
    f1[ii]->SetParLimits(0,5,200);
    f1[ii]->SetParLimits(2,1,20);
    h2[ii] = (TH2F*) _file ->Get(Form("top/DALI/h2_dali_E_T_%i",ii));
    h1[ii] = (TH1F*)h2[ii]->ProjectionY(Form("%i_py",ii),0,-1);
    can[ican]->cd(ipad);
    h1[ii]->Draw();
    h1[ii]->Fit(f1[ii]);
    fprintf(file_dali,"%i %f %f \n",ii,f1[ii]->GetParameter(1),f1[ii]->GetParameter(2));
    if(ipad==16){
      ican++;
      ipad = 1;
    }else ipad++;
  }//for
  fclose(file_dali);

}//time_dali

void time_QTC(){

  FILE* file_qtc = fopen("cutParameters/QTC_parameter_n1.txt","w");
  FILE* file_qtcl = fopen("cutParameters/QTC_parameter_n2.txt","w");
  TFile* _file = new TFile("data/Histo/Sn132/0.0_3.0_2.5_2.0_1.5_0.8_1.0_2.0_-1.0.root");//data/Histo/Sn132_15.0_2.0_2.0_2.0_2.0_0.8_1.0_2.0_-1.0.root");

  TH2F* h2[8];
  TH2F* h2l[8];
  TH1F* h1[8];
  TH1F* hl[8];
  TF1 * f1[8];
  TF1 * fl[8];
  TCanvas* can = new TCanvas("can","");
  can->Divide(3,3);
  TCanvas* canl = new TCanvas("canl","");
  canl->Divide(3,3);
  TCanvas* can2 = new TCanvas("can2","");
  can2->Divide(3,3);

  for(Int_t ii = 0; ii<8;ii++){
    cout<< ii <<endl;
    h2[ii] = (TH2F*) _file ->Get(Form("top/LaBr/h2_labr_E_T_MID_%i",ii));
    h1[ii] = (TH1F*)h2[ii]->ProjectionY(Form("%i_py",ii),35,-1);
    h1[ii]->GetXaxis()->SetRangeUser(-15,5);
    h2l[ii] = (TH2F*) _file ->Get(Form("top/LaBr/h2_labr_E_T_LOW_%i",ii));
    hl[ii] = (TH1F*)h2l[ii]->ProjectionY(Form("%i_pl",ii),17,-1);
    hl[ii]->GetXaxis()->SetRangeUser(-15,5);

    f1[ii] = new TF1(Form("f%i",ii),"gaus(0)+pol0(3)",-20,2);
    f1[ii]->SetParameters(1e2,-7,0.3,2);
    f1[ii]->SetParLimits(2,0.1,0.6);

    fl[ii] = new TF1(Form("f%i",ii),"gaus(0)+pol0(3)",-20,2);
    fl[ii]->SetParameters(1e2,-7,0.3,2);
    fl[ii]->SetParLimits(2,0.1,0.6);

    can->cd(ii+1);
    h1[ii]->Draw();
    h1[ii]->Fit(f1[ii],"R");
    fprintf(file_qtc,"%i %f %f \n",ii,f1[ii]->GetParameter(1),f1[ii]->GetParameter(2));
    canl->cd(ii+1);
    hl[ii]->Draw();
    hl[ii]->Fit(fl[ii],"R");
    fprintf(file_qtcl,"%i %f %f \n",ii,fl[ii]->GetParameter(1),fl[ii]->GetParameter(2));

    can2->cd(ii+1);
    h2[ii]->Draw("colz");
  }//for

  fclose(file_qtc);
  fclose(file_qtcl);

}
