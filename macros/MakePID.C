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
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"

#include "TVector3.h"
#include "TMath.h"

#include "signal.h"
#include "TSpectrum.h"


void MakePID(){

  FILE* aoq_in = fopen("cutParameters/PID/Sn132/peak_aoq_in.txt","w");
  FILE* aoq_out = fopen("cutParameters/PID/Sn132/peak_aoq_out.txt","w");
  FILE* zet_in = fopen("cutParameters/PID/Sn132/peak_zet_in.txt","w");
  FILE* zet_out = fopen("cutParameters/PID/Sn132/peak_zet_out.txt","w");

  TH1F* h1_AOQ[64][2];
  TH1F* h1_ZET[64][2];
  TCanvas* can_AOQ[4][2];
  TCanvas* can_ZET[4][2];
  for(Int_t ii = 0; ii<4;ii++){
    can_AOQ[ii][0] = new TCanvas(Form("can_AOQ_0_%i",ii), "AOQ 0");
    can_AOQ[ii][0] ->Divide(4,4);
    can_AOQ[ii][1] = new TCanvas(Form("can_AOQ_2_%i",ii), "AOQ 2");
    can_AOQ[ii][1] ->Divide(4,4);
    can_ZET[ii][0] = new TCanvas(Form("can_ZET_0_%i",ii), "ZET 0");
    can_ZET[ii][0] ->Divide(4,4);
    can_ZET[ii][1] = new TCanvas(Form("can_ZET_2_%i",ii), "ZET 2");
    can_ZET[ii][1] ->Divide(4,4);
  }

  Int_t ican = 0, ipad = 1, nentries = (int)3e5;//run = 9, 
  Float_t aoq = 2.64, zet = 50.;
  TF1* f1_aoq[64][2],*f1_zet[64][2];
  for(Int_t run = 9; run<64;run++){

    for(short ii = 0; ii < 2; ii++){     
      f1_aoq[run][ii] = new TF1(Form("f1_aoq_%i_%02d",ii,run),"gaus(0)+pol0(3)",2.63,2.665);
      f1_zet[run][ii] = new TF1(Form("f1_zet_%i_%02d",ii,run),"gaus(0)+pol0(3)",48,52);
      f1_zet[run][ii]->SetParLimits(0,10,1e6);
      f1_zet[run][ii]->SetParameter(1,12e3);
      f1_zet[run][ii]->SetParameter(1,zet);
      f1_zet[run][ii]->SetParLimits(1,zet-1.5,zet+4.5);
      f1_zet[run][ii]->SetParameter(2,0.5);
      f1_zet[run][ii]->SetParLimits(2,0,1.);

      f1_aoq[run][ii]->SetParLimits(0,10,1e6);
      f1_aoq[run][ii]->SetParameter(1,12e3);
      f1_aoq[run][ii]->SetParameter(1,aoq);
      f1_aoq[run][ii]->SetParLimits(1,aoq-0.1,aoq+0.1);
      f1_aoq[run][ii]->SetParameter(2,2e-3);
      f1_aoq[run][ii]->SetParLimits(2,0,5e-3);

    }//for (ii)

    TFile* _file = TFile::Open(Form("data/rootfiles/new/run%04d.root",run));
    TTree* tree = (TTree*)_file->Get("tree");
    h1_AOQ[run][0] = new TH1F(Form("h1_AOQ_%i_0",run),Form("run %i, AOQ 0", run),400,2.5,2.7);
    h1_AOQ[run][1] = new TH1F(Form("h1_AOQ_%i_2",run),Form("run %i, AOQ 2", run),400,2.5,2.7);
    h1_ZET[run][0] = new TH1F(Form("h1_ZET_%i_0",run),Form("run %i, ZET 0", run),200,45,55);
    h1_ZET[run][1] = new TH1F(Form("h1_ZET_%i_2",run),Form("run %i, ZET 2", run),200,45,55);

    tree->Draw(Form("AOQ[0]>>h1_AOQ_%i_0",run),"TMath::Abs(ZET[0]-50)<0.5","",nentries,0);
    tree->Draw(Form("ZET[0]>>h1_ZET_%i_0",run),"TMath::Abs(AOQ[0]-2.64)<0.005","",nentries,0);
    tree->Draw(Form("AOQ[2]>>h1_AOQ_%i_2",run),"TMath::Abs(ZET[0]-50)<0.5 && TMath::Abs(AOQ[0]-2.64)<0.005 && TMath::Abs(delta[2]-delta[4])<0.5","",nentries,0);
    tree->Draw(Form("ZET[2]>>h1_ZET_%i_2",run),"TMath::Abs(ZET[0]-50)<0.5 && TMath::Abs(AOQ[0]-2.64)<0.005 && TMath::Abs(AOQ[2]-2.64)<0.005","",nentries,0);

    can_AOQ[ican][0]->cd(ipad);
    h1_AOQ[run][0]->Draw();
    h1_AOQ[run][0]->Fit(f1_aoq[run][0]);
    can_AOQ[ican][1]->cd(ipad);
    h1_AOQ[run][1]->Draw();
    h1_AOQ[run][1]->Fit(f1_aoq[run][1]);

    can_ZET[ican][0]->cd(ipad);
    h1_ZET[run][0]->Draw();
    h1_ZET[run][0]->Fit(f1_zet[run][0]);
    can_ZET[ican][1]->cd(ipad);
    h1_ZET[run][1]->Draw();
    h1_ZET[run][1]->Fit(f1_zet[run][1]);

    if(ipad==16){
      ipad = 1;
      ican++;
    }else   ipad++;
    fprintf(aoq_in,"%i %f %f \n", run, f1_aoq[run][0]->GetParameter(1),f1_aoq[run][0]->GetParameter(2));
    fprintf(aoq_out,"%i %f %f \n", run, f1_aoq[run][1]->GetParameter(1),f1_aoq[run][1]->GetParameter(2));
    fprintf(zet_in,"%i %f %f \n", run, f1_zet[run][0]->GetParameter(1),f1_zet[run][0]->GetParameter(2));
    fprintf(zet_out,"%i %f %f \n", run, f1_zet[run][1]->GetParameter(1),f1_zet[run][1]->GetParameter(2));



    //delete _file;
    //delete tree;
  }

  fclose(aoq_in);
  fclose(aoq_out);
  fclose(zet_in);
  fclose(zet_out);

}
