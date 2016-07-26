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

void dali(short start, short end){
  TChain* chain = new TChain("tree");
  Char_t* path = "data/rootfiles/new_";
  for(short ii = start; ii < end; ii++)chain->AddFile(Form("%s/run%04d.root",path,ii));

  FILE *file;
  file = fopen("cutParameters/DALI_t_offset","w");

  TH2F* hadc[160];
  TH1F* h1[160];
  TCanvas* can[10];
  for(UInt_t ican=0;ican<10;ican++){
    can[ican] = new TCanvas(Form("can_%i",ican),Form("crystal %i - %i",ican*16,(ican+1)*16));
    can[ican]->Divide(4,4);
  }
  TCanvas* can_dummy = new TCanvas("can_dummy","");
  UInt_t ican = 0, ipad=1;
  for(UInt_t crystal = 0; crystal<160;crystal++){
    cout<< crystal <<endl;
    Float_t t_start=1.3e3, t_end=1.55e3;
    Int_t nbins=125;
    if(crystal <124){
      t_start = 1.4e3;
      t_end = 1.6e3;
      nbins = 100;
    }

    hadc[crystal]=new TH2F(Form("hadc_%i",crystal),Form("hadc %i",crystal),300,0,30e3,nbins,t_start,t_end);
    can_dummy->cd(0);
    chain->Draw(Form("DALINaI.fTimeOffseted-F7T+F8T-125./30./BETA[0]:DALINaI.fEnergy>>hadc_%i",crystal),
		Form("DALINaI.id==%i && TMath::Abs(AOQ[0]-2.64)<0.008 && TMath::Abs(ZET[0]-50)<0.5 && TMath::Abs(AOQ[3]-2.64)<0.005 && TMath::Abs(ZET[3]-50)<0.3",crystal),"colz");
    if(hadc[crystal]->GetEntries()<100){
      cout<< " no entries " <<endl;
      fprintf(file,"%i %f \n", crystal, -9999.);
      ipad++;
      if(ipad==17){
	ican++;
	ipad=1;
      }
      continue;
    }
    cout<< " found entries " <<endl;
    h1[crystal] = (TH1F*)hadc[crystal]->ProjectionY(Form("h1_%i_py",crystal),0,-1,"");
    TF1* f1 = new TF1("f1","gaus(0)+pol0(3)",500,2000);
    if(crystal<124){
      f1->SetParameter(1,1500);
      f1->SetParameter(2,15);
    }else{
      f1->SetParameter(1,1440);
      f1->SetParameter(2,9);
      f1->SetParameter(3,60);
    }
    h1[crystal]->Fit(f1,"","",t_start,t_end);
    fprintf(file,"%i %f\n", crystal,f1->GetParameter(1));
    can[ican]->cd(ipad);
    h1[crystal]->Draw();
    ipad++;
    if(ipad==17){
      ican++;
      ipad=1;
    }
    cout<< "done" <<endl;
    delete f1;
  }
  fclose(file);
}

void labr(short start, short end, short trig, short crystal = 0){
  TChain* chain = new TChain("tree");
  Char_t* path = "data/rootfiles/new2";
  for(short ii = start; ii < end; ii++)chain->AddFile(Form("%s/run%04d.root",path,ii));

  TCanvas* can[4];
  for(short ii = 0; ii<4; ii++){
    can[ii] = new TCanvas(Form("can_%i",ii),"");
    can[ii]->Divide(2,1);
  }
  Float_t mlength = 124.8;//cm
  Float_t clight = 29.979245;//  cm/ns
  Float_t length = mlength/clight;

  TH2F* h2[3];
  h2[0] = new TH2F("h2_0","QTCTime+F7T-F8T",200,0,2e4,400,0.8e3,1.2e3);
  h2[1] = new TH2F("h2_1","QTCTime-F7T+F8T",200,0,2e4,900,-1.e2,1e2);
  h2[2] = new TH2F("h2_2","QTCTime+F8T-F7T",200,0,2e4,500,-1.4e3,-0.9e3);
  h2[3] = new TH2F("h2_3","QTCTime-F8T+F7T",200,0,2e4,600,-0.3e3,0.3e3);/*  
  h2[4] = new TH2F("h2_4","QTCTime-F7T",100,0,1e4,1000,-1.9e3,-0.9e3);
  h2[5] = new TH2F("h2_5","QTCTime+F7T",100,0,1e4,1000,-1.0e3,0e3);
  h2[6] = new TH2F("h2_6","QTCTime+F8T",100,0,1e4,1000,-1.7e3,-700);
  h2[7] = new TH2F("h2_7","QTCTime-F8T",100,0,1e4,1000,-1.2e3,-200);  
									 */
  ifstream infile;
  infile.open("cutParameters/QTC_time.txt");
  Float_t corr[13];
  short ii = 0;
  while(infile.good()){
    std::string line, para;
    getline(infile, line);
    std::istringstream iss(line);
    iss >> para;
    corr[ii] = atof(para.c_str());
    cout<< corr[ii] <<endl;
    ii++;   
    if(ii==13)break;
  }
  ifstream infile2;
  infile2.open("cutParameters/QTC_time_pol.txt");
  Double_t pol[8];
  ii=0;
  while(infile2.good()){
    std::string line, para;
    getline(infile2, line);
    std::istringstream iss(line);
    iss >> para;
    pol[ii] =atof(para.c_str());
    cout<< pol[ii] <<endl;
    ii++;   
    if(ii==8)break;
  }
  Char_t*corrpol = Form("%3.30f+%3.30f*MIDGAIN[%i]+%3.30f*MIDGAIN[%i]*MIDGAIN[%i]+%3.30f*MIDGAIN[%i]*MIDGAIN[%i]*MIDGAIN[%i]+%3.30f*MIDGAIN[%i]*MIDGAIN[%i]*MIDGAIN[%i]*MIDGAIN[%i]+%3.30f*MIDGAIN[%i]*MIDGAIN[%i]*MIDGAIN[%i]*MIDGAIN[%i]*MIDGAIN[%i]+%3.30f*MIDGAIN[%i]*MIDGAIN[%i]*MIDGAIN[%i]*MIDGAIN[%i]*MIDGAIN[%i]*MIDGAIN[%i]+%3.34f*MIDGAIN[%i]*MIDGAIN[%i]*MIDGAIN[%i]*MIDGAIN[%i]*MIDGAIN[%i]*MIDGAIN[%i]*MIDGAIN[%i]",pol[0],pol[1],crystal,pol[2],crystal,crystal,pol[3],crystal,crystal,crystal,pol[4],crystal,crystal,crystal,crystal,pol[5],crystal,crystal,crystal,crystal,crystal,pol[6],crystal,crystal,crystal,crystal,crystal,crystal,pol[7],crystal,crystal,crystal,crystal,crystal,crystal,crystal);
  cout<< corrpol <<endl;
  Char_t* correction = Form("((%f*exp(%f*(MIDGAIN[%i]-%f)))+(%f*exp(%f*(MIDGAIN[%i]-%f))))+(%f*exp(%f*(MIDGAIN[%i]-%f)))-%f-%f*MIDGAIN[%i]-%f/(MIDGAIN[%i]-%f)",corr[0],corr[1],crystal,corr[2],corr[3],corr[4],crystal,corr[5],corr[6],corr[7],crystal,corr[8],corr[9],corr[10],crystal,corr[11],crystal,corr[12]);
  can[0]->cd(1);
  chain->Draw(Form("QTCTimeC:MIDGAIN>>h2_0",crystal, crystal),Form("TMath::Abs(AOQ[2]-2.64)<0.02 && TMath::Abs(ZET[2]-50)<0.5 && TMath::Abs(ZET[0]-50) < 0.5 && TMath::Abs(AOQ[0]-2.64)<0.01",trig),"colz");
  can[0]->cd(2);/*
  chain->Draw(Form("QTCTime[%i]-F7T+F8T-%f/BETA[0]-%s:MIDGAIN[%i]>>h2_1",crystal, length, correction, crystal),"TMath::Abs(AOQ[3]-2.635)<0.02 && TMath::Abs(ZET[3]-50)<0.5 && TMath::Abs(ZET[0]-50) < 0.5 && TMath::Abs(AOQ[0]-2.64)<0.01","colz");
  
  can[1]->cd(1);
  chain->Draw(Form("QTCTime[%i]-F8T+F7T-%f/BETA[0]:MIDGAIN[%i]>>h2_2",crystal,length, crystal),Form("TMath::Abs(AOQ[3]-2.635)<0.02 && TMath::Abs(ZET[3]-50)<0.5 && TMath::Abs(ZET[0]-50) < 0.5 && TMath::Abs(AOQ[0]-2.64)<0.01 && triggerbit == %i",trig),"colz");
  can[1]->cd(2);
  chain->Draw(Form("QTCTime[%i]-F8T+F7T-%f/BETA[0]-%s:MIDGAIN[%i]>>h2_3",crystal,length, corrpol, crystal),"TMath::Abs(AOQ[3]-2.635)<0.02 && TMath::Abs(ZET[3]-50)<0.5 && TMath::Abs(ZET[0]-50) < 0.5 && TMath::Abs(AOQ[0]-2.64)<0.01","colz");
  /*
  can[2]->cd(1);
  chain->Draw(Form("QTCTime[%i]-F7T:MIDGAIN[%i]>>h2_4",crystal, crystal),Form("TMath::Abs(AOQ[3]-2.635)<0.05 && TMath::Abs(ZET[3]-50)<0.5 && TMath::Abs(ZET[0]-50) < 0.5 && TMath::Abs(AOQ[0]-2.64)<0.01 && triggerbit == %i",trig),"colz");
  can[2]->cd(2);
  chain->Draw(Form("QTCTime[%i]+F7T:MIDGAIN[%i]>>h2_5",crystal, crystal),Form("TMath::Abs(AOQ[3]-2.635)<0.05 && TMath::Abs(ZET[3]-50)<0.5 && TMath::Abs(ZET[0]-50) < 0.5 && TMath::Abs(AOQ[0]-2.64)<0.01 && triggerbit == %i",trig),"colz");

  can[3]->cd(1);
  chain->Draw(Form("QTCTime[%i]-F8T:MIDGAIN[%i]>>h2_6",crystal, crystal),Form("TMath::Abs(AOQ[3]-2.635)<0.05 && TMath::Abs(ZET[3]-50)<0.5 && TMath::Abs(ZET[0]-50) < 0.5 && TMath::Abs(AOQ[0]-2.64)<0.01 && triggerbit == %i",trig),"colz");
  can[3]->cd(2);
  chain->Draw(Form("QTCTime[%i]+F8T:MIDGAIN[%i]>>h2_7",crystal, crystal),Form("TMath::Abs(AOQ[3]-2.635)<0.05 && TMath::Abs(ZET[3]-50)<0.5 && TMath::Abs(ZET[0]-50) < 0.5 && TMath::Abs(AOQ[0]-2.64)<0.01 && triggerbit == %i",trig),"colz");
  */

  TObjArray* dummy= new TObjArray();
  dummy->SetOwner(kTRUE);

  TH1F* h1_mean[3];
  TH1F* h1_sigma[3];
  short ican = 0, ipad = 1;
  /*  
  for(short ii = 0; ii < 3; ii++){
    TF1* f1 = new TF1(Form("f1_%i",ii),"landau",-1200,-900);
    h2[ii]->FitSlicesY(f1,0,-1,0,"QNRG4", dummy);

    h1_mean[ii]  = (TH1F*)dummy->FindObject(Form("h2_%i_1",ii));
    h1_sigma[ii] = (TH1F*)dummy->FindObject(Form("h2_%i_2",ii));
    
      can[ican]->cd(ipad);
      h1_mean[ii]->SetLineColor(ii+1);
      h1_mean[ii]->Draw("same");
    
      if(ii==0){
	can[4]->cd(1);
	h1_mean[ii]->SetLineColor(ii+1);
	h1_mean[ii]->Draw();
	can[0]->cd(1);
	h1_mean[ii]->Draw();

	can[4]->cd(2);
	h1_sigma[ii]->SetLineColor(ii+1);
	h1_sigma[ii]->Draw();
      }else {
	can[4]->cd(1);
	h1_mean[ii]->SetLineColor(ii+1);
	h1_mean[ii]->Draw();
	can[0]->cd(1);
	h1_mean[ii]->Draw("same");

	can[4]->cd(2);
	h1_sigma[ii]->SetLineColor(ii+1);
	h1_sigma[ii]->Draw("same");
      }
      if(ipad%2==0){
	ican++;
	ipad = 0;
      }
      ipad++;
  }
  */
}

void test_time(short start, short end){

  ifstream infile;
  infile.open("cutParameters/QTC_time.txt");
  Float_t parameter[7];
  short ii = 0;
  while(infile.good()){
    std::string line, para;
    getline(infile, line);
    std::istringstream iss(line);
    iss >> para;
    parameter[ii] = atof(para.c_str());
    cout<< parameter[ii] <<endl;
    ii++;   
  }


  TChain* chain = new TChain("tree");
  Char_t* path = "data/rootfiles/new";
  for(short ii = start; ii < end; ii++)chain->AddFile(Form("%s/run%04d.root",path,ii));

  TCanvas* can[2];
  for(short ii = 0; ii<2; ii++){
    can[ii] = new TCanvas(Form("can_%i",ii),"");
    can[ii]->Divide(2,3);
  }

  TH2F* h2[6];
  h2[0] = new TH2F("h2_0","QTCTime-F7T+F8T",100,0,1e4,60, 0,60);
  h2[1] = new TH2F("h2_1","QTCTime+F7T-F8T",100,0,1e4,120,-60,60);

  h2[2] = new TH2F("h2_2","QTCTime-F7T+F8T-beta",100,0,1e4,60, 0,60);
  h2[3] = new TH2F("h2_3","QTCTime-F7T+F8T+beta",100,0,1e4,60, 0,60);
  h2[4] = new TH2F("h2_4","QTCTime+F7T-F8T-beta",100,0,1e4,120,-60,60);
  h2[5] = new TH2F("h2_5","QTCTime+F7T-F8T+beta",100,0,1e4,120,-60,60);


  Char_t* corr = Form("exp(%f+%f*MIDGAIN)+exp(%f+%f*MIDGAIN)+expo(%f+%f*MIDGAIN)+%f",parameter[0],parameter[1],parameter[2],parameter[3],parameter[4],parameter[5],parameter[6]);

  can[0]->cd(1);
  chain->Draw("QTCTimeC+BigRIPSPlastic.fTime[1]-F8T:MIDGAIN>>h2_0","TMath::Abs(AOQ[2]-2.64)<0.02 && TMath::Abs(ZET[3]-50)<0.5 && TMath::Abs(ZET[0]-50) < 0.5 && TMath::Abs(AOQ[0]-2.64)<0.01","colz");
  can[0]->cd(2);
  // chain->Draw(Form("QTCTimeC+BigRIPSPlastic.fTime[1]-F8T-exp(%f+%f*MIDGAIN)-exp(%f+%f*MIDGAIN)-exp(%f+%f*MIDGAIN)-%f:MIDGAIN>>h2_1",parameter[0],parameter[1],parameter[2],parameter[3],parameter[4],parameter[5],parameter[6]),"TMath::Abs(AOQ[2]-2.64)<0.02 && TMath::Abs(ZET[3]-50)<0.5 && TMath::Abs(ZET[0]-50) < 0.5 && TMath::Abs(AOQ[0]-2.64)<0.01","colz");
  chain->Draw(Form("QTCTimeC-F7T-F8T-%s:MIDGAIN>>h2_1",corr),"TMath::Abs(AOQ[2]-2.64)<0.02 && TMath::Abs(ZET[3]-50)<0.5 && TMath::Abs(ZET[0]-50) < 0.5 && TMath::Abs(AOQ[0]-2.64)<0.01","colz");

  chain->Draw(Form("QTCTimeC+7T+F8T:MIDGAIN>>h2_1",corr),"TMath::Abs(AOQ[2]-2.64)<0.02 && TMath::Abs(ZET[3]-50)<0.5 && TMath::Abs(ZET[0]-50) < 0.5 && TMath::Abs(AOQ[0]-2.64)<0.01","colz");


}

void time(short start, short end, short trig){


  TChain* chain = new TChain("tree");
  Char_t* path = "data/rootfiles/new";
  for(short ii = start; ii < end; ii++)chain->AddFile(Form("%s/run%04d.root",path,ii));

  TCanvas* can[3];
  for(short ii = 0; ii<3; ii++){
    can[ii] = new TCanvas(Form("can_%i",ii),"");
    can[ii]->Divide(2,1);
  }

  TH2F* h2[6];
  h2[0] = new TH2F("h2_0","QTCTime-F7T+F8T",100,0,1e4,200,-100,100);
  h2[1] = new TH2F("h2_1","QTCTime+F7T-F8T",100,0,1e4,200,-100,100);
  
  h2[2] = new TH2F("h2_3","QTCTime+F8T",100,0,1e4,1000,0,1e3);
  h2[3] = new TH2F("h2_4","QTCTime-F8T",100,0,1e4,1000,-1e3,0);

  h2[4] = new TH2F("h2_5","QTCTime-F7T",100,0,1e4,1000,-1e3,0);
  h2[5] = new TH2F("h2_6","QTCTime+F7T",100,0,1e4,1000,-0,1e3);
  
  //h2[2] = new TH2F("h2_2","QTCTime",100,0,1e4,2000,-1000,1000);
  //  h2[3] = new TH2F("h2_3","QTCTime",100,0,1e4,1000,-500,500);

  can[0]->cd(1);
  chain->Draw("QTCTimeC[0]+BigRIPSPlastic.fTime-F8T:MIDGAIN[0]>>h2_0",Form("BigRIPSPlastic.fpl == 7 && TMath::Abs(AOQ[2]-2.64)<0.02 && TMath::Abs(ZET[3]-50)<0.5 && TMath::Abs(ZET[0]-50) < 0.5 && TMath::Abs(AOQ[0]-2.64)<0.01 && triggerbit == %i",trig),"colz");
  can[0]->cd(2);
  chain->Draw("QTCTimeC[0]-BigRIPSPlastic.fTime+F8T:MIDGAIN[0]>>h2_1",Form("BigRIPSPlastic.fpl == 7 && TMath::Abs(AOQ[2]-2.64)<0.02 && TMath::Abs(ZET[3]-50)<0.5 && TMath::Abs(ZET[0]-50) < 0.5 && TMath::Abs(AOQ[0]-2.64)<0.01 && triggerbit == %i",trig),"colz");
  /*  
  can[1]->cd(1);
  chain->Draw("QTCTimeC+F8T:MIDGAIN>>h2_3","TMath::Abs(AOQ[2]-2.64)<0.05 && TMath::Abs(ZET[3]-50)<0.5 && TMath::Abs(ZET[0]-50) < 0.5 && TMath::Abs(AOQ[0]-2.64)<0.01","colz");
  can[1]->cd(2);
  chain->Draw("QTCTimeC-F8T:MIDGAIN>>h2_4","TMath::Abs(AOQ[2]-2.64)<0.05 && TMath::Abs(ZET[3]-50)<0.5 && TMath::Abs(ZET[0]-50) < 0.5 && TMath::Abs(AOQ[0]-2.64)<0.01","colz");
  
  can[2]->cd(1);
  chain->Draw("QTCTimeC-BigRIPSPlastic.fTime:MIDGAIN>>h2_5","TMath::Abs(AOQ[2]-2.64)<0.05 && TMath::Abs(ZET[3]-50)<0.5 && TMath::Abs(ZET[0]-50) < 0.5 && TMath::Abs(AOQ[0]-2.64)<0.01 && BigRIPSPlastic.fpl == 7","colz");
  can[2]->cd(2);
  chain->Draw("QTCTimeC+BigRIPSPlastic.fTime:MIDGAIN>>h2_6","TMath::Abs(AOQ[2]-2.64)<0.05 && TMath::Abs(ZET[3]-50)<0.5 && TMath::Abs(ZET[0]-50) < 0.5 && TMath::Abs(AOQ[0]-2.64)<0.01 && BigRIPSPlastic.fpl == 7","colz");
  /*
  can[1]->cd(1);
  chain->Draw("QTCTimeC[0]:MIDGAIN[0]>>h2_2","TMath::Abs(AOQ[2]-2.64)<0.02 && TMath::Abs(ZET[3]-50)<0.5 && TMath::Abs(ZET[0]-50) < 0.5 && TMath::Abs(AOQ[0]-2.64)<0.01","colz");
  */                                       //"TMath::Abs(AOQ[2]-2.64)<0.02 && TMath::Abs(ZET[3]-50)<0.5 && TMath::Abs(ZET[0]-50) < 0.5 && TMath::Abs(AOQ[0]-2.64)<0.01 && TMath::Abs(QTCTimeC)<50"
  /*  can[1]->cd(2);
  chain->Draw("QTCTimeC:MIDGAIN>>h2_3","TMath::Abs(AOQ[2]-2.64)<0.05 && TMath::Abs(ZET[3]-50)<0.5 && TMath::Abs(ZET[0]-50) < 0.5 && TMath::Abs(AOQ[0]-2.64)<0.01","colz");
  */
  TObjArray* dummy= new TObjArray();
  dummy->SetOwner(kTRUE);

  TH1F* h1_mean[3];
  TH1F* h1_sigma[3];
  for(short ii = 0; ii < 3; ii++){
    TF1* f1 = new TF1(Form("f1_%i",ii),"landau",-100,100);
    h2[ii]->FitSlicesY(f1,0,-1,0,"QNRG4", dummy);

    h1_mean[ii]  = (TH1F*)dummy->FindObject(Form("h2_%i_1",ii));
    h1_sigma[ii] = (TH1F*)dummy->FindObject(Form("h2_%i_2",ii));

    if(ii>0){
      can[2]->cd(1);
      h1_mean[ii]->SetLineColor(ii+1);
      h1_mean[ii]->Draw("same");
      can[2]->cd(2);
      h1_sigma[ii]->SetLineColor(ii+1);
      h1_sigma[ii]->Draw("same");
    }else {
      can[2]->cd(1);
      h1_mean[ii]->SetLineColor(ii+1);
      h1_mean[ii]->Draw();
      can[0]->cd(1);
      h1_mean[ii]->Draw("same");

      can[2]->cd(2);
      h1_sigma[ii]->SetLineColor(ii+1);
      h1_sigma[ii]->Draw();
    }
  }



}



void time_DALI(short start, short end,Int_t trig){


  TChain* chain = new TChain("tree");
  Char_t* path = "data/rootfiles/new";
  for(short ii = start; ii < end; ii++)chain->AddFile(Form("%s/run%04d.root",path,ii));

  TCanvas* can[3];
  for(short ii = 0; ii<3; ii++){
    can[ii] = new TCanvas(Form("can_%i",ii),"");
    can[ii]->Divide(2,1);
  }

  TH2F* h2[4];
  h2[0] = new TH2F("h2_0","DALI Time-F7T+F8T",100,0,1e4,80,-40,40);
  h2[1] = new TH2F("h2_1","DALI Time-F7T-F8T",100,0,1e4,80,-10,70);
  /*
  h2[2] = new TH2F("h2_3","DALI Time+F8T",100,0,1e4,1000,0,1e3);
  h2[3] = new TH2F("h2_4","DALI Time-F8T",100,0,1e4,1000,-1e3,0);

  h2[4] = new TH2F("h2_5","DALI Time-F7T",100,0,1e4,1000,-1e3,0);
  h2[5] = new TH2F("h2_6","DALI Time+F7T",100,0,1e4,1000,-0,1e3);
  */
  h2[2] = new TH2F("h2_2","DALI Time",100,0,1e4,500,-1000,0);
  h2[3] = new TH2F("h2_3","DALI Time",100,0,1e4,100,-1000,1000);

  can[0]->cd(1);
  chain->Draw("DALINaI.fTimeOffseted+F8T-30:DALINaI.fEnergy>>h2_0","TMath::Abs(AOQ[2]-2.64)<0.02 && TMath::Abs(ZET[3]-50)<0.5 && TMath::Abs(ZET[0]-50) < 0.5 && TMath::Abs(AOQ[0]-2.64)<0.01 && DALINaI.id < 124","colz");
  can[0]->cd(2);
  chain->Draw("DALINaI.fTimeOffseted-F8T+2*BigRIPSPlastic.fTime[1]:DALINaI.fEnergy>>h2_1","TMath::Abs(AOQ[2]-2.64)<0.02 && TMath::Abs(ZET[3]-50)<0.5 && TMath::Abs(ZET[0]-50) < 0.5 && TMath::Abs(AOQ[0]-2.64)<0.01 && DALINaI.id < 124", "colz");
  /*
  can[1]->cd(1);
  chain->Draw("DALI TimeC+F8T:MIDGAIN>>h2_3","TMath::Abs(AOQ[2]-2.64)<0.05 && TMath::Abs(ZET[3]-50)<0.5 && TMath::Abs(ZET[0]-50) < 0.5 && TMath::Abs(AOQ[0]-2.64)<0.01","colz");
  can[1]->cd(2);
  chain->Draw("DALI TimeC-F8T:MIDGAIN>>h2_4","TMath::Abs(AOQ[2]-2.64)<0.05 && TMath::Abs(ZET[3]-50)<0.5 && TMath::Abs(ZET[0]-50) < 0.5 && TMath::Abs(AOQ[0]-2.64)<0.01","colz");
  
  can[2]->cd(1);
  chain->Draw("DALI TimeC-BigRIPSPlastic.fTime:MIDGAIN>>h2_5","TMath::Abs(AOQ[2]-2.64)<0.05 && TMath::Abs(ZET[3]-50)<0.5 && TMath::Abs(ZET[0]-50) < 0.5 && TMath::Abs(AOQ[0]-2.64)<0.01 && BigRIPSPlastic.fpl == 7","colz");
  can[2]->cd(2);
  chain->Draw("DALI TimeC+BigRIPSPlastic.fTime:MIDGAIN>>h2_6","TMath::Abs(AOQ[2]-2.64)<0.05 && TMath::Abs(ZET[3]-50)<0.5 && TMath::Abs(ZET[0]-50) < 0.5 && TMath::Abs(AOQ[0]-2.64)<0.01 && BigRIPSPlastic.fpl == 7","colz");
  */
  can[1]->cd(1);
  chain->Draw("DALINaI.fTimeOffseted-F8T:DALINaI.fEnergy>>h2_2",Form("TMath::Abs(AOQ[2]-2.64)<0.02 && TMath::Abs(ZET[3]-50)<0.5 && TMath::Abs(ZET[0]-50) < 0.5 && TMath::Abs(AOQ[0]-2.64)<0.01 && triggerbit == %i && DALINaI.id < 124", trig),"colz");
  can[1]->cd(2);
  chain->Draw("DALINaI.fTimeOffseted:DALINaI.fEnergy>>h2_3",Form("TMath::Abs(AOQ[2]-2.64)<0.02 && TMath::Abs(ZET[3]-50)<0.5 && TMath::Abs(ZET[0]-50) < 0.5 && TMath::Abs(AOQ[0]-2.64)<0.01 && triggerbit == %i && DALINaI.id < 124",  trig),"colz");
  TObjArray* dummy= new TObjArray();
  dummy->SetOwner(kTRUE);

  TH1F* h1_mean[3];
  TH1F* h1_sigma[3];
  for(short ii = 0; ii < 3; ii++){
    h2[ii]->FitSlicesY(0,1,60,0,"QNRG4", dummy);

    h1_mean[ii]  = (TH1F*)dummy->FindObject(Form("h2_%i_1",ii));
    h1_sigma[ii] = (TH1F*)dummy->FindObject(Form("h2_%i_2",ii));

    if(ii>0){
      can[2]->cd(1);
      h1_mean[ii]->SetLineColor(ii+1);
      h1_mean[ii]->Draw("same");
      can[2]->cd(2);
      h1_sigma[ii]->SetLineColor(ii+1);
      h1_sigma[ii]->Draw("same");
    }else {
      can[2]->cd(1);
      h1_mean[ii]->SetLineColor(ii+1);
      h1_mean[ii]->Draw();
      can[0]->cd(1);
      h1_mean[ii]->Draw("same");

      can[2]->cd(2);
      h1_sigma[ii]->SetLineColor(ii+1);
      h1_sigma[ii]->Draw();
    }
  }

}
