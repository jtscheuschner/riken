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

void compare_labr(){

  TFile* _file[8];
  UInt_t array[8]={124,125,126,127,141,142,143,144};

  TH2F* h2[8];
  TH1F* h1[8][8];
  TF1*  f1[8][8];
  TCanvas* can=new TCanvas("can");
  can->Divide(4,2);
  int i=0;
  for(int i=0; i<8;i++){
    _file[i]= new TFile(Form("data/Histo/Sn132/test/%i_12.0_3.0_130.0_2.0_1.5_0.8_2.0_2.0_-5.0.root",array[i]));
  //_file[i]= new TFile(Form("data/Histo/Sn132/test/%i_12.0_3.0_130.0_2.0_1.5_0.8_2.0_2.0_-5.0.root",0));
  //_file[i]= new TFile("data/Histo/Sn132/12.0_3.0_130.0_2.0_1.5_0.8_2.0_2.0_-5.0.root");
    h2[i] = (TH2F*) _file[i] ->Get("top/LaBr/h2_labr_EdoppCor_MID");
    h2[i]->SetName(Form("%i",array[i]));

    for(int j=0; j<8;j++){
      h1[i][j]=(TH1F*)h2[i]->ProjectionY(Form("%i_%i",array[i],j),j+1,j+1);
      h1[i][j]->SetLineColor(i+1);
      can->cd(j+1);
      if(i==0)h1[i][j]->Draw();
      else    h1[i][j]->Draw("same");
      f1[i][j] = new TF1(Form("%i_%i",i,j),"gaus(0)+pol1(3)",2e3,5e3);
      f1[i][j]->SetParameters(20,4e3,350,9,-1e-3);
      h1[i][j]->Fit(f1[i][j],"R+B");

    }
  }
  for(int j=0; j<8;j++)
    for(int i=0; i<8;i++)
      //    for(int j=0; j<8;j++)
      cout<< i << " " << j << " " << f1[i][j]->GetParameter(0) << " " << f1[i][j]->GetParameter(2) <<endl;
}
