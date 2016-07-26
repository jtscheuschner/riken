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

void Slope(){

  static const short MaxDALI = 86;
  static const short MaxQTC = 8;
  //define output
  TFile* output = new TFile("data/rootfiles/analysis/Slope.root","RECREATE");

  vector<TH2F*> h2_ADC;
  vector<TH2F*> h2_QTC;

  vector<TH2F*> h2_ADC_corrected;
  vector<TH2F*> h2_QTC_corrected;

  h2_ADC = vector<TH2F*>(145);
  h2_QTC = vector<TH2F*>(8);

  h2_ADC_corrected = vector<TH2F*>(145);
  h2_QTC_corrected = vector<TH2F*>(8);

  for(short ii=0; ii<145; ii++){
    h2_ADC[ii] = new TH2F(Form("h2_ADC_%i",ii), Form("ADC %i",ii), 300,-1200,-900,100,0,1e4);
    h2_ADC_corrected[ii] = new TH2F(Form("h2_ADC_corrected_%i",ii), Form("ADC_corrected %i",ii), 300,-1200,-900,100,0,1e4);
  }
  for(short ii=0; ii<8;ii++){
    h2_QTC[ii] = new TH2F(Form("h2_QTC_%i",ii),Form("QTC %i",ii), 200,-100,100,100,0,1e4);
    h2_QTC_corrected[ii] = new TH2F(Form("h2_QTC_corrected_%i",ii),Form("QTC_corrected %i",ii), 200,-100,100,100,0,1e4);
  }
  cout<< "bla" <<endl;
  for(short irun=9; irun<64; irun++){//Loop over all trees

    //Initialize tree
    //    Char_t* filename = Form("data/rootfiles/new/run%04d.root",irun);
    //    cout<< "back  " << filename <<endl;
    TChain* tree = new TChain("tree");
    tree->AddFile(Form("data/rootfiles/new/run%04d.root",irun));
    
    tree->Draw("DALINaI.id");
    Long64_t Entry = tree->GetEntriesFast();
    cout<< "run " << irun << " number of events " << Entry <<endl;
   
    //Set branches and leaves just need Energy and time and timeoffset
    Int_t fDALINaI_, fADC_ID[MaxDALI];
    Double_t 
      fADC_Energy[MaxDALI],
      fADC_Time[MaxDALI],
      fADC_TimeOffset[MaxDALI],
      fQTC_Time[MaxQTC],
      fQTC_TimeC[MaxQTC],
      fQTC_Energy[MaxQTC];
    TBranch* bDALINAI;
    TBranch* bDALINAI_id;
    TBranch* bDALINAI_fEnergy;
    TBranch* bDALINAI_fTime;
    TBranch* bDALINAI_fTimeOffseted;

    TBranch* bQTC_Time;
    TBranch* bQTC_fTimeC;
    TBranch* bQTC_MID;
    cout<< "blobb" <<endl;
    tree->SetBranchAddress("DALINaI",&fDALINaI_, &bDALINAI);
    tree->SetBranchAddress("DALINaI.id",fADC_ID, &bDALINAI_id);
    tree->SetBranchAddress("DALINaI.fEnergy",fADC_Energy, &bDALINAI_fEnergy);
    tree->SetBranchAddress("DALINaI.fTime",fADC_Time, &bDALINAI_fTime);
    tree->SetBranchAddress("DALINaI.fTimeOffseted",fADC_TimeOffset, &bDALINAI_fTimeOffseted);

    tree->SetBranchAddress("QTCTime",fQTC_Time, &bQTC_fTimeC);
    tree->SetBranchAddress("QTCTimeC",fQTC_TimeC, &bQTC_fTimeC);
    tree->SetBranchAddress("MIDGAIN", fQTC_Energy, &bQTC_MID);
    //tree->SetBranchAddress("MIDGAINCTA", fQTC_Energy_CTA);
    //cout<< "blub" <<endl;
   
    Long64_t bla=0, bla2=-1;
    for(Long64_t ievent=0;ievent<Entry;ievent++){//Loop over all events
      cout<< "blobb" <<endl;    
      bla2 = tree->LoadTree(ievent);
      
      //bla = tree->GetEntry(ievent);
      //cout<< "blA" <<endl;

      //fill histograms
      for(short iDALI=0; iDALI<fDALINaI_/*kMaxDALINaI*/; iDALI++){
	cout<< "bloe" <<endl;
	h2_ADC[fADC_ID[iDALI]]->Fill(fADC_Time[iDALI],fADC_Energy[iDALI]);
	h2_ADC_corrected[fADC_ID[iDALI]]->Fill(fADC_Time[iDALI],fADC_Energy[iDALI]);
	cout<< "bloe" <<endl;
      }//for(DALI)
      for(short icrystal=0; icrystal<MaxQTC; icrystal++){
	h2_QTC[icrystal]->Fill(fQTC_Time[icrystal],fQTC_Energy[icrystal]);
	h2_QTC_corrected[icrystal]->Fill(fQTC_Time[icrystal],fQTC_Energy[icrystal]);
      }//for(crystal)


    }//for(event)
    //delete tree
  }//Loop over all trees




}//Slope
