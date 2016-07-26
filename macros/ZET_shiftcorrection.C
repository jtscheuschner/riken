#include <iostream>
#include <string>
#include <time.h>

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
#include "TCanvas.h"
/*
#include "TArtBigRIPSParameters.hh"
#include "TArtDALIParameters.hh"
#include "TArtCalibPID.hh"
#include "TArtCalibDALI.hh"
#include "TArtCalibPPAC.hh"
#include "TArtCalibPlastic.hh"
#include "TArtCalibFocalPlane.hh"
#include "TArtCalibCoin.hh"
#include "TArtEventInfo.hh"
#include "TArtPlastic.hh"
#include "TArtPPAC.hh"
#include "TArtFocalPlane.hh"

#include "TArtRecoPID.hh"
#include "TArtRecoRIPS.hh"
#include "TArtRecoTOF.hh"
#include "TArtRecoBeam.hh"
#include "TArtRecoPID.hh"
#include "TArtRecoRIPS.hh"
#include "TArtRecoTOF.hh"
#include "TArtRecoBeam.hh"
#include "TArtRIPS.hh"
#include "TArtTOF.hh"
#include "TArtBeam.hh"
#include "TArtIC.hh"
*/
#include "TVector3.h"
#include "TMath.h"

#include "signal.h"
#include "TSpectrum.h"
#include "TF1.h"

void ZET_shiftcorrection(short start = 9, short end = 63){

  FILE *ZETfile;
  if(start < 9)        ZETfile = fopen("cutParameters/PID/empty/peakpos.txt","w");
  else if(start < 100 )ZETfile = fopen("cutParameters/PID/Sn132/peakpos.txt","w");
  else                 ZETfile = fopen("cutParameters/PID/Sn128/peakpos.txt","w");
  TH1F* h1_ESq[63];
  TCanvas* can[4];
  for(short ii =0; ii<4;ii++){
    can[ii]=new TCanvas(Form("can_%i",ii),"can");
    can[ii]->Divide(4,4);
  }
  short ican=0, ipad=1;
  for(short run = start; run<end+1;run++){
    //if(run==118)continue;
   
    Char_t* filename;
    if( run < 9 )        filename = Form("data/rootfiles/test/run%04d_299.798_-159.333.root",run);
    else if( run < 100 ) filename = Form("data/rootfiles/test/run%04d_299.730_-158.540.root",run);
    else                 filename = Form("data/rootfiles/test/run%04d_299.778_-158.436.root",run);
    TChain *tree = new TChain("tree");
    tree->AddFile(filename);
    
    if(run<9)h1_ESq[run] = new TH1F(Form("h1_ESq_%04d",run),Form("ESq_%04d",run),200,600,1000);
    else h1_ESq[run] = new TH1F(Form("h1_ESq_%04d",run),Form("ESq_%04d",run),200,900,1300);
    can[ican]->cd(ipad);
    if(run < 100)
      tree->Draw(Form("ESQ>>h1_ESq_%04d",run),"TMath::Abs(BigRIPSBeam.aoq[4]-2.64)<0.01 && TMath::Abs(BigRIPSBeam.aoq[0]-2.64)<0.005 && TMath::Abs(BigRIPSBeam.zet[0]-50)<0.3","");
    else
      tree->Draw(Form("ESQ>>h1_ESq_%04d",run),"BigRIPSBeam.aoq[0]>2.555 && BigRIPSBeam.aoq[0]<2.565 && BigRIPSBeam.aoq[4]>2.555 && BigRIPSBeam.aoq[4]<2.565 && TMath::Abs(BigRIPSBeam.zet[0]-50)<0.3","");// && BigRIPSBeam.aoq[4]>2.555 && BigRIPSBeam.aoq[2]<2.565

    TSpectrum *s = new TSpectrum(2);
    Int_t npeaks =  s->Search(h1_ESq[run],1,"new");
    Float_t *xpeaks = s->GetPositionX();
    for(Int_t ii=0; ii<npeaks; ii++){
      TF1* f1_gaus = new TF1("f1_gaus","gaus",800,1300);
      h1_ESq[run]->Fit(f1_gaus);
      //cout<< "peak-position  "<< xpeaks[ii] <<endl;
      fprintf(ZETfile,"%i %f\n", run,f1_gaus->GetParameter(1));

    }

    if(ipad%16==0){
      ican++;
      ipad=1;
    }else ipad++;
    delete tree;
  }//for

  fclose(ZETfile);
  std::cout<< "out" <<endl;

  /*
    TF1* gaussian = new TF1("gaussian","gaus",900,1100);
    gaussian->SetParameters(
    h1_ESq->Fit("gaussian");
  */
  return;

}
