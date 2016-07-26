

#include <iostream>
#include <string>
#include <time.h>
#include <vector>

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
#include "TH2F.h"
#include "TF1.h"
#include "TCanvas.h"

#include "TVector3.h"
#include "TMath.h"

#include "signal.h"
#include "TSpectrum.h"


void Improve_PID_ZDS(short run, Float_t aoq){

  //Read file
  Char_t* path = "data/rootfiles/new/";
  TChain* tree = new TChain("tree");
  std::cout<< Form("%srun%04d.root",path,run) <<std::endl;
  tree->AddFile(Form("%srun%04d.root",path,run));

  //define cuts
  Char_t* cut_in_delta = "TMath::Abs(AOQ[0]-aoq)<0.002 && TMath::Abs(ZET[0]-50)<0.3 && TMath::Abs(delta[2]-delta[4])<0.5"; //secure incoming Sn and outgoing delta cut
  Char_t* aoq_out = "TMath::Abs(AOQ[2]-aoq)<0.01";
  Char_t* zet_out = "TMath::Abs(ZET[3]-50) <0.5";

  //draw aoq against Angle, X- & Y-Position, beta, delta(?)
  tree->Draw("AOQ[2]:F8A>>h2_aoq_F8A",Form("%s && %s", cut_in_delta,zet_out),"colz");

  //Fit projection

  //subtract fit

  //next variable

  //start over again till convergence

}
