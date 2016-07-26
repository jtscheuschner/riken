
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

void Improve_PID2(Char_t* path, Char_t* filename, short aoq, short var=0,
		  Float_t A_Lin = 0, Float_t A_Squ = 0,
		  Float_t X_Lin = 0, Float_t X_Squ = 0,
		  Float_t Y_Lin = 0, Float_t Y_Squ = 0,
		  Float_t D_Lin = 0, Float_t D_Squ = 0,

		  Float_t A_Lin1 = 0, Float_t A_Squ1 = 0,
		  Float_t X_Lin1 = 0, Float_t X_Squ1 = 0,
		  Float_t Y_Lin1 = 0, Float_t Y_Squ1 = 0,
		  Float_t D_Lin1 = 0, Float_t D_Squ1 = 0,

		  Float_t B_Lin = 0, Float_t B_Squ = 0
){



  TChain* tree = new TChain("tree");
  std::cout<< Form("%s%s",path,filename) <<std::endl;
  tree->AddFile(Form("%s%s",path,filename));


  short beta_;
  vector<Char_t*> names;
  names = vector<Char_t*> (9);
  Float_t corr[2][9] = { {A_Lin, X_Lin, Y_Lin,D_Lin, A_Lin1, X_Lin1, Y_Lin1, D_Lin, B_Lin},
		      {A_Squ, X_Squ, Y_Squ,D_Squ, A_Squ1, X_Squ1, Y_Squ1, D_Squ, B_Squ} };
  Float_t plane_[2], delta_[2];


  if(aoq == 0){
    names[0]="FA[0]",  names[1]="FX[0]",  names[2]="FY[0]",  names[3]="delta[0]";//First detector
    names[4]="FA[1]",  names[5]="FX[1]",  names[6]="FY[1]",  names[7]="delta[1]";//Second detector
    names[8]="BETA[0]";                                                       //Beta before ([0]) or after ([1]) the target
    beta_ = 0;
    plane_[0] = 0;
    plane_[1] = 1;
    delta_[0] = 0;
    delta_[1] = 1;
  }else if(aoq == 2){
    names[0]="FA[4]",  names[1]="FX[4]",  names[2]="FY[4]",  names[3]="delta[2]";//First detector
    names[4]="FA[5]",  names[5]="FX[5]",  names[6]="FY[5]",  names[7]="delta[4]";//Second detector
    names[8]="BETA[1]";                                                       //Beta before ([0]) or afte
    beta_ = 0;
    plane_[0] = 4;
    plane_[1] = 5;
    delta_[0] = 2;
    delta_[1] = 4;
  }else if(aoq == 5){
    names[0]="FA[3]",  names[1]="FX[3]",  names[2]="FY[3]",  names[3]="delta[2]";//First detector
    names[4]="FA[5]",  names[5]="FX[5]",  names[6]="FY[5]",  names[7]="delta[4]";//Second detector
    names[8]="BETA[1]";                                                      //Beta before ([0]) or after
    beta_ = 0;
    plane_[0] = 3;
    plane_[1] = 5;
    delta_[0] = 2;
    delta_[1] = 4;
  }

  Char_t* LinCorr = "",* SquCorr = "",*cut = "";
  if(aoq != 0) cut = "ZET[0]>49.5 && ZET[0] < 50.5 && TMath::Abs(AOQ[0]-2.56)<0.01";
  for(short ii = 0; ii < 9; ii++){
    LinCorr = Form("%s+%f*%s"   ,LinCorr,corr[0][ii],names[ii]);
    SquCorr = Form("%s+%f*%s*%s",SquCorr,corr[1][ii],names[ii],names[ii]);
  }
  TH2F* h2_aoq;
  if(var == 8) h2_aoq = new TH2F("h2_aoq",Form("AOQ[%i] vs %s",aoq,names[var]),200, 0.5,0.6,200,2.54,2.58);
  else if(var == 3 || var == 7)h2_aoq = new TH2F("h2_aoq",Form("AOQ[%i] vs %s",aoq,names[var]),200,-2,2, 200,2.54,2.58);
  else  h2_aoq = new TH2F("h2_aoq",Form("AOQ[%i] vs %s",aoq, names[var]),200,-10,10,200,2.54,2.58);
  tree->Draw(Form("AOQ[%i]+%s+%s:%s>>h2_aoq",aoq,LinCorr,SquCorr,names[var]),cut,"colz");



  tree->Draw(Form("AOQ[2]+%s+%s:ZET[3]>>h50(100,45,55,200,2.5,2.7)",LinCorr,SquCorr),cut,"colz");

}
