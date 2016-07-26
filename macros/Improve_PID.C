
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

Double_t FindMean(TH1F* h1){//This function finds the mean value of the distribution
  TF1* f1 = new TF1("f1","gaus",-10,10);
  f1->SetParameter(1,2.);
  h1->Fit(f1,"N","",-10,10);
  return f1->GetParameter(1);
}

void Improve_PID_c(short max, Char_t* filename, short line = 2, Float_t AOQ=2.56, Float_t var = 0, Double_t corrLin[9]=0,Double_t corrQu[9]=0){//max - number of iteration, filename - obvious, line - which AOQ you use (0,2,5), aoq - the aoq value of your isotop
  //  short max = 2;
  /*
    This function looks at the variables 

    - angle
    - x-position
    - y-position
    - delta
    - beta

    to improve the AOQ resolution as this is a work in progress:
    it is not tested, debugged, and the result has to be checked. 
    The variables are supposed to be independent (this  should not be!!)
    The improvement are linear and afterwards quadratically

  */

  //Read Tree

  std::cout<<std::endl;
  std::cout<< "/******** Improve PID c "<<std::endl;
  std::cout<< "          this will improve the PID for AOQ[c]!  ********/" <<std::endl;
  std::cout<<std::endl;
  std::cout<< "----------------------------------------------------------" <<std::endl;
  std::cout<< " Did you adapt the variables to your purpose?? "<<std::endl;
  std::cout<<std::endl;

  Char_t//* filename = "run0102_301.400_-156.264.root",
    *path = "data/rootfiles/test/", *cuts;
  TChain* tree = new TChain("tree");
  std::cout<< Form("%s%s",path,filename) <<std::endl;
  tree->AddFile(Form("%s%s",path,filename));

  static const short variables = 9;// the number of variables you want to look at

  /*********************** 
	       the variables after this point to the termination comment have to be adapted corresponding to the isotop and particle identification channel
  ***********************/
  Int_t aoq = line;
  Float_t AOQ_value = AOQ;  //AOQ value for the isotop, it should before done that the desired value is close with choosing the correct timeoffset, this can be changed later again that the value is met again after the corrections

  std::vector<Char_t*> names;
  names = std::vector<Char_t*> (variables);
  Double_t correction[2][variables]={0x0};  //corrections to the AOQ linear(0) and quadratically(1) for the different variables
  for(short ii=0; ii < variables; ii++){
    correction[0][ii] = corrLin[ii];
    correction[1][ii] = corrQu[ii];
  }
  Double_t mean_value[variables]={0x0};     //Those corrections are needed because of the quadratic adaption

  if(aoq==0){
    names[0]="F5A",  names[1]="F5X",  names[2]="BigRIPSPPAC.fY",  names[3]="delta[1]";//First detector
    names[4]="F5A",  names[5]="F5X",  names[6]="BigRIPSPPAC.fY",  names[7]="delta[1]";//Second detector
    names[8]="BETA[0]";                                                       //Beta before ([0]) or after ([1]) the target
  /*
  for(short ii=0;ii<variables;ii++){//Find first the mean-value of each variable
    TH1F* h1 = new TH1F("h1","",100,-10,10);
    tree->Draw(Form("%s>>h1",names[ii]));
    mean_value[ii] = FindMean(h1);
    delete h1;
  }
  */
    cuts = /*"TMath::Abs(AOQ[0]-%f)<0.1";*/Form("TMath::Abs(%s-%f)<5 && TMath::Abs(%s-%f)<5 && TMath::Abs(%s-%f)<5 && TMath::Abs(%s-%f)<5 && TMath::Abs(%s-%f)<5 && TMath::Abs(%s-%f)<5 && TMath::Abs(%s-%f)<5 && TMath::Abs(%s-%f)<5 && TMath::Abs(%s-%f)<5 && TMath::Abs(AOQ[0]-%f)<0.1 && BigRIPSPPAC.fpl==5", 
							names[0],mean_value[0],names[1],mean_value[1],names[2],mean_value[2],
							names[3],mean_value[3],names[4],mean_value[4],names[5],mean_value[5],
							names[6],mean_value[6],names[7],mean_value[7],names[8],mean_value[8],AOQ_value);
  }else if(aoq==2){
    names[0]="F10A",  names[1]="F10X",  names[2]="F10Y",  names[3]="delta[2]";//First detector
    names[4]="F11A",  names[5]="F11X",  names[6]="F11Y",  names[7]="delta[4]";//Second detector
    names[8]="BETA[1]";                                                       //Beta before ([0]) or after ([1]) the target
  /*
  for(short ii=0;ii<variables;ii++){//Find first the mean-value of each variable
    TH1F* h1 = new TH1F("h1","",100,-10,10);
    tree->Draw(Form("%s>>h1",names[ii]));
    mean_value[ii] = FindMean(h1);
    delete h1;
  }
  */
    /*
    cuts = Form("TMath::Abs(%s-%f)<5 && TMath::Abs(%s-%f)<5 && TMath::Abs(%s-%f)<5 && TMath::Abs(%s-%f)<5 && TMath::Abs(%s-%f)<5 && TMath::Abs(%s-%f)<5 && TMath::Abs(%s-%f)<5 && TMath::Abs(%s-%f)<5 && TMath::Abs(%s-%f)<5 && TMath::Abs(AOQ[0]-%f)<0.1 && TMath::Abs(ZET[2]-53.7)<1", 
							names[0],mean_value[0],names[1],mean_value[1],names[2],mean_value[2],
							names[3],mean_value[3],names[4],mean_value[4],names[5],mean_value[5],
							names[6],mean_value[6],names[7],mean_value[7],names[8],mean_value[8],AOQ_value);
    */
    cuts = "TMath::Abs(AOQ[0]-2.56)<0.01 && TMath::Abs(AOQ[2]-2.56)<0.01 && TMath::Abs(ZET[2]-53.7)<1";
  }else {
    names[0]="F10A",  names[1]="F10X",  names[2]="F10Y",  names[3]="delta[2]";//First detector
    names[4]="F11A",  names[5]="F11X",  names[6]="F11Y",  names[7]="delta[4]";//Second detector
    names[8]="BETA[1]";                                                      //Beta before ([0]) or after ([1]) the target
  /*
  for(short ii=0;ii<variables;ii++){//Find first the mean-value of each variable
    TH1F* h1 = new TH1F("h1","",100,-10,10);
    tree->Draw(Form("%s>>h1",names[ii]));
    mean_value[ii] = FindMean(h1);
    delete h1;
  }
  */
    cuts = /*"TMath::Abs(AOQ[0]-%f)<0.1";*/Form("TMath::Abs(%s-%f)<5 && TMath::Abs(%s-%f)<5 && TMath::Abs(%s-%f)<5 && TMath::Abs(%s-%f)<5 && TMath::Abs(%s-%f)<5 && TMath::Abs(%s-%f)<5 && TMath::Abs(%s-%f)<5 && TMath::Abs(%s-%f)<5 && TMath::Abs(%s-%f)<5 && TMath::Abs(AOQ[0]-%f)<0.1 && TMath::Abs(ZET[2]-53.7)<1", 
							names[0],mean_value[0],names[1],mean_value[1],names[2],mean_value[2],
							names[3],mean_value[3],names[4],mean_value[4],names[5],mean_value[5],
							names[6],mean_value[6],names[7],mean_value[7],names[8],mean_value[8],AOQ_value);
	
  }
  TCanvas* can = new TCanvas("can","");
  tree->Draw("AOQ[2]:ZET[2]>>h2_bla(100,45,55,200,2.5,2.7)",cuts,"colz");
  
  TCanvas* can1 = new TCanvas("can1","variables");
  TCanvas* can2 = new TCanvas("can2","slices");
  TCanvas* can0 = new TCanvas("can0","AOQ");
  can1->Divide(3,3);
  can2->Divide(3,3);
  /*********************** 
	this should it be the rest goes automatically, besides the saving (should be done soon)
  ***********************/
  
  TObjArray* dummy = new TObjArray();//this is needed to delete properly the fitslices histograms
  dummy->SetOwner(kTRUE);
  

  for(short iteration = 0; iteration<max; iteration++){//how often do you want to look on all variables after the corrections have been changed?
    TH2F* h2[variables];
    TH1F* h2_Sl[variables];

    for(short ivar = 0; ivar<variables; ivar++){
      if(ivar!=var)continue;
      //if(ivar < 4)continue;
      Double_t correction_local[2][variables]={0x0};
      for(short jj = 0; jj<variables;jj++){
	for(short kk = 0; kk<2;kk++){
	  correction_local[kk][jj] = correction[kk][jj];
	  if(jj==ivar)correction_local[kk][jj]=0;
	}
      }

      if(ivar == 8 )                  h2[ivar] = new TH2F(Form("h2_%s",names[ivar]),names[ivar],200,0.525,0.535,100,2.55,2.57);//Beta
      else if(ivar == 7 || ivar == 3) h2[ivar] = new TH2F(Form("h2_%s",names[ivar]),names[ivar],100,-2,2,100,2.55,2.57);//delta
      else                            h2[ivar] = new TH2F(Form("h2_%s",names[ivar]),names[ivar],100,-10,10,100,2.55,2.57);//rest
      Char_t* Lin_Corr = Form("-%f*(%s-%f)-%f*(%s-%f)-%f*(%s-%f)-%f*(%s-%f)-%f*(%s-%f)-%f*(%s-%f)-%f*(%s-%f)-%f*(%s-%f)",
			      correction_local[0][0],names[0],mean_value[0],correction_local[0][1],names[1],mean_value[1],//First detector 
			      correction_local[0][2],names[2],mean_value[2],correction_local[0][3],names[3],mean_value[3],//First detector 
			      correction_local[0][4],names[4],mean_value[4],correction_local[0][5],names[5],mean_value[5],//Second detector
			      correction_local[0][6],names[6],mean_value[6],correction_local[0][7],names[7],mean_value[7],//Second Detector
			      correction_local[0][8],names[8],mean_value[8]);

      Char_t* Squ_Corr = Form("-%f*(%s-%f)*(%s-%f)-%f*(%s-%f)*(%s-%f)-%f*(%s-%f)*(%s-%f)-%f*(%s-%f)*(%s-%f)-%f*(%s-%f)*(%s-%f)-%f*(%s-%f)*(%s-%f)-%f*(%s-%f)*(%s-%f)-%f*(%s-%f)*(%s-%f)",
			      correction_local[1][0],names[0],mean_value[0],names[0],mean_value[0], correction_local[1][1],names[1],mean_value[1],names[1],mean_value[1],//First detector
			      correction_local[1][2],names[2],mean_value[2],names[2],mean_value[2], correction_local[1][3],names[3],mean_value[3],names[3],mean_value[3],//First detector
			      correction_local[1][4],names[4],mean_value[4],names[4],mean_value[4], correction_local[1][5],names[5],mean_value[5],names[5],mean_value[5],//Second detector
			      correction_local[1][6],names[6],mean_value[6],names[6],mean_value[6], correction_local[1][7],names[7],mean_value[7],names[7],mean_value[7],//Second Detector
			      correction_local[1][8],names[8],mean_value[8],names[8],mean_value[8]);
      //Drawing
      can1->cd(ivar+1);
      //std::cout<<Form("AOQ[%i]%s%s:%s>>h2_%s", line,Lin_Corr,Squ_Corr,names[ivar],names[ivar])<<std::endl;
      tree->Draw(Form("AOQ[%i]%s%s:%s>>h2_%s", line,Lin_Corr,Squ_Corr,names[ivar],names[ivar]),
		 cuts,"colz");
      //FitSlices
      TF1* gaus_A = new TF1("gaus_A","gaus",2.5,2.7);
      //gaus_A->SetParameters(AOQ_value, correction[ivar][0], correction[ivar][1]);
      h2[ivar]->FitSlicesY(gaus_A,0,-1,0,"QNRG3",dummy);
      //Drawing and Fitting 
      h2_Sl[ivar] = (TH1F*) dummy->FindObject(Form("h2_%s_1",names[ivar]));
      can2->cd(ivar+1);
      h2_Sl[ivar]->Draw();
      TF1* func = new TF1("func","pol2",-10,10);
      func->SetParameter(0,AOQ_value);
      func->SetParameter(1,correction[ivar][0]);
      //func->SetParLimits(1,-1e-3,1e-3);
      func->SetParameter(2,correction[ivar][1]);
      //func->SetParLimits(2,-1e-4,1e-4);
      h2_Sl[ivar]->Fit(func,"BM","",-10,10);//Q: quiet, B: user defined parameter settings, M: Improve fit results, N: Do not draw & do not store the graphics
      correction[ivar][0] = func->GetParameter(1);
      // if(TMath::Abs(correction[ivar][0]) == 1e-3)correction[ivar][0] =0;
      correction[ivar][1] = func->GetParameter(2);

    }//for variables

    //deleting objects which are not used keep those of the last iteration
    if( iteration+1 < max){
      /*for(short ivar=0;ivar<variables;ivar++)
	delete h2[ivar];*/
      dummy->Clear();
    }//if
  }//for iterations
      
  Char_t* Lin_Corr = Form("-%f*(%s-%f)-%f*(%s-%f)-%f*(%s-%f)-%f*(%s-%f)-%f*(%s-%f)-%f*(%s-%f)-%f*(%s-%f)-%f*(%s-%f)",
			  correction[0][0],names[0],mean_value[0],correction[1][0],names[1],mean_value[1],correction[2][0],names[2],mean_value[2],correction[3][0],names[3],mean_value[3],//First Detector
			  correction[4][0],names[4],mean_value[4],correction[5][0],names[5],mean_value[5],correction[6][0],names[6],mean_value[6],correction[7][0],names[7],mean_value[7],//Second Detector
			  correction[8][0],names[8],mean_value[8]);                                                                                                                       //Beta
  //  std::cout<< Lin_Corr <<std::endl;
  Char_t* Squ_Corr = Form("-%f*(%s-%f)*(%s-%f)-%f*(%s-%f)*(%s-%f)-%f*(%s-%f)*(%s-%f)-%f*(%s-%f)*(%s-%f)-%f*(%s-%f)*(%s-%f)-%f*(%s-%f)*(%s-%f)-%f*(%s-%f)*(%s-%f)-%f*(%s-%f)*(%s-%f)",
			  correction[0][0],names[0],mean_value[0],names[0],mean_value[0], correction[1][0],names[1],mean_value[1],names[1],mean_value[1],//First detector
			  correction[2][0],names[2],mean_value[2],names[2],mean_value[2], correction[3][0],names[3],mean_value[3],names[3],mean_value[3],//First detector
			  correction[4][0],names[4],mean_value[4],names[4],mean_value[4], correction[5][0],names[5],mean_value[5],names[5],mean_value[5],//Second detector
			  correction[6][0],names[6],mean_value[6],names[6],mean_value[6], correction[7][0],names[7],mean_value[7],names[7],mean_value[7],//Second Detector
			  correction[8][0],names[8],mean_value[8],names[8],mean_value[8]);
  //  std::cout<< Squ_Corr <<std::endl;
  /*
  cuts = Form("TMath::Abs(%s-%f)<10 && TMath::Abs(%s-%f)<10 && TMath::Abs(%s-%f)<10 && TMath::Abs(%s-%f)<10 && TMath::Abs(%s-%f)<10 && TMath::Abs(%s-%f)<10 && TMath::Abs(%s-%f)<10 && TMath::Abs(%s-%f)<10 && TMath::Abs(%s-%f)<5 && TMath::Abs(AOQ[0]-%f)<2.1", 
		      names[0],mean_value[0],names[1],mean_value[1],names[2],mean_value[2],
		      names[3],mean_value[3],names[4],mean_value[4],names[5],mean_value[5],
		      names[6],mean_value[6],names[7],mean_value[7],names[8],mean_value[8],AOQ_value);
  */
  TH1F*h1_AOQ = new TH1F("h1_AOQ","AOQ with corrections",500,2.5,2.7);
  can0->cd();
  //  std::cout<<Form("AOQ[%i]%s%s>>h1_AOQ", aoq,Lin_Corr,Squ_Corr) <<std::endl;
  tree->Draw(Form("AOQ[%i]%s%s>>h1_AOQ", aoq,Lin_Corr,Squ_Corr),"","colz");
  TF1* f1 = new TF1("f1","gaus(0)+pol1(3)",2.55,2.57);
  f1->SetParLimits(1,2.557,2.562);
  f1->SetParameter(2,3e-3);
  h1_AOQ->Fit(f1,"BM","",2.55,2.57);


}//Improve_c

void Improve_PID(short max, Char_t* filename){//AOQ[0]
  //  short max = 2;
  /*
    This function looks at the variables 

    - angle
    - x-position
    - y-position
    - delta
    - beta

    to improve the AOQ resolution as this is a work in progress:
    it is not tested, debugged, and the result has to be checked. 
    The variables are supposed to be independent (this  should not be!!)
    The improvement are linear and afterwards quadratically

  */

  //Read Tree
  Char_t//* filename = "run0102_301.400_-156.264.root",
    *path = "data/rootfiles/test/";
  TChain* tree = new TChain("tree");
  std::cout<< Form("%s%s",path,filename) <<std::endl;
  tree->AddFile(Form("%s%s",path,filename));

  //corrections to the AOQ linear(0) and quadratically(1) for the different variables
  Double_t angle[2]={0,0},     x_pos[2]={0,0},      y_pos[2]={0,0}, 		delta[2]={0,0},		  beta[2]={0,0};
  Float_t AOQ_value = 2.56;  
  //  TH2F* h2_AOQ = new TH2F("h2_AOQ","h2_AOQ",100,-10,10,100,2.55,2.57);
  //TH1F* h1_AOQ = new TH1F("h1_AOQ","h1_AOQ",100,2.55,2.57);
  //TH1F* h1_FitSlices_AOQ;

  TCanvas* can = new TCanvas("can", "variables");//for the 2d histo
  TCanvas* can_sl = new TCanvas("can_sl", "slices");//for the slices of the 2d histos
  TCanvas* can_aoq = new TCanvas("can_aoq", "AOQ");//in the end the new AOQ
  can->Divide(2,3);
  can_sl->Divide(2,3);
  TObjArray* dummy = new TObjArray();
  dummy->SetOwner(kTRUE);

  for(short iteration = 0; iteration<max; iteration++){//iterate to get best AOQ resolution

    //start with the angle
    std::cout<< "//****************** A **************// " <<std::endl;
    can->cd(1);
    TH2F* h2_AOQ_A = new TH2F("h2_AOQ_A","h2_AOQ_A",100,-10,10,100,2.55,2.57);
    tree->Draw(Form("AOQ[0]-BigRIPSPPAC.fX*%f-BigRIPSPPAC.fY*%f-delta[1]*%f-BETA[0]*%f-%f-BigRIPSPPAC.fX*BigRIPSPPAC.fX*%f-%f-BigRIPSPPAC.fY*BigRIPSPPAC.fY*%f-delta[1]*delta[1]*%f-BETA[0]*BETA[0]*%f:F5A>>h2_AOQ_A",
		    x_pos[0],y_pos[0],delta[0],beta[0],x_pos[1],y_pos[1],delta[1],beta[1]),
	       "BigRIPSPPAC.fpl==5","colz");
    std::cout<< " dummes A " << h2_AOQ_A <<std::endl;
    TF1* gaus_A = new TF1("gaus_A","gaus",2.5,2.7);
    gaus_A->SetParameters(AOQ_value, angle[0], angle[1]);
    h2_AOQ_A->FitSlicesY(gaus_A,0,-1,0,"QNRG3",dummy);
    TH1F* h2_AOQ_A_1 = (TH1F*) dummy->FindObject("h2_AOQ_A_1");
    std::cout<< " gans dummes A " << h2_AOQ_A_1 <<std::endl;
    can_sl->cd(1);
    h2_AOQ_A_1->Draw();
    TF1* funcangle = new TF1("funcangle","pol2",-10,10);
    funcangle->FixParameter(0,2.56);
    h2_AOQ_A_1->Fit(funcangle,"","",-10,10);//Q: quiet, B: user defined parameter settings, M: Improve fit results, N: Do not draw & do not store the graphics
    angle[0] = funcangle->GetParameter(1);
    angle[1] = funcangle->GetParameter(2);
    //break;
   
    //continue with X
    std::cout<< "//****************** X **************// " <<std::endl;
    can->cd(2);
    TH2F* h2_AOQ_X = new TH2F("h2_AOQ_X","h2_AOQ_X",100,-10,10,100,2.55,2.57);
    tree->Draw(Form("AOQ[0]-F5A*%f-BigRIPSPPAC.fY*%f-delta[1]*%f-BETA[0]*%f-%f-F5A*F5A*%f-%f-BigRIPSPPAC.fY*BigRIPSPPAC.fY*%f-delta[1]*delta[1]*%f-BETA[0]*BETA[0]*%f:BigRIPSPPAC.fX>>h2_AOQ_X",
		    angle[0],y_pos[0],delta[0],beta[0],angle[1],y_pos[1],delta[1],beta[1]),
	       "BigRIPSPPAC.fpl==5","colz"); 
    can_sl->cd(2);
    TF1* gaus_X = new TF1("gaus_X","gaus",2.5,2.7);
    gaus_X->SetParameters(AOQ_value, x_pos[0], x_pos[1]);
    h2_AOQ_X->FitSlicesY(gaus_X,0,-1,0,"QNRG3",dummy);
    TH1F* h2_AOQ_X_1 = (TH1F*) dummy->FindObject("h2_AOQ_X_1");
    h2_AOQ_X_1->Draw();
    TF1* funcX = new TF1("funcX","pol2",-10,10);
    funcX->FixParameter(0,2.56);
    h2_AOQ_X_1->Fit(funcX,"BM","",-10,10);
    x_pos[0] = funcX->GetParameter(1);
    x_pos[1] = funcX->GetParameter(2);
    //break;    
   
    //Y
    std::cout<< "//****************** Y **************// " <<std::endl;
    can->cd(3);
    TH2F* h2_AOQ_Y = new TH2F("h2_AOQ_Y","h2_AOQ_Y",100,-10,10,100,2.55,2.57);
    tree->Draw(Form("AOQ[0]-F5A*%f-BigRIPSPPAC.fX*%f-delta[1]*%f-BETA[0]*%f-%f-F5A*F5A*%f-BigRIPSPPAC.fX*BigRIPSPPAC.fX*%f-%f-delta[1]*delta[1]*%f-BETA[0]*BETA[0]*%f:BigRIPSPPAC.fY>>h2_AOQ_Y",
		    angle[0],x_pos[0],delta[0],beta[0],angle[1],x_pos[1],delta[1],beta[1]),
	       "BigRIPSPPAC.fpl==5","colz");
    can_sl->cd(3);
    TF1* gaus_Y = new TF1("gaus_Y","gaus",2.5,2.7);
    gaus_Y->SetParameters(AOQ_value, y_pos[0], y_pos[1]);
    h2_AOQ_Y->FitSlicesY(gaus_Y,0,-1,0,"QNRG3",dummy);
    TH1F* h2_AOQ_Y_1 = (TH1F*) dummy->FindObject("h2_AOQ_Y_1");
    h2_AOQ_Y_1->Draw();
    TF1* funcY = new TF1("funcY","pol2",-10,10);
    funcY->FixParameter(0,2.56);
    h2_AOQ_Y_1->Fit(funcY,"BM","",-10,10);
    y_pos[0] = funcY->GetParameter(1);
    y_pos[1] = funcY->GetParameter(2);
    // break;
    // at the moment not in use due to no further reasen as to be unsure if this is clever and it does not show significant dependence at least in run102 for Sn128 at F5
    //delta
    std::cout<< "//****************** D **************// " <<std::endl;
    can->cd(5);
    TH2F* h2_AOQ_D = new TH2F("h2_AOQ_D","h2_AOQ_D",100,-2,2,100,2.55,2.57);
    tree->Draw(Form("AOQ[0]-F5A*%f-BigRIPSPPAC.fX*%f-BigRIPSPPAC.fY*%f-BETA[0]*%f-%f-F5A*F5A*%f-BigRIPSPPAC.fX*BigRIPSPPAC.fX*%f-%f-BigRIPSPPAC.fY*BigRIPSPPAC.fY*%f-delta[1]*delta[1]*%f-BETA[0]*BETA[0]*%f:delta[0]>>h2_AOQ_D",angle[0],x_pos[0],y_pos[0],beta[0],angle[1],x_pos[1],y_pos[1],beta[1]),"BigRIPSPPAC.fpl==5","colz");
    h2_AOQ_D->FitSlicesY(0,0,-1,0,"QNRG3",dummy);

    can_sl->cd(5);
    TF1* funcdelta = new TF1("funcdelta","pol2",-10,10);
    funcdelta->FixParameter(0,2.56);
    TH1F* h2_AOQ_D_1 = (TH1F*) dummy->FindObject("h2_AOQ_D_1");
    h2_AOQ_D_1->Fit("pol2","QBMN","",-10,10);
    delta[0] = funcdelta->GetParameter(1);
    delta[1] = funcdelta->GetParameter(2);
    h2_AOQ_D_1->Draw();

    //BETA
    std::cout<< "//****************** B **************// " <<std::endl;
    can->cd(4);
    TH2F* h2_AOQ_B = new TH2F("h2_AOQ_B","h2_AOQ_B",100,0.54,0.6,100,2.55,2.57);
    tree->Draw(Form("AOQ[0]-F5A*%f-BigRIPSPPAC.fX*%f-BigRIPSPPAC.fY*%f-delta[1]*%f-F5A*F5A*%f-%f-BigRIPSPPAC.fX*BigRIPSPPAC.fX*%f-%f-BigRIPSPPAC.fY*BigRIPSPPAC.fY*%f-delta[1]*delta[1]*%f:BETA[0]>>h2_AOQ_B",
		    angle[0],x_pos[0],y_pos[0],delta[0],angle[1],x_pos[1],y_pos[1],delta[1]),
	       "BigRIPSPPAC.fpl==5","colz");
    TCanvas* can;
    can_sl->cd(4);
    TF1* gaus_B = new TF1("gaus_B","gaus",2.5,2.7);
    gaus_B->SetParameters(AOQ_value, beta[0], beta[1]);
    h2_AOQ_B->FitSlicesY(gaus_B,0,-1,0,"QNRG3",dummy);
    TH1F* h2_AOQ_B_1 = (TH1F*) dummy->FindObject("h2_AOQ_B_1");
    h2_AOQ_B_1->Draw();
    TF1* funcbeta = new TF1("funcbeta","pol2",-10,10);
    funcbeta->SetParameter(0,2.56);
    h2_AOQ_B_1->Fit("pol2","BM","",-10,10);
    beta[0] = funcbeta->GetParameter(1);
    beta[1] = funcbeta->GetParameter(2);

    //look at the AOQ resolution
    can_aoq->cd();
    TH1F* h1_AOQ_B = new TH1F("h1_AOQ_B","h1_AOQ_B",500,2.5,2.7);
    tree->Draw(Form("AOQ[0]-F5A*%f-BigRIPSPPAC.fX*%f-BigRIPSPPAC.fY*%f-delta[1]*%f-BETA[0]*%f-F5A*F5A*%f-%f-BigRIPSPPAC.fX*BigRIPSPPAC.fX*%f-%f-BigRIPSPPAC.fY*BigRIPSPPAC.fY*%f-delta[1]*delta[1]*%f-BETA[0]*BETA[0]*%f>>h1_AOQ_B",
		    angle[0],x_pos[0],y_pos[0],delta[0],beta[0],angle[1],x_pos[1],y_pos[1],delta[1],beta[1]),"BigRIPSPPAC.fpl==5","");


    TF1* funcaoq = new TF1("funcaoq","gaus(0)+pol1(3)",2.55,2.57);
    funcaoq->SetParLimits(0,1e3,1e7);
    funcaoq->SetParLimits(2,1e-5,8e-3);
    funcaoq->SetParLimits(1,2.555,2.565);
    h1_AOQ_B->Fit(funcaoq,"BM","",2.55,2.57);
    std::cout<< " --------------------- " << std::endl;
    std::cout<< "iteration: " << iteration+1 << "  sigma of AOQ: " << funcaoq->GetParameter(2) << std::endl;
    std::cout<<std::endl;
    //break;
    
    if( iteration+1 < max){

      delete funcangle;
      delete funcX;
      delete funcY;
      delete funcbeta;
      delete funcaoq;
      delete h2_AOQ_A;
      /*
      delete h2_AOQ_A_0;
      delete h2_AOQ_A_2;     
      delete h2_AOQ_A_chi2;    
      delete h2_AOQ_X_0;
      delete h2_AOQ_X_2;
      delete h2_AOQ_X_chi2;
      delete h2_AOQ_Y_0;
      delete h2_AOQ_Y_2;
      delete h2_AOQ_Y_chi2;
      delete h2_AOQ_B_0;
      delete h2_AOQ_B_2;
      delete h2_AOQ_B_chi2;
     
      delete h2_AOQ_A_1;
      delete h2_AOQ_X_1;
      delete h2_AOQ_Y_1;     
      delete h2_AOQ_B_1;   */
      delete h2_AOQ_X;
      delete h2_AOQ_Y;
      delete h2_AOQ_B;
      dummy->Clear();

    }
  }//for



}//improve PID


void Improve_PID_2(short max, Char_t* filename){//AOQ[2]
  //  short max = 2;
  /*
    This function looks at the variables 

    - angle
    - x-position
    - y-position
    - delta
    - beta

    to improve the AOQ resolution as this is a work in progress:
    it is not tested, debugged, and the result has to be checked. 
    The variables are supposed to be independent (this  should not be!!)
    The improvement are linear and afterwards quadratically

  */

  //Read Tree


  std::cout<< "/******** Improve PID 2 "<<std::endl;
  std::cout<< "          this will improve the PID for AOQ[2]!  ********/" <<std::endl;

  Char_t//* filename = "run0102_301.400_-156.264.root",
    *path = "data/rootfiles/test/";
  TChain* tree = new TChain("tree");
  std::cout<< Form("%s%s",path,filename) <<std::endl;
  tree->AddFile(Form("%s%s",path,filename));

  //corrections to the AOQ linear(0) and quadratically(1) for the different variables
  Double_t angle[2]={0,0},     x_pos[2]={0,0},      y_pos[2]={0,0}, 		delta[2]={0,0},		  beta[2]={0,0},
	    angle_11[2]={0,0},     x_pos_11[2]={0,0},      y_pos_11[2]={0,0},       delta2[2]={0,0}; 

	    /*********************** 
	       the variables after this point to the termination comment have to be adapted corresponding to the isotop and particle identification channel
	    ***********************/

  Float_t AOQ_value = 2.56;  
  static const short variables = 7;
 
  //F8
  TCanvas* can = new TCanvas("can", "variables - F10");//for the 2d histo
  TCanvas* can_sl = new TCanvas("can_sl", "slices - F10");//for the slices of the 2d histos

  TCanvas* can_11 = new TCanvas("can_11", "variable - F11");//for the 2d histo
  TCanvas* can_sl_11 = new TCanvas("can_sl_11", "slices - F11");//for the slices of the 2d histos

  TCanvas* can_aoq = new TCanvas("can_aoq", "AOQ");//in the end the new AOQ
  can->Divide(2,2);
  can_sl->Divide(2,2);

  can_11->Divide(2,2);
  can_sl_11->Divide(2,2);

  TObjArray* dummy = new TObjArray();
  dummy->SetOwner(kTRUE);
  if(max==0){//for testing
    can_aoq->cd();
    TH1F* h1_AOQ_B = new TH1F("h1_AOQ_B","h1_AOQ_B",600,2.4,2.7);
    tree->Draw(Form("AOQ[2]-F10A*%f-F10X*%f-F10Y*%f-delta[2]*%f-BETA[1]*%f-F10A*F10A*%f-F10X*F10X*%f-F10Y*F10Y*%f-delta[2]*delta[2]*%f-BETA[1]*BETA[1]*%f-F11A*%f-F11X*%f-F11Y*%f-delta[4]*%f-F11A*F11A*%f-F11X*F11X*%f-F11Y*F11Y*%f-delta[4]*delta[4]*%f>>h1_AOQ_B",
		    angle[0],x_pos[0],y_pos[0],delta[0],beta[0],  angle[1],x_pos[1],y_pos[1],delta[1],beta[1],      angle_11[0],x_pos_11[0],y_pos_11[0],delta2[0],   angle_11[1],x_pos_11[1],y_pos_11[1],delta2[1]),
	       "","colz");
  }
  /*
  Double_t mean_variables[variables]={0x0}
  for(short ii = 0;ii<variables;ii++){
    TH1F* h1=new TH1F("h1","h1",100,-10,10);
    tree->Draw(Form("%s>>h1"));
    mean_variables[ii] = FindMean(h1);
    delete h1;

  }
  */
  for(short iteration = 0; iteration<max; iteration++){//iterate to get best AOQ resolution
    /*
      Here starts a problem:

      The variables are not clustered around (0/AOQ[2]), but approximatly around (2/AOQ[2]). This offset has to be coorrected somehow.

      A solution would be to find out where the mean of the distribution is and then to adjust the variable to that. (before the loop for each variable and then to set the variable to the new value (everywhere))
    */


    //start with the angle
    std::cout<< "//****************** A-F10 **************// " <<std::endl;
    can->cd(1);
    TH2F* h2_AOQ_A = new TH2F("h2_AOQ_A","h2_AOQ_A",100,-10,10,100,2.55,2.57);
    tree->Draw(Form("AOQ[2]-F10X*%f-F10Y*%f-delta[2]*%f-BETA[1]*%f-F10X*F10X*%f-F10Y*F10Y*%f-delta[2]*delta[2]*%f-BETA[1]*BETA[1]*%f-F11A*%f-F11X*%f-F11Y*%f-delta[4]*%f-F11A*F11A*%f-F11X*F11X*%f-F11Y*F11Y*%f-delta[4]*delta[4]*%f:F10A>>h2_AOQ_A",
		    x_pos[0],y_pos[0],delta[0],beta[0],  x_pos[1],y_pos[1],delta[1],beta[1],      angle_11[0],x_pos_11[0],y_pos_11[0],delta2[0],   angle_11[1],x_pos_11[1],y_pos_11[1],delta2[1]),
	       "TMath::Abs(AOQ[2]-2.56)<0.01 && TMath::","colz");
    std::cout<< " dummes A " << h2_AOQ_A <<std::endl;
    TF1* gaus_A = new TF1("gaus_A","gaus",2.5,2.7);
    gaus_A->SetParameters(AOQ_value, angle[0], angle[1]);
    h2_AOQ_A->FitSlicesY(gaus_A,0,-1,0,"QNRG3",dummy);
    TH1F* h2_AOQ_A_1 = (TH1F*) dummy->FindObject("h2_AOQ_A_1");
    std::cout<< " gans dummes A " << h2_AOQ_A_1 <<std::endl;
    can_sl->cd(1);
    h2_AOQ_A_1->Draw();
    TF1* funcangle = new TF1("funcangle","pol2",-10,10);
    funcangle->FixParameter(0,2.56);
    h2_AOQ_A_1->Fit(funcangle,"","",-10,10);//Q: quiet, B: user defined parameter settings, M: Improve fit results, N: Do not draw & do not store the graphics
    angle[0] = funcangle->GetParameter(1);
    angle[1] = funcangle->GetParameter(2);
    //break;
   
    std::cout<< "//****************** A-F11 **************// " <<std::endl;
    can_11->cd(1);
    TH2F* h2_AOQ_11_A = new TH2F("h2_AOQ_11_A","h2_AOQ_11_A",100,-10,10,100,2.55,2.57);
    tree->Draw(Form("AOQ[2]-F10A*%f-F10X*%f-F10Y*%f-delta[2]*%f-BETA[1]*%f-F10A*F10A*%f-F10X*F10X*%f-F10Y*F10Y*%f-delta[2]*delta[2]*%f-BETA[1]*BETA[1]*%f-F11X*%f-F11Y*%f-delta[4]*%f-F11X*F11X*%f-F11Y*F11Y*%f-delta[4]*delta[4]*%f:F10A>>h2_AOQ_11_A",
		    angle[0],x_pos[0],y_pos[0],delta[0],beta[0],  angle[1],x_pos[1],y_pos[1],delta[1],beta[1],      x_pos_11[0],y_pos_11[0],delta2[0],   x_pos_11[1],y_pos_11[1],delta2[1]),
	       "TMath::Abs(AOQ[2]-2.56)<0.01","colz");
    TF1* gaus_11_A = new TF1("gaus_11_A","gaus",2.5,2.7);
    gaus_11_A->SetParameters(AOQ_value, angle[0], angle[1]);
    h2_AOQ_11_A->FitSlicesY(gaus_11_A,0,-1,0,"QNRG3",dummy);
    TH1F* h2_AOQ_11_A_1 = (TH1F*) dummy->FindObject("h2_AOQ_11_A_1");
    can_sl_11->cd(1);
    h2_AOQ_11_A_1->Draw();
    TF1* funcangle_11 = new TF1("funcangle_11","pol2",-10,10);
    funcangle_11->FixParameter(0,2.56);
    h2_AOQ_11_A_1->Fit(funcangle_11,"","",-10,10);//Q: quiet, B: user defined parameter settings, M: Improve fit results, N: Do not draw & do not store the graphics
    angle_11[0] = funcangle_11->GetParameter(1);
    angle_11[1] = funcangle_11->GetParameter(2);
    //break;

    //continue with X
    std::cout<< "//****************** X-F10 **************// " <<std::endl;
    can->cd(2);
    TH2F* h2_AOQ_X = new TH2F("h2_AOQ_X","h2_AOQ_X",100,-10,10,100,2.55,2.57);
    tree->Draw(Form("AOQ[2]-F10A*%f-F10Y*%f-delta[2]*%f-BETA[1]*%f-F10A*F10A*%f-F10Y*F10Y*%f-delta[2]*delta[2]*%f-BETA[1]*BETA[1]*%f-F11A*%f-F11X*%f-F11Y*%f-delta[4]*%f-F11A*F11A*%f-F11X*F11X*%f-F11Y*F11Y*%f-delta[4]*delta[4]*%f:F10A>>h2_AOQ_X",
		    angle[0],y_pos[0],delta[0],beta[0],  angle[1],y_pos[1],delta[1],beta[1],      angle_11[0],x_pos_11[0],y_pos_11[0],delta2[0],   angle_11[1],x_pos_11[1],y_pos_11[1],delta2[1]),
	       "TMath::Abs(AOQ[2]-2.56)<0.01","colz");
    can_sl->cd(2);
    TF1* gaus_X = new TF1("gaus_X","gaus",2.5,2.7);
    gaus_X->SetParameters(AOQ_value, x_pos[0], x_pos[1]);
    h2_AOQ_X->FitSlicesY(gaus_X,0,-1,0,"QNRG3",dummy);
    TH1F* h2_AOQ_X_1 = (TH1F*) dummy->FindObject("h2_AOQ_X_1");
    h2_AOQ_X_1->Draw();
    TF1* funcX = new TF1("funcX","pol2",-10,10);
    funcX->FixParameter(0,2.56);
    h2_AOQ_X_1->Fit(funcX,"BM","",-10,10);
    x_pos[0] = funcX->GetParameter(1);
    x_pos[1] = funcX->GetParameter(2);
    //break;    
    std::cout<< "//****************** X-F11 **************// " <<std::endl;
    can_11->cd(2);
    TH2F* h2_AOQ_11_X = new TH2F("h2_AOQ_11_X","h2_AOQ_11_X",100,-10,10,100,2.55,2.57);
    tree->Draw(Form("AOQ[2]-F10A*%f-F10X*%f-F10Y*%f-delta[2]*%f-BETA[1]*%f-F10A*F10A*%f-F10X*F10X*%f-F10Y*F10Y*%f-delta[2]*delta[2]*%f-BETA[1]*BETA[1]*%f-F11A*%f-F11Y*%f-delta[4]*%f-F11A*F11A*%f-F11Y*F11Y*%f-delta[4]*delta[4]*%f:F10A>>h2_AOQ_11_X",
		    angle[0],x_pos[0],y_pos[0],delta[0],beta[0],  angle[1],x_pos[1],y_pos[1],delta[1],beta[1],      angle_11[0],y_pos_11[0],delta2[0],   angle_11[1],y_pos_11[1],delta2[1]),
	       "TMath::Abs(AOQ[2]-2.56)<0.01","colz");
    can_sl_11->cd(2);
    TF1* gaus_11_X = new TF1("gaus_11_X","gaus",2.5,2.7);
    gaus_11_X->SetParameters(AOQ_value, x_pos[0], x_pos[1]);
    h2_AOQ_11_X->FitSlicesY(gaus_11_X,0,-1,0,"QNRG3",dummy);
    TH1F* h2_AOQ_11_X_1 = (TH1F*) dummy->FindObject("h2_AOQ_11_X_1");
    h2_AOQ_11_X_1->Draw();
    TF1* funcX_11 = new TF1("funcX_11","pol2",-10,10);
    funcX_11->FixParameter(0,2.56);
    h2_AOQ_11_X_1->Fit(funcX_11,"BM","",-10,10);
    x_pos_11[0] = funcX_11->GetParameter(1);
    x_pos_11[1] = funcX_11->GetParameter(2);
    //break;    
   
    //Y
    std::cout<< "//****************** Y F10 **************// " <<std::endl;
    can->cd(3);
    TH2F* h2_AOQ_Y = new TH2F("h2_AOQ_Y","h2_AOQ_Y",100,-10,10,100,2.55,2.57);
    tree->Draw(Form("AOQ[2]-F10A*%f-F10X*%f-delta[2]*%f-BETA[1]*%f-F10A*F10A*%f-F10X*F10X*%f-delta[2]*delta[2]*%f-BETA[1]*BETA[1]*%f-F11A*%f-F11X*%f-F11Y*%f-delta[4]*%f-F11A*F11A*%f-F11X*F11X*%f-F11Y*F11Y*%f-delta[4]*delta[4]*%f:F10A>>h2_AOQ_Y",
		    angle[0],x_pos[0],delta[0],beta[0],  angle[1],x_pos[1],delta[1],beta[1],      angle_11[0],x_pos_11[0],y_pos_11[0],delta2[0],   angle_11[1],x_pos_11[1],y_pos_11[1],delta2[1]),
	       "TMath::Abs(AOQ[2]-2.56)<0.01","colz");
    can_sl->cd(3);
    TF1* gaus_Y = new TF1("gaus_Y","gaus",2.5,2.7);
    gaus_Y->SetParameters(AOQ_value, y_pos[0], y_pos[1]);
    h2_AOQ_Y->FitSlicesY(gaus_Y,0,-1,0,"QNRG3",dummy);
    TH1F* h2_AOQ_Y_1 = (TH1F*) dummy->FindObject("h2_AOQ_Y_1");
    h2_AOQ_Y_1->Draw();
    TF1* funcY = new TF1("funcY","pol2",-10,10);
    funcY->FixParameter(0,2.56);
    h2_AOQ_Y_1->Fit(funcY,"BM","",-10,10);
    y_pos[0] = funcY->GetParameter(1);
    y_pos[1] = funcY->GetParameter(2);

    std::cout<< "//****************** Y F11 **************// " <<std::endl;
    can_11->cd(3);
    TH2F* h2_AOQ_11_Y = new TH2F("h2_AOQ_11_Y","h2_AOQ_11_Y",100,-10,10,100,2.55,2.57);
    tree->Draw(Form("AOQ[2]-F10A*%f-F10X*%f-F10Y*%f-delta[2]*%f-BETA[1]*%f-F10A*F10A*%f-F10X*F10X*%f-F10Y*F10Y*%f-delta[2]*delta[2]*%f-BETA[1]*BETA[1]*%f-F11A*%f-F11X*%f-delta[4]*%f-F11X*F11X*%f-F11Y*F11Y*%f-delta[4]*delta[4]*%f:F10A>>h2_AOQ_11_Y",
		    angle[0],x_pos[0],y_pos[0],delta[0],beta[0],  angle[1],x_pos[1],y_pos[1],delta[1],beta[1],      angle_11[0],x_pos_11[0],delta2[0],   angle_11[1],x_pos_11[1],delta2[1]),
	       "TMath::Abs(AOQ[2]-2.56)<0.01","colz");
    can_sl_11->cd(3);
    TF1* gaus_11_Y = new TF1("gaus_11_Y","gaus",2.5,2.7);
    gaus_11_Y->SetParameters(AOQ_value, y_pos[0], y_pos[1]);
    h2_AOQ_11_Y->FitSlicesY(gaus_11_Y,0,-1,0,"QNRG3",dummy);
    TH1F* h2_AOQ_11_Y_1 = (TH1F*) dummy->FindObject("h2_AOQ_11_Y_1");
    h2_AOQ_11_Y_1->Draw();
    TF1* funcY_11 = new TF1("funcY_11","pol2",-10,10);
    funcY_11->FixParameter(0,2.56);
    h2_AOQ_11_Y_1->Fit(funcY_11,"BM","",-10,10);
    y_pos_11[0] = funcY_11->GetParameter(1);
    y_pos_11[1] = funcY_11->GetParameter(2);
    // break;

    /***********  delta missing ***********/

    //BETA
    std::cout<< "//****************** B **************// " <<std::endl;
    can->cd(4);
    TH2F* h2_AOQ_B = new TH2F("h2_AOQ_B","h2_AOQ_B",100,-10,10,100,2.55,2.57);
    tree->Draw(Form("AOQ[2]-F10A*%f-F10X*%f-F10Y*%f-delta[2]*%f-F10A*F10A*%f-F10X*F10X*%f-F10Y*F10Y*%f-delta[2]*delta[2]*%f-BETA[1]*BETA[1]*%f-F11A*%f-F11X*%f-F11Y*%f-delta[4]*%f-F11A*F11A*%f-F11X*F11X*%f-F11Y*F11Y*%f-delta[4]*delta[4]*%f:F10A>>h2_AOQ_B",
		    angle[0],x_pos[0],y_pos[0],delta[0],  angle[1],x_pos[1],y_pos[1],delta[1],      angle_11[0],x_pos_11[0],y_pos_11[0],delta2[0],   angle_11[1],x_pos_11[1],y_pos_11[1],delta2[1]),
	       "TMath::Abs(AOQ[2]-2.56)<0.01","colz");
    TCanvas* can;
    can_sl->cd(4);
    TF1* gaus_B = new TF1("gaus_B","gaus",2.5,2.7);
    gaus_B->SetParameters(AOQ_value, beta[0], beta[1]);
    h2_AOQ_B->FitSlicesY(gaus_B,0,-1,0,"QNRG3",dummy);
    TH1F* h2_AOQ_B_1 = (TH1F*) dummy->FindObject("h2_AOQ_B_1");
    h2_AOQ_B_1->Draw();
    TF1* funcbeta = new TF1("funcbeta","pol2",-10,10);
    funcbeta->FixParameter(0,2.56);
    h2_AOQ_B_1->Fit("pol2","BM","",-10,10);
    beta[0] = funcbeta->GetParameter(1);
    beta[1] = funcbeta->GetParameter(2);

    //look at the AOQ resolution
    can_aoq->cd();
    TH1F* h1_AOQ_B = new TH1F("h1_AOQ_B","h1_AOQ_B",600,2.4,2.7);
    tree->Draw(Form("AOQ[2]-F10A*%f-F10X*%f-F10Y*%f-delta[2]*%f-BETA[1]*%f-F10A*F10A*%f-F10X*F10X*%f-F10Y*F10Y*%f-delta[2]*delta[2]*%f-BETA[1]*BETA[1]*%f-F11A*%f-F11X*%f-F11Y*%f-delta[4]*%f-F11A*F11A*%f-F11X*F11X*%f-F11Y*F11Y*%f-delta[4]*delta[4]*%f>>h1_AOQ_B",
		    angle[0],x_pos[0],y_pos[0],delta[0],beta[0],  angle[1],x_pos[1],y_pos[1],delta[1],beta[1],      angle_11[0],x_pos_11[0],y_pos_11[0],delta2[0],   angle_11[1],x_pos_11[1],y_pos_11[1],delta2[1]),
	       "","colz");

    TF1* funcaoq = new TF1("funcaoq","gaus(0)+pol1(3)",2.55,2.57);
    funcaoq->SetParLimits(0,1e3,1e7);
    funcaoq->SetParLimits(2,1e-5,8e-3);
    funcaoq->SetParLimits(1,2.555,2.565);
    h1_AOQ_B->Fit(funcaoq,"BM","",2.55,2.57);
    std::cout<< " --------------------- " << std::endl;
    std::cout<< "iteration: " << iteration+1 << "  sigma of AOQ: " << funcaoq->GetParameter(2) << std::endl;
    std::cout<<std::endl;
    //break;
    
    if( iteration+1 < max){

      delete funcangle;
      delete funcX;
      delete funcY;

      delete funcangle_11;
      delete funcX_11;
      delete funcY_11;

      delete funcbeta;
      delete funcaoq;

      /*
      delete h2_AOQ_A_0;
      delete h2_AOQ_A_2;     
      delete h2_AOQ_A_chi2;    
      delete h2_AOQ_X_0;
      delete h2_AOQ_X_2;
      delete h2_AOQ_X_chi2;
      delete h2_AOQ_Y_0;
      delete h2_AOQ_Y_2;
      delete h2_AOQ_Y_chi2;
      delete h2_AOQ_B_0;
      delete h2_AOQ_B_2;
      delete h2_AOQ_B_chi2;
      delete h2_AOQ_A_1;
      delete h2_AOQ_X_1;
      delete h2_AOQ_Y_1;     
      delete h2_AOQ_B_1;
      */
      delete h2_AOQ_A;
      delete h2_AOQ_X;
      delete h2_AOQ_Y;
      delete h2_AOQ_B;
     
      delete h2_AOQ_11_A;
      delete h2_AOQ_11_X;
      delete h2_AOQ_11_Y;
     
      dummy->Clear();

      delete h1_AOQ_B;

    }
  }//for



}//improve PID AOQ[2]

