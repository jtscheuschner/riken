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


void PID(Bool_t Sn132){

  FILE *aoqfile[4];
  FILE *zetfile[3];

  TFile* outfile = new TFile("cutParameters/PID_Sn128.root","RECREATE");
  std::vector< std::vector<TH1F*> >h1_aoq;
  std::vector< std::vector<TH1F*> > h1_zet;

  std::vector< std::vector<TF1*> >f1_aoq;
  std::vector< std::vector<TF1*> >f1_zet;

  const short num_can = 51;
  TCanvas* can_aoq[num_can];
  TCanvas* can_zet[num_can];
  for(short jj = 0; jj<num_can;jj++){
    can_aoq[jj] = new TCanvas(Form("can_aoq_%02d",jj),Form("can_aoq run %02d",jj+9));
    can_aoq[jj]->Divide(2,2);
    can_zet[jj] = new TCanvas(Form("can_zet_%02d",jj),Form("can_zet run %02d",jj+9));
    can_zet[jj]->Divide(2,2);
  }
  h1_aoq = std::vector< std::vector<TH1F*> >(4);
  h1_zet = std::vector< std::vector<TH1F*> >(3);
  f1_aoq = std::vector< std::vector<TF1*> > (4);
  f1_zet = std::vector< std::vector<TF1*> > (3);
  for(short ii = 0; ii < 4; ii++){
    h1_aoq[ii] = std::vector<TH1F*>(56);
    f1_aoq[ii] = std::vector<TF1*>(56);
    
    if(ii < 3){
      h1_zet[ii] = std::vector<TH1F*>(56);
      f1_zet[ii] = std::vector<TF1*>(56);
      aoqfile[ii] = fopen(Form("cutParameters/peakpos_AOQ_%i.txt",ii),"w");
    }else aoqfile[ii] = fopen(Form("cutParameters/peakpos_AOQ_%i.txt",5),"w");

  }
  zetfile[0] = fopen(Form("cutParameters/peakpos_ZET_%i.txt",0),"w");
  zetfile[1] = fopen(Form("cutParameters/peakpos_ZET_%i.txt",3),"w");
  zetfile[2] = fopen(Form("cutParameters/peakpos_ZET_%i.txt",5),"w");
 
  short ican=0;
  short start = 9; 
  short   end = 64;
  Float_t correction[3][16]={0x0};

  if(!Sn132){
    start = 101;
    end   = 152;//151

    ifstream infile;
    string dummy, line;
    infile.open("cutParameters/PID/Sn128/AOQ2_correction.txt");
    short ii = 0;
    while(infile.good()){
      getline(infile, line);
      std::istringstream iss(line);
      iss >> dummy >> correction[ii][0] >>  correction[ii][1] >>  correction[ii][2] >>  correction[ii][3] >>  correction[ii][4] >>  correction[ii][5] >>  correction[ii][6] >>  correction[ii][7] >>  correction[ii][8] >>  correction[ii][9] >>  correction[ii][10] >>  correction[ii][11] >>  correction[ii][12] >>  correction[ii][13] >>  correction[ii][14] >>  correction[ii][15];
      ii++;
    }//while
    
    for(short ii =0; ii<3;ii++)
      cout<< dummy << "  " << correction[ii][0] << "  " << correction[ii][1] << "  " <<  correction[ii][2] << "  " <<  correction[ii][3] << "  " <<  correction[ii][4] << "  " <<  correction[ii][5] << "  " <<  correction[ii][6] << "  " <<  correction[ii][7] << "  " <<  correction[ii][8] << "  " <<  correction[ii][9] << "  " <<  correction[ii][10] << "  " <<  correction[ii][11] << "  " <<  correction[ii][12] << "  " <<  correction[ii][13] << "  " <<  correction[ii][14] << "  " <<  correction[ii][15] <<endl;
  }//if


  Float_t aoq = 0;
  if(Sn132) aoq = 2.64;
  else      aoq = 2.56;
  cout<< "aoq = " << aoq <<endl;
  for(short run = start; run<end; run++){
    TChain* tree = new TChain("tree");
    tree->AddFile(Form("data/rootfiles/new/run%04d.root",run));

    h1_zet[0][run-start] = new TH1F(Form("h1_zet_%i_%02d",0,run),Form("ZET[%i] run %02d",0,run),500,45,55);
    f1_zet[0][run-start] = new TF1(Form("f1_zet_%i_%02d",0,run),"gaus",49.5,51);
    h1_zet[1][run-start] = new TH1F(Form("h1_zet_%i_%02d",3,run),Form("ZET[%i] run %02d",3,run),500,45,55);
    f1_zet[1][run-start] = new TF1(Form("f1_zet_%i_%02d",3,run),"gaus",49.,51);
    h1_zet[2][run-start] = new TH1F(Form("h1_zet_%i_%02d",5,run),Form("ZET[%i] run %02d",5,run),500,45,55);
    f1_zet[2][run-start] = new TF1(Form("f1_zet_%i_%02d",5,run),"gaus",49.,53);

    for(short ii = 0; ii < 4; ii++){     

      if(ii<3){	
	h1_aoq[ii][run-start] = new TH1F(Form("h1_aoq_%i_%02d",ii,run),Form("AOQ[%i] run %02d",ii,run),600,2.5,2.7);
	if(Sn132)f1_aoq[ii][run-start] = new TF1(Form("f1_aoq_%i_%02d",ii,run),"gaus",2.63,2.665);
	else f1_aoq[ii][run-start] = new TF1(Form("f1_aoq_%i_%02d",ii,run),"gaus",2.55,2.585);
     
      f1_zet[ii][run-start]->SetParameter(2,1);
      f1_zet[ii][run-start]->SetParameter(1,50);

      }else{
	h1_aoq[ii][run-start] = new TH1F(Form("h1_aoq_%i_%02d",5,run),Form("AOQ[%i] run %02d",5,run),600,2.5,2.7);
	if(Sn132)f1_aoq[ii][run-start] = new TF1(Form("f1_aoq_%i_%02d",5,run),"gaus",2.63,2.665);
	else f1_aoq[ii][run-start] = new TF1(Form("f1_aoq_%i_%02d",5,run),"gaus",2.55,2.585);
      }//else
      f1_aoq[ii][run-start]->SetParLimits(0,100,1e6);
      f1_aoq[ii][run-start]->SetParameter(1,12e3);
      f1_aoq[ii][run-start]->SetParameter(1,aoq);
      f1_aoq[ii][run-start]->SetParLimits(1,aoq-0.1,aoq+0.1);
      f1_aoq[ii][run-start]->SetParameter(2,2e-3);
      f1_aoq[ii][run-start]->SetParLimits(2,0,5e-3);

    }//for (ii)
    //First BigRIPS than ZDS

    //AOQ
    const Char_t* aoq_0 = "",* aoq_2 = "",* aoq_5 = "";
  
    aoq_0 = "AOQ[0]";
    aoq_2 = "AOQ[2]";
    aoq_5 = "AOQ[5]";
    Bool_t b_correction = kFALSE;
    if(b_correction){
      aoq_0 = Form("%s+%f*FA[0]+%f*FA[0]*FA[0]+%f*FX[0]+%f*FX[0]*FX[0]+%f*FY[0]+%f*FY[0]*FY[0]+delta[0]*%f+delta[0]*delta[0]*%f+%f*FA[1]+%f*FA[1]*FA[1]+%f*FX[1]+%f*FX[1]*FX[1]+%f*FY[1]+%f*FY[1]*FY[1]+delta[1]*%f+delta[1]*delta[1]*%f",
		   aoq_0,
		   correction[0][0],correction[0][1],correction[0][2],correction[0][3],correction[0][4],correction[0][5],correction[0][6],correction[0][7],//First part of the BRHO
		   correction[0][8],correction[0][9],correction[0][10],correction[0][11],correction[0][12],correction[0][13],correction[0][14],correction[0][15]);//Second part of the BRHO

      aoq_2 = Form("%s+%f*FA[4]+%f*FA[4]*FA[4]+%f*FX[4]+%f*FX[4]*FX[4]+%f*FY[4]+%f*FY[4]*FY[4]+delta[2]*%f+delta[2]*delta[2]*%f+%f*FA[5]+%f*FA[5]*FA[5]+%f*FX[5]+%f*FX[5]*FX[5]+%f*FY[5]+%f*FY[5]*FY[5]+delta[4]*%f+delta[4]*delta[4]*%f",
		   aoq_2
		   ,correction[1][0],correction[1][1],correction[1][2],correction[1][3],correction[1][4],correction[1][5],correction[1][6],correction[1][7]//First part of the BRHO
		   ,correction[1][8],correction[1][9],correction[1][10],correction[1][11],correction[1][12],correction[1][13],correction[1][14],correction[1][15]//Second part of the BRHO
		   );

      aoq_5 = Form("%s+%f*FA[3]+%f*FA[3]*FA[3]+%f*FX[3]+%f*FX[3]*FX[3]+%f*FY[3]+%f*FY[3]*FY[3]+delta[2]*%f+delta[2]*delta[2]*%f+%f*FA[5]+%f*FA[5]*FA[5]+%f*FX[5]+%f*FX[5]*FX[5]+%f*FY[5]+%f*FY[5]*FY[5]+delta[4]*%f+delta[4]*delta[4]*%f",
		   aoq_5,
		   correction[2][0],correction[2][1],correction[2][2],correction[2][3],correction[2][4],correction[2][5],correction[2][6],correction[2][7],//First part of the BRHO
		   correction[2][8],correction[2][9],correction[2][10],correction[2][11],correction[2][12],correction[2][13],correction[2][14],correction[2][15]);//Second part of the BRHO
    }
    

    //cout<< aoq_2 <<endl;
    //cout<< aoq_5 <<endl;

    //Draw and Fit

    //AOQ
    can_aoq[ican]->cd(1);
    tree->Draw(Form("%s>>h1_aoq_0_%02d",aoq_0,run),"TMath::Abs(ZET[0]-50)<0.5","",200000,0);
    //h1_aoq[0][run-start]->Fit(Form("f1_aoq_0_%02d",run),"R","",aoq-0.1,aoq+0.1);
    can_aoq[ican]->cd(2);
    tree->Draw(Form("AOQ[1]>>h1_aoq_1_%02d",run)  ,Form("TMath::Abs(AOQ[0]-%f)<0.005 && TMath::Abs(ZET[0]-50)<0.5 && TMath::Abs(ZET[3]-50)<0.5",aoq),"",200000,0);
    //h1_aoq[1][run-start]->Fit(Form("f1_aoq_1_%02d",run),"R");
    can_aoq[ican]->cd(3);
    tree->Draw(Form("%s>>h1_aoq_2_%02d",aoq_2,run),Form("TMath::Abs(AOQ[0]-%f)<0.005 && TMath::Abs(ZET[0]-50)<0.5 && TMath::Abs(ZET[3]-50)<0.5",aoq),"",200000,0);
    //h1_aoq[2][run-start]->Fit(Form("f1_aoq_2_%02d",run),"R");
    can_aoq[ican]->cd(4);
    tree->Draw(Form("%s>>h1_aoq_5_%02d",aoq_5,run),Form("TMath::Abs(AOQ[0]-%f)<0.005 && TMath::Abs(ZET[0]-50)<0.5 && TMath::Abs(ZET[3]-50)<0.5",aoq),"",200000,0);
    //    h1_aoq[3][run-start]->Fit(Form("f1_aoq_5_%02d",run),"R");

    h1_aoq[3][run-start]->Fit(Form("f1_aoq_5_%02d",run),"R");
    can_aoq[ican]->cd(3);
    h1_aoq[2][run-start]->Fit(Form("f1_aoq_2_%02d",run),"R");
    can_aoq[ican]->cd(2);
    h1_aoq[1][run-start]->Fit(Form("f1_aoq_1_%02d",run),"R");
    can_aoq[ican]->cd(1);
    h1_aoq[0][run-start]->Fit(Form("f1_aoq_0_%02d",run),"R","",aoq-0.1,aoq+0.1);
    //ZET
    can_zet[ican]->cd(1);
    tree->Draw(Form("ZET[0]>>h1_zet_0_%02d",run),"TMath::Abs(AOQ[0]-aoq)<0.005","",200000,0);
    h1_zet[0][run-start]->Fit(Form("f1_zet_0_%02d",run),"R");
    can_zet[ican]->cd(2);
    tree->Draw(Form("ZET[3]>>h1_zet_3_%02d",run),"TMath::Abs(AOQ[0]-aoq)<0.005 && TMath::Abs(ZET[0]-50)<0.5 && TMath::Abs(AOQ[2]-aoq)<0.005","",200000,0);
    h1_zet[1][run-start]->Fit(Form("f1_zet_3_%02d",run),"R");
    can_zet[ican]->cd(3);
    tree->Draw(Form("ZET[5]>>h1_zet_5_%02d",run),"TMath::Abs(AOQ[0]-aoq)<0.005 && TMath::Abs(ZET[0]-50)<0.5 && TMath::Abs(AOQ[2]-aoq)<0.005","",200000,0);
    h1_zet[2][run-start]->Fit(Form("f1_zet_5_%02d",run),"R");
    for(short jj = 0; jj<4; jj++){   
              fprintf(aoqfile[jj],"%i %f %f \n",run,f1_aoq[jj][run-start]->GetParameter(1),f1_aoq[jj][run-start]->GetParameter(2));
      if(jj<3)fprintf(zetfile[jj],"%i %f %f \n",run,f1_zet[jj][run-start]->GetParameter(1),f1_zet[jj][run-start]->GetParameter(2));
    }
    ican++;
  }//for(run)

  for(short jj = 0; jj<4; jj++){
    fclose(aoqfile[jj]);
    if(jj<3) fclose(zetfile[jj]);
  }
    
  outfile->Write();

}//void
