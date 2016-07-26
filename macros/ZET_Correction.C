#define ZET_Correction_cxx
#include "ZET_Correction.h"
#include <TH2.h>
#include "TH2F.h"
#include <TStyle.h>
#include <TCanvas.h>
#include "TF1.h"
#include "TClonesArray.h"

#include<fstream>
#include<stdlib.h>
#include "TString.h"
#include <iostream>
#include "Riostream.h"
#include <string>
#include <time.h>
#include <sstream>
#include "TLine.h"
#include "TArtIC.hh"
#include "TArtICPara.hh"
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
#include <TRandom3.h>
#include "TSystem.h"

#include "TArtStoreManager.hh"
#include "TArtEventStore.hh"


void ZET_Correction::Loop(){

  fChain->SetBranchStatus("*",0);
  fChain->SetBranchStatus("BigRIPSRIPS",1);
  fChain->SetBranchStatus("BigRIPSRIPS.downstream_fpl",1);
  fChain->SetBranchStatus("BigRIPSRIPS.upstream_fpl",1);
  fChain->SetBranchStatus("BigRIPSRIPS.angle",1);
  fChain->SetBranchStatus("BigRIPSPPAC",1);
  fChain->SetBranchStatus("BigRIPSPPAC.fpl",1);
  fChain->SetBranchStatus("BigRIPSPPAC.fX",1);
  fChain->SetBranchStatus("BigRIPSPPAC.fY",1);
  fChain->SetBranchStatus("AOQ");
  fChain->SetBranchStatus("ZET");
  fChain->SetBranchStatus("delta");
  fChain->SetBranchStatus("BETA");
  fChain->SetBranchStatus("BigRIPSPlastic");
  fChain->SetBranchStatus("BigRIPSPlastic.fTime");
  fChain->SetBranchStatus("BigRIPSPlastic.fpl");
  fChain->SetBranchStatus("BigRIPSIC.fpl");
  fChain->SetBranchStatus("BigRIPSIC.fADC[12]");

  Histo();
  Long64_t nentries = fChain->GetEntriesFast();
  if (fChain == 0) return;
   Long64_t nbytes = 0, nb = 0;

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      if( !Incoming() )continue;
      if( !Outgoing() )continue;
      GetVar();
      Plot();//plot the histos
      //fitsl();
   }//for(events)


}//loop

Bool_t ZET_Correction::Incoming(){

  if(TMath::Abs(AOQ[0]-2.64)>0.006)return kFALSE;
  //  else if(TMath::abs(zet[0]-50)>0.5)return kFALSE;

  return kTRUE;

}

Bool_t ZET_Correction::Outgoing(){

  //if(TMath::Abs(aoq[2]-2.64)>0.006)return kFALSE;
  // else if(TMath::Abs(zet[2]-50)>0.5)return kFALSE;

  return kTRUE;
}

void ZET_Correction::GetVar(){//get the needed variables

  for(short irips = 0; irips < kMaxBigRIPSRIPS; irips++){//getting the variables

    if(BigRIPSRIPS_upstream_fpl[irips] == 5){
      fvar[0][0] = BigRIPSRIPS_angle[irips];
    }else if(BigRIPSRIPS_downstream_fpl[irips] == 11){
      fvar[1][0] = BigRIPSRIPS_angle[irips];
    }else continue;

  }//for(rips)
  for(short ippac = 0; ippac < kMaxBigRIPSPPAC; ippac++){//getting the variables
    if(BigRIPSPPAC_fpl[ippac] == 5){
      fvar[0][1] = BigRIPSPPAC_fX[ippac];
      fvar[0][2] = BigRIPSPPAC_fY[ippac];
    }else if(BigRIPSPPAC_fpl[ippac] == 11){
      fvar[1][1] = BigRIPSPPAC_fX[ippac];
      fvar[1][2] = BigRIPSPPAC_fY[ippac];
    }else continue;
      
  }//ippac

  fvar[0][3] = BETA[0];
  fvar[1][3] = BETA[2];

}//getvar()

void ZET_Correction::FitSl(){

  TH1F* h1[2][kMaxICChannel][4];
  TF1*  f1[2][kMaxICChannel][4];
  TObjArray* dummy = new TObjArray();
  dummy->SetOwner(kTRUE);
  for(short iic=0;iic<kMaxBigRIPSIC;iic++){
    for(short ichannel = 0; ichannel<kMaxICChannel; ichannel++){
      fh2_a[iic][ichannel]->FitSlicesY(0,0,-1,0,"QNRG4",dummy);
      h1[iic][ichannel][0] = (TH1F*)dummy->FindObject(Form("fh2_a_%i_%i_1",iic,ichannel));
      fh2_x[iic][ichannel]->FitSlicesY(0,0,-1,0,"QNRG5",dummy);
      h1[iic][ichannel][1] = (TH1F*)dummy->FindObject(Form("fh2_x_%i_%i_1",iic,ichannel));
      fh2_y[iic][ichannel]->FitSlicesY(0,0,-1,0,"QNRG4",dummy);
      h1[iic][ichannel][2] = (TH1F*)dummy->FindObject(Form("fh2_y_%i_%i_1",iic,ichannel));
      fh2_b[iic][ichannel]->FitSlicesY(0,0,-1,0,"QNRG4",dummy);
      h1[iic][ichannel][3] = (TH1F*)dummy->FindObject(Form("fh2_b_%i_%i_1",iic,ichannel));
      for(short ivar = 3; ivar < 4; ivar++){
	cout<< h1[iic][ichannel][ivar] <<endl;
	f1[iic][ichannel][ivar] = new TF1(Form("f1_%i_%i_%i",iic,ichannel,ivar),"pol1",-80,80);
	f1[iic][ichannel][ivar]->SetParameters(10e3,-5);
	h1[iic][ichannel][ivar]->Fit(f1[iic][ichannel][ivar],"","",-80,80);
      }
    }//ichannel
  }//iic
  

}

void ZET_Correction::Histo(){//initialize for histograms

  for(short iic=0;iic<kMaxBigRIPSIC;iic++){
    for(short ichannel = 0; ichannel<kMaxICChannel; ichannel++){
      /*
      fh2_a[iic][ichannel] = new TH2F(Form("fh2_a_%i_%i",iic,ichannel),Form("a vs ic %i ch %i", iic, ichannel),200,-10,10,200,45,55);
      fh2_x[iic][ichannel] = new TH2F(Form("fh2_x_%i_%i",iic,ichannel),Form("x vs ic %i ch %i", iic, ichannel),400,-100,100,200,45,55);
      fh2_y[iic][ichannel] = new TH2F(Form("fh2_y_%i_%i",iic,ichannel),Form("y vs ic %i ch %i", iic, ichannel),200,-10,10,200,45,55);
      fh2_b[iic][ichannel] = new TH2F(Form("fh2_b_%i_%i",iic,ichannel),Form("beta vs ic %i ch %i", iic, ichannel),200,0.5,0.6,200,45,55);
      */
      fh2_a[iic][ichannel] = new TH2F(Form("fh2_a_%i_%i",iic,ichannel),Form("a vs ic %i ch %i", iic, ichannel),200,-10,10,150,0,15e3);
      fh2_x[iic][ichannel] = new TH2F(Form("fh2_x_%i_%i",iic,ichannel),Form("x vs ic %i ch %i", iic, ichannel),400,-100,100,150,0,15e3);
      fh2_y[iic][ichannel] = new TH2F(Form("fh2_y_%i_%i",iic,ichannel),Form("y vs ic %i ch %i", iic, ichannel),200,-10,10,150,0,15e3);
      fh2_b[iic][ichannel] = new TH2F(Form("fh2_b_%i_%i",iic,ichannel),Form("beta vs ic %i ch %i", iic, ichannel),200,0.5,0.6,150,0,15e3);
      
    }//for(channel)
    fh2_z[iic] = new TH2F(Form("fh2_z_%i",iic),Form("z vs aoq %i",iic),200,2.5,2.7,200,45,55);
  }//for(iic)
}//histo()

void ZET_Correction::Plot(){

 for(short iic=0;iic<kMaxBigRIPSIC;iic++){
   for(short ichannel = 0; ichannel<kMaxICChannel; ichannel++){
     /*
      fh2_a[iic][ichannel]->Fill(fvar[iic][0],ZET[iic]);
      fh2_x[iic][ichannel]->Fill(fvar[iic][1],ZET[iic]);
      fh2_y[iic][ichannel]->Fill(fvar[iic][2],ZET[iic]);
      fh2_b[iic][ichannel]->Fill(fvar[iic][3],ZET[iic]);
     */
     fh2_a[iic][ichannel]->Fill(fvar[iic][0],BigRIPSIC_fADC[iic][ichannel]);
     fh2_x[iic][ichannel]->Fill(fvar[iic][1],BigRIPSIC_fADC[iic][ichannel]);
     fh2_y[iic][ichannel]->Fill(fvar[iic][2],BigRIPSIC_fADC[iic][ichannel]);
     fh2_b[iic][ichannel]->Fill(fvar[iic][3],BigRIPSIC_fADC[iic][ichannel]);
      
    }//ichannel
    //    fh2_z[iic]->Fill(AOQ[iic],zet(iic));
 }//iic

}

void ZET_Correction::Draw(){

  TCanvas* can[2][4];
  for(short iic=0;iic<kMaxBigRIPSIC;iic++){
    for(short ivar = 0; ivar<4; ivar++){

      can[iic][ivar] = new TCanvas(Form("can_%i_%i",iic,ivar), Form("ic %i variable %i",iic,ivar));
      can[iic][ivar]->Divide(3,2);
    }
  }

  for(short iic=0;iic<kMaxBigRIPSIC;iic++){
    for(short ichannel = 0; ichannel<kMaxICChannel; ichannel++){
      can[iic][0]->cd(ichannel+1);
      fh2_a[iic][ichannel]->Draw("colz");
      can[iic][1]->cd(ichannel+1);
      fh2_x[iic][ichannel]->Draw("colz");
      can[iic][2]->cd(ichannel+1);
      fh2_y[iic][ichannel]->Draw("colz");
      can[iic][3]->cd(ichannel+1);
      fh2_b[iic][ichannel]->Draw("colz");


    }
  }
}


Double_t ZET_Correction::zet(short ic){

  TArtICPara *para;
  GetIC(ic);
  Float_t E_offset = ESq_offset(9);
  Float_t EnergySqSum = 1;
  for(short ichannel = 0; ichannel < kMaxICChannel; ichannel++){
    EnergySqSum = EnergySqSum*BigRIPSIC_fADC[ic][ichannel];
  }

  EnergySqSum = para->GetCh2MeV(0) + para->GetCh2MeV(1)*EnergySqSum;
  Float_t de_v = TMath::Log(fIC->GetIonPair()*BETA[ic]*BETA[ic])-TMath::Log(1-BETA[ic]*BETA[ic])-BETA[ic]*BETA[ic];
  Float_t Z = fIC->GetZCoef(0)*TMath::Sqrt((EnergySqSum+E_offset)/de_v)*BETA[ic]+fIC->GetZCoef(1);

  return Z;

}

void ZET_Correction::GetIC(short ic){

  gStyle->SetOptStat(111111);
  gSystem->Load("libXMLParser.so");
  gSystem->Load("libanacore.so");
  gSystem->Load("libanabrips.so");
  gROOT->cd();
  //      cout<<$TARTSYS<<endl;                                                                                                                                                                                                          
  gSystem->AddIncludePath("-I$TARTSYS/include/");


  TArtStoreManager * sman = TArtStoreManager::Instance();
  TArtRawEventObject *rawevent = (TArtRawEventObject *)sman->FindDataContainer("RawEvent");

  // Create BigRIPSParameters to get Plastics, PPACs, ICs and FocalPlanes parameters from ".xml" files                                                                                                                                   
  TArtBigRIPSParameters *para = TArtBigRIPSParameters::Instance();
  para -> LoadParameter("db/BigRIPSPPAC.xml");
  para -> LoadParameter("db/BigRIPSPlastic.xml");
  para -> LoadParameter("db/BigRIPSIC.xml");
  para -> LoadParameter("db/FocalPlane.xml");
  para -> LoadParameter("db/BigRIPSTKE.xml");


  TArtRIPS *rips;
  TArtTOF *tof; 
  const char * icname; 
  TArtBeam* beam;
  Bool_t br;
  TArtRecoPID* recopid = new TArtRecoPID();
  if(ic == 0){
    rips = recopid -> DefineNewRIPS(5,7,"matrix/mat2.mat","D5");
    tof = recopid -> DefineNewTOF("F3pl","F7pl",270/*Offset(9,0)*/,5);
    icname = "F7IC";
    br = kTRUE;
  }else{ 
    rips = recopid -> DefineNewRIPS(10,11,"matrix/F10F11_LargeAccAchr_132Sn.mat","D8");
    tof  = recopid -> DefineNewTOF("F8pl","F11pl-1",-158.344/*Offset(9,2)*/,9);
    icname = "F11IC";
    br = kFALSE;
  }
  beam = recopid ->DefineNewBeam(rips, tof, (Char_t*)icname);

  TClonesArray * ic_array = (TClonesArray *)sman->FindDataContainer("BigRIPSIC");
  char name[128];
  Int_t nbeam=1;//GetNumBeam();                                                                                         beam->SetDetectorName(name);

  beam->SetNumRIPS(1);
  beam->SetRIPSName(*(rips->GetDetectorName()));
  beam->SetTOFName(*(tof->GetDetectorName()));
  beam->SetICName(icname);

  TArtIC* ic_;// = NULL;
  Int_t num_ic = ic_array->GetEntries();
  for(Int_t j=0;j<num_ic;j++){
    TArtIC *tmp_ic = (TArtIC *)ic_array->At(j);
    TString icname = *(tmp_ic->GetDetectorName());
    if(icname ==  *(beam->GetICName())){
      //cout<< icname << endl;                                                                                                                                                                                                                 
      ic_ = tmp_ic;
      break;
    }
  }
  if(NULL == ic_){
    TArtCore::Error(__FILE__,"no ic: %s",beam->GetICName()->Data());
  }
  else{
   
    /*    if(Br)fIC2=ic;
	  else*/ fIC=ic_;
  }

}

Float_t ZET_Correction::ESq_offset(Int_t run){

  Int_t mean = 964, peak[2] ={0x0};
  Float_t offset =0;
  string tmp;
  ifstream infile;
  std::string irun, peak_found, line;
  if(run<9){
    mean = 867.05;
    infile.open("cutParameters/PID/empty/peakpos_ESq.txt");
  }else  if(run>90){
    infile.open("cutParameters/PID/Sn128/peakpos_ESq.txt");
    mean = 971;
  }else       infile.open("cutParameters/PID/Sn132/peakpos_ESq.txt");

  while(infile.good()) {
    getline(infile, line);

    std::istringstream iss(line);
    iss >> irun >> peak_found;
    //cout<< irun << endl;                                                                                                                                                                                                                     
    if(atoi(irun.c_str())==run){
      offset = mean - atof(peak_found.c_str());
      break;
    }

  }
  cout<< "ESq - offset " << offset <<endl;

  return offset;
}

Float_t ZET_Correction::Offset(Int_t run, Int_t aoq){

  Float_t time =0;
  string tmp;
  ifstream infile;
  std::string _aoq, offset, line;
  if(run<9)
    infile.open("cutParameters/PID/empty/timeoff.txt");
  else  if(run>90)
    infile.open("cutParameters/PID/Sn128/timeoff.txt");
  else
    infile.open("cutParameters/PID/Sn132/timeoff.txt");

  while(infile.good()) {
    getline(infile, line);

    std::istringstream iss(line);
    iss >> _aoq >> offset;
    //cout<< irun << endl;                                                                                                                                                                                                                     
    if(atoi(_aoq.c_str())==aoq){
      time = atof(offset.c_str());
      break;
    }

  }
  cout<< "time - offset for " << aoq << "  is  " << offset <<endl;

  return time;

}
