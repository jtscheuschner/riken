#define analyse_DALI_cxx
#include "analyse_DALI.h"
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


//#include "addback.C"


void analyse_DALI::Loop(short end, short trig, Float_t _sigma, Bool_t pieter){

  if(_sigma <= 0){
    cout<< " sigma smaller or equal to '0' please check your parameters!!" <<endl;
    return;
  }
  Sigma = _sigma;

  trigg = trig;
  if(trig <1 || trig >4) cout<< " **** All trigger are accepted!!  ****" <<endl;
  else if(trig==0) cout<< " **** Only trigger "<< trig << " (Beam) is accepted!!  ****" <<endl; 
  else if(trig==1) cout<< " **** Only trigger "<< trig << " (LaBr x Beam) is accepted!!  ****" <<endl; 
  else if(trig==2) cout<< " **** Only trigger "<< trig << " (DALI x Beam) is accepted!!  ****" <<endl; 
  else if(trig==3) cout<< " **** Only trigger "<< trig << " (DS(DALI x Beam) ) is accepted!!  ****" <<endl; 
  else if(trig==4) cout<< " **** Only trigger "<< trig << " (SumDALI x Beam) is accepted!!  ****" <<endl; 
  selected_events = 0;
  Bool_t b_tke = kTRUE, b_beta = kTRUE, b_beam = kTRUE;
  output = new TFile(Form("data/rootfiles/analysis/test_trig_%i.root",trig),"RECREATE");
  short run = 9;
  Histo();
  for(short run=40;run<end;run++){  

    reset();
    set(run, pieter);

    LoadData(0,run);

    if (fChain == 0) return;
    Long64_t nentries = fChain->GetEntriesFast();

    cout << "run " << run << " number of events " << nentries <<endl;

    Long64_t nbytes = 0, nb = 0;

    for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      //reset parameters run by run
      reset_run();
      F7T = BigRIPSPlastic_fTime[1];
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      /*
      if(b_tke)Tke();
      if(b_beta)Beta();
      if(b_beam)BeamProfile();
      */
      //Event selection
      if( !trigger() )continue;
      if( !PID(run) )continue;

      if(b_beam)BeamProfile();
      selected_events++;
      ADC();
      QTC();
      /*
      if(b_tke)Tke();
      if(b_beta)Beta();
      if(b_beam)BeamProfile();
      */
    }//for(entry)
  }//for(run)
  output->Write();
  cout<< " Selected Events "<< selected_events <<endl;
}

//*************  Gamma's *************//
void analyse_DALI::ADC(){

  //short hits[150]={0x0}; 
  for(short iDALI=0; iDALI<DALINaI_/*kMaxDALINaI*/; iDALI++){
    //    cout<< iDALI << "  " << DALINaI_id[iDALI] << endl;    

    //if(DALINaI_fEnergy[iDALI] > 100)hits[ DALINaI_id[iDALI] ];
    h2_HitsEvent->Fill(DALINaI_id[iDALI],selected_events);
    if( !ADC_time(iDALI) )continue;    

    if(DALINaI_id[iDALI] <124)       
      h2_DALI_EvsT[DALINaI_id[iDALI]]->Fill(DALINaI_fTimeOffseted[iDALI]+F8T+30,DALINaI_fEnergy[iDALI]);    
    else  h2_DALI_EvsT[DALINaI_id[iDALI]]->Fill(DALINaI_fTimeOffseted[iDALI]+F8T+95,DALINaI_fEnergy[iDALI]);    
    if(DALINaI_id[iDALI] >124)    h2_DALI_E->Fill(DALINaI_id[iDALI],DALINaI_fDoppCorEnergy[iDALI]);
   

  }//for(DALI)
  /*
  for(short ii = 0; ii<150;ii++){//test for Multihit in one crystal: non observed: maybe change readout system
    if(hits[ii] > 1)
      cout<< "crystal with Multihits: " << DALINaI_id[ii] <<endl;
  }
  */
}// void ADC

void analyse_DALI::QTC(){

  for(short icrystal = 0; icrystal<8; icrystal++){//loop over all crystal
    // if( !QTC_time(icrystal) )continue;
    Double_t LaBr_time = QTCTimeC[icrystal]+F7T-F8T;
    for(short ichannel = 0; ichannel < 3;ichannel++){//loop over all channels, each crystal has 3
      /*
      Float_t correction = +48 - 17e3/LaBr_GAIN[ichannel][icrystal];
      if(TMath::Abs(LaBr_time)<5){
      */
      h2_LaBr_EvsT[ichannel][icrystal]->Fill(LaBr_time,LaBr_GAIN[ichannel][icrystal]);
      if( !QTC_time(icrystal, ichannel,LaBr_time) )continue;
      //      if( !QTC_time(icrystal) )continue;
      //if(TMath::Abs(LaBr_time)<10)
      h2_LaBr_E[ichannel]->Fill(icrystal,LaBr_GAINC[ichannel][icrystal]);
      //}
     
    }//for(channel)
  }//for(crystal)

}

//************* Addback for DALI          **************//
/*
void analyse_DALI::Addback(){



  for(short iDALI_st=0; iDALI_st<DALINaI_; iDALI_st++){
    if( !ADC_time(iDALI_st) )continue;
    if(DALINaI_id[iDALI_st]>124)continue;

    for(short iDALI_sec=0; iDALI_sec<DALINaI_; iDALI_sec++){
      if( iDALI_sec==iDALI_st)continue;//exclude same detector
      if( DALINaI_id[iDALI_sec]>124)continue;
      if( Distance(iDALI_st, iDALI_sec) > Max_Dist)continue;
      if( !ADC_time(iDALI_sec) )continue;




    }//for(sec)
  }//for(first)

}//Addback
*/
//************* TKE & BETA and other stuff *************//
void analyse_DALI::Tke(){//Look into the TKE and TKE depending on other variables
  
  h2_TKE_beta->Fill(BETA[1],TKE);
  h2_TKE_A->Fill(F11A,TKE);
  h2_TKE_X->Fill(F11X,TKE);
  h2_TKE_Y->Fill(F11Y,TKE);
  h2_TKE_delta->Fill(delta[3],TKE);
  h2_TKE_AOQ->Fill(AOQ[2],TKE);
  h2_TKE_ZET->Fill(ZET[3],TKE);
}//void TKE

void analyse_DALI::Beta(){//Look into beta this variable can be useful to decrease the background in the crystals (DALI and LaBr)



  for(short ii=0; ii < kMaxBigRIPSPPAC; ii++){
    if(BigRIPSPPAC_fpl[ii]==7){
      F7X = BigRIPSPPAC_fX[ii];
      F7Y = BigRIPSPPAC_fY[ii];
    }else
    if(BigRIPSPPAC_fpl[ii]==8){
      F8X = BigRIPSPPAC_fX[ii];
      F8Y = BigRIPSPPAC_fY[ii];
    }else if(BigRIPSPPAC_fpl[ii]==9){
      F9X = BigRIPSPPAC_fX[ii];
      F9Y = BigRIPSPPAC_fY[ii];
    }
  }//for

  h2_BETA_TOF->Fill(TOF[3],BETA[1]);
  h2_BETA_AOQ->Fill(AOQ[2],BETA[1]);
  h2_BETA_A->Fill(F7A,BETA[1]);
  h2_BETA_X->Fill(F7X,BETA[1]);
  //  h2_BETA_Y->Fill(F7Y,BETA[1]);
  h2_BETA_delta->Fill(delta[2],BETA[1]);
  h2_BETA_beta->Fill(BETA[0],BETA[1]);
}//void beta

void analyse_DALI::BeamProfile(){//Look into the beamprofile before and after the target, x-y-pos, angular distribution etc.
  //this might help to get rid of ions scattered on the thick chamber-walls

  h2_Beam_F7_XY->Fill(F7X,F7Y);
  h2_Beam_F8_XY->Fill(F8X,F8Y);
  h2_Beam_F9_XY->Fill(F9X,F9Y);

  h2_Beam_AX->Fill(F9X,F9A);
  h2_Beam_AA->Fill(F7A,F9A);

  h2_Beam_delta_delta->Fill(DELTA[4],DELTA[2]);
  h2_Beam_XX->Fill(F7X,F9X);

}//void BeamProfile

//************* selection criteria *************//

/********

The time cut needs a slope-correction

********/
//Gamma criteria
Bool_t analyse_DALI::ADC_time(short dali){//cut on the time-signal of the detector
  //this has to be done more properly but for the first time being this should be good enough
  //the timing should be in some correlation with the time of flight of the beam to the target from F7 to F8

  if(DALINaI_id[dali] < 124){
    if(TMath::Abs(DALINaI_fTimeOffseted[dali]+F8T+30)>4)return kFALSE;
  }else if(TMath::Abs(DALINaI_fTimeOffseted[dali]+F8T+97)>4)return kFALSE;   

  return kTRUE;
}

Bool_t analyse_DALI::QTC_time(short crystal, short channel, Double_t LaBr_time){//time cut on the LaBr crystals read-out with the QTC's
  
  Double_t correction = tri_exp(LaBr_GAIN[1][crystal]);//tripple exponential
  if(channel == 1){
    if(TMath::Abs(LaBr_time-correction)>1)return kFALSE;
  }//else if(channel == 2){} else if(channel == 1){}
  return kTRUE;
}
//Incoming and outgoing PID
Bool_t analyse_DALI::PID(short run){//Do the incoming and outgoing PID with AOQ and ZET the criteria are set in set()

  //create a 2D cut with an ellipse around the mean values using the sigma
  Double_t dist_aoq_sq = (AOQ[0]-AOQ_ms[0][0])*(AOQ[0]-AOQ_ms[0][0])/(AOQ_ms[1][0]*AOQ_ms[1][0]);
  Double_t dist_zet_sq = (ZET[0]-ZET_ms[0][0])*(ZET[0]-ZET_ms[0][0])/(ZET_ms[1][0]*ZET_ms[1][0]);
  Double_t dist_sq = dist_aoq_sq+dist_zet_sq;
  h2_PID_before_0->Fill(AOQ[0],ZET[0]);
  
  if(dist_sq>Sigma)return kFALSE; //incoming
  h2_PID_before_2->Fill(AOQ[2],ZET[3]);
  /** outgoing **/
  dist_aoq_sq = (AOQ[1]-AOQ_ms[0][1])*(AOQ[1]-AOQ_ms[0][1])/(AOQ_ms[1][1]*AOQ_ms[1][1]);
  if(dist_aoq_sq>Sigma)return kFALSE;//AOQ 8-10

  dist_aoq_sq = (AOQ[2]-AOQ_ms[0][2])*(AOQ[2]-AOQ_ms[0][2])/(AOQ_ms[1][2]*AOQ_ms[1][2]); //10-11
  dist_zet_sq = (ZET[3]-ZET_ms[0][1])*(ZET[3]-ZET_ms[0][1])/(ZET_ms[1][1]*ZET_ms[1][1]); // 8-11
  dist_sq = dist_aoq_sq+dist_zet_sq;
  //if(dist_sq>Sigma)return kFALSE; //Best 
  
  dist_aoq_sq = (AOQ[5]-AOQ_ms[0][3])*(AOQ[5]-AOQ_ms[0][3])/(AOQ_ms[1][3]*AOQ_ms[1][3]); //8-11
  //dist_zet_sq = (ZET[3]-ZET_ms[0][1])*(ZET[3]-ZET_ms[0][2])/(ZET_ms[1][3]*ZET_ms[1][3]); //8-11??
  dist_sq = dist_aoq_sq+dist_zet_sq;
  //if(dist_sq>Sigma)return kFALSE; //

  /** delta cut on outgoing **/ 
  h2_PID_delta->Fill(delta[2],delta[4]);
  if( TMath::Abs(delta[2]-delta[4]) > 0.5) return kFALSE;

  h2_PID_after_0->Fill(AOQ[0],ZET[0]);
  h2_PID_after_2->Fill(AOQ[2],ZET[3]);

  //maybe add delta-cut
  return kTRUE;

}//PID()

Bool_t analyse_DALI::trigger(){

  if(trigg<0)    return kTRUE;
  else if(trigg==triggerbit) return kTRUE;
  
  return kFALSE;
    

}

//*************  Define histograms *************//
void analyse_DALI::Histo(){

  /** PID **/
  h2_PID_before_0 = new TH2F("h2_PID_before_0", "PID BigRIPS",1000,2.5,2.7, 500,45,55);
  h2_PID_after_0 = new TH2F("h2_PID_after_0", "PID BigRIPS after PID",1000,2.5,2.7, 500,45,55);
  h2_PID_after_2 = new TH2F("h2_PID_after_2", "PID ZeroDegree after PID",1000,2.5,2.7, 500,45,55);
  h2_PID_before_2 = new TH2F("h2_PID_before_2", "PID ZeroDegree",1000,2.5,2.7, 500,45,55);

  h2_PID_delta = new TH2F("h2_PID_delta", "delta[2] vs delta[4]",100,-10,10,100,-10,10);
  /** LaBr **/
  for(short ii = 0; ii<3;ii++){
    h2_LaBr_E[ii] = new TH2F(Form("h2_LaBr_E_%i",ii),"E vs Channel",9,-0.5,8.5,2500,0,5e4);  
    for(short ik=0;ik<8;ik++){
      h2_LaBr_EvsT[ii][ik] = new TH2F(Form("h2_LaBr_EvsT_%i_%i",ii,ik),Form("E vs T for channel %i crystal %i",ii,ik),400,-1e2,1e2,1000,500,5e4);
    }
  }   

  /** DALI **/
  h2_DALI_E = new TH2F("h2_DALI_E","E vs Channel DALI",150,-0.5,149.5,500,0,1e4);
  for(short iDALI=0;iDALI<150;iDALI++){    
    h2_DALI_EvsT[iDALI]=new TH2F(Form("h2_DALI_EvsT_%03d",iDALI),Form("E vs T for crystal %i",iDALI),200,-100,100,500,0,4e4);
  }

  //Control Multihits
  h2_HitsEvent = new TH2F("h2_hitsEvent", "hits per crystal and per event",150,-0.5,149.5,200000,0,2e6);
  /** TKE **/
  h2_TKE_beta = new TH2F("h2_TKE_beta","TKE vs #beta",200,0.4,0.6,1000,0,5e3);
  h2_TKE_A = new TH2F("h2_TKE_A","TKE vs angular",200,-20,20,1000,0,5e3);
  h2_TKE_X = new TH2F("h2_TKE_X","TKE vs X",200,-20,20,1000,0,5e3);
  h2_TKE_Y = new TH2F("h2_TKE_Y","TKE vs Y",200,-20,20,1000,0,5e3);
  h2_TKE_delta = new TH2F("h2_TKE_delta","TKE vs delta",200,-20,20,1000,0,5e3);
  h2_TKE_AOQ = new TH2F("h2_TKE_AOQ","TKE vs AOQ",200,2.50,2.70,1000,0,5e3);
  h2_TKE_ZET = new TH2F("h2_TKE_ZET","TKE vs ZET",200,45,55,1000,0,5e3);

  /** BETA **/
  h2_BETA_TOF = new TH2F("h2_BETA_TOF","#beta vs TOF",200,-1e3,1e3,200,0.4,0.6);
  h2_BETA_AOQ = new TH2F("h2_BETA_AOQ","#beta vs A/Q",200,2.5,2.7,200,0.4,0.6);
  h2_BETA_A = new TH2F("h2_BETA_A","#beta vs A",200,-20,20,200,0.4,0.6);
  h2_BETA_X = new TH2F("h2_BETA_X","#beta vs X",200,-20,20,200,0.4,0.6);
  h2_BETA_Y = new TH2F("h2_BETA_Y","#beta vs Y",200,-20,20,200,0.4,0.6);
  h2_BETA_delta = new TH2F("h2_BETA_delta","#beta vs delta",200,-2,2,200,0.4,0.6);
  h2_BETA_beta = new TH2F("h2_BETA_beta","#beta vs #beta[0]",200,0.4,0.6,200,0.4,0.6);

  /** BeamProfile **/
  h2_Beam_F9_XY = new TH2F("h2_Beam_F9_XY","Beam F9 x-y",200,-20,20,200,-20,20);
  h2_Beam_F8_XY = new TH2F("h2_Beam_F8_XY","Beam F8 x-y",200,-20,20,200,-20,20);
  h2_Beam_F7_XY = new TH2F("h2_Beam_F7_XY","Beam F7 x-y",200,-20,20,200,-20,20);
  h2_Beam_AX = new TH2F("h2_Beam_AX","Beam F9 A-x",200,-20,20,200,-20,20);

  h2_Beam_AA = new TH2F("h2_Beam_AA","Beam F7-F9 A-A",200,-20,20,200,-20,20);
  h2_Beam_delta_delta = new TH2F("h2_Beam_delta_delta","Beam F7-F9 delta-delta",200,-5,5,200,-5,5);
  h2_Beam_XX = new TH2F("h2_Beam_XX","Beam F7-F9 x-x",200,-20,20,200,-20,20);

}//Histo()

//************* Setter and Reseter *************//
void analyse_DALI::reset(){
 
  for(short ii = 0; ii<2;ii++){
    for(short ik = 0; ik<4; ik++){
      AOQ_ms[ii][ik] = -1;
      ZET_ms[ii][ik] = -1;
    }//for(ik)
  }//for(ii)
  delete fChain;

}

void analyse_DALI::set(short run, Bool_t pieter){//Set the initial variables which maybe are depending on the run
  

  ReadPIDCuts(run);
  if(pieter){
    for(short ii = 1; ii<4; ii++){
      AOQ_ms[0][ii]=2.67346938775510212;
      ZET_ms[0][ii]=49;
    }
  }

}


void analyse_DALI::reset_run(){

  F7X = -1000;
  F8X = -1000;
  F9X = -1000;

  F7Y = -1000;
  F8Y = -1000;
  F9Y = -1000;

}

void analyse_DALI::ReadPIDCuts(short run){//read the evaluated PID mean- and sigma-values (found with PID.C)

  Char_t* path;
  if(run < 100) path = "cutParameters/PID/Sn132";
  else path = "cutParameters/PID/Sn128";
  std::string crun, mean, sigma, line;
  for(short ii = 0; ii<4;ii++){
    ifstream aoqfile, zetfile;
    aoqfile.open(Form("%s/peakpos_AOQ_%i.txt",path,ii));
    while(aoqfile.good()){
      getline(aoqfile,line);
      std::istringstream iss(line);
      iss >> crun >> mean >> sigma;
      if(atoi(crun.c_str()) == run){
	AOQ_ms[0][ii]=atof(mean.c_str());
	AOQ_ms[1][ii]=atof(sigma.c_str());
	//	cout<< 	AOQ_ms[1][ii] << "  " << AOQ_ms[1][ii] <<endl;
	break;
      }
    }//while
    zetfile.open(Form("%s/peakpos_ZET_%i.txt",path,ii));
    while(zetfile.good()){
      getline(zetfile,line);
      std::istringstream iss(line);
      iss >> crun >> mean >> sigma;
      if(atoi(crun.c_str()) == run){
	ZET_ms[0][ii]=atof(mean.c_str());
	ZET_ms[1][ii]=atof(sigma.c_str());
	//	cout<< ZET_ms[0][ii] << "  " << ZET_ms[1][ii] << endl;
	break;
      }
    }//while
  }//for()

}//ReadPIDCuts()


//***************** Define Functions

Double_t analyse_DALI::tri_exp(Double_t En){

  /* TF1* slopecorrection = new TF1("slopecorrection","expo(0)+expo(2)+expo(4)+[6]",100,1e4);
    slopecorrection->SetParameters(7.69801e+00,
				    1.06976e-07,
				    3.12043e+00,
				    -2.95502e-04,
				    3.95422e+00,
				    -2.77578e-03,
				    -2.19232e+03);
  */
  Double_t val = 
    TMath::Exp(7.69801e+00 +  1.06976e-07 *En)+
    TMath::Exp(3.12043e+00 + -2.95502e-04 *En)+
    TMath::Exp(3.95422e+00 + -2.77578e-03 *En)+
    -2.19232e+03;

  return val;
}
