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

void analyse::main(short end){

  //TFile* output = new TFile("data/rootfiles/analysis/test.root","REACREATE");
  output = new TFile("data/rootfiles/analysis/test.root","REACREATE");
  Char_t* path = "data/rootfiles/new";
  /** PID **/
  h2_PID_before_0 = new TH2F("h2_PID_before_0", "PID BigRIPS",1000,2.5,2.7, 500,45,55);
  h2_PID_after_0 = new TH2F("h2_PID_after_0", "PID BigRIPS after PID",1000,2.5,2.7, 500,45,55);
  h2_PID_after_2 = new TH2F("h2_PID_after_2", "PID ZeroDegree after PID",1000,2.5,2.7, 500,45,55);
  h2_PID_before_2 = new TH2F("h2_PID_before_2", "PID ZeroDegree",1000,2.5,2.7, 500,45,55);

  /** LaBr **/
  for(short ii = 0; ii<3;ii++){
    h2_LaBr_E[ii] = new TH2F(Form("h2_LaBr_E_%i",ii),"E vs Channel",26,-0.5,25.5,500,0,4e4);  
    for(short ik=0;ik<8;ik++){
      h2_LaBr_EvsT[ii][ik] = new TH2F(Form("h2_LaBr_EvsT_%i_%i",ii,ik),Form("E vs T for channel %i crystal %i",ii,ik),100,100.5,200.5,500,0,4e4);
    }
  }   

  /** DALI **/
  h2_DALI_E = new TH2F("h2_DALI_E","E vs Channel DALI",180,-0.5,179.5,500,0,1e3);
  for(short iDALI=0;iDALI<180;iDALI++){
    
    h2_DALI_EvsT[iDALI]=new TH2F(Form("h2_DALI_EvsT_%03d",iDALI),Form("E vs T for crystal %i",iDALI),100,-1e3,1e3,500,0,4e4);
   }

  /** Loop over all Events **/
  for(UInt_t run = 10; run < end; run++){
    cout<< (run-9)/(end-9.) <<endl;
    fChain = new TChain("tree");
    fChain->AddFile(Form("%s/run%04d.root",path,run));
    cout<< Form("%s/run%04d.root",path,run) <<endl;
    LoadTree(run);
    UInt_t nEvents = fChain->GetEntriesFast();
    set(run);
    nb = 0;
    for(UInt_t iEvent = 0; iEvent < 1e6/*nEvents*/; iEvent++ ){//loop over all events
      if(iEvent%100000==0)cout<< iEvent << endl;     
      Long64_t ientry = fChain->LoadTree(iEvent);
      nb = fChain->GetEntry(ientry);//LoadTree(iEvent);//Load event
     
      // if(!PID())continue;
      //LaBr();
      // DALI();
      for(short iDALI=0; iDALI<180/*DALINaI_ kMaxDALINaI*/; iDALI++){
	//    cout<< iDALI << "  " << DALINaI_id[iDALI] << endl;
	h2_DALI_E->Fill(iDALI,DALINaI_fDoppCorEnergy[iDALI]);
	h2_DALI_EvsT[iDALI]->Fill(DALINaI_fTime[iDALI],DALINaI_fEnergy[iDALI]);

	if(DALINaI_id[iDALI] >121)
	  for(short iCrystal = 0; iCrystal<8; iCrystal++)h2_LaBr_EvsT[1][iCrystal]->Fill(DALINaI_id[iDALI],LaBr_GAINC[1][iCrystal]);
      }//for(DALI)
  
    }//for(tree)

    //    reset();
  }//for(runs)

  output->Write();

}

/********* Analyse the gamma's *********/ 
void analyse::LaBr(){
  /*
  for(short iChannel = 0; iChannel<3; iChannel++){//Loop over all LaBr channels, each LaBr has three channels, Low, Middle and High-Gain

    for(short iCrystal = 0; iCrystal<8; iCrystal++){//Loop over all LaBr crystals
      h2_LaBr_E[iChannel]->Fill(iCrystal, LaBr_GAINC[iChannel][iCrystal]);
      for(short iTDC=0; iTDC<kMaxDALINaI; iTDC++){
	if(DALINaI_id[iTDC]==??)h2_LaBr_EvsT[iChannel][iCrystal]->Fill(DALINaI_fTDC[iTDC], LaBr_GAINC[iChannel][iCrystal]);
      }//iTDC
    }//crystal
  }//channel
  */
  

  /************************
	  change this!!!!!!!!!!
maye into the loop over all dalis and then if(channel > bla) then
  ************************/
  
}//void

void analyse::DALI(){

  //    cout<< 0 << "  " << DALINaI_id[0] << endl;
  for(short iDALI=0; iDALI<180/*DALINaI_ kMaxDALINaI*/; iDALI++){
    //    cout<< iDALI << "  " << DALINaI_id[iDALI] << endl;
    h2_DALI_E->Fill(iDALI,DALINaI_fDoppCorEnergy[iDALI]);
    h2_DALI_EvsT[iDALI]->Fill(DALINaI_fTime[iDALI],DALINaI_fEnergy[iDALI]);

    if(DALINaI_id[iDALI] >121)
      for(short iCrystal = 0; iCrystal<8; iCrystal++)h2_LaBr_EvsT[1][iCrystal]->Fill(DALINaI_id[iDALI],LaBr_GAINC[1][iCrystal]);
  }//for(DALI)
}//void DALI()

/********* TKE and BETA *********/


/********* Functions to cut *********/
Bool_t analyse::PID(){

  //create a 2D cut with an ellipse around the mean values using the sigma
  Double_t dist_aoq_sq = (AOQ[0]-AOQ_ms[0][0])*(AOQ[0]-AOQ_ms[0][0])/(AOQ_ms[1][0]*AOQ_ms[1][0]);
  Double_t dist_zet_sq = (ZET[0]-ZET_ms[0][0])*(ZET[0]-ZET_ms[0][0])/(ZET_ms[1][0]*ZET_ms[1][0]);
  Double_t dist_sq = dist_aoq_sq+dist_zet_sq;
  h2_PID_before_0->Fill(AOQ[0],ZET[0]);
  h2_PID_before_2->Fill(AOQ[2],ZET[3]);
  if(dist_sq>Sigma)return kFALSE; //incoming
  
  /** outgoing **/
  dist_aoq_sq = (AOQ[1]-AOQ_ms[0][1])*(AOQ[1]-AOQ_ms[0][1])/(AOQ_ms[1][1]*AOQ_ms[1][1]);
  //if(dist_aoq_sq>Sigma)return kFALSE;//AOQ 8-10

  dist_aoq_sq = (AOQ[2]-AOQ_ms[0][2])*(AOQ[2]-AOQ_ms[0][2])/(AOQ_ms[1][2]*AOQ_ms[1][2]); //10-11
  dist_zet_sq = (ZET[3]-ZET_ms[0][1])*(ZET[3]-ZET_ms[0][1])/(ZET_ms[1][1]*ZET_ms[1][1]); // 8-10
  dist_sq = dist_aoq_sq+dist_zet_sq;
  //if(dist_sq>Sigma)return kFALSE; //Best 
  
  dist_aoq_sq = (AOQ[5]-AOQ_ms[0][3])*(AOQ[5]-AOQ_ms[0][3])/(AOQ_ms[1][3]*AOQ_ms[1][3]); //8-11
  dist_zet_sq = (ZET[5]-ZET_ms[0][2])*(ZET[5]-ZET_ms[0][2])/(ZET_ms[1][3]*ZET_ms[1][3]); //8-11??
  dist_sq = dist_aoq_sq+dist_zet_sq;
  //if(dist_sq>Sigma)return kFALSE; // 
  h2_PID_after_0->Fill(AOQ[0],ZET[0]);
  h2_PID_after_2->Fill(AOQ[2],ZET[3]);
  return kTRUE;

}
/********* Reset and Set the variables *********/
void analyse::reset(){
 
  for(short ii = 0; ii<2;ii++){
    for(short ik = 0; ik<4; ik++){
      AOQ_ms[ii][ik] = -1;
      ZET_ms[ii][ik] = -1;
    }//for(ik)
  }//for(ii)
  delete fChain;

}

void analyse::set(short run){

  Sigma = 3;
  ReadPIDCuts(run);

}

void analyse::ReadPIDCuts(short run){

  Char_t* path = "cutParameters/PID";
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
	break;
      }
    }//while


  }//for()

}



/*Keep it as last funtion*/
void analyse::LoadTree(short run){
  
  //fChain->SetBranchAddress("DALINaI", &DALINaI_, &b_DALINaI_);
   fChain->SetBranchAddress("DALINaI.fUniqueID", DALINaI_fUniqueID, &b_DALINaI_fUniqueID);
   fChain->SetBranchAddress("DALINaI.fBits", DALINaI_fBits, &b_DALINaI_fBits);
   fChain->SetBranchAddress("DALINaI.id", DALINaI_id, &b_DALINaI_id);
   fChain->SetBranchAddress("DALINaI.fpl", DALINaI_fpl, &b_DALINaI_fpl);
   fChain->SetBranchAddress("DALINaI.name", DALINaI_name, &b_DALINaI_name);
   fChain->SetBranchAddress("DALINaI.fDataState", DALINaI_fDataState, &b_DALINaI_fDataState);
   fChain->SetBranchAddress("DALINaI.fADC", DALINaI_fADC, &b_DALINaI_fADC);
   fChain->SetBranchAddress("DALINaI.fTDC", DALINaI_fTDC, &b_DALINaI_fTDC);
   fChain->SetBranchAddress("DALINaI.layer", DALINaI_layer, &b_DALINaI_layer);
   fChain->SetBranchAddress("DALINaI.theta", DALINaI_theta, &b_DALINaI_theta);
   /*  
       fChain->SetBranchAddress("DALINaI.fXPos", DALINaI_fXPos, &b_DALINaI_fXPos);
   fChain->SetBranchAddress("DALINaI.fYPos", DALINaI_fYPos, &b_DALINaI_fYPos);
   fChain->SetBranchAddress("DALINaI.fZPos", DALINaI_fZPos, &b_DALINaI_fZPos);
   */
   fChain->SetBranchAddress("DALINaI.costheta", DALINaI_costheta, &b_DALINaI_costheta);
   fChain->SetBranchAddress("DALINaI.fEnergy", DALINaI_fEnergy, &b_DALINaI_fEnergy);
   fChain->SetBranchAddress("DALINaI.fDoppCorEnergy", DALINaI_fDoppCorEnergy, &b_DALINaI_fDoppCorEnergy);
   /*
   fChain->SetBranchAddress("DALINaI.fTime", DALINaI_fTime, &b_DALINaI_fTime);
   fChain->SetBranchAddress("DALINaI.fTimeOffseted", DALINaI_fTimeOffseted, &b_DALINaI_fTimeOffseted);
   fChain->SetBranchAddress("DALINaI.fEnergyWithoutT", DALINaI_fEnergyWithoutT, &b_DALINaI_fEnergyWithoutT);  
   fChain->SetBranchAddress("DALINaI.fTimeTrueEnergy", DALINaI_fTimeTrueEnergy, &b_DALINaI_fTimeTrueEnergy);
   fChain->SetBranchAddress("DALINaI.fTimeTrueDoppCorEnergy", DALINaI_fTimeTrueDoppCorEnergy, &b_DALINaI_fTimeTrueDoppCorEnergy);
   fChain->SetBranchAddress("DALINaI.fTimeTrueDoppCorVertexEnergy", DALINaI_fTimeTrueDoppCorVertexEnergy, &b_DALINaI_fTimeTrueDoppCorVertexEnergy);
   fChain->SetBranchAddress("DALINaI.fTimeTrueTime", DALINaI_fTimeTrueTime, &b_DALINaI_fTimeTrueTime);
   fChain->SetBranchAddress("DALINaI.fTimeTrueTimeOffseted", DALINaI_fTimeTrueTimeOffseted, &b_DALINaI_fTimeTrueTimeOffseted);
   */
   fChain->SetBranchAddress("triggerbit", &triggerbit, &b_triggerbit);
   fChain->SetBranchAddress("neve", &neve, &b_neve);
   /*
   fChain->SetBranchAddress("tgtx", &tgtx, &b_tgtx);
   fChain->SetBranchAddress("tgty", &tgty, &b_tgty);
   fChain->SetBranchAddress("tgta", &tgta, &b_tgta);
   fChain->SetBranchAddress("tgtb", &tgtb, &b_tgtb);
   fChain->SetBranchAddress("F10X", &F10X, &b_F10X);
   fChain->SetBranchAddress("F10Y", &F10Y, &b_F10Y);
   fChain->SetBranchAddress("F10A", &F10A, &b_F10A);
   fChain->SetBranchAddress("F10B", &F10B, &b_F10B);
   fChain->SetBranchAddress("F11X", &F11X, &b_F11X);
   fChain->SetBranchAddress("F11Y", &F11Y, &b_F11Y);
   fChain->SetBranchAddress("F11A", &F11A, &b_F11A);
   fChain->SetBranchAddress("F11B", &F11B, &b_F11B);
   fChain->SetBranchAddress("F3PLA_QL", &F3PLA_QL, &b_F3PLA_QL);
   fChain->SetBranchAddress("F3PLA_QR", &F3PLA_QR, &b_F3PLA_QR);
   fChain->SetBranchAddress("F7PLA_QL", &F7PLA_QL, &b_F7PLA_QL);
   fChain->SetBranchAddress("F7PLA_QR", &F7PLA_QR, &b_F7PLA_QR);
   fChain->SetBranchAddress("F8PLA_QL", &F8PLA_QL, &b_F8PLA_QL);
   fChain->SetBranchAddress("F8PLA_QR", &F8PLA_QR, &b_F8PLA_QR);
   fChain->SetBranchAddress("F11PLA_QL", &F11PLA_QL, &b_F11PLA_QL);
   fChain->SetBranchAddress("F11PLA_QR", &F11PLA_QR, &b_F11PLA_QR);
   */
   fChain->SetBranchAddress("DELTA", DELTA, &b_DELTA);
   fChain->SetBranchAddress("ANGLE", ANGLE, &b_ANGLE);
   fChain->SetBranchAddress("BRHO", BRHO, &b_BRHO);
   fChain->SetBranchAddress("TOF", TOF, &b_TOF);
   fChain->SetBranchAddress("BETA", BETA, &b_BETA);
   fChain->SetBranchAddress("AOQ", AOQ, &b_AOQ);
   fChain->SetBranchAddress("ZET", ZET, &b_ZET);
   fChain->SetBranchAddress("TKE", &TKE, &b_TKE);
   /*
   fChain->SetBranchAddress("ASUM", &ASUM, &b_ASUM);
   fChain->SetBranchAddress("dalimultwotime", &dalimultwotime, &b_dalimultwotime);
   fChain->SetBranchAddress("dalimult", &dalimult, &b_dalimult);
   fChain->SetBranchAddress("dalitimetruemult", &dalitimetruemult, &b_dalitimetruemult);
   fChain->SetBranchAddress("dalimultthres", &dalimultthres, &b_dalimultthres);
   fChain->SetBranchAddress("dalitimetruemultthres", &dalitimetruemultthres, &b_dalitimetruemultthres);
   */
   fChain->SetBranchAddress("LOWGAIN", LaBr_GAIN[0], &b_LOWGAIN);
   fChain->SetBranchAddress("LOWGAINC", LaBr_GAINC[0], &b_LOWGAINC);
   fChain->SetBranchAddress("MIDGAIN", LaBr_GAIN[1], &b_MIDGAIN);
   fChain->SetBranchAddress("MIDGAINC", LaBr_GAINC[1], &b_MIDGAINC);
   fChain->SetBranchAddress("HIGAIN", LaBr_GAIN[2], &b_HIGAIN);
   fChain->SetBranchAddress("HIGAINC", LaBr_GAINC[2], &b_HIGAINC);
     
}
