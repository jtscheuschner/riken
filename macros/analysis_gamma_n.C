#define analysis_gamma_cxx
#include "analysis_gamma.h"
#include <iostream>
#include "Riostream.h"
#include <string>
#include <time.h>
#include <sstream>
#include "TThread.h"
#include "TF1.h"
#include "snprintf.h"

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void analysis_gamma::Loop(UInt_t start, UInt_t end, Char_t* savename, Float_t addback_dist){
//   In a ROOT session, you can do:
//      Root > .L analysis_gamma.C
//      Root > analysis_gamma t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch


  faddback_dist = addback_dist;
  ReadPIDCut();
  if(end<10)      savename = Form("%sempty/%2.1f",savename, faddback_dist);
  else if(end<100)savename = Form("%sSn132/%2.1f",savename, faddback_dist);
  else if(end<200)savename = Form("%sSn128/%2.1f",savename, faddback_dist);
  else            savename = Form("%scalib/%04d_%2.1f.root",savename,start,faddback_dist);
  if(end<200){
    for(int i = 0; i<8; i++)  
      savename = Form("%s_%1.1f",savename, nsigma[i]);
    savename = Form("%s.root",savename);
  }  
  cout<< "savename: " << savename <<endl;
  
  Histo(savename);
  for(UInt_t run = start; run<end; run++){
    cout<<"run  "<< run <<endl;
    Tree(run);
    Read(run);
    if (fChain == 0) return;

    Long64_t nentries = fChain->GetEntriesFast();

    Long64_t nbytes = 0, nb = 0;
    if(run<200){
      for (Long64_t jentry=0; jentry<nentries;jentry++) {
	Long64_t ientry = LoadTree(jentry);
	if (ientry < 0) break;
	nb = fChain->GetEntry(jentry);   nbytes += nb;
	// if (Cut(ientry) < 0) continue;

	if(!PID_cut())continue;
	LaBr();
	DALI();
	Gamma_Gamma();
      }//for(event
    }else{
      for (Long64_t jentry=0; jentry<nentries;jentry++) {
	Long64_t ientry = LoadTree(jentry);
	if (ientry < 0) break;
	nb = fChain->GetEntry(jentry);   nbytes += nb;
      
	LaBr();
	DALI();
      }    
      TFile* f = new TFile(Form("data/rootfiles/new_/run%04d_1.root",run));
      if(f->IsOpen()){
	TTree* tree;
	f->GetObject("tree",tree);
	Init(tree);
	run--;
      }
    }
  }//for(run
  foutfile->Write();
  foutfile->Close();
}//Loop

void analysis_gamma::LaBr(){
  fGamma_LaBr=0;//Reset for gamma-gamma correlation
  Double_t theta = TMath::Cos(30./180*TMath::Pi());
  for(UInt_t ii= 0; ii<kMaxLaBr;ii++){
    fEnergy_LaBr[ii] = 0;
    MIDGAINC[ii] = DALIDoppCor(MIDGAIN[ii],theta,BETA[2]);
    LOWGAINC[ii] = DALIDoppCor(LOWGAIN[ii],theta,BETA[2]);
    QTCTimeC[ii] = QTCTimeC[ii] - fQTC_mean[0][ii];
  }
  for(UInt_t ilabr=0; ilabr<kMaxLaBr;ilabr++){
    h2_labr_E_T[0][ilabr]->Fill(LOWGAIN[ilabr],QTCTimeC[ilabr]);
    h2_labr_E_T[1][ilabr]->Fill(MIDGAIN[ilabr],QTCTimeC[ilabr]);

    if(!Gamma_t_cut(0,ilabr))continue;
    h2_labr_AOQ->Fill(AOQ[2],LOWGAINC[ilabr]);
    h2_labr_E_LOW->Fill(ilabr,LOWGAIN[ilabr]);
    h2_labr_E_MID->Fill(ilabr,MIDGAIN[ilabr]);
    h2_labr_EdoppCor->Fill(ilabr,LOWGAINC[ilabr]);
    h2_labr_EdoppCor_MID->Fill(ilabr,MIDGAINC[ilabr]);
    fEnergy_LaBr[fGamma_LaBr]=LOWGAINC[ilabr];
    fGamma_LaBr++;
  }//for
}//LaBr

void analysis_gamma::DALI(){
  for(UInt_t ii = 0; ii<kMaxDALINaI;ii++)//reset E for Addback
    fEnergy_DALI[ii]=0;

  for(UInt_t idali=0; idali < DALINaI_/*kMaxDALINaI*/; idali++){
    if(DALINaI_id[idali]>150)continue;
    fTime_DALI[idali] = DALI_T(idali);
    fDALI_Energy[idali] = DALIEnergy(DALINaI_fADC[idali],DALINaI_id[idali]);
    h2_dali_adc_E->Fill(DALINaI_id[idali], DALINaI_fADC[idali]);
    h2_dali_E_T[DALINaI_id[idali]]->Fill(fDALI_Energy[idali],fTime_DALI[idali]);

    if(!Gamma_t_cut(1,idali))continue;
    fDALI_DoppCorEnergy[idali] = DALIDoppCor(fDALI_Energy[idali],DALINaI_costheta[idali],BETA[2]);
    h2_dali_EdoppCor->Fill(DALINaI_id[idali],fDALI_DoppCorEnergy[idali]);
    h2_dali_E       ->Fill(DALINaI_id[idali],fDALI_Energy[idali]);
    if(DALINaI_id[idali]<124){
      h2_dali_EdoppCor_layer->Fill(DALINaI_layer[idali],fDALI_DoppCorEnergy[idali]);
    }else{//labr-adc
      h2_dali_EdoppCor_layer->Fill(0.,fDALI_DoppCorEnergy[idali]);
    }
    if(DALINaI_id[idali]>123)continue;//exclude LaBr
    //if(DALINaI_layer[idali]>7)continue;//exclude backward
    fEnergy_DALI[idali]=fDALI_Energy[idali];//Save the Energies for addback in separate array
  }//for()
  DALIAddback(faddback_dist);//addback distance in cm
}

void analysis_gamma::DALIAddback(Float_t dist){
  Int_t *index; index = new Int_t[kMaxDALINaI];//Sort the array starting with the biggest one
  TMath::Sort((Int_t)kMaxDALINaI,fEnergy_DALI,index,kTRUE);//it is only index sorting so be careful elements will not be switched within the array!!
  Bool_t used[kMaxDALINaI] = {0x0};//check if crystal has already be used
  fGamma_DALI=0;//counting how many gammas have been registered
  for(UInt_t idali=0; idali < kMaxDALINaI; idali++){
    if(fEnergy_DALI[index[idali]]<150)continue;//check if already at the end
    if(!Gamma_t_cut(1,index[idali]))continue;
    if(used[idali])continue;//check if used
    fAddbackE_DALI[fGamma_DALI] = fEnergy_DALI[index[idali]];
    used[idali] = 1;

    for(UInt_t idalisec=0; idalisec < kMaxDALINaI; idalisec++){
      if(fEnergy_DALI[index[idalisec]]<150)continue;//check energy
      if(used[idalisec])continue;// check that it is not used already
      if(!Gamma_t_cut(1,index[idalisec]))continue;
      if(DALIDist(DALINaI_id[index[idali]],DALINaI_id[index[idalisec]])<dist){//check the distance
	fAddbackE_DALI[fGamma_DALI] = fAddbackE_DALI[fGamma_DALI] + fEnergy_DALI[index[idalisec]];
	used[idalisec] = 1;
      }//if
    }//for    
    fAddbackE_DALI[fGamma_DALI] = DALIDoppCor(fAddbackE_DALI[fGamma_DALI],DALINaI_costheta[index[idali]],BETA[2]);
    h1_dali_Addback->Fill(fAddbackE_DALI[fGamma_DALI]);
    fGamma_DALI++;
  }//for

}//DALIAddback

Double_t analysis_gamma::DALIDist(UInt_t id_1, UInt_t id_2){//calculates dist between two dalis
  Double_t distx = (DALINaI_fX[id_1]-DALINaI_fX[id_2]);//*(DALINaI_fX[0id_1]]-DALINaI_fX[0id_2]])+
  Double_t disty = (DALINaI_fY[id_1]-DALINaI_fY[id_2]);//*(DALINaI_fY[0id_1]]-DALINaI_fY[0id_2]])+
  Double_t distz = (DALINaI_fZ[id_1]-DALINaI_fZ[id_2]);//*(DALINaI_fZ[0id_1]]-DALINaI_fZ[0id_2]]);
  Double_t dist  = distx*distx+disty*disty+distz*distz;
  return TMath::Sqrt(dist);
}//dalidist

Double_t analysis_gamma::DALIEnergy(Int_t adc, Int_t id){
  /*if(adc>0 && adc < 4e3 && id > 0)*/
  return (adc*fDALI_calibpara[id][1]+fDALI_calibpara[id][0]);
  return adc;
  //else                return -1000.;
}

Double_t analysis_gamma::DALIDoppCor(Double_t energy, Double_t costheta, Double_t beta){
  Double_t doppCorE = -1;
  doppCorE = energy * (1- beta *costheta)/TMath::Sqrt(1-beta*beta);
  return doppCorE;
  return energy;
}

Double_t analysis_gamma::DALI_T(Int_t id){
  //for(short ii = 0; ii<150;ii++)
  
  Double_t time = -9999;
  Float_t length = 124.8;//cm
  Double_t clight= 29.979245;// cm/ns
  
  time = DALINaI_fTimeOffseted[id]-length/clight/BETA[0]-DALI_offset[id][0]-1400;
  return time;
}

void analysis_gamma::Gamma_Gamma(){
  Int_t *index; index = new Int_t[kMaxDALINaI];//Sort the array starting with the biggest one
  TMath::Sort((Int_t)kMaxDALINaI,fEnergy_DALI,index,kTRUE);//it is only index sorting so be careful elements will not be switched within the array!!

  Float_t energy_dali=0, energy_labr=0;
  for(Int_t idali = 0; idali < fGamma_DALI; idali++){
    for(Int_t ilabr = 0; ilabr < fGamma_LaBr; ilabr++){
      h2_dali_labr_addback->Fill(fAddbackE_DALI[idali],fEnergy_LaBr[ilabr]);
      if( Gamma_t_cut(1,index[idali]) )
	h2_dali_labr->Fill(fEnergy_DALI[index[idali]],fEnergy_LaBr[ilabr]);
    }
    for(Int_t idali2 = 0; idali2 < fGamma_DALI; idali2++){
      h2_dali_dali->Fill(fAddbackE_DALI[idali2],fAddbackE_DALI[idali]);
      if( Gamma_t_cut(1,index[idali]) &&  Gamma_t_cut(1,index[idali2]) )
      h2_dali_dali_without->Fill(fEnergy_DALI[index[idali]],fEnergy_DALI[index[idali2]]);
    }
    energy_dali+=fAddbackE_DALI[idali];
  }
  for(Int_t ilabr = 0; ilabr < fGamma_LaBr; ilabr++){
    energy_labr+=fEnergy_LaBr[ilabr];
    for(Int_t ilabr2 = ilabr+1; ilabr2 < fGamma_LaBr; ilabr2++)
      h2_labr_labr->Fill(fEnergy_LaBr[ilabr2],fEnergy_LaBr[ilabr]);
  }
  h1_sum->Fill(energy_labr+energy_dali);
}//Gamma_Gamma

//******* Cuts on the Gamma's
Bool_t analysis_gamma::Gamma_t_cut(Bool_t dali, UInt_t id){
  if(dali){
    //if(DALINaI_id[id]<124){
      if(TMath::Abs(fTime_DALI[id]/DALI_offset[id][1])>nsigma[5])return kFALSE;
      /*}else{
      if(TMath::Abs(fTime_DALI[id])>50)return kFALSE;
      }*/
  }else
    //if(TMath::Abs(QTCTimeC[id]/fQTC_sigma[0][id])>nsigma[6])return kFALSE;
    if(LOWGAIN[id]<5e3){
      if(TMath::Abs(QTCTimeC[id]/fQTC_sigma[0][id])>nsigma[6])return kFALSE;
    }else
      if(TMath::Abs((QTCTimeC[id]-fQTC_mean[1][id])/fQTC_sigma[1][id])>nsigma[6])return kFALSE;
  
  return kTRUE;
}

//******* Cuts on the Event/PID
Bool_t analysis_gamma::PID_cut(){
  /*
  AOQ_mean[0]  = _aoq,  AOQ_mean[1] = _aoq;
  ZET_mean[0]  = 50,    ZET_mean[1] = 50;//change to read-in
  ZET_sigma[0] = 0.2,   ZET_sigma[1]= 0.2;
  AOQ_sigma[0] = 0.004, AOQ_sigma[1]= 0.007;
  */

  AOQ[0] = AOQ[0]+_aoq-AOQ_mean[0];
  AOQ[1] = AOQ[1]+_aoq-AOQ_mean[1];
  AOQ[2] = AOQ[2]+_aoq-AOQ_mean[2];
  AOQ[3] = AOQ[3]+_aoq-AOQ_mean[3];
  ZET[0] = ZET[0]+_zet -ZET_mean[0];
  ZET[1] = ZET[1]+_zet -ZET_mean[1];
  ZET[2] = ZET[2]+_zet -ZET_mean[2];
  ZET[3] = ZET[3]+_zet -ZET_mean[3];

  h2_Z_B->Fill(BETA[2],ZET[2]);

  if( (TRIGGER != nsigma[7] && nsigma[7]>-1) || (nsigma[7]==-5 && TRIGGER<2) )return kFALSE;

  h2_PID_in->Fill(AOQ[0],ZET[0]);
  h2_PID_in1->Fill(BigRIPSBeam_aoq[1],BigRIPSBeam_zet[1]);
  h2_PID_in2->Fill(BigRIPSBeam_aoq[2],BigRIPSBeam_zet[2]);
  if(!PID_in())return kFALSE;
  h2_PID_in0->Fill(AOQ[0],ZET[0]);
  h2_PID_out1->Fill(AOQ[1],ZET[1]);
  h2_PID_out2->Fill(AOQ[2],ZET[2]);
  h2_PID_out3->Fill(AOQ[3],ZET[3]);
  h2_PID_out->Fill(AOQ[2],ZET[3]);
  if(!PID_out())return kFALSE;
  h2_AOQ_out_a1->Fill(AOQ[3],ZET[3]);
  h2_AOQ_out_a2->Fill(AOQ[2],ZET[3]);
  h2_x_y_before->Fill(FX[3],FY[3]);
  h2_x_y_after ->Fill(FX[4],FY[4]);
  h2_a_a       ->Fill(FA[3],FA[4]);
  h2_beta_beta ->Fill(BETA[0],BETA[2]);
  return kTRUE;
}

Bool_t analysis_gamma::PID_out(){
  Float_t aoq = _aoq, zet = _zet;
  //Float_t aoq = 2.67543e+00, zet = 49.;
  Double_t dist_aoq = (AOQ[3]-aoq)*(AOQ[3]-aoq)/(nsigma[2]*nsigma[2]*AOQ_sigma[3]*AOQ_sigma[3]);
  Double_t dist_zet = (ZET[3]-zet)*(ZET[3]-zet)/(nsigma[3]*nsigma[3]*ZET_sigma[3]*ZET_sigma[3]);

  //if(dist_zet > n) return kFALSE;
  //h2_PID_out_trig->Fill(AOQ[2],TRIGGER);
  if(TMath::Sqrt(dist_zet)>1)return kFALSE; 
  if(TMath::Abs(BigRIPSRIPS_delta[2]-BigRIPSRIPS_delta[4]-1.33930e-01)>nsigma[4])return kFALSE;
  h2_PID_out_trig->Fill(AOQ[3],TRIGGER);
  if(TMath::Sqrt(dist_aoq+dist_zet)>1)return kFALSE;
  h2_delta->Fill(BigRIPSRIPS_delta[2],BigRIPSRIPS_delta[2]-BigRIPSRIPS_delta[4]);

  /*
  if(TMath::Abs(ZET[2]-ZET_mean[1])>n*ZET_sigma[1])return kFALSE; 
  if(TMath::Abs(AOQ[2]-AOQ_mean[1])>n*AOQ_sigma[1])return kFALSE;
  h2_delta->Fill(BigRIPSRIPS_delta[2],BigRIPSRIPS_delta[4]);
  if(TMath::Abs(BigRIPSRIPS_delta[2]-BigRIPSRIPS_delta[4])>0.5)return kFALSE;
  */
 
  return kTRUE;  
}

Bool_t analysis_gamma::PID_in(){
  Float_t aoq = _aoq, zet = _zet;
  Double_t dist_aoq = (AOQ[0]-aoq)*(AOQ[0]-aoq)/(nsigma[0]*nsigma[0]*AOQ_sigma[0]*AOQ_sigma[0]);
  Double_t dist_zet = (ZET[0]-zet)*(ZET[0]-zet)/(nsigma[1]*nsigma[1]*ZET_sigma[0]*ZET_sigma[0]);
  if(TMath::Sqrt(dist_aoq+dist_zet)>1)return kFALSE;
  /*  
  if(TMath::Abs(AOQ[0]-AOQ_mean[0])>n*AOQ_sigma[0])return kFALSE;
  if(TMath::Abs(ZET[0]-ZET_mean[0])>n*ZET_sigma[0])return kFALSE;
  */
  return kTRUE;
}

void analysis_gamma::Histo(Char_t* savename){//initialize histos and savings
  foutfile = new TFile(Form("%s",savename),"RECREATE");
  ftop = foutfile->mkdir("top");
  ftop->cd();

  //*****dali
  fdali = ftop->mkdir("DALI");
  fdali->cd();
  h2_dali_E_T = vector<TH2F*>(kMaxDALICrystal);
  h2_dali_EdoppCor_layer = new TH2F("h2_dali_EdoppCor_layer", "Dopp corrected Energy vs layer",13,-0.5,12.5,4000,0,4e4);
  h2_dali_E        = new TH2F("h2_dali_E", "Energy vs id",150,-0.5,149.5,4000,0,4e4);
  h2_dali_EdoppCor = new TH2F("h2_dali_EdoppCor","Dopp Correceted Energy vs Crystal-id",150,-0.5,149.5,1000,0,5e4);
  for(UInt_t ii=0;ii<kMaxDALICrystal; ii++){
    h2_dali_E_T[ii] = new TH2F(Form("h2_dali_E_T_%i",ii),Form("Energy vs Time in Crystal %i",ii),200,0,4e4,600,-3e2,3e2);
  }
  h1_dali_Addback = new TH1F("h1_dali_Addback","Addback DALI",1000,0,4e4);
  h2_dali_adc_E = new TH2F("h2_dali_adc_E","",150,-0.5,149.5,1500,0,3e3);
  //*****labr
  flabr = ftop->mkdir("LaBr");
  flabr->cd();
  h2_labr_E_T    = vector<vector<TH2F*> >(2);
  h2_labr_E_T[0] = vector<TH2F*>(kMaxLaBr);
  h2_labr_E_T[1] = vector<TH2F*>(kMaxLaBr);
  for(UInt_t ii = 0; ii < kMaxLaBr; ii++){
    h2_labr_E_T[0][ii] = new TH2F(Form("h2_labr_E_T_LOW_%i",ii), Form("LOWEnergy vs Time in Crystal %i",ii),200,0,4e4,200,-50,50);
    h2_labr_E_T[1][ii] = new TH2F(Form("h2_labr_E_T_MID_%i",ii),Form("MIDEnergy vs Time in Crystal %i",ii),2000,0,2e4,200,-50,50);
  }//for
  h2_labr_EdoppCor = new TH2F("h2_labr_EdoppCor","Dopp corrected Energy vs Crystal-id",9,-0.5,8.5,1000,0,5e4);
  h2_labr_EdoppCor_MID = new TH2F("h2_labr_EdoppCor_MID","Dopp corrected Energy vs Crystal-id",9,-0.5,8.5,150,0,15e3);
  h2_labr_E_LOW = new TH2F("h2_labr_E_LOW","Energy vs Crystal-id",9,-0.5,8.5,1000,0,5e4);
  h2_labr_E_MID = new TH2F("h2_labr_E_MID","Energy vs Crystal-id",9,-0.5,8.5,150,0,15e3);
  h2_labr_AOQ = new TH2F("h2_labr_AOQ","AOQ vs LaBr E", 200,2.5,2.7,200,0,2e4);
  //*****corr
  fcorr = ftop->mkdir("Corr - gamma");
  fcorr->cd();
  h2_dali_dali = new TH2F("h2_dali_dali","Corr dali - dali",400,0,4e4,400,0,4e4);
  h2_dali_dali_without = new TH2F("h2_dali_dali_without","Corr dali - dali no addback",400,0,4e4,400,0,4e4);
  h2_dali_labr = new TH2F("h2_dali_labr","Corr dali - labr",400,0,4e4,400,0,4e4);
  h2_dali_labr_addback = new TH2F("h2_dali_labr_addback","Corr dali - dali",400,0,4e4,400,0,4e4);
  h2_labr_labr = new TH2F("h2_labr_labr","Corr labr - labr",400,0,4e4,400,0,4e4);
  h1_sum = new TH1F("h1_sum","sum energy of all #gamma's",400,0,4e4);

  //*****PID
  fpid = ftop->mkdir("PID");
  fpid->cd();
  h2_Z_B =   = new TH2F ("h2_PID_in"  , "Incoming A/Q vs ZET",  400,0.5,0.6,200,45,55);
  h2_PID_in  = new TH2F ("h2_PID_in"  , "Incoming A/Q vs ZET",  400,2.5,2.7,200,45,55);
  h2_PID_in0 = new TH2F ("h2_PID_in0" , "Incoming A/Q0 vs ZET", 400,2.5,2.7,200,45,55);
  h2_PID_in1 = new TH2F ("h2_PID_in1" , "Incoming A/Q1 vs ZET", 400,2.5,2.7,200,45,55);
  h2_PID_in2 = new TH2F ("h2_PID_in2" , "Incoming A/Q2 vs ZET", 400,2.5,2.7,200,45,55);
  h2_PID_out1 = new TH2F("h2_PID_out1", "Outgoing A/Q1 vs ZET", 400,2.5,2.7,200,45,55);
  h2_PID_out2 = new TH2F("h2_PID_out2", "Outgoing A/Q2 vs ZET", 400,2.5,2.7,200,45,55);
  h2_PID_out3 = new TH2F("h2_PID_out3", "Outgoing A/Q2 vs ZET", 400,2.5,2.7,200,45,55);
  h2_PID_out  = new TH2F("h2_PID_out" , "Outgoing A/Q vs ZET",  400,2.5,2.7,200,45,55);

  h2_AOQ_out_a1 = new TH2F("h2_AOQ_out_a1", "Outgoing A/Q vs ZET 3", 400,2.5,2.7,200,45,55);
  h2_AOQ_out_a2 = new TH2F("h2_AOQ_out_a2", "Outgoing A/Q vs ZET 2", 400,2.5,2.7,200,45,55);
  h2_PID_out_trig = new TH2F("h2_AOQ_out_trig","Outgoing A/Q vs trig", 400, 2.5,2.7,5,-0.5,4.5);
  h2_delta = new TH2F("h2_delta", "Delta F10 - F11",100,-6,6,200,-5,5);

  h2_x_y_before = new TH2F("h2_x_y_before" , "X-Y pos at F8",  200,-20,20, 200,-20,20);
  h2_x_y_after  = new TH2F("h2_x_y_after"  , "X-Y pos at F10", 200,-20,20, 200,-20,20);
  h2_a_a        = new TH2F("h2_a_a", "angle F10 vs F8"       , 200,-20,20, 200,-20,20);
  h2_beta_beta  = new TH2F("h2_b_b","Outgoing b vs incoming b",200,0.5,0.6,200,0.5,0.6);

}
void analysis_gamma::Read(UInt_t run){
  ReadDALI(run);
  ReadPID(run);
  ReadQTC();
}

void analysis_gamma::ReadQTC(){
  TString path = "cutParameters/QTC_parameter.txt";
  ifstream infile;
  std::string _crystal,_mean,_sigma, line;
  UInt_t ii=0;
  infile.open(path);
  while(infile.good()) {
    getline(infile, line);
    std::istringstream iss(line);
    iss >> _crystal >> _mean >> _sigma;
    fQTC_mean[0][ii] =atof(_mean.c_str());
    fQTC_sigma[0][ii]=atof(_sigma.c_str());
    ii++;
  }//while

  path = "cutParameters/QTC_parameter2.txt";
  ifstream infile2;
  ii=0;
  infile2.open(path);
  while(infile2.good()) {
    getline(infile2, line);
    std::istringstream iss(line);
    iss >> _crystal >> _mean >> _sigma;
    fQTC_mean[1][ii] =atof(_mean.c_str());
    fQTC_sigma[1][ii]=atof(_sigma.c_str());
    ii++;
  }//while
}

void analysis_gamma::ReadDALI(UInt_t run){
  TString path = "";
  if(run<9)          path = "cutParameters/Sn132_target_DALI_E.csv";
  else if(run < 100 )path = "cutParameters/Sn132_target_DALI_E.csv";
  else               path = "cutParameters/Sn132_target_DALI_E.csv";
  //cout<<path<<endl;
  ifstream infile;
  std::string slope, offset, line;
  UInt_t ii=1;
  infile.open(path);
  while(infile.good()) {
    getline(infile, line);
    std::istringstream iss(line);
    iss >> offset >> slope;
    fDALI_calibpara[ii][0]=atof(offset.c_str());
    fDALI_calibpara[ii][1]=atof(slope.c_str());
    ii++;
    if(ii==150)break;
  }//while
  ReadDALI_T();
  ReadDALIPos();
}

void analysis_gamma::ReadDALI_T(){

  ifstream infile;
  std::string crystal, offset, sigma, line;
  UInt_t ii=1;
  infile.open("cutParameters/time_dali2.txt");
  while(infile.good()) {
    getline(infile, line);
    std::istringstream iss(line);
    iss >> crystal >> offset >> sigma;
    ii = atoi(crystal.c_str());
    DALI_offset[ii][0] = atof(offset.c_str());
    //DALI_offset[ii][1] = atof(sigma.c_str());
    ii++;
    if(ii==150)break;
  }//while
/*
  ifstream infile2;
  infile2.open("cutParameters/time_dali2.txt");
  while(infile2.good()) {
    getline(infile, line);
    std::istringstream iss(line);
    iss >> crystal >> offset >> sigma;
    ii = atoi(crystal.c_str());
    DALI_offset[ii][0] = DALI_offset[ii][0]+atof(offset.c_str());
    DALI_offset[ii][1] = atof(sigma.c_str());
    ii++;
    if(ii==150)break;
  }//while
*/
}//ReadDALI

void analysis_gamma::ReadDALIPos(){
  TString path = "cutParameters/DALI_Pos.csv";
  string tmp;
  ifstream infile;
  std::string xx, yy, zz, line;
  UInt_t ii=1;
  infile.open(path);
  while(infile.good()) {
    getline(infile, line);
    std::istringstream iss(line);
    iss >> xx >> yy >> zz;
    DALINaI_fX[ii]=atof(xx.c_str());
    DALINaI_fY[ii]=atof(yy.c_str());
    DALINaI_fZ[ii]=atof(zz.c_str());
    ii++;
    if(ii==124)break;
  }//while
}//ReadDALI


void analysis_gamma::ReadPID(UInt_t run){

  TString path = "cutParameters/peakpos_";/*
  if(run<9)path ="cutParameters/PID/empty/peak_";
  else if(run<100)path ="cutParameters/PID/Sn132/peak_";
  else path ="cutParameters/PID/Sn128/peak_";
		    
  Char_t* name[4] = {"aoq_in","aoq_out","zet_in","zet_out"};
				       */
  TString name[9] = {"AOQ_0","AOQ_1","AOQ_2","AOQ_3","ZET_0","ZET_1","ZET_2","ZET_3",".txt"};
  std::string _run, _mean, _sigma, line;

  for(int ii = 0; ii<4;ii++){
    UInt_t irun = 9999;
    ifstream infile;
    TString bla = path+name[ii]+name[8];
    infile.open(bla);
    //infile.open(Form("%s%s.txt",path,name[ii]));
    while(infile.good()) {
      //cout<<"file aoq "<< path << name <<endl;
      getline(infile, line);
      std::istringstream iss(line);
      iss >> _run >> _mean >> _sigma;
      irun = atoi(_run.c_str());
      if(irun==run){
	AOQ_mean[ii]  = atof(_mean.c_str());
	AOQ_sigma[ii] = atof(_sigma.c_str());
	break;
      }      
    }//while
    
    irun = -1;
    ifstream infile2;
    bla =  path+name[ii+4]+name[8];
    infile2.open(bla);
    //infile2.open(Form("%s%s.txt",path,name[ii+4]));
    while(infile2.good()) {
      getline(infile2, line);
      //cout<< line <<endl;
      std::istringstream iss(line);
      iss >> _run >> _mean >> _sigma;
      irun = atoi(_run.c_str());
      if(irun==run){
	ZET_mean[ii]  = atof(_mean.c_str());
	ZET_sigma[ii] = atof(_sigma.c_str());
	break;      
      }      
    }//while
    //cout<< "run "<< run << " AOQ" << ii << " " << AOQ_mean[ii] << "+-" << AOQ_sigma[ii] << " ZET" << ii << " " << ZET_mean[ii] << "+-" << ZET_sigma[ii]<<endl;

  }//for

}//ReadPID

void analysis_gamma::ReadPIDCut(){
  TString path = "cutParameters/PID_nsigmacut.txt";
  ifstream infile;
  std::string _nsigma, _comment, line;
  UInt_t ii=0;
  infile.open(path);
  while(infile.good()) {
    getline(infile, line);
    std::istringstream iss(line);
    iss >> _nsigma >> _comment;
    nsigma[ii]=atof(_nsigma.c_str());
    ii++;
  }//while

}

//******* Initialize the tree
void analysis_gamma::SetBranch(){//Switch not needed Branches
  /*fChain->SetBranchStatus("DALINaI.fX",0);
  fChain->SetBranchStatus("DALINaI.fY",0);
  fChain->SetBranchStatus("DALINaI.fZ",0);
  */
  fChain->SetBranchStatus("BigRIPSPPAC*",0);
  fChain->SetBranchStatus("BigRIPSPlastic*",0);
  fChain->SetBranchStatus("BigRIPSIC*",0);
  fChain->SetBranchStatus("BigRIPSFocal*",0);
  fChain->SetBranchStatus("BigRIPSTOF*",0);
  fChain->SetBranchStatus("BigRIPSRIPS*",0);
  fChain->SetBranchStatus("BigRIPSRIPS.delta",1);
  fChain->SetBranchStatus("BigRIPSBeam*",0);
  fChain->SetBranchStatus("BigRIPSBeam.aoq",1);
  fChain->SetBranchStatus("BigRIPSBeam.zet",1);
  
  /*
  fChain->SetBranchStatus("*",1);
  fChain->SetBranchStatus("*",1);
  fChain->SetBranchStatus("*",1);
  fChain->SetBranchStatus("*",1);
  fChain->SetBranchStatus("*",1);
  fChain->SetBranchStatus("*",1);
  fChain->SetBranchStatus("*",1);
  fChain->SetBranchStatus("*",1);
  fChain->SetBranchStatus("*",1);
  fChain->SetBranchStatus("*",1);
  fChain->SetBranchStatus("*",1);
  fChain->SetBranchStatus("*",1);
  fChain->SetBranchStatus("*",1);
  fChain->SetBranchStatus("*",1);
  fChain->SetBranchStatus("*",1);
				*/
}

void analysis_gamma::Tree(UInt_t run){
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
  //if(fChain != 0 ) delete fChain->GetCurrentFile();
  TTree* tree = 0;
  if (tree == 0) {
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(Form("data/rootfiles/new_/run%04d.root",run));
    if (!f || !f->IsOpen()) {
      f = new TFile(Form("data/rootfiles/new_/run%04d.root",run));
    }
    f->GetObject("tree",tree);
    
  }
  Init(tree);
  SetBranch();
  _aoq=2.56;
  if(run<100)_aoq=2.64;
  _zet=50.;
}


void analysis_gamma::MakePID(UInt_t start, UInt_t end,short aoq, short zet){
  FILE* _outaoq;
  _outaoq = fopen(Form("cutParameters/peakpos_AOQ_%i.txt",aoq),"w");
  FILE* _outzet;
  _outzet = fopen(Form("cutParameters/peakpos_ZET_%i.txt",zet),"w");
  TH1F* h1_AOQ[152], *h1_ZET[152];
  TF1* f[2][152];
  TCanvas* can[2][10];
  for(short ii = 0; ii < 10; ii++){
    can[0][ii]=new TCanvas(Form("can0_%i",ii));
    can[0][ii]->Divide(4,4);
    can[1][ii]=new TCanvas(Form("can1_%i",ii));
    can[1][ii]->Divide(4,4);
  }
  UInt_t ican = 0, ipad = 1;
  for(UInt_t run = start; run<end;run++){
    Tree(run);
    h1_AOQ[run] = new TH1F(Form("h1_AOQ_%i",run),"AOQ 0",200,2.5,2.7);
    h1_ZET[run] = new TH1F(Form("h1_ZET_%i",run),"ZET 0",200,45,55);

    fChain->Draw(Form("AOQ[%i]>>h1_AOQ_%i",aoq,run),Form("(AOQ[0]-%f)<0.01 && (AOQ[0]-%f)>-0.01 && (ZET[0]-50)<0.5 && (ZET[0]-50)>-0.5 && (ZET[3]-50)<0.5 && (ZET[3]-50)>-0.5",_aoq,_aoq),"");
    fChain->Draw(Form("ZET[%i]>>h1_ZET_%i",zet,run),Form("(AOQ[0]-%f)<0.01 && (AOQ[0]-%f)>-0.01 && (ZET[0]-50)<0.5 && (ZET[0]-50)>-0.5 && (AOQ[3]-%f)<0.01 && (AOQ[3]-%f)>-0.01",_aoq,_aoq,_aoq,_aoq),"");

    //return;
    f[0][run] = new TF1(Form("f1_%i",run),"gaus",_aoq-0.01,_aoq+0.01);
    f[1][run] = new TF1(Form("f2_%i",run),"gaus",45,55);
    can[0][ican]->cd(ipad);
    h1_AOQ[run]->Draw();
    h1_AOQ[run]->Fit(f[0][run],"R");
    can[1][ican]->cd(ipad);
    h1_ZET[run]->Draw();
    h1_ZET[run]->Fit(f[1][run],"R");
    fprintf(_outaoq,"%i %f %f\n",run, f[0][run]->GetParameter(1), f[0][run]->GetParameter(2));
    fprintf(_outzet,"%i %f %f\n",run, f[1][run]->GetParameter(1), f[1][run]->GetParameter(2));
    if(ipad==16){
      ican++;
      ipad=1;
    }else ipad++;
    //delete h1_AOQ;
    //delete h1_AOQ[1];
    //delete h1_ZET;
    //delete h1_ZET[1];
  }
  fclose(_outaoq);
  fclose(_outzet);
}
