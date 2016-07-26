#include "MakeTree_jojo_Sn128.h"
#include <iostream>
#include "Riostream.h"
#include <string>
#include <time.h>
#include <sstream>
#include "TThread.h"

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TObject.h>
#include <TRandom3.h>

#include "TArtStoreManager.hh"
#include "TArtEventStore.hh"
#include "TClonesArray.h"
#include "TSystem.h"
#include "TTree.h"
#include "TFile.h"

#include "TArtBigRIPSParameters.hh"
#include "TArtDALIParameters.hh"
#include "TArtDALINaIPara.hh"
#include "TArtDALINaI.hh"
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

#include "TVector3.h"
#include "TVector.h"
#include "TMath.h"

#include "signal.h"
/*
//132Sn Empty RUN
#define brho35 6.1721 
#define brho57 5.9055 
#define brho89 5.7055 
#define brhoTE 5.6420 //run391
*/

//132Sn Filled RUN
#define brho35 6.1721 
#define brho57 5.9091
#define brho89 5.1481 
#define brho911 5.140 //guest
#define brhoTE 5.1349 //run31

/*
//128Sn
#define brho35 5.9970 
#define brho57 5.7298 
#define brho89 4.9686 
#define brho911 4.966
#define brhoTE 4.9655 //run381
*/
using namespace std;

//********** main task to call all other function initialize main-loop and so on.
void MakeTree::main(short run, Bool_t test){

  stoploop = false;
  fman = TArtStoreManager::Instance();

  fstore = new TArtEventStore();
  fstore->SetInterrupt(&stoploop); 
  if(run<9)        fstore->Open(Form("data/ridf/Sn/empty_target/labr14Sn%04d.ridf",run));
  else if(run<100) fstore->Open(Form("data/ridf/Sn/Sn132_LHe/labr14Sn%04d.ridf",run));
  else if(run<200) fstore->Open(Form("data/ridf/Sn/Sn128_LHe/labr14Sn%04d.ridf",run));
  else if(run<300) fstore->Open(Form("data/ridf/ggdaq04/labr14Ni%04d.ridf",run-200)); 
  else if(run<400) fstore->Open(Form("data/ridf/ggdaq04/calibration%04d.ridf",run-300)); 
  else             fstore->Open(Form("data/ridf/ggdaq04/calib2XO%04d.ridf",run-400)); 

  Read(run);// read the parameterfiles for corrections and QTC, selfintroduced things
  if( test )fout = new TFile(Form("data/rootfiles/test/run%04d_%3.3f_%3.3f.root",run,ftofoff[0],ftofoff[5]),"RECREATE");
  else fout = new TFile(Form("data/rootfiles/new2/run%04d.root",run),"RECREATE");

  Reset();//Reset all tree and global variables, add treebranches/leaves  
  fTree = new TTree("tree","tree");
  Init(run);//initialize the runvariables such as beam, tof, etc
  AddTree(run);
  Long64_t nev = 0;
  while(fstore->GetNextEvent()){
    //while(estore->GetNextEvent() && neve < end){
    //cout<<nev<<endl;
    nev++;
    if(nev%50000==0){
      cout<< "event:  " << nev <<endl;
      if(test)break;
    }
    if(run>200)
      Calibration(run);
    else{
      Reconstruct(run);//Reconstruct the PID
      Gamma();//Reconstruct the gammas timing(slope-correction), energy for all systems
    }
    fTree -> Fill();
    frawevent->Clear();
    
  }//eventloop
  fout -> Write();
  fout -> Close();
  delete fman;
  fpara->Delete();
  fdpara->Delete();
  delete frecopid;

}
//**********make the timing of the gamma-detectors better
void MakeTree::Gamma(){
  //Float_t clight = 3.0;
  //ftof_F8_T = F8T - 1.198 * clight/fbeam_br_38->GetBeta();//BETA[0]; // time to come from F8 to the target. fc: speed of light
  Float_t mlength = 119.8;//124.8;//cm
  ftof_F8_T = mlength / (clight*BETA[0]); // time to come from F8 to the target. fc: speed of light
  DALI();
  LaBr();
}

//**********DALI
void MakeTree::DALI(){
  /*  for(Int_t i=0; i<150; i++){
    Double_t nan = TMath::QuietNaN();
    DALIfID[i] = -1;      
    DALIfLayer[i] = -1;     
    DALIfx[i] = nan;     
    DALIfy[i] = nan;     
    DALIfz[i] = nan;     
    DALIfd[i] = nan;     
    DALIfTheta[i] = -9999;     
    DALIfTOF[i] = -9999;     
    DALIfRawEnergy[i] = -1;     
    DALIfRawTDC[i] = -1;     
    DALIfEnergy[i] = -1;     
    DALIfDoppCorEnergy[i] = -1;    
    DALIfTOFGEnergy[i] = -1;     
    DALIfTOFGDoppCorEnergy[i] = -1;
  }*/
  fdalicalib ->ClearData();
  fdalicalib -> SetPlTime(fpla->FindPlastic("F7pl")->GetTime());
  //fdalicalib -> SetPlTime(-F8T-ftof_F8_T);
  //fdalicalib -> SetTRef6(tref6);
  //fdalicalib -> SetTRef8(tref8);
  //fdalicalib -> SetMultiplicityThreshold(0);
  //Add above to remove F8plastic tof.
  //fdalicalib -> SetBeta(BETA[0]);
  fdalicalib -> ReconstructData();
  /*for(Int_t i = 0 ; i < fdali_array->GetEntriesFast() ; i++){
    TArtDALINaI * dali = (TArtDALINaI *)fdali_array -> At(i);
    Int_t id = (Int_t)dali -> GetID() - 1;
    
    if(id<150){
      DALIfID[id] = id;
      
      DALIfTOF[id] = test - fADCoffset_t[id];
      DALIfEnergy[id] = (Double_t)dali->GetEnergy();
      DALIfLayer[id] = (Int_t)dali -> GetLayer();
      DALIfx[id] = (Double_t)dali->GetXPos();
      DALIfy[id] = (Double_t)dali->GetYPos();
      DALIfz[id] = (Double_t)dali->GetZPos();
      //DALIfd[id] = sqrt(DALIfx[id]*DALIfx[id]+DALIfy[id]*DALIfy[id]+DALIfz[id]*DALIfz[id]);
      DALIfTheta[id] = (Double_t)dali->GetCosTheta();
      
      //DALIfRawEnergy[id] = (Double_t)dali->GetRawADC();
      //DALIfRawTDC[id] = (Double_t)dali->GetRawTDC();
      //DALIfDoppCorEnergy[id] = (Double_t)dali->GetDoppCorEnergy();
      //DALIfTOFGEnergy[id] = (Double_t)dali->GetTOFGEnergy();
      //DALIfTOFGDoppCorEnergy[id] = (Double_t)dali->GetTOFGDoppCorEnergy();
    }
  }
  
  DALIfMulti = dalicalib -> GetMulti();
  DALIfTOFGMulti = dalicalib -> GetTOFGMulti();
  DALIfEGMulti  = dalicalib -> GetEGMulti();
  DALIfTOFEGMulti  = dalicalib -> GetTOFEGMulti();
  */
}
//**********LaBr
void MakeTree::LaBr(){
  
  TRandom3 rndm;
  for( int i=0; i<fn_labr; i++ ){
    LOWGAIN[i] = -9999;
    MIDGAIN[i] = -9999;
    HIGAIN[i] = -9999;
    QTCTime[i] = -9999;
  }
  Int_t tref = -9000;
  Int_t qtc_dummy[2][130];
  for( int i=0; i<130; i++ ){
    for( int j=0; j<2; j++ ){
      qtc_dummy[j][i] = -9999;
    }
  }
  for(int i=0; i<frawevent -> GetNumSeg(); i++){
    TArtRawSegmentObject *seg = frawevent -> GetSegment(i);
    Int_t fpl = seg -> GetFP();
    if(fpl==8){
      for(int j=0; j < seg -> GetNumData(); j++){
	TArtRawDataObject *d = seg -> GetData(j);
	Int_t geo  = d -> GetGeo();
	Int_t ch   = d -> GetCh();
	Int_t val  = d -> GetVal();
	Int_t edge = d -> GetEdge();
	if(geo==9) qtc_dummy[edge][ch] = val + (rndm.Rndm()- 0.5);
	//if(ch == fQTC[1][1] )cout<<qtc_dummy[edge][ch]<< " qtc_dummy  "<<ch<<endl;
	if(geo ==9 && ch==127 && edge == 0) tref = val;					
      }//for
    }//if
  }//for
  
  for(short ilabr=0; ilabr<fn_labr; ilabr++){
    HIGAIN[ilabr]  = (qtc_dummy[1][ fQTC[0][ilabr] ] - qtc_dummy[0][ fQTC[0][ilabr] ]);//*fQTC_slope[0][ilabr] + fQTC_offset[0][ilabr];
    MIDGAIN[ilabr] = (qtc_dummy[1][ fQTC[1][ilabr] ] - qtc_dummy[0][ fQTC[1][ilabr] ])*fQTC_slope[1][ilabr] + fQTC_offset[1][ilabr];
    LOWGAIN[ilabr] = (qtc_dummy[1][ fQTC[2][ilabr] ] - qtc_dummy[0][ fQTC[2][ilabr] ])*fQTC_slope[2][ilabr] + fQTC_offset[2][ilabr];
    QTCTime[ilabr] = (qtc_dummy[0][fQTC[1][ilabr]] - tref)*ftime - fQTCTimeoffset[ilabr]- ftof_F8_T;

    HIGAINC[ilabr]  = HIGAIN[ilabr]  * (1 - BETA[0] * TMath::Cos(ftheta_labr))/TMath::Sqrt(1-BETA[0] * BETA[0]);
    MIDGAINC[ilabr] = MIDGAIN[ilabr] * (1 - BETA[0] * TMath::Cos(ftheta_labr))/TMath::Sqrt(1-BETA[0] * BETA[0]);
    LOWGAINC[ilabr] = LOWGAIN[ilabr] * (1 - BETA[0] * TMath::Cos(ftheta_labr))/TMath::Sqrt(1-BETA[0] * BETA[0]);
    Double_t corr_time = QTC_Time(qtc_dummy[1][ fQTC[1][ilabr] ] - qtc_dummy[0][ fQTC[1][ilabr] ], LOWGAIN[ilabr],ilabr);
    QTCTimeM[ilabr] = QTCTime[ilabr] - corr_time;
    corr_time = QTC_Time(qtc_dummy[1][ fQTC[1][ilabr] ] - qtc_dummy[0][ fQTC[1][ilabr] ], LOWGAIN[ilabr],ilabr);
    QTCTimeL[ilabr] = QTCTime[ilabr] - corr_time;
  }//for
}

Double_t MakeTree::QTC_Time(Double_t e_mid, Double_t e_low, UInt_t crystal){
 
  Double_t corr = 0.;
  /*if(e_mid < 10e3){
    corr = TMath::Exp(fQTC_slopecorr[crystal][0]*fQTC_slopecorr[crystal][1])+
      fQTC_slopecorr[crystal][2]+fQTC_slopecorr[crystal][3]*e_mid +
      fQTC_slopecorr[crystal][4]/e_mid;
  }else{
    corr = fQTC_slopecorr_low[crystal][0]+
      fQTC_slopecorr_low[crystal][1]*e_low;
  }
*/  
  if(e_mid>0){
    corr = (am[crystal][0]+
	    am[crystal][1]*TMath::Log(TMath::Abs(e_mid-am[crystal][2]))+
	    am[crystal][3]*TMath::Power(TMath::Log(TMath::Abs(e_mid-am[crystal][2])),2)+
	    am[crystal][4]*TMath::Power(TMath::Log(TMath::Abs(e_mid-am[crystal][2])),3));

    corr = corr + (bm[crystal][0]+
		   bm[crystal][1]*e_mid+
		   bm[crystal][2]*e_mid*e_mid+
		   bm[crystal][3]*e_mid*e_mid*e_mid+
		   bm[crystal][4]*e_mid*e_mid*e_mid*e_mid+
		   bm[crystal][5]*e_mid*e_mid*e_mid*e_mid*e_mid+
		   bm[crystal][6]*e_mid*e_mid*e_mid*e_mid*e_mid*e_mid);
  }else{
    corr = (al[crystal][0]+
	    al[crystal][1]*TMath::Log(TMath::Abs(e_mid-al[crystal][2]))+
	    al[crystal][3]*TMath::Power(TMath::Log(TMath::Abs(e_mid-al[crystal][2])),2)+
	    al[crystal][4]*TMath::Power(TMath::Log(TMath::Abs(e_mid-al[crystal][2])),3));

    corr = corr + (bl[crystal][0]+
		   bl[crystal][1]*e_mid+
		   bl[crystal][2]*e_mid*e_mid+
		   bl[crystal][3]*e_mid*e_mid*e_mid+
		   bl[crystal][4]*e_mid*e_mid*e_mid*e_mid+
		   bl[crystal][5]*e_mid*e_mid*e_mid*e_mid*e_mid+
		   bl[crystal][6]*e_mid*e_mid*e_mid*e_mid*e_mid*e_mid);
  }
  return corr;
  
}

//**********Reset the variables after each event (not needed?)
void MakeTree::Reset(){
  for(short ii=0;ii<3;ii++)
    for(short jj=0;jj<6;jj++)
      fVar[ii][jj]=-9999;

  for(short ii=0;ii<9;ii++)
    fVarout3[ii]=-9999;
}
//**********Get the Variables for the PID and the improvements
void MakeTree::GetVar(){
  BETA[0] = fbeam_br -> GetBeta();
  BETA[1] = fbeam_zd_810 -> GetBeta();
  BETA[2] = fbeam_zd_1011 -> GetBeta();
  BETA[3] = fbeam_zd_1011_test -> GetBeta();
  //cout<< "beta: " << BETA[0] << "  " << BETA[1] << "  " << BETA[2] <<endl;
  TArtFocalPlane* tfpl;
  Int_t plane_[6] = {3,5,7,8,10,11};
  for(short ii=0;ii<6;ii++){
    tfpl = ffpl->FindFocalPlane(plane_[ii]);
    if(tfpl){TVectorD* vec=tfpl->GetOptVector(); 
      FX[ii]=(*vec)(0); 
      FA[ii]=(*vec)(1);
      FY[ii]=(*vec)(2);
      FB[ii]=(*vec)(3);
      delta[ii]=(*vec)(4);
    }else cout<< "not found"<<endl;
  }
  F7T = fpla -> FindPlastic("F7pl") -> GetTime();
  F8T = fpla -> FindPlastic("F8pl") -> GetTime();

  //BR 37                    //ZD 810                     //ZD 1011                              		     
  fVar[0][0] = FA[0];        fVar[1][0] = FA[3];          fVar[2][0] = FA[4];
  fVar[0][1] = FX[0];	     fVar[1][1] = FX[3];	  fVar[2][1] = FX[4];
  fVar[0][2] = FY[0];	     fVar[1][1] = FY[3];	  fVar[2][2] = FY[4];
  fVar[0][3] = FA[1];	     fVar[1][3] = FA[4];	  fVar[2][3] = FA[5];
  fVar[0][4] = FX[1];	     fVar[1][4] = FX[4];	  fVar[2][4] = FX[5];
  fVar[0][5] = FY[1];	     fVar[1][5] = FY[4];	  fVar[2][5] = FY[5];

  //ZD 8-11
  fVarout3[0] = FA[3];
  fVarout3[1] = FX[3];
  fVarout3[2] = FY[3];
  fVarout3[3] = FA[4];
  fVarout3[4] = FX[4];
  fVarout3[5] = FY[4];
  fVarout3[6] = FA[5];
  fVarout3[7] = FX[5];
  fVarout3[8] = FY[5];
  fVarout3[9] = BETA[3];

  //tfpl->Delete();
}

//**********Make the PID
void MakeTree::PID(){
  GetVar();
  Float_t aoq[3] = {fbeam_br ->GetAoQ(),  fbeam_zd_810->GetAoQ(), fbeam_zd_1011 ->GetAoQ()};
  for(short iaoq = 0; iaoq<3; iaoq++){
    for(short icorr = 0; icorr < 6; icorr++)
      aoq[iaoq] = aoq[iaoq] - fCorr_aoq[iaoq][1][icorr]*(fVar[iaoq][icorr]/*-fCorr_aoq[iaoq][0][icorr]*/);
  }
  AOQ[0] = aoq[0];
  AOQ[1] = aoq[1];
  AOQ[2] = aoq[2];
  PID_out3();

  ZET[0] = fbeam_br -> GetZet();
  Double_t de_v2 = TMath::Log(fIC2->GetIonPair()*BETA[0]*BETA[0])-TMath::Log(1-BETA[0]*BETA[0])-BETA[0]*BETA[0];
  ZET[1] = fIC2->GetZCoef(0)*TMath::Sqrt((fIC2->GetEnergySqSum()+E_offset_br)/de_v2)*BETA[0]+fIC2->GetZCoef(1);
  ZET[2] = fbeam_zd_1011 -> GetZet();// - (1.17871e+04*BETA[2]+1.10883e+04*BETA[2]*BETA[2]+1.01191e-13*BETA[2]*BETA[2];
  Double_t de_v = TMath::Log(fIC->GetIonPair()*BETA[3]*BETA[3])-TMath::Log(1-BETA[3]*BETA[3])-BETA[3]*BETA[3];
  ZET[3] = fIC->GetZCoef(0)*TMath::Sqrt((fIC->GetEnergySqSum()+fESqoff)/de_v)*BETA[2]+fIC->GetZCoef(1);
  de_v = TMath::Log(fIC->GetIonPair()*BETA[2]*BETA[2])-TMath::Log(1-BETA[2]*BETA[2])-BETA[2]*BETA[2];
  ZET[4] = fIC->GetZCoef(0)*TMath::Sqrt((fIC->GetEnergySqSum()+fESqoff)/de_v)*BETA[3]+fIC->GetZCoef(1);
  ESQ = fIC->GetEnergySqSum();

  // TKE & ASUM 
  TKE = -9999; 
  ASUM = -9999;
	  
  for(int i=0;i<frawevent -> GetNumSeg();i++){
    TArtRawSegmentObject *seg = frawevent -> GetSegment(i);
    Int_t fpl = seg -> GetFP();
    if(fpl==63){
      for(int j=0; j < seg -> GetNumData(); j++){
	TArtRawDataObject *d = seg -> GetData(j);
	Int_t ch = d -> GetCh();
	Int_t geo = d -> GetGeo();
	if(ch==8 && geo==3) TKE = d -> GetVal();
	if(ch==9 && geo==3) ASUM = d -> GetVal();
      }
    }
  }

}
//**********make PID outgoing (test)
void MakeTree::PID_out3(){
  Float_t aoq = fbeam_zd_1011_test ->GetAoQ();
  for(short icorr = 0; icorr < 10; icorr++)
    aoq = aoq - fCorr_aoqout3[1][icorr]*(fVarout3[icorr]/*-fCorr_aoqout3[0][icorr]*/);    

  AOQ[3] = aoq;

}
//**********Initialize the reconstruction and Improve Z-resolution
void MakeTree::Reconstruct(short run){
  Plastic();
  Trigger();
  fpid->ClearData();
  fpid->ReconstructData();
  GetIC(0);
  GetIC(1);
  IC_F11();
  
  // for skip F3 ppac tracking --->
  fpid->ClearData();
  fpid->ReconstructData();    
  if(run>52){
    TArtFocalPlane *fpl = ffpl->FindFocalPlane(3);
    if( fpl ){
      TVectorD *vec=fpl->GetOptVector();
      Double_t vec_tmp[4]={0,0,0,0};
      vec->SetElements(vec_tmp);
    }
    //  fpl->Delete();
  }
  // <-----

  
  frecopid->ClearData();
  frecopid->ReconstructData();

  PID();
}
//**********Improve IC outgoing
void MakeTree::IC_F11(){
  //  GetIC(frips8to10,ftof8to11,"F11IC", fbeam_zd_810,kFALSE);
  TArtFocalPlane *fpl11; TVectorD *vec11;
  Float_t F11X, F11A, F11Y, F11B;
  fpl11 = ffpl -> FindFocalPlane(11);
  if( fpl11 ){vec11 = fpl11 -> GetOptVector();
    F11X = (*vec11)(0); F11A = (*vec11)(1); F11Y = (*vec11)(2); F11B = (*vec11)(3);
  }
  Float_t adc = fIC->GetRawADC(0);
  Float_t corr_a[6] = {5.34924, -6.20299, -5.5168,2.4,-3.4,-5.4094};
  Float_t corr_y[6] = {301.436, -24.1181, 207.03, 53.265, 9.27019, 3.45404};
  Float_t corr_x[6] = {-4.6392, -7.50558, -10.2507, -1.05718, -9.00642, -16.6056};
  for(short ichannel = 0; ichannel<6;ichannel++){
    adc = adc - corr_a[ichannel]*F11A;
    adc = adc - corr_y[ichannel]*F11Y;
	    
    if(TMath::Abs(F11X)<20 && ichannel == 0)
      fIC->SetRawADC(ichannel,adc+800);
    else
      fIC->SetRawADC(ichannel,adc);    
  }
  //  fpl11->Delete();
  //  vec11->Delete();
}

//**********Get IC's to modify them
void MakeTree::GetIC(Bool_t Br){

  TArtRIPS *rips; TArtTOF *tof; const char * icname; TArtBeam* beam;
  if(!Br){
    rips = frips10to11;
    tof = ftof8to11[1];
    icname = "F11IC";
    beam = fbeam_zd_1011;
  }else{//return;
    rips = frips3to5;
    tof = ftof3to7;
    icname = "F7IC";
    beam = fbeam_br2;
  }
  TArtStoreManager * sman = TArtStoreManager::Instance();
  TClonesArray * ic_array = (TClonesArray *)sman->FindDataContainer("BigRIPSIC");
  char name[128];
  Int_t nbeam=1;
  beam->SetDetectorName(name);

  beam->SetNumRIPS(1);
  beam->SetRIPSName(*(rips->GetDetectorName()));
  beam->SetTOFName(*(tof->GetDetectorName()));
  beam->SetICName(icname);

  TArtIC* ic = NULL;
  Int_t num_ic = ic_array->GetEntries();
  for(Int_t j=0;j<num_ic;j++){
    TArtIC *tmp_ic = (TArtIC *)ic_array->At(j);
    TString icname = *(tmp_ic->GetDetectorName());
    if(icname ==  *(beam->GetICName())){
      ic = tmp_ic; 
      break;
    }
  }
  if(NULL == ic){
    TArtCore::Error(__FILE__,"no ic: %s",beam->GetICName()->Data());
  }
  else{
    if(Br)fIC2=ic;
    else fIC=ic;
  }
}

//**********Reading functions for additional parameters
void MakeTree::ReadPID(short run){
  Char_t* path="";
  Char_t *name[4]={"AOQ_Correction_in.txt", "AOQ_Correction_out1.txt", "AOQ_Correction_out2.txt", "AOQ_Correction_out3.txt"};
  string tmp;
  ifstream infile[3], infile3;
  std::string slope, offset, line, nothing;
  if(run<9)          path = "cutParameters/PID/empty/";
  else if(run < 100 )path = "cutParameters/PID/Sn132/";
  else               path = "cutParameters/PID/Sn128/";
  for(short ifile = 0; ifile < 3; ifile++){
    infile[ifile].open(Form("%s%s",path, name[ifile]));
    /* short ivar = 0;
    while(infile[ifile].good()) {
    */
    for(short ivar = 0; ivar < 6; ivar++){
      getline(infile[ifile], line);
      std::istringstream iss(line);
      iss >> offset;// >> slope >> nothing;
      //fCorr_aoq[ifile][0][ivar] = atof(offset.c_str());
      fCorr_aoq[ifile][1][ivar] = atof(offset.c_str());//atof(slope.c_str());
      //ivar++;
      cout<< "Corr aoq " << ifile << "  is and var " << ivar << " is " << fCorr_aoq[ifile][0][ivar] << "  " << fCorr_aoq[ifile][1][ivar] <<endl;
      //if(ivar == 6)break;
    }//while
  }//for

  infile3.open(Form("%s%s",path, name[3]) );
  //while(infile3.good()) {
  if(!infile3.good())cout<< "baeh" <<endl;
  for(short ivar = 0; ivar<10;ivar++){
    getline(infile3, line);
    std::istringstream iss(line);
    iss >> offset;// >> slope >> nothing;
    //fCorr_aoqout3[0][ivar] = atof(offset.c_str());
    fCorr_aoqout3[1][ivar] = atof(offset.c_str());//atof(slope.c_str());
    //ivar++;
    cout<< "Corr aoq " << 3 << "  and var " << ivar << " is " << fCorr_aoqout3[0][ivar] << "  " << fCorr_aoqout3[1][ivar] <<endl;
    if(ivar == 10)break;
  }

}
//**********
void MakeTree::Read(short run){
  ReadTimeoff(run);
  ReadPID(run);
  ReadIC(run);
  ReadQTC();
  ReadDALI();
  ReadNakatsuka();
}
//**********

void MakeTree::ReadNakatsuka(){

  //Nakatsuka
  string tmp;
  ifstream infile1, infile2;
  std::string line;
  std::string b0mt,b1mt,b2mt,b3mt,b4mt,b5mt,b6mt;
  std::string b0lt,b1lt,b2lt,b3lt,b4lt,b5lt,b6lt;
  infile1.open("cutParameters/Nakatsuka1.txt");  
  infile2.open("cutParameters/Nakatsuka2.txt");
  Int_t i=0,crystal=0;

  while(infile1.good()){
    getline(infile1,line);
    std::istringstream iss(line);
    cout<<line<<endl;
    iss >> b0mt >> b1mt >> b2mt >> b3mt >> b4mt >> b5mt >> b6mt;
    am[i][0] = atof(b0mt.c_str());am[i][1] = atof(b1mt.c_str());am[i][2] = atof(b2mt.c_str());am[i][3] = atof(b3mt.c_str());
    am[i][4] = atof(b4mt.c_str());am[i][5] = atof(b5mt.c_str());am[i][6] = atof(b6mt.c_str());
    i++;
    if(i==8)continue;
  }i=0;
  while(infile2.good()){
    getline(infile2,line);
    cout<<line<<endl;
    std::istringstream iss(line);
    iss >> b0lt >> b1lt >> b2lt >> b3lt >> b4lt >> b5lt >> b6lt;
    bm[i][0] = atof(b0lt.c_str());bm[i][1] = atof(b1lt.c_str());bm[i][2] = atof(b2lt.c_str());bm[i][3] = atof(b3lt.c_str());bm[i][4] = atof(b4lt.c_str());
    bm[i][5] = atof(b5lt.c_str());bm[i][6] = atof(b6lt.c_str());
    i++;
    if(i==8)continue;
  }

  infile1.open("cutParameters/lowgainslew.txt");  
  infile2.open("cutParameters/lowgainslew2.txt");
  i=0,crystal=0;

  while(infile1.good()){
    getline(infile1,line);
    std::istringstream iss(line);
    cout<<line<<endl;
    iss >> b0mt >> b1mt >> b2mt >> b3mt >> b4mt >> b5mt >> b6mt;
    al[i][0] = atof(b0mt.c_str());al[i][1] = atof(b1mt.c_str());al[i][2] = atof(b2mt.c_str());al[i][3] = atof(b3mt.c_str());
    al[i][4] = atof(b4mt.c_str());al[i][5] = atof(b5mt.c_str());al[i][6] = atof(b6mt.c_str());
    i++;
    if(i==8)continue;
  }i=0;
  while(infile2.good()){
    getline(infile2,line);
    cout<<line<<endl;
    std::istringstream iss(line);
    iss >> b0lt >> b1lt >> b2lt >> b3lt >> b4lt >> b5lt >> b6lt;
    bl[i][0] = atof(b0lt.c_str());bl[i][1] = atof(b1lt.c_str());bl[i][2] = atof(b2lt.c_str());bl[i][3] = atof(b3lt.c_str());bl[i][4] = atof(b4lt.c_str());
    bl[i][5] = atof(b5lt.c_str());bl[i][6] = atof(b6lt.c_str());
    i++;
    if(i==8)continue;
  }

  for(int i =0; i<fn_labr;i++){
    cout<<i<<endl;
    for(int j =0; j<7;j++)
      cout<< bm[i][j] <<"  "<< bl[i][j] <<"  ";
    cout<<endl;
    for(int j=0; j<5;j++)
      cout<< am[i][j] <<"  "<< al[i][j] <<"  ";
    cout<<endl;
  }

}

//**********
void MakeTree::ReadQTC(){
  
  ftime = 0.09766;
  ftheta_labr = 30./180. * TMath::Pi();
  clight=29.979245;//3.00;
  Float_t qtcslope[3][fn_labr]={1.5194, 1.3509, 1.54495,1.3224,1.51863,1.54124,1.31455,1.43441,
				11.2806,11.1504,12.3476,10.983,12.1714,11.2647,10.1906,10.8408,
				45.6024,46.4729,51.1454,44.8537,47.2729,48.067, 41.5915,43.65};
  Float_t qtcoffset[3][fn_labr]={-2221.76,-1986.75,-2289.73,-2053.94,-2250.28,-2214.31,-1911.24,-2082.49,
				 -16752.7,-16238.2,-17523.1,-16555,  -17515.1,-15875.3,-15031.1,-14837.9,
				 -66695.4,-66298,-70456.5,-66355.7,-67140.3,-69184.8,-61608.1,-59968.6};
  UInt_t hig[fn_labr] = {28, 29, 30, 1, 31, 4, 3, 2};
  UInt_t mig[fn_labr] = {24, 25, 26, 6, 27, 9, 8, 7};
  UInt_t log[fn_labr] = {20, 21, 22,11, 23,14,13,12};//?logic maybe 10,13,12,11??
  Float_t QTCoffset [fn_labr] = {-1037., -1040., -1034., -1030., -1042., -1031., -1031., -1031.};
  Float_t QTCoffset2[fn_labr] = {4.12967e+02,4.16005e+02,4.09709e+02,4.05414e+02,4.17909e+02,4.06827e+02,4.06555e+02,4.06756e+02};//0x0};
  for(short ii=0; ii<8;ii++){
    fQTCTimeoffset[ii]=QTCoffset[ii]+QTCoffset2[ii];
    fQTC[0][ii] = hig[ii];
    fQTC[1][ii] = mig[ii];
    fQTC[2][ii] = log[ii];
    for(short jj=0;jj<3;jj++){
      fQTC_slope[jj][ii] = qtcslope[jj][ii];
      fQTC_offset[jj][ii] = qtcoffset[jj][ii];
    }
  }
  
  //expo+pol1+[5]/x
  string tmp;
  ifstream infile;
  std::string dummy, _parameter, line;
   
  infile.open("cutParameters/QTC2.txt");
  short ii=0,crystal=0;
  while(infile.good()){
    getline(infile, line);
    std::istringstream iss(line);
    iss >> _parameter;// >>;//     iss >> _parameter;
    fQTC_slopecorr[crystal][ii]=atof(_parameter.c_str());     
    cout<< "QTC  " << crystal << "  is  " << fQTC_slopecorr[crystal][ii];
    
    ii++;
    if(ii==5){
      ii=0;
      crystal++;
      cout<<endl;
      if(crystal==8)break;
    }
  }
/*
  //exp0+pol0   
  string tmp;
  ifstream infile;
  std::string dummy, _parameter, line;
   
  infile.open("cutParameters/QTC.txt");
  short ii=0;
  while(infile.good()){
    getline(infile, line);
    std::istringstream iss(line);
    iss >> _parameter;
    fQTC_slopecorr[ii]=atof(_parameter.c_str());     
      cout<< "QTC - slope-correction " << ii << "  is  " << fQTC_slopecorr[ii]<<endl;
      ii++;
    if(ii==3)break;
  }
 
  //pol7 
  string tmp;
  ifstream infile;
  std::string dummy, _parameter, line;
   
  infile.open("cutParameters/QTC.txt");
  short ii=0;
  while(infile.good()){
    getline(infile, line);
    std::istringstream iss(line);
    iss >> _parameter;//     iss >> _parameter;
    fQTC_slopecorr[ii]=atof(_parameter.c_str());     
    cout<< "QTC - slope-correction " << ii << "  is  " << fQTC_slopecorr[ii]<<endl;
    
    ii++;
    if(ii==8)break;
  }
  */
  Double_t _Parameter[2][8] = {-21.1423,-20.9035,-20.569,-20.7755,-19.4699,-18.8751,-19.7029, -20.2248,
			       -0.000179662, -0.000172842, -0.000205383,-0.000243128,-0.000316844,-0.000296179,-0.000302134,-0.000219983};
  /*  Double_t _Parameter[2][8] = {-20.,-21.2585,-20.4934,-21.0495,-20.2589,-19.3573,-17.6622,-21.3258,
      -2e-4,-0.000168967,-0.000209036,-0.000236572,-0.00029798,-0.00027749,-0.00034866,-1.84131e-4};*/
  for(short ii =0;ii<2;ii++)
    for(short jj = 0; jj <fn_labr;jj++){
      fQTC_slopecorr_low[jj][ii] = _Parameter[ii][jj];
      cout<< "QTC - slope-correction low " << ii <<"  "<< jj << "  is  " << fQTC_slopecorr_low[jj][ii]<<endl;
    }
  /*
    vector<std::string>_Parameter;
  _Parameter = vector<std::string>(8);
  infile.open("cutParameters/QTC_low.txt");
  ii=0;
  cout<< "tolle datei!" <<endl;
  while(infile.good()) {
    cout<< "cutParameters/QTC_low.txt" << endl;
    getline(infile, line);
    std::istringstream iss(line);
    iss >> _parameter[0] >> _Parameter[1] >> _Parameter[2] >> _Parameter[3] >> _Parameter[4] >> _Parameter[5] >> _Parameter[6] >> _Parameter[7];
    for(short icrystal = 0; icrystal<fn_labr;icrystal++){
      fQTC_slopecorr_low[icrystal][ii]=atof(_Parameter[icrystal].c_str());     
      cout<< "QTC - slope-correction low " << ii << "  is  " << fQTC_slopecorr_low[icrystal][ii]<<endl;
    }
    ii++;
    if(ii==4)break;
  }
  */
}
//**********
void MakeTree::ReadDALI(){
  string tmp;
  ifstream infile;
  std::string _crystal, offset, line;

  infile.open("cutParameters/DALI_t_offset");
  while(infile.good()) {
    getline(infile, line);
    std::istringstream iss(line);
    iss >> _crystal >> offset;
    fADCoffset_t[atoi(_crystal.c_str())]=atof(offset.c_str());

    //    cout<< "ADC - offset for crystal " << _crystal << "  is  " << fADCoffset_t[atoi(_crystal.c_str())] <<endl;
    if(atoi(_crystal.c_str())==149)break;
  }
}
//**********
void MakeTree::ReadIC(short run){
  string tmp;
  ifstream infile;
  std::string _run, offset, line;
  Float_t soll = 0;
  if(run<9)        
    infile.open("cutParameters/PID/empty/peakpos.txt");
  else  if(run>90)
    infile.open("cutParameters/PID/Sn128/peakpos.txt");
  else
    infile.open("cutParameters/PID/Sn132/peakpos.txt");

  while(infile.good()) {
    getline(infile, line);
    std::istringstream iss(line);
    iss >> _run >> offset;
    if(atoi(_run.c_str())==0){
      soll = atof(offset.c_str());
    }else if(atoi(_run.c_str())==run){
      fESqoff = atof(offset.c_str());
      fESqoff = soll - fESqoff;
      cout<< "ESq - offset for run " << atoi(_run.c_str()) << "  is  " << fESqoff <<endl;
      //break;
    }
  }
  
  //fESqoff = 0;
  Double_t _corr_a[8] = {5.34924, -6.20299, -5.5168,2.4,-3.4,-5.4094};
  Double_t _corr_y[8] = {301.436, -24.1181, 207.03, 53.265, 9.27019, 3.45404};
  Double_t _corr_x[8] = {-4.6392, -7.50558, -10.2507, -1.05718, -9.00642, -16.6056};
  for(short ii=0; ii<6;ii++){
    corr_a[ii] = _corr_a[ii];
    corr_x[ii] = _corr_x[ii];
    corr_y[ii] = _corr_y[ii];
  }
}
//**********
void MakeTree::ReadTimeoff(short run){
  string tmp;
  ifstream infile;
  std::string _aoq, offset, line;
  if(run<9)        
    infile.open("cutParameters/PID/empty/timeoff.txt");
  else  if(run>101)
    infile.open("cutParameters/PID/Sn128/timeoff.txt");
  else  if(run==101)
    infile.open("cutParameters/PID/Sn128/timeoff101.txt");
  else
    infile.open("cutParameters/PID/Sn132/timeoff.txt");
  short aoq=0;
  while(infile.good()) {
    getline(infile, line);
    std::istringstream iss(line);
    iss >> _aoq >> offset;
    if(atoi(_aoq.c_str())==aoq)
      ftofoff[aoq] = atof(offset.c_str());
    cout<< "time - offset for " << aoq << "  is  " << ftofoff[aoq] <<endl;
    aoq++;
    if(aoq == 8)break;
  }  
}
//**********
void MakeTree::LoadPara(){
  fpara -> LoadParameter("db/BigRIPSPPAC.xml");
  fpara -> LoadParameter("db/BigRIPSPlastic.xml");
  fpara -> LoadParameter("db/BigRIPSIC.xml");
  fpara -> LoadParameter("db/FocalPlane.xml");
  fpara -> LoadParameter("db/BigRIPSTKE.xml");
}

//**********Set the timeoffset and window for the TOF
void MakeTree::Plastic(){
  //**** Before doing something first get the time of the plastics
  //Time of the hits in Plastic at F3,F7,F8,F11
  Int_t PlaTDCMAP[8] = {16,17,18,19,20,21,22,23};
  Int_t PlaQMAP[8] = {0,1,2,3,4,5,6,7};
  Int_t PlaTDCRAW[8][100]={0x0};
  Int_t PlaQRAW[8][100];
  Int_t tref63 = -1000;
  for(short ii = 0; ii < fpl_ch; ii++)fplMulti[ii] = 0;
  for(int i = 0; i < frawevent -> GetNumSeg(); i++){
    TArtRawSegmentObject *seg = frawevent -> GetSegment(i);
    Int_t fpl = seg -> GetFP();
    if(fpl==63){
      for(int j=0; j < seg -> GetNumData(); j++){
	TArtRawDataObject *d = seg -> GetData(j);
	Int_t geo  = d -> GetGeo();
	if(geo==6){
	  Int_t ch   = d -> GetCh();
	  Int_t edge = d -> GetEdge();
	  Int_t val  = (Int_t) d -> GetVal();
	  if(edge==0 && ch==31){
	    tref63 = val;
	  }
	  for(Int_t k=0; k<fpl_ch; k++){
	    if(edge==0 && ch == PlaTDCMAP[k]){		      
	      fplMulti[k] = fplMulti[k]+1;
	      PlaTDCRAW[k][j] = val;		 
	    }
	  }
	}
      }
    }
  }

  int *index0; index0 = new int[100]; TMath::Sort(100, PlaTDCRAW[0], index0, true);
  int *index1; index1 = new int[100]; TMath::Sort(100, PlaTDCRAW[1], index1, true);
  int *index2; index2 = new int[100]; TMath::Sort(100, PlaTDCRAW[2], index2, true);
  int *index3; index3 = new int[100]; TMath::Sort(100, PlaTDCRAW[3], index3, true);
  int *index4; index4 = new int[100]; TMath::Sort(100, PlaTDCRAW[4], index4, true);
  int *index5; index5 = new int[100]; TMath::Sort(100, PlaTDCRAW[5], index5, true);
  int *index6; index6 = new int[100]; TMath::Sort(100, PlaTDCRAW[6], index6, true);
  int *index7; index7 = new int[100]; TMath::Sort(100, PlaTDCRAW[7], index7, true);

  Double_t fpltcal = 0.024413;
    for(short ii=0; ii<fpl_ch; ii++){
    F3Pla_TR[ii] =-9999;//(PlaTDCRAW[0][index0[ii]] - tref63)*fpltcal;
    F3Pla_TL[ii] =-9999;//(PlaTDCRAW[1][index1[ii]] - tref63)*fpltcal;
    F7Pla_TR[ii] =-9999;//(PlaTDCRAW[2][index2[ii]] - tref63)*fpltcal;
    F7Pla_TL[ii] =-9999;//(PlaTDCRAW[3][index3[ii]] - tref63)*fpltcal;
    F8Pla_TR[ii] =-9999;//(PlaTDCRAW[4][index4[ii]] - tref63)*fpltcal;
    F8Pla_TL[ii] =-9999;//(PlaTDCRAW[5][index5[ii]] - tref63)*fpltcal;
    F11Pla_TR[ii]=-9999;//(PlaTDCRAW[6][index6[ii]] - tref63)*fpltcal;
    F11Pla_TL[ii]=-9999;//(PlaTDCRAW[7][index7[ii]] - tref63)*fpltcal;
  }
  for(short ii=0; ii<fpl_ch; ii++){

    if(PlaTDCRAW[0][index0[ii]] > 0){
      F3Pla_TR[ii] =(PlaTDCRAW[0][index0[ii]] - tref63)*fpltcal;
    }else F3Pla_TR[ii] = -9999;
    if(PlaTDCRAW[1][index1[ii]] > 0){ 
      F3Pla_TL[ii] =(PlaTDCRAW[1][index1[ii]] - tref63)*fpltcal;
    }else F3Pla_TL[ii] = -9999;
    if(PlaTDCRAW[2][index2[ii]] > 0){
      F7Pla_TR[ii] =(PlaTDCRAW[2][index2[ii]] - tref63)*fpltcal;
    }else F7Pla_TR[ii] = -9999;
    if(PlaTDCRAW[3][index3[ii]] > 0){
      F7Pla_TL[ii] =(PlaTDCRAW[3][index3[ii]] - tref63)*fpltcal;
    }else F7Pla_TL[ii] = -9999;
    if(PlaTDCRAW[4][index4[ii]] > 0){ 
      F8Pla_TR[ii] =(PlaTDCRAW[4][index4[ii]] - tref63)*fpltcal;
    }else F8Pla_TR[ii] = -9999;
    if(PlaTDCRAW[5][index5[ii]] > 0){ 
      F8Pla_TL[ii] =(PlaTDCRAW[5][index5[ii]] - tref63)*fpltcal;
    }else F8Pla_TL[ii] = -9999;
    if(PlaTDCRAW[6][index6[ii]] > 0){
      F11Pla_TR[ii]=(PlaTDCRAW[6][index6[ii]] - tref63)*fpltcal;
    }else F11Pla_TR[ii] = -9999;
    if(PlaTDCRAW[7][index7[ii]] > 0){
      F11Pla_TL[ii]=(PlaTDCRAW[7][index7[ii]] - tref63)*fpltcal;
    }else F11Pla_TL[ii] = -9999;
  }
  delete [] index0;
  delete [] index1;
  delete [] index2;
  delete [] index3;
  delete [] index4;
  delete [] index5;
  delete [] index6;
  delete [] index7;

  //Set the plastic time
  Float_t F3Pla_Time = -9999, F7Pla_Time = -9999, F8Pla_Time = -9999, F11Pla_Time = -9999;

  for(short ii = 0; ii<4;ii++)gate[ii]=0;
  TArtPlastic *pla;

  for(short ii=0;ii<fpl_ch; ii++){
    if( Gated(F3Pla_TL[ii], 3) ){
      if(Gated(F3Pla_TR[ii], 3) ){
	F3Pla_Time = (F3Pla_TR[ii]+F3Pla_TL[ii])/2;
	gate[0]++;
      }else{
	F3Pla_Time = F3Pla_TL[ii];
	gate[0]++;
      }
    }else if(Gated(F3Pla_TR[ii], 3))	F3Pla_Time = F3Pla_TR[ii];
    pla = fpla->FindPlastic((char*)"F3pl");
    pla ->SetTime(F3Pla_Time);
    
    if( Gated(F7Pla_TL[ii],7) ){
      if(Gated(F7Pla_TR[ii],7) ){
	F7Pla_Time = (F7Pla_TR[ii]+F7Pla_TL[ii])/2;
	gate[1]++;
      }else{
	F7Pla_Time = F7Pla_TL[ii];
	gate[1]++;
      }      
    }else if(Gated(F7Pla_TR[ii], 3))	F7Pla_Time = F7Pla_TR[ii];
    pla = fpla->FindPlastic((char*)"F7pl");
    pla ->SetTime(F7Pla_Time);
    
    if( Gated(F8Pla_TL[ii],8) ){
      if(Gated(F8Pla_TR[ii],8) ){
	F8Pla_Time = (F8Pla_TR[ii]+F8Pla_TL[ii])/2;
	gate[2]++;
      }else{
	F8Pla_Time = F8Pla_TL[ii];
	gate[2]++;
      }
    }else if(Gated(F8Pla_TR[ii], 3))	F8Pla_Time = F8Pla_TR[ii];
    pla = fpla->FindPlastic((char*)"F8pl");
    pla ->SetTime(F7Pla_Time);
    
    if( Gated(F11Pla_TL[ii],11) ){
      if(Gated(F11Pla_TR[ii],11) ){
	F11Pla_Time = (F11Pla_TR[ii]+F11Pla_TL[ii])/2;
	gate[3]++;		
      }else{
	F11Pla_Time = F11Pla_TL[ii];
	gate[3]++;
      }
    }else if(Gated(F11Pla_TR[ii], 3))	F11Pla_Time = F7Pla_TR[ii];    
    pla = fpla->FindPlastic((char*)"F11pl-1");
    pla ->SetTime(F11Pla_Time);

  }  
}//void Plastic()

//**********Gates for the TOF
Bool_t MakeTree::Gated(Float_t stamp, short fpl){

  if(fpl == 3){
    if( TMath::Abs(stamp)<5 )return kTRUE;
  }else  if(fpl == 7){
    if( TMath::Abs(stamp)<5 )return kTRUE;
  }else  if(fpl == 8){
    if( TMath::Abs(stamp)<50 )return kTRUE;
  }else  if(fpl == 11){
    if( TMath::Abs(stamp)<5 )return kTRUE;
  }
  else return kFALSE;

}

//**********Give the trigger a number
void MakeTree::Trigger(){

  // Coincidence Register
  TRIGGER = -9999;
  Int_t ftriggerbit = -9999;
  for(int i=0;i<frawevent -> GetNumSeg();i++){
    TArtRawSegmentObject *seg = frawevent -> GetSegment(i);
    Int_t fpl = seg -> GetFP();
    Int_t detector = seg -> GetDetector();
    if(fpl==63 && detector==10){
      for(int j=0; j < seg -> GetNumData(); j++){
	TArtRawDataObject *d = seg -> GetData(j);
	ftriggerbit = d -> GetVal();
      }
    }
  }
  for(Int_t i = 0; i < 7; i++){
    if ((ftriggerbit >> i) & 0x1)
      TRIGGER = i+1;
  }
}

//**********Initialize the arrays of all Detectors
void MakeTree::Array(){
  fppac_array = 
    (TClonesArray *)fman->FindDataContainer("BigRIPSPPAC");
  fTree->Branch(fppac_array->GetName(),&fppac_array);
  fpla_array = 
    (TClonesArray *)fman->FindDataContainer("BigRIPSPlastic");
  fTree->Branch(fpla_array->GetName(),&fpla_array);

  fic_array = 
    (TClonesArray *)fman->FindDataContainer("BigRIPSIC");
  fTree->Branch(fic_array->GetName(),&fic_array);
  ffpl_array = 
    (TClonesArray *)fman->FindDataContainer("BigRIPSFocalPlane");
  fTree->Branch(ffpl_array->GetName(),&ffpl_array);

  //Dali data
  fdali_array=
    (TClonesArray *)fman->FindDataContainer("DALINaI");
  fTree->Branch(fdali_array->GetName(),&fdali_array);
  //PID data
  frips_array = 
    (TClonesArray *)fman->FindDataContainer("BigRIPSRIPS");
  fTree->Branch(frips_array->GetName(),&frips_array); 
  ftof_array  = 
    (TClonesArray *)fman->FindDataContainer("BigRIPSTOF");
  fTree->Branch(ftof_array->GetName(),&ftof_array);   
  fbeam_array = 
    (TClonesArray *)fman->FindDataContainer("BigRIPSBeam");	
  fTree->Branch(fbeam_array->GetName(),&fbeam_array);
}

//**********Initialize for each run the managers
void MakeTree::Init(short run){
  gStyle->SetOptStat(111111);
  gSystem->Load("libXMLParser.so");
  gSystem->Load("libanacore.so");
  gSystem->Load("libanabrips.so");
  gROOT->cd();
  gSystem->AddIncludePath("-I$TARTSYS/include/");
  frawevent = (TArtRawEventObject *)fman->FindDataContainer("RawEvent");
  fpara = TArtBigRIPSParameters::Instance();
  LoadPara();

  fdpara = TArtDALIParameters::Instance();
  fdpara -> LoadParameter("db/DALI.xml");
  fdalicalib = new TArtCalibDALI();  // Create CalibPID to get and calibrate raw data 

  fpid = new TArtCalibPID();
  fppac = fpid -> GetCalibPPAC();
  fpla = fpid -> GetCalibPlastic();
  ffpl = fpid -> GetCalibFocalPlane();

  // Create RecoPID to get calibrated data and to reconstruct TOF, AoQ, Z 
  frecopid = new TArtRecoPID();
  // Definition of observables we want to reconstruct
  frips3to5 = frecopid ->    DefineNewRIPS(3,5,  "matrix/mat1.mat","D3");//brho35);//"D3"); // F3 - F5
  frips5to7 = frecopid ->    DefineNewRIPS(5,7,  "matrix/mat2.mat","D5");//brho57);//"D5"); // F5 - F7
  frips8to10 = frecopid ->   DefineNewRIPS(8,10, "matrix/F8F10_LargeAccAchr.mat","D7");
  //  frips8to9 = frecopid ->    DefineNewRIPS(8,9,  "matrix/F8F9_LargeAccAchr.mat", brho89);
  frips8to9 = frecopid ->    DefineNewRIPS(8,9,"matrix/F8F9_LargeAccAchr.mat",brho89); // F8 - F10  
  frips9to11 = frecopid ->   DefineNewRIPS(9,11, "matrix/F10F11_LargeAccAchr_132Sn.mat",brho911);//.mat.modifiedPieter","D8"); // F8 - F10  
  frips10to11 = frecopid ->  DefineNewRIPS(10,11,"matrix/F10F11_LargeAccAchr_132Sn.mat","D8");//brhoTE);//"D8");
  ftof3to7  = frecopid ->    DefineNewTOF("F3pl","F7pl",ftofoff[0],5); // F3 - F7
  ftof3to8  = frecopid ->    DefineNewTOF("F3pl","F8pl",ftofoff[1],5); // F3 - F8
  ftof3to7_2  = frecopid ->  DefineNewTOF("F3pl","F7pl",ftofoff[2],5); // F3 - F7
  ftof8to11[0] = frecopid -> DefineNewTOF("F8pl","F11pl-1",ftofoff[3],9); // F8 - F11
  ftof8to11[1] = frecopid -> DefineNewTOF("F8pl","F11pl-1",ftofoff[4],9); // F8 - F11
  ftof8to11[2] = frecopid -> DefineNewTOF("F8pl","F11pl-1",ftofoff[5],9); // F8 - F11

  ftof8to11[3] = frecopid -> DefineNewTOF("F8pl","F11pl-1",ftofoff[6],9); // F8 - F11
  ftof8to11[4] = frecopid -> DefineNewTOF("F8pl","F11pl-1",ftofoff[7],9); // F8 - F11

  // Reconstruction of IC observables for ID
  fbeam_br = frecopid -> DefineNewBeam(frips3to5,frips5to7,ftof3to7,"F7IC");
  fbeam_br_38 = frecopid -> DefineNewBeam(frips3to5,frips5to7,ftof3to8,"F7IC");
  fbeam_br2 = frecopid -> DefineNewBeam(frips3to5,frips5to7,ftof3to7_2,"F7IC");
  fbeam_zd_810 = frecopid -> DefineNewBeam(frips8to10,ftof8to11[0],"F11IC");
  fbeam_zd_1011 = frecopid -> DefineNewBeam(frips10to11,ftof8to11[1],"F11IC");
  fbeam_zd_1011_test = frecopid -> DefineNewBeam(frips8to10,frips10to11,ftof8to11[2],"F11IC");
  fbeam_zd_8_9 = frecopid -> DefineNewBeam(frips8to9,ftof8to11[3],"F11IC");
  fbeam_zd_9_11 = frecopid -> DefineNewBeam(frips9to11,ftof8to11[4],"F11IC");
  Array();
}

//**********function to exit loop at keyboard interrupt. 
void MakeTree::stop_interrupt(){
  printf("keyboard interrupt\n");
  stoploop = true;
}

void MakeTree::Calibration(UInt_t run){

  //******************labr  
  TRandom3 rndm;
  for( int i=0; i<fn_labr; i++ ){
    LOWGAIN[i] = -9999;
    MIDGAIN[i] = -9999;
    HIGAIN[i] = -9999;
    QTCTime[i] = -9999;
  }
  Int_t tref = -9000;
  Int_t qtc_dummy[2][130];
  for( int i=0; i<130; i++ ){
    for( int j=0; j<2; j++ ){
      qtc_dummy[j][i] = -9999;
    }
  }
  for(int i=0; i<frawevent -> GetNumSeg(); i++){
    TArtRawSegmentObject *seg = frawevent -> GetSegment(i);
    Int_t fpl = seg -> GetFP();
    if(fpl==8){
      for(int j=0; j < seg -> GetNumData(); j++){
	TArtRawDataObject *d = seg -> GetData(j);
	Int_t geo  = d -> GetGeo();
	Int_t ch   = d -> GetCh();
	Int_t val  = d -> GetVal();
	Int_t edge = d -> GetEdge();
	if(geo==9) qtc_dummy[edge][ch] = val + (rndm.Rndm()- 0.5);
	//if(ch == fQTC[1][1] )cout<<qtc_dummy[edge][ch]<< " qtc_dummy  "<<ch<<endl;
	if(geo ==9 && ch==127 && edge == 0) tref = val;					
      }//for
    }//if
  }//for
  
  for(short ilabr=0; ilabr<fn_labr; ilabr++){
    HIGAIN[ilabr]  = (qtc_dummy[1][ fQTC[0][ilabr] ] - qtc_dummy[0][ fQTC[0][ilabr] ]);
    MIDGAIN[ilabr] = (qtc_dummy[1][ fQTC[1][ilabr] ] - qtc_dummy[0][ fQTC[1][ilabr] ]);
    LOWGAIN[ilabr] = (qtc_dummy[1][ fQTC[2][ilabr] ] - qtc_dummy[0][ fQTC[2][ilabr] ]);
    QTCTime[ilabr] = (qtc_dummy[0][ fQTC[1][ilabr] ] - tref)*ftime - fQTCTimeoffset[ilabr];
  }
  //*************DALI

  fdalicalib->ClearData();
  fdalicalib->ReconstructData();
}//Calibration


//**********Add the variables to the tree
void MakeTree::AddTree(UInt_t run){
  if(run<200){
    fTree->Branch("AOQ",AOQ,"AOQ[4]/D");
    fTree->Branch("ZET",ZET,"ZET[5]/D");
    fTree->Branch("FX",FX,"FX[6]/D");
    fTree->Branch("FY",FY,"FY[6]/D");
    fTree->Branch("FA",FA,"FA[6]/D");
    fTree->Branch("FB",FB,"FB[6]/D");
    fTree->Branch("delta",delta,"delta[6]/D");
    fTree->Branch("F7T",&F7T,"F7T/D");
    fTree->Branch("F8T",&F8T,"F8T/D");
    fTree->Branch("BETA",BETA,"BETA[4]/D");
    fTree->Branch("gate",gate,"gate[4]/I");
    fTree->Branch("TRIGGER",&TRIGGER,"TRIGGER/I");        
    fTree->Branch("fplMulti",fplMulti,Form("fplMulti[%i]/D",fpl_ch));
    fTree->Branch("ESQ",&ESQ,"ESQ/F");

    //PID
    //Plastic
    fTree->Branch("F3Pla_TR", F3Pla_TR, "F3Pla_TR[8]/F"); 
    fTree->Branch("F3Pla_TL", F3Pla_TL, "F3Pla_TL[8]/F"); 
    fTree->Branch("F7Pla_TR", F7Pla_TR, "F7Pla_TR[8]/F"); 
    fTree->Branch("F7Pla_TL", F7Pla_TL, "F7Pla_TL[8]/F"); 
    fTree->Branch("F8Pla_TR", F8Pla_TR, "F8Pla_TR[8]/F"); 
    fTree->Branch("F8Pla_TL", F8Pla_TL, "F8Pla_TL[8]/F"); 
    fTree->Branch("F11Pla_TR",F11Pla_TR,"F11Pla_TR[8]/F");
    fTree->Branch("F11Pla_TL",F11Pla_TL,"F11Pla_TL[8]/F");

    fTree->Branch("TKE",&TKE,"TKE/D");
    fTree->Branch("ASUM",&ASUM,"ASUM/D");
	
  }
  //*******Gammas
  
  //LaBr
  //Raw

  fTree->Branch("HIGAIN",HIGAIN,"HIGAIN[8]/F");        
  fTree->Branch("MIDGAIN",MIDGAIN,"MIDGAIN[8]/F");        
  fTree->Branch("LOWGAIN",LOWGAIN,"LOWGAIN[8]/F");        
  fTree->Branch("QTCTime",QTCTime,"QTCTime[8]/F");        
  //Doppler corrected and slope-correction
  fTree->Branch("HIGAINC",HIGAINC,"HIGAINC[8]/F");        
  fTree->Branch("MIDGAINC",MIDGAINC,"MIDGAINC[8]/F");        
  fTree->Branch("LOWGAINC",LOWGAINC,"LOWGAINC[8]/F");        
  fTree->Branch("QTCTimeM",QTCTimeM,"QTCTimeM[8]/F");        
  fTree->Branch("QTCTimeL",QTCTimeL,"QTCTimeL[8]/F");        

  //DALI
  /*
  fTree->Branch("DALIfID",DALIfID,"DALIfID[150]/I");	
  fTree->Branch("DALIfLayer",DALIfLayer,"DALIfLayer[150]/I");	
  fTree->Branch("DALIfx",DALIfx,"DALIfx[150]/F");	
  fTree->Branch("DALIfy",DALIfy,"DALIfy[150]/F");	
  fTree->Branch("DALIfz",DALIfz,"DALIfz[150]/F");
  fTree->Branch("DALIfd",DALIfd,"DALIfd[150]/F"); 
  fTree->Branch("DALIfTheta",DALIfTheta,"DALIfTheta[150]/F");	      
  fTree->Branch("DALIfTOF",DALIfTOF,"DALIfTOF[150]/F");	      
  fTree->Branch("DALIfRawEnergy",DALIfRawEnergy,"DALIfRawEnergy[150]/F");	      
  fTree->Branch("DALIfRawTDC",DALIfRawTDC,"DALIfRawTDC[150]/F");	      
  fTree->Branch("DALIfEnergy",DALIfEnergy,"DALIfEnergy[150]/F");           
  fTree->Branch("DALIfDoppCorEnergy",DALIfDoppCorEnergy,"DALIfDoppCorEnergy[150]/F");    
  fTree->Branch("DALIfTOFGEnergy",DALIfTOFGEnergy,"DALIfTOFGEnergy[150]/F");       
  fTree->Branch("DALIfTOFGDoppCorEnergy",DALIfTOFGDoppCorEnergy,"DALIfTOFGDoppCorEnergy[150]/F");
  */


}  //******end
