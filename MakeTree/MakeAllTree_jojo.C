#include "MakeAllTree_jojo.h"
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
#define brhoTE 5.1349 //run31

/*
//128Sn
#define brho35 5.9970 
#define brho57 5.7298 
#define brho89 1
#define brhoTE 1 //run381
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
  else             fstore->Open(Form("data/ridf/Sn/Sn128_LHe/labr14Sn%04d.ridf",run));
  Read(run);// read the parameterfiles for corrections and QTC, selfintroduced things
  if( test )fout = new TFile(Form("data/rootfiles/test/run%04d_%3.3f_%3.3f.root",run,ftofoff[0],ftofoff[3]),"RECREATE");
  else fout = new TFile(Form("data/rootfiles/new/run%04d.root",run),"RECREATE");

  Reset();//Reset all tree and global variables, add treebranches/leaves  
  fTree = new TTree("tree","tree");
  Init(run);//initialize the runvariables such as beam, tof, etc
  AddTree();
  Long64_t nev = 0;
  while(fstore->GetNextEvent()){
    //while(estore->GetNextEvent() && neve < end){
    nev++;
    if(nev%50000==0){
      cout<< "event:  " << nev <<endl;
      if(test)break;
    }
    //cout << " beta 1:  " << fbeam_br->GetBeta() << endl;
    Reconstruct(run);//Reconstruct the PID
    //Gamma();//Reconstruct the gammas timing(slope-correction), energy for all systems
    fTree -> Fill();
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
  Float_t clight = 3.0;
  ftof_F8_T = F8T - 1.198 * clight/BETA[3]; // time to come from F8 to the target. fc: speed of light
  DALI();
  LaBr();
}
//**********DALI
void MakeTree::DALI(){
  for(Int_t i=0; i<150; i++){
    Double_t nan = TMath::QuietNaN();
    DALIfID[i] = -1;      
    DALIfLayer[i] = -1;     
    DALIfx[i] = nan;     
    DALIfy[i] = nan;     
    DALIfz[i] = nan;     
    DALIfd[i] = nan;     
    DALIfTheta[i] = nan;     
    DALIfTOF[i] = nan;     
    DALIfRawEnergy[i] = -1;     
    DALIfRawTDC[i] = -1;     
    DALIfEnergy[i] = -1;     
    DALIfDoppCorEnergy[i] = -1;    
    DALIfTOFGEnergy[i] = -1;     
    DALIfTOFGDoppCorEnergy[i] = -1;
  }
  
  fdalicalib -> SetPlTime(ftof_F8_T);
  //fdalicalib -> SetTRef6(tref6);
  //fdalicalib -> SetTRef8(tref8);
  //fdalicalib -> SetMultiplicityThreshold(0);
  //cout << "fPlTOF in MAKE = " << tofoffg  << endl;
  //Add above to remove F8plastic tof.
  fdalicalib -> ReconstructData();
  
  for(Int_t i = 0 ; i < dali_array->GetEntriesFast() ; i++){
    TArtDALINaI * dali = (TArtDALINaI *)dali_array -> At(i);
    Int_t id = (Int_t)dali -> GetID() - 1;
    DALIfID[id] = (Int_t)dali->GetID() -1;
    DALIfLayer[id] = (Int_t)dali -> GetLayer();
    DALIfx[id] = (Double_t)dali->GetXPos();
    DALIfy[id] = (Double_t)dali->GetYPos();
    DALIfz[id] = (Double_t)dali->GetZPos();
    DALIfTheta[id] = (Double_t)dali->GetTheta();
    DALIfTOF[id] = (Double_t)dali->GetTOF();
    DALIfRawEnergy[id] = (Double_t)dali->GetRawADC();
    DALIfRawTDC[id] = (Double_t)dali->GetRawTDC();
    DALIfEnergy[id] = (Double_t)dali->GetEnergy();
    DALIfDoppCorEnergy[id] = (Double_t)dali->GetDoppCorEnergy();
    DALIfTOFGEnergy[id] = (Double_t)dali->GetTOFGEnergy();
    DALIfTOFGDoppCorEnergy[id] = (Double_t)dali->GetTOFGDoppCorEnergy(); 
  }
  
  DALIfMulti = dalicalib -> GetMulti();
  DALIfTOFGMulti = dalicalib -> GetTOFGMulti();
  DALIfEGMulti  = dalicalib -> GetEGMulti();
  DALIfTOFEGMulti  = dalicalib -> GetTOFEGMulti();
  
}
//**********LaBr
void MakeTree::LaBr(){

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
  ZET[2] = fbeam_zd_1011 -> GetZet();
  Double_t de_v = TMath::Log(fIC->GetIonPair()*BETA[3]*BETA[3])-TMath::Log(1-BETA[3]*BETA[3])-BETA[3]*BETA[3];
  ZET[3] = fIC->GetZCoef(0)*TMath::Sqrt((fIC->GetEnergySqSum()+fESqoff)/de_v)*BETA[2]+fIC->GetZCoef(1);
  de_v = TMath::Log(fIC->GetIonPair()*BETA[2]*BETA[2])-TMath::Log(1-BETA[2]*BETA[2])-BETA[2]*BETA[2];
  ZET[4] = fIC->GetZCoef(0)*TMath::Sqrt((fIC->GetEnergySqSum()+fESqoff)/de_v)*BETA[3]+fIC->GetZCoef(1);
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
  if(run>52){
    TArtFocalPlane *fpl = ffpl->FindFocalPlane(3);
    if( fpl ){
      TVectorD *vec=fpl->GetOptVector();
      Double_t vec_tmp[4]={0,0,0,0};
      vec->SetElements(vec_tmp);
    }
  }
  // <-----
 
  fpid->ClearData();
  fpid->ReconstructData();  
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
}

//**********Get IC's to modify them
void MakeTree::GetIC(Bool_t Br){

  TArtRIPS *rips; TArtTOF *tof; const char * icname; TArtBeam* beam;
  if(!Br){
    rips = frips10to11;
    tof = ftof8to11[1];
    icname = "F11IC";
    beam = fbeam_zd_1011;
  }else{
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
      iss >> offset >> slope >> nothing;
      fCorr_aoq[ifile][0][ivar] = atof(offset.c_str());
      fCorr_aoq[ifile][1][ivar] = atof(slope.c_str());
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
    iss >> offset >> slope >> nothing;
    fCorr_aoqout3[0][ivar] = atof(offset.c_str());
    fCorr_aoqout3[1][ivar] = atof(slope.c_str());
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
}
//**********
void MakeTree::ReadQTC(){
  
  Float_t QTCoffset[8] = {-1037., -1040., -1034., -1030., -1042., -1031., -1031., -1031.};
  for(short ii=0; ii<8;ii++)
    fQTCoffset[ii]=QTCoffset[ii];
}
//**********
void MakeTree::ReadIC(short run){
  string tmp;
  ifstream infile;
  std::string _run, offset, line;

  if(run<9)        
    infile.open("cutParameters/PID/empty/peakpos_ESq.txt");
  else  if(run>90)
    infile.open("cutParameters/PID/Sn128/peakpos_ESq.txt");
  else
    infile.open("cutParameters/PID/Sn132/peakpos_ESq.txt");

  while(infile.good()) {
    getline(infile, line);
    std::istringstream iss(line);
    iss >> _run >> offset;
    if(atoi(_run.c_str())==run){
      fESqoff = atof(offset.c_str());
      cout<< "ESq - offset for run " << run << "  is  " << fESqoff <<endl;
    }
  }
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
  else  if(run>90)
    infile.open("cutParameters/PID/Sn128/timeoff.txt");
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
    if(aoq == 4)break;
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
  Double_t F3Pla_TR[fpl_ch],F3Pla_TL[fpl_ch],
    F7Pla_TR[fpl_ch], F7Pla_TL[fpl_ch],
    F8Pla_TR[fpl_ch], F8Pla_TL[fpl_ch],
    F11Pla_TR[fpl_ch],F11Pla_TL[fpl_ch];
  for(short ii=0; ii<fpl_ch; ii++){
    F3Pla_TR[ii] =0;//(PlaTDCRAW[0][index0[ii]] - tref63)*fpltcal;
    F3Pla_TL[ii] =0;//(PlaTDCRAW[1][index1[ii]] - tref63)*fpltcal;
    F7Pla_TR[ii] =0;//(PlaTDCRAW[2][index2[ii]] - tref63)*fpltcal;
    F7Pla_TL[ii] =0;//(PlaTDCRAW[3][index3[ii]] - tref63)*fpltcal;
    F8Pla_TR[ii] =0;//(PlaTDCRAW[4][index4[ii]] - tref63)*fpltcal;
    F8Pla_TL[ii] =0;//(PlaTDCRAW[5][index5[ii]] - tref63)*fpltcal;
    F11Pla_TR[ii]=0;//(PlaTDCRAW[6][index6[ii]] - tref63)*fpltcal;
    F11Pla_TL[ii]=0;//(PlaTDCRAW[7][index7[ii]] - tref63)*fpltcal;
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
	F3Pla_Time = (F3Pla_TR[ii]+F3Pla_TL[ii]);
	gate[0]++;
      }else{
	F3Pla_Time = F3Pla_TL[ii]/2;
	gate[0]++;
      }
      pla = fpla->FindPlastic((char*)"F3pl");
      pla ->SetTime(F3Pla_Time);
    }
    if( Gated(F7Pla_TL[ii],7) ){
      if(Gated(F7Pla_TR[ii],7) ){
	F7Pla_Time = (F7Pla_TR[ii]+F7Pla_TL[ii]);
	gate[1]++;
      }else{
	F7Pla_Time = F7Pla_TL[ii];
	gate[1]++;
      }
      pla = fpla->FindPlastic((char*)"F7pl");
      pla ->SetTime(F7Pla_Time);
    }
    if( Gated(F8Pla_TL[ii],8) ){
      if(Gated(F8Pla_TR[ii],8) ){
	F8Pla_Time = (F8Pla_TR[ii]+F8Pla_TL[ii]);
	gate[2]++;
      }else{
	F8Pla_Time = F8Pla_TL[ii];
	gate[2]++;
      }
      pla = fpla->FindPlastic((char*)"F7pl");
      pla ->SetTime(F7Pla_Time);
    }
    if( Gated(F11Pla_TL[ii],11) ){
      if(Gated(F11Pla_TR[ii],11) ){
	F11Pla_Time = (F11Pla_TR[ii]+F11Pla_TL[ii]);
	gate[3]++;		
      }else{
	F11Pla_Time = F11Pla_TL[ii];
	gate[3]++;
      }
      pla = fpla->FindPlastic((char*)"F11pl-1");
      pla ->SetTime(F11Pla_Time);
    }
  }

}

//**********Gates for the TOF
Bool_t MakeTree::Gated(Float_t stamp, short fpl){

  if(fpl == 3){
    if( TMath::Abs(stamp+560)<20 )return kTRUE;
  }else  if(fpl == 7){
    if( TMath::Abs(stamp+600)<20 )return kTRUE;
  }else  if(fpl == 8){
    if( TMath::Abs(stamp+610)<20 )return kTRUE;
  }else  if(fpl == 11){
    if( TMath::Abs(stamp+225)<20 )return kTRUE;
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

  // Create BigRIPSParameters to get Plastics, PPACs, ICs and FocalPlanes parameters from ".xml" files
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
  frips3to5 = frecopid -> DefineNewRIPS(3,5,"matrix/mat1.mat","D3");//brho35);//"D3"); // F3 - F5
  frips5to7 = frecopid -> DefineNewRIPS(5,7,"matrix/mat2.mat","D5");//brho57);//"D5"); // F5 - F7
  frips8to10 = frecopid -> DefineNewRIPS(8,10,"matrix/F8F10_LargeAccAchr.mat","D7"); // F8 - F10  
  frips8to9 = frecopid -> DefineNewRIPS(8,9,"matrix/F8F9_LargeAccAchr.mat",brho89); // F8 - F10  
  frips10to11 = frecopid -> DefineNewRIPS(10,11,"matrix/F10F11_LargeAccAchr_132Sn.mat","D8");//brhoTE);//"D8");
  ftof3to7  = frecopid -> DefineNewTOF("F3pl","F7pl",ftofoff[0],5); // F3 - F7
  ftof8to11[0] = frecopid -> DefineNewTOF("F8pl","F11pl-1",ftofoff[1],9); // F8 - F11
  ftof8to11[1] = frecopid -> DefineNewTOF("F8pl","F11pl-1",ftofoff[2],9); // F8 - F11
  ftof8to11[2] = frecopid -> DefineNewTOF("F8pl","F11pl-1",ftofoff[3],9); // F8 - F11

  // Reconstruction of IC observables for ID
  fbeam_br = frecopid -> DefineNewBeam(frips3to5,frips5to7,ftof3to7,"F7IC");
  fbeam_br2 = frecopid -> DefineNewBeam(frips3to5,frips5to7,ftof3to7,"F7IC");
  fbeam_zd_810 = frecopid -> DefineNewBeam(frips8to10,ftof8to11[0],"F11IC");
  fbeam_zd_1011 = frecopid -> DefineNewBeam(frips10to11,ftof8to11[1],"F11IC");
  fbeam_zd_1011_test = frecopid -> DefineNewBeam(frips8to10,frips10to11,ftof8to11[2],"F11IC");

  Array();
}

//**********function to exit loop at keyboard interrupt. 
void MakeTree::stop_interrupt(){
  printf("keyboard interrupt\n");
  stoploop = true;
}

//**********Add the variables to the tree
void MakeTree::AddTree(){
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
  fTree->Branch("fplMulti",fplMulti,Form("fplMulti[%s]/I",fpl_ch));	

  fTree->Branch("QTCTime",QTCTime,Form("QTCTime[%i]/F",fpl_ch));	
  fTree->Branch("QTCTimeC",QTCTimeC,Form("QTCTime[%i]/F",fn_labr));
  fTree->Branch("LOWGAIN",LOWGAIN,Form("LOWGAIN[%i]/F",fpl_ch));	
  fTree->Branch("LOWGAINC",LOWGAINC,Form("LOWGAIN[%i]/F",fn_labr));
  fTree->Branch("MIDGAIN",MIDGAIN,Form("MIDGAIN[%i]/F",fn_labr));	
  fTree->Branch("MIDGAINC",MIDGAINC,Form("MIDGAINC[%i]/F",fn_labr));	
  fTree->Branch("HIGAIN",HIGAIN,Form("HIGAIN[%i]/F",fpl_ch));	
  fTree->Branch("HIGAINC",HIGAINC,Form("HIGAINC[%i]/F",fn_labr));
	
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
}
