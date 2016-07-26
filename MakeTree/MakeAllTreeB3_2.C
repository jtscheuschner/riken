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

//Sn-128
//#define tofoffbr  301.// TOF3-7 offset VME
//#define tofoffzd -156.5 // TOF8-11 offset VME

//Sn-132
// BigRIPS parameters
//#define tofoffbr  299.8// TOF3-7 offset VME
//#define tofoffzd  -158.1 // TOF7-11 offset VME

/*
//132Sn Empty RUN
#define brho35 6.1721 
#define brho57 5.9055 
#define brho89 5.7055 
#define brhoTE 5.6420 //run391
*/
/*
//132Sn Filled RUN
#define brho35 6.1721 
#define brho57 5.9091
#define brho89 5.1481 
#define brhoTE 5.1349 //run31

**/
//128Sn
#define brho35 5.9970 
#define brho57 5.7298 
#define brho89 1
#define brhoTE 1 //run381

using namespace std;

TArtIC* fIC;
TArtIC* fIC2;


// function to exit loop at keyboard interrupt. 
bool stoploop = false;
void stop_interrupt(){
  printf("keyboard interrupt\n");
  stoploop = true;
}

Bool_t gated(Float_t stamp, short fpl){

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

void GetIC(TArtRIPS *rips, TArtTOF *tof, const char * icname, TArtBeam* beam,Bool_t Br){

  TArtStoreManager * sman = TArtStoreManager::Instance();
  TClonesArray * ic_array = (TClonesArray *)sman->FindDataContainer("BigRIPSIC");
  char name[128];
  Int_t nbeam=1;//GetNumBeam();
  //new ((*fBeamArray)[nbeam]) TArtBeam();
  //  TArtBeam * beam = (TArtBeam *)beam;

  //  sprintf(name,"Beam_rips%s_tof%s_ic%s",rips->GetDetectorName()->Data(),tof->GetDetectorName()->Data(),icname);
  //TArtCore::Info(__FILE__,"define %s",name);
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
      //cout<< icname << endl;
      ic = tmp_ic; 
      break;
    }
  }
  if(NULL == ic){
    TArtCore::Error(__FILE__,"no ic: %s",beam->GetICName()->Data());
  }
  else{
    //    fUpstreamRIPSArrayBuffer.push_back(rips);
    //fDownstreamRIPSArrayBuffer.push_back(NULL);
    //fTOFArrayBuffer.push_back(tof);
    if(Br)fIC2=ic;
    else fIC=ic;
  }

}

//Calibration parameters
Double_t QTCTimeOffset[]=
  {-1037., -1040., -1034., -1030., -1042., -1031., -1031., -1031.};

Float_t ESq_offset(Int_t run){

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

  /*for(Int_t ii = 0; ii<55; ii++){
    for(Int_t jj=0; jj<2;jj++){
      infile>>peak[jj];
      cout<< peak[jj] <<endl;
    }
    if(peak[0] == run){
      offset = mean - peak[1];
      break;
    }  
  }
*/
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

Float_t Offset(Int_t run, Int_t aoq){

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

void MakeAllTreeB3_2(Int_t run, Bool_t test){//-158.39


  const Int_t fpl_ch = 8;
  Float_t E_offset=ESq_offset(run);
  Int_t E_offset_br=0;

  Float_t tofoffbr  = Offset(run,0);
  Float_t tofoffzd1 = Offset(run,1);
  Float_t tofoffzd2 = Offset(run,2);
  Float_t tofoffzd3 = Offset(run,3);
  Int_t TimeWindow[2][2][4];/*
  TimeWindow[0][0][0] = 3000; TimeWindow[0][1][0] = 3500;//F3Left
  TimeWindow[1][0][0] = 3000; TimeWindow[1][1][0] = 3500;//F3right
  TimeWindow[0][0][1] = 3200; TimeWindow[0][1][1] = 4000;//F7
  TimeWindow[1][0][1] = 3000; TimeWindow[1][1][1] = 4000;//F7
  TimeWindow[0][0][2] = TimeWindow[0][1][2] = 
  TimeWindow[1][0][1] = TimeWindow[1][1][2] = ;
  TimeWindow[0][0][3] = TimeWindow[0][1][3] = ;//F11
  TimeWindow[1][0][3] = TimeWindow[1][1][3] = 
			    */

  TRandom3 rndm;
	gStyle->SetOptStat(111111);
	gSystem->Load("libXMLParser.so");
	gSystem->Load("libanacore.so");
	gSystem->Load("libanabrips.so");
	gROOT->cd();
	//	cout<<$TARTSYS<<endl;
	gSystem->AddIncludePath("-I$TARTSYS/include/");

	// Create StoreManager both for calibration "TArtCalib..." and treatment "TArtReco..."
	TArtStoreManager * sman = TArtStoreManager::Instance();

	// Create EventStore to control the loop and get the EventInfo
	TArtEventStore *estore = new TArtEventStore();
	estore->SetInterrupt(&stoploop); 
	if(run<9)        estore->Open(Form("data/ridf/Sn/empty_target/labr14Sn%04d.ridf",run));
	else if(run<100) estore->Open(Form("data/ridf/Sn/Sn132_LHe/labr14Sn%04d.ridf",run));
	else{  
	  estore->Open(Form("data/ridf/Sn/Sn128_LHe/labr14Sn%04d.ridf",run));
	  E_offset_br = -14;
	  /*  tofoffbr = 301.526;
	      tofoffzd = -156.264;
	      tofoffzd2 = -156.264;*/
	}
	//	estore->Open(Form("ridf/ggdaq04/calibration%04d.ridf",run));
	//estore->Open(Form("ridf/bigrips/run%04d.ridf",run));
	//estore->Open(0);

	TArtRawEventObject *rawevent = (TArtRawEventObject *)sman->FindDataContainer("RawEvent");
	
	// Create BigRIPSParameters to get Plastics, PPACs, ICs and FocalPlanes parameters from ".xml" files
	TArtBigRIPSParameters *para = TArtBigRIPSParameters::Instance();
	para -> LoadParameter("db/BigRIPSPPAC.xml");
	para -> LoadParameter("db/BigRIPSPlastic.xml");
	para -> LoadParameter("db/BigRIPSIC.xml");
	para -> LoadParameter("db/FocalPlane.xml");
	para -> LoadParameter("db/BigRIPSTKE.xml");

	//Create DALI Parameters
	TArtDALIParameters *dpara = TArtDALIParameters::Instance();
	dpara -> LoadParameter("db/DALI.xml");
	TArtCalibDALI *dalicalib = new TArtCalibDALI();

	// Create CalibPID to get and calibrate raw data 
	TArtCalibPID *cpid = new TArtCalibPID();
	TArtCalibPPAC *cppac = cpid -> GetCalibPPAC();
	TArtCalibPlastic *cpla = cpid -> GetCalibPlastic();
	TArtCalibFocalPlane *cfpl = cpid -> GetCalibFocalPlane();
	TArtCalibIC* icCalib;

	// Create RecoPID to get calibrated data and to reconstruct TOF, AoQ, Z 
	TArtRecoPID *recopid = new TArtRecoPID();

	// Definition of observables we want to reconstruct
	TArtRIPS *rips3to5 = recopid -> DefineNewRIPS(3,5,"matrix/mat1.mat","D3");//brho35);//"D3"); // F3 - F5
	TArtRIPS *rips5to7 = recopid -> DefineNewRIPS(5,7,"matrix/mat2.mat","D5");//brho57);//"D5"); // F5 - F7
	TArtRIPS *rips8to10 = recopid -> DefineNewRIPS(8,10,"matrix/F8F10_LargeAccAchr.mat","D7"); // F8 - F10  
	TArtRIPS *rips8to9 = recopid -> DefineNewRIPS(8,9,"matrix/F8F9_LargeAccAchr.mat",brho89); // F8 - F10  
	TArtRIPS *rips10to11 = recopid -> DefineNewRIPS(10,11,"matrix/F10F11_LargeAccAchr_132Sn.mat","D8");//brhoTE);//"D8"); // F9 - F11
	//TArtRIPS *rips10to11 = recopid -> DefineNewRIPS(10,11,"matrix/F10F11_LargeAccAchr.mat",brhoTE);//"D8"); // F9 - F11

	// Reconstruction of TOF DefineNewTOF(fisrt plane, second plane, time offset)
	TArtTOF *tof3to7  = recopid -> DefineNewTOF("F3pl","F7pl",tofoffbr,5); // F3 - F7
	TArtTOF *tof8to11_1 = recopid -> DefineNewTOF("F8pl","F11pl-1",tofoffzd1,9); // F8 - F11
	TArtTOF *tof8to11_2 = recopid -> DefineNewTOF("F8pl","F11pl-1",tofoffzd2,9); // F8 - F11
	TArtTOF *tof8to11_3 = recopid -> DefineNewTOF("F8pl","F11pl-1",tofoffzd3,9); // F8 - F11
	//TArtTOF *tof8to11_test = recopid -> DefineNewTOF("F8pl","F11pl-1",tofoffzd2,9); // F8 - F11

	// Reconstruction of IC observables for ID
	TArtBeam *beam_br = recopid -> DefineNewBeam(rips3to5,rips5to7,tof3to7,"F7IC");
	TArtBeam *beam_br2 = recopid -> DefineNewBeam(rips3to5,rips5to7,tof3to7,"F7IC");
	//TArtBeam *beam_zd_89 = recopid -> DefineNewBeam(rips8to9,tof8to11_1,"F11IC");
	TArtBeam *beam_zd_810 = recopid -> DefineNewBeam(rips8to10,tof8to11_1,"F11IC");
	TArtBeam *beam_zd_1011 = recopid -> DefineNewBeam(rips10to11,tof8to11_2,"F11IC");
	//test
	//TArtBeam *beam_zd_811_test  = recopid -> DefineNewBeam(rips8to9,rips10to11,tof8to11,"F11IC");
	TArtBeam *beam_zd_1011_test = recopid -> DefineNewBeam(rips8to10,rips10to11,tof8to11_3,"F11IC");
	//TArtBeam *beam_zd_810_test  = recopid -> DefineNewBeam(rips8to9,rips8to10,tof8to11,"F11IC");
	/*		
	TClonesArray *rips_array = recopid -> GetRIPSArray();
	TClonesArray *tof_array  = recopid -> GetTOFArray();
	TClonesArray *beam_array = recopid -> GetBeamArray();
	*/
	//Define Out put tree	
	//TFile *fout = new TFile(Form("data/rootfiles/test/run%04d_%3.3f_%3.3f.root",run,tofoffbr,tofoffzd),"RECREATE");
	TFile *fout;
	if( test )fout = new TFile(Form("data/rootfiles/test/run%04d_%3.3f_%3.3f_x.root",run,tofoffbr,tofoffzd1),"RECREATE");
	else fout = new TFile(Form("data/rootfiles/new/run%04d.root",run),"RECREATE");
	//	TFile *fout = new TFile(Form("rootfiles/calibration/calibration%04d.root",run),"RECREATE");
	TTree *tree = new TTree("tree","tree");

	// define data nodes which are supposed to be dumped to tree 
	//EventInfo is importand for the fBit information to know the trigger!
//	TClonesArray * info_array = (TClonesArray *)sman->FindDataContainer("EventInfo");
//	tree->Branch(info_array->GetName(),&info_array);

	TClonesArray * ppac_array = 
		(TClonesArray *)sman->FindDataContainer("BigRIPSPPAC");
	tree->Branch(ppac_array->GetName(),&ppac_array);

	TClonesArray * pla_array = 
		(TClonesArray *)sman->FindDataContainer("BigRIPSPlastic");
	tree->Branch(pla_array->GetName(),&pla_array);

	TClonesArray * ic_array = 
		(TClonesArray *)sman->FindDataContainer("BigRIPSIC");
	tree->Branch(ic_array->GetName(),&ic_array);

	TClonesArray * fpl_array = 
		(TClonesArray *)sman->FindDataContainer("BigRIPSFocalPlane");
	tree->Branch(fpl_array->GetName(),&fpl_array);

	//Dali data
	TClonesArray * dali_array=
		(TClonesArray *)sman->FindDataContainer("DALINaI");
	tree->Branch(dali_array->GetName(),&dali_array);
	//PID data
	TClonesArray *rips_array = 
		(TClonesArray *)sman->FindDataContainer("BigRIPSRIPS");
	tree->Branch(rips_array->GetName(),&rips_array); 

	TClonesArray *tof_array  = 
		(TClonesArray *)sman->FindDataContainer("BigRIPSTOF");
	tree->Branch(tof_array->GetName(),&tof_array);   

	TClonesArray *beam_array = 
		(TClonesArray *)sman->FindDataContainer("BigRIPSBeam");	
	tree->Branch(beam_array->GetName(),&beam_array); 

	//_________________________________________________________________________
	//Making new branches
	// Coincidence Registerg
	Int_t neve = 0;
	Int_t triggerbit = 0;

	tree->Branch("triggerbit", &triggerbit,"triggerbit/I");
	tree->Branch("neve", &neve,"neve/I");

	//%%%%%%%%%%%%%%%%%%%%%%
	//F8 target position && PPACs
	Double_t tgtx = -9999;
	Double_t tgty = -9999;
	Double_t tgta = -9999;
	Double_t tgtb = -9999;

	tree->Branch("tgtx",&tgtx,"tgtx/D");
	tree->Branch("tgty",&tgty,"tgty/D");
	tree->Branch("tgta",&tgta,"tgta/D");
	tree->Branch("tgtb",&tgtb,"tgtb/D");

	Double_t FX[6] = {0x0},FY[6] = {0x0}, FA[6] = {0x0}, FB[6] = {0x0}, F8T, F7T;
	Double_t delta[6] = {0x0};

	tree->Branch("FX",FX,"FX[6]/D");
	tree->Branch("FY",FY,"FY[6]/D");
	tree->Branch("FA",FA,"FA[6]/D");
	tree->Branch("FB",FB,"FB[6]/D");
  
	tree->Branch("F8T",&F8T,"F8T/D");
	tree->Branch("F7T",&F8T,"F7T/D");

	tree->Branch("delta",delta,"delta[6]/D");
	Double_t F10B = -9999;
	Double_t F11B = -9999;
	tree->Branch("F10B",&F10B,"F10B/D"); 
	tree->Branch("F11B",&F11B,"F11B/D");

	//%%%%%%%%%%%%%%%%%%%%%%%%
	// Plastic
	Double_t F3PLA_QL;  tree->Branch("F3PLA_QL",&F3PLA_QL,"F3PLA_QL/D");
	Double_t F3PLA_QR;  tree->Branch("F3PLA_QR",&F3PLA_QR,"F3PLA_QR/D");

	Double_t F7PLA_QL;  tree->Branch("F7PLA_QL",&F7PLA_QL,"F7PLA_QL/D");
	Double_t F7PLA_QR;  tree->Branch("F7PLA_QR",&F7PLA_QR,"F7PLA_QR/D");

	Double_t F8PLA_QL;  tree->Branch("F8PLA_QL",&F8PLA_QL,"F8PLA_QL/D");
	Double_t F8PLA_QR;  tree->Branch("F8PLA_QR",&F8PLA_QR,"F8PLA_QR/D");

	Double_t F11PLA_QL;  tree->Branch("F11PLA_QL",&F11PLA_QL,"F11PLA_QL/D");
	Double_t F11PLA_QR;  tree->Branch("F11PLA_QR",&F11PLA_QR,"F11PLA_QR/D");

	//%%%%%%%%%%%%%%%%%%%%%%%%
	//PID  Plot
	Double_t DELTA[3]; 
	Double_t ANGLE[3]; 
	Double_t BRHO[3];  
	Double_t TOF[3];   
	Double_t BETA[5];  
	Double_t AOQ[5];   
	Double_t ZET[4];   
	Double_t ESq=-9999;
	//Double_t delta[2];
	Double_t DIST = 9999;
	Int_t F3Pla_TR[fpl_ch] = {0x0};
	Int_t F3Pla_TL[fpl_ch] = {0x0};
	Int_t F7Pla_TR[fpl_ch] = {0x0};
	Int_t F7Pla_TL[fpl_ch] = {0x0};
	Int_t F8Pla_TR[fpl_ch] = {0x0};
	Int_t F8Pla_TL[fpl_ch] = {0x0};
	//Int_t F10Pla_TR[fpl_ch] = {0x0};
	//Int_t F10Pla_TL[fpl_ch] = {0x0};
	Int_t F11Pla_TR[fpl_ch] = {0x0};
	Int_t F11Pla_TL[fpl_ch] = {0x0};

	Int_t F3Pla_QR[fpl_ch] = {0x0};
	Int_t F3Pla_QL[fpl_ch] = {0x0};
	Int_t F7Pla_QR[fpl_ch] = {0x0};
	Int_t F7Pla_QL[fpl_ch] = {0x0};
	Int_t F8Pla_QR[fpl_ch] = {0x0};
	Int_t F8Pla_QL[fpl_ch] = {0x0};
	Int_t F11Pla_QR[fpl_ch] = {0x0};
	Int_t F11Pla_QL[fpl_ch] = {0x0};
	//tree->Branch("delta",delta,"delta[2]/D");
	tree->Branch("Distanceto",&DIST,"DIST/D");

	tree->Branch("DELTA",DELTA,"DELTA[3]/D");
	tree->Branch("ANGLE",ANGLE,"ANGLE[3]/D");
	tree->Branch("BRHO",BRHO,"BRHO[3]/D");
	tree->Branch("TOF",TOF,"TOF[3]/D");
	tree->Branch("BETA",BETA,"BETA[5]/D");
	tree->Branch("AOQ",AOQ,"AOQ[5]/D");
	tree->Branch("ZET",ZET,"ZET[4]/D");
	tree->Branch("ESq",&ESq,"ESq/D");
	tree->Branch("F3Pla_TR",F3Pla_TR,Form("F3Pla_TR[%i]/I",fpl_ch));
	tree->Branch("F3Pla_TL",F3Pla_TL,Form("F3Pla_TL[%i]/I",fpl_ch));
	tree->Branch("F7Pla_TR",F7Pla_TR,Form("F7Pla_TR[%i]/I",fpl_ch));
	tree->Branch("F7Pla_TL",F7Pla_TL,Form("F7Pla_TL[%i]/I",fpl_ch));
	tree->Branch("F8Pla_TR",F8Pla_TR,Form("F8Pla_TR[%i]/I",fpl_ch));
	tree->Branch("F8Pla_TL",F8Pla_TL,Form("F8Pla_TL[%i]/I",fpl_ch));
	//	tree->Branch("F10Pla_TR",F10Pla_TR,Form("F10Pla_TR[%i]/I",fpl_ch));
	//tree->Branch("F10Pla_TL",F10Pla_TL,Form("F10Pla_TL[%i]/I",fpl_ch));
	tree->Branch("F11Pla_TR",F11Pla_TR,Form("F11Pla_TR[%i]/I",fpl_ch));
	tree->Branch("F11Pla_TL",F11Pla_TL,Form("F11Pla_TL[%i]/I",fpl_ch));

	tree->Branch("F3Pla_QR",F3Pla_QR,Form("F3Pla_QR[%i]/I",fpl_ch));
	tree->Branch("F3Pla_QL",F3Pla_QL,Form("F3Pla_QL[%i]/I",fpl_ch));
	tree->Branch("F7Pla_QR",F7Pla_QR,Form("F7Pla_QR[%i]/I",fpl_ch));
	tree->Branch("F7Pla_QL",F7Pla_QL,Form("F7Pla_QL[%i]/I",fpl_ch));
	tree->Branch("F8Pla_QR",F8Pla_QR,Form("F8Pla_QR[%i]/I",fpl_ch));
	tree->Branch("F8Pla_QL",F8Pla_QL,Form("F8Pla_QL[%i]/I",fpl_ch));
	tree->Branch("F11Pla_QR",F11Pla_QR,Form("F11Pla_QR[%i]/I",fpl_ch));
	tree->Branch("F11Pla_QL",F11Pla_QL,Form("F11Pla_QL[%i]/I",fpl_ch));

	Int_t  fplMulti[fpl_ch] = {0x0};  tree->Branch("fplMulti",fplMulti,Form("fplMulti[%i]/I",fpl_ch));
	Int_t gate[4] = {0x0}; tree->Branch("gate",gate,"gate[4]/I");
	//%%%%%%%%%%%%%%%%%%%%%%%%
	//TKE & Analog sum
	Double_t TKE = -9999;   
	Double_t ASUM = -9999;   

	tree->Branch("TKE",&TKE,"TKE/D");
	tree->Branch("ASUM",&ASUM,"ASUM/D");

	//%%%%%%%%%%%%%%%%%%%%%%
	//DALI

	Int_t dalimult = 0;
	Int_t dalitimetruemult = 0;


	
	Int_t dalimultwotime = 0;
	Int_t dalimultthres = 0;
	Int_t dalitimetruemultthres = 0;
	tree->Branch("dalimultwotime",&dalimultwotime,"dalimultwotime/I");
	tree->Branch("dalimultthres",&dalimultthres,"dalimultthres/I");
	tree->Branch("dalitimetruemultthres",&dalitimetruemultthres,"dalitimetruemultthres/I");
	
	tree->Branch("dalimult",&dalimult,"dalimult/I");
	tree->Branch("dalitimetruemult",&dalitimetruemult,"dalitimetruemult/I");



	//%%%%%%%%%%%%%%%%%%%%%%
	//LaBr QTC Array
	Double_t LOWGAIN[8]; 
	Double_t LOWGAINC[8]; 
	Double_t MIDGAIN[8]; 
	Double_t MIDGAINC[8]; 
	Double_t MIDGAINCT[8]; 
	Double_t HIGAIN[8];  
	Double_t HIGAINC[8];  
	Double_t QTCTime[8];  
	Double_t QTCTimeC[8];
	Double_t MIDGAINCTA;

	MIDGAINCTA = 0.;

	tree->Branch("LOWGAIN",LOWGAIN,"LOWGAIN[8]/D");
	tree->Branch("LOWGAINC",LOWGAINC,"LOWGAINC[8]/D");
	tree->Branch("MIDGAIN",MIDGAIN,"MIDGAIN[8]/D");
	tree->Branch("MIDGAINC",MIDGAINC,"MIDGAINC[8]/D");
	tree->Branch("MIDGAINCT",MIDGAINCT,"MIDGAINCT[8]/D");
	tree->Branch("HIGAIN",HIGAIN,"HIGAIN[8]/D");
	tree->Branch("HIGAINC",HIGAINC,"HIGAINC[8]/D");
	tree->Branch("QTCTime",QTCTime,"QTCtime[8]/D");
	tree->Branch("QTCTimeC",QTCTimeC,"QTCtimeC[8]/D");
	tree->Branch("MIDGAINCTA", &MIDGAINCTA, "MIDGAINCTA/D");
	Float_t end = 1e15;
	if(test) end = 0.1e6;

	//while(estore->GetNextEvent()){
	while(estore->GetNextEvent() && neve < end){
	  if(neve%50000==0 && neve != 0) std::cout << "event: " << neve << std::endl;

	  //**** Before doing something first get the time of the plastics
	  //Time of the hits in Plastic at F3,F7,F8,F11
	  Int_t PlaTDCMAP[8] = {16,17,18,19,20,21,22,23};
	  Int_t PlaQMAP[8] = {0,1,2,3,4,5,6,7};
	  Int_t PlaTDCRAW[8][100]={0x0};
	  Int_t PlaQRAW[8][100];
	  Int_t tref63 = -1000;
	  for(short ii = 0; ii < fpl_ch; ii++)fplMulti[ii] = 0;
	  for(int i = 0; i < rawevent -> GetNumSeg(); i++){
	    TArtRawSegmentObject *seg = rawevent -> GetSegment(i);
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
		      //if(val>0){
		      fplMulti[k] = fplMulti[k]+1;
		      PlaTDCRAW[k][j] = val;
		      // }else  PlaTDCRAW[k][j] = -1000;
		    }
		  }
		}/*else if(geo==3){
		  Int_t ch   = d -> GetCh();
		  Int_t edge = d -> GetEdge();
		  Int_t val  = (Int_t)d -> GetVal();
		  for(Int_t k=0; k<fpl_ch; k++)
		    if(edge==0 && ch == PlaQMAP[k])
		      PlaQRAW[k][j] = val;
		}*/
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
	  //int *index8; index8 = new int[100]; TMath::Sort(100, PlaTDCRAW[8], index8, true);
	  //int *index9; index9 = new int[100]; TMath::Sort(100, PlaTDCRAW[9], index9, true);
	  Double_t fpltcal = 0.024413;

	  for(short ii=0; ii<fpl_ch; ii++){
	    F3Pla_TR[ii] =0;//(PlaTDCRAW[0][index0[ii]] - tref63)*fpltcal;
	    F3Pla_TL[ii] =0;//(PlaTDCRAW[1][index1[ii]] - tref63)*fpltcal;
	    F7Pla_TR[ii] =0;//(PlaTDCRAW[2][index2[ii]] - tref63)*fpltcal;
	    F7Pla_TL[ii] =0;//(PlaTDCRAW[3][index3[ii]] - tref63)*fpltcal;
	    F8Pla_TR[ii] =0;//(PlaTDCRAW[4][index4[ii]] - tref63)*fpltcal;
	    F8Pla_TL[ii] =0;//(PlaTDCRAW[5][index5[ii]] - tref63)*fpltcal;
	    F11Pla_TR[ii]=0;//(PlaTDCRAW[6][index6[ii]] - tref63)*fpltcal;
	    F11Pla_TL[ii]=0;//(PlaTDCRAW[7][index7[ii]] - tref63)*fpltcal;
	    //	    F11Pla_TR[ii]=(PlaTDCRAW[8][index8[ii]] - tref63)*fpltcal;
	    //F11Pla_TL[ii]=(PlaTDCRAW[9][index9[ii]] - tref63)*fpltcal;
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
	    //	    F11Pla_TR[ii]=(PlaTDCRAW[8][index8[ii]] - tref63)*fpltcal;
	    //F11Pla_TL[ii]=(PlaTDCRAW[9][index9[ii]] - tref63)*fpltcal;
	  }
	  delete [] index0;
	  delete [] index1;
	  delete [] index2;
	  delete [] index3;
	  delete [] index4;
	  delete [] index5;
	  delete [] index6;
	  delete [] index7;
	  //delete [] index8;
	  //delete [] index9;
	  //**** Set

	  //Set the plastic time
	  Float_t F3Pla_Time = -9999, F7Pla_Time = -9999, F8Pla_Time = -9999, F11Pla_Time = -9999;
	  for(short ii = 0; ii<4;ii++)gate[ii]=0;
	  TArtPlastic *pla;

	  for(short ii=0;ii<fpl_ch; ii++){
	    if( gated(F3Pla_TL[ii], 3) ){
	      if(gated(F3Pla_TR[ii], 3) ){
		F3Pla_Time = (F3Pla_TR[ii]+F3Pla_TL[ii]);
		gate[0]++;
	      }else{
		F3Pla_Time = F3Pla_TL[ii]/2;
		gate[0]++;
	      }
	      pla = cpla->FindPlastic((char*)"F3pl");
	      pla ->SetTime(F3Pla_Time);
	    }
	    if( gated(F7Pla_TL[ii],7) ){
	      if(gated(F7Pla_TR[ii],7) ){
		F7Pla_Time = (F7Pla_TR[ii]+F7Pla_TL[ii]);
		gate[1]++;
	      }else{
		F7Pla_Time = F7Pla_TL[ii];
		gate[1]++;
	      }
	      pla = cpla->FindPlastic((char*)"F7pl");
	      pla ->SetTime(F7Pla_Time);
	    }
	    if( gated(F8Pla_TL[ii],8) ){
	      if(gated(F8Pla_TR[ii],8) ){
		F8Pla_Time = (F8Pla_TR[ii]+F8Pla_TL[ii]);
		gate[2]++;
	      }else{
		F8Pla_Time = F8Pla_TL[ii];
		gate[2]++;
	      }
	      pla = cpla->FindPlastic((char*)"F7pl");
	      pla ->SetTime(F7Pla_Time);
	    }
	    if( gated(F11Pla_TL[ii],11) ){
	      if(gated(F11Pla_TR[ii],11) ){
		F11Pla_Time = (F11Pla_TR[ii]+F11Pla_TL[ii]);
		gate[3]++;		
	      }else{
		F11Pla_Time = F11Pla_TL[ii];
		gate[3]++;
	      }
	      pla = cpla->FindPlastic((char*)"F11pl-1");
	      pla ->SetTime(F11Pla_Time);
	    }
	  }
	  //if(fplMulti[0]>1)continue;

	  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	  //Plastic ADC
	  F3PLA_QL=-9999; F3PLA_QR=-9999; 
	  
	  F7PLA_QL=-9999; F7PLA_QR=-9999;
	  
	  F8PLA_QL=-9999; F8PLA_QR=-9999;
	  
	  F11PLA_QL=-9999; F11PLA_QR=-9999; 
	  
	  //	  TArtPlastic *pla;
	  //F3 Plastic
	  pla = cpla->FindPlastic((char*)"F3pl");
	  if( pla ){
	    F3PLA_QL = pla->GetQLRaw(); F3PLA_QR = pla->GetQRRaw();
	    pla->SetTime(F3Pla_Time);
	  }
	  
	  // F7 Plastic
	  pla = cpla->FindPlastic((char*)"F7pl");
	  if( pla ){
	    F7PLA_QL = pla->GetQLRaw(); F7PLA_QR = pla->GetQRRaw();
	    pla->SetTime(F7Pla_Time);
	  }
	  
	  // F8 Plastic
	  pla = cpla->FindPlastic((char*)"F8pl");
	  if( pla ){
	    F8PLA_QL = pla->GetQLRaw(); F8PLA_QR = pla->GetQRRaw();
	    pla->SetTime(F8Pla_Time);
	  }
	  
	  // F11 Plastic
	  pla = cpla->FindPlastic((char*)"F11pl-1");
	  if( pla ){
	    F11PLA_QL = pla->GetQLRaw(); F11PLA_QR = pla->GetQRRaw();
	    pla->SetTime(F11Pla_Time);
	  }
	  

	  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	  // Coincidence Register
	  triggerbit = -9999;
	  //coin->ClearData();
	  //coin->LoadData();
	  //Int_t ftriggerbit = triggerinfo -> GetTriggerBit();
	  //cout <<  ftriggerbit << endl;
	  Int_t ftriggerbit = -9999;
	  for(int i=0;i<rawevent -> GetNumSeg();i++){
	    TArtRawSegmentObject *seg = rawevent -> GetSegment(i);
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
	      triggerbit = i+1;
	  }
	  
	  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	  //Making the BigRIPS tree calibration
	  cpid -> ClearData();
	  cpid -> ReconstructData();
	  
	  TArtFocalPlane* tfpl;
	  //F3X = -999; F5X = -999; F7X = -999; F9X = -999; F3A = -999;F5A = -999;F7A = -999;F9A = -999; 
	  Int_t plane_[6] = {3,5,7,8,10,11};
	  for(short ii=0;ii<6;ii++){
	    tfpl = cfpl->FindFocalPlane(plane_[ii]);
	    //	    FX[ii] = -999.; FY[ii] = -999.; FA[ii] = -999.; FB = -999.;
	    if(tfpl){TVectorD* vec=tfpl->GetOptVector(); 
	      FX[ii]=(*vec)(0); 
	      FA[ii]=(*vec)(1);
	      FY[ii]=(*vec)(2);
	      FB[ii]=(*vec)(3);
	      delta[ii]=(*vec)(4);
	    }else cout<< "not found"<<endl;
	  }
	    /*	  } 
	  if(tfpl){TVectorD* vec=tfpl->GetOptVector(); F3X=(*vec)(0); F3A=(*vec)(1);}else cout<< "not found"<<endl;
	  tfpl = cfpl->FindFocalPlane(5); 
	  if(tfpl){TVectorD* vec=tfpl->GetOptVector(); F5X=(*vec)(0); F5A=(*vec)(1);}
	  tfpl = cfpl->FindFocalPlane(7); 
	  if(tfpl){TVectorD* vec=tfpl->GetOptVector(); F7X=(*vec)(0); F7A=(*vec)(1);}
	  tfpl = cfpl->FindFocalPlane(9); 
	  if(tfpl){TVectorD* vec=tfpl->GetOptVector(); F9X=(*vec)(0); F9A=(*vec)(1);}
	    */
	  F8T = cpla -> FindPlastic("F8pl") -> GetTime();
	  F7T = F7Pla_Time;//cpla -> FindPlastic("F7pl") -> GetTime();

	  // for skip F3 ppac tracking --->
	  if(run>52){
	    TArtFocalPlane *fpl = cfpl->FindFocalPlane(3);
	    if( fpl ){
	      TVectorD *vec=fpl->GetOptVector();
	      Double_t vec_tmp[4]={0,0,0,0};
	      vec->SetElements(vec_tmp);
	    }
	  }
	  // <-----
	  
	  
	  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	  //Reconstructiong the PID
	  recopid->ClearData();
	  recopid->ReconstructData();
	  
	  GetIC(rips8to10,tof8to11_1,"F11IC", beam_zd_810,kFALSE);
	  TArtFocalPlane *fpl11; TVectorD *vec11;
	  Float_t F11X, F11A, F11Y, F11B;
	  fpl11 = cfpl -> FindFocalPlane(11);
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
	
	  //fIC->ReconstructData();
	  /*
	    Float_t fSqSum = 1;
	    Int_t fired = 0;
	    for(short ich=0; ich<12;ich++){
	    Float_t adc = fIC->GetRawADC(ich) - fIC->GetPedestal(ich);
	    if(adc>0){
	    if(ich==0)
	    adc = 1000;
	    fSqSum = fSqSum * adc;
	    fired++;
	    }
	    }
	    if(fired > 0)fSqSum = TMath::Power(fSqSum,1./fired);
	    fIC->SetEnergySqSum(para->GetCh2MeV(0) + para->GetCh2MeV(1)*fSqSum);
	  */  
	  
	  recopid->ClearData();
	  recopid->ReconstructData();
  	  
	  for( int i=0; i<3; i++ ){
	    DELTA[i] = -9999;
	    ANGLE[i] = -9999;
	    BRHO[i] = -9999;
	    TOF[i] = -9999;
	    BETA[i] = -9999;
	  }
	  for( int i=0; i<4; i++ ){
	    ZET[i] = -9999;
	    AOQ[i] = -9999;
	  }
	    AOQ[5] = -9999;

	  for(int i=0; i<rips_array -> GetEntries(); i++){
	    TArtRIPS *myrips = (TArtRIPS*)rips_array -> At(i);
	    //std::cout << "rips " << i << " : " << myrips->GetDelta() << std::endl;
	    DELTA[i] = myrips->GetDelta();
	    ANGLE[i] = myrips->GetAngle();
	    BRHO[i] = myrips->GetBrho();
	  }
	  for(int i=0; i<tof_array -> GetEntries(); i++){
	    TArtTOF *mytof = (TArtTOF*)tof_array -> At(i);
	    //std::cout << "tof " << i << " : " << mytof->GetTOF() << std::endl;
	    TOF[i] = mytof -> GetTOF();
	  }
	  
	  // delta[0] = beam_zd_810 ->GetDelta();
	  // delta[1] = beam_zd_1011 ->GetDelta();
	  DIST = beam_br->IsOnHeap();
	  BETA[0] = beam_br -> GetBeta();
	  BETA[1] = beam_zd_810 -> GetBeta();
	  BETA[2] = beam_zd_1011 -> GetBeta();
	  BETA[3] = beam_zd_1011_test -> GetBeta();

	  AOQ[0] = beam_br -> GetAoQ();
	  AOQ[1] = beam_zd_810 -> GetAoQ();
	  
	  AOQ[2] = beam_zd_1011 -> GetAoQ();//change so that the matrixelements are fitting see pid.C for values and add dependency on TOF then fix tof-offset and check everything again and ... then fix ZET same for Sn132 and 128 this is already done with the matrix Sn132 for the Sn132
	  AOQ[3] = beam_zd_1011_test -> GetAoQ();
	  //AOQ[4] = beam_zd_89->GetAoQ();// not a good quantity at all

	  if(run<100)ZET[0] = beam_br -> GetZet();
	  else {

	    GetIC(rips3to5,tof3to7,"F7IC", beam_br2,kTRUE);
	    Double_t de_v2 = TMath::Log(fIC2->GetIonPair()*BETA[0]*BETA[0])-TMath::Log(1-BETA[0]*BETA[0])-BETA[0]*BETA[0];
	    ZET[0] = fIC2->GetZCoef(0)*TMath::Sqrt((fIC2->GetEnergySqSum()+E_offset_br)/de_v2)*BETA[0]+fIC2->GetZCoef(1);
	  }
	  ZET[1] = beam_zd_810 -> GetZet();
	  //ZET[2] = ZET[1] - (BETA[1]*BETA[1]*1291.43-BETA[1]*878.207+170.106)+38;
	  
	  //Int_t ionpair = 4866;
	  //GetIC(rips8to10,tof8to11_1,"F11IC", beam_zd_810,kFALSE);
	  Double_t de_v = TMath::Log(fIC->GetIonPair()*BETA[2]*BETA[2])-TMath::Log(1-BETA[2]*BETA[2])-BETA[2]*BETA[2];
	  //E_offset = ESq_offset(run);
	  ZET[2] = fIC->GetZCoef(0)*TMath::Sqrt((fIC->GetEnergySqSum()+E_offset)/de_v)*BETA[2]+fIC->GetZCoef(1);
	  ESq = fIC->GetEnergySqSum();
	  
	  de_v = TMath::Log(fIC->GetIonPair()*BETA[3]*BETA[3])-TMath::Log(1-BETA[3]*BETA[3])-BETA[3]*BETA[3];
	  ZET[3] = fIC->GetZCoef(0)*TMath::Sqrt((fIC->GetEnergySqSum()+E_offset)/de_v)*BETA[3]+fIC->GetZCoef(1);

	  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%		


	  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	  //BigRIPS F8 target profile 
	  tgtx = -9999;
	  tgty = -9999;
	  tgta = -9999; 
	  tgtb = -9999;
	  
	  TArtFocalPlane *fpl8; TVectorD *vec8;
	  fpl8 = cfpl -> FindFocalPlane(8);
	  if( fpl8 ){vec8 = fpl8 -> GetOptVector();
	    tgtx = (*vec8)(0); tgta = (*vec8)(1); tgty = (*vec8)(2); tgtb = (*vec8)(3);
	  }
	  /*
	  TArtFocalPlane *fpl10; TVectorD *vec10;
	  fpl10 = cfpl -> FindFocalPlane(10);
	  if( fpl10 ){vec10 = fpl10 -> GetOptVector();
	    F10X = (*vec10)(0); F10A = (*vec10)(1); F10Y = (*vec10)(2); F10B = (*vec10)(3);
	  }
	  */
	  //ZET[3] = fIC->GetZCoef(0)*TMath::Sqrt((fIC->GetEnergySqSum()+E_offset)/de_v)*BETA[3]+fIC->GetZCoef(1);

	  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	  // TKE & ASUM 
	  TKE = -9999; 
	  ASUM = -9999;
	  
	  for(int i=0;i<rawevent -> GetNumSeg();i++){
	    TArtRawSegmentObject *seg = rawevent -> GetSegment(i);
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
	  
	  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	  // LaBr QTC Array
	  for( int i=0; i<8; i++ ){
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
	  for(int i=0; i<rawevent -> GetNumSeg(); i++){
	    TArtRawSegmentObject *seg = rawevent -> GetSegment(i);
	    Int_t fpl = seg -> GetFP();
	    if(fpl==8){
	      for(int j=0; j < seg -> GetNumData(); j++){
		TArtRawDataObject *d = seg -> GetData(j);
		Int_t geo  = d -> GetGeo();
		Int_t ch   = d -> GetCh();
		Int_t val  = d -> GetVal();
		Int_t edge = d -> GetEdge();
		if(geo==9) qtc_dummy[edge][ch] = val + (rndm.Rndm()- 0.5);
		if(geo ==9 && ch==127 && edge == 0) tref = val;					
	      }
	    }
	  }
	  //cout << "ch " << ch << "edge " << edge << "val " << endl; 
	  //cout << "1 " << qtc_dummy[1][1] << " " << qtc_dummy[0][1] << endl;
	  //cout << "4 " << qtc_dummy[1][4] << " " << qtc_dummy[0][4] << endl;
	  Double_t beta = BETA[0];
	  Double_t theta = 30.0 * TMath::Pi() / 180.;
	  //h- gain
	  HIGAIN[0] = qtc_dummy[1][28]-qtc_dummy[0][28];
	  HIGAIN[0] = HIGAIN[0] * 1.5194 -2221.76;
	  
	  HIGAIN[1]=qtc_dummy[1][29]-qtc_dummy[0][29];
	  HIGAIN[1] = HIGAIN[1] * 1.3509 - 1986.75;
	  
	  HIGAIN[2]=qtc_dummy[1][30]-qtc_dummy[0][30];
	  HIGAIN[2] = HIGAIN[2] * 1.54495 - 2289.73;
	  
	  HIGAIN[3]=qtc_dummy[1][1] -qtc_dummy[0][1];
	  HIGAIN[3] = HIGAIN[3] * 1.3224 - 2053.94;
	  
	  HIGAIN[4]=qtc_dummy[1][31]-qtc_dummy[0][31];
	  HIGAIN[4] = HIGAIN[4] * 1.51863 - 2250.28;
	  
	  HIGAIN[5]=qtc_dummy[1][4] -qtc_dummy[0][4];
	  HIGAIN[5] = HIGAIN[5] * 1.54124 - 2214.31;
	  
	  HIGAIN[6]=qtc_dummy[1][3] -qtc_dummy[0][3];
	  HIGAIN[6] = HIGAIN[6] * 1.31455 - 1911.24;
	  
	  HIGAIN[7]=qtc_dummy[1][2] -qtc_dummy[0][2];
	  HIGAIN[7] = HIGAIN[7] * 1.43441 - 2082.49;
	  
	  for(Int_t i=0; i < 8; i++){
	    HIGAINC[i] = HIGAIN[i] * (1 - beta * cos(theta))/sqrt(1-beta * beta);
	  };
	  
	  //for(int i=0;i<8;i++){cout<<QTC_h[i];cout<<endl;}
	  //m- gain
	  MIDGAIN[0]=qtc_dummy[1][24]-qtc_dummy[0][24];
	  MIDGAIN[0] = MIDGAIN[0] * 11.2806 - 16752.7;
	  
	  MIDGAIN[1]=qtc_dummy[1][25]-qtc_dummy[0][25];
	  MIDGAIN[1] = MIDGAIN[1] * 11.1504 - 16238.2;
	  
	  MIDGAIN[2]=qtc_dummy[1][26]-qtc_dummy[0][26];
	  MIDGAIN[2] = MIDGAIN[2] * 12.3476 - 17523.1;
	  
	  MIDGAIN[3]=qtc_dummy[1][6] -qtc_dummy[0][6];
	  MIDGAIN[3] = MIDGAIN[3] * 10.983 - 16555;
	  
	  MIDGAIN[4]=qtc_dummy[1][27]-qtc_dummy[0][27];
	  MIDGAIN[4] = MIDGAIN[4] * 12.1714 - 17515.1;
	  
	  MIDGAIN[5]=qtc_dummy[1][9] -qtc_dummy[0][9];
	  MIDGAIN[5] = MIDGAIN[5] * 11.2647 - 15875.3;
	  
	  MIDGAIN[6]=qtc_dummy[1][8] -qtc_dummy[0][8];
	  MIDGAIN[6] = MIDGAIN[6] * 10.1906 - 15031.1;
	  
	  MIDGAIN[7]=qtc_dummy[1][7] -qtc_dummy[0][7];
	  MIDGAIN[7] = MIDGAIN[7] * 10.8408 - 14837.9;
	  
	  for(Int_t i=0; i < 8; i++){
	    MIDGAINC[i] = MIDGAIN[i] * (1 - beta * cos(theta))/sqrt(1-beta * beta);
	  };
	  
	  //for(int i=0;i<8;i++){cout<<QTC_m[i];cout<<endl;}
	  //l- gain
	  LOWGAIN[0]=qtc_dummy[1][20]-qtc_dummy[0][20];
	  LOWGAIN[1]=qtc_dummy[1][21]-qtc_dummy[0][21];
	  LOWGAIN[2]=qtc_dummy[1][22]-qtc_dummy[0][22];
	  LOWGAIN[3]=qtc_dummy[1][11]-qtc_dummy[0][11];
	  LOWGAIN[4]=qtc_dummy[1][23]-qtc_dummy[0][23];
	  LOWGAIN[5]=qtc_dummy[1][14]-qtc_dummy[0][14];
	  LOWGAIN[6]=qtc_dummy[1][13]-qtc_dummy[0][13];
	  LOWGAIN[7]=qtc_dummy[1][12]-qtc_dummy[0][12];
	  
	  QTCTime[0] = qtc_dummy[0][24] - tref;
	  QTCTime[1] = qtc_dummy[0][25] - tref;
	  QTCTime[2] = qtc_dummy[0][26] - tref;
	  QTCTime[3] = qtc_dummy[0][6] - tref;
	  QTCTime[4] = qtc_dummy[0][27] - tref;
	  QTCTime[5] = qtc_dummy[0][9] - tref;
	  QTCTime[6] = qtc_dummy[0][8] - tref;
	  QTCTime[7] = qtc_dummy[0][7] - tref;
	  for(Int_t i=0;i<8;i++){
	    QTCTime[i] *= 0.09766;
	    QTCTimeC[i] = QTCTime[i] - QTCTimeOffset[i];
	  }
	  

	  for(Int_t i=0; i < 8; i++){
	    if(QTCTimeC[i] > -15. && QTCTimeC[i] < 15.){
	      MIDGAINCT[i] = MIDGAINC[i];
	      MIDGAINCTA += MIDGAIN[i];
	    }else{
	      MIDGAINCT[i] = -999.;
	    }
	  }

	  if(MIDGAINCTA > 0.){
	    MIDGAINCTA = MIDGAINCTA * (1 - beta * cos(theta))/sqrt(1-beta * beta);	  
	  }else{
	    MIDGAINCTA = -999.;
	  }

	  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	  //Making DALI
	  dalicalib -> ClearData();
	  dalicalib -> SetPlTime(F7Pla_Time);//cpla -> FindPlastic("F7pl") -> GetTime());
	  /*cout << "beta[2] " << BETA[0]<<endl;
	  dalicalib -> SetBeta(BETA[0]);
	  cout << "beta dali " << dalicalib->GetBeta() << endl;
	  */
	  //Add above to remove F8plastic tof.
	  dalicalib -> ReconstructData();
	  cout<<"ID "<< dpara->GetID() <<" QCal "<< dpara->GetQCal() <<" QPed "<<  dpara->GetQPed() <<endl
	  dalimultwotime = dalicalib -> GetMultWithoutT();
	  dalimultthres = dalicalib -> GetMultThres();
	  dalitimetruemultthres = dalicalib -> GetTimeTrueMultThres();
	  
	  dalimult = dalicalib -> GetMult();
	  dalitimetruemult = dalicalib -> GetTimeTrueMult();
	  
	  tree -> Fill();
	  neve++;
	}
	//cout<<"Writing the tree."<<endl;

	fout -> Write();
	fout -> Close();
	
	delete sman;
	sman = 0;
	para -> Delete();
	para = 0;
	dpara -> Delete();
	dpara = 0;
	delete recopid;
	recopid = 0;
}
