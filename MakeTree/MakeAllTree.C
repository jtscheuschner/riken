#include <iostream>
#include <string>
#include <time.h>

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

#include "TVector3.h"
#include "TMath.h"

#include "signal.h"

// BigRIPS parameters
#define tofoffbr  299.8// TOF3-7 offset VME
#define tofoffzd  -158.1 // TOF7-11 offset VME
//not needed knows already 
#define brho89 5.6435
/*
#define brho35 6.1667 
#define brho57 5.9024 

#define brhoTE 5.6435 //run185
/*
//128Sn Filled RUN
#define tofoffbr 301.// TOF3-7 offset VME
#define tofoffzd -156.5 // TOF8-11 offset VME

#define brho35 5.9930 
#define brho57 5.7312 
#define brho89 4.9686 */
#define brhoTE 4.9655 //run381

using namespace std;

// function to exit loop at keyboard interrupt. 
bool stoploop = false;
void stop_interrupt(){
  printf("keyboard interrupt\n");
  stoploop = true;
}

void MakeAllTree(Int_t run){

	gStyle->SetOptStat(111111);
	gSystem->Load("libXMLParser.so");
	gSystem->Load("libanacore.so");
	gSystem->Load("libanabrips.so");
	gROOT->cd();
	gSystem->AddIncludePath("-I$TARTSYS/include/");

	// Create StoreManager both for calibration "TArtCalib..." and treatment "TArtReco..."
	TArtStoreManager * sman = TArtStoreManager::Instance();

	// Create EventStore to control the loop and get the EventInfo
	TArtEventStore *estore = new TArtEventStore();
	estore->SetInterrupt(&stoploop); 
	estore->Open(Form("data/ridf/Sn/Sn132_LHe/labr14Sn%04d.ridf",run));
	//	estore->Open(Form("ridf/ggdaq04/calibration%04d.ridf",run));
	TArtRawEventObject *rawevent = (TArtRawEventObject *)sman->FindDataContainer("RawEvent");

	//Coincidence register LUPO
	TArtCalibCoin *coin = new TArtCalibCoin();
	TClonesArray * triggerbit_array = (TClonesArray *)estore -> GetEventInfoArray();
	TArtEventInfo *triggerinfo = (TArtEventInfo *)triggerbit_array -> At(0);
	
	// Create BigRIPSParameters to get Plastics, PPACs, ICs and FocalPlanes parameters from ".xml" files
	TArtBigRIPSParameters *para = TArtBigRIPSParameters::Instance();
	para -> LoadParameter("db/BigRIPSPPAC.xml");
	para -> LoadParameter("db/Plastic/BigRIPSPlastic.xml");
	para -> LoadParameter("db/BigRIPSIC.xml");
	para -> LoadParameter("db/FocalPlane.xml");

	//CREATE DALI PARAMETERS
	TArtDALIParameters *dpara = TArtDALIParameters::Instance();
	dpara -> LoadParameter("db/DALI.xml");

	TArtCalibDALI *dalicalib = new TArtCalibDALI();

	// Create CalibPID to get and calibrate raw data 
	TArtCalibPID *cpid = new TArtCalibPID();
	TArtCalibPPAC *cppac = cpid -> GetCalibPPAC();
	TArtCalibPlastic *cplastic = cpid -> GetCalibPlastic();
	TArtCalibFocalPlane *cfpl = cpid -> GetCalibFocalPlane();


	// Create RecoPID to get calibrated data and to reconstruct TOF, AoQ, Z 
	TArtRecoPID *recopid = new TArtRecoPID();

	// Definition of observables we want to reconstruct
	TArtRIPS *rips3to5 = recopid -> DefineNewRIPS(3,5,"matrix/mat1.mat","D3"); // F3 - F5
	TArtRIPS *rips5to7 = recopid -> DefineNewRIPS(5,7,"matrix/mat2.mat","D5"); // F5 - F7
	TArtRIPS *rips8to10 = recopid -> DefineNewRIPS(8,10,"matrix/F8F10_LargeAccAchr.mat","D7"); // F9 - F11  
	TArtRIPS *rips8to9 = recopid -> DefineNewRIPS(8,9,"matrix/F8F9_LargeAccAchr.mat",brho89); // F8 - F10  
	TArtRIPS *rips10to11 = recopid -> DefineNewRIPS(10,11,"matrix/F10F11_LargeAccAchr.mat","D8"); // F9 - F11  

	// Reconstruction of TOF DefineNewTOF(fisrt plane, second plane, time offset)
	TArtTOF *tof3to7  = recopid -> DefineNewTOF("F3pl","F7pl",tofoffbr,5); // F3 - F7
	TArtTOF *tof8to11 = recopid -> DefineNewTOF("F8pl","F11pl-1",tofoffzd,9); // F8 - F11

	// Reconstruction of IC observables for ID
	TArtBeam *beam_br = recopid -> DefineNewBeam(rips3to5,rips5to7,tof3to7,"F7IC");
	TArtBeam *beam_zd_89 = recopid -> DefineNewBeam(rips8to9,tof8to11,"F11IC");
	TArtBeam *beam_zd_810 = recopid -> DefineNewBeam(rips8to10,tof8to11,"F11IC");
	TArtBeam *beam_zd_1011 = recopid -> DefineNewBeam(rips10to11,tof8to11,"F11IC");

	TClonesArray *rips_array = recopid -> GetRIPSArray();
	TClonesArray *tof_array  = recopid -> GetTOFArray();
	TClonesArray *beam_array = recopid -> GetBeamArray();

	//Define Out put tree	
	TFile *fout = new TFile(Form("data/rootfiles/test/run%04d.root",run),"RECREATE");
	//TFile *fout = new TFile(Form("rootfiles/Sn128_LHe/run%04d.root",run),"RECREATE");
	TTree *tree = new TTree("tree","tree");

	tree->Branch(rips_array->GetName(),&rips_array); 
 	tree->Branch(tof_array->GetName(),&tof_array);   
	tree->Branch(beam_array->GetName(),&beam_array);
	// define data nodes which are supposed to be dumped to tree 
	//EventInfo is importand for the fBit information to know the trigger!
	TClonesArray * info_array = (TClonesArray *)sman->FindDataContainer("EventInfo");
	tree->Branch(info_array->GetName(),&info_array);

	TClonesArray * ppac_array = 
		(TClonesArray *)sman->FindDataContainer("BigRIPSPPAC");
	tree->Branch(ppac_array->GetName(),&ppac_array);

	TClonesArray * pla_array = 
		(TClonesArray *)sman->FindDataContainer("BigRIPSPlastic");
	tree->Branch(pla_array->GetName(),&pla_array);

	TClonesArray * ic_array = 
		(TClonesArray *)sman->FindDataContainer("BigRIPSIC");
	tree->Branch(ic_array->GetName(),&ic_array);

	//TClonesArray * fpl_array = 
	//  (TClonesArray *)sman->FindDataContainer("BigRIPSFocalPlane");
	//tree->Branch(fpl_array->GetName(),&fpl_array);
	
	//Dali data
	TClonesArray * dali_array=
		(TClonesArray *)sman->FindDataContainer("DALINaI");
	tree->Branch(dali_array->GetName(),&dali_array);
  
	//PID data
	/*TClonesArray *rips_array = 
	  (TClonesArray *)sman->FindDataContainer("BigRIPSRIPS");
	tree->Branch(rips_array->GetName(),&rips_array); 

	TClonesArray *tof_array  = 
		(TClonesArray *)sman->FindDataContainer("BigRIPSTOF");
	tree->Branch(tof_array->GetName(),&tof_array);   

	TClonesArray *beam_array = 
		(TClonesArray *)sman->FindDataContainer("BigRIPSBeam");	

	*/
	//_________________________________________________________________________
	//Making new branches
	// Coincidence Register
	Int_t neve = 0;
	Int_t triggerbit = 0;

	tree->Branch("triggerbit", &triggerbit,"triggerbit/I");
	tree->Branch("neve", &neve,"neve/I");

	//%%%%%%%%%%%%%%%%%%%%%%
	//F8 target position
	Double_t tgtx = -9999;
	Double_t tgty = -9999;
	Double_t tgta = -9999;
	Double_t tgtb = -9999;

	tree->Branch("tgtx",&tgtx,"tgtx/D");
	tree->Branch("tgty",&tgty,"tgty/D");
	tree->Branch("tgta",&tgta,"tgta/D");
	tree->Branch("tgtb",&tgtb,"tgtb/D");

	Double_t F3X, F5X, F7X, F9X, 
	  F3A,F5A,F7A,F9A, F8T;
	Double_t delta[4];

	tree->Branch("F3X",&F3X,"F3X/D");
	tree->Branch("F3A",&F3A,"F3A/D");
  
	tree->Branch("F5X",&F5X,"F5X/D");
	tree->Branch("F5A",&F5A,"F5A/D");

	tree->Branch("F7X",&F7X,"F7X/D");
	tree->Branch("F7A",&F7A,"F7A/D");

	tree->Branch("F9X",&F9X,"F9X/D");
	tree->Branch("F9A",&F9A,"F9A/D");

	tree->Branch("F8T",&F8T,"F8T/D");

	tree->Branch("delta",delta,"delta[4]/D");

	Double_t F10X = -9999;
	Double_t F10Y = -9999;
	Double_t F10A = -9999;
	Double_t F10B = -9999;

	tree->Branch("F10X",&F10X,"F10X/D");
	tree->Branch("F10Y",&F10Y,"F10Y/D");
	tree->Branch("F10A",&F10A,"F10A/D");
	tree->Branch("F10B",&F10B,"F10B/D");

	Double_t F11X = -9999;
	Double_t F11Y = -9999;
	Double_t F11A = -9999;
	Double_t F11B = -9999;
                                             
	tree->Branch("F11X",&F11X,"F11X/D");
	tree->Branch("F11Y",&F11Y,"F11Y/D");
	tree->Branch("F11A",&F11A,"F11A/D");
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
	Double_t BETA[3];  
	Double_t AOQ[3];   
	Double_t ZET[4];   

	tree->Branch("DELTA",DELTA,"DELTA[3]/D");
	tree->Branch("ANGLE",ANGLE,"ANGLE[3]/D");
	tree->Branch("BRHO",BRHO,"BRHO[3]/D");
	tree->Branch("TOF",TOF,"TOF[3]/D");
	tree->Branch("BETA",BETA,"BETA[3]/D");
	tree->Branch("AOQ",AOQ,"AOQ[3]/D");
	tree->Branch("ZET",ZET,"ZET[4]/D");

	//%%%%%%%%%%%%%%%%%%%%%%%%
	//TKE & Analog sum
	Double_t TKE = -9999;   
	Double_t ASUM = -9999;   

	tree->Branch("TKE",&TKE,"TKE/D");
	tree->Branch("ASUM",&ASUM,"ASUM/D");

	//%%%%%%%%%%%%%%%%%%%%%%
	//DALI
	Int_t dalimultwotime = 0;
	Int_t dalimult = 0;
	Int_t dalitimetruemult = 0;
	Int_t dalimultthres = 0;
	Int_t dalitimetruemultthres = 0;

	tree->Branch("dalimultwotime",&dalimultwotime,"dalimultwotime/I");
	tree->Branch("dalimult",&dalimult,"dalimult/I");
	tree->Branch("dalitimetruemult",&dalitimetruemult,"dalitimetruemult/I");
	tree->Branch("dalimultthres",&dalimultthres,"dalimultthres/I");
	tree->Branch("dalitimetruemultthres",&dalitimetruemultthres,"dalitimetruemultthres/I");

	//%%%%%%%%%%%%%%%%%%%%%%
	//LaBr QTC Array
	Double_t LOWGAIN[8]; 
	Double_t LOWGAINC[8]; 
	Double_t MIDGAIN[8]; 
	Double_t MIDGAINC[8]; 
	Double_t HIGAIN[8];  
	Double_t HIGAINC[8];  
	Double_t QTCTime[8];  

	tree->Branch("LOWGAIN",LOWGAIN,"LOWGAIN[8]/D");
	tree->Branch("LOWGAINC",LOWGAINC,"LOWGAINC[8]/D");
	tree->Branch("MIDGAIN",MIDGAIN,"MIDGAIN[8]/D");
	tree->Branch("MIDGAINC",MIDGAINC,"MIDGAINC[8]/D");
	tree->Branch("HIGAIN",HIGAIN,"HIGAIN[8]/D");
	tree->Branch("HIGAINC",HIGAINC,"HIGAINC[8]/D");
	tree->Branch("QTCTime",QTCTime,"QTCtime[8]/D");

	while( estore->GetNextEvent() ){
	  //while(estore->GetNextEvent() && neve < 100000){
		if(neve%100000==0)
			std::cout << "event: " << neve << std::endl;
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		// Coincidence Register
		triggerbit = -9999;
		//coin->ClearData();
		//coin->LoadData();
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
		F3X = -999; F5X = -999; F7X = -999; F9X = -999; F3A = -999;F5A = -999;F7A = -999;F9A = -999; 
		tfpl = cfpl->FindFocalPlane(3); 
		if(tfpl){TVectorD* vec=tfpl->GetOptVector(); F3X=(*vec)(0); F3A=(*vec)(1);}else cout<< "not found"<<endl;
		tfpl = cfpl->FindFocalPlane(5); 
		if(tfpl){TVectorD* vec=tfpl->GetOptVector(); F5X=(*vec)(0); F5A=(*vec)(1);}
		tfpl = cfpl->FindFocalPlane(7); 
		if(tfpl){TVectorD* vec=tfpl->GetOptVector(); F7X=(*vec)(0); F7A=(*vec)(1);}
		tfpl = cfpl->FindFocalPlane(9); 
		if(tfpl){TVectorD* vec=tfpl->GetOptVector(); F9X=(*vec)(0); F9A=(*vec)(1);}
		F8T = cplastic -> FindPlastic("F8pl") -> GetTime();
		//Delta 
		delta[0] = rips3to5->GetDelta();
		delta[1] = rips5to7->GetDelta();
		delta[2] = rips8to9->GetDelta();
		delta[3] = rips10to11->GetDelta();


		// for skip F3 ppac tracking --->
		TArtFocalPlane *fpl = cfpl->FindFocalPlane(3);
		if( fpl ){
			TVectorD *vec=fpl->GetOptVector();
			Double_t vec_tmp[4]={0,0,0,0};
			vec->SetElements(vec_tmp);
		}
		// <-----


                //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		//Reconstructiong the PID
		recopid->ClearData();
		recopid->ReconstructData();

		for( int i=0; i<3; i++ ){
			DELTA[i] = -9999;
			ANGLE[i] = -9999;
			BRHO[i] = -9999;
			TOF[i] = -9999;
			BETA[i] = -9999;
			AOQ[i] = -9999;
		}
		for( int i=0; i<4; i++ ){
			ZET[i] = -9999;
		}

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

		BETA[0] = beam_br -> GetBeta();
		BETA[1] = beam_zd_810 -> GetBeta();
		AOQ[0] = beam_br -> GetAoQ();
		AOQ[1] = beam_zd_810 -> GetAoQ();
		AOQ[2] = beam_zd_1011 -> GetAoQ();
		ZET[0] = beam_br -> GetZet();
		ZET[1] = beam_zd_810 -> GetZet();
		ZET[2] = ZET[1] + 0.3;//??
		//ZET[2] = ZET[1] - (BETA[1]*BETA[1]*1291.43-BETA[1]*878.207+170.106)+38;
		ZET[3] = ZET[2] * 1.5352 - 20.292;

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
		TArtFocalPlane *fpl10; TVectorD *vec10;
		fpl10 = cfpl -> FindFocalPlane(10);
		if( fpl10 ){vec10 = fpl10 -> GetOptVector();
			F10X = (*vec10)(0); F10A = (*vec10)(1); F10Y = (*vec10)(2); F10B = (*vec10)(3);
		}
		TArtFocalPlane *fpl11; TVectorD *vec11;
		fpl11 = cfpl -> FindFocalPlane(11);
		if( fpl11 ){vec11 = fpl11 -> GetOptVector();
			F11X = (*vec11)(0); F11A = (*vec11)(1); F11Y = (*vec11)(2); F11B = (*vec11)(3);
		} 

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
		Int_t tref = -9000;
	
		for( int i=0; i<8; i++ ){
			LOWGAIN[i] = -9999;
			MIDGAIN[i] = -9999;
			HIGAIN[i] = -9999;
		}

		Int_t qtc_dummy[2][40];
		for( int i=0; i<40; i++ ){
			for( int j=0; i<2; i++ ){
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
						if(geo==9) qtc_dummy[edge][ch] = val;
						if(geo ==9 && ch==127 && edge == 0) tref = val;					
					}
				}
			}
			Double_t beta = 0.5709654;
			Double_t theta = 30.0 * 3.141592 / 180.;
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
			} 

		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		//Making DALI
		dalicalib -> ClearData();
		dalicalib -> SetPlTime(cplastic -> FindPlastic("F7pl") -> GetTime());
		//Add above to remove F8plastic tof.
		dalicalib -> ReconstructData();

		dalimultwotime = dalicalib -> GetMultWithoutT();
		dalimult = dalicalib -> GetMult();
		dalitimetruemult = dalicalib -> GetTimeTrueMult();
		dalimultthres = dalicalib -> GetMultThres();
		dalitimetruemultthres = dalicalib -> GetTimeTrueMultThres();

		tree -> Fill();
		neve++;
	}
	cout<<"Writing the tree."<<endl;

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
