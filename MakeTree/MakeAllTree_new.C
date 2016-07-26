#include <iostream>
#include <string>
#include <time.h>

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TObject.h>

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

#include "signal.h"

// BigRIPS parameters
#define tofoffbr  302// TOF3-7 offset VME
#define tofoffzd  -749 // TOF7-11 offset VME

#define brho35 6.9181 
#define brho57 6.4551
#define brho89 2.285
#define brhoTE 2.348 //run185

using namespace std;

// function to exit loop at keyboard interrupt. 
bool stoploop = false;
void stop_interrupt(){
  printf("keyboard interrupt\n");
  stoploop = true;
}

void MakeAllTree_new(Int_t run){

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
	estore->Open(Form("ridf/ggdaq04/calibration%04d.ridf",run));

	//Coincidence register LUPO
	TArtCalibCoin *coin = new TArtCalibCoin();
	TClonesArray * triggerbit_array = (TClonesArray *)estore -> GetEventInfoArray();
	TArtEventInfo *triggerinfo = (TArtEventInfo *)triggerbit_array -> At(0);
	
	// Create BigRIPSParameters to get Plastics, PPACs, ICs and FocalPlanes parameters from ".xml" files
	TArtBigRIPSParameters *para = TArtBigRIPSParameters::Instance();
	para -> LoadParameter("db/BigRIPSPPAC.xml");
	para -> LoadParameter("db/BigRIPSPlastic.xml");
	para -> LoadParameter("db/BigRIPSIC.xml");
	para -> LoadParameter("db/FocalPlane.xml");

	//Create DALI Parameters
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
	TArtRIPS *rips3to5 = recopid -> DefineNewRIPS(3,5,"matrix/mat1.mat",brho35); // F3 - F5
	TArtRIPS *rips5to7 = recopid -> DefineNewRIPS(5,7,"matrix/mat2.mat",brho57); // F5 - F7
	TArtRIPS *rips8to10 = recopid -> DefineNewRIPS(8,10,"matrix/F8F10_LargeAccAchr.mat",brho89); // F9 - F11  
	TArtRIPS *rips10to11 = recopid -> DefineNewRIPS(10,11,"matrix/F10F11_LargeAccAchr.mat",brhoTE); // F9 - F11  

	// Reconstruction of TOF DefineNewTOF(fisrt plane, second plane, time offset)
	TArtTOF *tof3to7  = recopid -> DefineNewTOF("F3pl","F7pl",tofoffbr,5); // F3 - F7
	TArtTOF *tof7to11 = recopid -> DefineNewTOF("F7pl","F11pl-1",tofoffzd,9); // F8 - F11

	// Reconstruction of IC observables for ID
	TArtBeam *beam_br = recopid -> DefineNewBeam(rips3to5,rips5to7,tof3to7,"F7IC");
	TArtBeam *beam_zd_810 = recopid -> DefineNewBeam(rips8to10,tof7to11,"F11IC");
	TArtBeam *beam_zd_1011 = recopid -> DefineNewBeam(rips10to11,tof7to11,"F11IC");

	TClonesArray *rips_array = recopid -> GetRIPSArray();
	TClonesArray *tof_array  = recopid -> GetTOFArray();
	TClonesArray *beam_array = recopid -> GetBeamArray();

	//Define Out put tree	
	TFile *fout = new TFile(Form("run%04d.root",run),"RECREATE");
	TTree *tree = new TTree("tree","tree");
/*
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
*/
	//Dali data
	TClonesArray * dali_array=
		(TClonesArray *)sman->FindDataContainer("DALINaI");
	tree->Branch(dali_array->GetName(),&dali_array);
/*  
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

	while(estore->GetNextEvent() && neve < 100000){
		if(neve%10000==0)
			std::cout << "event: " << neve << std::endl;
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		// Coincidence Register
		triggerbit = -9999;
		coin->ClearData();
		coin->LoadData();
		Int_t ftriggerbit = triggerinfo -> GetTriggerBit();
		for(Int_t i = 0; i < 7; i++){
			if ((ftriggerbit >> i) & 0x1)
				triggerbit = i+1;
		}
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		//Making the BigRIPS tree calibration
		cpid -> ClearData();
		cpid -> ReconstructData();

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
		ZET[2] = ZET[1] - (BETA[1]*BETA[1]*1291.43-BETA[1]*878.207+170.106)+38;
		ZET[3] = ZET[2] * 1.5352 - 20.292;

                //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		//BigRIPS F8 target profile 
		tgtx = -9999;
		tgty = -9999;
		tgta = -9999; 
		tgtb = -9999;

		TArtFocalPlane *fpl; TVectorD *vec;
		fpl = cfpl -> FindFocalPlane(8);
		if( fpl ){vec = fpl -> GetOptVector();
			tgtx = (*vec)(0); tgta = (*vec)(1); tgty = (*vec)(2); tgtb = (*vec)(3);
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
