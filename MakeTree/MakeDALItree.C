
#include "TArtStoreManager.hh"
#include "TArtEventStore.hh"
#include "TArtCore.hh"

#include <iostream>
#include <string>
#include <time.h>

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TStyle.h>
#include <TSystem.h>

#include "TArtCalibPID.hh"
#include "TArtCalibPPAC.hh"
#include "TArtCalibPlastic.hh"
#include "TArtCalibIC.hh"
#include "TArtCalibFocalPlane.hh"
#include "TArtEventInfo.hh"
#include "TArtEventStore.hh"
#include "TArtPPAC.hh"
#include "TArtPPACPara.hh"
#include "TArtPlastic.hh"
#include "TArtIC.hh"
#include "TArtFocalPlane.hh"
#include "TArtBigRIPSParameters.hh"

#include "TArtStoreManager.hh"
#include "TClonesArray.h"

#include "TArtBigRIPSParameters.hh"
#include "TArtDALIParameters.hh"

#include "TArtCalibPID.hh"
#include "TArtCalibPlastic.hh"
#include "TArtCalibDALI.hh"

#include "TArtPlastic.hh"
#include "TArtDALINaI.hh"
#include "TArtDALINaIPara.hh"

void MakeDALItree(Int_t run){

	//TArtStoreManager * sman = TArtStoreManager::Instance();
	TArtStoreManager * sman = new TArtStoreManager();// if compile
	TArtEventStore *estore = new TArtEventStore();

	//TArtCalibCoin *coin = new TArtCalibCoin();

	//estore->Open(1);
	//if(run < 10)
	  estore -> Open(Form("ridf/ggdaq04/calibration%04d.ridf",run));
	  //else estore ->Open(Form("ridf/calibration00%i.ridf",run));
		//TArtDALIParameters *para = TArtDALIParameters::Instance();
		TArtDALIParameters *para = new TArtDALIParameters();// if compile
		para -> LoadParameter((char*)"db/141026/withTartget/DALI.xml");

		TFile* f = new TFile(Form("rootfile/calibration/calibration%04d.root",run),"RECREATE");
		TArtCalibDALI *calib_dali = new TArtCalibDALI();
                /*calib_dali->ReconstructData();
                TClonesArray * dali_array = (TClonesArray *)sman->FindDataContainer("DALINaI");*/

		TTree* tree = new TTree("eventtree","Event Tree");
                
		cout << "tree created" << endl;
		
		Double_t fTRaw[180] = {-9999};
		Double_t fTCal[180] = {-9999};
		Double_t fQRaw[180] = {-9999};
		Double_t fQCal[180] = {-9999};
		Double_t fQDopcal[180] = {-9999};

		Int_t fID[180] = {-9999};
		Double_t fx[180] = {-9999};
		Double_t fy[180] = {-9999};
		Double_t fz[180] = {-9999};
		Double_t ftheta[180] = {-9999};

		Int_t flayer[180] = {-9999};
		Int_t fneve = 0;
		Int_t multi = 0;
	
		//Int_t ftriggerbit = 0;
		//Int_t Triggerbit = 0;

		tree->Branch("fTRaw", fTRaw,"fTRaw[180]/D");
		tree->Branch("fTCal", fTCal,"fTCal[180]/D");
		tree->Branch("fQRaw", fQRaw,"fQRaw[180]/D");
		tree->Branch("fQCal", fQCal,"fQCal[180]/D");
		tree->Branch("fQDopcal", fQDopcal,"fQDopcal[180]/D");
		
		tree->Branch("fID", fID,"fID[180]/I");
		tree->Branch("fx", fx,"fx[180]/D");
		tree->Branch("fy", fy,"fy[180]/D");
		tree->Branch("fz", fz,"fz[180]/D");
		tree->Branch("ftheta", ftheta,"ftheta[180]/D");

		tree->Branch("flayer", flayer,"flayer[180]/I");
		tree->Branch("fneve", &fneve,"fneve/I");
		tree->Branch("multi", &multi,"multi/I");

		//tree->Branch("Triggerbit", &Triggerbit,"Triggerbit/I");
		
		/*		
		TArtCalibDALI *calib_dali = new TArtCalibDALI();
		calib_dali->ReconstructData();
		TClonesArray * dali_array = (TClonesArray *)sman->FindDataContainer("DALINaI");*/
		//calib_dali->ReconstructData();		

		//TClonesArray * triggerbit_array = (TClonesArray *)estore ->GetEventInfoArray();
		//TArtEventInfo *triggerinfo = triggerbit_array -> At(0);

		int neve =0;
		//while( estore -> GetNextEvent() ){
		while( estore->GetNextEvent() && neve < 1e6){
		  if(neve%(Int_t)1e5==0){ cout << "neve: " << neve << endl; }

		  calib_dali->ReconstructData();
		  TClonesArray * dali_array = (TClonesArray *)sman->FindDataContainer("DALINaI");
		  //estore->GetNextEvent();
		  /*for(Int_t jj = 0; jj < rawevent->GetNumSeg(); jj++){
		    TArtRawSegmentObject *seg = rawevent->GetSegment(i);
		    Int_t detid = seg -> GetDetector();
			/*------------------------------------------------/
			coin->ClearData();
			coin->LoadData();
			ftriggerbit = triggerinfo -> GetTriggerBit();
			for(Int_t i = 0; i < 7; i++){
				if ((ftriggerbit >> i) & 0x1)
					Triggerbit = i+1;
			}
			/------------------------------------------------*/
		    multi = calib_dali ->GetNumNaI();
		    //	if(calib_dali->GetNumNaI() > 0)
		    //	  cout<<"Number of Dalis "<<calib_dali->GetNumNaI()<<endl;


			for(Int_t i = 0 ; i < multi /*calib_dali -> GetNumNaI() */; i++){				

				TArtDALINaI *nai = (TArtDALINaI*)calib_dali -> GetNaI(i);
				TArtDALINaIPara *naip = (TArtDALINaIPara*)calib_dali -> GetNaIPara(i);

				Int_t id = (Int_t)nai -> GetID() - 1;
				flayer[i] = (Int_t)nai -> GetLayer();

				fTRaw[i] = (Double_t)nai -> GetRawTDC();
				//if(fTRaw[id] > 0)cout<<"TRAW "<<fTRaw[id]<<endl;
				fTCal[i] = (Double_t)nai -> GetTime();
				fQRaw[i] = (Double_t)nai -> GetRawADC();
				fQCal[i] = (Double_t)nai -> GetEnergy();
				fQDopcal[i] = (Double_t)nai -> GetDoppCorEnergy();

				fID[i] = (Int_t)nai -> GetID();
				fx[i] = (Double_t)nai -> GetXPos();
				fy[i] = (Double_t)nai -> GetYPos();
				fz[i] = (Double_t)nai -> GetZPos(); 
				ftheta[i] = (Double_t)nai -> GetTheta(); 

			}

				fneve = neve;
				tree->Fill();
		
				// initialize			
				for(Int_t i = 0; i<180; i++){
					fTRaw[i] = -9999;
					fTCal[i] = -9999;
					fQRaw[i] = -9999;
					fQCal[i] = -9999;
					fQDopcal[i] = -9999;

					fID[i] = -9999;
					fx[i] = -9999;
					fy[i] = -9999;
					fz[i] = -9999;
					ftheta[i] = -9999;
					
					flayer[i] = -9999;
					fneve = 0;
					multi = 0;
				}
				//------------------------------------------------
				estore -> ClearData();
				calib_dali -> ClearData();
				++neve;
				//------------------------------------------------
				//}
		}
			f->Write();
			//for ridf loop		
			//	}
}
