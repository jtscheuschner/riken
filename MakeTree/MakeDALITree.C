#include "TMath.h"
#include <TROOT.h>
#include <TFile.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TObject.h>
#include <TRandom3.h>
#include <iostream>
#include <string>
#include <time.h>
#include "TCanvas.h"

#include "TArtStoreManager.hh"
#include "TArtEventStore.hh"
#include "TArtDALIParameters.hh"
#include "TArtCalibDALI.hh"
#include "TSystem.h"
#include "TTree.h"
#include "TFile.h"
#include "TClonesArray.h"
#include <sstream>
#include "signal.h"

// function to exit loop at keyboard interrupt. 
bool stoploop = false;
void stop_interrupt()
{
  printf("keyboard interrupt\n");
  stoploop = true;
}

void MakeDALITree(UInt_t run){
  gSystem->Load("libXMLParser.so");
  gSystem->Load("libanaroot.so");

  TArtStoreManager * sman = TArtStoreManager::Instance();
  TArtEventStore *estore = new TArtEventStore();
  estore->SetInterrupt(&stoploop);
  //estore->Open(0);
  Char_t* name = "calibration";
  UInt_t run_=run;
  if(run>200){
    name = "calib2XO";
    run_=run-200;
  }else if(run>100){
    name = "labr14Ni";
    run_=run-100;
  }
  estore -> Open(Form("data/ridf/ggdaq04/%s%04d.ridf",name,run_));

  //  estore -> Open(Form("ridf/ggdaq04/labr14Ni%04d.ridf",run));

  TArtDALIParameters *dpara = TArtDALIParameters::Instance();
  dpara->LoadParameter("db/DALI.xml");
  TArtCalibDALI *dalicalib= new TArtCalibDALI();

  //check if the tree is readable, otherwise change here to array and save once in a while during the loop and delete the tree, start again.
  //TFile *fout[100];
  //TTree *tree;
   
  Char_t* path = "data/Histo/calib/";
  Char_t* rename;
  if(run<34)
    rename = "before_empty/run";
  else if(run<39)
    rename = "before_Sn132/run";
  else if(run<43)
    rename = "before_Sn128/run";
  else if(run<46)
    rename = "during_Sn128/run";
  else if(run<52)
    rename = "after_Sn128/run";
  else 
    rename = "after_exp/run";
  TFile *fout = new TFile(Form("%s%s%04d.root",path,rename,run),"RECREATE");
  TTree *tree = new TTree("tree","tree");
  
  // define data nodes which are supposed to be dumped to tree 
  TClonesArray * info_array = (TClonesArray *)sman->FindDataContainer("EventInfo");
  TArtRawEventObject* frawevent = (TArtRawEventObject *)sman->FindDataContainer("RawEvent");
  std::cout<<info_array->GetName()<<std::endl;
  //std::cout<<info_array->GetDate()<<std::endl;
  tree->Branch(info_array->GetName(),&info_array);

  TClonesArray * dali_array=
     (TClonesArray *)sman->FindDataContainer("DALINaI");
  tree->Branch(dali_array->GetName(),&dali_array);

  Int_t MIDGAIN[8] = {0x0};
  Int_t LOWGAIN[8] = {0x0};
  Int_t HIGAIN[8]  = {0x0};

  tree->Branch("MIDGAIN",MIDGAIN,"MIDGAIN[8]/I");
  tree->Branch("LOWGAIN",LOWGAIN,"LOWGAIN[8]/I");
  tree->Branch("HIGAIN" ,HIGAIN ,"HIGAIN[8]/I");

  int neve = 0;
  Int_t ii = -1;
  UInt_t  fQTC[3][8] = {28, 29, 30, 1, 31, 4, 3, 2,
			24, 25, 26, 6, 27, 9, 8, 7,
			20, 21, 22,11, 23,14,13,12};
  TRandom3 rndm;

  while( estore->GetNextEvent() ){
    //while(neve<9000000){
    //	  estore->GetNextEvent();
    if(neve%500000==0){
      std::cout << "event: " << neve << std::endl;
      //not tested
/*
      if(neve%(Int_t)1e7==0){
	if(neve!=0){
	  fout[ii]->Write();
	  delete tree;
	}    
	ii++;
	if(ii > 99){
	  std::cout<< "file array is not big enough please increase the array size"<<std::endl;
	  return;
	}
	tree = new TTree("tree","tree");
	fout[ii] = new TFile(Form("rootfiles/calibration/calibration%04d_%i.root",run,ii),"RECREATE");
	}*/
    }

    dalicalib->ClearData();
    dalicalib->ReconstructData();

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
  
    for(short ilabr=0; ilabr<8; ilabr++){
      HIGAIN[ilabr]  = (qtc_dummy[1][ fQTC[0][ilabr] ] - qtc_dummy[0][ fQTC[0][ilabr] ]);
      MIDGAIN[ilabr] = (qtc_dummy[1][ fQTC[1][ilabr] ] - qtc_dummy[0][ fQTC[1][ilabr] ]);
      LOWGAIN[ilabr] = (qtc_dummy[1][ fQTC[2][ilabr] ] - qtc_dummy[0][ fQTC[2][ilabr] ]);
    }
    tree->Fill();
    neve ++;
  }
  fout->Write();
  fout->Close();

 

}

