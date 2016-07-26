#include "signal.h"

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
#include <TRandom3.h>
// function to exit loop at keyboard interrupt. 
bool stoploop = false;
void stop_interrupt()
{
  printf("keyboard interrupt\n");
  stoploop = true;
}

void MakeDaliTree(char *infile, char *outfile="dalia.root"){
  //  void MakeDaliTree(char *outfile="dalia.root"){
  gSystem->Load("libXMLParser.so");
  gSystem->Load("libanaroot.so");

  TArtStoreManager * sman = TArtStoreManager::Instance();
  TArtEventStore *estore = new TArtEventStore();
  TArtRawEventObject* frawevent = (TArtRawEventObject *)sman->FindDataContainer("RawEvent");

  estore->SetInterrupt(&stoploop);
  estore->Open(infile);
  //estore->Open();

  TArtDALIParameters *dpara = TArtDALIParameters::Instance();
  dpara->LoadParameter("db/DALI.xml");
  TArtCalibDALI *dalicalib= new TArtCalibDALI();

  TFile *fout = new TFile(outfile,"RECREATE");
  TTree *tree = new TTree("tree","tree");

  // define data nodes which are supposed to be dumped to tree 
  TClonesArray * info_array = (TClonesArray *)sman->FindDataContainer("EventInfo");
  std::cout<<info_array->GetName()<<std::endl;
  tree->Branch(info_array->GetName(),&info_array);

  TClonesArray * dali_array=
     (TClonesArray *)sman->FindDataContainer("DALINaI");
  tree->Branch(dali_array->GetName(),&dali_array);

  Int_t dalimult = 0;
  const UInt_t fn_labr = 8;
  Int_t dalitimetruemult = 0;
  Float_t HIGAIN[fn_labr]={0x0},MIDGAIN[fn_labr]={0x0},LOWGAIN[fn_labr]={0x0}, QTCTime[fn_labr]={0x0};
  Float_t ftime = 0.09766;
  UInt_t fQTC[3][fn_labr] = {28, 29, 30, 1, 31, 4, 3, 2,
		       24, 25, 26, 6, 27, 9, 8, 7,
		       20, 21, 22,11, 23,14,13,12};//?logic maybe 10,13,12,11??

  tree->Branch("dalimult",&dalimult,"dalimult/I");
  tree->Branch("dalitimetruemult",&dalitimetruemult,"dalitimetruemult/I");

  int neve = 0;
  //  while(estore->GetNextEvent() && neve<5e6){
  while(estore->GetNextEvent() ){
    if(neve%100000==0)
      std::cout << "event: " << neve << std::endl;
    //if(neve>1e7)break;
    dalicalib->ClearData();
    dalicalib->ReconstructData();

    dalimult = dalicalib->GetMult();
    dalitimetruemult = dalicalib->GetTimeTrueMult();
  
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
      QTCTime[ilabr] = (qtc_dummy[0][fQTC[1][ilabr]] - tref)*ftime;
    }
    tree->Fill();
    neve ++;
  }
  fout->Write();
  fout->Close();

}

