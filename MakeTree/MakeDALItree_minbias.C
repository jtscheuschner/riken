#include "TArtStoreManager.hh"
#include "TArtEventStore.hh"
#include "TArtDALIParameters.hh"
#include "TArtCalibDALI.hh"
#include "TSystem.h"
#include "TTree.h"
#include "TFile.h"
#include "TClonesArray.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
//#include "segidlist.h"

// function to exit loop at keyboard interrupt. 
/*
bool stoploop = false;
void stop_interrupt()
{
  printf("keyboard interrupt\n");
  stoploop = true;
}
*/
void MakeDALItree_minbias(Int_t run){

  gSystem->Load("libXMLParser.so");
  gSystem->Load("libanacore.so");

  TArtEventStore *estore = new TArtEventStore();
  TArtRawEventObject *rawevent = estore->GetRawEventObject();
  estore -> Open(Form("ridf/ggdaq04/calibration%04d.ridf",run));

  TFile* f = new TFile(Form("rootfiles/calibration/calibration%04d.root",run),"RECREATE");
  TTree* tree = new TTree("eventtree","Event Tree");
  Double_t fTRaw[180] = {-9999};
  Double_t fQRaw[180] = {-9999};
  Int_t fneve = -9999;
  tree->Branch("fTRaw", fTRaw,"fTRaw[180]/D");
  tree->Branch("fQRaw", fQRaw,"fQRaw[180]/D");
  tree->Branch("fneve", &fneve,"fneve/I");

/*
  TArtDALIParameters *dpara = TArtDALIParameters::Instance();
  dpara->LoadParameter("db/DALI.xml");
  TArtCalibDALI *dalicalib= new TArtCalibDALI();
*/
//  TCanvas *c1 = new TCanvas("c1","ADC",700,500);
//  c1->Draw();
//  TH2F* h1 = new TH2F("module","module",50,-0.5,49.5,50,-0.5,49.5);
//  h1->Draw();
//  c1->Modified();
//  c1->Update();

  TCanvas *cDALI = new TCanvas("DALI2","ADC all",1400,800);
  cDALI->Divide(3,3);
  TH2I *adc1=new TH2I("adc1","adc(V785) geo=10",33,-0.5,33.5,512,0,8096); 
  TH2I *adc2=new TH2I("adc2","adc(V785) geo=11",33,-0.5,33.5,512,0,8096); 
  TH2I *adc3=new TH2I("adc3","adc(V785) geo=12",33,-0.5,33.5,512,0,8096); 
  TH2I *adc4=new TH2I("adc4","adc(V785) geo=13",33,-0.5,33.5,512,0,8096); 
  TH2I *adc5=new TH2I("adc5","adc(V785) geo=14",33,-0.5,33.5,512,0,8096); 
  TH2I *adc6=new TH2I("adc6","adc(V785) geo=15",33,-0.5,33.5,512,0,8096); 
  TH2I *adc7=new TH2I("adc7","adc(V785) geo=16",33,-0.5,33.5,512,0,8096); 
  TH2I *adc8=new TH2I("adc8","adc(V785) geo=17",33,-0.5,33.5,512,0,8096); 
  
  cDALI->cd(1);
  adc1->Draw("colz");
  cDALI->cd(2);
  adc2->Draw("colz");
  cDALI->cd(3);
  adc3->Draw("colz");
  cDALI->cd(4);
  adc4->Draw("colz");
  cDALI->cd(5);
  adc5->Draw("colz");
  cDALI->cd(6);
  adc6->Draw("colz");
  cDALI->cd(7);
  adc7->Draw("colz");
  cDALI->cd(8);
  adc8->Draw("colz");

  cDALI->Modified();
  cDALI->Update();

  //dali calib
  //TTree *tree = new TTree("tree","tree");
  // define data nodes which are supposed to be dumped to tree 
  //TClonesArray * info_array = (TClonesArray *)sman->FindDataContainer("EventInfo");
  //std::cout<<info_array->GetName()<<std::endl;
  //tree->Branch(info_array->GetName(),&info_array);
  //TClonesArray * dali_array=
  //   (TClonesArray *)sman->FindDataContainer("DALINaI");
  //tree->Branch(dali_array->GetName(),&dali_array);

  Int_t neve = 0;
  while(estore -> GetNextEvent()){
    for(Int_t i=0; i < rawevent -> GetNumSeg(); i++){
      TArtRawSegmentObject *seg = rawevent->GetSegment(i);
      Int_t detid = seg -> GetDetector();
      Int_t modid = seg -> GetModule();
      //h1->Fill(detid,modid);
      //add by Shiga
      //Int_t device = seg -> GetDevice();
      //Int_t fp = seg -> GetFP();
      //Int_t detector = seg -> GetDetector();
      //Int_t module = seg -> GetModule();
      //if(DALI==device&&DALIA==detector){
	//cout << "    seg:"<< i <<" dev:"<< device 
	//     << " fp:"<<fp<< " det:"<<detector<< " mod:"<<module
	//     << " #data=" << seg->GetNumData() << endl;
	for(Int_t j=0; j < seg -> GetNumData(); j++){
	  TArtRawDataObject *d = seg -> GetData(j);
	  Int_t geo = d -> GetGeo();
	  Int_t ch = d -> GetCh();
	  unsigned int val = d -> GetVal();
	  //cout << "       geo:" << geo 
	  //     << " ch:" << ch << " val:" << val << endl;
	  //ntp->Fill((float)geo,(float)ch,(float)val);

          //fill data
          if(geo==10)adc1 -> Fill(ch,val);
	  if(geo==10){//ADC
		  for(Int_t k = 0; k < 180; k++){
			  if(ch == k) fQRaw[k] = val;
		  }
	  };
	  if(geo==11)adc2->Fill(ch,val);
	  if(geo==11){//ADC
		  for(Int_t k = 0; k < 180; k++){
			  if(ch == k) fQRaw[k] = val;
		  }
	  };
	  if(geo==12)adc3->Fill(ch,val);
	  if(geo==12){//ADC
		  for(Int_t k = 0; k < 180; k++){
			  if(ch == k) fQRaw[k] = val;
		  }
	  };
	  if(geo==13)adc4->Fill(ch,val);
	  if(geo==13){//ADC
		  for(Int_t k = 0; k < 180; k++){
			  if(ch == k) fQRaw[k] = val;
		  }
	  };
	  if(geo==14)adc5->Fill(ch,val);
	  if(geo==14){//ADC
		  for(Int_t k = 0; k < 180; k++){
			  if(ch == k) fQRaw[k] = val;
		  }
	  };
	  if(geo==15)adc6->Fill(ch,val);
	  if(geo==15){//ADC
		  for(Int_t k = 0; k < 180; k++){
			  if(ch == k) fQRaw[k] = val;
		  }
	  };
	  if(geo==16)adc7->Fill(ch,val);
	  if(geo==16){//ADC
		  for(Int_t k = 0; k < 180; k++){
			  if(ch == k) fQRaw[k] = val;
		  }
	  };
	  if(geo==17)adc8->Fill(ch,val);
	  if(geo==17){//ADC
		  for(Int_t k = 0; k < 180; k++){
			  if(ch == k) fQRaw[k] = val;
		  }
	  };
	}
      //} 
    }
    fneve = neve;
    tree->Fill();
    rawevent->Clear();
    neve ++;

    if(neve%100==0){
      //h1->Draw();
      //c1->Update();
      cDALI->cd(1);
      adc1->Draw("colz");
      cDALI->cd(2);
      adc2->Draw("colz");
      cDALI->cd(3);
      adc3->Draw("colz");
      cDALI->cd(4);
      adc4->Draw("colz");
      cDALI->cd(5);
      adc5->Draw("colz");
      cDALI->cd(6);
      adc6->Draw("colz");
      cDALI->cd(7);
      adc7->Draw("colz");
      cDALI->Update();
    }
  }
 f -> Write();
}
