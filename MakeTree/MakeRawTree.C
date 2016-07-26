#include "TArtStoreManager.hh"
#include "TArtEventStore.hh"
#include "TArtRawEventObject.hh"

#include "TSystem.h"
#include "TTree.h"
#include "TFile.h"

void MakeRawTree(char *infile, char *outfile="raw.root"){

  gSystem->Load("libanaroot.so");

  TArtStoreManager *sman = TArtStoreManager::Instance();
  
  TArtEventStore *estore = new TArtEventStore();
  estore->Open(infile);
  //estore->Open();  //for online

  TArtRawEventObject *rawevent = estore->GetRawEventObject();

  TFile *fout = new TFile(outfile,"RECREATE");
  TTree *tree = new TTree("tree","tree");
  
  tree->Branch("rawdata",&rawevent);

  int neve = 0;
 
  cout<<"Starting loop."<<endl;

  while(estore->GetNextEvent()&&neve<2000){
    if(neve%1000==0) cout << "Event:" << neve << endl;

    tree->Fill();
    estore->ClearData();
    neve ++;
  }

  fout->Write();
  fout->Close();
}
