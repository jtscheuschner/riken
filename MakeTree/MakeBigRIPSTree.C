#include "TArtStoreManager.hh"
#include "TArtEventStore.hh"
#include "TArtBigRIPSParameters.hh"
#include "TArtDALIParameters.hh"
#include "TArtCalibPID.hh"
#include "TArtCalibDALI.hh"
#include "TArtCalibPPAC.hh"
#include "TArtCalibPlastic.hh"
#include "TArtEventInfo.hh"
#include "TArtPlastic.hh"
#include "TArtPPAC.hh"
#include "TArtRecoPID.hh"
#include "TArtRecoRIPS.hh"
#include "TArtRecoTOF.hh"
#include "TArtRecoBeam.hh"

#include "TSystem.h"
#include "TTree.h"
#include "TFile.h"

#include "TClonesArray.h"


#include "TVector3.h"

#include "signal.h"

using namespace std;

// function to exit loop at keyboard interrupt. 
bool stoploop = false;
void stop_interrupt(){
  printf("keyboard interrupt\n");
  stoploop = true;
}

void MakeBigRIPSTree(char *infile, char *outfile="bigrips.root"){

  //  signal(SIGINT,stop_interrupt); // CTRL + C , interrupt

  // Create StoreManager both for calibration "TArtCalib..." and treatment "TArtReco..."
  //------------------------------------------------------------------------------------
  TArtStoreManager * sman = TArtStoreManager::Instance();
  
  // Create EventStore to control the loop and get the EventInfo
  //------------------------------------------------------------
  TArtEventStore *estore = new TArtEventStore();
  estore->SetInterrupt(&stoploop); 
  estore->Open(infile);
  std::cout<<"estore ->"<< infile <<std::endl;
  
  // Create BigRIPSParameters to get Plastics, PPACs, ICs and FocalPlanes parameters from ".xml" files
  //--------------------------------------------------------------------------------------------------	
  TArtBigRIPSParameters *para = TArtBigRIPSParameters::Instance();
  para->LoadParameter("db/BigRIPSPPAC.xml");
  para->LoadParameter("db/BigRIPSPlastic.xml");
  para->LoadParameter("db/BigRIPSIC.xml");
  para->LoadParameter("db/FocalPlane.xml");
  para->SetFocusPosOffset(8,138.5);

  // Create CalibPID to get and calibrate raw data ( CalibPID -> 
  //[CalibPPAC , CalibIC, CalibPlastic , CalibFocalPlane] 
  TArtCalibPID *brcalib     = new TArtCalibPID();
  TArtCalibPPAC *ppaccalib  = brcalib->GetCalibPPAC();
  TArtCalibPlastic *plasticcalib = brcalib->GetCalibPlastic();

  // Create RecoPID to get calibrated data and to reconstruct TOF, AoQ, Z, ... (RecoPID -> 
  //[ RecoTOF , RecoRIPS , RecoBeam] )
  TArtRecoPID *recopid = new TArtRecoPID();

  para->PrintListOfPPACPara();
  //return;
 
  // Definition of observables we want to reconstruct
  std::cout << "Defining bigrips parameters" << std::endl;
  TArtRIPS *rips3to5 = recopid->DefineNewRIPS(3,5,"matrix/mat1.mat",7.527); // F3 - F5
  TArtRIPS *rips5to7 = recopid->DefineNewRIPS(5,7,"matrix/mat2.mat",7.296); // F5 - F7
  TArtRIPS *rips8to9 = recopid->DefineNewRIPS(8,9,"matrix/F8F9_LargeAccAchr.mat",7.154); // F8 - F9
  TArtRIPS *rips9to11 = recopid->DefineNewRIPS(9,11,"matrix/F9F11_LargeAccAchr.mat",7.1529); // F9 - F11  
  // Reconstruction of TOF DefineNewTOF(fisrt plane, second plane, time offset)
  TArtTOF *tof3to7  = recopid->DefineNewTOF("F3pl","F7pl",234.7,5); // F3 - F7
  TArtTOF *tof8to11 = recopid->DefineNewTOF("F8pl","F11pl-1",199.4,9); // F8 - F11
  // Reconstruction of IC observables for ID
  //TArtBeam *beam_br_37 = recopid->DefineNewBeam(rips3to5,rips5to7,tof3to7,"F7IC");
  TArtBeam *beam_br_35 = recopid->DefineNewBeam(rips3to5,rips5to7,tof3to7,"F7IC");
  TArtBeam *beam_br_57 = recopid->DefineNewBeam(rips5to7,tof3to7,"F7IC");   
  TArtBeam *beam_zd_89 = recopid->DefineNewBeam(rips8to9,tof8to11,"F11IC");
  TArtBeam *beam_zd_911 = recopid->DefineNewBeam(rips9to11,tof8to11,"F11IC");

  //to get trigger pattern
  TArtEventInfo *evtinfo=new TArtEventInfo();
  int trg=-777;
  cout<<"trigger "<<trg<<endl;

  TFile *fout = new TFile(outfile,"RECREATE");
  TTree *tree = new TTree("tree","tree");

  // define data nodes which are supposed to be dumped to tree 
  TClonesArray * info_array = (TClonesArray *)sman->FindDataContainer("EventInfo");
  std::cout<<info_array->GetName()<<std::endl;
  tree->Branch(info_array->GetName(),&info_array);
  
  TClonesArray * ppac_array = 
    (TClonesArray *)sman->FindDataContainer("BigRIPSPPAC");
  std::cout<<ppac_array->GetName()<<std::endl;
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
    
  //PID reconstructed data:
  TClonesArray *rips_array = 
   (TClonesArray *)sman->FindDataContainer("BigRIPSRIPS");
  std::cout<<rips_array->GetName()<<std::endl;     
  tree->Branch(rips_array->GetName(),&rips_array); 
  
  TClonesArray *tof_array  = 
    (TClonesArray *)sman->FindDataContainer("BigRIPSTOF");
  std::cout<<tof_array->GetName()<<std::endl;     
  tree->Branch(tof_array->GetName(),&tof_array);   
  
  TClonesArray *beam_array = 
    (TClonesArray *)sman->FindDataContainer("BigRIPSBeam");	
  std::cout<<beam_array->GetName()<<std::endl;     
  tree->Branch(beam_array->GetName(),&beam_array); 

  int neve = 0;
  while(estore->GetNextEvent()&& neve<1000){
    if(neve%1000==0)
      std::cout << "event: " << neve << std::endl;
      
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //Making the BigRIPS tree calibration
    brcalib->ClearData();
    brcalib->ReconstructData();
    //Reconstructiong the PID
    recopid->ClearData();
    recopid->ReconstructData();
    
    tree->Fill();
    neve++;
  }
  cout<<"Writing the tree."<<endl;
  
  fout->Write();
  fout->Close();
}
