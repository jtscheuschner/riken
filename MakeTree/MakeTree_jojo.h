#ifndef MakeTree_jojo_h
#define MakeTree_jojo_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include "TArtStoreManager.hh"
#include "TArtEventStore.hh"
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


class MakeTree {


public :

  const static Int_t fpl_ch = 8;
  const static Int_t fn_labr = 8;
  const static Int_t E_offset_br = -14;
   
  //Define function
  void main(short run, Bool_t test);
  void Reset();
  void Init(short run);
  void Reconstruct(short run);
  void Plastic();
  void GetVar();
  void IC_F11();
  void PID();
  void PID_out3();
  void Gamma();
  void LaBr();
  void DALI();
  void LoadPara();
  void Array();
  void Read(short run);
  void ReadTimeoff(short run);
  void ReadPID(short run);
  void ReadIC(short run);
  void ReadQTC();
  void ReadDALI();
  void ReadNakatsuka();
  void AddTree();
  void Trigger();
  //  void IC(Bool_t br);
  void GetIC(Bool_t Br);
  void stop_interrupt();
  //Float_t ESq_offset(short run);
  Bool_t Gated(Float_t stamp, short fpl);
  Double_t QTC_Time(Double_t e_mid, Double_t e_low, UInt_t crystal);
  void Calibration(UInt_t run);

  //global variables
  Float_t ftof_F8_T;
  Float_t ftofoff[6];
  Float_t fQTCTimeoffset[fn_labr];
  Float_t corr_a[6];
  Float_t corr_y[6];
  Float_t corr_x[6];
  Float_t fESqoff;
  Float_t fVar[3][6];
  Float_t fVarout3[10];
  Float_t fCorr_aoq[3][2][6];
  Float_t fCorr_aoqout3[2][10];
  Float_t fQTC_slope[3][fn_labr];
  Float_t fQTC_offset[3][fn_labr];
  Float_t fADCoffset_t[150];
  Float_t fQTC_slopecorr[fn_labr][5];
  Float_t fQTC_slopecorr_low[fn_labr][4];
  Float_t ESQ;
  Int_t fQTC[3][fn_labr];
  TTree* fTree;

  Float_t am[8][5];
  Float_t bm[8][7];
  Float_t al[8][5];
  Float_t bl[8][7];

  Float_t ftime;
  Float_t ftheta_labr;
  Float_t clight;

  //additional Tree variables
  Double_t AOQ[4];
  Double_t ZET[5];
  Double_t FX[6];
  Double_t FY[6];
  Double_t FA[6];
  Double_t FB[6];
  Double_t delta[6];
  Double_t F7T;
  Double_t F8T;
  Double_t BETA[4];
  short gate[4];
  Int_t TRIGGER;        
  Int_t fplMulti[fpl_ch];
  Float_t  HIGAIN[fn_labr];
  Float_t MIDGAIN[fn_labr];
  Float_t LOWGAIN[fn_labr];
  Float_t QTCTime[fn_labr];
  Float_t  HIGAINC[fn_labr];
  Float_t MIDGAINC[fn_labr];
  Float_t LOWGAINC[fn_labr];
  Float_t QTCTimeM[fn_labr];
  Float_t QTCTimeL[fn_labr];
  Int_t DALIfID[150];		      
  Int_t DALIfLayer[150];	      
  Double_t DALIfx[150];		      
  Double_t DALIfy[150];		      
  Double_t DALIfz[150];		      
  Double_t DALIfd[150];		      
  Double_t DALIfTheta[150];	      
  Double_t DALIfTOF[150];	      
  Double_t DALIfRawEnergy[150];	      
  Double_t DALIfRawTDC[150];	      
  Double_t DALIfEnergy[150];           
  Double_t DALIfDoppCorEnergy[150];    
  Double_t DALIfTOFGEnergy[150];       
  Double_t DALIfTOFGDoppCorEnergy[150];

  Double_t F3Pla_TR[fpl_ch]; 
  Double_t F3Pla_TL[fpl_ch]; 
  Double_t F7Pla_TR[fpl_ch]; 
  Double_t F7Pla_TL[fpl_ch]; 
  Double_t F8Pla_TR[fpl_ch]; 
  Double_t F8Pla_TL[fpl_ch]; 
  Double_t F11Pla_TR[fpl_ch];
  Double_t F11Pla_TL[fpl_ch];

  Double_t TKE;
  Double_t ASUM;
  //Anaroot stuff
  bool stoploop;
  TArtStoreManager * fman;
  TArtEventStore *fstore;
  TArtBigRIPSParameters *fpara;
  TArtDALIParameters *fdpara;
  TArtRawEventObject *frawevent;
  TArtRecoPID *frecopid;

  TArtCalibDALI* fdalicalib;
  TArtCalibPID *fpid;
  TArtCalibPPAC *fppac;
  TArtPlastic *fplastic;
  TArtCalibPlastic *fpla;
  TArtCalibFocalPlane *ffpl;

  TArtRIPS *frips3to5;
  TArtRIPS *frips5to7;
  TArtRIPS *frips8to10; 
  TArtRIPS *frips8to9;
  TArtRIPS *frips10to11;

  TArtTOF *ftof3to7;
  TArtTOF *ftof3to7_2;
  TArtTOF *ftof3to8;
  TArtTOF *ftof8to11[3];

  TArtIC* fIC;
  TArtIC* fIC2;

  TArtBeam *fbeam_br;
  TArtBeam *fbeam_br2;
  TArtBeam *fbeam_br_38;
  TArtBeam *fbeam_zd_810;
  TArtBeam *fbeam_zd_1011;
  TArtBeam *fbeam_zd_1011_test;

  TClonesArray* fppac_array;
  TClonesArray* fpla_array;
  TClonesArray* fic_array;
  TClonesArray* ffpl_array;
  TClonesArray* fdali_array;
  TClonesArray* frips_array;
  TClonesArray* ftof_array;
  TClonesArray* fbeam_array;
  TClonesArray *rips_array;
  TClonesArray *tof_array;
  TClonesArray *beam_array;

  TFile* fout;

  MakeTree();
  virtual ~MakeTree();
};
#endif 

MakeTree::MakeTree() : fTree(0){
}
MakeTree::~MakeTree(){

}

