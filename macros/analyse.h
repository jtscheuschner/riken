#ifndef analyse_h
#define analyse_h


#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

//global variables start with capital letter
//counting numbers start with "i"

class analyse{

 private:
  UInt_t maxEvents;

 public:
  //main function
  void main(short end);

  //reading and initialize variables and trees
  void ReadPIDCuts(short run);
  void LoadTree(short run);

  //function to analyse data
  void LaBr();
  void DALI();

  //function to cut
  Bool_t PID();
  //reset and delete all variables/tree
  void reset();
  void set(short run);//set all variables extern

  //needed global variables
  TFile* output;
  //cutvariables
  //PID
  Float_t AOQ_ms[2][4]; //mean and sigma
  Float_t ZET_ms[2][4];
  Float_t Sigma;
  //Histograms

  //DALI
  TH2F* h2_DALI_E;
  TH2F* h2_DALI_EvsT[180];
  //LaBr
  TH2F* h2_LaBr_E[3];
  TH2F* h2_LaBr_EvsT[3][8];
  //TKE
  //PID
  TH2F* h2_PID_before_0;
  TH2F* h2_PID_after_0;
  TH2F* h2_PID_before_2;
  TH2F* h2_PID_after_2;
  //treevariables
  TChain* fChain;
  Long64_t nb;
   //leaf types
  const static Int_t kMaxDALINaI = 81;
   Int_t           DALINaI_;
   UInt_t          DALINaI_fUniqueID[kMaxDALINaI];   //[DALINaI_]
   UInt_t          DALINaI_fBits[kMaxDALINaI];   //[DALINaI_]
   Int_t           DALINaI_id[kMaxDALINaI];   //[DALINaI_]
   Int_t           DALINaI_fpl[kMaxDALINaI];   //[DALINaI_]
   TString         DALINaI_name[kMaxDALINaI];
   Int_t           DALINaI_fDataState[kMaxDALINaI];   //[DALINaI_]
   Int_t           DALINaI_fADC[kMaxDALINaI];   //[DALINaI_]
   Int_t           DALINaI_fTDC[kMaxDALINaI];   //[DALINaI_]
   Int_t           DALINaI_layer[kMaxDALINaI];   //[DALINaI_]
   Double_t        DALINaI_theta[kMaxDALINaI];   //[DALINaI_]
   Double_t        DALINaI_fXPos[kMaxDALINaI];   //[DALINaI_]
   Double_t        DALINaI_fYPos[kMaxDALINaI];   //[DALINaI_]
   Double_t        DALINaI_fZPos[kMaxDALINaI];   //[DALINaI_]
   Double_t        DALINaI_costheta[kMaxDALINaI];   //[DALINaI_]
   Double_t        DALINaI_fEnergy[kMaxDALINaI];   //[DALINaI_]
   Double_t        DALINaI_fDoppCorEnergy[kMaxDALINaI];   //[DALINaI_]
   Double_t        DALINaI_fEnergyWithoutT[kMaxDALINaI];   //[DALINaI_]
   Double_t        DALINaI_fTime[kMaxDALINaI];   //[DALINaI_]
   Double_t        DALINaI_fTimeOffseted[kMaxDALINaI];   //[DALINaI_]
   Double_t        DALINaI_fTimeTrueEnergy[kMaxDALINaI];   //[DALINaI_]
   Double_t        DALINaI_fTimeTrueDoppCorEnergy[kMaxDALINaI];   //[DALINaI_]
   Double_t        DALINaI_fTimeTrueDoppCorVertexEnergy[kMaxDALINaI];   //[DALINaI_]
   Double_t        DALINaI_fTimeTrueTime[kMaxDALINaI];   //[DALINaI_]
   Double_t        DALINaI_fTimeTrueTimeOffseted[kMaxDALINaI];   //[DALINaI_]
   Int_t           triggerbit;
   Int_t           neve;
   /*
   Double_t        tgtx;
   Double_t        tgty;
   Double_t        tgta;
   Double_t        tgtb;
   Double_t        F10X;
   Double_t        F10Y;
   Double_t        F10A;
   Double_t        F10B;
   Double_t        F11X;
   Double_t        F11Y;
   Double_t        F11A;
   Double_t        F11B;
   Double_t        F3PLA_QL;
   Double_t        F3PLA_QR;
   Double_t        F7PLA_QL;
   Double_t        F7PLA_QR;
   Double_t        F8PLA_QL;
   Double_t        F8PLA_QR;
   Double_t        F11PLA_QL;
   Double_t        F11PLA_QR;
   */
   Double_t        DELTA[3];
   Double_t        ANGLE[3];
   Double_t        BRHO[3];
   Double_t        TOF[3];
   Double_t        BETA[3];
   Double_t        AOQ[5];
   Double_t        ZET[5];
   Double_t        TKE;
   /*
   Double_t        ASUM;
   Int_t           dalimultwotime;
   Int_t           dalimult;
   Int_t           dalitimetruemult;
   Int_t           dalimultthres;
   Int_t           dalitimetruemultthres;
   */
   Double_t        LaBr_GAIN[3][8];
   Double_t        LaBr_GAINC[3][8];
   /*
   Double_t        MIDGAIN[8];
   Double_t        MIDGAINC[8];
   Double_t        HIGAIN[8];
   Double_t        HIGAINC[8];
   */
   // List of branches
   TBranch        *b_DALINaI_;   //!
   TBranch        *b_DALINaI_fUniqueID;   //!
   TBranch        *b_DALINaI_fBits;   //!
   TBranch        *b_DALINaI_id;   //!
   TBranch        *b_DALINaI_fpl;   //!
   TBranch        *b_DALINaI_name;   //!
   TBranch        *b_DALINaI_fDataState;   //!
   TBranch        *b_DALINaI_fADC;   //!
   TBranch        *b_DALINaI_fTDC;   //!
   TBranch        *b_DALINaI_layer;   //!
   TBranch        *b_DALINaI_theta;   //!
   TBranch        *b_DALINaI_fXPos;   //!
   TBranch        *b_DALINaI_fYPos;   //!
   TBranch        *b_DALINaI_fZPos;   //!
   TBranch        *b_DALINaI_costheta;   //!
   TBranch        *b_DALINaI_fEnergy;   //!
   TBranch        *b_DALINaI_fDoppCorEnergy;   //!
   TBranch        *b_DALINaI_fEnergyWithoutT;   //!
   TBranch        *b_DALINaI_fTime;   //!
   TBranch        *b_DALINaI_fTimeOffseted;   //!
   TBranch        *b_DALINaI_fTimeTrueEnergy;   //!
   TBranch        *b_DALINaI_fTimeTrueDoppCorEnergy;   //!
   TBranch        *b_DALINaI_fTimeTrueDoppCorVertexEnergy;   //!
   TBranch        *b_DALINaI_fTimeTrueTime;   //!
   TBranch        *b_DALINaI_fTimeTrueTimeOffseted;   //!
   TBranch        *b_triggerbit;   //!
   TBranch        *b_neve;   //!
   /*
   TBranch        *b_tgtx;   //!
   TBranch        *b_tgty;   //!
   TBranch        *b_tgta;   //!
   TBranch        *b_tgtb;   //!
   TBranch        *b_F10X;   //!
   TBranch        *b_F10Y;   //!
   TBranch        *b_F10A;   //!
   TBranch        *b_F10B;   //!
   TBranch        *b_F11X;   //!
   TBranch        *b_F11Y;   //!
   TBranch        *b_F11A;   //!
   TBranch        *b_F11B;   //!
   TBranch        *b_F3PLA_QL;   //!
   TBranch        *b_F3PLA_QR;   //!
   TBranch        *b_F7PLA_QL;   //!
   TBranch        *b_F7PLA_QR;   //!
   TBranch        *b_F8PLA_QL;   //!
   TBranch        *b_F8PLA_QR;   //!
   TBranch        *b_F11PLA_QL;   //!
   TBranch        *b_F11PLA_QR;   //!
   */
   TBranch        *b_DELTA;   //!
   TBranch        *b_ANGLE;   //!
   TBranch        *b_BRHO;   //!
   TBranch        *b_TOF;   //!
   TBranch        *b_BETA;   //!
   TBranch        *b_AOQ;   //!
   TBranch        *b_ZET;   //!
   TBranch        *b_TKE;   //!
   /*
   TBranch        *b_ASUM;   //!
   TBranch        *b_dalimultwotime;   //!
   TBranch        *b_dalimult;   //!
   TBranch        *b_dalitimetruemult;   //!
   TBranch        *b_dalimultthres;   //!
   TBranch        *b_dalitimetruemultthres;   //!
   */
   TBranch        *b_LOWGAIN;   //!
   TBranch        *b_LOWGAINC;   //!
   TBranch        *b_MIDGAIN;   //!
   TBranch        *b_MIDGAINC;   //!
   TBranch        *b_HIGAIN;   //!
   TBranch        *b_HIGAINC;   //!


};//class

#endif
