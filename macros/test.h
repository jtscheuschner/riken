//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed May  6 09:46:33 2015 by ROOT version 5.34/25
// from TTree tree/tree
// found on file: data/rootfiles/test/run0102_301.388_-156.264.root
//////////////////////////////////////////////////////////

#ifndef test_h
#define test_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2F.h>

// Header file for the classes stored in the TTree if any.
#include "/home/jtscheuschner/.local/root/core/cont/inc/TClonesArray.h"
#include "/home/jtscheuschner/.local/root/core/base/inc/TObject.h"
#include "/home/jtscheuschner/.local/anaroot/include/TArtDataObject.hh"
#include "/home/jtscheuschner/.local/anaroot/include/TArtRIPS.hh"
#include "/home/jtscheuschner/.local/anaroot/include/TArtTOF.hh"
#include "/home/jtscheuschner/.local/anaroot/include/TArtBeam.hh"
#include "/home/jtscheuschner/.local/root/core/base/inc/TNamed.h"
#include "/home/jtscheuschner/.local/anaroot/include/TArtEventInfo.hh"
#include "/home/jtscheuschner/.local/anaroot/include/TArtPPAC.hh"
#include "/home/jtscheuschner/.local/anaroot/include/TArtPlastic.hh"
#include "/home/jtscheuschner/.local/anaroot/include/TArtIC.hh"
#include "/home/jtscheuschner/.local/anaroot/include/TArtFocalPlane.hh"


// Fixed size dimensions of array or collections stored in the TTree if any.
   const Int_t kMaxBigRIPSPPAC = 36;
   const Int_t kMaxBigRIPSPlastic = 4;
   const Int_t kMaxBigRIPSIC = 2;
   const Int_t kMaxBigRIPSFocalPlane = 14;
   const Int_t kMaxDALINaI = 76;
   const Int_t kMaxBigRIPSRIPS = 5;
   const Int_t kMaxBigRIPSTOF = 2;
   const Int_t kMaxBigRIPSBeam = 7;
   const Int_t variables = 9;

//functions
//  void Fitter(TH2F* h2, Int_t var, Bool_t clear);

//global variables
  Double_t correction[2][variables];

class test {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           BigRIPSPPAC_;
   UInt_t          BigRIPSPPAC_fUniqueID[kMaxBigRIPSPPAC];   //[BigRIPSPPAC_]
   UInt_t          BigRIPSPPAC_fBits[kMaxBigRIPSPPAC];   //[BigRIPSPPAC_]
   Int_t           BigRIPSPPAC_id[kMaxBigRIPSPPAC];   //[BigRIPSPPAC_]
   Int_t           BigRIPSPPAC_fpl[kMaxBigRIPSPPAC];   //[BigRIPSPPAC_]
   TString         BigRIPSPPAC_name[kMaxBigRIPSPPAC];
   Int_t           BigRIPSPPAC_fDataState[kMaxBigRIPSPPAC];   //[BigRIPSPPAC_]
   Double_t        BigRIPSPPAC_xzpos[kMaxBigRIPSPPAC];   //[BigRIPSPPAC_]
   Double_t        BigRIPSPPAC_yzpos[kMaxBigRIPSPPAC];   //[BigRIPSPPAC_]
   Int_t           BigRIPSPPAC_fTX1Raw[kMaxBigRIPSPPAC];   //[BigRIPSPPAC_]
   Int_t           BigRIPSPPAC_fTX2Raw[kMaxBigRIPSPPAC];   //[BigRIPSPPAC_]
   Int_t           BigRIPSPPAC_fTY1Raw[kMaxBigRIPSPPAC];   //[BigRIPSPPAC_]
   Int_t           BigRIPSPPAC_fTY2Raw[kMaxBigRIPSPPAC];   //[BigRIPSPPAC_]
   Int_t           BigRIPSPPAC_fTARaw[kMaxBigRIPSPPAC];   //[BigRIPSPPAC_]
   Int_t           BigRIPSPPAC_fQX1Raw[kMaxBigRIPSPPAC];   //[BigRIPSPPAC_]
   Int_t           BigRIPSPPAC_fQX2Raw[kMaxBigRIPSPPAC];   //[BigRIPSPPAC_]
   Int_t           BigRIPSPPAC_fQY1Raw[kMaxBigRIPSPPAC];   //[BigRIPSPPAC_]
   Int_t           BigRIPSPPAC_fQY2Raw[kMaxBigRIPSPPAC];   //[BigRIPSPPAC_]
   Int_t           BigRIPSPPAC_fQARaw[kMaxBigRIPSPPAC];   //[BigRIPSPPAC_]
   Double_t        BigRIPSPPAC_fTX1[kMaxBigRIPSPPAC];   //[BigRIPSPPAC_]
   Double_t        BigRIPSPPAC_fTX2[kMaxBigRIPSPPAC];   //[BigRIPSPPAC_]
   Double_t        BigRIPSPPAC_fTY1[kMaxBigRIPSPPAC];   //[BigRIPSPPAC_]
   Double_t        BigRIPSPPAC_fTY2[kMaxBigRIPSPPAC];   //[BigRIPSPPAC_]
   Double_t        BigRIPSPPAC_fTA[kMaxBigRIPSPPAC];   //[BigRIPSPPAC_]
   Double_t        BigRIPSPPAC_fTSumX[kMaxBigRIPSPPAC];   //[BigRIPSPPAC_]
   Double_t        BigRIPSPPAC_fTSumY[kMaxBigRIPSPPAC];   //[BigRIPSPPAC_]
   Double_t        BigRIPSPPAC_fTDiffX[kMaxBigRIPSPPAC];   //[BigRIPSPPAC_]
   Double_t        BigRIPSPPAC_fTDiffY[kMaxBigRIPSPPAC];   //[BigRIPSPPAC_]
   Double_t        BigRIPSPPAC_fX[kMaxBigRIPSPPAC];   //[BigRIPSPPAC_]
   Double_t        BigRIPSPPAC_fY[kMaxBigRIPSPPAC];   //[BigRIPSPPAC_]
   Bool_t          BigRIPSPPAC_fFiredX[kMaxBigRIPSPPAC];   //[BigRIPSPPAC_]
   Bool_t          BigRIPSPPAC_fFiredY[kMaxBigRIPSPPAC];   //[BigRIPSPPAC_]
   Int_t           BigRIPSPlastic_;
   UInt_t          BigRIPSPlastic_fUniqueID[kMaxBigRIPSPlastic];   //[BigRIPSPlastic_]
   UInt_t          BigRIPSPlastic_fBits[kMaxBigRIPSPlastic];   //[BigRIPSPlastic_]
   Int_t           BigRIPSPlastic_id[kMaxBigRIPSPlastic];   //[BigRIPSPlastic_]
   Int_t           BigRIPSPlastic_fpl[kMaxBigRIPSPlastic];   //[BigRIPSPlastic_]
   TString         BigRIPSPlastic_name[kMaxBigRIPSPlastic];
   Int_t           BigRIPSPlastic_fDataState[kMaxBigRIPSPlastic];   //[BigRIPSPlastic_]
   Double_t        BigRIPSPlastic_zposition[kMaxBigRIPSPlastic];   //[BigRIPSPlastic_]
   Double_t        BigRIPSPlastic_zoffset[kMaxBigRIPSPlastic];   //[BigRIPSPlastic_]
   Int_t           BigRIPSPlastic_fTLRaw[kMaxBigRIPSPlastic];   //[BigRIPSPlastic_]
   Int_t           BigRIPSPlastic_fTRRaw[kMaxBigRIPSPlastic];   //[BigRIPSPlastic_]
   Int_t           BigRIPSPlastic_fQLRaw[kMaxBigRIPSPlastic];   //[BigRIPSPlastic_]
   Int_t           BigRIPSPlastic_fQRRaw[kMaxBigRIPSPlastic];   //[BigRIPSPlastic_]
   Double_t        BigRIPSPlastic_fTime[kMaxBigRIPSPlastic];   //[BigRIPSPlastic_]
   Double_t        BigRIPSPlastic_fTimeL[kMaxBigRIPSPlastic];   //[BigRIPSPlastic_]
   Double_t        BigRIPSPlastic_fTimeR[kMaxBigRIPSPlastic];   //[BigRIPSPlastic_]
   Double_t        BigRIPSPlastic_fTimeLSlew[kMaxBigRIPSPlastic];   //[BigRIPSPlastic_]
   Double_t        BigRIPSPlastic_fTimeRSlew[kMaxBigRIPSPlastic];   //[BigRIPSPlastic_]
   Double_t        BigRIPSPlastic_fTimeSlew[kMaxBigRIPSPlastic];   //[BigRIPSPlastic_]
   Int_t           BigRIPSIC_;
   UInt_t          BigRIPSIC_fUniqueID[kMaxBigRIPSIC];   //[BigRIPSIC_]
   UInt_t          BigRIPSIC_fBits[kMaxBigRIPSIC];   //[BigRIPSIC_]
   Int_t           BigRIPSIC_id[kMaxBigRIPSIC];   //[BigRIPSIC_]
   Int_t           BigRIPSIC_fpl[kMaxBigRIPSIC];   //[BigRIPSIC_]
   TString         BigRIPSIC_name[kMaxBigRIPSIC];
   Int_t           BigRIPSIC_fDataState[kMaxBigRIPSIC];   //[BigRIPSIC_]
   Double_t        BigRIPSIC_zcoef[kMaxBigRIPSIC][2];   //[BigRIPSIC_]
   Double_t        BigRIPSIC_ionpair[kMaxBigRIPSIC];   //[BigRIPSIC_]
   Int_t           BigRIPSIC_nhitchannel[kMaxBigRIPSIC];   //[BigRIPSIC_]
   Int_t           BigRIPSIC_fADC[kMaxBigRIPSIC][12];   //[BigRIPSIC_]
   Double_t        BigRIPSIC_fEnergy[kMaxBigRIPSIC][12];   //[BigRIPSIC_]
   Double_t        BigRIPSIC_fRawADCSqSum[kMaxBigRIPSIC];   //[BigRIPSIC_]
   Double_t        BigRIPSIC_fRawADCAvSum[kMaxBigRIPSIC];   //[BigRIPSIC_]
   Double_t        BigRIPSIC_fCalMeVSqSum[kMaxBigRIPSIC];   //[BigRIPSIC_]
   Double_t        BigRIPSIC_fCalMeVAvSum[kMaxBigRIPSIC];   //[BigRIPSIC_]
   Int_t           BigRIPSFocalPlane_;
   UInt_t          BigRIPSFocalPlane_fUniqueID[kMaxBigRIPSFocalPlane];   //[BigRIPSFocalPlane_]
   UInt_t          BigRIPSFocalPlane_fBits[kMaxBigRIPSFocalPlane];   //[BigRIPSFocalPlane_]
   Int_t           BigRIPSFocalPlane_id[kMaxBigRIPSFocalPlane];   //[BigRIPSFocalPlane_]
   Int_t           BigRIPSFocalPlane_fpl[kMaxBigRIPSFocalPlane];   //[BigRIPSFocalPlane_]
   TString         BigRIPSFocalPlane_name[kMaxBigRIPSFocalPlane];
   Int_t           BigRIPSFocalPlane_fDataState[kMaxBigRIPSFocalPlane];   //[BigRIPSFocalPlane_]
   TVectorT<double> BigRIPSFocalPlane_opt_vector[kMaxBigRIPSFocalPlane];
   Double_t        BigRIPSFocalPlane_X[kMaxBigRIPSFocalPlane];   //[BigRIPSFocalPlane_]
   Double_t        BigRIPSFocalPlane_A[kMaxBigRIPSFocalPlane];   //[BigRIPSFocalPlane_]
   Double_t        BigRIPSFocalPlane_Y[kMaxBigRIPSFocalPlane];   //[BigRIPSFocalPlane_]
   Double_t        BigRIPSFocalPlane_B[kMaxBigRIPSFocalPlane];   //[BigRIPSFocalPlane_]
   Int_t           BigRIPSFocalPlane_nfired_ppacx[kMaxBigRIPSFocalPlane];   //[BigRIPSFocalPlane_]
   Int_t           BigRIPSFocalPlane_nfired_ppacy[kMaxBigRIPSFocalPlane];   //[BigRIPSFocalPlane_]
   Double_t        BigRIPSFocalPlane_zpos[kMaxBigRIPSFocalPlane];   //[BigRIPSFocalPlane_]
   Double_t        BigRIPSFocalPlane_zpos_offset[kMaxBigRIPSFocalPlane];   //[BigRIPSFocalPlane_]
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
   Int_t           BigRIPSRIPS_;
   UInt_t          BigRIPSRIPS_fUniqueID[kMaxBigRIPSRIPS];   //[BigRIPSRIPS_]
   UInt_t          BigRIPSRIPS_fBits[kMaxBigRIPSRIPS];   //[BigRIPSRIPS_]
   Int_t           BigRIPSRIPS_id[kMaxBigRIPSRIPS];   //[BigRIPSRIPS_]
   Int_t           BigRIPSRIPS_fpl[kMaxBigRIPSRIPS];   //[BigRIPSRIPS_]
   TString         BigRIPSRIPS_name[kMaxBigRIPSRIPS];
   Int_t           BigRIPSRIPS_fDataState[kMaxBigRIPSRIPS];   //[BigRIPSRIPS_]
   Int_t           BigRIPSRIPS_upstream_fpl[kMaxBigRIPSRIPS];   //[BigRIPSRIPS_]
   Int_t           BigRIPSRIPS_downstream_fpl[kMaxBigRIPSRIPS];   //[BigRIPSRIPS_]
   Double_t        BigRIPSRIPS_center_brho[kMaxBigRIPSRIPS];   //[BigRIPSRIPS_]
   Double_t        BigRIPSRIPS_brho[kMaxBigRIPSRIPS];   //[BigRIPSRIPS_]
   Double_t        BigRIPSRIPS_length[kMaxBigRIPSRIPS];   //[BigRIPSRIPS_]
   TMatrixT<double> BigRIPSRIPS_matrix[kMaxBigRIPSRIPS];
   Double_t        BigRIPSRIPS_delta[kMaxBigRIPSRIPS];   //[BigRIPSRIPS_]
   Double_t        BigRIPSRIPS_angle[kMaxBigRIPSRIPS];   //[BigRIPSRIPS_]
   TString         BigRIPSRIPS_dipolename[kMaxBigRIPSRIPS];
   Int_t           BigRIPSTOF_;
   UInt_t          BigRIPSTOF_fUniqueID[kMaxBigRIPSTOF];   //[BigRIPSTOF_]
   UInt_t          BigRIPSTOF_fBits[kMaxBigRIPSTOF];   //[BigRIPSTOF_]
   Int_t           BigRIPSTOF_id[kMaxBigRIPSTOF];   //[BigRIPSTOF_]
   Int_t           BigRIPSTOF_fpl[kMaxBigRIPSTOF];   //[BigRIPSTOF_]
   TString         BigRIPSTOF_name[kMaxBigRIPSTOF];
   Int_t           BigRIPSTOF_fDataState[kMaxBigRIPSTOF];   //[BigRIPSTOF_]
   Double_t        BigRIPSTOF_tof[kMaxBigRIPSTOF];   //[BigRIPSTOF_]
   Double_t        BigRIPSTOF_clight[kMaxBigRIPSTOF];   //[BigRIPSTOF_]
   Double_t        BigRIPSTOF_length[kMaxBigRIPSTOF];   //[BigRIPSTOF_]
   Double_t        BigRIPSTOF_ulength[kMaxBigRIPSTOF];   //[BigRIPSTOF_]
   Double_t        BigRIPSTOF_dlength[kMaxBigRIPSTOF];   //[BigRIPSTOF_]
   TString         BigRIPSTOF_upstream_plname[kMaxBigRIPSTOF];
   TString         BigRIPSTOF_downstream_plname[kMaxBigRIPSTOF];
   Int_t           BigRIPSTOF_upstream_plfpl[kMaxBigRIPSTOF];   //[BigRIPSTOF_]
   Int_t           BigRIPSTOF_downstream_plfpl[kMaxBigRIPSTOF];   //[BigRIPSTOF_]
   Double_t        BigRIPSTOF_time_offset[kMaxBigRIPSTOF];   //[BigRIPSTOF_]
   Int_t           BigRIPSBeam_;
   UInt_t          BigRIPSBeam_fUniqueID[kMaxBigRIPSBeam];   //[BigRIPSBeam_]
   UInt_t          BigRIPSBeam_fBits[kMaxBigRIPSBeam];   //[BigRIPSBeam_]
   Int_t           BigRIPSBeam_id[kMaxBigRIPSBeam];   //[BigRIPSBeam_]
   Int_t           BigRIPSBeam_fpl[kMaxBigRIPSBeam];   //[BigRIPSBeam_]
   TString         BigRIPSBeam_name[kMaxBigRIPSBeam];
   Int_t           BigRIPSBeam_fDataState[kMaxBigRIPSBeam];   //[BigRIPSBeam_]
   Double_t        BigRIPSBeam_brho[kMaxBigRIPSBeam];   //[BigRIPSBeam_]
   Double_t        BigRIPSBeam_aoq[kMaxBigRIPSBeam];   //[BigRIPSBeam_]
   Double_t        BigRIPSBeam_zet[kMaxBigRIPSBeam];   //[BigRIPSBeam_]
   Double_t        BigRIPSBeam_beta[kMaxBigRIPSBeam];   //[BigRIPSBeam_]
   Double_t        BigRIPSBeam_clight[kMaxBigRIPSBeam];   //[BigRIPSBeam_]
   Double_t        BigRIPSBeam_mnucleon[kMaxBigRIPSBeam];   //[BigRIPSBeam_]
   Int_t           BigRIPSBeam_nrips[kMaxBigRIPSBeam];   //[BigRIPSBeam_]
   TString         BigRIPSBeam_ripsname[2][kMaxBigRIPSBeam];
   TString         BigRIPSBeam_tofname[kMaxBigRIPSBeam];
   TString         BigRIPSBeam_icname[kMaxBigRIPSBeam];
   Int_t           triggerbit;
   Int_t           neve;
   Double_t        tgtx;
   Double_t        tgty;
   Double_t        tgta;
   Double_t        tgtb;
   Double_t        FX[6];
   Double_t        FY[6];
   Double_t        FA[6];
   Double_t        FB[6];
   Double_t        F8T;
   Double_t        delta[6];
   Double_t        F10B;
   Double_t        F11B;
   Double_t        F3PLA_QL;
   Double_t        F3PLA_QR;
   Double_t        F7PLA_QL;
   Double_t        F7PLA_QR;
   Double_t        F8PLA_QL;
   Double_t        F8PLA_QR;
   Double_t        F11PLA_QL;
   Double_t        F11PLA_QR;
   Double_t        DELTA[3];
   Double_t        ANGLE[3];
   Double_t        BRHO[3];
   Double_t        TOF[3];
   Double_t        BETA[4];
   Double_t        AOQ[7];
   Double_t        ZET[7];
   Double_t        ESq;
   Double_t        TKE;
   Double_t        ASUM;
   Int_t           dalimult;
   Int_t           dalitimetruemult;
   Double_t        LOWGAIN[8];
   Double_t        LOWGAINC[8];
   Double_t        MIDGAIN[8];
   Double_t        MIDGAINC[8];
   Double_t        MIDGAINCT[8];
   Double_t        HIGAIN[8];
   Double_t        HIGAINC[8];
   Double_t        QTCTime[8];
   Double_t        QTCTimeC[8];
   Double_t        MIDGAINCTA;

   // List of branches
   TBranch        *b_BigRIPSPPAC_;   //!
   TBranch        *b_BigRIPSPPAC_fUniqueID;   //!
   TBranch        *b_BigRIPSPPAC_fBits;   //!
   TBranch        *b_BigRIPSPPAC_id;   //!
   TBranch        *b_BigRIPSPPAC_fpl;   //!
   TBranch        *b_BigRIPSPPAC_name;   //!
   TBranch        *b_BigRIPSPPAC_fDataState;   //!
   TBranch        *b_BigRIPSPPAC_xzpos;   //!
   TBranch        *b_BigRIPSPPAC_yzpos;   //!
   TBranch        *b_BigRIPSPPAC_fTX1Raw;   //!
   TBranch        *b_BigRIPSPPAC_fTX2Raw;   //!
   TBranch        *b_BigRIPSPPAC_fTY1Raw;   //!
   TBranch        *b_BigRIPSPPAC_fTY2Raw;   //!
   TBranch        *b_BigRIPSPPAC_fTARaw;   //!
   TBranch        *b_BigRIPSPPAC_fQX1Raw;   //!
   TBranch        *b_BigRIPSPPAC_fQX2Raw;   //!
   TBranch        *b_BigRIPSPPAC_fQY1Raw;   //!
   TBranch        *b_BigRIPSPPAC_fQY2Raw;   //!
   TBranch        *b_BigRIPSPPAC_fQARaw;   //!
   TBranch        *b_BigRIPSPPAC_fTX1;   //!
   TBranch        *b_BigRIPSPPAC_fTX2;   //!
   TBranch        *b_BigRIPSPPAC_fTY1;   //!
   TBranch        *b_BigRIPSPPAC_fTY2;   //!
   TBranch        *b_BigRIPSPPAC_fTA;   //!
   TBranch        *b_BigRIPSPPAC_fTSumX;   //!
   TBranch        *b_BigRIPSPPAC_fTSumY;   //!
   TBranch        *b_BigRIPSPPAC_fTDiffX;   //!
   TBranch        *b_BigRIPSPPAC_fTDiffY;   //!
   TBranch        *b_BigRIPSPPAC_fX;   //!
   TBranch        *b_BigRIPSPPAC_fY;   //!
   TBranch        *b_BigRIPSPPAC_fFiredX;   //!
   TBranch        *b_BigRIPSPPAC_fFiredY;   //!
   TBranch        *b_BigRIPSPlastic_;   //!
   TBranch        *b_BigRIPSPlastic_fUniqueID;   //!
   TBranch        *b_BigRIPSPlastic_fBits;   //!
   TBranch        *b_BigRIPSPlastic_id;   //!
   TBranch        *b_BigRIPSPlastic_fpl;   //!
   TBranch        *b_BigRIPSPlastic_name;   //!
   TBranch        *b_BigRIPSPlastic_fDataState;   //!
   TBranch        *b_BigRIPSPlastic_zposition;   //!
   TBranch        *b_BigRIPSPlastic_zoffset;   //!
   TBranch        *b_BigRIPSPlastic_fTLRaw;   //!
   TBranch        *b_BigRIPSPlastic_fTRRaw;   //!
   TBranch        *b_BigRIPSPlastic_fQLRaw;   //!
   TBranch        *b_BigRIPSPlastic_fQRRaw;   //!
   TBranch        *b_BigRIPSPlastic_fTime;   //!
   TBranch        *b_BigRIPSPlastic_fTimeL;   //!
   TBranch        *b_BigRIPSPlastic_fTimeR;   //!
   TBranch        *b_BigRIPSPlastic_fTimeLSlew;   //!
   TBranch        *b_BigRIPSPlastic_fTimeRSlew;   //!
   TBranch        *b_BigRIPSPlastic_fTimeSlew;   //!
   TBranch        *b_BigRIPSIC_;   //!
   TBranch        *b_BigRIPSIC_fUniqueID;   //!
   TBranch        *b_BigRIPSIC_fBits;   //!
   TBranch        *b_BigRIPSIC_id;   //!
   TBranch        *b_BigRIPSIC_fpl;   //!
   TBranch        *b_BigRIPSIC_name;   //!
   TBranch        *b_BigRIPSIC_fDataState;   //!
   TBranch        *b_BigRIPSIC_zcoef;   //!
   TBranch        *b_BigRIPSIC_ionpair;   //!
   TBranch        *b_BigRIPSIC_nhitchannel;   //!
   TBranch        *b_BigRIPSIC_fADC;   //!
   TBranch        *b_BigRIPSIC_fEnergy;   //!
   TBranch        *b_BigRIPSIC_fRawADCSqSum;   //!
   TBranch        *b_BigRIPSIC_fRawADCAvSum;   //!
   TBranch        *b_BigRIPSIC_fCalMeVSqSum;   //!
   TBranch        *b_BigRIPSIC_fCalMeVAvSum;   //!
   TBranch        *b_BigRIPSFocalPlane_;   //!
   TBranch        *b_BigRIPSFocalPlane_fUniqueID;   //!
   TBranch        *b_BigRIPSFocalPlane_fBits;   //!
   TBranch        *b_BigRIPSFocalPlane_id;   //!
   TBranch        *b_BigRIPSFocalPlane_fpl;   //!
   TBranch        *b_BigRIPSFocalPlane_name;   //!
   TBranch        *b_BigRIPSFocalPlane_fDataState;   //!
   TBranch        *b_BigRIPSFocalPlane_opt_vector;   //!
   TBranch        *b_BigRIPSFocalPlane_X;   //!
   TBranch        *b_BigRIPSFocalPlane_A;   //!
   TBranch        *b_BigRIPSFocalPlane_Y;   //!
   TBranch        *b_BigRIPSFocalPlane_B;   //!
   TBranch        *b_BigRIPSFocalPlane_nfired_ppacx;   //!
   TBranch        *b_BigRIPSFocalPlane_nfired_ppacy;   //!
   TBranch        *b_BigRIPSFocalPlane_zpos;   //!
   TBranch        *b_BigRIPSFocalPlane_zpos_offset;   //!
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
   TBranch        *b_BigRIPSRIPS_;   //!
   TBranch        *b_BigRIPSRIPS_fUniqueID;   //!
   TBranch        *b_BigRIPSRIPS_fBits;   //!
   TBranch        *b_BigRIPSRIPS_id;   //!
   TBranch        *b_BigRIPSRIPS_fpl;   //!
   TBranch        *b_BigRIPSRIPS_name;   //!
   TBranch        *b_BigRIPSRIPS_fDataState;   //!
   TBranch        *b_BigRIPSRIPS_upstream_fpl;   //!
   TBranch        *b_BigRIPSRIPS_downstream_fpl;   //!
   TBranch        *b_BigRIPSRIPS_center_brho;   //!
   TBranch        *b_BigRIPSRIPS_brho;   //!
   TBranch        *b_BigRIPSRIPS_length;   //!
   TBranch        *b_BigRIPSRIPS_matrix;   //!
   TBranch        *b_BigRIPSRIPS_delta;   //!
   TBranch        *b_BigRIPSRIPS_angle;   //!
   TBranch        *b_BigRIPSRIPS_dipolename;   //!
   TBranch        *b_BigRIPSTOF_;   //!
   TBranch        *b_BigRIPSTOF_fUniqueID;   //!
   TBranch        *b_BigRIPSTOF_fBits;   //!
   TBranch        *b_BigRIPSTOF_id;   //!
   TBranch        *b_BigRIPSTOF_fpl;   //!
   TBranch        *b_BigRIPSTOF_name;   //!
   TBranch        *b_BigRIPSTOF_fDataState;   //!
   TBranch        *b_BigRIPSTOF_tof;   //!
   TBranch        *b_BigRIPSTOF_clight;   //!
   TBranch        *b_BigRIPSTOF_length;   //!
   TBranch        *b_BigRIPSTOF_ulength;   //!
   TBranch        *b_BigRIPSTOF_dlength;   //!
   TBranch        *b_BigRIPSTOF_upstream_plname;   //!
   TBranch        *b_BigRIPSTOF_downstream_plname;   //!
   TBranch        *b_BigRIPSTOF_upstream_plfpl;   //!
   TBranch        *b_BigRIPSTOF_downstream_plfpl;   //!
   TBranch        *b_BigRIPSTOF_time_offset;   //!
   TBranch        *b_BigRIPSBeam_;   //!
   TBranch        *b_BigRIPSBeam_fUniqueID;   //!
   TBranch        *b_BigRIPSBeam_fBits;   //!
   TBranch        *b_BigRIPSBeam_id;   //!
   TBranch        *b_BigRIPSBeam_fpl;   //!
   TBranch        *b_BigRIPSBeam_name;   //!
   TBranch        *b_BigRIPSBeam_fDataState;   //!
   TBranch        *b_BigRIPSBeam_brho;   //!
   TBranch        *b_BigRIPSBeam_aoq;   //!
   TBranch        *b_BigRIPSBeam_zet;   //!
   TBranch        *b_BigRIPSBeam_beta;   //!
   TBranch        *b_BigRIPSBeam_clight;   //!
   TBranch        *b_BigRIPSBeam_mnucleon;   //!
   TBranch        *b_BigRIPSBeam_nrips;   //!
   TBranch        *b_BigRIPSBeam_ripsname;   //!
   TBranch        *b_BigRIPSBeam_tofname;   //!
   TBranch        *b_BigRIPSBeam_icname;   //!
   TBranch        *b_triggerbit;   //!
   TBranch        *b_neve;   //!
   TBranch        *b_tgtx;   //!
   TBranch        *b_tgty;   //!
   TBranch        *b_tgta;   //!
   TBranch        *b_tgtb;   //!
   TBranch        *b_FX;   //!
   TBranch        *b_FY;   //!
   TBranch        *b_FA;   //!
   TBranch        *b_FB;   //!
   TBranch        *b_F8T;   //!
   TBranch        *b_delta;   //!
   TBranch        *b_F10B;   //!
   TBranch        *b_F11B;   //!
   TBranch        *b_F3PLA_QL;   //!
   TBranch        *b_F3PLA_QR;   //!
   TBranch        *b_F7PLA_QL;   //!
   TBranch        *b_F7PLA_QR;   //!
   TBranch        *b_F8PLA_QL;   //!
   TBranch        *b_F8PLA_QR;   //!
   TBranch        *b_F11PLA_QL;   //!
   TBranch        *b_F11PLA_QR;   //!
   TBranch        *b_DELTA;   //!
   TBranch        *b_ANGLE;   //!
   TBranch        *b_BRHO;   //!
   TBranch        *b_TOF;   //!
   TBranch        *b_BETA;   //!
   TBranch        *b_AOQ;   //!
   TBranch        *b_ZET;   //!
   TBranch        *b_ESq;   //!
   TBranch        *b_TKE;   //!
   TBranch        *b_ASUM;   //!
   TBranch        *b_dalimult;   //!
   TBranch        *b_dalitimetruemult;   //!
   TBranch        *b_LOWGAIN;   //!
   TBranch        *b_LOWGAINC;   //!
   TBranch        *b_MIDGAIN;   //!
   TBranch        *b_MIDGAINC;   //!
   TBranch        *b_MIDGAINCT;   //!
   TBranch        *b_HIGAIN;   //!
   TBranch        *b_HIGAINC;   //!
   TBranch        *b_QTCtime;   //!
   TBranch        *b_QTCtimeC;   //!
   TBranch        *b_MIDGAINCTA;   //!

   test(TTree *tree=0);
   virtual ~test();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(short aoq, short max);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   void LoadData(TTree *tree,short run);//Load the tree variables
};

#endif

#ifdef test_cxx
test::test(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
/*   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("data/rootfiles/test/run0102_301.388_-156.264.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("data/rootfiles/test/run0102_301.388_-156.264.root");
      }
      f->GetObject("tree",tree);

   }
   Init(tree);*/
}

test::~test()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t test::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t test::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}


void test::LoadData(TTree* tree, short run){
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
     TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(Form("data/rootfiles/test/run0102_301.388_-156.264.root",run));
      if (!f || !f->IsOpen()) {
	f = new TFile(Form("data/rootfiles/test/run0102_301.388_-156.264.root",run));
      }
      f->GetObject("tree",tree);

   }
   Init(tree);
   
}


void test::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("BigRIPSPPAC", &BigRIPSPPAC_, &b_BigRIPSPPAC_);
   fChain->SetBranchAddress("BigRIPSPPAC.fUniqueID", BigRIPSPPAC_fUniqueID, &b_BigRIPSPPAC_fUniqueID);
   fChain->SetBranchAddress("BigRIPSPPAC.fBits", BigRIPSPPAC_fBits, &b_BigRIPSPPAC_fBits);
   fChain->SetBranchAddress("BigRIPSPPAC.id", BigRIPSPPAC_id, &b_BigRIPSPPAC_id);
   fChain->SetBranchAddress("BigRIPSPPAC.fpl", BigRIPSPPAC_fpl, &b_BigRIPSPPAC_fpl);
   fChain->SetBranchAddress("BigRIPSPPAC.name", BigRIPSPPAC_name, &b_BigRIPSPPAC_name);
   fChain->SetBranchAddress("BigRIPSPPAC.fDataState", BigRIPSPPAC_fDataState, &b_BigRIPSPPAC_fDataState);
   fChain->SetBranchAddress("BigRIPSPPAC.xzpos", BigRIPSPPAC_xzpos, &b_BigRIPSPPAC_xzpos);
   fChain->SetBranchAddress("BigRIPSPPAC.yzpos", BigRIPSPPAC_yzpos, &b_BigRIPSPPAC_yzpos);
   fChain->SetBranchAddress("BigRIPSPPAC.fTX1Raw", BigRIPSPPAC_fTX1Raw, &b_BigRIPSPPAC_fTX1Raw);
   fChain->SetBranchAddress("BigRIPSPPAC.fTX2Raw", BigRIPSPPAC_fTX2Raw, &b_BigRIPSPPAC_fTX2Raw);
   fChain->SetBranchAddress("BigRIPSPPAC.fTY1Raw", BigRIPSPPAC_fTY1Raw, &b_BigRIPSPPAC_fTY1Raw);
   fChain->SetBranchAddress("BigRIPSPPAC.fTY2Raw", BigRIPSPPAC_fTY2Raw, &b_BigRIPSPPAC_fTY2Raw);
   fChain->SetBranchAddress("BigRIPSPPAC.fTARaw", BigRIPSPPAC_fTARaw, &b_BigRIPSPPAC_fTARaw);
   fChain->SetBranchAddress("BigRIPSPPAC.fQX1Raw", BigRIPSPPAC_fQX1Raw, &b_BigRIPSPPAC_fQX1Raw);
   fChain->SetBranchAddress("BigRIPSPPAC.fQX2Raw", BigRIPSPPAC_fQX2Raw, &b_BigRIPSPPAC_fQX2Raw);
   fChain->SetBranchAddress("BigRIPSPPAC.fQY1Raw", BigRIPSPPAC_fQY1Raw, &b_BigRIPSPPAC_fQY1Raw);
   fChain->SetBranchAddress("BigRIPSPPAC.fQY2Raw", BigRIPSPPAC_fQY2Raw, &b_BigRIPSPPAC_fQY2Raw);
   fChain->SetBranchAddress("BigRIPSPPAC.fQARaw", BigRIPSPPAC_fQARaw, &b_BigRIPSPPAC_fQARaw);
   fChain->SetBranchAddress("BigRIPSPPAC.fTX1", BigRIPSPPAC_fTX1, &b_BigRIPSPPAC_fTX1);
   fChain->SetBranchAddress("BigRIPSPPAC.fTX2", BigRIPSPPAC_fTX2, &b_BigRIPSPPAC_fTX2);
   fChain->SetBranchAddress("BigRIPSPPAC.fTY1", BigRIPSPPAC_fTY1, &b_BigRIPSPPAC_fTY1);
   fChain->SetBranchAddress("BigRIPSPPAC.fTY2", BigRIPSPPAC_fTY2, &b_BigRIPSPPAC_fTY2);
   fChain->SetBranchAddress("BigRIPSPPAC.fTA", BigRIPSPPAC_fTA, &b_BigRIPSPPAC_fTA);
   fChain->SetBranchAddress("BigRIPSPPAC.fTSumX", BigRIPSPPAC_fTSumX, &b_BigRIPSPPAC_fTSumX);
   fChain->SetBranchAddress("BigRIPSPPAC.fTSumY", BigRIPSPPAC_fTSumY, &b_BigRIPSPPAC_fTSumY);
   fChain->SetBranchAddress("BigRIPSPPAC.fTDiffX", BigRIPSPPAC_fTDiffX, &b_BigRIPSPPAC_fTDiffX);
   fChain->SetBranchAddress("BigRIPSPPAC.fTDiffY", BigRIPSPPAC_fTDiffY, &b_BigRIPSPPAC_fTDiffY);
   fChain->SetBranchAddress("BigRIPSPPAC.fX", BigRIPSPPAC_fX, &b_BigRIPSPPAC_fX);
   fChain->SetBranchAddress("BigRIPSPPAC.fY", BigRIPSPPAC_fY, &b_BigRIPSPPAC_fY);
   fChain->SetBranchAddress("BigRIPSPPAC.fFiredX", BigRIPSPPAC_fFiredX, &b_BigRIPSPPAC_fFiredX);
   fChain->SetBranchAddress("BigRIPSPPAC.fFiredY", BigRIPSPPAC_fFiredY, &b_BigRIPSPPAC_fFiredY);
   fChain->SetBranchAddress("BigRIPSPlastic", &BigRIPSPlastic_, &b_BigRIPSPlastic_);
   fChain->SetBranchAddress("BigRIPSPlastic.fUniqueID", BigRIPSPlastic_fUniqueID, &b_BigRIPSPlastic_fUniqueID);
   fChain->SetBranchAddress("BigRIPSPlastic.fBits", BigRIPSPlastic_fBits, &b_BigRIPSPlastic_fBits);
   fChain->SetBranchAddress("BigRIPSPlastic.id", BigRIPSPlastic_id, &b_BigRIPSPlastic_id);
   fChain->SetBranchAddress("BigRIPSPlastic.fpl", BigRIPSPlastic_fpl, &b_BigRIPSPlastic_fpl);
   fChain->SetBranchAddress("BigRIPSPlastic.name", BigRIPSPlastic_name, &b_BigRIPSPlastic_name);
   fChain->SetBranchAddress("BigRIPSPlastic.fDataState", BigRIPSPlastic_fDataState, &b_BigRIPSPlastic_fDataState);
   fChain->SetBranchAddress("BigRIPSPlastic.zposition", BigRIPSPlastic_zposition, &b_BigRIPSPlastic_zposition);
   fChain->SetBranchAddress("BigRIPSPlastic.zoffset", BigRIPSPlastic_zoffset, &b_BigRIPSPlastic_zoffset);
   fChain->SetBranchAddress("BigRIPSPlastic.fTLRaw", BigRIPSPlastic_fTLRaw, &b_BigRIPSPlastic_fTLRaw);
   fChain->SetBranchAddress("BigRIPSPlastic.fTRRaw", BigRIPSPlastic_fTRRaw, &b_BigRIPSPlastic_fTRRaw);
   fChain->SetBranchAddress("BigRIPSPlastic.fQLRaw", BigRIPSPlastic_fQLRaw, &b_BigRIPSPlastic_fQLRaw);
   fChain->SetBranchAddress("BigRIPSPlastic.fQRRaw", BigRIPSPlastic_fQRRaw, &b_BigRIPSPlastic_fQRRaw);
   fChain->SetBranchAddress("BigRIPSPlastic.fTime", BigRIPSPlastic_fTime, &b_BigRIPSPlastic_fTime);
   fChain->SetBranchAddress("BigRIPSPlastic.fTimeL", BigRIPSPlastic_fTimeL, &b_BigRIPSPlastic_fTimeL);
   fChain->SetBranchAddress("BigRIPSPlastic.fTimeR", BigRIPSPlastic_fTimeR, &b_BigRIPSPlastic_fTimeR);
   fChain->SetBranchAddress("BigRIPSPlastic.fTimeLSlew", BigRIPSPlastic_fTimeLSlew, &b_BigRIPSPlastic_fTimeLSlew);
   fChain->SetBranchAddress("BigRIPSPlastic.fTimeRSlew", BigRIPSPlastic_fTimeRSlew, &b_BigRIPSPlastic_fTimeRSlew);
   fChain->SetBranchAddress("BigRIPSPlastic.fTimeSlew", BigRIPSPlastic_fTimeSlew, &b_BigRIPSPlastic_fTimeSlew);
   fChain->SetBranchAddress("BigRIPSIC", &BigRIPSIC_, &b_BigRIPSIC_);
   fChain->SetBranchAddress("BigRIPSIC.fUniqueID", BigRIPSIC_fUniqueID, &b_BigRIPSIC_fUniqueID);
   fChain->SetBranchAddress("BigRIPSIC.fBits", BigRIPSIC_fBits, &b_BigRIPSIC_fBits);
   fChain->SetBranchAddress("BigRIPSIC.id", BigRIPSIC_id, &b_BigRIPSIC_id);
   fChain->SetBranchAddress("BigRIPSIC.fpl", BigRIPSIC_fpl, &b_BigRIPSIC_fpl);
   fChain->SetBranchAddress("BigRIPSIC.name", BigRIPSIC_name, &b_BigRIPSIC_name);
   fChain->SetBranchAddress("BigRIPSIC.fDataState", BigRIPSIC_fDataState, &b_BigRIPSIC_fDataState);
   fChain->SetBranchAddress("BigRIPSIC.zcoef[2]", BigRIPSIC_zcoef, &b_BigRIPSIC_zcoef);
   fChain->SetBranchAddress("BigRIPSIC.ionpair", BigRIPSIC_ionpair, &b_BigRIPSIC_ionpair);
   fChain->SetBranchAddress("BigRIPSIC.nhitchannel", BigRIPSIC_nhitchannel, &b_BigRIPSIC_nhitchannel);
   fChain->SetBranchAddress("BigRIPSIC.fADC[12]", BigRIPSIC_fADC, &b_BigRIPSIC_fADC);
   fChain->SetBranchAddress("BigRIPSIC.fEnergy[12]", BigRIPSIC_fEnergy, &b_BigRIPSIC_fEnergy);
   fChain->SetBranchAddress("BigRIPSIC.fRawADCSqSum", BigRIPSIC_fRawADCSqSum, &b_BigRIPSIC_fRawADCSqSum);
   fChain->SetBranchAddress("BigRIPSIC.fRawADCAvSum", BigRIPSIC_fRawADCAvSum, &b_BigRIPSIC_fRawADCAvSum);
   fChain->SetBranchAddress("BigRIPSIC.fCalMeVSqSum", BigRIPSIC_fCalMeVSqSum, &b_BigRIPSIC_fCalMeVSqSum);
   fChain->SetBranchAddress("BigRIPSIC.fCalMeVAvSum", BigRIPSIC_fCalMeVAvSum, &b_BigRIPSIC_fCalMeVAvSum);
   fChain->SetBranchAddress("BigRIPSFocalPlane", &BigRIPSFocalPlane_, &b_BigRIPSFocalPlane_);
   fChain->SetBranchAddress("BigRIPSFocalPlane.fUniqueID", BigRIPSFocalPlane_fUniqueID, &b_BigRIPSFocalPlane_fUniqueID);
   fChain->SetBranchAddress("BigRIPSFocalPlane.fBits", BigRIPSFocalPlane_fBits, &b_BigRIPSFocalPlane_fBits);
   fChain->SetBranchAddress("BigRIPSFocalPlane.id", BigRIPSFocalPlane_id, &b_BigRIPSFocalPlane_id);
   fChain->SetBranchAddress("BigRIPSFocalPlane.fpl", BigRIPSFocalPlane_fpl, &b_BigRIPSFocalPlane_fpl);
   fChain->SetBranchAddress("BigRIPSFocalPlane.name", BigRIPSFocalPlane_name, &b_BigRIPSFocalPlane_name);
   fChain->SetBranchAddress("BigRIPSFocalPlane.fDataState", BigRIPSFocalPlane_fDataState, &b_BigRIPSFocalPlane_fDataState);
   fChain->SetBranchAddress("BigRIPSFocalPlane.opt_vector", BigRIPSFocalPlane_opt_vector, &b_BigRIPSFocalPlane_opt_vector);
   fChain->SetBranchAddress("BigRIPSFocalPlane.X", BigRIPSFocalPlane_X, &b_BigRIPSFocalPlane_X);
   fChain->SetBranchAddress("BigRIPSFocalPlane.A", BigRIPSFocalPlane_A, &b_BigRIPSFocalPlane_A);
   fChain->SetBranchAddress("BigRIPSFocalPlane.Y", BigRIPSFocalPlane_Y, &b_BigRIPSFocalPlane_Y);
   fChain->SetBranchAddress("BigRIPSFocalPlane.B", BigRIPSFocalPlane_B, &b_BigRIPSFocalPlane_B);
   fChain->SetBranchAddress("BigRIPSFocalPlane.nfired_ppacx", BigRIPSFocalPlane_nfired_ppacx, &b_BigRIPSFocalPlane_nfired_ppacx);
   fChain->SetBranchAddress("BigRIPSFocalPlane.nfired_ppacy", BigRIPSFocalPlane_nfired_ppacy, &b_BigRIPSFocalPlane_nfired_ppacy);
   fChain->SetBranchAddress("BigRIPSFocalPlane.zpos", BigRIPSFocalPlane_zpos, &b_BigRIPSFocalPlane_zpos);
   fChain->SetBranchAddress("BigRIPSFocalPlane.zpos_offset", BigRIPSFocalPlane_zpos_offset, &b_BigRIPSFocalPlane_zpos_offset);
   fChain->SetBranchAddress("DALINaI", &DALINaI_, &b_DALINaI_);
   fChain->SetBranchAddress("DALINaI.fUniqueID", DALINaI_fUniqueID, &b_DALINaI_fUniqueID);
   fChain->SetBranchAddress("DALINaI.fBits", DALINaI_fBits, &b_DALINaI_fBits);
   fChain->SetBranchAddress("DALINaI.id", DALINaI_id, &b_DALINaI_id);
   fChain->SetBranchAddress("DALINaI.fpl", DALINaI_fpl, &b_DALINaI_fpl);
   fChain->SetBranchAddress("DALINaI.name", DALINaI_name, &b_DALINaI_name);
   fChain->SetBranchAddress("DALINaI.fDataState", DALINaI_fDataState, &b_DALINaI_fDataState);
   fChain->SetBranchAddress("DALINaI.fADC", DALINaI_fADC, &b_DALINaI_fADC);
   fChain->SetBranchAddress("DALINaI.fTDC", DALINaI_fTDC, &b_DALINaI_fTDC);
   fChain->SetBranchAddress("DALINaI.layer", DALINaI_layer, &b_DALINaI_layer);
   fChain->SetBranchAddress("DALINaI.theta", DALINaI_theta, &b_DALINaI_theta);
   fChain->SetBranchAddress("DALINaI.fXPos", DALINaI_fXPos, &b_DALINaI_fXPos);
   fChain->SetBranchAddress("DALINaI.fYPos", DALINaI_fYPos, &b_DALINaI_fYPos);
   fChain->SetBranchAddress("DALINaI.fZPos", DALINaI_fZPos, &b_DALINaI_fZPos);
   fChain->SetBranchAddress("DALINaI.costheta", DALINaI_costheta, &b_DALINaI_costheta);
   fChain->SetBranchAddress("DALINaI.fEnergy", DALINaI_fEnergy, &b_DALINaI_fEnergy);
   fChain->SetBranchAddress("DALINaI.fDoppCorEnergy", DALINaI_fDoppCorEnergy, &b_DALINaI_fDoppCorEnergy);
   fChain->SetBranchAddress("DALINaI.fEnergyWithoutT", DALINaI_fEnergyWithoutT, &b_DALINaI_fEnergyWithoutT);
   fChain->SetBranchAddress("DALINaI.fTime", DALINaI_fTime, &b_DALINaI_fTime);
   fChain->SetBranchAddress("DALINaI.fTimeOffseted", DALINaI_fTimeOffseted, &b_DALINaI_fTimeOffseted);
   fChain->SetBranchAddress("DALINaI.fTimeTrueEnergy", DALINaI_fTimeTrueEnergy, &b_DALINaI_fTimeTrueEnergy);
   fChain->SetBranchAddress("DALINaI.fTimeTrueDoppCorEnergy", DALINaI_fTimeTrueDoppCorEnergy, &b_DALINaI_fTimeTrueDoppCorEnergy);
   fChain->SetBranchAddress("DALINaI.fTimeTrueDoppCorVertexEnergy", DALINaI_fTimeTrueDoppCorVertexEnergy, &b_DALINaI_fTimeTrueDoppCorVertexEnergy);
   fChain->SetBranchAddress("DALINaI.fTimeTrueTime", DALINaI_fTimeTrueTime, &b_DALINaI_fTimeTrueTime);
   fChain->SetBranchAddress("DALINaI.fTimeTrueTimeOffseted", DALINaI_fTimeTrueTimeOffseted, &b_DALINaI_fTimeTrueTimeOffseted);
   fChain->SetBranchAddress("BigRIPSRIPS", &BigRIPSRIPS_, &b_BigRIPSRIPS_);
   fChain->SetBranchAddress("BigRIPSRIPS.fUniqueID", BigRIPSRIPS_fUniqueID, &b_BigRIPSRIPS_fUniqueID);
   fChain->SetBranchAddress("BigRIPSRIPS.fBits", BigRIPSRIPS_fBits, &b_BigRIPSRIPS_fBits);
   fChain->SetBranchAddress("BigRIPSRIPS.id", BigRIPSRIPS_id, &b_BigRIPSRIPS_id);
   fChain->SetBranchAddress("BigRIPSRIPS.fpl", BigRIPSRIPS_fpl, &b_BigRIPSRIPS_fpl);
   fChain->SetBranchAddress("BigRIPSRIPS.name", BigRIPSRIPS_name, &b_BigRIPSRIPS_name);
   fChain->SetBranchAddress("BigRIPSRIPS.fDataState", BigRIPSRIPS_fDataState, &b_BigRIPSRIPS_fDataState);
   fChain->SetBranchAddress("BigRIPSRIPS.upstream_fpl", BigRIPSRIPS_upstream_fpl, &b_BigRIPSRIPS_upstream_fpl);
   fChain->SetBranchAddress("BigRIPSRIPS.downstream_fpl", BigRIPSRIPS_downstream_fpl, &b_BigRIPSRIPS_downstream_fpl);
   fChain->SetBranchAddress("BigRIPSRIPS.center_brho", BigRIPSRIPS_center_brho, &b_BigRIPSRIPS_center_brho);
   fChain->SetBranchAddress("BigRIPSRIPS.brho", BigRIPSRIPS_brho, &b_BigRIPSRIPS_brho);
   fChain->SetBranchAddress("BigRIPSRIPS.length", BigRIPSRIPS_length, &b_BigRIPSRIPS_length);
   fChain->SetBranchAddress("BigRIPSRIPS.matrix", BigRIPSRIPS_matrix, &b_BigRIPSRIPS_matrix);
   fChain->SetBranchAddress("BigRIPSRIPS.delta", BigRIPSRIPS_delta, &b_BigRIPSRIPS_delta);
   fChain->SetBranchAddress("BigRIPSRIPS.angle", BigRIPSRIPS_angle, &b_BigRIPSRIPS_angle);
   fChain->SetBranchAddress("BigRIPSRIPS.dipolename", BigRIPSRIPS_dipolename, &b_BigRIPSRIPS_dipolename);
   fChain->SetBranchAddress("BigRIPSTOF", &BigRIPSTOF_, &b_BigRIPSTOF_);
   fChain->SetBranchAddress("BigRIPSTOF.fUniqueID", BigRIPSTOF_fUniqueID, &b_BigRIPSTOF_fUniqueID);
   fChain->SetBranchAddress("BigRIPSTOF.fBits", BigRIPSTOF_fBits, &b_BigRIPSTOF_fBits);
   fChain->SetBranchAddress("BigRIPSTOF.id", BigRIPSTOF_id, &b_BigRIPSTOF_id);
   fChain->SetBranchAddress("BigRIPSTOF.fpl", BigRIPSTOF_fpl, &b_BigRIPSTOF_fpl);
   fChain->SetBranchAddress("BigRIPSTOF.name", BigRIPSTOF_name, &b_BigRIPSTOF_name);
   fChain->SetBranchAddress("BigRIPSTOF.fDataState", BigRIPSTOF_fDataState, &b_BigRIPSTOF_fDataState);
   fChain->SetBranchAddress("BigRIPSTOF.tof", BigRIPSTOF_tof, &b_BigRIPSTOF_tof);
   fChain->SetBranchAddress("BigRIPSTOF.clight", BigRIPSTOF_clight, &b_BigRIPSTOF_clight);
   fChain->SetBranchAddress("BigRIPSTOF.length", BigRIPSTOF_length, &b_BigRIPSTOF_length);
   fChain->SetBranchAddress("BigRIPSTOF.ulength", BigRIPSTOF_ulength, &b_BigRIPSTOF_ulength);
   fChain->SetBranchAddress("BigRIPSTOF.dlength", BigRIPSTOF_dlength, &b_BigRIPSTOF_dlength);
   fChain->SetBranchAddress("BigRIPSTOF.upstream_plname", BigRIPSTOF_upstream_plname, &b_BigRIPSTOF_upstream_plname);
   fChain->SetBranchAddress("BigRIPSTOF.downstream_plname", BigRIPSTOF_downstream_plname, &b_BigRIPSTOF_downstream_plname);
   fChain->SetBranchAddress("BigRIPSTOF.upstream_plfpl", BigRIPSTOF_upstream_plfpl, &b_BigRIPSTOF_upstream_plfpl);
   fChain->SetBranchAddress("BigRIPSTOF.downstream_plfpl", BigRIPSTOF_downstream_plfpl, &b_BigRIPSTOF_downstream_plfpl);
   fChain->SetBranchAddress("BigRIPSTOF.time_offset", BigRIPSTOF_time_offset, &b_BigRIPSTOF_time_offset);
   fChain->SetBranchAddress("BigRIPSBeam", &BigRIPSBeam_, &b_BigRIPSBeam_);
   fChain->SetBranchAddress("BigRIPSBeam.fUniqueID", BigRIPSBeam_fUniqueID, &b_BigRIPSBeam_fUniqueID);
   fChain->SetBranchAddress("BigRIPSBeam.fBits", BigRIPSBeam_fBits, &b_BigRIPSBeam_fBits);
   fChain->SetBranchAddress("BigRIPSBeam.id", BigRIPSBeam_id, &b_BigRIPSBeam_id);
   fChain->SetBranchAddress("BigRIPSBeam.fpl", BigRIPSBeam_fpl, &b_BigRIPSBeam_fpl);
   fChain->SetBranchAddress("BigRIPSBeam.name", BigRIPSBeam_name, &b_BigRIPSBeam_name);
   fChain->SetBranchAddress("BigRIPSBeam.fDataState", BigRIPSBeam_fDataState, &b_BigRIPSBeam_fDataState);
   fChain->SetBranchAddress("BigRIPSBeam.brho", BigRIPSBeam_brho, &b_BigRIPSBeam_brho);
   fChain->SetBranchAddress("BigRIPSBeam.aoq", BigRIPSBeam_aoq, &b_BigRIPSBeam_aoq);
   fChain->SetBranchAddress("BigRIPSBeam.zet", BigRIPSBeam_zet, &b_BigRIPSBeam_zet);
   fChain->SetBranchAddress("BigRIPSBeam.beta", BigRIPSBeam_beta, &b_BigRIPSBeam_beta);
   fChain->SetBranchAddress("BigRIPSBeam.clight", BigRIPSBeam_clight, &b_BigRIPSBeam_clight);
   fChain->SetBranchAddress("BigRIPSBeam.mnucleon", BigRIPSBeam_mnucleon, &b_BigRIPSBeam_mnucleon);
   fChain->SetBranchAddress("BigRIPSBeam.nrips", BigRIPSBeam_nrips, &b_BigRIPSBeam_nrips);
   fChain->SetBranchAddress("BigRIPSBeam.ripsname[2]", BigRIPSBeam_ripsname, &b_BigRIPSBeam_ripsname);
   fChain->SetBranchAddress("BigRIPSBeam.tofname", BigRIPSBeam_tofname, &b_BigRIPSBeam_tofname);
   fChain->SetBranchAddress("BigRIPSBeam.icname", BigRIPSBeam_icname, &b_BigRIPSBeam_icname);
   fChain->SetBranchAddress("triggerbit", &triggerbit, &b_triggerbit);
   fChain->SetBranchAddress("neve", &neve, &b_neve);
   fChain->SetBranchAddress("tgtx", &tgtx, &b_tgtx);
   fChain->SetBranchAddress("tgty", &tgty, &b_tgty);
   fChain->SetBranchAddress("tgta", &tgta, &b_tgta);
   fChain->SetBranchAddress("tgtb", &tgtb, &b_tgtb);
   fChain->SetBranchAddress("FX", FX, &b_FX);
   fChain->SetBranchAddress("FY", FY, &b_FY);
   fChain->SetBranchAddress("FA", FA, &b_FA);
   fChain->SetBranchAddress("FB", FB, &b_FB);
   fChain->SetBranchAddress("F8T", &F8T, &b_F8T);
   fChain->SetBranchAddress("delta", delta, &b_delta);
   fChain->SetBranchAddress("F10B", &F10B, &b_F10B);
   fChain->SetBranchAddress("F11B", &F11B, &b_F11B);
   fChain->SetBranchAddress("F3PLA_QL", &F3PLA_QL, &b_F3PLA_QL);
   fChain->SetBranchAddress("F3PLA_QR", &F3PLA_QR, &b_F3PLA_QR);
   fChain->SetBranchAddress("F7PLA_QL", &F7PLA_QL, &b_F7PLA_QL);
   fChain->SetBranchAddress("F7PLA_QR", &F7PLA_QR, &b_F7PLA_QR);
   fChain->SetBranchAddress("F8PLA_QL", &F8PLA_QL, &b_F8PLA_QL);
   fChain->SetBranchAddress("F8PLA_QR", &F8PLA_QR, &b_F8PLA_QR);
   fChain->SetBranchAddress("F11PLA_QL", &F11PLA_QL, &b_F11PLA_QL);
   fChain->SetBranchAddress("F11PLA_QR", &F11PLA_QR, &b_F11PLA_QR);
   fChain->SetBranchAddress("DELTA", DELTA, &b_DELTA);
   fChain->SetBranchAddress("ANGLE", ANGLE, &b_ANGLE);
   fChain->SetBranchAddress("BRHO", BRHO, &b_BRHO);
   fChain->SetBranchAddress("TOF", TOF, &b_TOF);
   fChain->SetBranchAddress("BETA", BETA, &b_BETA);
   fChain->SetBranchAddress("AOQ", AOQ, &b_AOQ);
   fChain->SetBranchAddress("ZET", ZET, &b_ZET);
   fChain->SetBranchAddress("ESq", &ESq, &b_ESq);
   fChain->SetBranchAddress("TKE", &TKE, &b_TKE);
   fChain->SetBranchAddress("ASUM", &ASUM, &b_ASUM);
   fChain->SetBranchAddress("dalimult", &dalimult, &b_dalimult);
   fChain->SetBranchAddress("dalitimetruemult", &dalitimetruemult, &b_dalitimetruemult);
   fChain->SetBranchAddress("LOWGAIN", LOWGAIN, &b_LOWGAIN);
   fChain->SetBranchAddress("LOWGAINC", LOWGAINC, &b_LOWGAINC);
   fChain->SetBranchAddress("MIDGAIN", MIDGAIN, &b_MIDGAIN);
   fChain->SetBranchAddress("MIDGAINC", MIDGAINC, &b_MIDGAINC);
   fChain->SetBranchAddress("MIDGAINCT", MIDGAINCT, &b_MIDGAINCT);
   fChain->SetBranchAddress("HIGAIN", HIGAIN, &b_HIGAIN);
   fChain->SetBranchAddress("HIGAINC", HIGAINC, &b_HIGAINC);
   fChain->SetBranchAddress("QTCTime", QTCTime, &b_QTCtime);
   fChain->SetBranchAddress("QTCTimeC", QTCTimeC, &b_QTCtimeC);
   fChain->SetBranchAddress("MIDGAINCTA", &MIDGAINCTA, &b_MIDGAINCTA);
   Notify();
}

Bool_t test::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void test::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t test::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef test_cxx
