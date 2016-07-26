#define PID_2_cxx
#include "PID_2.h"
#include <TH2.h>
#include "TH2F.h"
#include <TStyle.h>
#include <TCanvas.h>
#include "TF1.h"
#include "TClonesArray.h"

#include<fstream>
#include<stdlib.h>
#include "TString.h"
#include <iostream>
#include "Riostream.h"
#include <string>
#include <time.h>
#include <sstream>
#include "TLine.h"


void PID_2::Loop(Bool_t beta, short var, short aoq){
//   In a ROOT session, you can do:
//      Root > .L PID_2.C
//      Root > PID_2 t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
  
  fChain->SetBranchStatus("*",0);
  fChain->SetBranchStatus("BigRIPSRIPS",1);
  fChain->SetBranchStatus("BigRIPSRIPS.downstream_fpl",1);
  fChain->SetBranchStatus("BigRIPSRIPS.upstream_fpl",1);
  fChain->SetBranchStatus("BigRIPSRIPS.angle",1);
  fChain->SetBranchStatus("BigRIPSPPAC",1);
  fChain->SetBranchStatus("BigRIPSPPAC.fpl",1);
  fChain->SetBranchStatus("BigRIPSPPAC.fX",1);
  fChain->SetBranchStatus("BigRIPSPPAC.fY",1);
  fChain->SetBranchStatus("AOQ");
  fChain->SetBranchStatus("ZET");
  fChain->SetBranchStatus("delta");
  fChain->SetBranchStatus("BETA");
  fChain->SetBranchStatus("BigRIPSPlastic");
  fChain->SetBranchStatus("BigRIPSPlastic.fTime");
  fChain->SetBranchStatus("BigRIPSPlastic.fpl");


  //Reset();
  //Histo(beta);
  Read();
  TH2F* h2;
  TH2F* h2_AOQ_ZET = new TH2F("fh2_AOQ_ZET",Form("aoq vs ZET %i",aoq),200,40,55,800,2.3,2.7);
  TH2F* fh2_timevs_evnumber = new TH2F("fh2_timevs_evnumber","Events vs time",10000,-0.5,1e4-0.5,800,400,2800);
  if( beta ){
    h2 = new TH2F(Form("fh2_aoq_%02d",aoq),Form("aoq %i vs var",3),1000,0.5,0.6,800,2.5,2.7);
    h2_z = new TH2F("h2_z","zet vs var",1000,0.5,0.6,200,45,55);
  }else{
    h2 = new TH2F(Form("fh2_aoq_%02d",aoq),Form("aoq %i vs var",3),400,-40,40,600,2.5,2.8);
    h2_z = new TH2F("h2_z","zet vs var",400,-20,20,200,45,55);
  }
  Long64_t nentries = fChain->GetEntriesFast();
  if (fChain == 0) return;
   Long64_t nbytes = 0, nb = 0;


   short plane = (short)var/3;
   short var_3 = var%3;
   cout<< "plane: " << plane << "  var_3: " << var_3<<endl;

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     reset();
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      if( jentry < 1e4 )
	for(short ii = 0; ii <kMaxBigRIPSPlastic; ii++)
	  if(BigRIPSPlastic_fpl[ii] == 7)
	    fh2_timevs_evnumber->Fill(jentry,BigRIPSPlastic_fTime[ii]);
      //      hitsperevent(jentry);
      //global cut
      if( !globalcut())continue;
      
      for(short irips = 0; irips < kMaxBigRIPSRIPS; irips++){//getting the variables

	if(BigRIPSRIPS_upstream_fpl[irips] == 3){
	  fin_var[0] = BigRIPSRIPS_angle[irips];
	  //	  angle_b_F3 = BigRIPSRIPS_angle[irips];
	}else if(BigRIPSRIPS_upstream_fpl[irips] == 5){
	  fin_var[1] = BigRIPSRIPS_angle[irips];
	}else if(BigRIPSRIPS_upstream_fpl[irips] == 8){
	  if(fvar[0][0]==-100){
	    fout_var[0][0] = BigRIPSRIPS_angle[irips];
	    fvar[0][0] = BigRIPSRIPS_angle[irips];
	  }
	  //	  if( outgoingcut() )
	  //h2->Fill(fvar[plane][var_3],AOQ[3]);//Correction_3(3, var));

	}

	/*if(BigRIPSRIPS_downstream_fpl[irips] == 9){
	  fout_var[0][0][0] = BigRIPSRIPS_angle[irips];
	  }else*/
	if(BigRIPSRIPS_downstream_fpl[irips] == 10){
	  fout_var[0][1] = BigRIPSRIPS_angle[irips];
	  fout_var[1][0] = BigRIPSRIPS_angle[irips];
	  //fout_var[2][0] = BigRIPSRIPS_angle[irips];
	  fvar[1][0] = BigRIPSRIPS_angle[irips];
	}else if(BigRIPSRIPS_downstream_fpl[irips] == 11){
	  fout_var[1][1] = BigRIPSRIPS_angle[irips];
	  //fout_var[2][1] = BigRIPSRIPS_angle[irips];
	  fvar[2][0] = BigRIPSRIPS_angle[irips];
	}else continue;

      }//for(RIPS)
      /*debugging
	continue;
	for(short irips = 0; irips<8;irips++)
	if(BigRIPSRIPS_upstream_fpl[irips] == 8 && fvar[0][0]!=BigRIPSRIPS_angle[irips])cout<<"yes"<<endl;
      */
      for(short ippac = 0; ippac < kMaxBigRIPSPPAC; ippac++){//getting the variables

	if(BigRIPSPPAC_fpl[ippac] == 3){
	  fin_var[2] = BigRIPSPPAC_fX[ippac];
	  fin_var[4] = BigRIPSPPAC_fY[ippac];
	}else if(BigRIPSPPAC_fpl[ippac] == 5){
	  fin_var[3] = BigRIPSPPAC_fX[ippac];
	  fin_var[5] = BigRIPSPPAC_fY[ippac];
	}else if(BigRIPSPPAC_fpl[ippac] == 8){
	  fout_var[0][2] = BigRIPSPPAC_fX[ippac];
	  fout_var[0][4] = BigRIPSPPAC_fY[ippac];
	  fvar[0][1] = BigRIPSPPAC_fX[ippac];
	  fvar[0][2] = BigRIPSPPAC_fY[ippac];
	}else if(BigRIPSPPAC_fpl[ippac] == 10){
	  fout_var[0][3] = BigRIPSPPAC_fX[ippac];
	  fout_var[0][5] = BigRIPSPPAC_fY[ippac];
	  fout_var[1][2] = BigRIPSPPAC_fX[ippac];
	  fout_var[1][4] = BigRIPSPPAC_fY[ippac];

	  fvar[1][1] = BigRIPSPPAC_fX[ippac];
	  fvar[1][2] = BigRIPSPPAC_fY[ippac];
	}else if(BigRIPSPPAC_fpl[ippac] == 11){
	  fout_var[1][3] = BigRIPSPPAC_fX[ippac];
	  fout_var[1][5] = BigRIPSPPAC_fY[ippac];

	  fvar[2][1] = BigRIPSPPAC_fX[ippac];
	  fvar[2][2] = BigRIPSPPAC_fY[ippac];
	}else continue;
      
      }//ippac

      //h2->Fill(fin_var[var],AOQ[0]);
      //continue;

      //      h2->Fill(fout_var[0][var],Correction(1,var));
      //continue;
      if(aoq==0)h2_AOQ_ZET->Fill(zet(aoq, var),Correction(aoq,100));
      else if( TMath::Abs(delta[2] - delta[4] )< 0.3 ){
	if(aoq==3)h2_AOQ_ZET->Fill(zet(aoq, var),Correction_3(3,100));
	else h2_AOQ_ZET->Fill(zet(aoq, var),Correction(aoq,100));
      }
      if(aoq==0){
	if(beta)h2->Fill(BETA[aoq], Correction(aoq,var));
	else h2->Fill(fin_var[var], Correction(aoq,var));
      } else if( outgoingcut() ){
	if(aoq!=3){
	  if(beta)h2->Fill(BETA[aoq], Correction(aoq,var));
	  else h2->Fill(fout_var[aoq-1][var], Correction(aoq,var));
	}else if(aoq==3){
	  if( beta )h2->Fill(BETA[3],           Correction_3(3, var));
	  else      h2->Fill(fvar[plane][var_3],AOQ[3]);//Correction_3(3, var));
	}
      }
      if( aoq != 0 && outgoingcut() )
	if( beta ) h2_z->Fill(BETA[aoq],zet(aoq, var));
	else if( aoq <  3) h2_z->Fill(fout_var[aoq-1][var],zet(aoq, var));
	else h2_z->Fill(fvar[plane][var_3],zet(aoq, var));
      else if(aoq == 0)
	if( beta ) h2_z->Fill(BETA[aoq],zet(aoq, var));
	else       h2_z->Fill(fin_var[var],zet(aoq, var));


/*****
      if(var<6){
	cout<<"hi"<<endl;
	fh2_aoq[0]->Fill(fin_var[var],Correction(0, var));
	if( !outgoingcut() )continue;
	cout<<"hi0"<<endl;
	fh2_aoq[1]->Fill(fout_var[0][var],Correction(1, var));
	cout<<"hi1"<<endl;
	fh2_aoq[2]->Fill(fout_var[0][var],Correction(2, var));
	cout<<"hi2"<<endl;
      }
      if( !outgoingcut() )continue;
      fh2_aoq[3]->Fill(fvar[var/3][var%3],Correction_3(3, var));
      //      cout<<"hi3"<<endl;
*****/
      
   }//for
   
   FitZET(h2_z,aoq,var);
   FitSlices(h2,aoq,var);
   cout<< "Drawing" <<endl;
   TCanvas* can = new TCanvas("aoq","aoq");
   h2->Draw("colz");   
   fh1_aoq_1[aoq]->Draw("same");
   return;
   can->Divide(2,2);
   for(short ii = 3; ii<4; ii++){
     can->cd(ii+1);
     fh2_aoq[ii]->Draw("colz");
     fh1_aoq_1[ii]->Draw("same");
   }
}//Loop()



void PID_2::Histo(Bool_t beta){

  Float_t minX = -10, maxX = 10;
  if(beta){
    minX=0.5;
    maxX=0.6;
  }
  //for(short jj = 0; jj<4; jj++)
  //  fh2_aoq[jj] = new TH2F(Form("fh2_aoq_%02d",jj),Form("aoq %i vs var",jj),200,minX,maxX,400,2.6,2.7);
  


}

void PID_2::Reset(){

  for(short jj = 0; jj<4; jj++){
    fh2_aoq[jj]->Delete();
    fh1_aoq_1[jj]->Delete();
  }
  h2_z->Delete();
  fdummy->Clear();
}


void PID_2::Read(){
  ReadIncoming();
  ReadOutgoing_1();
  ReadOutgoing_3();
  ReadZET();
}

void PID_2::ReadIncoming(){
 
  ifstream aoqfile;
  std::string off, lin, quadr, line;
  aoqfile.open(Form("%s/AOQ_Correction_in.txt",path));
  short ii =0;
  while(aoqfile.good()){
    getline(aoqfile,line);
    std::istringstream iss(line);
    iss >> off >> lin >> quadr;
    fcorr_Lin[0][ii][0]=atof(off.c_str());
    fcorr_Lin[0][ii][1]=atof(lin.c_str());
    fcorr_Lin[0][ii][2]=atof(quadr.c_str());
    ii++;
    if(ii == 7)break;
  }
  

}


void PID_2::ReadZET(){
  
  ifstream aoqfile;
  std::string off, lin, quadr, line;
  aoqfile.open(Form("%s/ZET_Correction.txt",path));
  short ii =0, jj=0;
  while(aoqfile.good()){
    getline(aoqfile,line);
    std::istringstream iss(line);
    iss >> off >> lin >> quadr;
    cout<< off << " " << lin << " " << quadr <<endl;
    fzet_corr_Lin[jj][ii][0]=atof(off.c_str());
    fzet_corr_Lin[jj][ii][1]=atof(lin.c_str());
    fzet_corr_Lin[jj][ii][2]=atof(quadr.c_str());
    ii++;
    if(ii == 4){
      ii=0;
      jj++;
    }
    if( jj == 2 ) break;
  }
  

}

void PID_2::ReadOutgoing_1(){
  
  ifstream aoqfile;
  std::string off, lin, quadr, line;
  for(short jj =1; jj<3;jj++){
    aoqfile.open(Form("%s/AOQ_Correction_out%i.txt",path,jj));
    short ii =0;
    while(aoqfile.good()){
      getline(aoqfile,line);
      std::istringstream iss(line);
      iss >> off >> lin >> quadr;
      fcorr_Lin[jj-1][ii][0]=atof(off.c_str());
      fcorr_Lin[jj-1][ii][1]=atof(lin.c_str());
      fcorr_Lin[jj-1][ii][2]=atof(quadr.c_str());
      ii++;
      if(ii == 7)break;
    }
  }//for(ii)


}

void PID_2::ReadOutgoing_3(){
  
  std::string var1,var2,var3, line;
  ifstream aoqfile;
  aoqfile.open(Form("%s/AOQ_Correction_out3.txt",path));
  short plane = 0, var = 0;
  while(aoqfile.good()){

    getline(aoqfile,line);
    std::istringstream iss(line);
    iss >> var1 >> var2 >> var3;
    fcorr_Lin_3[plane][var][0]=atof(var1.c_str());
    fcorr_Lin_3[plane][var][1]=atof(var2.c_str());
    fcorr_Lin_3[plane][var][2]=atof(var3.c_str());
    var++;
    if(var == 1 && plane == 3)break;
    if(var == 3){
      plane++;
      var = 0;
    }
  }

}

Bool_t PID_2::globalcut(){

  if( TMath::Abs(ZET[0]-50)>0.5 )return kFALSE;
  else if( TMath::Abs(AOQ[0]-f_aoq)>0.005)return kFALSE;

  return kTRUE;

}

Bool_t PID_2::outgoingcut(){

  if( TMath::Abs(delta[2] - delta[4] )> 0.3 ) return kFALSE;
  //  if( TMath::Abs(ZET[2]-50)>0.5 )return kFALSE;
  if( TMath::Abs(AOQ[2]-f_aoq)>0.008)return kFALSE;
  if( TMath::Abs(AOQ[3]-f_aoq)>0.008)return kFALSE;

  return kTRUE;

}

void PID_2::FitSlices(TH2F* h2, short aoq, short var){

  cout<< "Fitting"<<endl;
  TObjArray* dummy = new TObjArray();
  dummy->SetOwner(kTRUE);

  TF1* gausfunc = new TF1("gausfunc","gaus",f_aoq-0.02,f_aoq+0.02);//(0)+pol1(3)
  h2->FitSlicesY(gausfunc,0,-1,0,"QNRG4",dummy);
  TH1F* h1 = (TH1F*) dummy->FindObject(Form("fh2_aoq_%02d_1",aoq));
  TF1* Linfunc = new TF1("Linfunc",Form("[1]*(x+[0])+%f",f_aoq),-10,10);
  //TF1* Linfunc = new TF1("Linfunc","[1]*(x+[0])^2+[2]+2.64",-10,10);
  Linfunc->SetParameters(0,0);
  //  Linfunc->SetParLimits(0,-1e-3,1e-3);
  //  Linfunc->SetParLimits(1,-1e-3,1e-3);
  h1->Fit(Linfunc);
  if(aoq!=3){
    fcorr_Lin[aoq][var][1] = Linfunc->GetParameter(1);
    fcorr_Lin[aoq][var][0] = Linfunc->GetParameter(0);
  }else{
    short plane = (short)var/3;
    short var_ = var%3;
    fcorr_Lin_3[plane][var_][1] = Linfunc->GetParameter(1);
    fcorr_Lin_3[plane][var_][0] = Linfunc->GetParameter(0);
  }
  fh1_aoq_1[aoq] = h1;
}//FitSlices


void PID_2::FitZET(TH2F* h2, short aoq, short var){

  cout<< "Fitting"<<endl;
  TObjArray* dummy = new TObjArray();
  dummy->SetOwner(kTRUE);

  TF1* gausfunc = new TF1("gausfunc","gaus",49,51);//(0)+pol1(3)
  h2->FitSlicesY(gausfunc,0,-1,0,"QNRG4",dummy);
  TH1F* h1 = (TH1F*) dummy->FindObject("h2_z_1");
  //TF1* Linfunc = new TF1("Linfunc","[1]*(x+[0])+2.64",-10,10);
  TF1* Linfunc = new TF1("Linfunc","[1]*(x+[0])+50.",-10,10);
  Linfunc->SetParameters(0,0);
  //  Linfunc->SetParLimits(0,-1e-3,1e-3);
  //  Linfunc->SetParLimits(1,-1e-3,1e-3);
  h1->Fit(Linfunc);
  //  fzet_corr_Lin[aoq][var][2] = Linfunc->GetParameter(2);
  fzet_corr_Lin[aoq][var][1] = Linfunc->GetParameter(1);
  fzet_corr_Lin[aoq][var][0] = Linfunc->GetParameter(0);
  
  fh1_zet_1[aoq] = h1;
}//FitSlices



Double_t PID_2::Correction(short aoq, short var){

  Double_t naoq = AOQ[aoq];
  if(aoq>0){
    for(short jj = 0; jj < 6; jj++){
      if( jj == var ) continue;
      naoq = naoq - (fout_var[aoq-1][jj]+fcorr_Lin[aoq][jj][0])*fcorr_Lin[aoq][jj][1];
    }
  }else{
    for(short jj = 0; jj < 6; jj++){
      if( jj == var ) continue;
      naoq = naoq - (fin_var[jj]+fcorr_Lin[aoq][jj][0])*fcorr_Lin[aoq][jj][1];
    }
  }
  if( var != 6 )
    naoq = naoq - (BETA[aoq]+fcorr_Lin[aoq][6][0])*fcorr_Lin[aoq][6][1];

  return naoq;

}


Double_t PID_2::Correction_3(short aoq, short var){

  Double_t naoq = AOQ[aoq];
  
  for(short iplane = 0; iplane < 3; iplane++){
    for(short ivar = 0; ivar < 3; ivar++){
      if( iplane == var/3 && ivar ==  var%3) continue;
      naoq = naoq - (fvar[iplane][ivar]+fcorr_Lin_3[iplane][ivar][0])*fcorr_Lin_3[iplane][ivar][1];
      /*naoq = naoq -( (fvar[iplane][ivar] + fcorr_Lin_3[iplane][ivar][0]) *
		     (fvar[iplane][ivar] + fcorr_Lin_3[iplane][ivar][0]) *
		     fcorr_Lin_3[iplane][ivar][1]                        +
		     fcorr_Lin_3[iplane][ivar][2]);
      */
    }//for(iplane)
  }//for(ivar)
  
  if( var != 9 )
    naoq = naoq - (BETA[3]*fcorr_Lin_3[4][0][0])*fcorr_Lin_3[4][0][1];


  return naoq;

}//Correction_3

void PID_2::hitsperevent(Long64_t entry){

  for(short ii = 0; ii <kMaxBigRIPSPlastic; ii++)
    if(BigRIPSPlastic_fpl[ii] == 7)
      fh2_timevs_evnumber->Fill(entry,BigRIPSPlastic_fTime[ii]);



}

Double_t PID_2::zet(short aoq, short var){

  //  TF1* Linfunc = new TF1("Linfunc",Form("[1]*(x+[0])+[2]+%f", 50),-10,10);
  Double_t nzet = ZET[aoq];
  Double_t bla = 0;
  short kk =0;
  if(aoq==2)aoq=1;
  for(short jj = 1; jj < 6; jj=jj+2){//even first plane odd, second plane (F7 / F11)
    //cout<< fzet_corr_Lin[aoq][kk][0] << " " << fzet_corr_Lin[aoq][kk][1] <<endl;
    if( jj == var ) continue;
    else if(aoq > 0)  bla = fout_var[1][jj];//aoq==2 -1=1 -- last plane before MUSIC
    else              bla = fin_var[jj];
    nzet = nzet - fzet_corr_Lin[aoq][kk][1]*(bla+fzet_corr_Lin[aoq][kk][0]);//* (bla+fzet_corr_Lin[aoq][jj][0]);
    kk++;
  }
  if( var != 6 )
    nzet = nzet - fzet_corr_Lin[aoq][3][1]*(BETA[aoq]+fzet_corr_Lin[aoq][3][0]);//*(BETA[aoq]+fzet_corr_Lin[aoq][3][0]);
  //  cout<< fzet_corr_Lin[aoq][3][0] << " " << fzet_corr_Lin[aoq][3][1] <<endl;
  
  return nzet;

}//zet

void PID_2::reset(){//reset all the variables after each event

  for(short ii = 0; ii<7;ii++){

    fin_var[ii]=-100;
    for(short jj =0; jj<2; jj++){
      fout_var[jj][ii]=-100;
      fout_var[jj][ii]=-100;
    }
  }
  for(short ii = 0;ii<3;ii++)
    for(short jj = 0;jj<3;jj++)
      fvar[ii][jj] = -100;

  }
