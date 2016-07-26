void bla(Char_t* aoq, Char_t* var, Int_t nbin = 500, Float_t xmin = -10., Float_t xmax = 10. ){

  TChain* tree = new TChain("tree");
  tree->AddFile("data/rootfiles/test/run0009_299.793_-158.540.root");
  //TFile *file = TFile::Open("data/rootfiles/test/run0009_299.793_-158.540.root");
  TH2F* h2 = new TH2F("h2","h2",nbin,xmin,xmax,300,2.5,2.7);
  Char_t* cut = "TMath::Abs(AOQ[0]-2.64)<0.008 && TMath::Abs(ZET[0]-50)<0.4 && TMath::Abs(delta[0]-delta[2])<0.5 && TMath::Abs(AOQ[3]-2.64)<0.01 && TMath::Abs(ZET[2]-49.5)<0.6";
  tree->Draw(Form("%s:%s>>h2",aoq,var),cut,"colz");
  TObjArray* dummy = new TObjArray();
  dummy->SetOwner(kTRUE);
  h2->FitSlicesY(0,0,-1,0,"QNRG4",dummy);
  TH1F* h1 = (TH1F*)dummy->FindObject("h2_1");
  h1->Draw("same");
  TF1* f1 = new TF1("f1","[0]+[1]*x",xmin,xmax);
  f1->SetLineColor(kBlue);
  h1->Fit(f1,"","",xmin,xmax);

}
