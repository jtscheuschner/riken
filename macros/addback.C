

Float_t dopp_Corr(Float_t beta, Float_t cosTheta, Float_t energy){

  Float_t gamma = 1./TMath::Sqrt(1-beta*beta);
  Float_t doppcorr = energy * (1-beta*cosTheta)*gamma;

  return doppcorr;

}

Float_t distance(Float_t xpos1, Float_t ypos1, Float_t zpos1,Float_t xpos2, Float_t ypos2, Float_t zpos2){

  Float_t dist = (xpos1-xpos2)*(xpos1-xpos2)+
    (ypos1-ypos2)*(ypos1-ypos2)+
    (zpos1-zpos2)*(zpos1-zpos2);

  return TMath::Sqrt(dist);

}


Float_t addback(){


}
