//#include "macros/PID.h"

//gROOT->ProcessLine("macros/PID.C+");

void run(short maxIt, short var){

  PID pid(0);




  pid.Loop(maxIt,var);

  return;

  cout << " hi" <<endl;
  const UInt_t VAR = 9;

  for(short it = 0; it < maxIt; it++){
    for(short ivar = 0; ivar < VAR; ivar++){
      pid.Loop(it,ivar);
      //pid.Fitting(pid.fh2_aoq_0,0,ivar);
      //pid.Fitting(pid.fh2_aoq_1,1,ivar);
      //pid.Fitting(pid.fh2_aoq_2,2,ivar);
      //pid.Fitting(pid.fh2_aoq_3,3,ivar);
      if(it < maxIt -1 || ivar < VAR-1)pid.Reset();
    }//for(ivar)


  }//for(it)
  
  pid.AoQ();
  pid.Write("test.root");
  

}//run
