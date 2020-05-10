#include "BSHelper.cxx"
//
TFile * LoadJetResults( TString name, TString runnum){
  if( !runnum.IsNull() ) runnum = "_"+runnum;
  //auto fname = "../Results/"+name+"/AnalysisResults_"+name+runnum+".root";
  auto fname = "../results/"+name+"/AnalysisResults_"+name+runnum+".root";
	cout<<fname<<endl;
  return  LoadRoot( fname ) ;
}
//__________________________________________________________
TObject*  LoadJetResultList( TFile *fh, TString clistname){
  auto l = fh->Get( clistname );
  if(!l) ErrorExit("No list "+clistname);
  return l;
}


//__________________________________________________________
TObject*  LoadJetResultList( TString fname, TString clistname, TString runnum=""){
  auto f = LoadJetResults( fname, runnum );
  return LoadJetResultList( f, clistname );
}

