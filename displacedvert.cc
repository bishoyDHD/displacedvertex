#include "readROOT.h"

int main(int argc, char* argv[]){
  readROOT* anaroot=new readROOT();
  TString infile="PATH_TO_FILE";
  // Set output file for the histograms
  TChain *mychain = new TChain("events");
  int i=1;
  while(i<argc){
    infile=argv[i];
    mychain->Add(infile);
    i++;
  }
  anaroot->beginROOT();
  anaroot->readFile(mychain);
  anaroot->writeROOT();

  delete anaroot;

  return 0;
}
