#ifndef readROOT_H
#define readROOT_H 1

#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TVector3.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TTreeReader.h>
#include <TLorentzVector.h>
#include <TTreeReaderArray.h>

struct vec4{
  TLorentzVector particle4Vec;
  double vertX,vertY,vertZ;
  double vertR;
};
class readROOT{
public:
  readROOT();
  virtual ~readROOT();

  void beginROOT();
  void readFile(TChain*);
  void tracking();
  void writeROOT();

private:
  double r;
  // Create/Write a ROOT file
  TFile* ofile;

  TH1D *partEta;
  TH1D *matchedPartEta;
  TH1D* partMom;
  TH1D* matchedPartMom;
  TH1D* partPhi;
  TH1D* matchedPartPhi;

  TH2D* partPEta;
  TH2D* matchedPartPEta;
  TH2D* partPhiEta;
  TH2D* matchedPartPhiEta;

  TH1D *matchedPartTrackDeltaEta;
  TH1D *matchedPartTrackDeltaPhi;
  TH1D *matchedPartTrackDeltaR;
  TH1D *matchedPartTrackDeltaMom;

  // Define some histograms for our efficiencies
  TH1D *TrackEff_Eta;
  TH1D *TrackEff_Mom;
  TH1D *TrackEff_Phi;

  // 2D Efficiencies
  TH2D* TrackEff_PEta;
  TH2D* TrackEff_PhiEta;

  // All charged particle histos
  TH1D *h1posChargedEta[4];
  TH1D *h1posChargedTheta[4];
  TH1D *h1posChargedPhi[4];
  TH1D *h1posChargedP[4];
  TH1D *h1pospathLen[4];
  TH1D *h1posvertBegin[4];
  TH1D *h1posvertEnd[4];
  //negative
  TH1D *h1negChargedEta[4];
  TH1D *h1negChargedTheta[4];
  TH1D *h1negChargedPhi[4];
  TH1D *h1negChargedP[4];
  TH1D *h1negpathLen[4];
  TH1D *h1negvertBegin[4];
  TH1D *h1negvertEnd[4];

  // TRACKING DETECTOR HISTOS
  TH2D* h2VertBHits[4],*h2SiBHitsXY[4],*h2MPGDBHits[4];

  //Reco Heavy Meson
  TH1D *h1InvM[4];
  TH1D *h1VertX[4];
  TH1D *h1VertY[4];
  TH1D *h1VertZ[4];
  TH1D *h1VertR[4];

  //2-D plots
  TH2D* h2pTvy[4];
  TH2D* h2pvEta[4];
  TH2D* h2pTvalpha[4];
  TH2D* h2vertXY[4];
  TH2D *h2VertREta[4];
  // daughter particles p v. eta
  TH2D* h2pospvEta[4];
  TH2D* h2negpvEta[4];

  //4-vector of particles
  TLorentzVector M_D0[4], M_neg[4], M_pos[4];
  TVector3 prtlpos[4],prtlneg[4];
  std::vector<vec4> posPrtl,negPrtl;
  std::vector<vec4> kaon4Vec,pion4Vec;
  vec4 parentLV;
};
#endif
