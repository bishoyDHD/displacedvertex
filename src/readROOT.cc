#include "readROOT.h"

readROOT::readROOT(){

}
readROOT::~readROOT(){}

void readROOT::beginROOT(){
  // Set output file for the histograms
  ofile = new TFile("../outputs/anafiles/minQ2_1D0recoAnalysis.root","RECREATE");
}
void readROOT::writeROOT(){
  ofile->Write(); // Write histograms to file
  ofile->Close(); // Close output file
}
void readROOT::readFile(TChain* pchain){
  // Initialize reader
  TTreeReader tree_reader(pchain);

  // Get Particle Information
  TTreeReaderArray<int> partGenStat(tree_reader, "MCParticles.generatorStatus");
  TTreeReaderArray<float> simMomX(tree_reader, "MCParticles.momentum.x");
  TTreeReaderArray<float> simMomY(tree_reader, "MCParticles.momentum.y");
  TTreeReaderArray<float> simMomZ(tree_reader, "MCParticles.momentum.z");
  TTreeReaderArray<int> partPdg(tree_reader, "MCParticles.PDG");

  // Get Central Track Vertex (ctv) Information
  TTreeReaderArray<float> ctvChi2(tree_reader, "CentralTrackVertices.chi2");
  TTreeReaderArray<int> ctvNDF(tree_reader, "CentralTrackVertices.ndf");
  TTreeReaderArray<float> ctvErrX(tree_reader, "CentralTrackVertices.positionError.xx");
  TTreeReaderArray<float> ctvErrY(tree_reader, "CentralTrackVertices.positionError.yy");
  TTreeReaderArray<float> ctvErrZ(tree_reader, "CentralTrackVertices.positionError.zz");
  TTreeReaderArray<float> ctvPosX(tree_reader, "CentralTrackVertices.position.x");
  TTreeReaderArray<float> ctvPosY(tree_reader, "CentralTrackVertices.position.y");
  TTreeReaderArray<float> ctvPosZ(tree_reader, "CentralTrackVertices.position.z");
  // Get Central Track Segment Points Information
  TTreeReaderArray<unsigned long> ctsegSurface(tree_reader, "_CentralTrackSegments_points.surface");
  TTreeReaderArray<float> ctsegTheta(tree_reader, "_CentralTrackSegments_points.theta");
  TTreeReaderArray<float> ctsegPhi(tree_reader, "_CentralTrackSegments_points.phi");
  TTreeReaderArray<float> ctsegTime(tree_reader, "_CentralTrackSegments_points.time");
  TTreeReaderArray<float> ctsegLen(tree_reader, "_CentralTrackSegments_points.pathlength");
  TTreeReaderArray<float> ctsegMomX(tree_reader, "_CentralTrackSegments_points.momentum.x");
  TTreeReaderArray<float> ctsegMomY(tree_reader, "_CentralTrackSegments_points.momentum.y");
  TTreeReaderArray<float> ctsegMomZ(tree_reader, "_CentralTrackSegments_points.momentum.z");
  TTreeReaderArray<float> ctsegPosX(tree_reader, "_CentralTrackSegments_points.position.x");
  TTreeReaderArray<float> ctsegPosY(tree_reader, "_CentralTrackSegments_points.position.y");
  TTreeReaderArray<float> ctsegPosZ(tree_reader, "_CentralTrackSegments_points.position.z");

  // Get Central CKF Trajectory Information
  TTreeReaderArray<unsigned int> ckftrajErrX(tree_reader, "CentralCKFTrajectories.nStates");
  TTreeReaderArray<unsigned int> ckftrajErrY(tree_reader, "CentralCKFTrajectories.nMeasurements");
  TTreeReaderArray<unsigned int> ckftrajErrZ(tree_reader, "CentralCKFTrajectories.nOutliers");
  TTreeReaderArray<unsigned int> ckftrajPosX(tree_reader, "CentralCKFTrajectories.nHoles");
  TTreeReaderArray<unsigned int> ckftrajPosY(tree_reader, "CentralCKFTrajectories.nSharedHits");
  // Get Central CKF Track Parameter Information
  TTreeReaderArray<float> ckftparamLoca(tree_reader, "CentralCKFTrackParameters.loc.a");
  TTreeReaderArray<float> ckftparamLocb(tree_reader, "CentralCKFTrackParameters.loc.b");
  TTreeReaderArray<float> ckftparamTheta(tree_reader, "CentralCKFTrackParameters.theta");
  TTreeReaderArray<float> ckftparamPhi(tree_reader, "CentralCKFTrackParameters.phi");
  TTreeReaderArray<float> ckftparamqOverP(tree_reader, "CentralCKFTrackParameters.qOverP");

  // Get Reconstructed Track Information
  TTreeReaderArray<float> trackE(tree_reader, "ReconstructedChargedParticles.energy");
  TTreeReaderArray<float> trackCharge(tree_reader, "ReconstructedChargedParticles.charge");
  TTreeReaderArray<int> trackPDG(tree_reader, "ReconstructedChargedParticles.PDG");
  TTreeReaderArray<float> trackMomX(tree_reader, "ReconstructedChargedParticles.momentum.x");
  TTreeReaderArray<float> trackMomY(tree_reader, "ReconstructedChargedParticles.momentum.y");
  TTreeReaderArray<float> trackMomZ(tree_reader, "ReconstructedChargedParticles.momentum.z");
  TTreeReaderArray<float> trackRefX(tree_reader, "ReconstructedChargedParticles.referencePoint.x");
  TTreeReaderArray<float> trackRefY(tree_reader, "ReconstructedChargedParticles.referencePoint.y");
  TTreeReaderArray<float> trackRefZ(tree_reader, "ReconstructedChargedParticles.referencePoint.z");
  TTreeReaderArray<unsigned int> trackBegin(tree_reader, "ReconstructedChargedParticles.tracks_begin");
  TTreeReaderArray<unsigned int> trackEnd(tree_reader, "ReconstructedChargedParticles.tracks_end");

  // Get Associations Between MCParticles and ReconstructedChargedParticles
  TTreeReaderArray<unsigned int> recoAssoc(tree_reader, "ReconstructedChargedParticleAssociations.recID");
  TTreeReaderArray<unsigned int> simuAssoc(tree_reader, "ReconstructedChargedParticleAssociations.simID");
    
  // Define Histograms
  TH1D* partEta = new TH1D("partEta","#eta of Thrown Charged Particles; #eta", 120, -6, 6);
  TH1D* matchedPartEta = new TH1D("matchedPartEta","#eta of Thrown Charged Particles That Have Matching Track; #eta", 120, -6, 6);
  TH1D* partMom = new TH1D("partMom", "Momentum of Thrown Charged Particles (truth); P(GeV/c)", 150, 0, 150);
  TH1D* matchedPartMom = new TH1D("matchedPartMom", "Momentum of Thrown Charged Particles (truth), with matching track; P(GeV/c)", 150, 0, 150);
  TH1D* partPhi = new TH1D("partPhi", "#phi of Thrown Charged Particles (truth); #phi(rad)", 320, -3.2, 3.2);
  TH1D* matchedPartPhi = new TH1D("matchedPartPhi", "#phi of Thrown Charged Particles (truth), with matching track; #phi(rad)", 320, -3.2, 3.2);

  TH2D* partPEta = new TH2D("partPEta", "P vs #eta of Thrown Charged Particles; P(GeV/c); #eta", 150, 0, 150, 120, -6, 6);
  TH2D* matchedPartPEta = new TH2D("matchedPartPEta", "P vs #eta of Thrown Charged Particles, with matching track; P(GeV/c); #eta", 150, 0, 150, 120, -6, 6);
  TH2D* partPhiEta = new TH2D("partPhiEta", "#phi vs #eta of Thrown Charged Particles; #phi(rad); #eta", 160, -3.2, 3.2, 120, -6, 6);
  TH2D* matchedPartPhiEta = new TH2D("matchedPartPhiEta", "#phi vs #eta of Thrown Charged Particles; #phi(rad); #eta", 160, -3.2, 3.2, 120, -6, 6);
    
  TH1D *matchedPartTrackDeltaEta = new TH1D("matchedPartTrackDeltaEta","#Delta#eta Between Matching Thrown and Reconstructed Charged Particle; #Delta#eta", 100, -0.25, 0.25);
  TH1D *matchedPartTrackDeltaPhi = new TH1D("matchedPartTrackDeltaPhi","#Detla #phi Between Matching Thrown and Reconstructed Charged Particle; #Delta#phi", 200, -0.2, 0.2);
  TH1D *matchedPartTrackDeltaR = new TH1D("matchedPartTrackDeltaR","#Delta R Between Matching Thrown and Reconstructed Charged Particle; #Delta R", 300, 0, 0.3);
  TH1D *matchedPartTrackDeltaMom = new TH1D("matchedPartTrackDeltaMom","#Delta P Between Matching Thrown and Reconstructed Charged Particle; #Delta P", 200, -10, 10);
    
  // Define some histograms for our efficiencies
  TH1D *TrackEff_Eta = new TH1D("TrackEff_Eta", "Tracking efficiency as fn of #eta; #eta; Eff(%)", 120, -6, 6); 
  TH1D *TrackEff_Mom = new TH1D("TrackEff_Mom", "Tracking efficiency as fn of P; P(GeV/c); Eff(%)", 150, 0, 150); 
  TH1D *TrackEff_Phi = new TH1D("TrackEff_Phi", "Tracking efficiency as fn of #phi; #phi(rad); Eff(%)", 320, -3.2, 3.2);

  // 2D Efficiencies
  TH2D* TrackEff_PEta = new TH2D("TrackEff_PEta", "Tracking efficiency as fn of P and #eta; P(GeV/c); #eta", 150, 0, 150, 120, -6, 6);
  TH2D* TrackEff_PhiEta = new TH2D("TrackEff_PhiEta", "Tracking efficiency as fn of #phi and #eta; #phi(rad); #eta", 160, -3.2, 3.2, 120, -6, 6);

  // All charged particle histos
  std::string nameHist;
  for(int n=0; n<4; n++){
    // Positively charged particle
    nameHist="h1posChargedEta_"+std::to_string(n);
    h1posChargedEta[n] = new TH1D(nameHist.c_str(), "#eta of positive charged particles; #eta", 120, -6, 6);
    nameHist="h1posChargedPhi_"+std::to_string(n);
    h1posChargedPhi[n] = new TH1D(nameHist.c_str(), "#phi of positive charged particles; #phi (rad)", 120, -3.2, 3.2);
    nameHist="h1posChargedTheta_"+std::to_string(n);
    h1posChargedTheta[n] = new TH1D(nameHist.c_str(), "#theta of positive charged particles; #theta (rad)", 120, -3.2, 3.2);
    nameHist="h1posChargedP_"+std::to_string(n);
    h1posChargedP[n] = new TH1D(nameHist.c_str(), "P of positive charged particles; P(GeV/c)", 150, 0, 150);
    // Negatively charged particle
    nameHist="h1negChargedEta_"+std::to_string(n);
    h1negChargedEta[n] = new TH1D(nameHist.c_str(), "#eta of negavtive charged particles; #eta", 120, -6, 6);
    nameHist="h1negChargedPhi_"+std::to_string(n);
    h1negChargedPhi[n] = new TH1D(nameHist.c_str(), "#phi of negative charged particles; #phi (rad)", 120, -3.2, 3.2);
    nameHist="h1negChargedTheta_"+std::to_string(n);
    h1negChargedTheta[n] = new TH1D(nameHist.c_str(), "#theta of negative charged particles; #theta (rad)", 120, -3.2, 3.2);
    nameHist="h1negChargedP_"+std::to_string(n);
    h1negChargedP[n] = new TH1D(nameHist.c_str(), "P of negative charged particles; P(GeV/c)", 150, 0, 150);
    nameHist="h1InvM_"+std::to_string(n);
    h1InvM[n] = new TH1D(nameHist.c_str(), "Invariant Mass GeV/c^{2}; M[GeV/c^{2}]", 150, 0, 15);
    nameHist="h1VertX_"+std::to_string(n);
    h1VertX[n] = new TH1D(nameHist.c_str(), "Reconstructed Vertex X; x[mm]", 100, -5, 5);
    nameHist="h1VertY_"+std::to_string(n);
    h1VertY[n] = new TH1D(nameHist.c_str(), "Reconstructed Vertex Y; y[mm]", 100, -5, 5);
    nameHist="h1VertZ_"+std::to_string(n);
    h1VertZ[n] = new TH1D(nameHist.c_str(), "Reconstructed Vertex Z; z[mm]", 100, -150, 150);
    nameHist="h1VertR_"+std::to_string(n);
    h1VertR[n] = new TH1D(nameHist.c_str(), "Reconstructed Vertex R; r[mm]", 300, 0, 150);

    //2-D profile Histos
    nameHist="h2vertXY_"+std::to_string(n);
    h2vertXY[n]=new TH2D(nameHist.c_str(),"",100,-5,5,100,-5,5);
    nameHist="h2pTvy_"+std::to_string(n);
    h2pTvy[n]=new TH2D(nameHist.c_str(),"",60,-2,4,100,0,10);
    nameHist="h2pvEta_"+std::to_string(n);
    h2pvEta[n]=new TH2D(nameHist.c_str(),"",60,-2,4,100,0,20);
    //Daughter p v. Eta profiles
    nameHist="h2pospvEta_"+std::to_string(n);
    h2pospvEta[n]=new TH2D(nameHist.c_str(),"",60,-2,4,100,0,20);
    nameHist="h2negpvEta_"+std::to_string(n);
    h2negpvEta[n]=new TH2D(nameHist.c_str(),"",60,-2,4,100,0,20);
  }
  
  while(tree_reader.Next()) { // Loop over events
    // Empty the container at the start of each event
    posPrtl.clear(); negPrtl.clear();
    kaon4Vec.clear(); pion4Vec.clear();
    for(unsigned int j=0; j<trackCharge.GetSize(); j++){ // Loop over charged particles
      M_D0[0].SetPxPyPzE(trackMomX[j],trackMomY[j],trackMomZ[j],trackE[j]);
      h1InvM[0]->Fill(M_D0[0].M());
      h1vertX[0]->Fill(ctvPosX[j]);
      h1vertY[0]->Fill(ctvPosY[j]);
      h1vertZ[0]->Fill(ctvPosZ[j]);
      r=std::sqrt(std::pow(ctvPosX[j],2)+std::pow(ctvPosY[j],2)+std::pow(ctvPosZ[j],2));
      h2pTvy[0]->Fill(M_D0[0].Rapidity(),M_D0[0].Pt());
      h2pvEta[0]->Fill(M_D0[0].PseudoRapidity(),M_D0[0].P());
      if(trackCharge[j]==1){//positively charged particle
        parentLV.particle4Vec.SetPxPyPzE(trackMomX[j],trackMomY[j],trackMomZ[j],trackE[j]);
        posPrtl.push_back(parentLV);
        h2pospvEta[0]->Fill(parentLV.particle4Vec.Rapidity(),parentLV.particle4Vec.P());
        if(trackPDG[j]==321){//K+ 4-vector
          parentLV.particle4Vec.SetPxPyPzE(trackMomX[j],trackMomY[j],trackMomZ[j],trackE[j]);
          kaon4Vec.push_back(parentLV);
          h2pospvEta[1]->Fill(parentLV.particle4Vec.Rapidity(),parentLV.particle4Vec.P());
        }//end K+
      } //end positively charged particle if-loop
      if(trackCharge[j]==-1){//negatively charged particle
        parentLV.particle4Vec.SetPxPyPzE(trackMomX[j],trackMomY[j],trackMomZ[j],trackE[j]);
        negPrtl.push_back(parentLV);
        h2negpvEta[0]->Fill(parentLV.particle4Vec.Rapidity(),parentLV.particle4Vec.P());
        if(trackPDG[j]==-211){//pi- 4-vector
          parentLV.particle4Vec.SetPxPyPzE(trackMomX[j],trackMomY[j],trackMomZ[j],trackE[j]);
          pion4Vec.push_back(parentLV);
          h2negpvEta[1]->Fill(parentLV.particle4Vec.Rapidity(),parentLV.particle4Vec.P());
        }//end pi-
      } //end negatively charged particle if-loop
    } // End loop over all charged particles
    //Check that neither pos nor neg particle vectors are empty
    if(posPrtl.size()>0 && negPrtl.size()>0){
      for(int p=0; p<posPrtl.size(); p++){
        for(int n=0; n<negPrtl.size(); n++){
          M_D0[1]=posPrtl[p].particle4Vec+negPrtl[n].particle4Vec;
          h1InvM[1]->Fill(M_D0[1].M());
          h2pTvy[1]->Fill(M_D0[1].Rapidity(),M_D0[1].Pt());
          h2pvEta[1]->Fill(M_D0[1].PseudoRapidity(),M_D0[1].P());
        }
      }
    } //end of charged particle 4-vector sum
    //Check that neither kaon nor pion particle vectors are empty
    if(kaon4Vec.size()>0 && pion4Vec.size()>0){
      for(int p=0; p<kaon4Vec.size(); p++){
        h2pospvEta[2]->Fill(kaon4Vec[p].particle4Vec.Rapidity(),kaon4Vec[p].particle4Vec.P());
        for(int n=0; n<pion4Vec.size(); n++){
          M_D0[2]=kaon4Vec[p].particle4Vec+pion4Vec[n].particle4Vec;
          h1InvM[2]->Fill(M_D0[2].M());
          h2pTvy[2]->Fill(M_D0[2].Rapidity(),M_D0[2].Pt());
          h2pvEta[2]->Fill(M_D0[2].PseudoRapidity(),M_D0[2].P());
        }
      }
    } //end of invarian k+ + pi- particle 4-vector sum
/*
	if(partGenStat[i] == 1){ // Select stable thrown particles
	    int pdg = TMath::Abs(partPdg[i]);
	    if(pdg == 11 || pdg == 13 || pdg == 211 || pdg == 321 || pdg == 2212){ // Look at charged particles (electrons, muons, pions, kaons, protons)
		TVector3 trueMom(partMomX[i],partMomY[i],partMomZ[i]);

		float trueEta = trueMom.PseudoRapidity();
		float truePhi = trueMom.Phi();
	    
		partEta->Fill(trueEta);
		partPhi->Fill(truePhi);
		partMom->Fill(trueMom.Mag());
		partPEta->Fill(trueMom.Mag(), trueEta);
		partPhiEta->Fill(truePhi, trueEta);

		// Loop over associations to find matching ReconstructedChargedParticle
		for(unsigned int j=0; j<simuAssoc.GetSize(); j++){
		    if(simuAssoc[j] == i){ // Find association index matching the index of the thrown particle we are looking at
			TVector3 recMom(trackMomX[recoAssoc[j]],trackMomY[recoAssoc[j]],trackMomZ[recoAssoc[j]]); // recoAssoc[j] is the index of the matched ReconstructedChargedParticle

			// Check the distance between the thrown and reconstructed particle
			float deltaEta = trueEta - recMom.PseudoRapidity();
			float deltaPhi = TVector2::Phi_mpi_pi(truePhi - recMom.Phi());
			float deltaR = TMath::Sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);
			float deltaMom = ((trueMom.Mag()) - (recMom.Mag()));

			matchedPartTrackDeltaEta->Fill(deltaEta);
			matchedPartTrackDeltaPhi->Fill(deltaPhi);
			matchedPartTrackDeltaR->Fill(deltaR);
			matchedPartTrackDeltaMom->Fill(deltaMom);

			matchedPartEta->Fill(trueEta); // Plot the thrown eta if a matched ReconstructedChargedParticle was found
			matchedPartPhi->Fill(truePhi);
			matchedPartMom->Fill(trueMom.Mag());

			matchedPartPEta->Fill(trueMom.Mag(), trueEta);
			matchedPartPhiEta->Fill(truePhi, trueEta);
		      }
		  }// End loop over associations
	      } // End PDG check
	  } // End stable particles condition
      } // End loop over thrown particles
    // Loop over all charged particles and fill some histograms of kinematics quantities
    for(unsigned int k=0; k<trackMomX.GetSize(); k++){ // Loop over all charged particles, thrown or not
      
      TVector3 CPartMom(trackMomX[k], trackMomY[k], trackMomZ[k]);

      float CPartEta = CPartMom.PseudoRapidity();
      float CPartPhi = CPartMom.Phi();
      ChargedEta->Fill(CPartEta);
      ChargedPhi->Fill(CPartPhi);
      ChargedP->Fill(CPartMom.Mag());
*/      
  } // End loop over events

/*
  // Take the ratio of the histograms above to get our efficiency plots
  TrackEff_Eta->Divide(matchedPartEta, partEta, 1, 1, "b");
  TrackEff_Mom->Divide(matchedPartMom, partMom, 1, 1, "b");
  TrackEff_Phi->Divide(matchedPartPhi, partPhi, 1, 1, "b");
  TrackEff_PEta->Divide(matchedPartPEta, partPEta, 1, 1, "b");
  TrackEff_PhiEta->Divide(matchedPartPhiEta, partPhiEta, 1, 1, "b");
  
*/
}
