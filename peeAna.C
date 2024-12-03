#include "commonUtility.h"
#include <TChain.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TMath.h>
#include <TGraph2D.h>
#include <TRandom2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TF2.h>
#include <TH1.h>
#include <TGraphErrors.h>
#include <TMarker.h>
#include <Math/Functor.h>
#include <TPolyLine3D.h>
#include <Math/Vector3D.h>
#include <Fit/Fitter.h>
#include <cassert>
#include <TMinuit.h>
#include <TH3D.h>
#include <iostream>
#include <vector>
#include <TVector3.h>
#include <TPolyLine3D.h>
#include <TLorentzVector.h>


std::vector<double> xData, yData, zData, eData;

int displayCount=0;

int maxEvents = 1000;

float protonMass = 9.3827e-01;
float eleMass =  0.000511;

//void test(TString infile="reconstructed_data_N2100_neutron_theta_0_0.2mard_25GeV_OnlyHcal_info.root") {
//void peeAna(TString infile="podio_output_Pentaquark_hepmc_output_20241112_50GeV_10evts.hepmc.edm4hep.root") {
// void peeAna(TString infile="podio_output_Pentaquark_hepmc_output_20241113_50GeV_10000evts.hepmc.edm4hep.root") { 1k E = 50 GeV
void peeAna(TString infile="podio_files/podio_Pentaquark_hepmc_output_20241126E100-200GeV_eta1.9-7_100000evts_ip6_hidiv_275x18.root") {
  // void peeAna(TString infile="podio_pentaquark_6_newSource_2024.12.01.root") {
  
    const int kElse = 0;
    const int kProton = 1;
    const int kElectron = 2;
    const int kPositron = 3;
  
  // p bins for resolution study.  Must be four bins
    double pBin[5] = {0,20,50,100,200};
    TH1D* pBinHist = new TH1D("hpBinHist"," ;p (GeV);",4, pBin);
    
    bool kEvtDis = 1 ;
  
  // Set output file for the histograms
  // TFile *ofile = TFile::Open(Form("output_%s", infile.Data()),"RECREATE");
  
  // Set up input file chain
    TChain *mychain = new TChain("events");
    mychain->Add(infile);
  
  // Initialize reader
    TTreeReader tree_reader(mychain);
  
  // Get Particle Information
  //  TTreeReaderArray<float> schit_e(tree_reader, "HcalFarForwardZDCSubcellHits.energy");
  //  TTreeReaderArray<float> schit_x(tree_reader, "HcalFarForwardZDCSubcellHits.position.x");
  //  TTreeReaderArray<float> schit_y(tree_reader, "HcalFarForwardZDCSubcellHits.position.y");
  //  TTreeReaderArray<float> schit_z(tree_reader, "HcalFarForwardZDCSubcellHits.position.z");
  
    TTreeReaderArray<int> genp_pdg(tree_reader, "MCParticles.PDG");
    TTreeReaderArray<int> genp_status(tree_reader, "MCParticles.generatorStatus");
    TTreeReaderArray<float> genp_px(tree_reader, "MCParticles.momentum.x");
    TTreeReaderArray<float> genp_py(tree_reader, "MCParticles.momentum.y");
    TTreeReaderArray<float> genp_pz(tree_reader, "MCParticles.momentum.z");
  
    TTreeReaderArray<int>  recop_simID(tree_reader, "ReconstructedTruthSeededChargedParticleAssociations.simID");
    TTreeReaderArray<float> recop_px(tree_reader, "ReconstructedTruthSeededChargedParticles.momentum.x");
    TTreeReaderArray<float> recop_py(tree_reader, "ReconstructedTruthSeededChargedParticles.momentum.y");
    TTreeReaderArray<float> recop_pz(tree_reader, "ReconstructedTruthSeededChargedParticles.momentum.z");
  
    /*  TTreeReaderArray<int>  recop_simID(tree_reader, "ReconstructedParticleAssociations.simID");
    TTreeReaderArray<float> recop_px(tree_reader, "ReconstructedParticles.momentum.x");
    TTreeReaderArray<float> recop_py(tree_reader, "ReconstructedParticles.momentum.y");
    TTreeReaderArray<float> recop_pz(tree_reader, "ReconstructedParticles.momentum.z");
    */
    
    
    double schit_err[1500];
    double dbschit_x[1500];
    double dbschit_y[1500];
    double dbschit_z[1500];
  
  //  TTreeReaderArray<ULong_t> schit_subcellid(tree_reader, "HcalFarForwardZDCSubcellHits.cellID");
  
  // RECO particle histograms
  
    TH1D *hRecop_pt = new TH1D("recop_pt",";p_{T} (GeV)",100,0,30);
  
    TH1D *hRecop_phi = new TH1D ("recop_phi",";#phi",100,3.141592, -3.141592);
  
    TH1D *hRecop_p = new TH1D("recop_p",";p^{RECO} (GeV)",100,0,200);
    TH1D *hRecop_p_match[4];
    TH1D *hRecop_eta = new TH1D("recop_eta",";#eta^{RECO}",100,0,6);
    TH1D *hRecop_eta_match[4];
    TH1D *hMatch_dp[4][4];  // [PID][p bins]
    TH1D* hResP[4]; // PID
    for ( int i = 0 ; i<4 ; i++)  { // kElse = 0; kProton = 1; kElectron = 2; kPositron = 3
        hRecop_p_match[i] = (TH1D*)hRecop_p->Clone(Form("hRecop_p_match%d",i));
        hRecop_eta_match[i] = (TH1D*)hRecop_eta->Clone(Form("hRecop_eta_match%d",i));
        for ( int jp =0 ; jp<4; jp++) {
            hMatch_dp[i][jp] = new TH1D(Form("hmatch_dp_%di_%djp",i,jp), ";dp/p;",100,-0.5,0.5);
        }
        hResP[i] = (TH1D*)pBinHist->Clone(Form("hRes_d%d",i));
    }
  
  
  
  
  
  // Gen Particle histograms
    TH1D *hGen_p[4];
    TH1D *hGen_p_recoMatched[4];
    TH1D *hGen_eta[4];
    TH1D *hGen_eta_recoMatched[4];
    TH1D *hGen_phi[4];
    TH1D *hGen_phi_recoMatched[4];
  
    for ( int id=0 ; id<=3 ; id++) {
        hGen_p[id] = new TH1D(Form("hGen_p_id%d",id),"; p^{GEN} (GeV);",100,0,300);
        hGen_p_recoMatched[id] = (TH1D*)hGen_p[id]->Clone(Form("hGen_p_recoMatched_id%d",id));
        hGen_eta[id] = new TH1D(Form("hGen_eta_id%d",id),"; #eta;",100,0,7);
        hGen_eta_recoMatched[id] = (TH1D*)hGen_eta[id]->Clone(Form("hGen_eta_recoMatched_id%d",id));
        hGen_phi[id] = new TH1D(Form("hGen_phi_id%d",id),"; #phi;",100,-3.2,3.2);
        hGen_phi_recoMatched[id] = (TH1D*)hGen_phi[id]->Clone(Form("hGen_phi_recoMatched_id%d",id));
    }
  
    TH1D *hProton_pt = new TH1D("hProton_pt",";p_{T} (GeV)",100,0,30);
    TH1D *hProton_eta = new TH1D ("hProton_eta",";#eta",100,0,7);
    TH1D *hProton_phi = new TH1D ("hProton_phi",";#phi",100,3.141592, -3.141592);
  
    TH1D *hProtonReco_pt = (TH1D*)hProton_pt->Clone("hProtonReco_pt");
    TH1D *hProtonReco_eta = (TH1D*)hProton_eta->Clone("hProtonReco_eta");
    TH2D *hProtonReco_EtaPt = new TH2D("hProtonReco_EtaPt",";#eta;pT (GeV)",100,0,7,100,0,30);
  
    TH1D *hEle_pt = (TH1D*)hProton_pt->Clone("hEle_pt");
    TH1D *hEle_eta = (TH1D*)hProton_eta->Clone("hEle_eta");
    TH1D *hEle_phi = (TH1D*)hProton_phi->Clone("hEle_phi");
    
    TH1D *hPos_pt = (TH1D*)hProton_pt->Clone("hPos_pt");
    TH1D *hPos_eta = (TH1D*)hProton_eta->Clone("hPos_eta");
    TH1D *hPos_phi = (TH1D*)hProton_phi->Clone("hPos_phi");
  
    TH1D *hProton_pt_match = (TH1D*)hProton_pt->Clone("hProton_pt_match");
    TH1D *hProton_eta_match = (TH1D*)hProton_eta->Clone("hProton_eta_match");
    TH1D *hProton_phi_match = (TH1D*)hProton_phi->Clone("hProton_phi_match");
  
    TH1D *hEle_pt_match = (TH1D*)hEle_pt->Clone("hEle_pt_match");
    TH1D *hEle_eta_match = (TH1D*)hEle_eta->Clone("hEle_eta_match");
    TH1D *hEle_phi_match = (TH1D*)hEle_phi->Clone("hEle_phi_match");
    
    TH1D *hPos_pt_match = (TH1D*)hPos_pt->Clone("hPos_pt_match");
    TH1D *hPos_eta_match = (TH1D*)hPos_eta->Clone("hPos_eta_match");
    TH1D *hPos_phi_match = (TH1D*)hPos_phi->Clone("hPos_phi_match");
  
    TH1D *h_matchDr = new TH1D("h_matchdr",";#DeltaR;",100,0,3);
    TH2D *h_matchDr_eta = new TH2D("h_matchdr_eta",";#DeltaR;#eta",50,0,3,50,0,7);

    TH1D *hPc_mass_reco = new TH1D ("hPc_mass_reco","; m (GeV)",50,0,6);
    TH1D *hPc_mass_gen = (TH1D*)hPc_mass_reco->Clone("hPc_mass_gen");

    TH1D *hJpsi_mass_reco = new TH1D ("hJpsi_mass_reco","; m (GeV)",50,0,6);
    TH1D *hJpsi_mass_gen = (TH1D*)hJpsi_mass_reco->Clone("hJpsi_mass_gen");


    TH2D *hJpsi_mass_eta_reco = new TH2D ("hJpsi_mass_eta_reco",";J/#psi #eta^{Reco};J/#psi mass^{Reco} (GeV)",10,0,5,50,2.5,3.5);
    TH2D *hPc_mass_eta_reco = new TH2D ("hPc_mass_eta_reco",";P_{c} #eta^{Reco};J/P_{c} mass^{Reco} (GeV)",10,0,5,100,4,5);

  // pT vs eta
    TH2D* hProton_gen2d =  new TH2D("hProton_gen2d",";#eta;p (GeV);",60,0,6,40,0,200);
    TH2D* hProton_gen2d_recoMatched = (TH2D*)hProton_gen2d->Clone("hProton_gen2d_recoMatched");

    TH2D* hEle_gen2d =  new TH2D("hEle_gen2d",";#eta;p_{T} (GeV);",60,0,6,40,0,20);
    TH2D* hEle_recoMatchGen2d = (TH2D*)hEle_gen2d->Clone("hEle_recoMatchGen2d");

  // p vs eta
    TH2D* hProton_gen2d_peta =  new TH2D("hProton_gen2d_peta",";#eta;p (GeV);",50,0,10,100,0,100);
    TH2D* hProton_recoMatchGen2d_peta = (TH2D*)hProton_gen2d_peta->Clone("hProton_recoMatchGen2d_peta");

    TH2D* hEle_gen2d_peta =  new TH2D("hEle_gen2d_peta",";#eta;p (GeV);",50,0,10,100,0,100);
    TH2D* hEle_recoMatchGen2d_peta = (TH2D*)hEle_gen2d_peta->Clone("hEle_recoMatchGen2d_peta");

  
    handsomeTH1(hProton_pt);  handsomeTH1(hProton_eta); handsomeTH1(hProton_phi);
    handsomeTH1(hProtonReco_pt);  handsomeTH1(hProtonReco_eta);  handsomeTH1(hProtonReco_EtaPt);
    handsomeTH1(hEle_pt); handsomeTH1(hEle_eta); handsomeTH1(hEle_phi);
    handsomeTH1(hPos_pt); handsomeTH1(hPos_eta); handsomeTH1(hPos_phi);

    handsomeTH1(hProton_pt_match,2);  handsomeTH1(hProton_eta_match,2); handsomeTH1(hProton_phi_match,2);
    handsomeTH1(hEle_pt_match,2); handsomeTH1(hEle_eta_match,2); handsomeTH1(hEle_phi_match,2);
    handsomeTH1(hPos_pt_match,2); handsomeTH1(hPos_eta_match,2); handsomeTH1(hPos_phi_match,2);

    
    int evtCount = 0;
  
  
    while(tree_reader.Next()) { // Loop over events
        evtCount++;
        if (evtCount > maxEvents)
        break;
	
	TLorentzVector proton_4vec_reco;
        TLorentzVector ele_4vec_reco;
        TLorentzVector pos_4vec_reco;
        TLorentzVector pc_4vec_reco;
        TLorentzVector jpsi_4vec_reco;
    
        TLorentzVector proton_4vec_gen;
        TLorentzVector ele_4vec_gen;
        TLorentzVector pos_4vec_gen;
        TLorentzVector pc_4vec_gen;
        TLorentzVector jpsi_4vec_gen;

	// Fill the RECO histograms:
        for (unsigned int j=0; j<recop_px.GetSize(); j++) {
      
            TVector3 recop3V(recop_px[j], recop_py[j], recop_pz[j]);
      
            float recoPt = sqrt(recop_px[j]*recop_px[j] + recop_py[j]*recop_py[j]);
      
            hRecop_p->Fill(recop3V.Mag());   // RECO p
            hRecop_eta->Fill(recop3V.PseudoRapidity());  // RECO eta
      
      
      // Match with Generated particles:
            double minDr = 100;
            int iGenOfMinDr = -1;
            for(unsigned int igen=0; igen<genp_pdg.GetSize(); igen++){ // Loop over thrown particles
	        if (genp_status[igen] != 1)
		  continue;
                TVector3 iGen3V(genp_px[igen], genp_py[igen], genp_pz[igen]);
		//     cout << " PseudoRapidity() = " << iGen3V.PseudoRapidity() << endl;
		//     cout << " phi: = " << iGen3V.Phi() << endl;
                double theDr = sqrt( pow(recop3V.PseudoRapidity()- iGen3V.PseudoRapidity(),2) + pow(recop3V.Phi()- iGen3V.Phi(),2));
		
                if ( theDr <  minDr) {
		  minDr = theDr;
		  iGenOfMinDr = igen;
                }
            }
            if ( minDr > 0.3 ) { // matching failed
                minDr = 100;
                iGenOfMinDr = -1;
            }
      
            int thePid = kElse;
            if ( iGenOfMinDr < 0 ) {   // If it's not matached
                thePid = kElse;
            }
            else if ((genp_status[iGenOfMinDr]==1)&&(genp_pdg[iGenOfMinDr]==2212))
	      thePid = kProton;
            else if ((genp_status[iGenOfMinDr]==1)&&(genp_pdg[iGenOfMinDr]==11))
	      thePid = kElectron;
            else if ((genp_status[iGenOfMinDr]==1)&&(genp_pdg[iGenOfMinDr]==-11))
	      thePid = kPositron;
            else
	      thePid = kElse;
	    
            TVector3 genMatch3V(0,0,0);
            if ( thePid != kElse) {
	      genMatch3V = TVector3(genp_px[iGenOfMinDr], genp_py[iGenOfMinDr], genp_pz[iGenOfMinDr]);
                h_matchDr->Fill(minDr);
                h_matchDr_eta->Fill( minDr, genMatch3V.Eta());
            }
            int thePBin = pBinHist->FindBin( genMatch3V.Mag()) - 1; // GEN p
            hRecop_p_match[thePid]->Fill( recop3V.Mag());
            hRecop_eta_match[thePid]->Fill( recop3V.PseudoRapidity());
            if ( (thePBin>=0) && (thePBin<4) ) {
	      hMatch_dp[thePid][thePBin]->Fill ( recop3V.Mag()/genMatch3V.Mag() - 1);
            }

	    if ( thePid == kProton)
	      proton_4vec_reco.SetXYZM(recop_px[j], recop_py[j], recop_pz[j], protonMass);
	    if ( thePid == kElectron)
	      ele_4vec_reco.SetXYZM(recop_px[j], recop_py[j], recop_pz[j], eleMass);
	    if ( thePid == kPositron)
	      pos_4vec_reco.SetXYZM(recop_px[j], recop_py[j], recop_pz[j], eleMass);
        }
	
	jpsi_4vec_reco =  ele_4vec_reco + pos_4vec_reco;
	pc_4vec_reco =  jpsi_4vec_reco + proton_4vec_reco;
	
	
        for(unsigned int i=0; i<genp_pdg.GetSize(); i++)  { // Loop over thrown particles
	  if (genp_status[i] != 1)
            continue;
	  
	  int thePid = 0;
            if ( genp_pdg[i]==2212) 
                thePid = kProton; 
            else if ( genp_pdg[i]==11) 
                thePid = kElectron; 
            else if ( (int)genp_pdg[i]==-11) 
                thePid = kPositron; 
            else  
                thePid = 0; 
      
            TVector3 gen3V(genp_px[i], genp_py[i], genp_pz[i]);
            hGen_p[thePid]->Fill(gen3V.Mag());
            hGen_eta[thePid]->Fill(gen3V.Eta());
            hGen_phi[thePid]->Fill(gen3V.Phi());
      
            float iGenEta = gen3V.PseudoRapidity();
            float iGenPhi = gen3V.Phi();
        
            //ele_4vec_reco.SetXYZM(recop_px[j], recop_py[j], recop_pz[j], eleMass);
            
      // if ( genp_pdg[i] == 2212) {
      // hProton_pt->Fill(iGenPt);
      // hProton_eta->Fill(iGenEta);
      // hProton_phi->Fill(iGenPhi);
      // proton_4vec_gen.SetXYZM(genp_px[i], genp_py[i], genp_pz[i], protonMass);
      
      // hProton_gen2d->Fill(iGenEta,iGenPt);
      // hProton_gen2d_peta->Fill(iGenEta,iGenP);
      //      }
      //      else if ( genp_pdg[i] == 11) {
      // hEle_pt->Fill(iGenPt);
      // hEle_eta->Fill(iGenEta);
      // hEle_phi->Fill(iGenPhi);
      // ele_4vec_gen.SetXYZM(genp_px[i], genp_py[i], genp_pz[i], eleMass);
      
      // hEle_gen2d->Fill(iGenEta,iGenPt);
      // hEle_gen2d_peta->Fill(iGenEta,iGenP);
      //      }
      //      else if ( genp_pdg[i] == -11) {
      // hPos_pt->Fill(iGenPt);
      // hPos_eta->Fill(iGenEta);
      // hPos_phi->Fill(iGenPhi);
      // pos_4vec_gen.SetXYZM(genp_px[i], genp_py[i], genp_pz[i], eleMass);
      //      }
      
      
            for (unsigned int j=0; j<recop_px.GetSize(); j++) {
	
                if (recop_simID[j] == i) { // Reconstructed by the algorihtm
	  
	  
	  //	  TVector3 genMatReco3V(recop_px[j], recop_py[j], recop_pz[j]);
                    TVector3 genMatReco3V(recop_px[j], recop_py[j], recop_pz[j]);
                    float genMatRecoPt = sqrt(recop_px[j]*recop_px[j] + recop_py[j]*recop_py[j]);
                    float genMatRecoEta = genMatReco3V.PseudoRapidity();
                    float genMatRecoPhi = genMatReco3V.Phi();
                    float matchDeta = genMatRecoEta - iGenEta;
                    float matchDphi = genMatRecoPhi - iGenPhi;
                    float matchDr = sqrt  ( matchDphi*matchDphi + matchDeta*matchDeta );

	  // h_matchDr->Fill(matchDr);
	  // h_matchDr_eta->Fill( matchDr, iGenEta);
	  // if ( genp_pdg[i] == 2212) {
	  //   hProton_pt_match->Fill(iGenPt);
	  //   hProton_eta_match->Fill(iGenEta);
	  //   hProton_phi_match->Fill(iGenPhi);
	  //   proton_4vec_reco.SetXYZM(recop_px[j], recop_py[j], recop_pz[j], protonMass);

	  //   hProton_recoMatchGen2d->Fill(iGenEta,iGenPt);
	  //   hProton_recoMatchGen2d_peta->Fill(iGenEta,iGenP);

	  //   hProtonReco_eta->Fill(proton_4vec_reco.Eta());
	  //   hProtonReco_pt->Fill(proton_4vec_reco.Pt());
	  //   hProtonReco_EtaPt->Fill( proton_4vec_reco.Eta(), proton_4vec_reco.Pt());
	  // }
	  // else if ( genp_pdg[i] == 11) {
	  //   hEle_pt_match->Fill(iGenPt);
	  //   hEle_eta_match->Fill(iGenEta);
	  //   hEle_phi_match->Fill(iGenPhi);
	  //   ele_4vec_reco.SetXYZM(recop_px[j], recop_py[j], recop_pz[j], eleMass);

	  //   hEle_recoMatchGen2d->Fill(iGenEta,iGenPt);
	  //   hEle_recoMatchGen2d_peta->Fill(iGenEta,iGenP);
	  // }
	  // else if ( genp_pdg[i] == -11) {
	  //   hPos_pt_match->Fill(iGenPt);
	  //   hPos_eta_match->Fill(iGenEta);
	  //   hPos_phi_match->Fill(iGenPhi);
	  //   pos_4vec_reco.SetXYZM(recop_px[j], recop_py[j], recop_pz[j], eleMass);
	  // }
	  
                }
	
	
            }
      // for(unsigned int j=0; j<recop_type.GetSize(); j++) 
      // { 		              
      //     if(recop_[j] == i) // Find association index matching the index of the thrown particle we are looking at 	
      //     {
      //         TVector3 recMom(trackMomX[recoAssoc[j]],trackMomY[recoAssoc[j]],trackMomZ[recoAssoc[j]]); // recoAssoc[j] is the index of the matched ReconstructedChargedParticle                       
      //     }
      // }
        }
        if (proton_4vec_gen.X() !=0 && ele_4vec_gen.X() !=0 && pos_4vec_gen.X() !=0 ) {
	  //            pc_4vec_gen = proton_4vec_gen + ele_4vec_gen +pos_4vec_gen;
            //hPc_mass_gen->Fill (pc_4vec_gen.M());
        }


      
        if ( ele_4vec_gen.X() !=0 && pos_4vec_gen.X() !=0 ) {
	  //            jpsi_4vec_gen =  ele_4vec_gen +pos_4vec_gen;
	  //            hJpsi_mass_gen->Fill (jpsi_4vec_gen.M());
        }

	jpsi_4vec_reco = ele_4vec_reco +pos_4vec_reco;
	hJpsi_mass_reco->Fill (jpsi_4vec_reco.M());
	hJpsi_mass_eta_reco->Fill (jpsi_4vec_reco.Eta(), jpsi_4vec_reco.M());
	pc_4vec_reco = proton_4vec_reco + ele_4vec_reco +pos_4vec_reco;
	hPc_mass_reco->Fill (pc_4vec_reco.M());
	hPc_mass_eta_reco->Fill (pc_4vec_reco.Eta(), pc_4vec_reco.M());
	
    
    }



    TCanvas* cMom = new TCanvas("cMom","",1200,600);
    cMom->Divide(2,1);
    cMom->cd(1);
    hRecop_p->Sumw2();
    handsomeTH1(hRecop_p,1);
    handsomeTH1(hRecop_p_match[kProton],2);
    handsomeTH1(hRecop_p_match[kElectron],4);
    handsomeTH1(hRecop_p_match[kPositron],6);
    handsomeTH1(hRecop_p_match[kElse],8);

    easyRange(hRecop_p);
    hRecop_p->SetLineWidth(4);
    hRecop_p->Draw();
    for ( int i=0 ; i<4 ; i++) {
        hRecop_p_match[i]->SetLineWidth(3);
        hRecop_p_match[i]->Draw("same hist");
    }
    auto legend1 = new TLegend(0.4671329,0.6909621,0.9972028,0.8556851,NULL,"brNDC");
    easyLeg(legend1, "");
    legend1->AddEntry(hRecop_p,"Reconstructed particles","pe");
    legend1->AddEntry(hRecop_p_match[kProton],"matched to p");
    legend1->AddEntry(hRecop_p_match[kElectron],"matched to e^{-}");
    legend1->AddEntry(hRecop_p_match[kPositron],"matched to e^{+}");
    legend1->AddEntry(hRecop_p_match[kElse],"else");
    legend1->Draw();
    
    cMom->cd(2);
    handsomeTH1(hRecop_eta,1);
    handsomeTH1(hRecop_eta_match[kProton],2);
    handsomeTH1(hRecop_eta_match[kElectron],4);
    handsomeTH1(hRecop_eta_match[kPositron],6);
    handsomeTH1(hRecop_eta_match[kElse],8);
    easyRange(hRecop_eta);
    hRecop_eta->SetLineWidth(4);
    hRecop_eta->Draw();
    for ( int i=0 ; i<4 ; i++) {
        hRecop_eta_match[i]->SetLineWidth(4);
        hRecop_eta_match[i]->Draw("same hist");
    }
    

    
    TCanvas* cRes = new TCanvas("cRes","",1200,900);
    cRes->Divide(4,3);
    for ( int ipid= 1 ; ipid <4 ; ipid++) {
        for ( int ip = 0 ; ip<4 ; ip++ ) {
            cRes->cd((ipid-1)*4 + ip+1);
            hMatch_dp[ipid][ip]->Draw();
            
            TFitResultPtr fitResult =  hMatch_dp[ipid][ip]->Fit("gaus", "S");
            TF1 *fitFunc = hMatch_dp[ipid][ip]->GetFunction("gaus");
            if (fitResult.Get()) {
                hResP[ipid]->SetBinContent( ip+1, fitFunc->GetParameter(2));
                hResP[ipid]->SetBinError( ip+1, fitFunc->GetParError(2));
            }
        }
    }

    TCanvas* cRes2 = new TCanvas("cRes2","",500,500);
    hResP[kProton]->SetAxisRange(0,0.3,"Y");
    handsomeTH1(hResP[kProton],2);
    handsomeTH1(hResP[kElectron],4);
    handsomeTH1(hResP[kPositron],6);
 
    hResP[kProton]->Draw();
    hResP[kElectron]->Draw("same");
    hResP[kPositron]->Draw("same");
    

        
    TCanvas* c0 = new TCanvas("c0","",1200,800);
    c0->Divide(3,2);
    c0->cd(1);
    hRecop_pt->Draw();
    auto legend0 = new TLegend(0.2744468,0.5504032,0.9244663,0.8004032,NULL,"brNDC");
    easyLeg(legend0, "Reconstructd objects");
    legend0->Draw();
    
    c0->cd(2);
    hRecop_eta->Draw();
    c0->cd(3);
    hRecop_phi->Draw();

    
    TCanvas* c1 = new TCanvas("c1","",1000,1000);
    c1->Divide(3,3);
    for ( int kpid=1 ; kpid<=3 ; kpid++) {
        c1->cd(1+3*(kpid-1));
        hGen_p[kpid]->Draw("hist");
        //hProton_pt_match->Draw("same");
        c1->cd(2+3*(kpid-1));
        hGen_eta[kpid]->Draw("hist");
        //hProton_eta_match->Draw("same");
        c1->cd(3+3*(kpid-1));
        hGen_phi[kpid]->Draw("hist");
        //hProton_phi_match->Draw("same");
    }


  
    TH1D* hProton_pt_eff = (TH1D*)hProton_pt_match->Clone("hProton_pt_eff");
    TH1D* hEle_pt_eff = (TH1D*)hEle_pt_match->Clone("hEle_pt_eff");
    TH1D* hPos_pt_eff = (TH1D*)hPos_pt_match->Clone("hPos_pt_eff");
    hProton_pt_eff->Divide(hProton_pt);
    hEle_pt_eff->Divide(hEle_pt);
    hPos_pt_eff->Divide(hPos_pt);
  
    handsomeTH1(hProton_pt_eff,1);
    hProton_pt_eff->SetAxisRange(0,1.2,"Y");
    hEle_pt_eff->SetAxisRange(0,1.2,"Y");
    hPos_pt_eff->SetAxisRange(0,1.2,"Y");


  
    TCanvas* c2 = new TCanvas("c2","",1200,400);
    c2->Divide(3,1);
    c2->cd(1);
    hProton_pt_eff->Draw();
    c2->cd(2);
    hEle_pt_eff->Draw();
    c2->cd(3);
    hPos_pt_eff->Draw();

    TCanvas* c2d = new TCanvas("c2d","",800,800);
    c2d->Divide(2,2);
    c2d->cd(1);
    //    TH2D* hProton_eff2d =  (TH2D*)hProton_recoMatchGen2d->Clone("hProton_recoMatchGen2d");
    //    hProton_eff2d->Divide(hProton_gen2d);
    //    handsomeTH2(hProton_eff2d);
    //    hProton_eff2d->SetStats(0);
    //    hProton_eff2d->Draw("colz");

    c2d->cd(2);
    TH2D* hEle_eff2d =  (TH2D*)hEle_recoMatchGen2d->Clone("hEle_recoMatchGen2d");
    hEle_eff2d->Divide(hEle_gen2d);
    handsomeTH2(hEle_eff2d);
    hEle_eff2d->SetStats(0);
    hEle_eff2d->Draw("colz");

    c2d->cd(3);
    TH2D* hProton_eff2d_peta =  (TH2D*)hProton_recoMatchGen2d_peta->Clone("hProton_recoMatchGen2d_peta");
    hProton_eff2d_peta->Divide(hProton_gen2d_peta);
    handsomeTH2(hProton_eff2d_peta);
    hProton_eff2d_peta->SetStats(0);
    hProton_eff2d_peta->Draw("colz");

    c2d->cd(4);
    TH2D* hEle_eff2d_peta =  (TH2D*)hEle_recoMatchGen2d_peta->Clone("hEle_recoMatchGen2d_peta");
    hEle_eff2d_peta->Divide(hEle_gen2d_peta);
    handsomeTH2(hEle_eff2d_peta);
    hEle_eff2d_peta->SetStats(0);
    hEle_eff2d_peta->Draw("colz");


  
    TCanvas* c3 = new TCanvas("c3","",800,400);
    c3->Divide(2,1);
    c3->cd(1);
    h_matchDr->Draw();
    c3->cd(2);
    h_matchDr_eta->Draw("colz");

    TCanvas* c4 = new TCanvas("c4","",800,400);
    c4->Divide(2,1);
    c4->cd(1);
    handsomeTH1(hPc_mass_gen,1);
    handsomeTH1(hPc_mass_reco,2);
    hPc_mass_reco->SetMarkerStyle(20);
    hPc_mass_reco->SetMarkerSize(1);
    hPc_mass_reco->Draw();
    hPc_mass_gen->Draw("same");
    c4->cd(2);
    handsomeTH1(hJpsi_mass_gen,1);
    handsomeTH1(hJpsi_mass_reco,2);
  //  hJpsi_mass_gen->Draw();
    hJpsi_mass_reco->SetMarkerStyle(20);
    hJpsi_mass_reco->SetMarkerSize(1);
    hJpsi_mass_reco->Draw();

    TCanvas* c5 = new TCanvas("c5","",800,400);
    c5->Divide(2,1);
    c5->cd(1);
    hJpsi_mass_eta_reco->Draw("colz");
    c5->cd(2);
    hPc_mass_eta_reco->Draw("colz");

    TCanvas* c6 = new TCanvas("c6","",1200,800);
    c6->Divide(3,2);
    c6->cd(1);
    hProtonReco_pt->Draw();
    gPad->SetLogy();
    c6->cd(2);
    hProtonReco_eta->Draw();
    gPad->SetLogy();
    c6->cd(3);
    hProtonReco_EtaPt->Draw("colz");

  
  
  
}

