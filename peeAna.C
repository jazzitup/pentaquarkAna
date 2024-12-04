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

int maxEvents = 100000;

float protonMass = 9.3827e-01;
float eleMass =  0.000511;

//void test(TString infile="reconstructed_data_N2100_neutron_theta_0_0.2mard_25GeV_OnlyHcal_info.root") {
//void peeAna(TString infile="podio_output_Pentaquark_hepmc_output_20241112_50GeV_10evts.hepmc.edm4hep.root") {
// void peeAna(TString infile="podio_output_Pentaquark_hepmc_output_20241113_50GeV_10000evts.hepmc.edm4hep.root") { 1k E = 50 GeV
void peeAna(TString infile="podio_files/Pentaquark_hepmc_output_20241202_p275.0GeV_e18.0GeV_two_body_kinematics_eta1.9-8_100000evts_ip6_hidiv_275x18.root") {
 // void peeAna(TString infile="podio_files/Pentaquark_hepmc_output_20241202_p275.0GeV_e18.0GeV_two_body_kinematics_eta4-8_10000evts_ip6_hidiv_275x18.root") { 
//void peeAna(TString infile="podio.root") {
  
    const int kElse = 0;
    const int kProton = 1;
    const int kElectron = 2;
    const int kPositron = 3;

    bool isGenMatched[500];
        
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

    TH2D *hRecop_peta = new TH2D("hRecop_peta",";#eta^{RECO};p^{RECO} (GeV)",100,1,7,100,0,300);
    handsomeTH1(hRecop_peta);
    TH2D *hRecop_pteta = new TH2D("hRecop_pteta",";#eta^{RECO};p_{T}^{RECO} (GeV)",100,1,7,100,0,50);
    handsomeTH1(hRecop_pteta);
    
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
        hResP[i]->SetYTitle("Momentum resolution(#Deltap/p)"); 
    }
  
  
  
  
  
  // Gen Particle histograms
    TH1D *hGen_p[4];
    TH1D *hGen_p_recoMatched[4];
    TH1D *hGen_eta[4];
    TH1D *hGen_eta_recoMatched[4];
    TH1D *hGen_phi[4];
    TH1D *hGen_phi_recoMatched[4];
    TH1D *hEff_p[4];
    TH1D *hEff_eta[4];
    TH1D *hEff_phi[4];

    TH2D *hGen_peta[4];
    TH2D *hGen_peta_recoMatched[4];
    TH2D *hEff_peta[4];

    TH2D *hGenP_RecoP[4];
    TH2D *hGenEta_RecoEta[4];
    
    
    for ( int id=0 ; id<=3 ; id++) {
        hGen_p[id] = new TH1D(Form("hGen_p_id%d",id),"; p^{GEN} (GeV);",100,0,200);
        hGen_p_recoMatched[id] = (TH1D*)hGen_p[id]->Clone(Form("hGen_p_recoMatched_id%d",id));
        hGen_eta[id] = new TH1D(Form("hGen_eta_id%d",id),"; #eta^{GEN};",100,1,7);
        hGen_eta_recoMatched[id] = (TH1D*)hGen_eta[id]->Clone(Form("hGen_eta_recoMatched_id%d",id));
        hGen_phi[id] = new TH1D(Form("hGen_phi_id%d",id),"; #phi^{GEN};",50,-3.2,3.2);
        hGen_phi_recoMatched[id] = (TH1D*)hGen_phi[id]->Clone(Form("hGen_phi_recoMatched_id%d",id));

        hGen_peta[id] = new TH2D(Form("hgen_petaid%d",id),"; #eta^{GEN}; p^{GEN} (GeV)",60,1,5.5,50,0,300);
        hGen_peta_recoMatched[id] = (TH2D*)hGen_peta[id]->Clone(Form("hgen_petaid%d_recomatched",id));
        handsomeTH1(hGen_peta_recoMatched[id]);

        hGenP_RecoP[id] = new TH2D(Form("hGenP_RecoP_id%d",id),";p^{GEN} (GeV); p^{RECO} (GeV)",50,0,200,50,0,200);
        hGenEta_RecoEta[id] = new TH2D(Form("hGeneta_RecoEta_id%d",id),";#eta^{GEN} (GeV); #eta^{RECO}", 50, 0, 6, 50, 0, 6);
    }
  
    TH1D *hProton_pt = new TH1D("hProton_pt",";p_{T} (GeV)",100,0,30);
    TH1D *hProton_eta = new TH1D ("hProton_eta",";#eta",100,0,7);
    TH1D *hProton_phi = new TH1D ("hProton_phi",";#phi",100,3.141592, -3.141592);
  
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
  
    TH1D *h_matchDr = new TH1D("h_matchdr",";#DeltaR;",200,0,0.5);
    TH2D *h_matchDr_eta = new TH2D("h_matchdr_eta",";#DeltaR;#eta",200,0,0.5,140,0,7);

    TH1D *hPc_mass_reco = new TH1D ("hPc_mass_reco","; m (GeV)",50,2,5);
    TH1D *hPc_mass_gen = (TH1D*)hPc_mass_reco->Clone("hPc_mass_gen");

    TH1D *hJpsi_mass_reco = new TH1D ("hJpsi_mass_reco","; m (GeV)",50,2,5);
    TH1D *hJpsi_mass_gen = (TH1D*)hJpsi_mass_reco->Clone("hJpsi_mass_gen");
    
    TH1D *hJpsi_eta_gen = new TH1D("hJpsi_eta_gen",";#eta^{J/#psi};",50,0,7);
    TH1D *hJpsi_eta_reco = (TH1D*)hJpsi_eta_gen->Clone("hJpsi_eta_reco");
    TH1D *hPc_eta_gen = new TH1D("hJpsi_pc_gen",";#eta^{P_{c}};",50,0,7);
    TH1D *hPc_eta_reco = (TH1D*)hPc_eta_gen->Clone("hPc_eta_reco");
    


    TH2D *hJpsi_mass_eta_reco = new TH2D ("hJpsi_mass_eta_reco",";J/#psi #eta^{Reco};J/#psi mass^{Reco} (GeV)",10,0,5,50,2.5,3.5);
    TH2D *hPc_mass_eta_reco = new TH2D ("hPc_mass_eta_reco",";P_{c} #eta^{Reco};J/P_{c} mass^{Reco} (GeV)",10,0,5,100,4,5);
    handsomeTH1(hJpsi_mass_eta_reco);
        
  // pT vs eta
    TH2D* hProton_gen2d =  new TH2D("hProton_gen2d",";#eta;p (GeV);",60,0,6,40,0,200);
    TH2D* hProton_gen2d_recoMatched = (TH2D*)hProton_gen2d->Clone("hProton_gen2d_recoMatched");

    TH2D* hEle_gen2d =  new TH2D("hEle_gen2d",";#eta;p_{T} (GeV);",60,0,6,40,0,20);
    TH2D* hEle_recoMatchGen2d = (TH2D*)hEle_gen2d->Clone("hEle_recoMatchGen2d");

  
  
    handsomeTH1(hProton_pt);  handsomeTH1(hProton_eta); handsomeTH1(hProton_phi);
    
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
	

        TLorentzVector proton_4vec_gen;
        TLorentzVector ele_4vec_gen;
        TLorentzVector pos_4vec_gen;
        TLorentzVector pc_4vec_gen;
        TLorentzVector jpsi_4vec_gen;

	    TLorentzVector proton_4vec_reco;
        TLorentzVector ele_4vec_reco;
        TLorentzVector pos_4vec_reco;
        TLorentzVector pc_4vec_reco;
        TLorentzVector jpsi_4vec_reco;


        // Clear the isGenMatched before matching the GEN with RECO
        for (unsigned int ireco=0; ireco<recop_px.GetSize(); ireco++) {  
            isGenMatched[ireco] = false; 
        }
        
        // LOOP over GEN particles 
        for(unsigned int i=0; i<genp_pdg.GetSize(); i++)  { // Loop over thrown particles

            if (genp_status[i] != 1)
                continue;
            TVector3 gen3V(genp_px[i], genp_py[i], genp_pz[i]);
            
            int thePid = kElse;
            if ( genp_pdg[i]==2212)  { 
                thePid = kProton; 
                proton_4vec_gen.SetXYZM(gen3V.Px(), gen3V.Py(), gen3V.Pz(), protonMass);
            }
            else if ( genp_pdg[i]==11) {
                thePid = kElectron; 
                ele_4vec_gen.SetXYZM(gen3V.Px(), gen3V.Py(), gen3V.Pz(), eleMass);
            }
            else if ( (int)genp_pdg[i]==-11)  {
                thePid = kPositron; 
                pos_4vec_gen.SetXYZM(gen3V.Px(), gen3V.Py(), gen3V.Pz(), eleMass);
            }
            else  {
                thePid = kElse; 
            }
            
            hGen_p[thePid]->Fill(gen3V.Mag());
            hGen_eta[thePid]->Fill(gen3V.Eta());
            hGen_phi[thePid]->Fill(gen3V.Phi());
            hGen_peta[thePid]->Fill(gen3V.Eta(), gen3V.Mag());

            // Match with RECO:
            double minDr = 0.3;
            int iRecoOfMinDr = -1;
            
            // LOOP over RECO collection 
            for (unsigned int ireco=0; ireco<recop_px.GetSize(); ireco++) {
                if (isGenMatched[ireco] == true) 
                    continue;
                TVector3 recop3V(recop_px[ireco], recop_py[ireco], recop_pz[ireco]);
        		double theDr = sqrt( pow(recop3V.PseudoRapidity()- gen3V.PseudoRapidity(),2) + pow(recop3V.Phi()- gen3V.Phi(),2));
		
                if ( theDr <  minDr) {
        		  minDr = theDr;
        		  iRecoOfMinDr = ireco;
                }
            }
            
            if (iRecoOfMinDr >= 0 ) {  // If matched
                isGenMatched[iRecoOfMinDr] = true;
                TVector3 matchedReco3V(recop_px[iRecoOfMinDr], recop_py[iRecoOfMinDr], recop_pz[iRecoOfMinDr]);
                h_matchDr->Fill(minDr);
                h_matchDr_eta->Fill( minDr, gen3V.Eta());            
                hGenP_RecoP[thePid]->Fill( gen3V.Mag(), matchedReco3V.Mag());
                hGenEta_RecoEta[thePid]->Fill( gen3V.Eta(), matchedReco3V.Eta());

                hGen_p_recoMatched[thePid]->Fill   (gen3V.Mag());
                hGen_eta_recoMatched[thePid]->Fill (gen3V.Eta());
                hGen_phi_recoMatched[thePid]->Fill (gen3V.Phi());

                hGen_peta_recoMatched[thePid]->Fill( gen3V.Eta(), gen3V.Mag());           
                hRecop_p_match[thePid]->Fill( matchedReco3V.Mag());
                hRecop_eta_match[thePid]->Fill( matchedReco3V.PseudoRapidity());

                // Resolution:
                int thePBin = pBinHist->FindBin(gen3V.Mag()) - 1; // GEN p binning   
                if ( (thePBin>=0) && (thePBin<4) ) {
                  hMatch_dp[thePid][thePBin]->Fill ( matchedReco3V.Mag()/gen3V.Mag() - 1);
                }
            

                // Now assign the p,e,ebar
                if ( thePid == kProton)
                  proton_4vec_reco.SetXYZM(matchedReco3V.Px(), matchedReco3V.Py(), matchedReco3V.Pz(), protonMass);
                if ( thePid == kElectron)
                  ele_4vec_reco.SetXYZM(matchedReco3V.Px(), matchedReco3V.Py(), matchedReco3V.Pz(), eleMass);
                if ( thePid == kPositron)
                  pos_4vec_reco.SetXYZM(matchedReco3V.Px(), matchedReco3V.Py(), matchedReco3V.Pz(), eleMass);
            }

            
        }
    
        if ( (ele_4vec_gen.M()>0) && (pos_4vec_gen.M()>0) )  { // If all particles are generated
            jpsi_4vec_gen = ele_4vec_gen +pos_4vec_gen;
            hJpsi_eta_gen->Fill (jpsi_4vec_gen.Eta());
            

        }
        
        if (  (proton_4vec_gen.M()>0) && (ele_4vec_gen.M()>0) && (pos_4vec_gen.M()>0) )  { 
            pc_4vec_gen = proton_4vec_gen + ele_4vec_gen +pos_4vec_gen;
            hPc_eta_gen->Fill ( pc_4vec_gen.Eta());
        }

        if (  (proton_4vec_reco.M()>0) && (ele_4vec_reco.M()>0) && (pos_4vec_reco.M()>0) )  { // If all particles are reconstructed
            pc_4vec_reco =  jpsi_4vec_reco + proton_4vec_reco;
                
            jpsi_4vec_reco = ele_4vec_reco +pos_4vec_reco;
            hJpsi_mass_reco->Fill (jpsi_4vec_reco.M());
            hJpsi_eta_reco->Fill(jpsi_4vec_reco.Eta());
            hJpsi_mass_eta_reco->Fill (jpsi_4vec_reco.Eta(), jpsi_4vec_reco.M());
            pc_4vec_reco = proton_4vec_reco + ele_4vec_reco +pos_4vec_reco;
            hPc_mass_reco->Fill (pc_4vec_reco.M());
            hPc_eta_reco->Fill(pc_4vec_reco.Eta());
            hPc_mass_eta_reco->Fill (pc_4vec_reco.Eta(), pc_4vec_reco.M());
        }

        
        // Loop over RECO collection 
        for (unsigned int j=0; j<recop_px.GetSize(); j++) {
            TVector3 recop3V(recop_px[j], recop_py[j], recop_pz[j]);
            float recoPt = sqrt(recop_px[j]*recop_px[j] + recop_py[j]*recop_py[j]);
            hRecop_p->Fill(recop3V.Mag());   // RECO p
            hRecop_eta->Fill(recop3V.PseudoRapidity());  // RECO eta
            hRecop_peta->Fill( recop3V.PseudoRapidity(), recop3V.Mag());
            hRecop_pteta->Fill( recop3V.PseudoRapidity(), recoPt);
        }
               
    
        	
    
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

    
    


    TCanvas* cResp = new TCanvas("cResp","",1200,800);
    cResp->Divide(3,2);
    for ( int i = 1 ; i<=3 ; i++) { 
        cResp->cd(i);
        hGenP_RecoP[i]->Draw("colz");
        cResp->cd(i+3);
        hGenEta_RecoEta[i]->Draw("colz");
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
    
    auto legend2 =  new TLegend(0.4959839,0.6547368,0.9839357,0.8378947,NULL,"brNDC");
    easyLeg(legend2, "#frac{#Deltap}{p}");
    legend2->AddEntry(hResP[kProton],"proton","pe");
    legend2->AddEntry(hResP[kElectron],"electron","pe");
    legend2->AddEntry(hResP[kPositron],"positron","pe");
    legend2->Draw();
    
    

        
    TCanvas* c0 = new TCanvas("c0","",1200,400);
    c0->Divide(3,1);
    for ( int kpid=1 ; kpid<=3 ; kpid++) {
        c0->cd(kpid);
        hGen_peta[kpid]->Draw("colz");
    }
    
  
    
    TCanvas* c1 = new TCanvas("c1","",1000,1000);
    c1->Divide(3,3);
    for ( int kpid=1 ; kpid<=3 ; kpid++) {
        c1->cd(kpid);
        handsomeTH1(hGen_p[kpid],1);
        handsomeTH1(hGen_p_recoMatched[kpid],2);
        hGen_p[kpid]->Draw("hist");
        hGen_p_recoMatched[kpid]->Draw("same");
        
        //hProton_pt_match->Draw("same");
        c1->cd(kpid+3);
        handsomeTH1(hGen_eta[kpid],1);
        handsomeTH1(hGen_eta_recoMatched[kpid],2);
        hGen_eta[kpid]->Draw("hist");
        hGen_eta_recoMatched[kpid]->Draw("same");
        
        //hProton_eta_match->Draw("same");
        c1->cd(kpid+6);
        handsomeTH1(hGen_phi[kpid],1);
        handsomeTH1(hGen_phi_recoMatched[kpid],2);
        hGen_phi[kpid]->Draw("hist");
        hGen_phi_recoMatched[kpid]->Draw("same");
    }

     // efficiency 
    TCanvas* cEff = new TCanvas("cEff","",1000,1000);
    cEff->Divide(3,3);

    for ( int kpid=1 ; kpid<=3 ; kpid++) {
        cEff->cd(kpid);
        hEff_p[kpid] = (TH1D*)hGen_p_recoMatched[kpid]->Clone(Form("hEff_p_kpid%d",kpid));
        hEff_p[kpid]->Divide(hGen_p[kpid]);
        hEff_p[kpid]->SetAxisRange(0,1.2,"Y");
        hEff_p[kpid]->Draw();
        
    
        cEff->cd(kpid+3);
        hEff_eta[kpid] = (TH1D*)hGen_eta_recoMatched[kpid]->Clone(Form("hEff_eta_kpid%d",kpid));
        hEff_eta[kpid]->Divide(hGen_eta[kpid]);
        easyRange(hEff_eta[kpid]);
        hEff_eta[kpid]->SetAxisRange(0,1.2,"Y");
        hEff_eta[kpid]->Draw();
        
        //hProton_eta_match->Draw("same");
        cEff->cd(kpid+6);
        hEff_phi[kpid] = (TH1D*)hGen_phi_recoMatched[kpid]->Clone(Form("hEff_phi_kpid%d",kpid));
        hEff_phi[kpid]->Divide(hGen_phi[kpid]);
        easyRange(hEff_phi[kpid]);
        hEff_phi[kpid]->SetAxisRange(0,1.2,"Y");
        hEff_phi[kpid]->Draw();
    
    }

    cEff->cd(1); 
    addText("Reco. efficiency",100,1.05);
    addText("proton",80,0.9);
    cEff->cd(2); 
    addText("electron",80,0.9);
    cEff->cd(3); 
    addText("positron",80,0.9);
        

    



    TCanvas* c2d = new TCanvas("c2d","",1200,400);
    c2d->Divide(3,1);
    for ( int kpid=1 ; kpid<=3 ; kpid++) {
        c2d->cd(kpid);
        hEff_peta[kpid] = (TH2D*)hGen_peta_recoMatched[kpid]->Clone(Form("hGen_peta_recoMatched_kpid%d",kpid));
        hEff_peta[kpid]->Divide(hGen_peta[kpid]);
        hEff_peta[kpid]->Draw("colz");
        auto legendM = new TLegend(0.1671329,0.6909621,0.8972028,0.856851,NULL,"brNDC");
        if (kpid == kProton) easyLeg(legendM, "Trk. efficiency of Proton");
        if (kpid == kElectron) easyLeg(legendM, "Trk. efficiency of Electron");
        if (kpid == kPositron) easyLeg(legendM, "Trk. efficiency of Positron");
        legendM->Draw();
    
    }    

    TCanvas* cGen2d = new TCanvas("cGen2d","",800,400);
    cGen2d->Divide(2,1);
    cGen2d->cd(1);
    hRecop_peta->Draw("colz");
    auto legend11 = new TLegend(0.1671329,0.7509621,0.8972028,0.856851,NULL,"brNDC");
    easyLeg(legend11, "All tracks");
    legend11->Draw();
 
    cGen2d->cd(2);
    hRecop_pteta->Draw("colz");
    
    TCanvas* c3 = new TCanvas("c3","",800,400);
    c3->Divide(2,1);
    c3->cd(1);
    handsomeTH1(h_matchDr);
    easyRange(h_matchDr, 100);
    h_matchDr->Draw();
    gPad->SetLogy();
    c3->cd(2);
    h_matchDr_eta->Draw("colz");

    TCanvas* c4 = new TCanvas("c4","",500,500);
    c4->cd(1);
    handsomeTH1(hJpsi_mass_gen,1);
    handsomeTH1(hJpsi_mass_reco,1);

    
    handsomeTH1(hJpsi_mass_reco,1);
    handsomeTH1(hPc_mass_reco,2);
    easyRange(hJpsi_mass_reco);
     hJpsi_mass_reco->SetMarkerSize(1);
     hPc_mass_reco->SetMarkerSize(1);
    
    hJpsi_mass_reco->Draw();
    hPc_mass_reco->Draw("same");
    auto legendM = new TLegend(0.3671329,0.6909621,0.8972028,0.8556851,NULL,"brNDC");
    easyLeg(legendM, "");
    legendM->AddEntry(hJpsi_mass_reco,"e^{+} + e^{-} #leftarrow J/#psi ","pe");
    legendM->AddEntry(hRecop_p_match[kProton],"#it{p} + e^{+} + e^{-} #leftarrow P_{c}","pe");
    legendM->Draw();
    
    
    TCanvas* c5 = new TCanvas("c5","",800,400);
    c5->Divide(2,1);
    c5->cd(1);
    hJpsi_mass_eta_reco->Draw("colz");
    c5->cd(2);
    hPc_mass_eta_reco->Draw("colz");  

    TCanvas* c6 = new TCanvas("c6","",800,400);
    c6->Divide(2,1);
    c6->cd(1);
    TH1D* hJpsi_eta_eff = (TH1D*)hJpsi_eta_reco->Clone("hjpsi_eta_eff");
    for ( int i = 1 ; i<= hJpsi_eta_gen->GetNbinsX() ; i++) 
        hJpsi_eta_gen->SetBinError(i, 0);
    hJpsi_eta_eff->Divide(hJpsi_eta_gen);
    hJpsi_eta_eff->SetAxisRange(0,1.4,"Y");
    handsomeTH1(hJpsi_eta_eff);
    hJpsi_eta_eff->Draw();
    jumSun(0,1,7,1);
    auto legend6 = new TLegend(0.3671329,0.6909621,0.8972028,0.8556851,NULL,"brNDC");
    easyLeg(legend6, "Reco Efficiency");
    legend6->AddEntry(hJpsi_eta_eff,"e^{+} + e^{-} #leftarrow J/#psi ","pe");
    legend6->Draw();
 
    c6->cd(2);
    TH1D* hPc_eta_eff = (TH1D*)hPc_eta_reco->Clone("hPc_eta_eff");
    for ( int i = 1 ; i<= hPc_eta_gen->GetNbinsX() ; i++) 
        hPc_eta_gen->SetBinError(i, 0);
    hPc_eta_eff->Divide(hPc_eta_gen);
    hPc_eta_eff->SetAxisRange(0,1.4,"Y");
    handsomeTH1(hPc_eta_eff);
    hPc_eta_eff->Draw();
    jumSun(0,1,7,1);
    auto legend7 = new TLegend(0.3671329,0.6909621,0.8972028,0.8556851,NULL,"brNDC");
    easyLeg(legend7, "Reco Efficiency");
    legend7->AddEntry(hPc_eta_eff,"#it{p} + e^{+} + e^{-} #leftarrow P_{c}","pe");

    legend7->Draw();
 
    
  
}

