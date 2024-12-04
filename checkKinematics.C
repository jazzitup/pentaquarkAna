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

//void checkKinematics(TString infile="podio_files/Pentaquark_hepmc_output_20241202_p275.0GeV_e18.0GeV_two_body_kinematics_eta1.9-8_100000evts_ip6_hidiv_275x18.root") {
 // void checkKinematics(TString infile="podio_files/Pentaquark_hepmc_output_20241202_p275.0GeV_e18.0GeV_two_body_kinematics_eta4-8_10000evts_ip6_hidiv_275x18.root") { 
void checkKinematics(TString infile="podio.root") {
  
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


    TH2D* hPEta_ele_proton24 = new TH2D("hpeta_ele_proton24",";#eta;p (GeV)",30,0,6,200,0,200);    
    TH2D* hPEta_ele_proton45 = new TH2D("hpeta_ele_proton45",";#eta;p (GeV)",30,0,6,200,0,200);    
    
    
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
	    
        }
	
	if ( (proton_4vec_gen.Eta() > 2) &&(proton_4vec_gen.Eta() <4)) {
	  hPEta_ele_proton24->Fill(ele_4vec_gen.Eta(), ele_4vec_gen.P());
	  hPEta_ele_proton24->Fill(pos_4vec_gen.Eta(), pos_4vec_gen.P());
	}
	if ( (proton_4vec_gen.Eta() > 4) &&(proton_4vec_gen.Eta() <5)) {
	  hPEta_ele_proton45->Fill(ele_4vec_gen.Eta(), ele_4vec_gen.P());
	  hPEta_ele_proton45->Fill(pos_4vec_gen.Eta(), pos_4vec_gen.P());
	}
    }
    
    TCanvas* c0 = new TCanvas("c0","",800,400);
    c0->Divide(2,1);
    c0->cd(1);
    hPEta_ele_proton24->Draw("colz");
    c0->cd(2);
    hPEta_ele_proton45->Draw("colz");

  
}

