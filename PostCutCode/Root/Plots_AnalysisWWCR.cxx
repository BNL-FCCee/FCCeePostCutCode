#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <iostream>
#include <TGraph.h>

void Plots_AnalysisWWCR() {

    TFile *test_CR = TFile::Open("/usatlas/u/kleibensperger/FCCeePostCutCode/run/test.root", "READ");
    TFile *test_noCR = TFile::Open("/usatlas/u/kleibensperger/FCCeePostCutCode/run/test_noCR.root", "READ");
    TFile *test_SKI = TFile::Open("/usatlas/u/kleibensperger/FCCeePostCutCode/run/test_SKI.root", "READ");
    TFile *test_SKII = TFile::Open("/usatlas/u/kleibensperger/FCCeePostCutCode/run/test_SKII.root", "READ");

    if (!test_CR || test_CR->IsZombie()) {
        std::cerr << "Error: Cannot open test_CR!" << std::endl;
        return;
    } 

    if (!test_noCR || test_noCR->IsZombie()) {
        std::cerr << "Error: Cannot open test_noCR!" << std::endl;
        return;
    } 

    if (!test_SKI || test_SKI->IsZombie()) {
        std::cerr << "Error: Cannot open test_SKI!" << std::endl;
        return;
    } 

    if (!test_SKII || test_SKII->IsZombie()) {
        std::cerr << "Error: Cannot open test_SKII!" << std::endl;
        return;
    } 

    // *************************** CR Plots ***************************
    TH1F *h_W1_truth_mass = (TH1F*)test_CR->Get("Nominal/h_W1_truth_mass_WW");
    TH1F *h_W2_truth_mass = (TH1F*)test_CR->Get("Nominal/h_W2_truth_mass_WW");
    TH1F *h_W1_mass = (TH1F*)test_CR->Get("Nominal/h_W1_mass_WW");
    TH1F *h_W2_mass = (TH1F*)test_CR->Get("Nominal/h_W2_mass_WW");

    TH1F *h_W1_truth_p = (TH1F*)test_CR->Get("Nominal/h_W1_truth_p_WW");
    TH1F *h_W2_truth_p = (TH1F*)test_CR->Get("Nominal/h_W2_truth_p_WW");
    TH1F *h_W1_p = (TH1F*)test_CR->Get("Nominal/h_W1_p_WW");
    TH1F *h_W2_p = (TH1F*)test_CR->Get("Nominal/h_W2_p_WW");

    TH1F *h_W1_truth_e = (TH1F*)test_CR->Get("Nominal/h_W1_truth_e_WW");
    TH1F *h_W2_truth_e = (TH1F*)test_CR->Get("Nominal/h_W2_truth_e_WW");
    TH1F *h_W1_e = (TH1F*)test_CR->Get("Nominal/h_W1_e_WW");
    TH1F *h_W2_e = (TH1F*)test_CR->Get("Nominal/h_W2_e_WW");

    TH1F *h_theta_truth_cs = (TH1F*)test_CR->Get("Nominal/h_theta_truth_cs_WW");
    TH1F *h_theta_truth_ud = (TH1F*)test_CR->Get("Nominal/h_theta_truth_ud_WW");
    TH1F *h_theta_c_l0 = (TH1F*)test_CR->Get("Nominal/h_theta_c_l0_WW");
    TH1F *h_theta_l1_l2 = (TH1F*)test_CR->Get("Nominal/h_theta_l1_l2_WW");

    TH1F *h_delta_theta_truth_cs = (TH1F*)test_CR->Get("Nominal/h_delta_theta_truth_cs_WW");
    TH1F *h_delta_theta_truth_ud = (TH1F*)test_CR->Get("Nominal/h_delta_theta_truth_ud_WW");
    TH1F *h_delta_theta_c_l0 = (TH1F*)test_CR->Get("Nominal/h_delta_theta_c_l0_WW");
    TH1F *h_delta_theta_l1_l2 = (TH1F*)test_CR->Get("Nominal/h_delta_theta_l1_l2_WW");

    TH1F *h_delta_eta_truth_cs = (TH1F*)test_CR->Get("Nominal/h_delta_eta_truth_cs_WW");
    TH1F *h_delta_eta_truth_ud = (TH1F*)test_CR->Get("Nominal/h_delta_eta_truth_ud_WW");
    TH1F *h_delta_eta_c_l0 = (TH1F*)test_CR->Get("Nominal/h_delta_eta_c_l0_WW");
    TH1F *h_delta_eta_l1_l2 = (TH1F*)test_CR->Get("Nominal/h_delta_eta_l1_l2_WW");

    TH1F *h_delta_phi_truth_cs = (TH1F*)test_CR->Get("Nominal/h_phi_truth_cs_WW");
    TH1F *h_delta_phi_truth_ud = (TH1F*)test_CR->Get("Nominal/h_phi_truth_ud_WW");
    TH1F *h_delta_phi_c_l0 = (TH1F*)test_CR->Get("Nominal/h_phi_c_l0_WW");
    TH1F *h_delta_phi_l1_l2 = (TH1F*)test_CR->Get("Nominal/h_phi_l1_l2_WW");

    TH1F *h_cos_phi_truth_cs = (TH1F*)test_CR->Get("Nominal/h_cos_phi_truth_cs_WW");
    TH1F *h_cos_phi_truth_ud = (TH1F*)test_CR->Get("Nominal/h_cos_phi_truth_ud_WW");
    TH1F *h_cos_phi_c_l0 = (TH1F*)test_CR->Get("Nominal/h_cos_phi_c_l0_WW");
    TH1F *h_cos_phi_l1_l2 = (TH1F*)test_CR->Get("Nominal/h_cos_phi_l1_l2_WW");

    TH1F *h_eec_truth_cs = (TH1F*)test_CR->Get("Nominal/h_eec_truth_cs_WW");
    TH1F *h_eec_truth_ud = (TH1F*)test_CR->Get("Nominal/h_eec_truth_ud_WW");
    TH1F *h_eec_c_l0 = (TH1F*)test_CR->Get("Nominal/h_eec_c_l0_WW");
    TH1F *h_eec_l1_l2 = (TH1F*)test_CR->Get("Nominal/h_eec_l1_l2_WW");

    // *************************** No CR Plots ***************************
    TH1F *h_W1_truth_mass_noCR = (TH1F*)test_noCR->Get("Nominal/h_W1_truth_mass_WW");
    TH1F *h_W2_truth_mass_noCR = (TH1F*)test_noCR->Get("Nominal/h_W2_truth_mass_WW");
    TH1F *h_W1_mass_noCR = (TH1F*)test_noCR->Get("Nominal/h_W1_mass_WW");
    TH1F *h_W2_mass_noCR = (TH1F*)test_noCR->Get("Nominal/h_W2_mass_WW");

    TH1F *h_W1_truth_p_noCR = (TH1F*)test_noCR->Get("Nominal/h_W1_truth_p_WW");
    TH1F *h_W2_truth_p_noCR = (TH1F*)test_noCR->Get("Nominal/h_W2_truth_p_WW");
    TH1F *h_W1_p_noCR = (TH1F*)test_noCR->Get("Nominal/h_W1_p_WW");
    TH1F *h_W2_p_noCR = (TH1F*)test_noCR->Get("Nominal/h_W2_p_WW");

    TH1F *h_W1_truth_e_noCR = (TH1F*)test_noCR->Get("Nominal/h_W1_truth_e_WW");
    TH1F *h_W2_truth_e_noCR = (TH1F*)test_noCR->Get("Nominal/h_W2_truth_e_WW");
    TH1F *h_W1_e_noCR = (TH1F*)test_noCR->Get("Nominal/h_W1_e_WW");
    TH1F *h_W2_e_noCR = (TH1F*)test_noCR->Get("Nominal/h_W2_e_WW");

    TH1F *h_theta_truth_cs_noCR = (TH1F*)test_noCR->Get("Nominal/h_theta_truth_cs_WW");
    TH1F *h_theta_truth_ud_noCR = (TH1F*)test_noCR->Get("Nominal/h_theta_truth_ud_WW");
    TH1F *h_theta_c_l0_noCR = (TH1F*)test_noCR->Get("Nominal/h_theta_c_l0_WW");
    TH1F *h_theta_l1_l2_noCR = (TH1F*)test_noCR->Get("Nominal/h_theta_l1_l2_WW");

    TH1F *h_delta_theta_truth_cs_noCR = (TH1F*)test_noCR->Get("Nominal/h_delta_theta_truth_cs_WW");
    TH1F *h_delta_theta_truth_ud_noCR = (TH1F*)test_noCR->Get("Nominal/h_delta_theta_truth_ud_WW");
    TH1F *h_delta_theta_c_l0_noCR = (TH1F*)test_noCR->Get("Nominal/h_delta_theta_c_l0_WW");
    TH1F *h_delta_theta_l1_l2_noCR = (TH1F*)test_noCR->Get("Nominal/h_delta_theta_l1_l2_WW");

    TH1F *h_delta_eta_truth_cs_noCR = (TH1F*)test_noCR->Get("Nominal/h_delta_eta_truth_cs_WW");
    TH1F *h_delta_eta_truth_ud_noCR = (TH1F*)test_noCR->Get("Nominal/h_delta_eta_truth_ud_WW");
    TH1F *h_delta_eta_c_l0_noCR = (TH1F*)test_noCR->Get("Nominal/h_delta_eta_c_l0_WW");
    TH1F *h_delta_eta_l1_l2_noCR = (TH1F*)test_noCR->Get("Nominal/h_delta_eta_l1_l2_WW");

    TH1F *h_delta_phi_truth_cs_noCR = (TH1F*)test_noCR->Get("Nominal/h_phi_truth_cs_WW");
    TH1F *h_delta_phi_truth_ud_noCR = (TH1F*)test_noCR->Get("Nominal/h_phi_truth_ud_WW");
    TH1F *h_delta_phi_c_l0_noCR = (TH1F*)test_noCR->Get("Nominal/h_phi_c_l0_WW");
    TH1F *h_delta_phi_l1_l2_noCR = (TH1F*)test_noCR->Get("Nominal/h_phi_l1_l2_WW");

    TH1F *h_cos_phi_truth_cs_noCR = (TH1F*)test_noCR->Get("Nominal/h_cos_phi_truth_cs_WW");
    TH1F *h_cos_phi_truth_ud_noCR = (TH1F*)test_noCR->Get("Nominal/h_cos_phi_truth_ud_WW");
    TH1F *h_cos_phi_c_l0_noCR = (TH1F*)test_noCR->Get("Nominal/h_cos_phi_c_l0_WW");
    TH1F *h_cos_phi_l1_l2_noCR = (TH1F*)test_noCR->Get("Nominal/h_cos_phi_l1_l2_WW");

    TH1F *h_eec_truth_cs_noCR = (TH1F*)test_noCR->Get("Nominal/h_eec_truth_cs_WW");
    TH1F *h_eec_truth_ud_noCR = (TH1F*)test_noCR->Get("Nominal/h_eec_truth_ud_WW");
    TH1F *h_eec_c_l0_noCR = (TH1F*)test_noCR->Get("Nominal/h_eec_c_l0_WW");
    TH1F *h_eec_l1_l2_noCR = (TH1F*)test_noCR->Get("Nominal/h_eec_l1_l2_WW");

    TFile *Plots_WWCR = new TFile("/usatlas/u/kleibensperger/FCCeePostCutCode/run/WWCR_Plots.root", "Recreate");

    // *************************** THETA DISTRIBUTIONS ***************************
    TCanvas* c_theta_truth = new TCanvas("c_theta_truth", "Truth Histograms", 1500, 600);
    c_theta_truth->Divide(2, 1);

    c_theta_truth->cd(1);
    h_theta_truth_cs->Draw("HIST P");
    h_theta_truth_cs->SetMarkerStyle(20);
    h_theta_truth_cs->SetMarkerColor(kBlue);
    h_theta_truth_cs->SetMarkerSize(0.8);
    h_theta_truth_cs->GetXaxis()->SetTitle("#theta");
    h_theta_truth_cs->SetTitle("#theta: W1 #rightarrow cs");
    h_theta_truth_cs->Rebin(2);
    h_theta_truth_cs->Scale(1.0 / h_theta_truth_cs->Integral());

    h_theta_truth_cs_noCR->Draw("HIST P same");
    h_theta_truth_cs_noCR->SetMarkerStyle(21);
    h_theta_truth_cs_noCR->SetMarkerColor(kRed);
    h_theta_truth_cs_noCR->SetMarkerSize(0.8);
    h_theta_truth_cs_noCR->Rebin(2);
    h_theta_truth_cs_noCR->Scale(1.0 / h_theta_truth_cs_noCR->Integral());

    c_theta_truth->cd(2);
    h_theta_truth_ud->Draw("HIST P");
    h_theta_truth_ud->SetMarkerStyle(20);
    h_theta_truth_ud->SetMarkerColor(kBlue);
    h_theta_truth_ud->SetMarkerSize(0.7);
    h_theta_truth_ud->GetXaxis()->SetTitle("#theta");
    h_theta_truth_ud->SetTitle("#theta: W2 #rightarrow ud");
    h_theta_truth_ud->Rebin(2);
    h_theta_truth_ud->Scale(1.0 / h_theta_truth_ud->Integral());

    h_theta_truth_ud_noCR->Draw("HIST P same");
    h_theta_truth_ud_noCR->SetMarkerStyle(21);
    h_theta_truth_ud_noCR->SetMarkerColor(kRed);
    h_theta_truth_ud_noCR->SetMarkerSize(0.7);
    h_theta_truth_ud_noCR->Rebin(2);
    h_theta_truth_ud_noCR->Scale(1.0 / h_theta_truth_ud_noCR->Integral());


    TCanvas* c_theta_reco = new TCanvas("c_theta_reco", "Reconstructed Histograms", 1500, 600);
    c_theta_reco->Divide(2, 1);

    c_theta_reco->cd(1);
    h_theta_c_l0->Draw("HIST P");
    h_theta_c_l0->SetMarkerStyle(20);
    h_theta_c_l0->SetMarkerColor(kBlue);
    h_theta_c_l0->SetMarkerSize(0.8);
    h_theta_c_l0->GetXaxis()->SetTitle("#theta");
    h_theta_c_l0->SetTitle("#theta: W1 #rightarrow c #q_0");
    h_theta_c_l0->Rebin(2);
    h_theta_c_l0->Scale(1.0 / h_theta_c_l0->Integral());

    h_theta_c_l0_noCR->Draw("HIST P same");
    h_theta_c_l0_noCR->SetMarkerStyle(21);
    h_theta_c_l0_noCR->SetMarkerColor(kRed);
    h_theta_c_l0_noCR->SetMarkerSize(0.7);
    h_theta_c_l0_noCR->Rebin(2);
    h_theta_c_l0_noCR->Scale(1.0 / h_theta_c_l0_noCR->Integral());

    c_theta_reco->cd(2);
    h_theta_l1_l2->Draw("HIST P");
    h_theta_l1_l2->SetMarkerStyle(20);
    h_theta_l1_l2->SetMarkerColor(kBlue);
    h_theta_l1_l2->SetMarkerSize(0.8);
    h_theta_l1_l2->GetXaxis()->SetTitle("#theta");
    h_theta_l1_l2->SetTitle("#theta: W2 #rightarrow #q_1 #q_2");
    h_theta_l1_l2->Rebin(2);
    h_theta_l1_l2->Scale(1.0 / h_theta_l1_l2->Integral());

    h_theta_l1_l2_noCR->Draw("HIST P same");
    h_theta_l1_l2_noCR->SetMarkerStyle(21);
    h_theta_l1_l2_noCR->SetMarkerColor(kRed);
    h_theta_l1_l2_noCR->SetMarkerSize(0.7);
    h_theta_l1_l2_noCR->Rebin(2);
    h_theta_l1_l2_noCR->Scale(1.0 / h_theta_l1_l2_noCR->Integral());

    // // *************************** DELTA PHI DISTRIBUTIONS ***************************
    // TCanvas* c_delta_phi_truth = new TCanvas("c_delta_phi_truth", "Truth Histograms", 1500, 600);
    // c_delta_phi_truth->Divide(2, 1);

    // c_delta_phi_truth->cd(1);
    // h_delta_phi_truth_cs->Draw("HIST");
    // h_delta_phi_truth_cs->SetLineColor(kBlue);
    // h_delta_phi_truth_cs->GetXaxis()->SetTitle("#Delta#phi");
    // h_delta_phi_truth_cs->SetTitle("#Delta#phi: cs");

    // h_delta_phi_truth_cs_noCR->Draw("HIST same");
    // h_delta_phi_truth_cs_noCR->SetLineColor(kRed);

    // TLegend *legend_w1_truth_delta_phi = new TLegend(0.6, 0.7, 0.8, 0.8);
    // legend_w1_truth_delta_phi->SetTextSize(0.04);
    // legend_w1_truth_delta_phi->SetBorderSize(0);
    // legend_w1_truth_delta_phi->SetFillStyle(0); // Transparent background

    // legend_w1_truth_delta_phi->AddEntry(h_delta_phi_truth_cs, "#Delta#phi: cs (CR)", "l");
    // legend_w1_truth_delta_phi->AddEntry(h_delta_phi_truth_cs_noCR, "#Delta#phi: cs (no CR)", "l");
    // legend_w1_truth_delta_phi->Draw();

    // c_delta_phi_truth->cd(2);
    // h_delta_phi_truth_ud->Draw("HIST");
    // h_delta_phi_truth_ud->SetLineColor(kBlue);
    // h_delta_phi_truth_ud->GetXaxis()->SetTitle("#Delta#phi");
    // h_delta_phi_truth_ud->SetTitle("#Delta#phi: ud");

    // h_delta_phi_truth_ud_noCR->Draw("HIST same");
    // h_delta_phi_truth_ud_noCR->SetLineColor(kRed);

    // TLegend *legend_w2_truth_delta_phi = new TLegend(0.6, 0.7, 0.8, 0.8);
    // legend_w2_truth_delta_phi->SetTextSize(0.04);
    // legend_w2_truth_delta_phi->SetBorderSize(0);
    // legend_w2_truth_delta_phi->SetFillStyle(0); // Transparent background

    // legend_w2_truth_delta_phi->AddEntry(h_delta_phi_truth_ud, "#Delta#phi: ud (CR)", "l");
    // legend_w2_truth_delta_phi->AddEntry(h_delta_phi_truth_ud_noCR, "#Delta#phi: ud (no CR)", "l");
    // legend_w2_truth_delta_phi->Draw();


    // TCanvas* c_delta_phi_reco = new TCanvas("c_delta_phi_reco", "Reco Histograms", 1500, 600);
    // c_delta_phi_reco->Divide(2, 1);

    // c_delta_phi_reco->cd(1);
    // h_delta_phi_c_l0->Draw("HIST");
    // h_delta_phi_c_l0->SetLineColor(kBlue);
    // h_delta_phi_c_l0->GetXaxis()->SetTitle("#Delta#phi");
    // h_delta_phi_c_l0->SetTitle("#Delta#phi: W1 #rightarrow c q0");

    // h_delta_phi_c_l0_noCR->Draw("HIST same");
    // h_delta_phi_c_l0_noCR->SetLineColor(kRed);

    // TLegend *legend_w1_reco_delta_phi = new TLegend(0.6, 0.7, 0.8, 0.8);
    // legend_w1_reco_delta_phi->SetTextSize(0.04);
    // legend_w1_reco_delta_phi->SetBorderSize(0);
    // legend_w1_reco_delta_phi->SetFillStyle(0); // Transparent background

    // legend_w1_reco_delta_phi->AddEntry(h_delta_phi_c_l0, "#Delta#phi (CR)", "l");
    // legend_w1_reco_delta_phi->AddEntry(h_delta_phi_c_l0_noCR, "#Delta#phi (no CR)", "l");
    // legend_w1_reco_delta_phi->Draw();

    // c_delta_phi_reco->cd(2);
    // h_delta_phi_l1_l2->Draw("HIST");
    // h_delta_phi_l1_l2->SetLineColor(kBlue);
    // h_delta_phi_l1_l2->GetXaxis()->SetTitle("#Delta#phi");
    // h_delta_phi_l1_l2->SetTitle("#Delta#phi: W2 #rightarrow q1 q2");

    // h_delta_phi_l1_l2_noCR->Draw("HIST same");
    // h_delta_phi_l1_l2_noCR->SetLineColor(kRed);

    // TLegend *legend_w2_reco_delta_phi = new TLegend(0.6, 0.7, 0.8, 0.8);
    // legend_w2_reco_delta_phi->SetTextSize(0.04);
    // legend_w2_reco_delta_phi->SetBorderSize(0);
    // legend_w2_reco_delta_phi->SetFillStyle(0); // Transparent background

    // legend_w2_reco_delta_phi->AddEntry(h_delta_phi_l1_l2, "#Delta#phi (CR)", "l");
    // legend_w2_reco_delta_phi->AddEntry(h_delta_phi_l1_l2_noCR, "#Delta#phi (no CR)", "l");
    // legend_w2_reco_delta_phi->Draw();

    // // ****************************** COS(PHI) DISTRIBUTIONS ******************************
    // TCanvas* c_cos_phi_truth = new TCanvas("c_cos_phi_truth", "Truth Histograms", 1500, 600);
    // c_cos_phi_truth->Divide(2, 1);

    // c_cos_phi_truth->cd(1);
    // h_cos_phi_truth_cs->Draw("HIST");
    // h_cos_phi_truth_cs->SetLineColor(kBlue);
    // h_cos_phi_truth_cs->GetXaxis()->SetTitle("cos(#phi)");
    // h_cos_phi_truth_cs->SetTitle("cos(#phi): cs");

    // h_cos_phi_truth_cs_noCR->Draw("HIST same");
    // h_cos_phi_truth_cs_noCR->SetLineColor(kRed);

    // TLegend *legend_w1_truth_cos_phi = new TLegend(0.6, 0.7, 0.8, 0.8);
    // legend_w1_truth_cos_phi->SetTextSize(0.04);
    // legend_w1_truth_cos_phi->SetBorderSize(0);
    // legend_w1_truth_cos_phi->SetFillStyle(0); // Transparent background

    // legend_w1_truth_cos_phi->AddEntry(h_cos_phi_truth_cs, "cos(#phi): cs (CR)", "l");
    // legend_w1_truth_cos_phi->AddEntry(h_cos_phi_truth_cs_noCR, "cos(#phi): cs (no CR)", "l");
    // legend_w1_truth_cos_phi->Draw();

    // c_cos_phi_truth->cd(2);
    // h_cos_phi_truth_ud->Draw("HIST");
    // h_cos_phi_truth_ud->SetLineColor(kBlue);
    // h_cos_phi_truth_ud->GetXaxis()->SetTitle("cos(#phi)");
    // h_cos_phi_truth_ud->SetTitle("cos(#phi): ud");

    // h_cos_phi_truth_ud_noCR->Draw("HIST same");
    // h_cos_phi_truth_ud_noCR->SetLineColor(kRed);

    // TLegend *legend_w2_truth_cos_phi = new TLegend(0.6, 0.7, 0.8, 0.8);
    // legend_w2_truth_cos_phi->SetTextSize(0.04);
    // legend_w2_truth_cos_phi->SetBorderSize(0);
    // legend_w2_truth_cos_phi->SetFillStyle(0); // Transparent background

    // legend_w2_truth_cos_phi->AddEntry(h_delta_phi_truth_ud, "cos(#phi): ud (CR)", "l");
    // legend_w2_truth_cos_phi->AddEntry(h_delta_phi_truth_ud_noCR, "cos(#phi): ud (no CR)", "l");
    // legend_w2_truth_cos_phi->Draw();


    // TCanvas* c_cos_phi_reco = new TCanvas("c_cos_phi_reco", "Reco Histograms", 1500, 600);
    // c_cos_phi_reco->Divide(2, 1);

    // c_cos_phi_reco->cd(1);
    // h_cos_phi_c_l0->Draw("HIST");
    // h_cos_phi_c_l0->SetLineColor(kBlue);
    // h_cos_phi_c_l0->GetXaxis()->SetTitle("cos(#phi)");
    // h_cos_phi_c_l0->SetTitle("cos(#phi): W1 #rightarrow c q0");

    // h_cos_phi_c_l0_noCR->Draw("HIST same");
    // h_cos_phi_c_l0_noCR->SetLineColor(kRed);

    // TLegend *legend_w1_reco_cos_phi = new TLegend(0.6, 0.7, 0.8, 0.8);
    // legend_w1_reco_cos_phi->SetTextSize(0.04);
    // legend_w1_reco_cos_phi->SetBorderSize(0);
    // legend_w1_reco_cos_phi->SetFillStyle(0); // Transparent background

    // legend_w1_reco_cos_phi->AddEntry(h_cos_phi_c_l0, "cos(#phi) (CR)", "l");
    // legend_w1_reco_cos_phi->AddEntry(h_cos_phi_c_l0_noCR, "cos(#phi) (no CR)", "l");
    // legend_w1_reco_cos_phi->Draw();

    // c_cos_phi_reco->cd(2);
    // h_cos_phi_l1_l2->Draw("HIST");
    // h_cos_phi_l1_l2->SetLineColor(kBlue);
    // h_cos_phi_l1_l2->GetXaxis()->SetTitle("cos(#phi)");
    // h_cos_phi_l1_l2->SetTitle("cos(#phi): W2 #rightarrow q1 q2");

    // h_cos_phi_l1_l2_noCR->Draw("HIST same");
    // h_cos_phi_l1_l2_noCR->SetLineColor(kRed);

    // TLegend *legend_w2_cos_delta_phi = new TLegend(0.6, 0.7, 0.8, 0.8);
    // legend_w2_cos_delta_phi->SetTextSize(0.04);
    // legend_w2_cos_delta_phi->SetBorderSize(0);
    // legend_w2_cos_delta_phi->SetFillStyle(0); // Transparent background

    // legend_w2_cos_delta_phi->AddEntry(h_delta_phi_l1_l2, "cos(#phi) (CR)", "l");
    // legend_w2_cos_delta_phi->AddEntry(h_delta_phi_l1_l2_noCR, "cos(#phi) (no CR)", "l");
    // legend_w2_cos_delta_phi->Draw();

    // // ****************************** E-E CORRELATORS DISTRIBUTIONS ******************************
    // TCanvas* c_eec_truth = new TCanvas("c_eec_truth", "Truth Histograms", 1500, 600);
    // c_eec_truth->Divide(2, 1);

    // c_eec_truth->cd(1);
    // h_eec_truth_cs->Draw("HIST");
    // h_eec_truth_cs->SetLineColor(kBlue);
    // h_eec_truth_cs->GetXaxis()->SetTitle("1 - cos(#phi)");
    // h_eec_truth_cs->SetTitle("e-e correlator: cs");
    // h_eec_truth_cs->Rebin(2);
    // h_eec_truth_cs->Scale(1.0 / h_eec_truth_cs->Integral());

    // h_eec_truth_cs_noCR->Draw("HIST same");
    // h_eec_truth_cs_noCR->SetLineColor(kRed);
    // h_eec_truth_cs_noCR->Rebin(2);
    // h_eec_truth_cs_noCR->Scale(1.0 / h_eec_truth_cs_noCR->Integral());

    // TLegend *legend_w1_truth_eec = new TLegend(0.6, 0.7, 0.8, 0.8);
    // legend_w1_truth_eec->SetTextSize(0.04);
    // legend_w1_truth_eec->SetBorderSize(0);
    // legend_w1_truth_eec->SetFillStyle(0); // Transparent background

    // legend_w1_truth_eec->AddEntry(h_eec_truth_cs, "w/ CR", "l");
    // legend_w1_truth_eec->AddEntry(h_eec_truth_cs_noCR, "no CR", "l");
    // legend_w1_truth_eec->Draw();

    // c_eec_truth->cd(2);
    // h_eec_truth_ud->Draw("HIST");
    // h_eec_truth_ud->SetLineColor(kBlue);
    // h_eec_truth_ud->GetXaxis()->SetTitle("1 - cos(#phi)");
    // h_eec_truth_ud->SetTitle("e-e correlator: ud");
    // h_eec_truth_ud->Rebin(2);
    // h_eec_truth_ud->Scale(1.0 / h_eec_truth_ud->Integral());

    // h_eec_truth_ud_noCR->Draw("HIST same");
    // h_eec_truth_ud_noCR->SetLineColor(kRed);
    // h_eec_truth_ud_noCR->Rebin(2);
    // h_eec_truth_ud_noCR->Scale(1.0 / h_eec_truth_ud_noCR->Integral());

    // TLegend *legend_w2_truth_eec = new TLegend(0.6, 0.7, 0.8, 0.8);
    // legend_w2_truth_eec->SetTextSize(0.04);
    // legend_w2_truth_eec->SetBorderSize(0);
    // legend_w2_truth_eec->SetFillStyle(0); // Transparent background

    // legend_w2_truth_eec->AddEntry(h_delta_phi_truth_ud, "w/ CR", "l");
    // legend_w2_truth_eec->AddEntry(h_delta_phi_truth_ud_noCR, "no CR", "l");
    // legend_w2_truth_eec->Draw();

    // TCanvas* c_eec_truth_ratio = new TCanvas("c_eec_truth_ratio", "Truth Ratio", 1500, 600);

    // c_eec_truth_ratio->cd();
    // TH1D* h_CR_noCR_cs = (TH1D*)h_eec_truth_cs->Clone("h_CR_noCR_cs");
    // h_CR_noCR_cs->Divide(h_eec_truth_cs_noCR);

    // TH1D* h_CR_noCR_ud = (TH1D*)h_eec_truth_ud->Clone("h_CR_noCR_ud");
    // h_CR_noCR_ud->Divide(h_eec_truth_ud_noCR);

    // h_CR_noCR_cs->Draw("HIST");
    // h_CR_noCR_cs->SetLineColor(kRed);
    // h_CR_noCR_cs->GetXaxis()->SetTitle("1 - cos(#phi)");
    // h_CR_noCR_cs->GetYaxis()->SetTitle("Ratio: no CR / CR");

    // h_CR_noCR_ud->Draw("HIST same");
    // h_CR_noCR_ud->SetLineColor(kBlue);

    // TLegend *legend_CR_noCR_truth_ratio = new TLegend(0.6, 0.7, 0.8, 0.8);
    // legend_CR_noCR_truth_ratio->SetTextSize(0.04);
    // legend_CR_noCR_truth_ratio->SetBorderSize(0);
    // legend_CR_noCR_truth_ratio->SetFillStyle(0); // Transparent background

    // legend_CR_noCR_truth_ratio->AddEntry(h_CR_noCR_cs, "w/ CR", "l");
    // legend_CR_noCR_truth_ratio->AddEntry(h_CR_noCR_ud, "no CR", "l");
    // legend_CR_noCR_truth_ratio->Draw();



    // TCanvas* c_eec_reco = new TCanvas("c_eec_reco", "Reco Histograms", 1500, 600);
    // c_eec_reco->Divide(2, 1);

    // c_eec_reco->cd(1);
    // h_eec_c_l0->Draw("HIST");
    // h_eec_c_l0->SetLineColor(kBlue);
    // h_eec_c_l0->GetXaxis()->SetTitle("1 - cos(#phi)");
    // h_eec_c_l0->SetTitle("e-e correlator: W1 #rightarrow c q0");
    // h_eec_c_l0->Rebin(2);
    // h_eec_c_l0->Scale(1.0 / h_eec_c_l0->Integral());

    // h_eec_c_l0_noCR->Draw("HIST same");
    // h_eec_c_l0_noCR->SetLineColor(kRed);
    // h_eec_c_l0_noCR->Rebin(2);
    // h_eec_c_l0_noCR->Scale(1.0 / h_eec_c_l0_noCR->Integral());

    // TLegend *legend_w1_reco_eec = new TLegend(0.6, 0.7, 0.8, 0.8);
    // legend_w1_reco_eec->SetTextSize(0.04);
    // legend_w1_reco_eec->SetBorderSize(0);
    // legend_w1_reco_eec->SetFillStyle(0); // Transparent background

    // legend_w1_reco_eec->AddEntry(h_eec_c_l0, "w/ CR", "l");
    // legend_w1_reco_eec->AddEntry(h_eec_c_l0_noCR, "no CR", "l");
    // legend_w1_reco_eec->Draw();

    // c_eec_reco->cd(2);
    // h_eec_l1_l2->Draw("HIST");
    // h_eec_l1_l2->SetLineColor(kBlue);
    // h_eec_l1_l2->GetXaxis()->SetTitle("1 - cos(#phi)");
    // h_eec_l1_l2->SetTitle("e-e correlator: W2 #rightarrow q1 q2");
    // h_eec_l1_l2->Rebin(2);
    // h_eec_l1_l2->Scale(1.0 / h_eec_l1_l2->Integral());

    // h_eec_l1_l2_noCR->Draw("HIST same");
    // h_eec_l1_l2_noCR->SetLineColor(kRed);
    // h_eec_l1_l2_noCR->Rebin(2);
    // h_eec_l1_l2_noCR->Scale(1.0 / h_eec_l1_l2_noCR->Integral());

    // TLegend *legend_w2_reco_eec = new TLegend(0.6, 0.7, 0.8, 0.8);
    // legend_w2_reco_eec->SetTextSize(0.04);
    // legend_w2_reco_eec->SetBorderSize(0);
    // legend_w2_reco_eec->SetFillStyle(0); // Transparent background

    // legend_w2_reco_eec->AddEntry(h_eec_l1_l2, "w/ CR", "l");
    // legend_w2_reco_eec->AddEntry(h_eec_l1_l2_noCR, "no CR", "l");
    // legend_w2_reco_eec->Draw();

    // TCanvas* c_eec_reco_ratio = new TCanvas("c_eec_reco_ratio", "Reco Ratio", 1500, 600);

    // c_eec_reco_ratio->cd();
    // TH1D* h_CR_noCR_c_l0 = (TH1D*)h_eec_c_l0->Clone("h_CR_noCR_c_l0");
    // h_CR_noCR_c_l0->Divide(h_eec_c_l0_noCR);

    // TH1D* h_CR_noCR_l1_l2 = (TH1D*)h_eec_l1_l2->Clone("h_CR_noCR_l1_l2");
    // h_CR_noCR_l1_l2->Divide(h_eec_l1_l2_noCR);

    // h_CR_noCR_c_l0->Draw("HIST");
    // h_CR_noCR_c_l0->SetLineColor(kRed);
    // h_CR_noCR_c_l0->GetXaxis()->SetTitle("1 - cos(#phi)");
    // h_CR_noCR_c_l0->GetYaxis()->SetTitle("Ratio: no CR / CR");

    // h_CR_noCR_l1_l2->Draw("HIST same");
    // h_CR_noCR_l1_l2->SetLineColor(kBlue);

    // TLegend *legend_CR_noCR_reco_ratio = new TLegend(0.6, 0.7, 0.8, 0.8);
    // legend_CR_noCR_reco_ratio->SetTextSize(0.04);
    // legend_CR_noCR_reco_ratio->SetBorderSize(0);
    // legend_CR_noCR_reco_ratio->SetFillStyle(0); // Transparent background

    // legend_CR_noCR_reco_ratio->AddEntry(h_CR_noCR_c_l0, "w/ CR", "l");
    // legend_CR_noCR_reco_ratio->AddEntry(h_CR_noCR_l1_l2, "no CR", "l");
    // legend_CR_noCR_reco_ratio->Draw();

    c_theta_truth->Write();
    c_theta_reco->Write();
    // c_delta_phi_truth->Write();
    // c_delta_phi_reco->Write();
    // c_cos_phi_truth->Write();
    // c_cos_phi_reco->Write();
    // c_eec_truth->Write();
    // c_eec_truth_ratio->Write();
    // c_eec_reco->Write();
    // c_eec_reco_ratio->Write();

    Plots_WWCR->Write();
    Plots_WWCR->Close();
}