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

    TH1F *h_cos_phi_truth_cs_noCR = (TH1F*)test_noCR->Get("Nominal/h_cos_phi_truth_cs_WW");
    TH1F *h_cos_phi_truth_ud_noCR = (TH1F*)test_noCR->Get("Nominal/h_cos_phi_truth_ud_WW");
    TH1F *h_cos_phi_c_l0_noCR = (TH1F*)test_noCR->Get("Nominal/h_cos_phi_c_l0_WW");
    TH1F *h_cos_phi_l1_l2_noCR = (TH1F*)test_noCR->Get("Nominal/h_cos_phi_l1_l2_WW");

    TH1F *h_eec_truth_cs_noCR = (TH1F*)test_noCR->Get("Nominal/h_eec_truth_cs_WW");
    TH1F *h_eec_truth_ud_noCR = (TH1F*)test_noCR->Get("Nominal/h_eec_truth_ud_WW");
    TH1F *h_eec_c_l0_noCR = (TH1F*)test_noCR->Get("Nominal/h_eec_c_l0_WW");
    TH1F *h_eec_l1_l2_noCR = (TH1F*)test_noCR->Get("Nominal/h_eec_l1_l2_WW");

    // *************************** SK-I CR Plots ***************************
    TH1F *h_W1_truth_mass_SKI = (TH1F*)test_SKI->Get("Nominal/h_W1_truth_mass_WW");
    TH1F *h_W2_truth_mass_SKI = (TH1F*)test_SKI->Get("Nominal/h_W2_truth_mass_WW");
    TH1F *h_W1_mass_SKI = (TH1F*)test_SKI->Get("Nominal/h_W1_mass_WW");
    TH1F *h_W2_mass_SKI = (TH1F*)test_SKI->Get("Nominal/h_W2_mass_WW");

    TH1F *h_W1_truth_p_SKI = (TH1F*)test_SKI->Get("Nominal/h_W1_truth_p_WW");
    TH1F *h_W2_truth_p_SKI = (TH1F*)test_SKI->Get("Nominal/h_W2_truth_p_WW");
    TH1F *h_W1_p_SKI = (TH1F*)test_SKI->Get("Nominal/h_W1_p_WW");
    TH1F *h_W2_p_SKI = (TH1F*)test_SKI->Get("Nominal/h_W2_p_WW");

    TH1F *h_W1_truth_e_SKI = (TH1F*)test_SKI->Get("Nominal/h_W1_truth_e_WW");
    TH1F *h_W2_truth_e_SKI = (TH1F*)test_SKI->Get("Nominal/h_W2_truth_e_WW");
    TH1F *h_W1_e_SKI = (TH1F*)test_SKI->Get("Nominal/h_W1_e_WW");
    TH1F *h_W2_e_SKI = (TH1F*)test_SKI->Get("Nominal/h_W2_e_WW");

    TH1F *h_cos_phi_truth_cs_SKI = (TH1F*)test_SKI->Get("Nominal/h_cos_phi_truth_cs_WW");
    TH1F *h_cos_phi_truth_ud_SKI = (TH1F*)test_SKI->Get("Nominal/h_cos_phi_truth_ud_WW");
    TH1F *h_cos_phi_c_l0_SKI = (TH1F*)test_SKI->Get("Nominal/h_cos_phi_c_l0_WW");
    TH1F *h_cos_phi_l1_l2_SKI = (TH1F*)test_SKI->Get("Nominal/h_cos_phi_l1_l2_WW");

    TH1F *h_eec_truth_cs_SKI = (TH1F*)test_SKI->Get("Nominal/h_eec_truth_cs_WW");
    TH1F *h_eec_truth_ud_SKI = (TH1F*)test_SKI->Get("Nominal/h_eec_truth_ud_WW");
    TH1F *h_eec_c_l0_SKI = (TH1F*)test_SKI->Get("Nominal/h_eec_c_l0_WW");
    TH1F *h_eec_l1_l2_SKI = (TH1F*)test_SKI->Get("Nominal/h_eec_l1_l2_WW");

    // *************************** SK-II CR Plots ***************************
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

    TH1F *h_cos_phi_truth_cs = (TH1F*)test_CR->Get("Nominal/h_cos_phi_truth_cs_WW");
    TH1F *h_cos_phi_truth_ud = (TH1F*)test_CR->Get("Nominal/h_cos_phi_truth_ud_WW");
    TH1F *h_cos_phi_c_l0 = (TH1F*)test_CR->Get("Nominal/h_cos_phi_c_l0_WW");
    TH1F *h_cos_phi_l1_l2 = (TH1F*)test_CR->Get("Nominal/h_cos_phi_l1_l2_WW");

    TH1F *h_eec_truth_cs = (TH1F*)test_CR->Get("Nominal/h_eec_truth_cs_WW");
    TH1F *h_eec_truth_ud = (TH1F*)test_CR->Get("Nominal/h_eec_truth_ud_WW");
    TH1F *h_eec_c_l0 = (TH1F*)test_CR->Get("Nominal/h_eec_c_l0_WW");
    TH1F *h_eec_l1_l2 = (TH1F*)test_CR->Get("Nominal/h_eec_l1_l2_WW");


    TFile *Plots_WWCR = new TFile("/usatlas/u/kleibensperger/FCCeePostCutCode/run/WWCR_Plots.root", "Recreate");

    // *************************** MASS PLOTS ***************************
    TCanvas* c_W1_W2_mass = new TCanvas("c_W1_W2_mass", "Truth Histograms", 1500, 600);
    c_W1_W2_mass->Divide(2, 1);
    c_W1_W2_mass->SetWindowSize(1500, 600);      // Window size in pixels
    c_W1_W2_mass->SetCanvasSize(1500, 600);      // Drawable area in pixels

    c_W1_W2_mass->cd(1);
    h_W1_truth_mass->Draw();
    h_W1_truth_mass->SetLineColor(kBlue);
    h_W1_truth_mass->GetXaxis()->SetTitle("GeV");

    h_W1_mass->Draw("same");
    h_W1_mass->SetLineColor(kRed);

    TLegend *legend_w1_mass = new TLegend(0.6, 0.7, 0.8, 0.8);
    legend_w1_mass->SetTextSize(0.04);
    legend_w1_mass->SetBorderSize(0);
    legend_w1_mass->SetFillStyle(0); // Transparent background

    legend_w1_mass->AddEntry(h_W1_truth_mass, "w1_truth_mass", "l");
    legend_w1_mass->AddEntry(h_W1_mass, "w1_mass", "l");
    legend_w1_mass->Draw();

    c_W1_W2_mass->cd(2);
    h_W2_truth_mass->Draw();
    h_W2_truth_mass->SetLineColor(kBlue);
    h_W2_truth_mass->GetXaxis()->SetTitle("GeV");

    h_W2_mass->Draw("same");
    h_W2_mass->SetLineColor(kRed);

    TLegend *legend_w2_mass = new TLegend(0.6, 0.7, 0.8, 0.8);
    legend_w2_mass->SetTextSize(0.04);
    legend_w2_mass->SetBorderSize(0);
    legend_w2_mass->SetFillStyle(0); // Transparent background

    legend_w2_mass->AddEntry(h_W2_truth_mass, "w2_truth_mass", "l");
    legend_w2_mass->AddEntry(h_W2_mass, "w2_mass", "l");
    legend_w2_mass->Draw();

    // *************************** COMPARATIVE MOMENTA PLOTS ***************************
    TCanvas* c_W1_W2_p = new TCanvas("c_W1_W2_p", "Truth Histograms", 1500, 600);
    c_W1_W2_p->Divide(2, 1);

    c_W1_W2_p->cd(1);
    h_W1_truth_p->Draw();
    h_W1_truth_p->SetLineColor(kBlue);
    h_W1_truth_p->GetXaxis()->SetTitle("GeV");

    h_W1_p->Draw("same");
    h_W1_p->SetLineColor(kRed);

    TLegend *legend_w1_p = new TLegend(0.6, 0.7, 0.8, 0.8);
    legend_w1_p->SetTextSize(0.04);
    legend_w1_p->SetBorderSize(0);
    legend_w1_p->SetFillStyle(0); // Transparent background

    legend_w1_p->AddEntry(h_W1_truth_p, "w1_truth_p", "l");
    legend_w1_p->AddEntry(h_W1_p, "w1_p", "l");
    legend_w1_p->Draw();

    c_W1_W2_p->cd(2);
    h_W2_truth_p->Draw();
    h_W2_truth_p->SetLineColor(kBlue);
    h_W2_truth_p->GetXaxis()->SetTitle("GeV");

    h_W2_p->Draw("same");
    h_W2_p->SetLineColor(kRed);

    TLegend *legend_w2_p = new TLegend(0.6, 0.7, 0.8, 0.8);
    legend_w2_p->SetTextSize(0.04);
    legend_w2_p->SetBorderSize(0);
    legend_w2_p->SetFillStyle(0); // Transparent background

    legend_w2_p->AddEntry(h_W2_truth_p, "w2_truth_p", "l");
    legend_w2_p->AddEntry(h_W2_p, "w2_p", "l");
    legend_w2_p->Draw();

    // *************************** COMPARATIVE ENERGY PLOTS ***************************
    TCanvas* c_W1_W2_e = new TCanvas("c_W1_W2_e", "Truth Histograms", 1500, 600);
    c_W1_W2_e->Divide(2, 1);

    c_W1_W2_e->cd(1);
    h_W1_truth_e->Draw();
    h_W1_truth_e->SetLineColor(kBlue);
    h_W1_truth_e->GetXaxis()->SetTitle("GeV");

    h_W1_e->Draw("same");
    h_W1_e->SetLineColor(kRed);

    TLegend *legend_w1_e = new TLegend(0.6, 0.7, 0.8, 0.8);
    legend_w1_e->SetTextSize(0.04);
    legend_w1_e->SetBorderSize(0);
    legend_w1_e->SetFillStyle(0); // Transparent background

    legend_w1_e->AddEntry(h_W1_truth_e, "w1_truth_e", "l");
    legend_w1_e->AddEntry(h_W1_e, "w1_e", "l");
    legend_w1_e->Draw();

    c_W1_W2_e->cd(2);
    h_W2_truth_e->Draw();
    h_W2_truth_e->SetLineColor(kBlue);
    h_W2_truth_e->GetXaxis()->SetTitle("GeV");

    h_W2_e->Draw("same");
    h_W2_e->SetLineColor(kRed);

    TLegend *legend_w2_e = new TLegend(0.6, 0.7, 0.8, 0.8);
    legend_w2_e->SetTextSize(0.04);
    legend_w2_e->SetBorderSize(0);
    legend_w2_e->SetFillStyle(0); // Transparent background

    legend_w2_e->AddEntry(h_W2_truth_e, "w2_truth_e", "l");
    legend_w2_e->AddEntry(h_W2_e, "w2_e", "l");
    legend_w2_e->Draw();

    // *************************** COMPARATIVE ANGULAR DISTRIBUTIONS PLOTS ***************************
    TCanvas* c_angular_dist = new TCanvas("c_angular_dist", "Truth Histograms", 1500, 600);
    c_angular_dist->Divide(2, 1);

    c_angular_dist->cd(1);
    h_cos_phi_truth_cs->Draw();
    h_cos_phi_truth_cs->SetLineColor(kBlue);
    h_cos_phi_truth_cs->GetXaxis()->SetTitle("cos(#phi)");

    h_cos_phi_c_l0->Draw("same");
    h_cos_phi_c_l0->SetLineColor(kRed);

    TLegend *legend_w1_cos_phi = new TLegend(0.6, 0.7, 0.8, 0.8);
    legend_w1_cos_phi->SetTextSize(0.04);
    legend_w1_cos_phi->SetBorderSize(0);
    legend_w1_cos_phi->SetFillStyle(0); // Transparent background

    legend_w1_cos_phi->AddEntry(h_cos_phi_truth_cs, "truth cos_phi_cs", "l");
    legend_w1_cos_phi->AddEntry(h_cos_phi_c_l0, "truth cos_phi_c_l0", "l");
    legend_w1_cos_phi->Draw();

    c_angular_dist->cd(2);
    h_cos_phi_truth_ud->Draw();
    h_cos_phi_truth_ud->SetLineColor(kBlue);
    h_cos_phi_truth_ud->GetXaxis()->SetTitle("cos(#phi)");

    h_cos_phi_l1_l2->Draw("same");
    h_cos_phi_l1_l2->SetLineColor(kRed);

    TLegend *legend_w2_cos_phi = new TLegend(0.6, 0.7, 0.8, 0.8);
    legend_w2_cos_phi->SetTextSize(0.04);
    legend_w2_cos_phi->SetBorderSize(0);
    legend_w2_cos_phi->SetFillStyle(0); // Transparent background

    legend_w2_cos_phi->AddEntry(h_cos_phi_truth_ud, "truth cos_phi_ud", "l");
    legend_w2_cos_phi->AddEntry(h_cos_phi_l1_l2, "truth cos_phi_l1_l2", "l");
    legend_w2_cos_phi->Draw();

    // *************************** Energy-Energy Correlators ***************************
    TCanvas* c_truth_ee_correlators = new TCanvas("c_truth_ee_correlators", "Truth Histograms", 1500, 600);

    c_truth_ee_correlators->cd();
    h_eec_truth_cs->Draw();
    h_eec_truth_cs->SetLineColor(kBlue);
    h_eec_truth_cs->GetXaxis()->SetTitle("cos(#phi)");

    h_eec_truth_ud->Draw("same");
    h_eec_truth_ud->SetLineColor(kRed);

    TLegend *legend_eec_truth = new TLegend(0.6, 0.7, 0.8, 0.8);
    legend_eec_truth->SetTextSize(0.04);
    legend_eec_truth->SetBorderSize(0);
    legend_eec_truth->SetFillStyle(0); // Transparent background

    legend_eec_truth->AddEntry(h_eec_truth_cs, "W1 #rightarrow cs", "l");
    legend_eec_truth->AddEntry(h_eec_truth_ud, "W1 #rightarrow ud", "l");
    legend_eec_truth->Draw();

    TCanvas* c_ee_correlators = new TCanvas("c_ee_correlators", "Truth Histograms", 1500, 600);

    c_ee_correlators->cd();
    h_eec_c_l0->Draw();
    h_eec_c_l0->SetLineColor(kBlue);
    h_eec_c_l0->GetXaxis()->SetTitle("cos(#phi)");

    h_eec_l1_l2->Draw("same");
    h_eec_l1_l2->SetLineColor(kRed);

    TLegend *legend_eec = new TLegend(0.6, 0.7, 0.8, 0.8);
    legend_eec->SetTextSize(0.04);
    legend_eec->SetBorderSize(0);
    legend_eec->SetFillStyle(0); // Transparent background

    legend_eec->AddEntry(h_eec_c_l0, "truth cos_phi_ud", "l");
    legend_eec->AddEntry(h_eec_l1_l2, "truth cos_phi_l1_l2", "l");
    legend_eec->Draw();

    c_W1_W2_mass->Write();
    c_W1_W2_p->Write();
    c_W1_W2_e->Write();
    c_angular_dist->Write();
    c_truth_ee_correlators->Write();
    c_ee_correlators->Write();

    Plots_WWCR->Write();
    Plots_WWCR->Close();
}