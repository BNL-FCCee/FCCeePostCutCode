void plot() {
    // Turn off stats boxes
    gStyle->SetOptStat(0);

    // Open files and get trees
    TFile *fNom = TFile::Open("test.root");
    TFile *fCR = TFile::Open("test_noCR.root");
    TTree *tNom = (TTree*)fNom->Get("analysis");
    TTree *tCR = (TTree*)fCR->Get("analysis");

     // Custom binning to match z = (1 - cos(θ)) / 2 axis
    int nBins = 100;
    double xMin = 0.01, xMax = 3.1;
    double xMin_lin = 0.01, xMax_lin = 180;
    double xMin_lin_2 = -1.0, xMax_lin_2 = 1.0;

    // Linear bins
    double linearBinWidth = (xMax_lin - xMin_lin) / nBins;
    double* linearBins = new double[nBins + 1];
    for (int i = 0; i <= nBins; ++i)
        linearBins[i] = xMin_lin + i * linearBinWidth;

    double linearBinWidth_2 = (xMax_lin_2 - xMin_lin_2) / nBins;
    double* linearBins_2 = new double[nBins + 1];
    for (int i = 0; i <= nBins; ++i)
        linearBins_2[i] = xMin_lin_2 + i * linearBinWidth_2;


    // Logarithmic Bins
    double logMin = std::log10(xMin);
    double logMax = std::log10(xMax);
    double binWidth = (logMax - logMin) / nBins;
    double* logBins = new double[nBins + 1];
    for (int i = 0; i <= nBins; ++i) logBins[i] = std::pow(10, logMin + i * binWidth);

    const int totalBins = 100;
    const double lastEdgeBeforeMerge = 100.0;
    const double finalEdge = 180.0;

    const int binsBeforeLast = totalBins - 1; // 99 bins before last
    std::vector<double> binEdges;

    // Uniform binning from 0 to 105 GeV in 99 bins
    for (int i = 0; i <= binsBeforeLast; ++i) {
        binEdges.push_back(i * (lastEdgeBeforeMerge / binsBeforeLast));
    }

    // Add final bin edge at 180 GeV
    binEdges.push_back(finalEdge);

    // Convert to array if needed (e.g., for ROOT)
    double* binsArray = new double[binEdges.size()];
    std::copy(binEdges.begin(), binEdges.end(), binsArray);

    std::cout << "Last few bin edges:\n";
    for (size_t i = binEdges.size() - 5; i < binEdges.size(); ++i) {
        std::cout << binEdges[i] << " ";
    }
    std::cout << std::endl;

    // Create histograms
    // Histograms
    TH1F *h_p_jet_const_1_Nom = new TH1F("h_p_jet_const_1_Nom", "Jet constituent momentum", totalBins, binsArray);
    TH1F *h_p_jet_const_1_CR = new TH1F("h_p_jet_const_1_CR", "", totalBins, binsArray);

    TH1F *h_p_jet_const_2_Nom = new TH1F("h_p_jet_const_2_Nom", "Jet constituent momentum", totalBins, binsArray);
    TH1F *h_p_jet_const_2_CR = new TH1F("h_p_jet_const_2_CR", "", totalBins, binsArray);

    TH1F *h_theta_Nom = new TH1F("h_theta_Nom", "#Delta#theta;#theta;Normalized Events", nBins, logBins);
    TH1F *h_theta_CR = new TH1F("h_theta_CR",  "", nBins, logBins);

    TH1F *h_phi_Nom = new TH1F("h_phi_Nom", "#Delta#phi;Normalized Events", nBins, logBins);
    TH1F *h_phi_CR = new TH1F("h_phi_CR",  "", nBins, logBins);

    TH1F *h_chi_Nom = new TH1F("h_chi_Nom", "#Delta#chi;Normalized Events", nBins, logBins);
    TH1F *h_chi_CR = new TH1F("h_chi_CR",  "", nBins, logBins);

    TH1F *h_eec_theta_Nom = new TH1F("h_eec_theta_Nom", "Change Title;Normalized Events", nBins, logBins);
    TH1F *h_eec_theta_CR = new TH1F("h_eec_theta_CR",  "", nBins, logBins);

    TH1D* h_cos_theta_Nom = new TH1D("h_cos_theta_Nom", "Cos(theta) - Nominal", nBins, linearBins_2);
    TH1D* h_cos_theta_CR = new TH1D("h_cos_theta_CR", "Cos(theta) - noCR", nBins, linearBins_2);

    TH1D* h_eec_Nom = new TH1D("h_eec_Nom", "Cos(theta) - Nominal", nBins, linearBins_2);
    TH1D* h_eec_CR = new TH1D("h_eec_CR", "Cos(theta) - noCR", nBins, linearBins_2);

    TH1F* h_p_reco_cJet_Nom = new TH1F("h_p_reco_cJet_Nom", "Momentum", nBins, linearBins);
    TH1F* h_p_reco_cJet_CR = new TH1F("h_p_reco_cJet_CR", "", nBins, linearBins);

    TH1F* h_p_reco_lJet0_Nom = new TH1F("h_p_reco_lJet0_Nom", "Momentum", nBins, linearBins);
    TH1F* h_p_reco_lJet0_CR = new TH1F("h_p_reco_lJet0_CR", "", nBins, linearBins);

    TH1F* h_p_reco_lJet1_Nom = new TH1F("h_p_reco_lJet1_Nom", "Momentum", nBins, linearBins);
    TH1F* h_p_reco_lJet1_CR = new TH1F("h_p_reco_lJet1_CR", "", nBins, linearBins);

    TH1F* h_p_reco_lJet2_Nom = new TH1F("h_p_reco_lJet2_Nom", "Momentum", nBins, linearBins);
    TH1F* h_p_reco_lJet2_CR = new TH1F("h_p_reco_lJet2_CR", "", nBins, linearBins);

    // Fill
    tNom->Draw("b_p_jet_const_1 >> h_p_jet_const_1_Nom", "", "goff");
    tCR->Draw("b_p_jet_const_1 >> h_p_jet_const_1_CR", "", "goff");

    tNom->Draw("b_p_jet_const_2 >> h_p_jet_const_2_Nom", "", "goff");
    tCR->Draw("b_p_jet_const_2 >> h_p_jet_const_2_CR", "", "goff");

    tNom->Draw("b_theta_corr >> h_theta_Nom", "", "goff");
    tCR->Draw("b_theta_corr >> h_theta_CR",  "", "goff");

    tNom->Draw("b_phi_corr >> h_phi_Nom", "", "goff");
    tCR->Draw("b_phi_corr >> h_phi_CR",  "", "goff");

    tNom->Draw("b_chi_corr >> h_chi_Nom", "", "goff");
    tCR->Draw("b_chi_corr >> h_chi_CR",  "", "goff");

    tNom->Draw("b_p_reco_cJet >> h_p_reco_cJet_Nom", "", "goff");
    tCR->Draw("b_p_reco_cJet >> h_p_reco_cJet_CR", "", "goff");

    tNom->Draw("b_p_reco_lJet0 >> h_p_reco_lJet0_Nom", "", "goff");
    tCR->Draw("b_p_reco_lJet0 >> h_p_reco_lJet0_CR", "", "goff");

    tNom->Draw("b_p_reco_lJet1 >> h_p_reco_lJet1_Nom", "", "goff");
    tCR->Draw("b_p_reco_lJet1 >> h_p_reco_lJet1_CR", "", "goff");

    tNom->Draw("b_p_reco_lJet2 >> h_p_reco_lJet2_Nom", "", "goff");
    tCR->Draw("b_p_reco_lJet2 >> h_p_reco_lJet2_CR", "", "goff");

    // Normalize both histograms by their integrals
    h_p_jet_const_1_Nom->Scale(1.0 / h_p_jet_const_1_Nom->Integral());
    h_p_jet_const_1_CR->Scale(1.0 / h_p_jet_const_1_CR->Integral());
    
    h_p_jet_const_2_Nom->Scale(1.0 / h_p_jet_const_2_Nom->Integral());
    h_p_jet_const_2_CR->Scale(1.0 / h_p_jet_const_2_CR->Integral());

    h_theta_Nom->Scale(1.0 / h_theta_Nom->Integral());
    h_theta_CR->Scale(1.0 / h_theta_CR->Integral());

    h_phi_Nom->Scale(1.0 / h_phi_Nom->Integral());
    h_phi_CR->Scale(1.0 / h_phi_CR->Integral());

    h_chi_Nom->Scale(1.0 / h_chi_Nom->Integral());
    h_chi_CR->Scale(1.0 / h_chi_CR->Integral());

    h_p_reco_cJet_Nom->Scale(1.0 / h_p_reco_cJet_Nom->Integral());
    h_p_reco_cJet_CR->Scale(1.0 / h_p_reco_cJet_CR->Integral());

    h_p_reco_lJet0_Nom->Scale(1.0 / h_p_reco_lJet0_Nom->Integral());
    h_p_reco_lJet0_CR->Scale(1.0 / h_p_reco_lJet0_CR->Integral());

    h_p_reco_lJet1_Nom->Scale(1.0 / h_p_reco_lJet1_Nom->Integral());
    h_p_reco_lJet1_CR->Scale(1.0 / h_p_reco_lJet1_CR->Integral());

    h_p_reco_lJet2_Nom->Scale(1.0 / h_p_reco_lJet2_Nom->Integral());
    h_p_reco_lJet2_CR->Scale(1.0 / h_p_reco_lJet2_CR->Integral());

     // Ratio to Nominal
    TH1F *r_p_jet_const_1_CR = (TH1F*)h_p_jet_const_1_CR->Clone("r_p_jet_const_1_CR");   r_p_jet_const_1_CR->Divide(h_p_jet_const_1_Nom);
    TH1F *r_p_jet_const_2_CR = (TH1F*)h_p_jet_const_2_CR->Clone("r_p_jet_const_2_CR");   r_p_jet_const_2_CR->Divide(h_p_jet_const_2_Nom);
    TH1F *r_theta_CR = (TH1F*)h_theta_CR->Clone("r_theta_CR");   r_theta_CR->Divide(h_theta_Nom);
    TH1F *r_phi_CR = (TH1F*)h_phi_CR->Clone("r_phi_CR");   r_phi_CR->Divide(h_phi_Nom);
    TH1F *r_chi_CR = (TH1F*)h_chi_CR->Clone("r_chi_CR");   r_chi_CR->Divide(h_chi_Nom);
    TH1D *r_eec_CR = (TH1D*)h_eec_CR->Clone("r_eec_CR");     r_eec_CR->Divide(h_eec_Nom);
    TH1F *r_p_reco_cJet_CR = (TH1F*)h_p_reco_cJet_CR->Clone("r_p_reco_cJet_CR");   r_p_reco_cJet_CR->Divide(h_p_reco_cJet_Nom);
    TH1F *r_p_reco_lJet0_CR = (TH1F*)h_p_reco_lJet0_CR->Clone("r_p_reco_lJet0_CR");   r_p_reco_lJet0_CR->Divide(h_p_reco_lJet0_Nom);
    TH1F *r_p_reco_lJet1_CR = (TH1F*)h_p_reco_lJet1_CR->Clone("r_p_reco_lJet1_CR");   r_p_reco_lJet1_CR->Divide(h_p_reco_lJet1_Nom);
    TH1F *r_p_reco_lJet2_CR = (TH1F*)h_p_reco_lJet2_CR->Clone("r_p_reco_lJet2_CR");   r_p_reco_lJet2_CR->Divide(h_p_reco_lJet2_Nom);
    
    // **************************************** Jet Constituent Plots ****************************************
    TCanvas *c_p_jet_const_1 = new TCanvas("c_p_jet_const_1", "comparison", 1000, 800);
    c_p_jet_const_1->Divide(1, 2);

    // Top pad (60% height): from y = 0.4 to y = 1.0
    TPad *topPad_p_jetConst_1 = new TPad("topPad_p_jetConst_1", "Top Pad", 0.0, 0.4, 1.0, 1.0);
    // topPad->SetTopMargin(0.1);    // Optional
    topPad_p_jetConst_1->SetBottomMargin(0.025);
    topPad_p_jetConst_1->Draw();

    // Bottom pad (40% height): from y = 0.0 to y = 0.4
    TPad *bottomPad_p_jetConst_1 = new TPad("bottomPad_p_jetConst_1", "Bottom Pad", 0.0, 0.0, 1.0, 0.4);
    bottomPad_p_jetConst_1->SetTopMargin(0.025);    // Small top margin
    bottomPad_p_jetConst_1->SetBottomMargin(0.2); // More space for x-axis label
    bottomPad_p_jetConst_1->Draw();

    // Upper plot
    topPad_p_jetConst_1->cd();
    // gPad->SetLogx();
    gPad->SetLogy();

    h_p_jet_const_1_Nom->SetTitle("Momenta- Jet Constituents");
    h_p_jet_const_1_Nom->GetXaxis()->SetLabelSize(0);
    h_p_jet_const_1_Nom->GetYaxis()->SetTitle("normalized events");
    h_p_jet_const_1_Nom->GetYaxis()->CenterTitle(true);
    h_p_jet_const_1_Nom->GetXaxis()->SetTitleSize(0.08);
    h_p_jet_const_1_Nom->GetXaxis()->SetLabelSize(0);
    h_p_jet_const_1_Nom->GetYaxis()->SetTitleOffset(1.0);

    h_p_jet_const_1_Nom->SetLineColor(kBlue); h_p_jet_const_1_Nom->SetLineWidth(3);
    h_p_jet_const_1_CR->SetLineColor(kRed);   h_p_jet_const_1_CR->SetLineWidth(3);

    h_p_jet_const_1_Nom->Draw("HIST");
    h_p_jet_const_1_CR->Draw("HIST SAME");

    auto leg_p_jetConst_1 = new TLegend(0.75, 0.75, 0.95, 0.88);
    leg_p_jetConst_1->AddEntry(h_p_jet_const_1_Nom,  "Nominal", "l");
    leg_p_jetConst_1->AddEntry(h_p_jet_const_1_CR,   "noCR",   "l");
    leg_p_jetConst_1->Draw();

    // Ratio plot
    bottomPad_p_jetConst_1->cd();
    // gPad->SetLogx();
    // r_p_jet_const_1_CR->SetMinimum(0.85);  r_p_jet_const_1_CR->SetMaximum(1.1);
    r_p_jet_const_1_CR->SetLineColor(kRed);
    r_p_jet_const_1_CR->Draw("HIST");
    r_p_jet_const_1_CR->SetLineWidth(3);

    r_p_jet_const_1_CR->GetXaxis()->SetTitle("momenta (GeV/c)"); // Axis titles
    r_p_jet_const_1_CR->GetYaxis()->SetTitle("ratio");
    r_p_jet_const_1_CR->SetTitle("");
    r_p_jet_const_1_CR->GetXaxis()->CenterTitle(true); // Center them
    r_p_jet_const_1_CR->GetYaxis()->CenterTitle(true);
    r_p_jet_const_1_CR->GetXaxis()->SetTitleSize(0.055); // Set the title size
    r_p_jet_const_1_CR->GetYaxis()->SetTitleSize(0.055);
    r_p_jet_const_1_CR->GetXaxis()->SetLabelSize(0.05); // Set the label size
    r_p_jet_const_1_CR->GetYaxis()->SetLabelSize(0.05);
    r_p_jet_const_1_CR->GetXaxis()->SetTitleOffset(1.0); // Set the title offset to the axis
    r_p_jet_const_1_CR->GetYaxis()->SetTitleOffset(0.55);


    TCanvas *c_p_jet_const_2 = new TCanvas("c_p_jet_const_2", "comparison", 1000, 800);
    c_p_jet_const_2->Divide(1, 2);

    // Top pad (60% height): from y = 0.4 to y = 1.0
    TPad *topPad_p_jetConst_2 = new TPad("topPad_p_jetConst_2", "Top Pad", 0.0, 0.4, 1.0, 1.0);
    topPad_p_jetConst_2->SetBottomMargin(0.025);
    topPad_p_jetConst_2->Draw();

    // Bottom pad (40% height): from y = 0.0 to y = 0.4
    TPad *bottomPad_p_jetConst_2 = new TPad("bottomPad_p_jetConst_2", "Bottom Pad", 0.0, 0.0, 1.0, 0.4);
    bottomPad_p_jetConst_2->SetTopMargin(0.025);    // Small top margin
    bottomPad_p_jetConst_2->SetBottomMargin(0.2); // More space for x-axis label
    bottomPad_p_jetConst_2->Draw();

    // Upper plot
    topPad_p_jetConst_2->cd();
    // gPad->SetLogx();
    gPad->SetLogy();

    h_p_jet_const_2_Nom->SetTitle("Momenta- Jet Constituents");
    h_p_jet_const_2_Nom->GetXaxis()->SetLabelSize(0);
    h_p_jet_const_2_Nom->GetYaxis()->SetTitle("normalized events");
    h_p_jet_const_2_Nom->GetYaxis()->CenterTitle(true);
    h_p_jet_const_2_Nom->GetXaxis()->SetTitleSize(0.08);
    h_p_jet_const_2_Nom->GetXaxis()->SetLabelSize(0);
    h_p_jet_const_2_Nom->GetYaxis()->SetTitleOffset(1.0);

    h_p_jet_const_2_Nom->SetLineColor(kBlue); h_p_jet_const_2_Nom->SetLineWidth(3);
    h_p_jet_const_2_CR->SetLineColor(kRed);   h_p_jet_const_2_CR->SetLineWidth(3);

    h_p_jet_const_2_Nom->Draw("HIST");
    h_p_jet_const_2_CR->Draw("HIST SAME");

    auto leg_p_jetConst_2 = new TLegend(0.75, 0.75, 0.95, 0.88);
    leg_p_jetConst_2->AddEntry(h_p_jet_const_2_Nom,  "Nominal", "l");
    leg_p_jetConst_2->AddEntry(h_p_jet_const_2_CR,   "noCR",   "l");
    leg_p_jetConst_2->Draw();

    // Ratio plot
    bottomPad_p_jetConst_2->cd();
    // gPad->SetLogx();
    // r_p_jet_const_1_CR->SetMinimum(0.85);  r_p_jet_const_1_CR->SetMaximum(1.1);
    r_p_jet_const_2_CR->SetLineColor(kRed);
    r_p_jet_const_2_CR->Draw("HIST");
    r_p_jet_const_2_CR->SetLineWidth(3);

    r_p_jet_const_2_CR->GetXaxis()->SetTitle("momenta (GeV/c)"); // Axis titles
    r_p_jet_const_2_CR->GetYaxis()->SetTitle("ratio");
    r_p_jet_const_2_CR->SetTitle("");
    r_p_jet_const_2_CR->GetXaxis()->CenterTitle(true); // Center them
    r_p_jet_const_2_CR->GetYaxis()->CenterTitle(true);
    r_p_jet_const_2_CR->GetXaxis()->SetTitleSize(0.055); // Set the title size
    r_p_jet_const_2_CR->GetYaxis()->SetTitleSize(0.055);
    r_p_jet_const_2_CR->GetXaxis()->SetLabelSize(0.05); // Set the label size
    r_p_jet_const_2_CR->GetYaxis()->SetLabelSize(0.05);
    r_p_jet_const_2_CR->GetXaxis()->SetTitleOffset(1.0); // Set the title offset to the axis
    r_p_jet_const_2_CR->GetYaxis()->SetTitleOffset(0.55);

    TCanvas *c_theta = new TCanvas("c_theta", "comparison", 1000, 800);
    c_theta->Divide(1, 2);

    // Top pad (60% height): from y = 0.4 to y = 1.0
    TPad *topPad = new TPad("topPad", "Top Pad", 0.0, 0.4, 1.0, 1.0);
    // topPad->SetTopMargin(0.1);    // Optional
    topPad->SetBottomMargin(0.025);
    topPad->Draw();

    // Bottom pad (40% height): from y = 0.0 to y = 0.4
    TPad *bottomPad = new TPad("bottomPad", "Bottom Pad", 0.0, 0.0, 1.0, 0.4);
    bottomPad->SetTopMargin(0.025);    // Small top margin
    bottomPad->SetBottomMargin(0.2); // More space for x-axis label
    bottomPad->Draw();

    // Upper plot
    topPad->cd();
    gPad->SetLogx();
    gPad->SetLogy();

    h_theta_Nom->SetTitle("b_theta_corr");
    h_theta_Nom->GetXaxis()->SetLabelSize(0);
    h_theta_Nom->GetYaxis()->SetTitle("normalized events");
    h_theta_Nom->GetYaxis()->CenterTitle(true);
    h_theta_Nom->GetXaxis()->SetTitleSize(0.08);
    h_theta_Nom->GetXaxis()->SetLabelSize(0.09);
    h_theta_Nom->GetYaxis()->SetTitleOffset(1.0);

    h_theta_Nom->SetLineColor(kBlack); h_theta_Nom->SetLineWidth(3);
    h_theta_CR->SetLineColor(kRed);   h_theta_CR->SetLineWidth(3);

    h_theta_Nom->Draw("HIST");
    h_theta_CR->Draw("HIST SAME");

    auto leg_theta = new TLegend(0.15, 0.3, 0.35, 0.45);
    leg_theta->AddEntry(h_theta_Nom,  "Nominal", "l");
    leg_theta->AddEntry(h_theta_CR,   "noCR",   "l");
    leg_theta->Draw();

    // Ratio plot
    bottomPad->cd();
    gPad->SetLogx();
    r_theta_CR->SetMinimum(0.85);  r_theta_CR->SetMaximum(1.1);
    r_theta_CR->SetLineColor(kRed);
    r_theta_CR->Draw("HIST");
    r_theta_CR->SetLineWidth(3);

    r_theta_CR->GetXaxis()->SetTitle("#theta (radians)"); // Axis titles
    r_theta_CR->GetYaxis()->SetTitle("ratio");
    r_theta_CR->GetXaxis()->CenterTitle(true); // Center them
    r_theta_CR->GetYaxis()->CenterTitle(true);
    r_theta_CR->GetXaxis()->SetTitleSize(0.08); // Set the title size
    r_theta_CR->GetYaxis()->SetTitleSize(0.055);
    r_theta_CR->GetXaxis()->SetLabelSize(0.05); // Set the label size
    r_theta_CR->GetYaxis()->SetLabelSize(0.05);
    r_theta_CR->GetXaxis()->SetTitleOffset(1.0); // Set the title offset to the axis
    r_theta_CR->GetYaxis()->SetTitleOffset(0.55);


    // Canvas 2
    TCanvas *c_phi = new TCanvas("c_phi", "comparison", 1000, 800);
    c_phi->Divide(1, 2);

    // Top pad (60% height): from y = 0.4 to y = 1.0
    TPad *topPad_phi = new TPad("topPad_phi", "Top Pad", 0.0, 0.4, 1.0, 1.0);
    // topPad->SetTopMargin(0.1);    // Optional
    topPad_phi->SetBottomMargin(0.025);
    topPad_phi->Draw();

    // Bottom pad (40% height): from y = 0.0 to y = 0.4
    TPad *bottomPad_phi = new TPad("bottomPad_phi", "Bottom Pad", 0.0, 0.0, 1.0, 0.4);
    bottomPad_phi->SetTopMargin(0.025);    // Small top margin
    bottomPad_phi->SetBottomMargin(0.2); // More space for x-axis label
    bottomPad_phi->Draw();

    // Upper plot
    topPad_phi->cd();
    gPad->SetLogx();
    gPad->SetLogy();

    h_phi_Nom->SetTitle("b_phi_corr");
    h_phi_Nom->GetXaxis()->SetLabelSize(0);
    h_phi_Nom->GetYaxis()->SetTitle("normalized events");
    h_phi_Nom->GetYaxis()->CenterTitle(true);
    h_phi_Nom->GetXaxis()->SetTitleSize(0.08);
    h_phi_Nom->GetYaxis()->SetTitleOffset(1.0);

    h_phi_Nom->SetLineColor(kBlack); h_phi_Nom->SetLineWidth(3);
    h_phi_CR->SetLineColor(kRed);   h_phi_CR->SetLineWidth(3);

    h_phi_Nom->Draw("HIST");
    h_phi_CR->Draw("HIST SAME");

    auto leg_phi = new TLegend(0.15, 0.75, 0.45, 0.88);
    leg_phi->AddEntry(h_phi_Nom,  "Nominal", "l");
    leg_phi->AddEntry(h_phi_CR,   "noCR",   "l");
    leg_phi->Draw();

    // Ratio plot
    bottomPad_phi->cd();
    gPad->SetLogx();

    r_phi_CR->SetMinimum(0.85);  r_phi_CR->SetMaximum(1.1);
    r_phi_CR->SetLineColor(kRed);
    r_phi_CR->Draw("HIST");
    r_phi_CR->SetLineWidth(3);

    r_phi_CR->GetXaxis()->SetTitle("#phi (radians)"); // Axis titles
    r_phi_CR->GetYaxis()->SetTitle("ratio");
    r_phi_CR->GetXaxis()->CenterTitle(true); // Center them
    r_phi_CR->GetYaxis()->CenterTitle(true);
    r_phi_CR->GetXaxis()->SetTitleSize(0.07); // Set the title size
    r_phi_CR->GetYaxis()->SetTitleSize(0.055);
    r_phi_CR->GetXaxis()->SetLabelSize(0.05); // Set the label size
    r_phi_CR->GetYaxis()->SetLabelSize(0.05);
    r_phi_CR->GetXaxis()->SetTitleOffset(1.0); // Set the title offset to the axis
    r_phi_CR->GetYaxis()->SetTitleOffset(0.55);


    // Canvas 3
    TCanvas *c_chi = new TCanvas("c_chi", "comparison", 1000, 800);
    c_chi->Divide(1, 2);

    // Top pad (60% height): from y = 0.4 to y = 1.0
    TPad *topPad_chi = new TPad("topPad_chi", "Top Pad", 0.0, 0.4, 1.0, 1.0);
    // topPad->SetTopMargin(0.1);    // Optional
    topPad_chi->SetBottomMargin(0.025);
    topPad_chi->Draw();

    // Bottom pad (40% height): from y = 0.0 to y = 0.4
    TPad *bottomPad_chi = new TPad("bottomPad_chi", "Bottom Pad", 0.0, 0.0, 1.0, 0.4);
    bottomPad_chi->SetTopMargin(0.025);    // Small top margin
    bottomPad_chi->SetBottomMargin(0.2); // More space for x-axis label
    bottomPad_chi->Draw();

    // Upper plot
    topPad_chi->cd();
    gPad->SetLogx();
    gPad->SetLogy();

    h_chi_Nom->SetTitle("b_chi_corr");
    h_chi_Nom->GetXaxis()->SetLabelSize(0);
    h_chi_Nom->GetYaxis()->SetTitle("normalized events");
    h_chi_Nom->GetYaxis()->CenterTitle(true);
    h_chi_Nom->GetXaxis()->SetTitleSize(0.08);
    h_chi_Nom->GetYaxis()->SetTitleOffset(1.0);

    h_chi_Nom->SetLineColor(kBlack); h_chi_Nom->SetLineWidth(3);
    h_chi_CR->SetLineColor(kRed);   h_chi_CR->SetLineWidth(3);

    h_chi_Nom->Draw("HIST");
    h_chi_CR->Draw("HIST SAME");

    auto leg_chi = new TLegend(0.15, 0.75, 0.45, 0.88);
    leg_chi->AddEntry(h_chi_Nom,  "Nominal", "l");
    leg_chi->AddEntry(h_chi_CR,   "noCR",   "l");
    leg_chi->Draw();

    // Ratio plot
    bottomPad_chi->cd();
    gPad->SetLogx();

    r_chi_CR->SetMinimum(0.85);  r_chi_CR->SetMaximum(1.1);
    r_chi_CR->SetLineColor(kRed);
    r_chi_CR->Draw("HIST");
    r_chi_CR->SetLineWidth(3);

    r_chi_CR->GetXaxis()->SetTitle("#chi (radians)"); // Axis titles
    r_chi_CR->GetYaxis()->SetTitle("ratio");
    r_chi_CR->GetXaxis()->CenterTitle(true); // Center them
    r_chi_CR->GetYaxis()->CenterTitle(true);
    r_chi_CR->GetXaxis()->SetTitleSize(0.07); // Set the title size
    r_chi_CR->GetYaxis()->SetTitleSize(0.055);
    r_chi_CR->GetXaxis()->SetLabelSize(0.05); // Set the label size
    r_chi_CR->GetYaxis()->SetLabelSize(0.05);
    r_chi_CR->GetXaxis()->SetTitleOffset(1.0); // Set the title offset to the axis
    r_chi_CR->GetYaxis()->SetTitleOffset(0.55);


    c_p_jet_const_1->SaveAs("b_p_jet_const_comparison_full_all.png");
    c_p_jet_const_2->SaveAs("b_p_jet2_const_comparison_full_all.png");
    c_theta->SaveAs("b_theta_corr_comparison_full_all.png");
    c_phi->SaveAs("b_phi_corr_comparison_full_all.png");
    c_chi->SaveAs("b_chi_corr_comparison_full_all.png");

}

