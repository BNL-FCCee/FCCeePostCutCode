void plot() {
    gStyle->SetOptStat(0);

    // Custom binning to match z = (1 - cos(θ)) / 2 axis
    int nBins = 100;
    double xMin = 0.01, xMax = 3.1;
    double logMin = std::log10(xMin);
    double logMax = std::log10(xMax);
    double binWidth = (logMax - logMin) / nBins;
    double* logBins = new double[nBins + 1];
    for (int i = 0; i <= nBins; ++i) logBins[i] = std::pow(10, logMin + i * binWidth);

    // const int nLowBins = 20;
    // const int nHighBins = 20;
    // const int nMidBins = 10;
    // const int nBins = nLowBins + nMidBins + nHighBins;
    // double* bins = new double[nBins + 1];

    // const double lowStart = 1e-4, lowEnd = 1e-2;
    // const double highStart = 1 - 1e-2, highEnd = 1 - 1e-4;
    // const double midStart = lowEnd, midEnd = highStart;

    // double logLowStart = TMath::Log10(lowStart);
    // double logLowEnd = TMath::Log10(lowEnd);
    // for (int i = 0; i <= nLowBins; ++i)
    //     bins[i] = TMath::Power(10, logLowStart + i * (logLowEnd - logLowStart) / nLowBins);

    // for (int i = 1; i <= nMidBins; ++i)
    //     bins[nLowBins + i] = midStart + i * (midEnd - midStart) / nMidBins;

    // double logHighStart = TMath::Log10(1 - highStart);
    // double logHighEnd = TMath::Log10(1 - highEnd);
    // for (int i = 0; i <= nHighBins; ++i)
    //     bins[nLowBins + nMidBins + i] = 1 - TMath::Power(10, logHighStart + i * (logHighEnd - logHighStart) / nHighBins);

    // Load trees
    TFile *fNom = TFile::Open("test.root");
    TFile *fCR = TFile::Open("test_noCR.root");
    TFile *fSKI = TFile::Open("test_SKI.root");
    TFile *fSKII = TFile::Open("test_SKII.root");

    TTree *tNom = (TTree*)fNom->Get("analysis");
    TTree *tCR = (TTree*)fCR->Get("analysis");
    TTree *tSKI = (TTree*)fSKI->Get("analysis");
    TTree *tSKII = (TTree*)fSKII->Get("analysis");

    // Histograms
    TH1F *h_theta_Nom = new TH1F("h_theta_Nom", "#Delta#theta;#theta;Normalized Events", nBins, logBins);
    TH1F *h_theta_CR = new TH1F("h_theta_CR",  "", nBins, logBins);
    TH1F *h_theta_SKI = new TH1F("h_theta_SKI", "", nBins, logBins);
    TH1F *h_theta_SKII = new TH1F("h_theta_SKII", "", nBins, logBins);

    TH1F *h_phi_Nom = new TH1F("h_phi_Nom", "#Delta#phi;Normalized Events", nBins, logBins);
    TH1F *h_phi_CR = new TH1F("h_phi_CR",  "", nBins, logBins);
    TH1F *h_phi_SKI = new TH1F("h_phi_SKI", "", nBins, logBins);
    TH1F *h_phi_SKII = new TH1F("h_phi_SKII", "", nBins, logBins);

    TH1F *h_chi_Nom = new TH1F("h_chi_Nom", "#Delta#chi;Normalized Events", nBins, logBins);
    TH1F *h_chi_CR = new TH1F("h_chi_CR",  "", nBins, logBins);
    TH1F *h_chi_SKI = new TH1F("h_chi_SKI", "", nBins, logBins);
    TH1F *h_chi_SKII = new TH1F("h_chi_SKII", "", nBins, logBins);

    TH1F *h_eec_theta_Nom = new TH1F("h_eec_theta_Nom", "Change Title;Normalized Events", nBins, logBins);
    TH1F *h_eec_theta_CR = new TH1F("h_eec_theta_CR",  "", nBins, logBins);
    TH1F *h_eec_theta_SKI = new TH1F("h_eec_theta_SKI", "", nBins, logBins);
    TH1F *h_eec_theta_SKII = new TH1F("h_eec_theta_SKII", "", nBins, logBins);

    TH1D* h_cos_theta_Nom = new TH1D("h_cos_theta_Nom", "Cos(theta) - Nominal", 100, -1.0, 1.0);
    TH1D* h_cos_theta_CR = new TH1D("h_cos_theta_CR", "Cos(theta) - noCR", 100, -1.0, 1.0);
    TH1D* h_cos_theta_SKI = new TH1D("h_cos_theta_SKI", "Cos(theta) - SKI", 100, -1.0, 1.0);
    TH1D* h_cos_theta_SKII = new TH1D("h_cos_theta_SKII", "Cos(theta) - SKII", 100, -1.0, 1.0);

    // Fill
    tNom->Draw("b_theta_corr >> h_theta_Nom", "", "goff");
    tCR->Draw("b_theta_corr >> h_theta_CR",  "", "goff");
    tSKI->Draw("b_theta_corr >> h_theta_SKI", "", "goff");
    tSKII->Draw("b_theta_corr >> h_theta_SKII", "", "goff");

    tNom->Draw("b_phi_corr >> h_phi_Nom", "", "goff");
    tCR->Draw("b_phi_corr >> h_phi_CR",  "", "goff");
    tSKI->Draw("b_phi_corr >> h_phi_SKI", "", "goff");
    tSKII->Draw("b_phi_corr >> h_phi_SKII", "", "goff");

    tNom->Draw("b_chi_corr >> h_chi_Nom", "", "goff");
    tCR->Draw("b_chi_corr >> h_chi_CR",  "", "goff");
    tSKI->Draw("b_chi_corr >> h_chi_SKI", "", "goff");
    tSKII->Draw("b_chi_corr >> h_chi_SKII", "", "goff");

    // Normalize
    h_theta_Nom->Scale(1.0 / h_theta_Nom->Integral());
    h_theta_CR->Scale(1.0 / h_theta_CR->Integral());
    h_theta_SKI->Scale(1.0 / h_theta_SKI->Integral());
    h_theta_SKII->Scale(1.0 / h_theta_SKII->Integral());

    h_phi_Nom->Scale(1.0 / h_phi_Nom->Integral());
    h_phi_CR->Scale(1.0 / h_phi_CR->Integral());
    h_phi_SKI->Scale(1.0 / h_phi_SKI->Integral());
    h_phi_SKII->Scale(1.0 / h_phi_SKII->Integral());

    h_chi_Nom->Scale(1.0 / h_chi_Nom->Integral());
    h_chi_CR->Scale(1.0 / h_chi_CR->Integral());
    h_chi_SKI->Scale(1.0 / h_chi_SKI->Integral());
    h_chi_SKII->Scale(1.0 / h_chi_SKII->Integral());

    // Ratio to Nominal
    TH1F *r_theta_CR = (TH1F*)h_theta_CR->Clone("r_theta_CR");   r_theta_CR->Divide(h_theta_Nom);
    TH1F *r_theta_SKI = (TH1F*)h_theta_SKI->Clone("r_theta_SKI"); r_theta_SKI->Divide(h_theta_Nom);
    TH1F *r_theta_SKII = (TH1F*)h_theta_SKII->Clone("r_theta_SKII"); r_theta_SKII->Divide(h_theta_Nom);

    TH1F *r_phi_CR = (TH1F*)h_phi_CR->Clone("r_phi_CR");   r_phi_CR->Divide(h_phi_Nom);
    TH1F *r_phi_SKI = (TH1F*)h_phi_SKI->Clone("r_phi_SKI"); r_phi_SKI->Divide(h_phi_Nom);
    TH1F *r_phi_SKII = (TH1F*)h_phi_SKII->Clone("r_phi_SKII"); r_phi_SKII->Divide(h_phi_Nom);

    TH1F *r_chi_CR = (TH1F*)h_chi_CR->Clone("r_chi_CR");   r_chi_CR->Divide(h_chi_Nom);
    TH1F *r_chi_SKI = (TH1F*)h_chi_SKI->Clone("r_chi_SKI"); r_chi_SKI->Divide(h_chi_Nom);
    TH1F *r_chi_SKII = (TH1F*)h_chi_SKII->Clone("r_chi_SKII"); r_chi_SKII->Divide(h_chi_Nom);

    // Canvas 1
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
    h_theta_Nom->GetYaxis()->SetTitleOffset(1.0);

    h_theta_Nom->SetLineColor(kBlack); h_theta_Nom->SetLineWidth(3);
    h_theta_CR->SetLineColor(kRed);   h_theta_CR->SetLineWidth(3);
    h_theta_SKI->SetLineColor(kBlue); h_theta_SKI->SetLineWidth(3);
    h_theta_SKII->SetLineColor(kGreen+2); h_theta_SKII->SetLineWidth(3);

    h_theta_Nom->Draw("HIST");
    h_theta_CR->Draw("HIST SAME");
    h_theta_SKI->Draw("HIST SAME");
    h_theta_SKII->Draw("HIST SAME");

    auto leg_theta = new TLegend(0.15, 0.75, 0.45, 0.88);
    leg_theta->AddEntry(h_theta_Nom,  "Nominal", "l");
    leg_theta->AddEntry(h_theta_CR,   "noCR",   "l");
    leg_theta->AddEntry(h_theta_SKI,  "SKI",  "l");
    leg_theta->AddEntry(h_theta_SKII, "SKII", "l");
    leg_theta->Draw();

    // Ratio plot
    bottomPad->cd();
    gPad->SetLogx();
    r_theta_CR->SetMinimum(0.85);  r_theta_CR->SetMaximum(1.1);
    r_theta_CR->SetLineColor(kRed);
    // r_theta_CR->SetTitle("Ratio to Nominal;z = (1 - cos(#theta))/2;Ratio");
    r_theta_CR->Draw("HIST");
    r_theta_CR->SetLineWidth(3);
    r_theta_SKI->SetLineColor(kBlue);
    r_theta_SKI->SetLineWidth(3);
    r_theta_SKI->Draw("HIST SAME");
    r_theta_SKII->SetLineColor(kGreen+2);
    r_theta_SKII->SetLineWidth(3);
    r_theta_SKII->Draw("HIST SAME");

    r_theta_CR->GetXaxis()->SetTitle("#theta (degrees)"); // Axis titles
    r_theta_CR->GetYaxis()->SetTitle("ratio");
    r_theta_CR->GetXaxis()->CenterTitle(true); // Center them
    r_theta_CR->GetYaxis()->CenterTitle(true);
    r_theta_CR->GetXaxis()->SetTitleSize(0.07); // Set the title size
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
    h_phi_SKI->SetLineColor(kBlue); h_phi_SKI->SetLineWidth(3);
    h_phi_SKII->SetLineColor(kGreen+2); h_phi_SKII->SetLineWidth(3);

    h_phi_Nom->Draw("HIST");
    h_phi_CR->Draw("HIST SAME");
    h_phi_SKI->Draw("HIST SAME");
    h_phi_SKII->Draw("HIST SAME");

    auto leg_phi = new TLegend(0.15, 0.75, 0.45, 0.88);
    leg_phi->AddEntry(h_phi_Nom,  "Nominal", "l");
    leg_phi->AddEntry(h_phi_CR,   "noCR",   "l");
    leg_phi->AddEntry(h_phi_SKI,  "SKI",  "l");
    leg_phi->AddEntry(h_phi_SKII, "SKII", "l");
    leg_phi->Draw();

    // Ratio plot
    bottomPad_phi->cd();
    gPad->SetLogx();

    r_phi_CR->SetMinimum(0.85);  r_phi_CR->SetMaximum(1.1);
    r_phi_CR->SetLineColor(kRed);
    r_phi_CR->Draw("HIST");
    r_phi_CR->SetLineWidth(3);
    r_phi_SKI->SetLineColor(kBlue);
    r_phi_SKI->SetLineWidth(3);
    r_phi_SKI->Draw("HIST SAME");
    r_phi_SKII->SetLineColor(kGreen+2);
    r_phi_SKII->SetLineWidth(3);
    r_phi_SKII->Draw("HIST SAME");

    r_phi_CR->GetXaxis()->SetTitle("#phi (degrees)"); // Axis titles
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
    h_chi_SKI->SetLineColor(kBlue); h_chi_SKI->SetLineWidth(3);
    h_chi_SKII->SetLineColor(kGreen+2); h_chi_SKII->SetLineWidth(3);

    h_chi_Nom->Draw("HIST");
    h_chi_CR->Draw("HIST SAME");
    h_chi_SKI->Draw("HIST SAME");
    h_chi_SKII->Draw("HIST SAME");

    auto leg_chi = new TLegend(0.15, 0.75, 0.45, 0.88);
    leg_chi->AddEntry(h_chi_Nom,  "Nominal", "l");
    leg_chi->AddEntry(h_chi_CR,   "noCR",   "l");
    leg_chi->AddEntry(h_chi_SKI,  "SKI",  "l");
    leg_chi->AddEntry(h_chi_SKII, "SKII", "l");
    leg_chi->Draw();

    // Ratio plot
    bottomPad_chi->cd();
    gPad->SetLogx();

    r_chi_CR->SetMinimum(0.85);  r_chi_CR->SetMaximum(1.1);
    r_chi_CR->SetLineColor(kRed);
    r_chi_CR->Draw("HIST");
    r_chi_CR->SetLineWidth(3);
    r_chi_SKI->SetLineColor(kBlue);
    r_chi_SKI->SetLineWidth(3);
    r_chi_SKI->Draw("HIST SAME");
    r_chi_SKII->SetLineColor(kGreen+2);
    r_chi_SKII->SetLineWidth(3);
    r_chi_SKII->Draw("HIST SAME");

    r_chi_CR->GetXaxis()->SetTitle("#chi (degrees)"); // Axis titles
    r_chi_CR->GetYaxis()->SetTitle("ratio");
    r_chi_CR->GetXaxis()->CenterTitle(true); // Center them
    r_chi_CR->GetYaxis()->CenterTitle(true);
    r_chi_CR->GetXaxis()->SetTitleSize(0.07); // Set the title size
    r_chi_CR->GetYaxis()->SetTitleSize(0.055);
    r_chi_CR->GetXaxis()->SetLabelSize(0.05); // Set the label size
    r_chi_CR->GetYaxis()->SetLabelSize(0.05);
    r_chi_CR->GetXaxis()->SetTitleOffset(1.0); // Set the title offset to the axis
    r_chi_CR->GetYaxis()->SetTitleOffset(0.55);


    // **************************** Cosine Plots ****************************
    int nbins_2 = h_theta_Nom->GetNbinsX();
    for (int i = 1; i <= nbins_2; ++i) {
        double theta_Nom = h_theta_Nom->GetBinCenter(i);
        double binWidth_Nom = h_theta_Nom->GetBinWidth(i);
        double cos_theta_Nom = std::cos(theta_Nom);

        double theta_CR = h_theta_CR->GetBinCenter(i);
        double binWidth_CR = h_theta_CR->GetBinWidth(i);
        double cos_theta_CR = std::cos(theta_CR);

        double theta_SKI = h_theta_SKI->GetBinCenter(i);
        double binWidth_SKI = h_theta_SKI->GetBinWidth(i);
        double cos_theta_SKI = std::cos(theta_SKI);

        double theta_SKII = h_theta_SKII->GetBinCenter(i);
        double binWidt_SKII = h_theta_SKII->GetBinWidth(i);
        double cos_theta_SKII = std::cos(theta_SKII);

        // Nominal
        double weight = h_theta_Nom->GetBinContent(i) / binWidth_Nom;
        h_cos_theta_Nom->Fill(cos_theta_Nom, weight);

        // noCR
        weight = h_theta_CR->GetBinContent(i) / binWidth_CR;
        h_cos_theta_CR->Fill(cos_theta_CR, weight);

        // SKI
        weight = h_theta_SKI->GetBinContent(i) / binWidth_SKI;
        h_cos_theta_SKI->Fill(cos_theta_SKI, weight);

        // SKII
        weight = h_theta_SKII->GetBinContent(i) / binWidt_SKII;
        h_cos_theta_SKII->Fill(cos_theta_SKII, weight);

    }

    TCanvas *c_cos_theta = new TCanvas("c_cos_theta", "Cos(theta) Comparison", 1000, 800);
    c_cos_theta->cd();

    h_cos_theta_Nom->SetLineColor(kBlack); h_cos_theta_Nom->SetLineWidth(3);
    h_cos_theta_CR->SetLineColor(kRed); h_cos_theta_CR->SetLineWidth(3);
    h_cos_theta_SKI->SetLineColor(kBlue); h_cos_theta_SKI->SetLineWidth(3);
    h_cos_theta_SKII->SetLineColor(kGreen+2); h_cos_theta_SKII->SetLineWidth(3);

    h_cos_theta_Nom->SetTitle("cos(#theta) comparison");
    h_cos_theta_Nom->GetXaxis()->SetTitle("cos(#theta)");
    h_cos_theta_Nom->GetYaxis()->SetTitle("Events");

    h_cos_theta_Nom->Draw("HIST");
    h_cos_theta_CR->Draw("HIST SAME");
    h_cos_theta_SKI->Draw("HIST SAME");
    h_cos_theta_SKII->Draw("HIST SAME");

    auto leg_cos = new TLegend(0.15, 0.75, 0.45, 0.88);
    leg_cos->AddEntry(h_cos_theta_Nom,  "Nominal", "l");
    leg_cos->AddEntry(h_cos_theta_CR,   "noCR",   "l");
    leg_cos->AddEntry(h_cos_theta_SKI,  "SKI",  "l");
    leg_cos->AddEntry(h_cos_theta_SKII, "SKII", "l");
    leg_cos->Draw();

    // **************************** EE Correlation Plots ****************************



    c_theta->SaveAs("b_theta_corr_comparison_full_all.png");
    c_cos_theta->SaveAs("b_cos_theta_corr_comparison_full_all.png");

    c_phi->SaveAs("b_phi_corr_comparison_full_all.png");
    c_chi->SaveAs("b_chi_corr_comparison_full_all.png");
}

