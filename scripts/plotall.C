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
    TFile *fCR = TFile::Open("testnoCR.root");
    TFile *fSKI = TFile::Open("testnoSKI.root");
    TFile *fSKII = TFile::Open("testnoSKII.root");

    TTree *tNom = (TTree*)fNom->Get("analysis");
    TTree *tCR = (TTree*)fCR->Get("analysis");
    TTree *tSKI = (TTree*)fSKI->Get("analysis");
    TTree *tSKII = (TTree*)fSKII->Get("analysis");

    // Histograms
    TH1F *hNom = new TH1F("hNom", "delta_theta;Normalized Events", nBins, logBins);
    TH1F *hCR = new TH1F("hCR",  "", nBins, logBins);
    TH1F *hSKI = new TH1F("hSKI", "", nBins, logBins);
    TH1F *hSKII = new TH1F("hSKII", "", nBins, logBins);

    // Fill
    tNom->Draw("b_theta_corr >> hNom", "", "goff");
    tCR->Draw("b_theta_corr >> hCR",  "", "goff");
    tSKI->Draw("b_theta_corr >> hSKI", "", "goff");
    tSKII->Draw("b_theta_corr >> hSKII", "", "goff");

    // Normalize
    hNom->Scale(1.0 / hNom->Integral());
    hCR->Scale(1.0 / hCR->Integral());
    hSKI->Scale(1.0 / hSKI->Integral());
    hSKII->Scale(1.0 / hSKII->Integral());

    // Ratio to Nominal
    TH1F *rCR = (TH1F*)hCR->Clone("rCR");   rCR->Divide(hNom);
    TH1F *rSKI = (TH1F*)hSKI->Clone("rSKI"); rSKI->Divide(hNom);
    TH1F *rSKII = (TH1F*)hSKII->Clone("rSKII"); rSKII->Divide(hNom);

    // Canvas
    TCanvas *c = new TCanvas("c", "comparison", 800, 800);
    c->Divide(1, 2);

    // Upper plot
    c->cd(1);
    gPad->SetLogx();
    gPad->SetLogy();
    hNom->SetLineColor(kBlack); hNom->SetLineWidth(3);
    hCR->SetLineColor(kRed);   hCR->SetLineWidth(3);
    hSKI->SetLineColor(kBlue); hSKI->SetLineWidth(3);
    hSKII->SetLineColor(kGreen+2); hSKII->SetLineWidth(3);

    hNom->Draw("HIST");
    hCR->Draw("HIST SAME");
    hSKI->Draw("HIST SAME");
    hSKII->Draw("HIST SAME");

    auto leg = new TLegend(0.15, 0.75, 0.45, 0.88);
    leg->AddEntry(hNom,  "Nominal", "l");
    leg->AddEntry(hCR,   "noCR",   "l");
    leg->AddEntry(hSKI,  "SKI",  "l");
    leg->AddEntry(hSKII, "SKII", "l");
    leg->Draw();

    // Ratio plot
    c->cd(2);
    gPad->SetLogx();
    rCR->SetMinimum(0.85);  rCR->SetMaximum(1.1);
    rCR->SetLineColor(kRed);
    // rCR->SetTitle("Ratio to Nominal;z = (1 - cos(#theta))/2;Ratio");
    rCR->Draw("HIST");
    rCR->SetLineWidth(3);
    rSKI->SetLineColor(kBlue);
    rSKI->SetLineWidth(3);
    rSKI->Draw("HIST SAME");
    rSKII->SetLineColor(kGreen+2);
    rSKII->SetLineWidth(3);
    rSKII->Draw("HIST SAME");

    c->SaveAs("b_theta_corr_comparison_full_all.png");
}

