void plot() {
    // Turn off stats boxes
    gStyle->SetOptStat(0);

    // Open files and get trees
    TFile *f0 = TFile::Open("testnoCR.root");
    TFile *f1 = TFile::Open("test.root");
    TTree *t0 = (TTree*)f0->Get("analysis");
    TTree *t1 = (TTree*)f1->Get("analysis");

    // Define log bins
    int nBins = 100;
    double xMin = 0.01, xMax = 3.2;
    double logMin = std::log10(xMin);
    double logMax = std::log10(xMax);
    double binWidth = (logMax - logMin) / nBins;
    double* logBins = new double[nBins + 1];
    for (int i = 0; i <= nBins; ++i) logBins[i] = std::pow(10, logMin + i * binWidth);

    // Create histograms
    TH1F* h0 = new TH1F("h0", "b_theta_corr from testnoCR", nBins, logBins);
    TH1F* h1 = new TH1F("h1", "b_theta_corr from test", nBins, logBins);

    // Fill histograms
    t0->Draw("b_theta_corr>>h0", "", "goff");
    t1->Draw("b_theta_corr>>h1", "", "goff");

    // Normalize both histograms by their integrals
    h0->Scale(1.0 / h0->Integral());
    h1->Scale(1.0 / h1->Integral());

    // Styling
    h0->SetLineColor(kRed);   h0->SetLineWidth(2);
    h1->SetLineColor(kBlue);  h1->SetLineWidth(2);

    // Create canvas with 2 pads
    TCanvas* c = new TCanvas("c", "Comparison with Ratio", 800, 800);
    c->Divide(1, 2);

    // Top pad: overlay histograms
    c->cd(1);
    gPad->SetLogx();
    gPad->SetLogy();
    h0->SetTitle("b_theta_corr;b_theta_corr;Normalized Events");
    h0->Draw("hist");
    h1->Draw("hist same");

    // Legend moved to top-left
    TLegend* leg = new TLegend(0.15, 0.75, 0.45, 0.88);
    leg->AddEntry(h0, "testnoCR.root (normalized)", "l");
    leg->AddEntry(h1, "test.root (normalized)", "l");
    leg->Draw();

    // Bottom pad: ratio h1/h0
    c->cd(2);
    gPad->SetLogx();

    TH1F* h_ratio = (TH1F*)h1->Clone("h_ratio");
    h_ratio->Divide(h0);
    h_ratio->SetTitle("Ratio (test / testnoCR);b_theta_corr;Ratio");
    h_ratio->SetLineColor(kBlack);
    h_ratio->SetMinimum(0.9);
    h_ratio->SetMaximum(1.1);
    h_ratio->Draw("hist");

    c->Update();
    c->SaveAs("b_theta_corr_comparison.png");
}

