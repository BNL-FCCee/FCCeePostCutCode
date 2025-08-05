void plot() {
    // Turn off stats boxes
    gStyle->SetOptStat(0);

    // Open files and get trees
    TFile *fNom = TFile::Open("/atlasgpfs01/usatlas/data/kleibensperger/WWCR_100k_events/test.root");
    TFile *fCR = TFile::Open("/atlasgpfs01/usatlas/data/kleibensperger/WWCR_100k_events/test_noCR.root");
    TFile *fSKI_kl5 = TFile::Open("/atlasgpfs01/usatlas/data/kleibensperger/WWCR_100k_events/test_SKI_kl5.root");
    TFile *fSKI_kl50 = TFile::Open("/atlasgpfs01/usatlas/data/kleibensperger/WWCR_100k_events/test_SKI_kl50.root");
    TFile *fSKI_kl100 = TFile::Open("/atlasgpfs01/usatlas/data/kleibensperger/WWCR_100k_events/test_SKI_kl100.root");

    TTree *tNom = (TTree*)fNom->Get("analysis");
    TTree *tCR = (TTree*)fCR->Get("analysis");
    TTree *tSKI_kl5 = (TTree*)fSKI_kl5->Get("analysis");
    TTree *tSKI_kl50 = (TTree*)fSKI_kl50->Get("analysis");
    TTree *tSKI_kl100 = (TTree*)fSKI_kl100->Get("analysis");

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

    TH1F *h_delR_Nom = new TH1F("h_delR_Nom", "#deltaR- Nominal", nBins, linearBins);
    TH1F *h_delR_CR = new TH1F("h_delR_CR", "", nBins, linearBins);

    TH1F *h_eec_theta_Nom = new TH1F("h_eec_theta_Nom", "Change Title;Normalized Events", nBins, logBins);
    TH1F *h_eec_theta_CR = new TH1F("h_eec_theta_CR",  "", nBins, logBins);

    TH1D* h_cos_theta_Nom = new TH1D("h_cos_theta_Nom", "Cos(theta) - Nominal", nBins, linearBins_2);
    TH1D* h_cos_theta_CR = new TH1D("h_cos_theta_CR", "Cos(theta) - noCR", nBins, linearBins_2);

    TH1D* h_eec_Nom = new TH1D("h_eec_Nom", "Cos(theta) - Nominal", nBins, linearBins_2);
    TH1D* h_eec_CR = new TH1D("h_eec_CR", "Cos(theta) - noCR", nBins, linearBins_2);

    // ************************ Reconstructed W Bosons ************************
    TH1F* h_mass_W1_Nom = new TH1F("h_mass_W1_Nom", "Mass", nBins, linearBins);
    TH1F* h_mass_W1_CR = new TH1F("h_mass_W1_CR", "", nBins, linearBins);

    TH1F* h_mass_W2_Nom = new TH1F("h_mass_W2_Nom", "Mass", nBins, linearBins);
    TH1F* h_mass_W2_CR = new TH1F("h_mass_W2_CR", "", nBins, linearBins);

    TH1F* h_p_W1_Nom = new TH1F("h_p_W1_Nom", "Momentum", nBins, linearBins);
    TH1F* h_p_W1_CR = new TH1F("h_p_W1_CR", "", nBins, linearBins);

    TH1F* h_p_W2_Nom = new TH1F("h_p_W2_Nom", "Momentum", nBins, linearBins);
    TH1F* h_p_W2_CR = new TH1F("h_p_W2_CR", "", nBins, linearBins);

    TH1F* h_e_W1_Nom = new TH1F("h_e_W1_Nom", "Energy", nBins, linearBins);
    TH1F* h_e_W1_CR = new TH1F("h_e_W1_CR", "", nBins, linearBins);

    TH1F* h_e_W2_Nom = new TH1F("h_e_W2_Nom", "Energy", nBins, linearBins);
    TH1F* h_e_W2_CR = new TH1F("h_e_W2_CR", "", nBins, linearBins);

    // ************************ Reconstructed Jets ************************
    TH1F* h_p_reco_Jet1_Nom = new TH1F("h_p_reco_Jet1_Nom", "Momentum", nBins, linearBins);
    TH1F* h_p_reco_Jet1_CR = new TH1F("h_p_reco_Jet1_CR", "", nBins, linearBins);
    TH1F* h_p_reco_Jet2_Nom = new TH1F("h_p_reco_Jet2_Nom", "Momentum", nBins, linearBins);
    TH1F* h_p_reco_Jet2_CR = new TH1F("h_p_reco_Jet2_CR", "", nBins, linearBins);
    TH1F* h_p_reco_Jet3_Nom = new TH1F("h_p_reco_Jet3_Nom", "Momentum", nBins, linearBins);
    TH1F* h_p_reco_Jet3_CR = new TH1F("h_p_reco_Jet3_CR", "", nBins, linearBins);
    TH1F* h_p_reco_Jet4_Nom = new TH1F("h_p_reco_Jet4_Nom", "Momentum", nBins, linearBins);
    TH1F* h_p_reco_Jet4_CR = new TH1F("h_p_reco_Jet4_CR", "", nBins, linearBins);

    TH1F* h_theta_Jet1_Jet2_Nom = new TH1F("h_theta_Jet1_Jet2_Nom", "#theta", nBins, logBins);
    TH1F* h_theta_Jet1_Jet2_CR = new TH1F("h_theta_Jet1_Jet2_CR", "", nBins, logBins);
    TH1F* h_theta_Jet3_Jet4_Nom = new TH1F("h_theta_Jet3_Jet4_Nom", "#theta", nBins, logBins);
    TH1F* h_theta_Jet3_Jet4_CR = new TH1F("h_theta_Jet3_Jet4_CR", "", nBins, logBins);

    TH1F* h_phi_Jet1_Jet2_Nom = new TH1F("h_phi_Jet1_Jet2_Nom", "#theta", nBins, logBins);
    TH1F* h_phi_Jet1_Jet2_CR = new TH1F("h_phi_Jet1_Jet2_CR", "", nBins, logBins);
    TH1F* h_phi_Jet3_Jet4_Nom = new TH1F("h_phi_Jet3_Jet4_Nom", "#theta", nBins, logBins);
    TH1F* h_phi_Jet3_Jet4_CR = new TH1F("h_phi_Jet3_Jet4_CR", "", nBins, logBins);

    TH1F* h_chi_Jet1_Jet2_Nom = new TH1F("h_chi_Jet1_Jet2_Nom", "#theta", nBins, logBins);
    TH1F* h_chi_Jet1_Jet2_CR = new TH1F("h_chi_Jet1_Jet2_CR", "", nBins, logBins);
    TH1F* h_chi_Jet3_Jet4_Nom = new TH1F("h_chi_Jet3_Jet4_Nom", "#theta", nBins, logBins);
    TH1F* h_chi_Jet3_Jet4_CR = new TH1F("h_chi_Jet3_Jet4_CR", "", nBins, logBins);

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

    tNom->Draw("b_delR_corr >> h_delR_Nom", "", "goff");
    tCR->Draw("b_delR_corr >> h_delR_CR", "", "goff");

    // ***************** W Bosons *****************
    tNom->Draw("b_mass_W1 >> h_mass_W1_Nom", "", "goff");
    tCR->Draw("b_mass_W1 >> h_mass_W1_CR", "", "goff");
    tNom->Draw("b_mass_W2 >> h_mass_W2_Nom", "", "goff");
    tCR->Draw("b_mass_W2 >> h_mass_W2_CR", "", "goff");

    tNom->Draw("b_p_W1 >> h_p_W1_Nom", "", "goff");
    tCR->Draw("b_p_W1 >> h_p_W1_CR", "", "goff");
    tNom->Draw("b_p_W2 >> h_p_W2_Nom", "", "goff");
    tCR->Draw("b_p_W2 >> h_p_W2_CR", "", "goff");

    tNom->Draw("b_e_W1 >> h_e_W1_Nom", "", "goff");
    tCR->Draw("b_e_W1 >> h_e_W1_CR", "", "goff");
    tNom->Draw("b_e_W2 >> h_e_W2_Nom", "", "goff");
    tCR->Draw("b_e_W2 >> h_e_W2_CR", "", "goff");

    // ***************** Reconstructed Jets *****************
    tNom->Draw("b_p_reco_Jet1 >> h_p_reco_Jet1_Nom", "", "goff");
    tCR->Draw("b_p_reco_Jet1 >> h_p_reco_Jet1_CR", "", "goff");
    tNom->Draw("b_p_reco_Jet2 >> h_p_reco_Jet2_Nom", "", "goff");
    tCR->Draw("b_p_reco_Jet2 >> h_p_reco_Jet2_CR", "", "goff");
    tNom->Draw("b_p_reco_Jet3 >> h_p_reco_Jet3_Nom", "", "goff");
    tCR->Draw("b_p_reco_Jet3 >> h_p_reco_Jet3_CR", "", "goff");
    tNom->Draw("b_p_reco_Jet4 >> h_p_reco_Jet4_Nom", "", "goff");
    tCR->Draw("b_p_reco_Jet4 >> h_p_reco_Jet4_CR", "", "goff");

    tNom->Draw("b_theta_Jet1_Jet2 >> h_theta_Jet1_Jet2_Nom", "", "goff");
    tCR->Draw("b_theta_Jet1_Jet2 >> h_theta_Jet1_Jet2_CR", "", "goff");
    tNom->Draw("b_theta_Jet3_Jet4 >> h_theta_Jet3_Jet4_Nom", "", "goff");
    tCR->Draw("b_theta_Jet3_Jet4 >> h_theta_Jet3_Jet4_CR", "", "goff");

    tNom->Draw("b_phi_Jet1_Jet2 >> h_phi_Jet1_Jet2_Nom", "", "goff");
    tCR->Draw("b_phi_Jet1_Jet2 >> h_phi_Jet1_Jet2_CR", "", "goff");
    tNom->Draw("b_phi_Jet3_Jet4 >> h_phi_Jet3_Jet4_Nom", "", "goff");
    tCR->Draw("b_phi_Jet3_Jet4 >> h_phi_Jet3_Jet4_CR", "", "goff");

    tNom->Draw("b_chi_Jet1_Jet2 >> h_chi_Jet1_Jet2_Nom", "", "goff");
    tCR->Draw("b_chi_Jet1_Jet2 >> h_chi_Jet1_Jet2_CR", "", "goff");
    tNom->Draw("b_chi_Jet3_Jet4 >> h_chi_Jet3_Jet4_Nom", "", "goff");
    tCR->Draw("b_chi_Jet3_Jet4 >> h_chi_Jet3_Jet4_CR", "", "goff");

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

    h_delR_Nom->Scale(1.0 / h_delR_Nom->Integral());
    h_delR_CR->Scale(1.0 / h_delR_CR->Integral());

    // ***************** W Bosons *****************
    h_mass_W1_Nom->Scale(1.0 / h_mass_W1_Nom->Integral());
    h_mass_W1_CR->Scale(1.0 / h_mass_W1_CR->Integral());
    h_mass_W2_Nom->Scale(1.0 / h_mass_W2_Nom->Integral());
    h_mass_W2_CR->Scale(1.0 / h_mass_W2_CR->Integral());

    h_p_W1_Nom->Scale(1.0 / h_p_W1_Nom->Integral());
    h_p_W1_CR->Scale(1.0 / h_p_W1_CR->Integral());
    h_p_W2_Nom->Scale(1.0 / h_p_W2_Nom->Integral());
    h_p_W2_CR->Scale(1.0 / h_p_W2_CR->Integral());

    h_e_W1_Nom->Scale(1.0 / h_e_W1_Nom->Integral());
    h_e_W1_CR->Scale(1.0 / h_e_W1_CR->Integral());
    h_e_W2_Nom->Scale(1.0 / h_e_W2_Nom->Integral());
    h_e_W2_CR->Scale(1.0 / h_e_W2_CR->Integral());

    // ***************** Reconstructed Jets *****************
    h_p_reco_Jet1_Nom->Scale(1.0 / h_p_reco_Jet1_Nom->Integral());
    h_p_reco_Jet1_CR->Scale(1.0 / h_p_reco_Jet1_CR->Integral());
    h_p_reco_Jet2_Nom->Scale(1.0 / h_p_reco_Jet2_Nom->Integral());
    h_p_reco_Jet2_CR->Scale(1.0 / h_p_reco_Jet2_CR->Integral());
    h_p_reco_Jet3_Nom->Scale(1.0 / h_p_reco_Jet3_Nom->Integral());
    h_p_reco_Jet3_CR->Scale(1.0 / h_p_reco_Jet3_CR->Integral());
    h_p_reco_Jet4_Nom->Scale(1.0 / h_p_reco_Jet4_Nom->Integral());
    h_p_reco_Jet4_CR->Scale(1.0 / h_p_reco_Jet4_CR->Integral());

    h_theta_Jet1_Jet2_Nom->Scale(1.0 / h_theta_Jet1_Jet2_Nom->Integral());
    h_theta_Jet1_Jet2_CR->Scale(1.0 / h_theta_Jet1_Jet2_CR->Integral());
    h_theta_Jet3_Jet4_Nom->Scale(1.0 / h_theta_Jet3_Jet4_Nom->Integral());
    h_theta_Jet3_Jet4_CR->Scale(1.0 / h_theta_Jet3_Jet4_CR->Integral());

    h_phi_Jet1_Jet2_Nom->Scale(1.0 / h_phi_Jet1_Jet2_Nom->Integral());
    h_phi_Jet1_Jet2_CR->Scale(1.0 / h_phi_Jet1_Jet2_CR->Integral());
    h_phi_Jet3_Jet4_Nom->Scale(1.0 / h_phi_Jet3_Jet4_Nom->Integral());
    h_phi_Jet3_Jet4_CR->Scale(1.0 / h_phi_Jet3_Jet4_CR->Integral());

    h_chi_Jet1_Jet2_Nom->Scale(1.0 / h_chi_Jet1_Jet2_Nom->Integral());
    h_chi_Jet1_Jet2_CR->Scale(1.0 / h_chi_Jet1_Jet2_CR->Integral());
    h_chi_Jet3_Jet4_Nom->Scale(1.0 / h_chi_Jet3_Jet4_Nom->Integral());
    h_chi_Jet3_Jet4_CR->Scale(1.0 / h_chi_Jet3_Jet4_CR->Integral());

     // Ratio to Nominal
    TH1F *r_p_jet_const_1_CR = (TH1F*)h_p_jet_const_1_CR->Clone("r_p_jet_const_1_CR");   r_p_jet_const_1_CR->Divide(h_p_jet_const_1_Nom);
    TH1F *r_p_jet_const_2_CR = (TH1F*)h_p_jet_const_2_CR->Clone("r_p_jet_const_2_CR");   r_p_jet_const_2_CR->Divide(h_p_jet_const_2_Nom);
    TH1F *r_theta_CR = (TH1F*)h_theta_CR->Clone("r_theta_CR");   r_theta_CR->Divide(h_theta_Nom);
    TH1F *r_phi_CR = (TH1F*)h_phi_CR->Clone("r_phi_CR");   r_phi_CR->Divide(h_phi_Nom);
    TH1F *r_chi_CR = (TH1F*)h_chi_CR->Clone("r_chi_CR");   r_chi_CR->Divide(h_chi_Nom);
    TH1F *r_delR_CR = (TH1F*)h_delR_CR->Clone("r_delR_CR");   r_delR_CR->Divide(h_delR_Nom);
    TH1D *r_eec_CR = (TH1D*)h_eec_CR->Clone("r_eec_CR");     r_eec_CR->Divide(h_eec_Nom);

    TH1F *r_mass_W1_CR = (TH1F*)h_mass_W1_CR->Clone("r_mass_W1_CR");  r_mass_W1_CR->Divide(h_mass_W1_Nom);
    TH1F *r_mass_W2_CR = (TH1F*)h_mass_W2_CR->Clone("r_mass_W2_CR");  r_mass_W2_CR->Divide(h_mass_W2_Nom);
    TH1F *r_p_W1_CR = (TH1F*)h_p_W1_CR->Clone("r_p_W1_CR");  r_p_W1_CR->Divide(h_p_W1_Nom);
    TH1F *r_p_W2_CR = (TH1F*)h_p_W2_CR->Clone("r_p_W2_CR");  r_p_W2_CR->Divide(h_p_W2_Nom);
    TH1F *r_e_W1_CR = (TH1F*)h_e_W1_CR->Clone("r_e_W1_CR");  r_e_W1_CR->Divide(h_e_W1_Nom);
    TH1F *r_e_W2_CR = (TH1F*)h_e_W2_CR->Clone("r_e_W2_CR");  r_e_W2_CR->Divide(h_e_W2_Nom);

    TH1F *r_p_reco_Jet1_CR = (TH1F*)h_p_reco_Jet1_CR->Clone("r_p_reco_Jet1_CR");  r_p_reco_Jet1_CR->Divide(h_p_reco_Jet1_Nom);
    TH1F *r_p_reco_Jet2_CR = (TH1F*)h_p_reco_Jet2_CR->Clone("r_p_reco_Jet2_CR");  r_p_reco_Jet2_CR->Divide(h_p_reco_Jet2_Nom);
    TH1F *r_p_reco_Jet3_CR = (TH1F*)h_p_reco_Jet3_CR->Clone("r_p_reco_Jet3_CR");  r_p_reco_Jet3_CR->Divide(h_p_reco_Jet3_Nom);
    TH1F *r_p_reco_Jet4_CR = (TH1F*)h_p_reco_Jet4_CR->Clone("r_p_reco_Jet4_CR");  r_p_reco_Jet4_CR->Divide(h_p_reco_Jet4_Nom);
    TH1F *r_theta_Jet1_Jet2_CR = (TH1F*)h_theta_Jet1_Jet2_CR->Clone("r_theta_Jet1_Jet2_CR"); r_theta_Jet1_Jet2_CR->Divide(h_theta_Jet1_Jet2_Nom);
    TH1F *r_theta_Jet3_Jet4_CR = (TH1F*)h_theta_Jet3_Jet4_CR->Clone("r_theta_Jet3_Jet4_CR"); r_theta_Jet3_Jet4_CR->Divide(h_theta_Jet3_Jet4_Nom);
    TH1F *r_phi_Jet1_Jet2_CR = (TH1F*)h_phi_Jet1_Jet2_CR->Clone("r_phi_Jet1_Jet2_CR"); r_phi_Jet1_Jet2_CR->Divide(h_phi_Jet1_Jet2_Nom);
    TH1F *r_phi_Jet3_Jet4_CR = (TH1F*)h_phi_Jet3_Jet4_CR->Clone("r_phi_Jet3_Jet4_CR"); r_phi_Jet3_Jet4_CR->Divide(h_phi_Jet3_Jet4_Nom);
    TH1F *r_chi_Jet1_Jet2_CR = (TH1F*)h_chi_Jet1_Jet2_CR->Clone("r_chi_Jet1_Jet2_CR"); r_chi_Jet1_Jet2_CR->Divide(h_chi_Jet1_Jet2_Nom);
    TH1F *r_chi_Jet3_Jet4_CR = (TH1F*)h_chi_Jet3_Jet4_CR->Clone("r_chi_Jet3_Jet4_CR"); r_chi_Jet3_Jet4_CR->Divide(h_chi_Jet3_Jet4_Nom);
    
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

    auto leg_p_jetConst_1 = new TLegend(0.725, 0.75, 0.88, 0.88);
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

    auto leg_p_jetConst_2 = new TLegend(0.725, 0.75, 0.88, 0.88);
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

    auto leg_theta = new TLegend(0.15, 0.3, 0.25, 0.4);
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

    auto leg_phi = new TLegend(0.15, 0.75, 0.25, 0.88);
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

    auto leg_chi = new TLegend(0.15, 0.75, 0.25, 0.88);
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


    TCanvas *c_delR = new TCanvas("c_delR", "comparison", 1000, 800);
    c_delR->Divide(1, 2);

    // Top pad (60% height): from y = 0.4 to y = 1.0
    TPad *topPad_delR = new TPad("topPad_delR", "Top Pad", 0.0, 0.4, 1.0, 1.0);
    // topPad->SetTopMargin(0.1);    // Optional
    topPad_delR->SetBottomMargin(0.025);
    topPad_delR->Draw();

    // Bottom pad (40% height): from y = 0.0 to y = 0.4
    TPad *bottomPad_delR = new TPad("bottomPad_delR", "Bottom Pad", 0.0, 0.0, 1.0, 0.4);
    bottomPad_delR->SetTopMargin(0.025);    // Small top margin
    bottomPad_delR->SetBottomMargin(0.2); // More space for x-axis label
    bottomPad_delR->Draw();

    // Upper plot
    topPad_delR->cd();
    gPad->SetLogx();
    gPad->SetLogy();

    h_delR_Nom->SetTitle("b_delR_corr");
    h_delR_Nom->GetXaxis()->SetLabelSize(0);
    h_delR_Nom->GetYaxis()->SetTitle("normalized events");
    h_delR_Nom->GetYaxis()->CenterTitle(true);
    h_delR_Nom->GetXaxis()->SetTitleSize(0.08);
    h_delR_Nom->GetYaxis()->SetTitleOffset(1.0);

    h_delR_Nom->SetLineColor(kBlack); h_delR_Nom->SetLineWidth(3);
    h_delR_CR->SetLineColor(kRed); h_delR_CR->SetLineWidth(3);

    h_delR_Nom->Draw("HIST");
    h_delR_CR->Draw("HIST SAME");

    auto leg_delR = new TLegend(0.15, 0.75, 0.25, 0.88);
    leg_delR->AddEntry(h_delR_Nom,  "#deltaR- Nom", "l");
    leg_delR->AddEntry(h_delR_CR,   "#deltaR- no CR",   "l");
    leg_delR->Draw();

    // Ratio plot
    bottomPad_delR->cd();
    gPad->SetLogx();

    r_delR_CR->SetMinimum(0.85);  r_delR_CR->SetMaximum(1.1);
    r_delR_CR->SetLineColor(kRed);
    r_delR_CR->Draw("HIST");
    r_delR_CR->SetLineWidth(3);

    r_delR_CR->GetXaxis()->SetTitle("#deltaR (radians)"); // Axis titles
    r_delR_CR->GetYaxis()->SetTitle("ratio");
    r_delR_CR->GetXaxis()->CenterTitle(true); // Center them
    r_delR_CR->GetYaxis()->CenterTitle(true);
    r_delR_CR->GetXaxis()->SetTitleSize(0.07); // Set the title size
    r_delR_CR->GetYaxis()->SetTitleSize(0.055);
    r_delR_CR->GetXaxis()->SetLabelSize(0.05); // Set the label size
    r_delR_CR->GetYaxis()->SetLabelSize(0.05);
    r_delR_CR->GetXaxis()->SetTitleOffset(1.0); // Set the title offset to the axis
    r_delR_CR->GetYaxis()->SetTitleOffset(0.55);

    // **************************************** Reconstructed W Boson Plots ****************************************
    TCanvas *c_mass_W = new TCanvas("c_mass_W", "comparison", 1000, 800);
    c_mass_W->Divide(1, 2);

    // Top pad (60% height): from y = 0.4 to y = 1.0
    TPad *topPad_mass_W = new TPad("topPad_mass_W", "Top Pad", 0.0, 0.4, 1.0, 1.0);
    // topPad->SetTopMargin(0.1);    // Optional
    topPad_mass_W->SetBottomMargin(0.025);
    topPad_mass_W->Draw();

    // Bottom pad (40% height): from y = 0.0 to y = 0.4
    TPad *bottomPad_mass_W = new TPad("bottomPad_mass_W", "Bottom Pad", 0.0, 0.0, 1.0, 0.4);
    bottomPad_mass_W->SetTopMargin(0.025);    // Small top margin
    bottomPad_mass_W->SetBottomMargin(0.2); // More space for x-axis label
    bottomPad_mass_W->Draw();

    // Upper plot
    topPad_mass_W->cd();
    // gPad->SetLogx();
    // gPad->SetLogy();

    h_mass_W1_Nom->SetTitle("Mass- W Bosons");
    h_mass_W1_Nom->GetXaxis()->SetLabelSize(0);
    h_mass_W1_Nom->GetYaxis()->SetTitle("normalized events");
    h_mass_W1_Nom->GetYaxis()->CenterTitle(true);
    h_mass_W1_Nom->GetXaxis()->SetTitleSize(0.08);
    h_mass_W1_Nom->GetXaxis()->SetLabelSize(0);
    h_mass_W1_Nom->GetYaxis()->SetTitleOffset(1.0);

    h_mass_W1_Nom->SetLineColor(kBlack); h_mass_W1_Nom->SetLineWidth(3);
    h_mass_W1_CR->SetLineColor(kRed);   h_mass_W1_CR->SetLineWidth(3);
    h_mass_W2_Nom->SetLineColor(kBlue); h_mass_W2_Nom->SetLineWidth(3);
    h_mass_W2_CR->SetLineColor(kGreen);   h_mass_W2_CR->SetLineWidth(3);

    h_mass_W1_Nom->Draw("HIST");
    h_mass_W1_CR->Draw("HIST SAME");
    h_mass_W2_Nom->Draw("HIST SAME");
    h_mass_W2_CR->Draw("HIST SAME");

    auto leg_W_mass = new TLegend(0.7, 0.75, 0.88, 0.88);
    leg_W_mass->AddEntry(h_mass_W1_Nom,  "W1 Mass- Nom", "l");
    leg_W_mass->AddEntry(h_mass_W1_CR,   "W1 Mass- No CR",   "l");
    leg_W_mass->AddEntry(h_mass_W2_Nom,  "W2 Mass- Nom", "l");
    leg_W_mass->AddEntry(h_mass_W2_CR,   "W2 Mass- No CR",   "l");
    leg_W_mass->Draw();

    // Ratio plot
    bottomPad_mass_W->cd();
    // gPad->SetLogx();
    // r_p_jet_const_1_CR->SetMinimum(0.85);  r_p_jet_const_1_CR->SetMaximum(1.1);
    r_mass_W1_CR->SetLineColor(kRed);
    r_mass_W1_CR->Draw("HIST");
    r_mass_W1_CR->SetLineWidth(3);
    r_mass_W2_CR->SetLineColor(kBlue);
    r_mass_W2_CR->Draw("HIST SAME");
    r_mass_W2_CR->SetLineWidth(3);

    r_mass_W1_CR->GetXaxis()->SetTitle("Momentum (GeV/c)"); // Axis titles
    r_mass_W1_CR->GetYaxis()->SetTitle("ratio");
    r_mass_W1_CR->SetTitle("");
    r_mass_W1_CR->GetXaxis()->CenterTitle(true); // Center them
    r_mass_W1_CR->GetYaxis()->CenterTitle(true);
    r_mass_W1_CR->GetXaxis()->SetTitleSize(0.055); // Set the title size
    r_mass_W1_CR->GetYaxis()->SetTitleSize(0.055);
    r_mass_W1_CR->GetXaxis()->SetLabelSize(0.05); // Set the label size
    r_mass_W1_CR->GetYaxis()->SetLabelSize(0.05);
    r_mass_W1_CR->GetXaxis()->SetTitleOffset(1.0); // Set the title offset to the axis
    r_mass_W1_CR->GetYaxis()->SetTitleOffset(0.55);

    auto leg_ratio_W_mass = new TLegend(0.15, 0.75, 0.35, 0.9);
    leg_ratio_W_mass->AddEntry(r_mass_W1_CR,  "W1 Mass Nominal", "l");
    leg_ratio_W_mass->AddEntry(r_mass_W2_CR,   "W2 Mass noCR",   "l");
    leg_ratio_W_mass->Draw();

    
    TCanvas *c_p_W = new TCanvas("c_p_W", "comparison", 1000, 800);
    c_p_W->Divide(1, 2);

    // Top pad (60% height): from y = 0.4 to y = 1.0
    TPad *topPad_p_W1 = new TPad("topPad_p_W1", "Top Pad", 0.0, 0.4, 1.0, 1.0);
    // topPad->SetTopMargin(0.1);    // Optional
    topPad_p_W1->SetBottomMargin(0.025);
    topPad_p_W1->Draw();

    // Bottom pad (40% height): from y = 0.0 to y = 0.4
    TPad *bottomPad_p_W1 = new TPad("bottomPad_p_W1", "Bottom Pad", 0.0, 0.0, 1.0, 0.4);
    bottomPad_p_W1->SetTopMargin(0.025);    // Small top margin
    bottomPad_p_W1->SetBottomMargin(0.2); // More space for x-axis label
    bottomPad_p_W1->Draw();

    // Upper plot
    topPad_p_W1->cd();
    // gPad->SetLogx();
    // gPad->SetLogy();

    h_p_W1_Nom->SetTitle("Momenta- W Bosons");
    h_p_W1_Nom->GetXaxis()->SetLabelSize(0);
    h_p_W1_Nom->GetYaxis()->SetTitle("normalized events");
    h_p_W1_Nom->GetYaxis()->CenterTitle(true);
    h_p_W1_Nom->GetXaxis()->SetTitleSize(0.08);
    h_p_W1_Nom->GetXaxis()->SetLabelSize(0);
    h_p_W1_Nom->GetYaxis()->SetTitleOffset(1.0);

    h_p_W1_Nom->SetLineColor(kBlack); h_p_W1_Nom->SetLineWidth(3);
    h_p_W1_CR->SetLineColor(kRed);   h_p_W1_CR->SetLineWidth(3);
    h_p_W2_Nom->SetLineColor(kBlue); h_p_W2_Nom->SetLineWidth(3);
    h_p_W2_CR->SetLineColor(kGreen);   h_p_W2_CR->SetLineWidth(3);

    h_p_W1_Nom->Draw("HIST");
    h_p_W1_CR->Draw("HIST SAME");
    h_p_W2_Nom->Draw("HIST SAME");
    h_p_W2_CR->Draw("HIST SAME");

    auto leg_p_W = new TLegend(0.7, 0.75, 0.88, 0.88);
    leg_p_W->AddEntry(h_p_W1_Nom,  "W1 Nom", "l");
    leg_p_W->AddEntry(h_p_W1_CR,   "W1 No CR",   "l");
    leg_p_W->AddEntry(h_p_W2_Nom,  "W2 Nom", "l");
    leg_p_W->AddEntry(h_p_W2_CR,   "W2 No CR",   "l");
    leg_p_W->Draw();

    // Ratio plot
    bottomPad_p_W1->cd();
    // gPad->SetLogx();
    r_p_W1_CR->SetMinimum(-0.01);  r_p_W1_CR->SetMaximum(3.5);
    r_p_W1_CR->SetLineColor(kRed);
    r_p_W1_CR->Draw("HIST");
    r_p_W1_CR->SetLineWidth(3);
    r_p_W2_CR->SetLineColor(kBlue);
    r_p_W2_CR->Draw("HIST SAME");
    r_p_W2_CR->SetLineWidth(3);

    r_p_W1_CR->GetXaxis()->SetTitle("momenta (GeV/c)"); // Axis titles
    r_p_W1_CR->GetYaxis()->SetTitle("ratio");
    r_p_W1_CR->SetTitle("");
    r_p_W1_CR->GetXaxis()->CenterTitle(true); // Center them
    r_p_W1_CR->GetYaxis()->CenterTitle(true);
    r_p_W1_CR->GetXaxis()->SetTitleSize(0.055); // Set the title size
    r_p_W1_CR->GetYaxis()->SetTitleSize(0.055);
    r_p_W1_CR->GetXaxis()->SetLabelSize(0.05); // Set the label size
    r_p_W1_CR->GetYaxis()->SetLabelSize(0.05);
    r_p_W1_CR->GetXaxis()->SetTitleOffset(1.0); // Set the title offset to the axis
    r_p_W1_CR->GetYaxis()->SetTitleOffset(0.55);

    auto leg_ratio_p_mass = new TLegend(0.75, 0.75, 0.88, 0.9);
    leg_ratio_p_mass->AddEntry(r_p_W1_CR,  "W1 Momenta", "l");
    leg_ratio_p_mass->AddEntry(r_p_W2_CR,   "W2 Momenta",   "l");
    leg_ratio_p_mass->Draw();


    TCanvas *c_e_W = new TCanvas("c_e_W", "comparison", 1000, 800);
    c_e_W->Divide(1, 2);

    // Top pad (60% height): from y = 0.4 to y = 1.0
    TPad *topPad_e_W1 = new TPad("topPad_e_W1", "Top Pad", 0.0, 0.4, 1.0, 1.0);
    // topPad->SetTopMargin(0.1);    // Optional
    topPad_e_W1->SetBottomMargin(0.025);
    topPad_e_W1->Draw();

    // Bottom pad (40% height): from y = 0.0 to y = 0.4
    TPad *bottomPad_e_W1 = new TPad("bottomPad_e_W1", "Bottom Pad", 0.0, 0.0, 1.0, 0.4);
    bottomPad_e_W1->SetTopMargin(0.025);    // Small top margin
    bottomPad_e_W1->SetBottomMargin(0.2); // More space for x-axis label
    bottomPad_e_W1->Draw();

    // Upper plot
    topPad_e_W1->cd();
    // gPad->SetLogx();
    // gPad->SetLogy();

    h_e_W1_Nom->SetTitle("Momenta- W Bosons");
    h_e_W1_Nom->GetXaxis()->SetLabelSize(0);
    h_e_W1_Nom->GetYaxis()->SetTitle("normalized events");
    h_e_W1_Nom->GetYaxis()->CenterTitle(true);
    h_e_W1_Nom->GetXaxis()->SetTitleSize(0.08);
    h_e_W1_Nom->GetXaxis()->SetLabelSize(0);
    h_e_W1_Nom->GetYaxis()->SetTitleOffset(1.0);

    h_e_W1_Nom->SetLineColor(kBlack); h_e_W1_Nom->SetLineWidth(3);
    h_e_W1_CR->SetLineColor(kRed);   h_e_W1_CR->SetLineWidth(3);
    h_e_W2_Nom->SetLineColor(kBlue); h_e_W2_Nom->SetLineWidth(3);
    h_e_W2_CR->SetLineColor(kGreen);   h_e_W2_CR->SetLineWidth(3);

    h_e_W1_Nom->Draw("HIST");
    h_e_W1_CR->Draw("HIST SAME");
    h_e_W2_Nom->Draw("HIST SAME");
    h_e_W2_CR->Draw("HIST SAME");

    auto leg_e_W = new TLegend(0.15, 0.7, 0.35, 0.88);
    leg_e_W->AddEntry(h_e_W1_Nom,  "W1- Nom", "l");
    leg_e_W->AddEntry(h_e_W1_CR,   "W1- CR",   "l");
    leg_e_W->AddEntry(h_e_W2_Nom,  "W2- Nom", "l");
    leg_e_W->AddEntry(h_e_W2_CR,   "W2- CR",   "l");
    leg_e_W->Draw();

    // Ratio plot
    bottomPad_e_W1->cd();
    // gPad->SetLogx();
    // r_p_jet_const_1_CR->SetMinimum(0.85);  r_p_jet_const_1_CR->SetMaximum(1.1);
    r_e_W1_CR->SetLineColor(kRed);
    r_e_W1_CR->Draw("HIST");
    r_e_W1_CR->SetLineWidth(3);
    r_e_W2_CR->SetLineColor(kBlue);
    r_e_W2_CR->Draw("HIST SAME");
    r_e_W2_CR->SetLineWidth(3);

    r_e_W1_CR->GetXaxis()->SetTitle("momenta (GeV/c)"); // Axis titles
    r_e_W1_CR->GetYaxis()->SetTitle("ratio");
    r_e_W1_CR->SetTitle("");
    r_e_W1_CR->GetXaxis()->CenterTitle(true); // Center them
    r_e_W1_CR->GetYaxis()->CenterTitle(true);
    r_e_W1_CR->GetXaxis()->SetTitleSize(0.055); // Set the title size
    r_e_W1_CR->GetYaxis()->SetTitleSize(0.055);
    r_e_W1_CR->GetXaxis()->SetLabelSize(0.05); // Set the label size
    r_e_W1_CR->GetYaxis()->SetLabelSize(0.05);
    r_e_W1_CR->GetXaxis()->SetTitleOffset(1.0); // Set the title offset to the axis
    r_e_W1_CR->GetYaxis()->SetTitleOffset(0.55);

    auto leg_ratio_e_W = new TLegend(0.725, 0.75, 0.88, 0.9);
    leg_ratio_e_W->AddEntry(r_e_W1_CR,  "W1 p Ratio", "l");
    leg_ratio_e_W->AddEntry(r_e_W2_CR,   "W2 p Ratio",   "l");
    leg_ratio_e_W->Draw();


    // **************************************** Reconstructed Jet Plots ****************************************
    TCanvas *c_jet_momentum = new TCanvas("c_jet_momentum", "comparison", 1000, 800);
    c_jet_momentum->Divide(1, 2);

    // Top pad (60% height): from y = 0.4 to y = 1.0
    TPad *topPad_jet_momentum = new TPad("topPad_jet_momentum", "Top Pad", 0.0, 0.4, 1.0, 1.0);
    topPad_jet_momentum->SetBottomMargin(0.025);
    topPad_jet_momentum->Draw();

    // Bottom pad (40% height): from y = 0.0 to y = 0.4
    TPad *bottomPad_jet_momentum = new TPad("bottomPad_jet_momentum", "Bottom Pad", 0.0, 0.0, 1.0, 0.4);
    bottomPad_jet_momentum->SetTopMargin(0.025);    // Small top margin
    bottomPad_jet_momentum->SetBottomMargin(0.2); // More space for x-axis label
    bottomPad_jet_momentum->Draw();

    // Upper plot
    topPad_jet_momentum->cd();
    // gPad->SetLogx();
    // gPad->SetLogy();

    h_p_reco_Jet1_Nom->SetTitle("Momenta- Jets");
    h_p_reco_Jet1_Nom->GetXaxis()->SetLabelSize(0);
    h_p_reco_Jet1_Nom->GetYaxis()->SetTitle("normalized events");
    h_p_reco_Jet1_Nom->GetYaxis()->CenterTitle(true);
    h_p_reco_Jet1_Nom->GetXaxis()->SetTitleSize(0.08);
    h_p_reco_Jet1_Nom->GetXaxis()->SetLabelSize(0);
    h_p_reco_Jet1_Nom->GetYaxis()->SetTitleOffset(1.0);

    h_p_reco_Jet1_Nom->SetLineColor(kBlack); h_p_reco_Jet1_Nom->SetLineWidth(3);
    h_p_reco_Jet1_CR->SetLineColor(kRed); h_p_reco_Jet1_CR->SetLineWidth(3);
    h_p_reco_Jet2_Nom->SetLineColor(kBlue); h_p_reco_Jet2_Nom->SetLineWidth(3);
    h_p_reco_Jet2_CR->SetLineColor(kGreen); h_p_reco_Jet2_CR->SetLineWidth(3);

    h_p_reco_Jet3_Nom->SetLineColor(kOrange+7); h_p_reco_Jet3_Nom->SetLineWidth(3);
    h_p_reco_Jet3_CR->SetLineColor(kPink+7); h_p_reco_Jet3_CR->SetLineWidth(3);
    h_p_reco_Jet4_Nom->SetLineColor(kAzure+7); h_p_reco_Jet4_Nom->SetLineWidth(3);
    h_p_reco_Jet4_CR->SetLineColor(kGreen+3); h_p_reco_Jet4_CR->SetLineWidth(3);

    h_p_reco_Jet1_Nom->Draw("HIST");
    h_p_reco_Jet1_CR->Draw("HIST SAME");
    h_p_reco_Jet2_Nom->Draw("HIST SAME");
    h_p_reco_Jet2_CR->Draw("HIST SAME");
    h_p_reco_Jet3_Nom->Draw("HIST SAME");
    h_p_reco_Jet3_CR->Draw("HIST SAME");
    h_p_reco_Jet4_Nom->Draw("HIST SAME");
    h_p_reco_Jet4_CR->Draw("HIST SAME");

    auto leg_reco_jets = new TLegend(0.725, 0.48, 0.88, 0.88);
    leg_reco_jets->AddEntry(h_p_reco_Jet1_Nom,  "Jet1- Nom", "l");
    leg_reco_jets->AddEntry(h_p_reco_Jet1_CR,   "Jet1- CR",   "l");
    leg_reco_jets->AddEntry(h_p_reco_Jet2_Nom,  "Jet2- Nom", "l");
    leg_reco_jets->AddEntry(h_p_reco_Jet2_CR,   "Jet2- CR",   "l");
    leg_reco_jets->AddEntry(h_p_reco_Jet3_Nom,  "Jet3- Nom", "l");
    leg_reco_jets->AddEntry(h_p_reco_Jet3_CR,   "Jet3- CR",   "l");
    leg_reco_jets->AddEntry(h_p_reco_Jet4_Nom,  "Jet4- Nom", "l");
    leg_reco_jets->AddEntry(h_p_reco_Jet4_CR,   "Jet4- CR",   "l");
    leg_reco_jets->Draw();

    // Ratio plot
    bottomPad_jet_momentum->cd();
    // gPad->SetLogx();
    r_p_reco_Jet1_CR->SetMinimum(0.0);  r_p_jet_const_1_CR->SetMaximum(3.5);
    r_p_reco_Jet1_CR->SetLineColor(kBlack);
    r_p_reco_Jet1_CR->Draw("HIST");
    r_p_reco_Jet1_CR->SetLineWidth(3);
    r_p_reco_Jet2_CR->SetLineColor(kRed);
    r_p_reco_Jet2_CR->Draw("HIST SAME");
    r_p_reco_Jet2_CR->SetLineWidth(3);
    r_p_reco_Jet3_CR->SetLineColor(kBlue);
    r_p_reco_Jet3_CR->Draw("HIST SAME");
    r_p_reco_Jet3_CR->SetLineWidth(3);
    r_p_reco_Jet4_CR->SetLineColor(kGreen);
    r_p_reco_Jet4_CR->Draw("HIST SAME");
    r_p_reco_Jet4_CR->SetLineWidth(3);

    r_p_reco_Jet1_CR->GetXaxis()->SetTitle("Momenta (GeV/c)"); // Axis titles
    r_p_reco_Jet1_CR->GetYaxis()->SetTitle("ratio");
    r_p_reco_Jet1_CR->SetTitle("");
    r_p_reco_Jet1_CR->GetXaxis()->CenterTitle(true); // Center them
    r_p_reco_Jet1_CR->GetYaxis()->CenterTitle(true);
    r_p_reco_Jet1_CR->GetXaxis()->SetTitleSize(0.055); // Set the title size
    r_p_reco_Jet1_CR->GetYaxis()->SetTitleSize(0.055);
    r_p_reco_Jet1_CR->GetXaxis()->SetLabelSize(0.05); // Set the label size
    r_p_reco_Jet1_CR->GetYaxis()->SetLabelSize(0.05);
    r_p_reco_Jet1_CR->GetXaxis()->SetTitleOffset(1.0); // Set the title offset to the axis
    r_p_reco_Jet1_CR->GetYaxis()->SetTitleOffset(0.55);

    auto leg_ratio_reco_jets = new TLegend(0.68, 0.68, 0.88, 0.88);
    leg_ratio_reco_jets->AddEntry(r_p_reco_Jet1_CR, "Jet1- Momenta Ratio", "l");
    leg_ratio_reco_jets->AddEntry(r_p_reco_Jet2_CR, "Jet2- Momenta Ratio",   "l");
    leg_ratio_reco_jets->AddEntry(r_p_reco_Jet3_CR, "Jet3- Momenta Ratio", "l");
    leg_ratio_reco_jets->AddEntry(r_p_reco_Jet4_CR, "Jet4- Momenta Ratio",   "l");
    leg_ratio_reco_jets->Draw();


    TCanvas *c_theta_reco_jet = new TCanvas("c_theta_reco_jet", "comparison", 1000, 800);
    c_theta_reco_jet->Divide(1, 2);

    // Top pad (60% height): from y = 0.4 to y = 1.0
    TPad *topPad_theta_reco_jet = new TPad("topPad_theta_reco_jet", "Top Pad", 0.0, 0.4, 1.0, 1.0);
    topPad_theta_reco_jet->SetBottomMargin(0.025);
    topPad_theta_reco_jet->Draw();

    // Bottom pad (40% height): from y = 0.0 to y = 0.4
    TPad *bottomPad_theta_reco_jet = new TPad("bottomPad_theta_reco_jet", "Bottom Pad", 0.0, 0.0, 1.0, 0.4);
    bottomPad_theta_reco_jet->SetTopMargin(0.025);    // Small top margin
    bottomPad_theta_reco_jet->SetBottomMargin(0.2); // More space for x-axis label
    bottomPad_theta_reco_jet->Draw();

    // Upper plot
    topPad_theta_reco_jet->cd();
    gPad->SetLogx();
    gPad->SetLogy();

    h_theta_Jet1_Jet2_Nom->SetTitle("#theta- Reconstructed Jets");
    h_theta_Jet1_Jet2_Nom->GetXaxis()->SetLabelSize(0);
    h_theta_Jet1_Jet2_Nom->GetYaxis()->SetTitle("normalized events");
    h_theta_Jet1_Jet2_Nom->GetYaxis()->CenterTitle(true);
    h_theta_Jet1_Jet2_Nom->GetXaxis()->SetTitleSize(0.08);
    h_theta_Jet1_Jet2_Nom->GetXaxis()->SetLabelSize(0);
    h_theta_Jet1_Jet2_Nom->GetYaxis()->SetTitleOffset(1.0);

    h_theta_Jet1_Jet2_Nom->SetLineColor(kBlack); h_theta_Jet1_Jet2_Nom->SetLineWidth(3);
    h_theta_Jet1_Jet2_CR->SetLineColor(kRed); h_theta_Jet1_Jet2_CR->SetLineWidth(3);
    h_theta_Jet3_Jet4_Nom->SetLineColor(kBlue); h_theta_Jet3_Jet4_Nom->SetLineWidth(3);
    h_theta_Jet3_Jet4_CR->SetLineColor(kGreen); h_theta_Jet3_Jet4_CR->SetLineWidth(3);

    h_theta_Jet1_Jet2_Nom->Draw("HIST");
    h_theta_Jet1_Jet2_CR->Draw("HIST SAME");
    h_theta_Jet3_Jet4_Nom->Draw("HIST SAME");
    h_theta_Jet3_Jet4_CR->Draw("HIST SAME");

    auto leg_theta_reco_jets = new TLegend(0.15, 0.73, 0.35, 0.88);
    leg_theta_reco_jets->AddEntry(h_theta_Jet1_Jet2_Nom,  "W1 Jets- Nom", "l");
    leg_theta_reco_jets->AddEntry(h_theta_Jet1_Jet2_CR,   "W1 Jets- CR",   "l");
    leg_theta_reco_jets->AddEntry(h_theta_Jet3_Jet4_Nom,  "W2 Jets- Nom", "l");
    leg_theta_reco_jets->AddEntry(h_theta_Jet3_Jet4_CR,   "W2 Jets- CR",   "l");
    leg_theta_reco_jets->Draw();

    // Ratio plot
    bottomPad_theta_reco_jet->cd();
    gPad->SetLogx();
    // r_p_jet_const_1_CR->SetMinimum(0.85);  r_p_jet_const_1_CR->SetMaximum(1.1);
    r_theta_Jet1_Jet2_CR->SetLineColor(kBlack);
    r_theta_Jet1_Jet2_CR->Draw("HIST");
    r_theta_Jet1_Jet2_CR->SetLineWidth(3);
    r_theta_Jet3_Jet4_CR->SetLineColor(kRed);
    r_theta_Jet3_Jet4_CR->Draw("HIST SAME");
    r_theta_Jet3_Jet4_CR->SetLineWidth(3);

    r_theta_Jet1_Jet2_CR->GetXaxis()->SetTitle("#theta (radians)"); // Axis titles
    r_theta_Jet1_Jet2_CR->GetYaxis()->SetTitle("ratio");
    r_theta_Jet1_Jet2_CR->SetTitle("");
    r_theta_Jet1_Jet2_CR->GetXaxis()->CenterTitle(true); // Center them
    r_theta_Jet1_Jet2_CR->GetYaxis()->CenterTitle(true);
    r_theta_Jet1_Jet2_CR->GetXaxis()->SetTitleSize(0.055); // Set the title size
    r_theta_Jet1_Jet2_CR->GetYaxis()->SetTitleSize(0.055);
    r_theta_Jet1_Jet2_CR->GetXaxis()->SetLabelSize(0.05); // Set the label size
    r_theta_Jet1_Jet2_CR->GetYaxis()->SetLabelSize(0.05);
    r_theta_Jet1_Jet2_CR->GetXaxis()->SetTitleOffset(1.0); // Set the title offset to the axis
    r_theta_Jet1_Jet2_CR->GetYaxis()->SetTitleOffset(0.55);

    auto leg_ratio_theta_reco_jets = new TLegend(0.15, 0.73, 0.35, 0.88);
    leg_ratio_theta_reco_jets->AddEntry(r_theta_Jet1_Jet2_CR,  "#theta W1 Jets", "l");
    leg_ratio_theta_reco_jets->AddEntry(r_theta_Jet3_Jet4_CR,   "#theta W2 Jets",   "l");
    leg_ratio_theta_reco_jets->Draw();


    TCanvas *c_phi_reco_jet = new TCanvas("c_phi_reco_jet", "comparison", 1000, 800);
    c_phi_reco_jet->Divide(1, 2);

    // Top pad (60% height): from y = 0.4 to y = 1.0
    TPad *topPad_phi_reco_jet = new TPad("topPad_phi_reco_jet", "Top Pad", 0.0, 0.4, 1.0, 1.0);
    topPad_phi_reco_jet->SetBottomMargin(0.025);
    topPad_phi_reco_jet->Draw();

    // Bottom pad (40% height): from y = 0.0 to y = 0.4
    TPad *bottomPad_phi_reco_jet = new TPad("bottomPad_phi_reco_jet", "Bottom Pad", 0.0, 0.0, 1.0, 0.4);
    bottomPad_phi_reco_jet->SetTopMargin(0.025);    // Small top margin
    bottomPad_phi_reco_jet->SetBottomMargin(0.2); // More space for x-axis label
    bottomPad_phi_reco_jet->Draw();

    // Upper plot
    topPad_phi_reco_jet->cd();
    gPad->SetLogx();
    gPad->SetLogy();

    h_phi_Jet1_Jet2_Nom->SetTitle("#phi- Reconstructed Jets");
    h_phi_Jet1_Jet2_Nom->GetXaxis()->SetLabelSize(0);
    h_phi_Jet1_Jet2_Nom->GetYaxis()->SetTitle("normalized events");
    h_phi_Jet1_Jet2_Nom->GetYaxis()->CenterTitle(true);
    h_phi_Jet1_Jet2_Nom->GetXaxis()->SetTitleSize(0.08);
    h_phi_Jet1_Jet2_Nom->GetXaxis()->SetLabelSize(0);
    h_phi_Jet1_Jet2_Nom->GetYaxis()->SetTitleOffset(1.0);

    h_phi_Jet1_Jet2_Nom->SetLineColor(kBlack); h_phi_Jet1_Jet2_Nom->SetLineWidth(3);
    h_phi_Jet1_Jet2_CR->SetLineColor(kRed); h_phi_Jet1_Jet2_CR->SetLineWidth(3);
    h_phi_Jet3_Jet4_Nom->SetLineColor(kBlue); h_phi_Jet3_Jet4_Nom->SetLineWidth(3);
    h_phi_Jet3_Jet4_CR->SetLineColor(kGreen); h_phi_Jet3_Jet4_CR->SetLineWidth(3);

    h_phi_Jet1_Jet2_Nom->Draw("HIST");
    h_phi_Jet1_Jet2_CR->Draw("HIST SAME");
    h_phi_Jet3_Jet4_Nom->Draw("HIST SAME");
    h_phi_Jet3_Jet4_CR->Draw("HIST SAME");

    auto leg_phi_reco_jets = new TLegend(0.15, 0.73, 0.35, 0.88);
    leg_phi_reco_jets->AddEntry(h_phi_Jet1_Jet2_Nom,  "W1 Jets- Nom", "l");
    leg_phi_reco_jets->AddEntry(h_phi_Jet1_Jet2_CR,   "W1 Jets- CR",   "l");
    leg_phi_reco_jets->AddEntry(h_phi_Jet3_Jet4_Nom,  "W2 Jets- Nom", "l");
    leg_phi_reco_jets->AddEntry(h_phi_Jet3_Jet4_CR,   "W2 Jets- CR",   "l");
    leg_phi_reco_jets->Draw();

    // Ratio plot
    bottomPad_phi_reco_jet->cd();
    gPad->SetLogx();
    // r_p_jet_const_1_CR->SetMinimum(0.85);  r_p_jet_const_1_CR->SetMaximum(1.1);
    r_phi_Jet1_Jet2_CR->SetLineColor(kBlack);
    r_phi_Jet1_Jet2_CR->Draw("HIST");
    r_phi_Jet1_Jet2_CR->SetLineWidth(3);
    r_phi_Jet3_Jet4_CR->SetLineColor(kRed);
    r_phi_Jet3_Jet4_CR->Draw("HIST SAME");
    r_phi_Jet3_Jet4_CR->SetLineWidth(3);

    r_phi_Jet1_Jet2_CR->GetXaxis()->SetTitle("#phi (radians)"); // Axis titles
    r_phi_Jet1_Jet2_CR->GetYaxis()->SetTitle("ratio");
    r_phi_Jet1_Jet2_CR->SetTitle("");
    r_phi_Jet1_Jet2_CR->GetXaxis()->CenterTitle(true); // Center them
    r_phi_Jet1_Jet2_CR->GetYaxis()->CenterTitle(true);
    r_phi_Jet1_Jet2_CR->GetXaxis()->SetTitleSize(0.055); // Set the title size
    r_phi_Jet1_Jet2_CR->GetYaxis()->SetTitleSize(0.055);
    r_phi_Jet1_Jet2_CR->GetXaxis()->SetLabelSize(0.05); // Set the label size
    r_phi_Jet1_Jet2_CR->GetYaxis()->SetLabelSize(0.05);
    r_phi_Jet1_Jet2_CR->GetXaxis()->SetTitleOffset(1.0); // Set the title offset to the axis
    r_phi_Jet1_Jet2_CR->GetYaxis()->SetTitleOffset(0.55);

    auto leg_ratio_phi_reco_jets = new TLegend(0.725, 0.73, 0.88, 0.88);
    leg_ratio_phi_reco_jets->AddEntry(r_phi_Jet1_Jet2_CR,  "#phi W1 Jets", "l");
    leg_ratio_phi_reco_jets->AddEntry(r_phi_Jet3_Jet4_CR,   "#phi W2 Jets",   "l");
    leg_ratio_phi_reco_jets->Draw();


    TCanvas *c_chi_reco_jet = new TCanvas("c_chi_reco_jet", "comparison", 1000, 800);
    c_chi_reco_jet->Divide(1, 2);

    // Top pad (60% height): from y = 0.4 to y = 1.0
    TPad *topPad_chi_reco_jet = new TPad("topPad_chi_reco_jet", "Top Pad", 0.0, 0.4, 1.0, 1.0);
    topPad_chi_reco_jet->SetBottomMargin(0.025);
    topPad_chi_reco_jet->Draw();

    // Bottom pad (40% height): from y = 0.0 to y = 0.4
    TPad *bottomPad_chi_reco_jet = new TPad("bottomPad_chi_reco_jet", "Bottom Pad", 0.0, 0.0, 1.0, 0.4);
    bottomPad_chi_reco_jet->SetTopMargin(0.025);    // Small top margin
    bottomPad_chi_reco_jet->SetBottomMargin(0.2); // More space for x-axis label
    bottomPad_chi_reco_jet->Draw();

    // Upper plot
    topPad_chi_reco_jet->cd();
    gPad->SetLogx();
    gPad->SetLogy();

    h_chi_Jet1_Jet2_Nom->SetTitle("#chi- Reconstructed Jets");
    h_chi_Jet1_Jet2_Nom->GetXaxis()->SetLabelSize(0);
    h_chi_Jet1_Jet2_Nom->GetYaxis()->SetTitle("normalized events");
    h_chi_Jet1_Jet2_Nom->GetYaxis()->CenterTitle(true);
    h_chi_Jet1_Jet2_Nom->GetXaxis()->SetTitleSize(0.08);
    h_chi_Jet1_Jet2_Nom->GetXaxis()->SetLabelSize(0);
    h_chi_Jet1_Jet2_Nom->GetYaxis()->SetTitleOffset(1.0);

    h_chi_Jet1_Jet2_Nom->SetLineColor(kBlack); h_chi_Jet1_Jet2_Nom->SetLineWidth(3);
    h_chi_Jet1_Jet2_CR->SetLineColor(kRed); h_chi_Jet1_Jet2_CR->SetLineWidth(3);
    h_chi_Jet3_Jet4_Nom->SetLineColor(kBlue); h_chi_Jet3_Jet4_Nom->SetLineWidth(3);
    h_chi_Jet3_Jet4_CR->SetLineColor(kGreen); h_chi_Jet3_Jet4_CR->SetLineWidth(3);

    h_chi_Jet1_Jet2_Nom->Draw("HIST");
    h_chi_Jet1_Jet2_CR->Draw("HIST SAME");
    h_chi_Jet3_Jet4_Nom->Draw("HIST SAME");
    h_chi_Jet3_Jet4_CR->Draw("HIST SAME");

    auto leg_chi_reco_jets = new TLegend(0.15, 0.7, 0.35, 0.88);
    leg_chi_reco_jets->AddEntry(h_chi_Jet1_Jet2_Nom,  "#chi W1 Jets- Nom", "l");
    leg_chi_reco_jets->AddEntry(h_chi_Jet1_Jet2_CR,   "#chi W1 Jets- no CR", "l");
    leg_chi_reco_jets->AddEntry(h_chi_Jet3_Jet4_Nom,  "#chi W2 Jets- Nom", "l");
    leg_chi_reco_jets->AddEntry(h_chi_Jet3_Jet4_CR,   "#chi W1 Jets- no CR", "l");
    leg_chi_reco_jets->Draw();

    // Ratio plot
    bottomPad_chi_reco_jet->cd();
    gPad->SetLogx();
    // r_p_jet_const_1_CR->SetMinimum(0.85);  r_p_jet_const_1_CR->SetMaximum(1.1);
    r_chi_Jet1_Jet2_CR->SetLineColor(kBlack);
    r_chi_Jet1_Jet2_CR->Draw("HIST");
    r_chi_Jet1_Jet2_CR->SetLineWidth(3);
    r_chi_Jet3_Jet4_CR->SetLineColor(kRed);
    r_chi_Jet3_Jet4_CR->Draw("HIST SAME");
    r_chi_Jet3_Jet4_CR->SetLineWidth(3);

    r_chi_Jet1_Jet2_CR->GetXaxis()->SetTitle("#chi (radians)"); // Axis titles
    r_chi_Jet1_Jet2_CR->GetYaxis()->SetTitle("ratio");
    r_chi_Jet1_Jet2_CR->SetTitle("");
    r_chi_Jet1_Jet2_CR->GetXaxis()->CenterTitle(true); // Center them
    r_chi_Jet1_Jet2_CR->GetYaxis()->CenterTitle(true);
    r_chi_Jet1_Jet2_CR->GetXaxis()->SetTitleSize(0.055); // Set the title size
    r_chi_Jet1_Jet2_CR->GetYaxis()->SetTitleSize(0.055);
    r_chi_Jet1_Jet2_CR->GetXaxis()->SetLabelSize(0.05); // Set the label size
    r_chi_Jet1_Jet2_CR->GetYaxis()->SetLabelSize(0.05);
    r_chi_Jet1_Jet2_CR->GetXaxis()->SetTitleOffset(1.0); // Set the title offset to the axis
    r_chi_Jet1_Jet2_CR->GetYaxis()->SetTitleOffset(0.55);

    auto leg_ratio_chi_reco_jets = new TLegend(0.15, 0.75, 0.3, 0.9);
    leg_ratio_chi_reco_jets->AddEntry(r_chi_Jet1_Jet2_CR, "ratio- #chi W1 Jets", "l");
    leg_ratio_chi_reco_jets->AddEntry(r_chi_Jet3_Jet4_CR, "ratio- #chi W2 Jets", "l");
    leg_ratio_chi_reco_jets->Draw();


    c_p_jet_const_1->SaveAs("Jet_const_Momenta_1.png");
    c_p_jet_const_2->SaveAs("Jet_const_Momenta_2.png");
    c_theta->SaveAs("Jet_const_theta.png");
    c_phi->SaveAs("Jet_const_phi.png");
    c_chi->SaveAs("Jet_const_chi.png");
    c_delR->SaveAs("Jet_const_delR.png");

    c_mass_W->SaveAs("W_boson_Mass.png");
    c_p_W->SaveAs("W_boson_Momenta.png");
    c_e_W->SaveAs("W_boson_Energy.png");

    c_jet_momentum->SaveAs("Reco_Jet_Momenta.png");
    c_theta_reco_jet->SaveAs("Reco_Jet_theta.png");
    c_phi_reco_jet->SaveAs("Reco_Jet_phi.png");
    c_chi_reco_jet->SaveAs("Reco_Jet_chi.png");
}