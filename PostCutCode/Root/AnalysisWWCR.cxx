// Local includes
#include <cstdlib>  // for std::exit
#include "PostCutCode/AnalysisWWCR.h"
#include "PostCutCode/TreeContainer.h"
#include <iomanip>
#include <TFile.h>
#include <fstream>
#include <TRandom.h>
#include "TLorentzVector.h"
// #include "Root"
#include <nlohmann/json.hpp>

#include "TH1F.h"
#include <set>

using namespace std;

// Base class 
AnalysisWWCR::AnalysisWWCR():
    AnalysisBase(),
    m_debug(true) 
{}
AnalysisWWCR::~AnalysisWWCR()
{}

struct QuarkCounter {
        int n_c = 0;
        int n_s = 0;
        int n_u = 0;
        int n_d = 0;
        int n_b = 0;
    };

QuarkCounter countQuarks(const ROOT::VecOps::RVec<int>* jet_truth) {
    QuarkCounter counter;

    for (size_t i = 0; i < jet_truth->size(); ++i) {
        int pdg_code = std::abs(jet_truth->at(i));

        if (pdg_code == 5) counter.n_b++;
        if (pdg_code == 4) counter.n_c++;
        if (pdg_code == 3) counter.n_s++;
        if (pdg_code == 2) counter.n_u++;
        if (pdg_code == 1) counter.n_d++;
    }
    return counter;
}

// void CosPhi_Angle(const TLorentzVector& jet1, const TLorentzVector& jet2, TH1F* hist_global_theta, TH1F *hist_delta_theta, TH1F *hist_cos_theta, 
//                   TH1F* hist_delta_eta, TH1F* hist_delta_phi, TH1F* hist_cosphi, TH1F* hist_eec_theta, TH1F* hist_eec_phi) {

//     double global_theta = jet1.Angle(jet2.Vect());
//     double delta_theta = jet1.Theta() - jet2.Theta();
//     double cos_theta = cos(global_theta);
    
//     double delta_eta = jet1.Eta() - jet2.Eta();
//     double delta_phi = jet1.Phi() - jet2.Phi();
//     double cos_phi = cos(delta_phi);

//     double ee_correlator_theta = 0.5 * (1 - cos_theta);
//     double ee_correlator_phi = 0.5 * (1 - cos_phi);

//     // double degrees = radians * 180.0 / TMath::Pi();

//     double global_theta_deg = global_theta * 180.0 / TMath::Pi();
//     double delta_theta_deg = delta_theta * 180.0 / TMath::Pi();
//     double delta_eta_deg = delta_eta * 180 / TMath::Pi();
//     double delta_phi_deg = delta_phi * 180.0 / TMath::Pi();

//     hist_global_theta->Fill(global_theta_deg);
//     hist_delta_theta->Fill(delta_theta_deg);
//     hist_cos_theta->Fill(cos_theta);
//     hist_delta_eta->Fill(delta_eta_deg);
//     hist_delta_phi->Fill(delta_phi_deg);
//     hist_cosphi->Fill(cos_phi);
//     hist_eec_theta->Fill(ee_correlator_theta);
//     hist_eec_phi->Fill(ee_correlator_phi);
// }

void AnalysisWWCR::run() { 

    std::vector<std::string> cutFlowMap {"All Events", "DecayCuts", "NoNaNFlavScore", "Has 4 Jets", "leptonCut", "d123Cut", "d34Cut", "Events have a c, s, and 2 l tagged jets"};
     // using json = nlohmann::json;
    std::ifstream f(MDC::GetInstance()->getSOWJSONFile());
    nlohmann::json data = nlohmann::json::parse(f);

    // override the sum of weights, if it is inside the extra files that we built by hand
    std::ifstream customF(MDC::GetInstance()->getCustomSOWJSONFile());
    nlohmann::json customData = nlohmann::json::parse(customF);

     auto sName = MDC::GetInstance()->getSampleName();

    double norm_weight = (double)data[sName]["crossSection"]/(double)data[sName]["sumOfWeights"];
    
    // Get the histograms
    auto countingHist = m_histContainer->getCountingHist();
    auto cutFlowHist = m_histContainer->get1DHist("cutFlowHist", cutFlowMap.size(), 0, 8, cutFlowMap);

    // ******************* SUBSTRCTURE HISTOGRAMS *******************
    auto h_subJet_1_p = m_histContainer->get1DHist("h_subJet_1_p", 500, 0, 150);
    auto h_subJet_2_p = m_histContainer->get1DHist("h_subJet_2_p", 500, 0, 150);

    auto h_chi_subjet = m_histContainer->get1DHist("h_chi_subjet", 500, 0, TMath::Pi());
    auto h_cos_chi_subjet = m_histContainer->get1DHist("h_cos_chi_subjet", 500, -1, 1);
    auto h_eec_chi_subjet = m_histContainer->get1DHist("h_eec_chi_subjet", 500, -1, 1);

    auto h_theta_subjet = m_histContainer->get1DHist("h_theta_subjet", 500, 0, TMath::Pi());
    auto h_cos_theta_subjet = m_histContainer->get1DHist("h_cos_theta_subjet", 500, -1, 1);
    auto h_eec_theta_subjet = m_histContainer->get1DHist("h_eec_theta_subjet", 500, -1, 1);

    auto h_eta_subjet = m_histContainer->get1DHist("h_eta_subjet", 500, 0, TMath::Pi());

    auto h_phi_subjet = m_histContainer->get1DHist("h_phi_subjet", 500, -TMath::Pi(), TMath::Pi());
    auto h_cos_phi_subjet = m_histContainer->get1DHist("h_cos_phi_subjet", 500, -1, 1);
    auto h_eec_phi_subjet = m_histContainer->get1DHist("h_eec_phi_subjet", 500, -1, 1);
    
    // ******************************** TRUTH HISTOGRAMS ********************************
    auto h_W1_truth_mass = m_histContainer->get1DHist("h_W1_truth_mass", 500, 0, 150);
    auto h_W2_truth_mass = m_histContainer->get1DHist("h_W2_truth_mass", 500, 0, 150);

    auto h_W1_truth_p = m_histContainer->get1DHist("h_W1_truth_p", 500, 0, 150);
    auto h_W2_truth_p = m_histContainer->get1DHist("h_W2_truth_p", 500, 0, 150);

    auto h_W1_truth_e = m_histContainer->get1DHist("h_W1_truth_e", 500, 0, 150);
    auto h_W2_truth_e = m_histContainer->get1DHist("h_W2_truth_e", 500, 0, 150);

    auto h_chi_truth_cs = m_histContainer->get1DHist("h_chi_truth_cs", 500, 0, TMath::Pi());
    auto h_chi_truth_ud = m_histContainer->get1DHist("h_chi_truth_ud", 500, 0, TMath::Pi());
    auto h_cos_chi_truth_cs = m_histContainer->get1DHist("h_cos_chi_truth_cs", 500, -1, 1);
    auto h_cos_chi_truth_ud = m_histContainer->get1DHist("h_cos_chi_truth_ud", 500, -1, 1);
    auto h_eec_chi_truth_cs = m_histContainer->get1DHist("h_eec_chi_truth_cs", 500, -1, 1);
    auto h_eec_chi_truth_ud = m_histContainer->get1DHist("h_eec_chi_truth_ud", 500, -1, 1);

    auto h_theta_truth_cs = m_histContainer->get1DHist("h_theta_truth_cs", 500, 0, TMath::Pi());
    auto h_theta_truth_ud = m_histContainer->get1DHist("h_theta_truth_ud", 500, 0, TMath::Pi());
    auto h_cos_theta_truth_cs = m_histContainer->get1DHist("h_cos_theta_truth_cs", 500, -1, 1);
    auto h_cos_theta_truth_ud = m_histContainer->get1DHist("h_cos_theta_truth_ud", 500, -1, 1);
    auto h_eec_theta_truth_cs = m_histContainer->get1DHist("h_eec_theta_truth_cs", 500, -1, 1);
    auto h_eec_theta_truth_ud = m_histContainer->get1DHist("h_eec_theta_truth_ud", 500, -1, 1);

    auto h_eta_truth_cs = m_histContainer->get1DHist("h_eta_truth_cs", 500, 0, TMath::Pi());
    auto h_eta_truth_ud = m_histContainer->get1DHist("h_eta_truth_ud", 500, 0, TMath::Pi());

    auto h_phi_truth_cs = m_histContainer->get1DHist("h_phi_truth_cs", 500, -TMath::Pi(), TMath::Pi());
    auto h_phi_truth_ud = m_histContainer->get1DHist("h_phi_truth_ud", 500, -TMath::Pi(), TMath::Pi());
    auto h_cos_phi_truth_cs = m_histContainer->get1DHist("h_cos_phi_truth_cs", 500, -1, 1);
    auto h_cos_phi_truth_ud = m_histContainer->get1DHist("h_cos_phi_truth_ud", 500, -1, 1);
    auto h_eec_phi_truth_cs = m_histContainer->get1DHist("h_eec_phi_truth_cs", 500, -1, 1);
    auto h_eec_phi_truth_ud = m_histContainer->get1DHist("h_eec_phi_truth_ud", 500, -1, 1);

    // ******************* RECOJET HISTOGRAMS *******************
    auto h_W1_mass = m_histContainer->get1DHist("h_W1_mass", 500, 0, 150);
    auto h_W2_mass = m_histContainer->get1DHist("h_W2_mass", 500, 0, 150);
    auto h_chi2 = m_histContainer->get1DHist("h_chi2", 500, 0, 150);

    auto h_cJet_p = m_histContainer->get1DHist("cJet_p", 500, 0, 150);
    auto h_lJet0_p = m_histContainer->get1DHist("lJet0_p", 500, 0, 150);
    auto h_lJet1_p = m_histContainer->get1DHist("lJet1_p", 500, 0, 150);
    auto h_lJet2_p = m_histContainer->get1DHist("lJet2_p", 500, 0, 150);

    auto h_W1_p = m_histContainer->get1DHist("h_W1_p", 500, 0, 150);
    auto h_W2_p = m_histContainer->get1DHist("h_W2_p", 500, 0, 150);

    auto h_W1_e = m_histContainer->get1DHist("h_W1_e", 500, 0, 150);
    auto h_W2_e = m_histContainer->get1DHist("h_W2_e", 500, 0, 150);

    auto h_chi_c_l0 = m_histContainer->get1DHist("h_chi_c_l0", 500, 0, TMath::Pi());
    auto h_chi_l1_l2 = m_histContainer->get1DHist("h_chi_l1_l2", 500, 0, TMath::Pi());
    auto h_cos_chi_c_l0 = m_histContainer->get1DHist("h_cos_chi_c_l0", 500, -1, 1);
    auto h_cos_chi_l1_l2 = m_histContainer->get1DHist("h_cos_chi_l1_l2", 500, -1, 1);
    auto h_eec_chi_c_l0 = m_histContainer->get1DHist("h_eec_chi_c_l0", 500, -1, 1);
    auto h_eec_chi_l1_l2 = m_histContainer->get1DHist("h_eec_chi_l1_l2", 500, -1, 1);

    auto h_theta_c_l0 = m_histContainer->get1DHist("h_theta_c_l0", 500, 0, TMath::Pi());
    auto h_theta_l1_l2 = m_histContainer->get1DHist("h_theta_l1_l2", 500, 0, TMath::Pi());
    auto h_cos_theta_c_l0 = m_histContainer->get1DHist("h_cos_theta_c_l0", 500, -1, 1);
    auto h_cos_theta_l1_l2 = m_histContainer->get1DHist("h_cos_theta_l1_l2", 500, -1, 1);
    auto h_eec_theta_c_l0 = m_histContainer->get1DHist("h_eec_theta_c_l0", 500, -1, 1);
    auto h_eec_theta_l1_l2 = m_histContainer->get1DHist("h_eec_theta_l1_l2", 500, -1, 1);

    auto h_eta_c_l0 = m_histContainer->get1DHist("h_eta_c_l0", 500, -TMath::Pi(), TMath::Pi());
    auto h_eta_l1_l2 = m_histContainer->get1DHist("h_eta_l1_l2", 500, -TMath::Pi(), TMath::Pi());

    auto h_phi_c_l0 = m_histContainer->get1DHist("h_phi_c_l0", 500, 0, TMath::Pi());
    auto h_phi_l1_l2 = m_histContainer->get1DHist("h_phi_l1_l2", 500, 0, TMath::Pi());
    auto h_cos_phi_c_l0 = m_histContainer->get1DHist("h_cos_phi_c_l0", 500, -1, 1);
    auto h_cos_phi_l1_l2 = m_histContainer->get1DHist("h_cos_phi_l1_l2", 500, -1, 1);
    auto h_eec_phi_c_l0 = m_histContainer->get1DHist("h_eec_phi_c_l0", 500, -1, 1);
    auto h_eec_phi_l1_l2 = m_histContainer->get1DHist("h_eec_phi_l1_l2", 500, -1, 1);
    
    //Make new tree
    TTree* t = new TTree("analysis","my analysis tree");
    my_tree = (TTree*) t;
    //define the output branches 
     
    my_tree->Branch("b_ee_corr",&ee_corr);
    my_tree->Branch("b_theta_corr",&theta_corr);
    my_tree->Branch("b_phi_corr",&phi_corr);
    my_tree->Branch("b_chi_corr",&chi_corr);
    my_tree->Branch("b_w",&mc_weight);
    // Get the trees
    auto treeCont = std::make_shared<TreeContainer>();
 
    // Get max events to run on
    int nEntries = treeCont->getEventsToRun();

    // Connect branches to trees
    auto tree = treeCont->getTree();

    // calls the branch_name (var_type branch {tree, "branch_name"};)
    varMember<int> event_njet {tree, "event_njet"};
    varMember<ulong> event_nmu {tree, "event_nmu"};
    varMember<ulong> event_nel {tree, "event_nel"};
    varMember<float> muons_p {tree, "muons_p"};
    varMember<float> elecrons_p {tree, "electrons_p"};
    varMember<float> d_12 {tree, "d_12"};
    varMember<float> d_23 {tree, "d_23"};
    varMember<float> d_34 {tree, "d_34"};

    // reconstructed jets (for flavor tagger score)
    varMember<ROOT::VecOps::RVec<float>> recojet_isB {tree, "recojet_isB"};
    varMember<ROOT::VecOps::RVec<float>> recojet_isC {tree, "recojet_isC"};
    varMember<ROOT::VecOps::RVec<float>> recojet_isS {tree, "recojet_isS"};
    varMember<ROOT::VecOps::RVec<float>> recojet_isU {tree, "recojet_isU"};
    varMember<ROOT::VecOps::RVec<float>> recojet_isD {tree, "recojet_isD"};
    varMember<ROOT::VecOps::RVec<float>> recojet_isG {tree, "recojet_isG"};
    varMember<ROOT::VecOps::RVec<float>> recojet_isTAU {tree, "recojet_isTAU"};

    // truth information
    varMember<ROOT::VecOps::RVec<int>> truth_Wm_Daugthers_pdg {tree, "truth_Wm_Daugthers_pdg"};
    varMember<ROOT::VecOps::RVec<int>> truth_Wp_Daugthers_pdg {tree, "truth_Wp_Daugthers_pdg"};

    // For the new input file
    varMember<ROOT::VecOps::RVec<float>> jet_px {tree, "jet_px_corr"};
    varMember<ROOT::VecOps::RVec<float>> jet_py {tree, "jet_py_corr"};
    varMember<ROOT::VecOps::RVec<float>> jet_pz {tree, "jet_pz_corr"};
    varMember<ROOT::VecOps::RVec<float>> jet_e {tree, "jet_e_corr"};

    ROOT::VecOps::RVec<int>* Wm_jet_truth = truth_Wm_Daugthers_pdg();
    ROOT::VecOps::RVec<int>* Wp_jet_truth = truth_Wp_Daugthers_pdg();

    varMember<ROOT::VecOps::RVec<float>> truth_Wm_p {tree, "truth_Wm_Daugthers_p"};
    varMember<ROOT::VecOps::RVec<float>> truth_Wm_e {tree, "truth_Wm_Daugthers_e"};
    varMember<ROOT::VecOps::RVec<float>> truth_Wm_theta {tree, "truth_Wm_Daugthers_theta"};
    varMember<ROOT::VecOps::RVec<float>> truth_Wm_phi {tree, "truth_Wm_Daugthers_phi"};

    varMember<ROOT::VecOps::RVec<float>> truth_Wp_p {tree, "truth_Wp_Daugthers_p"};
    varMember<ROOT::VecOps::RVec<float>> truth_Wp_e {tree, "truth_Wp_Daugthers_e"};
    varMember<ROOT::VecOps::RVec<float>> truth_Wp_theta {tree, "truth_Wp_Daugthers_theta"};
    varMember<ROOT::VecOps::RVec<float>> truth_Wp_phi {tree, "truth_Wp_Daugthers_phi"};

    // Jet Constituents
    varMember<ROOT::VecOps::RVec<int>> jetconstituents_kt4 {tree, "jetconstituents_kt4"};
    varMember<ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>>> jetconstituents_kt4_e {tree, "jetconstituents_kt4_e"};
    varMember<ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>>> jetconstituents_kt4_p {tree, "jetconstituents_kt4_p"};
    varMember<ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>>> jetconstituents_kt4_theta {tree, "jetconstituents_kt4_theta"};
    varMember<ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>>> jetconstituents_kt4_phi {tree, "jetconstituents_kt4_phi"};

    // Increment for CutFlow
    int NEvents = 0;
    int NokFlav = 0;
    int NjetCut = 0;
    int NleptonCut = 0;
    int eventNum = 1;
    int nFlavScore = 0; 
    int NdCutd123 = 0;
    int NdCutd34 = 0;
    int NdecayCuts = 0;

    const int maxPrint = 10;
    int nPrinted = 0;
    bool do_debug = true; // Set to false to turn off all debug printing

    // *************************** SAME AS BELOW, LOOKING OVER THE MACRO STEP BY STEP ***************************
    for(int i = 0; i < nEntries; i++)
    {
        treeCont->getEntry(i);
        countingHist->Fill(1);
        NEvents++;

        if(i % 10000 == 0) std::cout<<"Done i: "<<i<<" out of "<<nEntries<<std::endl;

        if (do_debug && nPrinted < maxPrint) {
            std::cout << "========== Event " << eventNum << " ==========\n";

            for (size_t i = 0; i < Wm_jet_truth->size(); ++i) {
                int pdg_code = Wm_jet_truth->at(i);
                std::cout << "W_1 daughter " << i << " PDG code: " << pdg_code << std::endl;
            }

            for (size_t i = 0; i < Wp_jet_truth->size(); ++i) {
                int pdg_code = Wp_jet_truth->at(i);
                std::cout << "W_2 daughter " << i << " PDG code: " << pdg_code << std::endl;
            }

            std::cout << "  " << std::endl;

            nPrinted++;
        }

        eventNum++;

        if(event_njet() != 4) continue;
        NjetCut++;

        bool flage_toss = false;
        if (recojet_isB.size() == 0){
            flage_toss = true;
        }
        if (flage_toss) continue;
        // define vectors 
        std::vector<double> vec_ee_corr;
        std::vector<double> vec_theta_corr;
        std::vector<double> vec_phi_corr;
        std::vector<double> vec_chi_corr;

        float jet0_scoreQ = std::max({recojet_isU.at(0), recojet_isD.at(0), recojet_isS.at(0)});
        float jet1_scoreQ = std::max({recojet_isU.at(1), recojet_isD.at(1), recojet_isS.at(1)});
        float jet2_scoreQ = std::max({recojet_isU.at(2), recojet_isD.at(2), recojet_isS.at(2)});
        float jet3_scoreQ = std::max({recojet_isU.at(3), recojet_isD.at(3), recojet_isS.at(3)});

        std::vector<float> j0_flav {recojet_isB.at(0), recojet_isC.at(0), jet0_scoreQ, recojet_isG.at(0), recojet_isTAU.at(0)};
        std::vector<float> j1_flav {recojet_isB.at(1), recojet_isC.at(1), jet1_scoreQ, recojet_isG.at(1), recojet_isTAU.at(1)};
        std::vector<float> j2_flav {recojet_isB.at(2), recojet_isC.at(2), jet2_scoreQ, recojet_isG.at(2), recojet_isTAU.at(2)};
        std::vector<float> j3_flav {recojet_isB.at(3), recojet_isC.at(3), jet3_scoreQ, recojet_isG.at(3), recojet_isTAU.at(3)};
       
       
        for (float j0_f : j0_flav) {
            if (std::isnan(j0_f)) {
                flage_toss = true;
                continue;
            }
        }
        for (float j1_f : j1_flav) {
            if (std::isnan(j1_f)) {
                flage_toss = true;
                continue;
            }
        }
        for (float j2_f : j2_flav) {
            if (std::isnan(j2_f)) {
                flage_toss = true;
                continue;
            }
        }
        for (float j3_f : j3_flav) {
            if (std::isnan(j3_f)) {
                flage_toss = true;
                continue;
            }
        }

        NokFlav++;

        // cut on the number of electrons
        // if (event_nmu() > 2.) continue;
        // if (event_nel() > 2.) continue;
        // if (muons_p() > 20.) continue;
        // if (elecrons_p() > 20.) continue;

        NleptonCut++;

        // if (d_12()<=15000. || d_12()>=58000.) continue;//missing!
        // if ((d_23()<=400.) || (d_23()>=18000.))continue;
        // NdCutd123++;
        // if ((d_34()<=100.) || (d_34()>=6000.))continue;
        // NdCutd34++; 

        std::map<int,std::vector<float>> jetFlavScores;
        jetFlavScores[0] = j0_flav;
        jetFlavScores[1] = j1_flav;
        jetFlavScores[2] = j2_flav;
        jetFlavScores[3] = j3_flav;

        std::vector<std::string> flavLabels = {"b", "c", "s", "l", "g", "tau"};

        if (do_debug && nPrinted < maxPrint) {
            for (const auto& [jetIdx, flavVec] : jetFlavScores) {
                std::cout << "Jet " << jetIdx << " flavor scores:\n";
                for (size_t i = 0; i < flavVec.size(); ++i) {
                    std::cout << "  " << flavLabels[i] << ": ";
                    if (std::abs(flavVec[i]) < 1e-6) {
                        std::cout << std::scientific << std::setprecision(2) << flavVec[i];
                    } else {
                        std::cout << std::fixed << std::setprecision(6) << flavVec[i];
                    }
                    std::cout << "\n";
                }
                std::cout << std::endl;
            }
        }

        // Get the max flavor scores
        auto j0_MaxScoreIt = std::max_element(j0_flav.begin(), j0_flav.end());
        auto j1_MaxScoreIt = std::max_element(j1_flav.begin(), j1_flav.end());
        auto j2_MaxScoreIt = std::max_element(j2_flav.begin(), j2_flav.end());
        auto j3_MaxScoreIt = std::max_element(j3_flav.begin(), j3_flav.end());

        int j0_maxScoreIdx = std::distance(j0_flav.begin(), j0_MaxScoreIt);
        int j1_maxScoreIdx = std::distance(j1_flav.begin(), j1_MaxScoreIt);
        int j2_maxScoreIdx = std::distance(j2_flav.begin(), j2_MaxScoreIt);
        int j3_maxScoreIdx = std::distance(j3_flav.begin(), j3_MaxScoreIt);

        if (do_debug && nPrinted < maxPrint) {
            std::cout << "Jet 0: " << flavLabels[j0_maxScoreIdx] << " = " << *j0_MaxScoreIt << "\n";
            std::cout << "Jet 1: " << flavLabels[j1_maxScoreIdx] << " = " << *j1_MaxScoreIt << "\n";
            std::cout << "Jet 2: " << flavLabels[j2_maxScoreIdx] << " = " << *j2_MaxScoreIt << "\n";
            std::cout << "Jet 3: " << flavLabels[j3_maxScoreIdx] << " = " << *j3_MaxScoreIt << "\n";
            std::cout << "      " << std::endl;
        }

        std::array<int, 4> maxScoreIdx {j0_maxScoreIdx,j1_maxScoreIdx,j2_maxScoreIdx,j3_maxScoreIdx};

       if (nPrinted < maxPrint) {
            std::cout << "maxScoreIdx: ";
            for (int ScoreIdx : maxScoreIdx) {
                std::cout << ScoreIdx << " ";
            }
            std::cout << std::endl;
            std::cout << "      " << std::endl;
        }

        std::map<int,std::vector<int>> jetFlavMaxScore;
        for (std::size_t i = 0; i < maxScoreIdx.size(); ++i){
            jetFlavMaxScore[maxScoreIdx[i]].push_back(i);
        }

        int n_c = 0;
        int n_l = 0;
        // int n_b = 0;

        bool has_invalid_flavor = false;

        for (int i = 0; i < 4; ++i) {
            if (maxScoreIdx[i] == 1) n_c++;
            else if (maxScoreIdx[i] == 2) n_l++;
        }

        if (has_invalid_flavor) continue;

        if (!(n_c == 1 && n_l == 3)) continue;

        nFlavScore++;

        // Start Pairing the jets in W1 and W2
        int cJet = -1;
        std::vector<int> lJets;

        for (int i = 0; i < 4; ++i) {
            if (maxScoreIdx[i] == 1) cJet = i;
            else if (maxScoreIdx[i] == 2) lJets.push_back(i); 
        }

        // Prevent jet reuse: ensure all jet indices are unique
        std::set<int> uniqueJets = {cJet, lJets[0], lJets[1], lJets[2]};

        if (cJet == -1 || lJets.size() != 3 || uniqueJets.size() < 4) {
            if (cJet == lJets[0] && cJet == lJets[1] && cJet == lJets[2] && 
                lJets[0] == lJets[1] && lJets[0] == lJets[2] && lJets[1] == lJets[2]) {

            std::cerr << "Error: Jet reuse detected or invalid jet assignment!" << std::endl;
            continue;

            }
        }

        // W pair 1: c + s and u+d
        std::pair<int, int> W1_pair = {cJet, lJets[0]};
        std::pair<int, int> W2_pair = {lJets[1], lJets[2]};

        if (do_debug && nPrinted < maxPrint) {
            std::cout << "W1 pair: Jet " << W1_pair.first << " (c), Jet " << W1_pair.second << " (s)" << std::endl;
            std::cout << "W2 pair: Jet " << W2_pair.first << " (l), Jet " << W2_pair.second << " (l)" << std::endl;
            std::cout << "      " << std::endl;
        }

        // ****************************** CALCULATIONS USING JET CONSTITUENTS ******************************
        std::vector<TLorentzVector> jetConstituents; // flattens the overall vector and allows us to open it up
        // jetconstituents_kt4.at(i): number of jet const. belogning to jet i. 
        for (int i = 0; i < jetconstituents_kt4.size(); ++i){
            // std::cout << "instance i: " << i<< std::endl;
            // std::cout << "jetconstituents_kt4_p.at(i): " << jetconstituents_kt4_p.at(i)<< std::endl;

            if (jetconstituents_kt4_p.at(i).size() != jetconstituents_kt4.at(i)){
                    std::cerr << "ERROR: Input mismatch — jet size doesn't match!" << std::endl;
                    std::exit(EXIT_FAILURE); 
                }
                for (int k = 0; k < jetconstituents_kt4_p.at(i).size(); ++k){
                    // std::cout << "jetconstituents_kt4_p.at(i).at(k): " << jetconstituents_kt4_p.at(i).at(k)<< std::endl;
                    float p = jetconstituents_kt4_p.at(i).at(k);
                    float theta = jetconstituents_kt4_theta.at(i).at(k);
                    float phi = jetconstituents_kt4_phi.at(i).at(k);
                    float e = jetconstituents_kt4_e.at(i).at(k);


                    float px = p * sin(theta) * cos(phi);
                    float py = p * sin(theta) * sin(phi);
                    float pz = p * cos(theta);

                    TLorentzVector vec;
                    vec.SetPxPyPzE(px, py, pz, e);
                    jetConstituents.push_back(vec);
                }
        }

        for (size_t m = 0; m < jetConstituents.size(); ++m){
            for (size_t n = m + 1; n < jetConstituents.size(); ++n){

                TLorentzVector subJet_1 = jetConstituents.at(m);
                TLorentzVector subJet_2 = jetConstituents.at(n);

                h_subJet_1_p->Fill(subJet_1.P());
                h_subJet_2_p->Fill(subJet_2.P());

                double jetSub_chi_1 = subJet_1.Angle(subJet_2.Vect());
                double jetSub_cos_chi_1 = cos(jetSub_chi_1);
                double eec_jetSub_chi_1 = 0.5 * (1 - jetSub_cos_chi_1);

                double jetSub_theta_2 = std::fabs(subJet_1.Theta() - subJet_2.Theta());
                double jetSub_cos_theta_2 = cos(jetSub_theta_2);
                double eec_jetSub_theta_2 = 0.5 * (1 - jetSub_cos_theta_2);

                double jetSub_eta_3 = std::fabs(subJet_1.Eta() - subJet_2.Eta());

                double jetSub_phi_4 = subJet_1.Phi() - subJet_2.Phi();
                double jetSub_cos_phi_4 = cos(jetSub_phi_4);
                double eec_jetSub_phi_4 = 0.5 * (1 - jetSub_cos_phi_4);

                vec_ee_corr.push_back(eec_jetSub_theta_2);
                vec_theta_corr.push_back(jetSub_theta_2);
                vec_phi_corr.push_back(jetSub_phi_4);
                vec_chi_corr.push_back(jetSub_chi_1);

                h_chi_subjet->Fill(jetSub_chi_1);
                h_cos_chi_subjet->Fill(jetSub_cos_chi_1);
                h_eec_chi_subjet->Fill(eec_jetSub_chi_1);
                
                h_theta_subjet->Fill(jetSub_theta_2);
                h_cos_theta_subjet->Fill(jetSub_cos_theta_2);
                h_eec_theta_subjet->Fill(eec_jetSub_theta_2);

                h_eta_subjet->Fill(jetSub_eta_3);

                h_phi_subjet->Fill(jetSub_phi_4);
                h_cos_phi_subjet->Fill(jetSub_cos_phi_4);
                h_eec_phi_subjet->Fill(eec_jetSub_phi_4);
            }
        }

        // ****************************** CALCULATIONS USING TRUTH QUARKS ******************************
        std::vector<TLorentzVector> truthJets_Wm, truthJets_Wp;

        for (size_t i = 0; i < truth_Wm_p()->size(); ++i) {
            float p = truth_Wm_p()->at(i);
            float e = truth_Wm_e()->at(i);
            float theta = truth_Wm_theta()->at(i);
            float phi = truth_Wm_phi()->at(i);

            float px = p * sin(theta) * cos(phi);
            float py = p * sin(theta) * sin(phi);
            float pz = p * cos(theta);
            
            TLorentzVector lv;

            lv.SetPxPyPzE(px, py, pz, e);
            truthJets_Wm.push_back(lv);
        }

        for (size_t i = 0; i < truth_Wp_p()->size(); ++i) {
            float p = truth_Wp_p()->at(i);
            float e = truth_Wp_e()->at(i);
            float theta = truth_Wp_theta()->at(i);
            float phi = truth_Wp_phi()->at(i);

            float px = p * sin(theta) * cos(phi);
            float py = p * sin(theta) * sin(phi);
            float pz = p * cos(theta);
            
            TLorentzVector lv;

            lv.SetPxPyPzE(px, py, pz, e);
            truthJets_Wp.push_back(lv);
        }

        if (truthJets_Wm.size() >= 2 && truthJets_Wp.size() >= 2) {
            TLorentzVector W1_jet1 = truthJets_Wm[0];
            TLorentzVector W1_jet2 = truthJets_Wm[1];
            TLorentzVector W2_jet1 = truthJets_Wp[0];
            TLorentzVector W2_jet2 = truthJets_Wp[1];

            TLorentzVector W1_truthJet = W1_jet1 + W1_jet2;
            TLorentzVector W2_truthJet = W2_jet1 + W2_jet2;

            h_W1_truth_mass->Fill(W1_truthJet.M());
            h_W2_truth_mass->Fill(W2_truthJet.M());

            h_W1_truth_p->Fill(W1_truthJet.P());
            h_W2_truth_p->Fill(W2_truthJet.P());

            h_W1_truth_e->Fill(W1_truthJet.E());
            h_W2_truth_e->Fill(W2_truthJet.E());

            double chi_truth_cs = W1_jet1.Angle(W1_jet2.Vect()); // Is this the proper variable name for the angle?
            double chi_truth_ud = W2_jet1.Angle(W2_jet2.Vect()); // If not, please note the correct varaible name

            double theta_truth_cs = std::fabs(W1_jet1.Theta() - W1_jet2.Theta());
            double theta_truth_ud = std::fabs(W2_jet1.Theta() - W2_jet2.Theta());

            double eta_truth_cs = std::fabs(W1_jet1.Eta() - W1_jet2.Eta());
            double eta_truth_ud = std::fabs(W2_jet1.Eta() - W2_jet2.Eta());
            
            double phi_truth_cs = W1_jet1.Phi() - W1_jet2.Phi();
            double phi_truth_ud = W2_jet1.Phi() - W2_jet2.Phi();

            // ************** COSINE MEASUREMENTS **************
            double cos_chi_truth_cs = cos(chi_truth_cs);
            double cos_chi_truth_ud = cos(chi_truth_ud);

            double cos_theta_truth_cs = cos(theta_truth_cs);
            double cos_theta_truth_ud = cos(theta_truth_ud);

            double cos_phi_truth_cs = cos(phi_truth_cs);
            double cos_phi_truth_ud = cos(phi_truth_ud);



            // if (nPrinted < maxPrint) {
            //         std::cout << "Truth Quarks" << std::endl;
            //         std::cout << "chi: " << chi_truth_cs << " " << "cos(chi): " << cos_chi_truth_cs << std::endl;
            //         std::cout << "chi: " << chi_truth_ud << " " << "cos(chi): " << cos_chi_truth_ud << std::endl;
            //         std::cout << "theta: " << theta_truth_cs << " " << "cos(chi): " << cos_theta_truth_cs << std::endl;
            //         std::cout << "theta: " << theta_truth_ud << " " << "cos(chi): " << cos_theta_truth_ud << std::endl;
            //         std::cout << "phi: " << phi_truth_cs << " " << "cos(chi): " << cos_phi_truth_cs << std::endl;
            //         std::cout << "phi: " << phi_truth_ud << " " << "cos(chi): " << cos_phi_truth_ud << std::endl;
            //         std::cout << "  " << std::endl;
            // }

            // ************** EE CORRELATIONS **************
            double ee_correlation_chi_truth_cs = 0.5 * (1 - cos_chi_truth_cs);
            double ee_correlation_chi_truth_ud = 0.5 * (1 - cos_chi_truth_ud);

            double ee_correlation_theta_truth_cs = 0.5 * (1 - cos_theta_truth_cs);
            double ee_correlation_theta_truth_ud = 0.5 * (1 - cos_theta_truth_ud);

            double ee_correlation_phi_truth_cs = 0.5 * (1 - cos_phi_truth_cs);
            double ee_correlation_phi_truth_ud = 0.5 * (1 - cos_phi_truth_ud);

            // *********************** HISTOGRAMS ***********************
            h_chi_truth_cs->Fill(chi_truth_cs);
            h_chi_truth_ud->Fill(chi_truth_ud);
            h_theta_truth_cs->Fill(theta_truth_cs);
            h_theta_truth_ud->Fill(theta_truth_ud);
            h_eta_truth_cs->Fill(eta_truth_cs);
            h_eta_truth_ud->Fill(eta_truth_ud);
            h_phi_truth_cs->Fill(phi_truth_cs);
            h_phi_truth_ud->Fill(phi_truth_ud);

            h_cos_chi_truth_cs->Fill(cos_chi_truth_cs);
            h_cos_chi_truth_ud->Fill(cos_chi_truth_ud);
            h_cos_theta_truth_cs->Fill(cos_theta_truth_cs);
            h_cos_theta_truth_ud->Fill(cos_theta_truth_ud);
            h_cos_phi_truth_cs->Fill(cos_phi_truth_cs);
            h_cos_phi_truth_ud->Fill(cos_phi_truth_ud);

            h_eec_chi_truth_cs->Fill(ee_correlation_chi_truth_cs);
            h_eec_chi_truth_ud->Fill(ee_correlation_chi_truth_ud);
            h_eec_theta_truth_cs->Fill(ee_correlation_theta_truth_cs);
            h_eec_theta_truth_ud->Fill(ee_correlation_theta_truth_ud);
            h_eec_phi_truth_cs->Fill(ee_correlation_phi_truth_cs);
            h_eec_phi_truth_ud->Fill(ee_correlation_phi_truth_ud);
        }

        // ************************** CALCULATIONS FOR MASS AND SUCH BEGIN HERE (USING RECO-JETS) ************************** 
        TLorentzVector cTag_Jet, lTag_Jet_0, lTag_Jet_1, lTag_Jet_2;

        cTag_Jet.SetPxPyPzE(jet_px.at(cJet), jet_py.at(cJet), jet_pz.at(cJet), jet_e.at(cJet));
        lTag_Jet_0.SetPxPyPzE(jet_px.at(lJets[0]), jet_py.at(lJets[0]), jet_pz.at(lJets[0]), jet_e.at(lJets[0]));

        lTag_Jet_1.SetPxPyPzE(jet_px.at(lJets[1]), jet_py.at(lJets[1]), jet_pz.at(lJets[1]), jet_e.at(lJets[1]));
        lTag_Jet_2.SetPxPyPzE(jet_px.at(lJets[2]), jet_py.at(lJets[2]), jet_pz.at(lJets[2]), jet_e.at(lJets[2]));
        
        const double m_W_true = 80.379;

        TLorentzVector W1_option1 = cTag_Jet + lTag_Jet_0;
        TLorentzVector W2_option1 = lTag_Jet_1 + lTag_Jet_2;

        TLorentzVector W1_option2 = cTag_Jet + lTag_Jet_1;
        TLorentzVector W2_option2 = lTag_Jet_0 + lTag_Jet_2;

        TLorentzVector W1_option3 = cTag_Jet + lTag_Jet_2;
        TLorentzVector W2_option3 = lTag_Jet_0 + lTag_Jet_1;

        double chi2_option1 = 
            (pow(W1_option1.M() - m_W_true, 2) + pow(W2_option1.M() - m_W_true, 2)) / (m_W_true);

        double chi2_option2 = 
            (pow(W1_option2.M() - m_W_true, 2) + pow(W2_option2.M() - m_W_true, 2)) / (m_W_true);

        double chi2_option3 = 
            (pow(W1_option3.M() - m_W_true, 2) + pow(W2_option3.M() - m_W_true, 2)) / (m_W_true);

        // for the cos(phi) angle
        TLorentzVector W1_j1, W1_j2, W2_j1, W2_j2;

        if (chi2_option1 <= chi2_option2 && chi2_option1 <= chi2_option3) {
            h_W1_mass->Fill(W1_option1.M());
            h_W2_mass->Fill(W2_option1.M());
            h_chi2->Fill(chi2_option1);

            h_W1_p->Fill(W1_option1.P());
            h_W2_p->Fill(W2_option1.P());
            
            h_W1_e->Fill(W1_option1.E());
            h_W2_e->Fill(W2_option1.E());

            h_cJet_p->Fill(cTag_Jet.P());
            h_lJet0_p->Fill(lTag_Jet_0.P());
            h_lJet1_p->Fill(lTag_Jet_1.P());
            h_lJet2_p->Fill(lTag_Jet_2.P());

            W1_j1 = cTag_Jet; 
            W1_j2 = lTag_Jet_0;

            W2_j1 = lTag_Jet_1;   
            W2_j2 = lTag_Jet_2;

        } else if (chi2_option2 <= chi2_option1 && chi2_option2 <= chi2_option3) {
            h_W1_mass->Fill(W1_option2.M());
            h_W2_mass->Fill(W2_option2.M());
            h_chi2->Fill(chi2_option2);

            h_W1_p->Fill(W1_option2.P());
            h_W2_p->Fill(W2_option2.P());

            h_W1_e->Fill(W1_option2.E());
            h_W2_e->Fill(W2_option2.E());

            h_cJet_p->Fill(cTag_Jet.P());
            h_lJet0_p->Fill(lTag_Jet_0.P());
            h_lJet1_p->Fill(lTag_Jet_1.P());
            h_lJet2_p->Fill(lTag_Jet_2.P());

            W1_j1 = cTag_Jet; 
            W1_j2 = lTag_Jet_1;

            W2_j1 = lTag_Jet_0;   
            W2_j2 = lTag_Jet_2;

        } else {
            h_W1_mass->Fill(W1_option3.M());
            h_W2_mass->Fill(W2_option3.M());
            h_chi2->Fill(chi2_option3);

            h_W1_p->Fill(W1_option3.P());
            h_W2_p->Fill(W2_option3.P());

            h_W1_e->Fill(W1_option3.E());
            h_W2_e->Fill(W2_option3.E());

            h_cJet_p->Fill(cTag_Jet.P());
            h_lJet0_p->Fill(lTag_Jet_0.P());
            h_lJet1_p->Fill(lTag_Jet_1.P());
            h_lJet2_p->Fill(lTag_Jet_2.P());

            W1_j1 = cTag_Jet; 
            W1_j2 = lTag_Jet_2;

            W2_j1 = lTag_Jet_1;   
            W2_j2 = lTag_Jet_0;
        }

        mc_weight = norm_weight;
        ee_corr = vec_ee_corr;
        theta_corr = vec_theta_corr;
        phi_corr= vec_phi_corr;
        chi_corr=vec_chi_corr;
        my_tree->Fill();

        // ************* CHI *************
        double chi_c_l0 = W1_j1.Angle(W1_j2.Vect());
        double chi_l1_l2 = W2_j1.Angle(W2_j2.Vect());

        double cos_chi_c_l0 = cos(chi_c_l0);
        double cos_chi_l1_l2 = cos(chi_l1_l2);

        double ee_correlation_chi_c_l0 = 0.5 * (1 - cos_chi_c_l0);
        double ee_correlation_chi_l1_l2 = 0.5 * (1 - cos_chi_l1_l2);

        // ************* THETA *************
        double theta_c_l0 = std::fabs(W1_j1.Theta() - W1_j2.Theta());
        double theta_l1_l2 = std::fabs(W2_j1.Theta() - W2_j2.Theta());

        double cos_theta_c_l0 = cos(theta_c_l0);
        double cos_theta_l1_l2 = cos(theta_l1_l2);

        double ee_correlation_theta_c_l0 = 0.5 * (1 - cos_theta_c_l0);
        double ee_correlation_theta_l1_l2 = 0.5 * (1 - cos_theta_l1_l2);

        // ************* ETA *************
        double eta_c_l0 = std::fabs(W1_j1.Eta() - W1_j2.Eta());
        double eta_l1_l2 = std::fabs(W2_j1.Eta() - W2_j2.Eta());
        
        // ************* PHI *************
        double phi_c_l0 = W1_j1.Phi() - W1_j2.Phi();
        double phi_l1_l2 = W2_j1.Phi() - W2_j2.Phi();

        double cos_phi_c_l0 = cos(phi_c_l0);
        double cos_phi_l1_l2 = cos(phi_l1_l2);

        double ee_correlation_phi_c_l0 = 0.5 * (1 - cos_phi_c_l0);
        double ee_correlation_phi_l1_l2 = 0.5 * (1 - cos_phi_l1_l2);

        // *********************** HISTOGRAMS ***********************
        h_chi_c_l0->Fill(chi_c_l0);
        h_chi_l1_l2->Fill(chi_l1_l2);
        h_theta_c_l0->Fill(theta_c_l0);
        h_theta_l1_l2->Fill(theta_l1_l2);
        h_eta_c_l0->Fill(eta_c_l0);
        h_eta_l1_l2->Fill(eta_l1_l2);
        h_phi_c_l0->Fill(phi_c_l0);
        h_phi_l1_l2->Fill(phi_l1_l2);

        h_cos_chi_c_l0->Fill(cos_chi_c_l0);
        h_cos_chi_l1_l2->Fill(cos_chi_l1_l2);
        h_cos_theta_c_l0->Fill(cos_theta_c_l0);
        h_cos_theta_l1_l2->Fill(cos_theta_l1_l2);
        h_cos_phi_c_l0->Fill(cos_phi_c_l0);
        h_cos_phi_l1_l2->Fill(cos_phi_l1_l2);

        h_eec_chi_c_l0->Fill(ee_correlation_chi_c_l0);
        h_eec_chi_l1_l2->Fill(ee_correlation_chi_l1_l2);
        h_eec_theta_c_l0->Fill(ee_correlation_theta_c_l0);
        h_eec_theta_l1_l2->Fill(ee_correlation_theta_l1_l2);
        h_eec_phi_c_l0->Fill(ee_correlation_phi_c_l0);
        h_eec_phi_l1_l2->Fill(ee_correlation_phi_l1_l2);


        // cutflow histograms
        cutFlowHist->SetBinContent(1, NEvents);
        cutFlowHist->SetBinContent(2, NdecayCuts);
        cutFlowHist->SetBinContent(3, NokFlav);
        cutFlowHist->SetBinContent(4, NjetCut);
        cutFlowHist->SetBinContent(5, NleptonCut);
        cutFlowHist->SetBinContent(6, NdCutd123);
        cutFlowHist->SetBinContent(7, NdCutd34);
        cutFlowHist->SetBinContent(8, nFlavScore);
    }

    std::cout << "      " << std::endl;
    std::cout << "-------------------- Outputs --------------------" << std::endl;
    std::cout << "Number of events: " << NEvents << std::endl;
    std::cout << "Number of events passing W decay cuts: " << NdecayCuts << std::endl;
    std::cout << "Number of events w/ 4 jets: " << NjetCut << std::endl;
    std::cout << "Number of Leptons Cut: " << NleptonCut << std::endl;
    std::cout << "Number of events with 1 c-tagged, 1 s-tagged, and 2 light-tagged jets: " << nFlavScore << std::endl;
    std::cout << "      " << std::endl;

    std::cout << "      " << std::endl;
    std::cout << "Let there be data :)" << std::endl;

    // end of macro
}