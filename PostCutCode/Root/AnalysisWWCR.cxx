// Local includes
#include "PostCutCode/AnalysisWWCR.h"
#include "PostCutCode/TreeContainer.h"
#include <iomanip>
#include <TFile.h>
#include <fstream>
#include <TRandom.h>
#include "TLorentzVector.h"

#include "TH1F.h"
#include "TCanvas.h"


using namespace std;

// Base class 
AnalysisWWCR::AnalysisWWCR():
    AnalysisBase(),
    m_debug(true) 
{}
AnalysisWWCR::~AnalysisWWCR()
{}


void AnalysisWWCR::run()
{ 

    std::vector<std::string> cutFlowMap {"All Events", "NoNaNFlavScore", "Has 4 Jets", "leptonCut", "Events have a c, s, and 2 l tagged jets"};

    // Get the histograms
    auto countingHist = m_histContainer->getCountingHist();
    // auto truth_W_e =  m_histContainer->get1DHist("truth_W_e", 100, 0, 200);
    auto cutFlowHist = m_histContainer->get1DHist("cutFlowHist", cutFlowMap.size(), 0, 5, cutFlowMap);
    auto h_W1_mass = m_histContainer->get1DHist("h_W1_mass", 300, 0, 150);
    auto h_W2_mass = m_histContainer->get1DHist("h_W2_mass", 300, 0, 150);
    auto h_chi2 = m_histContainer->get1DHist("h_chi2", 300, 0, 150);

    // auto h_W1_truth_mass = m_histContainer->get1DHist("h_W1_truth_mass", 300, 0, 150);
    // auto h_W2_truth_mass = m_histContainer->get1DHist("h_W2_truth_mass", 300, 0, 150);
    // auto h_chi2_truth = m_histContainer->get1DHist("h_chi2_truth", 300, 0, 150);
    // auto h_W2_truth_mass = m_histContainer->get1DHist("h_W2_truth_mass", 300, 0, 150);

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

    // varMember<ROOT::VecOps::RVec<float>>  truth_Wp_HS_e {tree, "truth_Wp_HS_e"};

    // reconstructed jets (for flavor tagger score)
    varMember<ROOT::VecOps::RVec<float>> recojet_isB {tree, "recojet_isB"};
    varMember<ROOT::VecOps::RVec<float>> recojet_isC {tree, "recojet_isC"};
    varMember<ROOT::VecOps::RVec<float>> recojet_isS {tree, "recojet_isS"};
    varMember<ROOT::VecOps::RVec<float>> recojet_isU {tree, "recojet_isU"};
    varMember<ROOT::VecOps::RVec<float>> recojet_isD {tree, "recojet_isD"};
    varMember<ROOT::VecOps::RVec<float>> recojet_isG {tree, "recojet_isG"};
    varMember<ROOT::VecOps::RVec<float>> recojet_isTAU {tree, "recojet_isTAU"};
    // truth jet flav (from the Monte Carlo simulation4)
    varMember<ROOT::VecOps::RVec<int>> truth_flav {tree, "jets_truth"};

    // Mass of W boson
    // float W_mass = 80.3692;

    // For the new input file
    varMember<ROOT::VecOps::RVec<float>> jet_px {tree, "jet_px_corr"};
    varMember<ROOT::VecOps::RVec<float>> jet_py {tree, "jet_py_corr"};
    varMember<ROOT::VecOps::RVec<float>> jet_pz {tree, "jet_pz_corr"};
    varMember<ROOT::VecOps::RVec<float>> jet_e {tree, "jet_e_corr"};

    // Increment for CutFlow
    int NEvents = 0;
    int NokFlav = 0;
    int NjetCut = 0;
    int NleptonCut = 0;
    int eventNum = 1;
    int nFlavScore = 0; 

    // int NEvents_truth = 0;
    // int nFlavScore_truth = 0; 

    const int maxPrint = 10;
    int nPrinted = 0;
    bool do_debug = true; // Set to false to turn off all debug printing

    // Main event loop 
    for(int i = 0; i < nEntries; i++)
    {
        treeCont->getEntry(i);
        // Just to store how many events were run over
        countingHist->Fill(1);
        NEvents++;

        // ROOT::VecOps::RVec<int>* jet_truth = truth_flav();

        if(i % 10000 == 0) std::cout<<"Done i: "<<i<<" out of "<<nEntries<<std::endl;

        // Adds the requirement that an event has 4 jets
        if(event_njet() != 4) continue;
        NjetCut++;

        ROOT::VecOps::RVec<int>* jet_truth = truth_flav();

        if (jet_truth->size() < 4) {
            std::cout << "Warning: jets_truth has fewer than 4 entries!" << std::endl;
            continue;
        }

        int jet0_truth = jet_truth->at(0);
        int jet1_truth = jet_truth->at(1);
        int jet2_truth = jet_truth->at(2);
        int jet3_truth = jet_truth->at(3);

        NEvents_truth++;

        bool flage_toss = false;
        if (recojet_isB.size() == 0){
            flage_toss = true;
        }
        if (flage_toss) continue;

        // Find the tagger score
        // both D and U are consisdered "light flavors" and the tagger struggles to differentiate the two
        float jet0_scoreQ = recojet_isU.at(0) < recojet_isD.at(0)? recojet_isD.at(0) :recojet_isU.at(0);
        float jet1_scoreQ = recojet_isU.at(1) < recojet_isD.at(1)? recojet_isD.at(1) :recojet_isU.at(1);
        float jet2_scoreQ = recojet_isU.at(2) < recojet_isD.at(2)? recojet_isD.at(2) :recojet_isU.at(2);
        float jet3_scoreQ = recojet_isU.at(3) < recojet_isD.at(3)? recojet_isD.at(3) :recojet_isU.at(3);

        std::vector<float> j0_flav {recojet_isB.at(0), recojet_isC.at(0), recojet_isS.at(0), jet0_scoreQ, recojet_isG.at(0), recojet_isTAU.at(0)};
        std::vector<float> j1_flav {recojet_isB.at(1), recojet_isC.at(1), recojet_isS.at(1), jet1_scoreQ, recojet_isG.at(1), recojet_isTAU.at(1)};
        std::vector<float> j2_flav {recojet_isB.at(2), recojet_isC.at(2), recojet_isS.at(2), jet2_scoreQ, recojet_isG.at(2), recojet_isTAU.at(2)};
        std::vector<float> j3_flav {recojet_isB.at(3), recojet_isC.at(3), recojet_isS.at(3), jet3_scoreQ, recojet_isG.at(3), recojet_isTAU.at(3)};

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

        // // cut on the number of electrons
        // if (event_nmu() > 2.) continue;
        // if (event_nel() > 2.) continue;
        // if (muons_p() > 20.) continue;
        // if (elecrons_p() > 20.) continue;

        // NleptonCut++;

        std::map<int,std::vector<float>> jetFlavScores;
        jetFlavScores[0] = j0_flav;
        jetFlavScores[1] = j1_flav;
        jetFlavScores[2] = j2_flav;
        jetFlavScores[3] = j3_flav;

        std::vector<std::string> flavLabels = {"b", "c", "s", "l", "g", "tau"};

        if (do_debug && nPrinted < maxPrint) {
            std::cout << "========== Event " << eventNum << " ==========\n";

            std::cout << "Truth Flavors: " << std::endl;
            std::cout << "Jet 0 flavor = " << jet0_truth << std::endl;
            std::cout << "Jet 1 flavor = " << jet1_truth << std::endl;
            std::cout << "Jet 2 flavor = " << jet2_truth << std::endl;
            std::cout << "Jet 3 flavor = " << jet3_truth << std::endl;
            std::cout << "  " << std::endl;

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

            nPrinted++; // Count only after printing all blocks for this event
        }

        eventNum++;

        // ************* THIS CALCULATES THE MASS WITH THE RECONSTRUCTED JETS *************
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
        int n_s = 0;
        int n_l = 0;

        bool has_invalid_flavor = false;

        for (int i = 0; i < 4; ++i) {
            if (maxScoreIdx[i] == 1) n_c++; // c-tag
            else if (maxScoreIdx[i] == 2) n_s++; // s-tag
            else if (maxScoreIdx[i] == 3) n_l++; // light-tag
        }

        if (has_invalid_flavor) continue;

        if (!(n_c == 1 && n_s == 1 && n_l == 2)) continue;

        nFlavScore++; // Add this as a 3rd bin in cutflow

        // Start Pairing the jets in W1 and W2
        int cJet = -1;
        int sJet = -1;
        std::vector<int> lJets;

        for (int i = 0; i < 4; ++i) {
            if (maxScoreIdx[i] == 1) cJet = i;
            else if (maxScoreIdx[i] == 2) sJet = i;
            else if (maxScoreIdx[i] == 3) lJets.push_back(i); 
        }

        if (cJet == -1 || sJet == -1 || lJets.size() != 2) {
            std::cerr << "Error: Incorrect number of flavored jets for pairing!" << std::endl;
            continue;
        }

        // W pair 1: c + s
        std::pair<int, int> W1_pair = {cJet, sJet};

        // W pair 2: l + l
        std::pair<int, int> W2_pair = {lJets[0], lJets[1]};

        if (do_debug && nPrinted < maxPrint) {
            std::cout << "W1 pair: Jet " << W1_pair.first << " (c), Jet " << W1_pair.second << " (s)" << std::endl;
            std::cout << "W2 pair: Jet " << W2_pair.first << " (l), Jet " << W2_pair.second << " (l)" << std::endl;
            std::cout << "      " << std::endl;
        }

        // calculate the mass 
        TLorentzVector cTag_Jet, sTag_Jet, lTag_Jet_0, lTag_Jet_1;

        cTag_Jet.SetPxPyPzE(jet_px.at(cJet), jet_py.at(cJet), jet_py.at(cJet), jet_e.at(cJet));
        sTag_Jet.SetPxPyPzE(jet_px.at(sJet), jet_py.at(sJet), jet_py.at(sJet), jet_e.at(sJet));

        lTag_Jet_0.SetPxPyPzE(jet_px.at(lJets[0]), jet_py.at(lJets[0]), jet_py.at(lJets[0]), jet_e.at(lJets[0]));
        lTag_Jet_1.SetPxPyPzE(jet_px.at(lJets[1]), jet_py.at(lJets[1]), jet_py.at(lJets[1]), jet_e.at(lJets[1]));

        // ************* THIS CALCULATES THE MASS WITH THE TRUTH JETS USING THE CHI^2 METHOD *************

        const double m_W_true = 80.379;
        // const double sigma_W = 10.0; <-- should this be the mass of W in the chi^2 calculations

        TLorentzVector W1_option1 = cTag_Jet + sTag_Jet;
        TLorentzVector W2_option1 = lTag_Jet_0 + lTag_Jet_1;

        TLorentzVector W1_option2 = cTag_Jet + lTag_Jet_0;
        TLorentzVector W2_option2 = sTag_Jet + lTag_Jet_1;

        TLorentzVector W1_option3 = cTag_Jet + lTag_Jet_1;
        TLorentzVector W2_option3 = sTag_Jet + lTag_Jet_0;

        // double chi2_option1 = 
        //     (pow(W1_option1.M() - m_W_true, 2) + pow(W2_option1.M() - m_W_true, 2));

        // double chi2_option2 = 
        //     (pow(W1_option2.M() - m_W_true, 2) + pow(W2_option2.M() - m_W_true, 2));

        double chi2_option1 = 
            (pow(W1_option1.M() - m_W_true, 2) + pow(W2_option1.M() - m_W_true, 2)) / (m_W_true);

        double chi2_option2 = 
            (pow(W1_option2.M() - m_W_true, 2) + pow(W2_option2.M() - m_W_true, 2)) / (m_W_true);

        double chi2_option3 = 
            (pow(W1_option3.M() - m_W_true, 2) + pow(W2_option3.M() - m_W_true, 2)) / (m_W_true);

        // double chi2_option1 = 
        //     (pow(W1_option1.M() - m_W_true, 2) + pow(W2_option1.M() - m_W_true, 2)) / (sigma_W * sigma_W);

        // double chi2_option2 = 
        //     (pow(W1_option2.M() - m_W_true, 2) + pow(W2_option2.M() - m_W_true, 2)) / (sigma_W * sigma_W);

        if (chi2_option1 <= chi2_option2) {
            // Option 1 is chosen
            h_W1_mass->Fill(W1_option1.M());
            h_W2_mass->Fill(W2_option1.M());
            h_chi2->Fill(chi2_option1);
        } else if (chi2_option2 <= chi2_option3) {
            // Option 2 is chosen
            h_W1_mass->Fill(W1_option2.M());
            h_W2_mass->Fill(W2_option2.M());
            h_chi2->Fill(chi2_option2);
        } else {
            // Option 2 is chosen
            h_W1_mass->Fill(W1_option3.M());
            h_W2_mass->Fill(W2_option3.M());
            h_chi2->Fill(chi2_option3);
        }

        // TLorentzVector W1 = cTag_Jet + sTag_Jet;
        // TLorentzVector W2 = lTag_Jet_0 + lTag_Jet_1;

        // h_W1_mass->Fill(W1.M());
        // h_W2_mass->Fill(W2.M());

        // if (do_debug && nPrinted < maxPrint) {
        //     std::cout << "W1 mass: " << W1.M() << " GeV" << std::endl;
        //     std::cout << "W2 mass: " << W2.M() << " GeV" << std::endl;
        //     std::cout << "  " << std::endl;
        // }



        // // **************************** THIS CALCULATES THE MASS WITH THE TRUTH JETS ****************************
        // // Pair the truth jets
        // int nTruth_c = 0;
        // int nTruth_s = 0;
        // int nTruth_u = 0;
        // int nTruth_d = 0;

        // bool truth_has_invalid_flavor = false;

        // for (int i = 0; i < 4; ++i) {
        //     int flav = jet_truth->at(i);

        //     // PDG flavor codes for quarks and anti-quarks:
        //     if (flav == 4 || flav == -4) nTruth_c++;
        //     else if (flav == 3 || flav == -3) nTruth_s++;
        //     else if (flav == 2 || flav == -2) nTruth_u++;
        //     else if (flav == 1 || flav == -1) nTruth_d++;
        //     else {
        //         truth_has_invalid_flavor = true;
        //         break;  // stop counting if invalid flavor found
        //     }
        // }

        // if (truth_has_invalid_flavor) continue;

        // if (!(nTruth_c == 1 && nTruth_s == 1 && nTruth_u == 1 && nTruth_d == 1)) continue;

        // nFlavScore_truth++;

        // int cJet_truth = -1;
        // int sJet_truth = -1;
        // int uJet_truth = -1;
        // int dJet_truth = -1;

        // for (int i = 0; i < 4; ++i) {
        //     int flav = jet_truth->at(i);

        //     if (flav == 4 || flav == -4) cJet_truth = i;
        //     else if (flav == 3 || flav == -3) sJet_truth = i;
        //     else if (flav == 2 || flav == -2) uJet_truth = i;
        //     else if (flav == 1 || flav == -1) dJet_truth = i;
        // }

        // if (cJet_truth == -1 || sJet_truth == -1 || uJet_truth == -1 || dJet_truth == -1) {
        //     std::cerr << "Error: Incorrect number of flavored jets for pairing!" << std::endl;
        //     continue;
        // }

        // // ************* LOOKS AT EVENTS WITH THE PAIRS (U, D, S, C) AND (U, D, C, C) *************

        // // if (!((nTruth_c == 1 && nTruth_s == 1 && nTruth_u == 1 && nTruth_d == 1) ||
        // // (nTruth_c == 2 && nTruth_u == 1 && nTruth_d == 1)))
        // // continue;

        // // nFlavScore_truth++;

        // // int cJet1_truth = -1;
        // // int cJet2_truth = -1;
        // // int sJet_truth = -1;
        // // int uJet_truth = -1;
        // // int dJet_truth = -1;

        // // int c_count = 0;

        // // for (int i = 0; i < 4; ++i) {
        // //     int flav = jet_truth->at(i);

        // //     if (flav == 4 || flav == -4) {
        // //         if (c_count == 0) cJet1_truth = i;
        // //         else if (c_count == 1) cJet2_truth = i;
        // //         c_count++;
        // //     }
        // //     else if (flav == 3 || flav == -3) sJet_truth = i;
        // //     else if (flav == 2 || flav == -2) uJet_truth = i;
        // //     else if (flav == 1 || flav == -1) dJet_truth = i;
        // // }

        // // if (!(
        // //         (c_count == 1 && sJet_truth != -1 && uJet_truth != -1 && dJet_truth != -1)  // case 1
        // //         ||
        // //         (c_count == 2 && uJet_truth != -1 && dJet_truth != -1) // case 2
        // //     )) 
        // // {
        // //     std::cerr << "Error: Incorrect number of flavored jets for pairing!" << std::endl;
        // //     continue;
        // // }

        // // ************************************ END OF EDITS ************************************


        // // W pair 1: c + s
        // std::pair<int, int> W1_pair_truth = {cJet_truth, sJet_truth};

        // // W pair 2: u + d
        // std::pair<int, int> W2_pair_truth = {uJet_truth, dJet_truth};

        // if (do_debug && nPrinted < maxPrint) {
        //     std::cout << "W1 pair: Jet " << W1_pair_truth.first << " (c), Jet " << W1_pair_truth.second << " (s)" << std::endl;
        //     std::cout << "W2 pair: Jet " << W2_pair_truth.first << " (u), Jet " << W2_pair_truth.second << " (d)" << std::endl;
        //     std::cout << "      " << std::endl;
        // }

        // // calculate the mass 
        // TLorentzVector cTag_truthJet, sTag_truthJet, uTag_truthJet, dTag_truthJet;

        // cTag_truthJet.SetPxPyPzE(jet_px.at(cJet_truth), jet_py.at(cJet_truth), jet_py.at(cJet_truth), jet_e.at(cJet_truth));
        // sTag_truthJet.SetPxPyPzE(jet_px.at(sJet_truth), jet_py.at(sJet_truth), jet_py.at(sJet_truth), jet_e.at(sJet_truth));
        // uTag_truthJet.SetPxPyPzE(jet_px.at(uJet_truth), jet_py.at(uJet_truth), jet_py.at(uJet_truth), jet_e.at(uJet_truth));
        // dTag_truthJet.SetPxPyPzE(jet_px.at(dJet_truth), jet_py.at(dJet_truth), jet_py.at(dJet_truth), jet_e.at(dJet_truth));


        // // TLorentzVector W1_truth = cTag_truthJet + sTag_truthJet;
        // // TLorentzVector W2_truth = uTag_truthJet + dTag_truthJet;

        // // h_W1_truth_mass->Fill(W1_truth.M());
        // // h_W2_truth_mass->Fill(W2_truth.M());

        // // if (do_debug && nPrinted < maxPrint) {
        // //     std::cout << "W1 mass: " << W1_truth.M() << " GeV" << std::endl;
        // //     std::cout << "W2 mass: " << W2_truth.M() << " GeV" << std::endl;
        // //     std::cout << "  " << std::endl;
        // // }


        // // ************* THIS CALCULATES THE MASS WITH THE TRUTH JETS USING THE CHI^2 METHOD *************

        // const double m_W_true = 80.379;
        // // const double sigma_W = 10.0; <-- should this be the mass of W in the chi^2 calculations

        // TLorentzVector W1_option1 = cTag_truthJet + sTag_truthJet;
        // TLorentzVector W2_option1 = uTag_truthJet + dTag_truthJet;

        // TLorentzVector W1_option2 = cTag_truthJet + dTag_truthJet;
        // TLorentzVector W2_option2 = uTag_truthJet + sTag_truthJet;

        // TLorentzVector W1_option3 = cTag_truthJet + uTag_truthJet;
        // TLorentzVector W2_option3 = sTag_truthJet + dTag_truthJet;

        // // double chi2_option1 = 
        // //     (pow(W1_option1.M() - m_W_true, 2) + pow(W2_option1.M() - m_W_true, 2));

        // // double chi2_option2 = 
        // //     (pow(W1_option2.M() - m_W_true, 2) + pow(W2_option2.M() - m_W_true, 2));

        // double chi2_option1 = 
        //     (pow(W1_option1.M() - m_W_true, 2) + pow(W2_option1.M() - m_W_true, 2)) / (m_W_true);

        // double chi2_option2 = 
        //     (pow(W1_option2.M() - m_W_true, 2) + pow(W2_option2.M() - m_W_true, 2)) / (m_W_true);

        // double chi2_option3 = 
        //     (pow(W1_option3.M() - m_W_true, 2) + pow(W2_option3.M() - m_W_true, 2)) / (m_W_true);

        // // double chi2_option1 = 
        // //     (pow(W1_option1.M() - m_W_true, 2) + pow(W2_option1.M() - m_W_true, 2)) / (sigma_W * sigma_W);

        // // double chi2_option2 = 
        // //     (pow(W1_option2.M() - m_W_true, 2) + pow(W2_option2.M() - m_W_true, 2)) / (sigma_W * sigma_W);

        // if (chi2_option1 <= chi2_option2) {
        //     // Option 1 is chosen
        //     h_W1_truth_mass->Fill(W1_option1.M());
        //     h_W2_truth_mass->Fill(W2_option1.M());
        //     h_chi2_truth->Fill(chi2_option1);
        // } else if (chi2_option2 <= chi2_option3) {
        //     // Option 2 is chosen
        //     h_W1_truth_mass->Fill(W1_option2.M());
        //     h_W2_truth_mass->Fill(W2_option2.M());
        //     h_chi2_truth->Fill(chi2_option2);
        // } else {
        //     // Option 2 is chosen
        //     h_W1_truth_mass->Fill(W1_option3.M());
        //     h_W2_truth_mass->Fill(W2_option3.M());
        //     h_chi2_truth->Fill(chi2_option3);
        // }

        // // ************************************ END OF EDITS ************************************




        // Fill the histograms
        // truth_W_e->Fill(truth_Wp_HS_e.at(0));

        // cutflow histograms
        cutFlowHist->SetBinContent(1, NEvents);
        cutFlowHist->SetBinContent(2, NokFlav);
        cutFlowHist->SetBinContent(3, NjetCut);
        cutFlowHist->SetBinContent(4, NleptonCut);
        cutFlowHist->SetBinContent(5, nFlavScore);
    }

    std::cout << "      " << std::endl;
    std::cout << "-------------------- Outputs --------------------" << std::endl;
    std::cout << "Number of events: " << NEvents << std::endl;
    std::cout << "Number of events w/ 4 jets: " << NjetCut << std::endl;
    std::cout << "Number of Leptons Cut: " << NleptonCut << std::endl;
    std::cout << "Number of events with 1 c-tagged, 1 s-tagged, and 2 light-tagged jets: " << nFlavScore << std::endl;
    std::cout << "      " << std::endl;
    // std::cout << "----------------- Truth Events ------------------" << std::endl;
    // std::cout << "Number of truth events: " << NEvents_truth << std::endl;
    // std::cout << "Number of truth events with 1 c-tagged, 1 s-tagged, and 2 light-tagged jets: " << nFlavScore_truth << std::endl;

    std::cout << "      " << std::endl;
    std::cout << "Let there be data!" << std::endl;


    // end of macro
}


