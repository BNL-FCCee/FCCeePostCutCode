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

// struct Reco_wPair {
//     TLorentzVector W1;
//     TLorentzVector W2;
//     double chi2;
// };

// Reco_wPair Chi2_wMass(const TLorentzVector& cTag_Jet,
//                       const TLorentzVector& lTag_Jet_0,
//                       const TLorentzVector& lTag_Jet_1,
//                       const TLorentzVector& lTag_Jet_2,
//                       double m_W_true = 80.379) {

//     // ************************* CHANGE // COMMENT OUT WHEN RUNNING WITH DIFFERENT WW DECAY OPTIONS *******************************
//     std::vector<std::pair<TLorentzVector, TLorentzVector>> pairings = {
//         {cTag_Jet + lTag_Jet_0, lTag_Jet_1 + lTag_Jet_2},
//         {cTag_Jet + lTag_Jet_1, lTag_Jet_0 + lTag_Jet_2},
//         {cTag_Jet + lTag_Jet_2, lTag_Jet_0 + lTag_Jet_1}

//         // {lTag_Jet_1 + lTag_Jet_2, lTag_Jet_1 + lTag_Jet_2},
//     };

//     Reco_wPair best_result;
//     best_result.chi2 = std::numeric_limits<double>::max();

//     for (const auto& [W1, W2] : pairings) {
//         double chi2 = (pow(W1.M() - m_W_true, 2) + pow(W2.M() - m_W_true, 2)) / (m_W_true);

//         if (chi2 < best_result.chi2) {
//             best_result = {W1, W2, chi2};
//         }
//     }

//     return best_result;
// }

void CosPhi_Angle(const TLorentzVector& jet1, const TLorentzVector& jet2, TH1F* hist) {
    double phi = jet1.Angle(jet2.Vect());
    double cos_phi= cos(phi);

    hist->Fill(cos_phi);
}

void AnalysisWWCR::run()
{ 

    std::vector<std::string> cutFlowMap {"All Events", "DecayCuts", "NoNaNFlavScore", "Has 4 Jets", "leptonCut", "d123Cut", "d34Cut", "Events have a c, s, and 2 l tagged jets"};

    // Get the histograms
    auto countingHist = m_histContainer->getCountingHist();
    auto cutFlowHist = m_histContainer->get1DHist("cutFlowHist", cutFlowMap.size(), 0, 8, cutFlowMap);

    // ******************* TRUTH HISTOGRAMS *******************
    auto h_W1_truth_mass = m_histContainer->get1DHist("h_W1_truth_mass", 300, 0, 150);
    auto h_W2_truth_mass = m_histContainer->get1DHist("h_W2_truth_mass", 300, 0, 150);

    auto h_W1_truth_p = m_histContainer->get1DHist("h_W1_truth_p", 300, 0, 150);
    auto h_W2_truth_p = m_histContainer->get1DHist("h_W2_truth_p", 300, 0, 150);

    auto h_W1_truth_e = m_histContainer->get1DHist("h_W1_truth_e", 300, 0, 150);
    auto h_W2_truth_e = m_histContainer->get1DHist("h_W2_truth_e", 300, 0, 150);

    auto h_cos_phi_truth_cs = m_histContainer->get1DHist("h_cos_phi_truth_cs", 150, -1, 1);
    auto h_cos_phi_truth_ud = m_histContainer->get1DHist("h_cos_phi_truth_ud", 150, -1, 1);

    // ******************* RECOJET HISTOGRAMS *******************
    auto h_W1_mass = m_histContainer->get1DHist("h_W1_mass", 300, 0, 150);
    auto h_W2_mass = m_histContainer->get1DHist("h_W2_mass", 300, 0, 150);
    auto h_chi2 = m_histContainer->get1DHist("h_chi2", 300, 0, 150);

    auto h_cJet_pt = m_histContainer->get1DHist("cJet_pt", 300, 0, 150);
    auto h_lJet0_pt = m_histContainer->get1DHist("lJet0_pt", 300, 0, 150);
    auto h_lJet1_pt = m_histContainer->get1DHist("lJet1_pt", 300, 0, 150);
    auto h_lJet2_pt = m_histContainer->get1DHist("lJet2_pt", 300, 0, 150);

    auto h_cJet_p = m_histContainer->get1DHist("cJet_p", 300, 0, 150);
    auto h_lJet0_p = m_histContainer->get1DHist("lJet0_p", 300, 0, 150);
    auto h_lJet1_p = m_histContainer->get1DHist("lJet1_p", 300, 0, 150);
    auto h_lJet2_p = m_histContainer->get1DHist("lJet2_p", 300, 0, 150);

    auto h_W1_p = m_histContainer->get1DHist("h_W1_p", 150, 0, 150);
    auto h_W2_p = m_histContainer->get1DHist("h_W2_p", 150, 0, 150);

    auto h_W1_e = m_histContainer->get1DHist("h_W1_e", 300, 0, 150);
    auto h_W2_e = m_histContainer->get1DHist("h_W2_e", 300, 0, 150);

    auto h_cos_phi_c_s = m_histContainer->get1DHist("h_cos_phi_c_s", 150, -1, 1);
    auto h_cos_phi_l1_l2 = m_histContainer->get1DHist("h_cos_phi_l1_l2", 150, -1, 1);
    // auto h_cos_phi_c_l1 = m_histContainer->get1DHist("h_cos_phi_c_l1", 300, -1, 1);  
    // auto h_cos_phi_c_l2 = m_histContainer->get1DHist("h_cos_phi_c_l2", 300, -1, 1); 
    // auto h_cos_phi_s_l1 = m_histContainer->get1DHist("h_cos_phi_s_l1", 300, -1, 1);  
    // auto h_cos_phi_s_l2 = m_histContainer->get1DHist("h_cos_phi_s_l2", 300, -1, 1);

    // auto h_W1_cos_phi_l1_l2 = m_histContainer->get1DHist("h_W1_cos_phi_l1_l2", 150, -1, 1); 
    // auto h_W2_cos_phi_l1_l2 = m_histContainer->get1DHist("h_W2_cos_phi_l1_l2", 150, -1, 1); 

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

    const int maxPrint = 25;
    int nPrinted = 0;
    bool do_debug = true; // Set to false to turn off all debug printing

    // *************************** SAME AS BELOW, LOOKING OVER THE MACRO STEP BY STEP ***************************
    for(int i = 0; i < nEntries; i++)
    {
        treeCont->getEntry(i);
        // Just to store how many events were run over
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

            nPrinted++; // Count only after printing all blocks for this event
        }

        eventNum++;

        QuarkCounter W1_quarks = countQuarks(Wm_jet_truth);
        QuarkCounter W2_quarks = countQuarks(Wp_jet_truth);

        bool W1_is_cs = ( (W1_quarks.n_c == 1 && W1_quarks.n_s == 1) );
        // bool W1_is_ud = ( (W1_quarks.n_u == 1 && W1_quarks.n_d == 1) );

        // bool W2_is_cs = ( (W2_quarks.n_c == 1 && W2_quarks.n_s == 1) );
        bool W2_is_ud = ( (W2_quarks.n_u == 1 && W2_quarks.n_d == 1) );

        if (! ( 
            (W1_is_cs && W2_is_ud) 
            // (W1_is_ud && W2_is_ud)
            // (W1_is_cs && W2_is_cs)
        )) continue;

        NdecayCuts++;
    
        // Adds the requirement that an event has 4 jets
        if(event_njet() != 4) continue;
        NjetCut++;

        // ****************************** CALCULATIONS USING TRUTH VARIABLES ******************************
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

            CosPhi_Angle(W1_jet1, W1_jet2, h_cos_phi_truth_cs);
            CosPhi_Angle(W2_jet1, W2_jet2, h_cos_phi_truth_ud);
        }

        // *************************** END OF CALCULATIONS USING TRUTH VARIABLES ***************************

        bool flage_toss = false;
        if (recojet_isB.size() == 0){
            flage_toss = true;
        }
        if (flage_toss) continue;

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
        if (event_nmu() > 2.) continue;
        if (event_nel() > 2.) continue;
        if (muons_p() > 20.) continue;
        if (elecrons_p() > 20.) continue;

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

        // if (cJet == -1 || lJets.size() != 3) {
        //     std::cerr << "Error: Incorrect number of flavored jets for pairing!" << std::endl;
        //     continue;

        //     // continue;
        // }

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

        // calculate the mass 
        TLorentzVector cTag_Jet, lTag_Jet_0, lTag_Jet_1, lTag_Jet_2;

        cTag_Jet.SetPxPyPzE(jet_px.at(cJet), jet_py.at(cJet), jet_pz.at(cJet), jet_e.at(cJet));
        lTag_Jet_0.SetPxPyPzE(jet_px.at(lJets[0]), jet_py.at(lJets[0]), jet_pz.at(lJets[0]), jet_e.at(lJets[0]));

        lTag_Jet_1.SetPxPyPzE(jet_px.at(lJets[1]), jet_py.at(lJets[1]), jet_pz.at(lJets[1]), jet_e.at(lJets[1]));
        lTag_Jet_2.SetPxPyPzE(jet_px.at(lJets[2]), jet_py.at(lJets[2]), jet_pz.at(lJets[2]), jet_e.at(lJets[2]));
        

        // ************************** CALCULATIONS FOR MASS AND SUCH BEGIN HERE **************************
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

        if (chi2_option1 <= chi2_option2 && chi2_option1 <= chi2_option3) {
            h_W1_mass->Fill(W1_option1.M());
            h_W2_mass->Fill(W2_option1.M());
            h_chi2->Fill(chi2_option1);

            h_W1_p->Fill(W1_option1.P());
            h_W2_p->Fill(W2_option1.P());
            
            h_W1_e->Fill(W1_option1.E());
            h_W2_e->Fill(W2_option1.E());

            h_cJet_pt->Fill(cTag_Jet.Pt());
            h_lJet0_pt->Fill(lTag_Jet_0.Pt());
            h_lJet1_pt->Fill(lTag_Jet_1.Pt());
            h_lJet2_pt->Fill(lTag_Jet_2.Pt());

            h_cJet_p->Fill(cTag_Jet.P());
            h_lJet0_p->Fill(lTag_Jet_0.P());
            h_lJet1_p->Fill(lTag_Jet_1.P());
            h_lJet2_p->Fill(lTag_Jet_2.P());

            CosPhi_Angle(cTag_Jet, lTag_Jet_0, h_cos_phi_c_s);
            CosPhi_Angle(lTag_Jet_1, lTag_Jet_2, h_cos_phi_l1_l2);

        } else if (chi2_option2 <= chi2_option1 && chi2_option2 <= chi2_option3) {
            h_W1_mass->Fill(W1_option2.M());
            h_W2_mass->Fill(W2_option2.M());
            h_chi2->Fill(chi2_option2);

            h_W1_p->Fill(W1_option2.P());
            h_W2_p->Fill(W2_option2.P());

            h_W1_e->Fill(W1_option2.E());
            h_W2_e->Fill(W2_option2.E());

            h_cJet_pt->Fill(cTag_Jet.Pt());
            h_lJet0_pt->Fill(lTag_Jet_0.Pt());
            h_lJet1_pt->Fill(lTag_Jet_1.Pt());
            h_lJet2_pt->Fill(lTag_Jet_2.Pt());

            h_cJet_p->Fill(cTag_Jet.P());
            h_lJet0_p->Fill(lTag_Jet_0.P());
            h_lJet1_p->Fill(lTag_Jet_1.P());
            h_lJet2_p->Fill(lTag_Jet_2.P());

            CosPhi_Angle(cTag_Jet, lTag_Jet_0, h_cos_phi_c_s);
            CosPhi_Angle(lTag_Jet_1, lTag_Jet_2, h_cos_phi_l1_l2);

        } else {
            h_W1_mass->Fill(W1_option3.M());
            h_W2_mass->Fill(W2_option3.M());
            h_chi2->Fill(chi2_option3);

            h_W1_p->Fill(W1_option3.P());
            h_W2_p->Fill(W2_option3.P());

            h_W1_e->Fill(W1_option3.E());
            h_W2_e->Fill(W2_option3.E());

            h_cJet_pt->Fill(cTag_Jet.Pt());
            h_lJet0_pt->Fill(lTag_Jet_0.Pt());
            h_lJet1_pt->Fill(lTag_Jet_1.Pt());
            h_lJet2_pt->Fill(lTag_Jet_2.Pt());

            h_cJet_p->Fill(cTag_Jet.P());
            h_lJet0_p->Fill(lTag_Jet_0.P());
            h_lJet1_p->Fill(lTag_Jet_1.P());
            h_lJet2_p->Fill(lTag_Jet_2.P());

            CosPhi_Angle(cTag_Jet, lTag_Jet_0, h_cos_phi_c_s);
            CosPhi_Angle(lTag_Jet_1, lTag_Jet_2, h_cos_phi_l1_l2);
        }

        // Reco_wPair result = Chi2_wMass(cTag_Jet, lTag_Jet_0, lTag_Jet_1, lTag_Jet_2);
        // Reco_wPair result = Chi2_wMass(lTag_Jet_1, lTag_Jet_1, lTag_Jet_1, lTag_Jet_2);

        // h_W1_mass->Fill(result.W1.M());
        // h_W2_mass->Fill(result.W2.M());
        // h_chi2->Fill(result.chi2);

        // h_W1_p->Fill(W1_option1.P());
        // h_W2_p->Fill(W2_option1.P());

        // Calculates the transverse and total momentum for the jets
        // h_cJet_pt->Fill(cTag_Jet.Pt());
        // h_lJet0_pt->Fill(lTag_Jet_0.Pt());
        // h_lJet1_pt->Fill(lTag_Jet_1.Pt());
        // h_lJet2_pt->Fill(lTag_Jet_2.Pt());

        // h_cJet_p->Fill(cTag_Jet.P());
        // h_lJet0_p->Fill(lTag_Jet_0.P());
        // h_lJet1_p->Fill(lTag_Jet_1.P());
        // h_lJet2_p->Fill(lTag_Jet_2.P());

        // ************************* CHANGE // COMMENT OUT WHEN RUNNING WITH DIFFERENT WW DECAY OPTIONS *******************************
        // CosPhi_Angle(cTag_Jet, lTag_Jet_0, h_cos_phi_c_s);
        // CosPhi_Angle(lTag_Jet_1, lTag_Jet_2, h_cos_phi_l1_l2);
        // CosPhi_Angle(cTag_Jet, lTag_Jet_1, h_cos_phi_c_l1);
        // CosPhi_Angle(cTag_Jet, lTag_Jet_2, h_cos_phi_c_l2);
        // CosPhi_Angle(lTag_Jet_0, lTag_Jet_1, h_cos_phi_s_l1);
        // CosPhi_Angle(lTag_Jet_0, lTag_Jet_2, h_cos_phi_s_l2);

        // CosPhi_Angle(lTag_Jet_1, lTag_Jet_2, h_W1_cos_phi_l1_l2);
        // CosPhi_Angle(lTag_Jet_1, lTag_Jet_2, h_W2_cos_phi_l1_l2);

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
    std::cout << "Let there be data ;)" << std::endl;


    // end of macro
}


