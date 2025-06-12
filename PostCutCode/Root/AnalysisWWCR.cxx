// Local includes
#include "PostCutCode/AnalysisWWCR.h"
#include "PostCutCode/TreeContainer.h"
#include <iomanip>
#include <TFile.h>
#include <fstream>
#include <TRandom.h>

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

    std::vector<std::string> cutFlowMap {"All Events", "Has 4 Jets"};

    // Get the histograms
    auto countingHist = m_histContainer->getCountingHist();
    auto truth_W_e =  m_histContainer->get1DHist("truth_W_e", 100, 0, 200);
    auto cutFlowHist = m_histContainer->get1DHist("cutFlowHist", cutFlowMap.size(), 0, 2, cutFlowMap);
    // Get the trees
    auto treeCont = std::make_shared<TreeContainer>();
 
    // Get max events to run on
    int nEntries = treeCont->getEventsToRun();

    // Connect branches to trees
    auto tree = treeCont->getTree();

    // calls the branch_name (var_type branch {tree, "branch_name"};)
    varMember<int> event_njet {tree, "event_njet"};
    varMember<ROOT::VecOps::RVec<float>>  truth_Wp_HS_e {tree, "truth_Wp_HS_e"};

    // reconstructed jets (for flavor tagger score)
    varMember<ROOT::VecOps::RVec<float>> recojet_isB {tree, "recojet_isB"};
    varMember<ROOT::VecOps::RVec<float>> recojet_isC {tree, "recojet_isC"};
    varMember<ROOT::VecOps::RVec<float>> recojet_isS {tree, "recojet_isS"};
    varMember<ROOT::VecOps::RVec<float>> recojet_isU {tree, "recojet_isU"};
    varMember<ROOT::VecOps::RVec<float>> recojet_isD {tree, "recojet_isD"};
    varMember<ROOT::VecOps::RVec<float>> recojet_isG {tree, "recojet_isG"};
    varMember<ROOT::VecOps::RVec<float>> recojet_isTAU {tree, "recojet_isTAU"};

    // Increment for CutFlow
    int NEvents = 0;
    int NjetCut = 0;
    int eventNum = 1;
    // int n_c = 0, n_s = 0, n_l = 0;
    // int n_c = 0;
    // int n_s = 0;
    // int n_l = 0;
    int nFlavScore = 0; 

    // Main event loop 
    for(int i = 0; i < nEntries; i++)
    {
        treeCont->getEntry(i);
        // Just to store how many events were run over
        countingHist->Fill(1);
        NEvents++;

        if(i % 10000 == 0) std::cout<<"Done i: "<<i<<" out of "<<nEntries<<std::endl;

        // Adds the requirement that an event has 4 jets
        if(event_njet() != 4) continue;
        NjetCut++;

        // Old section (copy into if statement)
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

        std::map<int,std::vector<float>> jetFlavScores;
        jetFlavScores[0] = j0_flav;
        jetFlavScores[1] = j1_flav;
        jetFlavScores[2] = j2_flav;
        jetFlavScores[3] = j3_flav;

        std::vector<std::string> flavLabels = {"b", "c", "s", "l", "g", "tau"};

        std::cout << "========== Event " << eventNum << " ==========\n";

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

        eventNum++;

        // Get the max flavor scores
        auto j0_MaxScoreIt = std::max_element(j0_flav.begin(), j0_flav.end());
        auto j1_MaxScoreIt = std::max_element(j1_flav.begin(), j1_flav.end());
        auto j2_MaxScoreIt = std::max_element(j2_flav.begin(), j2_flav.end());
        auto j3_MaxScoreIt = std::max_element(j3_flav.begin(), j3_flav.end());

        int j0_maxScoreIdx = std::distance(j0_flav.begin(), j0_MaxScoreIt);
        int j1_maxScoreIdx = std::distance(j1_flav.begin(), j1_MaxScoreIt);
        int j2_maxScoreIdx = std::distance(j2_flav.begin(), j2_MaxScoreIt);
        int j3_maxScoreIdx = std::distance(j3_flav.begin(), j3_MaxScoreIt);

        std::cout << "Jet 0: " << flavLabels[j0_maxScoreIdx] << " = " << *j0_MaxScoreIt << "\n";
        std::cout << "Jet 1: " << flavLabels[j1_maxScoreIdx] << " = " << *j1_MaxScoreIt << "\n";
        std::cout << "Jet 2: " << flavLabels[j2_maxScoreIdx] << " = " << *j2_MaxScoreIt << "\n";
        std::cout << "Jet 3: " << flavLabels[j3_maxScoreIdx] << " = " << *j3_MaxScoreIt << "\n";

        std::cout << "      " << std::endl;

        std::array<int, 4> maxScoreIdx {j0_maxScoreIdx,j1_maxScoreIdx,j2_maxScoreIdx,j3_maxScoreIdx};
        if (m_debug) {
            std::cout << "maxScoreIdx: ";
                for (float ScoreIdx : maxScoreIdx) { // change back to float if need be
                    std::cout << ScoreIdx << " ";
                }
                std::cout << std::endl;
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

        // if (m_debug) {
        //     std::cout << "n_c: " << n_c << ", n_s: " << n_s << ", n_l: " << n_l << std::endl;
        // }

        if (has_invalid_flavor) continue;

        if (!(n_c == 1 && n_s == 1 && n_l == 2)) continue;

        nFlavScore++; // Add this as a 3 bin in cutflow

        

        // Fill the histograms
        truth_W_e->Fill(truth_Wp_HS_e.at(0));

        cutFlowHist->SetBinContent(1, NEvents);
        cutFlowHist->SetBinContent(2, NjetCut);
    }

    std::cout << "      " << std::endl;
    std::cout << "-------------------- Outputs --------------------" << std::endl;
    std::cout << "Number of events: " << NEvents << std::endl;
    std::cout << "Number of events w/ 4 jets: " << NjetCut << std::endl;
    std::cout << "Number of events with 1 c-tag, 1 s-tag, and 2 light-tag jets: " << nFlavScore << std::endl;

    std::cout << "      " << std::endl;
    std::cout << "Let there be data!" << std::endl;

    // end of macro
}


