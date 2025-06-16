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

    std::vector<std::string> cutFlowMap {"All Events", "Has 4 Jets", "Events have a c, s, and 2 l tagged jets"};

    // Get the histograms
    auto countingHist = m_histContainer->getCountingHist();
    // auto truth_W_e =  m_histContainer->get1DHist("truth_W_e", 100, 0, 200);
    auto cutFlowHist = m_histContainer->get1DHist("cutFlowHist", cutFlowMap.size(), 0, 3, cutFlowMap);
    // Get the trees
    auto treeCont = std::make_shared<TreeContainer>();
 
    // Get max events to run on
    int nEntries = treeCont->getEventsToRun();

    // Connect branches to trees
    auto tree = treeCont->getTree();

    // calls the branch_name (var_type branch {tree, "branch_name"};)
    varMember<int> event_njet {tree, "event_njet"};
    // varMember<ROOT::VecOps::RVec<float>>  truth_Wp_HS_e {tree, "truth_Wp_HS_e"};

    // reconstructed jets (for flavor tagger score)
    varMember<ROOT::VecOps::RVec<float>> recojet_isB {tree, "recojet_isB"};
    varMember<ROOT::VecOps::RVec<float>> recojet_isC {tree, "recojet_isC"};
    varMember<ROOT::VecOps::RVec<float>> recojet_isS {tree, "recojet_isS"};
    varMember<ROOT::VecOps::RVec<float>> recojet_isU {tree, "recojet_isU"};
    varMember<ROOT::VecOps::RVec<float>> recojet_isD {tree, "recojet_isD"};
    varMember<ROOT::VecOps::RVec<float>> recojet_isG {tree, "recojet_isG"};
    varMember<ROOT::VecOps::RVec<float>> recojet_isTAU {tree, "recojet_isTAU"};

    // Mass of W boson
    // float W_mass = 80.3692;

    // momentum
    varMember<ROOT::VecOps::RVec<float>> jet_px {tree, "recojet_px"};
    varMember<ROOT::VecOps::RVec<float>> jet_py {tree, "recojet_py"};
    varMember<ROOT::VecOps::RVec<float>> jet_pz {tree, "recojet_pz"};
    varMember<ROOT::VecOps::RVec<float>> jet_e {tree, "recojet_e"};

    // TFile* outFile_WWCR = new TFile("test.root", "Recreate");

    TH1F *h_W1_mass = new TH1F("h_W1_mass", "W1 Mass", 50, 0, 150);
    TH1F *h_W2_mass = new TH1F("h_W2_mass", "W2 Mass", 50, 0, 150);

    // Increment for CutFlow
    int NEvents = 0;
    int NjetCut = 0;
    int eventNum = 1;
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

        std::cout << "      " << std::endl;  

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

        // Optional: print the pairs
        std::cout << "W1 pair: Jet " << W1_pair.first << " (c), Jet " << W1_pair.second << " (s)" << std::endl;
        std::cout << "W2 pair: Jet " << W2_pair.first << " (l), Jet " << W2_pair.second << " (l)" << std::endl;
        std::cout << "      " << std::endl;

        // calculate the mass 
        TLorentzVector cTag_Jet, sTag_Jet, lTag_Jet_0, lTag_Jet_1;

        cTag_Jet.SetPxPyPzE(jet_px.at(cJet), jet_py.at(cJet), jet_py.at(cJet), jet_e.at(cJet));
        sTag_Jet.SetPxPyPzE(jet_px.at(sJet), jet_py.at(sJet), jet_py.at(sJet), jet_e.at(sJet));

        lTag_Jet_0.SetPxPyPzE(jet_px.at(lJets[0]), jet_py.at(lJets[0]), jet_py.at(lJets[0]), jet_e.at(lJets[0]));
        lTag_Jet_1.SetPxPyPzE(jet_px.at(lJets[1]), jet_py.at(lJets[1]), jet_py.at(lJets[1]), jet_e.at(lJets[1]));

        TLorentzVector W1 = cTag_Jet + sTag_Jet;
        TLorentzVector W2 = lTag_Jet_0 + lTag_Jet_1;

        h_W1_mass->Fill(W1.M());
        h_W2_mass->Fill(W2.M());

        std::cout << "W1 mass: " << W1.M() << " GeV" << std::endl;
        std::cout << "W2 mass: " << W2.M() << " GeV" << std::endl;

        std::cout << "  " << std::endl;



        // Fill the histograms
        // truth_W_e->Fill(truth_Wp_HS_e.at(0));

        cutFlowHist->SetBinContent(1, NEvents);
        cutFlowHist->SetBinContent(2, NjetCut);
        cutFlowHist->SetBinContent(3, nFlavScore);
    }

    std::cout << "      " << std::endl;
    std::cout << "-------------------- Outputs --------------------" << std::endl;
    std::cout << "Number of events: " << NEvents << std::endl;
    std::cout << "Number of events w/ 4 jets: " << NjetCut << std::endl;
    std::cout << "Number of events with 1 c-tagged, 1 s-tagged, and 2 light-tagged jets: " << nFlavScore << std::endl;

    std::cout << "      " << std::endl;
    std::cout << "Let there be data!" << std::endl;

    // // Plot the histograms and save them to "test.root"
    // TCanvas* c_Wmass = new TCanvas("c_Wmass", "W Boson Masses", 1000, 500);
    // c_Wmass->Divide(1, 1);

    // // Draw W1 mass
    // c_Wmass->cd(1);
    // h_W1_mass->SetLineColor(kBlue);
    // h_W1_mass->Draw();

    // h_W2_mass->SetLineColor(kRed);
    // h_W2_mass->Draw("same");

    // outFile_WWCR->cd();

    // c_Wmass->Write;

    // outFile_WWCR->Write();
    // outFile_WWCR->Close();



    // end of macro
}


