// Local includes
#include "PostCutCode/AnalysisZHAllHad.h"
#include "PostCutCode/TreeContainer.h"
#include <iomanip>
#include <TFile.h>
#include <fstream>
#include <set>
#include "TLorentzVector.h"
#include <TRandom.h>

// DEFINE sel_visMtheta, sel_d123, sel_d34! WHICH SEL TO PICK! CAREFUL!
using namespace std;

// Base class 
// std::vector<unsigned int> AnalysisZHAllHad::missing_pair(std::vector<int> BLA){

    
// }

AnalysisZHAllHad::AnalysisZHAllHad():
    AnalysisBase(),
    m_debug(false)
{}
AnalysisZHAllHad::~AnalysisZHAllHad()
{}


void AnalysisZHAllHad::run()
{       
    gRandom->SetSeed(42); 
    float Z_mass = 91.1876;
    float H_mass = 125.11;
    float W_mass = 80.377;
    std::map<int, std::string> char_jet_flav = {
        {1, "q"},
        {2, "q"},
        {3, "s"},
        {4, "c"},
        {5, "b"}
    };
    std::vector<std::string> flavourJets {"b", "c", "s", "g", "q"};
    std::vector<std::string> flavours{"B","C","S","D","U","G","TAU"};
    std::vector<std::string> flavourCategory {"B", "C", "S","Q", "G","TAU"};
    std::vector<std::string> fitCategory {"LowHss","MidHss","HiHss","LowbbHbb","LowccHbb","LowssHbb","LowqqHbb","LowbbHcc","LowccHcc","LowssHcc","LowqqHcc","LowbbHgg","LowccHgg","LowssHgg","LowqqHgg","MidbbHbb","MidccHbb","MidssHbb","MidqqHbb","MidbbHcc","MidccHcc","MidssHcc","MidqqHcc","MidbbHgg","MidccHgg","MidssHgg","MidqqHgg","HibbHbb","HiccHbb","HissHbb","HiqqHbb","HibbHcc","HiccHcc","HissHcc","HiqqHcc","HibbHgg","HiccHgg","HissHgg","HiqqHgg","Incl"};
    std::vector<std::string> cutFlowMap {"NoCut","njet=4","leptonCut","KineCut", "d123Cut", "d34Cut","Pairing","jjMassCut","ZHmassCut","YieldsFit" };
    std::vector<std::string> fitCategorySimple  {"LowHbb","LowHcc","LowHss","LowHgg","MidHbb","MidHcc","MidHss","MidHgg","HiHbb","HiHcc", "HiHss","HiHgg"};

    
    // Get the histograms
    auto histo_DetVars = m_histContainer->histo_DetVarsScoreSmear(flavourJets);

    auto scoreMapHist = m_histContainer->get1DHist("scoreMap_1D", flavourCategory.size(), 0, flavourCategory.size(), flavourCategory);
    auto scoreMapFitCatHist = m_histContainer->get1DHist("scoreMapFitCategory_1D", fitCategorySimple.size(), 0, fitCategorySimple.size(), fitCategorySimple);
    auto CountsMapHist = m_histContainer->get1DHist("CountsMap_1D", flavourCategory.size(), 0, flavourCategory.size(), flavourCategory);
    auto CountsFitCatHist = m_histContainer->get1DHist("CountsFitCategory_1D", fitCategorySimple.size(), 0, fitCategorySimple.size(), fitCategorySimple);
    auto obsHist = m_histContainer->getObsHistinFitCategory(fitCategory, 250, 0., 250.,250, 0., 250);
//     auto obsHist = m_histContainer->getObsHistinFitCategory(fitCategory, 90,90.,180.,70,50.,120.);
//     auto obsHist = m_histContainer->getObsHistinFitCategory(fitCategory, 250, 0., 250.);
    auto countingHist = m_histContainer->getCountingHist();
    auto cutFlowHist =  m_histContainer->get1DHist("CutFlow", cutFlowMap.size(), 0, 12,cutFlowMap);

    auto Bscore =  m_histContainer->get1DHist("Higgs_Bscore", 100, 0, 2);
    auto Cscore =  m_histContainer->get1DHist("Higgs_Cscore", 100, 0, 2);
    auto Sscore =  m_histContainer->get1DHist("Higgs_Sscore", 100, 0, 2);
    auto Qscore =  m_histContainer->get1DHist("Higgs_Qscore", 100, 0, 2);
    auto Gscore =  m_histContainer->get1DHist("Higgs_Gscore", 100, 0, 2);
    
    auto Z_Bscore =  m_histContainer->get1DHist("Z_Bscore", 100, 0, 2);
    auto Z_Cscore =  m_histContainer->get1DHist("Z_Cscore", 100, 0, 2);
    auto Z_Sscore =  m_histContainer->get1DHist("Z_Sscore", 100, 0, 2);
    auto Z_Qscore =  m_histContainer->get1DHist("Z_Qscore", 100, 0, 2);
    auto Z_Gscore =  m_histContainer->get1DHist("Z_Gscore", 100, 0, 2);
//     auto TAUscore =  m_histContainer->get1DHist("TAUscore", 100, 0, 2);
        // auto Uscore =  m_histContainer->get1DHist("Uscore", 100, 0, 2);


    auto mZ =  m_histContainer->get1DHist("mZ",250, 0., 250.);
    auto mH =  m_histContainer->get1DHist("mH", 250, 0., 250.);
    auto mH_corr = m_histContainer->get1DHist("mH_corr", 250, 0., 250.);

    auto ZZ_cut =  m_histContainer->get1DHist("ZZ_cut",250, 0., 250.);
    auto WW_cut =  m_histContainer->get1DHist("WW_cut", 250, 0., 250.);

    auto h_vis_M =  m_histContainer->get1DHist("vis_M",250, 0., 250.);
    auto h_vis_E =  m_histContainer->get1DHist("vis_e", 250, 0., 250.);
    auto h_vis_theta =  m_histContainer->get1DHist("vis_theta", 80, 0., 4.);

    auto h_muons_p =  m_histContainer->get1DHist("muons_p", 125, 0., 250.);
    auto h_electrons_p =  m_histContainer->get1DHist("electrons_p", 125, 0., 250.);
    auto h_event_nmu =  m_histContainer->get1DHist("event_nmu", 10, 0., 10.);
    auto h_event_nel =  m_histContainer->get1DHist("event_nel", 10, 0., 10.);

    auto h_d_12 =  m_histContainer->get1DHist("d_12", 1000,14000, 60000.);
    auto h_d_23 =  m_histContainer->get1DHist("d_23", 1000,200, 20000.);
    auto h_d_34 =  m_histContainer->get1DHist("d_34", 1000,50, 6500.);

        // Extra hist for easy fill up
    auto Incl_obsHist = obsHist["Incl"];
    // H->bb
    auto Low_bbZ_Hbb_obsHist = obsHist["LowbbHbb"];
    auto Low_ccZ_Hbb_obsHist = obsHist["LowccHbb"];
    auto Low_ssZ_Hbb_obsHist = obsHist["LowssHbb"];
    auto Low_qqZ_Hbb_obsHist = obsHist["LowqqHbb"];

    auto Mid_bbZ_Hbb_obsHist = obsHist["MidbbHbb"];
    auto Mid_ccZ_Hbb_obsHist = obsHist["MidccHbb"];
    auto Mid_ssZ_Hbb_obsHist = obsHist["MidssHbb"];
    auto Mid_qqZ_Hbb_obsHist = obsHist["MidqqHbb"];

    auto Hi_bbZ_Hbb_obsHist = obsHist["HibbHbb"];
    auto Hi_ccZ_Hbb_obsHist = obsHist["HiccHbb"];
    auto Hi_ssZ_Hbb_obsHist = obsHist["HissHbb"];
    auto Hi_qqZ_Hbb_obsHist = obsHist["HiqqHbb"];
    // H->cc
    auto Low_bbZ_Hcc_obsHist = obsHist["LowbbHcc"];
    auto Low_ccZ_Hcc_obsHist = obsHist["LowccHcc"];
    auto Low_ssZ_Hcc_obsHist = obsHist["LowssHcc"];
    auto Low_qqZ_Hcc_obsHist = obsHist["LowqqHcc"];

    auto Mid_bbZ_Hcc_obsHist = obsHist["MidbbHcc"];
    auto Mid_ccZ_Hcc_obsHist = obsHist["MidccHcc"];
    auto Mid_ssZ_Hcc_obsHist = obsHist["MidssHcc"];
    auto Mid_qqZ_Hcc_obsHist = obsHist["MidqqHcc"];

    auto Hi_bbZ_Hcc_obsHist = obsHist["HibbHcc"];
    auto Hi_ccZ_Hcc_obsHist = obsHist["HiccHcc"];
    auto Hi_ssZ_Hcc_obsHist = obsHist["HissHcc"];
    auto Hi_qqZ_Hcc_obsHist = obsHist["HiqqHcc"];
    //H->ss
    auto LowS_obsHist = obsHist["LowHss"];
    auto MidS_obsHist = obsHist["MidHss"];
    auto HiS_obsHist = obsHist["HiHss"];

//     // H->qq
//     auto Low_bbZ_Hqq_obsHist = obsHist["LowbbHqq"];
//     auto Low_ccZ_Hqq_obsHist = obsHist["LowccHqq"];
//     auto Low_ssZ_Hqq_obsHist = obsHist["LowssHqq"];
//     auto Low_qqZ_Hqq_obsHist = obsHist["LowqqHqq"];

//     auto Mid_bbZ_Hqq_obsHist = obsHist["MidbbHqq"];
//     auto Mid_ccZ_Hqq_obsHist = obsHist["MidccHqq"];
//     auto Mid_ssZ_Hqq_obsHist = obsHist["MidssHqq"];
//     auto Mid_qqZ_Hqq_obsHist = obsHist["MidqqHqq"];

//     auto Hi_bbZ_Hqq_obsHist = obsHist["HibbHqq"];
//     auto Hi_ccZ_Hqq_obsHist = obsHist["HiccHqq"];
//     auto Hi_ssZ_Hqq_obsHist = obsHist["HissHqq"];
//     auto Hi_qqZ_Hqq_obsHist = obsHist["HiqqHqq"];

    // H->gg
    auto Low_bbZ_Hgg_obsHist = obsHist["LowbbHgg"];
    auto Low_ccZ_Hgg_obsHist = obsHist["LowccHgg"];
    auto Low_ssZ_Hgg_obsHist = obsHist["LowssHgg"];
    auto Low_qqZ_Hgg_obsHist = obsHist["LowqqHgg"];

    auto Mid_bbZ_Hgg_obsHist = obsHist["MidbbHgg"];
    auto Mid_ccZ_Hgg_obsHist = obsHist["MidccHgg"];
    auto Mid_ssZ_Hgg_obsHist = obsHist["MidssHgg"];
    auto Mid_qqZ_Hgg_obsHist = obsHist["MidqqHgg"];

    auto Hi_bbZ_Hgg_obsHist = obsHist["HibbHgg"];
    auto Hi_ccZ_Hgg_obsHist = obsHist["HiccHgg"];
    auto Hi_ssZ_Hgg_obsHist = obsHist["HissHgg"];
    auto Hi_qqZ_Hgg_obsHist = obsHist["HiqqHgg"];
    

    // Get the trees
    auto treeCont = std::make_shared<TreeContainer>();
    std::cout<<"sampleName: "<<MDC::GetInstance()->getSampleName()<<" events: "<<treeCont->getEntries()<<std::endl;

    // Get max events to run on
    int nEntries = treeCont->getEventsToRun();
  

    // Connect branches to trees
    auto tree = treeCont->getTree();
  
    // assuming you have vectors as input (should also save option to run w/o vectors)
  
    // flavor scores
    varMember<ROOT::VecOps::RVec<float>> recojet_isB {tree, "recojet_isB"};
    varMember<ROOT::VecOps::RVec<float>> recojet_isC {tree, "recojet_isC"};
    varMember<ROOT::VecOps::RVec<float>> recojet_isG {tree, "recojet_isG"};
    varMember<ROOT::VecOps::RVec<float>> recojet_isU {tree, "recojet_isU"};
    varMember<ROOT::VecOps::RVec<float>> recojet_isTAU {tree, "recojet_isTAU"};
    varMember<ROOT::VecOps::RVec<float>> recojet_isD {tree, "recojet_isD"};
    varMember<ROOT::VecOps::RVec<float>> recojet_isS {tree, "recojet_isS"};
    //add truth jet flav
    varMember<ROOT::VecOps::RVec<int>> truth_flav {tree, "jets_truth"};
    
  
    // corrected momentum
    varMember<ROOT::VecOps::RVec<float>> jet_px {tree, "jet_px_corr"};
    varMember<ROOT::VecOps::RVec<float>> jet_py {tree, "jet_py_corr"}; 
    varMember<ROOT::VecOps::RVec<float>> jet_pz {tree, "jet_pz_corr"};
    varMember<ROOT::VecOps::RVec<float>> jet_e {tree, "jet_e_corr"}; 

    varMember<int> event_njet {tree, "event_njet"};
    varMember<ulong> event_nmu {tree, "event_nmu"};
    varMember<ulong> event_nel {tree, "event_nel"};
    varMember<float> muons_p {tree, "muons_p"};
    varMember<float> electrons_p {tree, "electrons_p"};
    varMember<double> vis_E {tree, "vis_E"};
    varMember<double> vis_M {tree, "vis_M"};
    varMember<double> vis_theta {tree, "vis_theta"};
    varMember<float> d_12 {tree, "d_12"};
    varMember<float> d_23 {tree, "d_23"};
    varMember<float> d_34 {tree, "d_34"};

    // varMember<double> costhetainv {tree, "costhetainv"};
    int NjetCut = 0;
    int NleptonCut = 0;
    int NkineCut = 0;
    int NdCutd123 = 0;
    int NdCutd34 = 0;
    int NafterPairing = 0;
    int NjjMassCut = 0;
    int NafterSel = 0;
    int Nfit = 0;
    int NEventsInt = 0;
    int BlikeEvents = 0;
    int ClikeEvents = 0;
    int SlikeEvents = 0;
    int GlikeEvents = 0;
    int QlikeEvents = 0;
//     int TAUlikeEvents = 0;

    std::array<int, 3> BlikeEvents_cat {0, 0, 0};
    std::array<int, 3> ClikeEvents_cat {0, 0, 0};
    std::array<int, 3> SlikeEvents_cat {0, 0, 0};
    std::array<int, 3> GlikeEvents_cat {0, 0, 0};
    std::array<int, 3> QlikeEvents_cat {0, 0, 0};
//     std::array<int, 3> TAUlikeEvents_cat {0, 0, 0};

//     int debug = 1;
  
    // Loop over the trees here
    for(int i = 0; i < nEntries; i++)
    {

        treeCont->getEntry(i);
        // Just to store how many events were run over
        countingHist->Fill(1);
        NEventsInt++;
        int H_flav = -1;
        int Z_flav = -1;

        if(i % 1000 == 0) std::cout<<"Done i: "<<i<<" out of "<<nEntries<<std::endl;
        // Add some basic selection        
        if(event_njet() != 4) continue; //Require EXACTLY 4 jet! This SHOULD be the case!
        NjetCut++;
        //numebr of electrons and muons cut! 2!
        // reject events with leptons with > 20 GeV leptons
        if(event_nmu() > 2.) continue;
        if(event_nel() > 2.) continue;
        if(muons_p() >= 20.) continue;
        if(electrons_p() >= 20.) continue;
        NleptonCut++;

        //Add this  to the options, as in which selection to choose
        if (vis_M() <= 150.) continue;
        if (vis_E() <= 150.) continue;
        if (vis_theta()<= 0.15 || vis_theta()>=3.0 ) continue; //typo! rerun!
        NkineCut++;
        // change the cuts! update! 
        if (d_12()<=15000. || d_12()>=58000.) continue;//missing!
        if ((d_23()<=400.) || (d_23()>=18000.))continue;
        NdCutd123++;
        if ((d_34()<=100.) || (d_34()>=6000.))continue;
        NdCutd34++;
        // ALL GOOD UP TO HERE!!!!!!
        // Step 1: check efficiency ~ 80%
        // increment counter of events that pass the cuts!
       //make lorentz vectors
//        std::cout<<"Passed intial selection!" << std::endl; 
       std::vector<TLorentzVector> LVjets;    
       std::vector<float> NEWrecojet_isB;
       std::vector<float> NEWrecojet_isC;
       std::vector<float> NEWrecojet_isS;
       for (size_t lv = 0; lv < 4; ++lv) {  
            TLorentzVector LVjet(jet_px.at(lv), jet_py.at(lv), jet_pz.at(lv), jet_e.at(lv));
            LVjets.push_back(LVjet);
            std::string jet_flav = char_jet_flav[abs(truth_flav.at(lv))];
//             std::cout<<"Jet: "<< lv <<std::endl; 
            // if (m_debug) std::cout<<"Truth Flav: "<< abs(truth_flav.at(lv)) <<std::endl;
            std::vector<int> set = {1, 2, 3, 4, 5};
            if ((std::find(set.begin(), set.end(), abs(truth_flav.at(lv))) == set.end())){
                // if (m_debug) std::cout<<"EMPTY truth_flav.at(lv): "<< truth_flav.at(lv) <<std::endl;
                NEWrecojet_isB.push_back(recojet_isB.at(lv));
                NEWrecojet_isC.push_back(recojet_isC.at(lv));
                NEWrecojet_isS.push_back(recojet_isS.at(lv));
            }
            else{
                double b_score;
                double c_score;
                double s_score;
                histo_DetVars[jet_flav]->GetRandom3(b_score, c_score,s_score);
                NEWrecojet_isB.push_back(b_score);
                NEWrecojet_isC.push_back(c_score);
                NEWrecojet_isS.push_back(s_score);
            } 
        }
      
       // Flav scores of each jet
       std::array<float,7> j0_flav {NEWrecojet_isB.at(0), NEWrecojet_isC.at(0), NEWrecojet_isS.at(0), recojet_isD.at(0), recojet_isU.at(0), recojet_isG.at(0), recojet_isTAU.at(0)};
       std::array<float,7> j1_flav {NEWrecojet_isB.at(1), NEWrecojet_isC.at(1), NEWrecojet_isS.at(1), recojet_isD.at(1), recojet_isU.at(1), recojet_isG.at(1), recojet_isTAU.at(1)};
       std::array<float,7> j2_flav {NEWrecojet_isB.at(2), NEWrecojet_isC.at(2), NEWrecojet_isS.at(2), recojet_isD.at(2), recojet_isU.at(2), recojet_isG.at(2), recojet_isTAU.at(2)};
       std::array<float,7> j3_flav {NEWrecojet_isB.at(3), NEWrecojet_isC.at(3), NEWrecojet_isS.at(3), recojet_isD.at(3), recojet_isU.at(3), recojet_isG.at(3), recojet_isTAU.at(3)};

        // - look for max score of jet 
       std::array<float,7>::iterator j0_maxScore;
       std::array<float,7>::iterator j1_maxScore;
       std::array<float,7>::iterator j2_maxScore; 
       std::array<float,7>::iterator j3_maxScore; 

       j0_maxScore = std::max_element(j0_flav.begin(), j0_flav.end());
       j1_maxScore = std::max_element(j1_flav.begin(), j1_flav.end());
       j2_maxScore = std::max_element(j2_flav.begin(), j2_flav.end());
       j3_maxScore = std::max_element(j3_flav.begin(), j3_flav.end());
       //Can simplify with j<>_maxScore
       int j0_maxScoreIdx = std::distance(j0_flav.begin(),std::max_element(j0_flav.begin(), j0_flav.end()));
       int j1_maxScoreIdx = std::distance(j1_flav.begin(),std::max_element(j1_flav.begin(), j1_flav.end()));
       int j2_maxScoreIdx = std::distance(j2_flav.begin(),std::max_element(j2_flav.begin(), j2_flav.end()));
       int j3_maxScoreIdx = std::distance(j3_flav.begin(),std::max_element(j3_flav.begin(), j3_flav.end()));

       std::map<int, std::array<float,7>> flavMap;
       flavMap[0]=j0_flav;
       flavMap[1]=j1_flav;
       flavMap[2]=j2_flav;
       flavMap[3]=j3_flav;

       std::array<int, 4> maxScoreIdx {j0_maxScoreIdx,j1_maxScoreIdx,j2_maxScoreIdx,j3_maxScoreIdx};

       if (m_debug) {
            std::cout<<"j0_flav: "<<" Max score: "<< *j0_maxScore<< " maxScoreIdx: "<< maxScoreIdx[0]<<std::endl;
            std::cout<<"j1_flav: "<<" Max score: "<< *j1_maxScore<< " maxScoreIdx: "<< maxScoreIdx[1]<<std::endl;
            std::cout<<"j2_flav: "<<" Max score: "<< *j2_maxScore<< " maxScoreIdx: "<< maxScoreIdx[2]<<std::endl;
            std::cout<<"j3_flav: "<<" Max score: "<< *j3_maxScore<< " maxScoreIdx: "<< maxScoreIdx[3]<<std::endl;
       }
       if (m_debug) {
            std::cout << "maxScoreIdx: ";
                for (float ScoreIdx : maxScoreIdx) {
                    std::cout << ScoreIdx << " ";
                }
                std::cout << std::endl;
       }     
       if (m_debug) {
            std::cout << "Jet score for j0: ";
                for (float j0_f : j0_flav) {
                    std::cout << j0_f << " ";
                }
                std::cout << std::endl;

            std::cout << "Jet score for j1: ";
                for (float j1_f : j1_flav) {
                    std::cout << j1_f << " ";
                }
                std::cout << std::endl;
            std::cout << "Jet score for j2: ";
                for (float j2_f : j2_flav) {
                    std::cout << j2_f << " ";
                }
                std::cout << std::endl;
            std::cout << "Jet score for j3: ";
                for (float j3_f : j3_flav) {
                    std::cout << j3_f << " ";
                }
                std::cout << std::endl;}
              //Find the pairs

        //Check how many pairs you get from the simple method  
        //Pair the jets, make sth very simple 
        std::map<int, std::vector<int>> tmp_jet_pair;
        std::map<int, std::vector<int>> jet_pair;
        std::map<int, int> jet_pair_flav;
        std::map<int, int> tmp_jet_pair_flav;
        std::vector<int> z_idx;
        std::vector<int> h_idx;
        std::map<int, std::vector<int>> flavOccurance;
        int hz_check = 0 ;

        for (int i = 0; i < 4; ++i) {
            flavOccurance[maxScoreIdx[i]].push_back(i);
        } 
        //Make tmp flav pairs
        int pair = 0;    
        for (const auto& entry : flavOccurance) {
            if (entry.second.size() > 1) {  
                tmp_jet_pair_flav[pair]=entry.first;
                tmp_jet_pair[pair]=entry.second; 
                pair++;
                if (m_debug) std::cout << "Flav " << flavours[entry.first] << " found at indices: ";
                if (m_debug) {
                    for (size_t index : entry.second) {
                        std::cout << index << " ";
                    }
                    std::cout << std::endl;
                } 
            }
        } 
        // if (m_debug) std::cout << "jet_pair: " <<tmp_jet_pair.size()<<std::endl;
        if (tmp_jet_pair.size() == 2){
            //So given that there are only 4 jets reconstructed, I am pretty sure that this should work.
            if (m_debug) std::cout << "OLE! You have two paird! You can continue with your mission!"<<std::endl;
            jet_pair[0]= tmp_jet_pair[0];
            jet_pair[1]= tmp_jet_pair[1];
            jet_pair_flav[0]=tmp_jet_pair_flav[0];
            jet_pair_flav[1]=tmp_jet_pair_flav[1];
        }
        else if (tmp_jet_pair.size() == 3){
            std::cout << "THIS IS IMPOSSIBLE! I THINK???"<<std::endl;
            break;
        }
        else if (flavOccurance.size()==4){
                if (m_debug) std::cout << "No same flavour pair found. Ignore for now!" <<std::endl;
                continue;}
        else if (tmp_jet_pair.size() == 1){
            if (m_debug) std::cout << "Only 1 jet pair found!" <<std::endl;
            // if (m_debug) std::cout << "tmp_jet_pair[0].size(): " <<tmp_jet_pair[0].size()<<std::endl;
            if (flavOccurance.size()==3){
                if (m_debug) std::cout << "2 same flav, 3 diff flav -> 3 unique!"<<std::endl;
                //Need to sum up each flav and see what max?
                jet_pair[0] = tmp_jet_pair[0];
                jet_pair_flav[0]=tmp_jet_pair_flav[0];
                auto found_pair = jet_pair[0];
                //Finding jets that are not paired 
                std::vector<int> idx_jet ={0, 1, 2, 3};
                std::sort(found_pair.begin(), found_pair.end());
                std::sort(idx_jet.begin(), idx_jet.end());
                std::vector<int> missing_pair;
                std::set_symmetric_difference(found_pair.begin(), found_pair.end(), idx_jet.begin(), idx_jet.end(),
                    std::back_inserter(missing_pair));
                //Maybe this should be a a function! Outside this mess defined... Because you will need to call it again...
                float max_score = 0;
                int new_fl = -1;
                for (int fl = 0; fl < 7; ++fl) {
                    if (flavMap[missing_pair[0]][fl]+flavMap[missing_pair[1]][fl]>max_score){
                        new_fl=fl;
                        max_score=flavMap[missing_pair[0]][fl]+flavMap[missing_pair[1]][fl];
                        if (m_debug) std::cout << "Updating max score: " << max_score << ", flav is: "<<flavours[fl] <<std::endl;
                        };
                }
                jet_pair[1]=missing_pair;
                jet_pair_flav[1]=new_fl;
                if (m_debug) std::cout << "New pair found! With flavour " << flavours[jet_pair_flav[1]] << " found at indices: ";
                if (m_debug) {
                    for (size_t index : jet_pair[1]) {
                        std::cout << index << " ";
                    }
                    std::cout << std::endl;
                    }
            }
            else if (flavOccurance.size()==2){
                if (m_debug) std::cout << "3 same flav, 1 diff flav -> 2 unique!"<<std::endl;
                std::vector<int> new_pair;
                float max_sc = 0;
                for (size_t j1_in : tmp_jet_pair[0]){
                    for (size_t j2_in : tmp_jet_pair[0]){
                        if ( (j1_in!=j2_in)  && (j2_in>j1_in) && (flavMap[j1_in][tmp_jet_pair_flav[0]]+flavMap[j2_in][tmp_jet_pair_flav[0]] > max_sc) ){
                            max_sc = flavMap[j1_in][tmp_jet_pair_flav[0]]+flavMap[j2_in][tmp_jet_pair_flav[0]];
                            new_pair={};
                            new_pair.push_back(j1_in);
                            new_pair.push_back(j2_in);
                            if (m_debug) std::cout << "Updating max score: " << max_sc << ", with jet: "<<j1_in << " and "<< j2_in <<std::endl;

                        }
                    }
                }
                jet_pair[0]=new_pair;
                jet_pair_flav[0]=tmp_jet_pair_flav[0];
                if (m_debug) std::cout << "New pair found! With flavour " << flavours[jet_pair_flav[0]] << " found at indices: ";
                if (m_debug) {
                    for (size_t index : jet_pair[0]) {
                        std::cout << index << " ";
                    }
                    std::cout << std::endl;
                    }
                //Need to sum up each flav and see what max?
                auto found_pair = jet_pair[0];
                std::vector<int> id_jet ={0, 1, 2, 3};
                std::sort(found_pair.begin(), found_pair.end());
                std::sort(id_jet.begin(), id_jet.end());
                std::vector<int> missing_p;
                std::set_symmetric_difference(found_pair.begin(), found_pair.end(), id_jet.begin(), id_jet.end(),
                    std::back_inserter(missing_p));
                //Maybe this should be a a function! Outside this mess defined... Because you will need to call it again...
                float max_score = 0;
                int new_fl = -1;
                for (int fla = 0; fla < 7; ++fla) {
                    if (flavMap[missing_p[0]][fla]+flavMap[missing_p[1]][fla]>max_score){
                        new_fl=fla;
                        max_score=flavMap[missing_p[0]][fla]+flavMap[missing_p[1]][fla];
                        if (m_debug) std::cout << "Updating max score: " << max_score << ", flav is: "<<flavours[fla] <<std::endl;
                        };
                }
                jet_pair[1]=missing_p;
                jet_pair_flav[1]=new_fl;
                if (m_debug)  std::cout << "New pair found! With flavour " << flavours[jet_pair_flav[1]] << " found at indices: ";
                if (m_debug) {
                    for (size_t index : jet_pair[1]) {
                        std::cout << index << " ";
                    }
                    std::cout << std::endl;
                    }
            }
            if (flavOccurance.size()==1 || (jet_pair_flav[1] == jet_pair_flav[0]) ){
                hz_check++;
                if (m_debug) std::cout << "4 same flav jets, you continue and decide the paris based on the H/Z chi^2!"<<std::endl;
                // std::cout << "Flav PAIR 1:" <<jet_pair_flav[0]<<std::endl;
                // std::cout << "Flav PAIR 2:" <<jet_pair_flav[1]<<std::endl;
                float chi_z=0;
                float chi_h=0; 
                std::vector<float> chi_comb; 
                std::vector<vector<int>> z_ijets;
                std::vector<vector<int>> h_ijets;
                std::vector<int>z_pair;
                std::vector<int> idx_j ={0, 1, 2, 3};
                for (const auto& j : idx_j) {
                    for (const auto& jj : idx_j)  {
                        if (j<jj){
                            z_pair = {j,jj};
                            std::sort(z_pair.begin(), z_pair.end());
                            std::sort(idx_j.begin(), idx_j.end());
                            std::vector<int> h_pair;
                            std::set_symmetric_difference(z_pair.begin(), z_pair.end(), idx_j.begin(), idx_j.end(),std::back_inserter(h_pair));
                            chi_z = std::pow((LVjets[z_pair[0]]+LVjets[z_pair[1]]).M()-Z_mass,2);
                            chi_h = std::pow((LVjets[h_pair[0]]+LVjets[h_pair[1]]).M()-H_mass,2);
                            z_ijets.push_back(z_pair);
                            h_ijets.push_back(h_pair);
                            chi_comb.push_back(chi_z+chi_h);
                            // if (m_debug) {
                            //     std::cout << "Z-H CHI METHOD: "<< chi_z+chi_h  << std::endl;
                            //     std::cout << "Z indx: " << z_pair[0] << z_pair[1] << ", H indx: " << h_pair[0]<<h_pair[1]<< std::endl;
                            //     std::cout << "H mass: " << (LVjets[h_pair[0]]+LVjets[h_pair[1]]).M() << " Z mass: " << (LVjets[z_pair[0]]+LVjets[z_pair[1]]).M() << std::endl;
                            //     }

                        }
                    }
                }
                 // Using std::min_element to find the minimum element
                auto chi_comb_min = std::min_element(chi_comb.begin(), chi_comb.end());
                size_t minIndex = std::distance(chi_comb.begin(), chi_comb_min);
                jet_pair[0]=z_ijets[minIndex];
                jet_pair[1]=h_ijets[minIndex];
                z_idx=jet_pair[0];
                h_idx=jet_pair[1];
                H_flav=jet_pair_flav[1];
                Z_flav=jet_pair_flav[0];
                }
                }
        
        
        if (hz_check==0){
            float chi_p1 = std::pow((LVjets[jet_pair[0][0]]+LVjets[jet_pair[0][1]]).M()-Z_mass,2);
            float chi_p2 = std::pow((LVjets[jet_pair[1][0]]+LVjets[jet_pair[1][1]]).M()-Z_mass,2);
            if (m_debug) {
                std::cout << "Z CHI METHOD, chi 1: " <<chi_p1 << " chi 2:" <<chi_p2<<std::endl;
                }
            if (chi_p1<chi_p2){
                // if (m_debug) std::cout << "Pair 0 is Z!" <<std::endl;
                z_idx=jet_pair[0];
                Z_flav=jet_pair_flav[0];
                h_idx=jet_pair[1];
                H_flav=jet_pair_flav[1];
                
            }
            else{
                // if (m_debug) std::cout << "Pair 0 is H!" <<std::endl;
                // if (m_debug) std::cout<<"Z_flav " <<flavours[jet_pair_flav[1]]<< std::endl;
                // if (m_debug) std::cout<<"H_flav " <<flavours[jet_pair_flav[0]]<< std::endl;
                z_idx=jet_pair[1];
                Z_flav=jet_pair_flav[1];
                h_idx=jet_pair[0];
                H_flav=jet_pair_flav[0];
            }
        }

        float m_zjj = (LVjets[z_idx[0]]+LVjets[z_idx[1]]).M();
        float m_hjj = (LVjets[h_idx[0]]+LVjets[h_idx[1]]).M();

        float H_flav_sc = flavMap[h_idx[0]][H_flav]+flavMap[h_idx[1]][H_flav];
        float Z_flav_sc = flavMap[z_idx[0]][Z_flav]+flavMap[z_idx[1]][Z_flav];
        if (m_debug) {
            std::cout <<"H: " << m_hjj <<" sc: " <<H_flav_sc<<std::endl;
            std::cout << "Z: "<< m_zjj<< " sc: "<<Z_flav_sc<<std::endl;
        }
        NafterPairing++;
        Incl_obsHist->Fill(m_zjj,m_hjj);
        //correcte m_hjj
        float m_hjj_corr = m_hjj + m_zjj - Z_mass;

        //A bit of selection 
        float WW_cuts = sqrt(pow( m_zjj-W_mass ,2) + pow( m_hjj-W_mass ,2));
        float ZZ_cuts = sqrt(pow( m_zjj-Z_mass ,2) + pow( m_hjj-Z_mass ,2));
        // if (m_debug) {
        //     std::cout << "WW_cuts: " <<sqrt(pow( m_zjj-W_mass ,2) + pow( m_hjj-W_mass ,2))<<std::endl;
        //     std::cout << "ZZ_cuts: " <<sqrt(pow( m_zjj-Z_mass ,2) + pow( m_hjj-Z_mass ,2))<<std::endl;}
        if(WW_cuts<=10) continue;
        if(ZZ_cuts<=10) continue;
        NjjMassCut++;
        if (50. >= m_zjj)continue; 
        if (m_zjj >= H_mass) continue; 
        //Fix
        if (m_hjj_corr<=Z_mass) continue;
        // For now remove the Htautau as a cat, only bkg....
        // if (H_flav == 6) continue;
        //After parining and cuts 
        //The flag e correction falg...
        float flag_ecorr = 0.0;
        for (size_t je = 0; je < 4; ++je) { 
            if (jet_e.at(je) > 240.0 || jet_e.at(je) < 0.0) {
                flag_ecorr += 1000.0;}
        }
        if (flag_ecorr>=1000.) continue;

        

        NafterSel++;
        ZZ_cut->Fill(ZZ_cuts);
        WW_cut->Fill(WW_cuts);
        mZ->Fill(m_zjj);
        mH->Fill(m_hjj);
        mH_corr->Fill(m_hjj_corr);
        h_vis_M->Fill(vis_M());
        h_vis_E->Fill(vis_E());
        h_vis_theta->Fill(vis_theta());
        h_muons_p->Fill(muons_p());
        h_electrons_p->Fill(electrons_p());
        h_event_nmu->Fill(event_nmu());
        h_event_nel->Fill(event_nel());
        h_d_12->Fill(d_12());
        h_d_23->Fill(d_23());
        h_d_34->Fill(d_34());
        
        if (Z_flav == 0){
            Z_Bscore->Fill(Z_flav_sc);
        }
        else if (Z_flav == 1){
            Z_Cscore->Fill(Z_flav_sc);
        }
        else if (Z_flav == 2){
            Z_Sscore->Fill(Z_flav_sc);
        }
        else if (Z_flav == 5){
            Z_Gscore->Fill(Z_flav_sc);
        }
        else if (Z_flav == 4 || Z_flav == 3){
            Z_Qscore->Fill(Z_flav_sc);
        }
        
        if (Z_flav == 5 || Z_flav == 6 || H_flav == 6 || H_flav == 4 || H_flav == 3) continue;
        Nfit++; 
        if (H_flav == 0){
            BlikeEvents++;
            // Hbb_obsHist->Fill(m_zjj,m_hjj);
            // std::cout << "B Flav" << std::endl;
            Bscore->Fill(H_flav_sc);
            if (H_flav_sc < 1.1)
            {
                // std::cout << "B Low Cat" << std::endl;
                if (Z_flav == 0){
//                     Low_bbZ_Hbb_obsHist->Fill(m_zjj,m_hjj);
                    Low_bbZ_Hbb_obsHist->Fill(m_hjj_corr,m_zjj);

                }
                else if (Z_flav == 1){
//                     Low_ccZ_Hbb_obsHist->Fill(m_zjj,m_hjj);
                    Low_ccZ_Hbb_obsHist->Fill(m_hjj_corr,m_zjj);


                }
                else if (Z_flav == 2){
//                     Low_ssZ_Hbb_obsHist->Fill(m_zjj,m_hjj);
                    Low_ssZ_Hbb_obsHist->Fill(m_hjj_corr,m_zjj);

                }
                else if (Z_flav == 4 || Z_flav == 3 ){
//                     Low_qqZ_Hbb_obsHist->Fill(m_zjj,m_hjj);
                    Low_qqZ_Hbb_obsHist->Fill(m_hjj_corr,m_zjj);

                }
                
                BlikeEvents_cat[0]++;
            }
            else if ((H_flav_sc >= 1.1) && (H_flav_sc <= 1.8))
            {
                // std::cout << "B Mid Cat" << std::endl;
                // MidB_obsHist->Fill(m_zjj,m_hjj);
                if (Z_flav == 0){
//                     Mid_bbZ_Hbb_obsHist->Fill(m_zjj,m_hjj);
                    Mid_bbZ_Hbb_obsHist->Fill(m_hjj_corr,m_zjj);


                }
                else if (Z_flav == 1){
//                     Mid_ccZ_Hbb_obsHist->Fill(m_zjj,m_hjj);
                    Mid_ccZ_Hbb_obsHist->Fill(m_hjj_corr,m_zjj);

                }
                else if (Z_flav == 2){
//                     Mid_ssZ_Hbb_obsHist->Fill(m_zjj,m_hjj);
                    Mid_ssZ_Hbb_obsHist->Fill(m_hjj_corr,m_zjj);

                }
                else if (Z_flav == 4 || Z_flav == 3 ){
//                     Mid_qqZ_Hbb_obsHist->Fill(m_zjj,m_hjj);
                    Mid_qqZ_Hbb_obsHist->Fill(m_hjj_corr,m_zjj);
                }
                BlikeEvents_cat[1]++;
            }
            else if (H_flav_sc > 1.8)
            {
                // std::cout << "B Hi Cat" << std::endl;
                // HiB_obsHist->Fill(m_zjj,m_hjj);
                if (Z_flav == 0){
//                     Hi_bbZ_Hbb_obsHist->Fill(m_zjj,m_hjj);
                    Hi_bbZ_Hbb_obsHist->Fill(m_hjj_corr,m_zjj);

                }
                else if (Z_flav == 1){
//                     Hi_ccZ_Hbb_obsHist->Fill(m_zjj,m_hjj);
                    Hi_ccZ_Hbb_obsHist->Fill(m_hjj_corr,m_zjj);

                }
                else if (Z_flav == 2){
//                     Hi_ssZ_Hbb_obsHist->Fill(m_zjj,m_hjj);
                    Hi_ssZ_Hbb_obsHist->Fill(m_hjj_corr,m_zjj);

                }
                else if (Z_flav == 4 || Z_flav == 3 ){
//                     Hi_qqZ_Hbb_obsHist->Fill(m_zjj,m_hjj);
                     Hi_qqZ_Hbb_obsHist->Fill(m_hjj_corr,m_zjj);
                }
                BlikeEvents_cat[2]++;
            }
        }
        else if (H_flav == 1){
            ClikeEvents++;
            Cscore->Fill(H_flav_sc);
            // Hcc_obsHist->Fill(m_zjj,m_hjj);
            if (H_flav_sc < 1.1)
            {
                // LowC_obsHist->Fill(m_zjj,m_hjj);
                if (Z_flav == 0){
//                     Low_bbZ_Hcc_obsHist->Fill(m_zjj,m_hjj);
                        Low_bbZ_Hcc_obsHist->Fill(m_hjj_corr,m_zjj);


                }
                else if (Z_flav == 1){
//                     Low_ccZ_Hcc_obsHist->Fill(m_zjj,m_hjj);
                    Low_ccZ_Hcc_obsHist->Fill(m_hjj_corr,m_zjj);

                }
                else if (Z_flav == 2){
//                     Low_ssZ_Hcc_obsHist->Fill(m_zjj,m_hjj);
                    Low_ssZ_Hcc_obsHist->Fill(m_hjj_corr,m_zjj);

                }
                else if (Z_flav == 4 || Z_flav == 3 ){
//                     Low_qqZ_Hcc_obsHist->Fill(m_zjj,m_hjj);
                    Low_qqZ_Hcc_obsHist->Fill(m_hjj_corr,m_zjj);
                }
                ClikeEvents_cat[0]++;
            }
            else if ((H_flav_sc >= 1.1) && (H_flav_sc <= 1.8))
            {
                // MidC_obsHist->Fill(m_zjj,m_hjj);
                if (Z_flav == 0){
//                     Mid_bbZ_Hcc_obsHist->Fill(m_zjj,m_hjj);
                    Mid_bbZ_Hcc_obsHist->Fill(m_hjj_corr,m_zjj);

                }
                else if (Z_flav == 1){
//                     Mid_ccZ_Hcc_obsHist->Fill(m_zjj,m_hjj);
                     Mid_ccZ_Hcc_obsHist->Fill(m_hjj_corr,m_zjj);

                }
                else if (Z_flav == 2){
//                     Mid_ssZ_Hcc_obsHist->Fill(m_zjj,m_hjj);
                    Mid_ssZ_Hcc_obsHist->Fill(m_hjj_corr,m_zjj);

                }
                else if (Z_flav == 4 || Z_flav == 3 ){
//                     Mid_qqZ_Hcc_obsHist->Fill(m_zjj,m_hjj);
                    Mid_qqZ_Hcc_obsHist->Fill(m_hjj_corr,m_zjj);
                }
                ClikeEvents_cat[1]++;
            }
            else if (H_flav_sc > 1.8)
            {
                // HiC_obsHist->Fill(m_zjj,m_hjj);
                if (Z_flav == 0){
//                     Hi_bbZ_Hcc_obsHist->Fill(m_zjj,m_hjj);
                    Hi_bbZ_Hcc_obsHist->Fill(m_hjj_corr,m_zjj);

                }
                else if (Z_flav == 1){
//                     Hi_ccZ_Hcc_obsHist->Fill(m_zjj,m_hjj);
                    Hi_ccZ_Hcc_obsHist->Fill(m_hjj_corr,m_zjj);

                }
                else if (Z_flav == 2){
//                     Hi_ssZ_Hcc_obsHist->Fill(m_zjj,m_hjj);
                    Hi_ssZ_Hcc_obsHist->Fill(m_hjj_corr,m_zjj);

                }
                else if (Z_flav == 4 || Z_flav == 3 ){
//                     Hi_qqZ_Hcc_obsHist->Fill(m_zjj,m_hjj);
                    Hi_qqZ_Hcc_obsHist->Fill(m_hjj_corr,m_zjj);
                }
                ClikeEvents_cat[2]++;
            }
        }
        else if (H_flav == 2){
            Sscore->Fill(H_flav_sc);
            SlikeEvents++;
            // Hss_obsHist->Fill(m_zjj,m_hjj);
            if (H_flav_sc < 0.8)
            {
//                 LowS_obsHist->Fill(m_zjj,m_hjj);
                LowS_obsHist->Fill(m_hjj_corr,m_zjj);
                SlikeEvents_cat[0]++;
            }
            else if ((H_flav_sc >= 0.8) && (H_flav_sc <= 1.4))
            {
//                 MidS_obsHist->Fill(m_zjj,m_hjj);
                MidS_obsHist->Fill(m_hjj_corr,m_zjj);
                SlikeEvents_cat[1]++;
            }
            else if (H_flav_sc > 1.4)
            {
//                 HiS_obsHist->Fill(m_zjj,m_hjj);
                HiS_obsHist->Fill(m_hjj_corr,m_zjj);
                SlikeEvents_cat[2]++;
            }
        }
        else if (H_flav == 5){
            GlikeEvents++;
            Gscore->Fill(H_flav_sc);
            // Hgg_obsHist->Fill(m_zjj,m_hjj);
            if (H_flav_sc < 1.1)
            {
                // LowG_obsHist->Fill(m_zjj,m_hjj);
                if (Z_flav == 0){
//                     Low_bbZ_Hgg_obsHist->Fill(m_zjj,m_hjj);
                    Low_bbZ_Hgg_obsHist->Fill(m_hjj_corr,m_zjj);

                }
                else if (Z_flav == 1){
//                     Low_ccZ_Hgg_obsHist->Fill(m_zjj,m_hjj);
                    Low_ccZ_Hgg_obsHist->Fill(m_hjj_corr,m_zjj);

                }
                else if (Z_flav == 2){
//                     Low_ssZ_Hgg_obsHist->Fill(m_zjj,m_hjj);
                    Low_ssZ_Hgg_obsHist->Fill(m_hjj_corr,m_zjj);

                }
                else if (Z_flav == 4 || Z_flav == 3 ){
//                     Low_qqZ_Hgg_obsHist->Fill(m_zjj,m_hjj);
                    Low_qqZ_Hgg_obsHist->Fill(m_hjj_corr,m_zjj);
                }
                GlikeEvents_cat[0]++;
            }
            else if ((H_flav_sc >= 1.1) && (H_flav_sc <= 1.8))
            {
                // MidG_obsHist->Fill(m_zjj,m_hjj);
                if (Z_flav == 0){
//                     Mid_bbZ_Hgg_obsHist->Fill(m_zjj,m_hjj);
                    Mid_bbZ_Hgg_obsHist->Fill(m_hjj_corr,m_zjj);

                }
                else if (Z_flav == 1){
//                     Mid_ccZ_Hgg_obsHist->Fill(m_zjj,m_hjj);
                    Mid_ccZ_Hgg_obsHist->Fill(m_hjj_corr,m_zjj);

                }
                else if (Z_flav == 2){
//                     Mid_ssZ_Hgg_obsHist->Fill(m_zjj,m_hjj);
                    Mid_ssZ_Hgg_obsHist->Fill(m_hjj_corr,m_zjj);

                }
                else if (Z_flav == 4 || Z_flav == 3 ){
//                     Mid_qqZ_Hgg_obsHist->Fill(m_zjj,m_hjj);
                    Mid_qqZ_Hgg_obsHist->Fill(m_hjj_corr,m_zjj);
                }
                GlikeEvents_cat[1]++;
            }
            else if (H_flav_sc > 1.8)
            {
                // HiG_obsHist->Fill(m_zjj,m_hjj);
                if (Z_flav == 0){
//                     Hi_bbZ_Hgg_obsHist->Fill(m_zjj,m_hjj);
                    Hi_bbZ_Hgg_obsHist->Fill(m_hjj_corr,m_zjj);

                }
                else if (Z_flav == 1){
//                     Hi_ccZ_Hgg_obsHist->Fill(m_zjj,m_hjj);
                    Hi_ccZ_Hgg_obsHist->Fill(m_hjj_corr,m_zjj);

                }
                else if (Z_flav == 2){
//                     Hi_ssZ_Hgg_obsHist->Fill(m_zjj,m_hjj);
                    Hi_ssZ_Hgg_obsHist->Fill(m_hjj_corr,m_zjj);

                }
                else if (Z_flav == 4 || Z_flav == 3 ){
//                     Hi_qqZ_Hgg_obsHist->Fill(m_zjj,m_hjj);
                    Hi_qqZ_Hgg_obsHist->Fill(m_hjj_corr,m_zjj);
                }
                GlikeEvents_cat[2]++;
            }
        }
        else if (H_flav == -1){
            std::cout << "Something went wrong! No Higgs found!"<<std::endl;
            break;
            }
        else{
            std::cout<<"Couldn't find best score"<<std::endl;
        }
        // }
        // Step 2: check efficiency ~ 60-70%, for each cat you fit you need to check eff 
                // else if (H_flav == 6){
        //     TAUlikeEvents++;
        //     TAUscore->Fill(H_flav_sc);
        //     // Htautau_obsHist->Fill(m_zjj,m_hjj);
        //     if (H_flav_sc < 1.1)
        //     {
        //         LowTAU_obsHist->Fill(m_zjj,m_hjj);
        //         TAUlikeEvents_cat[0]++;
        //     }
        //     else if ((H_flav_sc >= 1.1) && (H_flav_sc <= 1.9))
        //     {
        //         MidTAU_obsHist->Fill(m_zjj,m_hjj);
        //         TAUlikeEvents_cat[1]++;
        //     }
        //     else
        //     {
        //         HiTAU_obsHist->Fill(m_zjj,m_hjj);
        //         TAUlikeEvents_cat[2]++;
        //     }
        // }

    }
    scoreMapHist->SetBinContent(1, BlikeEvents*100./NafterSel);
    scoreMapHist->SetBinContent(2, ClikeEvents*100./NafterSel);
    scoreMapHist->SetBinContent(3, SlikeEvents*100./NafterSel);
    scoreMapHist->SetBinContent(4, QlikeEvents*100./NafterSel);
    scoreMapHist->SetBinContent(5, GlikeEvents*100./NafterSel);
    // scoreMapHist->SetBinContent(6, TAUlikeEvents*100./NafterSel);

    CountsMapHist->SetBinContent(1, BlikeEvents);
    CountsMapHist->SetBinContent(2, ClikeEvents);
    CountsMapHist->SetBinContent(3, SlikeEvents);
    CountsMapHist->SetBinContent(4, QlikeEvents);
    CountsMapHist->SetBinContent(5, GlikeEvents);
    // CountsMapHist->SetBinContent(6, TAUlikeEvents);

    
    //CutFlow

    cutFlowHist->SetBinContent(1,NEventsInt);
    cutFlowHist->SetBinContent(2,NjetCut);
    cutFlowHist->SetBinContent(3,NleptonCut);
    cutFlowHist->SetBinContent(4,NkineCut);
    cutFlowHist->SetBinContent(5,NdCutd123);
    cutFlowHist->SetBinContent(6,NdCutd34);
    cutFlowHist->SetBinContent(7,NafterPairing);
    cutFlowHist->SetBinContent(8,NjjMassCut);
    cutFlowHist->SetBinContent(9,NafterSel);
    cutFlowHist->SetBinContent(10,Nfit);

    for(int i = 0; i < 3; i++)
    {
        scoreMapFitCatHist->SetBinContent(1+i*5, BlikeEvents_cat[i]*100./NafterSel);
        scoreMapFitCatHist->SetBinContent(2+i*5, ClikeEvents_cat[i]*100./NafterSel);
        scoreMapFitCatHist->SetBinContent(3+i*5, SlikeEvents_cat[i]*100./NafterSel);
        scoreMapFitCatHist->SetBinContent(4+i*5, GlikeEvents_cat[i]*100./NafterSel);
        scoreMapFitCatHist->SetBinContent(5+i*5, QlikeEvents_cat[i]*100./NafterSel);
        // scoreMapFitCatHist->SetBinContent(6+i*6, TAUlikeEvents_cat[i]*100./NafterSel);

        CountsFitCatHist->SetBinContent(1+i*4, BlikeEvents_cat[i]);
        CountsFitCatHist->SetBinContent(2+i*4, ClikeEvents_cat[i]);
        CountsFitCatHist->SetBinContent(3+i*4, SlikeEvents_cat[i]);
        CountsFitCatHist->SetBinContent(4+i*4, GlikeEvents_cat[i]);
//         CountsFitCatHist->SetBinContent(5+i*5, QlikeEvents_cat[i]);
        // CountsFitCatHist->SetBinContent(6+i*6, TAUlikeEvents_cat[i]);
    }


    //Make histograms for the cats 

}


