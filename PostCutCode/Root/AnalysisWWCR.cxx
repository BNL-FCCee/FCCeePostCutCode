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
    // Get the histograms
    auto countingHist = m_histContainer->getCountingHist();
    auto truth_W_e =  m_histContainer->get1DHist("truth_W_e", 100, 0, 200);
    // Get the trees
    auto treeCont = std::make_shared<TreeContainer>();
 
    // Get max events to run on
    int nEntries = treeCont->getEventsToRun();

    // Connect branches to trees
    auto tree = treeCont->getTree();
    varMember<int> event_njet {tree, "event_njet"};
    varMember<float>  truth_Wp_HS_e {tree, "truth_Wp_HS_e"};

    // Loop over the trees here
    for(int i = 0; i < nEntries; i++)
    {
        treeCont->getEntry(i);
        // Just to store how many events were run over
        countingHist->Fill(1);
        if (m_debug) std::cout<<"IN WW CR!"<<std::endl;
        if(i % 10000 == 0) std::cout<<"Done i: "<<i<<" out of "<<nEntries<<std::endl;


        if(event_njet() != 4) continue; //Require EXACTLY 4 jet! This SHOULD be the case!
        if (m_debug) std::cout<<"HAS 4 JETS EXACTLY!"<<std::endl;
       
        truth_W_e->Fill(truth_Wp_HS_e());
    }


}


