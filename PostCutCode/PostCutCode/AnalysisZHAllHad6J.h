#ifndef _AnalysisZHAllHad6J_H
#define _AnalysisZHAllHad6J_H

// c++ includes
#include <map>
#include <vector>
#include <iostream>
#include <memory>
#include <algorithm>
//ROOT things

// local include
#include "PostCutCode/AnalysisBase.h"

class AnalysisZHAllHad6J: public AnalysisBase
{
   public:
      bool m_debug;
      AnalysisZHAllHad6J();
      virtual ~AnalysisZHAllHad6J();
      

      int  getMaxScoreFlav(const std::vector<int>& vec, const std::map<int, std::vector<float>>& mapFlav);
      void getMissingPair(const std::vector<int>& foundPair, std::vector<int>& missingPair);
      void printZcandidate(struct ZCandidate zCandidate);


      // std::vector<unsigned int> missing_pair(std::vector<int>);

      void run() override;

   protected:
      double mH_ZZ;
      double mZ1_jj;
      double mZ2_jj;
      double mZ3_jj;
      int Z1_j1FlavorIdx;
      int Z1_j2FlavorIdx;
      int Z2_j1FlavorIdx;
      int Z2_j2FlavorIdx;
      int Zstar_j1FlavorIdx;
      int Zstar_j2FlavorIdx;
      double Z1_FlavorScoreCombined;
      double Z2_FlavorScoreCombined;
      double Z3_FlavorScoreCombined;
      double DeltaR_HZZ_Z2;
      double DeltaPhi_HZZ_Z2;
      double Hscore;
      double Zscore;
      int Hflav;
      int Zflav;
      double ChiH;
      double ChiZ;
      double mc_weight;
      double v_d12;
      double v_d23;
      double v_d34;
      double v_d45;
      double v_d56;
      double v_vis_E;
      double v_vis_M;
      // double Jet0_E;
      // double Jet0_px;
      // double Jet0_py;
      // double Jet0_pz;
      // double Jet1_E;
      // double Jet1_px;
      // double Jet1_py;
      // double Jet1_pz;
      // double Jet2_E;
      // double Jet2_px;
      // double Jet2_py;
      // double Jet2_pz;
      // double Jet3_E;
      // double Jet3_px;
      // double Jet3_py;
      // double Jet3_pz;
      // double Jet4_E;
      // double Jet4_px;
      // double Jet4_py;
      // double Jet4_pz;
      // double Jet5_E;
      // double Jet5_px;
      // double Jet5_py;
      // double Jet5_pz;
      
};


#endif

