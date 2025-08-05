#ifndef _AnalysisWWCR_H
#define _AnalysisWWCR_H

// c++ includes
#include <map>
#include <vector>
#include <iostream>
#include <memory>
#include <algorithm>

// local include
#include "PostCutCode/AnalysisBase.h"

class AnalysisWWCR: public AnalysisBase
{
   public:
      bool m_debug;
      AnalysisWWCR();
      virtual ~AnalysisWWCR();

      void run() override;

   protected:
      std::vector<double> ee_corr;
      std::vector<double> theta_corr;
      std::vector<double> phi_corr;
      std::vector<double> chi_corr;
      std::vector<double> delR_corr;
      std::vector<double> p_jet_const_1;
      std::vector<double> p_jet_const_2;

      std::vector<double> p_reco_Jet1;
      std::vector<double> p_reco_Jet2;
      std::vector<double> p_reco_Jet3;
      std::vector<double> p_reco_Jet4;
      std::vector<double> theta_Jet1_Jet2;
      std::vector<double> theta_Jet3_Jet4;
      std::vector<double> theta_Jet1_Jet3;
      std::vector<double> theta_Jet2_Jet4;
      std::vector<double> theta_Jet1_Jet4;
      std::vector<double> theta_Jet2_Jet3;
      std::vector<double> phi_Jet1_Jet2;
      std::vector<double> phi_Jet3_Jet4;
      std::vector<double> phi_Jet1_Jet3;
      std::vector<double> phi_Jet2_Jet4;
      std::vector<double> phi_Jet1_Jet4;
      std::vector<double> phi_Jet2_Jet3;
      std::vector<double> chi_Jet1_Jet2;
      std::vector<double> chi_Jet3_Jet4;
      std::vector<double> chi_Jet1_Jet3;
      std::vector<double> chi_Jet2_Jet4;
      std::vector<double> chi_Jet1_Jet4;
      std::vector<double> chi_Jet2_Jet3;

      std::vector<double> mass_W1;
      std::vector<double> mass_W2;
      std::vector<double> p_W1;
      std::vector<double> p_W2;
      std::vector<double> e_W1;
      std::vector<double> e_W2;
      // std::vector<double> delta_R_W_boson;

      double mc_weight;

};


#endif

