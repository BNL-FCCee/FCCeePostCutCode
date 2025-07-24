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
      double mc_weight;

};


#endif

