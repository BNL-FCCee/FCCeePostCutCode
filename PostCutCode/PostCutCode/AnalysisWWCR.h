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


};


#endif

