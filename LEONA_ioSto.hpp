//
//  LEONA_ioSto.hpp
//  LEONA1.0
//
//  Created by Shuotao Diao on 10/24/22.
//

#ifndef LEONA_ioSto_hpp
#define LEONA_ioSto_hpp

#include <stdio.h>

#include "LEONA_dataStructure.hpp"

// including be, bi, Ce and Ci
secondStageRHSmap readStochasticMap(const std::string& stochasticPath);

// merge randomVector
secondStageRHSpoint merge_randomVector(const dataPoint& be_point, const dataPoint& bi_point, const dataPoint& Ce_point, const dataPoint& Ci_point);
#endif /* LEONA_ioSto_hpp */
