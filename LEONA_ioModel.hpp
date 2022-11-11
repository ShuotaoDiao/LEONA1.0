//
//  LEONA_ioModel.hpp
//  LEONA1.0
//
//  Created by Shuotao Diao on 10/24/22.
//

#ifndef LEONA_ioModel_hpp
#define LEONA_ioModel_hpp

#include <stdio.h>

#include "LEONA_dataStructure.hpp"

// convert into standard form
/* |De 0| |y| +  |Ce| x = |be|
   |Di I| |s|    |Ci|     |bi|
 y,s >=0
 */
/* D = |De 0|
       |Di I|
 */
std::vector<std::vector<double>> standard_D(const std::vector<std::vector<double>>& De, const std::vector<std::vector<double>>& Di);
/* C = |Ce|
       |Ci|
 */
std::vector<std::vector<double>> standard_C(const std::vector<std::vector<double>>& Ce,
                                            const std::vector<std::vector<double>>& Ci);
/* e = |be|
       |bi|
 */
std::vector<double> standard_e(const std::vector<double>& be, const std::vector<double>& bi);
// standard model parameters
standardTwoStageParameters readStandardTwoStageParameters(const std::string& parameterPath);

#endif /* LEONA_ioModel_hpp */
