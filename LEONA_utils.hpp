//
//  LEONA_utils.hpp
//  LEONA1.0
//
//  Created by Shuotao Diao on 10/24/22.
//

#ifndef LEONA_utils_hpp
#define LEONA_utils_hpp

#include <stdio.h>

#include "LEONA_solver.hpp"

validationResult twoStageLP_validation_outputResults(const std::string& folder_path, const std::vector<double>& x_candidate);


void interface_leona_knn(const std::string& folder_path,
                                    const std::string& validation_folder_path,
                                    int max_it,
                                    const std::vector<double>& x_init,
                                    const std::vector<double>& observed_predictor,
                                    int N_pre,
                                    int N_incre,
                                    double Dx,
                                    int m,
                                    double max_subgradient);

void interface_leona_kernel(const std::string& folder_path,
                            const std::string& validation_folder_path,
                            int max_it,
                            const std::vector<double>& x_init,
                            const std::vector<double>& observed_predictor,
                            int N_pre,
                            int N_incre,
                            int flag_kernel,
                            double Dx,
                            int m,
                            double max_subgradient);
#endif /* LEONA_utils_hpp */
