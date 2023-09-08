//
//  LEONA_solver.hpp
//  LEONA1.0
//
//  Created by Shuotao Diao on 10/24/22.
//

#ifndef LEONA_solver_hpp
#define LEONA_solver_hpp

#include <stdio.h>
#include <stdlib.h> // rand
#include <ctime>
#include <cmath>

#include "LEONA_dataStructure.hpp"
#include "LEONA_ioDB.hpp"
#include "LEONA_ioSto.hpp"
#include "LEONA_ioModel.hpp"


// projection for the first stage in the two stage linear programming
std::vector<double> leona_projection(const std::vector<double>& x, standardTwoStageParameters& model_parameters);

// LEONA kNN solver
solverOutput leona_knn(const std::string& folder_path,
                       int max_it,
                       const std::vector<double>& x_init,
                       const std::vector<double>& observed_predictor,
                       int N_pre,
                       int N_incre,
                       double Dx,
                       int m,
                       double max_subgradient);

// LEONA kernel solver
solverOutput leona_kernel(const std::string& folder_path,
                          int max_it,
                          const std::vector<double>& x_init,
                          const std::vector<double>& observed_predictor,
                          int N_pre,
                          int N_incre,
                          int flag_kernel,
                          double Dx,
                          int m,
                          double max_subgradient);

// output a sequence of decisions
solverOutput leona_knn2(const std::string& folder_path,
                       int max_it,
                       const std::vector<int>& it_pointer,
                       const std::vector<double>& x_init,
                       const std::vector<double>& observed_predictor,
                       int N_pre,
                       int N_incre,
                       double Dx,
                       int m,
                       double max_subgradient);


solverOutput leona_kernel2(const std::string& folder_path,
                          int max_it,
                          const std::vector<int>& it_pointer,
                          const std::vector<double>& x_init,
                          const std::vector<double>& observed_predictor,
                          int N_pre,
                          int N_incre,
                          int flag_kernel,
                          double Dx,
                          int m,
                          double max_subgradient);
#endif /* LEONA_solver_hpp */
