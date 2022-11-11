//
//  LEONA_solver.cpp
//  LEONA1.0
//
//  Created by Shuotao Diao on 10/24/22.
//

#include "LEONA_solver.hpp"

// declare global variables (need to be defined in the source file)
const double SOLVER_PRECISION_LOWER = -1e-6;
const double SOLVER_PRECISION_UPPER = 1e-6;
const double SOLVER_INF = 1e10;

// projection for the first stage in the two stage linear programming
std::vector<double> leona_projection(const std::vector<double>& x, standardTwoStageParameters& model_parameters) {
    // solve a quadratic programming
    IloEnv env;
    IloModel mod(env);
    IloNumVarArray x_temp(env,x.size(),-IloInfinity,IloInfinity,ILOFLOAT);
    mod.add(x_temp);
    IloExpr expr_obj(env);
    for (int x_index = 0; x_index < x.size(); ++x_index) {
        expr_obj += x_temp[x_index] * x_temp[x_index] - 2.0 * x_temp[x_index] * x[x_index];
    }
    IloObjective obj = IloMinimize(env,expr_obj); // objective function
    mod.add(obj);
    // constraints
    std::vector<IloExpr> exprs;
    for (int index_cons = 0; index_cons < model_parameters.A.getRowLength(); ++index_cons) {
        IloExpr expr(env);
        exprs.push_back(expr);
    }
    // A x
    for (int col_idx = 0; col_idx < model_parameters.A.getColLength(); ++col_idx) {
        int beg_idx = model_parameters.A.getCbeg(col_idx);
        for (int idx = beg_idx; idx < model_parameters.A.getClen(col_idx) + beg_idx; ++idx) {
            exprs[model_parameters.A.getRow(idx)] += model_parameters.A.getVal(idx) * x_temp[col_idx];
        }
    }
    // right hand side b
    for (int idx = 0; idx < model_parameters.b.getNzeroLen(); ++idx) {
        exprs[model_parameters.b.getLoc(idx)] -= model_parameters.b.getVal(idx);
    }
    // add constraints Ax <= b
    for (int idx = 0; idx < model_parameters.A.getRowLength(); ++idx) {
        mod.add(exprs[idx] <= 0);
    }
    // create cplex environment
    IloCplex cplex(env);
    cplex.extract(mod);
    cplex.setOut(env.getNullStream());
    cplex.solve();
    std::vector<double> x_proj;
    // obtain the projected point
    for (int x_index = 0; x_index < model_parameters.A.getColLength(); ++x_index) {
        x_proj.push_back(cplex.getValue(x_temp[x_index]));
        //std::cout << cplex.getValue(x_temp[x_index]) << std::endl;
    }
    env.end();
    // return results
    return x_proj;
}


// LEON solver
solverOutput leona_knn(const std::string& folder_path,
                       int max_it,
                       const std::vector<double>& x_init,
                       const std::vector<double>& observed_predictor,
                       int N_pre,
                       int N_incre,
                       double Dx,
                       int m,
                       double max_subgradient) {
    // timer
    std::clock_t time_start;
    time_start = std::clock();
    // current time
    //std::time_t currTime = std::time(nullptr);
    // STEP 1: INITIALIZATION
    //double beta = 0.5; 0.6 is okay
    double beta = 0.6; // 0 < beta < 1
    //int k = 1;
    int k_new = 1;
    //int N = 0;
    std::vector<double> distanceSet;
    std::vector<int> orderSet;
    std::vector<int> kNNSet;
    double knn_radius = 0;
    bool flag_be; // tell if be stochastic is generated
    bool flag_bi; // tell if bi stochastic is generated
    bool flag_Ce; // tell if Ce stochastic is generated
    bool flag_Ci; // tell if Ci stochastic is generated
    std::vector<secondStageRHSpoint> RHS_dataset;
    // create directory paths for database and model
    std::string be_DB_path = folder_path + "/be_DB.txt";
    std::string bi_DB_path = folder_path + "/bi_DB.txt";
    std::string Ce_DB_path = folder_path + "/Ce_DB.txt";
    std::string Ci_DB_path = folder_path + "/Ci_DB.txt";
    std::string model_path = folder_path + "/model.txt";
    std::string sto_path = folder_path + "/sto.txt";
    std::string resultsOutput_path = folder_path + "/computationalResults(leona_knn1.0).txt";
    // convert all the paths into constant chars
    const char* be_DB_path_const = be_DB_path.c_str();
    const char* bi_DB_path_const = bi_DB_path.c_str();
    const char* Ce_DB_path_const = Ce_DB_path.c_str();
    const char* Ci_DB_path_const = Ci_DB_path.c_str();
    // create stream object
    std::ifstream readFile_be(be_DB_path_const);
    std::ifstream readFile_bi(bi_DB_path_const);
    std::ifstream readFile_Ce(Ce_DB_path_const);
    std::ifstream readFile_Ci(Ci_DB_path_const);
    // create database
    std::vector<std::vector<dataPoint>> be_DB;
    std::vector<std::vector<dataPoint>> bi_DB;
    std::vector<std::vector<dataPoint>> Ce_DB;
    std::vector<std::vector<dataPoint>> Ci_DB;
    // read  be
    if (readFile_be.is_open()) {
        std::cout << "be_DB is found." << std::endl;
        readFile_be.close(); // close the file
        // read be database
        be_DB = readNonparametricDB(be_DB_path);
        flag_be = true;
    }
    else {
        readFile_be.close(); // close the file
        flag_be = false;
        std::cout << "be_DB is not found!" << std::endl;
    }
    // read bi
    if (readFile_bi.is_open()) {
        std::cout << "bi_DB is found." << std::endl;
        readFile_be.close(); // close the file
        // read bi database
        bi_DB = readNonparametricDB(bi_DB_path);
        flag_bi = true;
    }
    else {
        readFile_bi.close(); // close the file
        flag_bi = false;
        std::cout << "bi_DB is not found!" << std::endl;
    }
    // read Ce
    if (readFile_Ce.is_open()) {
        std::cout << "Ce_DB stochastic part is found." << std::endl;
        readFile_Ce.close(); // close the file
        // Ce database
        Ce_DB = readNonparametricDB(Ce_DB_path);
        flag_Ce = true;
    }
    else {
        readFile_Ce.close(); // close the file
        flag_Ce = false;
        std::cout << "Ce_DB is not found!" << std::endl;
    }
    // read Ci
    if (readFile_Ci.is_open()) {
        std::cout << "Ci_DB stochastic part is found." << std::endl;
        readFile_Ci.close(); // close the file
        // Ci database
        Ci_DB = readNonparametricDB(Ci_DB_path);
        flag_Ci = true;
    }
    else {
        readFile_Ci.close(); // close the file
        flag_Ci = false;
        std::cout << "Ci_DB is not found!" << std::endl;
    }
    // read model file
    // create model structure
    standardTwoStageParameters model_parameters = readStandardTwoStageParameters(model_path);
    // read sto file
    secondStageRHSmap RHSmap = readStochasticMap(sto_path);
    // STEP 2: SOLVING PROCESS (SD-kNN)
    // initialize feasibility cut collection
    std::vector<feasibilityCut> feasibility_cuts;
    // initialization of output file
    const char* writeFilePath = resultsOutput_path.c_str();
    std::fstream writeFile;
    writeFile.open(writeFilePath,std::fstream::app); // append results to the end of the file
    //
    // write initial setup
    std::cout << "*******************************************\n";
    writeFile << "*******************************************\n";
    std::cout << "LEONA kNN(v1.0) is initialized\n";
    writeFile << "LEONA kNN(v1.0) is initialized\n";
    std::cout << "Algorithmic Parameters\n";
    writeFile << "Algorithmic Parameters\n";
    std::cout << "Dx, m, max_subgradient, beta (kNN), N_pre, N_incre\n";
    writeFile << "Dx, m, max_subgradient, beta (kNN), N_pre, N_incre\n";
    std::cout << Dx << ", " << m << ", " << max_subgradient << ", " << beta << ", " << N_pre << ", " << N_incre << std::endl;
    writeFile << Dx << ", " << m << ", " << max_subgradient << ", " << beta << ", " << N_pre << ", " << N_incre << std::endl;
    writeFile << "Initial solution: ";
    long x_init_size = x_init.size();
    for (int index = 0; index < x_init_size - 1; ++index) {
        writeFile << x_init[index] << ", ";
    }
    writeFile << x_init[x_init_size - 1] << "\n";
    writeFile << "Max number of iterations: " << max_it << "\n";
    std::cout << "Problem Complexity\n";
    writeFile << "Problem Complexity\n";
    std::cout << "A_num_row, A_num_col\n";
    writeFile << "A_num_row, A_num_col\n";
    std::cout << model_parameters.A.getRowLength() << ", " << model_parameters.A.getRowLength() << std::endl;
    writeFile << model_parameters.A.getRowLength() << ", " << model_parameters.A.getRowLength() << std::endl;
    std::cout << "D_num_row, D_num_col (after converting into standard form)\n";
    writeFile << "D_num_row, D_num_col (after converting into standard form)\n";
    std::cout << model_parameters.D.getRowLength() << ", " << model_parameters.D.getColLength() << std::endl;
    writeFile << model_parameters.D.getRowLength() << ", " << model_parameters.D.getColLength() << std::endl;
    // set up initial incumbent solution
    //long A_rowsize = model_parameters.A.getRowLength();
    long A_colsize = model_parameters.A.getColLength();
    // output observed predictor
    std::cout << "Observed Predictor: ";
    writeFile << "Observed Predictor: ";
    for (int predictor_index = 0; predictor_index < observed_predictor.size() - 1; ++predictor_index) {
        std::cout << observed_predictor[predictor_index] << ", ";
        writeFile << observed_predictor[predictor_index] << ", ";
    }
    std::cout << observed_predictor[observed_predictor.size() - 1] << std::endl;
    writeFile << observed_predictor[observed_predictor.size() - 1] << std::endl;
    // initial variable
    std::vector<double> x_new(A_colsize,0.0);
    std::vector<double> x_leon(A_colsize,0.0);
    // projection
    std::vector<double> x_old = leona_projection(x_init, model_parameters);
    for (int x_index = 0; x_index < A_colsize; ++x_index) {
        x_new[x_index] = x_old[x_index];
        x_leon[x_index] = x_old[x_index];
    }
    // set up subproblem environment
    // set up the model
    IloEnv env_sub;
    IloModel mod_sub(env_sub);
    IloNumVarArray y(env_sub,model_parameters.D.getColLength(),model_parameters.y_lb,IloInfinity,ILOFLOAT);
    mod_sub.add(y);
    IloExpr expr_obj_sub(env_sub);
    for (int idx = 0; idx < model_parameters.d.getNzeroLen(); ++idx) {
        expr_obj_sub += model_parameters.d.getVal(idx) * y[model_parameters.d.getLoc(idx)];
    }
    IloObjective obj_sub = IloMinimize(env_sub,expr_obj_sub);
    mod_sub.add(obj_sub);
    // stndard form equality constraints [D Islack] [y s] + Cx = e
    IloRangeArray constraintsEquality_sub(env_sub);
    std::vector<IloExpr> exprs_eq_sub;
    for (int idx = 0; idx < model_parameters.D.getRowLength(); ++idx) {
        IloExpr expr(env_sub);
        exprs_eq_sub.push_back(expr);
    }
    // coefficients before y; Dy
    for (int col_idx = 0; col_idx < model_parameters.D.getColLength(); ++col_idx) {
        int beg_idx = model_parameters.D.getCbeg(col_idx);
        for (int idx = beg_idx; idx < model_parameters.D.getClen(col_idx) + beg_idx; ++idx) {
            exprs_eq_sub[model_parameters.D.getRow(idx)] += model_parameters.D.getVal(idx) * y[col_idx];
        }
    }
    // coefficients before x; Cx
    for (int col_idx = 0; col_idx < model_parameters.C.getColLength(); ++col_idx) {
        int beg_idx = model_parameters.C.getCbeg(col_idx);
        for (int idx = beg_idx; idx < model_parameters.C.getClen(col_idx) + beg_idx; ++idx) {
            exprs_eq_sub[model_parameters.C.getRow(idx)] += model_parameters.C.getVal(idx) * x_new[col_idx];
        }
    }
    // right hand side e
    for (int idx = 0; idx < model_parameters.e.getNzeroLen(); ++idx) {
        exprs_eq_sub[model_parameters.e.getLoc(idx)] -= model_parameters.e.getVal(idx);
    }
    // add the equality constraints
    for (int idx = 0; idx< model_parameters.D.getRowLength(); ++idx) {
        constraintsEquality_sub.add(exprs_eq_sub[idx] == 0);
    }
    mod_sub.add(constraintsEquality_sub);
    // set up cplex solver for the subproblem
    IloCplex cplex_sub(env_sub);
    cplex_sub.extract(mod_sub);
    cplex_sub.setOut(env_sub.getNullStream());
    // intermediate values for rhs update
    // deterministic part of the rhs_bounds
    std::vector<double> e_det(model_parameters.D.getRowLength(), 0.0);
    for (int idx = 0; idx < model_parameters.e.getNzeroLen(); ++idx) {
        e_det[model_parameters.e.getLoc(idx)] += model_parameters.e.getVal(idx);
    }
    std::vector<double> det_rhs_bounds = e_det;
    std::vector<double> rhs_bounds(model_parameters.D.getRowLength(), 0.0);
    // end settig up subproblem ==================================
    // main solving process
    // write main loop
    writeFile << "Main loop\n";
    // data start index
    int data_start_idx= 0;
    // sample size for NSQG update
    int N_nsqg = N_pre;
    int N_total = 0;
    // main loop (outer loop)
    for (int outer_it = 0; outer_it < max_it; ++outer_it) {
        // initialization reset x_leon
        for (int x_index = 0; x_index < A_colsize; ++x_index) {
            x_leon[x_index] = 0;
        }
        // max number iterates in the inner loop
        int maxInner_iterates = (outer_it + 1) * m;
        std::cout << "Outer Loop: Iteration " << outer_it + 1 << std::endl;
        writeFile << "##################################################\n";
        writeFile << "Outer Loop: Iteration " << outer_it + 1 << "\n";
        // step size
        double stepsize = 2.0 * Dx / (max_subgradient * sqrt((double) maxInner_iterates));
        // inner loop NSQG update
        for (int inner_it = 0; inner_it < maxInner_iterates; ++inner_it) {
            // kNN estimation
            // clear RHS dataset
            RHS_dataset.clear();
            // clear distance set
            distanceSet.clear();
            // clear order set
            orderSet.clear();
            k_new = (int) pow(N_nsqg, beta);
            for (int idx = 0; idx < N_nsqg; ++idx) {
                // obtain data points
                int idx_dataPoint = data_start_idx + idx;
                // obtain a new data point
                dataPoint be_datapoint;
                if (flag_be == true) {
                    be_datapoint = be_DB[0][idx_dataPoint];
                }
                dataPoint bi_datapoint;
                if (flag_bi == true) {
                    bi_datapoint = bi_DB[0][idx_dataPoint];
                }
                dataPoint Ce_datapoint;
                if (flag_Ce == true) {
                    Ce_datapoint = Ce_DB[0][idx_dataPoint];
                }
                dataPoint Ci_datapoint;
                if (flag_Ci == true) {
                    Ci_datapoint = Ci_DB[0][idx_dataPoint];
                }
                // merge all the datapoints
                secondStageRHSpoint RHS_datapoint = merge_randomVector(be_datapoint, bi_datapoint, Ce_datapoint, Ci_datapoint);
                RHS_dataset.push_back(RHS_datapoint);
                // calculate the squared distance
                double distance_squared = 0;
                for (int predictor_idx = 0; predictor_idx < RHS_datapoint.predictor.size(); ++predictor_idx) {
                    distance_squared += (RHS_datapoint.predictor[predictor_idx] - observed_predictor[predictor_idx]) * (RHS_datapoint.predictor[predictor_idx] - observed_predictor[predictor_idx]);
                }
                distanceSet.push_back(distance_squared);
                // store the new squared distance
                // sorting (like insert sorting)
                if (idx_dataPoint == 0) { // first iteration
                    orderSet.push_back(1);
                }
                else {// from left to right in increasing order
                    //int left_index = 0; // the index corresponds to the largest distance that is smaller than the current one
                    //double left_distance = -1;
                    // double indices used for tie-breaking
                    int right_index = -1; // the index corresponds to the smallest distance that is larger than  the current one
                    double right_distance = -1;
                    for (int index = 0; index < orderSet.size(); ++index) {
                        /*
                        if (distanceSet[index] < distance_squared) {
                            if (left_index == 0) {
                                left_distance = distanceSet[index];
                                left_index = orderSet[index];
                            }
                            else if (distanceSet[index] > left_distance) {
                                left_distance = distanceSet[index];
                                left_index = orderSet[index];
                            }
                        }
                         */
                        if (distanceSet[index] > distance_squared) {
                            if (right_index == -1) {
                                right_distance = distanceSet[index];
                                right_index = orderSet[index];
                            }
                            else if (distanceSet[index] < right_distance) {
                                right_distance = distanceSet[index];
                                right_index = orderSet[index];
                            }
                            else if (distanceSet[index] == right_distance && right_index > orderSet[index]) {
                                right_index = orderSet[index];
                            }
                        }
                    }
                    /*
                    if (flag_debug == true) {
                        std::cout << "Output double indices\n";
                        writeFile << "Output double indices\n";
                        std::cout << "left index: " << left_index << std::endl;
                        writeFile << "left index: " << left_index << std::endl;
                        std::cout << "right index: " << right_index << std::endl;
                        writeFile << "right index: " << right_index << std::endl;
                    }
                     */
                    // update the orderSet
                    for (int index = 0; index < orderSet.size(); ++index) {
                        if (right_index != -1 && orderSet[index] >= right_index) {
                            orderSet[index] = orderSet[index] + 1;
                        }
                        //if (left_index == 0) { // current one is the nearest neighbor
                        //    orderSet[index] = orderSet[index] + 1;
                        //}
                        //else if (orderSet[index] > left_index) {
                        //    orderSet[index] = orderSet[index] + 1;
                        //}
                    }
                    if (right_index == -1) {
                        orderSet.push_back((int) orderSet.size() + 1);
                    }
                    else {
                        orderSet.push_back(right_index);
                    }
                    /*
                    if (flag_debug == true) {
                        std::cout << "Updated Order in the scenario set\n";
                        writeFile << "Updated Order in the scenario set\n";
                        std::cout << "Index, Order, Distance (Squared)\n";
                        writeFile << "Index, Order, Distance (Squared)\n";
                        // update the kNN set
                        for (int index = 0; index < orderSet.size(); ++index) {
                            std::cout << index << ", "<< orderSet[index] << ", " << distanceSet[index];
                            writeFile << index << ", "<< orderSet[index] << ", " << distanceSet[index];
                            if (orderSet[index] <= k_new) {
                                std::cout << "*";
                                writeFile << "*";
                                kNNSet_new.push_back(index);
                            }
                            std::cout << std::endl;
                            writeFile << std::endl;
                        }
                    }
                    else {
                        // update the kNN set
                        for (int index = 0; index < orderSet.size(); ++index) {
                            if (orderSet[index] <= k_new) {
                                kNNSet_new.push_back(index);
                            }
                        }
                    }
                     */
                }
            }
            // clear kNN set
            kNNSet.clear();
            for (int index = 0; index < orderSet.size(); ++index) {
                if (orderSet[index] <= k_new) {
                    kNNSet.push_back(index);
                    if (orderSet[index] == k_new) {
                        knn_radius = distanceSet[index];
                    }
                }
            }
            // end kNN estimation
            // compute NSQG
            std::vector<double> nsqg(A_colsize,0.0);
            double one_over_k = 1.0 / ((double) k_new);
            det_rhs_bounds = e_det;
            // det C x
            for (int col_idx = 0; col_idx < model_parameters.C.getColLength(); ++col_idx) {
                int beg_idx = model_parameters.C.getCbeg(col_idx);
                for (int idx = beg_idx; idx < model_parameters.C.getClen(col_idx) + beg_idx; ++idx) {
                    det_rhs_bounds[model_parameters.C.getRow(idx)] -= model_parameters.C.getVal(idx) * x_old[col_idx];
                }
            }
            for (int knn_idx = 0; knn_idx < k_new; ++knn_idx) {
                // update the subproblem
                for (int idx_row = 0; idx_row < model_parameters.D.getRowLength(); ++idx_row) {
                    rhs_bounds[idx_row] = det_rhs_bounds[idx_row];
                } // update the deterministic part
                // update the stochastic parts of e
                for (int idx_be = 0; idx_be < RHS_dataset[kNNSet[knn_idx]].be.size(); ++idx_be) {
                    rhs_bounds[RHSmap.be_map[idx_be]] += RHS_dataset[kNNSet[knn_idx]].be[idx_be];
                }
                // right hand side (stochastic part) equality bi_(i) inequality
                for (int idx_bi = 0; idx_bi < RHS_dataset[kNNSet[knn_idx]].bi.size(); ++idx_bi) {
                    rhs_bounds[RHSmap.bi_map[idx_bi] + model_parameters.num_eq] += RHS_dataset[kNNSet[knn_idx]].bi[idx_bi];
                }
                // coefficients before x (stochastic part) equality (i.e., Cij * xj map: <i,j> )
                for (int idx_Ce = 0; idx_Ce < RHS_dataset[kNNSet[knn_idx]].Ce.size(); ++idx_Ce) {
                    rhs_bounds[RHSmap.Ce_map[idx_Ce].first] -= RHS_dataset[kNNSet[knn_idx]].Ce[idx_Ce] * x_old[RHSmap.Ce_map[idx_Ce].second];
                }
                // coefficients before x (stochastic part) inequality (location is behind equality constraints)
                for (int idx_Ci = 0; idx_Ci < RHS_dataset[kNNSet[knn_idx]].Ci.size(); ++idx_Ci) {
                    rhs_bounds[RHSmap.Ci_map[idx_Ci].first + model_parameters.num_eq] -= RHS_dataset[kNNSet[knn_idx]].Ci[idx_Ci] * x_old[RHSmap.Ci_map[idx_Ci].second];
                }
                
                // update the RHS
                for (int idx_row = 0; idx_row < rhs_bounds.size(); ++idx_row) {
                    constraintsEquality_sub[idx_row].setBounds(rhs_bounds[idx_row], rhs_bounds[idx_row]);
                    // reset rhs_bounds to 0
                    rhs_bounds[idx_row] = 0;
                }
                // calculate the dual multipliers
                IloBool flag_solve = cplex_sub.solve();
                if (flag_solve == IloTrue) {
                    IloNumArray dual_equality_sub(env_sub);
                    //double optimal_value = cplex_sub.getObjValue(); // get the optimal value
                    cplex_sub.getDuals(dual_equality_sub,constraintsEquality_sub);
                    std::vector<double> minus_dual_sub;
                    for (int dual_idx = 0; dual_idx < model_parameters.D.getRowLength(); ++dual_idx) {
                        minus_dual_sub.push_back((-1.0) * dual_equality_sub[dual_idx]);
                    }
                    // compute NSQG
                    std::vector<double> subgradient = model_parameters.C.fast_rightMultiply(minus_dual_sub);
                    // stochastic C
                    if (flag_Ce == true || flag_Ci == true) {
                        // equality
                        for (int idx_Ce = 0; idx_Ce < RHS_dataset[kNNSet[knn_idx]].Ce.size(); ++idx_Ce) {
                            subgradient[RHSmap.Ce_map[idx_Ce].second] += RHS_dataset[kNNSet[knn_idx]].Ce[idx_Ce] * minus_dual_sub[RHSmap.Ce_map[idx_Ce].first];
                        }
                        // inequality
                        for (int idx_Ci = 0; idx_Ci < RHS_dataset[kNNSet[knn_idx]].Ci.size(); ++idx_Ci) {
                            subgradient[RHSmap.Ci_map[idx_Ci].second] += RHS_dataset[kNNSet[knn_idx]].Ci[idx_Ci] * minus_dual_sub[RHSmap.Ci_map[idx_Ci].first + model_parameters.num_eq];
                        }
                    }
                    for (int idx = 0; idx < A_colsize; ++idx) {
                        nsqg[idx] += subgradient[idx];
                    }
                }
                else {
                    throw std::logic_error("Main solving process terminate. Subproblem is infeasible.\n");
                }
            } // end for (int knn_idx = 0; knn_idx < k_new; ++knn_idx)
            // averaging
            for (int idx = 0; idx < A_colsize; ++idx) {
                nsqg[idx] = nsqg[idx] * one_over_k;
            }
            // first stage
            for (int idx = 0; idx < model_parameters.c.getNzeroLen(); ++idx) {
                nsqg[model_parameters.c.getLoc(idx)] += model_parameters.c.getVal(idx);
            }
            // update solution
            for (int idx = 0; idx < A_colsize; ++idx) {
                x_new[idx] = x_old[idx] - nsqg[idx] * stepsize;
            }
            writeFile << "------------------------------------------------\n";
            writeFile << "Inner loop: Iteration  " << inner_it + 1 << "\n";
            writeFile << "Number of scenarios: " << N_nsqg << "\n";
            // write stepzie
            writeFile << "Stepsize: " << stepsize << std::endl;
            // write quasigradient
            writeFile << "NSQG: ";
            for (int idx = 0; idx < A_colsize - 1; ++idx) {
                writeFile << nsqg[idx] << ", ";
            }
            writeFile << nsqg[A_colsize - 1] << std::endl;
            // projection
            x_new = leona_projection(x_new, model_parameters);
            // update x_leon
            for (int idx = 0; idx < A_colsize; ++idx) {
                x_leon[idx] += x_new[idx];
            }
            // update x_old
            std::cout << "x (nsqg update): (";
            writeFile << "x (nsqg update): (";
            // update x_old and go to the next iteration
            for (int x_index = 0; x_index < A_colsize; ++x_index) {
                x_old[x_index] = x_new[x_index];
                std::cout << x_new[x_index] << " ";
                writeFile << x_new[x_index] << " ";
            }
            std::cout << ")" <<std::endl;
            writeFile << ")" <<std::endl;
            std::cout << "**************************" << std::endl;
            writeFile << "**************************" << std::endl;
            // update sample size and data start index
            N_total += N_nsqg;
            data_start_idx += N_nsqg;
            N_nsqg += N_incre;
        } // end inner loop NSQG updates
        double mq = maxInner_iterates;
        // averaging to get x_leon
        for (int idx = 0; idx < A_colsize; ++idx) {
            x_leon[idx] = x_leon[idx] / mq;
        }
    } // end main outer loop
    
    std::cout << "*******************************************\n";
    writeFile << "*******************************************\n";
    std::cout << "Total Number of Samples Used: " << N_total << std::endl;
    writeFile << "Total Number of Samples Used: " << N_total << std::endl;
    std::cout << "Output Solution: ";
    writeFile << "Output Solution: ";
    for (int index = 0; index < A_colsize-1; ++index) {
        std::cout << x_leon[index] << ", ";
        writeFile << x_leon[index] << ", ";
    }
    std::cout << x_leon[A_colsize - 1] << std::endl;
    writeFile << x_leon[A_colsize - 1] << std::endl;
    std::cout << "Computation Log: Finish Solving Process.\n";
    writeFile << "Computation Log: Finish Solving Process.\n";
    // write time elapsed
    double duration = (std::clock() - time_start ) / (double) CLOCKS_PER_SEC;
    writeFile << "Time elapsed(secs) : " << duration << "\n";
    writeFile << "*******************************************\n";
    
    // return results
    solverOutput res;
    res.x = x_leon;
    res.time_elapse = duration;
    res.num_it = max_it;
    res.N_total = N_total;
    return res;
}


// LEONA kernel solver 1: Naive Kernel, 2: Epanechnikov Kernel, 3: Quartic Kernel, 4: Gaussian Kernel
solverOutput leona_kernel(const std::string& folder_path,
                          int max_it,
                          const std::vector<double>& x_init,
                          const std::vector<double>& observed_predictor,
                          int N_pre,
                          int N_incre,
                          int flag_kernel,
                          double Dx,
                          int m,
                          double max_subgradient) {
    // timer
    std::clock_t time_start;
    time_start = std::clock();
    // current time
    //std::time_t currTime = std::time(nullptr);
    // STEP 1: INITIALIZATION
    long predictor_dimension = observed_predictor.size();
    //double C = 5; // bk19 kernel regression is sensitive to the bandwidth
    double C = 1;// bk19 for Gaussian
    //double C = 50; // pbaa99
    //double C = 10; // pbaa99 for Gaussian
    double beta = 0.5 / ((double) predictor_dimension); // beta is between 0 and 1/predictor_dimension
    bool flag_be; // tell if be stochastic is generated
    bool flag_bi; // tell if bi stochastic is generated
    bool flag_Ce; // tell if Ce stochastic is generated
    bool flag_Ci; // tell if Ci stochastic is generated
    std::vector<secondStageRHSpoint> RHS_dataset;
    // create directory paths for database and model
    std::string be_DB_path = folder_path + "/be_DB.txt";
    std::string bi_DB_path = folder_path + "/bi_DB.txt";
    std::string Ce_DB_path = folder_path + "/Ce_DB.txt";
    std::string Ci_DB_path = folder_path + "/Ci_DB.txt";
    std::string model_path = folder_path + "/model.txt";
    std::string sto_path = folder_path + "/sto.txt";
    std::string resultsOutput_path = folder_path + "/computationalResults(leona_kernel1.0).txt";
    // convert all the paths into constant chars
    const char* be_DB_path_const = be_DB_path.c_str();
    const char* bi_DB_path_const = bi_DB_path.c_str();
    const char* Ce_DB_path_const = Ce_DB_path.c_str();
    const char* Ci_DB_path_const = Ci_DB_path.c_str();
    // create stream object
    std::ifstream readFile_be(be_DB_path_const);
    std::ifstream readFile_bi(bi_DB_path_const);
    std::ifstream readFile_Ce(Ce_DB_path_const);
    std::ifstream readFile_Ci(Ci_DB_path_const);
    // create database
    std::vector<std::vector<dataPoint>> be_DB;
    std::vector<std::vector<dataPoint>> bi_DB;
    std::vector<std::vector<dataPoint>> Ce_DB;
    std::vector<std::vector<dataPoint>> Ci_DB;
    // read  be
    if (readFile_be.is_open()) {
        std::cout << "be_DB is found." << std::endl;
        readFile_be.close(); // close the file
        // read be database
        be_DB = readNonparametricDB(be_DB_path);
        flag_be = true;
    }
    else {
        readFile_be.close(); // close the file
        flag_be = false;
        std::cout << "be_DB is not found!" << std::endl;
    }
    // read bi
    if (readFile_bi.is_open()) {
        std::cout << "bi_DB is found." << std::endl;
        readFile_be.close(); // close the file
        // read bi database
        bi_DB = readNonparametricDB(bi_DB_path);
        flag_bi = true;
    }
    else {
        readFile_bi.close(); // close the file
        flag_bi = false;
        std::cout << "bi_DB is not found!" << std::endl;
    }
    // read Ce
    if (readFile_Ce.is_open()) {
        std::cout << "Ce_DB stochastic part is found." << std::endl;
        readFile_Ce.close(); // close the file
        // Ce database
        Ce_DB = readNonparametricDB(Ce_DB_path);
        flag_Ce = true;
    }
    else {
        readFile_Ce.close(); // close the file
        flag_Ce = false;
        std::cout << "Ce_DB is not found!" << std::endl;
    }
    // read Ci
    if (readFile_Ci.is_open()) {
        std::cout << "Ci_DB stochastic part is found." << std::endl;
        readFile_Ci.close(); // close the file
        // Ci database
        Ci_DB = readNonparametricDB(Ci_DB_path);
        flag_Ci = true;
    }
    else {
        readFile_Ci.close(); // close the file
        flag_Ci = false;
        std::cout << "Ci_DB is not found!" << std::endl;
    }
    // read model file
    // create model structure
    standardTwoStageParameters model_parameters = readStandardTwoStageParameters(model_path);
    // read sto file
    secondStageRHSmap RHSmap = readStochasticMap(sto_path);
    // STEP 2: SOLVING PROCESS (SD-kNN)
    // initialize feasibility cut collection
    std::vector<feasibilityCut> feasibility_cuts;
    // initialization of output file
    const char* writeFilePath = resultsOutput_path.c_str();
    std::fstream writeFile;
    writeFile.open(writeFilePath,std::fstream::app); // append results to the end of the file
    // write initial setup
    std::cout << "*******************************************\n";
    writeFile << "*******************************************\n";
    std::cout << "LEONA kernel(v1.0) is initialized\n";
    writeFile << "LEONA kernel(v1.0) is initialized\n";
    std::cout << "Algorithmic Parameters\n";
    writeFile << "Algorithmic Parameters\n";
    std::cout << "Dx, m, max_subgradient, beta (kernel), C, N_pre, N_incre\n";
    writeFile << "Dx, m, max_subgradient, beta (kernel), C, N_pre, N_incre\n";
    std::cout << Dx << ", " << m << ", " << max_subgradient << ", " << beta << ", " << C << ", " << N_pre << ", " << N_incre << std::endl;
    writeFile << Dx << ", " << m << ", " << max_subgradient << ", " << beta << ", " << C << ", " << N_pre << ", " << N_incre << std::endl;
    switch (flag_kernel) {
        case 1:
            std::cout << "Naive kernel is used\n";
            writeFile << "Naive kernel is used\n";
            break;
        case 2:
            std::cout << "Epanechnikov kernel is used.\n";
            writeFile << "Epanechnikov kernel is used.\n";
            break;
        case 3:
            std::cout << "Quartic kernel is used.\n";
            writeFile << "Quartic kernel is used.\n";
            break;
        case 4:
            std::cout << "Gaussian kernel is used.\n";
            writeFile << "Gaussian kernel is used.\n";
            break;
        default:
            throw std::logic_error("Kernel flag is incorrect! Solver is terminated.\n");
            break;
    }
    writeFile << "Initial solution: ";
    long x_init_size = x_init.size();
    for (int index = 0; index < x_init_size - 1; ++index) {
        writeFile << x_init[index] << ", ";
    }
    writeFile << x_init[x_init_size - 1] << "\n";
    writeFile << "Max number of iterations: " << max_it << "\n";
    std::cout << "Problem Complexity\n";
    writeFile << "Problem Complexity\n";
    std::cout << "A_num_row, A_num_col\n";
    writeFile << "A_num_row, A_num_col\n";
    std::cout << model_parameters.A.getRowLength() << ", " << model_parameters.A.getRowLength() << std::endl;
    writeFile << model_parameters.A.getRowLength() << ", " << model_parameters.A.getRowLength() << std::endl;
    std::cout << "D_num_row, D_num_col (after converting into standard form)\n";
    writeFile << "D_num_row, D_num_col (after converting into standard form)\n";
    std::cout << model_parameters.D.getRowLength() << ", " << model_parameters.D.getColLength() << std::endl;
    writeFile << model_parameters.D.getRowLength() << ", " << model_parameters.D.getColLength() << std::endl;
    // set up initial incumbent solution
    //long A_rowsize = model_parameters.A.getRowLength();
    long A_colsize = model_parameters.A.getColLength();
    // output observed predictor
    std::cout << "Observed Predictor: ";
    writeFile << "Observed Predictor: ";
    for (int predictor_index = 0; predictor_index < observed_predictor.size() - 1; ++predictor_index) {
        std::cout << observed_predictor[predictor_index] << ", ";
        writeFile << observed_predictor[predictor_index] << ", ";
    }
    std::cout << observed_predictor[observed_predictor.size() - 1] << std::endl;
    writeFile << observed_predictor[observed_predictor.size() - 1] << std::endl;
    // initial variable
    std::vector<double> x_new(A_colsize,0.0);
    std::vector<double> x_leon(A_colsize,0.0);
    // projection
    std::vector<double> x_old = leona_projection(x_init, model_parameters);
    for (int x_index = 0; x_index < A_colsize; ++x_index) {
        x_new[x_index] = x_old[x_index];
        x_leon[x_index] = x_old[x_index];
    }
    // set up subproblem environment
    // set up the model
    IloEnv env_sub;
    IloModel mod_sub(env_sub);
    IloNumVarArray y(env_sub,model_parameters.D.getColLength(),model_parameters.y_lb,IloInfinity,ILOFLOAT);
    mod_sub.add(y);
    IloExpr expr_obj_sub(env_sub);
    for (int idx = 0; idx < model_parameters.d.getNzeroLen(); ++idx) {
        expr_obj_sub += model_parameters.d.getVal(idx) * y[model_parameters.d.getLoc(idx)];
    }
    IloObjective obj_sub = IloMinimize(env_sub,expr_obj_sub);
    mod_sub.add(obj_sub);
    // stndard form equality constraints [D Islack] [y s] + Cx = e
    IloRangeArray constraintsEquality_sub(env_sub);
    std::vector<IloExpr> exprs_eq_sub;
    for (int idx = 0; idx < model_parameters.D.getRowLength(); ++idx) {
        IloExpr expr(env_sub);
        exprs_eq_sub.push_back(expr);
    }
    // coefficients before y; Dy
    for (int col_idx = 0; col_idx < model_parameters.D.getColLength(); ++col_idx) {
        int beg_idx = model_parameters.D.getCbeg(col_idx);
        for (int idx = beg_idx; idx < model_parameters.D.getClen(col_idx) + beg_idx; ++idx) {
            exprs_eq_sub[model_parameters.D.getRow(idx)] += model_parameters.D.getVal(idx) * y[col_idx];
        }
    }
    // coefficients before x; Cx
    for (int col_idx = 0; col_idx < model_parameters.C.getColLength(); ++col_idx) {
        int beg_idx = model_parameters.C.getCbeg(col_idx);
        for (int idx = beg_idx; idx < model_parameters.C.getClen(col_idx) + beg_idx; ++idx) {
            exprs_eq_sub[model_parameters.C.getRow(idx)] += model_parameters.C.getVal(idx) * x_new[col_idx];
        }
    }
    // right hand side e
    for (int idx = 0; idx < model_parameters.e.getNzeroLen(); ++idx) {
        exprs_eq_sub[model_parameters.e.getLoc(idx)] -= model_parameters.e.getVal(idx);
    }
    // add the equality constraints
    for (int idx = 0; idx< model_parameters.D.getRowLength(); ++idx) {
        constraintsEquality_sub.add(exprs_eq_sub[idx] == 0);
    }
    mod_sub.add(constraintsEquality_sub);
    // set up cplex solver for the subproblem
    IloCplex cplex_sub(env_sub);
    cplex_sub.extract(mod_sub);
    cplex_sub.setOut(env_sub.getNullStream());
    // intermediate values for rhs update
    // deterministic part of the rhs_bounds
    std::vector<double> e_det(model_parameters.D.getRowLength(), 0.0);
    for (int idx = 0; idx < model_parameters.e.getNzeroLen(); ++idx) {
        e_det[model_parameters.e.getLoc(idx)] += model_parameters.e.getVal(idx);
    }
    std::vector<double> det_rhs_bounds = e_det;
    std::vector<double> rhs_bounds(model_parameters.D.getRowLength(), 0.0);
    // end settig up subproblem ==================================
    // main solving process
    // write main loop
    writeFile << "Main loop\n";
    // data start index
    int data_start_idx= 0;
    // sample size for NSQG update
    int N_total = 0;
    int N_nsqg = N_pre;
    // main loop (outer loop)
    for (int outer_it = 0; outer_it < max_it; ++outer_it) {
        // initialization reset x_leon
        for (int x_index = 0; x_index < A_colsize; ++x_index) {
            x_leon[x_index] = 0;
        }
        // max number iterates in the inner loop
        int maxInner_iterates = (outer_it + 1) * m;
        std::cout << "Outer Loop: Iteration " << outer_it + 1 << std::endl;
        writeFile << "##################################################\n";
        writeFile << "Outer Loop: Iteration " << outer_it + 1 << "\n";
        // step size
        double stepsize = 2.0 * Dx / (max_subgradient * sqrt((double) maxInner_iterates));
        for (int inner_it = 0; inner_it < maxInner_iterates; ++inner_it) {
            // clear RHS dataset
            RHS_dataset.clear();
            // bandwidth
            double hN = C * pow(N_nsqg, -beta);
            std::vector<double> kernel_weight;
            double total_weight = 0;
            std::cout << "Bandwidth: " << hN << std::endl;
            writeFile << "Bandwidth: " << hN << std::endl;
            for (int idx = 0; idx < N_nsqg; ++idx) {
                // obtain data points
                int idx_dataPoint = data_start_idx + idx;
                // obtain a new data point
                dataPoint be_datapoint;
                if (flag_be == true) {
                    be_datapoint = be_DB[0][idx_dataPoint];
                }
                dataPoint bi_datapoint;
                if (flag_bi == true) {
                    bi_datapoint = bi_DB[0][idx_dataPoint];
                }
                dataPoint Ce_datapoint;
                if (flag_Ce == true) {
                    Ce_datapoint = Ce_DB[0][idx_dataPoint];
                }
                dataPoint Ci_datapoint;
                if (flag_Ci == true) {
                    Ci_datapoint = Ci_DB[0][idx_dataPoint];
                }
                // merge all the datapoints
                secondStageRHSpoint RHS_datapoint = merge_randomVector(be_datapoint, bi_datapoint, Ce_datapoint, Ci_datapoint);
                RHS_dataset.push_back(RHS_datapoint); // 
                // calculate the (w - w_obs) / h_N
                double distance_squared = 0;
                for (int predictor_idx = 0; predictor_idx < RHS_datapoint.predictor.size(); ++predictor_idx) {
                    distance_squared += (RHS_datapoint.predictor[predictor_idx] - observed_predictor[predictor_idx]) * (RHS_datapoint.predictor[predictor_idx] - observed_predictor[predictor_idx]);
                }
                // calculate weight
                // distance squared normalized by the bandwidth
                double normalized_distance_squared = distance_squared / hN / hN;
                switch (flag_kernel) {
                    case 1: // Naive kernel
                        if (normalized_distance_squared <= 1.0) {
                            kernel_weight.push_back(1.0);
                        }
                        else {
                            kernel_weight.push_back(0.0);
                        }
                        break;
                    case 2: // Epanechnikov kernel
                        if (normalized_distance_squared <= 1.0) {
                            double tmp_val = 1 - normalized_distance_squared;
                            kernel_weight.push_back(tmp_val);
                        }
                        else {
                            kernel_weight.push_back(0.0);
                        }
                        break;
                    case 3: // Quartic kernel
                        if (normalized_distance_squared <= 1.0) {
                            double tmp_val = 1 - normalized_distance_squared;
                            kernel_weight.push_back(tmp_val * tmp_val);
                        }
                        else {
                            kernel_weight.push_back(0.0);
                        }
                        break;
                    case 4: // Gaussian kernel
                        kernel_weight.push_back(exp((-0.5) * normalized_distance_squared));
                        //std::cout << "normalized_distance_squared = " << normalized_distance_squared << std::endl;
                        //std::cout << "e(-0.5 * normalized_distance_squared) = " << kernel_weight[idx] << std::endl;
                        break;
                    default:
                        break;
                }
                total_weight += kernel_weight[idx];
            } // end for (int idx = 0; idx < N_nsqg; ++idx)
            // compute NSQG
            std::vector<double> nsqg(A_colsize,0.0); // initialize nsqg
            det_rhs_bounds = e_det;
            // det C x
            for (int col_idx = 0; col_idx < model_parameters.C.getColLength(); ++col_idx) {
                int beg_idx = model_parameters.C.getCbeg(col_idx);
                for (int row_idx = beg_idx; row_idx < model_parameters.C.getClen(col_idx) + beg_idx; ++row_idx) {
                    det_rhs_bounds[model_parameters.C.getRow(row_idx)] -= model_parameters.C.getVal(row_idx) * x_old[col_idx];
                }
            }
            if (total_weight > SOLVER_PRECISION_UPPER) { // total weight is greater than 0
                for (int idx = 0; idx < N_nsqg; ++idx) {
                    // compute normalized weight
                    double normalized_weight = kernel_weight[idx] / total_weight;
                    //std::cout << "Debug: idx = " << idx << ", normalized_weight = " << normalized_weight << std::endl;
                    if (normalized_weight > SOLVER_PRECISION_UPPER) { // normalized_weight is significantly greater than 0
                        // update the subproblem
                        for (int idx_row = 0; idx_row < model_parameters.D.getRowLength(); ++idx_row) {
                            rhs_bounds[idx_row] = det_rhs_bounds[idx_row];
                        } // update the deterministic part
                        // update the stochastic parts of e
                        for (int idx_be = 0; idx_be < RHS_dataset[idx].be.size(); ++idx_be) {
                            rhs_bounds[RHSmap.be_map[idx_be]] += RHS_dataset[idx].be[idx_be];
                        }
                        // right hand side (stochastic part) equality bi_(i) inequality
                        for (int idx_bi = 0; idx_bi < RHS_dataset[idx].bi.size(); ++idx_bi) {
                            rhs_bounds[RHSmap.bi_map[idx_bi] + model_parameters.num_eq] += RHS_dataset[idx].bi[idx_bi];
                        }
                        // coefficients before x (stochastic part) equality (i.e., Cij * xj map: <i,j> )
                        for (int idx_Ce = 0; idx_Ce < RHS_dataset[idx].Ce.size(); ++idx_Ce) {
                            rhs_bounds[RHSmap.Ce_map[idx_Ce].first] -= RHS_dataset[idx].Ce[idx_Ce] * x_old[RHSmap.Ce_map[idx_Ce].second];
                        }
                        // coefficients before x (stochastic part) inequality (location is behind equality constraints)
                        for (int idx_Ci = 0; idx_Ci < RHS_dataset[idx].Ci.size(); ++idx_Ci) {
                            rhs_bounds[RHSmap.Ci_map[idx_Ci].first + model_parameters.num_eq] -= RHS_dataset[idx].Ci[idx_Ci] * x_old[RHSmap.Ci_map[idx_Ci].second];
                        }
                        
                        // update the RHS
                        for (int idx_row = 0; idx_row < rhs_bounds.size(); ++idx_row) {
                            constraintsEquality_sub[idx_row].setBounds(rhs_bounds[idx_row], rhs_bounds[idx_row]);
                            // reset rhs_bounds to 0
                            rhs_bounds[idx_row] = 0;
                        }
                        // calculate the dual multipliers
                        IloBool flag_solve = cplex_sub.solve();
                        if (flag_solve == IloTrue) {
                            IloNumArray dual_equality_sub(env_sub);
                            //double optimal_value = cplex_sub.getObjValue(); // get the optimal value
                            cplex_sub.getDuals(dual_equality_sub,constraintsEquality_sub);
                            std::vector<double> minus_dual_sub;
                            for (int dual_idx = 0; dual_idx < model_parameters.D.getRowLength(); ++dual_idx) {
                                minus_dual_sub.push_back((-1.0) * dual_equality_sub[dual_idx]);
                            }
                            // compute NSQG
                            std::vector<double> subgradient = model_parameters.C.fast_rightMultiply(minus_dual_sub);
                            // stochastic C
                            if (flag_Ce == true || flag_Ci == true) {
                                // equality
                                for (int idx_Ce = 0; idx_Ce < RHS_dataset[idx].Ce.size(); ++idx_Ce) {
                                    subgradient[RHSmap.Ce_map[idx_Ce].second] += RHS_dataset[idx].Ce[idx_Ce] * minus_dual_sub[RHSmap.Ce_map[idx_Ce].first];
                                }
                                // inequality
                                for (int idx_Ci = 0; idx_Ci < RHS_dataset[idx].Ci.size(); ++idx_Ci) {
                                    subgradient[RHSmap.Ci_map[idx_Ci].second] += RHS_dataset[idx].Ci[idx_Ci] * minus_dual_sub[RHSmap.Ci_map[idx_Ci].first + model_parameters.num_eq];
                                }
                            }
                            for (int nsqg_idx = 0; nsqg_idx < A_colsize; ++nsqg_idx) {
                                nsqg[nsqg_idx] += subgradient[nsqg_idx] * normalized_weight;
                            }
                        }
                        else {
                            throw std::logic_error("Main solving process terminate. Subproblem is infeasible.\n");
                        }
                    } // end if (normalized_weight > SOLVER_PRECISION_UPPER)
                } // end for (int idx = 0; idx < N_nsqg; ++idx)
            } // end if (total_weight > SOLVER_PRECISION_UPPER)
            writeFile << "------------------------------------------------\n";
            writeFile << "Inner loop: Iteration  " << inner_it + 1 << "\n";
            writeFile << "Number of scenarios: " << N_nsqg << "\n";
            // write stepzie
            writeFile << "Stepsize: " << stepsize << std::endl;
            // first stage
            for (int idx = 0; idx < model_parameters.c.getNzeroLen(); ++idx) {
                nsqg[model_parameters.c.getLoc(idx)] += model_parameters.c.getVal(idx);
            }
            if (total_weight < SOLVER_PRECISION_UPPER) {
                std::cout << "Total weight is equal to 0. The estimated solution is retained.\n" << std::endl;
                for (int idx = 0; idx < A_colsize; ++idx) {
                    x_new[idx] = x_old[idx];
                }
                // write quasigradient
                writeFile << "NSQG: ";
                for (int idx = 0; idx < A_colsize - 1; ++idx) {
                    writeFile << 0 << ", ";
                }
                writeFile << 0 << std::endl;
            }
            else {
                // update solution
                for (int idx = 0; idx < A_colsize; ++idx) {
                    x_new[idx] = x_old[idx] - nsqg[idx] * stepsize;
                }
                // write quasigradient
                writeFile << "NSQG: ";
                for (int idx = 0; idx < A_colsize - 1; ++idx) {
                    writeFile << nsqg[idx] << ", ";
                }
                writeFile << nsqg[A_colsize - 1] << std::endl;
            }
            // projection
            x_new = leona_projection(x_new, model_parameters);
            // update x_leon
            for (int idx = 0; idx < A_colsize; ++idx) {
                x_leon[idx] += x_new[idx];
            }
            // update x_old
            std::cout << "x (nsqg update): (";
            writeFile << "x (nsqg update): (";
            // update x_old and go to the next iteration
            for (int x_index = 0; x_index < A_colsize; ++x_index) {
                x_old[x_index] = x_new[x_index];
                std::cout << x_new[x_index] << " ";
                writeFile << x_new[x_index] << " ";
            }
            std::cout << ")" <<std::endl;
            writeFile << ")" <<std::endl;
            std::cout << "**************************" << std::endl;
            writeFile << "**************************" << std::endl;
            // update sample size and data start index
            N_total += N_nsqg;
            data_start_idx += N_nsqg;
            N_nsqg += N_incre;
        } // end for (int inner_it = 0; inner_it < maxInner_iterates; ++inner_it)
        double mq = maxInner_iterates;
        // averaging to get x_leon
        for (int idx = 0; idx < A_colsize; ++idx) {
            x_leon[idx] = x_leon[idx] / mq;
        }
    } // end outer loop
    std::cout << "*******************************************\n";
    writeFile << "*******************************************\n";
    std::cout << "Total Number of Samples Used: " << N_total << std::endl;
    writeFile << "Total Number of Samples Used: " << N_total << std::endl;
    std::cout << "Output Solution: ";
    writeFile << "Output Solution: ";
    for (int index = 0; index < A_colsize-1; ++index) {
        std::cout << x_leon[index] << ", ";
        writeFile << x_leon[index] << ", ";
    }
    std::cout << x_leon[A_colsize - 1] << std::endl;
    writeFile << x_leon[A_colsize - 1] << std::endl;
    std::cout << "Computation Log: Finish Solving Process.\n";
    writeFile << "Computation Log: Finish Solving Process.\n";
    // write time elapsed
    double duration = (std::clock() - time_start ) / (double) CLOCKS_PER_SEC;
    writeFile << "Time elapsed(secs) : " << duration << "\n";
    writeFile << "*******************************************\n";
    
    // return results
    solverOutput res;
    res.x = x_leon;
    res.time_elapse = duration;
    res.num_it = max_it;
    res.N_total = N_total;
    return res;
}
