//
//  main.cpp
//  LEONA1.0
//
//  Created by Shuotao Diao on 10/24/22.
//

#include <iostream>

#include "LEONA_utils.hpp"

// bk19_3

void new_bk19_knn2(int caseNumber, int dimension) {
    std::string folder_path = "/Users/sonny/Documents/numericalExperiment/SDkNN2/bk19_";
    folder_path = folder_path + std::to_string(dimension) + "/experiment6/case" + std::to_string(caseNumber);
    std::string validation_path = "/Users/sonny/Documents/numericalExperiment/SDkNN2/bk19/validationSet";
    std::vector<double> observed_predictor(dimension,0.5);
    observed_predictor[0] = -1.1401;//-1.1401, 0.3406, 1.3871
    observed_predictor[1] = 0.3406;
    observed_predictor[2] = 1.3871;
    std::vector<double> x_init(4,0.0);
    int m = 1;
    double max_subgradient = 190;
    double Dx = 80;
    int N_pre = 50;
    int N_incre = 1;
    //int maxOuterLoops[] = {5,10,15};
    //int maxOuterLoops[] = {1,2,3,4,5,6,7,8,9,10,11,12,13};
    int maxOuterLoops[] = {14,15,16,17,18,19,20};
    std::vector<int> it_pointer;
    //for (int idx = 0; idx < 13; ++idx) {
    //    it_pointer.push_back(maxOuterLoops[idx]);
    //}
    for (int idx = 0; idx < 7; ++idx) {
        it_pointer.push_back(maxOuterLoops[idx]);
    }
    //interface_leona_knn2(folder_path, validation_path, maxOuterLoops[12], it_pointer, x_init, observed_predictor, N_pre, N_incre, Dx, m, max_subgradient);
    interface_leona_knn2(folder_path, validation_path, maxOuterLoops[6], it_pointer, x_init, observed_predictor, N_pre, N_incre, Dx, m, max_subgradient);
}

void new_bk19_knn(int caseNumber, int dimension) {
    std::string folder_path = "/Users/sonny/Documents/numericalExperiment/SDkNN2/bk19_";
    folder_path = folder_path + std::to_string(dimension) + "/experiment6/case" + std::to_string(caseNumber);
    std::string validation_path = "/Users/sonny/Documents/numericalExperiment/SDkNN2/bk19/validationSet";
    std::vector<double> observed_predictor(dimension,0.5);
    observed_predictor[0] = -1.1401;//-1.1401, 0.3406, 1.3871
    observed_predictor[1] = 0.3406;
    observed_predictor[2] = 1.3871;
    std::vector<double> x_init(4,0.0);
    int m = 1;
    double max_subgradient = 190;
    double Dx = 80;
    int N_pre = 50;
    int N_incre = 1;
    //int maxOuterLoops[] = {5,10,15};
    int maxOuterLoops[] = {1,2,3,4,5,6,7,8,9,10,11,12,13};
    for (int idx = 12; idx < 13; ++idx) {
        interface_leona_knn(folder_path, validation_path, maxOuterLoops[idx], x_init, observed_predictor, N_pre, N_incre, Dx, m, max_subgradient);
    }
}

void new_bk19_kernel(int caseNumber, int dimension, int flag_kernel) {
    std::string folder_path = "/Users/sonny/Documents/numericalExperiment/SDkNN2/bk19_";
    folder_path = folder_path + std::to_string(dimension) + "/experiment6/case" + std::to_string(caseNumber);
    std::string validation_path = "/Users/sonny/Documents/numericalExperiment/SDkNN2/bk19/validationSet";
    std::vector<double> observed_predictor(dimension,0.5);
    observed_predictor[0] = -1.1401;//-1.1401, 0.3406, 1.3871
    observed_predictor[1] = 0.3406;
    observed_predictor[2] = 1.3871;
    std::vector<double> x_init(4,0.0);
    int m = 1;
    double max_subgradient = 190;
    double Dx = 80;
    int N_pre = 50;
    int N_incre = 1;
    int maxOuterLoops[] = {1,2,3,4,5,6,7,8,9,10,11,12,13};
    for (int idx = 12; idx < 13; ++idx) {
        interface_leona_kernel(folder_path, validation_path, maxOuterLoops[idx], x_init, observed_predictor, N_pre, N_incre, flag_kernel, Dx, m, max_subgradient);
    }
}

void new_bk19_kernel2(int caseNumber, int dimension, int flag_kernel) {
    std::string folder_path = "/Users/sonny/Documents/numericalExperiment/SDkNN2/bk19_";
    folder_path = folder_path + std::to_string(dimension) + "/experiment6/case" + std::to_string(caseNumber);
    std::string validation_path = "/Users/sonny/Documents/numericalExperiment/SDkNN2/bk19/validationSet";
    std::vector<double> observed_predictor(dimension,0.5);
    observed_predictor[0] = -1.1401;//-1.1401, 0.3406, 1.3871
    observed_predictor[1] = 0.3406;
    observed_predictor[2] = 1.3871;
    std::vector<double> x_init(4,0.0);
    int m = 1;
    double max_subgradient = 190;
    double Dx = 80;
    int N_pre = 50;
    int N_incre = 1;
    //int maxOuterLoops[] = {1,2,3,4,5,6,7,8,9,10,11,12,13};
    int maxOuterLoops[] = {14,15,16,17,18,19,20};
    std::vector<int> it_pointer;
    //for (int idx = 0; idx < 13; ++idx) {
    //    it_pointer.push_back(maxOuterLoops[idx]);
    //}
    for (int idx = 0; idx < 7; ++idx) {
        it_pointer.push_back(maxOuterLoops[idx]);
    }
    //interface_leona_kernel2(folder_path, validation_path, maxOuterLoops[12], it_pointer, x_init, observed_predictor, N_pre, N_incre, flag_kernel, Dx, m, max_subgradient);
    interface_leona_kernel2(folder_path, validation_path, maxOuterLoops[6], it_pointer, x_init, observed_predictor, N_pre, N_incre, flag_kernel, Dx, m, max_subgradient);
}


void bk19_3_knn() {
    std::string folder_path = "/Users/sonny/Documents/numericalExperiment/SDkNN2/bk19_3/experiment6/case1";
    std::string validation_path = "/Users/sonny/Documents/numericalExperiment/SDkNN2/bk19/validationSet";
    std::vector<double> observed_predictor(3,0.0);
    observed_predictor[0] = -1.1401;//-1.1401, 0.3406, 1.3871
    observed_predictor[1] = 0.3406;
    observed_predictor[2] = 1.3871;
    std::vector<double> x_init(4,0.0);
    int m = 1;
    double max_subgradient = 190;
    double Dx = 80;
    int N_pre = 50;
    int N_incre = 1;
    //int maxOuterLoops[] = {5,10,15};
    int maxOuterLoops[] = {1,2,3,4,5,6,7,8,9,10,11,12,13};
    for (int idx = 0; idx < 13; ++idx) {
        interface_leona_knn(folder_path, validation_path, maxOuterLoops[idx], x_init, observed_predictor, N_pre, N_incre, Dx, m, max_subgradient);
    }
}


void bk19_3_knn(int caseNumber) {
    std::string folder_path = "/Users/sonny/Documents/numericalExperiment/SDkNN2/bk19_3/experiment6/case" + std::to_string(caseNumber);
    std::string validation_path = "/Users/sonny/Documents/numericalExperiment/SDkNN2/bk19/validationSet";
    std::vector<double> observed_predictor(3,0.0);
    observed_predictor[0] = -1.1401;//-1.1401, 0.3406, 1.3871
    observed_predictor[1] = 0.3406;
    observed_predictor[2] = 1.3871;
    std::vector<double> x_init(4,0.0);
    int m = 1;
    double max_subgradient = 190;
    double Dx = 80;
    int N_pre = 50;
    int N_incre = 1;
    //int maxOuterLoops[] = {5,10,15};
    int maxOuterLoops[] = {1,2,3,4,5,6,7,8,9,10,11,12,13};
    for (int idx = 0; idx < 13; ++idx) {
        interface_leona_knn(folder_path, validation_path, maxOuterLoops[idx], x_init, observed_predictor, N_pre, N_incre, Dx, m, max_subgradient);
    }
}


void bk19_3_kernel() {
    std::string folder_path = "/Users/sonny/Documents/numericalExperiment/SDkNN2/bk19_3/experiment6/case1";
    std::string validation_path = "/Users/sonny/Documents/numericalExperiment/SDkNN2/bk19/validationSet";
    std::vector<double> observed_predictor(3,0.0);
    observed_predictor[0] = -1.1401;//-1.1401, 0.3406, 1.3871
    observed_predictor[1] = 0.3406;
    observed_predictor[2] = 1.3871;
    std::vector<double> x_init(4,0.0);
    int m = 1;
    double max_subgradient = 190;
    double Dx = 80;
    int N_pre = 50;
    int N_incre = 1;
    int flag_kernel = 4; // 1: Naive, 2: Epanechnikov, 3: Quartic, 4: Gaussian
    //int maxOuterLoops[] = {1,10,15};
    int maxOuterLoops[] = {1,2,3,4,5,6,7,8,9,10,11,12,13};
    for (int idx = 0; idx < 13; ++idx) {
        interface_leona_kernel(folder_path, validation_path, maxOuterLoops[idx], x_init, observed_predictor, N_pre, N_incre, flag_kernel, Dx, m, max_subgradient);
    }
}


void bk19_3_kernel(int caseNumber, int flag_kernel) {
    std::string folder_path = "/Users/sonny/Documents/numericalExperiment/SDkNN2/bk19_3/experiment6/case" + std::to_string(caseNumber);
    std::string validation_path = "/Users/sonny/Documents/numericalExperiment/SDkNN2/bk19/validationSet";
    std::vector<double> observed_predictor(3,0.0);
    observed_predictor[0] = -1.1401;//-1.1401, 0.3406, 1.3871
    observed_predictor[1] = 0.3406;
    observed_predictor[2] = 1.3871;
    std::vector<double> x_init(4,0.0);
    int m = 1;
    double max_subgradient = 190;
    double Dx = 80;
    int N_pre = 50;
    int N_incre = 1;
    int maxOuterLoops[] = {1,2,3,4,5,6,7,8,9,10,11,12,13};
    for (int idx = 0; idx < 13; ++idx) {
        interface_leona_kernel(folder_path, validation_path, maxOuterLoops[idx], x_init, observed_predictor, N_pre, N_incre, flag_kernel, Dx, m, max_subgradient);
    }
}


// bk19_5
void bk19_5_knn(int caseNumber) {
    std::string folder_path = "/Users/sonny/Documents/numericalExperiment/SDkNN2/bk19_5/experiment6/case" + std::to_string(caseNumber);
    std::string validation_path = "/Users/sonny/Documents/numericalExperiment/SDkNN2/bk19/validationSet";
    std::vector<double> observed_predictor(5,0.0);
    observed_predictor[0] = -1.1401;//-1.1401, 0.3406, 1.3871
    observed_predictor[1] = 0.3406;
    observed_predictor[2] = 1.3871;
    observed_predictor[3] = 0.5;
    observed_predictor[4] = 0.5;
    std::vector<double> x_init(4,0.0);
    int m = 1;
    double max_subgradient = 190;
    double Dx = 80;
    int N_pre = 50;
    int N_incre = 1;
    //int maxOuterLoops[] = {5,10,15};
    int maxOuterLoops[] = {1,2,3,4,5,6,7,8,9,10,11,12,13};
    for (int idx = 0; idx < 13; ++idx) {
        interface_leona_knn(folder_path, validation_path, maxOuterLoops[idx], x_init, observed_predictor, N_pre, N_incre, Dx, m, max_subgradient);
    }
}


void bk19_5_kernel(int caseNumber, int flag_kernel) {
    std::string folder_path = "/Users/sonny/Documents/numericalExperiment/SDkNN2/bk19_5/experiment6/case" + std::to_string(caseNumber);
    std::string validation_path = "/Users/sonny/Documents/numericalExperiment/SDkNN2/bk19/validationSet";
    std::vector<double> observed_predictor(5,0.0);
    observed_predictor[0] = -1.1401;//-1.1401, 0.3406, 1.3871
    observed_predictor[1] = 0.3406;
    observed_predictor[2] = 1.3871;
    observed_predictor[3] = 0.5;
    observed_predictor[4] = 0.5;
    std::vector<double> x_init(4,0.0);
    int m = 1;
    double max_subgradient = 190;
    double Dx = 80;
    int N_pre = 50;
    int N_incre = 1;
    int maxOuterLoops[] = {1,2,3,4,5,6,7,8,9,10,11,12,13};
    for (int idx = 0; idx < 13; ++idx) {
        interface_leona_kernel(folder_path, validation_path, maxOuterLoops[idx], x_init, observed_predictor, N_pre, N_incre, flag_kernel, Dx, m, max_subgradient);
    }
}

// bk19_7
void bk19_7_knn(int caseNumber) {
    std::string folder_path = "/Users/sonny/Documents/numericalExperiment/SDkNN2/bk19_7/experiment6/case" + std::to_string(caseNumber);
    std::string validation_path = "/Users/sonny/Documents/numericalExperiment/SDkNN2/bk19/validationSet";
    std::vector<double> observed_predictor(7,0.5);
    observed_predictor[0] = -1.1401;//-1.1401, 0.3406, 1.3871
    observed_predictor[1] = 0.3406;
    observed_predictor[2] = 1.3871;
    std::vector<double> x_init(4,0.0);
    int m = 1;
    double max_subgradient = 190;
    double Dx = 80;
    int N_pre = 50;
    int N_incre = 1;
    //int maxOuterLoops[] = {5,10,15};
    int maxOuterLoops[] = {1,2,3,4,5,6,7,8,9,10,11,12,13};
    for (int idx = 0; idx < 13; ++idx) {
        interface_leona_knn(folder_path, validation_path, maxOuterLoops[idx], x_init, observed_predictor, N_pre, N_incre, Dx, m, max_subgradient);
    }
}


void bk19_7_kernel(int caseNumber, int flag_kernel) {
    std::string folder_path = "/Users/sonny/Documents/numericalExperiment/SDkNN2/bk19_7/experiment6/case" + std::to_string(caseNumber);
    std::string validation_path = "/Users/sonny/Documents/numericalExperiment/SDkNN2/bk19/validationSet";
    std::vector<double> observed_predictor(7,0.5);
    observed_predictor[0] = -1.1401;//-1.1401, 0.3406, 1.3871
    observed_predictor[1] = 0.3406;
    observed_predictor[2] = 1.3871;
    std::vector<double> x_init(4,0.0);
    int m = 1;
    double max_subgradient = 190;
    double Dx = 80;
    int N_pre = 50;
    int N_incre = 1;
    int maxOuterLoops[] = {1,2,3,4,5,6,7,8,9,10,11,12,13};
    for (int idx = 0; idx < 13; ++idx) {
        interface_leona_kernel(folder_path, validation_path, maxOuterLoops[idx], x_init, observed_predictor, N_pre, N_incre, flag_kernel, Dx, m, max_subgradient);
    }
}

// p-baa99
void p_baa99_knn(int caseNumber) {
    //std::string folder_path = "/Users/sonny/Documents/numericalExperiment/SDkNN2/slp/baa99large/experiment8/case" + std::to_string(caseNumber);
    std::string folder_path = "/Users/sonny/Documents/numericalExperiment/SDkNN2/pbaa99/experiment2/case" + std::to_string(caseNumber);
    //std::string validation_path = "/Users/sonny/Documents/numericalExperiment/SDkNN2/slp/baa99large/experiment9/kNNValidation2";
    std::string validation_path = "/Users/sonny/Documents/numericalExperiment/SDkNN2/pbaa99/kNNValidation1";
    int num_predictor = 25;
    std::vector<double> observed_predictor;
    for (int idx = 0; idx < num_predictor; ++idx) {
        std::cout << 108 - ((double) idx) * (1.0 / 3.0)  << std::endl;
        observed_predictor.push_back(108 - ((double) idx) * (1.0 / 3.0));
    }
    std::vector<double> x_init(50,0.0);
    int m = 1;
    double max_subgradient = 1;
    double Dx = 0.5; // C = 1, max_subgradient and Dx is unknown
    int N_pre = 50;
    int N_incre = 1;
    //int maxOuterLoops[] = {5,10,15};
    int maxOuterLoops[] = {1,2,3,4,5,6,7,8,9,10,11,12,13};
    //int maxOuterLoops[] = {12,20};
    for (int idx = 11; idx < 12; ++idx) {
        interface_leona_knn(folder_path, validation_path, maxOuterLoops[idx], x_init, observed_predictor, N_pre, N_incre, Dx, m, max_subgradient);
    }
}

void p_baa99_kernel(int caseNumber, int flag_kernel) {
    //std::string folder_path = "/Users/sonny/Documents/numericalExperiment/SDkNN2/slp/baa99large/experiment9/case" + std::to_string(caseNumber);
    std::string folder_path = "/Users/sonny/Documents/numericalExperiment/SDkNN2/pbaa99/experiment2/case" + std::to_string(caseNumber);
    //std::string validation_path = "/Users/sonny/Documents/numericalExperiment/SDkNN2/slp/baa99large/experiment9/kNNValidation2";
    std::string validation_path = "/Users/sonny/Documents/numericalExperiment/SDkNN2/pbaa99/kNNValidation1";
    int num_predictor = 25;
    std::vector<double> observed_predictor;
    for (int idx = 0; idx < num_predictor; ++idx) {
        std::cout << 108 - ((double) idx) * (1.0 / 3.0)  << std::endl;
        observed_predictor.push_back(108 - ((double) idx) * (1.0 / 3.0));
    }
    std::vector<double> x_init(50,0.0);
    int m = 1;
    double max_subgradient = 1;
    double Dx = 0.5; // C = 1, max_subgradient and Dx is unknown
    int N_pre = 50;
    int N_incre = 1;
    //int maxOuterLoops[] = {5,10,15};
    int maxOuterLoops[] = {1,2,3,4,5,6,7,8,9,10,11,12,13};
    //int maxOuterLoops[] = {12,20};
    for (int idx = 11; idx < 12; ++idx) {
        interface_leona_kernel(folder_path, validation_path, maxOuterLoops[idx], x_init, observed_predictor, N_pre, N_incre, flag_kernel, Dx, m, max_subgradient);
    }
}
int main(int argc, const char * argv[]) {
    //bk19_3_knn();
    /*
    for (int caseNumber = 2; caseNumber < 11; ++caseNumber) {
        //bk19_3_knn(caseNumber);
        bk19_3_kernel(caseNumber, 4);
    }
     */
    //bk19_3_kernel();
    //bk19_5_knn(1);
    /*
    int dimension = 15;
    for (int caseNumber = 1; caseNumber < 31; ++caseNumber) {
        //new_bk19_knn(caseNumber, dimension);
        //new_bk19_kernel(caseNumber, dimension, 4);
        p_baa99_kernel(caseNumber, 4);
        //p_baa99_knn(caseNumber);
        //bk19_5_knn(caseNumber);
        //bk19_5_kernel(caseNumber, 4); // 1 N 2 E 3 Q 4 G
        //bk19_7_knn(caseNumber);
        //bk19_7_kernel(caseNumber, 4);
    }
     */
    for (int idx_dim = 2; idx_dim < 3; ++idx_dim) {
        int dimension = 5 + idx_dim * 2;
        for (int caseNumber = 15; caseNumber < 31; ++caseNumber) {
            //new_bk19_knn2(caseNumber, dimension);
            //new_bk19_kernel2(caseNumber, dimension,4);
            for (int kernelFlag = 1; kernelFlag < 3; ++kernelFlag) {
                new_bk19_kernel2(caseNumber, dimension, kernelFlag);
            }
        }
    }
    
    //bk19_5_kernel(1, 4);
    return 0;
}
