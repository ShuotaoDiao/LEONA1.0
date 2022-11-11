//
//  LEONA_ioDB.hpp
//  LEONA1.0
//
//  Created by Shuotao Diao on 10/24/22.
//

#ifndef LEONA_ioDB_hpp
#define LEONA_ioDB_hpp

#include <stdio.h>
#include <fstream>
#include <sstream>

#include "LEONA_dataStructure.hpp"

std::vector<std::vector<dataPoint>> readNonparametricDB(std::string readPath); // read database from a text file

void printNonparametricDB(const std::vector<std::vector<dataPoint>>& dataPointDB);
// print dataPoint
void printDataPoint(const dataPoint& dataPoint01);
void inputDBTest(); // test on input functions of nonparametric DB
#endif /* LEONA_ioDB_hpp */
