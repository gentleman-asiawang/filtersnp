//
// Created by 王亚洲 on 25-1-2.
//

#ifndef PARAMS_HPP
#define PARAMS_HPP

#include <string>
#include <set>

struct Params {
    std::string inputFilename;
    std::string outputFilename;
    std::string filetype;
    std::set<std::string> keep_fields = {"AC", "AF", "AN", "DP"};
};


#endif //PARAMS_HPP
