// Standary Library
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <string>
#include <iomanip>
#include <vector>

// Handlers
#include "randnormflt.h"
#include "internalcoordinate.h"


int main(int argc, char *argv[]) {
    using namespace itrnl;
    std::vector< std::string > tcrd =
        {"H\n",
         "O 1 0.9\n",
         "H 2 0.9 1 105.0"};

    Internalcoordinates icrds(tcrd);

    icrds.printdata();

    std::vector<int> seeds = {822650048,8938205,51381752,90742112};
    RandomReal reng(seeds,0.0,1.0f,std::string("uniform"));

    unsigned N(0);
    while (N<10) {
        t_iCoords ictest = icrds.generateRandomICoords(reng);
        std::string zmat;
        iCoordToZMat(ictest,zmat);
        std::cout << zmat << std::endl;
        ++N;
    }
    return 0;
};
