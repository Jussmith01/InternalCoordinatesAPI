// Standary Library
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <fstream>
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
         "O 1 0.9737\n",
         "O 2 1.4558 1 99.6811\n",
         "H 3 0.9737 2 99.6811 1 118.6258\n"};

    std::vector< std::string > rrng =
        {"B 0.3 1.0",
         "B 0.8 2.0",
         "B 0.3 1.0",
         "A 15.0 .0",
         "A 15.0 15.0",
         "D 120.0 60.0"};

    Internalcoordinates icrds(tcrd);

    icrds.getRandRng().setRandomRanges(rrng);

    icrds.printdata();

    std::vector<int> seeds = {822650048,8938205,51381752,90742112};
    RandomReal reng(seeds,0.0f,1.0f,std::string("uniform"));

    std::ofstream out("frames.xyz");
    out.setf( std::ios::fixed, std::ios::floatfield );

    unsigned N(0);
    while (N<360) {
        //t_iCoords ictest = icrds.generateRandomICoords(reng);
        //std::string zmat;
        //iCoordToZMat(ictest,zmat);
        //std::vector<glm::vec3> xyz;
        //itrnl::iCoordToXYZ(ictest,xyz);
        //std::cout << zmat << std::endl;

        t_iCoords ictest = icrds.getiCoords();
        ictest.dhls[0] = static_cast<float>(N) * 1.0;
        std::vector<glm::vec3> xyz;
        itrnl::iCoordToXYZ(ictest,xyz);

        out << xyz.size() << "\n\n";
        unsigned cnt(0);
        for (auto&& i : xyz) {
            out << ictest.type[cnt] << " " << std::setprecision(7) << i.x << " " << i.y << " " << i.z << std::endl;
            ++cnt;
        }

        ++N;
    }

    out.close();

    //std::vector<glm::vec3> xyz;
   // itrnl::iCoordToXYZ(icrds.getiCoords(),xyz);

    return 0;
};
