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

    std::vector<std::string> xyz =
        {std::string("C  C1  1.0000 2.0000 3.0000 0.5000")
        ,std::string("H  H1  1.0000 0.0000 0.0000 0.5000")
        ,std::string("H  H1  0.0000 2.0000 3.0000 0.25000")
        ,std::string("H  H2  0.0000 1.0000 2.0000 0.25000")
        ,std::string("He He2 0.0000 1.0000 0.0000 0.2000")};


    RandomCartesian rcrds(xyz);



    //std::vector<int> seeds = {822650048,8938205,51381752,90742112};
    //RandomReal reng(seeds,0.0f,1.0f,std::string("uniform"));

    //std::ofstream out("frames.xyz");
    //out.setf( std::ios::fixed, std::ios::floatfield );

    //unsigned N(0);
    //while (N<360) {
        //t_iCoords ictest = icrds.generateRandomICoords(reng);
        //std::string zmat;
        //iCoordToZMat(ictest,zmat);
        //std::vector<glm::vec3> xyz;
        //itrnl::iCoordToXYZ(ictest,xyz);
        //std::cout << zmat << std::endl;

        //t_iCoords ictest = icrds.getiCoords();
        //ictest.dhls[0] = static_cast<float>(N) * 1.0;
        //std::vector<glm::vec3> xyz;
        //itrnl::iCoordToXYZ(ictest,xyz);

        //out << xyz.size() << "\n\n";
        //unsigned cnt(0);
        //for (auto&& i : xyz) {
        //    out << ictest.type[cnt] << " " << std::setprecision(7) << i.x << " " << i.y << " " << i.z << std::endl;
        //    ++cnt;
        //}

        //++N;
    //}

    //out.close();

    //std::vector<glm::vec3> xyz;
   // itrnl::iCoordToXYZ(icrds.getiCoords(),xyz);

    return 0;
};
