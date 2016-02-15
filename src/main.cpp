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

    std::stringstream xyz;
    xyz << "O O 0.0000 0.0000 -0.6\n";
    xyz << "O O 0.0000 0.0000 0.60\n";
    xyz << "H H 0.0000 1.0000 -0.6\n";
    xyz << "H H 1.0000 0.0000 0.60\n";

    std::stringstream conn; // Connectivity
    conn << "3 1\n";
    conn << "2 1\n";
    conn << "2 4\n";

    std::stringstream rin;
    //rin << "B 1 2 0.5 0.5\n";
    //rin << "B 2 4 0.2 0.2\n";
    //rin << "A 2 1 3 0.0 90.0\n";
    rin << "D 3 1 2 4 90.0 90.0\n";

    RandomCartesian rcrds(xyz.str(),conn.str(),rin.str());

    std::vector<int> seeds = {822650048,8938205,51381752,90742112};
    RandomReal reng(seeds,std::string("uniform"));

    std::vector<glm::vec3> oxyz;
    rcrds.generateRandomCoordsSpherical(oxyz,reng);


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
