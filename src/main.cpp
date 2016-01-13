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
    std::vector< glm::ivec2 > bndidx =
        {glm::ivec2(1,2)
        ,glm::ivec2(1,3)
        ,glm::ivec2(1,4)
        ,glm::ivec2(3,6)
        ,glm::ivec2(4,7)
        ,glm::ivec2(4,8)
        ,glm::ivec2(4,5)
        ,glm::ivec2(5,9)
        ,glm::ivec2(5,10)};

    std::vector<glm::vec3> xyz =
        {glm::vec3(0.0000,0.5519,0.0000)
        ,glm::vec3(1.1805,0.8224,0.0000)
        ,glm::vec3(-0.9734,1.4951,0.0000)
        ,glm::vec3(-0.5804,-0.8589,0.0000)
        ,glm::vec3(0.3994,-1.9298,0.0000)
        ,glm::vec3(	-0.5187,2.3588,0.0000)
        ,glm::vec3(	-1.2368,-0.9540,0.8740)
        ,glm::vec3(-1.2368,-0.9540,-0.8740)
        ,glm::vec3(	1.0113,-1.8200,0.8075)
        ,glm::vec3(	1.0113,-1.8200,-0.8075)};

    std::vector<std::string> type({"C","O","O","C","N","H","H","H","H","H"});

    Internalcoordinates icrds(bndidx,xyz,type);

    //icrds.getRandRng().setRandomRanges(rrng);

    icrds.printdata();

    std::string zmat;
    iCoordToZMat(icrds.getiCoords(),zmat);

    std::cout << zmat << std::endl;

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
