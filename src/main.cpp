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
        {glm::ivec2(0,1)
        ,glm::ivec2(0,2)
        ,glm::ivec2(0,4)
        ,glm::ivec2(2,5)
        ,glm::ivec2(3,6)
        ,glm::ivec2(3,7)
        ,glm::ivec2(3,4)
        ,glm::ivec2(4,8)
        ,glm::ivec2(4,9)};

    std::vector<glm::vec3> xyz =
        {glm::vec3(0.0000,0.7671,0.0000)
        ,glm::vec3(0.0000	-0.7671	0.0000)
        ,glm::vec3(-1.4078	1.3718	0.0000)
        ,glm::vec3()
        ,glm::vec3()
        ,glm::vec3()
        ,glm::vec3()
        ,glm::vec3()
        ,glm::vec3()
        ,glm::vec3()
        ,glm::vec3()
        ,glm::vec3()
        ,glm::vec3()
        ,glm::vec3()};

    std::vector<std::string> type({
                             "C" /*1*/
                            ,"C" /*2*/
                            ,"C" /*3*/
                            ,"C" /*4*/
                            ,"H" /*5*/
                            ,"H" /*6*/
                            ,"H" /*7*/
                            ,"H" /*8*/
                            ,"H" /*9*/
                            ,"H" /*10*/
                            ,"H" /*11*/
                            ,"H" /*12*/
                            ,"H" /*13*/
                            ,"H" /*14*/ });

    Internalcoordinates icrds(bndidx,xyz,type);

    //icrds.getRandRng().setRandomRanges(rrng);

    icrds.printdata();

    std::string zmat;
    iCoordToZMat(icrds.getiCoords(),zmat);
    std::cout << zmat << std::endl;

    std::vector<glm::vec3> xyzf;
    iCoordToXYZ(icrds.getiCoords(),xyzf);
    std::cout << "COORDINATES:" << std::endl;
    for (auto && i : xyzf) {
        std::cout << "[" << i.x << "," << i.y << "," << i.z << "]" << std::endl;
    }

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
