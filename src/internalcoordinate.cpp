// Standary Library
#include <stdlib.h>
#include <stdio.h>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <cmath>
#include <iomanip>

// Namespace header
#include "internalcoordinate.h"

#define Rad 180.0/M_PI; /*radians to degrees*/

/****************************************

    itrnl::t_ICRandRng functions

*****************************************/
void itrnl::t_ICRandRng::setRandomRanges(std::vector< std::string > &rngin) {
    rset=true;

    std::regex patt_rrngb("([B])\\s+([-+]?[0-9]*\\.?[0-9]+)\\s+([-+]?[0-9]*\\.?[0-9]+)");
    std::regex patt_rrnga("([A])\\s+([-+]?[0-9]*\\.?[0-9]+)\\s+([-+]?[0-9]*\\.?[0-9]+)");
    std::regex patt_rrngd("([D])\\s+([-+]?[0-9]*\\.?[0-9]+)\\s+([-+]?[0-9]*\\.?[0-9]+)");

    for (auto && i : rngin) {
        std::smatch m;

        // Bond Matching
        if        (std::regex_search(i,m,patt_rrngb)) {
            rngb.push_back(std::pair<float,float>(atof(m.str(2).c_str()),atof(m.str(3).c_str())));
        } else if (std::regex_search(i,m,patt_rrnga)) {
            rnga.push_back(std::pair<float,float>(atof(m.str(2).c_str()),atof(m.str(3).c_str())));
        } else if (std::regex_search(i,m,patt_rrngd)) {
            rngd.push_back(std::pair<float,float>(atof(m.str(2).c_str()),atof(m.str(3).c_str())));
        } else {
            itrnlThrowException("A random range line does not match the expected syntax. Check the input.");
        }
    };
};

/****************************************

    itrnl::t_ICRandRng functions

*****************************************/
void itrnl::t_ICScanRng::setScanRanges(std::vector< std::string > &rngin) {
    sset=true;

    std::regex patt_rrngb("([B])\\s+([-+]?[0-9]*\\.?[0-9]+)\\s+([-+]?[0-9]*\\.?[0-9]+)");
    std::regex patt_rrnga("([A])\\s+([-+]?[0-9]*\\.?[0-9]+)\\s+([-+]?[0-9]*\\.?[0-9]+)");
    std::regex patt_rrngd("([D])\\s+([-+]?[0-9]*\\.?[0-9]+)\\s+([-+]?[0-9]*\\.?[0-9]+)");

    for (auto && i : rngin) {
        std::smatch m;

        // Bond Matching
        if        (std::regex_search(i,m,patt_rrngb)) {
            rngb.push_back(std::pair<float,float>(atof(m.str(2).c_str()),atof(m.str(3).c_str())));
        } else if (std::regex_search(i,m,patt_rrnga)) {
            rnga.push_back(std::pair<float,float>(atof(m.str(2).c_str()),atof(m.str(3).c_str())));
        } else if (std::regex_search(i,m,patt_rrngd)) {
            rngd.push_back(std::pair<float,float>(atof(m.str(2).c_str()),atof(m.str(3).c_str())));
        } else {
            itrnlThrowException("A scan range line does not match the expected syntax. Check the input.");
        }
    };
};

/****************************************

  itrnl::InternalCoordinates functions

*****************************************/
/*---------Store the bond index------------

------------------------------------------*/
void itrnl::Internalcoordinates::m_calculateBondIndex(const std::vector< glm::ivec2 > &mbond) {

    iic.bidx.reserve(mbond.size());
    for (auto && bnd : mbond) {
        int first(bnd.x);
        int second(bnd.y);

        // Some error checking
        if (first==second)
            itrnlThrowException("Self bonding detected. Check the input file.");
        if (first<0||second<0)
            itrnlThrowException("Negative atom index detected. Check the input file.");

        // Order the atoms by atom number
        if (first < second)
            iic.bidx.push_back(itrnl::CreateBondIndex(first,second));
        else
            iic.bidx.push_back(itrnl::CreateBondIndex(second,first));

        // Check that bond is unique
        for (int i=0; i<int(iic.bidx.size())-1; ++i) {
            if (itrnl::bndCompareeq(iic.bidx.back(),iic.bidx[i]))
                itrnlThrowException("Duplicate bond detected. Check the input file.");
        }
    }

    // Sort by the first element from lowest to highest
    std::sort(iic.bidx.begin(),iic.bidx.end(),itrnl::bndComparelt);
};

/*---------Calculate the Angle index-----------

Requires the bond index bidx to be populated.
----------------------------------------------*/
void itrnl::Internalcoordinates::m_calculateAngleIndex() {
    if (iic.bidx.empty()) itrnlThrowException("Cannot calculate angles; Bonds have not been defined");

    for (uint32_t i=0; i<iic.bidx.size(); ++i) {
        for (uint32_t j=i+1; j<iic.bidx.size(); ++j) {
            if (iic.bidx[j].v1==iic.bidx[i].v1)
                iic.aidx.push_back(itrnl::CreateAngleIndex(iic.bidx[i].v1,iic.bidx[i].v2,iic.bidx[j].v2));

            if (iic.bidx[j].v1>iic.bidx[i].v1) {
                break;
            }
        }

        for (uint32_t j=i+1; j<iic.bidx.size(); ++j) {
            if (iic.bidx[j].v1==iic.bidx[i].v2)
                iic.aidx.push_back(itrnl::CreateAngleIndex(iic.bidx[i].v2,iic.bidx[i].v1,iic.bidx[j].v2));

            if (iic.bidx[j].v1>iic.bidx[i].v2) {
                break;
            }
        }
    }
};

/*----------Calculate the Dihedral index------------


Requires the angle index aidx to be populated.
----------------------------------------------------*/
void itrnl::Internalcoordinates::m_calculateDihedralIndex() {
    //if (iic.aidx.empty()) throwException("Cannot calculate dihedrals; Angles have not been defined");

    for (uint32_t i=0; i<iic.aidx.size(); ++i) {
        if(iic.aidx[i].v2 > iic.aidx[i].v1) {
            for (uint32_t j=i+1; j<iic.aidx.size(); ++j) {
                if (iic.aidx[j].v1==iic.aidx[i].v2) {
                    if (iic.aidx[j].v2==iic.aidx[i].v1) {
                        iic.didx.push_back(itrnl::CreateDihedralIndex(iic.aidx[i].v3,iic.aidx[i].v1,iic.aidx[j].v1,iic.aidx[j].v3));
                    }

                    if (iic.aidx[j].v3==iic.aidx[i].v1) {
                        iic.didx.push_back(itrnl::CreateDihedralIndex(iic.aidx[i].v3,iic.aidx[i].v1,iic.aidx[j].v1,iic.aidx[j].v2));
                    }
                }
            }
        }

        if(iic.aidx[i].v3 > iic.aidx[i].v1) {
            for (uint32_t j=i+1; j<iic.aidx.size(); ++j) {
                if (iic.aidx[j].v1==iic.aidx[i].v3) {
                    if (iic.aidx[j].v2==iic.aidx[i].v1) {
                        iic.didx.push_back(itrnl::CreateDihedralIndex(iic.aidx[i].v2,iic.aidx[i].v1,iic.aidx[j].v1,iic.aidx[j].v3));
                    }

                    if (iic.aidx[j].v3==iic.aidx[i].v1) {
                        iic.didx.push_back(itrnl::CreateDihedralIndex(iic.aidx[i].v2,iic.aidx[i].v1,iic.aidx[j].v1,iic.aidx[j].v2));
                    }
                }
            }
        }
    }
};

/*---------Get the bond index------------

------------------------------------------*/
void itrnl::Internalcoordinates::m_getAtomTypes(const std::vector< std::string > &icoords) {

    std::regex patt_atom("([A-Za-z]+)", std::regex_constants::icase);
    iic.type.reserve(icoords.size());
    for ( unsigned i=0; i<icoords.size(); ++i ) {
        std::smatch m;
        std::regex_search(icoords[i],m,patt_atom);
        iic.type.push_back(m.str(0));
        //std::cout << m.str(0) << std::endl;
    }
};

/*---------Get the bond index------------

------------------------------------------*/
void itrnl::Internalcoordinates::m_getBondIndex(const std::vector< std::string > &icoords) {
    std::regex patt_bond("[A-Za-z]+\\s+(\\d+)\\s+(\\d+\\.\\d+)", std::regex_constants::icase);
    for ( unsigned i=0; i<icoords.size(); ++i ) {
        std::smatch m;
        std::regex_search(icoords[i],m,patt_bond);
        if (!m.str(1).empty()) {
            unsigned first = atoi(m.str(1).c_str());
            unsigned second = i+1;
            iic.bidx.push_back(itrnl::CreateBondIndex(first,second));
            iic.bnds.push_back(atof(m.str(2).c_str()));
            //std::cout << "BOND: " << iic.bidx.back().v1 << "->" << iic.bidx.back().v2 << " VAL: " << iic.bnds.back() << std::endl;
        }
    }
};

/*---------Get the Angle index-----------

----------------------------------------------*/
void itrnl::Internalcoordinates::m_getAngleIndex(const std::vector< std::string > &icoords) {
    std::regex patt_angle("[A-Za-z]+\\s+(\\d+)\\s+\\d+\\.\\d+\\s+(\\d+)\\s+(\\d+\\.\\d+)", std::regex_constants::icase);
    for ( unsigned i=0; i<icoords.size(); ++i ) {
        std::smatch m;
        std::regex_search(icoords[i],m,patt_angle);
        if (!m.str(2).empty()) {
            unsigned first = atoi(m.str(1).c_str());
            unsigned second = atoi(m.str(2).c_str());
            unsigned third = i+1;
            iic.aidx.push_back(itrnl::CreateAngleIndex(first,second,third));
            iic.angs.push_back(atof(m.str(3).c_str()));
            //std::cout << "ANGLE: " << iic.aidx.back().v1 << "-" << iic.aidx.back().v2 << "-" << iic.aidx.back().v3 << " VAL: " << iic.angs.back() << std::endl;
        }
    }
};

/*----------Get the Dihedral index------------

----------------------------------------------------*/
void itrnl::Internalcoordinates::m_getDihedralIndex(const std::vector< std::string > &icoords) {
    std::regex patt_dihed("[A-Za-z]+\\s+(\\d+)\\s+\\d+\\.\\d+\\s+(\\d+)\\s+\\d+\\.\\d+\\s+(\\d+)\\s+(\\d+\\.\\d+)", std::regex_constants::icase);
    for ( unsigned i=0; i<icoords.size(); ++i ) {
        std::smatch m;
        std::regex_search(icoords[i],m,patt_dihed);
        if (!m.str(3).empty()) {
            unsigned first = atoi(m.str(1).c_str());
            unsigned second = atoi(m.str(2).c_str());
            unsigned third = atoi(m.str(3).c_str());
            unsigned fourth = i+1;
            iic.didx.push_back(itrnl::CreateDihedralIndex(first,second,third,fourth));
            iic.dhls.push_back(atof(m.str(4).c_str()));
            //std::cout << "DIHLS: " << iic.didx.back().v1 << "-" << iic.didx.back().v2 << "-" << iic.didx.back().v3 << "-" << iic.didx.back().v4 << " VAL: " << iic.dhls.back() << std::endl;
        }
    }
};

/*----------Calculate the Bonds------------

Requires the bond index, bidx, to be populated.
------------------------------------------*/
void itrnl::Internalcoordinates::m_calculateBonds(const std::vector<glm::vec3> &xyz) {
    if (iic.bidx.empty()) itrnlThrowException("Cannot calculate internals without bonds");
    if (xyz.empty()) itrnlThrowException("Cannot calculate internals without coordinates");

    //std::cout << "Bonds: " << " \n";
    for (uint32_t i=0; i<iic.bidx.size(); ++i) {
        iic.bnds.push_back(glm::length(xyz[iic.bidx[i].v1]-xyz[iic.bidx[i].v2]));
    }
};

/*----------Calculate the Angles------------

Calculates the angle.  Returns empty string
if iic.aidx is not populated.
------------------------------------------*/
void itrnl::Internalcoordinates::m_calculateAngles(const std::vector<glm::vec3> &xyz) {
    if (iic.aidx.empty()) return;

    // Calculate angles
    glm::vec3 v1;
    glm::vec3 v2;

    for (uint32_t i=0; i<iic.aidx.size(); ++i) {
        v1 = glm::normalize(xyz[iic.aidx[i].v1]-xyz[iic.aidx[i].v2]);
        v2 = glm::normalize(xyz[iic.aidx[i].v1]-xyz[iic.aidx[i].v3]);

        iic.angs.push_back(glm::dot(v1,v2));
        iic.angs.back() = std::acos(iic.angs.back());
    }
};

/*--------Calculate the Dihedrals----------

Calculates the dihedrals. Returns empty string
if iic.didx is not populated.
------------------------------------------*/
void itrnl::Internalcoordinates::m_calculateDihedrals(const std::vector<glm::vec3> &xyz) {
    if (iic.didx.empty()) return;

    // Calculate Dihedrals
    glm::vec3 b1;
    glm::vec3 b2;
    glm::vec3 b3;

    glm::vec3 n1;
    glm::vec3 n2;

    for (uint32_t i=0; i<iic.didx.size(); ++i) {
        b1 = xyz[iic.didx[i].v1]-xyz[iic.didx[i].v2];
        b2 = xyz[iic.didx[i].v2]-xyz[iic.didx[i].v3];
        b3 = xyz[iic.didx[i].v3]-xyz[iic.didx[i].v4];

        n1 = glm::normalize(glm::cross( b1,b2));
        n2 = glm::normalize(glm::cross(-b3,b2));

        iic.dhls.push_back(glm::dot(n1,n2));
        iic.dhls.back() = std::acos(iic.dhls.back());
    }
};

/*---------Calculate Numeric IC-----------

Calculates internal coordinates of an xyz input
based on stored Internal Coordinate (IC) index.

------------------------------------------*/
void itrnl::Internalcoordinates::m_setInternalCoordinatesFromXYZ(const std::vector<glm::vec3> &xyz,const std::vector<std::string> &type) {
    try {
        // These Functions fill the t_iCoords iic class member variable
        // with the bonds, angles and dihedrals from the varibles
        // predefined indices
        m_calculateBonds(xyz);
        m_calculateAngles(xyz);
        m_calculateDihedrals(xyz);

        iic.type = type;

    } catch (std::string error) itrnlErrorcatch(error);
};

/*------Generate a random IC Struct--------

Generate a random IC structure based upon
initial. Store in the working IC vectors.

------------------------------------------*/
itrnl::t_iCoords itrnl::Internalcoordinates::generateRandomICoords(RandomReal &rnGen) {
    t_iCoords oic(iic); // Copy the classes initial coords

    if (rrg.isset()) {

        for (unsigned i=0; i<iic.bnds.size(); ++i) {
            float fst = rrg.getRngBnd(i).first;
            float sec = rrg.getRngBnd(i).second;

            rnGen.setRandomRange(iic.bnds[i]-fst,iic.bnds[i]+sec);
            rnGen.getRandom(oic.bnds[i]);
        }

        for (unsigned i=0; i<iic.angs.size(); ++i) {
            float fst = rrg.getRngAng(i).first;
            float sec = rrg.getRngAng(i).second;

            rnGen.setRandomRange(iic.angs[i]-fst,iic.angs[i]+sec);
            rnGen.getRandom(oic.angs[i]);
        }

        for (unsigned i=0; i<iic.dhls.size(); ++i) {
            float fst = rrg.getRngDhl(i).first;
            float sec = rrg.getRngDhl(i).second;

            rnGen.setRandomRange(iic.dhls[i]-fst,iic.dhls[i]+sec);
            rnGen.getRandom(oic.dhls[i]);
        }
    } else {
        itrnlThrowException("Random range class is not set!");
    }

    return oic;
};

/*------Generate a scan IC Struct--------

Generate a scan IC structure based upon
initial. Store in the working IC vectors.

------------------------------------------*/
itrnl::t_iCoords itrnl::Internalcoordinates::generateScanICoords() {
    t_iCoords oic(iic); // Copy the classes initial coords

    if (srg.isset()) {
        unsigned inc = srg.getCounter();

        for (unsigned i=0; i<iic.bnds.size(); ++i) {
            float fst = srg.getRngBnd(i).first;
            float sec = srg.getRngBnd(i).second;
            oic.bnds[i] = iic.bnds[i] - sec + inc * fst;
        }
        ;

        for (unsigned i=0; i<iic.angs.size(); ++i) {
            float fst = srg.getRngAng(i).first;
            float sec = srg.getRngAng(i).second;
            oic.angs[i] = iic.angs[i] - sec + inc * fst;
        }

        for (unsigned i=0; i<iic.dhls.size(); ++i) {
            float fst = srg.getRngDhl(i).first;
            float sec = srg.getRngDhl(i).second;
            oic.dhls[i] = iic.dhls[i] - sec + inc * fst;
        }
    } else {
        itrnlThrowException("Scan range class is not set!");
    }

    return oic;
};

/****************************************

     itrnl:: namespace functions

*****************************************/


/*-----------ICoords to Z-mat-------------

Convert a t_iCoords to a string containing
the zmatrix.

------------------------------------------*/
void itrnl::iCoordToZMat(const t_iCoords &ics,std::string &zmats) {
    std::vector< std::stringstream > zmat_line(ics.type.size());
    for (auto&& z : zmat_line)
        z.setf( std::ios::fixed, std::ios::floatfield );

    for (unsigned i=0; i<ics.type.size(); ++i) {
        zmat_line[i] << ics.type[i] << " ";
    }

    for (unsigned i=0; i<ics.bidx.size(); ++i) {
        zmat_line[ics.bidx[i].v2-1] << ics.bidx[i].v1 << " " << std::setprecision(7) << ics.bnds[i] << " ";
    }

    for (unsigned i=0; i<ics.aidx.size(); ++i) {
        zmat_line[ics.aidx[i].v3-1] << ics.aidx[i].v2 << " " << std::setprecision(7) << ics.angs[i] << " ";
    }

    for (unsigned i=0; i<ics.didx.size(); ++i) {
        zmat_line[ics.didx[i].v4-1] << ics.didx[i].v3 << " " << std::setprecision(7) << ics.dhls[i] << " ";
    }

    //std::cout << "|---ZMAT TEST---|\n";
    std::stringstream zmat;
    for (auto&& z : zmat_line)
        zmat << z.str() << std::endl;

    zmats = zmat.str();
};

/*--------Get CSV String From IC-----------

Returns a string of the Internal Coordinates
(IC) with a bond,angle,dihedral count pre-
pended.

------------------------------------------*/
std::string itrnl::getCsvICoordStr(const t_iCoords &ics,std::string units) {
    std::stringstream icstr;
    icstr << ics.bnds.size() << "," << ics.angs.size() << "," << ics.dhls.size() << ",";
    icstr.setf( std::ios::scientific, std::ios::floatfield );

    if (units.compare("radians")==0) {
        for (auto&& i : ics.bnds)
            icstr << std::setprecision(7) << i << ",";

        for (auto&& i : ics.angs)
            icstr << std::setprecision(7) << glm::radians(i) << ",";

        for (auto&& i : ics.dhls)
            icstr << std::setprecision(7) << glm::radians(i) << ",";
    } else {
        for (auto&& i : ics.bnds)
            icstr << std::setprecision(7) << i << ",";

        for (auto&& i : ics.angs)
            icstr << std::setprecision(7) << i << ",";

        for (auto&& i : ics.dhls)
            icstr << std::setprecision(7) << i << ",";
    }


    //std::string rtn(icstr.str());

    return icstr.str();
};

std::string v3ToStr(glm::vec3 &v) {
    std::stringstream ss;
    ss.setf( std::ios::fixed, std::ios::floatfield );
    ss << " [" << std::setprecision(7) << v.x << "," << v.y << "," << v.z << "]";
    return ss.str();
};

void printv3ToStr(std::string comment,glm::vec3 &v) {
    std::stringstream ss;
    ss.setf( std::ios::fixed, std::ios::floatfield );
    ss << " [" << std::setprecision(7) << v.x << "," << v.y << "," << v.z << "]";
    std::cout << comment << ss.str() << std::endl;
};
/*------------IC to Cartesian-------------

Calculates cartesian coordinates based on
internal coordinates

------------------------------------------*/
void itrnl::iCoordToXYZ(const t_iCoords &ics,std::vector<glm::vec3> &xyz) {
    xyz.clear();

    xyz.push_back(glm::vec3(0.0,0.0,0.0));
    xyz.push_back(glm::vec3(ics.bnds[0],0.0,0.0));

    //std::cout << ics.aidx[0].v2-1 << " " << ics.aidx[0].v1-1 << std::endl;

    if (ics.type.size() > 2) {
        glm::vec3 xyztmp = glm::normalize(glm::vec3(xyz[ics.aidx[0].v2-1] - xyz[ics.aidx[0].v1-1]));
        //printv3ToStr("Tvec:",xyztmp);
        if (xyztmp.x > 0)
            xyz.push_back(glm::rotate(xyztmp,-glm::radians(ics.angs[0]),glm::vec3(0.0,0.0,1.0)) * ics.bnds[1] + xyz[ics.aidx[0].v1-1]);
        else
            xyz.push_back(glm::rotate(xyztmp,-glm::radians(ics.angs[0]),glm::vec3(0.0,0.0,1.0)) * ics.bnds[1] + xyz[ics.aidx[0].v1-1]);

        for (unsigned i=3; i<ics.bnds.size()+1; ++i) {
            //std::cout << ics.didx[i-3].v3-1 << std::endl;
            //std::cout << ics.didx[i-3].v1-1 << std::endl;
            if (ics.angs[i-2] > 180.0 || ics.angs[i-2] < 0.0) {
                itrnlThrowException("Angle outside of bounds!");
            }

            glm::vec3 R10 = xyz[ics.didx[i-3].v3-1] - xyz[ics.didx[i-3].v2-1]; // vec[0] - vec[1]
            //printv3ToStr("R10:",R10);
            glm::vec3 R12 = xyz[ics.didx[i-3].v1-1] - xyz[ics.didx[i-3].v2-1]; // vec[0] - vec[1]
            //printv3ToStr("R12:",R12);
            glm::vec3 N = glm::normalize(glm::cross(R10,R12));
            //printv3ToStr("CROSS:",N);
            glm::vec3 rw = glm::rotate(-glm::normalize(R12),glm::radians(ics.angs[i-2]),N);
            //printv3ToStr("RW1:",rw);
            rw = glm::rotate(rw,glm::radians(ics.dhls[i-3]),glm::normalize(R12));
            //printv3ToStr("RW2:",rw);
            rw = ics.bnds[i-1] * rw + xyz[ics.didx[i-3].v1-1];
            //printv3ToStr("RW3:",rw);
            xyz.push_back(rw);
        }
    }
    //for (auto && i : xyz)
    //    std::cout << "XYZ: " << v3ToStr(i) << std::endl;
};
