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
    for ( unsigned i=0;i<icoords.size();++i ) {
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
    for ( unsigned i=0;i<icoords.size();++i ) {
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
    for ( unsigned i=0;i<icoords.size();++i ) {
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
    for ( unsigned i=0;i<icoords.size();++i ) {
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
void itrnl::Internalcoordinates::m_setInternalCoordinatesFromXYZ(const std::vector<glm::vec3> &xyz) {
    try {
        // These Functions fill the t_iCoords iic class member variable
        // with the bonds, angles and dihedrals from the varibles
        // predefined indices
        m_calculateBonds(xyz);
        m_calculateAngles(xyz);
        m_calculateDihedrals(xyz);

    } catch (std::string error) itrnlErrorcatch(error);
};

/*------Generate a random IC Struct--------

Generate a random IC structure based upon
initial. Store in the working IC vectors.

------------------------------------------*/
itrnl::t_iCoords itrnl::Internalcoordinates::generateRandomICoords(RandomReal &rnGen) {
    t_iCoords oic(iic);

    for (unsigned i=0; i<iic.bnds.size(); ++i) {
        rnGen.setRandomRange(iic.bnds[i]-0.1f,iic.bnds[i]+0.1f);
        rnGen.getRandom(oic.bnds[i]);
    }

    for (unsigned i=0; i<iic.angs.size(); ++i) {
        rnGen.setRandomRange(iic.angs[i]-5.0f,iic.angs[i]+5.0f);
        rnGen.getRandom(oic.angs[i]);
    }

    for (unsigned i=0; i<iic.dhls.size(); ++i) {
        rnGen.setRandomRange(iic.dhls[i]-10.0f,iic.dhls[i]+10.0f);
        rnGen.getRandom(oic.dhls[i]);
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
std::string itrnl::getCsvICoordStr(const t_iCoords &ics) {
    std::stringstream icstr;
    icstr << ics.bnds.size() << "," << ics.angs.size() << "," << ics.dhls.size() << ",";
    icstr.setf( std::ios::scientific, std::ios::floatfield );

    for (auto&& i : ics.bnds)
        icstr << std::setprecision(7) << i << ",";

    for (auto&& i : ics.angs)
        icstr << std::setprecision(7) << i << ",";

    for (auto&& i : ics.dhls)
        icstr << std::setprecision(7) << i << ",";

    //std::string rtn(icstr.str());

    return icstr.str();
};

/*------------IC to Cartesian-------------

Calculates cartesian coordinates based on
internal coordinates

------------------------------------------*/
//void itrnl::Internalcoordinates::m_calculateCartesianCoordinates(std::vector<glm::vec3> &xyz) {

//};
