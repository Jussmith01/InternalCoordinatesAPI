#ifndef INTERNALS_H
#define INTERNALS_H
/*----------------------------------------------
        Written by Justin Smith ~August 2015
        E-Mail Jussmith48@gmail.com
        Copyright the Roitberg research group
        Chemistry Department
        University of Florida
        Gainesville FL.
------------------------------------------------*/
//________________________________________________________________________//
//      *************************************************************     //
//                      Compiler Error Handling
//      *************************************************************     //
/* Check for sufficient compiler version */
#if defined(__GNUC__) || defined(__GNUG__)
#if (__GNUC__ < 4)
#error "Insufficient GNU compiler version to install this library-- 4.9 or greater required"
#endif
#if (__GNUC__ == 4 && __GNUC_MINOR__ < 9)
#error "Insufficient GNU compiler version to install this library-- 4.9 or greater required"
#endif
#else
#warning "Currently only GNU compilers are supported and tested, but go ahead if you know what you're doing."
#endif

/*----------------------------------------------
                 Fatal Error
    Aborts runtime with location and _error
    text if called.
-----------------------------------------------*/
#define itrnlFatalError(_error)                            \
{                                                     \
    std::stringstream _location,_message;             \
    _location << __FILE__ << ":" << __LINE__;         \
    _message << "Error -- "+_error.str();             \
    std::cerr << "File "                              \
              << _location.str() << "\n"              \
              << _message.str() << "\n"               \
              << "Aborting!" << "\n";                 \
    exit(EXIT_FAILURE);                               \
};

/*----------------------------------------------
             Standard Error
    Checks for if an input string is empty,
    if it is pass the string to FatalError.
-----------------------------------------------*/
#define itrnlErrorcatch(_errchk)                             \
{                                                            \
    if(!_errchk.empty())                                     \
    {                                                        \
        std::stringstream _error;                            \
        _error << std::string(_errchk);                      \
        itrnlFatalError(_error);                                  \
    }                                                        \
};

/*----------------------------------------------
               Throw Exception
    Pass the macro a char string with an error
    and throw a std::string exception.
-----------------------------------------------*/
#define itrnlThrowException(_errstr)                              \
{                                                            \
    if (!std::string(_errstr).empty())                       \
    {                                                        \
        throw std::string("In Function: ")                   \
            + std::string(__FUNCTION__)                      \
            + std::string("() -- ")                          \
            + std::string(_errstr);                          \
    }                                                        \
};

// Random
#include "randnormflt.h"
#include <fstream>
#include <regex>
#include <random>
#include <utility>

// GLM Vector Math
#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtx/rotate_vector.hpp>
#include <glm/gtx/vector_angle.hpp>

namespace itrnl {

/*----------Bond Index-------------

Consists of:
bndindex - 8 byte aligned
           container for two
           integer indicies.

CreateBondIndex() - A function
           that produces a
           bndindex type from
           two values.

bndCompare() - A function that
            can be used with a
            std::sort for
            sorting the bonds.
-----------------------------*/
// 8 byte alignment forced to prevent cache mis-alignment
struct bndindex {
    int v1,v2;
    float bs,bf,bi;

    bndindex () {};

    bndindex (int v1,int v2,float bs,float bf)
        : v1(v1),v2(v2),bs(bs),bf(bf) {};

    bndindex (int v1,int v2,float bs,float bf,float bi)
        : v1(v1),v2(v2),bs(bs),bf(bf),bi(bi) {};

} __attribute__ ((__aligned__(16)));

// Function for creating bond indexes
inline bndindex CreateBondIndex(const int v1,const int v2) {
    bndindex tmpidx;
    tmpidx.v1=v1;
    tmpidx.v2=v2;
    return tmpidx;
};

// Function for less than comparing bndindex types -- for use by std::sort
inline bool bndComparelt (bndindex i,bndindex j) {
    return ((i.v1<=j.v1) && (i.v2<j.v2));
}

// Function for equality comparing bndindex types
inline bool bndCompareeq (bndindex i,bndindex j) {
    return ((i.v1==j.v1) && (i.v2==j.v2));
}

/*---------Angle Index------------

Consists of:
angindex - 8 byte aligned
           container for 3
           integer indicies.

CreateAngleIndex() - A function
           that produces a
           angindex type from
           3 values.
----------------------------------*/
// 8 byte alignment forced to prevent cache mis-alignment
struct angindex {
    int v1,v2,v3;
    float as,af,ai;

    angindex () {};

    angindex (int v1,int v2,int v3,float as,float af)
        : v1(v1),v2(v2),v3(v3),as(as),af(af) {};

    angindex (int v1,int v2,int v3,float as,float af,float ai)
        : v1(v1),v2(v2),v3(v3),as(as),af(af),ai(ai) {};

} __attribute__ ((__aligned__(16)));

inline angindex CreateAngleIndex(const int v1,const int v2,const int v3) {
    angindex tmpidx;
    tmpidx.v1=v1;
    tmpidx.v2=v2;
    tmpidx.v3=v3;
    return tmpidx;
};

// Function for less than comparing bndindex types -- for use by std::sort
inline bool angComparelt (angindex i,angindex j) {
    return (i.v1<j.v1);
}

/*----------Dihedral Index------------

Consists of:
dhlindex - 8 byte aligned
           container for 4
           integer indicies.

CreateDihedralIndex() - A
           function that
           produces a dhlindex
           type from 4 values.
-------------------------------------*/
// 8 byte alignment forced to prevent cache mis-alignment
struct dhlindex {
    int v1,v2,v3,v4;
    float ds,df,di;

    dhlindex () {};

    dhlindex (int v1,int v2,int v3,int v4,float ds,float df)
        : v1(v1),v2(v2),v3(v3),v4(v4),ds(ds),df(df) {};

    dhlindex (int v1,int v2,int v3,int v4,float ds,float df,float di)
        : v1(v1),v2(v2),v3(v3),v4(v4),ds(ds),df(df),di(di) {};

} __attribute__ ((__aligned__(16)));

inline dhlindex CreateDihedralIndex(const int v1,const int v2,const int v3,const int v4) {
    dhlindex tmpidx;
    tmpidx.v1=v1;
    tmpidx.v2=v2;
    tmpidx.v3=v3;
    tmpidx.v4=v4;
    return tmpidx;
};

/*--------Internal Coordinates Type----------


This type stores the indexes for the bonds,
angles and dihedrals along with the values
of the coordinates and types of atoms in
the molecule.
----------------------------------------------*/
class t_iCoords {
public:
    std::vector<bndindex> bidx; // Bonding index
    std::vector<angindex> aidx; // Angle index
    std::vector<dhlindex> didx; // Dihedral index

    std::vector<std::string> type; // Atom Types

    std::vector<float> bnds; // Storage for bonds
    std::vector<float> angs; // Storage for angles
    std::vector<float> dhls; // Storage for dihedrals

    void clear () {
        bidx.clear();
        aidx.clear();
        didx.clear();
        type.clear();
        bnds.clear();
        angs.clear();
        dhls.clear();
    };
};

/*--------Internal Coordinates RandRng----------


This type stores the Ranges of the random
purturbations of the IC.

----------------------------------------------*/
class t_ICRandRng {

    bool rset;

    std::vector< std::pair<float,float> > rngb; // Range of bonds
    std::vector< std::pair<float,float> > rnga; // Range of angles
    std::vector< std::pair<float,float> > rngd; // Range of dihedrals

public:

    t_ICRandRng () :
        rset(false) {
    };

    // Returns true if the range values are set.
    bool isset() {
        return rset;
    };

    // Const Access Functions
    const std::pair<float,float>& getRngBnd(unsigned i) {
        return rngb[i];
    };
    const std::pair<float,float>& getRngAng(unsigned i) {
        return rnga[i];
    };
    const std::pair<float,float>& getRngDhl(unsigned i) {
        return rngd[i];
    };

    // Function for defining the random ranges
    void setRandomRanges (std::vector< std::string > &rngin);

    void clear () {
        rnga.clear();
        rngb.clear();
        rngd.clear();
    };
};

/*--------Internal Coordinates RandRng----------


This type stores the Ranges of the random
purturbations of the IC.

----------------------------------------------*/
class t_ICScanRng {

    bool sset;
    unsigned scnt;

    std::vector< std::pair<float,float> > rngb; // Range of bonds
    std::vector< std::pair<float,float> > rnga; // Range of angles
    std::vector< std::pair<float,float> > rngd; // Range of dihedrals

public:

    t_ICScanRng () :
        sset(false),scnt(0) {
    };

    // Returns true if the range values are set.
    bool isset() {
        return sset;
    };

    // Returns the scan counter and increment it by one.
    unsigned getCounter() {
        unsigned tmp(scnt);
        ++scnt;
        return tmp;
    };

    // Const Access Functions
    const std::pair<float,float>& getRngBnd(unsigned i) {
        return rngb[i];
    };
    const std::pair<float,float>& getRngAng(unsigned i) {
        return rnga[i];
    };
    const std::pair<float,float>& getRngDhl(unsigned i) {
        return rngd[i];
    };

    // Function for defining the scan ranges
    void setScanRanges (std::vector< std::string > &scnin);

    void clear () {
        rnga.clear();
        rngb.clear();
        rngd.clear();
    };
};

/*--------Internal Coordinates Class----------


This class stores the indexes for the bonds,
angles and dihedrals of the molecule.
----------------------------------------------*/
class Internalcoordinates {

    t_iCoords iic; // Initial Internal Coordinates
    t_ICRandRng rrg; // Random Range Container
    t_ICScanRng srg; // Scan Range Container


    /** Member Fucntions **/
    // Calculate the bonding index
    void m_calculateBondIndex(const std::vector< glm::ivec2 > &mbond);

    // Calculate the angle index
    void m_calculateAngleIndex();

    // Calculate the dihedral index
    void m_calculateDihedralIndex();

    // Get atom types
    void m_getAtomTypes(const std::vector< std::string > &icoords);

    // Get the bonding index
    void m_getBondIndex(const std::vector< std::string > &icoords);

    // Get the angle index
    void m_getAngleIndex(const std::vector< std::string > &icoords);

    // Get the dihedral index
    void m_getDihedralIndex(const std::vector< std::string > &icoords);

    // Calculate bond lengths
    void m_calculateBonds(const std::vector<glm::vec3> &xyz);

    // Calculate angles
    void m_calculateAngles(const std::vector<glm::vec3> &xyz);

    // Calculate dihedrals
    void m_calculateDihedrals(const std::vector<glm::vec3> &xyz);

    // Create and return Comma Separated Values Internal Coordinates string
    std::string m_createCSVICstring(const std::vector<glm::vec3> &xyz);

    // Calculate the values for the internal coords based on xyz input
    void m_setInternalCoordinatesFromXYZ(const std::vector<glm::vec3> &xyz,const std::vector<std::string> &type);

public:

    // Class index constructor
    Internalcoordinates (const std::vector< glm::ivec2 > &mbond) {
        try {
            /* Determing Internal Coords */
            m_calculateBondIndex(mbond);
            m_calculateAngleIndex();
            m_calculateDihedralIndex();

        } catch (std::string error) itrnlErrorcatch(error);
    };

    // Class index and initial iternals constructor
    Internalcoordinates (const std::vector< glm::ivec2 > &mbond,const std::vector<glm::vec3> &ixyz,const std::vector<std::string> &type) {
        try {
            /* Determing (IC) Internal Coords Index */
            m_calculateBondIndex(mbond);
            m_calculateAngleIndex();
            m_calculateDihedralIndex();

            /* Calculate and store the initial IC */
            m_setInternalCoordinatesFromXYZ(ixyz,type);

        } catch (std::string error) itrnlErrorcatch(error);
    };

    // Class index and initial iternals constructor
    Internalcoordinates (const std::vector< std::string > &icoords) {
        try {
            /* Determing (IC) Internal Coords Index */
            m_getAtomTypes(icoords);
            m_getBondIndex(icoords);
            m_getAngleIndex(icoords);
            m_getDihedralIndex(icoords);

        } catch (std::string error) itrnlErrorcatch(error);
    };

    // Calculate the CSV (Comma Separated Values) string of internal coords based on xyz input
    std::string calculateCSVInternalCoordinates(const std::vector<glm::vec3> &xyz);

    // Generate a random structure
    t_iCoords generateRandomICoords(RandomReal &rnGen);

    // Generate a random structure
    t_iCoords getInitialICoords();

    // Generate a random structure
    t_iCoords generateScanICoords();

    std::vector<std::string> getAtomTypes (void) {
        return iic.type;
    }

    // Data Printer
    void printdata() {
        std::cout << "Internal Coordinates Class Setup" << std::endl;
        std::cout << "Bonds: " << iic.bidx.size() << std::endl;
        std::cout << "Angles: " << iic.aidx.size() << std::endl;
        std::cout << "Dihedrals: " << iic.didx.size() << std::endl;
        std::cout << std::endl;
    };

    t_iCoords& getiCoords() {
        return iic;
    }

    t_ICRandRng& getRandRng() {
        return rrg;
    }

    t_ICScanRng& getScanRng() {
        return srg;
    }

    // Destructor
    ~Internalcoordinates() {
        iic.clear();
    };
};

// Produce a Z-mat string from input internal coordinates
extern void iCoordToZMat(const t_iCoords &ics,std::string &zmats);

// Calculate the CSV (Comma Separated Values) string of internal coords based on xyz input
extern std::string getCsvICoordStr(const t_iCoords &ics, std::string units="degrees");

// Convert internal coordinates to XYZ
extern void iCoordToXYZ(const t_iCoords &ics,std::vector<glm::vec3> &xyz);

/*---------Random Cartesian Class------------


----------------------------------------------*/
class RandomCartesian {

    std::vector<glm::vec3>   ixyz; // Initial Cartesian Coords
    std::vector<float>       irnd; // Random ranges for each atom

    std::vector<std::string> ityp; // Input Types
    std::vector<std::string> otyp; // Output types

    std::vector<bndindex> bidx; // Bond index and random or scan bounds
    std::vector<angindex> aidx; // Angle index and random or scan bounds
    std::vector<dhlindex> didx; // Dihedral index and random or scan bounds

    std::vector<unsigned> conn1; // These vectors define the connectivity of the molecule
    std::vector<unsigned> conn2; //

    // Parse Coords Input
    void m_parsecrdsin(const std::string &crdsin);

    // Parse Connectivity Input
    void m_parseconnin(const std::string &connin);

    // Parse Random Input
    void m_parserandin(const std::string &randin);

    // Bond tranform
    void m_bndtransform(std::vector<glm::vec3> &oxyz,std::mt19937& rgenerator);

    // Angle Transform
    void m_angtransform(std::vector<glm::vec3> &oxyz,std::mt19937& rgenerator);

    // Dihedral Transform
    void m_dhltransform(std::vector<glm::vec3> &oxyz,std::mt19937& rgenerator);

    // transform molecule based on index data
    void m_tranformviaidx(std::vector<glm::vec3> &oxyz,std::mt19937& rgenerator) {
        m_bndtransform(oxyz,rgenerator);
        m_angtransform(oxyz,rgenerator);
        m_dhltransform(oxyz,rgenerator);
    };

    bool m_searchforidx(unsigned idx
                       ,std::vector<unsigned> &carr);

    void m_searchconnectivity(unsigned catom
                             ,std::vector<unsigned> &bond);

    void m_determinebondconnectivity(const bndindex &bidx
                                    ,std::vector<unsigned> &sbond
                                    ,std::vector<unsigned> &dbond);

    void m_determineangleconnectivity(const angindex &aidx
                                     ,std::vector<unsigned> &sbond
                                     ,std::vector<unsigned> &dbond);

    void m_determinedihedralconnectivity(const dhlindex &didx
                                        ,std::vector<unsigned> &sbond
                                        ,std::vector<unsigned> &dbond);

public:

    // Class index constructor
    RandomCartesian (const std::string crdsin,const std::string connin,const std::string randin);

    // Generate a set of spherical random coordinates
    void generateRandomCoordsSpherical(std::vector<glm::vec3> &oxyz,std::mt19937& rgenerator);

    // Generate a set of boxed random coordinates
    void generateRandomCoordsBox(std::vector<glm::vec3> &oxyz,RandomReal &rnGen);

    // Generate a set of random coordinates based on a random force
    void generateRandomCoordsForce(std::vector<glm::vec3> &oxyz,RandomReal &rnGen);

    // Generate a set of random coords based on dist matrix
    void generateRandomCoordsDistmat(std::vector<glm::vec3> &oxyz,RandomReal &rnGen);

    unsigned getNa () {return ixyz.size();};

    // Get the input coords
    std::vector<glm::vec3> getixyz() {
        return ixyz;
    };

    // Set the input coords
    void setixyz(std::vector<glm::vec3> & xyz) {
        ixyz = xyz;
    };

    // Get the input types
    std::vector<std::string> getitype() {
        return ityp;
    };

    // Get the output types
    std::vector<std::string> getotype() {
        return otyp;
    };
};

/*---------Random Cartesian Class------------


----------------------------------------------*/
class RandomStructureNormalMode {

    std::vector<glm::vec3> ixyz; // Initial Cartesian Coords
    std::vector<std::vector<glm::vec3>> nm; // Normal modes
    std::vector<float> fc; // Force constant

    std::vector<std::string> ityp; // Input Types
    std::vector<std::string> otyp; // Output types

    std::ofstream peout;

    // Parse Coords Input
    void m_parsecrdsin(const std::string &crdsin);

    // Parse Connectivity Input
    void m_parsenormalmodes(const std::string &normmodein);

public:

    // Class index constructor
    RandomStructureNormalMode ( const std::string crdsin,const std::string normmodein );

    ~RandomStructureNormalMode () {peout.close();};

    // Generate a set of spherical random coordinates
    void generateRandomCoords(std::vector<glm::vec3> &oxyz,float temp,std::mt19937& rgenerator);

    unsigned getNa () {return ixyz.size();};


    // Get the input coords
    std::vector<glm::vec3> getixyz() {
        return ixyz;
    };

    // Set the input coords
    void setixyz(std::vector<glm::vec3> & xyz) {
        ixyz = xyz;
    };

    // Get the input types
    std::vector<std::string> getitype() {
        return ityp;
    };

    // Get the output types
    std::vector<std::string> getotype() {
        return otyp;
    };
};

/*---------IC Distance Matrix Class------------


----------------------------------------------*/
class ScanCartesian {

    std::vector<glm::vec3> ixyz; // Initial Cartesian Coords

    unsigned scanidx; // Used to determine which scan is being set
    std::vector<unsigned> scanrcnt; // Running Counter
    std::vector<unsigned> scantcnt; // Totals Counter

    bool complete;

    std::vector<std::string> ityp; // Input Types
    std::vector<std::string> otyp; // Output types

    std::vector<bndindex> bidx; // Bond index and random or scan bounds
    std::vector<angindex> aidx; // Angle index and random or scan bounds
    std::vector<dhlindex> didx; // Dihedral index and random or scan bounds

    std::vector<unsigned> conn1; // These vectors define the connectivity of the molecule
    std::vector<unsigned> conn2; //

    unsigned totalstrct;

    // Increment Scan Counter
    void incScanCounter () {
        for (int i = scanrcnt.size() - 1; i >= 0; --i) {
            if (scanrcnt[i] == scantcnt[i]-1) {
                scanrcnt[i] = 0;
            } else {
                ++scanrcnt[i];
                break;
            }
        }

        for (unsigned i = 0;i < scanrcnt.size();++i) {
            if ( scanrcnt[i] + 1 == scantcnt[i] ) {
                complete = true;
            } else {
                complete = false;
                break;
            }
        }
    };

    // Increment Scan Index
    void incScanIndex (unsigned &idx) {
        if (idx == scanrcnt.size() - 1) {
            idx = 0;
        } else {
            ++idx;
        }
    };

    // Parse Coords Input
    void m_parsecrdsin(const std::string &crdsin);

    // Parse Connectivity Input
    void m_parseconnin(const std::string &connin);

    // Parse Random Input
    void m_parsescanin(const std::string &randin);

    // Bond tranform
    void m_bndtransform(std::vector<glm::vec3> &oxyz);

    // Angle Transform
    void m_angtransform(std::vector<glm::vec3> &oxyz);

    // Dihedral Transform
    void m_dhltransform(std::vector<glm::vec3> &oxyz);

    // transform molecule based on index data
    void m_tranformviaidx(std::vector<glm::vec3> &oxyz) {
        m_bndtransform(oxyz);
        m_angtransform(oxyz);
        m_dhltransform(oxyz);
    };

    bool m_searchforidx(unsigned idx
                       ,std::vector<unsigned> &carr);

    void m_searchconnectivity(unsigned catom
                             ,std::vector<unsigned> &bond);

    void m_determinebondconnectivity(const bndindex &bidx
                                    ,std::vector<unsigned> &sbond
                                    ,std::vector<unsigned> &dbond);

    void m_determineangleconnectivity(const angindex &aidx
                                     ,std::vector<unsigned> &sbond
                                     ,std::vector<unsigned> &dbond);

    void m_determinedihedralconnectivity(const dhlindex &didx
                                        ,std::vector<unsigned> &sbond
                                        ,std::vector<unsigned> &dbond);

public:

    // Class index constructor
    ScanCartesian (const std::string crdsin,const std::string connin,const std::string scanin);

    // Generate a set of spherical random coordinates
    bool generateNextScanStructure(std::vector<glm::vec3> &oxyz);

    /* Return number of atoms */
    unsigned getNa () {return ixyz.size();};

    /* Return number of atoms */
    unsigned getScanCount () {return totalstrct;};

    // Get the input types
    std::vector<std::string> getitype() {
        return ityp;
    };

    // Get the output types
    std::vector<std::string> getotype() {
        return otyp;
    };
};

};
#endif
