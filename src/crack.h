#pragma once

#include "libmesh/point.h"
#include "util.h"

namespace geom {

/*! @brief A structure to edge crack of any orientation */
struct EdgeCrack {

    /*!
        * @brief Orientation of crack
        *
        * - 1 for orientation along horizontal axis (x-axis),
        * - -1 for orientation along vertical axis (y-axis),
        * - 0 for any crack which makes nonzero angle with the horizontal axis
        */
    int d_o;
    
    /*!
        * @brief Angle the edge crack makes with the horizontal axis
        *
        * - 0 if orientation = 1,
        * - \f$ \pi/2 \f$ if orientation = -1
        * - Value for orientation = 0
        */
    double d_theta;
    
    /*! @brief Total current length of crack */
    double d_l;
    
    /*! @brief Current length of crack on top side (right side) */
    double d_lt;
    
    /*! @brief Current length of crack on bottom side (left side) */
    double d_lb;
    
    /*! @brief Velocity of top (right) crack tip */
    libMesh::Point d_vt;
    
    /*! @brief Velocity of bottom (left) crack tip */
    libMesh::Point d_vb;
    
    /*! @brief Closest id of node to top (right) crack tip */
    int d_it;
    
    /*! @brief Closest id of node to bottom (left) crack tip */
    int d_ib;
    
    /*! @brief Top (right) crack tip location */
    libMesh::Point d_pt;
    
    /*! @brief Bottom (left) crack tip location */
    libMesh::Point d_pb;
    
    /*! @brief Top (right) crack tip location (old) */
    libMesh::Point d_oldPt;
    
    /*! @brief Bottom (left) crack tip location (old) */
    libMesh::Point d_oldPb;
    
    /*! @brief Top (right) crack tip location (initial) */
    libMesh::Point d_initPt;
    
    /*! @brief Bottom (left) crack tip location (initial) */
    libMesh::Point d_initPb;
    
    /*! @brief Track top point (right) of crack */
    bool d_trackt;
    
    /*! @brief Track bottom point (left) of crack */
    bool d_trackb;
    
    /*! @brief Activation time (for example when we abruptly slit the domain
        * at specified time) */
    double d_activationTime;
    
    /*! @brief Is crack activated */
    bool d_crackAcrivated;
    
    /*!
        * @brief Constructor
        */
    EdgeCrack()
        : d_o(1), d_theta(0.), d_l(0.), d_lt(0.), d_lb(0.), d_it(-1), d_ib(-1),
            d_vt(libMesh::Point()), d_vb(libMesh::Point()),
            d_pt(libMesh::Point()), d_pb(libMesh::Point()),
            d_oldPt(libMesh::Point()), d_oldPb(libMesh::Point()),
            d_initPt(libMesh::Point()), d_initPb(libMesh::Point()),
            d_trackt(false), d_trackb(false), d_activationTime(-1.),
            d_crackAcrivated(false) {};
    
    /*!
        * @brief Checks if point lies outside the crack
        * @param p Point p to check
        * @param o Orientation of crack
        * @param pb Bottom (left) point of crack line
        * @param pt Top (right) point of crack line
        * @param theta Angle that crack line makes with horizontal axis
        * @return true if lies outside
        */
    bool ptOutside(libMesh::Point p, int o, libMesh::Point pb, libMesh::Point pt,
                    double theta = 0.0) {
    
        if (o == -1) {
    
            // straight crack along y-axis
            return util::isLess(p(1), pb(1)) or
                    util::isGreater(p(1), pt(1));
    
        } else if (o == 1) {
    
            // straight crack along x-axis ===> pb == pl, pt == pr
            return util::isLess(p(0), pb(0)) or
                    util::isGreater(p(0), pt(0));
    
        } else if (o == 0) {
    
            // straight crack at theta angle with x-axis ===> pb == pl, pt == pr
        
            // 1. pb ==> (0,0), pt ==> pt - pb, p ==> p - pb
            // 2. Apply CW rotation to new pt and new p so that crack line
            // after rotation is simply along x-axis with left point at
            // origin and right point at transformed new pt
        
            libMesh::Point pmap = util::rotateCW2D(p - pb, theta);
            libMesh::Point ptmap = util::rotateCW2D(pt - pb, theta);
        
            return util::isLess(pmap(0), 0.0) or
                    util::isGreater(pmap(0), ptmap(0));
        }
    
        return true;
    };
    
    /*!
        * @brief Checks if point lies on left(top) or right(bottom) of crack
        * @param p Point p to check
        * @param pb Bottom (left) point of crack line
        * @param pt Top (right) point of crack line
        * @return true if lies on left(top)
        */
    bool ptLeftside(libMesh::Point p, libMesh::Point pb, libMesh::Point pt) {
    
        // algorithm is same for any orientation of crack
        //            pt
        //           o
        //          /
        //  p      /
        //   o    /
        //       /
        //      /
        //     /
        //    o
        //    pb
        double a = 0.5 * ((pt(0) - pb(0)) * (p(1) - pb(1)) -
                        (p(0) - pb(0)) * (pt(1) - pb(1)));
    
        // crack is closer to left nodes
        return !util::isLess(a, 0.0);
    };
    
    /*!
        * @brief Checks if point lies on left(top) or right(bottom) of crack
        * @param p Point p to check
        * @param pb Bottom (left) point of crack line
        * @param pt Top (right) point of crack line
        * @return true if lies on right(bottom)
        */
    bool ptRightside(libMesh::Point p, libMesh::Point pb, libMesh::Point pt) {
    
        return !this->ptLeftside(p, pb, pt);
    };
    
    /*!
        * @brief Returns the string containing information about the instance of
        * the object
        *
        * @param nt Number of tabs to append before each line of string
        * @param lvl Level of information sought (higher level means more
        * information)
        * @return string String containing information about this object
        * */
    std::string printStr(int nt = 0, int lvl = 0) const {
        auto tabS = util::io::getTabS(nt);
        std::ostringstream oss;
        oss << tabS << "------- geom::EdgeCrack --------" << std::endl << std::endl;
        oss << tabS << "Orientation = " << d_o << std::endl;
        oss << tabS << "Point 1 (b/l) = " << d_pb << std::endl;
        oss << tabS << "Point 2 (t/r) = " << d_pt << std::endl;
        oss << tabS << "Length = " << d_l << std::endl;
        oss << tabS << "Angle = " << d_theta << std::endl;
        oss << tabS << std::endl;
    
        return oss.str();
    };
    
    /*!
        * @brief Prints the information about the instance of the object
        *
        * @param nt Number of tabs to append before each line of string
        * @param lvl Level of information sought (higher level means more
        * information)
        * */
    void print(int nt = 0, int lvl = 0) const { std::cout << printStr(nt, lvl); };
}; // class EdgeCrack

} // namespace geom