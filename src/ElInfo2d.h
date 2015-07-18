/******************************************************************************
 *
 * AMDiS - Adaptive multidimensional simulations
 *
 * Copyright (C) 2013 Dresden University of Technology. All Rights Reserved.
 * Web: https://fusionforge.zih.tu-dresden.de/projects/amdis
 *
 * Authors: 
 * Simon Vey, Thomas Witkowski, Andreas Naumann, Simon Praetorius, et al.
 *
 * This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
 * WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 *
 *
 * This file is part of AMDiS
 *
 * See also license.opensource.txt in the distribution.
 * 
 ******************************************************************************/



/** \file ElInfo2d.h */

#ifndef AMDIS_ELINFO2D_H
#define AMDIS_ELINFO2D_H

#include <boost/numeric/mtl/mtl.hpp>
#include "ElInfo.h"

namespace AMDiS {

  /** \ingroup Traverse
   * \brief
   * ElInfo class for 2-dimensional elements (\ref Triangle).
   */
  class ElInfo2d : public ElInfo
  {
  public:
    /// Constructor. Calls ElInfo's protected Constructor.
    ElInfo2d(Mesh* aMesh);

    ///
    ~ElInfo2d();

    /// 2-dimensional realisation of ElInfo's fillElInfo method.
    virtual void fillElInfo(int ichild, const ElInfo *elinfo_old) override;

    /// 2-dimensional realisation of ElInfo's fillMacroInfo method.
    virtual void fillMacroInfo(const MacroElement*) override;

    /// 2-dimensional realisation of ElInfo's worldToCoord method.
    virtual int worldToCoord(const WorldVector<double>& w, DimVec<double>* l) const override;

    /// 2-dimensional realisation of ElInfo's calcGrdLambda method.
    virtual double calcGrdLambda(DimVec<WorldVector<double> >& grd_lam) override;

    /// 2-dimensional realisation of ElInfo's getNormal method.
    virtual double getNormal(int side, WorldVector<double> &normal) const override;

    /// 2-dimensional realisation of ElInfo's getElementNormal method.
    virtual double getElementNormal(WorldVector<double> &normal) const override;

    /// implements \ref Elnfo::getSubElemCoordsMat
    virtual mtl::dense2D<double>& getSubElemCoordsMat(int degree) const override;

  protected:
    /// Temp vectors for function \ref calcGrdLambda.
    WorldVector<double> e1, e2, normal;

    static double mat_d1_left_val[3][3];
    static mtl::dense2D<double> mat_d1_left;

    static double mat_d1_right_val[3][3];
    static mtl::dense2D<double> mat_d1_right;

    static double mat_d2_left_val[6][6];
    static mtl::dense2D<double> mat_d2_left;

    static double mat_d2_right_val[6][6];
    static mtl::dense2D<double> mat_d2_right;

    static double mat_d3_left_val[10][10];
    static mtl::dense2D<double> mat_d3_left;

    static double mat_d3_right_val[10][10];
    static mtl::dense2D<double> mat_d3_right;

    static double mat_d4_left_val[15][15];
    static mtl::dense2D<double> mat_d4_left;

    static double mat_d4_right_val[15][15];
    static mtl::dense2D<double> mat_d4_right;
  };

}

#endif // AMDIS_ELINFO2D_H
