/** \file ElInfo1d.h */

#pragma once

#include <boost/numeric/mtl/mtl.hpp>

#include "ElInfo.h"

namespace AMDiS 
{

  /** \ingroup Traverse
   * \brief
   * ElInfo class for 1-dimensional elements (\ref Line).
   */
  class ElInfo1d : public ElInfo
  {
  public:
    /// Constructor. Calls ElInfo's protected Constructor.
    ElInfo1d(Mesh* aMesh) 
      : ElInfo(aMesh) 
    {}

    /// 1-dimensional realisation of ElInfo's fillElInfo method.
    virtual void fillElInfo(int ichild, const ElInfo *elinfo_old) override;

    /// 1-dimensional realisation of ElInfo's fillMacroInfo method.
    virtual void fillMacroInfo(const MacroElement*) override;

    /// 1-dimensional realisation of ElInfo's worldToCoord method.
    virtual int worldToCoord(const WorldVector<double>& w, DimVec<double>& l) const override;

    /// 1-dimensional realisation of ElInfo's calcGrdLambda method.
    virtual double calcGrdLambda(DimVec<WorldVector<double> >& grd_lam) override;

    /// 1-dimensional realisation of ElInfo's getNormal method.
    virtual double getNormal(int side, WorldVector<double> &normal) const override;

    /// 1-dimensional realisation of ElInfo's getElementNormal method.
    virtual double getElementNormal(WorldVector<double> &normal) const override;

    int getSideOfNeighbour(int i) const
    { 
      return (i + 1) % 2; 
    }

    /// implements \ref Elnfo::getSubElemCoordsMat
    virtual mtl::dense2D<double>& getSubElemCoordsMat(int degree) const override;

  protected:
    static double mat_d1_val[2][2];

    static mtl::dense2D<double> mat_d1;

    static double mat_d1_left_val[2][2];

    static mtl::dense2D<double> mat_d1_left;

    static double mat_d1_right_val[2][2];

    static mtl::dense2D<double> mat_d1_right;

    static double mat_d2_val[3][3];

    static mtl::dense2D<double> mat_d2;

    static double mat_d2_left_val[3][3];

    static mtl::dense2D<double> mat_d2_left;

    static double mat_d2_right_val[3][3];

    static mtl::dense2D<double> mat_d2_right;
  };

} // end namespace AMDiS
