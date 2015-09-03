/** \file ElInfo3d.h */

#pragma once

#include <boost/numeric/mtl/mtl.hpp>

#include "ElInfo.h"

namespace AMDiS
{

  /** \ingroup Traverse
   * \brief
   * ElInfo class for 3-dimensional elements (\ref Tetrahedron).
   */
  class ElInfo3d : public ElInfo
  {
  public:
    /// Constructor. Calls ElInfo's protected Constructor.
    ElInfo3d(Mesh* aMesh)
      : ElInfo(aMesh)
    {
      tmpWorldVecs.resize(3);
      for (int i = 0; i < 3; i++)
        tmpWorldVecs[i].resize(3);
    }

    /// Assignment operator
    ElInfo3d& operator=(const ElInfo3d& rhs)
    {
      ElInfo::operator=(rhs);
      elType = rhs.elType;
      orientation = rhs.orientation;
      return *this;
    }

    /// get ElInfo's \ref orientation
    signed char getOrientation() const
    {
      return orientation;
    }

    /// set ElInfo's \ref orientation to o
    void setOrientation(signed char o)
    {
      orientation = o;
    }

    /// 3-dimensional realisation of ElInfo's fillElInfo method.
    virtual void fillElInfo(int ichild, const ElInfo* elinfo_old) override;

    /// 3-dimensional realisation of ElInfo's fillMacroInfo method.
    virtual void fillMacroInfo(const MacroElement*) override;

    /// 3-dimensional realisation of ElInfo's worldToCoord method.
    virtual int worldToCoord(const WorldVector<double>& w, DimVec<double>& l) const override;

    /// 3-dimensional realisation of ElInfo's calcGrdLambda method.
    virtual double calcGrdLambda(DimVec<WorldVector<double>>& grd_lam) override;

    /// 3-dimensional realisation of ElInfo's getNormal method.
    virtual double getNormal(int side, WorldVector<double>& normal) const override;


    virtual double getElementNormal(WorldVector<double>& normal) const override
    {
      ERROR_EXIT("ElementNormal not available in 3d!");
      return 0.0;
    }

    /// update ElInfo after refinement (of some neighbours). Only in 3d!
    void update();

    /// implements \ref Elnfo::getSubElemCoordsMat
    virtual mtl::dense2D<double>& getSubElemCoordsMat(int degree) const override;

  protected:
    /** \brief
     * +/- 1: sign of the determinant of the transformation to the reference
     * element with vertices (0,0,0), (1,1,1), (1,1,0), (1,0,0).
     */
    signed char orientation;

    /// Tmp vectors used for calculations in calcGrdLambda and getNormal().
    std::vector<std::vector<double>> tmpWorldVecs;

    static double mat_d1_left_val[4][4];
    static mtl::dense2D<double> mat_d1_left;

    static double mat_d1_l0_right_val[4][4];
    static mtl::dense2D<double> mat_d1_l0_right;

    static double mat_d1_l12_right_val[4][4];
    static mtl::dense2D<double> mat_d1_l12_right;



    static double mat_d4_left_val[35][35];
    static mtl::dense2D<double> mat_d4_left;

    static double mat_d4_l0_right_val[35][35];
    static mtl::dense2D<double> mat_d4_l0_right;

    static double mat_d4_l12_right_val[35][35];
    static mtl::dense2D<double> mat_d4_l12_right;
  };

} // end namespace AMDiS
