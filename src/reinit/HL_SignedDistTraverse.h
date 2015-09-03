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




#ifndef HL_SIGNEDDISTTRAVERSE
#define HL_SIGNEDDISTTRAVERSE

#include "ElInfo.h"
#include "FixVec.h"
#include "Traverse.h"
#include "ElementLevelSet.h"
#include "BoundaryElementDist.h"
#include "ElementUpdate.h"
#include "ElementUpdate_2d.h"
#include "ElementUpdate_3d.h"
#include "HL_SignedDist.h"
#include "VelocityExt.h"

namespace reinit
{

  using namespace AMDiS;

  class HL_SignedDistTraverse : public HL_SignedDist
  {
  public:

    HL_SignedDistTraverse(const char* name_,
                          int dim_,
                          bool doVelocityExt = false,
                          Flag velExtType_ = VEL_EXT)
      : HL_SignedDist(name_, dim_, doVelocityExt, velExtType_),
        sDOld_DOF(NULL),
        update_DOF(NULL),
        tol_reached(false),
        elVert(dim_),
        uhVal(dim_)
    {
      FUNCNAME("HL_SignedDistTraverse::HL_SignedDistTraverse");

      // ===== Read parameters from init file. =====
      Parameters::get(name + "->tolerance", tol);
      Parameters::get(name + "->maximal number of iteration steps", maxIt);
      Parameters::get(name + "->Gauss-Seidel iteration", GaussSeidelFlag);

      TEST_EXIT(tol > 0)("illegal tolerance !\n");

      // ---> for test purposes: initialization of counter variables
      calcUpdate_Cntr = 0;
      setUpdate_Cntr = 0;
      // ---> end: for test purposes
    }


    ~HL_SignedDistTraverse()
    {
      if (sDOld_DOF)
        delete sDOld_DOF;

      // ---> for test purposes: print result of update counting
      printUpdateCntr();
      // ---> end: for test purposes
    }

  protected:
    /**
     * Initializes the boundary: calculation of the distance of boundary
     * vertices to the interface.
     * Interface is given by lS_DOF and result is stored in
     * sD_DOF.
     */
    void initializeBoundary();

    /**
     * Calculates the distance function and stores result in sD_DOF.
     * Requirement: The boundary values are already set in sD_DOF.
     */
    void HL_updateIteration();

    /// Hopf-Lax element update on element elInfo.
    void HL_elementUpdate(ElInfo* elInfo);

    /// Calculate the update for the vertex xh of element elInfo.
    double calcElementUpdate(ElInfo* elInfo,
                             int nXh,
                             const DegreeOfFreedom* locInd);

    /**
     * In iteration loop: checks whether tolerance is reached.
     *    true  -  tolerance is reached
     *    false -  tolerance is not reached
     */
    bool checkTol();

    // ---> for test purposes: print result of update counting
    void printUpdateCntr()
    {
      /* cout << "\n"; */
      /* cout << "\tUpdate statistic: \n"; */
      /* cout << "\t-----------------\n"; */
      /* cout << "\tcalculated updates: " << calcUpdate_Cntr << "\n"; */
      /* cout << "\tset updates: " << setUpdate_Cntr << "\n"; */
      /* cout << "\n"; */
    }
    // ---> end: for test purposes

  protected:
    /**
     * DOF vector for the last iteration step.
     * Used during Jacobi iteration and to check whether tolerance is reached.
     */
    DOFVector<double>* sDOld_DOF;

    /**
     * Pointer to the DOF vector holding the values of the last iteration
     * step.
     *      Gauss-Seidel iteration  -  update_DOF == sD_DOF
     *      Jacobi iteration        -  update_DOF == sDOld_DOF
     */
    DOFVector<double>* update_DOF;

    /// Tolerance for Hopf-Lax update iteration loop.
    double tol;

    /// Flag showing whether tolerance tol is reached.
    bool tol_reached;

    /// Maximal number of mesh iterations for Hopf-Lax update.
    int maxIt;

    /**
     * Indicates whether Gauss-Seidel or Jacobi iteration is used.
     *     0 - Jacobi
     *   !=0 - Gauss-Seidel
     */
    int GaussSeidelFlag;

    std::vector<DegreeOfFreedom> locInd;

    FixVec<WorldVector<double> *, VERTEX> elVert;

    FixVec<double, VERTEX> uhVal;

    // ---> for test purposes: variables to count calculated and set updates
    int calcUpdate_Cntr;
    int setUpdate_Cntr;
    // ---> end: for test purposes

  };

} // end namespace reinit

using reinit::HL_SignedDistTraverse;

#endif  // HL_SIGNEDDISTTRAVERSE
