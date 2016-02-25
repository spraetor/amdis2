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




#ifndef HL_SIGNEDDIST
#define HL_SIGNEDDIST

#include "Global.h"
#include "AdaptInfo.h"
#include "DOFVector.h"
#include "ElementFunction.h"
#include "io/FileWriter.hpp"
#include "FixVec.h"
#include "Flag.h"
#include "Initfile.h"
#include "ElementLevelSet.h"
#include "BoundaryElementDist.h"
#include "BoundaryElementLevelSetDist.h"
#include "BoundaryElementTopDist.h"
#include "BoundaryElementEdgeDist.h"
#include "BoundaryElementNormalDist.h"
#include "ElementUpdate.h"
#include "ElementUpdate_2d.h"
#include "ElementUpdate_3d.h"
#include "VelocityExt.h"
#include "VelocityExtFromVelocityField.h"

namespace reinit
{

  //////////////////////////////////////////////////////////////////////////////
  //
  //   class HL_SignedDist:
  //   --------------------
  //   Holds functionality for
  //   - calculation of signed distance function for a level set function
  //     (calcSignedDistFct())
  //   - extension of velocity from an interface on complete mesh in direction
  //     normal to the interface, where the interface is given implicitly
  //     by a level set function
  //     (calcVelocityExt())
  //
  //////////////////////////////////////////////////////////////////////////////
  class HL_SignedDist
  {
  public:
    HL_SignedDist(const char* name_,
                  int dim_,
                  bool doVelocityExt = false,
                  AMDiS::Flag velExtType_ = VEL_EXT)
      : name(name_),
        adaptInfo(NULL),
        dim(dim_),
        lS_DOF(NULL),
        sD_DOF(NULL),
        bound_DOF(NULL),
        feSpace(NULL),
        elUpdate(NULL),
        bndElDist(NULL),
        elLS(NULL),
        phi(NULL),
        velExt(NULL),
        velExtType(velExtType_)
    {
      FUNCNAME("HL_SignedDist::HL_SignedDist");

      TEST_EXIT(dim == 2 || dim == 3)("only works for dimension 2 and 3 !\n");

      // ===== Read parameters from init file. =====
      AMDiS::Parameters::get(name + "->infinity value", inftyValue);

      TEST_EXIT(inftyValue > 1000)("illegal infinity value !\n");

      // ===== Create functionality for velocity extension. =====
      if (doVelocityExt)
      {
        if (velExtType.isSet(VEL_EXT))
          velExt = new VelocityExt(dim);
        else
          velExt = new VelocityExtFromVelocityField(dim);
      }
    }

    /**
     * Virtual destructor.
     */
    virtual ~HL_SignedDist()
    {
      if (elUpdate)
        delete elUpdate;
      if (bndElDist)
        delete bndElDist;

      if (elLS)
        delete elLS;
      if (phi)
        delete phi;

      if (bound_DOF)
        delete bound_DOF;

      delete velExt;
    };

    /**
     * Calculates the signed distance function for the interface given
     * implicitly by the zero level set of lS_DOF_. The result
     * is stored in sD_DOF_.
     *
     * Alternative to level set function as DOF vector:
     * If elFct != NULL, this ElementFunction is used as level set function.
     * In this case: only calculation of distance function (positive sign) !
     */
    void calcSignedDistFct(AMDiS::AdaptInfo& adaptInfo_,
                           const AMDiS::DOFVector<double>* lS_DOF_,
                           AMDiS::DOFVector<double>* sD_DOF_,
                           AMDiS::ElementFunction<double>* elFct = NULL);

    /**
     * Calculates the signed distance function for the interface given
     * implicitly by the zero level set of lS_DOF_. The result
     * is stored in lS_DOF_.
     */
    void calcSignedDistFct(AMDiS::AdaptInfo& adaptInfo_,
                           AMDiS::DOFVector<double>* lS_DOF_);

    /**
     * Calculates the extension of a velocity origVel_DOF_ from the interface
     * on the complete mesh in direction normal to the interface. The interface
     * is given implicitly as the zero level set of lS_DOF_. The result
     * is stored in vel_DOF_. If calcSDFct is true, the signed distance function
     * which is calculated during the extension of the velocity anyway is
     * stored in sD_DOF_.
     *
     * Alternative to level set function as DOF vector:
     * If elFct != NULL, this ElementFunction is used as level set function.
     * In this case: only calculation of distance function (positive sign) !
     */
    void calcVelocityExt(AMDiS::AdaptInfo& adaptInfo_,
                         AMDiS::DOFVector<double>* origVel_DOF_,
                         AMDiS::DOFVector<double>* vel_DOF_,
                         const AMDiS::DOFVector<double>* lS_DOF_,
                         AMDiS::DOFVector<double>* sD_DOF_,
                         bool calcSDFct,
                         AMDiS::ElementFunction<double>* elFct = NULL);

    /**
     * Calculates the extension of a velocity vel_DOF_ from the interface
     * on the complete mesh in direction normal to the interface. The interface
     * is given implicitly as the zero level set of lS_DOF_. The result
     * is stored in origVel_DOF_. If calcSDFct is true, the signed distance
     * function, which is calculated during the extension of the velocity
     * anyway, is stored in lS_DOF_.
     *
     * Alternative to level set function as DOF vector:
     * If elFct != NULL, this ElementFunction is used as level set function.
     * In this case: only calculation of distance function (positive sign) !
     */
    void calcVelocityExt(AMDiS::AdaptInfo& adaptInfo_,
                         AMDiS::DOFVector<double>* vel_DOF_,
                         AMDiS::DOFVector<double>* lS_DOF_,
                         bool calcSDFct,
                         AMDiS::ElementFunction<double>* elFct = NULL);

    /**
     * Calculates the extension of the velocity vectors in origVel_DOF_
     * from the interface on the complete mesh in direction normal to the
     * interface. The interface is given implicitly as the zero level set of
     * lS_DOF_. The result is stored in vel_DOF_.
     * If calcSDFct is true, the signed distance function,
     * which is calculated during the extension of the velocity anyway, is
     * stored in sD_DOF_.
     *
     * Alternative to level set function as DOF vector:
     * If elFct != NULL, this ElementFunction is used as level set function.
     * In this case: only calculation of distance function (positive sign) !
     */
    void calcVelocityExt(AMDiS::AdaptInfo& adaptInfo_,
                         std::vector<AMDiS::DOFVector<double> *> origVel_DOF_,
                         std::vector<AMDiS::DOFVector<double> *> vel_DOF_,
                         const AMDiS::DOFVector<double>* lS_DOF_,
                         AMDiS::DOFVector<double>* sD_DOF_,
                         bool calcSDFct,
                         AMDiS::ElementFunction<double>* elFct = NULL);

    /**
     * Calculates the extension of the velocity vectors in vel_DOF_
     * from the interface on the complete mesh in direction normal to the
     * interface. The interface is given implicitly as the zero level set
     * of lS_DOF_. The result is stored in vel_DOF_.
     * If calcSDFct is true, the signed distance function, which is calculated
     * during the extension of the velocity anyway, is stored in lS_DOF_.
     *
     * Alternative to level set function as DOF vector:
     * If elFct != NULL, this ElementFunction is used as level set function.
     * In this case: only calculation of distance function (positive sign) !
     */
    void calcVelocityExt(AMDiS::AdaptInfo& adaptInfo_,
                         std::vector<AMDiS::DOFVector<double> *> vel_DOF_,
                         AMDiS::DOFVector<double>* lS_DOF_,
                         bool calcSDFct,
                         AMDiS::ElementFunction<double>* elFct = NULL);

    /**
     * Calculates the normal velocity for the velocity field velField_ and its
     * extension from the interface on the complete mesh in direction normal
     * to the interface. The interface is given implicitly as the
     * zero level set of lS_DOF_. The result is stored in vel_DOF_.
     * If calcSDFct is true, the signed distance function
     * which is calculated during the extension of the velocity anyway is
     * stored in sD_DOF_.
     *
     * Alternative to level set function as DOF vector:
     * If elFct != NULL, this ElementFunction is used as level set function.
     * In this case: only calculation of distance function (positive sign) !
     */
    void calcVelocityExtFromVelocityField(
      AMDiS::AdaptInfo& adaptInfo_,
      std::vector<AMDiS::DOFVector<double> *>& velField_,
      AMDiS::DOFVector<double>* vel_DOF_,
      const AMDiS::DOFVector<double>* lS_DOF_,
      AMDiS::DOFVector<double>* sD_DOF_,
      bool calcSDFct,
      AMDiS::ElementFunction<double>* elFct = NULL);

    /**
     * Calculates the normal velocity for the velocity field velField_ and its
     * extension from the interface on the complete mesh in direction normal
     * to the interface. The interface is given implicitly as the
     * zero level set of lS_DOF_. The result is stored in vel_DOF_.
     * If calcSDFct is true, the signed distance function
     * which is calculated during the extension of the velocity anyway is
     * stored in lS_DOF_.
     *
     * Alternative to level set function as DOF vector:
     * If elFct != NULL, this ElementFunction is used as level set function.
     * In this case: only calculation of distance function (positive sign) !
     */
    void calcVelocityExtFromVelocityField(
      AMDiS::AdaptInfo& adaptInfo_,
      std::vector<AMDiS::DOFVector<double> *>& velField_,
      AMDiS::DOFVector<double>* vel_DOF_,
      AMDiS::DOFVector<double>* lS_DOF_,
      bool calcSDFct,
      AMDiS::ElementFunction<double>* elFct = NULL);

    /**
     * Print initial function: level set function defining the interface.
     */
    void printLevelSetFct()
    {
      AMDiS::FileWriter* fileWriter = new AMDiS::FileWriter(
        "SignedDist->level set fct output",
        feSpace->getMesh(),
        const_cast<AMDiS::DOFVector<double> *>(lS_DOF));
      fileWriter->writeFiles(*adaptInfo, false);

      delete fileWriter;
    };

    /**
     * Print signed distance function.
     */
    void printSignedDistFct()
    {
      AMDiS::FileWriter* fileWriter = new AMDiS::FileWriter(
        "SignedDist->result output",
        feSpace->getMesh(),
        sD_DOF);
      fileWriter->writeFiles(*adaptInfo, false);

      delete fileWriter;
    };

  protected:
    /**
     * Initialization.
     */
    virtual void initialize(AMDiS::ElementFunction<double>* elFct = NULL);

    /**
     * Initializes the boundary: calculation of the distance of boundary
     * vertices to the interface.
     * Interface is given by lS_DOF and result is stored in
     * sD_DOF.
     */
    virtual void initializeBoundary() = 0;

    /**
     * Calculates the distance function and stores result in sD_DOF.
     * Requirement: The boundary values are already set in sD_DOF.
     */
    virtual void HL_updateIteration() = 0;

    /**
     * Transforms the distance function into a signed distance function.
     * The sign is given by the level set function lS_DOF. The
     * signed distance function is stored in sD_DOF.
     */
    void setSign();

    /**
     * Print boundary initialization (initial function for Hopf-Lax iteration).
     */
    void printBoundInitFct()
    {
      AMDiS::FileWriter* fileWriter = new AMDiS::FileWriter("SignedDist->boundary initialization output",
          feSpace->getMesh(),
          sD_DOF);
      fileWriter->writeFiles(*adaptInfo, false);

      delete fileWriter;
    };

  public:
    /**
     * Flags to distinguish velocity extension types.
     */
    static const Flag VEL_EXT;
    static const Flag VEL_EXT_FROM_VEL_FIELD;

  protected:
    /**
     * Name of this instantiation of HL_SignedDist.
     */
    std::string name;

    /**
     * AdaptInfo.
     */
    AMDiS::AdaptInfo* adaptInfo;

    /**
     * Dimension.
     */
    int dim;

    /**
     * Level set function giving implicitely (as zero level set) the
     * interface for which the signed distance function is calculated.
     */
    const AMDiS::DOFVector<double>* lS_DOF;

    /**
     * DOF vector for the calculated signed distance function.
     * Also used during calculation.
     */
    AMDiS::DOFVector<double>* sD_DOF;

    /**
     * Marker for boundary vertices:
     *   0 - vertex is no boundary vertex
     *   1 - vertex is boundary vertex
     */
    AMDiS::DOFVector<double>* bound_DOF;

    /**
     * Finite element space.
     */
    const AMDiS::FiniteElemSpace* feSpace;

    /**
     * Initialization value "inifinity" for non-boundary vertices.
     */
    double inftyValue;

    /**
     * Pointer to ElementUpdate. Used for Hopf-Lax element update.
     */
    ElementUpdate* elUpdate;

    /**
     * Used for boundary vertex initialization: calculation of the distance
     * to the interface for all vertices of a boundary element.
     */
    BoundaryElementDist* bndElDist;

    /**
     * Holds level set function and functionalities for intersection point
     * calculation.
     */
    ElementLevelSet* elLS;

    /*
     * Level set function which implicitely gives the interface as its
     * zero level set.
     * This representation is needed for the use of class ElementLevelSet.
     */
    AMDiS::ElementFunction<double>* phi;

    /**
     * Object needed to extrapolate velocity from the interface.
     */
    VelocityExt* velExt;

    /**
     * Type of velocity extension method. Possible types:
     *     VEL_EXT                - object of class VelocityExt
     *     VEL_EXT_FROM_VEL_FIELD - object of class VelocityExtFromVelocityField
     */
    AMDiS::Flag velExtType;
  };

} // end namespace reinit

using reinit::HL_SignedDist;

#endif  // HL_SIGNEDDIST
