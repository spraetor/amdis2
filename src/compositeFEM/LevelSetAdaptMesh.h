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



/** \file LevelSetAdaptMesh.h */

#ifndef AMDIS_LEVELSETADAPTMESH_H
#define AMDIS_LEVELSETADAPTMESH_H

#include "AdaptInfo.h"
#include "ElInfo.h"
#include "Global.h"
#include "Mesh.h"
#include "ProblemStatBase.h"

#include "ElementLevelSet.h"

namespace compositeFEM
{

  using namespace AMDiS;

  class LevelSetAdaptMesh
  {
  public:
    /// Constructor.
    LevelSetAdaptMesh(ProblemStatBase* prob_,
                      Mesh* mesh_,
                      ElementLevelSet* elLS_,
                      double sizeInterf_,
                      double sizeNegLs_,
                      double sizePosLs_,
                      bool onlyInterfRef_ = false,
                      bool bandAtInterf_ = false,
                      double bandSize_ = -1.0)
      : prob(prob_),
        mesh(mesh_),
        elLS(elLS_),
        sizeInterf(sizeInterf_),
        sizeNegLs(sizeNegLs_),
        sizePosLs(sizePosLs_),
        onlyInterfRef(onlyInterfRef_),
        bandAtInterf(bandAtInterf_),
        bandSize(bandSize_),
        dim(mesh->getDim())
    {}

    /**
     * Calculates sizes sizeInterf, sizeNegLs, sizePosLs and bandSize
     * from number of refinements additional to initial global refinements
     * read from init file.
     */
    static void getElementSizesFromInit(const std::string probName,
                                        double& sizeInterf_,
                                        double& sizeNegLs_,
                                        double& sizePosLs_,
                                        double& bandSize_);

    /// Sets sizeInterf, sizeNegLs, sizePosLs and bandSize.
    void setElementSizes(double& sizeInterf_,
                         double& sizeNegLs_,
                         double& sizePosLs_,
                         double bandSize_ = -1.0)
    {
      FUNCNAME("LevelSetAdaptMesh::setElementSizes()");

      sizeInterf = sizeInterf_;
      sizeNegLs = sizeNegLs_;
      sizePosLs = sizePosLs_;
      bandSize = bandSize_;

      TEST_EXIT(bandSize > 0 || bandAtInterf == false)("illegal band size!\n");
    }

    /// Sets bandSize.
    void setBandSize(double& bandSize_)
    {
      FUNCNAME("LevelSetAdaptMesh::setBandSize()");

      TEST_EXIT(bandSize_ > 0)("illegal bandSize !\n");

      bandAtInterf = true;
      bandSize = bandSize_;
    }

    /// Adapts mesh to the appropriate element sizes.
    void adaptMesh(AdaptInfo* adaptInfo);

    /// Get sizeInterf.
    inline double getSizeInterf()
    {
      return sizeInterf;
    }

  protected:
    /// Element mark function in case no band around the interface is used.
    void markElement_noBand(ElInfo* elInfo);

    /// Element mark function in case a band around the interface is used.
    void markElement_Band(ElInfo* elInfo);

    /**
     * Checks whether element size elSize corresponds to the prescribed size
     * sizeWanted and sets a mark for refinement/coarsening if needed.
     */
    void markElement(ElInfo* elInfo, double elSize, double sizeWanted);

    /**
     * Determines the position of element with respect to the band around
     * the zero level set.
     */
    int calculateElementStatusBand(ElInfo* elInfo);

  protected:
    /// Problem for mesh adaption which gives refinement and coarsening manager.
    ProblemStatBase* prob;

    /// Mesh to be adapted.
    Mesh* mesh;

    /// Holds the level set function and functionalities on intersection point calculation.
    ElementLevelSet* elLS;

    /**
     * Mesh element size after mesh adaption.
     *    sizeInterf - elements cut by zero level set and - if bandAtInterf
     *                 is true - elements with at most distance bandSize to
     *                 the zero level set
     *    sizeNegLs  - elements in domain with negative level set function
     *                 values (and not within band at interface)
     *    sizePosLs  - elements in domain with positive level set function
     *                 values (and not within band at interface)
     */
    double sizeInterf;
    double sizeNegLs;
    double sizePosLs;

    /**
     * Indicates whether only interface elements or elements within a
     * band around the interface are refined (true) or refinement/coarsening
     * on the complete mesh is done (false).
     */
    bool onlyInterfRef;

    /**
     * Indicates whether band around interface/zero level set is used (true)
     * or not (false).
     * If band is used, elements with at most distance bandSize to the
     * zero level set are refined to size sizeInterf.
     */
    bool bandAtInterf;

    /// Defines bandwidth if band around interface is used.
    double bandSize;

    /**
     * Constants to describe the position of an element with respect to a band
     * around the interface.
     *    INTERIOR - in domain with negative level set function values and
     *               not in band
     *    IN_BAND  - at least on vertex of element has distance smaller than
     *               bandSize to the zero level set; the absolute value of
     *               the level set function is taken as the distance to
     *               the zero level set
     *    EXTERIOR - in domain with positive level set function values and
     *               not in band
     */
    static const int LS_ADAPT_MESH_INTERIOR = -1;
    static const int LS_ADAPT_MESH_IN_BAND = 0;
    static const int LS_ADAPT_MESH_EXTERIOR = 1;

    /// Dimension of the mesh.
    int dim;

    /// Variables used in adaptMesh(), markElement_Band() and  markElement_noBand().
    bool doRefine;
    bool doCoarsen;
    bool refinementMarkIsSet;
    bool coarseningMarkIsSet;
  };

}

using compositeFEM::LevelSetAdaptMesh;

#endif  // AMDIS_LEVELSETADAPTMESH_H
