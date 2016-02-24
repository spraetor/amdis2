/** \file DualTraverse.h */

#pragma once

#include "Traverse.h"
#include "Flag.h"
#include "AMDiS_fwd.h"

namespace AMDiS
{
  /** \brief
   * Stores the four pointers to element info structures, that are required for the
   * dual mesh traverse.
   */
  struct DualElInfo
  {
    ElInfo* rowElInfo;    ///< elInfo related to testfunction
    ElInfo* colElInfo;    ///< elInfo related to trialfunction
    ElInfo* smallElInfo;  ///< the smaller element of (rowElInfo, colElInfo) with refinementPath relative to largeElInfo
    ElInfo* largeElInfo;  ///< the larger element of (rowElInfo, colElInfo)
  };

  /// Parallel traversal of two meshes.
  class DualTraverse
  {
  public:
    DualTraverse()
      : fillSubElemMat(false),
        basisFcts(NULL)
    {}

    virtual ~DualTraverse() {}

    /// Start dual traversal
    bool traverseFirst(Mesh* mesh1,
                       Mesh* mesh2,
                       int level1,
                       int level2,
                       Flag flag1,
                       Flag flag2,
                       ElInfo** elInfo1,
                       ElInfo** elInfo2,
                       ElInfo** elInfoSmall,
                       ElInfo** elInfoLarge);

    /// Alternative use for starting dual traversal.
    inline bool traverseFirst(Mesh* mesh1, Mesh* mesh2,
                              int level1, int level2,
                              Flag flag1, Flag flag2,
                              DualElInfo& dualElInfo)
    {
      return traverseFirst(mesh1, mesh2, level1, level2, flag1, flag2,
                           &(dualElInfo.rowElInfo),
                           &(dualElInfo.colElInfo),
                           &(dualElInfo.smallElInfo),
                           &(dualElInfo.largeElInfo));
    }

    /// Get next ElInfo combination
    bool traverseNext(ElInfo** elInfoNext1,
                      ElInfo** elInfoNext2,
                      ElInfo** elInfoSmall,
                      ElInfo** elInfoLarge);

    /// Alternative use for getting the next elements in the dual traversal.
    inline bool traverseNext(DualElInfo& dualElInfo)
    {
      return traverseNext(&(dualElInfo.rowElInfo),
                          &(dualElInfo.colElInfo),
                          &(dualElInfo.smallElInfo),
                          &(dualElInfo.largeElInfo));
    }

    bool check(ElInfo** elInfo1,
               ElInfo** elInfo2,
               ElInfo** elInfoSmall,
               ElInfo** elInfoLarge)
    {
      prepareNextStep(elInfo1, elInfo2, elInfoSmall, elInfoLarge);
      return true;
    }

    virtual bool skipEl1(ElInfo* /*elInfo*/)
    {
      return false;
    }

    virtual bool skipEl2(ElInfo* /*elInfo*/)
    {
      return false;
    }

    inline void setFillSubElemMat(bool b, const BasisFunction* fcts)
    {
      fillSubElemMat = b;
      basisFcts = fcts;
    }

    /** \brief
     * Checks if the small element has an edge/face which is part of a given edge/face
     * of the large element. If this is the case, it returns the local number of the
     * small edge/face, and -1 otherwise.
     *
     * \param[in]  dualElInfo    Dual element info with large and small element infos.
     * \param[in]  largeFace     A local edge/face number on the large element.
     */
    static int getFace(DualElInfo* dualElInfo, int largeFace);

  protected:
    /** \brief
     * Determines smaller and larger element, determines which element(s) has to
     * be incremented in the next step
     */
    void prepareNextStep(ElInfo** elInfo1,
                         ElInfo** elInfo2,
                         ElInfo** elInfoSmall,
                         ElInfo** elInfoLarge);

    void fillSubElInfo(ElInfo* elInfo1,
                       ElInfo* elInfo2,
                       ElInfo* elInfoSmall,
                       ElInfo* elInfoLarge);

  protected:
    /// Stack for mesh 1
    TraverseStack stack1;

    /// Stack for mesh 2
    TraverseStack stack2;

    /** \brief
     * used to determine whether all small elements belonging to the large
     * element are traversed.
     */
    double rest;

    /// true if element 1 should be incremented (set in prepareNextStep())
    bool inc1;

    /// true if element 2 should be incremented (set in prepareNextStep())
    bool inc2;

    /// for level traverse of mesh 1
    int level1_;

    /// for level traverse of mesh 2
    int level2_;

    /// for leaf element level traverse of mesh 1
    bool callLeafElLevel1_;

    /// for leaf element level traverse of mesh 2
    bool callLeafElLevel2_;

    /** \brief
     * If true, during dual mesh traverse for the smaller element the transformation
     * matrix will be computed. This matrix defines the transformation mapping for
     * points defined on the larger element to the coordinates of the smaller element.
     */
    bool fillSubElemMat;

    /** \brief
     * If \ref fillSubElemMat is set to true, the corresponding transformation
     * matrices are computed. These depend on the basis functions that are used.
     */
    const BasisFunction* basisFcts;
  };

} // end namespace AMDiS
