/** \file MeshStructure.h */

#pragma once

#include <vector>
#include <stdint.h>

#include "AMDiS_fwd.h"
#include "Global.h"
#include "BoundaryObject.h"

#if (DEBUG != 0)
#define WITH_ELINDEX true
#else
#define WITH_ELINDEX false
#endif

namespace AMDiS
{
  class MeshStructure
  {
  public:
    MeshStructure()
      : currentIndex(0),
        currentCode(0),
        pos(0),
        currentElement(0),
        nElements(0),
        debugMode(false)
    {}

    void clear();


    /** \brief
     * Creates a mesh structure code from a mesh object by traversing it in
     * preorder.
     *
     * \param[in]  mesh           Mesh that is used to create the structure code.
     * \param[in]  macroElIndex   If the value is set to -1, the structure code is
     *                            created for the whole mesh. Otherwise, it is done
     *                            only for the macro element with this index.
     */
    void init(Mesh* mesh, int macroElIndex = -1);

    void init(BoundaryObject& bound, Element* mesh = NULL);

    void init(std::vector<uint64_t> const& initCode, int n)
    {
      code = initCode;
      nElements = n;
      reset();
    }

    /// Sets all position counters, that are used to traverse the code, to the
    /// starting position. The code itself is not changed.
    void reset();

    /// Returns whether the code is empty or not.
    bool empty()
    {
      return (nElements == 0);
    }

    void commit()
    {
      if (pos > 0)
        code.push_back(currentCode);
      reset();
    }

    bool skipBranch(MeshStructure* insert = NULL);

    ElInfo* skipBranch(ElInfo* elInfo, TraverseStack* stack);

    bool nextElement(MeshStructure* insert = NULL);

    int lookAhead(unsigned int n = 1);

    bool isLeafElement()
    {
      return (currentCode & 1) == 0;
    }

    /// Merges a mesh structure code with its own mesh structure code. The result
    /// overwrites the own mesh structure code.
    void merge(MeshStructure* struc)
    {
      MeshStructure temp(*this);
      merge(&temp, struc, this);
    }

    /** \brief
     * Fits a given mesh to the mesh structure code.
     *
     * \param debugMode     In debugMode, the whole mesh is fitted to the mesh
     *                      structure code. Otherwise, the mesh is fitted only on
     *                      the partition of the current process.
     * \param macroElIndex  If the mesh structure code represents only one macro
     *                      element, this can be denoted here by its index. In this
     *                      case, only the corresponding macro element will be
     *                      fitted to the code. Otherwise, this variable is
     *                      negative and the whole mesh will be adapted.
     */
    void fitMeshToStructure(Mesh* mesh,
                            RefinementManager* manager,
                            bool debugMode = false,
                            int macroElIndex = -1,
                            bool ignoreFinerMesh = false);

    /// Converts the mesh structure code to a string (for debugging).
    std::string toStr(bool resetCode = true);

    /// Prints the mesh structure code.
    void print(bool resetCode = true);

    /// Returns the mesh structure code.
    std::vector<uint64_t> const& getCode() const
    {
      return code;
    }

    int getNumElements() const
    {
      return nElements;
    }

    int getCurrentElement() const
    {
      return currentElement;
    }

    void setDebugMode(bool b)
    {
      debugMode = b;
    }

    /// Returns true, if the given mesh structure code is equal to this one.
    bool compare(MeshStructure& other);

    /** \brief
     * Creates a value array of a given DOFVector. This value array corresponds
     * to the mesh structure code of the element and thus can easily be used
     * to reconstruct the values of the DOFVector on the same element (e.g., after
     * the mesh and the value array has been redistributed in parallel
     * computations).
     *
     * \param[in]  macroElIndex  Index of the macro element for which the value
     *                           structure code must be created.
     * \param[in]  vec           DOFVector to be used for creating the value code.
     * \param[out] values        Resulting value structure code.
     */
    void getMeshStructureValues(int macroElIndex,
                                DOFVector<double> const* vec,
                                std::vector<double>& values,
                                bool withElIndex = WITH_ELINDEX);


    /** \brief
     * Uses a value structure code, e.g. created by \ref getMeshStructureValues,
     * to reconstruct the data of a DOFVector on a given macro element.
     *
     * \param[in]  macroElIndex  Macro element index the code is related to.
     * \param[out] vec           DOFVector that should be reconstructed.
     * \param[in]  values        Value structure code.
     */
    void setMeshStructureValues(int macroElIndex,
                                DOFVector<double>* vec,
                                std::vector<double> const& values,
                                bool withElIndex = WITH_ELINDEX);

    /// Insert a new element to the structure code. Is used by the init function.
    void insertElement(bool isLeaf);

  protected:

    void addAlongSide(Element* el,
                      GeoIndex subObj,
                      int ithObj,
                      int elType,
                      bool reverseOrder);

    /// Merges two mesh structure codes to one structure code.
    void merge(MeshStructure* structure1,
               MeshStructure* structure2,
               MeshStructure* result);

  protected:
    /// Mesh structure code.
    std::vector<uint64_t> code;

    int currentIndex;

    uint64_t currentCode;

    int pos;

    int currentElement;

    int nElements;

    /// If true, some output is printed to screen during mesh structure
    /// code generation.
    bool debugMode;

    static constexpr int structureSize = 64;

  };

} // endnmespace AMDiS
