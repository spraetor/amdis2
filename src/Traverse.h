/** \file Traverse.h */

/** \defgroup Traverse Traverse module
 * @{ <img src="traverse.png"> @}
 *
 * \brief
 * Contains classes used for mesh traversal.
 */

#pragma once

#include <vector>
#include <deque>

#include "Flag.h"
#include "Global.h"
#include "AMDiS_fwd.h"

namespace AMDiS 
{
  /** \ingroup Traverse
   * \brief
   * Mesh refinement and coarsening routines are examples of functions which 
   * need a non-recursive access to the mesh elements.
   * The implementation of the non-recursive mesh traversal routines uses a 
   * stack to save the tree path from a macro element to the current element. 
   * TraverseStack holds such information. Before calling the non-recursive mesh
   * traversal routines, such a stack must be allocated. The stack is 
   * initialized by each call to \ref traverseFirst().
   */
  class TraverseStack 
  {
  public:
    /// Creates an empty TraverseStack
    TraverseStack() 
      : limitedToMacroElement(-1),
        traverse_mel(NULL),
        stack_size(0),
        stack_used(0),
        save_stack_used(0),
        myThreadId(0),
        maxThreads(1)
    {}

    /// Destructor
    ~TraverseStack();

  public:
    /// Returns the first ElInfo of a non recursive traversal which fullfills the
    /// selection criterion specified by level and fill_flag and initializes the
    /// stack. After a call of traverseFirst the next Elements are traversed via 
    /// \ref traverseNext(). 
    ElInfo* traverseFirst(Mesh* mesh, int level, Flag fill_flag);

    /// Works in the same way as \ref traverseFirst defined above, but limits the
    /// traversal to one macro mesh only.
    ElInfo* traverseFirstOneMacro(Mesh* mesh, int macroIndex, int level, 
				                          Flag fill_flag);

    /// Returns the next ElInfo in a traversal initiated by \ref traverseFirst()
    ///  If NULL is returned, the traversal is finished.
    ElInfo* traverseNext(ElInfo* elinfo_old);

    /// Returns the neighbour-th neighbour of elInfoOld
    ElInfo* traverseNeighbour(ElInfo* elInfoOld, int neighbour);

    /// Returns the neighbour-th neighbour of elInfoOld
    ElInfo* traverseNeighbour3d(ElInfo* elInfoOld, int neighbour);

    /// Returns the neighbour-th neighbour of elInfoOld
    ElInfo* traverseNeighbour2d(ElInfo* elInfoOld, int neighbour);

    /// Not yet implemented
    ElInfo* traverseMultiGridLevel();

    /// Preorder traversal of all elements
    ElInfo* traverseEveryElementPreorder();

    /// Inorder traversal of all elements
    ElInfo* traverseEveryElementInorder();

    /// Postorder traversal of all elements
    ElInfo* traverseEveryElementPostorder();

    /// Only for 3d: Calls update of all ElInfo3d onjects in \ref elinfo_stack
    void update();

    void fillRefinementPath(ElInfo& elInfo, ElInfo const& upperElInfo);

    /// Is used for parallel mesh traverse.
    void setMyThreadId(int n) 
    {
      myThreadId = n;
    }

    /// Is used for parallel mesh traverse.
    void setMaxThreads(int n) 
    {
      maxThreads = n;
    }

    int getStackData(std::vector<ElInfo*>& elInfos, std::vector<int>& infos) 
    {
      elInfos = elinfo_stack;
      infos = info_stack;

      return stack_used;
    }

    /// Returns the elInfo object on the top of the stack.
    ElInfo* getElInfo()
    {
      FUNCNAME_DBG("TraverseStack::getElInfo()");

      if (stack_used < 0)
        return NULL;

      TEST_EXIT_DBG(elinfo_stack.size() > static_cast<size_t>(stack_used))
        ("Should not happen!\n");

      return elinfo_stack[stack_used];
    }

    Flag getTraverseFlag() const
    {
      return traverse_fill_flag;
    }

  private:
    /// Enlargement of the stack
    void enlargeTraverseStack();

    /// Used by \ref traverseFirst() \ref traverseNext()
    ElInfo* traverseLeafElement();

    /// Used by \ref traverseFirst() \ref traverseNext()
    ElInfo* traverseLeafElementLevel();

    /// Used by \ref traverseFirst() \ref traverseNext()
    ElInfo* traverseElementLevel();

    /// Avoids copy of a traverse stack. If copy should be possible
    /// the operator must be implemented (deep copy not flat copy!)
    void operator=(TraverseStack const& /*rhs*/) 
    {
      FUNCNAME("TraverseStack::operator=()");
      ERROR_EXIT("not implemented");
    }

    /// Avoids copy of a traverse stack. If copy should be possible
    /// the operator must be implemented (deep copy not flat copy!)
    TraverseStack(TraverseStack const&) 
    {
      FUNCNAME("TraverseStack::TraverseStack()");
      ERROR_EXIT("not implemented");
    }

  private:
    /// Iterator to the current MacroElement
    std::deque<MacroElement*>::const_iterator currentMacro;

    /// Mesh which is currently traversed
    Mesh* traverse_mesh;

    /// Traverse level. Used if CALL_EL_LEVEL or CALL_LEAF_EL_LEVEL are set in
    /// \ref traverse_fill_flag
    int traverse_level;

    /// Flags for traversal. Specifies which elements are visited and which
    /// information are filled to the ElInfo objects.
    Flag traverse_fill_flag;

    /// If -1, the whole mesh is traversed. Otherwise, only the macro element
    /// with the given element index is traversed.
    int limitedToMacroElement;
  
    /// current macro element
    MacroElement const* traverse_mel;

    /// Current size of the stack
    int stack_size;

    /// Used size of the stack
    int stack_used;

    ///
    std::vector<ElInfo*> elinfo_stack;

    /// Stores for all level, which children of this element was already 
    /// visited. If 0, no children were visited (this status occurs only
    /// in intermediate computations). If 1, only the first was vistied.
    /// If 2, both children of the element on this level were already 
    /// visited.
    std::vector<int> info_stack;

    ///
    MacroElement const* save_traverse_mel;

    ///
    std::vector<ElInfo*> save_elinfo_stack;

    ///
    std::vector<unsigned char> save_info_stack;

    ///
    int save_stack_used;

    ///
    int id;

    /// Is used for parallel mesh traverse. The thread with the id
    /// myThreadId is only allowed to access coarse elements, which id
    /// satisfies: elId % maxThreads = myThreadId.
    int myThreadId;

    /// Is used for parallel mesh traverse. The thread with the id
    /// myThreadId is only allowed to access coarse elements, which id
    /// satisfies: elId % maxThreads = myThreadId.
    int maxThreads;
  };

} // end namespace AMDiS
