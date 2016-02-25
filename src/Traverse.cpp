#include "Traverse.hpp"

#include "Debug.hpp"
#include "ElInfo.hpp"
#include "Element.hpp"
#include "FixVec.hpp"
#include "Flag.hpp"
#include "Line.hpp"
#include "MacroElement.hpp"
#include "Mesh.hpp"
#include "Parametric.hpp"
#include "Tetrahedron.hpp"
#include "Triangle.hpp"

namespace AMDiS
{

  TraverseStack::~TraverseStack()
  {
    for (size_t i = 0; i < elinfo_stack.size(); i++)
      delete elinfo_stack[i];

    for (size_t i = 0; i < save_elinfo_stack.size(); i++)
      delete save_elinfo_stack[i];
  }

  ElInfo* TraverseStack::traverseFirst(Mesh* mesh, int level, Flag fill_flag)
  {
    FUNCNAME("TraverseStack::traverseFirst()");

    TEST_EXIT_DBG(fill_flag.isSet(Mesh::CALL_REVERSE_MODE) == false ||
                  fill_flag.isSet(Mesh::CALL_EVERY_EL_PREORDER))
    ("Not yet implemented!\n");

    traverse_mesh = mesh;
    traverse_level = level;
    traverse_fill_flag = fill_flag;

    TEST_EXIT_DBG(mesh)("No mesh!\n");
    TEST_EXIT(traverse_mesh->getMacroElements().size() > 0)("Mesh is empty!\n");

    if (stack_size < 1)
      enlargeTraverseStack();

    Flag FILL_ANY = Mesh::getFillAnyFlag(mesh->getDim());

    for (int i = 0; i < stack_size; i++)
      elinfo_stack[i]->setFillFlag(fill_flag & FILL_ANY);

    elinfo_stack[0]->setMesh(mesh);
    elinfo_stack[1]->setMesh(mesh);

    if (fill_flag.isSet(Mesh::CALL_LEAF_EL_LEVEL))
    {
      TEST_EXIT_DBG(level >= 0)("invalid level: %d\n", level);
    }

    traverse_mel = NULL;
    stack_used = 0;

    return traverseNext(NULL);
  }


  ElInfo* TraverseStack::traverseFirstOneMacro(Mesh* mesh, int macroIndex, int level,
      Flag fill_flag)
  {
    FUNCNAME_DBG("TraverseStack::traverseFirstOneMacro()");

    TEST_EXIT_DBG(macroIndex >= 0)("Invalid macro element index!\n");
    TEST_EXIT_DBG(traverse_fill_flag.isSet(Mesh::CALL_MG_LEVEL) == false)
    ("Multigrid level traverse not supported for one macro element only!\n");

    limitedToMacroElement = macroIndex;
    return TraverseStack::traverseFirst(mesh, level, fill_flag);
  }


  ElInfo* TraverseStack::traverseNext(ElInfo* elinfo_old)
  {
    FUNCNAME("TraverseStack::traverseNext()");

    ElInfo* elinfo = NULL;
    Parametric* parametric = traverse_mesh->getParametric();

    if (stack_used)
    {
      if (parametric)
        elinfo_old = parametric->removeParametricInfo(elinfo_old);

      TEST_EXIT_DBG(elinfo_old == elinfo_stack[stack_used])("invalid old elinfo\n");
    }
    else
    {
      TEST_EXIT_DBG(elinfo_old == NULL)("invalid old elinfo != nil\n");
    }

    if (traverse_fill_flag.isSet(Mesh::CALL_LEAF_EL))
    {
      elinfo = traverseLeafElement();
    }
    else
    {
      if (traverse_fill_flag.isSet(Mesh::CALL_LEAF_EL_LEVEL))
        elinfo = traverseLeafElementLevel();
      else if (traverse_fill_flag.isSet(Mesh::CALL_EL_LEVEL))
        elinfo = traverseElementLevel();
      else if (traverse_fill_flag.isSet(Mesh::CALL_MG_LEVEL))
        elinfo = traverseMultiGridLevel();
      else
      {
        if (traverse_fill_flag.isSet(Mesh::CALL_EVERY_EL_PREORDER))
        {
          elinfo = traverseEveryElementPreorder();
        }
        else if (traverse_fill_flag.isSet(Mesh::CALL_EVERY_EL_INORDER))
          elinfo = traverseEveryElementInorder();
        else if (traverse_fill_flag.isSet(Mesh::CALL_EVERY_EL_POSTORDER))
          elinfo = traverseEveryElementPostorder();
        else
          ERROR_EXIT("invalid traverse_flag\n");
      }
    }

    if (elinfo)
    {
      if (parametric)
        elinfo = parametric->addParametricInfo(elinfo);
      elinfo->fillDetGrdLambda();
    }

    return elinfo;
  }


  void TraverseStack::enlargeTraverseStack()
  {
    int new_stack_size = stack_size + 10;

    elinfo_stack.resize(new_stack_size, NULL);

    // create new elinfos
    for (int i = stack_size; i < new_stack_size; i++)
    {
      TEST_EXIT_DBG(elinfo_stack[i] == NULL)("???\n");
      elinfo_stack[i] = traverse_mesh->createNewElInfo();
    }

    if (stack_size > 0)
      for (int i = stack_size; i < new_stack_size; i++)
        elinfo_stack[i]->setFillFlag(elinfo_stack[0]->getFillFlag());

    info_stack.resize(new_stack_size);
    save_elinfo_stack.resize(new_stack_size, NULL);

    // create new elinfos
    for (int i = stack_size; i < new_stack_size; i++)
    {
      TEST_EXIT_DBG(save_elinfo_stack[i] == NULL)("???\n");
      save_elinfo_stack[i] = traverse_mesh->createNewElInfo();
    }
    save_info_stack.resize(new_stack_size);

    stack_size = new_stack_size;
  }


  ElInfo* TraverseStack::traverseLeafElement()
  {
    FUNCNAME_DBG("TraverseStack::traverseLeafElement()");

    Element* el = NULL;

    if (stack_used == 0)     /* first call */
    {
      currentMacro = traverse_mesh->firstMacroElement();

      if (limitedToMacroElement >= 0)
      {
        while ((*currentMacro)->getIndex() != limitedToMacroElement &&
               currentMacro != traverse_mesh->endOfMacroElements())
          currentMacro++;

        TEST_EXIT_DBG(currentMacro != traverse_mesh->endOfMacroElements())
        ("Coult not find macro element with index %d!\n", limitedToMacroElement);

      }
      else
      {
        while (((*currentMacro)->getIndex() % maxThreads != myThreadId) &&
               currentMacro != traverse_mesh->endOfMacroElements())
          currentMacro++;
      }

      if (currentMacro == traverse_mesh->endOfMacroElements())
        return NULL;

      traverse_mel = *currentMacro;
      stack_used = 1;
      elinfo_stack[stack_used]->fillMacroInfo(traverse_mel);
      info_stack[stack_used] = 0;

      el = elinfo_stack[stack_used]->getElement();
      if (el == NULL || el->getFirstChild() == NULL)
        return elinfo_stack[stack_used];
    }
    else
    {
      el = elinfo_stack[stack_used]->getElement();

      /* go up in tree until we can go down again */
      while ((stack_used > 0) &&
             ((info_stack[stack_used] >= 2) || (el->getFirstChild() == NULL)))
      {
        stack_used--;
        el = elinfo_stack[stack_used]->getElement();
      }

      /* goto next macro element */
      if (stack_used < 1)
      {
        if (limitedToMacroElement >= 0)
          return NULL;

        do
        {
          currentMacro++;
        }
        while ((currentMacro != traverse_mesh->endOfMacroElements()) &&
               ((*currentMacro)->getIndex() % maxThreads != myThreadId));

        if (currentMacro == traverse_mesh->endOfMacroElements())
          return NULL;

        traverse_mel = *currentMacro;
        stack_used = 1;
        elinfo_stack[stack_used]->fillMacroInfo(traverse_mel);
        info_stack[stack_used] = 0;
        el = elinfo_stack[stack_used]->getElement();

        if (el == NULL || el->getFirstChild() == NULL)
          return elinfo_stack[stack_used];
      }
    }

    /* go down tree until leaf */
    while (el->getFirstChild())
    {
      if (stack_used >= stack_size - 1)
        enlargeTraverseStack();

      int i = info_stack[stack_used];
      el = const_cast<Element*>(((i == 0) ? el->getFirstChild() : el->getSecondChild()));
      info_stack[stack_used]++;
      elinfo_stack[stack_used + 1]->fillElInfo(i, elinfo_stack[stack_used]);
      stack_used++;

      TEST_EXIT_DBG(stack_used < stack_size)
      ("stack_size = %d too small, level = (%d, %d)\n",
       stack_size, elinfo_stack[stack_used]->getLevel());

      info_stack[stack_used] = 0;
    }

    return elinfo_stack[stack_used];
  }


  ElInfo* TraverseStack::traverseLeafElementLevel()
  {
    FUNCNAME("TraverseStack::traverseLeafElementLevel()");

    ERROR_EXIT("not yet implemented\n");

    return NULL;
  }


  ElInfo* TraverseStack::traverseElementLevel()
  {
    ElInfo* elInfo;
    do
    {
      elInfo = traverseEveryElementPreorder();
    }
    while (elInfo != NULL && elInfo->getLevel() != traverse_level);

    return elInfo;
  }


  ElInfo* TraverseStack::traverseMultiGridLevel()
  {
    FUNCNAME_DBG("TraverseStack::traverseMultiGridLevel()");

    if (stack_used == 0)     /* first call */
    {
      currentMacro = traverse_mesh->firstMacroElement();
      traverse_mel = *currentMacro;
      if (traverse_mel == NULL)
        return NULL;

      stack_used = 1;
      elinfo_stack[stack_used]->fillMacroInfo(traverse_mel);
      info_stack[stack_used] = 0;

      if ((elinfo_stack[stack_used]->getLevel() == traverse_level) ||
          (elinfo_stack[stack_used]->getLevel() < traverse_level &&
           elinfo_stack[stack_used]->getElement()->isLeaf()))
        return elinfo_stack[stack_used];
    }

    Element* el = elinfo_stack[stack_used]->getElement();

    /* go up in tree until we can go down again */
    while ((stack_used > 0) &&
           ((info_stack[stack_used] >= 2) || (el->getFirstChild()==NULL)))
    {
      stack_used--;
      el = elinfo_stack[stack_used]->getElement();
    }


    /* goto next macro element */
    if (stack_used < 1)
    {
      currentMacro++;
      if (currentMacro == traverse_mesh->endOfMacroElements())
        return NULL;

      traverse_mel = *currentMacro;
      stack_used = 1;
      elinfo_stack[stack_used]->fillMacroInfo(traverse_mel);
      info_stack[stack_used] = 0;

      if ((elinfo_stack[stack_used]->getLevel() == traverse_level) ||
          (elinfo_stack[stack_used]->getLevel() < traverse_level &&
           elinfo_stack[stack_used]->getElement()->isLeaf()))
        return elinfo_stack[stack_used];
    }


    /* go down tree */

    if (stack_used >= stack_size - 1)
      enlargeTraverseStack();

    int i = info_stack[stack_used];
    info_stack[stack_used]++;

    elinfo_stack[stack_used + 1]->fillElInfo(i, elinfo_stack[stack_used]);
    stack_used++;

    TEST_EXIT_DBG(stack_used < stack_size)
    ("stack_size=%d too small, level=%d\n",
     stack_size, elinfo_stack[stack_used]->getLevel());

    info_stack[stack_used] = 0;

    if ((elinfo_stack[stack_used]->getLevel() == traverse_level) ||
        (elinfo_stack[stack_used]->getLevel() < traverse_level &&
         elinfo_stack[stack_used]->getElement()->isLeaf()))
      return elinfo_stack[stack_used];

    return traverseMultiGridLevel();
  }


  ElInfo* TraverseStack::traverseEveryElementPreorder()
  {
    FUNCNAME_DBG("TraverseStack::traverseEveryElementPreorder()");

    if (stack_used == 0)     /* first call */
    {
      currentMacro = traverse_mesh->firstMacroElement();

      if (limitedToMacroElement >= 0)
      {
        while ((*currentMacro)->getIndex() != limitedToMacroElement &&
               currentMacro != traverse_mesh->endOfMacroElements())
          currentMacro++;

        TEST_EXIT_DBG(currentMacro != traverse_mesh->endOfMacroElements())
        ("Coult not find macro element with index %d!\n", limitedToMacroElement);
      }

      traverse_mel = *currentMacro;
      if (traverse_mel == NULL)
        return NULL;

      stack_used = 1;
      elinfo_stack[stack_used]->fillMacroInfo(traverse_mel);
      info_stack[stack_used] = 0;

      return elinfo_stack[stack_used];
    }

    Element* el = elinfo_stack[stack_used]->getElement();

    /* go up in tree until we can go down again */
    while (stack_used > 0 &&
           (info_stack[stack_used] >= 2 || el->getFirstChild() == NULL))
    {
      stack_used--;
      el = elinfo_stack[stack_used]->getElement();
    }


    /* goto next macro element */
    if (stack_used < 1)
    {
      if (limitedToMacroElement >= 0)
        return NULL;

      currentMacro++;
      if (currentMacro == traverse_mesh->endOfMacroElements())
        return NULL;
      traverse_mel = *currentMacro;

      stack_used = 1;
      elinfo_stack[stack_used]->fillMacroInfo(traverse_mel);
      info_stack[stack_used] = 0;

      return elinfo_stack[stack_used];
    }


    /* go down tree */

    if (stack_used >= stack_size - 1)
      enlargeTraverseStack();

    int fillIthChild = info_stack[stack_used];

    info_stack[stack_used]++;
    if (traverse_fill_flag.isSet(Mesh::CALL_REVERSE_MODE))
      fillIthChild = 1 - fillIthChild;

    elinfo_stack[stack_used + 1]->fillElInfo(fillIthChild, elinfo_stack[stack_used]);

    stack_used++;

    TEST_EXIT_DBG(stack_used < stack_size)
    ("stack_size = %d too small, level = %d\n",
     stack_size, elinfo_stack[stack_used]->getLevel());

    info_stack[stack_used] = 0;

    return elinfo_stack[stack_used];
  }


  ElInfo* TraverseStack::traverseEveryElementInorder()
  {
    FUNCNAME("TraverseStack::traverseEveryElementInorder");
    ERROR_EXIT("not yet implemented\n");
    return NULL;
  }


  ElInfo* TraverseStack::traverseEveryElementPostorder()
  {
    FUNCNAME_DBG("TraverseStack::traverseEveryElementPostorder()");

    if (stack_used == 0)     /* first call */
    {
      currentMacro = traverse_mesh->firstMacroElement();
      if (limitedToMacroElement >= 0)
      {
        while ((*currentMacro)->getIndex() != limitedToMacroElement &&
               currentMacro != traverse_mesh->endOfMacroElements())
          currentMacro++;

        TEST_EXIT_DBG(currentMacro != traverse_mesh->endOfMacroElements())
        ("Coult not find macro element with index %d!\n", limitedToMacroElement);
      }

      if (currentMacro == traverse_mesh->endOfMacroElements())
        return NULL;
      traverse_mel = *currentMacro;

      stack_used = 1;
      elinfo_stack[stack_used]->fillMacroInfo(traverse_mel);
      info_stack[stack_used] = 0;

      //return(elinfo_stack[stack_used]);
    }
    else     /* don't go up on first call */
    {
      Element* el = elinfo_stack[stack_used]->getElement();

      /* go up in tree until we can go down again */          /* postorder!!! */
      while (stack_used > 0 &&
             (info_stack[stack_used] >= 3 || el->getFirstChild() == NULL))
      {
        stack_used--;
        el = elinfo_stack[stack_used]->getElement();
      }


      /* goto next macro element */
      if (stack_used < 1)
      {
        if (limitedToMacroElement >= 0)
          return NULL;

        currentMacro++;
        if (currentMacro == traverse_mesh->endOfMacroElements())
          return NULL;
        traverse_mel = *currentMacro;

        stack_used = 1;
        elinfo_stack[stack_used]->fillMacroInfo(traverse_mel);
        info_stack[stack_used] = 0;

        /*    return(elinfo_stack+stack_used); */
      }
    }
    /* go down tree */

    while (elinfo_stack[stack_used]->getElement()->getFirstChild() &&
           info_stack[stack_used] < 2)
    {
      if (stack_used >= stack_size-1)
        enlargeTraverseStack();

      int i = info_stack[stack_used];
      info_stack[stack_used]++;
      elinfo_stack[stack_used + 1]->fillElInfo(i, elinfo_stack[stack_used]);
      stack_used++;
      info_stack[stack_used] = 0;
    }

    info_stack[stack_used]++;      /* postorder!!! */

    return elinfo_stack[stack_used];
  }


  ElInfo* TraverseStack::traverseNeighbour(ElInfo* elinfo_old, int neighbour)
  {
    int dim = elinfo_old->getMesh()->getDim();
    switch(dim)
    {
    case 1:
      ERROR_EXIT("invalid dim\n");
      break;
    case 2:
      return traverseNeighbour2d(elinfo_old, neighbour);
      break;
    case 3:
      return traverseNeighbour3d(elinfo_old, neighbour);
      break;
    default:
      ERROR_EXIT("invalid dim\n");
    }
    return NULL;
  }


  ElInfo* TraverseStack::traverseNeighbour3d(ElInfo* elinfo_old, int neighbour)
  {
    FUNCNAME("TraverseStack::traverseNeighbour3d()");

    Element* el2 = NULL;
    ElInfo* elinfo2 = NULL;
    int stack2_used = 0;
    int sav_neighbour = neighbour;

    // father.neigh[coarse_nb[i][j]] == child[i - 1].neigh[j]
    static constexpr int coarse_nb[3][3][4] =
    {
      {{-2, -2, -2, -2}, {-1, 2, 3, 1}, {-1, 3, 2, 0}},
      {{-2, -2, -2, -2}, {-1, 2, 3, 1}, {-1, 2, 3, 0}},
      {{-2, -2, -2, -2}, {-1, 2, 3, 1}, {-1, 2, 3, 0}}
    };

    TEST_EXIT_DBG(stack_used > 0)("no current element\n");

    Parametric* parametric = traverse_mesh->getParametric();
    if (parametric)
      elinfo_old = parametric->removeParametricInfo(elinfo_old);

    TEST_EXIT_DBG(elinfo_old == elinfo_stack[stack_used])("invalid old elinfo\n");
    TEST_EXIT_DBG(elinfo_stack[stack_used]->getFillFlag().isSet(Mesh::FILL_NEIGH))
    ("FILL_NEIGH not set");

    Element* el = elinfo_stack[stack_used]->getElement();
    int sav_index = el->getIndex();

    // First, goto to leaf level, if necessary ...
    if ((traverse_fill_flag & Mesh::CALL_LEAF_EL).isAnySet())
    {
      if (el->getChild(0) && neighbour < 2)
      {
        if (stack_used >= stack_size - 1)
          enlargeTraverseStack();
        int i = 1 - neighbour;

        elinfo_stack[stack_used + 1]->fillElInfo(i, elinfo_stack[stack_used]);
        if (traverse_fill_flag.isSet(Mesh::CALL_REVERSE_MODE))
          info_stack[stack_used] = (i == 0 ? 2 : 1);
        else
          info_stack[stack_used] = i + 1;
        stack_used++;
        info_stack[stack_used] = 0;
        neighbour = 3;
      }
    }

    /* save information about current element and its position in the tree */
    save_traverse_mel = traverse_mel;
    save_stack_used = stack_used;

    // === First phase (see 2D). ===

    int nb = neighbour;

    while (stack_used > 1)   /* go up in tree until we can go down again */
    {
      stack_used--;
      int typ = elinfo_stack[stack_used]->getType();
      int elIsIthChild = info_stack[stack_used];
      if (traverse_fill_flag.isSet(Mesh::CALL_REVERSE_MODE) && elIsIthChild != 0)
        elIsIthChild = (elIsIthChild == 1 ? 2 : 1);

      TEST_EXIT_DBG(!elinfo_stack[stack_used + 1]->getParent() ||
                    elinfo_stack[stack_used + 1]->getParent()->getChild(elIsIthChild - 1) ==
                    elinfo_stack[stack_used + 1]->getElement())
      ("Should not happen!\n");

      nb = coarse_nb[typ][elIsIthChild][nb];

      if (nb == -1)
        break;

      TEST_EXIT_DBG(nb >= 0)("Invalid coarse_nb %d!\n", nb);
    }

    for (int i = stack_used; i <= save_stack_used; i++)
    {
      save_info_stack[i] = info_stack[i];
      *(save_elinfo_stack[i]) = *(elinfo_stack[i]);
    }
    ElInfo* old_elinfo = save_elinfo_stack[save_stack_used];
    int opp_vertex = old_elinfo->getOppVertex(neighbour);


    if (nb >= 0)
    {
      // Go to macro element neighbour.

      int i = traverse_mel->getOppVertex(nb);

      traverse_mel = traverse_mel->getNeighbour(nb);
      if (traverse_mel == NULL)
        return NULL;

      if (nb < 2 && save_stack_used > 1)
      {
        // go down one level in OLD hierarchy
        stack2_used = 2;
      }
      else
      {
        stack2_used = 1;
      }

      elinfo2 = save_elinfo_stack[stack2_used];
      el2 = elinfo2->getElement();
      stack_used = 1;
      elinfo_stack[stack_used]->fillMacroInfo(traverse_mel);
      info_stack[stack_used] = 0;
      nb = i;
    }
    else
    {
      // Goto other child.

      stack2_used = stack_used + 1;
      if (save_stack_used > stack2_used)
      {
        // go down one level in OLD hierarchy
        stack2_used++;
      }

      elinfo2 = save_elinfo_stack[stack2_used];
      el2 = elinfo2->getElement();

      if (stack_used >= stack_size - 1)
        enlargeTraverseStack();

      int i = 2 - info_stack[stack_used];
      info_stack[stack_used] = i + 1;
      int fillIthChild = i;
      if (traverse_fill_flag.isSet(Mesh::CALL_REVERSE_MODE))
        fillIthChild = 1 - fillIthChild;
      elinfo_stack[stack_used + 1]->fillElInfo(fillIthChild, elinfo_stack[stack_used]);
      stack_used++;
      info_stack[stack_used] = 0;
      nb = 0;
    }


    // === Second phase. ===

    ElInfo* elinfo = elinfo_stack[stack_used];
    el = elinfo->getElement();

    while (el->getChild(0))
    {
      if (nb < 2)
      {
        // Go down one level in hierarchy.

        if (stack_used >= stack_size - 1)
          enlargeTraverseStack();


        if (traverse_fill_flag.isSet(Mesh::CALL_REVERSE_MODE))
        {
          int t = 2 - nb;
          info_stack[stack_used] = (t == 2 ? 1 : 2);
        }
        else
        {
          info_stack[stack_used] = 2 - nb;
        }

        int fillIthChild = 1 - nb;
        elinfo_stack[stack_used + 1]->fillElInfo(fillIthChild, elinfo_stack[stack_used]);

        stack_used++;
        info_stack[stack_used] = 0;
        elinfo = elinfo_stack[stack_used];
        el = elinfo->getElement();
        nb = 3;
      }

      if (save_stack_used > stack2_used)
      {
        // `refine' both el and el2.

        TEST_EXIT_DBG(el->getChild(0))
        ("Element %d has no children!\n", el->getIndex());

        int i = 0;
        if (el->getDof(0) == el2->getDof(0))
          i = save_info_stack[stack2_used] - 1;
        else if (el->getDof(1) == el2->getDof(0))
          i = 2 - save_info_stack[stack2_used];
        else
        {
          if (traverse_mesh->associated(el->getDof(0, 0), el2->getDof(0, 0)))
            i = save_info_stack[stack2_used] - 1;
          else if (traverse_mesh->associated(el->getDof(1, 0), el2->getDof(0, 0)))
            i = 2 - save_info_stack[stack2_used];
          else
          {
            ERROR_EXIT("No common refinement edge! %d\n", traverse_fill_flag.isSet(Mesh::CALL_REVERSE_MODE));
          }
        }

        int testChild = i;
        if (traverse_fill_flag.isSet(Mesh::CALL_REVERSE_MODE))
          testChild = 1 - testChild;

        if (el->getChild(0) &&
            (el->getChild(testChild)->getDof(1) == el->getDof(nb) ||
             traverse_mesh->associated(el->getChild(testChild)->getDof(1, 0), el->getDof(nb, 0))))
          nb = 1;
        else
          nb = 2;

        info_stack[stack_used] = i + 1;

        if (stack_used >= stack_size - 1)
          enlargeTraverseStack();

        int fillIthChild = i;
        if (traverse_fill_flag.isSet(Mesh::CALL_REVERSE_MODE))
          fillIthChild = 1 - fillIthChild;
        elinfo_stack[stack_used + 1]->fillElInfo(fillIthChild, elinfo_stack[stack_used]);

        stack_used++;
        info_stack[stack_used] = 0;

        elinfo = elinfo_stack[stack_used];
        el = elinfo->getElement();

        stack2_used++;
        elinfo2 = save_elinfo_stack[stack2_used];
        el2 = elinfo2->getElement();

        if (save_stack_used > stack2_used)
        {
          const DegreeOfFreedom* dof = el2->getDof(1);

          if (dof != el->getDof(1) && dof != el->getDof(2) &&
              !traverse_mesh->associated(dof[0], el->getDof(1, 0)) &&
              !traverse_mesh->associated(dof[0], el->getDof(2, 0)))
          {
            // go down one level in OLD hierarchy
            stack2_used++;
            elinfo2 = save_elinfo_stack[stack2_used];
            el2 = elinfo2->getElement();
          }
        }
      }
      else
      {
        // Now we're done...

        elinfo = elinfo_stack[stack_used];
        el = elinfo->getElement();

        break;
      }
    }


    if (elinfo->getNeighbour(opp_vertex) != old_elinfo->getElement())
    {
      MSG(" looking for neighbour %d of element %d at %p\n",
          neighbour, old_elinfo->getElement()->getIndex(), reinterpret_cast<void*>(old_elinfo->getElement()));
      MSG(" originally: neighbour %d of element %d at %p\n",
          sav_neighbour, sav_index, reinterpret_cast<void*>(old_elinfo->getElement()));
      MSG(" got element %d at %p with opp_vertex %d neigh %d\n",
          elinfo->getElement()->getIndex(), reinterpret_cast<void*>(elinfo->getElement()),
          opp_vertex, (elinfo->getNeighbour(opp_vertex))?(elinfo->getNeighbour(opp_vertex))->getIndex():-1);
      TEST_EXIT_DBG(elinfo->getNeighbour(opp_vertex) == old_elinfo->getElement())
      ("didn't succeed !?!?!?\n");
    }


    if ((traverse_fill_flag & Mesh::CALL_EVERY_EL_POSTORDER).isAnySet())
      info_stack[stack_used] = 3;
    else if ((traverse_fill_flag & Mesh::CALL_EVERY_EL_INORDER).isAnySet())
      info_stack[stack_used] = 1;  /* ??? */

    if (elinfo)
    {
      if (parametric)
        elinfo = parametric->addParametricInfo(elinfo);
      elinfo->fillDetGrdLambda();
    }

    return elinfo;
  }


  ElInfo* TraverseStack::traverseNeighbour2d(ElInfo* elinfo_old, int neighbour)
  {
    FUNCNAME("TraverseStack::traverseNeighbour2d()");

    //     Triangle *el2 = NULL;
    //     ElInfo *elinfo2 = NULL;
    int stack2_used = 0;
    int sav_neighbour = neighbour;

    // father.neigh[coarse_nb[i][j]] == child[i-1].neigh[j]
    // TODO: REMOVE STATIC
    static int coarse_nb[3][3] = {{-2, -2, -2}, {2, -1, 1}, {-1, 2, 0}};

    TEST_EXIT_DBG(stack_used > 0)("no current element");

    Parametric* parametric = traverse_mesh->getParametric();
    if (parametric)
      elinfo_old = parametric->removeParametricInfo(elinfo_old);

    TEST_EXIT_DBG(elinfo_old == elinfo_stack[stack_used])("invalid old elinfo");

    elinfo_stack[stack_used]->testFlag(Mesh::FILL_NEIGH);
    Triangle* el =
      dynamic_cast<Triangle*>(const_cast<Element*>(elinfo_stack[stack_used]->getElement()));
    int sav_index = el->getIndex();
    Triangle* sav_el = el;

    /* first, goto to leaf level, if necessary... */
    if (!(el->isLeaf()) && neighbour < 2)
    {

      if (stack_used >= static_cast<int>(elinfo_stack.size()) - 1)
        enlargeTraverseStack();

      // If we should search for neighbour 0, take second child, if for
      // neighbour 1, take the first child.
      int i = 1 - neighbour;

      elinfo_stack[stack_used + 1]->fillElInfo(i, elinfo_stack[stack_used]);
      if (traverse_fill_flag.isSet(Mesh::CALL_REVERSE_MODE))
        info_stack[stack_used] = (i == 0 ? 2 : 1);
      else
        info_stack[stack_used] = i + 1;

      stack_used++;
      neighbour = 2;
    }

    /* save information about current element and its position in the tree */
    save_traverse_mel = traverse_mel;
    save_stack_used = stack_used;
    for (int i = 0; i <= stack_used; i++)
    {
      save_info_stack[i] = info_stack[i];
      (*(save_elinfo_stack[i])) = (*(elinfo_stack[i]));
    }
    ElInfo* old_elinfo = save_elinfo_stack[stack_used];
    int opp_vertex = old_elinfo->getOppVertex(neighbour);


    /****************************************************************************/
    /* First phase: go up in tree until we can go down again.                   */
    /*                                                                          */
    /* During this first phase, nb is the neighbour index which points from an  */
    /* element of the OLD hierarchy branch to the new branch                    */
    /****************************************************************************/

    int nb = neighbour;

    while (stack_used > 1)
    {
      stack_used--;
      int elIsIthChild = info_stack[stack_used];

      if (traverse_fill_flag.isSet(Mesh::CALL_REVERSE_MODE) && elIsIthChild != 0)
        elIsIthChild = (elIsIthChild == 1 ? 2 : 1);

      TEST_EXIT_DBG(!elinfo_stack[stack_used + 1]->getParent() ||
                    elinfo_stack[stack_used + 1]->getParent()->getChild(elIsIthChild - 1) ==
                    elinfo_stack[stack_used + 1]->getElement())
      ("Should not happen!\n");

      nb = coarse_nb[elIsIthChild][nb];
      if (nb == -1)
        break;

      TEST_EXIT_DBG(nb >= 0)("invalid coarse_nb %d\n",nb);
    }


    /****************************************************************************/
    /* Now, goto neighbouring element at the local hierarchy entry              */
    /* This is either a macro element neighbour or the other child of parent.   */
    /* initialize nb for second phase (see below)                               */
    /****************************************************************************/

    if (nb >= 0)
    {
      // Go to macro element neighbour.

      if (nb < 2 && save_stack_used > 1)
        stack2_used = 2;           /* go down one level in OLD hierarchy */
      else
        stack2_used = 1;

      //       elinfo2 = save_elinfo_stack[stack2_used];
      //       el2 = dynamic_cast<Triangle*>(const_cast<Element*>(elinfo2->getElement()));

      int i = traverse_mel->getOppVertex(nb);
      traverse_mel = traverse_mel->getNeighbour(nb);
      if (traverse_mel == NULL)
        return NULL;
      nb = i;

      stack_used = 1;
      elinfo_stack[stack_used]->fillMacroInfo(const_cast<MacroElement*>(traverse_mel));
      info_stack[stack_used] = 0;
    }
    else
    {
      // Go to other child.

      stack2_used = stack_used + 1;
      if (save_stack_used > stack2_used)
        stack2_used++;               /* go down one level in OLD hierarchy */

      //       elinfo2 = save_elinfo_stack[stack2_used];
      //       el2 = dynamic_cast<Triangle*>(const_cast<Element*>(elinfo2->getElement()));

      if (stack_used >= stack_size - 1)
        enlargeTraverseStack();

      TEST_EXIT_DBG(info_stack[stack_used] == 1 || info_stack[stack_used] == 2)
      ("Should not happen!\n");

      int fillIthChild = -1;
      if (info_stack[stack_used] == 1)
      {
        info_stack[stack_used] = 2;
        fillIthChild = 1;
        nb = 0;
      }
      else
      {
        info_stack[stack_used] = 1;
        fillIthChild = 0;
        nb = 1;
      }

      if (traverse_fill_flag.isSet(Mesh::CALL_REVERSE_MODE))
      {
        fillIthChild = 1 - fillIthChild;
        nb = 1 - nb;
      }
      elinfo_stack[stack_used + 1]->fillElInfo(fillIthChild, elinfo_stack[stack_used]);
      stack_used++;
    }

    /****************************************************************************/
    /* Second phase: go down in a new hierarchy branch until leaf level.        */
    /* Now, nb is the neighbour index which points from an element of the       */
    /* new hierarchy branch to the OLD branch.                                  */
    /****************************************************************************/

    ElInfo* elinfo = elinfo_stack[stack_used];
    el = dynamic_cast<Triangle*>(const_cast<Element*>(elinfo->getElement()));

    while (el->getFirstChild())
    {
      if (nb < 2)
      {
        // Go down one level in hierarchy.

        if (stack_used >= stack_size - 1)
          enlargeTraverseStack();

        if (traverse_fill_flag.isSet(Mesh::CALL_REVERSE_MODE))
        {
          int t = 2 - nb;
          info_stack[stack_used] = (t == 2 ? 1 : 2);
        }
        else
        {
          info_stack[stack_used] = 2 - nb;
        }

        int fillIthChild = 1 - nb;
        elinfo_stack[stack_used + 1]->fillElInfo(fillIthChild, elinfo_stack[stack_used]);

        stack_used++;
        nb = 2;
      }

      if (save_stack_used > stack2_used)
      {
        // `refine' both el and el2.

        TEST_EXIT_DBG(el->getFirstChild())("invalid new refinement?");

        // Use child i, neighbour of el2->child[nb - 1].
        int i = 2 - save_info_stack[stack2_used];
        TEST_EXIT_DBG(i < 2)("invalid OLD refinement?");
        info_stack[stack_used] = i + 1;

        int fillIthChild = i;
        nb = i;

        if (traverse_fill_flag.isSet(Mesh::CALL_REVERSE_MODE))
        {
          fillIthChild = 1 - fillIthChild;
          nb = 1 - i;
        }

        elinfo_stack[stack_used + 1]->fillElInfo(fillIthChild, elinfo_stack[stack_used]);

        stack_used++;

        elinfo = elinfo_stack[stack_used];
        el = dynamic_cast<Triangle*>(const_cast<Element*>(elinfo->getElement()));

        stack2_used++;
        if (save_stack_used > stack2_used)
          stack2_used++;                /* go down one level in OLD hierarchy */

        // 	elinfo2 = save_elinfo_stack[stack2_used];
        // 	el2 = dynamic_cast<Triangle*>(const_cast<Element*>(elinfo2->getElement()));
      }
      else
      {
        // Now we're done...

        elinfo = elinfo_stack[stack_used];
        el = dynamic_cast<Triangle*>(const_cast<Element*>(elinfo->getElement()));
      }
    }

    if (elinfo->getNeighbour(opp_vertex) != old_elinfo->getElement())
    {
      MSG(" looking for neighbour %d of element %d at %8X\n",
          neighbour, old_elinfo->getElement()->getIndex(), old_elinfo->getElement());
      MSG(" originally: neighbour %d of element %d at %8X\n",
          sav_neighbour, sav_index, sav_el);
      MSG(" got element %d at %8X with opp_vertex %d neigh %d\n",
          elinfo->getElement()->getIndex(), elinfo->getElement(), opp_vertex,
          elinfo->getNeighbour(opp_vertex) ? elinfo->getNeighbour(opp_vertex)->getIndex() : -1);
      TEST_EXIT_DBG(elinfo->getNeighbour(opp_vertex) == old_elinfo->getElement())
      ("didn't succeed !?!?!?");
    }

    if (elinfo->getElement()->getFirstChild())
    {
      MSG(" looking for neighbour %d of element %d at %8X\n",
          neighbour, old_elinfo->getElement()->getIndex(), old_elinfo->getElement());
      MSG(" originally: neighbour %d of element %d at %8X\n",
          sav_neighbour, sav_index, sav_el);
      MSG(" got element %d at %8X with opp_vertex %d neigh %d\n",
          elinfo->getElement()->getIndex(), elinfo->getElement(), opp_vertex,
          elinfo->getNeighbour(opp_vertex)->getIndex());
      ERROR_EXIT("got no leaf element\n");
    }

    if (elinfo)
    {
      if (parametric)
        elinfo = parametric->addParametricInfo(elinfo);

      elinfo->fillDetGrdLambda();
    }

    return elinfo;
  }


  void TraverseStack::update()
  {
    FUNCNAME_DBG("TraverseStack::update()");

    TEST_EXIT_DBG(traverse_mesh->getDim() == 3)
    ("Update only in 3d, mesh is d = %d\n", traverse_mesh->getDim());

    for (int i = stack_used; i > 0; i--)
      dynamic_cast<ElInfo3d*>(elinfo_stack[i])->update();
  }


  void TraverseStack::fillRefinementPath(ElInfo& elInfo, const ElInfo& upperElInfo)
  {
    int levelDif = elinfo_stack[stack_used]->getLevel() - upperElInfo.getLevel();
    unsigned long rPath = 0;

    for (int i = 1; i <= levelDif; i++)
      if (elinfo_stack[stack_used - levelDif + i]->getIChild())
        rPath = rPath | (1 << (i - 1));

    elInfo.setRefinementPath(rPath);
    elInfo.setRefinementPathLength(levelDif);
  }

} // end namespace AMDiS
