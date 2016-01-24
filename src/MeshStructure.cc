#include "Debug.h"
#include "DOFVector.h"
#include "Element.h"
#include "ElementDofIterator.h"
#include "ElInfo.h"
#include "MeshStructure.h"
#include "MeshStructure_ED.h"
#include "Mesh.h"
#include "RefinementManager.h"
#include "Traverse.h"
#include "MacroElement.h"

using namespace std;

namespace AMDiS
{
  void MeshStructure::insertElement(bool isLeaf)
  {
    // overflow? -> next index
    if (pos >= structureSize)
    {
      code.push_back(currentCode);
      pos = 0;
      currentCode = 0;
    }

    // insert element in binary code
    if (!isLeaf)
    {
      uint64_t one = 1;
      currentCode += (one << pos);
    }

    pos++;
    nElements++;
  }


  void MeshStructure::clear()
  {
    currentCode = 0;
    code.resize(0);
    pos = 0;
    nElements = 0;
    currentElement = 0;
  }


  void MeshStructure::init(Mesh* mesh, int macroElIndex)
  {
    clear();

    TraverseStack stack;

    ElInfo* elInfo;
    if (macroElIndex == -1)
      elInfo = stack.traverseFirst(mesh, -1, Mesh::CALL_EVERY_EL_PREORDER);
    else
      elInfo = stack.traverseFirstOneMacro(mesh, macroElIndex, -1,
                                           Mesh::CALL_EVERY_EL_PREORDER);

    while (elInfo)
    {
      insertElement(elInfo->getElement()->isLeaf());
      elInfo = stack.traverseNext(elInfo);
    }

    commit();
  }


  void MeshStructure::init(BoundaryObject& bound, Element* element)
  {
    FUNCNAME("MeshStructure::init()");

    Element* el = (element == NULL) ? bound.el : element;

    TEST_EXIT_DBG(el)("No element!\n");

    clear();

    int s1 = el->getSubObjOfChild(0, bound.subObj, bound.ithObj, bound.elType);
    int s2 = el->getSubObjOfChild(1, bound.subObj, bound.ithObj, bound.elType);

    TEST_EXIT(s1 != -1 || s2 != -1)("This should not happen!\n");

    if (debugMode)
    {
      MSG("addAlondSide(%d, %d, %d, %d)\n",
          bound.elIndex, bound.ithObj, bound.elType, bound.reverseMode);
      MSG("Element is leaf: %d\n", el->isLeaf());
      MSG("s1 = %d    s2 = %d\n", s1, s2);
    }

    /*    switch (bound.subObj)
    {
      case EDGE:
    */
    if (!el->isLeaf())
    {
      if (s1 == -1)
        addAlongSide(el->getSecondChild(), bound.subObj, s2,
                     el->getChildType(bound.elType), bound.reverseMode);
      else if (s2 == -1)
        addAlongSide(el->getFirstChild(), bound.subObj, s1,
                     el->getChildType(bound.elType), bound.reverseMode);
      else
        addAlongSide(el, bound.subObj, bound.ithObj, bound.elType, bound.reverseMode);
    }
    /*	break;
    case FACE:
    addAlongSide(el, bound.subObj, bound.ithObj, bound.elType, bound.reverseMode);
    break;
      default:
    ERROR_EXIT("What is this?\n");
    }
      */

    commit();
  }


  void MeshStructure::addAlongSide(Element* el, GeoIndex subObj, int ithObj,
                                   int elType, bool reverseOrder)
  {
    FUNCNAME("MeshStructure::addAlongSide()");

    if (debugMode)
    {
      MSG("addAlondSide(%d, %d, %d, %d)\n",
          el->getIndex(), ithObj, elType, reverseOrder);
      MSG("Element is leaf: %d\n", el->isLeaf());
    }

    insertElement(el->isLeaf());

    if (!el->isLeaf())
    {
      int s1 = el->getSubObjOfChild(0, subObj, ithObj, elType);
      int s2 = el->getSubObjOfChild(1, subObj, ithObj, elType);

      if (debugMode)
      {
        MSG("Child index %d  %d\n",
            el->getFirstChild()->getIndex(),
            el->getSecondChild()->getIndex());
        MSG("s1 = %d    s2 = %d\n", s1, s2);
        MSG("   \n");
      }

      if (!reverseOrder)
      {
        if (s1 != -1)
          addAlongSide(el->getFirstChild(), subObj, s1,
                       el->getChildType(elType), reverseOrder);
        if (s2 != -1)
          addAlongSide(el->getSecondChild(), subObj, s2,
                       el->getChildType(elType), reverseOrder);
      }
      else
      {
        if (s2 != -1)
          addAlongSide(el->getSecondChild(), subObj, s2,
                       el->getChildType(elType), reverseOrder);
        if (s1 != -1)
          addAlongSide(el->getFirstChild(), subObj, s1,
                       el->getChildType(elType), reverseOrder);
      }
    }
  }


  void MeshStructure::reset()
  {
    currentIndex = 0;
    pos = 0;
    currentElement = 0;

    if (code.size() > 0)
      currentCode = code[0];
    else
      currentCode = 0;
  }


  bool MeshStructure::nextElement(MeshStructure* insert)
  {
    FUNCNAME_DBG("MeshStructure::nextElement()");

    if (insert)
      insert->insertElement(isLeafElement());

    pos++;
    currentElement++;

    if (currentElement >= nElements)
      return false;

    if (pos >= structureSize)
    {
      currentIndex++;
      TEST_EXIT_DBG(currentIndex < static_cast<int>(code.size()))
      ("End of structure reached!\n");
      pos = 0;
      currentCode = code[currentIndex];
    }
    else
    {
      currentCode >>= 1;
    }

    return true;
  }


  int MeshStructure::lookAhead(unsigned int n)
  {
    int returnValue = 0;

    int tmp_pos = pos;
    int tmp_currentElement = currentElement;
    int tmp_currentIndex = currentIndex;
    uint64_t tmp_currentCode = currentCode;

    for (unsigned int i = 0; i < n; i++)
    {
      if (nextElement() == false)
      {
        returnValue = -1;
        break;
      }
    }

    if (returnValue != -1)
      returnValue = static_cast<int>(!isLeafElement());

    pos = tmp_pos;
    currentElement = tmp_currentElement;
    currentIndex = tmp_currentIndex;
    currentCode = tmp_currentCode;

    return returnValue;
  }


  bool MeshStructure::skipBranch(MeshStructure* insert)
  {
    FUNCNAME_DBG("MeshStructure::skipBranch()");

    if (isLeafElement())
    {
      return nextElement(insert);
    }
    else
    {
      bool cont = nextElement(insert);
      cont = skipBranch(insert); // left branch
      TEST_EXIT_DBG(cont)("Invalid structure!\n");
      cont = skipBranch(insert); // righ branch
      return cont;
    }
  }


  void MeshStructure::merge(MeshStructure* structure1,
                            MeshStructure* structure2,
                            MeshStructure* result)
  {
    FUNCNAME_DBG("MeshStructure::merge()");

    result->clear();
    structure1->reset();
    structure2->reset();

    bool cont = true;
    while (cont)
    {
      bool cont1;
#if DEBUG != 0
      bool cont2;
#endif
      if (structure1->isLeafElement() == structure2->isLeafElement())
      {
        cont1 = structure1->nextElement(result);
#if DEBUG != 0
        cont2 = structure2->nextElement();
#endif
      }
      else
      {
        if (structure1->isLeafElement())
        {
          cont1 = structure1->nextElement();
#if DEBUG != 0
          cont2 = structure2->skipBranch(result);
#endif
        }
        else
        {
          cont1 = structure1->skipBranch(result);
#if DEBUG != 0
          cont2 = structure2->nextElement();
#endif
        }
      }
      TEST_EXIT_DBG(cont1 == cont2)("Structures don't match!\n");
      cont = cont1;
    }

    result->commit();
  }


  void MeshStructure::fitMeshToStructure(Mesh* mesh,
                                         RefinementManager* manager,
                                         bool debugMode,
                                         int macroElIndex,
                                         bool ignoreFinerMesh)
  {
    FUNCNAME("MeshStructure::fitMeshToStructure()");

    bool cont = true;

    // decorate leaf data
    reset();
    TraverseStack stack;
    ElInfo* elInfo = NULL;
    if (macroElIndex == -1)
      elInfo = stack.traverseFirst(mesh, -1, Mesh::CALL_EVERY_EL_PREORDER);
    else
      elInfo = stack.traverseFirstOneMacro(mesh, macroElIndex, -1, Mesh::CALL_EVERY_EL_PREORDER);


    while (elInfo)
    {
      Element* element = elInfo->getElement();

      TEST_EXIT(cont)("unexpected structure code end!\n");

      if (isLeafElement())
      {
        if (ignoreFinerMesh && !element->isLeaf())
        {
          int level = elInfo->getLevel();
          while (elInfo && level >= elInfo->getLevel())
            elInfo = stack.traverseNext(elInfo);
        }
        else
        {
          TEST_EXIT(element->isLeaf())
          ("Mesh is finer than strucutre code! (Element index: %d Macro element index: %d)\n",
           element->getIndex(), elInfo->getMacroElement()->getIndex());
        }
      }

      TEST_EXIT_DBG(element)("Should not happen!\n");

      if (element->isLeaf() && !isLeafElement())
      {
        MeshStructure* structure = new MeshStructure();
        cont = skipBranch(structure);
        structure->commit();

        MeshStructure_ED* elData = new MeshStructure_ED(element->getElementData());
        elData->setStructure(structure);
        element->setElementData(elData);
      }
      else
      {
        cont = nextElement();
      }

      if (elInfo)
        elInfo = stack.traverseNext(elInfo);
    }

    // refine mesh
    bool finished = true;

    do
    {
      finished = true;
      if (macroElIndex == -1)
        elInfo = stack.traverseFirst(mesh, -1, Mesh::CALL_LEAF_EL);
      else
        elInfo = stack.traverseFirstOneMacro(mesh, macroElIndex, -1, Mesh::CALL_LEAF_EL);
      while (elInfo)
      {
        Element* element = elInfo->getElement();
        if (element->getElementData(AMDIS_MESH_STRUCTURE) != NULL)
        {
          element->setMark(1);
          finished = false;
        }
        else
        {
          element->setMark(0);
        }
        elInfo = stack.traverseNext(elInfo);
      }

      if (!finished)
      {
#if (DEBUG != 0)
        int oldMeshIndex = mesh->getChangeIndex();
#endif

        if (macroElIndex == -1)
          manager->refineMesh(mesh);
        else
          manager->refineMacroElement(mesh, macroElIndex);

#if (DEBUG != 0)
        TEST_EXIT(oldMeshIndex != mesh->getChangeIndex())
        ("Mesh has not been changed by refinement procedure!\n");
#endif
      }
    }
    while (!finished);
  }


  string MeshStructure::toStr(bool resetCode)
  {
    std::stringstream oss;

    if (empty())
    {
      oss << "-" << std::endl;
    }
    else
    {
      if (resetCode)
        reset();

      bool cont = true;
      while (cont)
      {
        if (isLeafElement())
          oss << "0";
        else
          oss << "1";

        cont = nextElement();
      }
    }

    return oss.str();
  }


  void MeshStructure::print(bool resetCode)
  {
    FUNCNAME("MeshStructure::print()");

    string str = toStr(resetCode);

    if (str.length() < 255)
    {
      MSG("Mesh structure code: %s\n", str.c_str());
    }
    else
    {
#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
      std::cout << "[" << MPI::COMM_WORLD.Get_rank() << "]                Mesh structure code: " << str << "\n";
#else
      std::cout << "                Mesh structure code: " << str << "\n";
#endif
    }
  }


  bool MeshStructure::compare(MeshStructure& other)
  {
    return (other.getCode() == code);
  }


  void MeshStructure::getMeshStructureValues(int macroElIndex,
      DOFVector<double> const* vec,
      std::vector<double>& values,
      bool withElIndex)
  {
    FUNCNAME_DBG("MeshStructure::getMeshStructureValues()");

    TEST_EXIT_DBG(vec)("No DOFVector defined!\n");

    const FiniteElemSpace* feSpace = vec->getFeSpace();
    Mesh* mesh = feSpace->getMesh();
    int nVertexPreDofs = feSpace->getAdmin()->getNumberOfPreDofs(VERTEX);
    bool feSpaceHasNonVertexDofs = (feSpace->getBasisFcts()->getDegree() > 1);
    values.clear();

    // In debug mode we add the macro element index to the value code.
    if (withElIndex)
      values.push_back(static_cast<double>(macroElIndex));

    ElementDofIterator elDofIter(feSpace);
    TraverseStack stack;
    ElInfo* elInfo = stack.traverseFirstOneMacro(mesh, macroElIndex, -1,
                     Mesh::CALL_EVERY_EL_PREORDER);
    while (elInfo)
    {
      // For the macro element the mesh structure code stores all vertex values.
      if (elInfo->getLevel() == 0)
        for (int i = 0; i < mesh->getGeo(VERTEX); i++)
          values.push_back((*vec)[elInfo->getElement()->getDof(i, nVertexPreDofs)]);

      if (!elInfo->getElement()->isLeaf())
      {
        // If no leaf element store the value of the "new" DOF that is created
        // by bisectioning of this element.

        DegreeOfFreedom dof0 =
          elInfo->getElement()->getChild(0)->getDof(mesh->getDim(), nVertexPreDofs);
        values.push_back((*vec)[dof0]);
      }
      else
      {
        // If leaf element store all non vertex values of this element, thus
        // only relevant for higher order basis functions.

        if (feSpaceHasNonVertexDofs)
        {
          elDofIter.reset(elInfo->getElement());
          do
          {
            if (elDofIter.getPosIndex() != VERTEX)
              values.push_back((*vec)[elDofIter.getDof()]);
          }
          while (elDofIter.next());
        }
      }

      elInfo = stack.traverseNext(elInfo);
    }
  }


  void MeshStructure::setMeshStructureValues(int macroElIndex,
      DOFVector<double>* vec,
      const std::vector<double>& values,
      bool withElIndex)
  {
    FUNCNAME_DBG("MeshStructure::setMeshStructureValues()");

    TEST_EXIT_DBG(vec)("No DOFVector defined!\n");

    const FiniteElemSpace* feSpace = vec->getFeSpace();
    Mesh* mesh = feSpace->getMesh();
    bool feSpaceHasNonVertexDofs = (feSpace->getBasisFcts()->getDegree() > 1);
    int nVertexPreDofs = feSpace->getAdmin()->getNumberOfPreDofs(VERTEX);
    unsigned int counter = 0;

    if (withElIndex)
    {
      TEST_EXIT(static_cast<int>(values[0]) == macroElIndex)
      ("Value structure code was created for macro element index %d, but should be set to macro element index %d\n",
       static_cast<int>(values[0]), macroElIndex);
      counter++;
    }

    TEST_EXIT_DBG(static_cast<int>(values.size()) >= mesh->getGeo(VERTEX))
    ("Should not happen!\n");

    ElementDofIterator elDofIter(feSpace);
    TraverseStack stack;
    ElInfo* elInfo = stack.traverseFirstOneMacro(mesh, macroElIndex, -1,
                     Mesh::CALL_EVERY_EL_PREORDER);
    while (elInfo)
    {
      // For the macro element all vertex nodes are set first.
      if (elInfo->getLevel() == 0)
        for (int i = 0; i < mesh->getGeo(VERTEX); i++)
          (*vec)[elInfo->getElement()->getDof(i, nVertexPreDofs)] =
            values[counter++];

      if (!elInfo->getElement()->isLeaf())
      {
        // If no leaf element set the value of the "new" DOF that is created
        // by bisectioning of this element.
        TEST_EXIT_DBG(counter < values.size())("Should not happen!\n");

        (*vec)[elInfo->getElement()->getChild(0)->getDof(mesh->getDim(), nVertexPreDofs)] =
          values[counter++];
      }
      else
      {
        // On leaf elements set all non vertex values (thus DOFs of higher order
        // basis functions).

        if (feSpaceHasNonVertexDofs)
        {
          elDofIter.reset(elInfo->getElement());
          do
          {
            if (elDofIter.getPosIndex() != VERTEX)
              (*vec)[elDofIter.getDof()] = values[counter++];
          }
          while (elDofIter.next());
        }
      }

      elInfo = stack.traverseNext(elInfo);
    }

    TEST_EXIT_DBG(values.size() == counter)
    ("Should not happen! values size %d, counter %d\n", values.size(), counter);
  }

} // end namespace AMDiS
