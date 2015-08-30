/** \file Expressions.hh */

#include <Traverse.h>
#include <DualTraverse.h>

#include <traits/basic.hpp>

namespace AMDiS
{
  template <class M>
  inline Value_t<M>
  integrate(BaseTerm<M> const& term_base, Mesh* mesh_opt)
  {
    M term = term_base.sub();
    using TOut = Assign_t<Value_t<M>>;

    std::set<const FiniteElemSpace*> feSpaces;
    term.insertFeSpaces(feSpaces);

    TEST_EXIT(mesh_opt || !feSpaces.empty())
    ("The expression has no reference to a mesh. Pass a mesh explicitely as second argument to integrate(...)!\n");
    Mesh* mesh = mesh_opt ? mesh_opt : (*feSpaces.begin())->getMesh();

    int deg = term.getDegree();
    int dim = mesh->getDim();
    Quadrature* quad = Quadrature::provideQuadrature(dim, deg);

    TOut value;
    nullify(value);

    Flag traverseFlag = Mesh::CALL_LEAF_EL | Mesh::FILL_COORDS | 
                        Mesh::FILL_GRD_LAMBDA | Mesh::FILL_DET;
    TraverseStack stack;
    ElInfo* elInfo = stack.traverseFirst(mesh, -1, traverseFlag);
    while (elInfo)
    {
      term.initElement(elInfo, NULL, quad);
      TOut tmp;
      nullify(tmp);
      for (int iq = 0; iq < quad->getNumPoints(); iq++)
      {
        tmp += quad->getWeight(iq) * term[iq];
      }
      value += tmp * elInfo->getDet();

      elInfo = stack.traverseNext(elInfo);
    }

#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
    Parallel::mpi::globalAdd(value);
#endif

    return value;
  }

  namespace detail
  {
    template <class M, class Functor>
    inline Value_t<M>
    accumulate(BaseTerm<M> const& term_base, Functor f, Value_t<M> value0)
    {
      M term = term_base.sub();
      std::set<const FiniteElemSpace*> feSpaces;
      term.insertFeSpaces(feSpaces);

      TEST_EXIT(!feSpaces.empty())
      ("The expression must contain a DOFVector or FeSpace depended value!\n");
      const FiniteElemSpace* feSpace0 = *feSpaces.begin();
      Mesh* mesh = feSpace0->getMesh();

      const BasisFunction* basisFcts = feSpace0->getBasisFcts();
      int nBasisFcts = basisFcts->getNumber();

      DOFVector<bool> assigned(feSpace0, "assigned");
      assigned.set(false);

      std::vector<DegreeOfFreedom> localIndices(nBasisFcts);
      Flag traverseFlag = Mesh::CALL_LEAF_EL | Mesh::FILL_COORDS | Mesh::FILL_GRD_LAMBDA;
      TraverseStack stack;
      ElInfo* elInfo = stack.traverseFirst(mesh, -1, traverseFlag);

      while (elInfo)
      {
        term.initElement(elInfo, NULL, NULL, basisFcts);
        basisFcts->getLocalIndices(elInfo->getElement(), feSpace0->getAdmin(), localIndices);

        for (int i = 0; i < nBasisFcts; i++)
        {
          if (!assigned[localIndices[i]])
          {
            value0 = f(value0, term[i]);
            assigned[localIndices[i]] = true;
          }
        }
        elInfo = stack.traverseNext(elInfo);
      }

      return value0;
    }
    
    
    template <class T, class M>
    inline void transformDOF(BaseTerm<M> const& term_base, DOFVector<T>* result)
    {
      M term = term_base.sub();
      
      using TOut = Assign_t<Value_t<M>>;
      static_assert( traits::IsConvertible<TOut, T>::value,
        "ValueType of expression not convertible to ValueType of DOFVector!" );

      std::set<const FiniteElemSpace*> feSpaces;
      term.insertFeSpaces(feSpaces);

      Mesh* mesh = result->getFeSpace()->getMesh();
      //     if (feSpaces.size() > 0 && mesh != (*feSpaces.begin())->getMesh())
      //       return transformDOF_mm(term, result);


      DOFVector<TOut> temp(result->getFeSpace(), "temp");
      DOFVector<int> assigned(result->getFeSpace(), "assigned");

      const FiniteElemSpace* resultFeSpace = temp.getFeSpace();
      const BasisFunction* basisFcts = resultFeSpace->getBasisFcts();
      int nBasisFcts = basisFcts->getNumber();

      std::vector<DegreeOfFreedom> localIndices(nBasisFcts);
      Flag traverseFlag = Mesh::CALL_LEAF_EL | Mesh::FILL_COORDS | Mesh::FILL_GRD_LAMBDA;
      TraverseStack stack;
      ElInfo* elInfo = stack.traverseFirst(mesh, -1, traverseFlag);
      term.initElement(elInfo, NULL, NULL, basisFcts);

      TOut tmp(term[0]);
      nullify(tmp);
      assigned.set(0);
      temp.set(tmp);

      while (elInfo)
      {
        term.initElement(elInfo, NULL, NULL, basisFcts);
        basisFcts->getLocalIndices(elInfo->getElement(), resultFeSpace->getAdmin(), localIndices);

        for (int i = 0; i < nBasisFcts; i++)
        {
          temp[localIndices[i]] += term[i];
          assigned[localIndices[i]]++;
        }
        elInfo = stack.traverseNext(elInfo);
      }

#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
      Parallel::MeshDistributor::globalMeshDistributor->synchAddVector(temp);
      Parallel::MeshDistributor::globalMeshDistributor->synchAddVector(assigned);
#endif

      DOFIterator<TOut> tempIter(&temp, USED_DOFS);
      DOFIterator<T> resultIter(result, USED_DOFS);
      DOFIterator<int> assignedIter(&assigned, USED_DOFS);
      for (tempIter.reset(), resultIter.reset(), assignedIter.reset();
          !resultIter.end();
          ++tempIter, ++resultIter, ++assignedIter)
      {
        *resultIter = (*tempIter);
        *resultIter/= (*assignedIter);
      }
    }
    
  } // end namespace detail
} // end namespace AMDiS
