#pragma once

// AMDiS includes
#include "AMDiS_fwd.hpp"
#include "BasisFunction.hpp"
#include "DOFVectorBase.hpp"
#include "SubAssembler.hpp"
#include "expressions/BaseTerms.hpp"
#include "expressions/LazyOperatorTerm.hpp"

namespace AMDiS
{
  // a dummy name assigned by default to each DOFVector term
  namespace name
  {
    struct u {};
    template <char> struct v {};
  }

  /// Expressions that extracts the values of a DOFVector at QPs
  template <class Vector, class Name>
  struct ValueOf;

  template <class T, class Name>
  struct ValueOf<DOFVector<T>, Name>
    : public ShapedTerm_t<T, ValueOf<DOFVector<T>, Name>>,
      public LazyOperatorTermBase
  {
    using value_type = T;
    using id         = Name;

    constexpr ValueOf(DOFVectorBase<T> const& vector) : vecDV(&vector) {}
    constexpr ValueOf(DOFVectorBase<T> const* vector) : vecDV(vector) {}

    template <class List>
    void insertFeSpaces(List& feSpaces) const
    {
      feSpaces.insert(vecDV->getFeSpace());
    }

    int getDegree() const
    {
      return vecDV->getFeSpace()->getBasisFcts()->getDegree();
    }

    void initElement(ElInfo const* elInfo,
                     SubAssembler* subAssembler,
                     Quadrature* quad,
                     BasisFunction const* basisFct = NULL)
    {
      if (subAssembler)
        subAssembler->getVectorAtQPs(vecDV, elInfo, quad, values);
      else if (quad)
        vecDV->getVecAtQPs(elInfo, quad, NULL, values);
      else if (basisFct)
      {
        const BasisFunction* localBasisFct = vecDV->getFeSpace()->getBasisFcts();

        // get coefficients of DOFVector
        DenseVector<T> coeff((size_t)localBasisFct->getNumber());
        vecDV->getLocalVector(elInfo->getElement(), coeff);

        // eval basisfunctions of DOFVector at coords of given basisFct
        const int nBasisFct = basisFct->getNumber();
        values.change_dim(nBasisFct);
        for (int i = 0; i < nBasisFct; i++)
          values[i] = localBasisFct->evalUh(*basisFct->getCoords(i), coeff);
      }
    }

    /// eval at point with index iq, initialized in \ref initElement
    value_type evalAtIdx(int iq) const
    {
      return values[iq];
    }

    /// eval at point with coordinates x
    value_type operator()(WorldVector<double> const& x) const
    {
      return (*vecDV)(x);
    }

    std::string str() const
    {
      return std::string("value(") + vecDV->getName() + ")";
    }

  private:
    DOFVectorBase<T> const* vecDV;

    DenseVector<T> values;
  };


  // ---------------------------------------------------------------------------


  /// Expressions that extracts the values of a DOFVector at QPs
  template <class Vector, class Name, class = void>
  struct GradientOf;

  template <class T, class Name>
  struct GradientOf<DOFVector<T>, Name>
    : public ShapedTerm_t<Gradient_t<T>, GradientOf<DOFVector<T>, Name>>,
      public LazyOperatorTermBase
  {
    using value_type = Gradient_t<T>;
    using id         = Name;

    constexpr GradientOf(DOFVectorBase<T> const& vector) : vecDV(&vector) {}
    constexpr GradientOf(DOFVectorBase<T> const* vector) : vecDV(vector) {}

    template <class List>
    void insertFeSpaces(List& feSpaces) const
    {
      feSpaces.insert(vecDV->getFeSpace());
    }

    int getDegree() const
    {
      return vecDV->getFeSpace()->getBasisFcts()->getDegree();
    }

    void initElement(ElInfo const* elInfo,
                     SubAssembler* subAssembler,
                     Quadrature* quad,
                     BasisFunction const* basisFct = NULL)
    {
      if (subAssembler)
        subAssembler->getGradientsAtQPs(vecDV, elInfo, quad, values);
      else if (quad)
        vecDV->getGrdAtQPs(elInfo, quad, NULL, values);
      else if (basisFct)
      {
        BasisFunction const* localBasisFct = vecDV->getFeSpace()->getBasisFcts();

        // get coefficients of DOFVector
        DenseVector<T> coeff{(size_t)localBasisFct->getNumber()};
        vecDV->getLocalVector(elInfo->getElement(), coeff);

        // eval basisfunctions of DOFVector at coords of given basisFct
        const int nBasisFct = basisFct->getNumber();
        values.change_dim(nBasisFct);
        auto& grdLambda = elInfo->getGrdLambda();
        for (int i = 0; i < nBasisFct; i++)
          localBasisFct->evalGrdUh(*basisFct->getCoords(i), grdLambda, coeff, values[i]);
      }
    }

    value_type evalAtIdx(int iq) const
    {
      return values[iq];
    }

    value_type operator()(WorldVector<double> const& x) const
    {
      ElInfo* elInfo = NULL;
      DimVec<double> lambda(vecDV->getMesh()->getDim());
      bool found = const_cast<Mesh*>(vecDV->getMesh())->findElInfoAtPoint(x, elInfo, lambda, NULL, NULL, NULL);

      value_type result;
      if (found && elInfo != NULL)
      {
        auto& grdLambda = elInfo->getGrdLambda();
        BasisFunction const* localBasisFct = vecDV->getFeSpace()->getBasisFcts();
        // get coefficients of DOFVector
        DenseVector<T> coeff{(size_t)localBasisFct->getNumber()};
        vecDV->getLocalVector(elInfo->getElement(), coeff);

        localBasisFct->evalGrdUh(lambda, grdLambda, coeff, result);
      }

      return result;
    }

    std::string str() const
    {
      return std::string("grad(") + vecDV->getName() + ")";
    }

  private:
    DOFVectorBase<T> const* vecDV;

    DenseVector<value_type> values;
  };

} // end namespace AMDiS
