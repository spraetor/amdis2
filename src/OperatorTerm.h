/** \file OperatorTerm.h */

#pragma once

#include <set>

#include "AMDiS_fwd.h"
#include "SubAssembler.h"
#include <traits/traits.hpp>

namespace AMDiS
{
  /**
   * \ingroup Assembler
   *
   * \brief
   * Base class for ZeroOrderTerm, FirstOrderTerm and SecondOrderTerm.
   * OperatorTerms are the building blocks of an Operator. Each OperatorTerm
   * has its properties which are regarded, when constructing
   * an Assembler for the corresponding Operator.
   */
  class OperatorTerm
  {
  public:
    /// Constructs an OperatorTerm with initially no properties.
    /// degree_ is used to determine the degree of the needed quadrature
    /// for the assemblage.
    OperatorTerm(int deg)
      : degree(deg),
        dimOfWorld(Global::getGeo(WORLD))
    {}

    /// Destructor.
    virtual ~OperatorTerm() {}

    /// Virtual method. It's called by SubAssembler::initElement() for
    /// each OperatorTerm belonging to this SubAssembler. E.g., vectors
    /// and coordinates at quadrature points can be calculated here.
    void initElement(ElInfo const* elInfo, 
		     SubAssembler* subAssembler,
                     Quadrature* quad = NULL)
    {
      initImpl(elInfo, subAssembler, quad);
    }

    /// Returs \auxFeSpaces, the list of all aux fe spaces the operator makes
    /// use off.
    std::set<FiniteElemSpace const*>& getAuxFeSpaces()
    {
      return auxFeSpaces;
    }

    /// Specifies whether the matrix of the term is symmetric
    void setSymmetric(bool symm);

    /// Returns true, if the term is piecewise constant, returns false otherwise.
    bool isPWConst() const
    {
      return (degree == 0);
    }

    /// Returns true, if the term has a symmetric matrix, returns false otherwise.
    bool isSymmetric();

    /// Returns \ref degree.
    int getDegree() const
    {
      return degree;
    }

    /// Sets one component of the b vector to be one. See \ref bOne.
    void setB(int b)
    {
      bOne = b;
    }

    void setOperator(Operator* op)
    {
      operat = op;
    }

    /// Evaluation of the OperatorTerm at all quadrature points.
    void eval(int nPoints,
              DenseVector<double> const& uhAtQP,
              DenseVector<WorldVector<double>> const& grdUhAtQP,
              DenseVector<WorldMatrix<double>> const& D2UhAtQP,
              DenseVector<double>& result,
              double factor) const
    {
      evalImpl(nPoints, uhAtQP, grdUhAtQP, D2UhAtQP, result, factor);
    }

  private:
    // default behavior: init nothing
    virtual void initImpl(ElInfo const*, SubAssembler*, Quadrature*) { }

    // must be implemented by derived class
    virtual void evalImpl(int nPoints,
                          DenseVector<double> const& uhAtQP,
                          DenseVector<WorldVector<double>> const& grdUhAtQP,
                          DenseVector<WorldMatrix<double>> const& D2UhAtQP,
                          DenseVector<double>& result,
                          double factor) const = 0;

  protected:
    /// Stores the properties of this OperatorTerm
    Flag properties = 0;

    /// Polynomial degree of the term. Used to detemine the degree of the quadrature.
    int degree;

    /// Stores the dimension of the world.
    int dimOfWorld;

    /// List off all fe spaces, the operator term makes use off.
    std::set<FiniteElemSpace const*> auxFeSpaces;

    /// Pointer to the Operator this OperatorTerm belongs to.
    Operator* operat;

    /// In many cases, the vector b in the evaluation \f$ \Lambda \cdot b\f$ has
    /// zeros in all components expect one that is set to one. Using the function
    /// \ref lb is then unnecessary time consuming. Instead, this variable
    /// defines the component of the vector b to be one. The function \ref lb_one
    /// is used if this variable is not -1.
    int bOne = -1;

    /// Flag for piecewise constant terms
    static const Flag PW_CONST;

    /// Flag for symmetric terms
    static const Flag SYMMETRIC;
  };


  // forward declarations
  class ZeroOrderTerm;
  class FirstOrderTerm;
  class SecondOrderTerm;

  /// helper class to adopt the correct OperatorTerm based on the term order
  template <int Order>
  struct GetTerm
  {
    using type =  if_then_else<Order == 0, ZeroOrderTerm,
		  if_then_else<Order == 1, FirstOrderTerm,
		  if_then_else<Order == 2, SecondOrderTerm,
					   OperatorTerm >>>;
  };


  /// basic interface for OperatorTerms based on expressions
  template <class Term, int Order = -1>
  class GenericOperatorTerm : public GetTerm<Order>::type
  {
    using Super = typename GetTerm<Order>::type;

  public:
    /// Expression term stored as copy
    Term term;

    /// constructor
    /// adds all feSpaces provided by the expression term to auxFeSpaces liste
    template <class Term_, 
	      class = Requires_t< concepts::Compatible<Term, Term_> >>
    GenericOperatorTerm(Term_&& term_)
      : Super(term_.getDegree()),
        term{term_}
    {
      term.insertFeSpaces(this->auxFeSpaces);
#ifndef NDEBUG
      testAuxFeSpaces(this->auxFeSpaces);
#endif
    }

  private:
    // Implements OperatorTerm::initImpl().
    // calls initElement() on term
    virtual void initImpl(ElInfo const* elInfo,
                          SubAssembler* subAssembler,
                          Quadrature* quad) override
    {
      term.initElement(elInfo, subAssembler, quad, NULL);
    }

    /// test for only one mesh allowed in expressions
    template <class FeSpaceList>
    void testAuxFeSpaces(FeSpaceList const& auxFeSpaces)
    {
      FUNCNAME("GenericOperatorTerm::testAuxFeSpaces()");
      if (auxFeSpaces.size() > 0)
      {
        Mesh* mesh0 = (*begin(auxFeSpaces))->getMesh();
        for (auto const* feSpace : auxFeSpaces)
        {
          if (feSpace->getMesh() != mesh0)
          {
            ERROR_EXIT("Only one mesh allowed in expression.\n");
          }
        }
      }
    }
  };


  /// \brief Implementation of GenericOperatorTerm for the default Order parameter.
  /// This class is instantiated if no Order parameter is given and inherits the
  /// default implementation from the primary template (e.g. Order < -1)
  template <class Term>
  class GenericOperatorTerm<Term, -1> : public GenericOperatorTerm<Term, -2>
  {
    using Super = GenericOperatorTerm<Term, -2>;
    
  public:
    template <class Term_, 
	      class = Requires_t< concepts::Compatible<Term, Term_> >>
    GenericOperatorTerm(Term_&& term_)
      : Super(std::forward<Term_>(term_)) {}

  private:
    // Implements \ref OperatorTerm::eval().
    virtual void evalImpl(int nPoints,
                          DenseVector<double> const& uhAtQP,
                          DenseVector<WorldVector<double>> const& grdUhAtQP,
                          DenseVector<WorldMatrix<double>> const& D2UhAtQP,
                          DenseVector<double>& result,
                          double factor) const override {};
  };

} // end namespace AMDiS
