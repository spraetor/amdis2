#pragma once

#include <vector>

#include "AMDiS_fwd.hpp"
#include "FirstOrderTerm.hpp"
#include "Flag.hpp"
#include "OperatorTerm.hpp"
#include "SecondOrderTerm.hpp"
#include "SubAssembler.hpp"
#include "ZeroOrderTerm.hpp"
#include "expressions/TermGenerator.hpp"
#include "traits/category.hpp"
// #include "ElInfo.hpp"
// #include "FixVec.hpp"
// #include "MatrixVector.hpp"

namespace AMDiS
{
  /** \brief
   * An Operator holds all information needed to assemble the system matrix
   * and the right hand side. It consists of four OperatorTerm lists each storing
   * Terms of a specific order and type. You can define your own Operator by
   * creating an empty Operator and than adding OperatorTerms to it.
   * By calling \ref getElementMatrix() or \ref getElementVector() one can
   * initiate the assembling procedure. Therefore each Operator has its own
   * Assembler, especially created for this Operator, by the first call of
   * \ref getElementMatrix() or \ref getElementVector().
   */
  class Operator
  {
  public:
    /// Constructs an empty Operator of type operatorType for the given FiniteElemSpace.
    Operator(FiniteElemSpace const* rowFeSpace,
             FiniteElemSpace const* colFeSpace = NULL);

    Operator(Operator const&) = default;

    Operator(Operator&&) = default;

    /// Destructor.
    virtual ~Operator() {}

    /// Sets \ref optimized.
    void useOptimizedAssembler(bool opt)
    {
      optimized = opt;
    }

    /// Returns \ref optimized.
    bool isOptimized() const
    {
      return optimized;
    }

    /// Adds a ZeroOrderTerm to the Operator: < c * u, v >
    template <class C>
    Operator& addZeroOrderTerm(C const& c)
    {
      addZOTImpl(toTerm(c));
      return *this;
    }

    template <class C>
    Operator& addZOT(C const& c)
    {
      return addZeroOrderTerm(c);
    }


    /// Adds a FirstOrderTerm to the Operator: < 1 * b * u, v >
    template <class B>
    Operator& addFirstOrderTerm(B const& b, FirstOrderType type = GRD_PHI, int i = -1)
    {
      using ValueType = Value_t<traits::category<B>>;
      using Tag = typename traits::category<ValueType>::tag;
      addFOTImpl(Tag(), toTerm(b), type, i);
      return *this;
    }

    template <class B>
    Operator& addFOT(B const& b, FirstOrderType type = GRD_PHI, int i = -1)
    {
      return addFirstOrderTerm(b, type, i);
    }

    /// Adds a SecondOrderTerm to the Operator
    template <class A, bool symmetric = false>
    Operator& addSecondOrderTerm(A const& a, int i = -1, int j = -1,
                                 bool_<symmetric> /*s*/ = bool_<symmetric>())
    {
      using ValueType = Value_t<traits::category<A>>;
      using Tag = typename traits::category<ValueType>::tag;
      using Sym = if_then_else<std::is_same<Tag, tag::scalar>::value,
            bool_<true>, bool_<symmetric>>;
      addSOTImpl(Tag(), toTerm(a), i, j, Sym::value);
      return *this;
    }

    template <class A, bool symmetric = false>
    Operator& addSOT(A const& a, int i = -1, int j = -1,
		     bool_<symmetric> s = bool_<symmetric>())
    {
      return addSecondOrderTerm(a, i, j, s);
    }

    /// Calculates the element matrix for this ElInfo and adds it multiplied by
    /// factor to userMat.
    virtual void getElementMatrix(ElInfo const* elInfo,
                                  ElementMatrix& userMat,
                                  double factor = 1.0);

    /// Calculates the element vector for this ElInfo and adds it multiplied by
    /// factor to userVec.
    virtual void getElementVector(ElInfo const* elInfo,
                                  DenseVector<double>& userVec,
                                  double factor = 1.0);

    /// That function must be called after one assembling cycle has been finished.
    void finishAssembling();

    /// Returns \ref rowFeSpace.
    FiniteElemSpace const* getRowFeSpace() const
    {
      return rowFeSpace;
    }

    /// Returns \ref colFeSpace.
    FiniteElemSpace const* getColFeSpace() const
    {
      return colFeSpace;
    }

    /// Returns \ref auxFeSpaces.
    std::set<FiniteElemSpace const*>& getAuxFeSpaces()
    {
      return auxFeSpaces;
    }

    /// Sets \ref uhOld.
    void setUhOld(DOFVectorBase<double> const* uhOld);

    /// Returns \ref uhOld.
    DOFVectorBase<double> const* getUhOld() const
    {
      return uhOld;
    }

    /// Returns \ref assembler
    Assembler* getAssembler();

    /// Sets \ref assembler
    void setAssembler(Assembler* ass)
    {
      assembler.set(ass);
    }

    /// Sets \ref fillFlag, the flag used for this Operator while mesh traversal.
    void setFillFlag(Flag f)
    {
      fillFlag = f;
    }

    /// Sets \ref needDualTraverse.
    void setNeedDualTraverse(bool b)
    {
      needDualTraverse = b;
    }

    /// Returns \ref fillFlag
    Flag getFillFlag() const
    {
      return fillFlag;
    }

    /// Returns \ref needDualTraverse
    bool getNeedDualTraverse() const
    {
      return needDualTraverse;
    }

    /// Initialization of the needed SubAssemblers using the given quadratures.
    void initAssembler(Quadrature* quad2,
                       Quadrature* quad1GrdPsi,
                       Quadrature* quad1GrdPhi,
                       Quadrature* quad0);


    /// Calculates the needed quadrature degree for the given order.
    int getQuadratureDegree(int order, FirstOrderType firstOrderType = GRD_PHI);

    /// Evaluation of all terms in \ref zeroOrder.
    void evalZeroOrder(int nPoints,
                       DenseVector<double> const& uhAtQP,
                       DenseVector<WorldVector<double>> const& grdUhAtQP,
                       DenseVector<WorldMatrix<double>> const& D2UhAtQP,
                       DenseVector<double>& result,
                       double factor) const
    {
      for (OperatorTerm const* term : zeroOrder)
        term->eval(nPoints, uhAtQP, grdUhAtQP, D2UhAtQP, result, factor);
    }

    /// Evaluation of all terms in \ref firstOrderGrdPsi.
    void evalFirstOrderGrdPsi(int nPoints,
                              DenseVector<double> const& uhAtQP,
                              DenseVector<WorldVector<double>> const& grdUhAtQP,
                              DenseVector<WorldMatrix<double>> const& D2UhAtQP,
                              DenseVector<double>& result,
                              double factor) const
    {
      for (OperatorTerm const* term : firstOrderGrdPsi)
        term->eval(nPoints, uhAtQP, grdUhAtQP, D2UhAtQP, result, factor);
    }

    /// Evaluation of all terms in \ref firstOrderGrdPhi.
    void evalFirstOrderGrdPhi(int nPoints,
                              DenseVector<double> const& uhAtQP,
                              DenseVector<WorldVector<double>> const& grdUhAtQP,
                              DenseVector<WorldMatrix<double>> const& D2UhAtQP,
                              DenseVector<double>& result,
                              double factor) const
    {
      for (OperatorTerm const* term : firstOrderGrdPhi)
        term->eval(nPoints, uhAtQP, grdUhAtQP, D2UhAtQP, result, factor);
    }

    /// Evaluation of all terms in \ref secondOrder.
    void evalSecondOrder(int nPoints,
                         DenseVector<double> const& uhAtQP,
                         DenseVector<WorldVector<double>> const& grdUhAtQP,
                         DenseVector<WorldMatrix<double>> const& D2UhAtQP,
                         DenseVector<double>& result,
                         double factor) const
    {
      for (OperatorTerm const* term : secondOrder)
        term->eval(nPoints, uhAtQP, grdUhAtQP, D2UhAtQP, result, factor);
    }

    /// Weak evaluation of all terms in \ref secondOrder.
    void weakEvalSecondOrder(std::vector<WorldVector<double>> const& grdUhAtQP,
                             std::vector<WorldVector<double>>& result) const;

    /// Calls getLALt() for each term in \ref secondOrder and adds the results to LALt.
    void getLALt(ElInfo const* elInfo, std::vector<mtl::dense2D<double>>& LALt) const
    {
      for (OperatorTerm const* term : secondOrder)
        static_cast<SecondOrderTerm const*>(term)->getLALt(elInfo, LALt);
    }

    /// Calls getLb() for each term in \ref firstOrderGrdPsi and adds the
    /// results to Lb.
    void getLbGrdPsi(ElInfo const* elInfo,
                     std::vector<DenseVector<double>>& Lb) const
    {
      for (OperatorTerm const* term : firstOrderGrdPsi)
        static_cast<FirstOrderTerm const*>(term)->getLb(elInfo, Lb);
    }

    /// Calls getLb() for each term in \ref firstOrderGrdPhi and adds the
    /// results to Lb.
    void getLbGrdPhi(ElInfo const* elInfo,
                     std::vector<DenseVector<double>>& Lb) const
    {
      for (OperatorTerm const* term : firstOrderGrdPhi)
        static_cast<FirstOrderTerm const*>(term)->getLb(elInfo, Lb);
    }

    /// Calls getC() for each term in \ref zeroOrder and adds the results to c.
    void getC(ElInfo const* elInfo, int nPoints, DenseVector<double>& c)
    {
      for (OperatorTerm const* term : zeroOrder)
        static_cast<ZeroOrderTerm const*>(term)->getC(elInfo, nPoints, c);
    }

    /// Returns true, if there are second order terms. Returns false otherwise.
    bool secondOrderTerms() const
    {
      return !secondOrder.empty();
    }

    std::vector<OperatorTerm*> const& getSecondOrder() const
    {
      return secondOrder;
    }

    /// Returns true, if there are first order terms (grdPsi). Returns
    /// false otherwise.
    bool firstOrderTermsGrdPsi() const
    {
      return !firstOrderGrdPsi.empty();
    }

    std::vector<OperatorTerm*> const& getFirstOrderGrdPsi() const
    {
      return firstOrderGrdPsi;
    }

    /// Returns true, if there are first order terms (grdPhi). Returns
    /// false otherwise.
    bool firstOrderTermsGrdPhi() const
    {
      return !firstOrderGrdPhi.empty();
    }

    std::vector<OperatorTerm*> const& getFirstOrderGrdPhi() const
    {
      return firstOrderGrdPhi;
    }

    /// Returns true, if there are zero order terms. Returns false otherwise.
    bool zeroOrderTerms() const
    {
      return !zeroOrder.empty();
    }

    std::vector<OperatorTerm*> const& getZeroOrder() const
    {
      return zeroOrder;
    }

  protected:
    template <class C>
    void addZOTImpl(C const& c);

    template <class B>
    void addFOTImpl(tag::scalar, B const& b, FirstOrderType type, int i);

    template <class B>
    void addFOTImpl(tag::vector, B const& b, FirstOrderType type, int i);

    template <class A>
    void addSOTImpl(tag::scalar, A const& a, int i, int j, bool sym);

    template <class A>
    void addSOTImpl(tag::matrix, A const& a, int i, int j, bool sym);

  protected:
    /// FiniteElemSpace for matrix rows and element vector
    FiniteElemSpace const* rowFeSpace;

    /// FiniteElemSpace for matrix columns. Can be equal to \rowFeSpace.
    FiniteElemSpace const* colFeSpace;

    /// List of aux fe spaces, e.g., if a term is multiplied with a DOF vector
    std::set<FiniteElemSpace const*> auxFeSpaces;

    /// Number of rows in the element matrix
    int nRow;

    /// Number of columns in the element matrix
    int nCol;

    /// Flag for mesh traversal
    Flag fillFlag;

    /// If true, the operator needs a dual traverse over two meshes for assembling.
    bool needDualTraverse;

    /// Calculates the element matrix and/or the element vector. It is
    /// created especially for this Operator, when \ref getElementMatrix()
    /// or \ref getElementVector is called for the first time.
    ThreadPrivate<Assembler*> assembler;

    /// List of all second order terms
    std::vector<OperatorTerm*> secondOrder;

    /// List of all first order terms derived to psi
    std::vector<OperatorTerm*> firstOrderGrdPsi;

    /// List of all first order terms derived to phi
    std::vector<OperatorTerm*> firstOrderGrdPhi;

    /// List of all zero order terms
    std::vector<OperatorTerm*> zeroOrder;

    /// Pointer to the solution of the last timestep. Can be used for a more
    /// efficient assemblage of the element vector when the element matrix
    /// was already computed.
    DOFVectorBase<double> const* uhOld;

    /// Spezifies whether optimized assemblers are used or not.
    bool optimized;
  };

} // end namespace AMDiS

#include "Operator.hh"
