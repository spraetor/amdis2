/** \file Operator.h */

#pragma once

#include <vector>

#include "AMDiS_fwd.h"
// #include "FixVec.h"
#include "Flag.h"
// #include "MatrixVector.h"
// #include "ElInfo.h"
#include "SubAssembler.h"
#include "OperatorTerm.h"
#include "ZeroOrderTerm.h"
#include "FirstOrderTerm.h"
#include "SecondOrderTerm.h"

#include <traits/category.hpp>
#include <expressions/expr_traits.hpp>

namespace AMDiS 
{
  /** \brief
   * An Operator holds all information needed to assemble the system matrix
   * and the right hand side. It consists of four OperatorTerm lists each storing
   * Terms of a specific order and type. You can define your own Operator by 
   * creating an empty Operator and than adding OperatorTerms to it.
   * By calling \ref getElementMatrix() or \ref getElementVector() one can 
   * initiate the assembling procedure. Therefor each Operator has its own
   * Assembler, especially created for this Operator, by the first call of
   * \ref getElementMatrix() or \ref getElementVector(). 
   */
  class Operator
  {
  public:
    /// Constructs an empty Operator of type operatorType for the given FiniteElemSpace.
    Operator(const FiniteElemSpace *rowFeSpace,
             const FiniteElemSpace *colFeSpace = NULL);
    
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
    void addZeroOrderTerm(C const& c)
    {
      using Expr = traits::to_expr<C>;
      addZOTImpl(Expr::get(c));
    }

    /// Adds a FirstOrderTerm to the Operator: < 1 * b * u, v >
    template <class B>
    void addFirstOrderTerm(B const& b, FirstOrderType type = GRD_PHI, int i = -1)
    {
      typedef typename traits::category<B>::value_type    ValueType;
      typedef typename traits::category<ValueType>::tag   Tag;
      using Expr = traits::to_expr<B>;
      addFOTImpl(Tag(), Expr::get(b), type, i);
    }

    /// Adds a SecondOrderTerm to the Operator
    template <class A, bool symmetric = false>
    void addSecondOrderTerm(A const& a, int i = -1, int j = -1, 
                            bool_<symmetric> s = bool_<symmetric>())
    {
      typedef typename traits::category<A>::value_type    ValueType;
      typedef typename traits::category<ValueType>::tag   Tag;
      typedef if_then_else< std::is_same<Tag, tag::scalar>::value, 
          bool_<true>, bool_<symmetric> > Sym;
      
      using Expr = traits::to_expr<A>;
      addSOTImpl(Tag(), Expr::get(a), i, j, Sym::value);
    }

    /// Calculates the element matrix for this ElInfo and adds it multiplied by
    /// factor to userMat.
    virtual void getElementMatrix(const ElInfo *elInfo, 
                                  ElementMatrix& userMat, 
                                  double factor = 1.0);

    /// Calculates the element vector for this ElInfo and adds it multiplied by
    /// factor to userVec.
    virtual void getElementVector(const ElInfo *elInfo, 
                                  DenseVector<double>& userVec, 
                                  double factor = 1.0);

    /// That function must be called after one assembling cycle has been finished.
    void finishAssembling();

    /// Returns \ref rowFeSpace.
    const FiniteElemSpace *getRowFeSpace() const
    { 
      return rowFeSpace; 
    }

    /// Returns \ref colFeSpace.
    const FiniteElemSpace *getColFeSpace() const
    { 
      return colFeSpace; 
    }

    /// Returns \ref auxFeSpaces.
    std::set<const FiniteElemSpace*>& getAuxFeSpaces()
    {
      return auxFeSpaces;
    }

    /// Sets \ref uhOld.
    void setUhOld(const DOFVectorBase<double> *uhOld);

    /// Returns \ref uhOld.
    const DOFVectorBase<double> *getUhOld() const
    {
      return uhOld;
    }    

    /// Returns \ref assembler
    Assembler *getAssembler();

    /// Sets \ref assembler
    void setAssembler(Assembler *ass)
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
    void initAssembler(Quadrature *quad2,
                       Quadrature *quad1GrdPsi,
                       Quadrature *quad1GrdPhi,
                       Quadrature *quad0);


    /// Calculates the needed quadrature degree for the given order. 
    int getQuadratureDegree(int order, FirstOrderType firstOrderType = GRD_PHI);

    /// Evaluation of all terms in \ref zeroOrder. 
    void evalZeroOrder(int nPoints,		       
                       const DenseVector<double>& uhAtQP,
                       const DenseVector<WorldVector<double> >& grdUhAtQP,
                       const DenseVector<WorldMatrix<double> >& D2UhAtQP,
                       DenseVector<double>& result,
                       double factor) const
    {
      for (OperatorTerm const* term : zeroOrder)
        term->eval(nPoints, uhAtQP, grdUhAtQP, D2UhAtQP, result, factor);
    }

    /// Evaluation of all terms in \ref firstOrderGrdPsi. 
    void evalFirstOrderGrdPsi(int nPoints,
                              const DenseVector<double>& uhAtQP,
                              const DenseVector<WorldVector<double> >& grdUhAtQP,
                              const DenseVector<WorldMatrix<double> >& D2UhAtQP,
                              DenseVector<double>& result,
                              double factor) const
    {      
      for (OperatorTerm const* term : firstOrderGrdPsi)
        term->eval(nPoints, uhAtQP, grdUhAtQP, D2UhAtQP, result, factor);
    }

    /// Evaluation of all terms in \ref firstOrderGrdPhi. 
    void evalFirstOrderGrdPhi(int nPoints,
                              const DenseVector<double>& uhAtQP,
                              const DenseVector<WorldVector<double> >& grdUhAtQP,
                              const DenseVector<WorldMatrix<double> >& D2UhAtQP,
                              DenseVector<double>& result,
                              double factor) const
    {
      for (OperatorTerm const* term : firstOrderGrdPhi)
        term->eval(nPoints, uhAtQP, grdUhAtQP, D2UhAtQP, result, factor);
    }

    /// Evaluation of all terms in \ref secondOrder. 
    void evalSecondOrder(int nPoints,
                         const DenseVector<double>& uhAtQP,
                         const DenseVector<WorldVector<double> >& grdUhAtQP,
                         const DenseVector<WorldMatrix<double> >& D2UhAtQP,
                         DenseVector<double>& result,
                         double factor) const
    {      
      for (OperatorTerm const* term : secondOrder)
        term->eval(nPoints, uhAtQP, grdUhAtQP, D2UhAtQP, result, factor);
    }

    /// Weak evaluation of all terms in \ref secondOrder. 
    void weakEvalSecondOrder(const std::vector<WorldVector<double> > &grdUhAtQP,
                             std::vector<WorldVector<double> > &result) const;
  
    /// Calls getLALt() for each term in \ref secondOrder and adds the results to LALt.
    void getLALt(const ElInfo *elInfo, std::vector<mtl::dense2D<double> > &LALt) const
    {      
      for (OperatorTerm const* term : secondOrder)
        static_cast<SecondOrderTerm const*>(term)->getLALt(elInfo, LALt);
    }
  
    /// Calls getLb() for each term in \ref firstOrderGrdPsi and adds the 
    /// results to Lb.
    void getLbGrdPsi(const ElInfo *elInfo, 
                     std::vector<DenseVector<double> >& Lb) const
    {      
      for (OperatorTerm const* term : firstOrderGrdPsi)
        static_cast<FirstOrderTerm const*>(term)->getLb(elInfo, Lb);
    }

    /// Calls getLb() for each term in \ref firstOrderGrdPhi and adds the 
    /// results to Lb.
    void getLbGrdPhi(const ElInfo *elInfo, 
                     std::vector<DenseVector<double> >& Lb) const
    {      
      for (OperatorTerm const* term : firstOrderGrdPhi)
        static_cast<FirstOrderTerm const*>(term)->getLb(elInfo, Lb);
    }

    /// Calls getC() for each term in \ref zeroOrder and adds the results to c.
    void getC(const ElInfo *elInfo, int nPoints, DenseVector<double>& c)
    {      
      for (OperatorTerm const* term : zeroOrder)
        static_cast<ZeroOrderTerm const*>(term)->getC(elInfo, nPoints, c);
    }

    /// Returns true, if there are second order terms. Returns false otherwise.
    bool secondOrderTerms() const
    {
      return !secondOrder.empty();
    }

    /// Returns true, if there are first order terms (grdPsi). Returns 
    /// false otherwise.
    bool firstOrderTermsGrdPsi() const
    {
      return !firstOrderGrdPsi.empty();
    }

    /// Returns true, if there are first order terms (grdPhi). Returns 
    /// false otherwise.
    bool firstOrderTermsGrdPhi() const
    {
      return !firstOrderGrdPhi.empty();
    }

    /// Returns true, if there are zero order terms. Returns false otherwise.
    bool zeroOrderTerms() const
    {
      return !zeroOrder.empty();
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
    const FiniteElemSpace *rowFeSpace;

    /// FiniteElemSpace for matrix columns. Can be equal to \rowFeSpace.
    const FiniteElemSpace *colFeSpace;

    /// List of aux fe spaces, e.g., if a term is multiplied with a DOF vector
    std::set<const FiniteElemSpace*> auxFeSpaces;

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

    friend class Assembler;
    friend class SubAssembler;
    friend class ZeroOrderAssembler;
    friend class FirstOrderAssembler;
    friend class SecondOrderAssembler;
  };

} // end namespace AMDiS

#include "Operator.hh"
