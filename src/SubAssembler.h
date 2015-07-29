/** \file SubAssembler.h */

#pragma once

#include <map>
#include <vector>

#include <boost/any.hpp>

#include "AMDiS_fwd.h"
#include "FixVec.h"
#include "Flag.h"

namespace AMDiS 
{
  /**
   * \ingroup Assembler
   * 
   * \brief
   * Base class for SecondOrderAssembler, FirstOrderAssembler, 
   * ZeroOrderAssembler. The task of a SubAssembler is to assemble a list of 
   * terms of a special order and add their contributions to a DOFMatrix or a 
   * DOFVector. An Assembler can consist of up to four SubAssemblers: one 
   * SecondOrderAssembler for second order terms, one ZeroOrderAssembler for 
   * terms of order zero, and two FirstOrderAssemblers. One for terms with
   * derivatives of the basis functions connected to to row DOFs and one for 
   * those connected to the column DOFs.
   */
  class SubAssembler
  {
  public:
    /// Creates a SubAssembler belonging to assembler for the terms of given 
    /// order of Operator op. If order is equal to one, type spezifies what kind 
    /// of FirstOrderType are to assemble. During construction of a SubAssembler
    /// the needs and properties of the terms are considered.
    SubAssembler(Operator *op,
            		 Assembler *assembler,
            		 Quadrature *quadrat,
            		 int order, 
            		 bool optimized,
            		 FirstOrderType type = GRD_PHI);

    /// Calculates the element matrix for elInfo and adds it to mat. Memory for
    /// mat must be provided by the caller.
    void calculateElementMatrix(const ElInfo *elInfo,
				                        ElementMatrix& mat)
    {
      calculateElementMatrixImpl(elInfo, mat);
    }

    /// Calculates the element vector for elInfo and adds it to vec. Memory for
    /// vec must be provided by the caller.
    void calculateElementVector(const ElInfo *elInfo,
			                          ElementVector& vec)
    {
      calculateElementVectorImpl(elInfo, vec);
    }

    /// Called once for each ElInfo when \ref calculateElementMatrix() or
    /// \ref calculateElementVector() is called for the first time for this
    /// Element.
    void initElement(const ElInfo *smallElInfo,
            		     const ElInfo *largeElInfo = NULL,
            		     Quadrature *quad = NULL)
    {
      initImpl(smallElInfo, largeElInfo, quad);
    }
    
    /// Returns \ref terms
    std::vector<OperatorTerm*>* getTerms()
    { 
      return &terms; 
    }

    /// Returns \ref quadrature.
    Quadrature* getQuadrature()
    {
      return quadrature;
    }

    /// Returns \ref psiFast.
    const FastQuadrature* getPsiFast() const 
    { 
      return psiFast; 
    }

    // Returns \ref phiFast.
    const FastQuadrature* getPhiFast() const 
    { 
      return phiFast; 
    }

    /// Returns \ref name.
    std::string getName() const
    {
      return name;
    }

    /// Sets \ref quadrature to q.
    void setQuadrature(Quadrature* q) 
    {
      quadrature = q;
    }
  
    /// Creates a vector with the world coordinates of the quadrature points
    /// of \ref quadrature on the element of elInfo. 
    /// Used by \ref OperatorTerm::initElement().
    void getCoordsAtQPs(const ElInfo* elInfo, 
                  			Quadrature *quad,
                  			DenseVector<WorldVector<double> >& coordsAtQPs);

    /// DOFVector dv evaluated at quadrature points.
    /// Used by \ref OperatorTerm::initElement().
    template <class T>
    void getVectorAtQPs(DOFVectorBase<T>* dv, 
                  			const ElInfo* elInfo,
                  			Quadrature *quad,
                  			DenseVector<T>& vecAtQPs);
    
    /// Gradients of DOFVector dv evaluated at quadrature points.
    /// Used by \ref OperatorTerm::initElement().
    template <class T>
    void getGradientsAtQPs(DOFVectorBase<T>* dv,
                  			   const ElInfo* elInfo,
                  			   Quadrature *quad,
                  			   DenseVector<Gradient_t<T>>& grdAtQPs);

    /// The comp'th component of the derivative of DOFVector dv evaluated at
    /// quadrature points. Used by \ref OperatorTerm::initElement().
    /// Attention: not caching at the moment! Using cache if gradients for read 
    /// but not for write.
    template <class T>
    void getDerivativeAtQPs(DOFVectorBase<T>* dv,
                  			    const ElInfo* elInfo,
                  			    Quadrature *quad,
                  			    int comp,
                  			    DenseVector<T>& grdAtQPs);
    
  private:
    // must be implemented by derived class
    virtual void calculateElementMatrixImpl(const ElInfo *elInfo,
				                                    ElementMatrix& mat) = 0;

    // must be implemented by derived class.
    virtual void calculateElementVectorImpl(const ElInfo *elInfo,
					                                  ElementVector& vec) = 0;

    // calls initElement for all OperatorTerms
    virtual void initImpl(const ElInfo *smallElInfo, 
                  			  const ElInfo *largeElInfo,
                  			  Quadrature *quad);
    
    
  protected:
    /// Updates \ref psiFast and \ref phiFast.
    FastQuadrature *updateFastQuadrature(FastQuadrature *quadFast,
                              					 const BasisFunction *psi,
                              					 Flag updateFlag);
  
  protected:
    /// Problem dimension
    int dim;

    /// Row FiniteElemSpace.
    const FiniteElemSpace *rowFeSpace;

    /// Column FiniteElemSpace.
    const FiniteElemSpace *colFeSpace;

    /// Number of rows of the element matrix and length of the element
    /// vector. Is equal to the number of row basis functions
    int nRow;

    /// Number of columns of the element matrix. Is equal to the number
    /// of column basis functions
    int nCol;

    // TODO: try to remove boost::any
    /// Used for \ref getVectorAtQPs() and \ref getGradientsAtQPs().
    struct ValuesAtQPs {
      ValuesAtQPs()
        : valid(false), quad(NULL)
      {}
      ValuesAtQPs(boost::any values_, bool valid_, Quadrature* quad_=NULL)
        : values(values_), valid(valid_), quad(quad_)
      {}
      
      boost::any values; // used for DenseVector<T>
      bool valid;
      Quadrature* quad;
    };

    std::map<const DOFIndexedBase*, ValuesAtQPs* > cachedValuesAtQPs;
    std::map<const DOFIndexedBase*, ValuesAtQPs* > cachedGradientsAtQPs;
    
    /// Set and updated by \ref initElement() for each ElInfo. 
    /// coordsAtQPs[i] points to the coordinates of the i-th quadrature point.
    DenseVector<WorldVector<double> > cacheCoordsAtQPs;

    /// Used for \ref getCoordsAtQPs().
    bool coordsValid;

    /// Used for \ref getCoordsAtQP(). Stores the number of allocated 
    /// WorldVectors.
    int coordsNumAllocated;

    /// Quadrature object to be used for assembling.
    Quadrature *quadrature;

    /// FastQuadrature for row basis functions
    FastQuadrature *psiFast;

    /// FastQuadrature for column basis functions
    FastQuadrature *phiFast;

    /// Flag that specifies whether the element matrix is symmetric.
    bool symmetric;

    /// List of all terms with a contribution to this SubAssembler
    std::vector<OperatorTerm*> terms;

    ///
    bool opt;

    ///
    bool firstCall;

    /// Name of the assembler. Is used to print information about 
    /// used assembler.
    std::string name;

    friend class Assembler;
  };

} // end namespace AMDiS

#include "SubAssembler.hh"
