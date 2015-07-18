/******************************************************************************
 *
 * AMDiS - Adaptive multidimensional simulations
 *
 * Copyright (C) 2013 Dresden University of Technology. All Rights Reserved.
 * Web: https://fusionforge.zih.tu-dresden.de/projects/amdis
 *
 * Authors: 
 * Simon Vey, Thomas Witkowski, Andreas Naumann, Simon Praetorius, et al.
 *
 * This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
 * WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 *
 *
 * This file is part of AMDiS
 *
 * See also license.opensource.txt in the distribution.
 * 
 ******************************************************************************/



/** \file Lagrange.h */

#ifndef AMDIS_LAGRANGE_H
#define AMDIS_LAGRANGE_H

#include <list>
#include <boost/numeric/mtl/mtl.hpp>
#include "BasisFunction.h"
#include "FixVec.h"

namespace AMDiS {

#define MAX_DIM 3
#define MAX_DEGREE 4

  /** \ingroup FEMSpace
   * \brief
   * Lagrange basis functions. Sub class of BasisFunction
   */
  class Lagrange : public BasisFunction
  {
  public:
    /// Creator class used in the BasisFunctionCreatorMap.
    class Creator : public BasisFunctionCreator
    {
    public:
      Creator(int degree_) : degree(degree_) {}
      virtual ~Creator() {}
      
      /// Returns a new Lagrange object.
      BasisFunction* create() 
      { 
	return getLagrange(this->dim, degree); 
      }
      
    protected:
      int degree;
    };
    
  protected:
    /// Constructs lagrange basis functions with the given dim and degree.
    /// Constructor is protected to avoid multiple instantiation of identical
    /// basis functions. Use \ref getLagrange instead.
    Lagrange(int dim_, int degree_);

    /** \brief
     * destructor
     */
    virtual ~Lagrange();

  public:
    /// Returns a pointer to lagrange basis functions with the given dim and
    /// degree. Multiple instantiation of identical basis functions is avoided
    /// by rembering once created basis functions in \ref allBasFcts.
    static Lagrange* getLagrange(int dim, int degree);

    /// Implements BasisFunction::interpol
    void interpol(const ElInfo *, int, const int *, 
		  std::function<double(WorldVector<double>)>, 
		  mtl::dense_vector<double>&) const;

    /// Implements BasisFunction::interpol
    void interpol(const ElInfo *, int, 
		  const int *b_no,
		  std::function<WorldVector<double>(WorldVector<double>)>, 
		  mtl::dense_vector<WorldVector<double> >&) const;

    /// Returns the barycentric coordinates of the i-th basis function.
    DimVec<double> *getCoords(int i) const;

    /// Implements BasisFunction::getBound
    void getBound(const ElInfo*, BoundaryType *) const;

    /** \brief
     * Calculates the local vertex indices which are involved in evaluating
     * the nodeIndex-th DOF at the positionIndex-th part of type position 
     * (VERTEX/EDGE/FACE/CENTER). nodeIndex determines the permutation
     * of the involved vertices. So in 1d for lagrange4 there are two DOFs at
     * the CENTER (which is an edge in this case). Then vertices[0] = {0, 1} and
     * vertices[1] = {1, 0}. This allows to use the same local basis function 
     * for all DOFs at the same position.
     */
    static void setVertices(int dim, int degree, 
			    GeoIndex position, int positionIndex, int nodeIndex, 
			    int** vertices);

    /// Implements BasisFunction::refineInter
    void refineInter(DOFIndexed<double> *drv, RCNeighbourList* list, int n) 
    {
      if (refineInter_fct)
	(*refineInter_fct)(drv, list, n, this);
    }

    /// Implements BasisFunction::coarseRestrict
    void coarseRestr(DOFIndexed<double> *drv, RCNeighbourList* list, int n)
    {
      if (coarseRestr_fct)
	(*coarseRestr_fct)(drv, list, n, this);
    }
  
    /// Implements BasisFunction::coarseInter
    void coarseInter(DOFIndexed<double> *drv, RCNeighbourList* list, int n) 
    {
      if (coarseInter_fct)
	(*coarseInter_fct)(drv, list, n, this);
    }
  
    /// Implements BasisFunction::getLocalIndices().
    void getLocalIndices(const Element *el,
			 const DOFAdmin *admin,
			 std::vector<DegreeOfFreedom> &dofs) const;

    void getLocalDofPtrVec(const Element *el, 
			   const DOFAdmin *admin,
			   std::vector<const DegreeOfFreedom*>& vec) const;

    /// Implements BasisFunction::l2ScpFctBas
    void l2ScpFctBas(Quadrature* q,
		     std::function<double(WorldVector<double>)> f,
		     DOFVector<double>* fh);

    /// Implements BasisFunction::l2ScpFctBas
    void l2ScpFctBas(Quadrature* q,
		     std::function<WorldVector<double>(WorldVector<double>)> f,
		     DOFVector<WorldVector<double> >* fh);

    static void clear();
    
    /// Implements BasisFunction::isnodal
    bool isNodal() const
    {
      return true;
    }

  protected:
    /// sets the barycentric coordinates (stored in \ref bary) of the local 
    /// basis functions.
    void setBary();

    /// Recursive calculation of coordinates. Used by \ref setBary
    void createCoords(int* coordInd, int numCoords, int dimIndex, int rest, 
		      DimVec<double>* vec = NULL);

    /// Used by \ref setBary
    int** getIndexPermutations(int numIndices) const;

    /// Implements BasisFunction::setNDOF
    void setNDOF();

    /// Sets used function pointers
    void setFunctionPointer();

    /// Used by \ref getVec
    int* orderOfPositionIndices(const Element* el, 
				GeoIndex position, 
				int positionIndex) const;

    /// Calculates the number of DOFs needed for Lagrange of the given dim 
    /// and degree.
    static int getNumberOfDofs(int dim, int degree);

  private:
    /// barycentric coordinates of the locations of all basis functions
    std::vector<DimVec<double>* > *bary;

    /** \name static dim-degree-arrays
     * \{
     */
    static std::vector<DimVec<double>* > baryDimDegree[MAX_DIM + 1][MAX_DEGREE + 1];
    static DimVec<int>* ndofDimDegree[MAX_DIM + 1][MAX_DEGREE + 1];
    static int nBasFctsDimDegree[MAX_DIM + 1][MAX_DEGREE + 1];
    static std::vector<BasFctType*> phiDimDegree[MAX_DIM + 1][MAX_DEGREE + 1];
    static std::vector<GrdBasFctType*> grdPhiDimDegree[MAX_DIM + 1][MAX_DEGREE + 1];
    static std::vector<D2BasFctType*> D2PhiDimDegree[MAX_DIM + 1][MAX_DEGREE + 1];
    /** \} */

    /// List of all used BasisFunctions in the whole program. Avoids duplicate
    /// instantiation of identical BasisFunctions. 
    static std::list<Lagrange*> allBasFcts;


  protected:
    /// Pointer to the used refineInter function
    void (*refineInter_fct)(DOFIndexed<double> *, RCNeighbourList*, int, BasisFunction*);

    /** \name refineInter functions
     * \{
     */  
    static void  refineInter0(DOFIndexed<double> *, RCNeighbourList*, int, 
			      BasisFunction*);
    static void  refineInter1(DOFIndexed<double> *, RCNeighbourList*, int, 
			      BasisFunction*);
    static void  refineInter2_1d(DOFIndexed<double> *, RCNeighbourList*, int, 
				 BasisFunction*);
    static void  refineInter2_2d(DOFIndexed<double> *, RCNeighbourList*, int, 
				 BasisFunction*);
    static void  refineInter2_3d(DOFIndexed<double> *, RCNeighbourList*, int, 
				 BasisFunction*);
    static void  refineInter3_1d(DOFIndexed<double> *, RCNeighbourList*, int, 
				 BasisFunction*);
    static void  refineInter3_2d(DOFIndexed<double> *, RCNeighbourList*, int, 
				 BasisFunction*);
    static void  refineInter3_3d(DOFIndexed<double> *, RCNeighbourList*, int, 
				 BasisFunction*);
    static void  refineInter4_1d(DOFIndexed<double> *, RCNeighbourList*, int, 
				 BasisFunction*);
    static void  refineInter4_2d(DOFIndexed<double> *, RCNeighbourList*, int, 
				 BasisFunction*);
    static void  refineInter4_3d(DOFIndexed<double> *, RCNeighbourList*, int, 
				 BasisFunction*);
    /** \} */

    /// Pointer to the used coarseRestr function
    void (*coarseRestr_fct)(DOFIndexed<double> *, RCNeighbourList*, int, BasisFunction*);

    /** \name coarseRestr functions
     * \{
     */  
    static void  coarseRestr0(DOFIndexed<double> *, RCNeighbourList*, int, 
			      BasisFunction*);
    static void  coarseRestr1(DOFIndexed<double> *, RCNeighbourList*, int, 
			      BasisFunction*);
    static void  coarseRestr2_1d(DOFIndexed<double> *, RCNeighbourList*, int, 
				 BasisFunction*);
    static void  coarseRestr2_2d(DOFIndexed<double> *, RCNeighbourList*, int, 
				 BasisFunction*);
    static void  coarseRestr2_3d(DOFIndexed<double> *, RCNeighbourList*, int, 
				 BasisFunction*);
    static void  coarseRestr3_1d(DOFIndexed<double> *, RCNeighbourList*, int, 
				 BasisFunction*);
    static void  coarseRestr3_2d(DOFIndexed<double> *, RCNeighbourList*, int, 
				 BasisFunction*);
    static void  coarseRestr3_3d(DOFIndexed<double> *, RCNeighbourList*, int, 
				 BasisFunction*);
    static void  coarseRestr4_1d(DOFIndexed<double> *, RCNeighbourList*, int, 
				 BasisFunction*);
    static void  coarseRestr4_2d(DOFIndexed<double> *, RCNeighbourList*, int, 
				 BasisFunction*);
    static void  coarseRestr4_3d(DOFIndexed<double> *, RCNeighbourList*, int, 
				 BasisFunction*);
    /** \} */

    /// Pointer to the used coarseInter function
    void (*coarseInter_fct)(DOFIndexed<double> *, RCNeighbourList*, int, BasisFunction*);  

    /** \name coarseInter functions
     * \{
     */
    static void  coarseInter0(DOFIndexed<double> *, RCNeighbourList*, int, 
			      BasisFunction*);
    static void  coarseInter1(DOFIndexed<double> *, RCNeighbourList*, int, 
			      BasisFunction*);
    static void  coarseInter2_1d(DOFIndexed<double> *, RCNeighbourList*, int, 
				 BasisFunction*);
    static void  coarseInter2_2d(DOFIndexed<double> *, RCNeighbourList*, int, 
				 BasisFunction*);
    static void  coarseInter2_3d(DOFIndexed<double> *, RCNeighbourList*, int, 
				 BasisFunction*);
    static void  coarseInter3_1d(DOFIndexed<double> *, RCNeighbourList*, int, 
				 BasisFunction*);
    static void  coarseInter3_2d(DOFIndexed<double> *, RCNeighbourList*, int, 
				 BasisFunction*);
    static void  coarseInter3_3d(DOFIndexed<double> *, RCNeighbourList*, int, 
				 BasisFunction*);
    static void  coarseInter4_1d(DOFIndexed<double> *, RCNeighbourList*, int, 
				 BasisFunction*);
    static void  coarseInter4_2d(DOFIndexed<double> *, RCNeighbourList*, int, 
				 BasisFunction*);
    static void  coarseInter4_3d(DOFIndexed<double> *, RCNeighbourList*, int, 
				 BasisFunction*);
    /** \} */


    /// AbstractFunction which implements lagrange basis functions
    class Phi : public BasFctType
    {
    public:
      /// Constructs the local lagrange basis function for the given position,
      /// positionIndex and nodeIndex. owner_ is a pointer to the Lagrange 
      /// object this basis function belongs to.
      Phi(Lagrange* owner, GeoIndex position, int positionIndex, int nodeIndex);

      /// Destructor
      virtual ~Phi();

    private:
      /// vertices needed for evaluation of this function
      int* vertices;
    
      /// Pointer to the evaluating function
      double (*func)(const DimVec<double>& lambda, int* vert);

      /// Returns \ref func(lambda, vertices)
      double operator()(const DimVec<double>& lambda) const 
      {
	return func(lambda, vertices);
      }

      /** \name basis functions for different degrees
       * \{
       */

      // ====== Lagrange, degree = 0 =====================================
      // center
      static double phi0c(const DimVec<double>&, int*) 
      {
	return 1.0;
      }

      // ====== Lagrange, degree = 1 =====================================
      // vertex
      static double phi1v(const DimVec<double>& lambda, int* vertices) 
      {
	return lambda[vertices[0]]; 
      }

      // ====== Lagrange, degree = 2 =====================================
      // vertex
      static double phi2v(const DimVec<double>& lambda, int* vertices) 
      {
	return lambda[vertices[0]] * (2.0 * lambda[vertices[0]] - 1.0);
      }    

      // edge
      static double phi2e(const DimVec<double>& lambda, int* vertices) 
      {
	return (4.0 * lambda[vertices[0]] * lambda[vertices[1]]);
      }

      // ====== Lagrange, degree = 3 =====================================
      // vertex
      static double phi3v(const DimVec<double>& lambda, int* vertices) 
      {
	return (4.5 * (lambda[vertices[0]] - 1.0) * lambda[vertices[0]] + 1.0) * 
	  lambda[vertices[0]];
      }

      // edge
      static double phi3e(const DimVec<double>& lambda, int* vertices) 
      {
	return (13.5 * lambda[vertices[0]] - 4.5) * 
	  lambda[vertices[0]] * lambda[vertices[1]];
      }

      // face
      static double phi3f(const DimVec<double>& lambda, int* vertices) 
      {
	return 27.0 * lambda[vertices[0]] * lambda[vertices[1]] * 
	  lambda[vertices[2]];
      }

      // ====== Lagrange, degree = 4 ======================================
      // vertex
      static double phi4v(const DimVec<double>& lambda, int* vertices) 
      {
	return (((32.0 * lambda[vertices[0]] - 48.0) * lambda[vertices[0]] + 22.0)
	       * lambda[vertices[0]] - 3.0) * lambda[vertices[0]] / 3.0;
      }

      // edge
      static double phi4e0(const DimVec<double>& lambda, int* vertices) 
      {
	return ((128.0 * lambda[vertices[0]] - 96.0) * lambda[vertices[0]] + 16.0)
	  * lambda[vertices[0]] * lambda[vertices[1]] / 3.0;
      }

      static double phi4e1(const DimVec<double>& lambda, int* vertices) 
      {
	return (4.0 * lambda[vertices[0]] - 1.0) * lambda[vertices[0]] * 
	  (4.0 * lambda[vertices[1]] - 1.0) * lambda[vertices[1]] * 4.0;
      }

      // face
      static double phi4f(const DimVec<double>& lambda,  int* vertices) 
      {
	return (4.0 * lambda[vertices[0]] - 1.0) * lambda[vertices[0]] * 
	  lambda[vertices[1]] * lambda[vertices[2]] * 32.0;
      }

      // center
      static double phi4c(const DimVec<double>& lambda, int* vertices) 
      {
	return 256.0 * lambda[vertices[0]] * lambda[vertices[1]] * 
	  lambda[vertices[2]] * lambda[vertices[3]];
      }

    };
  
    /** \} */



    /// AbstractFunction which implements gradients of lagrange basis functions.
    /// See \ref Phi
    class GrdPhi : public GrdBasFctType
    {
    public:
      GrdPhi(Lagrange* owner, GeoIndex position, int positionIndex, int nodeIndex);

      virtual ~GrdPhi();
    private:
      int* vertices;

      void (*func)(const DimVec<double>& lambda, 
		   int* vertices_, 
		   mtl::dense_vector<double>& result);

      void operator()(const DimVec<double>& lambda, 
			     mtl::dense_vector<double>& result) const 
      {
	func(lambda, vertices, result);
      }

      // ====== Lagrange0 ================================================
      // center
      static void grdPhi0c(const DimVec<double>&, 
				  int*, 
				  mtl::dense_vector<double>& result) 
      {
	result = 0.0;
      }

      // ====== Lagrange1 ================================================
      // vertex
      static void grdPhi1v(const DimVec<double>&, 
				  int* vertices, 
				  mtl::dense_vector<double>& result) 
      {
	result = 0.0;
	result[vertices[0]] = 1.0;
      }

      // ====== Lagrange2 ================================================
      // vertex
      static void grdPhi2v(const DimVec<double>& lambda, 
				  int* vertices, 
				  mtl::dense_vector<double>& result) 
      {
	result = 0.0;
	result[vertices[0]] = 4.0 * lambda[vertices[0]] - 1.0;
      }

      // edge
      static void grdPhi2e(const DimVec<double>& lambda, 
				  int* vertices, 
				  mtl::dense_vector<double>& result) 
      {
	result = 0.0;
	result[vertices[0]] = 4.0 * lambda[vertices[1]];
	result[vertices[1]] = 4.0 * lambda[vertices[0]];
      }

      // ===== Lagrange3 ================================================
      // vertex
      static void grdPhi3v(const DimVec<double>& lambda,
				  int* vertices,
				  mtl::dense_vector<double>& result) 
      {
	result = 0.0;
	result[vertices[0]] = (13.5 * lambda[vertices[0]] - 9.0) * 
	  lambda[vertices[0]] + 1.0;
      }

      // edge
      static void grdPhi3e(const DimVec<double>& lambda,
				  int* vertices, 
				  mtl::dense_vector<double>& result)
      {
	result = 0.0;
	result[vertices[0]] = (27.0 * lambda[vertices[0]] - 4.5) * 
	  lambda[vertices[1]];
	result[vertices[1]] = (13.5 * lambda[vertices[0]] - 4.5) * 
	  lambda[vertices[0]];
      }

      // face
      static void grdPhi3f(const DimVec<double>& lambda, 
				  int* vertices, 
				  mtl::dense_vector<double>& result) 
      {
	result = 0.0;
	result[vertices[0]] = 27.0 * lambda[vertices[1]] * lambda[vertices[2]];
	result[vertices[1]] = 27.0 * lambda[vertices[0]] * lambda[vertices[2]];
	result[vertices[2]] = 27.0 * lambda[vertices[0]] * lambda[vertices[1]];
      }

    
      // ===== Lagrange4 ================================================
      // vertex
      static void grdPhi4v(const DimVec<double>& lambda, 
				  int* vertices, 
				  mtl::dense_vector<double>& result) 
      {
	result = 0.0;
	result[vertices[0]] = 
	  ((128.0 * lambda[vertices[0]] - 144.0) * lambda[vertices[0]] + 44.0) * 
	  lambda[vertices[0]] / 3.0 - 1.0;
      }
    
      // edge
      static void grdPhi4e0(const DimVec<double>& lambda,
				   int* vertices, 
				   mtl::dense_vector<double>& result)
      {
	result = 0.0;
	result[vertices[0]] = ((128.0 * lambda[vertices[0]] - 64.0) * 
			       lambda[vertices[0]] + 16.0 / 3.0) * lambda[vertices[1]];
	result[vertices[1]] = ((128.0 * lambda[vertices[0]] - 96.0) * 
			       lambda[vertices[0]] + 16.0)*lambda[vertices[0]] / 3.0;
      }

      static void grdPhi4e1(const DimVec<double>& lambda,
				   int* vertices,
				   mtl::dense_vector<double>& result)
      {
	result = 0.0;
	result[vertices[0]] = 4.0 * (8.0 * lambda[vertices[0]] - 1.0) * 
	  lambda[vertices[1]] * (4.0 * lambda[vertices[1]] - 1.0);
	result[vertices[1]] = 4.0 * lambda[vertices[0]] * 
	  (4.0 * lambda[vertices[0]] - 1.0) * (8.0 * lambda[vertices[1]] - 1.0);
      }

      // face
      static void grdPhi4f(const DimVec<double>& lambda, 
				  int* vertices, 
				  mtl::dense_vector<double>& result)
      {
	result = 0.0;
	result[vertices[0]] = 32.0 * (8.0 * lambda[vertices[0]] - 1.0) * 
	  lambda[vertices[1]] * lambda[vertices[2]];
	result[vertices[1]] = 32.0 * (4.0 * lambda[vertices[0]] - 1.0) * 
	  lambda[vertices[0]] * lambda[vertices[2]];
	result[vertices[2]] = 32.0 * (4.0 * lambda[vertices[0]] - 1.0) * 
	  lambda[vertices[0]] * lambda[vertices[1]];	
      }
    
      // center
      static void grdPhi4c(const DimVec<double>& lambda, 
				  int* vertices,
				  mtl::dense_vector<double>& result) 
      {
	result = 0.0;
	result[0] = 
	  256.0 * lambda[vertices[1]] * lambda[vertices[2]] * lambda[vertices[3]];
	result[1] = 
	  256.0 * lambda[vertices[0]] * lambda[vertices[2]] * lambda[vertices[3]];
	result[2] = 
	  256.0 * lambda[vertices[0]] * lambda[vertices[1]] * lambda[vertices[3]];
	result[3] = 
	  256.0 * lambda[vertices[0]] * lambda[vertices[1]] * lambda[vertices[2]];
      }
    };



    /// AbstractFunction which implements second derivatives of Lagrange basis
    /// functions. See \ref Phi
    class D2Phi : public D2BasFctType
    {
    public:
      D2Phi(Lagrange* owner, GeoIndex position, int positionIndex, int nodeIndex);

      virtual ~D2Phi();
    private:
      int* vertices;
    
      void (*func)(const DimVec<double>& lambda, int* vertices_, DimMat<double>& result);

      void operator()(const DimVec<double>& lambda, DimMat<double>& result) const {
	return func(lambda, vertices, result);
      }

      // ===== Lagrange0 ================================================
      // center
      static void D2Phi0c(const DimVec<double>&, int*, DimMat<double>& result) 
      {
	result.set(0.0);       
      }

      // ===== Lagrange1 ================================================
      // vertex
      static void D2Phi1v(const DimVec<double>&, int*, DimMat<double>& result) 
      {
	result.set(0.0);
      }

      // ===== Lagrange2 ================================================
      // vertex
      static void D2Phi2v(const DimVec<double>&, int* vertices, 
				 DimMat<double>& result) 
      {
	result.set(0.0);
	result[vertices[0]][vertices[0]] = 4.0;
      }

      // edge
      static void D2Phi2e(const DimVec<double>&, int* vertices, 
				 DimMat<double>& result) 
      {
	result.set(0.0);
	result[vertices[0]][vertices[1]] = 4.0;
	result[vertices[1]][vertices[0]] = 4.0;
      }


      // ===== Lagrange3 ================================================
      // vertex
      static void D2Phi3v(const DimVec<double>& lambda, int* vertices, 
				 DimMat<double>& result) 
      {
	result.set(0.0);
	result[vertices[0]][vertices[0]] = 27.0 * lambda[vertices[0]] - 9.0;
      }

      // edge
      static void D2Phi3e(const DimVec<double>& lambda, int* vertices, 
				 DimMat<double>& result) 
      {
	result.set(0.0);
	result[vertices[0]][vertices[0]] = 27.0 * lambda[vertices[1]];
	result[vertices[0]][vertices[1]] = 
	  result[vertices[1]][vertices[0]] = 27.0 * lambda[vertices[0]] - 4.5;
      }

      // face
      static void D2Phi3f(const DimVec<double>& lambda, int* vertices, 
				 DimMat<double>& result) 
      {
	result.set(0.0);
	result[vertices[0]][vertices[1]] = 
	  result[vertices[1]][vertices[0]] = 27.0 * lambda[vertices[2]];
	result[vertices[0]][vertices[2]] = 
	  result[vertices[2]][vertices[0]] = 27.0 * lambda[vertices[1]];
	result[vertices[1]][vertices[2]] = 
	  result[vertices[2]][vertices[1]] = 27.0 * lambda[vertices[0]];
      }


      // ===== Lagrange4 ================================================
      // vertex
      static void D2Phi4v(const DimVec<double>& lambda, 
				 int* vertices, 
				 DimMat<double>& result) {
	result.set(0.0);
	result[vertices[0]][vertices[0]] = 
	  (128.0 * lambda[vertices[0]] - 96.0) * lambda[vertices[0]] + 44.0 / 3.0;
      }

      // edge
      static void D2Phi4e0(const DimVec<double>& lambda, int* vertices, 
				  DimMat<double>& result) 
      {
	result.set(0.0);
	result[vertices[0]][vertices[0]] = 
	  (256.0 * lambda[vertices[0]] - 64.0) * lambda[vertices[1]];
	result[vertices[0]][vertices[1]] = 
	  result[vertices[1]][vertices[0]] =
	  (128.0 * lambda[vertices[0]] - 64.0) * lambda[vertices[0]] + 16.0 / 3.0;
      }

      static void D2Phi4e1(const DimVec<double>& lambda, int* vertices, 
				  DimMat<double>& result) 
      {
	result.set(0.0);
	result[vertices[0]][vertices[0]] = 
	  32.0 * lambda[vertices[1]] * (4.0 * lambda[vertices[1]] - 1.0);
	result[vertices[0]][vertices[1]] = 
	  result[vertices[1]][vertices[0]] = 
	  4.0 * (8.0 * lambda[vertices[0]] - 1.0) * (8.0 * lambda[vertices[1]] - 1.0);
	result[vertices[1]][vertices[1]] = 
	  32.0 * lambda[vertices[0]] * (4.0 * lambda[vertices[0]] - 1.0);
      }

      // face
      static void D2Phi4f(const DimVec<double>& lambda, int* vertices, 
				 DimMat<double>& result) 
      {
	result.set(0.0);
	result[vertices[0]][vertices[0]] = 
	  256.0 * lambda[vertices[1]] * lambda[vertices[2]];
	result[vertices[0]][vertices[1]] = 
	  result[vertices[1]][vertices[0]] = 
	  32.0 * (8.0 * lambda[vertices[0]] - 1.0) * lambda[vertices[2]];
	result[vertices[0]][vertices[2]] = 
	  result[vertices[2]][vertices[0]] = 
	  32.0 * (8.0 * lambda[vertices[0]] - 1.0) * lambda[vertices[1]];
	result[vertices[1]][vertices[2]] = 
	  result[vertices[2]][vertices[1]] = 
	  32.0 * (4.0 * lambda[vertices[0]] - 1.0) * lambda[vertices[0]];
      }

      // center
      static void D2Phi4c(const DimVec<double>& lambda, int* vertices, 
				 DimMat<double>& result) 
      {
	result.set(0.0);
	result[vertices[0]][vertices[1]] = 
	  result[vertices[1]][vertices[0]] = 
	  256.0 * lambda[vertices[2]] * lambda[vertices[3]];
	result[vertices[0]][vertices[2]] = 
	  result[vertices[2]][vertices[0]] = 
	  256.0 * lambda[vertices[1]] * lambda[vertices[3]];
	result[vertices[0]][vertices[3]] = 
	  result[vertices[3]][vertices[0]] = 
	  256.0 * lambda[vertices[1]] * lambda[vertices[2]];
	result[vertices[1]][vertices[2]] = 
	  result[vertices[2]][vertices[1]] = 
	  256.0 * lambda[vertices[0]] * lambda[vertices[3]];
	result[vertices[1]][vertices[3]] = 
	  result[vertices[3]][vertices[1]] = 
	  256.0 * lambda[vertices[0]] * lambda[vertices[2]];
	result[vertices[2]][vertices[3]] = 
	  result[vertices[3]][vertices[2]] = 
	  256.0 * lambda[vertices[0]] * lambda[vertices[1]];
      }
    };
  };

}

#endif // AMDIS_LAGRANGE_H
