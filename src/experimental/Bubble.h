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

// Created by Roman Weissflog & Philipp Schulz



/** \file Bubble.h */

#ifndef AMDIS_BUBBLE_H
#define AMDIS_BUBBLE_H

#include <list>
#include <boost/numeric/mtl/mtl.hpp>
#include "AbstractFunction.h"
#include "BasisFunction.h"
#include "FixVec.h"

namespace AMDiS {

  /** \ingroup FEMSpace
  * \brief
  * Lagrange basis functions plus Bubble function. Sub class of BasisFunction
  */

  class Bubble : public BasisFunction
  {
  public:
    /// Creator class used in the BasisFunctionCreatorMap.
    class Creator : public BasisFunctionCreator
    {
      public:
	virtual ~Creator() {}
		
	/// Returns a new Lagrange object.
	BasisFunction* create()
	{
	  return getBubble(this->dim, this->dim + 1); // has to be generalized in later versions of bubble functions
	}
    };
    
    protected:
      /// Constructs lagrange/bubble basis functions with the given dim and degree.
      /// Constructor is protected to avoid multiple instantiation of identical
      /// basis functions. Use \ref getBubble instead.
      Bubble(int dim_, int degree_);

      /** \brief
      * destructor
      */
      virtual ~Bubble();

    public:
      /// Returns a pointer to lagrange and bubble basis functions with the given dim and
      /// degree. Multiple instantiation of identical basis functions is avoided
      /// by rembering once created basis functions in \ref allBasFcts.
      static Bubble* getBubble(int dim, int degree);

      /// Implements BasisFunction::interpol
      void interpol(const ElInfo *, int, const int *, 
		    AbstractFunction<double, WorldVector<double> >*, 
		    mtl::dense_vector<double>&) const override;

      /// Implements BasisFunction::interpol
      void interpol(const ElInfo *, int, 
		  const int *b_no,
		  AbstractFunction<WorldVector<double>, WorldVector<double> >*, 
		  mtl::dense_vector<WorldVector<double> >&) const override;

      /// Returns the barycentric coordinates of the i-th basis function.
      DimVec<double> *getCoords(int i) const override;


      /// Implements BasisFunction::getBound
      void getBound(const ElInfo*, BoundaryType *) const override;


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
      inline void refineInter(DOFIndexed<double> *drv, RCNeighbourList* list, int n) override
      {
	if (refineInter_fct)
	  (*refineInter_fct)(drv, list, n, this);
      }

      /// Implements BasisFunction::coarseRestrict
      inline void coarseRestr(DOFIndexed<double> *drv, RCNeighbourList* list, int n) override
      {
	if (coarseRestr_fct)
	  (*coarseRestr_fct)(drv, list, n, this);
      }

      /// Implements BasisFunction::coarseInter
      inline void coarseInter(DOFIndexed<double> *drv, RCNeighbourList* list, int n) override
      {
	if (coarseInter_fct)
	  (*coarseInter_fct)(drv, list, n, this);
      }

      /// Implements BasisFunction::getLocalIndices().
      void getLocalIndices(const Element *el,
			  const DOFAdmin *admin,
			  std::vector<DegreeOfFreedom> &dofs) const override;

      /// Implements BasisFunction::getLocalDofPtrVec()
      /// Returns an vector filled with all DOFs per position
      void getLocalDofPtrVec(const Element *el, 
			  const DOFAdmin *admin,
			  std::vector<const DegreeOfFreedom*>& vec) const override;

      /// Implements BasisFunction::l2ScpFctBas
      void l2ScpFctBas(Quadrature* q,
		      AbstractFunction<double, WorldVector<double> >* f,
		      DOFVector<double>* fh) override;

      /// Implements BasisFunction::l2ScpFctBas
      void l2ScpFctBas(Quadrature* q,
		      AbstractFunction<WorldVector<double>, WorldVector<double> >* f,
		      DOFVector<WorldVector<double> >* fh) override;

      static void clear();

      /// Implements BasisFunction::isnodal
      bool isNodal() const override
      {
	return false;
      }


	
    protected:
      /// sets the barycentric coordinates (stored in \ref bary) of the local 
      /// basis functions.
      void setBary();		

      /// Implements BasisFunction::setNDOF
      void setNDOF() override;		

      /// Sets used function pointers
      void setFunctionPointer();

      /// Used by \ref getVec
      int* orderOfPositionIndices(const Element* el, 
				GeoIndex position, 
				int positionIndex) const override;

    private:
      /// barycentric coordinates of the locations of all basis functions
      std::vector<DimVec<double>* > *bary;

      /** \name static dim-degree-arrays
      * \{
      */
      static std::vector<DimVec<double>* > baryDimDegree;
      static DimVec<int>* ndofDimDegree;
      static int nBasFctsDimDegree;
      static std::vector<BasFctType*> phifunc;
      static std::vector<GrdBasFctType*> grdPhifunc;
      static std::vector<D2BasFctType*> D2Phifunc;
      /** \} */

      /// List of all used BasisFunctions in the whole program. Avoids duplicate
      /// instantiation of identical BasisFunctions. 
      static Bubble* Singleton;


    protected:
      /// Pointer to the used refineInter function
      void (*refineInter_fct)(DOFIndexed<double> *, RCNeighbourList*, int, BasisFunction*);

      static void  refineInter2_1d(DOFIndexed<double> *, RCNeighbourList*, int, 
				  BasisFunction*);
      static void  refineInter3_2d(DOFIndexed<double> *, RCNeighbourList*, int, 
				  BasisFunction*);
      static void  refineInter4_3d(DOFIndexed<double> *, RCNeighbourList*, int, 
				  BasisFunction*);

      /// Pointer to the used coarseRestr function
      void (*coarseRestr_fct)(DOFIndexed<double> *, RCNeighbourList*, int, BasisFunction*);

      static void  coarseRestr2_1d(DOFIndexed<double> *, RCNeighbourList*, int, 
				  BasisFunction*);
      static void  coarseRestr3_2d(DOFIndexed<double> *, RCNeighbourList*, int, 
				  BasisFunction*);
      static void  coarseRestr4_3d(DOFIndexed<double> *, RCNeighbourList*, int, 
				  BasisFunction*);

      /// Pointer to the used coarseInter function
      void (*coarseInter_fct)(DOFIndexed<double> *, RCNeighbourList*, int, BasisFunction*);  

      static void  coarseInter2_1d(DOFIndexed<double> *, RCNeighbourList*, int, 
				  BasisFunction*);
      static void  coarseInter3_2d(DOFIndexed<double> *, RCNeighbourList*, int, 
				  BasisFunction*);
      static void  coarseInter4_3d(DOFIndexed<double> *, RCNeighbourList*, int, 
				  BasisFunction*);


      /// AbstractFunction which implements lagrange/bubble basis functions
      class Phi : public BasFctType
      {
	public:
	  /// Constructs the local lagrange/bubble basis function for the given position,
	  /// positionIndex and nodeIndex. owner_ is a pointer to the Bubble 
	  /// object this basis function belongs to.
	  Phi(Bubble* owner, GeoIndex position, int positionIndex, int nodeIndex);

	  /// Destructor
	  virtual ~Phi();

	private:
	  /// vertices needed for evaluation of this function
	  int* vertices;

	  /// Pointer to the evaluating function
	  double (*func)(const DimVec<double>& lambda, int* vert);

	  /// Returns \ref func(lambda, vertices)
	  inline double operator()(const DimVec<double>& lambda) const override
	  {
	    return func(lambda, vertices);
	  }

	  // ====== Lagrange, degree = 1 =====================================
	  // vertex
	  inline static double phi1v(const DimVec<double>& lambda, int* vertices) 
	  {
	    return lambda[vertices[0]]; 
	  }

	  // ====== Bubble ===================================================
          // 1d
          inline static double phi2c(const DimVec<double>& lambda, int* vertices) 
          {
	    return (4.0 * lambda[vertices[0]] * lambda[vertices[1]]);
          }

          // 2d
	  inline static double phi3c(const DimVec<double>& lambda, int* vertices) 
	  {
	    return 27.0 * lambda[vertices[0]] * lambda[vertices[1]] * lambda[vertices[2]];
	  }

          // 3d
	  inline static double phi4c(const DimVec<double>& lambda, int* vertices) 
	  {
	    return 256.0 * lambda[vertices[0]] * lambda[vertices[1]] * 
	      lambda[vertices[2]] * lambda[vertices[3]];
	  }
      };


      /// AbstractFunction which implements gradients of Lagrange/Bubble basis functions.
      /// See \ref Phi
      class GrdPhi : public GrdBasFctType
      {	
	public:
	  GrdPhi(Bubble* owner, GeoIndex position, int positionIndex, int nodeIndex);

	  virtual ~GrdPhi();
	private:
	  int* vertices;
  
	  void (*func)(const DimVec<double>& lambda, 
			int* vertices_, 
			mtl::dense_vector<double>& result);

	  inline void operator()(const DimVec<double>& lambda, 
				mtl::dense_vector<double>& result) const override
	  {
	    func(lambda, vertices, result);
	  }

	// ====== Lagrange1 ================================================
	// vertex
	inline static void grdPhi1v(const DimVec<double>&, 
				  int* vertices, 
				  mtl::dense_vector<double>& result) 
	{
	  result = 0.0;
	  result[vertices[0]] = 1.0;
	}

	// ======= Bubble ==================================================
        inline static void grdPhi2c(const DimVec<double>& lambda, 
	  			  int* vertices, 
				  mtl::dense_vector<double>& result) 
        {
	  result = 0.0;
	  result[vertices[0]] = 4.0 * lambda[vertices[1]];
	  result[vertices[1]] = 4.0 * lambda[vertices[0]];
        }


	inline static void grdPhi3c(const DimVec<double>& lambda, 
				      int* vertices, 
				      mtl::dense_vector<double>& result) 
	{
	  result = 0.0;
	  result[vertices[0]] = 27.0 * lambda[vertices[1]] * lambda[vertices[2]];
	  result[vertices[1]] = 27.0 * lambda[vertices[0]] * lambda[vertices[2]];
	  result[vertices[2]] = 27.0 * lambda[vertices[0]] * lambda[vertices[1]];
	}

        inline static void grdPhi4c(const DimVec<double>& lambda, 
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


      /// AbstractFunction which implements second derivatives of Lagrange/Bubble basis
      /// functions. See \ref Phi
      class D2Phi : public D2BasFctType
      {
      public:
	D2Phi(Bubble* owner, GeoIndex position, int positionIndex, int nodeIndex);

	virtual ~D2Phi();

      private:
	int* vertices;

	void (*func)(const DimVec<double>& lambda, int* vertices_, DimMat<double>& result);

	inline void operator()(const DimVec<double>& lambda, DimMat<double>& result) const override
	{
	  return func(lambda, vertices, result);
	}

	// ===== Lagrange1 ================================================
	// vertex
	inline static void D2Phi1v(const DimVec<double>&, int*, DimMat<double>& result) 
	{
	  result.set(0.0);
	}


	// ===== Bubble ===================================================
        inline static void D2Phi2c(const DimVec<double>&, int* vertices, 
				 DimMat<double>& result) 
        {
	  result.set(0.0);
	  result[vertices[0]][vertices[1]] = 4.0;
	  result[vertices[1]][vertices[0]] = 4.0;
        }

	inline static void D2Phi3c(const DimVec<double>& lambda, int* vertices, 
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

        inline static void D2Phi4c(const DimVec<double>& lambda, int* vertices, 
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

#endif // AMDIS_Bubble_H
