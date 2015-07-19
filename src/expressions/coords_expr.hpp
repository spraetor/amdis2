/** \file coords_expr.hpp */

#pragma once

#include "AMDiS_fwd.h"
#include "Mesh.h"
#include "LazyOperatorTerm.h"

namespace AMDiS 
{    
  namespace expressions 
  {	
    /// Expression that representy the coordinate vector
    struct Coords : public LazyOperatorTermBase
    {
      typedef WorldVector<double> value_type;
      mutable mtl::dense_vector<WorldVector<double> > x;

      Coords() {}

      template <class List>
      void insertFeSpaces(List& feSpaces) const {}
      
      int getDegree() const { return 1; }

      template <class OT>
      void initElement(OT* ot, const ElInfo* elInfo,
		       SubAssembler* subAssembler, Quadrature *quad, 
		       const BasisFunction *basisFct = NULL)
      {
	if (subAssembler) {
	  subAssembler->getCoordsAtQPs(elInfo, quad, x); 
	}
	else if (quad) {
	  const int nPoints = quad->getNumPoints();
	
	  x.change_dim(nPoints);
	  for (int i = 0; i < nPoints; i++)
	    elInfo->coordToWorld(quad->getLambda(i), x[i]);
	}
	else if (basisFct) {
	  const int nBasisFct = basisFct->getNumber();
	
	  x.change_dim(nBasisFct);
	  for (int i = 0; i < nBasisFct; i++)
	    elInfo->coordToWorld(*basisFct->getCoords(i), x[i]);
	}
      }

      value_type operator()(const int& iq) const { return x[iq]; }
      
      std::string str() { return "X"; }
    };
    
    
    /// Expression that represents the Ith component of the coordinate vector
    template <int I_>
    struct Coord : public LazyOperatorTermBase
    {
      typedef double value_type;
      mutable mtl::dense_vector<WorldVector<double> > x;
      int I;

      Coord(int i = -1) : I(i >= 0 ? i : I_) { }

      template <class List>
      void insertFeSpaces(List& feSpaces) const {}
      
      int getDegree() const { return 1; }

      template <class OT>
      void initElement(OT* ot, const ElInfo* elInfo,
		       SubAssembler* subAssembler, Quadrature *quad, 
		       const BasisFunction *basisFct = NULL)
      {
	if (subAssembler) {
	  subAssembler->getCoordsAtQPs(elInfo, quad, x); 
	}
	else if (quad) {
	  const int nPoints = quad->getNumPoints();
	
	  x.change_dim(nPoints);
	  for (int i = 0; i < nPoints; i++)
	    elInfo->coordToWorld(quad->getLambda(i), x[i]);
	}
	else if (basisFct) {
	  const int nBasisFct = basisFct->getNumber();
	
	  x.change_dim(nBasisFct);
	  for (int i = 0; i < nBasisFct; i++)
	    elInfo->coordToWorld(*basisFct->getCoords(i), x[i]);
	}
      }

      double operator()(const int& iq) const { return x[iq][I]; }
      
      std::string str() { return std::string("X<") + std::to_string(I) + ">"; }
    };
    
    
    /// Expression that represents the element normal vector
    /** In the constructor the corresponding face must be given, of which 
	you want to calculate the normal vector. **/
    struct Normals : public LazyOperatorTermBase
    {
      typedef WorldVector<double> value_type;
      mutable WorldVector<double> normal;
      int boundary;

      Normals(int boundary_) : boundary(boundary_) {}

      template <class List>
      void insertFeSpaces(List& feSpaces) const {}
      
      int getDegree() const { return 1; }

      template <class OT>
      void initElement(OT* ot, const ElInfo* elInfo,
		       SubAssembler* subAssembler, Quadrature *quad, 
		       const BasisFunction *basisFct = NULL)
      {
	int dim = elInfo->getMesh()->getDim();
	for (int side = 0; side < dim+1; ++side) {
	  if (elInfo->getBoundary(side) == boundary) {
	    elInfo->getNormal(side, normal);
	    break;
	  }
	}
      }

      value_type operator()(const int& iq) const { return normal; }
      
      std::string str() { return "N"; }
    };
    
    
    /// Expression that represents the Ith component of an element normal vector
    /** In the constructor the corresponding face must be given, of which 
	you want to calculate the normal vector. **/
    struct Normal : public LazyOperatorTermBase
    {
      typedef double value_type;
      mutable WorldVector<double> normal;
      BoundaryType boundary;
      int I;

      Normal(BoundaryType boundary_, int I_) : boundary(boundary_), I(I_) {}

      template <class List>
      void insertFeSpaces(List& feSpaces) const {}
      
      int getDegree() const { return 1; }

      template <class OT>
      void initElement(OT* ot, const ElInfo* elInfo,
		       SubAssembler* subAssembler, Quadrature *quad, 
		       const BasisFunction *basisFct = NULL)
      {
	int dim = elInfo->getMesh()->getDim();
	for (int side = 0; side < dim+1; ++side) {
	  if (elInfo->getBoundary(side) == boundary) {
	    elInfo->getNormal(side, normal);
	    break;
	  }
	}
      }

      value_type operator()(const int& iq) const { return normal[I]; }
      
      std::string str() { return std::string("N<") + std::to_string(I) + ">"; }
    };  
    
    
    /// Expression that represents the surface normal vector
    struct ElementNormals : public LazyOperatorTermBase
    {
      typedef WorldVector<double> value_type;
      mutable WorldVector<double> elementNormal;

      ElementNormals() {}

      template <class List>
      void insertFeSpaces(List& feSpaces) const {}
      
      int getDegree() const { return 1; }

      template <class OT>
      void initElement(OT* ot, const ElInfo* elInfo,
		       SubAssembler* subAssembler, Quadrature *quad, 
		       const BasisFunction *basisFct = NULL)
      {
	elInfo->getElementNormal(elementNormal);
      }

      value_type operator()(const int& iq) const { return elementNormal; }
      
      std::string str() { return "M"; }
    };
    
    
    /// Expression that represents the Ith component of a surface normal vector
    struct ElementNormal : public LazyOperatorTermBase
    {
      typedef double value_type;
      mutable WorldVector<double> elementNormal;
      int I;

      ElementNormal(int I_) : I(I_) {}

      template <class List>
      void insertFeSpaces(List& feSpaces) const {}
      
      int getDegree() const { return 1; }

      template <class OT>
      void initElement(OT* ot, const ElInfo* elInfo,
		       SubAssembler* subAssembler, Quadrature *quad, 
		       const BasisFunction *basisFct = NULL)
      {
	elInfo->getElementNormal(elementNormal);
      }

      value_type operator()(const int& iq) const { return elementNormal[I]; }
      
      std::string str() { return std::string("M<") + std::to_string(I) + ">"; }
    };

  } // end namespace expressions
  

  inline expressions::Coords X() 
  { return expressions::Coords(); }
  
  template <int I>
  inline expressions::Coord<I> X() 
  { return expressions::Coord<I>(); }
  
  inline expressions::Coord<-1> X(int i) 
  { return expressions::Coord<-1>(i); }
  

  inline expressions::Normals N(int side) 
  { return expressions::Normals(side); }
  
  inline expressions::Normal N(int side, int I) 
  { return expressions::Normal(side, I); }
  

  inline expressions::ElementNormals M() 
  { return expressions::ElementNormals(); }
  
  inline expressions::ElementNormal M(int I) 
  { return expressions::ElementNormal(I); }


} // end namespace AMDiS
