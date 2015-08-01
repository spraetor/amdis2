#include "ElInfo1d.h"
#include "BasisFunction.h"
#include "Element.h"
#include "Line.h"
#include "Triangle.h"
#include "Tetrahedron.h"
#include "FiniteElemSpace.h"
#include "Flag.h"
#include "MacroElement.h"
#include "Mesh.h"
#include "Global.h"
#include "FixVec.h"
#include "DOFVector.h"

namespace AMDiS 
{
  
  double ElInfo1d::mat_d1_val[2][2] = {{1.0, 0.0}, 
				       {0.0, 1.0}};
  mtl::dense2D<double> ElInfo1d::mat_d1(mat_d1_val);

  double ElInfo1d::mat_d1_left_val[2][2] = {{1.0, 0.5}, 
					    {0.0, 0.5}};
  mtl::dense2D<double> ElInfo1d::mat_d1_left(mat_d1_left_val);

  double ElInfo1d::mat_d1_right_val[2][2] = {{0.5, 0.0}, 
					     {0.5, 1.0}};
  mtl::dense2D<double> ElInfo1d::mat_d1_right(mat_d1_right_val);

  double ElInfo1d::mat_d2_val[3][3] = {{1.0, 0.0, 0.0},
				       {0.0, 1.0, 0.0},
				       {0.0, 0.0, 1.0}};
  mtl::dense2D<double> ElInfo1d::mat_d2(mat_d2_val);

  double ElInfo1d::mat_d2_left_val[3][3] = {{1.0, 0.375, 0.0},
					    {0.0, 0.75, 1.0},
					    {0.0, -0.125, 0.0}};
  mtl::dense2D<double> ElInfo1d::mat_d2_left(mat_d2_left_val);
  
  double ElInfo1d::mat_d2_right_val[3][3] = {{0.0, -0.125, 0.0},
					     {1.0, 0.75, 0.0},
					     {0.0, 0.375, 1.0}};
  mtl::dense2D<double> ElInfo1d::mat_d2_right(mat_d2_right_val);
 

  void ElInfo1d::fillMacroInfo(const MacroElement * mel)
  {
    Element *nb;
    MacroElement *mnb;

    macroElement = const_cast<MacroElement*>(mel);
    element = const_cast<Element*>(mel->getElement());
    parent = NULL;
    level = 0;

    int vertices = mesh->getGeo(VERTEX);

    if (fillFlag.isSet(Mesh::FILL_COORDS) || fillFlag.isSet(Mesh::FILL_DET) ||
	fillFlag.isSet(Mesh::FILL_GRD_LAMBDA)) {
      
      for (int i = 0; i < vertices; i++)
	coord[i] = mel->coord[i];
    }

    if (fillFlag.isSet(Mesh::FILL_NEIGH) || fillFlag.isSet(Mesh::FILL_OPP_COORDS)) {
      WorldVector<double> oppC;
      
      int neighbours =  mesh->getGeo(NEIGH);
      for (int i = 0; i < neighbours; i++) {
	nb = NULL;
	if ((mnb = const_cast<MacroElement*>(mel->getNeighbour(i)))) {
	  if (fillFlag.isSet(Mesh::FILL_OPP_COORDS)) {
	    oppC = mnb->coord[i];
	  }

	  nb = const_cast<Element*>(mnb->getElement());

	  while (!(nb->isLeaf())) { // make nb nearest element
	    if (fillFlag.isSet(Mesh::FILL_OPP_COORDS)) {
	      if (nb->isNewCoordSet()) {
		oppC = *(nb->getNewCoord());
	      } else {
		oppC = (mel->coord[i] + oppC) * 0.5;
	      }
	    }
	    nb = const_cast<Element*>(nb->getChild(1-i));
	  }

	  if (fillFlag.isSet(Mesh::FILL_OPP_COORDS)) {
	    oppCoord[i] = oppC;
	  }
	}
	neighbour[i] = nb;
	oppVertex[i] = nb ? i : -1;
      }
    }

    if (fillFlag.isSet(Mesh::FILL_BOUND) ) {
      for (int i = 0; i < vertices; i++)
	boundary[i] = mel->getBoundary(i);

      for (int i = 0; i < element->getGeo(PROJECTION); i++)
	projection[i] = mel->getProjection(i);
    }
  }

  /****************************************************************************/
  /*  compute gradients of basis functions on element; return the absulute    */
  /*  value of the determinante from the transformation to the reference      */
  /*  element                                                                 */
  /****************************************************************************/
  double ElInfo1d::calcGrdLambda(DimVec<WorldVector<double> >& grd_lam)
  {
    FUNCNAME("ElInfo1d::calcGrdLambda()");

    testFlag(Mesh::FILL_COORDS);

    WorldVector<double> e;
    e = coord[1]; 
    e -= coord[0];
    double adet2 = e * e;

    if (adet2 < 1.0E-15) {
      MSG("det*det = %g\n", adet2);
      grd_lam[0] = grd_lam[1] = 0.0;
    } else {
      grd_lam[1] = e * (1.0 / adet2);
      grd_lam[0] = grd_lam[1] * -1.0;
    }

    return std::sqrt(adet2);
  }
  

  int ElInfo1d::worldToCoord(const WorldVector<double>& x,
				   DimVec<double>* lambda) const
  {
    FUNCNAME("ElInfo1d::worldToCoord()");

    double lmin;
    double a = coord[0][0];
    double length = (coord[1][0] - a);
    int dim = mesh->getDim();

    static DimVec<double> vec(dim);

    TEST_EXIT_DBG(lambda)("lambda must not be NULL\n");
    TEST_EXIT_DBG(dim == 1)("dim!=1\n");
    TEST_EXIT_DBG(dimOfWorld == dim)("not yet for DIM != DIM_OF_WORLD\n");

    if (math::abs(length) < DBL_TOL) {
      ERROR_EXIT("length = %le; abort\n", length);
      return 0;
    }

    (*lambda)[1] = (x[0] - a) / length;
    (*lambda)[0] = 1.0 - (*lambda)[1];

    int k = -1;
    lmin = 0.0;
    for (int i = 0; i <= dim; i++) {
      if ((*lambda)[i] < -1.E-5) {
	if ((*lambda)[i] < lmin) {
	  k = i;
	  lmin = (*lambda)[i];
	}
      }
    }

    return k;
  }

  /****************************************************************************/
  /*  calculate a facenormal of edge side of a triangle with coordinates      */
  /*  coord; return the absulute value of the determinant from the           */
  /*  transformation to the reference element                                 */
  /****************************************************************************/
  double ElInfo1d::getNormal(int side, WorldVector<double> &normal) const
  {
    FUNCNAME_DBG("ElInfo::getNormal()");
    normal = coord[side] - coord[(side + 1) % 2];
    double det = norm(normal);
    TEST_EXIT_DBG(det > 1.e-30)("det = 0 on side %d\n", side);
    normal *= 1.0 / det;

    return det;
  }


  /****************************************************************************/
  /*  calculate the normal of the element for dim of world = 2                */
  /*  return the absulute value of the determinant from the                   */
  /*  transformation to the reference element                                 */
  /****************************************************************************/
  double ElInfo1d::getElementNormal(WorldVector<double> &elementNormal) const
  {
    FUNCNAME_DBG("ElInfo::getElementNormal()");

    TEST_EXIT_DBG(dimOfWorld == 2)
      (" element normal only well defined for  DIM_OF_WORLD = DIM + 1 !!");

    elementNormal[0] = coord[1][1] - coord[0][1];
    elementNormal[1] = coord[0][0] - coord[1][0];

    double det = norm(elementNormal);

    TEST_EXIT_DBG(det > 1.e-30)("det = 0");

    elementNormal *= 1.0 / det;
    
    return det;
  }


  void ElInfo1d::fillElInfo(int ichild, const ElInfo *elInfoOld)
  {
    FUNCNAME_DBG("ElInfo1d::fillElInfo()");

    Element *nb;
    Element *elem = elInfoOld->getElement();

    TEST_EXIT_DBG(elem->getChild(0))("no children?\n");
    element = const_cast<Element*>(elem->getChild(ichild));

    TEST_EXIT_DBG(element)("missing child %d?\n", ichild);

    macroElement = elInfoOld->getMacroElement();
    fillFlag = elInfoOld->getFillFlag();
    parent = elem;
    level = elInfoOld->level + 1;
    iChild = ichild;

    int neighbours = mesh->getGeo(NEIGH);

    if (fillFlag.isSet(Mesh::FILL_COORDS) || fillFlag.isSet(Mesh::FILL_DET) ||
	fillFlag.isSet(Mesh::FILL_GRD_LAMBDA)) {

      const FixVec<WorldVector<double>, VERTEX> *old_coord = &(elInfoOld->coord);

      coord[ichild] = (*old_coord)[ichild];
      if (elem->isNewCoordSet())
	coord[1 - ichild] = *(elem->getNewCoord());
      else
	coord[1 - ichild] = ((*old_coord)[0] + (*old_coord)[1]) * 0.5;
    }

    if (fillFlag.isSet(Mesh::FILL_NEIGH) || fillFlag.isSet(Mesh::FILL_OPP_COORDS)) {
      WorldVector<double> oppC;
      
      TEST_EXIT_DBG(fillFlag.isSet(Mesh::FILL_COORDS))
	("FILL_OPP_COORDS only with FILL_COORDS\n");
      
      for (int i = 0; i < neighbours; i++) {
	if (i != ichild) {
	  nb = const_cast<Element*>(elem->getChild(1-ichild));  
	  if (fillFlag.isSet(Mesh::FILL_OPP_COORDS))
	    oppC = elInfoOld->coord[i];
	} else {
	  nb = const_cast<Element*>(elInfoOld->getNeighbour(i));

	  if (nb && fillFlag.isSet(Mesh::FILL_OPP_COORDS))
	    oppC = elInfoOld->oppCoord[i];
	}

	if (nb) {
	  while (nb->getChild(0)) {  // make nb nearest element
	    if (fillFlag.isSet(Mesh::FILL_OPP_COORDS)) {
	      if (nb->isNewCoordSet())
		oppC = *(nb->getNewCoord());
	      else
		oppC = (coord[i] + oppC) * 0.5;
	    }
	    nb = const_cast<Element*>(nb->getChild(1-i));
	  }

	  if (fillFlag.isSet(Mesh::FILL_OPP_COORDS))
	    oppCoord[i] = oppC;
	}
	neighbour[i] = nb;
	oppVertex[i] = nb ? i : -1;
      }
    }

    if (fillFlag.isSet(Mesh::FILL_BOUND)) {
      boundary[ichild] = elInfoOld->getBoundary(ichild);
      boundary[1 - ichild] = INTERIOR;

      if (elInfoOld->getProjection(0) && 
	  elInfoOld->getProjection(0)->getType() == VOLUME_PROJECTION)
	projection[0] = elInfoOld->getProjection(0);      
    }   
  }


  mtl::dense2D<double>& ElInfo1d::getSubElemCoordsMat(int degree) const
  {
    FUNCNAME("ElInfo1d::getSubElemCoordsMat()");

    using namespace mtl;

    if (subElemMatrices[degree].count(std::make_pair(refinementPathLength, refinementPath)) == 0) {
      switch (degree) {
      case 1:
	{
	  dense2D<double> mat(mat_d1);
	  dense2D<double> tmpMat(num_rows(mat), num_rows(mat));
	  
	  for (int i = 0; i < refinementPathLength; i++) {
	    if (refinementPath & (1 << i)) {
	      tmpMat = mat * mat_d1_right;
	      mat = tmpMat;
	    } else  {
	      tmpMat = mat * mat_d1_left;
	      mat = tmpMat;
	    }
	  }

	  subElemMatrices[1][std::make_pair(refinementPathLength, refinementPath)] = mat;
	}
	
	break;	
      case 2:
	{
	  dense2D<double> mat(mat_d2);
	  dense2D<double> tmpMat(num_rows(mat), num_rows(mat));
	  
	  for (int i = 0; i < refinementPathLength; i++) {
	    if (refinementPath & (1 << i)) {
	      tmpMat = mat * mat_d2_right;
	      mat = tmpMat;
	    } else  {
	      tmpMat = mat * mat_d2_left;
	      mat = tmpMat;
	    }
	  }

	  subElemMatrices[2][std::make_pair(refinementPathLength, refinementPath)] = mat;  
	}
	break;
      default:
	ERROR_EXIT("Not supported for basis function degree: %d\n", degree);
      }
    }

    return subElemMatrices[degree][std::make_pair(refinementPathLength, refinementPath)];
  }

} // enamesace AMDiS
