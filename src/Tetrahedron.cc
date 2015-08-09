#include "Tetrahedron.h"
#include "DOFAdmin.h"
#include "Mesh.h"
#include "CoarseningManager.h"
#include "FixVec.h"
#include "ElementDofIterator.h"
#include "Global.h"
#include "FiniteElemSpace.h"

using namespace std;

namespace AMDiS 
{
  constexpr unsigned char Tetrahedron::nChildEdge[3][2][2];
  constexpr unsigned char Tetrahedron::nChildFace[3][2][2];
  constexpr int Tetrahedron::childVertex[3][2][4];
  constexpr int Tetrahedron::childEdge[3][2][6];
  constexpr unsigned char Tetrahedron::adjacentChild[2][2];
  constexpr signed char Tetrahedron::childOrientation[3][2];
  constexpr unsigned char Tetrahedron::edgeOfDofs[4][4];
  constexpr int Tetrahedron::vertexOfEdge[6][2];
  constexpr int Tetrahedron::vertexOfFace[4][3];
  constexpr int Tetrahedron::sideOfChild[3][2][4];
  constexpr int Tetrahedron::edgeOfChild[3][2][6];
  constexpr int Tetrahedron::vertexOfParent[3][2][4];
  constexpr int Tetrahedron::edgeOfFace[4][3];


  bool Tetrahedron::hasSide(Element* sideElem) const
  {
    FUNCNAME("Tetrahedron::hasSide()");
    TEST_EXIT_DBG(sideElem->isTriangle())("called for sideElem-type != Triangle\n");
    ERROR_EXIT("not yet\n");
    return false;
  }


  int Tetrahedron::getVertexOfPosition(GeoIndex position,
                                       int positionIndex,
                                       int vertexIndex) const 
  {
    FUNCNAME("Triangle::getVertexOfPosition()");
    switch(position) {
    case VERTEX:
      return positionIndex;
      break;
    case EDGE:
      return vertexOfEdge[positionIndex][vertexIndex];
      break;
    case FACE:
      return vertexOfFace[positionIndex][vertexIndex];
      break;
    case CENTER:
      return vertexIndex;
      break;
    default:
      ERROR_EXIT("invalid position\n");
      return 0;
    }
  }


  void Tetrahedron::sortFaceIndices(int face, FixVec<int,WORLD> &vec) const
  {
    // TODO: REMOVE STATIC
    static MatrixOfFixVecs<FixVec<int,WORLD> > *sorted_3d = NULL;

    if (sorted_3d == NULL) {
      sorted_3d = new MatrixOfFixVecs<FixVec<int,WORLD> >(3, 4, 7, NO_INIT);
 
      (*sorted_3d)[0][0][0] = (*sorted_3d)[0][0][1] =
	(*sorted_3d)[0][0][2] = (*sorted_3d)[1][0][0] =
	(*sorted_3d)[1][0][1] = (*sorted_3d)[1][0][2] =
	(*sorted_3d)[1][1][0] = (*sorted_3d)[1][2][1] =
	(*sorted_3d)[1][3][0] = (*sorted_3d)[1][4][2] =
	(*sorted_3d)[1][5][1] = (*sorted_3d)[1][6][2] =
	(*sorted_3d)[2][0][0] = (*sorted_3d)[2][0][1] =
	(*sorted_3d)[2][0][2] = (*sorted_3d)[2][1][0] =
	(*sorted_3d)[2][2][1] = (*sorted_3d)[2][3][0] =
	(*sorted_3d)[2][4][2] = (*sorted_3d)[2][5][1] =
	(*sorted_3d)[2][6][2] = (*sorted_3d)[3][0][0] =
	(*sorted_3d)[3][0][1] = (*sorted_3d)[3][0][2] =
	(*sorted_3d)[3][1][0] = (*sorted_3d)[3][2][1] =
	(*sorted_3d)[3][3][0] = (*sorted_3d)[3][4][2] =
	(*sorted_3d)[3][5][1] = (*sorted_3d)[3][6][2] = 0;

      (*sorted_3d)[0][1][0] = (*sorted_3d)[0][2][1] =
	(*sorted_3d)[0][3][0] = (*sorted_3d)[0][4][2] =
	(*sorted_3d)[0][5][1] = (*sorted_3d)[0][6][2] =
	(*sorted_3d)[2][1][2] = (*sorted_3d)[2][2][0] =
	(*sorted_3d)[2][3][1] = (*sorted_3d)[2][4][1] =
	(*sorted_3d)[2][5][2] = (*sorted_3d)[2][6][0] =
	(*sorted_3d)[3][1][2] = (*sorted_3d)[3][2][0] =
	(*sorted_3d)[3][3][1] = (*sorted_3d)[3][4][1] =
	(*sorted_3d)[3][5][2] = (*sorted_3d)[3][6][0] = 1;

      (*sorted_3d)[0][1][2] = (*sorted_3d)[0][2][0] =
	(*sorted_3d)[0][3][1] = (*sorted_3d)[0][4][1] =
	(*sorted_3d)[0][5][2] = (*sorted_3d)[0][6][0] =
	(*sorted_3d)[1][1][2] = (*sorted_3d)[1][2][0] =
	(*sorted_3d)[1][3][1] = (*sorted_3d)[1][4][1] =
	(*sorted_3d)[1][5][2] = (*sorted_3d)[1][6][0] =
	(*sorted_3d)[3][1][1] = (*sorted_3d)[3][2][2] =
	(*sorted_3d)[3][3][2] = (*sorted_3d)[3][4][0] =
	(*sorted_3d)[3][5][0] = (*sorted_3d)[3][6][1] = 2;

      (*sorted_3d)[0][1][1] = (*sorted_3d)[0][2][2] =
	(*sorted_3d)[0][3][2] = (*sorted_3d)[0][4][0] =
	(*sorted_3d)[0][5][0] = (*sorted_3d)[0][6][1] =
	(*sorted_3d)[1][1][1] = (*sorted_3d)[1][2][2] =
	(*sorted_3d)[1][3][2] = (*sorted_3d)[1][4][0] =
	(*sorted_3d)[1][5][0] = (*sorted_3d)[1][6][1] =
	(*sorted_3d)[2][1][1] = (*sorted_3d)[2][2][2] =
	(*sorted_3d)[2][3][2] = (*sorted_3d)[2][4][0] =
	(*sorted_3d)[2][5][0] = (*sorted_3d)[2][6][1] = 3;
    }

    const int *vof = vertexOfFace[face];
    int no = 0;
    if (dof[vof[0]][0] < dof[vof[1]][0])
      no++;
    if (dof[vof[1]][0] < dof[vof[2]][0])
      no += 2;
    if (dof[vof[2]][0] < dof[vof[0]][0])
      no += 4;

    TEST_EXIT_DBG(no >= 1 && no <= 6)
      ("Cannot sort face indices of element %d at face %d\n", getIndex(), face);

    vec = (*sorted_3d)[face][no];
  }


  void Tetrahedron::getNodeDofs(const FiniteElemSpace* feSpace, 
				BoundaryObject bound,
				DofContainer& dofs,
				bool baseDofPtr) const
  {
    FUNCNAME("Tetrahedron::getNodeDofs()");

    switch (bound.subObj) {
    case VERTEX:
      {
	if (baseDofPtr) {
	  dofs.push_back(dof[bound.ithObj]);
	} else {
	  int n0 = feSpace->getAdmin()->getNumberOfPreDofs(VERTEX);
	  dofs.push_back(&(dof[bound.ithObj][n0]));
	}
      }
      break;
    case EDGE:
      getNodeDofsAtEdge(feSpace, bound, dofs, baseDofPtr);
      break;
    case FACE:
      getNodeDofsAtFace(feSpace, bound, dofs, baseDofPtr);
      break;      
    default:
      ERROR_EXIT("Should not happen!\n");
    }
  }


  void Tetrahedron::getNodeDofsAtFace(const FiniteElemSpace* feSpace, 
				      BoundaryObject bound,
				      DofContainer& dofs,
				      bool baseDofPtr) const
  {
    FUNCNAME("Tetrahedron::getNodeDofsAtFace()");

    if (!child[0]){
      return;
	}

    switch (bound.ithObj) {
    case 0:
      {
	BoundaryObject nextBound = bound;
	prepareNextBound(nextBound, 1);
	child[1]->getNodeDofs(feSpace, nextBound, dofs, baseDofPtr);
      }
      break;
    case 1:
      {
	BoundaryObject nextBound = bound;
	prepareNextBound(nextBound, 0);
	child[0]->getNodeDofs(feSpace, nextBound, dofs, baseDofPtr);
      }
      break;

    case 2:
    case 3:
      {
	int n0 = (baseDofPtr ? 0 : feSpace->getAdmin()->getNumberOfPreDofs(VERTEX));
	BoundaryObject nextBound0 = bound, nextBound1 = bound;
	prepareNextBound(nextBound0, 0);
	prepareNextBound(nextBound1, 1);
	nextBound0.reverseMode = false;
	nextBound1.reverseMode = false;
	
	bool addDof = true;
	for (ExcludeList::iterator it = bound.excludedSubstructures.begin(); 
	     it != bound.excludedSubstructures.end(); ++it)
	  if (it->first == EDGE && it->second == 0)
	    addDof = false;	
	
	if (bound.reverseMode) {
	  child[1]->getNodeDofs(feSpace, nextBound1, dofs, baseDofPtr);
	  
	  if (addDof)
	    dofs.push_back(&(child[0]->getDof()[3][n0]));
	  
	  child[0]->getNodeDofs(feSpace, nextBound0, dofs, baseDofPtr);
	} else {
	  child[0]->getNodeDofs(feSpace, nextBound0, dofs, baseDofPtr);
	  
	  if (addDof)
	    dofs.push_back(&(child[0]->getDof()[3][n0]));
	  
	  child[1]->getNodeDofs(feSpace, nextBound1, dofs, baseDofPtr);
	}
      }
      break;
    default:
      ERROR_EXIT("Should never happen: %d\n", bound.ithObj);
    }
  }


  void Tetrahedron::getNodeDofsAtEdge(const FiniteElemSpace* feSpace, 
				      BoundaryObject bound,
				      DofContainer& dofs,
				      bool baseDofPtr) const
  {
    FUNCNAME_DBG("Tetrahedron::getNodeDofsAtEdge()");


    if (!child[0]){
      return;
	}

    int n0 = (baseDofPtr ? 0 : feSpace->getAdmin()->getNumberOfPreDofs(VERTEX));
    BoundaryObject nextBound0 = bound, nextBound1 = bound;
    prepareNextBound(nextBound0, 0);
    prepareNextBound(nextBound1, 1);

    switch (bound.ithObj) {
    case 0:      
      nextBound0.reverseMode = false;
      nextBound1.reverseMode = false;

      if (bound.reverseMode) {
			child[1]->getNodeDofs(feSpace, nextBound0, dofs, baseDofPtr);
			dofs.push_back(&(child[0]->getDof()[3][n0]));
			child[0]->getNodeDofs(feSpace, nextBound1, dofs, baseDofPtr);

      } else {
			child[0]->getNodeDofs(feSpace, nextBound0, dofs, baseDofPtr);
			dofs.push_back(&(child[0]->getDof()[3][n0]));
			child[1]->getNodeDofs(feSpace, nextBound1, dofs, baseDofPtr);

      }

      break;
    case 5:
      TEST_EXIT_DBG(nextBound0.ithObj == nextBound1.ithObj)
	("Should not happen!\n");

      if (nextBound0.ithObj != -1)
	child[0]->getNodeDofs(feSpace, nextBound0, dofs, baseDofPtr);
      break;
    default:
      TEST_EXIT_DBG(nextBound0.ithObj == -1 || nextBound1.ithObj == -1)
	("This should not happen!\n");
      
      if (nextBound0.ithObj != -1)
	child[0]->getNodeDofs(feSpace, nextBound0, dofs, baseDofPtr);

      if (nextBound1.ithObj != -1)
	child[1]->getNodeDofs(feSpace, nextBound1, dofs, baseDofPtr);
    }

  }


  void Tetrahedron::getHigherOrderDofs(const FiniteElemSpace* feSpace,
				       BoundaryObject bound,
				       DofContainer& dofs,
				       bool baseDofPtr,
				       vector<GeoIndex>* dofGeoIndex) const
  {
    FUNCNAME("Tetrahedron::getHigherOrderDofs()");

    switch (bound.subObj) {
    case VERTEX:
      return;
      break;
    case EDGE:
      {
	// === Create boundary information objects for children elements. ===

	BoundaryObject nextBound0 = bound;
	prepareNextBound(nextBound0, 0);

	BoundaryObject nextBound1 = bound;
	prepareNextBound(nextBound1, 1);

	// === Check for boundary on children elements. ===

	if ((nextBound0.ithObj >= 0 || nextBound1.ithObj >= 0) &&  child[0]) {
	  // So, the edge is contained in at least on of the children and the
	  // element is also refined. Then we have go down further in refinement
	  // hierarchie.

	 	// Additonal Check if Boundary Edge is NOT Refinement Edge of BoundaryObject
	 	// AND the complete Edge is part of both childs
		// therefore only recursion on ONE child is necessary
	  if ( (bound.ithObj > 0) && (nextBound0.ithObj >= 0 && nextBound1.ithObj >= 0) ){
			child[0]->getHigherOrderDofs(feSpace, nextBound0, dofs,
						   baseDofPtr, dofGeoIndex);
	  } else {
	  	if (bound.reverseMode) {
	    	if (nextBound1.ithObj >= 0)
	    	  child[1]->getHigherOrderDofs(feSpace, nextBound1, dofs, 
						   baseDofPtr, dofGeoIndex);
	    	if (nextBound0.ithObj >= 0)
	    	  child[0]->getHigherOrderDofs(feSpace, nextBound0, dofs,
						   baseDofPtr, dofGeoIndex);
		  } else {
		    if (nextBound0.ithObj >= 0)
		      child[0]->getHigherOrderDofs(feSpace, nextBound0, dofs,
						   baseDofPtr, dofGeoIndex);
		    if (nextBound1.ithObj >= 0)
		      child[1]->getHigherOrderDofs(feSpace, nextBound1, dofs,
						   baseDofPtr, dofGeoIndex);
		  }
	  }

	} else {
	  // Either the edge is not contained in further refined children, or
	  // the element is not refined further on this edge. Then we can get
	  // all the DOFs on this edge.

	  ElementDofIterator elDofIter(feSpace, true);
	  elDofIter.reset(this);

	  if (baseDofPtr) {
	    do {
	      if (elDofIter.getCurrentPos() == 1 && 
		  elDofIter.getCurrentElementPos() == bound.ithObj) {
		dofs.push_back(elDofIter.getBaseDof());	
		if (dofGeoIndex != NULL)
		  dofGeoIndex->push_back(EDGE);
	      }
	    } while (elDofIter.nextStrict());
	  } else {
	    do {
	      if (elDofIter.getCurrentPos() == 1 && 
		  elDofIter.getCurrentElementPos() == bound.ithObj) {
				dofs.push_back(elDofIter.getDofPtr());	

				if (dofGeoIndex != NULL)
				  dofGeoIndex->push_back(EDGE);
	      }
	    } while (elDofIter.next());
	  }
	}
      }

      break;
    case FACE:
      {      
	if (child[0]) {
	  BoundaryObject nextBound0 = bound, nextBound1 = bound;
	  prepareNextBound(nextBound0, 0);
	  prepareNextBound(nextBound1, 1);

	  TEST_EXIT_DBG(nextBound0.ithObj != -1 || nextBound1.ithObj != 1)
	    ("No new face for child elements!\n");

	  TEST_EXIT_DBG(!bound.reverseMode)("Not yet implemented!\n");
	  
	  if (nextBound0.ithObj != -1)
	    child[0]->getHigherOrderDofs(feSpace, nextBound0, dofs,
					 baseDofPtr, dofGeoIndex);
	  
	  if (nextBound1.ithObj != -1)
	    child[1]->getHigherOrderDofs(feSpace, nextBound1, dofs,
					 baseDofPtr, dofGeoIndex);
	} else { 

	  // === On leaf elements iterate over all DOFs and collect the  === 
	  // === right ones.                                             ===

	  ElementDofIterator elDofIter(feSpace, true);
	  elDofIter.reset(this);

	  bool next = true;
	  do {
	    bool addDof = false;

	    switch (elDofIter.getCurrentPos()) {
	    case 0:
	      // VERTEX is never a higher order DOF
	      addDof = false;
	      break;
	    case 1:
	      // On EDGE, first check whether this edge is part of the FACE and
	      // then check if this EDGE is not inside of excludedSubstructures.

	      {
		int localEdge = elDofIter.getCurrentElementPos();
		if (edgeOfFace[bound.ithObj][0] == localEdge ||
		    edgeOfFace[bound.ithObj][1] == localEdge ||
		    edgeOfFace[bound.ithObj][2] == localEdge) {
		  addDof = true;

		  for (ExcludeList::iterator it = bound.excludedSubstructures.begin(); 
		       it != bound.excludedSubstructures.end(); ++it)
		    if (it->first == EDGE && it->second == localEdge)
		      addDof = false;	
		} else {
		  addDof = false;
		}
	      }

	      break;
	    case 2:
	      // On FACE, check for the right one.
	      if (elDofIter.getCurrentElementPos() == bound.ithObj)
		addDof = true;
	      else
		addDof = false;
	      break;
		case 3:
		  // On CENTER, dof can not be part on internal boundary
		addDof = false;
		  break;
	    default:
	      ERROR_EXIT("Should not happen!\n");
	    }


	    if (addDof) {
	      if (baseDofPtr)
			dofs.push_back(elDofIter.getBaseDof());	
	      else
			dofs.push_back(elDofIter.getDofPtr());	
	      
	      if (dofGeoIndex != NULL)
			dofGeoIndex->push_back(elDofIter.getPosIndex());

	    }

	    next = (baseDofPtr ? elDofIter.nextStrict() : elDofIter.next());	    
	  } while (next);
	}
      }
      break;
    default:
      ERROR_EXIT("Should not happen!\n");
    }

  }


  void Tetrahedron::getSubBoundary(BoundaryObject bound, 
				   vector<BoundaryObject> &subBound) const
  {}


  void Tetrahedron::prepareNextBound(BoundaryObject &bound, int ithChild) const
  {
    FUNCNAME("Tetrahedron::prepareNextBound()");

    switch (bound.subObj) {
    case FACE:      
      for (ExcludeList::iterator it = bound.excludedSubstructures.begin(); 
   	   it != bound.excludedSubstructures.end(); ++it)
	if (it->first == EDGE && it->second != -1)
	  it->second = edgeOfChild[bound.elType][ithChild][it->second];	
      
      bound.ithObj = sideOfChild[bound.elType][ithChild][bound.ithObj];
      bound.elType = (bound.elType + 1) % 3;

      break;
    case EDGE:
      bound.ithObj = edgeOfChild[bound.elType][ithChild][bound.ithObj];
      bound.elType = (bound.elType + 1) % 3;
      break;
    default:
      ERROR_EXIT("Should not happen!\n");
    }
  }

} // end namespace AMDiS
