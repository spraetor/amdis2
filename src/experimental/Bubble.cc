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


#include <stdio.h>
#include <algorithm>
#include <list>
#include "boost/lexical_cast.hpp"
#include "Mesh.h"
#include "Element.h"
#include "Bubble.h"
#include "DOFAdmin.h"
#include "RCNeighbourList.h"
#include "DOFVector.h"
#include "Traverse.h"
#include "Line.h"
#include "Triangle.h"
#include "Tetrahedron.h"
#include "Parametric.h"
#include "Debug.h"

namespace AMDiS {

  using namespace std;
  using boost::lexical_cast;

  std::vector<DimVec<double>* > Bubble::baryDimDegree; 	// lists with all basis-functions (important for later versions)
  DimVec<int>* Bubble::ndofDimDegree; 			// number of DOFs per element
  int Bubble::nBasFctsDimDegree;			// number of basis functions per element
  std::vector<BasFctType*> Bubble::phifunc; 		// basis functions
  std::vector<GrdBasFctType*> Bubble::grdPhifunc;	// first derivatives
  std::vector<D2BasFctType*> Bubble::D2Phifunc; 	// second derivatives

  Bubble* Bubble::Singleton = NULL;

  Bubble::Bubble(int dim, int degree)
	: BasisFunction(std::string("Bubble"), dim, degree)
  {
    // set name
    name += lexical_cast<std::string>(dim) + " " + lexical_cast<std::string>(degree);
	
    // set nDOF
    setNDOF();
	
    // set barycentric coordinates
   setBary();
		
   // set function pointer
   setFunctionPointer();		
  }


  Bubble::~Bubble()
  {
    for (int i = 0; i < static_cast<int>(bary->size()); i++)
    if ((*bary)[i]) {
      delete (*bary)[i];
      (*bary)[i] = NULL;
    }
  }


  // creates a new FE-Space of this instance
   Bubble* Bubble::getBubble(int dim, int degree) 
  {
    if (Singleton == NULL)  //if their is no instance
    {
      Singleton = new Bubble(dim, degree);
    }
    return Singleton;
  }


  void Bubble::clear()
  {    
    if (Singleton) {
      delete (Singleton);
    Singleton = NULL;
    }     
  }

  // set pointers to the used  basis functions and their first and
  // second derivatives.
  // set also the refinement and coarsening functions
  void Bubble::setFunctionPointer() 
  {
    if (static_cast<int>(phifunc.size()) == 0) {

      for(int i = 0; i < dim+1; i++){
	phifunc.push_back(new Phi(this, VERTEX, i, 0)); // Vertex
	grdPhifunc.push_back(new GrdPhi(this, VERTEX, i, 0));
	D2Phifunc.push_back(new D2Phi(this, VERTEX, i, 0));
      }
      phifunc.push_back(new Phi(this, CENTER, 0, 0)); // Bubble
      grdPhifunc.push_back(new GrdPhi(this, CENTER, 0, 0)); 
      D2Phifunc.push_back(new D2Phi(this, CENTER, 0, 0));

    }
    phi = &phifunc;
    grdPhi = &grdPhifunc;
    d2Phi = &D2Phifunc;

    // set refinement & coarse interpolation/restriction Function
    switch (dim) {
    case 1:
      refineInter_fct = refineInter2_1d;
      coarseRestr_fct = coarseRestr2_1d;
      coarseInter_fct = coarseInter2_1d;
      break;
    case 2:
      refineInter_fct = refineInter3_2d;
      coarseRestr_fct = coarseRestr3_2d;
      coarseInter_fct = coarseInter3_2d;
      break;
    case 3:
      refineInter_fct = refineInter4_3d;
      coarseRestr_fct = coarseRestr4_3d;
      coarseInter_fct = coarseInter4_3d;
      break;
    default:
      ERROR_EXIT("invalid dim!\n");
    }
  }


  
  void Bubble::setNDOF()
  {
    nBasFcts = dim + 2;
    nDOF = new DimVec<int>(dim, DEFAULT_VALUE, 0);
//     for (int i = 0; i < dim; i++)
      (*nDOF)[0] = 1; // CENTER
      (*nDOF)[1] = 1; // VERTEX
  } 

 
  DimVec<double> *Bubble::getCoords(int i) const
  {
    return (*bary)[i];
  }


  void Bubble::setVertices(int dim, int degree, 
	      GeoIndex position, int positionIndex, int nodeIndex, 
	      int** vertices)
  {
    FUNCNAME_DBG("Bubble::setVertices()");

    TEST_EXIT_DBG(*vertices == NULL)("vertices != NULL\n");

    int dimOfPosition = DIM_OF_INDEX(position, dim);

    *vertices = new int[dimOfPosition + 1];

    if (position == VERTEX) { // Vertex
      (*vertices)[0] = Global::getReferenceElement(dim)->getVertexOfPosition(position,
				positionIndex, 0);
    }
    else {  // Bubble
      for (int i = 0; i < dim+1; i++) {
	(*vertices)[i] = Global::getReferenceElement(dim)->getVertexOfPosition(position,
							  positionIndex, i);
      }
    }
  }


  Bubble::Phi::Phi(Bubble* owner, 
		GeoIndex position, 
		int positionIndex, 
		int nodeIndex)
		: vertices(NULL)
  {
    FUNCNAME("Bubble::Phi::Phi()");

    // get relevant vertices
    Bubble::setVertices(owner->getDim(), 
			owner->getDegree(), 
			position, 
			positionIndex, 
			nodeIndex, 
			&vertices);

    // set function pointer 
    switch (position) {
      case VERTEX:
	func = phi1v;
	break;
      case CENTER:
        switch (owner->getDim()) {
        case 1:
	  func = phi2c;
	  break;
        case 2:
          func = phi3c;
          break;
        case 3:
          func = phi4c;
          break;
	default:
	  ERROR_EXIT("invalid dim\n");
	  break;
	}
        break;
      default:
	ERROR_EXIT("invalid position\n");
    }
  }


  Bubble::Phi::~Phi()
  { 
    delete [] vertices; 
  }


  Bubble::GrdPhi::GrdPhi(Bubble* owner, 
		GeoIndex position, 
		int positionIndex, 
		int nodeIndex)
		: vertices(NULL)
  {
    // get relevant vertices
    Bubble::setVertices(owner->getDim(), 
			owner->getDegree(), 
			position, 
			positionIndex, 
			nodeIndex, 
			&vertices);

    // set function pointer 
    switch (position) {
      case VERTEX:
	func = grdPhi1v;
	break;
      case CENTER:
        switch (owner->getDim()) {
        case 1:
	  func = grdPhi2c;
	  break;
        case 2:
          func = grdPhi3c;
          break;
        case 3:
          func = grdPhi4c;
          break;
	default:
	  ERROR_EXIT("invalid dim\n");
	  break;
	}
        break;
      default:
	ERROR_EXIT("invalid position\n");
     }
  }

  Bubble::GrdPhi::~GrdPhi() 
  { 
    delete [] vertices; 
  }

  Bubble::D2Phi::D2Phi(Bubble* owner, 
			GeoIndex position, 
			int positionIndex, 
			int nodeIndex)
			: vertices(NULL)
  {
    // get relevant vertices
    Bubble::setVertices(owner->getDim(), 
		owner->getDegree(), 
		position, 
		positionIndex, 
		nodeIndex, 
		&vertices);

    // set function pointer 
    switch (position) {
      case VERTEX:
	func = D2Phi1v;
	break;
      case CENTER:
        switch (owner->getDim()) {
        case 1:
	  func = D2Phi2c;
	  break;
        case 2:
          func = D2Phi3c;
          break;
        case 3:
          func = D2Phi4c;
          break;
	default:
	  ERROR_EXIT("invalid dim\n");
	  break;
	}
        break;
      default:
	ERROR_EXIT("invalid position\n");
    }
  }

  Bubble::D2Phi::~D2Phi()
  { 
    delete [] vertices; 
  }

  // fills bary with the barycentric coordinates of the basis functions
  void Bubble::setBary()
  {
    bary = &baryDimDegree;
    if (static_cast<int>(bary->size()) == 0) { //nur ausführen, falls bary leer
      bary = new vector<DimVec<double>* >;

      for (int i = 0; i < dim+1; i++) {
	DimVec<double>* unit = new DimVec<double>(dim, DEFAULT_VALUE, 0.0);
	(*unit)[i] = 1.0;
	bary->push_back(unit); //coordinates of the Lagrange functions
      }
      bary->push_back(new DimVec<double>(dim, DEFAULT_VALUE,1.0/(dim + 1.0))); //coordinates of the bubble function
    }	
  }

  // position = {vertex, edge, face, center}
  // positionIndex i'th DOF at position
  int* Bubble::orderOfPositionIndices(const Element* el, 
				GeoIndex position, 
				int positionIndex) const
  {
    // returns 0 in all valid cases
    static int sortedVertex = 0;
    static int sortedCenter = 0;

    if (position == VERTEX)
      return &sortedVertex; 

    if (position == CENTER)
      return &sortedCenter;

    ERROR_EXIT("should not be reached\n");
    return NULL;
  }


  // assign each DOF to the corresponding boundary type, 
  // following a strict order (vertices DOFs, edge-DOFs, faces-DOFs, center-DOFs)
  // with the help of the variable offset
  void Bubble::getBound(const ElInfo* elInfo, BoundaryType* bound) const
  {
    elInfo->testFlag(Mesh::FILL_BOUND);

    // boundaries
    int index = 0;
    int offset = 0;
    BoundaryType boundaryType;
    for (int i = dim - 1; i > 0; i--) //EDGES
      offset += Global::getGeo(INDEX_OF_DIM(i, dim), dim);

    for (int i = 0; i < dim; i++) {
      int jto = offset + Global::getGeo(INDEX_OF_DIM(i, dim), dim); 
      for (int j = offset; j < jto; j++) {
	boundaryType = elInfo->getBoundary(j);
	int kto = (*nDOF)[INDEX_OF_DIM(i, dim)]; 
	for (int k = 0; k < kto; k++)
	  bound[index++] = boundaryType;	
      }
      offset -= Global::getGeo(INDEX_OF_DIM(i + 1, dim), dim);
    }

    // interior nodes in the center
    for (int i = 0; i < (*nDOF)[CENTER]; i++)
      bound[index++] = INTERIOR;    

    TEST_EXIT_DBG(index == nBasFcts)("found not enough boundarys. index=%d, nBasFcts=%d\n",index,nBasFcts);
  }


  // no = number of the basis fuctions where the interpolation coefficient is calculated
  // b_no = pointer to the indices of the basis-functions
  // f = function to be interpolated
  // rvec = vector with the coefficient values
  void Bubble::interpol(const ElInfo *elInfo, 
			int no, 
			const int *b_no, 
			AbstractFunction<double, WorldVector<double> > *f, 
			mtl::dense_vector<double> &rvec) const
  {
    FUNCNAME_DBG("Bubble::interpol()");
			
    WorldVector<double> x;

    elInfo->testFlag(Mesh::FILL_COORDS);

    double sum = 0;
    int BubbleIndex=-1;
				
    if (b_no) {
      TEST_EXIT_DBG(no >= 0 && no < getNumber())("Something is wrong!\n");
      for (int i = 0; i < no; i++) {
	if (b_no[i] < Global::getGeo(VERTEX, dim)) {
	  // interpolation in vertices:
	  rvec[i] = (*f)(elInfo->getCoord(b_no[i]));
	  sum += rvec[i];
	} else {
	  // interpolation at center:
	  elInfo->coordToWorld(*(*bary)[b_no[i]], x);
	  rvec[i] = (*f)(x);
	  BubbleIndex = i;
	}
      }
    }
    else {
      int vertices = Global::getGeo(VERTEX, dim);
      BubbleIndex = vertices;
      for (int i = 0; i < vertices; i++) {
	// interpolation in vertices:						
	rvec[i] = (*f)(elInfo->getCoord(i));
	sum += rvec[i];
      }
      // interpolation at center:
      elInfo->coordToWorld(*(*bary)[vertices], x);
      rvec[vertices] = (*f)(x);
    }  

    // correction of the center-interpolation
    if (BubbleIndex >= 0)
      rvec[BubbleIndex] -= sum/(dim + 1.0);
  }

  // no = number of the basis fuctions where the interpolation coefficient is calculated
  // b_no = pointer to the indices of the basis-functions
  // f = function to be interpolated
  // rvec = vector with the coefficient values
  void Bubble::interpol(const ElInfo *elInfo, 
				int no, 
				const int *b_no,
				AbstractFunction<WorldVector<double>, WorldVector<double> > *f, 
				mtl::dense_vector<WorldVector<double> > &rvec) const   //DONE
  {
    FUNCNAME_DBG("*Bubble::interpol()");

    WorldVector<double> x;

    elInfo->testFlag(Mesh::FILL_COORDS);

    WorldVector<double> sum(DEFAULT_VALUE, 0.0);
    int BubbleIndex=-1;
    int vertices = Global::getGeo(VERTEX, dim);

    if (b_no) {
      TEST_EXIT_DBG(no >= 0 && no < getNumber())("Something is wrong!\n");

      for (int i = 0; i < no; i++) {
	if (b_no[i] < Global::getGeo(VERTEX, dim)) {
	  // interpolation in vertices:
	  rvec[i] = (*f)(elInfo->getCoord(b_no[i]));
	  sum += rvec[i];
	} 
	else {
	  // interpolation at center:
	  elInfo->coordToWorld(*(*bary)[b_no[i]], x);
	  rvec[i] = (*f)(x);
	  BubbleIndex = i;
	}
      }
    } 
    else {
      BubbleIndex = vertices;
      for (int i = 0; i < vertices; i++) {
	// interpolation in vertices:
	rvec[i] = (*f)(elInfo->getCoord(i));
	sum += rvec[i];
      }
      // interpolation at center
      elInfo->coordToWorld(*(*bary)[vertices], x);
      rvec[vertices] = (*f)(x);
    } 
    // correction of the center-interpolation
    if (BubbleIndex >= 0)
      rvec[BubbleIndex] -= sum/(dim + 1.0) ;
  }

  // collects and order the DOFs per element
  // at first the DOF at the vertices and then the DOF in the center
  void Bubble::getLocalIndices(const Element *el, 
				const DOFAdmin *admin,
				std::vector<DegreeOfFreedom>& dofs) const
  {
    if (static_cast<int>(dofs.size()) < nBasFcts)
      dofs.resize(nBasFcts);

    const DegreeOfFreedom **dof = el->getDof();
    GeoIndex posIndex;
    int n0, node0, j=0;

    // DOFs at vertices
    posIndex = VERTEX;
    n0 = admin->getNumberOfPreDofs(posIndex);
    node0 = admin->getMesh()->getNode(posIndex);
    int num = Global::getGeo(posIndex, dim);
    for (int i = 0; i < num; node0++, i++) {
      const int *indi = orderOfPositionIndices(el, posIndex, i);
      dofs[j++] = dof[node0][n0 + indi[0]];
    }

    // DOFs at center
    posIndex = CENTER;
    n0 = admin->getNumberOfPreDofs(posIndex);
    node0 = admin->getMesh()->getNode(posIndex);
    const int *indi = orderOfPositionIndices(el, posIndex, 0);
    dofs[j] = dof[node0][n0 + indi[0]];
  }


  // returns a Vector filled with Pointer to the corresponding DOFs,
  // sorted in following way: vertices, edges, faces, center
  void Bubble::getLocalDofPtrVec(const Element *el, 
				const DOFAdmin *admin,
				std::vector<const DegreeOfFreedom*>& dofs) const
  {
    if (static_cast<int>(dofs.size()) < nBasFcts)
      dofs.resize(nBasFcts);

    const DegreeOfFreedom **dof = el->getDof();
    GeoIndex posIndex;
    int n0, node0, j=0;

    // DOFs at vertices
    posIndex = VERTEX;
    n0 = admin->getNumberOfPreDofs(posIndex);
    node0 = admin->getMesh()->getNode(posIndex);
    int num = Global::getGeo(posIndex, dim);
    for (int i = 0; i < num; node0++, i++)
    {
      const int *indi = orderOfPositionIndices(el, posIndex, i);
      dofs[j++] = &(dof[node0][n0 + indi[0]]);
     }

    // DOFs at center
    posIndex = CENTER;
    n0 = admin->getNumberOfPreDofs(posIndex);
    node0 = admin->getMesh()->getNode(posIndex);
    const int *indi = orderOfPositionIndices(el, posIndex, 0);
    dofs[j] = &(dof[node0][n0 + indi[0]]);
  }


  void Bubble::l2ScpFctBas(Quadrature *q,
			AbstractFunction<WorldVector<double>, WorldVector<double> >* f,
			DOFVector<WorldVector<double> >* fh)
  {
    ERROR_EXIT("not yet\n");
  }


  // calculates the L2-scalarproduct of given function f with all basis-functions by quadrature
  void Bubble::l2ScpFctBas(Quadrature *quad,
			AbstractFunction<double, WorldVector<double> >* f,
			DOFVector<double>* fh) 
  {
    FUNCNAME_DBG("Bubble::l2ScpFctBas()");

    TEST_EXIT_DBG(fh)("no DOF_REAL_VEC fh\n");
    TEST_EXIT_DBG(fh->getFeSpace())
	("no fe_space in DOF_REAL_VEC %s\n", fh->getName().c_str());
    TEST_EXIT_DBG(fh->getFeSpace()->getBasisFcts() == this)
	("wrong basis fcts for fh\n");
    TEST_EXIT_DBG(!fh->getFeSpace()->getMesh()->getParametric())
	("Not yet implemented!");

    if (!quad)
      quad = Quadrature::provideQuadrature(dim, 2 * degree - 2);

    const FastQuadrature *quad_fast = FastQuadrature::provideFastQuadrature(this, *quad, INIT_PHI);
    vector<double> wdetf_qp(quad->getNumPoints());
    int nPoints = quad->getNumPoints();
    DOFAdmin *admin = fh->getFeSpace()->getAdmin();
    WorldVector<double> x;
    vector<DegreeOfFreedom> dof;

    TraverseStack stack;
    ElInfo *elInfo = stack.traverseFirst(fh->getFeSpace()->getMesh(), -1, 
					Mesh::CALL_LEAF_EL | Mesh::FILL_COORDS);
    while (elInfo) {
      getLocalIndices(elInfo->getElement(), admin, dof);
      double det = elInfo->getDet();

      for (int iq = 0; iq < nPoints; iq++) {
	elInfo->coordToWorld(quad->getLambda(iq), x);
	wdetf_qp[iq] = quad->getWeight(iq) * det * ((*f)(x));
      }

      for (int j = 0; j < nBasFcts; j++) {
	double val = 0.0;
	for (int iq = 0; iq < nPoints; iq++)
	  val += quad_fast->getPhi(iq, j) * wdetf_qp[iq];

	(*fh)[dof[j]] += val;
      }

      elInfo = stack.traverseNext(elInfo);
    }
  }

  void Bubble::refineInter2_1d(DOFIndexed<double> *drv, 
				 RCNeighbourList* list, 
				 int n, 
				 BasisFunction* basFct)
  {
    if (n < 1) 
      return;
    
    Element *el = list->getElement(0);
    const DOFAdmin *admin = drv->getFeSpace()->getAdmin();
    vector<DegreeOfFreedom> pdof;
    basFct->getLocalIndices(el, admin, pdof);
  
    int node = drv->getFeSpace()->getMesh()->getNode(VERTEX);        
    int n0 = admin->getNumberOfPreDofs(VERTEX);
        
    // vertex-dof 1 of child[0]
    DegreeOfFreedom cdof1 = el->getChild(0)->getDof(node + 1, n0);
    (*drv)[cdof1] = 0.5*((*drv)[pdof[0]]+(*drv)[pdof[1]]) + (*drv)[pdof[2]];

    // vertex-dof 0 of child[1]
    DegreeOfFreedom cdof0 = el->getChild(1)->getDof(node, n0);
    (*drv)[cdof0] = 0.5*((*drv)[pdof[0]]+(*drv)[pdof[1]]) + (*drv)[pdof[2]];
    
    node = drv->getFeSpace()->getMesh()->getNode(CENTER);		       
    n0 = admin->getNumberOfPreDofs(CENTER);					

    // barycenter of child[0]	
    DegreeOfFreedom cdof2 = el->getChild(0)->getDof(node, n0); 
    (*drv)[cdof2] = 0.25*(*drv)[pdof[2]];

    // barycenter of child[1]
    cdof2 = el->getChild(1)->getDof(node, n0); 
    (*drv)[cdof2] = 0.25*(*drv)[pdof[2]];
  }

  // drv = global DOF-values
  // list = informations about element and the neighbours
  // n = number of neighbours
  void Bubble::refineInter3_2d(DOFIndexed<double> *drv, 
				RCNeighbourList* list, 
				int n, BasisFunction* basFct)
  {
    if (n < 1) 
      return;
			
    const DOFAdmin *admin = drv->getFeSpace()->getAdmin();
    Element *el = list->getElement(0);
    vector<DegreeOfFreedom> pdof;
    basFct->getLocalIndices(el, admin, pdof);

    int node = drv->getFeSpace()->getMesh()->getNode(VERTEX);		
    int n0 = admin->getNumberOfPreDofs(VERTEX);	

    // newest vertex of child[0] and child[1], newest vertex is node+2 
    DegreeOfFreedom cdof2 = el->getChild(0)->getDof(node + 2, n0);
    (*drv)[cdof2] = 0.5*((*drv)[pdof[0]]+(*drv)[pdof[1]]);

    node = drv->getFeSpace()->getMesh()->getNode(CENTER);		       
    n0 = admin->getNumberOfPreDofs(CENTER);					

    // barycenter of child[0]	
    DegreeOfFreedom cdof3 = el->getChild(0)->getDof(node, n0); 
    (*drv)[cdof3] = 0.75*(*drv)[pdof[3]];

    // barycenter of child[1]
    cdof3 = el->getChild(1)->getDof(node, n0); 
    (*drv)[cdof3] = 0.75*(*drv)[pdof[3]];

    if (n > 1) {
      // adjust the values at the barycenters of the neighbour's children
      el = list->getElement(1);
      basFct->getLocalIndices(el, admin, pdof);

      // barycenter of child[0]
      cdof3 = el->getChild(0)->getDof(node, n0); 
      (*drv)[cdof3] = 0.75*(*drv)[pdof[3]];

      // barycenter of child[1]								                 
      cdof3 = el->getChild(1)->getDof(node, n0); 
      (*drv)[cdof3] = 0.75*(*drv)[pdof[3]];

    }
  }

  void Bubble::refineInter4_3d(DOFIndexed<double> *drv, 
				 RCNeighbourList* list, 
				 int n, 
				 BasisFunction* basFct)
  {
    if (n < 1) 
      return;
    
    Element *el = list->getElement(0);
    const DOFAdmin *admin = drv->getFeSpace()->getAdmin();
    vector<DegreeOfFreedom> pdof;
    basFct->getLocalIndices(el, admin, pdof);
  
    int node = drv->getFeSpace()->getMesh()->getNode(VERTEX);        
    int n0 = admin->getNumberOfPreDofs(VERTEX);
    
    // newest vertex of child[0] and child[1], newest vertex is node+3 
    DegreeOfFreedom cdof3 = el->getChild(0)->getDof(node + 3, n0);
    (*drv)[cdof3] = 0.5*((*drv)[pdof[0]]+(*drv)[pdof[1]]);

    node = drv->getFeSpace()->getMesh()->getNode(CENTER);		       
    n0 = admin->getNumberOfPreDofs(CENTER);					

    // barycenter of child[0]	
    DegreeOfFreedom cdof4 = el->getChild(0)->getDof(node, n0); 
    (*drv)[cdof4] = 0.75*(*drv)[pdof[4]];

    // barycenter of child[1]
    cdof4 = el->getChild(1)->getDof(node, n0); 
    (*drv)[cdof4] = 0.75*(*drv)[pdof[4]];

    for (int i = 1; i < n; i++) {
      el = list->getElement(i);
      TEST_EXIT_DBG(el)("Should not happen!\n");

      basFct->getLocalIndices(el, admin, pdof);
      
      // barycenter of child[0]
      cdof4 = el->getChild(0)->getDof(node, n0); 
      (*drv)[cdof4] = 0.75*(*drv)[pdof[4]];

      // barycenter of child[1]								                 
      cdof4 = el->getChild(1)->getDof(node, n0); 
      (*drv)[cdof4] = 0.75*(*drv)[pdof[4]];
    }
  }

  void Bubble::coarseRestr2_1d(DOFIndexed<double> *drv, 
			       RCNeighbourList *list, 
			       int n, BasisFunction* basFct)
  {
    FUNCNAME("Bubble::coarseRestr2_1d()");
				
    TEST_EXIT(drv->getFeSpace())("No fe_space in dof_real_vec!\n");
    TEST_EXIT(drv->getFeSpace()->getBasisFcts())("No basis functions in fe_space!\n");

    ERROR_EXIT("not yet implemented!\n");
    if (n < 1)
      return;
#if 0
    const Element *el = list->getElement(0);
    const DOFAdmin *admin = drv->getFeSpace()->getAdmin();
#endif
    // ... //
  }

  // drv = global DOF-values
  // list = informations about element and the neighbours
  // n = number of neighbours
  // Calculation of the values can be traced in the documentation
  // NOTE: maybe the implementation is wrong!
  void Bubble::coarseRestr3_2d(DOFIndexed<double> *drv, 
			       RCNeighbourList *list, 
			       int n, BasisFunction* basFct)
  {
    FUNCNAME("Bubble::coarseRestr3_2d()");
				
    TEST_EXIT(drv->getFeSpace())("No fe_space in dof_real_vec!\n");
    TEST_EXIT(drv->getFeSpace()->getBasisFcts())("No basis functions in fe_space!\n");

    if (n < 1)
      return;

    const Element *el = list->getElement(0);
    const DOFAdmin *admin = drv->getFeSpace()->getAdmin();
    vector<DegreeOfFreedom> pdof(4), cdof0(4), cdof1(4);
    vector<double> cdof(9);
    basFct->getLocalIndices(el, admin, pdof);
    basFct->getLocalIndices(el->getChild(0), admin, cdof0);
    basFct->getLocalIndices(el->getChild(1), admin, cdof1);
				
    // Create a Vektor with the values of all Child-DOFs
    cdof[0] = (*drv)[cdof0[0]];
    cdof[1] = (*drv)[cdof0[1]];
    cdof[2] = (*drv)[cdof0[2]];
    cdof[3] = (*drv)[cdof0[3]];
    cdof[4] = (*drv)[cdof1[0]];
    cdof[5] = (*drv)[cdof1[3]];

    if (n==1){
      //no neighbours:

      // Calculate Parent DOF values via inv(A)*B
      (*drv)[pdof[0]] = 61.0/72.0*cdof[1] + 11.0/36.0*cdof[2]
		      + 87.0/160.0*cdof[3] - 11.0/72.0*cdof[4]
		      - 57.0/160.0*cdof[5];
      (*drv)[pdof[1]] = -11.0/72.0*cdof[1] + 11.0/36.0*cdof[2]
		      - 57.0/160.0*cdof[3] + 61.0/72.0*cdof[4]
		      + 87.0/160.0*cdof[5];
      (*drv)[pdof[2]] = cdof[0] + 7.0/72.0*cdof[1] - 7.0/36.0*cdof[2]
		      + 3.0/32.0*cdof[3] + 7.0/72.0*cdof[4]
		      + 3.0/32.0*cdof[5];
      (*drv)[pdof[3]] = -35.0/162.0*cdof[1] + 35.0/81.0*cdof[2]
		      + 7.0/24.0*cdof[3] - 35.0/162.0*cdof[4]
		      + 7.0/24.0*cdof[5];

    } 
    else {
      vector<DegreeOfFreedom> pdof1(4);
      // there is a neighbour:
      el = list->getElement(1);
      basFct->getLocalIndices(el, admin, pdof1);
      basFct->getLocalIndices(el->getChild(0), admin, cdof0);
      basFct->getLocalIndices(el->getChild(1), admin, cdof1);

      // Create a Vektor with the values of all Child-DOFs
      // cdof[0], cdof[1] and cdof[4] are the values before coarsening of the first element
      cdof[6] = (*drv)[cdof0[0]];
      cdof[7] = (*drv)[cdof0[3]];
      cdof[8] = (*drv)[cdof1[3]];

      // Calculate Parent DOF values via modified inv(A)*B
      (*drv)[pdof[0]] = 61.0/72.0 * cdof[1] + 11.0/36.0 * cdof[2]
		      + 87.0/320.0 * cdof[3] -11.0/72.0 * cdof[4] -57.0/320.0 * cdof[5] 
		      - 57.0/320.0 * cdof[7] + 87.0/320.0 * cdof[8];
      (*drv)[pdof[1]] = -11.0/72.0 * cdof[1] + 11.0/36.0 * cdof[2]
		      - 57.0/320.0 * cdof[3] + 61.0/72.0  * cdof[4] + 87.0/320.0 * cdof[5] 
		      + 87.0/320.0 * cdof[7] - 57.0/320.0 * cdof[8];  
      (*drv)[pdof[2]] =  cdof[0] + 7.0/72.0 * cdof[1] - 7.0/36.0 * cdof[2] 
		      + 51.0/512.0 * cdof[3] + 7.0/72.0 * cdof[4] + 51.0/512.0  * cdof[5] 
		      - 3.0/512.0 * cdof[7] - 3.0/512.0 * cdof[8];  
      (*drv)[pdof[3]] = -35.0/162.0 * cdof[1] + 35.0/81.0 * cdof[2]
		      + 259.0/768.0 * cdof[3] - 35.0/162.0  * cdof[4] + 259.0/768.0 * cdof[5] 
		      - 35.0/768.0 * cdof[7] - 35.0/768.0 * cdof[8];  
      (*drv)[pdof1[2]] = 7.0/72.0 * cdof[1] - 7.0/36.0 * cdof[2] 
			- 3.0/512.0  * cdof[3] + 7.0/72.0 * cdof[4] - 3.0/512.0 * cdof[5] 
			+  cdof[6] + 51.0/512.0  * cdof[7] + 51.0/512.0  * cdof[8];  
      (*drv)[pdof1[3]] = -35.0/162.0 * cdof[1] + 35.0/81.0* cdof[2]
			- 35.0/768.0 * cdof[3] - 35.0/162.0 * cdof[4] - 35.0/768.0  * cdof[5] 
			+ 259.0/768.0  * cdof[7] + 259.0/768.0  * cdof[8];    
    }
  }

  void Bubble::coarseRestr4_3d(DOFIndexed<double> *drv, 
			       RCNeighbourList *list, 
			       int n, BasisFunction* basFct)
  {
    FUNCNAME("Bubble::coarseRestr4_3d()");
				
    TEST_EXIT(drv->getFeSpace())("No fe_space in dof_real_vec!\n");
    TEST_EXIT(drv->getFeSpace()->getBasisFcts())("No basis functions in fe_space!\n");

    ERROR_EXIT("not yet implemented!\n");
    if (n < 1)
      return;
#if 0
    const Element *el = list->getElement(0);
    const DOFAdmin *admin = drv->getFeSpace()->getAdmin();
#endif
    // ... //
  }

  void Bubble::coarseInter2_1d(DOFIndexed<double> *drv, 
			       RCNeighbourList* list, 
			       int n, BasisFunction* basFct)
  {
    FUNCNAME("Bubble::coarseInter2_1d()");

    TEST_EXIT(drv->getFeSpace())("No fe_space in dof_real_vec!\n");
    TEST_EXIT(drv->getFeSpace()->getBasisFcts())("No basis functions in fe_space!\n");

    if (n < 1) 
      return;

    Element *el = list->getElement(0);
    const DOFAdmin *admin = drv->getFeSpace()->getAdmin();

    int node = drv->getFeSpace()->getMesh()->getNode(VERTEX);  
    int n0 = admin->getNumberOfPreDofs(VERTEX);					

    DegreeOfFreedom cdof0_0 = el->getChild(0)->getDof(node, n0);
    DegreeOfFreedom cdof1_0 = el->getChild(0)->getDof(node + 1, n0);
    DegreeOfFreedom cdof1_1 = el->getChild(1)->getDof(node + 1, n0);
    vector<DegreeOfFreedom> pdof;
    basFct->getLocalIndices(el, admin, pdof);

    (*drv)[pdof[2]] = (*drv)[cdof1_0] - 1.0/2.0*((*drv)[cdof0_0]+(*drv)[cdof1_1]);
  }

  // drv = global DOF-values
  // list = informations about element and the neighbours
  // n = number of neighbours
  void Bubble::coarseInter3_2d(DOFIndexed<double> *drv, 
			       RCNeighbourList* list, 
			       int n, BasisFunction* basFct)
  {
    FUNCNAME("Bubble::coarseInter3_2d()");

    TEST_EXIT(drv->getFeSpace())("No fe_space in dof_real_vec!\n");
    TEST_EXIT(drv->getFeSpace()->getBasisFcts())("No basis functions in fe_space!\n");

    if (n < 1) 
      return;

    Element *el = list->getElement(0);
    const DOFAdmin *admin = drv->getFeSpace()->getAdmin();

    int node = drv->getFeSpace()->getMesh()->getNode(VERTEX);  
    int n0 = admin->getNumberOfPreDofs(VERTEX);					

    DegreeOfFreedom cdof2 = el->getChild(0)->getDof(node + 2, n0);
    DegreeOfFreedom cdof1 = el->getChild(0)->getDof(node + 1, n0);
    DegreeOfFreedom cdof0 = el->getChild(1)->getDof(node, n0);
    vector<DegreeOfFreedom> pdof;
    basFct->getLocalIndices(el, admin, pdof);

    (*drv)[pdof[3]] = 2.0/3.0*(*drv)[cdof2] - 1.0/3.0*((*drv)[cdof1]+(*drv)[cdof0]);
				
    if (n > 1) {
      Element *el = list->getElement(1);
      cdof2 = el->getChild(0)->getDof(node + 2, n0);
      cdof1 = el->getChild(0)->getDof(node + 1, n0);
      cdof0 = el->getChild(1)->getDof(node, n0);
      basFct->getLocalIndices(el, admin, pdof);
      (*drv)[pdof[3]] = 2.0/3.0*(*drv)[cdof2] - 1.0/3.0*((*drv)[cdof1]+(*drv)[cdof0]);
    }				
  }

  void Bubble::coarseInter4_3d(DOFIndexed<double> *drv, 
			       RCNeighbourList* list, 
			       int n, BasisFunction* basFct)
  {
    FUNCNAME("Bubble::coarseInter4_3d()");

    TEST_EXIT(drv->getFeSpace())("No fe_space in dof_real_vec!\n");
    TEST_EXIT(drv->getFeSpace()->getBasisFcts())("No basis functions in fe_space!\n");

    if (n < 1) 
      return;

    Element *el = list->getElement(0);
    const DOFAdmin *admin = drv->getFeSpace()->getAdmin();

    int node = drv->getFeSpace()->getMesh()->getNode(VERTEX);  
    int n0 = admin->getNumberOfPreDofs(VERTEX);					
    
    DegreeOfFreedom cdof3_0 = el->getChild(0)->getDof(node + 3, n0);
    DegreeOfFreedom cdof0_0 = el->getChild(0)->getDof(node, n0);
    DegreeOfFreedom cdof0_1 = el->getChild(1)->getDof(node, n0);
    vector<DegreeOfFreedom> pdof;
    basFct->getLocalIndices(el, admin, pdof);

    (*drv)[pdof[4]] = 0.5*(*drv)[cdof3_0] - 0.25*((*drv)[cdof0_0]+(*drv)[cdof0_1]);
				
    for (int i = 1; i < n; i++) {
      Element *el = list->getElement(i);
      basFct->getLocalIndices(el, admin, pdof);
      
      cdof3_0 = el->getChild(0)->getDof(node + 3, n0);
      cdof0_0 = el->getChild(0)->getDof(node, n0);
      cdof0_1 = el->getChild(1)->getDof(node, n0);
      (*drv)[pdof[4]] = 0.5*(*drv)[cdof3_0] - 0.25*((*drv)[cdof0_0]+(*drv)[cdof0_1]);
    }				
  }
}
