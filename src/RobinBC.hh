
namespace AMDiS {

  template <class JExpr, class AlphaExpr>
  void init(JExpr const& j, AlphaExpr const& alpha)
  {
    int dim = rowFeSpace->getMesh()->getDim();

    // create barycentric coords for each vertex of each side
    const Element *refElement = Global::getReferenceElement(dim);
    coords = new VectorOfFixVecs<DimVec<double> >*[dim+1];

    // for all element sides
    for (int i = 0; i < dim + 1; i++) {
      coords[i] =
	new VectorOfFixVecs<DimVec<double> >(dim, dim, DEFAULT_VALUE, 
					     DimVec<double>(dim, DEFAULT_VALUE, 0.0));
      // for each vertex of the side
      for (int k = 0; k < dim; k++) {
	int index = refElement->getVertexOfPosition(INDEX_OF_DIM(dim - 1, dim), i, k);
	(*coords[i])[k][index] = 1.0;
      }
    }
    
    Operator *jOp = new Operator(rowFeSpace);
    Operator *alphaOp = new Operator(rowFeSpace, colFeSpace);
    
    jOp->addZeroOrderTerm(j);
    alphaOp->addZeroOrderTerm(alpha);
    
    neumannOperators = new DimVec<SurfaceOperator*>(dim, NO_INIT);
    robinOperators = new DimVec<SurfaceOperator*>(dim, NO_INIT);
  
    for (int i = 0; i < dim + 1; i++) {
      (*neumannOperators)[i] = new SurfaceOperator(jOp, *coords[i]);    
      (*robinOperators)[i] = new SurfaceOperator(alphaOp, *coords[i]);  
    }

    delete jOp;
    delete alphaOp;
  }
  

  template <class JExpr>
  void init(JExpr const& j)
  {
    int dim = rowFeSpace->getMesh()->getDim();

    // create barycentric coords for each vertex of each side
    const Element *refElement = Global::getReferenceElement(dim);
    coords = new VectorOfFixVecs<DimVec<double> >*[dim+1];

    // for all element sides
    for (int i = 0; i < dim + 1; i++) {
      coords[i] =
	new VectorOfFixVecs<DimVec<double> >(dim, dim, DEFAULT_VALUE, 
					     DimVec<double>(dim, DEFAULT_VALUE, 0.0));
      // for each vertex of the side
      for (int k = 0; k < dim; k++) {
	int index = refElement->getVertexOfPosition(INDEX_OF_DIM(dim - 1, dim), i, k);
	(*coords[i])[k][index] = 1.0;
      }
    }
    
    Operator *jOp = new Operator(rowFeSpace);
    jOp->addZeroOrderTerm(j);
    neumannOperators = new DimVec<SurfaceOperator*>(dim, NO_INIT);
  
    for (int i = 0; i < dim + 1; i++)
      (*neumannOperators)[i] = new SurfaceOperator(jOp, *coords[i]);    

    delete jOp;
  }

} // end namespace AMDiS
