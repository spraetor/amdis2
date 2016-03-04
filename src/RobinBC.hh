namespace AMDiS
{
  template <class JExpr, class AlphaExpr>
  void init(JExpr&& j_expr, AlphaExpr&& a_expr)
  {
    int dim = rowFeSpace->getMesh()->getDim();

    // create barycentric coords for each vertex of each side
    const Element* refElement = Global::getReferenceElement(dim);
    coords = new VectorOfFixVecs<DimVec<double>>* [dim+1];

    // for all element sides
    for (int i = 0; i < dim + 1; i++)
    {
      coords[i] =
        new VectorOfFixVecs<DimVec<double>>(dim, dim, DEFAULT_VALUE,
                                            DimVec<double>(dim, 0.0));
      // for each vertex of the side
      for (int k = 0; k < dim; k++)
      {
        int index = refElement->getVertexOfPosition(INDEX_OF_DIM(dim - 1, dim), i, k);
        (*coords[i])[k][index] = 1.0;
      }
    }

    Operator jOp(rowFeSpace);
    Operator alphaOp(rowFeSpace, colFeSpace);

    jOp.addZeroOrderTerm(std::forward<JExpr>(j_expr));
    alphaOp.addZeroOrderTerm(std::forward<AlphaExpr>(a_expr));

    neumannOperators = new DimVec<SurfaceOperator*>(dim);
    robinOperators = new DimVec<SurfaceOperator*>(dim);

    for (int i = 0; i < dim + 1; i++)
    {
      (*neumannOperators)[i] = new SurfaceOperator(jOp, *coords[i]);
      (*robinOperators)[i] = new SurfaceOperator(alphaOp, *coords[i]);
    }
  }


  template <class JExpr>
  void init(JExpr&& j_expr)
  {
    int dim = rowFeSpace->getMesh()->getDim();

    // create barycentric coords for each vertex of each side
    const Element* refElement = Global::getReferenceElement(dim);
    coords = new VectorOfFixVecs<DimVec<double>>* [dim+1];

    // for all element sides
    for (int i = 0; i < dim + 1; i++)
    {
      coords[i] =
        new VectorOfFixVecs<DimVec<double>>(dim, dim, DEFAULT_VALUE,
                                            DimVec<double>(dim, 0.0));
      // for each vertex of the side
      for (int k = 0; k < dim; k++)
      {
        int index = refElement->getVertexOfPosition(INDEX_OF_DIM(dim - 1, dim), i, k);
        (*coords[i])[k][index] = 1.0;
      }
    }

    Operator jOp(rowFeSpace);
    jOp->addZeroOrderTerm(std::forward<JExpr>(j_expr));
    neumannOperators = new DimVec<SurfaceOperator*>(dim);

    for (int i = 0; i < dim + 1; i++)
      (*neumannOperators)[i] = new SurfaceOperator(jOp, *coords[i]);
  }

} // end namespace AMDiS
