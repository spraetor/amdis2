namespace AMDiS
{
  template <class Expr>
  void ProblemStatSeq::addDirichletBC(BoundaryType type, int row, int col, Expr&& expr)
  {
    FUNCNAME("ProblemStat::addDirichletBC()");

    TEST_EXIT(row >= 0 && row < nComponents)("Wrong row number: %d\n", row);
    TEST_EXIT(col >= 0 && col < nComponents)("Wrong col number: %d\n", col);

    boundaryConditionSet = true;

    DirichletBC<Expr> *dirichletApply = 
      new DirichletBC<Expr>(type, std::forward<Expr>(expr), componentSpaces[row], componentSpaces[col], true);
    DirichletBC<Expr> *dirichletNotApply = 
      new DirichletBC<Expr>(type, std::forward<Expr>(expr), componentSpaces[row], componentSpaces[col], false);

    for (int i = 0; i < nComponents; i++)
      if (systemMatrix && (*systemMatrix)[row][i]) {
      	if (i == col)
      	  (*systemMatrix)[row][i]->getBoundaryManager()->addBoundaryCondition(dirichletApply);
      	else
      	  (*systemMatrix)[row][i]->getBoundaryManager()->addBoundaryCondition(dirichletNotApply);
      }	

    if (rhs)
      rhs->getDOFVector(row)->getBoundaryManager()->addBoundaryCondition(dirichletApply);
    if (solution)
      solution->getDOFVector(col)->getBoundaryManager()->addBoundaryCondition(dirichletApply);
  }

  
  template <class Expr>  
  void ProblemStatSeq::addNeumannBC(BoundaryType type, int row, int col, Expr&& expr)
  {
    boundaryConditionSet = true;

    NeumannBC *neumann = 
      new NeumannBC(type, std::forward<Expr>(expr), componentSpaces[row], componentSpaces[col]);

    if (rhs)
      rhs->getDOFVector(row)->getBoundaryManager()->addBoundaryCondition(neumann);
  }


  template <class ExprRhs, class ExprLhs>
  void ProblemStatSeq::addRobinBC(BoundaryType type, int row, int col, 
				  ExprRhs&& exprRhs, ExprLhs&& exprLhs)
  {
    boundaryConditionSet = true;

    RobinBC *robin = 
      new RobinBC(type, std::forward<ExprRhs>(exprRhs), std::forward<ExprLhs>(exprLhs), componentSpaces[row], componentSpaces[col]);

    if (systemMatrix && (*systemMatrix)[row][col])
      (*systemMatrix)[row][col]->getBoundaryManager()->addBoundaryCondition(robin);
    if (rhs)
      rhs->getDOFVector(row)->getBoundaryManager()->addBoundaryCondition(robin);
  }

} // end namespace AMDiS