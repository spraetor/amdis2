#include "DirichletBC.hpp"
#include "RobinBC.hpp"

namespace AMDiS
{
  template <class Term>
  void ProblemStatSeq::addDirichletBC(BoundaryType type, int row, int col, Term&& term)
  {
    FUNCNAME("ProblemStat::addDirichletBC()");

    TEST_EXIT(row >= 0 && row < nComponents)("Wrong row number: %d\n", row);
    TEST_EXIT(col >= 0 && col < nComponents)("Wrong col number: %d\n", col);

    boundaryConditionSet = true;

    DirichletBC<Term>* dirichletApply =
      new DirichletBC<Term>(type, std::forward<Term>(term), componentSpaces[row], componentSpaces[col], true);
    DirichletBC<Term>* dirichletNotApply =
      new DirichletBC<Term>(type, std::forward<Term>(term), componentSpaces[row], componentSpaces[col], false);

    for (int i = 0; i < nComponents; i++)
      if (systemMatrix && (*systemMatrix)[row][i])
      {
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


  template <class Term>
  void ProblemStatSeq::addNeumannBC(BoundaryType type, int row, int col, Term&& term)
  {
    boundaryConditionSet = true;

    NeumannBC* neumann =
      new NeumannBC(type, std::forward<Term>(term), componentSpaces[row], componentSpaces[col]);

    if (rhs)
      rhs->getDOFVector(row)->getBoundaryManager()->addBoundaryCondition(neumann);
  }


  template <class TermRhs, class TermLhs>
  void ProblemStatSeq::addRobinBC(BoundaryType type, int row, int col,
                                  TermRhs&& termRhs, TermLhs&& termLhs)
  {
    boundaryConditionSet = true;

    RobinBC* robin =
      new RobinBC(type, std::forward<TermRhs>(termRhs), std::forward<TermLhs>(termLhs), 
			componentSpaces[row], componentSpaces[col]);

    if (systemMatrix && (*systemMatrix)[row][col])
      (*systemMatrix)[row][col]->getBoundaryManager()->addBoundaryCondition(robin);
    if (rhs)
      rhs->getDOFVector(row)->getBoundaryManager()->addBoundaryCondition(robin);
  }

} // end namespace AMDiS