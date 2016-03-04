#include "time/RosenbrockStationary.hpp"

// AMDiS headers
#include "io/VtkWriter.hpp"
#include "ProblemStat.hpp"
#include "SystemVector.hpp"
// #include <Debug.hpp>

#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
#include "parallel/MeshDistributor.hpp"
#endif

namespace AMDiS
{
  void RosenbrockStationary::init()
  {
    stageSolution = new SystemVector(*solution);
    unVec = new SystemVector(*solution);
    timeRhsVec = new SystemVector(*solution);
    newUn = new SystemVector(*solution);
    tmp = new SystemVector(*solution);
    lowSol = new SystemVector(*solution);

    stageSolution->set(0.0);
    unVec->set(0.0);

    stageSolutions.resize(rm->getStages());
    for (int i = 0; i < rm->getStages(); i++)
    {
      stageSolutions[i] = new SystemVector(*solution);
      stageSolutions[i]->set(0.0);
    }
  }

  Flag RosenbrockStationary::oneIteration(AdaptInfo& adaptInfo, Flag toDo)
  {
    Flag flag = 0, markFlag = 0;


    if (toDo.isSet(MARK))
      markFlag = problem->markElements(adaptInfo);
    else
      markFlag = 3;

    // refine
    if (toDo.isSet(ADAPT) && markFlag.isSet(MESH_REFINED))
      flag = problem->refineMesh(adaptInfo);

    // coarsen
    if (toDo.isSet(ADAPT) && markFlag.isSet(MESH_COARSENED))
      flag |= problem->coarsenMesh(adaptInfo);

    if (toDo.isSet(BUILD) || toDo.isSet(SOLVE))
    {
      flag = stageIteration(adaptInfo, toDo, true, true);
      estimateTimeError(adaptInfo);
    }

    if (toDo.isSet(ESTIMATE))
      problem->estimate(adaptInfo);

    return flag;
  }


  Flag RosenbrockStationary::stageIteration(AdaptInfo& adaptInfo, Flag flag,
      bool asmMatrix, bool asmVector)
  {
    FUNCNAME("RosenbrockStationary::stageIteration()");

    TEST_EXIT(tauPtr)("No tau pointer defined in stationary problem!\n");

    if (first)
    {
      first = false;
      *unVec = *solution;
    }

    *newUn = *unVec;
    *lowSol = *unVec;
    for (int i = 0; i < rm->getStages(); i++)
    {
      stageTime = oldTime + rm->getAlphaI(i) * (*tauPtr);
      tauGammaI = rm->getGammaI(i) * (*tauPtr);

      // stage-solution: u_s(i) = u_old + sum_{j=0}^{i-1} a_ij*U_j
      *stageSolution = *unVec;
      for (int j = 0; j < i; j++)
      {
        *tmp = *(stageSolutions[j]);
        *tmp *= rm->getA(i, j);
        *stageSolution += *tmp;
      }

      // Dirichlet-BC implemented as additional algebraic equation u = g(x,t) on boundary
      // => U_i = -u_s(i) + g(x,t_s(i)) + tau*gamma_i* d_t(g)(t_old) on boundary
      // where u_s(i) = ith stage-solution, t_s(i) = ith stage-time
      for (unsigned int j = 0; j < boundaries.size(); j++)
      {
        boundaries[j].vec->interpol(boundaries[j].fct);
        *(boundaries[j].vec) -= *(stageSolution->getDOFVector(boundaries[j].col));

        if (boundaries[j].fctDt != NULL)
        {
          // time derivative of dirichlet bc is given
          DOFVector<double> tmpDt(getFeSpace(boundaries[j].col), "tmp");
          tmpDt.interpol(boundaries[j].fctDt);
          tmpDt *= tauGammaI;
          *(boundaries[j].vec) += tmpDt;
        }
      }

      // timeRhs: sum_{j=0}^{i-1} c_ij / tau * U_j
      timeRhsVec->set(0.0);
      for (int j = 0; j < i; j++)
      {
        *tmp = *(stageSolutions[j]);
        *tmp *= (rm->getC(i, j) / *tauPtr);
        *timeRhsVec += *tmp;
      }

      // assemble and solve stage equation
      ProblemStat::buildAfterCoarsen(adaptInfo, flag, (i == 0), true);
#if defined HAVE_PARALLEL_PETSC
      // TODO: Problems with reuse of Matrix with parallel PETSC-Solvers
      //       Thus, Rosenbrock not efficient but working (Rainer)
      ProblemStat::solve(adaptInfo, true , false);
#else
      ProblemStat::solve(adaptInfo, i == 0, i + 1 < rm->getStages());
#endif

      *(stageSolutions[i]) = *solution;

      *tmp = *solution;
      *tmp *= rm->getM1(i);
      *newUn += *tmp;

      *tmp = *solution;
      *tmp *= rm->getM2(i);
      *lowSol += *tmp;
    }

    stageTime = oldTime + (*tauPtr);

    Flag flag_;
    return flag_;
  }


  void RosenbrockStationary::estimateTimeError(AdaptInfo& adaptInfo)
  {
    for (int i = 0; i < nComponents; i++)
    {
      (*(lowSol->getDOFVector(i))) -= (*(newUn->getDOFVector(i)));
      adaptInfo.setTimeEstSum(lowSol->getDOFVector(i)->L2Norm(), i+componentShift);
    }
  }


  void RosenbrockStationary::acceptTimestep(AdaptInfo& adaptInfo)
  {
    *solution = *newUn;
    *unVec = *newUn;
    oldTime = adaptInfo.getTime();
  }


  void RosenbrockStationary::addOperator(Operator& op, int row, int col,
                                         double* factor, double* estFactor)
  {
    FUNCNAME("RosenbrockStationary::addOperator()");

    TEST_EXIT(op.getUhOld() == NULL)("UhOld is not allowed to be set!\n");

    op.setUhOld(stageSolution->getDOFVector(col));
    ProblemStat::addVectorOperator(op, row, factor, estFactor);
  }


  void RosenbrockStationary::addJacobianOperator(Operator& op, int row, int col,
      double* factor, double* estFactor)
  {
    FUNCNAME("RosenbrockStationary::addJacobianOperator()");

    TEST_EXIT(factor == NULL)("Not yet implemented!\n");
    TEST_EXIT(estFactor == NULL)("Not yet implemented!\n");

    ProblemStat::addMatrixOperator(op, row, col, &minusOne, &minusOne);
  }


  void RosenbrockStationary::addTimeOperator(int row, int col)
  {
    FUNCNAME("RosenbrockStationary::addTimeOperator()");

    TEST_EXIT(invTauGamma)("This should not happen!\n");

    Operator* op = new Operator(componentSpaces[row], componentSpaces[col]);
    op->addZeroOrderTerm(new Simple_ZOT);
    ProblemStat::addMatrixOperator(op, row, col, invTauGamma, invTauGamma);

    Operator* opRhs = new Operator(componentSpaces[row]);
    opRhs->addZeroOrderTerm(new VecAtQP_ZOT(timeRhsVec->getDOFVector(col)));
    ProblemStat::addVectorOperator(opRhs, row);
  }


  void RosenbrockStationary::addDirichletBC(BoundaryType type, int row, int col,
      AbstractFunction<double, WorldVector<double>>* fct)
  {
    DOFVector<double>* vec = new DOFVector<double>(componentSpaces[row], "vec");
    RosenbrockBoundary bound(fct, NULL, vec, row, col);
    boundaries.push_back(bound);

    ProblemStat::addDirichletBC(type, row, col, vec);
  }


  void RosenbrockStationary::addDirichletBC(BoundaryType type, int row, int col,
      AbstractFunction<double, WorldVector<double>>* fct,
      AbstractFunction<double, WorldVector<double>>* fctDt)
  {
    DOFVector<double>* vec = new DOFVector<double>(componentSpaces[col], "vec");
    RosenbrockBoundary bound(fct, fctDt, vec, row, col);
    boundaries.push_back(bound);

    ProblemStat::addDirichletBC(type, row, col, vec);
  }

}
