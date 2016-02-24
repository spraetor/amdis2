#include "AMDiS.h"

using namespace std;
using namespace AMDiS;

// ===========================================================================
// ===== function definitions ================================================
// ===========================================================================

/// RHS function
class F : public AbstractFunction<double, WorldVector<double> >
{
public:
  F(int degree) : AbstractFunction<double, WorldVector<double> >(degree) {}

  /// Implementation of AbstractFunction::operator().
  double operator()(const WorldVector<double>& x) const 
  {
    return -2.0 * x[0];
  }
};

/** \brief
 * Implementation of a scalar problem. In \ref buildBeforeCoarsen() parametric
 * coordinates for the vertices are filled in a DOFVector. This DOFVector is
 * used in \ref parametric to parametrize elements while mesh traversal. 
 */
class ParametricSphere : public ProblemStat
{
public:
  /// constructor
  ParametricSphere(const char *name) 
    : ProblemStat(name),
      parametric(NULL)
  {}

  /// destructor
  ~ParametricSphere() 
  {
    for (int i = 0; i < Global::getGeo(WORLD); i++)
      delete parametricCoords[i];

    delete parametric;
  }

  /// initialization of the base class and creation of \ref parametric.
  void initialize(Flag initFlag,
		  ProblemStat *adoptProblem = NULL,
		  Flag adoptFlag = INIT_NOTHING)
  {
    ProblemStat::initialize(initFlag, adoptProblem, adoptFlag);

    // ===== create projection =====
    WorldVector<double> ballCenter;
    ballCenter.set(0.0);
    new BallProject(1, 
		    VOLUME_PROJECTION, 
		    ballCenter, 
		    1.0);

    // ===== create coords vector =====
    for (int i = 0; i < Global::getGeo(WORLD); i++)
      parametricCoords[i] = new DOFVector<double>(this->getFeSpace(), 
						  "parametric coords");

    // ===== create parametric object =====
    parametric = new ParametricFirstOrder(&parametricCoords);

    // ===== enable parametric traverse =====
    this->getMesh()->setParametric(parametric);
  }


  /** \brief
   * Implementation of ProblemStatBase::buildBeforeCoarsen().
   */
  void buildBeforeCoarsen(AdaptInfo *adaptInfo, Flag flag) 
  {
    FUNCNAME("ParametricSphere::buildAfterCoarsen()");
    MSG("calculation of parametric coordinates\n");
    int preDOFs = this->getFeSpace()->getAdmin()->getNumberOfPreDofs(VERTEX);
    int dim = this->getMesh()->getDim();
    int dow = Global::getGeo(WORLD);
    WorldVector<double> vertexCoords;
    const DegreeOfFreedom **dof;
    DegreeOfFreedom vertexIndex;

    // ===== disable parametric traverse =====
    this->getMesh()->setParametric(NULL);

    TraverseStack stack;
    ElInfo *elInfo = NULL;
    elInfo = stack.traverseFirst(this->getMesh(), -1, 
				 Mesh::CALL_LEAF_EL | Mesh::FILL_COORDS);
    while (elInfo) {
      dof = elInfo->getElement()->getDof();
      for (int i = 0; i < dim + 1; i++) {
	vertexCoords = elInfo->getCoord(i);
	Projection::getProjection(1)->project(vertexCoords);
	for (int j = 0; j < dow; j++)
	  (*(parametricCoords[j]))[dof[i][preDOFs]] = vertexCoords[j];
      }
      elInfo = stack.traverseNext(elInfo);
    }
    
    // ===== enable parametric traverse =====
    this->getMesh()->setParametric(parametric);
  }


protected:
  /// DOFVector storing parametric coordinates.
  WorldVector<DOFVector<double>*> parametricCoords;

  /// Parametrizes vertex coordinates while mesh traversal.
  ParametricFirstOrder *parametric;
};

// ===========================================================================
// ===== main program ========================================================
// ===========================================================================

int main(int argc, char* argv[])
{
  FUNCNAME("parametric main");

  AMDiS::init(argc, argv);

  // ===== create and init the scalar problem ===== 
  ParametricSphere parametric("parametric");
  parametric.initialize(INIT_ALL);

  // === create adapt info ===
  AdaptInfo *adaptInfo = new AdaptInfo("parametric->adapt", parametric.getNumComponents());

  // === create adapt ===
  AdaptStationary *adapt = new AdaptStationary("parametric->adapt",
					       &parametric,
					       adaptInfo);
  
  // ===== create matrix operator =====
  Operator matrixOperator(parametric.getFeSpace());
  matrixOperator.addTerm(new Simple_SOT);
  parametric.addMatrixOperator(&matrixOperator, 0, 0);

  // ===== create rhs operator =====
  Operator rhsOperator(parametric.getFeSpace());

  int degree = parametric.getFeSpace()->getBasisFcts()->getDegree();

  rhsOperator.addTerm(new CoordsAtQP_ZOT(new F(degree)));
  parametric.addVectorOperator(&rhsOperator, 0);

  // ===== start adaption loop =====
  adapt->adapt();

  parametric.writeFiles(adaptInfo, true);

  AMDiS::finalize();
}


