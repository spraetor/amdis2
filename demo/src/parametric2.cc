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

struct TransformMesh : AbstractFunction<WorldVector<double>, WorldVector<double> >
{
  TransformMesh(double a_, double b_, double s_, double theta_) : a(a_), b(b_), s(s_), theta(theta_) {}
  
  WorldVector<double> operator()(const WorldVector<double>& x) const
  {
    WorldVector<double> x_bar, x_hat;
    x_bar[0] = a + (b-a)*pow(x[0], s);
    x_bar[1] = x[1];
    
    x_hat[0] = x_bar[0]*cos(theta*x_bar[1]);
    x_hat[1] = x_bar[0]*sin(theta*x_bar[1]);
    return x_hat;
  }
  
private:
  double a;
  double b;
  double s;
  double theta;
};

/** \brief
 * Implementation of a scalar problem. In \ref buildBeforeCoarsen() parametric
 * coordinates for the vertices are filled in a DOFVector. This DOFVector is
 * used in \ref parametric to parametrize elements while mesh traversal. 
 */
class HollowCylinder : public ProblemStat
{
public:
  /// constructor
  HollowCylinder(std::string name) 
    : ProblemStat(name),
      parametric(NULL),
      a(1.0),
      b(5.0),
      s(1.3),
      theta(m_pi/2.0)
  {
    Parameters::get(name + "->a", a);
    Parameters::get(name + "->b", b);
    Parameters::get(name + "->s", s);
    Parameters::get(name + "->theta", theta);
  }

  /// destructor
  ~HollowCylinder() 
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

    // ===== create coords vector =====
    for (int i = 0; i < Global::getGeo(WORLD); i++)
      parametricCoords[i] = new DOFVector<double>(this->getFeSpace(), 
						  "parametric coords "+boost::lexical_cast<std::string>(i));

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
    int dim = this->getMesh()->getDim();
    int dow = Global::getGeo(WORLD);
    WorldVector<double> vertexCoords, newCoords;

    DegreeOfFreedom dof;
    WorldVector<double> x;
  
    const FiniteElemSpace *feSpace = this->getFeSpace();
    const BasisFunction *basFcts = feSpace->getBasisFcts();
    int numBasFcts = basFcts->getNumber();
    std::vector<DegreeOfFreedom> localIndices(numBasFcts);
    DOFAdmin *admin = feSpace->getAdmin();
    
    // ===== disable parametric traverse =====
    this->getMesh()->setParametric(NULL);

    std::map<DegreeOfFreedom, bool> visited;

    TraverseStack stack;
    ElInfo *elInfo = stack.traverseFirst(this->getMesh(), -1, 
					Mesh::CALL_LEAF_EL | Mesh::FILL_COORDS);
    while (elInfo) {
      basFcts->getLocalIndices(elInfo->getElement(), admin, localIndices);
      for (int i = 0; i < dim + 1; i++) {
	dof = localIndices[i];
	vertexCoords = elInfo->getCoord(i);
	newCoords = TransformMesh(a,b,s,theta)(vertexCoords);
	if (!visited[dof]) {
	  for (int j = 0; j < dow; j++)
	    (*(parametricCoords[j]))[dof] = newCoords[j];

	  visited[dof] = true;
	}
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
  
  double a;
  double b;
  double s;
  double theta;
};

// ===========================================================================
// ===== main program ========================================================
// ===========================================================================

int main(int argc, char* argv[])
{
  FUNCNAME("parametric main");

  AMDiS::init(argc, argv);

  // ===== create and init the scalar problem ===== 
  HollowCylinder parametric("parametric");
  parametric.initialize(INIT_ALL);

  // === create adapt info ===
  AdaptInfo *adaptInfo = new AdaptInfo("parametric->adapt", parametric.getNumComponents());

  // === create adapt ===
  AdaptStationary *adapt = new AdaptStationary("parametric->adapt",
					       &parametric,
					       adaptInfo);
  
  // ===== create matrix operator =====
//   Operator matrixOperator(parametric.getFeSpace());
//   matrixOperator.addTerm(new Simple_ZOT);
//   parametric.addMatrixOperator(&matrixOperator, 0, 0);

  // ===== create rhs operator =====
//   Operator rhsOperator(parametric.getFeSpace());

//   int degree = parametric.getFeSpace()->getBasisFcts()->getDegree();

//   rhsOperator.addTerm(new CoordsAtQP_ZOT(new F(degree)));
//   parametric.addVectorOperator(&rhsOperator, 0);

  // ===== start adaption loop =====
  adapt->adapt();

  parametric.writeFiles(adaptInfo, true);

  AMDiS::finalize();
}


