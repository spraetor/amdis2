#include "AMDiS.h"

using namespace AMDiS;

class YRotation
{
public:
  static WorldVector<double>& rotate(WorldVector<double> &x, double angle) 
  {
    double x0 = x[0] * cos(angle) + x[2] * sin(angle);
    x[2] = -x[0] * sin(angle) + x[2] * cos(angle);
    x[0] = x0;
    return x;
  }
};

/// RHS function
class F : public AbstractFunction<double, WorldVector<double> >
{
public:
  F(int degree) 
    : AbstractFunction<double, WorldVector<double> >(degree),
      rotation(0.0)
  {}

  /// Implementation of AbstractFunction::operator().
  double operator()(const WorldVector<double>& x) const 
  {
    WorldVector<double> myX = x;
    YRotation::rotate(myX, -rotation);
    return -2.0 * myX[0];
  }

  void rotate(double r) 
  { 
    rotation += r; 
  }

private:
  double rotation;
};

class TorusProject : public Projection
{
public:
  /// Constructor.
  TorusProject(int id, 
	       ProjectionType type,
	       double radius1,
	       double radius2) 
    : Projection(id, type),
      radius1_(radius1),
      radius2_(radius2)
  {}

  /// Destructor.
  virtual ~TorusProject() {}

  /// Implementation of Projection::project();
  void project(WorldVector<double> &x) 
  {

    WorldVector<double> xPlane = x;
    xPlane[2] = 0.0;

    double norm = sqrt(xPlane*xPlane);
    TEST_EXIT(norm != 0.0)("can't project vector x\n");
    
    WorldVector<double> center = xPlane;
    center *= radius1_ / norm;

    x -= center;

    norm = sqrt(x*x);
    TEST_EXIT(norm != 0.0)("can't project vector x\n");
    x *= radius2_/norm;

    x += center;
  }

protected:
  double radius1_;

  double radius2_;
};


int main(int argc, char* argv[])
{
  FUNCNAME("torus main");

  AMDiS::init(argc, argv);

  // ===== create projection =====
  double r2 = (1.5 - 1.0 / sqrt(2.0)) / 2.0;
  double r1 = 1.0 / sqrt(2.0) + r2;

  new TorusProject(1, 
		   VOLUME_PROJECTION, 
		   r1,
		   r2);

  // ===== create and init the scalar problem ===== 
  ProblemStat torus("torus");
  torus.initialize(INIT_ALL);

  // === create adapt info ===
  AdaptInfo *adaptInfo = new AdaptInfo("torus->adapt", torus.getNumComponents());

  // === create adapt ===
  AdaptStationary *adapt = new AdaptStationary("torus->adapt",
					       &torus,
					       adaptInfo);
  
  // ===== create matrix operator =====
  Operator matrixOperator(torus.getFeSpace());
  matrixOperator.addTerm(new Simple_SOT);
  torus.addMatrixOperator(&matrixOperator, 0, 0);

  // ===== create rhs operator =====
  Operator rhsOperator(torus.getFeSpace());

  int degree = torus.getFeSpace()->getBasisFcts()->getDegree();

  F f(degree);
  rhsOperator.addTerm(new CoordsAtQP_ZOT(&f));
  torus.addVectorOperator(&rhsOperator, 0);

  // ===== start adaption loop =====
  adapt->adapt();

  torus.writeFiles(adaptInfo, true);


  double rotation = M_PI / 3.0;
  int dim = torus.getMesh()->getDim();
  int dow = Global::getGeo(WORLD);

  DegreeOfFreedom dof;
  WorldVector<double> x;
 
  const FiniteElemSpace *feSpace = torus.getFeSpace();
  const BasisFunction *basFcts = feSpace->getBasisFcts();
  int numBasFcts = basFcts->getNumber();
  std::vector<DegreeOfFreedom> localIndices(numBasFcts);
  DOFAdmin *admin = feSpace->getAdmin();
 
  WorldVector<DOFVector<double>*> parametricCoords;
  for (int i = 0; i < dow; i++)
    parametricCoords[i] = new DOFVector<double>(feSpace, "parametric coords");

  std::map<DegreeOfFreedom, bool> visited;

  TraverseStack stack;
  ElInfo *elInfo = stack.traverseFirst(torus.getMesh(), -1, 
				       Mesh::CALL_LEAF_EL | Mesh::FILL_COORDS);
  while (elInfo) {
    basFcts->getLocalIndices(elInfo->getElement(), admin, localIndices);
    for (int i = 0; i < dim + 1; i++) {
      dof = localIndices[i];
      x = elInfo->getCoord(i);
      YRotation::rotate(x, rotation);
      if (!visited[dof]) {
	for (int j = 0; j < dow; j++)
	  (*(parametricCoords[j]))[dof] = x[j];

	visited[dof] = true;
      }
    }
    elInfo = stack.traverseNext(elInfo);
  }

  ParametricFirstOrder parametric(&parametricCoords);
  torus.getMesh()->setParametric(&parametric);

  f.rotate(rotation);
  adaptInfo->reset();
  adapt->adapt();

  visited.clear();
  elInfo = stack.traverseFirst(torus.getMesh(), -1, 
			       Mesh::CALL_LEAF_EL | Mesh::FILL_COORDS);
  while (elInfo) {
    basFcts->getLocalIndices(elInfo->getElement(), admin, localIndices);
    for (int i = 0; i < dim + 1; i++) {
      dof = localIndices[i];
      x = elInfo->getCoord(i);
      YRotation::rotate(x, rotation);
      if (!visited[dof]) {
	for (int j = 0; j < dow; j++)
	  (*(parametricCoords[j]))[dof] = x[j];

	visited[dof] = true;
      }
    }
    elInfo = stack.traverseNext(elInfo);
  }

  f.rotate(rotation);
  adaptInfo->reset();
  adapt->adapt();

  torus.writeFiles(adaptInfo, true);

  for (int i = 0; i < dow; i++)
    delete parametricCoords[i];

  AMDiS::finalize();
}


