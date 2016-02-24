#pragma once

#include <list>
#include <vector>

#include <boost/numeric/mtl/mtl_fwd.hpp>

#include "AMDiS_fwd.hpp"
#include "BasisFunction.hpp"
#include "FixVec.hpp"
#include "Flag.hpp"

namespace AMDiS
{
  /**
   * \ingroup Assembler
   *
   * \brief
   * For the assemblage of the system matrix and right hand side vector of the
   * linear system, we have to compute integrals, for example:
   * \f[ \int_{\Omega} f(x)\varphi_i(x) dx \f]
   * For general data A, b, c, and f, these integrals can not be calculated
   * exactly. Quadrature formulas have to be used in order to calculate the
   * integrals approximately. Numerical integration in finite element methods is
   * done by looping over all grid elements and using a quadrature formula on
   * each element.
   */
  class Quadrature
  {
  protected:
    /// Avoids call of default constructor
    Quadrature();

    /** \brief
     * Constructs a Quadrature with name name_ of degree degree_ for dim dim_.
     * The Quadrature has n_points_ quadrature points with barycentric
     * coordinates lambda_ and weights w_. The constructor is protected because
     * access to a Quadrature should be done via \ref provideQuadrature.
     */
    Quadrature(const char* name_,
               int degree_,
               int dim_,
               int n_points_,
               VectorOfFixVecs<DimVec<double>>* lambda_,
               double* w_)
      : name(name_),
        degree(degree_),
        dim(dim_),
        n_points(n_points_),
        lambda(lambda_),
        w(w_)
    {}

  public:
    /// Copy constructor
    Quadrature(const Quadrature&);

    /// Destructor
    virtual ~Quadrature();

    /// Returns a Quadrature for dimension dim exact for degree degree.
    static Quadrature* provideQuadrature(int dim, int degree);

    /** \brief
     * Approximates an integral by the numerical quadrature described by quad;
     * f is a pointer to an AbstractFunction to be integrated, evaluated in
     * barycentric coordinates; the return value is
     * \f[ \sum_{k = 0}^{n_points-1} w[k] * (*f)(lambda[k]) \f]
     * For the approximation of \f$ \int_S f\f$ we have to multiply this value
     * with d!|S| for a simplex S; for a parametric simplex f should be a pointer
     * to a function which calculates
     * \f$ f(\lambda)|det DF_S(\hat{x}(\lambda))| \f$.
     */
    double integrateStdSimplex(std::function<double(DimVec<double>)> f);

    /// Returns \ref name
    std::string getName() const
    {
      return name;
    }

    /// Returns \ref n_points
    int getNumPoints() const
    {
      return n_points;
    }

    /// Returns \ref w[p]
    double getWeight(int p) const
    {
      return w[p];
    }

    /// Returns \ref w.
    double* getWeight() const
    {
      return w;
    }

    /// Returns \ref dim
    int getDim() const
    {
      return dim;
    }

    /// Returns \ref degree
    int getDegree() const
    {
      return degree;
    }

    /** \brief
     * Returns a pointer to a vector storing the values of a doubled valued
     * function at all quadrature points; f is that AbstractFunction
     * , evaluated in barycentric coordinates; if vec is not NULL, the values are
     * stored in this vector, otherwise the values are stored in some static
     * local vector, which is overwritten on the next call
     */
    const double* fAtQp(const std::function<double(DimVec<double>)>& f,
                        double* vec) const ;

    /** \brief
     * Returns a pointer to a vector storing the gradient (with respect to world
     * coordinates) of a double valued function at all quadrature points;
     * grdF is a pointer to a AbstractFunction, evaluated in barycentric
     * coordinates and returning a pointer to a WorldVector storing the gradient;
     * if vec is not NULL, the values are stored in this vector, otherwise the
     * values are stored in some local static vector, which is overwritten on the
     * next call
     */
    const WorldVector<double>* grdFAtQp(const std::function<WorldVector<double>(DimVec<double>)>& grdF,
                                        WorldVector<double>* vec) const;



    /// Returns \ref lambda[a][b] which is the b-th coordinate entry of the a-th
    /// quadrature point
    double getLambda(int a, int b) const
    {
      return (lambda ? (*lambda)[a][b] : 0.0);
    }

    /// Returns \ref lambda[a] which is a DimVec<double> containing the
    /// coordiantes of the a-th quadrature point
    const DimVec<double>& getLambda(int a) const
    {
      return (*lambda)[a];
    }

    /// Returns \ref lambda which is a VectorOfFixvecs<DimVec<double> >.
    VectorOfFixVecs<DimVec<double>>* getLambda() const
    {
      return lambda;
    }


  public:
    /// Maximal number of quadrature points for the different dimensions
    static constexpr int maxNQuadPoints[4] = {0, 10, 61, 64};

  protected:
    /// Name of this Quadrature
    std::string name;

    /// Quadrature is exact of this degree
    int degree;

    /// Quadrature for dimension dim
    int dim;

    /// Number of quadrature points
    int n_points;

    /// Vector of quadrature points given in barycentric coordinates
    VectorOfFixVecs<DimVec<double>>* lambda;

    /// Vector of quadrature weights
    double* w;

  protected:
    /// Initialisation of all static Quadrature objects which will be returned
    /// by \ref provideQuadrature()
    static void initStaticQuadratures();

    /** \name static quadratures, used weights, and barycentric coords
     * \{
     */
    static Quadrature** quad_nd[4];
    static Quadrature* quad_0d[1];
    static Quadrature* quad_1d[20];
    static Quadrature* quad_2d[18];
    static Quadrature* quad_3d[8];

    static VectorOfFixVecs<DimVec<double>>* x_0d;
    static double* w_0d;

    static VectorOfFixVecs<DimVec<double>>* x0_1d;
    static VectorOfFixVecs<DimVec<double>>* x1_1d;
    static VectorOfFixVecs<DimVec<double>>* x2_1d;
    static VectorOfFixVecs<DimVec<double>>* x3_1d;
    static VectorOfFixVecs<DimVec<double>>* x4_1d;
    static VectorOfFixVecs<DimVec<double>>* x5_1d;
    static VectorOfFixVecs<DimVec<double>>* x6_1d;
    static VectorOfFixVecs<DimVec<double>>* x7_1d;
    static VectorOfFixVecs<DimVec<double>>* x8_1d;
    static VectorOfFixVecs<DimVec<double>>* x9_1d;
    static double* w0_1d;
    static double* w1_1d;
    static double* w2_1d;
    static double* w3_1d;
    static double* w4_1d;
    static double* w5_1d;
    static double* w6_1d;
    static double* w7_1d;
    static double* w8_1d;
    static double* w9_1d;

    static VectorOfFixVecs<DimVec<double>>* x1_2d;
    static VectorOfFixVecs<DimVec<double>>* x2_2d;
    static VectorOfFixVecs<DimVec<double>>* x3_2d;
    static VectorOfFixVecs<DimVec<double>>* x4_2d;
    static VectorOfFixVecs<DimVec<double>>* x5_2d;
    static VectorOfFixVecs<DimVec<double>>* x7_2d;
    static VectorOfFixVecs<DimVec<double>>* x8_2d;
    static VectorOfFixVecs<DimVec<double>>* x9_2d;
    static VectorOfFixVecs<DimVec<double>>* x10_2d;
    static VectorOfFixVecs<DimVec<double>>* x11_2d;
    static VectorOfFixVecs<DimVec<double>>* x12_2d;
    static VectorOfFixVecs<DimVec<double>>* x17_2d;
    static double* w1_2d;
    static double* w2_2d;
    static double* w3_2d;
    static double* w4_2d;
    static double* w5_2d;
    static double* w7_2d;
    static double* w8_2d;
    static double* w9_2d;
    static double* w10_2d;
    static double* w11_2d;
    static double* w12_2d;
    static double* w17_2d;

    static VectorOfFixVecs<DimVec<double>>* x1_3d;
    static VectorOfFixVecs<DimVec<double>>* x2_3d;
    static VectorOfFixVecs<DimVec<double>>* x3_3d;
    static VectorOfFixVecs<DimVec<double>>* x4_3d;
    static VectorOfFixVecs<DimVec<double>>* x5_3d;
    static VectorOfFixVecs<DimVec<double>>* x7_3d;
    static double* w1_3d;
    static double* w2_3d;
    static double* w3_3d;
    static double* w4_3d;
    static double* w5_3d;
    static double* w7_3d;

    /** \} */
  };



  /// Pre-compute the values of all basis functions at all quadrature nodes;
  const Flag INIT_PHI=1;

  /// Pre-compute the gradients (with respect to the barycentric coordinates) of
  /// all basis functions at all quadrature nodes
  const Flag INIT_GRD_PHI=2;

  /// pre-compute all 2nd derivatives (with respect to the barycentric
  /// coordinates) of all basis functions at all quadrature nodes;
  const Flag INIT_D2_PHI=4;


  /**
   * \ingroup Integration
   *
   *\brief
   * Often numerical integration involves basis functions, such as the assembling
   * of the system matrix and right hand side, or the integration of finite
   * element functions. Since numerical quadrature involves only the values at
   * the quadrature points and the values of basis functions and its derivatives
   * are the same at these points for all elements of the grid, such routines can
   * be much more efficient, if they can use pre-computed values of the basis
   * functions at the quadrature points. In this case the basis functions do not
   * have to be evaluated for each quadrature point on every element newly.
   * Information that should be pre-computed can be specified by the following
   * symbolic constants:
   * \ref INIT_PHI, \ref INIT_GRD_PHI, \ref INIT_D2_PHI
   */
  class FastQuadrature
  {
  protected:
    /// Constructs a FastQuadrature for the given Quadrature, BasisFunction, and
    /// flag.
    FastQuadrature(BasisFunction* basFcts, Quadrature* quad, Flag flag)
      : init_flag(flag),
        phi(0, 0),
        D2Phi(NULL),
        quadrature(quad),
        basisFunctions(basFcts)
    {}

    /// Copy constructor
    FastQuadrature(const FastQuadrature&);

    /// Extended copy constructor
    FastQuadrature(const FastQuadrature&, const Flag);

    /// Destructor
    virtual ~FastQuadrature();

  public:
    /// Returns a FastQuadrature for the given BasisFunction, Quadrature, and flags.
    static FastQuadrature* provideFastQuadrature(const BasisFunction*,
        const Quadrature&,
        Flag);

    /// inits FastQuadrature like speciefied in flag
    void init(Flag init_flag);

    bool initialized(Flag flag)
    {
      if (flag == INIT_PHI)
        return (num_rows(phi) > 0);

      if (flag == INIT_GRD_PHI)
        return (!grdPhi.empty());

      if (flag == INIT_D2_PHI)
        return (D2Phi != NULL);

      ERROR_EXIT("invalid flag\n");
      return false;
    }

    /// Returns \ref quadrature
    const Quadrature* getQuadrature() const
    {
      return quadrature;
    }

    /// Returns \ref max_points
    int getMaxQuadPoints() const
    {
      return max_points;
    }

    /// Returns (*\ref D2Phi)[q][i][j][m]
    double getSecDer(int q, int i, int j, int m) const;

    /// Returns (*\ref D2Phi)[q]
    const VectorOfFixVecs<DimMat<double>>* getSecDer(int q) const;

    /// Returns (*\ref grdPhi)[q][i][j]
    double getGradient(int q, int i ,int j) const
    {
      return (!grdPhi.empty()) ? grdPhi[q][i][j] : 0.0;
    }

    /// Returns (*\ref grdPhi)[q]
    const std::vector<DenseVector<double>>& getGradient(int q) const
    {
      return grdPhi[q];
    }

    const DenseVector<double>& getGradient(int q, int i) const
    {
      return grdPhi[q][i];
    }

    const DenseMatrix<double>& getPhi() const
    {
      return phi;
    }

    /// Returns \ref phi[q][i]
    double getPhi(int q, int i) const
    {
      return phi[q][i];
    }

    /// Returns \ref quadrature ->integrateStdSimplex(f)
    double integrateStdSimplex(std::function<double(DimVec<double>)> f)
    {
      return quadrature->integrateStdSimplex(f);
    }

    /// Returns \ref quadrature ->getNumPoints()
    int getNumPoints() const
    {
      return quadrature->getNumPoints();
    }

    /// Returns \ref quadrature ->getWeight(p)
    double getWeight(int p) const
    {
      return quadrature->getWeight(p);
    }

    /// Returns \ref quadrature ->getDim()
    int getDim() const
    {
      return quadrature->getDim();
    }

    /// Returns \ref quadrature ->getDegree()
    int getDegree() const
    {
      return quadrature->getDegree();
    }

    /// Returns \ref quadrature ->grdFAtQp(f, vec)
    const WorldVector<double>
    * grdFAtQp(const std::function<WorldVector<double>(DimVec<double>)>& f,
               WorldVector<double>* vec) const
    {
      return quadrature->grdFAtQp(f, vec);
    }

    /// Returns \ref quadrature ->fAtQp(f, vec)
    const double* fAtQp(const std::function<double(DimVec<double>)>& f, double* vec) const
    {
      return quadrature->fAtQp(f, vec);
    }

    /// Returns \ref quadrature ->getLambda(a,b)
    double getLambda(int a,int b) const
    {
      return quadrature->getLambda(a,b);
    }

    /// Returns \ref quadrature ->getLambda(a)
    const DimVec<double>& getLambda(int a) const
    {
      return quadrature->getLambda(a);
    }

    /// Returns \ref basisFunctions
    BasisFunction* getBasisFunctions() const
    {
      return basisFunctions;
    }

  protected:
    /// Specifies which information should be pre-computed. Can be \ref INIT_PHI,
    /// \ref INIT_GRD_PHI, or \ref INIT_D2_PHI
    Flag init_flag;

    /** \brief
     * Matrix storing function values if the flag \ref INIT_PHI is set;
     * phi[i][j] stores the value \ref basisFunctions->phi[j]
     * (quadrature->lambda[i]), 0 <= j < basisFunctions->getNumber()  and
     * 0 <= i < n_points
     */
    DenseMatrix<double> phi;

    /** \brief
     * Matrix storing all gradients (with respect to the barycentric coordinates)
     * if the flag \ref INIT_GRD_PHI is set; grdPhi[i][j][k] stores the value
     * basisFunctions->grdPhi[j](quadrature->lambda[i])[k]
     * for 0 <= j < basisFunctions->getNumber(),
     * 0 <= i < . . . , n_points, and 0 <= k < DIM
     */
    std::vector<std::vector<DenseVector<double>>> grdPhi;

    /** \brief
     * Matrix storing all second derivatives (with respect to the barycentric
     * coordinates) if the flag \ref INIT_D2_PHI is set; D2Phi[i][j][k][l] stores
     * the value basisFunctions->D2Phi[j](quadrature->lambda[i])[k][l]
     * for 0 <= j < basisFunctions->getNumber(),
     * 0 <= i < n_points, and 0 <= k,l < DIM
     */
    MatrixOfFixVecs<DimMat<double>>* D2Phi;

    /// List of all used FastQuadratures
    static std::list<FastQuadrature*> fastQuadList;

    /// Maximal number of quadrature points for all yet initialised FastQuadrature
    /// objects. This value may change after a new initialisation of a
    /// FastQuadrature
    static int max_points;

    /// This FastQuadrature stores values for Quadrature quadrature
    Quadrature* quadrature;

    /// Values stored for basis functions basisFunctions
    BasisFunction* basisFunctions;

  };

} // end namespace AMDiS
