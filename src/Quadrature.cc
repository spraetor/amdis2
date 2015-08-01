#include <algorithm>

#include "Quadrature.h"
#include "FixVec.h"

// TODO: extend quadrature formulas in 3d, extend precision of formulas

using namespace std;

namespace AMDiS 
{
  constexpr int Quadrature::maxNQuadPoints[4];

  list<FastQuadrature*> FastQuadrature::fastQuadList;
  int FastQuadrature::max_points = 0;

  Quadrature::Quadrature(const Quadrature& q)
  {
    name = q.name;
    degree = q.degree;
    dim = q.dim;
    n_points = q.n_points;
  
    // copy barycentric coordinates
    lambda = new VectorOfFixVecs<DimVec<double> >(*(q.lambda));

    // copy weights
    w = new double[n_points];

    for (int i = 0; i < n_points; i++)
      w[i] = q.w[i];
  }

  /****************************************************************************/
  /*  initialize gradient values of a function f in local coordinates at the  */
  /*  quadrature points                                                       */
  /****************************************************************************/

  const WorldVector<double>
  *Quadrature::grdFAtQp(const std::function<WorldVector<double>(DimVec<double>)>& f, 
	                      WorldVector<double>* vec) const
  {
    static WorldVector<double> *quad_vec_d = NULL;
    static int size = 0;
    WorldVector<double> *val;
    WorldVector<double> grd;

    if (vec) {
      val = vec;
    } else {
      if (size < n_points) {
      	int new_size = std::max(maxNQuadPoints[dim], n_points);
      	if (quad_vec_d)
      	  delete [] quad_vec_d;
      	quad_vec_d = new WorldVector<double>[new_size];
      	size = new_size;
      }
      val = quad_vec_d;
    }

    int dow = Global::getGeo(WORLD);

    for (int i = 0; i < n_points; i++) {
      grd = f((*lambda)[i]);
      for (int j = 0; j < dow; j++)
        val[i][j] = grd[j];
    }
    
    return val;
  }

  const double *Quadrature::fAtQp(const std::function<double(DimVec<double>)>& f,
				  double *vec) const
  {
    static double *quad_vec = NULL;
    static int size = 0;
    double *val;
 
    if (vec) {
      val = vec;
    } else {
      if (size < n_points) {
        int new_size = std::max(maxNQuadPoints[dim], n_points);
        if (quad_vec)
          delete [] quad_vec;
        quad_vec = new double[new_size];
        size = new_size;
      }
      val = quad_vec;
    }

    for (int i = 0; i < n_points; i++)
      val[i] = f((*lambda)[i]);
    
    return const_cast<const double *>(val);
  }


  Quadrature **Quadrature::quad_nd[4];
  Quadrature *Quadrature::quad_0d[1];
  Quadrature *Quadrature::quad_1d[20];
  Quadrature *Quadrature::quad_2d[18];
  Quadrature *Quadrature::quad_3d[8];

  VectorOfFixVecs<DimVec<double> > *Quadrature::x_0d;
  double *Quadrature::w_0d;

  VectorOfFixVecs<DimVec<double> > *Quadrature::x0_1d = NULL;
  VectorOfFixVecs<DimVec<double> > *Quadrature::x1_1d;
  VectorOfFixVecs<DimVec<double> > *Quadrature::x2_1d;
  VectorOfFixVecs<DimVec<double> > *Quadrature::x3_1d;
  VectorOfFixVecs<DimVec<double> > *Quadrature::x4_1d;
  VectorOfFixVecs<DimVec<double> > *Quadrature::x5_1d;
  VectorOfFixVecs<DimVec<double> > *Quadrature::x6_1d;
  VectorOfFixVecs<DimVec<double> > *Quadrature::x7_1d;
  VectorOfFixVecs<DimVec<double> > *Quadrature::x8_1d;
  VectorOfFixVecs<DimVec<double> > *Quadrature::x9_1d;
  double *Quadrature::w0_1d;
  double *Quadrature::w1_1d;
  double *Quadrature::w2_1d;
  double *Quadrature::w3_1d;
  double *Quadrature::w4_1d;
  double *Quadrature::w5_1d;
  double *Quadrature::w6_1d;
  double *Quadrature::w7_1d;
  double *Quadrature::w8_1d;
  double *Quadrature::w9_1d;

  VectorOfFixVecs<DimVec<double> > *Quadrature::x1_2d;
  VectorOfFixVecs<DimVec<double> > *Quadrature::x2_2d;
  VectorOfFixVecs<DimVec<double> > *Quadrature::x3_2d;
  VectorOfFixVecs<DimVec<double> > *Quadrature::x4_2d;
  VectorOfFixVecs<DimVec<double> > *Quadrature::x5_2d;
  VectorOfFixVecs<DimVec<double> > *Quadrature::x7_2d;
  VectorOfFixVecs<DimVec<double> > *Quadrature::x8_2d;
  VectorOfFixVecs<DimVec<double> > *Quadrature::x9_2d;
  VectorOfFixVecs<DimVec<double> > *Quadrature::x10_2d;
  VectorOfFixVecs<DimVec<double> > *Quadrature::x11_2d;
  VectorOfFixVecs<DimVec<double> > *Quadrature::x12_2d;
  VectorOfFixVecs<DimVec<double> > *Quadrature::x17_2d;
  double *Quadrature::w1_2d;
  double *Quadrature::w2_2d;
  double *Quadrature::w3_2d;
  double *Quadrature::w4_2d;
  double *Quadrature::w5_2d;
  double *Quadrature::w7_2d;
  double *Quadrature::w8_2d;
  double *Quadrature::w9_2d;
  double *Quadrature::w10_2d;
  double *Quadrature::w11_2d;
  double *Quadrature::w12_2d;
  double *Quadrature::w17_2d;

  VectorOfFixVecs<DimVec<double> > *Quadrature::x1_3d;
  VectorOfFixVecs<DimVec<double> > *Quadrature::x2_3d;
  VectorOfFixVecs<DimVec<double> > *Quadrature::x3_3d;
  VectorOfFixVecs<DimVec<double> > *Quadrature::x4_3d;
  VectorOfFixVecs<DimVec<double> > *Quadrature::x5_3d;
  VectorOfFixVecs<DimVec<double> > *Quadrature::x7_3d;
  double *Quadrature::w1_3d;
  double *Quadrature::w2_3d;
  double *Quadrature::w3_3d;
  double *Quadrature::w4_3d;
  double *Quadrature::w5_3d;
  double *Quadrature::w7_3d;


  void Quadrature::initStaticQuadratures()
  {

    TEST_EXIT(x0_1d == NULL)("static quadratures already initialized\n");

#define  zero  0.0
#define  one   1.0
#define  half  0.5
#define  third 1.0/3.0
#define  quart 1.0/4.0

#define  StdVol 1.0

    x_0d = createAndInit(0, 1, 1.0);

    w_0d = createAndInitArray(1, StdVol * 1.0);

    Quadrature::quad_0d[0] = new Quadrature("0d", 0, 0, 1, Quadrature::x_0d, Quadrature::w_0d);

    /****************************************************************************/
    /*  1d quadrature formulas using 2 barycentric coordinates                  */
    /****************************************************************************/

#define MAX_QUAD_DEG_1d   19

    /****************************************************************************/
    /*  quadrature exact on P_1                                                 */
    /****************************************************************************/

    x0_1d = createAndInit(1, 1, 0.5, 0.5);
    w0_1d = createAndInitArray(1, StdVol * 1.0);

    /****************************************************************************/
    /*  quadrature exact on P_3                                                 */
    /****************************************************************************/

    x1_1d = createAndInit(1, 2, 
			  0.788675134594813, 0.211324865405187,
			  0.211324865405187, 0.788675134594813);
    w1_1d = createAndInitArray(2, StdVol * 0.5, StdVol * 0.5);

    /****************************************************************************/
    /*  quadrature exact on P_5                                                 */
    /****************************************************************************/


    x2_1d = createAndInit(1, 3,
			  0.887298334620741, 0.112701665379259,
			  0.500000000000000, 0.500000000000000,
			  0.112701665379259, 0.887298334620741);
    w2_1d = createAndInitArray(3, 
			       StdVol * 0.277777777777778,
			       StdVol * 0.444444444444444,
			       StdVol * 0.277777777777778);

    /****************************************************************************/
    /*  quadrature exact on P_7                                                 */
    /****************************************************************************/

    x3_1d = createAndInit(1, 4,
			  0.930568155797026, 0.069431844202973,
			  0.669990521792428, 0.330009478207572,
			  0.330009478207572, 0.669990521792428,
			  0.069431844202973, 0.930568155797026);
    w3_1d = createAndInitArray(4, 
			       StdVol * 0.173927422568727,
			       StdVol * 0.326072577431273,
			       StdVol * 0.326072577431273,
			       StdVol * 0.173927422568727);

    /****************************************************************************/
    /*  quadrature exact on P_9                                                 */
    /****************************************************************************/

    x4_1d = createAndInit(1, 5,
			  0.953089922969332, 0.046910077030668,
			  0.769234655052841, 0.230765344947159,
			  0.500000000000000, 0.500000000000000,
			  0.230765344947159, 0.769234655052841,
			  0.046910077030668, 0.953089922969332);
    w4_1d = createAndInitArray(5, 
			       StdVol * 0.118463442528095,
			       StdVol * 0.239314335249683,
			       StdVol * 0.284444444444444,
			       StdVol * 0.239314335249683,
			       StdVol * 0.118463442528095);

    /****************************************************************************/
    /*  quadrature exact on P_11                                                */
    /****************************************************************************/

    x5_1d = createAndInit(1, 6,
			  0.966234757101576, 0.033765242898424,
			  0.830604693233133, 0.169395306766867,
			  0.619309593041598, 0.380690406958402,
			  0.380690406958402, 0.619309593041598,
			  0.169395306766867, 0.830604693233133,
			  0.033765242898424, 0.966234757101576);
    w5_1d = createAndInitArray(6, 
			       StdVol * 0.085662246189585,
			       StdVol * 0.180380786524069,
			       StdVol * 0.233956967286345,
			       StdVol * 0.233956967286345,
			       StdVol * 0.180380786524069,
			       StdVol * 0.085662246189585);

    /****************************************************************************/
    /*  quadrature exact on P_13                                                */
    /****************************************************************************/

    x6_1d = createAndInit(1, 7,
			  0.974553956171380, 0.025446043828620,
			  0.870765592799697, 0.129234407200303,
			  0.702922575688699, 0.297077424311301,
			  0.500000000000000, 0.500000000000000,
			  0.297077424311301, 0.702922575688699,
			  0.129234407200303, 0.870765592799697,
			  0.025446043828620, 0.974553956171380);
    w6_1d = createAndInitArray(7, 
			       StdVol * 0.064742483084435,
			       StdVol * 0.139852695744614,
			       StdVol * 0.190915025252559,
			       StdVol * 0.208979591836735,
			       StdVol * 0.190915025252559,
			       StdVol * 0.139852695744614,
			       StdVol * 0.064742483084435);

    /****************************************************************************/
    /*  quadrature exact on P_15                                                */
    /****************************************************************************/

    x7_1d = createAndInit(1, 8, 
			  0.980144928248768, 0.019855071751232,
			  0.898333238706813, 0.101666761293187,
			  0.762766204958164, 0.237233795041836,
			  0.591717321247825, 0.408282678752175,
			  0.408282678752175, 0.591717321247825,
			  0.237233795041836, 0.762766204958164,
			  0.101666761293187, 0.898333238706813,
			  0.019855071751232, 0.980144928248768);
    w7_1d = createAndInitArray(8, 
			       StdVol * 0.050614268145188,
			       StdVol * 0.111190517226687,
			       StdVol * 0.156853322938943,
			       StdVol * 0.181341891689181,
			       StdVol * 0.181341891689181,
			       StdVol * 0.156853322938943,
			       StdVol * 0.111190517226687,
			       StdVol * 0.050614268145188);

    /****************************************************************************/
    /*  quadrature exact on P_17                                                */
    /****************************************************************************/

    x8_1d = createAndInit(1, 9,
			  0.984080119753813, 0.015919880246187,
			  0.918015553663318, 0.081984446336682,
			  0.806685716350295, 0.193314283649705,
			  0.662126711701905, 0.337873288298095,
			  0.500000000000000, 0.500000000000000,
			  0.337873288298095, 0.662126711701905,
			  0.193314283649705, 0.806685716350295,
			  0.081984446336682, 0.918015553663318,
			  0.015919880246187, 0.984080119753813);
    w8_1d = createAndInitArray(9, 
			       StdVol * 0.040637194180787,
			       StdVol * 0.090324080347429,
			       StdVol * 0.130305348201467,
			       StdVol * 0.156173538520001,
			       StdVol * 0.165119677500630,
			       StdVol * 0.156173538520001,
			       StdVol * 0.130305348201467,
			       StdVol * 0.090324080347429,
			       StdVol * 0.040637194180787);

    /****************************************************************************/
    /*  quadrature exact on P_19                                                */
    /****************************************************************************/

    x9_1d = createAndInit(1, 10,
			  0.986953264258586, 0.013046735741414,
			  0.932531683344493, 0.067468316655508,
			  0.839704784149512, 0.160295215850488,
			  0.716697697064623, 0.283302302935377,
			  0.574437169490815, 0.425562830509185,
			  0.425562830509185, 0.574437169490815,
			  0.283302302935377, 0.716697697064623,
			  0.160295215850488, 0.839704784149512,
			  0.067468316655508, 0.932531683344493,
			  0.013046735741414, 0.986953264258586);
    w9_1d = createAndInitArray(10, 
			       StdVol * 0.033335672154344,
			       StdVol * 0.074725674575291,
			       StdVol * 0.109543181257991,
			       StdVol * 0.134633359654998,
			       StdVol * 0.147762112357376,
			       StdVol * 0.147762112357376,
			       StdVol * 0.134633359654998,
			       StdVol * 0.109543181257991,
			       StdVol * 0.074725674575291,
			       StdVol * 0.033335672154344);

    Quadrature::quad_1d[0]= new Quadrature("1d-Gauss: P_1", 1, 1, 1, Quadrature::x0_1d, Quadrature::w0_1d); /* P_0   */
    Quadrature::quad_1d[1]= new Quadrature("1d-Gauss: P_1", 1, 1, 1, Quadrature::x0_1d, Quadrature::w0_1d); /* P_1   */
    Quadrature::quad_1d[2]= new Quadrature("1d-Gauss: P_3", 3, 1, 2, Quadrature::x1_1d, Quadrature::w1_1d); /* P_2   */
    Quadrature::quad_1d[3]= new Quadrature("1d-Gauss: P_3", 3, 1, 2, Quadrature::x1_1d, Quadrature::w1_1d); /* P_3   */
    Quadrature::quad_1d[4]= new Quadrature("1d-Gauss: P_5", 5, 1, 3, Quadrature::x2_1d, Quadrature::w2_1d); /* P_4   */
    Quadrature::quad_1d[5]= new Quadrature("1d-Gauss: P_5", 5, 1, 3, Quadrature::x2_1d, Quadrature::w2_1d); /* P_5   */
    Quadrature::quad_1d[6]= new Quadrature("1d-Gauss: P_7", 7, 1, 4, Quadrature::x3_1d, Quadrature::w3_1d); /* P_6   */
    Quadrature::quad_1d[7]= new Quadrature("1d-Gauss: P_7", 7, 1, 4, Quadrature::x3_1d, Quadrature::w3_1d); /* P_7   */
    Quadrature::quad_1d[8]= new Quadrature("1d-Gauss: P_9", 9, 1, 5, Quadrature::x4_1d, Quadrature::w4_1d); /* P_8   */
    Quadrature::quad_1d[9]= new Quadrature("1d-Gauss: P_9", 9, 1, 5, Quadrature::x4_1d, Quadrature::w4_1d); /* P_9   */
    Quadrature::quad_1d[10]= new Quadrature("1d-Gauss: P_11", 11, 1, 6, Quadrature::x5_1d, Quadrature::w5_1d); /* P_10  */
    Quadrature::quad_1d[11]= new Quadrature("1d-Gauss: P_11", 11, 1, 6, Quadrature::x5_1d, Quadrature::w5_1d); /* P_11  */
    Quadrature::quad_1d[12]= new Quadrature("1d-Gauss: P_13", 13, 1, 7, Quadrature::x6_1d, Quadrature::w6_1d); /* P_12  */
    Quadrature::quad_1d[13]= new Quadrature("1d-Gauss: P_13", 13, 1, 7, Quadrature::x6_1d, Quadrature::w6_1d); /* P_13  */
    Quadrature::quad_1d[14]= new Quadrature("1d-Gauss: P_15", 15, 1, 8, Quadrature::x7_1d, Quadrature::w7_1d); /* P_14  */
    Quadrature::quad_1d[15]= new Quadrature("1d-Gauss: P_15", 15, 1, 8, Quadrature::x7_1d, Quadrature::w7_1d); /* P_15  */
    Quadrature::quad_1d[16]= new Quadrature("1d-Gauss: P_17", 17, 1, 9, Quadrature::x8_1d, Quadrature::w8_1d); /* P_16  */
    Quadrature::quad_1d[17]= new Quadrature("1d-Gauss: P_17", 17, 1, 9, Quadrature::x8_1d, Quadrature::w8_1d); /* P_17  */
    Quadrature::quad_1d[18]= new Quadrature("1d-Gauss: P_19", 19, 1, 10, Quadrature::x9_1d, Quadrature::w9_1d); /* P_18 */
    Quadrature::quad_1d[19]= new Quadrature("1d-Gauss: P_19", 19, 1, 10, Quadrature::x9_1d, Quadrature::w9_1d);/* P_19 */

#undef StdVol
    /****************************************************************************/
    /****************************************************************************/
    /****************************************************************************/
    /*  2d quadrature formulas using 3 barycentric coordinates                  */
    /****************************************************************************/
    /****************************************************************************/
    /****************************************************************************/

#define CYCLE(c1,c2,c3) c1, c2, c3, c2, c3, c1, c3, c1, c2

#define ALL_COMB(c1,c2,c3)  CYCLE(c1,c2,c3), CYCLE(c1,c3,c2)
#define W_CYCLE(w1)         w1, w1, w1
#define W_ALL_COMB(w1)      W_CYCLE(w1), W_CYCLE(w1)

#define MAX_QUAD_DEG_2d   17

#define StdVol 0.5

    /****************************************************************************/
    /*  quadrature exact on P 1                                                 */
    /****************************************************************************/

#define N1  1

#define c1  1.0/3.0
#define w1  StdVol*1.0

    x1_2d = createAndInit(2, 1, c1, c1, c1);
    w1_2d = createAndInitArray(N1, w1);

#undef c1
#undef w1

    /****************************************************************************/
    /*  quadrature exact on P 2                                                 */
    /* Stroud, A.H.: Approximate calculation of multiple integrals              */
    /* Prentice-Hall Series in Automatic Computation. (1971)                    */
    /* optimal number of points: 3, number of points: 3                         */
    /* interior points, completly symmetric in barycentric coordinates          */
    /****************************************************************************/

#define N2  3

#define c1  2.0/3.0
#define c2  1.0/6.0
#define w1  StdVol/3.0

    x2_2d = createAndInit(2, 3, CYCLE(c1, c2, c2));
    w2_2d = createAndInitArray(3, W_CYCLE(w1));

#undef c1
#undef c2
#undef w1

    /****************************************************************************/
    /*  quadrature exact on P_3                                                 */
    /****************************************************************************/

#define N3  6

#define c1  0.0
#define c2  1.0/2.0
#define c3  4.0/6.0
#define c4  1.0/6.0
#define w1  StdVol*1.0/30.0
#define w2  StdVol*9.0/30.0

    x3_2d = createAndInit(2, N3, CYCLE(c1, c2, c2), CYCLE(c3, c4, c4));
    w3_2d = createAndInitArray(N3, W_CYCLE(w1), W_CYCLE(w2));

#undef c1
#undef c2
#undef c3
#undef c4
#undef w1
#undef w2

    /****************************************************************************/
    /*  quadrature exact on P 4                                                 */
    /* Dunavant, D.A.: High degree efficient symmetrical Gaussian quadrature    */
    /* rules for the triangle. Int. J. Numer. Methods Eng. 21, 1129-1148 (1985) */
    /* nearly optimal number of (interior) points, positive wheights  (PI)      */
    /* number of points: 6, optimal number of points: 6                         */
    /****************************************************************************/

#define N4  6

#define c1  0.816847572980459
#define c2  0.091576213509771
#define c3  0.108103018168070
#define c4  0.445948490915965
#define w1  StdVol*0.109951743655322
#define w2  StdVol*0.223381589678011

    x4_2d = createAndInit(2, 6, CYCLE(c1, c2, c2), CYCLE(c3, c4, c4));
    w4_2d = createAndInitArray(6, W_CYCLE(w1), W_CYCLE(w2));

#undef c1
#undef c2
#undef c3
#undef c4
#undef w1
#undef w2


    /****************************************************************************/
    /*  quadrature exact on P 5                                                 */
    /****************************************************************************/

    /****************************************************************************/
    /* Dunavant, D.A.: High degree efficient symmetrical Gaussian quadrature    */
    /* rules for the triangle. Int. J. Numer. Methods Eng. 21, 1129-1148 (1985) */
    /* nealy optimal number of (interior) points, positive wheights  (PI)       */
    /* number of points: 7, optimal number of points: 7                         */
    /****************************************************************************/

#define N5   7

#define c1   1.0/3.0
#define c2   0.797426985353087
#define c3   0.101286507323456
#define c4   0.059715871789770
#define c5   0.470142064105115
#define w1   StdVol*0.225000000000000
#define w2   StdVol*0.125939180544827
#define w3   StdVol*0.132394152788506

    x5_2d = createAndInit(2, 7, c1, c1, c1, CYCLE(c2, c3, c3), CYCLE(c4, c5, c5));
    w5_2d = createAndInitArray(N5, w1, W_CYCLE(w2), W_CYCLE(w3));

#undef c1
#undef c2
#undef c3
#undef c4
#undef c5
#undef w1
#undef w2
#undef w3

    /****************************************************************************/
    /*  quadrature exact on P 6: only 12 point rule available in the literature */
    /*  ->  use quadrature exact on P 7 with 12 points                          */
    /****************************************************************************/

    /****************************************************************************/
    /*  quadrature exact on P 7                                                 */
    /****************************************************************************/

    /****************************************************************************/
    /* Gatermann, Karin: The construction of symmetric cubature formulas for    */
    /* the square and the triangle. Computing 40, No.3, 229-240 (1988)          */
    /* optimal number of points: 12, number of points: 12                       */
    /* only interior points, not completly symmetric in barycentric coordinates */
    /****************************************************************************/

#define N7   12

#define c1   0.06238226509439084
#define c2   0.06751786707392436
#define c3   0.87009986783168480
#define c4   0.05522545665692000
#define c5   0.32150249385201560
#define c6   0.62327204949106440
#define c7   0.03432430294509488
#define c8   0.66094919618679800
#define c9   0.30472650086810720
#define c10  0.5158423343536001
#define c11  0.2777161669764050
#define c12  0.2064414986699949
#define w1   0.02651702815743450
#define w2   0.04388140871444811
#define w3   0.02877504278497528
#define w4   0.06749318700980879

    x7_2d = createAndInit(2, 12,
			  CYCLE(c1, c2, c3), CYCLE(c4, c5, c6),
			  CYCLE(c7, c8, c9), CYCLE(c10, c11, c12));
    w7_2d = createAndInitArray(N7, W_CYCLE(w1), W_CYCLE(w2),
			       W_CYCLE(w3), W_CYCLE(w4));

#undef c1
#undef c2
#undef c3
#undef c4
#undef c5
#undef c6
#undef c7
#undef c8
#undef c9
#undef c10
#undef c11
#undef c12
#undef w1
#undef w2
#undef w3
#undef w4


    /****************************************************************************/
    /*  quadrature exact on P 8                                                 */
    /* Dunavant, D.A.: High degree efficient symmetrical Gaussian quadrature    */
    /* rules for the triangle. Int. J. Numer. Methods Eng. 21, 1129-1148 (1985) */
    /* nealy optimal number of (interior) points, positive wheights  (PI)       */
    /* number of points: 16, optimal number of points: 15                       */
    /* only interior points, completly symmetric in barycentric coordinates     */
    /****************************************************************************/

#define N8  16

#define c1   1.0/3.0
#define c2   0.081414823414554
#define c3   0.459292588292723
#define c4   0.658861384496480
#define c5   0.170569307751760
#define c6   0.898905543365938
#define c7   0.050547228317031
#define c8   0.008394777409958
#define c9   0.263112829634638
#define c10  0.728492392955404
#define w1   StdVol*0.144315607677787
#define w2   StdVol*0.095091634267285
#define w3   StdVol*0.103217370534718
#define w4   StdVol*0.032458497623198
#define w5   StdVol*0.027230314174435

    x8_2d = createAndInit(2, 16,
			  c1, c1, c1, 
			  CYCLE(c2, c3, c3),
			  CYCLE(c4, c5, c5),
			  CYCLE(c6, c7, c7),
			  ALL_COMB(c8, c9, c10));
    w8_2d = createAndInitArray(N8, w1, W_CYCLE(w2), W_CYCLE(w3),
			       W_CYCLE(w4), W_ALL_COMB(w5));
				  

#undef c1
#undef c2
#undef c3
#undef c4
#undef c5
#undef c6
#undef c7
#undef c8
#undef c9
#undef c10
#undef w1
#undef w2
#undef w3
#undef w4
#undef w5

    /****************************************************************************/
    /*  quadrature exact on P 9                                                 */
    /* Dunavant, D.A.: High degree efficient symmetrical Gaussian quadrature    */
    /* rules for the triangle. Int. J. Numer. Methods Eng. 21, 1129-1148 (1985) */
    /* nealy optimal number of (interior) points, positive wheights  (PI)       */
    /* optimal number of points: ?, number of points: 19                        */
    /* only interior points, completly symmetric in barycentric coordinates     */
    /****************************************************************************/

#define N9  19

#define c1   1.0/3.0
#define c2   0.020634961602525
#define c3   0.489682519198738
#define c4   0.125820817014127
#define c5   0.437089591492937
#define c6   0.623592928761935
#define c7   0.188203535619033
#define c8   0.910540973211095
#define c9   0.044729513394453
#define c10  0.036838412054736
#define c11  0.221962989160766
#define c12  0.741198598784498
#define w1   StdVol*0.097135796282799
#define w2   StdVol*0.031334700227139
#define w3   StdVol*0.077827541004774
#define w4   StdVol*0.079647738927210
#define w5   StdVol*0.025577675658698
#define w6   StdVol*0.043283539377289

    x9_2d = createAndInit(2, 19,
			  c1, c1, c1,
			  CYCLE(c2, c3, c3),
			  CYCLE(c4, c5, c5),
			  CYCLE(c6, c7, c7),
			  CYCLE(c8, c9, c9),
			  ALL_COMB(c10, c11, c12));
    w9_2d = createAndInitArray(N9, w1, 
			       W_CYCLE(w2),
			       W_CYCLE(w3),
			       W_CYCLE(w4),
			       W_CYCLE(w5),
			       W_ALL_COMB(w6));

#undef c1
#undef c2
#undef c3
#undef c4
#undef c5
#undef c6
#undef c7
#undef c8
#undef c9
#undef c10
#undef c11
#undef c12
#undef w1
#undef w2
#undef w3
#undef w4
#undef w5
#undef w6

    /****************************************************************************/
    /*  quadrature exact on P 10                                                */
    /* Dunavant, D.A.: High degree efficient symmetrical Gaussian quadrature    */
    /* rules for the triangle. Int. J. Numer. Methods Eng. 21, 1129-1148 (1985) */
    /* nealy optimal number of (interior) points, positive wheights  (PI)       */
    /* optimal number of points: ?, number of points: 25                        */
    /* only interior points, completly symmetric in barycentric coordinates     */
    /****************************************************************************/

#define N10 25

#define c1   1.0/3.0

#define c2   0.028844733232685
#define c3   0.485577633383657

#define c4   0.781036849029926
#define c5   0.109481575485037

#define c6   0.141707219414880
#define c7   0.307939838764121
#define c8   0.550352941820999

#define c9   0.025003534762686
#define c10  0.246672560639903
#define c11  0.728323904597411

#define c12  0.009540815400299
#define c13  0.066803251012200
#define c14  0.923655933587500

#define w1   StdVol*0.090817990382754
#define w2   StdVol*0.036725957756467
#define w3   StdVol*0.045321059435528
#define w4   StdVol*0.072757916845420
#define w5   StdVol*0.028327242531057
#define w6   StdVol*0.009421666963733

    x10_2d = createAndInit(2, 25,
			   c1, c1, c1,
			   CYCLE(c2, c3, c3),
			   CYCLE(c4, c5, c5),
			   ALL_COMB(c6, c7, c8),
			   ALL_COMB(c9, c10, c11),
			   ALL_COMB(c12, c13, c14));
    w10_2d = createAndInitArray(N10, w1, 
				W_CYCLE(w2),
				W_CYCLE(w3),
				W_ALL_COMB(w4),
				W_ALL_COMB(w5),
				W_ALL_COMB(w6));

#undef c1
#undef c2
#undef c3
#undef c4
#undef c5
#undef c6
#undef c7
#undef c8
#undef c9
#undef c10
#undef c11
#undef c12
#undef c13
#undef c14
#undef w1
#undef w2
#undef w3
#undef w4
#undef w5
#undef w6

    /****************************************************************************/
    /*  quadrature exact on P 11                                                */
    /* Dunavant, D.A.: High degree efficient symmetrical Gaussian quadrature    */
    /* rules for the triangle. Int. J. Numer. Methods Eng. 21, 1129-1148 (1985) */
    /* nealy optimal number of (interior) points, positive wheights  (PI)       */
    /* optimal number of points: ?, number of points: 27                        */
    /* only interior points, completly symmetric in barycentric coordinates     */
    /****************************************************************************/

#define N11 27

#define c1  -0.069222096541517
#define c2   0.534611048270758

#define c3   0.202061394068290
#define c4   0.398969302965855

#define c5   0.593380199137435
#define c6   0.203309900431282

#define c7   0.761298175434837
#define c8   0.119350912282581

#define c9   0.935270103777448
#define c10  0.032364948111276

#define c11  0.050178138310495
#define c12  0.356620648261293
#define c13  0.593201213428213

#define c14  0.021022016536166
#define c15  0.171488980304042
#define c16  0.807489003159792

#define w1   StdVol*0.000927006328961
#define w2   StdVol*0.077149534914813
#define w3   StdVol*0.059322977380774
#define w4   StdVol*0.036184540503418
#define w5   StdVol*0.013659731002678
#define w6   StdVol*0.052337111962204
#define w7   StdVol*0.020707659639141

    x11_2d = createAndInit(2, 27,
			   CYCLE(c1,c2,c2),
			   CYCLE(c3,c4,c4),
			   CYCLE(c5,c6,c6),
			   CYCLE(c7,c8,c8),
			   CYCLE(c9,c10,c10),
			   ALL_COMB(c11,c12,c13),
			   ALL_COMB(c14,c15,c16));
    w11_2d = createAndInitArray(N11, W_CYCLE(w1),
				W_CYCLE(w2),
				W_CYCLE(w3),
				W_CYCLE(w4),
				W_CYCLE(w5),
				W_ALL_COMB(w6),
				W_ALL_COMB(w7));

#undef c1
#undef c2
#undef c3
#undef c4
#undef c5
#undef c6
#undef c7
#undef c8
#undef c9
#undef c10
#undef c11
#undef c12
#undef c13
#undef c14
#undef c15
#undef c16
#undef w1
#undef w2
#undef w3
#undef w4
#undef w5
#undef w6
#undef w7

    /****************************************************************************/
    /*  quadrature exact on P 12                                                */
    /* Dunavant, D.A.: High degree efficient symmetrical Gaussian quadrature    */
    /* rules for the triangle. Int. J. Numer. Methods Eng. 21, 1129-1148 (1985) */
    /* nealy optimal number of (interior) points, positive wheights  (PI)       */
    /* optimal number of points: 2, number of points: 25                        */
    /* only interior points, completly symmetric in barycentric coordinates     */
    /****************************************************************************/

#define N12 33

#define c1   0.023565220452390
#define c2   0.488217389773805

#define c3   0.120551215411079
#define c4   0.439724392294460

#define c5   0.457579229975768
#define c6   0.271210385012116

#define c7   0.744847708916828
#define c8   0.127576145541586
	       
#define c9   0.957365299093579
#define c10  0.021317350453210

#define c11  0.115343494534698
#define c12  0.275713269685514
#define c13  0.608943235779788

#define c14  0.022838332222257
#define c15  0.281325580989940
#define c16  0.695836086787803

#define c17  0.025734050548330
#define c18  0.116251915907597
#define c19  0.858014033544073

#define w1   StdVol*0.025731066440455
#define w2   StdVol*0.043692544538038
#define w3   StdVol*0.062858224217885
#define w4   StdVol*0.034796112930709
#define w5   StdVol*0.006166261051559
#define w6   StdVol*0.040371557766381
#define w7   StdVol*0.022356773202303
#define w8   StdVol*0.017316231108659

    x12_2d = createAndInit(2, 33,
			   CYCLE(c1,c2,c2),
			   CYCLE(c3,c4,c4),
			   CYCLE(c5,c6,c6),
			   CYCLE(c7,c8,c8),
			   CYCLE(c9,c10,c10),
			   ALL_COMB(c11,c12,c13),
			   ALL_COMB(c14,c15,c16),
			   ALL_COMB(c17,c18,c19));
    w12_2d = createAndInitArray(N12, W_CYCLE(w1),
				W_CYCLE(w2),
				W_CYCLE(w3),
				W_CYCLE(w4),
				W_CYCLE(w5),
				W_ALL_COMB(w6),
				W_ALL_COMB(w7),
				W_ALL_COMB(w8));

#undef c1
#undef c2
#undef c3
#undef c4
#undef c5
#undef c6
#undef c7
#undef c8
#undef c9
#undef c10
#undef c11
#undef c12
#undef c13
#undef c14
#undef c15
#undef c16
#undef c17
#undef c18
#undef c19
#undef w1
#undef w2
#undef w3
#undef w4
#undef w5
#undef w6
#undef w7
#undef w8

    /****************************************************************************/
    /*  quadrature exact on P 17                                                */
    /* Dunavant, D.A.: High degree efficient symmetrical Gaussian quadrature    */
    /* rules for the triangle. Int. J. Numer. Methods Eng. 21, 1129-1148 (1985) */
    /* nealy optimal number of (interior) points, positive wheights  (PI)       */
    /* optimal number of points: ?, number of points: 61                        */
    /* only interior points, completly symmetric in barycentric coordinates     */
    /****************************************************************************/

#define N17 61

#define c1   1.0/3.0
	       
#define c2   0.005658918886452
#define c3   0.497170540556774
	       
#define c4   0.035647354750751
#define c5   0.482176322624625
	       
#define c6   0.099520061958437
#define c7   0.450239969020782

#define c8   0.199467521245206
#define c9   0.400266239377397

#define c10  0.495717464058095
#define c11  0.252141267970953
	       
#define c12  0.675905990683077
#define c13  0.162047004658461

#define c14  0.848248235478508
#define c15  0.075875882260746

#define c16  0.968690546064356
#define c17  0.015654726967822

#define c18  0.010186928826919
#define c19  0.334319867363658
#define c20  0.655493203809423

#define c21  0.135440871671036
#define c22  0.292221537796944
#define c23  0.572337590532020

#define c24  0.054423924290583
#define c25  0.319574885423190
#define c26  0.626001190286228

#define c27  0.012868560833637
#define c28  0.190704224192292
#define c29  0.796427214974071

#define c30  0.067165782413524
#define c31  0.180483211648746
#define c32  0.752351005937729

#define c33  0.014663182224828
#define c34  0.080711313679564
#define c35  0.904625504095608

#define w1   StdVol*0.033437199290803
#define w2   StdVol*0.005093415440507
#define w3   StdVol*0.014670864527638
#define w4   StdVol*0.024350878353672
#define w5   StdVol*0.031107550868969
#define w6   StdVol*0.031257111218620
#define w7   StdVol*0.024815654339665
#define w8   StdVol*0.014056073070557
#define w9   StdVol*0.003194676173779
#define w10  StdVol*0.008119655318993
#define w11  StdVol*0.026805742283163
#define w12  StdVol*0.018459993210822
#define w13  StdVol*0.008476868534328
#define w14  StdVol*0.018292796770025
#define w15  StdVol*0.006665632004165

    x17_2d = createAndInit(2, 61,
			   c1, c1, c1,
			   CYCLE(c2,c3,c3),
			   CYCLE(c4,c5,c5),
			   CYCLE(c6,c7,c7),
			   CYCLE(c8,c9,c9),
			   CYCLE(c10,c11,c11),
			   CYCLE(c12,c13,c13),
			   CYCLE(c14,c15,c15),
			   CYCLE(c16,c17,c17),
			   ALL_COMB(c18,c19,c20),
			   ALL_COMB(c21,c22,c23),
			   ALL_COMB(c24,c25,c26),
			   ALL_COMB(c27,c28,c29),
			   ALL_COMB(c30,c31,c32),
			   ALL_COMB(c33,c34,c35));
    w17_2d = createAndInitArray(N17, w1, 
				W_CYCLE(w2),
				W_CYCLE(w3),
				W_CYCLE(w4),
				W_CYCLE(w5),
				W_CYCLE(w6),
				W_CYCLE(w7),
				W_CYCLE(w8),
				W_CYCLE(w9),
				W_ALL_COMB(w10),
				W_ALL_COMB(w11),
				W_ALL_COMB(w12),
				W_ALL_COMB(w13),
				W_ALL_COMB(w14),
				W_ALL_COMB(w15));

#undef c1
#undef c2
#undef c3
#undef c4
#undef c5
#undef c6
#undef c7
#undef c8
#undef c9
#undef c10
#undef c11
#undef c12
#undef c13
#undef c14
#undef c15
#undef c16
#undef c17
#undef c18
#undef c19
#undef c20
#undef c21
#undef c22
#undef c23
#undef c24
#undef c25
#undef c26
#undef c27
#undef c28
#undef c29
#undef c30
#undef c31
#undef c32
#undef c33
#undef c34
#undef c35

#undef w1
#undef w2
#undef w3
#undef w4
#undef w5
#undef w6
#undef w7
#undef w8
#undef w9
#undef w10
#undef w11
#undef w12
#undef w13
#undef w14
#undef w15

    Quadrature::quad_2d[0] = new Quadrature("2d-P_1", 1, 2, N1, Quadrature::x1_2d, Quadrature::w1_2d);   /* P 0  */
    Quadrature::quad_2d[1] = new Quadrature("2d-P_1", 1, 2, N1, Quadrature::x1_2d, Quadrature::w1_2d);   /* P 1  */
    Quadrature::quad_2d[2] = new Quadrature("2d  Stroud: P_2", 2, 2, N2, Quadrature::x2_2d, Quadrature::w2_2d);   /* P 2  */
    Quadrature::quad_2d[3] = new Quadrature("2d  Stroud: P_3", 3, 2, N3, Quadrature::x3_2d, Quadrature::w3_2d);   /* P 3  */
    Quadrature::quad_2d[4] = new Quadrature("2d  Dunavant: P_4", 4, 2, N4, Quadrature::x4_2d, Quadrature::w4_2d);   /* P 4  */
    Quadrature::quad_2d[5] = new Quadrature("2d  Dunavant: P_5", 5, 2, N5, Quadrature::x5_2d, Quadrature::w5_2d);   /* P 5  */
    Quadrature::quad_2d[6] = new Quadrature("2d  Gattermann: P_7", 7, 2, N7, Quadrature::x7_2d, Quadrature::w7_2d);   /* P 6  */
    Quadrature::quad_2d[7] = new Quadrature("2d  Gattermann: P_7", 7, 2, N7, Quadrature::x7_2d, Quadrature::w7_2d);   /* P 7  */
    Quadrature::quad_2d[8] = new Quadrature("2d  Dunavant: P_8", 8, 2, N8, Quadrature::x8_2d, Quadrature::w8_2d);   /* P 8  */
    Quadrature::quad_2d[9] = new Quadrature("2d  Dunavant: P_9", 9, 2, N9, Quadrature::x9_2d, Quadrature::w9_2d);   /* P 9  */
    Quadrature::quad_2d[10] = new Quadrature("2d  Dunavant: P_10", 10, 2, N10, Quadrature::x10_2d, Quadrature::w10_2d);/* P 10 */
    Quadrature::quad_2d[11] = new Quadrature("2d  Dunavant: P_11", 11, 2, N11, Quadrature::x11_2d, Quadrature::w11_2d);/* P 11 */
    Quadrature::quad_2d[12] = new Quadrature("2d  Dunavant: P_12", 12, 2, N12, Quadrature::x12_2d, Quadrature::w12_2d);/* P 12 */
    Quadrature::quad_2d[13] = new Quadrature("2d  Dunavant: P_17", 17, 2, N17, Quadrature::x17_2d, Quadrature::w17_2d);/* P 13 */
    Quadrature::quad_2d[14] = new Quadrature("2d  Dunavant: P_17", 17, 2, N17, Quadrature::x17_2d, Quadrature::w17_2d);/* P 14 */
    Quadrature::quad_2d[15] = new Quadrature("2d  Dunavant: P_17", 17, 2, N17, Quadrature::x17_2d, Quadrature::w17_2d);/* P 15 */
    Quadrature::quad_2d[16] = new Quadrature("2d  Dunavant: P_17", 17, 2, N17, Quadrature::x17_2d, Quadrature::w17_2d);/* P 16 */
    Quadrature::quad_2d[17] = new Quadrature("2d  Dunavant: P_17", 17, 2, N17, Quadrature::x17_2d, Quadrature::w17_2d); /* P 17 */


#undef StdVol
#undef N1
#undef N2
#undef N3
#undef N4
#undef N5
#undef N6
#undef N7
#undef N8
#undef N9
#undef N10
#undef N11
#undef N12
#undef N17

    /****************************************************************************/
    /*  3d quadrature formulas using 4 barycentric coordinates                  */
    /****************************************************************************/

#define MAX_QUAD_DEG_3d   7

#define StdVol (1.0/6.0)

    /****************************************************************************/
    /*  quadrature exact on P_1                                                 */
    /****************************************************************************/

    x1_3d = createAndInit(3, 1,
			  quart, quart, quart, quart);
    w1_3d = createAndInitArray(1, StdVol*one);



    /****************************************************************************/
    /*  Quad quadrature exact on P_2                                           */
    /****************************************************************************/

#define c14   0.585410196624969
#define c15   0.138196601125011

    x2_3d = createAndInit(3, 4,
			  c14, c15, c15, c15,
			  c15, c14, c15, c15,
			  c15, c15, c14, c15,
			  c15, c15, c15, c14);
    w2_3d = createAndInitArray(4, StdVol*quart, StdVol*quart,
			       StdVol*quart, StdVol*quart);

    /****************************************************************************/
    /*  quadrature exact on P_3                                                 */
    /****************************************************************************/

#define w8  1.0/40.0
#define w9  9.0/40.0

    x3_3d = createAndInit(3, 8,
			  one,  zero,  zero,  zero,
			  zero,   one,  zero,  zero,
			  zero,  zero,   one,  zero,
			  zero,  zero,  zero,   one,
			  zero, third, third, third,
			  third, zero, third, third,
			  third, third, zero, third,
			  third, third, third, zero);
    w3_3d = createAndInitArray(8,  StdVol*w8, StdVol*w8, StdVol*w8, StdVol*w8,
			       StdVol*w9, StdVol*w9, StdVol*w9, StdVol*w9);

    /****************************************************************************/
    /*  quadrature exact on P_4                                                 */
    /****************************************************************************/

#define c18  0.091971078
#define c19  0.724086765
#define c20  0.319793627
#define c21  0.040619116
#define c22  0.056350832
#define c23  0.443649167

#define w10  0.118518515999999990
#define w11  0.071937078000000002
#define w12  0.069068201999999995
#define w13  0.052910052911999995

    x4_3d = createAndInit(3, 15,
			  quart, quart, quart, quart,
			  c18, c18, c18, c19,
			  c18, c18, c19, c18,
			  c18, c19, c18, c18,
			  c19, c18, c18, c18,
			  c20, c20, c20, c21,
			  c20, c20, c21, c20,
			  c20, c21, c20, c20,
			  c21, c20, c20, c20,
			  c22, c22, c23, c23,
			  c22, c23, c22, c23,
			  c23, c22, c22, c23,
			  c22, c23, c23, c22,
			  c23, c22, c23, c22,
			  c23, c23, c22, c22);
    w4_3d = createAndInitArray(15, StdVol*w10, 
			       StdVol*w11, StdVol*w11,
			       StdVol*w11, StdVol*w11, 
			       StdVol*w12, StdVol*w12,
			       StdVol*w12, StdVol*w12,
			       StdVol*w13, StdVol*w13, StdVol*w13,
			       StdVol*w13, StdVol*w13, StdVol*w13);

    /****************************************************************************/
    /*  quadrature exact on P_5                                                 */
    /****************************************************************************/

    x5_3d = createAndInit(3, 15,
			  0.250000000000000, 0.250000000000000, 0.250000000000000, 0.250000000000000,
			  0.091971078052723, 0.091971078052723, 0.091971078052723, 0.724086765841831,
			  0.724086765841831, 0.091971078052723, 0.091971078052723, 0.091971078052723,
			  0.091971078052723, 0.724086765841831, 0.091971078052723, 0.091971078052723,
			  0.091971078052723, 0.091971078052723, 0.724086765841831, 0.091971078052723,
			  0.319793627829630, 0.319793627829630, 0.319793627829630, 0.040619116511110,
			  0.040619116511110, 0.319793627829630, 0.319793627829630, 0.319793627829630,
			  0.319793627829630, 0.040619116511110, 0.319793627829630, 0.319793627829630,
			  0.319793627829630, 0.319793627829630, 0.040619116511110, 0.319793627829630,
			  0.443649167310371, 0.056350832689629, 0.056350832689629, 0.443649167310371,
			  0.056350832689629, 0.443649167310371, 0.056350832689629, 0.443649167310371,
			  0.056350832689629, 0.056350832689629, 0.443649167310371, 0.443649167310371,
			  0.443649167310371, 0.056350832689629, 0.443649167310371, 0.056350832689629,
			  0.443649167310371, 0.443649167310371, 0.056350832689629, 0.056350832689629,
			  0.056350832689629, 0.443649167310371, 0.443649167310371, 0.056350832689629);
    w5_3d = createAndInitArray(15, StdVol*0.118518518518519,
			       StdVol*0.071937083779019,
			       StdVol*0.071937083779019,
			       StdVol*0.071937083779019,
			       StdVol*0.071937083779019,
			       StdVol*0.069068207226272,
			       StdVol*0.069068207226272,
			       StdVol*0.069068207226272,
			       StdVol*0.069068207226272,
			       StdVol*0.052910052910053,
			       StdVol*0.052910052910053,
			       StdVol*0.052910052910053,
			       StdVol*0.052910052910053,
			       StdVol*0.052910052910053,
			       StdVol*0.052910052910053);

    /****************************************************************************/
    /*  quadrature exact on P_7                                                 */
    /****************************************************************************/

    x7_3d = createAndInit(3, 64,
			  0.0485005494, 0.0543346112, 0.0622918076, 0.8348730318,
			  0.0485005494, 0.0543346112, 0.2960729005, 0.6010919389,
			  0.0485005494, 0.0543346112, 0.6010919389, 0.2960729005,
			  0.0485005494, 0.0543346112, 0.8348730300, 0.0622918093,
			  0.0485005494, 0.2634159753, 0.0477749033, 0.6403085720,
			  0.0485005494, 0.2634159753, 0.2270740686, 0.4610094066,
			  0.0485005494, 0.2634159753, 0.4610094066, 0.2270740686,
			  0.0485005494, 0.2634159753, 0.6403085706, 0.0477749047,
			  0.0485005494, 0.5552859758, 0.0275098315, 0.3687036433,
			  0.0485005494, 0.5552859758, 0.1307542021, 0.2654592727,
			  0.0485005494, 0.5552859758, 0.2654592727, 0.1307542021,
			  0.0485005494, 0.5552859758, 0.3687036425, 0.0275098323,
			  0.0485005494, 0.8185180165, 0.0092331459, 0.1237482881,
			  0.0485005494, 0.8185180165, 0.0438851337, 0.0890963004,
			  0.0485005494, 0.8185180165, 0.0890963004, 0.0438851337,
			  0.0485005494, 0.8185180165, 0.1237482879, 0.0092331462,
			  0.2386007376, 0.0434790928, 0.0498465199, 0.6680736497,
			  0.2386007376, 0.0434790928, 0.2369204606, 0.4809997090,
			  0.2386007376, 0.0434790928, 0.4809997090, 0.2369204606,
			  0.2386007376, 0.0434790928, 0.6680736482, 0.0498465214,
			  0.2386007376, 0.2107880664, 0.0382299497, 0.5123812464,
			  0.2386007376, 0.2107880664, 0.1817069135, 0.3689042825,
			  0.2386007376, 0.2107880664, 0.3689042825, 0.1817069135,
			  0.2386007376, 0.2107880664, 0.5123812453, 0.0382299508,
			  0.2386007376, 0.4443453248, 0.0220136390, 0.2950402987,
			  0.2386007376, 0.4443453248, 0.1046308045, 0.2124231331,
			  0.2386007376, 0.4443453248, 0.2124231331, 0.1046308045,
			  0.2386007376, 0.4443453248, 0.2950402980, 0.0220136396,
			  0.2386007376, 0.6549862048, 0.0073884546, 0.0990246030,
			  0.2386007376, 0.6549862048, 0.0351173176, 0.0712957400,
			  0.2386007376, 0.6549862048, 0.0712957400, 0.0351173176,
			  0.2386007376, 0.6549862048, 0.0990246028, 0.0073884548,
			  0.5170472951, 0.0275786260, 0.0316174612, 0.4237566177,
			  0.5170472951, 0.0275786260, 0.1502777622, 0.3050963168,
			  0.5170472951, 0.0275786260, 0.3050963168, 0.1502777622,
			  0.5170472951, 0.0275786260, 0.4237566168, 0.0316174621,
			  0.5170472951, 0.1337020823, 0.0242491141, 0.3250015085,
			  0.5170472951, 0.1337020823, 0.1152560157, 0.2339946069,
			  0.5170472951, 0.1337020823, 0.2339946069, 0.1152560157,
			  0.5170472951, 0.1337020823, 0.3250015078, 0.0242491148,
			  0.5170472951, 0.2818465779, 0.0139631689, 0.1871429581,
			  0.5170472951, 0.2818465779, 0.0663669280, 0.1347391990,
			  0.5170472951, 0.2818465779, 0.1347391990, 0.0663669280,
			  0.5170472951, 0.2818465779, 0.1871429577, 0.0139631693,
			  0.5170472951, 0.4154553004, 0.0046864691, 0.0628109354,
			  0.5170472951, 0.4154553004, 0.0222747832, 0.0452226213,
			  0.5170472951, 0.4154553004, 0.0452226213, 0.0222747832,
			  0.5170472951, 0.4154553004, 0.0628109352, 0.0046864693,
			  0.7958514179, 0.0116577407, 0.0133649937, 0.1791258477,
			  0.7958514179, 0.0116577407, 0.0635238021, 0.1289670393,
			  0.7958514179, 0.0116577407, 0.1289670393, 0.0635238021,
			  0.7958514179, 0.0116577407, 0.1791258473, 0.0133649941,
			  0.7958514179, 0.0565171087, 0.0102503252, 0.1373811482,
			  0.7958514179, 0.0565171087, 0.0487197855, 0.0989116879,
			  0.7958514179, 0.0565171087, 0.0989116879, 0.0487197855,
			  0.7958514179, 0.0565171087, 0.1373811479, 0.0102503255,
			  0.7958514179, 0.1191391593, 0.0059023608, 0.0791070620,
			  0.7958514179, 0.1191391593, 0.0280539153, 0.0569555075,
			  0.7958514179, 0.1191391593, 0.0569555075, 0.0280539153,
			  0.7958514179, 0.1191391593, 0.0791070618, 0.0059023610,
			  0.7958514179, 0.1756168040, 0.0019810139, 0.0265507642,
			  0.7958514179, 0.1756168040, 0.0094157572, 0.0191160209,
			  0.7958514179, 0.1756168040, 0.0191160209, 0.0094157572,
			  0.7958514179, 0.1756168040, 0.0265507642, 0.0019810140);
    w7_3d = createAndInitArray(64, StdVol*0.0156807540, StdVol*0.0293976870,
			       StdVol*0.0293976870, StdVol*0.0156807540,
			       StdVol*0.0235447608, StdVol*0.0441408300,
			       StdVol*0.0441408300, StdVol*0.0235447608,
			       StdVol*0.0150258564, StdVol*0.0281699100,
			       StdVol*0.0281699100, StdVol*0.0150258564,
			       StdVol*0.0036082374, StdVol*0.0067645878,
			       StdVol*0.0067645878, StdVol*0.0036082374,
			       StdVol*0.0202865376, StdVol*0.0380324358,
			       StdVol*0.0380324358, StdVol*0.0202865376,
			       StdVol*0.0304603764, StdVol*0.0571059660,
			       StdVol*0.0571059660, StdVol*0.0304603764,
			       StdVol*0.0194392824, StdVol*0.0364440336,
			       StdVol*0.0364440336, StdVol*0.0194392824,
			       StdVol*0.0046680564, StdVol*0.0087514968,
			       StdVol*0.0087514968, StdVol*0.0046680564,
			       StdVol*0.0097055322, StdVol*0.0181955664,
			       StdVol*0.0181955664, StdVol*0.0097055322,
			       StdVol*0.0145729242, StdVol*0.0273207684,
			       StdVol*0.0273207684, StdVol*0.0145729242,
			       StdVol*0.0093001866, StdVol*0.0174356394,
			       StdVol*0.0174356394, StdVol*0.0093001866,
			       StdVol*0.0022333026, StdVol*0.0041869110,
			       StdVol*0.0041869110, StdVol*0.0022333026,
			       StdVol*0.0014639124, StdVol*0.0027444882,
			       StdVol*0.0027444882, StdVol*0.0014639124,
			       StdVol*0.0021980748, StdVol*0.0041208678,
			       StdVol*0.0041208678, StdVol*0.0021980748,
			       StdVol*0.0014027730, StdVol*0.0026298660,
			       StdVol*0.0026298660, StdVol*0.0014027730,
			       StdVol*0.0003368550, StdVol*0.0006315234,
			       StdVol*0.0006315234, StdVol*0.0003368550);

    /****************************************************************************/
    /*  build a vector of Quad' quadrature formulars. For quadrature of degree */
    /*  use that of degree (only on function evaluation also)                   */
    /****************************************************************************/

    Quadrature::quad_3d[0] = new Quadrature("3d Stroud: P_1", 1, 3,  1, Quadrature::x1_3d, Quadrature::w1_3d);   /* P_0  */
    Quadrature::quad_3d[1] = new Quadrature("3d Stroud: P_1", 1, 3,  1, Quadrature::x1_3d, Quadrature::w1_3d);   /* P_1  */
    Quadrature::quad_3d[2] = new Quadrature("3d Stroud: P_2", 2, 3,  4, Quadrature::x2_3d, Quadrature::w2_3d);   /* P_2  */
    Quadrature::quad_3d[3] = new Quadrature("3d Stroud: P_3", 3, 3,  8, Quadrature::x3_3d, Quadrature::w3_3d);   /* P_3  */
    Quadrature::quad_3d[4] = new Quadrature("3d ???: P_5", 5, 3, 15, Quadrature::x5_3d, Quadrature::w5_3d);   /* P_4  */
    Quadrature::quad_3d[5] = new Quadrature("3d ???: P_5", 5, 3, 15, Quadrature::x5_3d, Quadrature::w5_3d);   /* P_5  */
    Quadrature::quad_3d[6] = new Quadrature("3d ???: P_7", 7, 3, 64, Quadrature::x7_3d, Quadrature::w7_3d);   /* P_6  */
    Quadrature::quad_3d[7] = new Quadrature("3d ???: P_7", 7, 3, 64, Quadrature::x7_3d, Quadrature::w7_3d);   /* P_7  */


#undef StdVol

    /****************************************************************************/
    /*  integration in different dimensions                                     */
    /****************************************************************************/

    Quadrature::quad_nd[0] = Quadrature::quad_0d;
    Quadrature::quad_nd[1] = Quadrature::quad_1d;
    Quadrature::quad_nd[2] = Quadrature::quad_2d;
    Quadrature::quad_nd[3] = Quadrature::quad_3d;
  }


  Quadrature* Quadrature::provideQuadrature(int dim_, int degree_)
  {
    FUNCNAME("Quadrature::provideQuadrature()");

    switch (dim_) {
    case 0:
      degree_ = 0;
      break;
    case 1:
      degree_ = std::min(degree_, MAX_QUAD_DEG_1d);
      break;
    case 2:
      degree_ = std::min(degree_, MAX_QUAD_DEG_2d);
      break;
    case 3:
      degree_ = std::min(degree_, MAX_QUAD_DEG_3d);
      break;
    default:
      ERROR_EXIT("invalid dim\n");
    }
    
    if (x0_1d == NULL) 
      initStaticQuadratures();

    return (quad_nd[dim_][degree_]);
  }


  Quadrature::~Quadrature()
  {
    if (lambda)
      delete lambda;
    
    delete [] w;
  }


  double Quadrature::integrateStdSimplex(std::function<double(DimVec<double>)> f)
  {
    double result = 0.0;
    // calculate weighted sum over all quadrature-points
    for (int i = 0; i < n_points; i++)
      result += w[i] * f((*lambda)[i]);

    return result;
  }


  FastQuadrature* FastQuadrature::provideFastQuadrature(const BasisFunction* bas_fcts,
                                          							const Quadrature& quad, 
                                          							Flag init_flag)
  {
    FastQuadrature *quad_fast = NULL;

// #pragma omp critical
    {
      list<FastQuadrature*>::iterator fast = fastQuadList.begin(); 
      for (; fast != fastQuadList.end(); fast++)
      	if ((*fast)->basisFunctions == bas_fcts && 
      	    (*fast)->quadrature == &quad)  
      	  break;
      
      if (fast != fastQuadList.end() && 
      	  ((*fast)->init_flag & init_flag) == init_flag) {
      	quad_fast = *fast;
      } else {
      	if (fast == fastQuadList.end()) {
      	  quad_fast = 
      	    new FastQuadrature(const_cast<BasisFunction*>(bas_fcts), 
      			       const_cast<Quadrature*>(&quad), 0);
	  
      	  fastQuadList.push_front(quad_fast);
      	  
      	  max_points = std::max(max_points, quad.getNumPoints());
      	} else {
      	  quad_fast = (*fast);
      	}
      }
      
      quad_fast->init(init_flag);  
    }

    return quad_fast;
  }


  void FastQuadrature::init(Flag init_flag)
  {
    int dim = quadrature->getDim();
    int nPoints = quadrature->getNumPoints();
    int nBasFcts = basisFunctions->getNumber();

    DimVec<double> lambda(dim);

    // ----- initialize phi ---------------------------------------------

    if (num_rows(phi) == 0 && init_flag.isSet(INIT_PHI)) {
      phi.change_dim(nPoints, nBasFcts);

      // fill memory
      for (int i = 0; i< nPoints; i++) {
      	lambda = quadrature->getLambda(i);
      	for (int j = 0; j < nBasFcts; j++)
      	  phi[i][j] = (*(basisFunctions->getPhi(j)))(lambda);
      }
    
      // update flag
      init_flag |= INIT_PHI;
    }

    // initialize grd_phi
    if (grdPhi.empty() && init_flag.isSet(INIT_GRD_PHI)) {

      // allocate memory
      grdPhi.resize(nPoints);

      // fill memory
      for (int i = 0; i< nPoints; i++) {
      	grdPhi[i].resize(nBasFcts);
      	lambda = quadrature->getLambda(i);
      
      	for (int j = 0; j < nBasFcts; j++) {
      	  grdPhi[i][j].change_dim(dim + 1);
      	  (*(basisFunctions->getGrdPhi(j)))(lambda, grdPhi[i][j]);
      	}
      }
    
      // update flag
      init_flag |= INIT_GRD_PHI;
    }

    // initialize D2_phi

    if (!D2Phi && init_flag.isSet(INIT_D2_PHI)) {

      // allocate memory
      D2Phi = new MatrixOfFixVecs<DimMat<double> >(dim, nPoints, nBasFcts, NO_INIT);

      // fill memory
      for (int i = 0; i < nPoints; i++) {
      	lambda = quadrature->getLambda(i);
      	for (int j = 0; j < nBasFcts; j++)
      	  (*(basisFunctions->getD2Phi(j)))(lambda, (*(D2Phi))[i][j]);
      }
    
      // update flag
      init_flag |= INIT_D2_PHI;
    }
  }


  FastQuadrature::FastQuadrature(const FastQuadrature& fastQuad)
  {
    FUNCNAME("FastQuadrature::FastQuadrature()");

    TEST_EXIT(quadrature)("no quadrature!\n");

    int dim = quadrature->getDim();

    if (max_points == 0)
      max_points = Quadrature::maxNQuadPoints[dim];

    init_flag = fastQuad.init_flag;
    basisFunctions = fastQuad.basisFunctions;
    quadrature = fastQuad.quadrature;
  
    int nPoints = quadrature->getNumPoints();
    int nBasFcts = basisFunctions->getNumber();

    if (num_rows(fastQuad.phi) > 0) {
      phi.change_dim(nPoints, nBasFcts);
      phi = fastQuad.phi;
    }

    if (!fastQuad.grdPhi.empty()) {
      grdPhi.resize(nPoints);
      for (int i = 0; i < nPoints; i++) {
      	grdPhi[i].resize(nBasFcts);
      	for (int j = 0; j < nBasFcts; j++)
      	  grdPhi[i][j] = fastQuad.grdPhi[i][j];
      }
    }

    if (fastQuad.D2Phi) {
      D2Phi = new MatrixOfFixVecs<DimMat<double> >(dim, nPoints, nBasFcts, NO_INIT);
      for (int i = 0; i < nPoints; i++)
    	for (int j = 0; j < nBasFcts; j++)
    	  (*D2Phi)[i][j] = (*(fastQuad.D2Phi))[i][j];
    }
  }


  FastQuadrature::~FastQuadrature()
  {
    delete D2Phi;
  }


  double FastQuadrature::getSecDer(int q,int i ,int j, int m) const 
  {
    return (D2Phi) ? (*D2Phi)[q][i][j][m] : 0.0;
  }


  const VectorOfFixVecs<DimMat<double> > *FastQuadrature::getSecDer(int q) const 
  {
    return D2Phi ? (&((*D2Phi)[q])) : NULL;
  }
  
} // end namespace AMDiS
