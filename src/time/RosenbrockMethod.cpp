#include "time/RosenbrockMethod.hpp"

namespace AMDiS
{

  void RosenbrockMethod::createData()
  {
    a.resize(stages);
    for (int i = 0; i < stages; i++)
    {
      a[i].resize(stages);
      for (int j = 0; j < stages; j++)
        a[i][j] = 0.0;
    }

    c.resize(stages);
    for (int i = 0; i < stages; i++)
    {
      c[i].resize(stages);
      for (int j = 0; j < stages; j++)
        c[i][j] = 0.0;
    }

    m1.resize(stages);
    m2.resize(stages);
    for (int i = 0; i < stages; i++)
    {
      m1[i] = 0.0;
      m2[i] = 0.0;
    }
  }


  Ros2::Ros2()
  {
    order = 2;
    stages = 2;
    gamma = 1.707106781186547;

    createData();

    // b(2) = 0.5
    // b(1) = 1-b(2) = 0.5
    // alpha(2,1) = 1/(2*b(2)) = 1
    // gamma = 1+1/sqrt(2) = 1.707...
    // gamma(2,1) = -gamma/b(2) = -3.41421356237309
    // b_(1) = 1
    // b_(2) = 0

    a[0][0] = 0.0;
    a[1][0] = 5.857864376269050e-01;
    a[1][1] = 1.0;

    c[0][0] = gamma;
    c[1][0] = -1.171572875253810e+00;
    c[1][1] = -gamma;

    m1[0] = 8.786796564403575e-01;
    m1[1] = 2.928932188134525e-01;

    m2[0] = 5.857864376269050e-01;
    m2[1] = 0.0;

    MSG("Rosenbrock scheme Ros2\n");
  }


  Rowda3::Rowda3()
  {
    order = 3;
    stages = 3;
    gamma = 4.358665215084590e-01;

    createData();

    a[0][0] = 0.0;
    a[1][0] = 1.605996252195329e+00;
    a[1][1] = 0.7;
    a[2][0] = 1.605996252195329e+00;
    a[2][2] = 0.7;

    c[0][0] =  gamma;
    c[1][0] =  8.874044410657833e-01;
    c[1][1] =  0.604455284065559;
    c[2][0] =  2.398747971635036e+01;
    c[2][1] =  5.263722371562129e+00;
    c[2][2] =  6.37978879934488;

    m1[0] =  2.236727045296590e+00;
    m1[1] =  2.250067730969644e+00;
    m1[2] = -2.092514044390320e-01;

    m2[0] = 2.059356167645940e+00;
    m2[1] = 1.694014319346528e-01;

    MSG("Rosenbrock scheme Rowda3\n");
  }


  Ros3p::Ros3p()
  {
    order = 3;
    stages = 3;
    gamma = 7.886751345948129e-01;

    createData();

    a[0][0] = 0.0;
    a[1][0] = 1.267949192431123;
    a[1][1] = 1.0;
    a[2][0] = 1.267949192431123;
    a[2][1] = 0.0;
    a[2][2] = 1.0;

    c[0][0] =  gamma;
    c[1][0] = -1.607695154586736;
    c[1][1] = -0.2113248654051871;
    c[2][0] = -3.464101615137755;
    c[2][1] = -1.732050807568877;
    c[2][2] = -1.077350269189626;

    m1[0] = 2.0;
    m1[1] = 0.5773502691896258;
    m1[2] = 0.4226497308103742;

    m2[0] = 2.113248654051871;
    m2[1] = 1.0;
    m2[2] = 0.4226497308103742;

    MSG("Rosenbrock scheme Ros3p\n");
  }


  Rodasp::Rodasp()
  {
    order = 4;
    stages = 6;
    gamma = 0.25;

    createData();

    a[0][0] = 0.0;
    a[1][0] =  3.0;
    a[1][1] =  0.75;
    a[2][0] =  1.831036793486759e+00;
    a[2][1] =  4.955183967433795e-01;
    a[2][2] =  0.21;
    a[3][0] =  2.304376582692669e+00;
    a[3][1] = -5.249275245743001e-02;
    a[3][2] = -1.176798761832782e+00;
    a[3][3] =  0.63;
    a[4][0] = -7.170454962423025e+00;
    a[4][1] = -4.741636671481786e+00;
    a[4][2] = -1.631002631330971e+01;
    a[4][3] = -1.062004044111401e+00;
    a[4][4] =  1.0;
    a[5][0] = -7.170454962423025e+00;
    a[5][1] = -4.741636671481785e+00;
    a[5][2] = -1.631002631330971e+01;
    a[5][3] = -1.062004044111401e+00;
    a[5][4] =  1.0;
    a[5][5] =  1.0;

    c[0][0] =  gamma;
    c[1][0] = -1.200000000000000e+01;
    c[1][1] = -0.5;
    c[2][0] = -8.791795173947035e+00;
    c[2][1] = -2.207865586973518e+00;
    c[2][2] = -0.023504;
    c[3][0] =  1.081793056857153e+01;
    c[3][1] =  6.780270611428266e+00;
    c[3][2] =  1.953485944642410e+01;
    c[3][3] = -0.0362;
    c[4][0] =  3.419095006749677e+01;
    c[4][1] =  1.549671153725963e+01;
    c[4][2] =  5.474760875964130e+01;
    c[4][3] =  1.416005392148534e+01;
    c[4][4] =  0.0;
    c[5][0] =  3.462605830930533e+01;
    c[5][1] =  1.530084976114473e+01;
    c[5][2] =  5.699955578662667e+01;
    c[5][3] =  1.840807009793095e+01;
    c[5][4] = -5.714285714285717e+00;
    c[5][5] =  0.0;

    m1[0] = -7.170454962423026e+00;
    m1[1] = -4.741636671481786e+00;
    m1[2] = -1.631002631330971e+01;
    m1[3] = -1.062004044111401e+00;
    m1[4] =  1.0;
    m1[5] =  1.0;

    m2[0] = -7.170454962423026e+00;
    m2[1] = -4.741636671481786e+00;
    m2[2] = -1.631002631330971e+01;
    m2[3] = -1.062004044111401e+00;
    m2[4] =  1.0;

    MSG("Rosenbrock scheme Rodasp\n");
  }


  ROSI2P1::ROSI2P1()
  {
    order = 3;
    stages = 4;
    gamma = 4.35866521508459e-01;

    createData();

    a[0][0] =  0.0;
    a[1][0] =  1.147140180139521;
    a[1][1] =  0.5;
    a[2][0] =  1.785764587181959;
    a[2][1] =  0.442124760965983;
    a[2][2] =  0.75;
    a[3][0] =  2.506239510951673;
    a[3][1] =  4.558210876565182;
    a[3][2] = -1.373615544906449;
    a[3][3] =  1.0;

    c[0][0] =  gamma;
    c[1][0] = -2.631861185781065;
    c[1][1] = -0.064133478491541;
    c[2][0] = -3.011310475541004;
    c[2][1] =  0.334203214637756;
    c[2][2] = -0.145563307177156;
    c[3][0] = -3.743590594301783;
    c[3][1] = -1.089941238157158;
    c[3][2] = -1.718365430214442;
    c[3][3] = -0.135847884055848;

    m1[0] =  2.833375148827832;
    m1[1] =  3.953417886999603;
    m1[2] = -1.215227714218472;
    m1[3] =  1.165417447459307;

    m2[0] =  2.747785798103605;
    m2[1] =  1.770380793635233;
    m2[2] =  0.257316038155499;
    m2[3] =  0.343556220548095;

    MSG("Rosenbrock scheme ROSI2P1\n");
  }


  ROSI2P2::ROSI2P2()
  {
    order = 3;
    stages = 4;
    gamma = 4.3586652150845900e-01;

    createData();

    a[0][0] =  0.0;
    a[1][0] =  1.147140180139521;
    a[1][1] =  0.5;
    a[2][0] =  2.807348188211369;
    a[2][1] =  3.486932172067671;
    a[2][2] =  1.0;
    a[3][0] =  a[2][0];
    a[3][1] =  a[2][1];
    a[3][2] =  0.0;
    a[3][3] =  1.0;

    c[0][0] =  gamma;
    c[1][0] = -2.631861185781065;
    c[1][1] = -0.064133478491541;
    c[2][0] =  4.976389977276388;
    c[2][1] =  6.181041021340408;
    c[2][2] =  1.208496649176010;
    c[3][0] = -1.761050184345382;
    c[3][1] = -6.545972652439727;
    c[3][2] = -0.539706236424999;
    c[3][3] =  0.0;

    m1[0] =  2.807348188211369;
    m1[1] =  3.486932172067672;
    m1[2] =  0.0;
    m1[3] =  1.0;

    m2[0] =  0.420084258522926;
    m2[1] = -5.943299341711317;
    m2[2] =  0.360559439940373;
    m2[3] = -3.343738912424890;

    MSG("Rosenbrock scheme ROSI2P2\n");
  }

  ROSI2Pw::ROSI2Pw()
  {
    order = 3;
    stages = 4;
    gamma = 4.3586652150845e-01;

    createData();

    a[0][0] =  0.0;
    a[1][0] =  2.0;
    a[1][1] =  0.871733043016918;
    a[2][0] =  1.63034046718537;
    a[2][1] = -0.0903698030239437;
    a[2][2] =  0.75;
    a[3][0] =  3.28685750853379;
    a[3][1] =  15.4503745896552;
    a[3][2] = -15.0445558288081;
    a[3][3] =  1.0;

    c[0][0] = gamma;
    c[1][0] = -4.58856072055827;
    c[1][1] = -gamma;
    c[2][0] = -4.56739138878308;
    c[2][1] = -0.0683107605436897;
    c[2][2] = -0.418867127163069;
    c[3][0] = -2.99296365068231;
    c[3][1] = -42.6472578877056;
    c[3][2] =  43.6510246478391;
    c[3][3] =  0.0;

    m1[0] =  3.28685750853379;
    m1[1] =  15.4503745896552;
    m1[2] = -15.0445558288081;
    m1[3] =  1.0;

    m2[0] =  3.3904377475357;
    m2[1] =  5.865078412615;
    m2[2] = -4.96246395067988;
    m2[3] =  0.260825116305751;

    MSG("Rosenbrock scheme ROSI2Pw\n");
  }


  ROSI2PW::ROSI2PW()
  {
    order = 3;
    stages = 4;
    gamma = 4.3586652150845e-01;

    createData();

    a[0][0] =  0.0;
    a[1][0] =  2.0;
    a[1][1] =  0.871733043016918;
    a[2][0] = -5.501959790112310;
    a[2][1] = -1.833986596704078;
    a[2][2] = -1.59874671679705;
    a[3][0] =  4.338560720558247;
    a[3][1] =  1.147140180139553;
    a[3][2] = -0.059559354637499;
    a[3][3] =  1.0;

    c[0][0] = gamma;
    c[1][0] = -4.58856072055827;
    c[1][1] = -gamma;
    c[2][0] =  48.3965596116246;
    c[2][1] =  16.132186537208;
    c[2][2] =  6.56544000523295;
    c[3][0] = -2.31911621201727;
    c[3][1] = -1.14714018013959;
    c[3][2] = -0.0745075552140216;
    c[3][3] =  0.0;

    m1[0] =  4.33856072055832;
    m1[1] =  1.14714018013958;
    m1[2] = -0.0595593546375007;
    m1[3] =  1.0;

    m2[0] =  3.31383144448529;
    m2[1] =  1.14714018013955;
    m2[2] =  0.00847038665840751;
    m2[3] =  0.260825116305751;

    MSG("Rosenbrock scheme ROSI2PW\n");
  }


  Ros3Pw::Ros3Pw()
  {
    // index1 | index2 | pdes | R(infty) | stiffly acc.
    // -------+--------+------+----------+-------------
    //    1   |   0    |   1  |  0.73    |   0
    order = 3;
    stages = 3;
    gamma = 7.8867513459481287e-01;

    createData();

    a[0][0] = 0.0;
    a[1][0] = 2.0;
    a[1][1] = 1.57735026918963;
    a[2][0] = 0.633974596215561;
    a[2][1] = 0.0;
    a[2][2] = 0.5;

    c[0][0] =  gamma;
    c[1][0] = -2.53589838486225;
    c[1][1] = -gamma;
    c[2][0] = -1.62740473580836;
    c[2][1] = -0.274519052838329;
    c[2][2] = -0.0528312163512967;

    m1[0] =  1.63397459621556;
    m1[1] =  0.294228634059948;
    m1[2] =  1.07179676972449;

    m2[0] =  1.99444650053487  ;
    m2[1] =  0.654700538379252;
    m2[2] =  1.07179676972449;

    MSG("Rosenbrock scheme Ros3Pw\n");
  }



  Ros34PW2::Ros34PW2()
  {
    // index1 | index2 | pdes | R(infty) | stiffly acc.
    // -------+--------+------+----------+-------------
    //    1   |   0    |   1  |  0.0     |   1
    order = 3;
    stages = 4;
    gamma = 4.3586652150845e-01;

    createData();

    a[0][0] =  0.0;
    a[1][0] =  2.0;
    a[1][1] =  0.871733043016918;
    a[2][0] =  1.41921731745576;
    a[2][1] = -0.25923221167297;
    a[2][2] =  0.731579957788852;
    a[3][0] =  4.18476048231916;
    a[3][1] = -0.285192017355496;
    a[3][2] =  2.29428036027904;
    a[3][3] =  1.0;

    c[0][0] =  gamma;
    c[1][0] = -4.58856072055809;
    c[1][1] = -gamma;
    c[2][0] = -4.18476048231916;
    c[2][1] =  0.285192017355496;
    c[2][2] = -0.413333376233886;
    c[3][0] = -6.36817920012836;
    c[3][1] = -6.79562094446684;
    c[3][2] =  2.87009860433106;
    c[3][3] =  0.0;

    m1[0] =  4.18476048231916;
    m1[1] = -0.285192017355496;
    m1[2] =  2.29428036027904;
    m1[3] =  1.0;

    m2[0] =  3.90701053467119;
    m2[1] =  1.1180478778205;
    m2[2] =  0.521650232611491;
    m2[3] =  0.5;

    MSG("Rosenbrock scheme Ros34PW2\n");
  }
}