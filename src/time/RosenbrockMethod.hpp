/** \file RosenbrochMethod.h */

#pragma once

#include "AMDiS_fwd.hpp"
#include "CreatorInterface.hpp"

namespace AMDiS
{

  class RosenbrockMethod
  {
  public:
    RosenbrockMethod() {}

    virtual ~RosenbrockMethod() {}

    inline int getOrder()
    {
      return order;
    }

    inline int getStages()
    {
      return stages;
    }

    inline double getGamma()
    {
      return gamma;
    }

    inline double getA(int i, int j)
    {
      return a[i][j];
    }

    inline double getC(int i, int j)
    {
      return c[i][j];
    }

    inline double getM1(int i)
    {
      return m1[i];
    }

    inline double getM2(int i)
    {
      return m2[i];
    }

    inline double getAlphaI(int i)
    {
      return a[i][i];
    }

    inline double getGammaI(int i)
    {
      return c[i][i];
    }

  protected:
    void createData();

  protected:
    int order;

    int stages;

    double gamma;

    std::vector<std::vector<double>> a;

    std::vector<std::vector<double>> c;

    std::vector<double> m1;

    std::vector<double> m2;
  };

  class RosenbrockMethodCreator : public CreatorInterface<RosenbrockMethod>
  {
  public:
    RosenbrockMethodCreator()
    {}
  };


  class Ros2 : public RosenbrockMethod
  {
  public:
    class Creator : public RosenbrockMethodCreator
    {
    public:
      Creator() : RosenbrockMethodCreator() {}

      ~Creator() {}

      RosenbrockMethod* create()
      {
        return new Ros2();
      }
    };

    Ros2();

    ~Ros2() {}
  };


  class Rowda3 : public RosenbrockMethod
  {
  public:
    class Creator : public RosenbrockMethodCreator
    {
    public:
      Creator() : RosenbrockMethodCreator() {}

      ~Creator() {}

      RosenbrockMethod* create()
      {
        return new Rowda3();
      }
    };

    Rowda3();

    ~Rowda3() {}
  };


  class Ros3p : public RosenbrockMethod
  {
  public:
    class Creator : public RosenbrockMethodCreator
    {
    public:
      Creator() : RosenbrockMethodCreator() {}

      ~Creator() {}

      RosenbrockMethod* create()
      {
        return new Ros3p();
      }
    };

    Ros3p();

    ~Ros3p() {}
  };


  class Rodasp : public RosenbrockMethod
  {
  public:
    class Creator : public RosenbrockMethodCreator
    {
    public:
      Creator() : RosenbrockMethodCreator() {}

      ~Creator() {}

      RosenbrockMethod* create()
      {
        return new Rodasp();
      }
    };

    Rodasp();

    ~Rodasp() {}
  };


  class ROSI2P1 : public RosenbrockMethod
  {
  public:
    class Creator : public RosenbrockMethodCreator
    {
    public:
      Creator() : RosenbrockMethodCreator() {}

      ~Creator() {}

      RosenbrockMethod* create()
      {
        return new ROSI2P1();
      }
    };

    ROSI2P1();

    ~ROSI2P1() {}
  };


  class ROSI2P2 : public RosenbrockMethod
  {
  public:
    class Creator : public RosenbrockMethodCreator
    {
    public:
      Creator() : RosenbrockMethodCreator() {}

      ~Creator() {}

      RosenbrockMethod* create()
      {
        return new ROSI2P2();
      }
    };

    ROSI2P2();

    ~ROSI2P2() {}
  };


  class ROSI2Pw : public RosenbrockMethod
  {
  public:
    class Creator : public RosenbrockMethodCreator
    {
    public:
      Creator() : RosenbrockMethodCreator() {}

      ~Creator() {}

      RosenbrockMethod* create()
      {
        return new ROSI2Pw();
      }
    };

    ROSI2Pw();

    ~ROSI2Pw() {}
  };


  class ROSI2PW : public RosenbrockMethod
  {
  public:
    class Creator : public RosenbrockMethodCreator
    {
    public:
      Creator() : RosenbrockMethodCreator() {}

      ~Creator() {}

      RosenbrockMethod* create()
      {
        return new ROSI2PW();
      }
    };

    ROSI2PW();

    ~ROSI2PW() {}
  };


  class Ros3Pw : public RosenbrockMethod
  {
  public:
    class Creator : public RosenbrockMethodCreator
    {
    public:
      Creator() : RosenbrockMethodCreator() {}

      ~Creator() {}

      RosenbrockMethod* create()
      {
        return new Ros3Pw();
      }
    };

    Ros3Pw();

    ~Ros3Pw() {}
  };


  class Ros34PW2 : public RosenbrockMethod
  {
  public:
    class Creator : public RosenbrockMethodCreator
    {
    public:
      Creator() : RosenbrockMethodCreator() {}

      ~Creator() {}

      RosenbrockMethod* create()
      {
        return new Ros34PW2();
      }
    };

    Ros34PW2();

    ~Ros34PW2() {}
  };
}
