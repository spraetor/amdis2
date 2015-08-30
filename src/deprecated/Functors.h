/******************************************************************************
 *
 * AMDiS - Adaptive multidimensional simulations
 *
 * Copyright (C) 2013 Dresden University of Technology. All Rights Reserved.
 * Web: https://fusionforge.zih.tu-dresden.de/projects/amdis
 *
 * Authors:
 * Simon Vey, Thomas Witkowski, Andreas Naumann, Simon Praetorius, et al.
 *
 * This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
 * WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 *
 *
 * This file is part of AMDiS
 *
 * See also license.opensource.txt in the distribution.
 *
 ******************************************************************************/

#ifndef AMDIS_FUNCTORS_H
#define AMDIS_FUNCTORS_H

#include "AbstractFunction.h"
#include <boost/math/special_functions/cbrt.hpp>
#include <boost/math/special_functions/pow.hpp>

// NOTE: should be removed or replaced by functors in operation/functors.hpp

namespace AMDiS
{

  template<typename T>
  struct Id : public AbstractFunction<T,T>
  {
    Id(int degree = 1) : AbstractFunction<T,T>(degree) {}
    T operator()(const T& v) const
    {
      return v;
    }
  };

  template<typename T1,typename T2>
  struct Const : public AbstractFunction<T1,T2>
  {
    Const(T1 val_) : AbstractFunction<T1,T2>(0), val(val_) {}
    T1 operator()(const T2& v) const
    {
      return val;
    }
  private:
    T1 val;
  };

  template<typename T=double>
  struct Factor : public AbstractFunction<T,T>
  {
    Factor(double fac_, int degree = 1) : AbstractFunction<T,T>(degree), fac(fac_) {}
    T operator()(const T& x) const
    {
      return fac*x;
    }
  private:
    double fac;
  };

  template<typename T=double>
  struct Add : public BinaryAbstractFunction<T,T,T>
  {
    Add(int degree = 1) : BinaryAbstractFunction<T,T,T>(degree) {}
    T operator()(const T& v1, const T& v2) const
    {
      return v1+v2;
    }
  };

  template<typename T=double>
  struct AddFactor : public BinaryAbstractFunction<T,T,T>
  {
    AddFactor(double factor_ = 1.0, int degree = 1) : BinaryAbstractFunction<T,T,T>(degree), factor(factor_) {}
    T operator()(const T& v1, const T& v2) const
    {
      return v1 + factor*v2;
    }
  private:
    double factor;
  };

  template<typename T=double>
  struct Subtract : public BinaryAbstractFunction<T,T,T>
  {
    Subtract(int degree = 1) : BinaryAbstractFunction<T,T,T>(degree) {}
    T operator()(const T& v1, const T& v2) const
    {
      return v1-v2;
    }
  };

  template<typename T=double>
  struct AddScal : public AbstractFunction<T,T>
  {
    AddScal(T scal_, int degree = 1) : AbstractFunction<T,T>(degree), scal(scal_) {}
    T operator()(const T& v) const
    {
      return v+scal;
    }
  private:
    T scal;
  };

  template<typename T>
  struct Mult : public BinaryAbstractFunction<T,T,T>
  {
    Mult(int degree = 2) : BinaryAbstractFunction<T,T,T>(degree) {}
    T operator()(const T& v1, const T& v2) const
    {
      return v1*v2;
    }
  };

  template<typename T1, typename T2, typename T3>
  struct Mult2 : public BinaryAbstractFunction<T1,T2,T3>
  {
    Mult2(int degree = 2) : BinaryAbstractFunction<T1,T2,T3>(degree) {}
    T1 operator()(const T2& v1, const T3& v2) const
    {
      return v1*v2;
    }
  };

  template<typename T>
  struct MultScal : public BinaryAbstractFunction<T,T,T>
  {
    MultScal(T scal_, int degree = 2) : BinaryAbstractFunction<T,T,T>(degree), scal(scal_) {}
    T operator()(const T& v1, const T& v2) const
    {
      return v1*v2*scal;
    }
  private:
    T scal;
  };

  template<typename T>
  struct Max : public BinaryAbstractFunction<T,T,T>
  {
    Max(int degree = 1) : BinaryAbstractFunction<T,T,T>(degree) {}
    T operator()(const T& v1, const T& v2) const
    {
      return std::max(v1,v2);
    }
  };

  template<typename T>
  struct Min : public BinaryAbstractFunction<T,T,T>
  {
    Min(int degree = 1) : BinaryAbstractFunction<T,T,T>(degree) {}
    T operator()(const T& v1, const T& v2) const
    {
      return std::min(v1,v2);
    }
  };

  template<typename T>
  struct Diff : public BinaryAbstractFunction<T,T,T>
  {
    Diff(int degree = 1) : BinaryAbstractFunction<T,T,T>(degree) {}
    T operator()(const T& v1, const T& v2) const
    {
      return std::abs(v1-v2);
    }
  };

  template<typename T=double>
  struct Abs : public AbstractFunction<T,T>
  {
    Abs(int degree = 1) : AbstractFunction<T,T>(degree) {}
    T operator()(const T& v) const
    {
      return std::abs(v);
    }
  };

  template<typename T=double>
  struct Signum : public AbstractFunction<T,T>
  {
    Signum() : AbstractFunction<T,T>(0) {}
    T operator()(const T& v) const
    {
      return (v>0.0?1.0:(v<0.0?-1.0:0.0));
    }
  };

  template<typename T=double>
  struct Sqr : public AbstractFunction<T,T>
  {
    Sqr(int degree = 2) : AbstractFunction<T, T>(degree) {}
    T operator()(const T& v) const
    {
      return sqr(v);
    }
  };

  template<typename T=double>
  struct Sqrt : public AbstractFunction<T,T>
  {
    Sqrt(int degree = 4) : AbstractFunction<T,T>(degree) {}
    T operator()(const T& v) const
    {
      return std::sqrt(v);
    }
  };

  namespace detail
  {
    template<int p, typename T, typename Enabled = void>
    struct Pow
    {
      typedef typename traits::mult_type<T, typename Pow<p-1,T>::result_type>::type result_type;
      static result_type eval(const T& v)
      {
        return v*Pow<p-1,T>::eval(v);
      }
    };

    template<int p, typename T>
    struct Pow<p, T, typename boost::enable_if_c<
      boost::is_same<T, typename traits::mult_type<T, T>::type>::value&&
    (p > 1)
    >::type >
    {
      typedef T result_type;
      static T eval(const T& v)
    {
      return boost::math::pow<p>(v);
    }
    };

    template<typename T>
    struct Pow<1,T>
    {
      typedef T result_type;
      static result_type eval(const T& v)
      {
        return v;
      }
    };

    template<typename T>
    struct Pow<0, T>
    {
      typedef double result_type;
      static result_type eval(const T& v)
      {
        return 1.0;
      }
    };
  }

  template<int p, typename T=double>
  struct Pow : public AbstractFunction<typename detail::Pow<p,T>::result_type, T>
  {
    typedef typename detail::Pow<p,T>::result_type result_type;
    Pow(double factor_=1.0, int degree = p) : AbstractFunction<result_type,T>(degree), factor(factor_) {}
    result_type operator()(const T& v) const
    {
      return factor * detail::Pow<p,T>::eval(v);
    }
  private:
    double factor;
  };

  template<typename TIn, typename TOut = typename traits::mult_type<TIn, TIn>::type>
  struct Norm2 : public AbstractFunction<TOut, TIn>
  {
    Norm2(int degree = 4) : AbstractFunction<TOut, TIn>(degree) {}
    TOut operator()(const TIn& v) const
    {
      return std::sqrt(v*v);
    }
  };

  template<typename TIn, typename TOut = typename traits::mult_type<TIn, TIn>::type>
  struct Norm2Sqr : public AbstractFunction<TOut, TIn>
  {
    Norm2Sqr(int degree = 2) : AbstractFunction<TOut, TIn>(degree) {}
    TOut operator()(const TIn& v) const
    {
      return v*v;
    }
  };

  template<typename T>
  struct Norm2_comp2 : public BinaryAbstractFunction<T,T,T>
  {
    Norm2_comp2(int degree = 4) : BinaryAbstractFunction<T,T,T>(degree) {}
    T operator()(const T& v1, const T& v2) const
    {
      return std::sqrt(sqr(v1)+sqr(v2));
    }
  };

  template<typename T>
  struct Norm2Sqr_comp2 : public BinaryAbstractFunction<T,T,T>
  {
    Norm2Sqr_comp2(int degree = 2) : BinaryAbstractFunction<T,T,T>(degree) {}
    T operator()(const T& v1, const T& v2) const
    {
      return sqr(v1)+sqr(v2);
    }
  };

  template<typename T>
  struct Norm2_comp3 : public TertiaryAbstractFunction<T,T,T,T>
  {
    Norm2_comp3(int degree = 4) : TertiaryAbstractFunction<T,T,T,T>(degree) {}
    T operator()(const T& v1, const T& v2, const T& v3) const
    {
      return std::sqrt(sqr(v1)+sqr(v2)+sqr(v3));
    }
  };

  template<typename T>
  struct Norm2Sqr_comp3 : public TertiaryAbstractFunction<T,T,T,T>
  {
    Norm2Sqr_comp3(int degree = 2) : TertiaryAbstractFunction<T,T,T,T>(degree) {}
    T operator()(const T& v1, const T& v2, const T& v3) const
    {
      return sqr(v1)+sqr(v2)+sqr(v3);
    }
  };

  template<typename T>
  struct L1Diff : public BinaryAbstractFunction<T,T,T>
  {
    T operator()(const T& v1, const T& v2) const
    {
      return std::abs(v1-v2);
    }
  };

  template<typename TOut, typename T=TOut>
  struct L2Diff : public BinaryAbstractFunction<TOut,T,T>
  {
    TOut operator()(const T& v1, const T& v2) const
    {
      return Norm2<TOut, T>()(v1-v2);
    }
  };

  template<typename T>
  struct Vec1WorldVec : public AbstractFunction<WorldVector<T>,T>
  {
    WorldVector<T> operator()(const T& v0) const
    {
      WorldVector<T> result;
      result[0]=v0;
      return result;
    }
  };
  template<typename T>
  struct Vec2WorldVec : public BinaryAbstractFunction<WorldVector<T>,T,T>
  {
    WorldVector<T> operator()(const T& v0, const T& v1) const
    {
      WorldVector<T> result;
      result[0]=v0;
      result[1]=v1;
      return result;
    }
  };
  template<typename T>
  struct Vec3WorldVec : public TertiaryAbstractFunction<WorldVector<T>,T,T,T>
  {
    WorldVector<T> operator()(const T& v0, const T& v1, const T& v2) const
    {
      WorldVector<T> result;
      result[0]=v0;
      result[1]=v1;
      result[2]=v2;
      return result;
    }
  };
  template<int c, typename T=double>
  struct Component : public AbstractFunction<T, WorldVector<T>>
  {
    Component(int degree = 1) : AbstractFunction<T, WorldVector<T>>(degree) {}
    T operator()(const WorldVector<T>& x) const
    {
      return x[c];
    }
  };
  template<typename T=double>
  struct Component2 : public AbstractFunction<T, WorldVector<T>>
  {
    Component2(int comp, int degree = 1) : AbstractFunction<T, WorldVector<T>>(degree), c_(comp) {}
    T operator()(const WorldVector<T>& x) const
    {
      return x[c_];
    }
  private:
    int c_;
  };

  struct FadeOut : public TertiaryAbstractFunction<double, double, double ,double>
  {
    double operator()(const double& v, const double& dist, const double& mean) const
    {
      return dist*mean+(1.0-dist)*v;
    }
  };

  struct Random : public AbstractFunction<double, WorldVector<double>>
  {
    Random(double mean_, double amplitude_) : mean(mean_), amplitude(amplitude_)
    {
      std::srand(time(0));
    }

    double operator()(const WorldVector<double>& x) const
    {
      return mean + 2.0*amplitude * ((std::rand() / static_cast<double>(RAND_MAX)) - 0.5);
    }

  private:
    double mean;
    double amplitude;
  };

}

#endif // AMDIS_FUNCTORS_H

