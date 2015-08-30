/** \file value_expr.hpp */

#pragma once

#include "LazyOperatorTerm.h"

namespace AMDiS
{
  namespace expressions
  {
    /// Expression that encapsulates a runtime value
    template <class T>
    struct RValue : public LazyOperatorTermBase
    {
      typedef T value_type;
      T value;
      RValue(const T& value_) : value(value_) {}

      inline value_type operator()(const int& iq) const
      {
        return value;
      }

      std::string str() const
      {
        return std::to_string(value);
      }
    };


    /// Expression that encapsulates a compiletime value
    template <int V>
    struct CValue : public LazyOperatorTermBase
    {
      typedef int value_type;
      inline value_type operator()(const int& iq) const
      {
        return V;
      }

      std::string str() const
      {
        return std::string("[") + std::to_string(V) + "]";
      }
    };


    /// Expression that encapsulates a compiletime vector
    template <int Size, int V0 = 0, int V1 = 0, int V2 = 0>
    struct CVector : public LazyOperatorTermBase
    {
      typedef Vector<int> value_type;
      value_type V;

      CVector() : V(Size)
      {
        int values[3] = {V0, V1, V2};
        V.setValues(values);
      }

      inline value_type operator()(const int& iq) const
      {
        return V;
      }

      std::string str() const
      {
        return std::string("[V(") + std::to_string(Size) + ")]";
      }
    };


    /// Expression that encapsulates a compiletime vector
    template <int V0 = 0, int V1 = 0, int V2 = 0>
    struct WVector : public LazyOperatorTermBase
    {
      typedef WorldVector<int> value_type;
      value_type V;

      WVector()
      {
        V[0] = V0;
        if (size(V) > 1) V[1] = V1;
        if (size(V) > 2) V[2] = V2;
      }

      inline value_type operator()(const int& iq) const
      {
        return V;
      }

      std::string str() const
      {
        return std::string("[V(") + std::to_string(size(V)) + ")]";
      }
    };

    template <int I = -1>
    struct E : public LazyOperatorTermBase
    {
      typedef WorldVector<int> value_type;
      value_type V;

      E(int I_ = -1)
      {
        assert((I >= 0 && I_ < 0) || (I < 0 && I_ >= 0));
        V = 0;
        V[std::max(I, I_)] = 1;
      }

      inline value_type operator()(const int& iq) const
      {
        return V;
      }

      std::string str() const
      {
        return std::string("[V(") + std::to_string(size(V)) + ")]";
      }
    };

    template<> struct E<0> : WVector<1> {};
    template<> struct E<1> : WVector<0,1> {};
    template<> struct E<2> : WVector<0,0,1> {};


    /// Expression that encapsulates a compiletime matrix
    template <int Rows, int Cols, int V0 = 0, int V1 = 0, int V2 = 0,
              int V3 = 0, int V4 = 0, int V5 = 0,
              int V6 = 0, int V7 = 0, int V8 = 0>
    struct CMatrix : public LazyOperatorTermBase
    {
      typedef Matrix<int> value_type;
      value_type M;

      CMatrix() : M(Rows, Cols)
      {
        int values[9] = {V0, V1, V2, V3, V4, V5, V6, V7, V8};
        M.setValues(values);
      }

      inline value_type operator()(const int& iq) const
      {
        return M;
      }

      std::string str() const
      {
        return std::string("[M(") + std::to_string(Rows) + ","+ std::to_string(Cols) + ")]";
      }
    };


    /// Expression that encapsulates a compiletime matrix
    template<int V0 = 0, int V1 = 0, int V2 = 0,
             int V3 = 0, int V4 = 0, int V5 = 0,
             int V6 = 0, int V7 = 0, int V8 = 0>
    struct WMatrix : public LazyOperatorTermBase
    {
      typedef WorldMatrix<int> value_type;
      value_type M;

      WMatrix()
      {
        M[0][0] = V0;
        if (Global::getGeo(WORLD) == 2)
        {
          M[0][1] = V1;
          M[1][0] = V2;
          M[1][1] = V3;
        }
        else if (Global::getGeo(WORLD) == 3)
        {
          M[0][1] = V1;
          M[0][2] = V2;
          M[1][0] = V3;
          M[1][1] = V4;
          M[1][2] = V5;
          M[2][0] = V6;
          M[2][1] = V7;
          M[2][2] = V8;
        }
      }

      inline value_type operator()(const int& iq) const
      {
        return M;
      }

      std::string str() const
      {
        return std::string("[M(") + std::to_string(Global::getGeo(WORLD)) + ","+ std::to_string(Global::getGeo(WORLD)) + ")]";
      }
    };


    template <int N = -1>
    struct Eye : public LazyOperatorTermBase
    {
      typedef Matrix<int> value_type;
      size_t n;
      value_type M;

      Eye(int n_) : n(n_ < 0 ? N : n_), M(n,n)
      {
        M = 0;
        for (size_t i = 0; i < n; ++i)
          M[i][i] = 1;
      }

      inline value_type operator()(const int& iq) const
      {
        return M;
      }

      std::string str() const
      {
        return std::string("[M(") + std::to_string(n) + ","+ std::to_string(n) + ")]";
      }
    };

    template <>
    struct Eye<-2> : public LazyOperatorTermBase
    {
      typedef WorldMatrix<int> value_type;
      value_type M;

      Eye()
      {
        M = 0;
        for (size_t i = 0; i < (size_t)Global::getGeo(WORLD); ++i)
          M[i][i] = 1;
      }

      inline value_type operator()(const int& iq) const
      {
        return M;
      }

      std::string str() const
      {
        return std::string("[M(") + std::to_string(Global::getGeo(WORLD)) + ","+ std::to_string(Global::getGeo(WORLD)) + ")]";
      }
    };

    template<> struct Eye<1> : CMatrix<1,1, 1> {};
    template<> struct Eye<2> : CMatrix<2,2, 1,0, 0,1> {};
    template<> struct Eye<3> : CMatrix<3,3, 1,0,0, 0,1,0, 0,0,1> {};

  } // end namespace expressions


  namespace traits
  {
    template <class T>
    struct remove_all_qualifiers
    {
      typedef typename boost::remove_cv
      <typename boost::remove_pointer
      <typename boost::remove_reference<T>::type>::type
      >::type type;
    };


    template <class T>
    struct pure_value
    {
      typedef typename remove_all_qualifiers<T>::type type;
      static type eval(const T& t)
      {
        return t;
      }
    };

    template <class T>
    struct pure_value<T&>
    {
      typedef typename remove_all_qualifiers<T>::type type;
      static type eval(T& t)
      {
        return t;
      }
    };

    template <class T>
    struct pure_value<T*>
    {
      typedef typename remove_all_qualifiers<T>::type type;
      static type eval(T* t)
      {
        return pure_value<T>::eval(*t);
      }
    };

    template <class T>
    struct pure_value<const T>
    {
      typedef typename remove_all_qualifiers<T>::type type;
      static type eval(const T& t)
      {
        return pure_value<T>::eval(t);
      }
    };

  } // end namespace traits

  namespace expressions
  {
    /// Expression that points to a value given by reference
    template <class T>
    struct Reference : public LazyOperatorTermBase
    {
      typedef typename traits::pure_value<T>::type value_type;
      const value_type& value;
      Reference(const T& value_) : value(value_) {}
      Reference(const T* value_) : value(*value_) {}

      template <class List>
      void insertFeSpaces(List& feSpaces) const {}

      int getDegree() const
      {
        return 0;
      }

      void initElement(const ElInfo* elInfo,
                       SubAssembler* subAssembler, Quadrature* quad,
                       const BasisFunction* basisFct = NULL) {}

      inline value_type operator()(const int& iq) const
      {
        return traits::pure_value<T>::eval(value);
      }

      std::string str() const
      {
        return std::string("&(") + std::to_string(traits::pure_value<T>::eval(value)) + ")";
      }
    };

  } // end namespace expressions


  // generator functions
  // _____________________________________________________________________________

  template <class T>
  inline expressions::RValue<T>
  constant_(T value)
  {
    return {value};
  }

  template <int I>
  inline expressions::CValue<I>
  constant_()
  {
    return {};
  }

  template <class T>
  inline expressions::Reference<T>
  ref_(T& value)
  {
    return {value};
  }

  template <class T>
  inline expressions::Reference<T>
  ref_(T* value)
  {
    return {value};
  }


  // unit vectors
  template <int I>
  inline expressions::E<I>
  E_()
  {
    return {};
  }

  inline expressions::E<-1>
  E_(int i)
  {
    return{i};
  }

  // unit matrix
  template <int N>
  inline expressions::Eye<N>
  eye_()
  {
    return {};    // NxN unit-matrix
  }

  inline expressions::Eye<-1>
  eye_(int n)
  {
    return {n};    // nxn unit-matrix
  }

  inline expressions::Eye<-2>
  eye_()
  {
    return {};    // dow x dow unit-matrix
  }


} // end namespace AMDiS
