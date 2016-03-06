#pragma once

// AMDiS includes
#include "expressions/FunctorTerm.hpp"
#include "expressions/TermConcepts.hpp"
#include "matrix_vector/MatrixVectorOperations.hpp"


/// Macro that generates a unary (vector) functor.
/**
 *  \p NAME    Name of the class.
 *  \p DEGREE  Expression in 'd0' that gives the polynomial degree, with
 *             'd0' the degree of the argument passed to the functor.
 *  \p FCT     Name of a unary c++-function that represents the functor.
 */
#define AMDIS_MAKE_VECTOR_FUNCTOR( NAME, DEGREE, FCT )    \
    struct NAME                                           \
    {                                                     \
      constexpr int getDegree(int d0) const               \
      {                                                   \
        return DEGREE ;                                   \
      }                                                   \
      template <class T>                                  \
      auto operator()(T const& t) const RETURNS                \
      (                                                   \
        FCT( t )                         \
      )                                                   \
    };


namespace AMDiS
{
  namespace functors
  {

    AMDIS_MAKE_VECTOR_FUNCTOR( TwoNorm,  2*d0+1, two_norm       )
    AMDIS_MAKE_VECTOR_FUNCTOR( OneNorm,      d0, one_norm       )
    AMDIS_MAKE_VECTOR_FUNCTOR( UnaryDot,   2*d0, unary_dot      )
    AMDIS_MAKE_VECTOR_FUNCTOR( Sum,          d0, sum            )
    AMDIS_MAKE_VECTOR_FUNCTOR( Mean,         d0, mean           )
    AMDIS_MAKE_VECTOR_FUNCTOR( Prod,   d0*d0*d0, prod           )
    AMDIS_MAKE_VECTOR_FUNCTOR( Max,          d0, AMDiS::max     )
    AMDIS_MAKE_VECTOR_FUNCTOR( AbsMax,       d0, AMDiS::abs_max )
    AMDIS_MAKE_VECTOR_FUNCTOR( Min,          d0, AMDiS::min     )
    AMDIS_MAKE_VECTOR_FUNCTOR( AbsMin,       d0, AMDiS::abs_min )


    /// Convert a vector to a diagonal matrix and extract the diagonal of
    /// a matrix and store it in a vector.
    struct Diagonal
    {
      constexpr int getDegree(int d0) const
      {
        return d0;
      }

    private:
      // Extract the diagonal of a matrix
      template <class M, class V = detail::MatrixToVector<M>>
      static typename V::type eval(M const& m, tag::matrix)
      {
        using value_type = Value_t<M>;
        typename V::type vector(num_rows(m), value_type{0});
        for (size_t i = 0; i < num_rows(m); ++i)
          vector(i) = m(i,i);
        return vector;
      }

      // Generate a diagonal matrix
      template <class V, class M = detail::VectorToMatrix<V>>
      static typename M::type eval(V const& v, tag::vector)
      {
        using value_type = Value_t<V>;
        typename M::type matrix(size(v), size(v), value_type{0});
        for (size_t i = 0; i < size(v); ++i)
          matrix(i,i) = v(i);
        return matrix;
      }

    public:
      template <class T>
      auto operator()(T&& t) const RETURNS
      (
        Diagonal::eval(std::forward<T>(t),
		       typename traits::category<Decay_t<T>>::tag())
      )
    };

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    struct Dot
    {
      constexpr int getDegree(int d0, int d1) const
      {
        return d0+d1;
      }

      template <class T1, class T2>
      auto operator()(T1&& t1, T2&& t2) const RETURNS
      (
        dot(std::forward<T1>(t1), std::forward<T2>(t2))
      )
    };

    struct Cross
    {
      constexpr int getDegree(int d0, int d1) const
      {
        return d0+d1;
      }

      template <class Vec1, class Vec2,
	class = Requires_t<and_<traits::is_vector<Vec1>, traits::is_vector<Vec2>>>>
      Assign_t<Vec1> operator()(Vec1 const& v1, Vec2 const& v2) const
      {
        using size_type = Size_t<traits::category<Vec1>>;
        Assign_t<Vec1> result(size(v1));

        TEST_EXIT_DBG( size(v1) == 3 && size(v1) == size(v2) )
        ("cross: inkompatible sizes!\n");

        for (size_type i = 0; i < size(v1); ++i)
        {
          size_type k = (i+1) % 3, l = (i+2) % 3;
          result(i) = v1(k) * v2(l) - v1(l) * v2(k);
        }
        return result;
      }
    };

    struct Outer
    {
      constexpr int getDegree(int d0, int d1) const
    {
      return d0+d1;
    }

    template <class Vec1, class Vec2,
      class = Requires_t<and_<traits::is_vector<Vec1>, traits::is_vector<Vec2>>>>
    VectorToMatrix_t<Vec1> operator()(Vec1 const& v1, Vec2 const& v2) const
    {
      using size_type1 = Size_t<traits::category<Vec1>>;
      using size_type2 = Size_t<traits::category<Vec2>>;
      VectorToMatrix_t<Vec1> result(size(v1), size(v2));

      for (size_type1 i = 0; i < size(v1); ++i)
      {
        for (size_type2 j = 0; j < size(v2); ++j)
        {
          result(i,j) = v1(i) * v2(j);
        }
      }
      return result;
    }
        };

    struct Distance
    {
      constexpr int getDegree(int d0, int d1) const
      {
	return d0+d1+1;
      }

      template <class T1, class T2>
      auto operator()(T1&& t1, T2&& t2) const RETURNS
      (
	distance(std::forward<T1>(t1), std::forward<T2>(t2))
      )
    };
  }


  // ---------------------------------------------------------------------------


  template <class T>
  FunctorTerm<functors::OneNorm, T>
  inline one_norm(VectorTerm<T> const& t)
  {
    return {t.sub()};
  }

  template <class T>
  FunctorTerm<functors::TwoNorm, T>
  inline two_norm(VectorTerm<T> const& t)
  {
    return {t.sub()};
  }

  template <class T>
  FunctorTerm<functors::TwoNorm, T>
  inline frobenius_norm(MatrixTerm<T> const& t)
  {
    return {t.sub()};
  }

  template <class T>
  FunctorTerm<functors::TwoNorm, T>
  inline norm(VectorTerm<T> const& t)
  {
    return {t.sub()};
  }

  template <class T>
  FunctorTerm<functors::TwoNorm, T>
  inline norm(MatrixTerm<T> const& t)
  {
    return {t.sub()};
  }

  template <class T>
  FunctorTerm<functors::UnaryDot, T>
  inline unary_dot(VectorTerm<T> const& t)
  {
    return {t.sub()};
  }

  template <class T>
  FunctorTerm<functors::Diagonal, T>
  inline diag(VectorTerm<T> const& t)
  {
    return {t.sub()};
  }

  template <class T>
  FunctorTerm<functors::Diagonal, T>
  inline diag(MatrixTerm<T> const& t)
  {
    return {t.sub()};
  }


  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


  template <class T1, class T2>
  requires::Term<FunctorTerm<functors::Dot, ToTerm_t<T1>, ToTerm_t<T2>>, T1, T2>
  inline dot(T1&& t1, T2&& t2)
  {
    return {toTerm(std::forward<T1>(t1)), toTerm(std::forward<T2>(t2))};
  }

  template <class T1, class T2>
  requires::Term<FunctorTerm<functors::Cross, ToTerm_t<T1>, ToTerm_t<T2>>, T1, T2>
  inline cross(T1&& t1, T2&& t2)
  {
    return {toTerm(std::forward<T1>(t1)), toTerm(std::forward<T2>(t2))};
  }

  template <class T1, class T2>
  requires::Term<FunctorTerm<functors::Outer, ToTerm_t<T1>, ToTerm_t<T2>>, T1, T2>
  inline outer(T1&& t1, T2&& t2)
  {
    return {toTerm(std::forward<T1>(t1)), toTerm(std::forward<T2>(t2))};
  }

  template <class T1, class T2>
  requires::Term<FunctorTerm<functors::Distance, ToTerm_t<T1>, ToTerm_t<T2>>, T1, T2>
  inline distance(T1&& t1, T2&& t2)
  {
    return {toTerm(std::forward<T1>(t1)), toTerm(std::forward<T2>(t2))};
  }

} // end namespace AMDiS
