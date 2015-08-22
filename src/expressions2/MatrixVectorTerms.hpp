/** \file MatrixVectorTerms.hpp */

#pragma once

#include <matrix_vector/MatrixVectorOperations.hpp>
#include "BinaryTerms.hpp"
#include "TermConcepts.hpp"

namespace AMDiS
{  
  namespace functors
  {
    struct TwoNorm
    {
      constexpr int getDegree(int d) const { return 2*d+1; }
      template <class T>
      auto operator()(T&& t) const RETURNS( two_norm(std::forward<T>(t)) )
    };
    
    struct OneNorm
    {
      constexpr int getDegree(int d) const { return d; }
      template <class T>
      auto operator()(T&& t) const RETURNS( one_norm(std::forward<T>(t)) )
    };
    
    struct UnaryDot
    {
      constexpr int getDegree(int d) const { return 2*d; }
      template <class T>
      auto operator()(T&& t) const RETURNS( unary_dot(std::forward<T>(t)) )
    };
    
    struct Max
    {
      constexpr int getDegree(int d) const { return d; }
      template <class T>
      auto operator()(T&& t) const RETURNS( AMDiS::max(std::forward<T>(t)) )
    };
    
    struct AbsMax
    {
      constexpr int getDegree(int d) const { return d; }
      template <class T>
      auto operator()(T&& t) const RETURNS( AMDiS::abs_max(std::forward<T>(t)) )
    };
    
    struct Min
    {
      constexpr int getDegree(int d) const { return d; }
      template <class T>
      auto operator()(T&& t) const RETURNS( AMDiS::min(std::forward<T>(t)) )
    };
    
    struct AbsMin
    {
      constexpr int getDegree(int d) const { return d; }
      template <class T>
      auto operator()(T&& t) const RETURNS( AMDiS::abs_min(std::forward<T>(t)) )
    };
    
    struct Sum
    {
      constexpr int getDegree(int d) const { return d; }
      template <class T>
      auto operator()(T&& t) const RETURNS( sum(std::forward<T>(t)) )
    };
    
    struct Mean
    {
      constexpr int getDegree(int d) const { return d; }
      template <class T>
      auto operator()(T&& t) const RETURNS( mean(std::forward<T>(t)) )
    };
    
    struct Prod
    {
      constexpr int getDegree(int d) const { return d*d*d; } // maximal degree for dow=3
      template <class T>
      auto operator()(T&& t) const RETURNS( prod(std::forward<T>(t)) )
    };
    
    struct Diagonal
    {
      constexpr int getDegree(int d0) const { return d0; }
      
    private:
      template <class M, class V = detail::MatrixToVector<M> >
      static typename V::type eval(M const& m, tag::matrix)
      {
        using value_type = Value_t<M>;
        typename V::type vector(num_rows(m), value_type(0)); 
        for (size_t i = 0; i < num_rows(m); ++i)
          vector(i) = m(i,i);
        return vector;
      }
      
      template <class V, class M = detail::VectorToMatrix<V> >
      static typename M::type eval(V const& v, tag::vector)
      {
        using value_type = Value_t<V>;
        typename M::type matrix(size(v), size(v), value_type(0)); 
        for (size_t i = 0; i < size(v); ++i)
          matrix(i,i) = v(i);
        return matrix;
      }
      
    public:      
      template <class T>
      auto operator()(T&& t) const RETURNS
      ( 
        Diagonal::eval(std::forward<T>(t), typename traits::category<Decay_t<T>>::tag()) 
      )
    };
    
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    
    struct Dot
    {
      constexpr int getDegree(int d0, int d1) const { return d0+d1; }
      template <class T1, class T2>
      auto operator()(T1&& t1, T2&& t2) const RETURNS
      ( 
        dot(std::forward<T1>(t1), std::forward<T2>(t2)) 
      )
    };
    
    struct Cross
    {
      constexpr int getDegree(int d0, int d1) const { return d0+d1; }
      
      template <class Vec1, class Vec2, 
        class = Requires_t<and_<traits::is_vector<Vec1>, traits::is_vector<Vec2>>> >
      Assign_t<Vec1> operator()(Vec1 const& v1, Vec2 const& v2) const
      {
        using size_type = Size_t<traits::category<Vec1>>;
        Assign_t<Vec1> result(size(v1));
        
        TEST_EXIT_DBG( size(v1) == 3 && size(v1) == size(v2) )
          ("cross: inkompatible sizes!\n");

        for (size_type i = 0; i < size(v1); ++i) {
          size_type k = (i+1) % 3, l = (i+2) % 3;
          result(i) = v1(k) * v2(l) - v1(l) * v2(k);
        }
        return result;
      }
    };
    
    struct Outer
    {
      constexpr int getDegree(int d0, int d1) const { return d0+d1; }
      
      template <class Vec1, class Vec2, 
        class = Requires_t<and_<traits::is_vector<Vec1>, traits::is_vector<Vec2>>> >
      VectorToMatrix_t<Vec1> operator()(Vec1 const& v1, Vec2 const& v2) const
      {
        using size_type1 = Size_t<traits::category<Vec1>>;
        using size_type2 = Size_t<traits::category<Vec2>>;
        VectorToMatrix_t<Vec1> result(size(v1), size(v2));

        for (size_type1 i = 0; i < size(v1); ++i) {
          for (size_type2 j = 0; j < size(v2); ++j) {
            result(i,j) = v1(i) * v2(j);
          }
        }
        return result;
      }
    };
    
    struct Distance
    {
      constexpr int getDegree(int d0, int d1) const { return d0+d1+1; }
      template <class T1, class T2>
      auto operator()(T1&& t1, T2&& t2) const RETURNS
      ( 
        distance(std::forward<T1>(t1), std::forward<T2>(t2)) 
      )
    };
  }
  
  
  // ---------------------------------------------------------------------------
  
      
  template <class T>
    UnaryTerm<T, functors::OneNorm>
  inline one_norm(VectorTerm<T> const& t) { return {t.sub()}; }
  
  template <class T>
    UnaryTerm<T, functors::TwoNorm>
  inline two_norm(VectorTerm<T> const& t) { return {t.sub()}; }
  
  template <class T>
    UnaryTerm<T, functors::TwoNorm>
  inline frobenius_norm(MatrixTerm<T> const& t) { return {t.sub()}; }
  
  template <class T>
    UnaryTerm<T, functors::TwoNorm>
  inline norm(VectorTerm<T> const& t) { return {t.sub()}; }
  
  template <class T>
    UnaryTerm<T, functors::TwoNorm>
  inline norm(MatrixTerm<T> const& t) { return {t.sub()}; }
  
  template <class T>
    UnaryTerm<T, functors::UnaryDot>
  inline unary_dot(VectorTerm<T> const& t) { return {t.sub()}; }
  
  template <class T>
    UnaryTerm<T, functors::Max>
  inline max(BaseTerm<T> const& t) { return {t.sub()}; }
  
  template <class T>
    UnaryTerm<T, functors::AbsMax>
  inline abs_max(BaseTerm<T> const& t) { return {t.sub()}; }
  
  template <class T>
    UnaryTerm<T, functors::Min>
  inline min(BaseTerm<T> const& t) { return {t.sub()}; }
  
  template <class T>
    UnaryTerm<T, functors::AbsMin>
  inline min(BaseTerm<T> const& t) { return {t.sub()}; }
  
  template <class T>
    UnaryTerm<T, functors::Sum>
  inline sum(BaseTerm<T> const& t) { return {t.sub()}; }
  
  template <class T>
    UnaryTerm<T, functors::Mean>
  inline mean(BaseTerm<T> const& t) { return {t.sub()}; }
  
  template <class T>
    UnaryTerm<T, functors::Prod>
  inline prod(BaseTerm<T> const& t) { return {t.sub()}; }
  
  template <class T>
    UnaryTerm<T, functors::Diagonal>
  inline diag(VectorTerm<T> const& t) { return {t.sub()}; }
  
  template <class T>
    UnaryTerm<T, functors::Diagonal>
  inline diag(MatrixTerm<T> const& t) { return {t.sub()}; }
  
  
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
  
  template <class T1, class T2>
    requires::Term< BinaryTerm<ToTerm_t<T1>, ToTerm_t<T2>, functors::Dot>, T1, T2 >  
  inline dot(T1&& t1, T2&& t2)
  {
    return {toTerm(std::forward<T1>(t1)), toTerm(std::forward<T2>(t2))};
  }
  
  template <class T1, class T2>
    requires::Term< BinaryTerm<ToTerm_t<T1>, ToTerm_t<T2>, functors::Cross>, T1, T2 >  
  inline cross(T1&& t1, T2&& t2)
  {
    return {toTerm(std::forward<T1>(t1)), toTerm(std::forward<T2>(t2))};
  }
  
  template <class T1, class T2>
    requires::Term< BinaryTerm<ToTerm_t<T1>, ToTerm_t<T2>, functors::Outer>, T1, T2 >  
  inline outer(T1&& t1, T2&& t2)
  {
    return {toTerm(std::forward<T1>(t1)), toTerm(std::forward<T2>(t2))};
  }
  
  template <class T1, class T2>
    requires::Term< BinaryTerm<ToTerm_t<T1>, ToTerm_t<T2>, functors::Distance>, T1, T2 >  
  inline distance(T1&& t1, T2&& t2)
  {
    return {toTerm(std::forward<T1>(t1)), toTerm(std::forward<T2>(t2))};
  }

} // end namespace AMDiS
