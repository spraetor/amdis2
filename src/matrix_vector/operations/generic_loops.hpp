#pragma once

#include "traits/at.hpp"

namespace AMDiS {
  namespace meta {

    // _________________________________________________________________________
    // generic loops
    
    /// generic loop for arrays, using the traits function at(...) for element access
    template <long I, long N> 
    struct LOOP 
    {
      /// for_each: op(vec) (elementwise, with shift)
      template <class A, class Op>
      static void for_each(A& a, Op op, size_t shift)
      {
	op(at(a, I+shift));  LOOP<I+1,N>::for_each(a, op, shift);
      }
      
      /// for_each: op(vec) (elementwise, without shift)
      template <class A, class Op>
      static void for_each(A& a, Op op)
      {
	op(at(a, I));  LOOP<I+1,N>::for_each(a, op);
      }
      
      /// transform: vec1 = op(vec2) (elementwise, with shift)
      template <class A, class B, class Op>
      static void transform(A const& a, B& b, Op op, size_t shift)
      {
	at(b, I+shift) = op(at(a, I+shift));  LOOP<I+1,N>::transform(a, b, op, shift);
      }
      
      /// transform: vec1 = op(vec2) (elementwise, without shift)
      template <class A, class B, class Op>
      static void transform(A const& a, B& b, Op op)
      {
	at(b, I) = op(at(a, I));  LOOP<I+1,N>::transform(a, b, op);
      }
      
      /// assign: vec1 = vec2, or vec1 += vec2, or vec1 *= vec2
      /// (compound assignment operator based on Assigner template)
      template <class A, class B, class Assigner>
      static void assign(A& a, B const& b, Assigner assigner, size_t shift)
      {
	Assigner::apply(at(a, I+shift), at(b, I+shift));
	LOOP<I+1,N>::assign(a, b, assigner, shift);
      }
      
      template <class A, class B, class Assigner>
      static void assign(A& a, B const& b, Assigner assigner)
      {
	Assigner::apply(at(a, I), at(b, I));
	LOOP<I+1,N>::assign(a, b, assigner);
      }
      
      /// accumulate: sum_i{ op(vec_i) }, prod_i{ op(vec_i) }
      /// (outer operation is a binary operator, inner operation is a unary operator)
      template <class A, class T, class Op, class BinaryOp>
      static T accumulate(A const& a, T init, Op op, BinaryOp binary_op, size_t shift)
      {
	return binary_op(LOOP<I+1,N>::accumulate(a, init, op, binary_op, shift), 
			 op(at(a, I+shift)));
      }
      
      template <class A, class T, class Op, class BinaryOp>
      static T accumulate(A const& a, T init, Op op, BinaryOp binary_op)
      {
	return binary_op(LOOP<I+1,N>::accumulate(a, init, op, binary_op), op(at(a, I)));
      }
      
      template <class A, class T, class Functor>
      static void accumulate(A const& a, T& init, Functor f)
      {
	Functor::update(init, at(a, I));
	LOOP<I+1,N>::accumulate(a, init, f);
      }
      
      /// inner_product: sum_i{ b_op(vec1_i, vec2_i) }, prod_i{ b_op(vec1_i, vec2_i) }
      /// (outer operation is a binary operator, inner operation is a binary operator)
      template <class A, class B, class T, class BinaryOp1, class BinaryOp2>
      static T inner_product(A const& a, B const& b, T init, BinaryOp1 binary_op1, BinaryOp2 binary_op2, size_t shift)
      {
	return binary_op1(LOOP<I+1,N>::inner_product(a, b, init, binary_op1, binary_op2, shift), 
			  binary_op2(at(a, I+shift), at(b, I+shift)));
      }
      
      template <class A, class B, class T, class BinaryOp1, class BinaryOp2>
      static T inner_product(A const& a, B const& b, T init, BinaryOp1 binary_op1, BinaryOp2 binary_op2)
      {
	return binary_op1(LOOP<I+1,N>::inner_product(a, b, init, binary_op1, binary_op2), 
			  binary_op2(at(a, I), at(b, I)));
      }
      
      // inner product, using one functor that includes both binary operators
      template <class A, class B, class T, class Functor>
      static void inner_product(A const& a, B const& b, T& init, Functor f)
      {
	Functor::update(init, at(a, I), at(b, I));
	LOOP<I+1,N>::inner_product(a, b, init, f);
      }
    };
    
    /// \cond HIDDEN_SYMBOLS
    // Abbruchbedingung I==N
    template <long N> 
    struct LOOP<N, N>
    { 
      template <class A, class Op>
      static void for_each(A&, Op, size_t = 0) {}
      
      template <class A, class B, class Op>
      static void transform(A const&, B&, Op, size_t = 0) {}
      
      template <class A, class B, class Assigner>
      static void assign(A&, B const&, Assigner, size_t = 0) {}
      
      template <class A, class T, class Op, class BinaryOp>
      static T accumulate(A const&, T init, Op, BinaryOp, size_t = 0)
      {
	return init;
      }
      
      template <class A, class T, class Functor>
      static void accumulate(A const&, T&, Functor) { }
      
      template <class A, class B, class T, class BinaryOp1, class BinaryOp2>
      static T inner_product(A const&, B const&, T init, BinaryOp1, BinaryOp2, size_t = 0)
      {
	return init;
      }
      
      template <class A, class B, class T, class Functor>
      static void inner_product(A const&, B const&, T&, Functor) { }
    };
    /// \endcond
    
    
    
    
    
    
    /// generic loop for expressions, using the operator() for element access
    template <long I, long N> 
    struct FOR 
    { 
      /// for_each: op(vec) (elementwise)
      template<class A, class Op>
      static void for_each(A& a, Op op, size_t shift)
      {
	op(a(I+shift));  FOR<I+1,N>::for_each(a, op, shift);
      }
      
      // specialization for pointer types
      template<class T, class Op>
      static void for_each(T* a, Op op, size_t shift)
      {
	op(a[I+shift]);  FOR<I+1,N>::for_each(a, op, shift);
      }
      
      template<class A, class Op>
      static void for_each(A& a, Op op)
      {
	op(a(I));  FOR<I+1,N>::for_each(a, op);
      }
      
      // specialization for pointer types
      template<class T, class Op>
      static void for_each(T* a, Op op)
      {
	op(a[I]);  FOR<I+1,N>::for_each(a, op);
      }
      
      /// transform: vec1 = op(vec2) (elementwise)
      template<class A, class B, class Op>
      static void transform(A const& a, B& b, Op op, size_t shift)
      {
	b(I+shift) = op(a(I+shift));
	FOR<I+1,N>::transform(a, b, op, shift);
      }
      
      template<class A, class B, class Op>
      static void transform(A const& a, B& b, Op op)
      {
	b(I) = op(a(I)); 
	FOR<I+1,N>::transform(a, b, op);
      }
      
      /// assign: vec1 = vec2, or vec1 += vec2, or vec1 *= vec2
      /// (compound assignment operator based on Assigner template)
      template<class A, class B, class Assigner>
      static void assign(A& a, B const& b, Assigner assigner, size_t shift)
      {
	Assigner::apply(a(I+shift), b(I+shift));
	FOR<I+1,N>::assign(a, b, assigner, shift);
      }
      
      template<class A, class B, class Assigner>
      static void assign(A& a, B const& b, Assigner assigner)
      {
	Assigner::apply(a(I), b(I));
	FOR<I+1,N>::assign(a, b, assigner);
      }
      
      /// accumulate: sum_i{ op(vec_i) }, prod_i{ op(vec_i) }
      /// (outer operation is a binary operator, inner operation is a unary operator)
      template<class A, class T, class Op, class BinaryOp>
      static T accumulate(A const& a, T init, Op op, BinaryOp binary_op, size_t shift)
      {
	return binary_op(FOR<I+1,N>::accumulate(a, init, op, binary_op, shift), 
			 op(a(I+shift)));
      }
      
      template<class A, class T, class Op, class BinaryOp>
      static T accumulate(A const* a, T init, Op op, BinaryOp binary_op, size_t shift)
      {
	return binary_op(FOR<I+1,N>::accumulate(a, init, op, binary_op, shift), 
			 op(a[I+shift]));
      }
      
      template<class A, class T, class Op, class BinaryOp>
      static T accumulate(A const& a, T init, Op op, BinaryOp binary_op)
      {
	return binary_op(FOR<I+1,N>::accumulate(a, init, op, binary_op), op(a(I)));
      }
      
      // specialization for pointer types
      template<class A, class T, class Op, class BinaryOp>
      static T accumulate(A const* a, T init, Op op, BinaryOp binary_op)
      {
	return binary_op(FOR<I+1,N>::accumulate(a, init, op, binary_op), op(a[I]));
      }
      
      template<class A, class T, class Functor>
      static void accumulate(A const& a, T& init, Functor f)
      {
	Functor::update(init, a(I));
	FOR<I+1,N>::accumulate(a, init, f);
      }
      
      // specialization for pointer types
      template<class A, class T, class Functor>
      static void accumulate(A const* a, T& init, Functor f)
      {
	Functor::update(init, a[I]);
	FOR<I+1,N>::accumulate(a, init, f);
      }
      
      /// inner_product: sum_i{ b_op(vec1_i, vec2_i) }, prod_i{ b_op(vec1_i, vec2_i) }
      /// (outer operation is a binary operator, inner operation is a binary operator)
      template<class A, class B, class T, class BinaryOp1, class BinaryOp2>
      static T inner_product(A const& a, B const& b, T init, BinaryOp1 binary_op1, BinaryOp2 binary_op2, size_t shift)
      {
	return binary_op1(FOR<I+1,N>::inner_product(a, b, init, binary_op1, binary_op2, shift), 
			  binary_op2(a(I+shift), b(I+shift)));
      }
      
      template<class A, class B, class T, class BinaryOp1, class BinaryOp2>
      static T inner_product(A const& a, B const& b, T init, BinaryOp1 binary_op1, BinaryOp2 binary_op2)
      {
	return binary_op1(FOR<I+1,N>::inner_product(a, b, init, binary_op1, binary_op2), 
			  binary_op2(a(I), b(I)));
      }
      
      // inner product, using one functor that includes both binary operators
      template<class A, class B, class T, class Functor>
      static void inner_product(A const& a, B const& b, T& init, Functor f)
      {
	Functor::update(init, a(I), b(I));
	FOR<I+1,N>::inner_product(a, b, init, f);
      }
    };

    /// \cond HIDDEN_SYMBOLS
    // Abbruchbedingung I==N
    template<long N> 
    struct FOR<N, N> : LOOP<N,N> {};
    
//     { 
//       template<class A, class Op>
//       static void for_each(A&, Op, size_t = 0) {}
//       
//       template<class A, class B, class Op>
//       static void transform(A const&, B&, Op, size_t = 0) {}
//       
//       template<class A, class B, class Assigner>
//       static void assign(A&, B const&, Assigner, size_t = 0) {}
//       
//       template<class A, class T, class Op, class BinaryOp>
//       static T accumulate(A const&, T init, Op, BinaryOp, size_t = 0)
//       {
// 	return init;
//       }
//       
//       template<class A, class T, class Functor>
//       static void accumulate(A const&, T&, Functor) { }
//       
//       template<class A, class B, class T, class BinaryOp1, class BinaryOp2>
//       static T inner_product(A const&, B const&, T init, BinaryOp1, BinaryOp2, size_t = 0)
//       {
// 	return init;
//       }
//       
//       template<class A, class B, class T, class Functor>
//       static void inner_product(A const&, B const&, T&, Functor) { }
//     };
    /// \endcond
    
  } // end namespace meta
} // end namespace AMDiS
