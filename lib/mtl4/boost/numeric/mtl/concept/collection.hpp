// Software License for MTL
// 
// Copyright (c) 2007 The Trustees of Indiana University.
//               2008 Dresden University of Technology and the Trustees of Indiana University.
//               2010 SimuNova UG (haftungsbeschränkt), www.simunova.com.
// All rights reserved.
// Authors: Peter Gottschling and Andrew Lumsdaine
// 
// This file is part of the Matrix Template Library
// 
// See also license.mtl.txt in the distribution.

#ifndef MTL_COLLECTION_INCLUDE
#define MTL_COLLECTION_INCLUDE

#include <boost/type_traits/remove_const.hpp>
#include <boost/numeric/mtl/mtl_fwd.hpp>
#include <vector>

#ifdef __GXX_CONCEPTS__
#  include <concepts>
#else 
#  include <boost/numeric/linear_algebra/pseudo_concept.hpp>
#  include <boost/numeric/mtl/concept/std_concept.hpp>
#endif

#include <boost/numeric/mtl/utility/transposed_orientation.hpp>

namespace mtl {

/** @addtogroup Concepts
 *  @{
 */

#ifdef __GXX_CONCEPTS__
    auto concept Collection<typename T>
    {
	typename value_type;
	typename const_reference;
	typename size_type;
    };
#else
    /// Concept Collection
    template <typename T>
    struct Collection
    {
	/// Associated type: elements in the collection
	typedef associated_type value_type;

	/// Associated type: return type of const element access (however implemented)
	typedef associated_type const_reference;

	/// Associated type: size type used for indexing in collection
	typedef associated_type size_type;
    };
#endif


#ifdef __GXX_CONCEPTS__
    auto concept MutableCollection<typename T>
      : Collection<T>
    {
	typename reference;
    }
#else
    /// Concept MutableCollection
    template <typename T>
    struct MutableCollection
	: public Collection<T>
    {
	/// Associated type: return type of non-const element access (however implemented)
	typedef associated_type reference;
    };
#endif


#ifdef __GXX_CONCEPTS__
    concept ConstantSizeCollection<typename T>
      : Collection<T>
    {};
#else
    /// Concept ConstantSizeCollection: size parameters of collection are completely given at compile time
    /* Which parameters determine collection size depends on type of collection, e.g. different for vector and matrix
       \par Refinement of:
       - Collection < T >
    */
    template <typename T>
    struct ConstantSizeCollection
	: Collection<T>
    {};
#endif


#ifdef __GXX_CONCEPTS__
    auto concept AlgebraicCollection<typename T>
      : Collection<T>
    {
	size_type num_rows(T);
	size_type num_cols(T);
	size_type size(T);
    };
#else
    /// Concept AlgebraicCollection: common requirements of matrices, vectors, and scalars in computations
    /** For more design clarity we consider them all as matrices (as Matlab does) and we regard 
	Scalar and Vector as special cases (see there).  However, the implementation of vectors
	is supposed to differ from the ones of matrices in order to provide efficient computations and storage.
        \par Refinement of:
	- Collection < T >
	\par Notation:
	- X is a type that models AlgebraicCollection
	- x is an object of type X
	\par Valid expressions:
	- Number of rows: \n num_rows(x) \n Return Type: size_type
	- Number of columns: \n num_cols(x) \n Return Type: size_type
	- Number of elements: \n size(x) \n Return Type: size_type
	  \n Sematics: num_rows(x) * num_cols(x) (but possibly faster implemented)
    */
    template <typename T>
    struct AlgebraicCollection
	: public Collection<T>
    {};
#endif


#ifdef __GXX_CONCEPTS__
    auto concept ConstantSizeAlgebraicCollection<typename T>
      : AlgebraicCollection<T>,
        ConstantSizeCollection<T>
    {
#if 0
	// Is there a way to require static members????
	static Collection<T>::size_type T::static_num_rows;
	static Collection<T>::size_type T::static_num_cols;
	static Collection<T>::size_type T::static_size;
#endif 
    };
#else
    /// Concept ConstantSizeAlgebraicCollection: extension of AlgebraicCollection with meta-functions
    /** This concept is used for algebraic collections with sizes known at compile time. 
	The motivation is that if the size of the collection is
	is small, arithmetic operations can be unrolled at compile time.

        \par Refinement of:
	- Collection < T >
	\par Notation:
	- X is a type that models ConstantSizeAlgebraicCollection
	- x is an object of type X
	\par Valid expressions:
	- Number of rows: \n static_num_rows<X>::value
	- Number of columns: \n static_num_cols<X>::value
	- Number of elements: \n static_size<X>::value
	  \n Sematics: static_num_rows<X>::value * static_size<X>::value
	\note
	-# For more design clarity we consider them all as matrices (as Matlab does) and we regard 
	   Scalar and Vector as special cases (see there).  However, the implementation of vectors
	   is supposed to differ from the ones of matrices in order to provide efficient computations and storage.

    */
    template <typename T>
    struct ConstantSizeAlgebraicCollection
      : public AlgebraicCollection<T>,
        public ConstantSizeCollection<T>
    {
	/// Associated type: meta-function for number of rows
	typedef associated_type static_num_rows;
	/// Associated type: meta-function for number of columns
	typedef associated_type static_num_cols;
	/// Associated type: meta-function for number of elements
	typedef associated_type static_size;
    };
#endif



#ifdef __GXX_CONCEPTS__
    auto concept TraversableCollection<typename Tag, typename C>
      : Collection<C>
    {
#if 0
	// This might be impossible to declare with concepts
	typename cursor_type;

	cursor_type begin<Tag>(const C& c);
	cursor_type   end<Tag>(const C& c);
#endif

	// Maybe we switch to this syntax for the sake of concepts
	typename cursor_type;

	cursor_type begin(const C& c, Tag);
	cursor_type   end(const C& c, Tag);
    }
#else
    /// Concept TraversableCollection: collections that can be traversed by cursor or iterator
    template <typename Tag, typename C>
    struct TraversableCollection
	: public Collection<C>
    {
	/// Associated type: return type of tagged begin and end function
	typedef associated_type cursor_type;

	/// Tagged free function that returns a cursor or iterator at the begin of an interval 
	/** The interval is specified by the Tag, i.e. the function is called begin<Tag>(c); */
	cursor_type begin(const C& c);

	/// Tagged free function that returns a cursor or iterator at the end of an interval 
	/** The interval is specified by the Tag, i.e. the function is called end<Tag>(c);  */
	cursor_type end(const C& c);
    };
#endif


#ifdef __GXX_CONCEPTS__
    auto concept TraversableMutableCollection<typename Tag, typename C>
      : MutableCollection<C>
    {
#if 0
	typename cursor_type;

	cursor_type begin<Tag>(C& c);
	cursor_type   end<Tag>(C& c);
#endif

	// Maybe we switch to this syntax for the sake of concepts
	typename cursor_type;

	cursor_type begin(C& c, Tag);
	cursor_type   end(C& c, Tag);
    }
#else
    /// Concept TraversableMutableCollection: collections that can be traversed by (mutable) iterator
    template <typename Tag, typename C>
    struct TraversableMutableCollection
	: public MutableCollection<C>
    {
	/// Associated type: return type of tagged begin function
	typedef associated_type cursor_type;

	/// Tagged free function that returns a cursor or iterator at the begin of an interval 
	/** The interval is specified by the Tag, i.e. the function is called begin<Tag>(c); */
	cursor_type begin(const C& c);

	/// Tagged free function that returns a cursor or iterator at the end of an interval 
	/** The interval is specified by the Tag, i.e. the function is called end<Tag>(c);  */
	cursor_type end(const C& c);
    };
#endif




#ifdef __GXX_CONCEPTS__
#if 0
    concept CategorizedType<typename T>
    {
	typedef associated_type type;
    };
#endif
#endif


#ifdef __GXX_CONCEPTS__
    concept OrientedCollection<typename T>
      : Collection<T>
    {
	typename orientation;

    };
#else
    /// Concept OrientedCollection: collections with concept-awareness in terms of associated type
    /** Concept-awareness is given for matrices as well as for vectors consistent to the unification in
	AlgebraicCollection. The orientation of vectors determines whether it is a row or a column vector.
	The orientation of matrices only characterizes the internal representation and has no semantic consequences.
        \par Refinement of:
	- Collection < T >
	\par Associated type:
	- orientation
    */
    template <typename T>
    struct OrientedCollection
      : public Collection<T>
    {
	/// Associated type for orientation; by default identical with member type
	typedef typename T::orientation orientation;
    };
#endif







// ============================================
// Concept maps (and emulations as type traits)
// ============================================

#ifdef __GXX_CONCEPTS__
    // Needs redefinition in refinements (at least in conceptg++)
    // -> as a consequence we define it directly there
#else
    template <typename Value, typename Parameters>
    struct Collection<mtl::matrix::dense2D<Value, Parameters> >
    {
	typedef Value            value_type;
	typedef const Value&     const_reference;
        // Alleged ambiguity with mtl::tag::dense2D on MSVC
        typedef typename mtl::matrix::dense2D<Value, Parameters>::size_type size_type;
    };
#endif


#ifdef __GXX_CONCEPTS__

#else
    template <typename Value, std::size_t Mask, typename Parameters>
    struct Collection<mtl::matrix::morton_dense<Value, Mask, Parameters> >
    {
	typedef Value            value_type;
	typedef const Value&     const_reference;
	typedef typename Parameters::size_type size_type;
    };
#endif


#ifdef __GXX_CONCEPTS__
    template <typename Value, typename Parameters>
    concept_map Collection<mtl::matrix::compressed2D<Value, Parameters> >
    {
	typedef Value            value_type;
	typedef Value            const_reference;
	typedef typename mtl::matrix::compressed2D<Value, Parameters>::size_type size_type;
    };
#else
    template <typename Value, typename Parameters>
    struct Collection<mtl::matrix::compressed2D<Value, Parameters> >
    {
	typedef Value            	       value_type;
	typedef Value            	       const_reference;
	typedef typename Parameters::size_type size_type;
    };

#endif


#ifdef __GXX_CONCEPTS__
#else
    template <typename Value, typename Parameters>
    struct Collection<mtl::matrix::coordinate2D<Value, Parameters> >
    {
	typedef Value            	       value_type;
	typedef Value            	       const_reference;
	typedef typename Parameters::size_type size_type;
    };

#endif


#ifdef __GXX_CONCEPTS__
#else
    template <typename Value, typename Parameters>
    struct Collection<mtl::matrix::sparse_banded<Value, Parameters> >
    {
	typedef Value            	       value_type;
	typedef Value            	       const_reference;
	typedef typename Parameters::size_type size_type;
    };

#endif


#ifdef __GXX_CONCEPTS__
    template <typename Vector>
    concept_map Collection<mtl::matrix::multi_vector<Vector> >
    {
	typedef typename mtl::matrix::multi_vector<Vector>::value_type       value_type;
	typedef typename mtl::matrix::multi_vector<Vector>::const_reference  const_reference;
	typedef typename mtl::matrix::multi_vector<Vector>::size_type        size_type;
    };
#else
    template <typename Vector>
    struct Collection<mtl::matrix::multi_vector<Vector> >
    {
	typedef typename mtl::matrix::multi_vector<Vector>::value_type       value_type;
	typedef typename mtl::matrix::multi_vector<Vector>::const_reference  const_reference;
	typedef typename mtl::matrix::multi_vector<Vector>::size_type        size_type;
    };

#endif




#ifdef __GXX_CONCEPTS__
    template <typename Vector>
    concept_map Collection<mtl::matrix::multi_vector_range<Vector> >
    {
	typedef typename mtl::matrix::multi_vector_range<Vector>::value_type       value_type;
	typedef typename mtl::matrix::multi_vector_range<Vector>::const_reference  const_reference;
	typedef typename mtl::matrix::multi_vector_range<Vector>::size_type        size_type;
    };
#else
    template <typename Vector>
    struct Collection<mtl::matrix::multi_vector_range<Vector> >
    {
	typedef typename mtl::matrix::multi_vector_range<Vector>::value_type       value_type;
	typedef typename mtl::matrix::multi_vector_range<Vector>::const_reference  const_reference;
	typedef typename mtl::matrix::multi_vector_range<Vector>::size_type        size_type;
    };

#endif

#ifdef __GXX_CONCEPTS__
#else
    template <typename Value>
    struct Collection<matrix::element_structure<Value> >
    {
	typedef Value            value_type;
      // typedef Value            const_reference;
      // typedef typename mtl::matrix::compressed2D<Value, Parameters>::size_type size_type;
    };

#endif

#ifdef __GXX_CONCEPTS__
#else
    template <typename Value, typename Parameters>
    struct Collection<matrix::ell_matrix<Value, Parameters> >
    {
	typedef Value            	       value_type;
	typedef Value            	       const_reference;
	typedef typename Parameters::size_type size_type;
    };

#endif


#ifdef __GXX_CONCEPTS__
    template <typename Scaling, typename Coll>
    concept_map Collection<matrix::scaled_view<Scaling, Coll> >
    {
	typedef typename matrix::scaled_view<Scaling, Coll>::value_type        value_type;
	typedef typename matrix::scaled_view<Scaling, Coll>::const_reference   const_reference;
	typedef typename matrix::scaled_view<Scaling, Coll>::size_type         size_type;
    };
#else
    template <typename Scaling, typename Coll>
    struct Collection<matrix::scaled_view<Scaling, Coll> >
    {
	typedef typename matrix::scaled_view<Scaling, Coll>::value_type        value_type;
	typedef typename matrix::scaled_view<Scaling, Coll>::const_reference   const_reference;
	typedef typename matrix::scaled_view<Scaling, Coll>::size_type         size_type;
    };
#endif

// added by Hui Li
#ifdef __GXX_CONCEPTS__
    template <typename Coll, typename RScaling>
    concept_map Collection<matrix::rscaled_view<Coll,RScaling> >
    {
	typedef typename matrix::rscaled_view<Coll,RScaling>::value_type        value_type;
	typedef typename matrix::rscaled_view<Coll,RScaling>::const_reference   const_reference;
	typedef typename matrix::rscaled_view<Coll,RScaling>::size_type         size_type;
    };
#else
    template <typename Coll, typename RScaling>
    struct Collection<matrix::rscaled_view<Coll,RScaling> >
    {
	typedef typename matrix::rscaled_view<Coll,RScaling>::value_type        value_type;
	typedef typename matrix::rscaled_view<Coll,RScaling>::const_reference   const_reference;
	typedef typename matrix::rscaled_view<Coll,RScaling>::size_type         size_type;
    };
#endif
	


#ifdef __GXX_CONCEPTS__
    // template <typename Functor, typename Coll>
    // concept_map Collection< vector::map_view<Functor, Coll> >
    // {
    // 	typedef typename vector::map_view<Functor, Coll>::value_type        value_type;
    // 	typedef typename vector::map_view<Functor, Coll>::const_reference   const_reference;
    // 	typedef typename vector::map_view<Functor, Coll>::size_type         size_type;
    // };
#else
    template <typename Functor, typename Coll>
    struct Collection< vector::map_view<Functor, Coll> >
    {
	// typedef typename Functor::result_type              value_type;
	// typedef typename Functor::result_type              const_reference;
	// typedef typename Collection<Coll>::size_type     size_type;
	
	typedef typename vector::map_view<Functor, Coll>::value_type        value_type;
	typedef typename vector::map_view<Functor, Coll>::const_reference   const_reference;
	typedef typename vector::map_view<Functor, Coll>::size_type         size_type;
    };
#endif

	
#ifdef __GXX_CONCEPTS__

#else
    template <typename Scaling, typename Coll>
    struct Collection<vector::scaled_view<Scaling, Coll> >
      : Collection< vector::map_view<tfunctor::rscale<typename Collection<Coll>::value_type, Scaling>, Coll> >
    {
	// typedef typename vector::scaled_view<Scaling, Coll>::value_type        value_type;
	// typedef typename vector::scaled_view<Scaling, Coll>::const_reference   const_reference;
	// typedef typename vector::scaled_view<Scaling, Coll>::size_type         size_type;
    };
#endif



// added by Hui Li
#ifdef __GXX_CONCEPTS__
    template <typename Coll, typename RScaling>
    concept_map Collection<vector::rscaled_view<Coll,RScaling> >
    {
	typedef typename vector::rscaled_view<Coll,RScaling>::value_type        value_type;
	typedef typename vector::rscaled_view<Coll,RScaling>::const_reference   const_reference;
	typedef typename vector::rscaled_view<Coll,RScaling>::size_type         size_type;
    };
#else
    template <typename Coll, typename RScaling>
    struct Collection<vector::rscaled_view<Coll,RScaling> >
      : Collection< vector::map_view<tfunctor::rscale<typename Collection<Coll>::value_type, RScaling>, Coll> >
    {
	// typedef typename vector::rscaled_view<Coll,RScaling>::value_type        value_type;
	// typedef typename vector::rscaled_view<Coll,RScaling>::const_reference   const_reference;
	// typedef typename vector::rscaled_view<Coll,RScaling>::size_type         size_type;
    };
#endif


#ifdef __GXX_CONCEPTS__
    template <typename Coll>
    concept_map Collection<vector::conj_view<Coll> >
    {
	typedef typename vector::conj_view<Coll>::value_type        value_type;
	typedef typename vector::conj_view<Coll>::const_reference   const_reference;
	typedef typename vector::conj_view<Coll>::size_type         size_type;
    };
#else
    template <typename Coll>
    struct Collection<vector::conj_view<Coll> >
    {
	typedef typename vector::conj_view<Coll>::value_type        value_type;
	typedef typename vector::conj_view<Coll>::const_reference   const_reference;
	typedef typename vector::conj_view<Coll>::size_type         size_type;
    };
#endif

#ifdef __GXX_CONCEPTS__
#else
    template <typename Coll>
    struct Collection<vector::real_view<Coll> >
    {
	typedef typename vector::real_view<Coll>::value_type        value_type;
	typedef typename vector::real_view<Coll>::const_reference   const_reference;
	typedef typename vector::real_view<Coll>::size_type         size_type;
    };
#endif

#ifdef __GXX_CONCEPTS__
#else
    template <typename Coll>
    struct Collection<vector::imag_view<Coll> >
    {
	typedef typename vector::imag_view<Coll>::value_type        value_type;
	typedef typename vector::imag_view<Coll>::const_reference   const_reference;
	typedef typename vector::imag_view<Coll>::size_type         size_type;
    };
#endif

#ifdef __GXX_CONCEPTS__
    template <typename Coll>
    concept_map Collection<vector::negate_view<Coll> >
    {
	typedef typename vector::negate_view<Coll>::value_type        value_type;
	typedef typename vector::negate_view<Coll>::const_reference   const_reference;
	typedef typename vector::negate_view<Coll>::size_type         size_type;
    };
#else
    template <typename Coll>
    struct Collection<vector::negate_view<Coll> >
    {
	typedef typename vector::negate_view<Coll>::value_type        value_type;
	typedef typename vector::negate_view<Coll>::const_reference   const_reference;
	typedef typename vector::negate_view<Coll>::size_type         size_type;
    };
#endif

#ifdef __GXX_CONCEPTS__
#else
    template <class E1, class E2, typename SFunctor>
    struct Collection< vector::vec_scal_aop_expr<E1, E2, SFunctor> >
    {
	typedef typename Collection<E1>::value_type                   value_type;
	typedef value_type                                            const_reference;
	typedef typename Collection<E1>::size_type                    size_type; 
    };
#endif


#ifdef __GXX_CONCEPTS__
#else
    template <class E1, class E2>
    struct Collection<vector::vec_scal_asgn_expr<E1, E2> >
      : Collection< vector::vec_scal_aop_expr<E1, E2, mtl::sfunctor::assign<typename Collection<E1>::value_type, E2> > >
    {};
#endif


#ifdef __GXX_CONCEPTS__
#else
    template <unsigned BSize, typename Vector>
    struct Collection<vector::unrolled1<BSize, Vector> >
    {
	typedef typename Collection<Vector>::value_type value_type;
	typedef value_type                              const_reference;
	typedef typename Collection<Vector>::size_type  size_type;
    };
#endif




#ifdef __GXX_CONCEPTS__
    template <typename Coll>
    concept_map Collection<matrix::conj_view<Coll> >
    {
	typedef typename matrix::conj_view<Coll>::value_type        value_type;
	typedef typename matrix::conj_view<Coll>::const_reference   const_reference;
	typedef typename matrix::conj_view<Coll>::size_type         size_type;
    };
#else
    template <typename Coll>
    struct Collection<matrix::conj_view<Coll> >
    {
	typedef typename matrix::conj_view<Coll>::value_type        value_type;
	typedef typename matrix::conj_view<Coll>::const_reference   const_reference;
	typedef typename matrix::conj_view<Coll>::size_type         size_type;
    };
#endif




#ifdef __GXX_CONCEPTS__
    template <typename Coll>
    concept_map Collection<matrix::imag_view<Coll> >
    {
	typedef typename matrix::imag_view<Coll>::value_type        value_type;
	typedef typename matrix::imag_view<Coll>::const_reference   const_reference;
	typedef typename matrix::imag_view<Coll>::size_type         size_type;
    };
#else
    template <typename Coll>
    struct Collection<matrix::imag_view<Coll> >
    {
	typedef typename matrix::imag_view<Coll>::value_type        value_type;
	typedef typename matrix::imag_view<Coll>::const_reference   const_reference;
	typedef typename matrix::imag_view<Coll>::size_type         size_type;
    };
#endif




#ifdef __GXX_CONCEPTS__
    template <typename Coll>
    concept_map Collection<matrix::negate_view<Coll> >
    {
	typedef typename matrix::negate_view<Coll>::value_type        value_type;
	typedef typename matrix::negate_view<Coll>::const_reference   const_reference;
	typedef typename matrix::negate_view<Coll>::size_type         size_type;
    };
#else
    template <typename Coll>
    struct Collection<matrix::negate_view<Coll> >
    {
	typedef typename matrix::negate_view<Coll>::value_type        value_type;
	typedef typename matrix::negate_view<Coll>::const_reference   const_reference;
	typedef typename matrix::negate_view<Coll>::size_type         size_type;
    };
#endif




#ifdef __GXX_CONCEPTS__
    template <typename Coll>
    concept_map Collection<matrix::real_view<Coll> >
    {
	typedef typename matrix::real_view<Coll>::value_type        value_type;
	typedef typename matrix::real_view<Coll>::const_reference   const_reference;
	typedef typename matrix::real_view<Coll>::size_type         size_type;
    };
#else
    template <typename Coll>
    struct Collection<matrix::real_view<Coll> >
    {
	typedef typename matrix::real_view<Coll>::value_type        value_type;
	typedef typename matrix::real_view<Coll>::const_reference   const_reference;
	typedef typename matrix::real_view<Coll>::size_type         size_type;
    };
#endif


#ifdef __GXX_CONCEPTS__
    template <typename Functor>
    concept_map Collection<matrix::implicit_dense<Functor> >
    {
	typedef typename matrix::implicit_dense<Functor>::value_type        value_type;
	typedef typename matrix::implicit_dense<Functor>::const_reference   const_reference;
	typedef typename matrix::implicit_dense<Functor>::size_type         size_type;
    };
#else
    template <typename Functor>
    struct Collection<matrix::implicit_dense<Functor> >
    {
	typedef typename matrix::implicit_dense<Functor>::value_type        value_type;
	typedef typename matrix::implicit_dense<Functor>::const_reference   const_reference;
	typedef typename matrix::implicit_dense<Functor>::size_type         size_type;
    };
#endif


#ifdef __GXX_CONCEPTS__
    template <typename Value>
    concept_map Collection<matrix::ones_matrix<Value> >
    {
	typedef typename matrix::ones_matrix<Value>::value_type        value_type;
	typedef typename matrix::ones_matrix<Value>::const_reference   const_reference;
	typedef typename matrix::ones_matrix<Value>::size_type         size_type;
    };
#else
    template <typename Value>
    struct Collection<matrix::ones_matrix<Value> >
    {
	typedef typename matrix::ones_matrix<Value>::value_type        value_type;
	typedef typename matrix::ones_matrix<Value>::const_reference   const_reference;
	typedef typename matrix::ones_matrix<Value>::size_type         size_type;
    };
#endif


#ifdef __GXX_CONCEPTS__
    template <typename Value>
    concept_map Collection<matrix::hilbert_matrix<Value> >
    {
	typedef typename matrix::hilbert_matrix<Value>::value_type        value_type;
	typedef typename matrix::hilbert_matrix<Value>::const_reference   const_reference;
	typedef typename matrix::hilbert_matrix<Value>::size_type         size_type;
    };
#else
    template <typename Value>
    struct Collection<matrix::hilbert_matrix<Value> >
    {
	typedef typename matrix::hilbert_matrix<Value>::value_type        value_type;
	typedef typename matrix::hilbert_matrix<Value>::const_reference   const_reference;
	typedef typename matrix::hilbert_matrix<Value>::size_type         size_type;
    };
#endif


#ifdef __GXX_CONCEPTS__
    template <typename Vector1, typename Vector2>
    concept_map Collection<matrix::outer_product_matrix<Vector1, Vector2> >
    {
	typedef typename matrix::outer_product_matrix<Vector1, Vector2>::value_type        value_type;
	typedef typename matrix::outer_product_matrix<Vector1, Vector2>::const_reference   const_reference;
	typedef typename matrix::outer_product_matrix<Vector1, Vector2>::size_type         size_type;
    };
#else
    template <typename Vector1, typename Vector2>
    struct Collection<matrix::outer_product_matrix<Vector1, Vector2> >
    {
	typedef typename matrix::outer_product_matrix<Vector1, Vector2>::value_type        value_type;
	typedef typename matrix::outer_product_matrix<Vector1, Vector2>::const_reference   const_reference;
	typedef typename matrix::outer_product_matrix<Vector1, Vector2>::size_type         size_type;
    };
#endif




#ifdef __GXX_CONCEPTS__
    template <typename Matrix>
    concept_map Collection<mtl::matrix::transposed_view<Matrix> >
    {
	typedef typename mtl::matrix::transposed_view<Matrix>::value_type        value_type;
	typedef typename mtl::matrix::transposed_view<Matrix>::const_reference   const_reference;
	typedef typename mtl::matrix::transposed_view<Matrix>::size_type         size_type;
    };
#else
    template <typename Matrix>
    struct Collection<mtl::matrix::transposed_view<Matrix> >
    {
	typedef typename mtl::matrix::transposed_view<Matrix>::value_type        value_type;
	typedef typename mtl::matrix::transposed_view<Matrix>::const_reference   const_reference;
	typedef typename mtl::matrix::transposed_view<Matrix>::size_type         size_type;
    };
#endif


#ifdef __GXX_CONCEPTS__
    template <typename Matrix>
    concept_map Collection<matrix::hermitian_view<Matrix> >
    {
	typedef typename Collection<Matrix>::value_type        value_type;
	typedef typename Collection<Matrix>::const_reference   const_reference;
	typedef typename Collection<Matrix>::size_type         size_type;
    };
#else
    template <typename Matrix>
    struct Collection<matrix::hermitian_view<Matrix> >
    {
	typedef typename Collection<Matrix>::value_type        value_type;
	typedef typename Collection<Matrix>::const_reference   const_reference;
	typedef typename Collection<Matrix>::size_type         size_type;
    };
#endif

#ifdef __GXX_CONCEPTS__
    template <typename Coll>
    concept_map Collection< matrix::banded_view<Coll> >
    {
	typedef typename matrix::banded_view<Coll>::value_type        value_type;
	typedef typename matrix::banded_view<Coll>::const_reference   const_reference;
	typedef typename matrix::banded_view<Coll>::size_type         size_type;
    };
#else
    template <typename Coll>
    struct Collection< matrix::banded_view<Coll> >
    {
	typedef typename matrix::banded_view<Coll>::value_type        value_type;
	typedef typename matrix::banded_view<Coll>::const_reference   const_reference;
	typedef typename matrix::banded_view<Coll>::size_type         size_type;
    };
#endif

#ifdef __GXX_CONCEPTS__
#else
    template <typename Matrix>
    struct Collection< matrix::indirect<Matrix> >
    {
	typedef typename Collection<Matrix>::value_type               value_type;
	typedef value_type                                            const_reference; // maybe change later
	typedef typename Collection<Matrix>::size_type                size_type;
    };
#endif


#ifdef __GXX_CONCEPTS__

#else
    template <typename Matrix, typename Tag, int level>
    struct Collection<traits::detail::sub_matrix_cursor<Matrix, Tag, level> >
    {
	typedef typename Collection<Matrix>::value_type               value_type;
	typedef typename Collection<Matrix>::const_reference          const_reference;
	typedef typename Collection<Matrix>::size_type                size_type;
    };
#endif

#ifdef __GXX_CONCEPTS__

#else
    template <typename E1, typename E2>
    struct Collection<matrix::mat_mat_times_expr<E1, E2> >
    {
	typedef typename Collection<E1>::value_type                   ft;
	typedef typename Collection<E2>::value_type                   st;

	typedef typename Multiplicable<ft, st>::result_type           value_type;
	typedef const value_type&                                     const_reference;
	typedef typename Collection<E1>::size_type                    size_type;
    };
#endif

#ifdef __GXX_CONCEPTS__

#else
    template <typename E1, typename E2>
    struct Collection<matrix::mat_mat_plus_expr<E1, E2> >
    {
	typedef typename Collection<E1>::value_type                   ft;
	typedef typename Collection<E2>::value_type                   st;

	typedef typename Addable<ft, st>::result_type                 value_type;
	typedef const value_type&                                     const_reference;
	typedef typename Collection<E1>::size_type                    size_type;
    };
#endif

#ifdef __GXX_CONCEPTS__

#else
    template <typename E1, typename E2>
    struct Collection<matrix::mv_mv_plus_expr<E1, E2> >
      : Collection<matrix::mat_mat_plus_expr<E1, E2> >
    {};
#endif

#ifdef __GXX_CONCECPTS__

#else
    template< typename Matrix, typename Vector>
    struct Collection<mtl::mat_cvec_times_expr<Matrix, Vector> > 
    {
	typedef typename Collection< Matrix >::value_type value_type;
	typedef const value_type& const_reference;
	typedef typename Collection< Matrix >::size_type size_type;
    };
#endif


#ifdef __GXX_CONCEPTS__

#else
    template <typename Value, typename Parameters>
    struct Collection<mtl::vector::dense_vector<Value, Parameters> >
    {
	typedef Value            value_type;
	typedef const Value&     const_reference;
	typedef typename Parameters::size_type size_type;
    };
#endif

#ifdef __GXX_CONCEPTS__

#else
    template <typename Value, typename Parameters>
    struct Collection<mtl::vector::strided_vector_ref<Value, Parameters> >
    {
	typedef typename boost::remove_const<Value>::type            value_type;
	typedef const Value&                                         const_reference;
	typedef typename mtl::vector::strided_vector_ref<Value, Parameters>::size_type size_type;
    };
#endif

#ifdef __GXX_CONCEPTS__
#else
    template <typename Value, typename Parameters>
    struct Collection<mtl::vector::sparse_vector<Value, Parameters> >
    {
	typedef Value            value_type;
	typedef value_type       const_reference;
	typedef typename Parameters::size_type size_type;
    };
#endif


#ifdef __GXX_CONCEPTS__
#else
    template <class E1, class E2, typename SFunctor> 
    struct Collection<mtl::vector::vec_vec_ele_prod_expr<E1, E2, SFunctor> >
    {
	typedef typename Multiplicable<typename Collection<E1>::value_type,
				       typename Collection<E2>::value_type>::result_type value_type;
	typedef const value_type&     const_reference;
	typedef typename Collection<E1>::size_type  size_type;
    };
#endif

#ifdef __GXX_CONCEPTS__
#else
    template <class E1, class E2, typename SFunctor> 
    struct Collection<mtl::vector::vec_vec_pmop_expr<E1, E2, SFunctor> >
    {
	typedef typename Addable<typename Collection<E1>::value_type,
				 typename Collection<E2>::value_type>::result_type value_type;
	typedef const value_type&     const_reference;
	typedef typename Collection<E1>::size_type  size_type;
    };
#endif

#ifdef __GXX_CONCEPTS__
#else
    template <class E1, class E2, typename SFunctor> 
    struct Collection<mtl::vector::vec_vec_op_expr<E1, E2, SFunctor> >
    {
	typedef typename SFunctor::result_type value_type;
	typedef const value_type&     const_reference;
	typedef typename Collection<E1>::size_type  size_type;
    };
#endif

#ifdef __GXX_CONCEPTS__
#else
    template <class E1, class E2, typename SFunctor> 
    struct Collection<mtl::vector::vec_vec_aop_expr<E1, E2, SFunctor> >
    {
	typedef typename Collection<E1>::value_type value_type;
	typedef const value_type&                   const_reference;
	typedef typename Collection<E1>::size_type  size_type;
    };
#endif

#ifdef __GXX_CONCEPTS__
#else
    template <typename Matrix, typename VectorIn>
    struct Collection<mtl::vector::mat_cvec_multiplier<Matrix, VectorIn> >
    {
	typedef typename Multiplicable<typename Collection<Matrix>::value_type,
				       typename Collection<VectorIn>::value_type>::result_type value_type;
	typedef const value_type&                   const_reference;
	typedef typename Collection<Matrix>::size_type  size_type;
    };
#endif


#ifdef __GXX_CONCEPTS__

#else
    // Dunno if this is really a good idea
    template <typename T>
    struct Collection<T const>
      : Collection<T>
    {};
#endif


#ifdef __GXX_CONCEPTS__

#else
    template <typename Value>
    struct Collection<std::vector<Value> >
    {
	typedef typename std::vector<Value>::value_type      value_type;
	typedef typename std::vector<Value>::const_reference const_reference;
	typedef typename std::vector<Value>::size_type       size_type;
    };
#endif

#ifdef __GXX_CONCEPTS__
#else
    template <typename Vector, typename Functor>
    struct Collection<mtl::vector::lazy_reduction<Vector, Functor> >
    {
	typedef std::size_t size_type;
    };
#endif



#ifdef __GXX_CONCEPTS__
    template <typename Value, typename Parameters>
    concept_map MutableCollection<mtl::matrix::dense2D<Value, Parameters> >
    {
	typedef Value            value_type;
	typedef const Value&     const_reference;
	typedef typename mtl::matrix::dense2D<Value, Parameters>::size_type size_type;

	typedef Value&           reference;
    };
#else
    template <typename Value, typename Parameters>
    struct MutableCollection<mtl::matrix::dense2D<Value, Parameters> >
	: public Collection<mtl::matrix::dense2D<Value, Parameters> >
    {
	typedef Value&           reference;
    };
#endif

#ifdef __GXX_CONCEPTS__

    template <typename Value, unsigned long Mask, typename Parameters>
    concept_map MutableCollection<mtl::matrix::morton_dense<Value, Mask, Parameters> >
    {
	typedef Value            value_type;
	typedef const Value&     const_reference;
	typedef typename mtl::matrix::morton_dense<Value, Mask, Parameters>::size_type size_type;

	typedef Value&           reference;
    };

#else

    template <typename Value, unsigned long Mask, typename Parameters>
    struct MutableCollection<mtl::matrix::morton_dense<Value, Mask, Parameters> >
	: public Collection<mtl::matrix::morton_dense<Value, Mask, Parameters> >
    {
	typedef Value&           reference;
    };

#endif


#ifdef __GXX_CONCEPTS__
    template <typename Value, typename Parameters>
    concept_map MutableCollection<mtl::vector::strided_vector_ref<Value, Parameters> >
    {
	typedef typename boost::remove_const<Value>::type            value_type;
	typedef const Value&     const_reference;
	typedef typename mtl::vector::strided_vector_ref<Value, Parameters>::size_type size_type;

	typedef Value&           reference;
    };
#else
    template <typename Value, typename Parameters>
    struct MutableCollection<mtl::vector::strided_vector_ref<Value, Parameters> >
	: public Collection<mtl::vector::strided_vector_ref<Value, Parameters> >
    {
	typedef Value&           reference;
    };
#endif



#ifdef __GXX_CONCEPTS__
    template <typename Value>
    concept_map MutableCollection<<std::vector<Value> >
    {
	typedef typename std::vector<Value>::value_type      value_type;
	typedef typename std::vector<Value>::const_reference const_reference;
	typedef typename std::vector<Value>::size_type       size_type;

	typedef typename std::vector<Value>::reference       reference;
    };
#else
    template <typename Value>
    struct MutableCollection<std::vector<Value> >
	: public Collection<std::vector<Value> >
    {
	typedef typename std::vector<Value>::reference       reference;
    };
#endif

#ifdef __GXX_CONCEPTS__
#else
    template <typename T>
    struct OrientedCollection<const T>
      : OrientedCollection<T>
    {};
#endif



#ifdef __GXX_CONCEPTS__
    template <typename Value, typename Parameters>
    concept_map OrientedCollection<mtl::matrix::dense2D<Value, Parameters> >
    {
	typedef Value            value_type;
	typedef const Value&     const_reference;
	typedef typename mtl::matrix::dense2D<Value, Parameters>::size_type size_type;

	typedef typename mtl::matrix::dense2D<Value, Parameters>::orientation   orientation;
    };
#else
    template <typename Value, typename Parameters>
    struct OrientedCollection<mtl::matrix::dense2D<Value, Parameters> >
	: public Collection<mtl::matrix::dense2D<Value, Parameters> >
    {
	typedef typename mtl::matrix::dense2D<Value, Parameters>::orientation   orientation;
    };
#endif

#ifdef __GXX_CONCEPTS__

    template <typename Value, unsigned long Mask, typename Parameters>
    concept_map OrientedCollection<mtl::matrix::morton_dense<Value, Mask, Parameters> >
    {
	typedef Value            value_type;
	typedef const Value&     const_reference;
	typedef typename mtl::matrix::morton_dense<Value, Mask, Parameters>::size_type size_type;

	typedef typename mtl::matrix::morton_dense<Value, Mask, Parameters>::orientation   orientation;
    };

#else

    template <typename Value, unsigned long Mask, typename Parameters>
    struct OrientedCollection<mtl::matrix::morton_dense<Value, Mask, Parameters> >
	: public Collection<mtl::matrix::morton_dense<Value, Mask, Parameters> >
    {
	typedef typename mtl::matrix::morton_dense<Value, Mask, Parameters>::orientation   orientation;
    };

#endif

#ifdef __GXX_CONCEPTS__
    template <typename Value, typename Parameters>
    concept_map OrientedCollection<mtl::matrix::compressed2D<Value, Parameters> >
    {
	typedef Value            value_type;
	typedef const Value&     const_reference;
	typedef typename mtl::matrix::compressed2D<Value, Parameters>::size_type size_type;

	typedef typename mtl::matrix::compressed2D<Value, Parameters>::orientation   orientation;
    };
#else
    template <typename Value, typename Parameters>
    struct OrientedCollection<mtl::matrix::compressed2D<Value, Parameters> >
	: public Collection<mtl::matrix::compressed2D<Value, Parameters> >
    {
	typedef typename mtl::matrix::compressed2D<Value, Parameters>::orientation   orientation;
    };
#endif

#ifdef __GXX_CONCEPTS__
#else
    template <typename Matrix>
    struct OrientedCollection<mtl::matrix::indirect<Matrix> >
	: public Collection<mtl::matrix::indirect<Matrix> >
    {
	typedef row_major   orientation;
    };
#endif

#ifdef __GXX_CONCEPTS__
    template <typename Value, typename Parameters>
    concept_map OrientedCollection<mtl::vector::dense_vector<Value, Parameters> >
    {
	typedef Value            value_type;
	typedef const Value&     const_reference;
	typedef typename mtl::vector::dense_vector<Value, Parameters>::size_type size_type;

	typedef typename mtl::vector::dense_vector<Value, Parameters>::orientation   orientation;
    };
#else
    template <typename Value, typename Parameters>
    struct OrientedCollection<mtl::vector::dense_vector<Value, Parameters> >
	: public Collection<mtl::vector::dense_vector<Value, Parameters> >
    {
	typedef typename mtl::vector::dense_vector<Value, Parameters>::orientation   orientation;
    };
#endif

#ifdef __GXX_CONCEPTS__
    template <typename Value, typename Parameters>
    concept_map OrientedCollection<mtl::vector::strided_vector_ref<Value, Parameters> >
    {
	typedef typename boost::remove_const<Value>::type            value_type;
	typedef const Value&     const_reference;
	typedef typename mtl::vector::strided_vector_ref<Value, Parameters>::size_type size_type;

	typedef typename mtl::vector::strided_vector_ref<Value, Parameters>::orientation   orientation;
    };
#else
    template <typename Value, typename Parameters>
    struct OrientedCollection<mtl::vector::strided_vector_ref<Value, Parameters> >
	: public Collection<mtl::vector::strided_vector_ref<Value, Parameters> >
    {
	typedef typename mtl::vector::strided_vector_ref<Value, Parameters>::orientation   orientation;
    };
#endif



#ifdef __GXX_CONCEPTS__
    template <typename Value>
    concept_map OrientedCollection<std::vector<Value> >
    {
	typedef Value            value_type;
	typedef const Value&     const_reference;
	typedef typename std::vector<Value>::size_type size_type;

	typedef mtl::tag::col_major   orientation;
    };
#else
    template <typename Value>
    struct OrientedCollection<std::vector<Value> >
	: public Collection<std::vector<Value> >
    {
	typedef mtl::tag::col_major   orientation;
    };
#endif


#ifdef __GXX_CONCEPTS__
    template <typename Scaling, typename Coll>
    concept_map OrientedCollection< matrix::scaled_view<Scaling, Coll> >
    {
	typedef typename matrix::scaled_view<Scaling, Coll>::value_type        value_type;
	typedef typename matrix::scaled_view<Scaling, Coll>::const_reference   const_reference;
	typedef typename matrix::scaled_view<Scaling, Coll>::size_type         size_type;

	typedef typename OrientedCollection<Coll>::orientation       orientation;
    };
#else
    template <typename Scaling, typename Coll>
    struct OrientedCollection< matrix::scaled_view<Scaling, Coll> >
	: public Collection< matrix::scaled_view<Scaling, Coll> >
    {
	typedef typename OrientedCollection<Coll>::orientation       orientation;
    };
#endif

// added by Hui Li
#ifdef __GXX_CONCEPTS__
    template <typename Coll, typename RScaling>
    concept_map OrientedCollection< matrix::rscaled_view<Coll,RScaling> >
    {
	typedef typename matrix::rscaled_view<Coll,RScaling>::value_type        value_type;
	typedef typename matrix::rscaled_view<Coll,RScaling>::const_reference   const_reference;
	typedef typename matrix::rscaled_view<Coll,RScaling>::size_type         size_type;

	typedef typename OrientedCollection<Coll>::orientation       orientation;
    };
#else
    template <typename Coll, typename RScaling>
    struct OrientedCollection< matrix::rscaled_view<Coll,RScaling> >
	: public Collection< matrix::rscaled_view<Coll,RScaling> >
    {
	typedef typename OrientedCollection<Coll>::orientation       orientation;
    };
#endif

	
#ifdef __GXX_CONCEPTS__
    template <typename Scaling, typename Coll>
    concept_map OrientedCollection< mtl::vector::scaled_view<Scaling, Coll> >
    {
	typedef typename mtl::vector::scaled_view<Scaling, Coll>::value_type        value_type;
	typedef typename mtl::vector::scaled_view<Scaling, Coll>::const_reference   const_reference;
	typedef typename mtl::vector::scaled_view<Scaling, Coll>::size_type         size_type;

	typedef typename OrientedCollection<Coll>::orientation       orientation;
    };
#else
    template <typename Scaling, typename Coll>
    struct OrientedCollection< mtl::vector::scaled_view<Scaling, Coll> >
	: public Collection< mtl::vector::scaled_view<Scaling, Coll> >
    {
	typedef typename OrientedCollection<Coll>::orientation       orientation;
    };
#endif

//added by Hui Li
#ifdef __GXX_CONCEPTS__
    template <typename Coll, typename RScaling>
    concept_map OrientedCollection< mtl::vector::rscaled_view<Coll,RScaling> >
    {
	typedef typename mtl::vector::rscaled_view<Coll,RScaling>::value_type        value_type;
	typedef typename mtl::vector::rscaled_view<Coll,RScaling>::const_reference   const_reference;
	typedef typename mtl::vector::rscaled_view<Coll,RScaling>::size_type         size_type;

	typedef typename OrientedCollection<Coll>::orientation       orientation;
    };
#else
    template <typename Coll, typename RScaling>
    struct OrientedCollection< mtl::vector::rscaled_view<Coll,RScaling> >
	: public Collection< mtl::vector::rscaled_view<Coll,RScaling> >
    {
	typedef typename OrientedCollection<Coll>::orientation       orientation;
    };
#endif

	
#ifdef __GXX_CONCEPTS__
    template <typename Coll>
    concept_map OrientedCollection<matrix::conj_view<Coll> >
    {
	typedef typename matrix::conj_view<Coll>::value_type         value_type;
	typedef typename matrix::conj_view<Coll>::const_reference    const_reference;
	typedef typename matrix::conj_view<Coll>::size_type          size_type;

	typedef typename OrientedCollection<Coll>::orientation       orientation;
    };
#else
    template <typename Coll>
    struct OrientedCollection<matrix::conj_view<Coll> >
	: public Collection<matrix::conj_view<Coll> >
    {
	typedef typename OrientedCollection<Coll>::orientation       orientation;
    };
#endif


#ifdef __GXX_CONCEPTS__
    template <typename Coll>
    concept_map OrientedCollection<mtl::vector::conj_view<Coll> >
    {
	typedef typename mtl::vector::conj_view<Coll>::value_type        value_type;
	typedef typename mtl::vector::conj_view<Coll>::const_reference   const_reference;
	typedef typename mtl::vector::conj_view<Coll>::size_type         size_type;

	typedef typename OrientedCollection<Coll>::orientation       orientation;
    };
#else
    template <typename Coll>
    struct OrientedCollection<mtl::vector::conj_view<Coll> >
	: public Collection<mtl::vector::conj_view<Coll> >
    {
	typedef typename OrientedCollection<Coll>::orientation       orientation;
    };
#endif


#ifdef __GXX_CONCEPTS__
    template <typename Functor, typename Coll>
    concept_map OrientedCollection< mtl::vector::map_view<Functor, Coll> >
    {
	typedef typename mtl::vector::map_view<Functor, Coll>::value_type        value_type;
	typedef typename mtl::vector::map_view<Functor, Coll>::const_reference   const_reference;
	typedef typename mtl::vector::map_view<Functor, Coll>::size_type         size_type;

	typedef typename OrientedCollection<Coll>::orientation       orientation;
    };
#else
    template <typename Functor, typename Coll>
    struct OrientedCollection< mtl::vector::map_view<Functor, Coll> >
	: public Collection< mtl::vector::map_view<Functor, Coll> >
    {
	typedef typename OrientedCollection<Coll>::orientation       orientation;
    };
#endif


#ifdef __GXX_CONCEPTS__
    template <typename Coll>
    concept_map OrientedCollection<mtl::matrix::transposed_view<Coll> >
    {
	typedef typename mtl::matrix::transposed_view<Coll>::value_type        value_type;
	typedef typename mtl::matrix::transposed_view<Coll>::const_reference   const_reference;
	typedef typename mtl::matrix::transposed_view<Coll>::size_type         size_type;

	typedef typename mtl::traits::transposed_orientation<typename OrientedCollection<Coll>::orientation>::type   orientation;
    };
#else
    template <typename Coll>
    struct OrientedCollection<mtl::matrix::transposed_view<Coll> >
	: public Collection<mtl::matrix::transposed_view<Coll> >
    {
	typedef typename mtl::traits::transposed_orientation<typename OrientedCollection<Coll>::orientation>::type   orientation;
    };
#endif

#ifdef __GXX_CONCEPTS__
    template <typename Coll>
    concept_map OrientedCollection<matrix::hermitian_view<Coll> >
    {
	typedef typename matrix::hermitian_view<Coll>::value_type        value_type;
	typedef typename matrix::hermitian_view<Coll>::const_reference   const_reference;
	typedef typename matrix::hermitian_view<Coll>::size_type         size_type;

	typedef typename mtl::traits::transposed_orientation<typename OrientedCollection<Coll>::orientation>::type   orientation;
    };
#else
    template <typename Coll>
    struct OrientedCollection<matrix::hermitian_view<Coll> >
      : public Collection<matrix::hermitian_view<Coll> >
    {
	typedef typename mtl::traits::transposed_orientation<typename OrientedCollection<Coll>::orientation>::type   orientation;
    };
#endif



/*@}*/ // end of group Concepts

} // namespace mtl

#endif // MTL_COLLECTION_INCLUDE
