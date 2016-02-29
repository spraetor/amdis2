#pragma once

#include "traits/basic.hpp"
#include "operations/functors.hpp"
#include "matrix_vector/ExprConcepts.hpp"
#include "expressions/ConstantTerms.hpp"

namespace AMDiS
{  
  namespace functors
  {
    template <class T>
    struct Component : FunctorBase
    {
      int getDegree(int d0, int /* d1 */) { return d0; }
      
      auto operator()(T const& vector, size_t idx) const RETURNS 
      ( 
        vector(idx)
      )
      
      auto operator()(T const& matrix, size_t r, size_t c) const RETURNS 
      ( 
        matrix(r, c) 
      )
    };
  }
  
  
  template <class Model>
  struct ComponentViewBase
  {  
  };
  
  // forward declaration
  template <class F, class Term1, class... Terms>
  struct FunctorTerm;
  
  // For voctor containers a component access methods is provided
  template <class Container, class Model>
  class ComponentView<Container, Model, 
                      Requires_t< concepts::VectorExpression<Container> > >
    : private ComponentViewBase<Model>
  {
    using size_type = Size_t<Container>;
    
    using ComponentFunctor = functors::Component<Container>;
    using IndexTerm        = RTConstant<size_type>;
    
    Model& self()
    {
      return static_cast<Model&>(*this);
    }
    Model const& self() const
    {
      return static_cast<Model const&>(*this);
    }
  public:    
    /// generate a component_view
    FunctorTerm<ComponentFunctor, Model, IndexTerm>
    operator[](size_type i) const
    {
      return {this->self(), IndexTerm(i)};
    }
    
    /// generate a component_view
    FunctorTerm<ComponentFunctor, Model, IndexTerm>
    operator()(size_type i) const
    {
      return {this->self(), IndexTerm(i)};
    }
  };
  
  // For matrix containers a component access methods is provided
  template <class Container, class Model>
  class ComponentView<Container, Model, 
                      Requires_t< concepts::MatrixExpression<Container> > >
    : private ComponentViewBase<Model>
  {
    using size_type = Size_t<Container>;
    
    using ComponentFunctor = functors::Component<Container>;
    using IndexTerm        = RTConstant<size_type>;
    
    Model& self()
    {
      return static_cast<Model&>(*this);
    }
    Model const& self() const
    {
      return static_cast<Model const&>(*this);
    }
  public:
    /// generate a component_view
    FunctorTerm<ComponentFunctor, Model, IndexTerm, IndexTerm>
    operator()(size_type r, size_type c) const
    {
      return {this->self(), IndexTerm(r), IndexTerm(c)};
    }
  };
  
} // end namespace AMDiS
