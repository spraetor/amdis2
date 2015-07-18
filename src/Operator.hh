#include "Expressions.h"

namespace AMDiS 
{
  template <class C>
  void addZOTAux(C const& c)
  {
    GenericZeroOrderTerm<C>* term = new GenericZeroOrderTerm<C>(c)
    zeroOrder.push_back(term);
    term->operat = this;
    c.insertFeSpaces(auxFeSpaces);
  }
  
  
  template <class B>
  void addFOTAux(tag::scalar, B const& b, FirstOrderType type, int i)
  {
    GenericOperatorTerm<B, 1>* term;
    if (i >= 0)
      term = new GenericFirstOrderTerm_i<-1,B>(b, i)
    else
      term = new GenericFirstOrderTerm_1<B>(b);
    
    if (type == GRD_PSI)
      firstOrderGrdPsi.push_back(term);
    else
      firstOrderGrdPhi.push_back(term);
    
    term->operat = this;
    b.insertFeSpaces(auxFeSpaces);
  }
  
  
  template <class B>
  void addFOTAux(tag::vector, B const& b, FirstOrderType type, int i)
  {
    GenericOperatorTerm<B, 1>* term = new GenericFirstOrderTerm_b<B>(b);
    
    if (type == GRD_PSI)
      firstOrderGrdPsi.push_back(term);
    else
      firstOrderGrdPhi.push_back(term);
    
    term->operat = this;
    b.insertFeSpaces(auxFeSpaces);
  }
  
  
  template <class A>
  void addSOTAux(tag::scalar, A const& a, int i, int j, bool sym)
  {
    GenericOperatorTerm<A, 2>* term;
    if (i >= 0 && j >= 0)
      term = new GenericSecondOrderTerm_ij<-1,-1,A>(a, i, j)
    else
      term = new GenericSecondOrderTerm_1<A>(a)
    
    secondOrder.push_back(term);
    term->operat = this;
    a.insertFeSpaces(auxFeSpaces);
  }
  
  
  template <class A>
  void addSOTAux(tag::matrix, A const& a, int i, int j, bool sym)
  {
    GenericOperatorTerm<A, 2>* term;
    if (sym)
      term = new GenericSecondOrderTerm_A<A, true>(a);
    else
      term = new GenericSecondOrderTerm_A<A, false>(a);
    
    secondOrder.push_back(term);
    term->operat = this;
    a.insertFeSpaces(auxFeSpaces);
  }
  
} // end namespace AMDiS
