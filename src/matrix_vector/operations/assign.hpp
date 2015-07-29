#pragma once

namespace AMDiS {
  namespace assign {
    
    // -------------------------------------------------------------------------
    // unary assigners
    
    /// assign_constant(v) --> v = value
    template <class T, class S=T>
    struct value
    {
      typedef T result_type;
      value(S val = 0) : val(val) {}
      
      T& operator()(T& v) const { return (v = val); }
      
    private:
      S val;
    };
        
    /// assign_constant(v) --> v = value
    template <class T, class S, S Val>
    struct ct_value : value<T, S>
    {
      ct_value() : value<T, S>(Val) {}
    };
    
    template <class T, class S=T>
    struct min_value : value<T, S>
    {
      min_value() : value<T, S>(std::numeric_limits<S>::min()) {}
    };
    
    template <class T, class S=T>
    struct max_value : value<T, S>
    {
      max_value() : value<T, S>(std::numeric_limits<S>::max()) {}
    };
    
    /// add_constant(v) --> v += value
    template <class T, class S=T>
    struct plus_value
    {
      typedef T result_type;
      plus_value(S val) : val(val) {}
      
      T& operator()(T& v) const { return (v += val); }
      
    private:
      S val;
    };
    
    /// minus_constant(v) --> v -= value
    template <class T, class S=T>
    struct minus_value
    {
      typedef T result_type;
      minus_value(S val) : val(val) {}
      
      T& operator()(T& v) const { return (v -= val); }
      
    private:
      S val;
    };
    
    /// mult_constant(v) --> v *= value
    template <class T, class S=T>
    struct mult_value
    {
      typedef T result_type;
      mult_value(S val) : val(val) {}
      
      T& operator()(T& v) const { return (v *= val); }
      
    private:
      S val;
    };
    
    /// div_constant(v) --> v /= value
    template <class T, class S=T>
    struct div_value
    {
      typedef T result_type;
      div_value(S val) : val(val) {}
      
      T& operator()(T& v) { return (v /= val); }
      
    private:
      S val;
    };

    // -------------------------------------------------------------------------
    /// binary assigners
    
    
    /// functor for operator=
    template <class T>
    struct assign
    {
      typedef T result_type;
      
      static T& apply(T& v, T const& v0) { return (v = v0); }
      T& operator()(T& v, T const& v0) const { return (v = v0); }
    };

    /// functor for operator+=
    template <class T>
    struct plus
    {
      typedef T result_type;
      
      static T& apply(T& v, T const& v0) { return (v += v0); }
      T& operator()(T& v, T const& v0) const { return (v += v0); }
    };

    /// functor for operator*=
    template <class T>
    struct multiplies
    {
      typedef T result_type;
      
      static T& apply(T& v, T const& v0) { return (v *= v0); }
      T& operator()(T& v, T const& v0) { return (v *= v0); }
    };

    /// functor for operator/=
    template <class T>
    struct divides
    {
      typedef T result_type;
      
      static T& apply(T& v, T const& v0) { return (v /= v0); }
      T& operator()(T& v, T const& v0) { return (v /= v0); }
    };

    /// functor for v = max(v, v0)
    template <class T>
    struct max
    {
      typedef T result_type;
      
      static T& apply(T& v, T const& v0) { return v = std::max(v, v0); }
      T& operator()(T& v, T const& v0) { return v = std::max(v, v0); }
    };

    /// functor for v = min(v, v0)
    template <class T>
    struct min
    {
      typedef T result_type;
      
      static T& apply(T& v, T const& v0) { return v = std::min(v, v0); }
      T& operator()(T& v, T const& v0) { return v = std::min(v, v0); }
    };

  } // end namespace assign
} // end namespace AMDiS
