/** \file Initfile.h */

#pragma once

#include <string>
#include <iostream>

#include <boost/lexical_cast.hpp> 
#include <boost/numeric/conversion/cast.hpp> 

#include <boost/tokenizer.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include <boost/type_traits.hpp>

// a parser for arithmetic expressions
#include "muParser.h"
#include "traits/basic.hpp"
#include "traits/scalar_types.hpp"
#include "traits/meta_basic.hpp"

#include "Math.h"

namespace AMDiS 
{
  namespace detail 
  {   
    /// convert string to intrinsic type
    template <class T, class Enable = void>
    struct Convert
    {
      static void eval(std::string valStr, T& value)
      {
        value = boost::lexical_cast<T>(valStr);
      }      
    };
    
    template <class T>
    struct Convert<T, typename std::enable_if<traits::is_arithmetic<T>::value >::type>
    {
      static void eval(const std::string valStr, T& value)
      {
        using boost::numeric_cast;
  
        mu::Parser parser;
        parser.DefineConst(_T("M_PI"), m_pi);
        parser.DefineConst(_T("M_E"), m_e);
  
        parser.SetExpr(valStr);
        value = numeric_cast<T>(parser.Eval());
      }
    };
    
    // convert string to vector
    template <class T>
    struct Convert<T, typename std::enable_if<IsVector<T>::value>::type>
    {
      static void eval(const std::string valStr, T& value)
      {
        using value_type = Value_t<T>;
        using Tokenizer = boost::tokenizer<boost::char_separator<char> >;
        
        std::vector<value_type> values;
        boost::char_separator<char> sep(",; ");
        Tokenizer tokens(valStr, sep);
        int i = 0;
        for (auto token : tokens) {
          value_type v;
          Convert<value_type>::eval(token, v);
          value[i] = v;
        }
      }
    };
    
    // convert string to vector
    template <class T>
    struct Convert<std::vector<T> >
    {
      static void eval(const std::string valStr, std::vector<T>& values)
      {
        using value_type = T;
        using Tokenizer = boost::tokenizer<boost::char_separator<char> >;
        
        boost::char_separator<char> sep(",; ");
        Tokenizer tokens(valStr, sep);
        for (auto token : tokens) {
          value_type v;
          Convert<value_type>::eval(token, v);
          values.push_back(v);
        }
      }
    };

  } // end namespace detail


  /// output-stream for std::list
  template<typename T>
  std::ostream& operator<<(std::ostream& o, const std::list< T >& l)
  {
    typename std::list< T >::const_iterator it = l.begin();
    o << "[";
    for (unsigned i = 0; it != l.end() && i < l.size(); i++) {
      o << *it << (i < l.size() - 1 ? ", " : "");
      ++it;
    }
    o << "]";
    return o;
  }


  /// output-stream for std::vector
  template<typename T>
  std::ostream& operator<<(std::ostream& o, const std::vector<T>& l)
  {
    typename std::vector<T>::const_iterator it = l.begin();
    o << "[";
    for (unsigned i = 0; it != l.end() && i < l.size(); i++) {
      o << *it << (i < l.size() - 1 ? ", " : "");
      ++it;
    }
    o << "]";
    return o;
  }

// _____________________________________________________________________________

  /** Basis data container as a map of tag on a value as strings. The container 
  *  throws an exception, if the tag was not found.
  */
  struct Initfile
  {
    /** initialize init-file from file with filename in, read data and save it 
    *  to singleton-map
    *  @param in: filename string
    */
    static void init(std::string in);


    /** Static get routine for getting parameter-values from init-file 
    *  initialized in init()-method.
    *  Cast the value to the desired type using std::stringstream.
    *  @param tag: The tag to look for
    *  @param value: The result.
    *  @param debugInfo: msgInfo for current parameter. (0..no printing,
    *    1..print missing parameter info, 2..print parameter value) [optional]
    */
    template <class T>
    static void get(const std::string tag, T& value, int debugInfo = -1)
    {
      initIntern();
      if (debugInfo == -1)
        debugInfo = singlett->getMsgInfo();
      else {
      	int swap(debugInfo);
      	debugInfo = singlett->getMsgInfo();
      	singlett->msgInfo=swap;
      }

      // TODO: use boost::optional instead
      // TODO: use convert method from above
      std::string valueStr = "-";
      valueStr = singlett->pt.get(tag, valueStr);
      
      if (valueStr != "-")
        detail::Convert<T>::eval(valueStr, value);

      if (debugInfo == 2) {
      	std::cout << "Parameter '" << tag << "'"
      		        << " initialized with: " << value << std::endl;
      }
      singlett->msgInfo = debugInfo;
    }


    /// update map tag->value_old to tag->value in singleton
    template <class T>
    static void set(const std::string tag, T& value, int debugInfo=  -1) 
    {
      initIntern();
      if (debugInfo == -1)
        debugInfo = singlett->getMsgInfo();

      singlett->pt.put(tag, value);
      
      // update msg parameters msgInfo, msgWait, paramInfo
      singlett->getInternalParameters();
      if (debugInfo == 2)
        std::cout << "Parameter '" << tag << "'"
                  << " set to: " << value << std::endl;
    }


    /// add map tag->value to data in singleton
    template <class T>
    static void add(const std::string tag, T& value, int debugInfo = -1) 
    {
      set(tag, value, debugInfo);
    }

    /// Returns specified info level
    static int getMsgInfo() 
    { 
      return (singlett != NULL) ? singlett->msgInfo : 0; 
    }


    /// Returns specified wait value
    static int getMsgWait() 
    { 
      return (singlett != NULL) ? singlett->msgWait : 0; 
    }


    /// Checks whether parameters are initialized. if not, call init()
    static bool initialized() 
    { 
      return (singlett != NULL); 
    }


    /// return pointer to singleton
    static Initfile *getSingleton() 
    { 
      return singlett; 
    }


    /// print all data in singleton to std::cout
    static void printParameters();


    /// clear data in singleton
    static void clearData()
    {
      initIntern();
      // singlett->clear();
    }


    /// save singlett-data to file with filename fn
    static void save(std::string fn)
    {
      using namespace boost::property_tree;
      initIntern();
      // json_parser::write_jason(fn, singlett->pt);
    }
	
protected:	
    Initfile() 
      : msgInfo(0), 
        msgWait(1), 
        paramInfo(1), 
        breakOnMissingTag(0) 
    {}


    static void initIntern() 
    {
      if (singlett == NULL)
        singlett = new Initfile;
    }

    /// pointer to the singleton that contains the data
    static Initfile* singlett;

    /** Fill the initfile from an input stream.
    * @param in: the stream to fill the data from.
    * Current dataformat: tag:value
    * Comment char: percent '%'
    * Include files: #include "filename" or #include <filename>
    */
    void read(std::string fn, bool force = false);

    /// Write data-map to initfile with filename fn
    void write(std::string fn);

    /// read parameters for msgInfo, msgWait, paramInfo
    void getInternalParameters();

    int msgInfo, msgWait, paramInfo, breakOnMissingTag;	
    
    /// boost:property_tree to read/store parameter values
    boost::property_tree::ptree pt;
  };
  
  using Parameters = Initfile;
  
#ifndef AMDIS_NO_EXTERN_INITFILE
extern template void Initfile::get(const std::string, int&, int);
extern template void Initfile::get(const std::string, double&, int);
extern template void Initfile::get(const std::string, std::string&, int);
#endif

} // end namespace AMDiS
