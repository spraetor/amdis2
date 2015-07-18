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

#ifndef INITFILE_H
#define INITFILE_H

#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <list>
#include <set>
#include <vector>
#include <stdexcept>
#include <iostream>
#include <typeinfo>
#include "FixVec.h"

#include "traits/size.hpp"

#include <boost/lexical_cast.hpp>
#include <boost/numeric/conversion/cast.hpp> 

#include <boost/type_traits.hpp>

// a parser for arithmetic expressions
#include "muParser.h"

namespace AMDiS {

  namespace detail {
  
    // Exceptions
    //_________________________________________________________________________________________
    
    /// Exception for wrong vector size when parsing for concrete vector type
    struct WrongVectorSize : std::runtime_error {
      WrongVectorSize(std::string m) : std::runtime_error(m) {}
    };


    /// Exception for no delimiter found in vector
    struct NoDelim : std::runtime_error {
      NoDelim(std::string m) : std::runtime_error(m) {}
    };


    /// Exception for begin and end brackets are different in vector
    struct WrongVectorFormat : std::runtime_error {
      WrongVectorFormat(std::string m) : std::runtime_error(m) {}
    };
    
    
    /// Exception for wrong value format
    template<typename T>
    struct WrongValueFormat : std::runtime_error 
    {  
      static std::string name(bool) { return "bool"; }
      static std::string name(double) { return "double"; }
      static std::string name(float) { return "float"; }
      static std::string name(int) { return "int"; }
      static std::string name(unsigned int) { return "unsigned int"; }

      template<typename G>
      static std::string name(G) 
      { 
        return std::string(typeid(G).name()); 
      }

      WrongValueFormat(std::string value)
      : std::runtime_error("cannot convert '" + value + "' into <" + name(T()) + ">")
      {}
    };
    

    /// Exception for bad arithmetic expression that can not be evaluated
    template<typename T>
    struct BadArithmeticExpression : std::runtime_error 
    {
      static std::string name(bool) { return "bool"; }
      static std::string name(double) { return "double"; }
      static std::string name(float) { return "float"; }
      static std::string name(int) { return "int"; }
      static std::string name(unsigned int) { return "unsigned int"; }

      template<typename G>
      static std::string name(G) 
      { 
        return std::string(typeid(G).name());
      }

      BadArithmeticExpression(std::string m, std::string value) 
      : std::runtime_error("cannot evaluate expression '" + value + "' into <" + name(T()) + ">\n"
			   "Parser message: '" + m + "'")
      {}
    };
      
// _____________________________________________________________________________

    /// return the delimiter or throw an exception if there is no known 
    /// delimiter in value
    inline size_t checkDelim(const std::string& value, const std::string& delims)
    {
      size_t pos(std::string::npos);
      for (size_t i = 0; i < delims.length(); i++) {
        pos = value.find(delims[i]);
        if (pos != std::string::npos)
          return i;
      }
      //      throw NoDelim("cannot detect the delimiter in " + value);
      return 0; 
    }


    /// convert string to string
    inline void convert(const std::string valStr, std::string& value) 
    {
      value = trim(valStr);
    }


    /// convert string to intrinsic type
    template<typename T> inline
    typename boost::enable_if
	     < boost::mpl::and_
	       < boost::is_pod<T>, 
	         boost::mpl::not_< boost::is_enum<T> > >, 
	       void
	     >::type
    convert(const std::string valStr, T& value)
    {
      using boost::lexical_cast;
      using boost::numeric_cast;

      mu::Parser parser;
      parser.DefineConst(_T("M_PI"), m_pi);
      parser.DefineConst(_T("M_E"), m_e);

      //      try {
        parser.SetExpr(valStr);
        value = numeric_cast< T >(parser.Eval());
/*       } catch (boost::bad_lexical_cast e) { */
/*         throw WrongValueFormat< T >(valStr); */
/*       } catch (boost::bad_numeric_cast e) { */
/*         throw WrongValueFormat< T >(valStr); */
/*       } catch (mu::Parser::exception_type &e) { */
/*         throw BadArithmeticExpression<T>(e.GetMsg(), valStr); */
/*       } */
    }


    template<typename T> inline
    typename boost::enable_if
	     < boost::is_enum<T>, 
	       void
	     >::type
    convert(const std::string valStr, T& value) 
    {
      int swap = 0;
      //      try {
        swap = boost::lexical_cast<int>(trim(valStr));
/*       } catch (boost::bad_lexical_cast e) { */
/*         throw WrongValueFormat< T >(valStr); */
/*       } */
      value = static_cast< T >(swap);
    }


    /// convert special enums
    inline void convert(const std::string valStr, Norm& value) 
    {
      std::string swapStr = boost::to_upper_copy(valStr);

      if (swapStr == "NO_NORM")
        value = static_cast< Norm >(NO_NORM);
      else if (swapStr == "H1_NORM")
        value = static_cast< Norm >(H1_NORM);
      else if (swapStr == "L2_NORM")
        value = static_cast< Norm >(L2_NORM);
      else {
        int swap = 0;
        convert(valStr, swap);
        value = static_cast< Norm >(swap);
      }
    }


    /// convert value of arbitrary type to string using stringstream and 
    /// operator<< for type
    template<typename T>
    inline void convert(const T value, std::string& valStr) 
    {
      std::stringstream ss;
      ss.precision(6);
      ss << value;
      valStr = ss.str();
    }


    /// convert WorldVector to string
    // TODO: allgemeine Funktion zum Schreiben von Vektoren implementieren
    //       unter Ausnutzung von Type-Traits
    template<typename T>
    inline void convert(const WorldVector<T>& c, std::string& valStr)
    {
      std::vector<T> temp_vec(size(c));
      for (unsigned i = 0; i < size(temp_vec); i++)
        temp_vec[i] = c[i];
      convert(temp_vec, valStr);
    }

    // forward declarations
    template< typename T >
    inline void convert(const std::string valStr, WorldVector<T>& c);

    template<typename T>
    inline void convert(const std::string valStr, std::list<T>& value);

    template<typename T>
    inline void convert(const std::string valStr, std::vector<T>& value);

    /** parse an container from tag tag. The Container must have the properties:
    * 	- type value_type
    * 	- member function push_back
    */
    template< typename Container >
    inline void getContainer(const std::string val_, Container& c)
    {
      // accepted brackets and delimiters for vector input
      std::string begBrackets= "{[(";
      std::string endBrackets= "}])";
      std::string delims= ";,";

      c.clear();
      std::string val = trim(val_);
      bool hasBrackets = true;
      size_t pos = begBrackets.find(val[0]);
      if (pos == std::string::npos)
        hasBrackets = false;
       if (hasBrackets && val[val.length() - 1] != endBrackets[pos])
        throw WrongVectorFormat("begin and end bracket are different in"
            " value '" + val + "'"); 
      size_t oldPos = (hasBrackets ? 1 : 0);
      size_t curDelim = 0;
      typedef typename Container::value_type ValueType;
      ValueType swap;
      try {
        curDelim = checkDelim(val, delims);
        pos = val.find(delims[curDelim], oldPos);
        while (pos != std::string::npos) {
          std::string curWord = val.substr(oldPos, pos - oldPos);
          oldPos = pos + 1;
          convert(curWord, swap);
          c.push_back(swap);
          pos = val.find(delims[curDelim], oldPos);
        }
        //last entry
        std::string curWord = val.substr(oldPos, val.length() - (hasBrackets ? 1 : 0) - oldPos);
        convert(curWord, swap);
        c.push_back(swap);
      } catch (NoDelim nd) {
        std::string curWord = val.substr(oldPos, val.length() - (hasBrackets ? 2 : 0));
        curWord = trim(curWord);
        if (curWord.length() > 0) {
          // container with one entry
          convert(curWord, swap);
          c.push_back(swap);
        }
      }
    }

    // TODO: Verallgemeinerung der convert-method für bel. Vektoren und
    //       hinzunahme der Möglichkeit Matrizen anzugeben

    /// convert string to WorldVector
    template< typename T >
    inline void convert(const std::string valStr, WorldVector<T>& c) 
    {
      std::vector<T> temp_vec;
      getContainer(valStr, temp_vec);
       if (static_cast<int>(temp_vec.size()) != c.getSize())
        throw WrongVectorSize("wrong number of entries for WorldVector"); 
  
      for (size_t i = 0; i < temp_vec.size(); i++)
        c[i] = temp_vec[i];
    }

    /// convert string to std::list using begBrackets, endBrackets and delims
    template<typename T>
    inline void convert(const std::string valStr, std::list<T>& value)
    {
      getContainer(valStr, value);
    }


    /// convert string to std::vector using begBrackets, endBrackets and delims
    template<typename T>
    inline void convert(const std::string valStr, std::vector<T>& value)
    {
      getContainer(valStr, value);
    }


  } // end namespace detail

// _____________________________________________________________________________

  /** The entry in an initfile. This helper class was constructed to allow calls 
  *  like val = data.get(tag) for arbitrary types of val. At current stage, only
  *  double and bool is supported
  */
  struct InitEntry {
    ///the value as string
    std::string valStr;

    /// initialize with value as string
    InitEntry(std::string v = "")
    : valStr(v) 
    {}

    /// cast string to type T
    template<typename T>
    operator T() const 
    { 
      T t; 
      detail::convert(valStr, t); 
      return t;
    }
  };


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
  struct Initfile : public std::map<std::string, std::string> 
  {
    typedef std::map< std::string, std::string > super;

    static const int TAG_NOT_FOUND = 1;
    static const int TAG_NOT_FOUND_BREAK = 2;

    /// Exceptions
    struct TagNotFound : std::invalid_argument {
      TagNotFound(std::string m) 
      : std::invalid_argument(m) 
      {}
    };


    struct TagNotFoundBreak : std::invalid_argument {
    // print 'tag not found' and exit
      TagNotFoundBreak(std::string m)
      : std::invalid_argument(m) 
      {}
    };


    /** initialize init-file from file with filename in, read data and save it 
    *  to singleton-map
    *  @param in: filename string
    */
    static void init(std::string in);

    static void init(int print, std::string filename, const char *flags = NULL) 
    {
        WARNING("Parameters::init(int,std::string,const char*) is depreciated. "
          "Use Parameters::init(std::string) instead!\n");
        init(filename);
    }


    /** Static get routine for getting parameter-values from init-file 
    *  initialized in init()-method.
    *  Cast the value to the desired type using std::stringstream.
    *  @param tag: The tag to look for
    *  @param value: The result.
    *  @param debugInfo: msgInfo for current parameter. (0..no printing,
    *    1..print missing parameter info, 2..print parameter value) [optional]
    */
    template<typename T>
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

      std::string valStr;
      try {
      int error_code = singlett->checkedGet(tag, valStr);
      if (error_code == 0) {
	valStr = trim(valStr);
	detail::convert(valStr, value);
      } 
      } catch(mu::ParserError& e) {
	std::string parser_error = "Could not parse: " + tag;
	throw std::runtime_error(parser_error);
      }
      /*
else if(error_code == TAG_NOT_FOUND_BREAK)
	throw TagNotFoundBreak("required tag '" + tag + "' not found");
      else if (error_code == TAG_NOT_FOUND) {
	if (debugInfo == 2)
	  std::cout << "there is no tag '" + tag + "'" << std::endl;
      } else
	throw std::runtime_error("unknown error_code (" + boost::lexical_cast<std::string>(error_code) + ") returned for tag '" + tag + "'");
      */

      if (debugInfo == 2) {
	std::cout << "Parameter '" << tag << "'"
		  << " initialized with: " << value << std::endl;
      }
      singlett->msgInfo = debugInfo;
    }


    /** Static get routine for getting a parameter-map from init-file 
    *  initialized in init()-method.
    *  Idea: 
    *     line in Initfile: <tag><key>: <value>
    *     is extracted in am map [key]->[value].
    *  @param tag: The tag withch labels the parameter map 
    *  @param debugInfo: msgInfo for current parameter. (0..no printing,
    *    2..print parameter value) [optional]
    */
    static void getParameterMap(const std::string tag, 
				std::map<std::string,std::string> &pm,
				int debugInfo = -1)
    {
	initIntern();
	for (Initfile::iterator it = singlett->begin(); it != singlett->end(); it++){	
	    std::string longTag= (*it).first ;
	    std::string value=(*it).second;
	    if(longTag.length()>tag.length() && 
	       longTag.compare(0,tag.length(),tag)==0){
		if (debugInfo == 2){
		  std::cout <<"tag "<< tag <<"\n";
		  std::cout <<"Extract Parameter map from "<<(*it).first   
			    << " => " <<(*it).second  << std::endl;
		}

		std::string key=trim(longTag.substr(tag.length()));

		if (debugInfo == 2){
		  std::cout <<"Parameter map "<< key << " => " 
			    << value << std::endl;
		}
		pm[key]=value;
	    }
	}	    
    }
    
    /// return InitEntry object for tag tag
    static InitEntry get(const std::string tag) 
    {
      InitEntry result;
      
      std::string valStr;
      int error_code = singlett->checkedGet(tag, valStr);
      if (error_code == 0) {
	valStr = trim(valStr);
	result = InitEntry(valStr);
      } 

#if 0
else if(error_code == TAG_NOT_FOUND_BREAK)
	throw TagNotFoundBreak("get(): required tag '" + tag + "' not found");
      else if (error_code == TAG_NOT_FOUND)
	throw TagNotFound("get(): there is no tag '" + tag + "'"); 
	// exception must be thrown, because an empty object would be return otherwise
      else
	throw std::runtime_error("get(): unknown error_code returned for tag '" + tag + "'");
#endif

      return result;
    }


    /// update map tag->value_old to tag->value in singleton
    template<typename T>
    static void set(const std::string tag, T& value, int debugInfo=  -1) 
    {
      initIntern();
      if (debugInfo == -1)
        debugInfo = singlett->getMsgInfo();

      std::string swap = "";
      detail::convert(value, swap);
      (*singlett)[trim(tag)] = swap;
      // update msg parameters msgInfo, msgWait, paramInfo
      singlett->getInternalParameters();
      if (debugInfo == 2)
        std::cout << "Parameter '" << tag << "'"
                  << " set to: " << value << std::endl;
    }


    /// add map tag->value to data in singleton
    template< typename T >
    static void add(const std::string tag, T& value, int debugInfo = -1) 
    {
      set(tag, value, debugInfo);
    }


    /// rescheduling parameter
    static void readArgv(std::string parameters, int debugInfo = 2);


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
      singlett->clear();
    }


    /// save singlett-data to file with filename fn
    static void save(std::string fn)
    {
      initIntern();
      singlett->write(fn);
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


    /// list of processed files
    static std::set< std::string > fn_include_list;


    /// pointer to the singleton that contains the data
    static Initfile* singlett;


    /// return the value of the given tag or throws an exception if the tag 
    /// does not exist
    int checkedGet(const std::string& tag, std::string &valStr) const
    {
      super::const_iterator it = find(tag);
      if (it == end()) {
        if (breakOnMissingTag == 0 || msgInfo <= 2)
	  return TAG_NOT_FOUND;
        else
	  return TAG_NOT_FOUND_BREAK;
      }
      valStr = it->second;
      return 0;
    }

    /// replace variables by its value defined as parameter previousely
    /// variable definition is simple parameter definition
    /// variable evaluation by ${variablename} or $variablename
    /// the last version only for variablenames without whitespaces
    std::string variableReplacement(const std::string& input) const;

    std::string variableEvaluation(const std::string& input) const;
    
    /** Fill the initfile from an input stream.
    * @param in: the stream to fill the data from.
    * Current dataformat: tag:value
    * Comment char: percent '%'
    * Include files: #include "filename" or #include <filename>
    */
    void read(std::istream& in);

    /// Fill the initfile from a file with filename fn
    void read(std::string fn, bool force = false);

    /** Write data-map to initfile
    * @param out: the stream to fill the data in.
    */
    void write(std::ostream& out);

    /// Write data-map to initfile with filename fn
    void write(std::string fn);

    /// read parameters for msgInfo, msgWait, paramInfo
    void getInternalParameters();

    int msgInfo, msgWait, paramInfo, breakOnMissingTag;	
  };
  
  typedef Initfile Parameters;

} // end namespace AMDiS
#endif
