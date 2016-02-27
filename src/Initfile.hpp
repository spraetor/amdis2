#pragma once

#include <string>
#include <iostream>

#include <boost/lexical_cast.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/tokenizer.hpp>

#include <boost/property_tree/ptree.hpp>
// #include <boost/type_traits.hpp>

#include "AMDiS_base.hpp"
#include "Log.hpp"
#include "Math.hpp"
#include "MatrixVector.hpp"
#include "traits/basic.hpp"
#include "traits/concepts_base.hpp"
#include "traits/meta_basic.hpp"
#include "traits/scalar_types.hpp"

namespace AMDiS
{
  namespace detail
  {
    double mu_parser_eval(std::string const& valStr);

    template <class T> inline T numeric_cast(double value)
    {
      return boost::numeric_cast< T >(value);
    }

    template <> inline bool numeric_cast<bool>(double value)
    {
      return value != 0.0;
    }

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
    struct Convert<T, Requires_t<concepts::Arithmetic<T>> >
    {
      static void eval(std::string valStr, T& value)
      {
        try {
            value = numeric_cast< T >(mu_parser_eval(valStr));
        } catch(...) {
            ERROR("Could not parse '%s' to '%s'\n", valStr.c_str(), typeid(T).name());
        }
      }
    };

    template <class T>
    struct Convert<T, Requires_t<std::is_enum<T>> >
    {
      static void eval(std::string valStr, T& value)
      {
        EnumParser<T>()(valStr, value);
      }
    };

    // convert string to vector
//     template <class T>
//     struct Convert<T, Requires_t<concepts::Vector<T>> >
//     {
//       static void eval(std::string valStr, T& values)
//       {
//         using value_type = Value_t<T>;
//         using Tokenizer = boost::tokenizer<boost::char_separator<char>>;
// 
//         boost::char_separator<char> sep(",; ");
//         Tokenizer tokens(valStr, sep);
//         int i = 0;
//         for (auto token : tokens)
//         {
//           value_type v;
//           Convert<value_type>::eval(token, v);
//           values[i] = v;
//         }
//       }
//     };

    template <class M, class S>
    struct Convert<VectorBase<M,S>>
    {
      using T = VectorBase<M,S>;
      static void eval(std::string valStr, T& values)
      {
        using value_type = Value_t<T>;
        using Tokenizer = boost::tokenizer<boost::char_separator<char>>;

        boost::char_separator<char> sep(",; ");
        Tokenizer tokens(valStr, sep);
        int i = 0;
        for (auto token : tokens)
        {
          value_type v;
          Convert<value_type>::eval(token, v);
          values[i] = v;
        }
      }
    };

    // convert string to vector
    template <class T, class Alloc>
    struct Convert<std::vector<T, Alloc>>
    {
      static void eval(std::string valStr, std::vector<T>& values)
      {
        using value_type = T;
        using Tokenizer = boost::tokenizer<boost::char_separator<char>>;

        boost::char_separator<char> sep(",; ");
        Tokenizer tokens(valStr, sep);
        for (auto token : tokens)
        {
          value_type v;
          Convert<value_type>::eval(token, v);
          values.push_back(v);
        }
      }
    };

  } // end namespace detail


  /// output-stream for std::list
  template <class T, class Alloc>
  std::ostream& operator<<(std::ostream& out, std::list<T,Alloc> const& l)
  {
    auto it = l.begin();
    out << "["; if (l.size() > 0) out << *it;
    for (; it != l.end(); ++it)
      out << ", " << *it;
    out << "]";
    return out;
  }


  /// output-stream for std::vector
  template <class T, class Alloc>
  std::ostream& operator<<(std::ostream& out, std::vector<T,Alloc> const& l)
  {
    auto it = l.begin();
    out << "["; if (l.size() > 0) out << *it;
    for (; it != l.end(); ++it)
      out << ", " << *it;
    out << "]";
    return out;
  }

  inline void replaceAll(std::string& str, std::string const& from, std::string const& to)
  {
    if (from.empty())
      return;
    size_t start_pos = 0;
    while ((start_pos = str.find(from, start_pos)) != std::string::npos)
    {
      str.replace(start_pos, from.length(), to);
      start_pos += to.length();
    }
  }

  // _____________________________________________________________________________

  /** Basis data container as a map of tag on a value as strings. The container
  *  throws an exception, if the tag was not found.
  */
  struct Initfile
  {
    using Self = Initfile;

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
    static void get(std::string tag, T& value, int debugInfo = -1)
    {
      initIntern();
      if (debugInfo == -1)
        debugInfo = singlett->getMsgInfo();
      else
      {
        int swap(debugInfo);
        debugInfo = singlett->getMsgInfo();
        singlett->msgInfo=swap;
      }

      using path = boost::property_tree::ptree::path_type;
      replaceAll(tag, "->", ">");
      auto tagPath = path(tag, '>');

      // TODO: use boost::optional instead
      // TODO: use convert method from above
      std::string valueStr = "-";
      valueStr = singlett->pt.get(tagPath, valueStr);

      if (valueStr != "-")
        detail::Convert<T>::eval(valueStr, value);

      if (debugInfo == 2)
      {
        std::cout << "Parameter '" << tag << "'"
                  << " initialized with: " << value << std::endl;
      }
      singlett->msgInfo = debugInfo;
    }


    template <class T, class S>
    static T get(std::string tag, S const& default_value)
    {
      T value = default_value;
      Self::get(tag, value);
      return value;
    }


    template <class T>
    static T get(std::string tag)
    {
      T value; nullify(value);
      Self::get(tag, value);
      return value;
    }


    /// update map tag->value_old to tag->value in singleton
    template <class T>
    static void set(std::string tag, T& value, int debugInfo=  -1)
    {
      initIntern();
      if (debugInfo == -1)
        debugInfo = singlett->getMsgInfo();

      using path = boost::property_tree::ptree::path_type;
      replaceAll(tag, "->", ">");
      auto tagPath = path(tag, '>');
      singlett->pt.put(tagPath, value);

      // update msg parameters msgInfo, msgWait, paramInfo
      singlett->getInternalParameters();
      if (debugInfo == 2)
        std::cout << "Parameter '" << tag << "'"
                  << " set to: " << value << std::endl;
    }


    /// add map tag->value to data in singleton
    template <class T>
    static void add(std::string tag, T& value, int debugInfo = -1)
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
    static Initfile* getSingleton()
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
    static void save(std::string /*fn*/)
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

    int msgInfo;
    int msgWait;
    int paramInfo;
    int breakOnMissingTag;

    /// boost:property_tree to read/store parameter values
    boost::property_tree::ptree pt;
  };

  using Parameters = Initfile;

#ifndef AMDIS_NO_EXTERN_INITFILE
  extern template void Initfile::get(std::string, int&, int);
  extern template void Initfile::get(std::string, double&, int);
  extern template void Initfile::get(std::string, std::string&, int);
#endif

} // end namespace AMDiS
