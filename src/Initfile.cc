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
#include "Initfile.h"

#include <string>
#include <stdexcept>
#include <iostream>
#include <sstream>
#include<boost/tokenizer.hpp>

using namespace std;

namespace AMDiS {

  /// the small parser for the initfile. see description of 
  /// read(Initfile&,istream&)
  struct Parser {
    Parser(const string& line) 
    {
      size_t pos = line.find(':');
      if (pos == string::npos) {
	throw runtime_error("cannot find the delimiter ':' in line "
			    "'" + line + "'");
      }
      name = line.substr(0, pos);
      value = line.substr(pos + 1, line.length() - (pos + 1));

      // remove everything after the %
      pos = value.find('%');
      if (pos != string::npos)
        value = value.substr(0, pos);      
    }
    string name;
    string value;
  };


  Initfile* Initfile::singlett = NULL;
  std::set<std::string> Initfile::fn_include_list;


  /// initialize singleton object an global parameters
  void Initfile::init(std::string in)
  {
    initIntern();
    fn_include_list.clear();
    singlett->read(in);	
    singlett->getInternalParameters();
  }


  /// Fill an initfile from a file with filename fn
  void Initfile::read(std::string fn, bool force)
  {
    // read file if its not parsed already
    if (fn_include_list.find(fn) == fn_include_list.end() || force) {
      std::ifstream inputFile;
      inputFile.open(fn.c_str(), std::ios::in);
      if (!inputFile.is_open())
	throw runtime_error("init-file '" + fn + "' cannot be opened for reading");

      fn_include_list.insert(fn);
      read(inputFile);
    }
  }


  /// Fill an initfile from an input stream
  void Initfile::read(istream& in) 
  {
    const unsigned line_length = 512;
    char swap[line_length];
    in.getline(swap, line_length);
    while (in.good() || in.gcount() > 0) {
      std::string whitespaces = " \t\r\f\n";
      std::string delimiter = "\r\n";
      std::string sw(swap);
      size_t pos0 = sw.find_first_not_of(whitespaces);

      if (pos0 != std::string::npos
          && sw[pos0] != '%' 
          && sw[pos0] != '#'
          && sw[pos0] != 0) {
        // parse line and extract map: tag->value
        Parser parser(sw);
        
        // add parameter to map after variable replacement
        std::string paramName = variableReplacement(trim(parser.name));
        std::string paramValue = variableReplacement(trim(parser.value));
        paramValue = variableEvaluation(paramValue);
      
        operator[](paramName) = paramValue;
        int info = 0;
        get("parameter information", info, 0);
        if (info >= 2)
          std::cout << "read [" << paramName << " => " << paramValue << "]\n";
      } else if (pos0 != std::string::npos &&
		 sw[pos0] == '#'  &&
		 static_cast<size_t>(sw.find("#include")) == pos0) {
        // include file by '#include "filename"' or '#include <filename>'
        bool forceRead = false;
        int posAdd = 1;
        if (sw[pos0 + std::string("#include").size()] == '!') {
          forceRead = true;
          posAdd++;
        }
        size_t pos = sw.find_first_not_of(whitespaces, 
          pos0 + std::string("#include").size() + posAdd);
        size_t epos = 0;
        std::string fn =  "";
        std::stringstream errorMsg;
        switch (char c= swap[pos++]) {
          case '<':
            c= '>';
          case '\"':
            delimiter += c;
            epos = sw.find_first_of(delimiter, pos);
            fn = sw.substr(pos, epos - pos);

            if (sw[epos]!=c) {
              errorMsg << "filename in #include not terminated by " << c;
	      throw std::runtime_error(errorMsg.str());
            }
            break;
	default:
	  throw std::runtime_error("no filename given for #include");
        }
        
        
        read(fn, forceRead);
      }

      in.getline(swap, line_length);
    }
  }


  std::string Initfile::variableReplacement(const std::string& input) const
  {
    std::string whitespaces = " \t\r\f";
    std::string allowedChars = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_";
    std::string inputSwap = input; 
    size_t posVar = inputSwap.find_first_of('$');
    while (posVar != string::npos) {
      size_t posVarBegin, posVarEnd;
      if (inputSwap[posVar+1] == '{') {						// ${var_name}
        posVarEnd = inputSwap.find_first_of('}',posVar + 2);
        posVarBegin = posVar + 1;
      } else if (inputSwap[posVar+1] != '(' && inputSwap[posVar+1] != '[') {
        posVarEnd = inputSwap.find_first_not_of(allowedChars, posVar + 1);
        posVarBegin = posVar;
      } else {
	posVar = inputSwap.find_first_of('$',posVar+1);
	continue;
      }
      
      std::string varName = inputSwap.substr(posVarBegin + 1 , posVarEnd - posVarBegin - 1);


      // if varname is found in parameter list then replace variable by value
      // otherwise throw tagNotFound exception
      std::string varParam;
      int error_code = checkedGet(varName, varParam);
      if (error_code > 2)
	throw TagNotFoundBreak("required tag '" + varName + "' for variable substitution not found");
      
      std::string replaceName = inputSwap.substr(posVar , posVarEnd - posVar + (posVarBegin - posVar));
      inputSwap.replace(inputSwap.find(replaceName), replaceName.length(), varParam);

      posVar = inputSwap.find_first_of('$',posVarBegin);
    }

    return inputSwap;
  }
  

  std::string Initfile::variableEvaluation(const std::string& input) const
  {
    std::string whitespaces = " \t\r\f";
    std::string inputSwap = input;
    size_t posVar = inputSwap.find("$(");
    while (posVar != string::npos) {
      size_t posVarBegin, posVarEnd;
      posVarEnd = inputSwap.find_first_of(')',posVar + 2);
      posVarBegin = posVar + 1;
      std::string varName = inputSwap.substr(posVarBegin + 1 , posVarEnd - posVarBegin - 1);

      double value = 0.0;
      detail::convert(varName, value); // string -> double (using muparser)
      detail::convert(value, varName); // double -> string

      std::string replaceName = inputSwap.substr(posVar , posVarEnd - posVar + (posVarBegin - posVar));
      inputSwap.replace(inputSwap.find(replaceName), replaceName.length(), varName);

      posVar = inputSwap.find("$(",posVarBegin);
    }

    return inputSwap;
  }

  void Initfile::readArgv(std::string parameters, int debugInfo)
  {
    initIntern();
    
    char seperator = ';';
    typedef boost::escaped_list_separator<char> TokenizerFunc;
    typedef boost::tokenizer<TokenizerFunc> Tokenizer;
    TokenizerFunc tokenizerFunc('\\', seperator, '\"');
    Tokenizer tok(parameters, tokenizerFunc);

    // split parameterstring by seperator ";"
    std::vector<std::string> val;
    for (Tokenizer::iterator cell = tok.begin(); cell != tok.end(); ++cell) {
      val.push_back(trim(*cell));
    }
    
    // split each parameter by ":"
    for (size_t i = 0; i < val.size(); i++) {
      int found = val[i].find_first_of(':');
      if (found != static_cast<int>(std::string::npos)) {
	std::string value = trim(val[i].substr(found+1));
	std::string key = trim(val[i].substr(0, found));
	set(key, value, debugInfo);
      }
    }
  }


  /// read standard values for output and information of parameter-values
  void Initfile::getInternalParameters() 
  {
    int val = 0;
    get("level of information", val, 0);
    msgInfo = val;

    val = 1;
    get("WAIT", val, 0);
    msgWait = val;

    val = 1;
    get("parameter information", val, 0);
    paramInfo = val;

    val = 0;
    get("break on missing tag", val, 0);
    breakOnMissingTag = val;

    if (msgInfo == 0)
      paramInfo = 0;
  }


  /// print all parameters to std::cout
  void Initfile::printParameters() 
  {
    initIntern();
    for (Initfile::iterator it = singlett->begin(); it != singlett->end(); it++)
      std::cout << (*it).first << " => " << (*it).second << std::endl;
  }


  /// Write data-map to initfile
  void Initfile::write(ostream& out)
  {
    for (Initfile::iterator it = begin() ; it!=end(); it++)
      out << (*it).first << ": " << (*it).second << std::endl;	
  }


  /// Write data-map to initfile
  void Initfile::write(std::string fn) 
  {
    std::ofstream outFile;
    outFile.open(fn.c_str(), std::ios::out);
//     if (!outFile.is_open())
//       throw runtime_error("init-file cannot be opened for writing");

    write(outFile);
  }
} // end namespace AMDiS
