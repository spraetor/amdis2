#pragma once

// std c++ headers
#include <algorithm>
#include <cctype>
#include <locale>
#include <string>

namespace AMDiS
{
  // convert all characters in a string to upper case
  inline std::string to_upper(std::string input)
  {
    for (auto & c : input) c = toupper(c);
    return input;
  }
  
  // convert all characters in a string to upper case
  inline std::string to_lower(std::string input)
  {
    for (auto & c : input) c = tolower(c);
    return input;
  }
  
  // trim a string from the left
  inline std::string& ltrim(std::string& str)
  {
    auto it =  std::find_if(str.begin(), str.end(), [](char ch)
    { 
      return !std::isspace<char>(ch, std::locale::classic()); 
    });
    str.erase(str.begin() , it);
    return str;   
  }

  // trim a string from the right
  inline std::string& rtrim(std::string& str)
  {
    auto it =  std::find_if(str.rbegin(), str.rend(), [](char ch)
    { 
      return !std::isspace<char>(ch, std::locale::classic()); 
    });
    str.erase(it.base(), str.end());
    return str;   
  }
  
  // trim a string from both sides
  inline std::string& trim(std::string& str)
  {
    return ltrim(rtrim(str));
  }
  
  // trim a (copy of the) string from both sides
  inline std::string trim_copy(std::string const& str)
  {
    auto s = str;
    return trim(s);
  }

} // end namspace AMDiS
