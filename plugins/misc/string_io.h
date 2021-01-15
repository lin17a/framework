// -*- C++ -*-
// author: afiq anuar
// short: a listing of free functions that contains common string and i/o operations

#ifndef FWK_STRING_IO_H
#define FWK_STRING_IO_H

#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>

/// string replacement - replaces the first occurence only
inline bool replace(std::string &str, const std::string &from, const std::string &to) {
  std::string::size_type start_pos = str.find(from);
  if (start_pos == std::string::npos)
    return false;
  str.replace(start_pos, from.length(), to);
  return true;
}



/// count occurences of substrings
/// credit https://stackoverflow.com/questions/22489073/counting-the-number-of-occurrences-of-a-string-within-a-string
inline int count_substring(const std::string &str, const std::string &sub) {
  int nSub = 0;
  std::string::size_type iSub = 0;
  while ((iSub = str.find(sub, iSub)) != std::string::npos) {
    ++nSub;
    iSub += sub.length();
  }

  return nSub;
}



/// number to string; to_string tend to give more precision than needed
template <typename Number> 
std::string to_str(Number num, const int prec = -1, const bool fixed = false) 
{
  std::ostringstream out_str;
  if (fixed)
    out_str << std::fixed;
  if (prec > 0)
    out_str << std::setprecision(prec);

  out_str << num; 
  return out_str.str(); 
}



/// tokenize a string i.e. split a string into multiple strings by a separator
/// returns a vector of string; vector contains original string if separator isn't present within it
std::vector<std::string> tokenize(const std::string &str, const std::string &sep = ",")
{
  std::vector<std::string> v_str = {};
  std::string::size_type isep = str.find(sep), iini = 0;
  while (isep != std::string::npos) {
    v_str.emplace_back( str.substr(iini, isep - iini) );
    iini = isep + sep.length();
    isep = str.find(sep, iini);
  }
  v_str.emplace_back( str.substr(iini, isep - iini) );

  return v_str;
}



/// strips a given substring from a string
/// i.e. runs replace() until the substring can no longer be found
std::string& strip(std::string &str, const std::string &sub = " ")
{
  if (sub.length() < 1)
    return str;

  while( replace(str, sub, "") )
    ;

  return str;
}


/// returns a logging function to a desired stream
/// works, but makes the compilation time MUCH longer when used with Group::iterate()...
/// and using it with 2 calls, 8 (same) arguments each makes the executable increase by 8MB!!!
/// however with 2 calls and 2 same arguments each it costs only 4kB, and short compilation time
/// likely just an issue with how it interacts with std::visit
auto logger(std::ostream &out) 
{
  return [&out] (const auto &arg, const auto &...args) {
    out << arg;
    ((out << " " << args), ...);
    out << "\n";
  };
}

#endif
