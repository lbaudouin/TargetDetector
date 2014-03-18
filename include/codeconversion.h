// Copyright 2014, Léo Baudouin
//
// This file is part of a free library: you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of
// the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details. You should have
// received a copy of the GNU General Public License along with
// this library. If not, see <http://www.gnu.org/licenses/>.

#ifndef CODECONVERSION_H
#define CODECONVERSION_H

#include <string>
#include <sstream>
#include <algorithm>

namespace {
  /** Convert binary value to decimal
   * @param bin is input binary string value
   **/
  inline int bin2dec(std::string bin)
  {
    int val = 0;
    int pow2 = 1;  
    for(int j=bin.size()-1;j>=0;j--){
      if(bin[j]=='1'){
	val += pow2;
      }
      pow2 *= 2;
    }
    return val;
  }

  /** Convert decimal value to binary string
   * @param val is input decimal value
   * @param nbBits is the number of bit to code the value
   **/
  inline std::string dec2bin(int val, int nbBits)
  {
    if(val>pow(2,nbBits)-1)
      return "";
    std::stringstream ss;
    while(val > 0) {
      ss << val % 2 == 0 ? '0' : '1';
      val /= 2;
    }
    std::string str = ss.str();
    std::reverse(str.begin(), str.end());
    return str;
  }

  /** Convert binary value to Gray code
   * @param bin is input binary string value
   **/
  inline std::string bin2gray(std::string bin)
  {
    if(bin.empty())
      return "0";
    int cnt = 0;
    std::stringstream gray;
    gray << bin[0];
    while(cnt<bin.length()-1){
      gray << ((bin[cnt] != bin[cnt+1])?"1":"0");
      cnt++;
    }
    return gray.str();
  }

  /** Convert Gray code to binary
   * @param gray is input Gray code string value
   **/
  inline std::string gray2bin(std::string gray)
  {
    if(gray.empty())
      return "0";
    int cnt = 0;
    std::stringstream bin;
    bin << gray[0];
    while(cnt<gray.length()-1){
      bin << ((bin.str()[cnt] != gray[cnt+1])?"1":"0");
      cnt++;
    }
    return bin.str();
  }

  /** Convert decimal to Gray code binary
   * @param val is input decimal value
   * @param nbBits is the number of bit to code the value
   **/
  inline std::string dec2gray(int val, int nbBits)
  {
    return bin2gray(dec2bin(val,nbBits));
  }

  /** Convert Gray code to decimal
   * @param gray is input Gray code string value
   **/
  inline int gray2dec(std::string gray)
  {
    return bin2dec(gray2bin(gray));
  }

  /** Get the parity of message
   * @param bin is input binary message
   * @return 0 if number of ones is even, 1 if odd
   **/
  inline int parity(std::string bin)
  {
    int nbOnes = 0;
    for(int i=0;i<bin.size();i++)
      if(bin[i]=='1')
	nbOnes++;
      return nbOnes%2;
  }
}

#endif // CODECONVERSION_H