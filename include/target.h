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

#ifndef TARGET_H
#define TARGET_H

#include <opencv2/opencv.hpp>

#include "codeconversion.h"

/** @class Target
 * Create an object to search targets in image defining the type and header
*/
class Target
{
public:
  /** @enum Type
   *  Define the type of target
   */
  enum Type
  {
    OneBlob = 0, /*!< One blob in the center */ 
    ThreeBlobs, /*!< Three center blobs */
    TwoRings /*!< Two rings of data */
  };
  
  /** Only position target objet (don't read code)
   * @arg target_type is the target @ref Type
   */
  Target(Type target_type) : m_type(target_type), m_readCode(false), m_orientation(-1)
  {
    initType();
    m_headerValue = -1;
    m_headerBits = -1;
  }
  
  /** Full information target objet
   * @arg target_type is the target @ref Type
   */
  Target(Type target_type, int header_bits, int message_bits, bool use_parity_bit, int header_value, bool use_gray_code = false) : 
	m_type(target_type), m_headerBits(header_bits), m_messageBits(message_bits), m_useParityBit(use_parity_bit), m_headerValue(header_value),
	m_useGrayCode(use_gray_code), m_nbPointPerBit(8), m_readCode(true), m_firstHeaderIndex(0), m_inverseParityBit(false), m_orientation(-1)
  {
    m_nbBits = m_headerBits + m_messageBits + (use_parity_bit?1:0);
        
    initType();
  }
   
  /** Get target type
   * @see Type
   * @return The Type of target
   **/
  Type type() const { return m_type; }
  
  /** Return true if gray code is use in header and message **/
  bool useGrayCode() const { return m_useGrayCode; }
  
  /** Return true is the message contains a parity bit **/
  bool useParityBit() const { return m_useParityBit; }
  
  /** Set if the parity bit need to be inversed **/
  void setInverseParityBit(bool inverse) { m_inverseParityBit = inverse; }
  
  /** Return true if the parity bit need to be inversed **/
  bool inverseParityBit() const { return m_inverseParityBit; }
  
  /** Return true if the message need to be read **/
  bool readCode() const { return m_readCode; }
  
  /** Set number of point per bit
   * @param nbPointPerBit is the number of point per bit
   **/
  void setNumberPointPerBit(int nbPointPerBit) { m_nbPointPerBit = nbPointPerBit; }
  
  /** Return the number of point per bit **/
  int numberPointPerBit() const { return m_nbPointPerBit; }
  
  /** Return the number of bits (header + message **/
  int nbBits() const { return m_nbBits; }
  
  /** Return the header value **/
  int header() const { return m_headerValue; }
  
  /** Return the message value **/
  int message() const { return m_messageValue; }
  
  /** Return the number of bits in the header **/
  int headerBits() const { return m_headerBits; }
  
  /** Return the number of bits in the message **/
  int messageBits() const { return m_messageBits; }
  
  /** Return the string formated header **/
  std::string headerString() const
  {
    if(m_useGrayCode)
      return dec2gray(m_headerValue,m_headerBits);
    else
      return dec2bin(m_headerValue,m_headerBits);
  }
    
  /** Return the string formated message **/
  std::string messageString() const
  {
    if(m_useGrayCode)
      return dec2gray(m_messageValue,m_messageBits);
    else
      return dec2bin(m_messageValue,m_messageBits);
  }
   
  /** Set message after reading
   * @param value is new message value
   * @param nbBits is number of bits for the message
   **/
  void setMessage(int value, int nbBits)
  {
    if(nbBits!=m_messageBits)
      throw "Wrong number of bits";
    if( value > (pow(nbBits,2)-1) )
      throw "Number of bits required for message to high";
    m_messageValue = value;
  }
  
  /** Set message after reading
   * @param value is new message value
   * @param useGrayCode set if Gray code is used
   **/
  void setMessage(std::string value, bool useGrayCode = false)
  {
    int nbBits = value.length();
    if(nbBits!=m_messageBits)
      throw "Wrong number of bits";
    m_useGrayCode = useGrayCode;
    if(m_useGrayCode)
      m_messageValue = gray2dec(value);
    else
      m_messageValue = bin2dec(value);
  }
  
  /** Set center of target
   * @param center is the input center of target
   **/
  void setCenter(cv::Point2f center) { m_center = center; }
  
  /** Return the center of target **/
  cv::Point2f center() const { return m_center; }
  
  /** cv::Point2f type casting operator.**/
  operator cv::Point2f () const { return center() ; } ;
  
  /** Return the orientation of target
   * angle is in [0:2/pi]
   * Return -1 if invalid
   **/
  double orientation() const { return m_orientation; }
  
  /** Set ellipse parameters
   * @param a is major axis
   * @param b is minor axis
   * @param orientation is ellipse orientation
   **/
  void setEllipseParameters(double a, double b, double orientation) { m_ellipseMajorAxis = a; m_ellipseMinorAxis = b; m_ellipseOrientation = orientation; }
  
  /** Get ellipse parameters
   * @param a is major axis
   * @param b is minor axis
   * @param orientation is ellipse orientation
   **/
  void ellipseParameters(double &a, double &b, double &orientation) const { a = m_ellipseMajorAxis; b = m_ellipseMinorAxis; orientation = m_ellipseOrientation; }
  
  /** Return first radius factor to read the header/message **/
  double radiusFactor() const { return m_radiusFactor; }
  
  /** Return second radius factor to read message **/
  double radiusFactor2() const { return m_radiusFactor2; }
  
  /** Set points and value used to extract header/message
   * @param points points position
   * @param pointsValue value read at point position
   **/
  void setPoints(std::vector<cv::Point2f> points, std::vector<int> pointsValue)
  {
    assert(points.size()==pointsValue.size());
    m_points = points;
    m_pointsValue = pointsValue;
  }
  
  /** Return number of points **/
  int nbPoints() const { return m_points.size(); }
  
  /** Get point position
   * @param index is the index of point in the vector
   **/
  cv::Point2f point(int index) const { return m_points.at(index); }
  
  /** Get point value
   * @param index is the index of point in the vector
   **/
  int pointValue(int index) const { return m_pointsValue.at(index); }
    
  /** Set index of the first header points
   * @param index is the index of the first header point
   **/
  int setFirstHeaderIndex(int index) {
    assert(index<m_points.size());
    m_firstHeaderIndex = index;
    cv::Point2f delta = m_points[index] - m_center;
    m_orientation = std::atan2( delta.y, delta.x );
  }
  
  /** Get the index of the first header point **/
  int firstHearderIndex() const { return m_firstHeaderIndex; }
  
  /** @warning Not read-only list of contours **/
  std::vector<cv::Point> contours;  
    
protected:
  /** Initialize the type of target **/
  void initType()
  {
    switch(m_type){
      case OneBlob:
	m_radiusFactor = 3.25;
	m_radiusFactor2 = -1;
	m_nbPointPerBit = 16;
	break;
      case ThreeBlobs:
	m_radiusFactor = 0.97*3.0;
	m_radiusFactor2 = -1;
	m_nbPointPerBit = 16;
	break;
      case TwoRings:
	m_radiusFactor = 0.97*5.0/6.0;
	m_radiusFactor2 = 0.97*7.0/6.0;
	m_nbPointPerBit = 16;
	break;
    }
  }
  
private:
  
  Type m_type; /*!< The target type */
   
  int m_nbBits; /*!< Number of bits in the code */
  int m_nbPointPerBit; /*!< Number of image point in a bit */
  bool m_useGrayCode; /*!< Use Gray code */
  bool m_useParityBit; /*!< Use parity bit */
  bool m_inverseParityBit; /*!< Inverse parity bit */
  
  int m_headerValue; /*!< Header value */
  int m_headerBits; /*!< Header number of bits */
  int m_messageValue; /*!< Message value */
  int m_messageBits; /*!< Message number of bits */
  
  int m_firstHeaderIndex; /*!< First header index */
  
  cv::Point2f m_center; /*!< Center of the target */
  double m_orientation; /*!< Orientation of the target */
  double m_ellipseMajorAxis; /*!< Major ellipse axis */
  double m_ellipseMinorAxis; /*!< Minor ellipse axis */
  double m_ellipseOrientation; /*!< Ellipse orientation */
  double m_radiusFactor; /*!< Factor to apply on the ellipse to be on the code ellipse */
  double m_radiusFactor2; /*!< Second factor to apply on the ellipse to be on the code ellipse */
  
  bool m_readCode; /*!< Read the target code */
  
  std::vector<cv::Point2f> m_points; /*!< Points used to extract code */
  std::vector<int> m_pointsValue; /*!< Point values at extracted position */
};

#endif // TARGET_H
