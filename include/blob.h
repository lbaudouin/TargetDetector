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

#ifndef BLOB_H
#define BLOB_H

#include <opencv2/opencv.hpp>

/** @struct Blob
 * Structure containing information about the blob and its children
 */
struct Blob
{
  int index; /*!< Index of blob during extraction */
  cv::Point2f center; /*!< Center of blob */
  double majorAxis; /*!< Ellipse major axis (a) */
  double minorAxis; /*!< Ellipse minor axis (b) */
  double orientation; /*!< Ellipse orientation */
  double area; /*!< External area of blob */
  double realArea; /*!< Real area (remove intrior blobs) */
  std::vector<cv::Point> contours; /*!< Contours of blob */
  std::vector<Blob> children; /*!< List of children */
  bool isValid; /*!< True if blob has more than 5 points (ellipse approx) */
};

#endif // BLOB_H
