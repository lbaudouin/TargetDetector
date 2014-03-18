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
