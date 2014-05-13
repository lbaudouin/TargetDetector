// Copyright 2014, Léo Baudouin
//
// This file is part of a free library: you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of
// the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// General Public License for more details. You should have
// received a copy of the GNU General Public License along with
// this library. If not, see <http://www.gnu.org/licenses/>.

#ifndef TARGETDETECTOR_H
#define TARGETDETECTOR_H

#include <vector>
#include <iostream>

#include <opencv2/opencv.hpp>
#include <opencv2/highgui/highgui.hpp>

#include "blob.h"

#include "target.h"


/** @class TargetDetector
 * This detector can find several type of target in an image.
**/
class TargetDetector
{
public:
  TargetDetector(int threshold = 125);
  ~TargetDetector();

  //Find grid
  bool findTargetsGrid(const cv::Mat &image, cv::Size size, std::vector<cv::Point2f> &centers, const Target &search, const std::vector<int> &values = std::vector<int>()) const;
  static bool findTargetsGrid(const cv::Mat &image, cv::Size size, std::vector<cv::Point2f> &centers, const Target &search, const int threshold, const std::vector<int> &values = std::vector<int>());

  //Find tragets
  std::vector<Target> track(const cv::Mat &image, const std::vector<Target> &searches) const;
  std::vector<Target> track(const cv::Mat &image, const Target &search) const;
  static std::vector<Target> track(const cv::Mat &image, const std::vector<Target> &searches, const int threshold);
  static std::vector<Target> track(const cv::Mat &image, const Target &search, const int threshold);

  //Draw targets and grid
  static void drawTarget(cv::Mat &image, const Target &target);
  static void drawTargets(cv::Mat &image, const std::vector<Target> &targets);
  static void drawTargetGrid(cv::Mat& image, const cv::Size &size, const std::vector< cv::Point2f >& centers, const bool &found);

  //Threshold
  void setThreshold(int threshold);
  int threshold() const;

  //Autothreshold
  bool autoThreshold(const cv::Mat &image, Target search, const int &step = 10);
  bool autoThreshold(const cv::Mat &image, const std::vector<Target> &search, const int &step = 10);
  bool autoThreshold(cv::VideoCapture &capture, Target search, const int &step = 10, const int &nbIterationMax = 100, const bool &display = false);
  bool autoThreshold(cv::VideoCapture &capture, const std::vector<Target> &search, const int &step = 10, const int &nbIterationMax = 100, const bool &display = false);

protected:
  //Find blobs
  std::vector<Blob> findBlobs(const cv::Mat &frameGray, const int &thresholdType = CV_THRESH_BINARY) const;
  Blob createBlob(int index, const std::vector<std::vector<cv::Point> > &contours, const std::vector<cv::Vec4i> &hierarchy, bool &ok) const;
  Blob createBlob(const std::vector<cv::Point> &external, bool &ok) const;

  //Check if blob is target
  std::vector<Target> checkBlobs(const std::vector<Blob> &blobs, const cv::Mat &frameGray, const Target &search) const;
  std::vector<Target> checkOneBlob(const Blob &blob, const cv::Mat &frameGray, const Target &search) const;
  std::vector<Target> checkThreeBlobs(const Blob &blob, const cv::Mat &frameGray, const Target &search) const;
  std::vector<Target> checkTwoRings(const Blob &blob, const cv::Mat &frameGray, const Target &search) const;

  bool isEllipse(const Blob &blob, const double &threshold = -1) const;
  bool isEllipseWithoutHole(const Blob &blob, const double &threshold = -1) const;
  bool isConcentricEllispe(const Blob &reference, const Blob &ellipse, const double &maxCenterDistance, const double &minFactor, const double &maxFactor) const;

  //Code extractor
  bool extractCode(const Blob &blob, const cv::Mat &frameGray, Target &target) const;
  bool extractCodeFromImage(const cv::Mat &frameGray, cv::Point2d center, double a, double b, double orientation, double angle, std::vector<int> &codeValues, std::vector<cv::Point2f> &codePoints, double shift = 0.0) const;

  bool alignPointsOnFirstBit(const std::vector<int> &header, std::vector<int> &code, std::vector<cv::Point2f> &points, const double &minScore, int *shift = NULL) const;
  
  void filterCode(std::vector<int> &code) const;
  std::string code2string(const std::vector<int> &code, const int &nbPointPerBit) const;

  int findHeader(const std::string &str, const std::string &headerString) const;

  bool readCode(Target &target, const std::string &codeString) const;

  //Display
  static void drawLineBetweenPoints(cv::Mat &image, const cv::Point2f &p1, const cv::Point2f &p2);

private:
  int m_binaryThreshold;
  double m_ellipseThreshold;
};

#endif // TARGETDETECTOR_H
