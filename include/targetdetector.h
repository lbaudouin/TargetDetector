#ifndef TARGETDETECTOR_H
#define TARGETDETECTOR_H

#define TARGET_DEBUG false

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
  
  std::vector<Target> track(const cv::Mat &image, std::vector<Target> searches);
  std::vector<Target> track(const cv::Mat &image, Target search);
  
  void drawTarget(cv::Mat &image, const Target &target);
  void drawTargets(cv::Mat &image, const std::vector<Target> &targets);
  
  void setThreshold(int threshold);
  bool autoThreshold(const cv::Mat &image, Target search, int step = 10);
  bool autoThreshold(cv::VideoCapture &capture, Target search, int step = 10, int nbIterationMax = 100);
    
protected:
  //Find blobs
  std::vector<Blob> findBlobs(const cv::Mat &frameGray, int thresholdType = CV_THRESH_BINARY);
  Blob createBlob(int index, const std::vector<std::vector<cv::Point> > &contours, const std::vector<cv::Vec4i> &hierarchy, bool &ok);
  Blob createBlob(std::vector<cv::Point> external, bool &ok);
  
  //Check if blob is target
  std::vector<Target> checkBlobs(const std::vector<Blob> &blobs, const cv::Mat &frameGray, Target search);
  std::vector<Target> checkOneBlob(const Blob &blob, const cv::Mat &frameGray, Target search);
  std::vector<Target> checkThreeBlobs(const Blob &blob, const cv::Mat &frameGray, Target search);
  std::vector<Target> checkTwoRings(const Blob &blob, const cv::Mat &frameGray, Target search);
  
  bool isEllipse(const Blob &blob, double threshold = -1);
  bool isEllipseWithoutHole(const Blob &blob, double threshold = -1);
  bool isConcentricEllispe(const Blob &reference,const Blob &ellipse, double maxCenterDistance, double minFactor, double maxFactor);
  
  //Code extractor
  bool extractCode(const Blob &blob, const cv::Mat &frameGray, Target &target);
  bool extractCodeFromImage(const cv::Mat &frameGray, cv::Point2d center, double a, double b, double orientation, double angle, std::vector<int> &codeValues, std::vector<cv::Point2f> &codePoints, double shift = 0.0);
  
  int alignPointsOnBit(std::vector<int> &code, std::vector<cv::Point2f> &points);
  double alignPointsOnBit(std::vector<int> &code, std::vector<cv::Point2f> &points, double angle);
  
  void filterCode(std::vector<int> &code);
  std::string code2string(std::vector<int> &code, int nbPointPerBit);
  
  int findHeader(std::string str, std::string headerString);
  
  bool readCode(Target &target, std::string codeString);
  
  //Display
  void drawLineBetweenPoints(cv::Mat &image, cv::Point2f p1, cv::Point2f p2);
    
private:
  int m_binaryThreshold;
  double m_ellipseThreshold;
};

#endif // TARGETDETECTOR_H
