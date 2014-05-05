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

#include "../include/targetdetector.h"

///////////////////////////////////// DEBUG FUNCTIONS ///////////////////////////////////////////

#ifdef DEBUG_MODE
#include "ostream"
std::vector<cv::Point2f> pts_images_check;
std::vector<int> pts_images_val;

std::ostream& operator<<(std::ostream& os, const std::vector<int>& v)
{
  for(unsigned int i=0; i<v.size(); i++){
    if(i%16==0 && i!=0)
      os << "-";
    os << v[i];
  }
  return os;
}

bool isEllipse(const Blob &blob, const double &threshold);
void drawBlobRecursively(cv::Mat &image, const Blob &blob, int depth = 0);
void drawEllipseRecursively(cv::Mat &image, const Blob &blob, int depth = 0);
#endif

//////////////////////////////////////////////////////////////////////////////////

/** Create TargetDetector object
 * @param threshold define the binary threshold of gray scale image
 */
TargetDetector::TargetDetector(int threshold) : m_binaryThreshold(threshold), m_ellipseThreshold(0.90)
{

}

/** Destroy TargetDetector object
 */
TargetDetector::~TargetDetector()
{

}

/** Set binary threshold
 * @param threshold define the binary threshold of gray scale image
 */
void TargetDetector::setThreshold(int threshold)
{
  m_binaryThreshold = threshold;
}

/** Get binary threshold **/
int TargetDetector::threshold() const
{
  return m_binaryThreshold;
}

/** Auto determine the optimal binary threshold
 * @param image input image
 * @param searches define the list of targets to search
 * @param step step between two threshold tests
 */
bool TargetDetector::autoThreshold(const cv::Mat &image, const std::vector<Target> &searches, const int &step)
{
  int min,max;

  //Try threshold from 0
  for(min=step;min<=255;min+=step){
    setThreshold(min);
    if(track(image,searches).size()>0)
      break;
  }
  if(min>=255)
    return false;

  //Try threshold from 255
  for(max=255-step;max>=0;max-=step){
    setThreshold(max);
    if(track(image,searches).size()>0)
      break;
  }
  if(max<0)
    return false;

  //Set and display result
  setThreshold(0.5*(min+max));
  std::cout << "New threshold = " << 0.5*(min+max) << std::endl;
  return true;
}

/** Auto determine the optimal binary threshold
 * @param image input image
 * @param search define the type of target to search
 * @param step step between two threshold tests
 */
bool TargetDetector::autoThreshold(const cv::Mat &image, Target search, const int &step)
{
  std::vector<Target> searches;
  searches.push_back(search);
  return autoThreshold(image, searches, step);
}

/** Auto determine the optimal binary threshold
 * @param capture video capture object to auto grab images
 * @param searches define types of target to search
 * @param step step between two threshold tests
 * @param nbIterationMax number of iteration before leaving the auto-thresholding function
 * @param display display video during autoThreshold process
 */
bool TargetDetector::autoThreshold(cv::VideoCapture &capture, const std::vector<Target> &searches, const int &step, const int &nbIterationMax, const bool &display)
{
  for(int i=0;i<nbIterationMax;i++){

    //Grab the image
    cv::Mat image;
    capture >> image;

    if(display){
      //Display a clone
      cv::Mat imageClone = image.clone();
      std::stringstream ss;
      ss << i << "/" << nbIterationMax;
      cv::putText(imageClone,ss.str(),cv::Point2f(5,25),cv::FONT_HERSHEY_SIMPLEX,1,CV_RGB(255,255,255),3);
      cv::putText(imageClone,ss.str(),cv::Point2f(5,25),cv::FONT_HERSHEY_SIMPLEX,1,CV_RGB(0,0,0),2);
      cv::imshow("autoThreshold",imageClone);
      char key = cv::waitKey(5);
      if(key!=-1)
	break;
    }

    if(autoThreshold(image,searches,step))
      break;
  }

  //Destroy display window
  if(display){
    cv::destroyWindow("autoThreshold");
  }
  return true;
}

/** Auto determine the optimal binary threshold
 * @param capture video capture object to auto grab images
 * @param search define the type of target to search
 * @param step step between two threshold tests
 * @param nbIterationMax number of iteration before leaving the auto-thresholding function
 * @param display display video during autoThreshold process
 */
bool TargetDetector::autoThreshold(cv::VideoCapture &capture, Target search, const int &step, const int &nbIterationMax, const bool &display)
{
  std::vector<Target> searches;
  searches.push_back(search);
  return autoThreshold(capture, searches, step, nbIterationMax, display);
}

/** Find a grid of targets in the image
 * @param image is the input image
 * @param size is the size of the grid
 * @param centers is the output vector of target centers
 * @param search is the target to search
 * @param values is an optional list of target values (use values order if provide, inscreasing order if not)
 **/
bool TargetDetector::findTargetsGrid(const cv::Mat &image, cv::Size size, std::vector<cv::Point2f> &centers, const Target &search, const std::vector<int> &values) const
{
  //Find targets
  std::vector<Target> targets = track(image,search);

  unsigned int nbTargetsNeeded = size.height*size.width;
  centers.clear();

  if(values.size()==nbTargetsNeeded){
    //Find each value
    for(unsigned int i=0; i<values.size(); i++){
      for(unsigned int j=0; j<targets.size(); j++){
	 if(targets[j].message()==values[i]){
	   centers.push_back (targets[j].center());
	   break;
	 }
      }
    }
    if(centers.size()!=values.size())
      return false;
  }else{
    std::sort(targets.begin(),targets.end(),TargetSorting());
    for(unsigned int i=0; i<targets.size(); i++){
      centers.push_back (targets[i].center());
    }
    if(targets.size() != nbTargetsNeeded)
      return false;
  }

  return true;
}

/** Track a list of targets in an image
 * @param image where to search target
 * @param searches kind of targets to search
 */
std::vector<Target> TargetDetector::track(const cv::Mat &image, const std::vector<Target> &searches) const
{
  //Convert from color to gray
  cv::Mat imageGray;
  if(image.channels()==1){
    imageGray = image.clone();
  }else{
    cv::cvtColor(image,imageGray,CV_BGR2GRAY);
  }

  //Extract blobs
  std::vector<Blob> blobs = findBlobs(imageGray);

#ifdef DEBUG_MODE
  cv::Mat color = image.clone();
  for(unsigned int i=0;i<blobs.size();i++)
    drawBlobRecursively(color,blobs[i]);
  cv::imshow("Blobs",color);
  pts_images_check.clear();
  pts_images_val.clear();
#endif

  //Call different method according to the type of target
  std::vector<Target> targets;
  for(unsigned int i=0;i<searches.size();i++){
    std::vector<Target> targetsTemp = checkBlobs(blobs,imageGray,searches[i]);
    if(targetsTemp.size()>0)
      targets.insert(targets.end(), targetsTemp.begin(), targetsTemp.end());
  }

#ifdef DEBUG_MODE
  cv::Mat ptsI = image.clone();
  for(unsigned int i=0;i<pts_images_check.size();i++){
    switch(pts_images_val[i]){
      case -1: cv::circle(ptsI,pts_images_check[i],1,CV_RGB(0,0,255),-1); break;
      case 0: cv::circle(ptsI,pts_images_check[i],1,CV_RGB(255,0,0),-1); break;
      case 1: cv::circle(ptsI,pts_images_check[i],1,CV_RGB(0,255,0),-1); break;
    }
  }
  cv::imshow("Pts",ptsI);
#endif

  return targets;
}

/** Track a target in an image
 * @param image where to search target
 * @param search kind of target to search
 */
std::vector<Target> TargetDetector::track(const cv::Mat& image, const Target &search) const
{
  std::vector<Target> searches;
  searches.push_back(search);
  return track(image,searches);
}

/** Find blobs in binary image
 * @param frameGray is input gray color image
 * @param thresholdType is the type of threshold used : CV_THRESH_BINARY (black blobs) or CV_THRESH_BINARY_INV (white blobs)
 **/
std::vector<Blob> TargetDetector::findBlobs(const cv::Mat &frameGray, const int &thresholdType) const
{
  //Segment image
  cv::Mat segmented;
  cv::threshold(frameGray, segmented, m_binaryThreshold, 255, CV_THRESH_BINARY_INV);
  //cv::Mat segmented = thresholdType==CV_THRESH_BINARY ? (frameGray > m_threshold) : (frameGray < m_threshold);

#ifdef DEBUG_MODE
  cv::imshow("Binary",segmented);
  cv::Mat segmentedClone = segmented.clone();
#endif

  std::vector<std::vector<cv::Point> > contours;
  std::vector<cv::Vec4i> hierarchy;

  cv::findContours(segmented, contours, hierarchy, CV_RETR_TREE, CV_CHAIN_APPROX_SIMPLE);

#ifdef DEBUG_MODE
  cv::Mat contoursImage;
  cv::cvtColor(segmentedClone,contoursImage,CV_GRAY2BGR);

  if(hierarchy.size()>0){
    int idx = 0;
    for( ; idx>=0; idx=hierarchy[idx][0]){
	cv::Scalar color( rand()&255, rand()&255, rand()&255 );
	cv::drawContours( contoursImage, contours, idx, color, CV_FILLED, 8, hierarchy );
    }
  }
  cv::imshow("Contours",contoursImage);
#endif

  std::vector<Blob> blobs;

  if( contours.empty() || hierarchy.empty() )
    return blobs;

  int index = 0;
  while(index >=0 ){
    bool ok;
    Blob blob = createBlob(index,contours,hierarchy,ok);
    if(ok)
      blobs.push_back(blob);
    index = hierarchy[index][0];
  }

#ifdef DEBUG_MODE
  cv::Mat ellipseImage;
  cv::cvtColor(segmentedClone,ellipseImage,CV_GRAY2BGR);

  for(unsigned int i=0;i<blobs.size();i++){
    if(blobs[i].isValid)
      cv::ellipse(ellipseImage, cv::RotatedRect(blobs[i].center,cv::Size(blobs[i].minorAxis*2.0,blobs[i].majorAxis*2.0),blobs[i].orientation), CV_RGB(255,0,0), 2 );
    for(unsigned int j=0;j<blobs[i].children.size();j++){
      if(blobs[i].children[j].isValid)
	cv::ellipse(ellipseImage, cv::RotatedRect(blobs[i].children[j].center,cv::Size(blobs[i].children[j].minorAxis*2.0,blobs[i].children[j].majorAxis*2.0),blobs[i].children[j].orientation), CV_RGB(0,255,0), 2 );
      for(unsigned int k=0;k<blobs[i].children[j].children.size();k++){
	if(blobs[i].children[j].children[k].isValid)
	  cv::ellipse(ellipseImage, cv::RotatedRect(blobs[i].children[j].children[k].center,cv::Size(blobs[i].children[j].children[k].minorAxis*2.0,blobs[i].children[j].children[k].majorAxis*2.0),blobs[i].children[j].children[k].orientation), CV_RGB(0,0,255), 2 );
      }
    }
  }

  cv::imshow("Ellipse",ellipseImage);
#endif

  return blobs;
}

/** Create blob recursively using hierarchy tree
 * @param index is the index of the blob
 * @param contours is the list of contours
 * @param hierarchy is the hierarchy tree, for each index: hierarchy[index] = [next contours at the same hierarchical level, previous contours at the same hierarchical level, the first child contour, parent contour]
 * @param ok is set to true if the blob is correct
 **/
Blob TargetDetector::createBlob(int index, const std::vector<std::vector<cv::Point> > &contours, const std::vector<cv::Vec4i> &hierarchy, bool &ok) const
{
  //hierarchy: next and previous contours at the same hierarchical level, the first child contour and the parent contour
  Blob blob = createBlob(contours[index],ok);

  if(!ok) return blob;

  blob.children.clear();

  int childIndex = hierarchy[index][2];
  while(childIndex >=0 ){
    bool childOk;
    Blob child = createBlob(childIndex,contours,hierarchy,childOk);
    if(childOk)
      blob.children.push_back(child);
    childIndex = hierarchy[childIndex][0];
  }

  //Compute area with hole
  blob.realArea = blob.area;
  for(unsigned int i=0;i<blob.children.size();i++)
    blob.realArea -= blob.children[i].area;

  blob.index = index;
  return blob;
}

/** Create blob using contours
 * @param external is the external contour
  * @param ok is set to true if the blob is correct
 **/
Blob TargetDetector::createBlob(const std::vector<cv::Point> &external, bool &ok) const
{
  Blob blob;
  blob.isValid = false;
  blob.area = cv::contourArea( external );
  blob.realArea = blob.area;

  if(blob.area<25){
    ok = false;
    return blob;
  }

  std::vector<cv::Point> hull;
  cv::convexHull(external,hull);

  if( hull.size() > 5 ){
    cv::RotatedRect ellipse = cv::fitEllipse( hull );

    if(ellipse.size.width<=0 || ellipse.size.height<=0){
      ok = false;
      std::cerr << "Bad ellipse" << std::endl;
      return blob;
    }

    blob.center = ellipse.center;
    blob.orientation = ellipse.angle;
    blob.majorAxis =  std::max(ellipse.size.width,ellipse.size.height) / 2.0;
    blob.minorAxis =  std::min(ellipse.size.width,ellipse.size.height) / 2.0;
    blob.contours = external;

    //Ellipse axis ratio too important
    if(blob.majorAxis>100.0*blob.minorAxis || blob.majorAxis<1.0){
      ok = false;
      return blob;
    }

    blob.isValid = true;
  }

  ok = true;
  return blob;
}

/** Check if blob is an ellipse using blob area
 * @param blob input blob pointer
 * @param threshold is the ellipse area threshold (area > threshold*a*b*Pi)
 **/
bool TargetDetector::isEllipse(const Blob &blob, const double &threshold) const
{
  if(threshold<=0)
    return (blob.area > m_ellipseThreshold*blob.majorAxis*blob.minorAxis*M_PI);
  else
    return (blob.area > threshold*blob.majorAxis*blob.minorAxis*M_PI);
}

/** Check if blob is an ellipse without hole using blob area
 * @param blob input blob pointer
 * @param threshold is the ellipse area threshold (area > threshold*a*b*Pi)
 **/
bool TargetDetector::isEllipseWithoutHole(const Blob &blob, const double &threshold) const
{
  //If blob is too small
  if(blob.realArea<25){
    return fabs(blob.majorAxis-blob.minorAxis) < 3;
  }

  if(threshold<=0)
    return (blob.realArea > m_ellipseThreshold*blob.majorAxis*blob.minorAxis*M_PI);
  else
    return (blob.realArea > threshold*blob.majorAxis*blob.minorAxis*M_PI);
}

/** Check if a blob is a scaled ellipse blob
 * @param reference is the reference ellipse blob
 * @param ellipse is the blob to test
 * @param maxCenterDistance is the maximum distance between centers
 * @param minFactor,maxFactor are min/max acceptable factors
 **/
bool TargetDetector::isConcentricEllispe(const Blob &reference, const Blob &ellipse, const double &maxCenterDistance, const double &minFactor, const double &maxFactor) const
{
  cv::Point2f referenceCenter = reference.center;
  cv::Point2f ellipseCenter = ellipse.center;
  double a1 = reference.majorAxis;
  double b1 = reference.minorAxis;
  double a2 = ellipse.majorAxis;
  double b2 = ellipse.minorAxis;

  return (cv::norm(referenceCenter-ellipseCenter)<maxCenterDistance && ( a1*minFactor<a2 && a2<a1*maxFactor ) && ( b1*minFactor<b2 && b2<b1*maxFactor ) );
}

/** Find center blob of target and check other blobs if needed
 * @param blobs is a list of blobs
 * @param frameGray is the grayscale image
 * @param search is the target to find
 **/
std::vector<Target> TargetDetector::checkBlobs(const std::vector<Blob> &blobs, const cv::Mat &frameGray, const Target &search) const
{
  std::vector<Target> targets;

  switch(search.type()){
    case Target::OneBlob:
      for(unsigned int i=0;i<blobs.size();i++){
	std::vector<Target> tmp = checkOneBlob(blobs[i],frameGray,search);
	if(tmp.size()>0){
	  targets.insert(targets.end(),tmp.begin(),tmp.end());
	}
      }
      break;
    case Target::ThreeBlobs:
      for(unsigned int i=0;i<blobs.size();i++){
	std::vector<Target> tmp = checkThreeBlobs(blobs[i],frameGray,search);
	if(tmp.size()>0){
	  targets.insert(targets.end(),tmp.begin(),tmp.end());
	}
      }
      break;
    case Target::TwoRings:
      for(unsigned int i=0;i<blobs.size();i++){
	std::vector<Target> tmp = checkTwoRings(blobs[i],frameGray,search);
	if(tmp.size()>0){
	  targets.insert(targets.end(),tmp.begin(),tmp.end());
	}
      }
      break;
  }

  return targets;
}

/** Check if a blob od children are the center blobs of OneBlob target type
 * @param blob is the blob to test
 * @param frameGray is the grascale image
 * @param search is the target researched
 **/
std::vector<Target> TargetDetector::checkOneBlob(const Blob &blob, const cv::Mat &frameGray, const Target &search) const
{
  std::vector<Target> targets;

  if(blob.isValid){
    if(isEllipseWithoutHole(blob)){
      //Create target and set center
      Target target = search;
      target.setCenter( blob.center );
      target.contours = blob.contours;

      //Read code if needed
      if(target.readCode()){
	if(extractCode(blob, frameGray, target)){
	  targets.push_back(target);
	}
      }else{
	targets.push_back(target);
      }
    }
  }

  for(unsigned int i=0;i<blob.children.size();i++){
    std::vector<Target> tmp = checkOneBlob(blob.children[i],frameGray,search);
    if(tmp.size()>0){
      targets.insert(targets.end(),tmp.begin(),tmp.end());
    }
  }

  return targets;
}

/** Check if a blob od children are the center blobs of ThreeBlobs target type
 * @param blob is the blob to test
 * @param frameGray is the grascale image
 * @param search is the target researched
 **/
std::vector<Target> TargetDetector::checkThreeBlobs(const Blob &blob, const cv::Mat &frameGray, const Target &search) const
{
  std::vector<Target> targets;

  if(blob.isValid){
    if(isEllipse(blob,0.8)){
      //Find internal circle
      for(unsigned int i=0;i<blob.children.size();i++){
	if(blob.children[i].index==blob.index || !blob.children[i].isValid)
	  continue;

	//if(isConcentricEllispe(blob, blob.children[i], blob.minorAxis/10.0, 1.5/3.0, 2.5/3.0)){
	if(isConcentricEllispe(blob, blob.children[i], blob.minorAxis/10.0, 1.0/3.0, 3.0/3.0)){
	  //Find central blob
	  for(unsigned int j=0;j<blob.children[i].children.size();j++){
	    if(blob.children[i].children[j].index==blob.children[i].index || blob.children[i].children[j].index==blob.index || !blob.children[i].children[j].isValid)
	      continue;

	    if(isEllipseWithoutHole(blob.children[i].children[j],0.8) && isConcentricEllispe(blob, blob.children[i].children[j], blob.minorAxis/10.0, 0.5/3.0, 1.5/3.0)){
	      //Create target and set center
	      Target target = search;
	      target.setCenter(  0.5* (blob.children[i].children[j].center + blob.center) );
	      target.contours = blob.contours;
	
	      //Read code if needed
	      if(target.readCode()){
		if(extractCode(blob, frameGray, target)){
		  targets.push_back(target);
		}
	      }else{
		targets.push_back(target);
	      }
	    }
	  }
	}
      }
    }
  }

  if(blob.children.size()>0){
    for(unsigned int i=0;i<blob.children.size();i++){
      std::vector<Target> tmp = checkThreeBlobs(blob.children[i],frameGray,search);
      if(tmp.size()>0){
	targets.insert(targets.end(),tmp.begin(),tmp.end());
      }
    }
  }
  return targets;
}

/** Check if a blob od children are the center blobs of TwoRings target type
 * @param blob is the blob to test
 * @param frameGray is the grascale image
 * @param search is the target researched
 **/
std::vector<Target> TargetDetector::checkTwoRings(const Blob &blob, const cv::Mat &frameGray, const Target &search) const
{
  std::vector<Target> targets;

  if(blob.isValid){
    if(isEllipse(blob)){
      //Find internal circle
      for(unsigned int i=0;i<blob.children.size();i++){
	if(blob.children[i].index==blob.index || !blob.children[i].isValid)
	  continue;
	if(isConcentricEllispe(blob, blob.children[i], blob.minorAxis/10.0, 4.5/6.0, 5.5/6.0)){
	  //Find central blob
	  for(unsigned int j=0;j<blob.children[i].children.size();j++){
	    //TODO, if message and header form a new circle check sub blob
	    for(unsigned int k=0;k<blob.children[i].children[j].children.size();k++){
	       if(blob.children[i].children[j].children[k].index==blob.children[i].index || blob.children[i].children[j].children[k].index==blob.index || !blob.children[i].children[j].children[k].isValid)
		continue;
	      if(isEllipseWithoutHole(blob.children[i].children[j].children[k],0.8) && isConcentricEllispe(blob, blob.children[i].children[j].children[k], blob.minorAxis/10.0, 0.5/6.0, 1.5/6.0)){
		//Create target and set center
		Target target = search;
		target.setCenter(  0.5* (blob.children[i].children[j].children[k].center + blob.center) );
		target.contours = blob.contours;

		//Read code if needed
		if(target.readCode()){
		  if(extractCode(blob, frameGray, target)){
		    targets.push_back(target);
		  }
		}else{
		  targets.push_back(target);
		}
	      }
	    }
	    //If normal
	    if(blob.children[i].children[j].index==blob.children[i].index || blob.children[i].children[j].index==blob.index || !blob.children[i].children[j].isValid)
	      continue;
	    if(isEllipseWithoutHole(blob.children[i].children[j],0.8) && isConcentricEllispe(blob, blob.children[i].children[j], blob.minorAxis/10.0, 0.5/6.0, 1.5/6.0)){
	      //Create target and set center
	      Target target = search;
	      target.setCenter(  0.5* (blob.children[i].children[j].center + blob.center) );
	      target.contours = blob.contours;

	      //Read code if needed
	      if(target.readCode()){
		if(extractCode(blob, frameGray, target)){
		  targets.push_back(target);
		}
	      }else{
		targets.push_back(target);
	      }
	    }
	  }
	}
      }
    }
  }

  if(blob.children.size()>0){
    for(unsigned int i=0;i<blob.children.size();i++){
      std::vector<Target> tmp = checkTwoRings(blob.children[i],frameGray,search);
      if(tmp.size()>0){
	targets.insert(targets.end(),tmp.begin(),tmp.end());
      }
    }
  }
  return targets;
}

/** Extract the bar code from the grayscale image using ellipse parameters
 * @param frameGray is input grayscale image
 * @param center,a,b,orientation are ellipse parameters
 * @param angle is the step between two points
 * @param codeValues are values of points
 * @param codePoints are positions of points
 * @param shift start the scanning at an other place
 **/
bool TargetDetector::extractCodeFromImage(const cv::Mat &frameGray, cv::Point2d center, double a, double b, double orientation, double angle, std::vector< int >& codeValues, std::vector< cv::Point2f >& codePoints, double shift) const
{
  int width = frameGray.size().width;
  int height =  frameGray.size().height;
  int min = 255;
  int max = 0;

  std::vector<int> codeInt;
  codePoints.clear();

  //Check if ellipse point is in the image and then read code
  double t = 360.0;
  while ( t >0.0 ){
    //Generate point on ellipse
    double xellipse_temp = (a/2.0)*cos(((t+shift)/360.0)*2.0*M_PI);
    double yellipse_temp = (b/2.0)*sin(((t+shift)/360.0)*2.0*M_PI);
    //Rotate and move point
    double xellipse = (xellipse_temp*cos(orientation))-(yellipse_temp*sin(orientation)) + center.x;
    double yellipse = (xellipse_temp*sin(orientation))+(yellipse_temp*cos(orientation)) + center.y;

    //Check if point is in the image
    if (((std::floor(xellipse) >= 0) && (std::floor(xellipse) <= width))  && ((std::floor(yellipse) >=0) && (std::floor(yellipse) <= height))){
      //Read value
      int val = frameGray.at<uchar>(std::floor(yellipse), std::floor(xellipse));

      max = std::max(max,val);
      min = std::min(min,val);

      codeInt.push_back(val);
      codePoints.push_back(cv::Point2f(xellipse,yellipse));

#ifdef DEBUG_MODE
      pts_images_check.push_back(cv::Point2f(xellipse,yellipse));
#endif

    }else{
      //Code outside of image
#ifdef DEBUG_MODE
      for(unsigned int i=0;i<codeInt.size();i++){
	pts_images_val.push_back(-1);
      }
#endif
      return false;
    }
    t -= angle;
  }

  //If not enough value range
  if(max-min<50){
#ifdef DEBUG_MODE
    for(unsigned int i=0;i<codeInt.size();i++){
      pts_images_val.push_back(-1);
    }
#endif
    return false;
  }

  codeValues.resize(codeInt.size());

  //Normalize values
  double threshold = 0.5;
  for(unsigned int i=0;i<codeInt.size();i++){
    double val = (((double)codeInt[i] - min)/(max - min));
    codeValues[i] = (val<threshold?1:0);
#ifdef DEBUG_MODE
    pts_images_val.push_back(codeValues[i]);
#endif
  }

  return true;
}

/** Extract code and read it
 * @param blob is the main blob of the target
 * @param frameGray is the inpur grayscale image
 * @param target is the input/output target
 * @return True is the extraction and reading are correct
 **/
bool TargetDetector::extractCode(const Blob &blob, const cv::Mat &frameGray, Target &target) const
{
  //Extract ellipse parameters
  double orientation = blob.orientation;
  orientation = -M_PI/2.0 + (( orientation / 360.0 ) * (2.0 * M_PI));
  double xcenter = blob.center.x;
  double ycenter = blob.center.y;
  double a = blob.majorAxis;
  double b = blob.minorAxis;
  cv::Point2d center = blob.center;

  //Check ellipse center
  int width = frameGray.size().width;
  int height =  frameGray.size().height;
  if(std::floor(xcenter) < 0 || std::floor(xcenter) >= width  || std::floor(ycenter) <0 || std::floor(ycenter) >= height){
    return false;
  }

  //Split ellipse
  int nbPointPerBit = target.numberPointPerBit();
  double codeAngle = 360.0/(target.nbBits()*nbPointPerBit);

  std::string header = target.headerString();
  std::vector<int> originalHeaderCode;
  for(unsigned int i=0;i<header.size();i++){
    for(int j=0;j<nbPointPerBit;j++){
      originalHeaderCode.push_back(header[i]=='0'?0:1);
    }
  }

  //Extract code from image
  if(target.type()==Target::OneBlob || target.type()==Target::ThreeBlobs){
    std::vector<int> code;
    std::vector<cv::Point2f> points;

    std::vector<double> scales;
    scales.push_back(1.0);
    scales.push_back(0.95);
    scales.push_back(1.05);
    scales.push_back(0.9);
    scales.push_back(1.1);


    for(unsigned int s=0;s<scales.size();s++){
      if(!extractCodeFromImage(frameGray, center, scales[s]*a*target.radiusFactor(), scales[s]*b*target.radiusFactor(), orientation, codeAngle, code, points)){
	continue;
      }

      if(!alignPointsOnFirstBit(originalHeaderCode,code,points,0.75))
	continue;
	
      //Filter the code
      //filterCode(code);

      //String code from value
      std::string codeString = code2string(code,nbPointPerBit);

      //Save points value
      std::vector<int> pointsValue;
      for(unsigned int i=0;i<points.size();i++){
	if( codeString[ (i/target.numberPointPerBit())%codeString.length() ]=='0' )
	  pointsValue.push_back(0);
	else
	  pointsValue.push_back(1);
      }
      target.setPoints( points, pointsValue );
      target.setEllipseParameters(a,b,orientation);

      if(readCode(target,codeString))
	return true;
      else
	continue;
    }

    return false;
  }

  if(target.type()==Target::TwoRings){
    std::vector<int> codeHeader, codeMessage;
    std::vector<cv::Point2f> pointsHeader, pointsMessage;

    int nbHeadetBits = target.headerBits();
    double angleHeader = 360.0/(nbHeadetBits*nbPointPerBit);
    if(!extractCodeFromImage(frameGray, center, a*target.radiusFactor(), b*target.radiusFactor(), orientation, angleHeader, codeHeader, pointsHeader))
      return false;

    //Filter the code
    //filterCode(codeHeader);

    int shift = 0;
    if(!alignPointsOnFirstBit(originalHeaderCode,codeHeader,pointsHeader,0.75,&shift))
	return false;

    //Compute rotation of header
    double rotationHeader = angleHeader * shift;

    std::string headerString = code2string(codeHeader,nbPointPerBit);
    int headerPosition = findHeader(headerString,target.headerString());

    if(headerPosition<0){
      return false;
    }

    int nbMessageBits = target.messageBits() + (target.useParityBit()?1:0);
    double angleMessage = 360.0/(nbMessageBits*nbPointPerBit);
    double rotationMessage = - (rotationHeader + headerPosition*angleHeader*nbPointPerBit);

    if(!extractCodeFromImage(frameGray, center, a*target.radiusFactor2(), b*target.radiusFactor2(), orientation, angleMessage, codeMessage, pointsMessage, rotationMessage))
      return false;

    //Filter the code
    filterCode(codeMessage);

    std::string messageString = code2string(codeMessage,nbPointPerBit);

    std::vector<int> pointsMessageValue;
    for(unsigned int i=0;i<pointsMessage.size();i++){
      if( messageString[ (i/target.numberPointPerBit())%messageString.length() ]=='0' )
	pointsMessageValue.push_back(0);
      else
	pointsMessageValue.push_back(1);
    }

    std::vector<int> pointsHeaderValue;
    for(unsigned int i=0;i<pointsHeader.size();i++){
      if( headerString[ (i/target.numberPointPerBit())%headerString.length() ]=='0' )
	pointsHeaderValue.push_back(0);
      else
	pointsHeaderValue.push_back(1);
    }

    pointsHeader.insert( pointsHeader.end(), pointsMessage.begin(), pointsMessage.end());
    pointsHeaderValue.insert( pointsHeaderValue.end(), pointsMessageValue.begin(), pointsMessageValue.end());

    target.setPoints( pointsHeader, pointsHeaderValue );

    target.setFirstHeaderIndex(0);
    target.setMessage( messageString );

    return true;
  }

  return false;
}

/** Convert code from value to string using number of value per bit
 * @param code is the input code (one number for each point)
 * @param nbPointPerBit is the number of points in a bit
 **/
std::string TargetDetector::code2string(const std::vector<int> &code, const int &nbPointPerBit) const
{
  std::string codeString;
  double val = 0;
  for(unsigned int i=0;i<code.size();i++){
    val += code[i];
    if(i!=0 && (i%nbPointPerBit==0 || i==code.size()-1)){
      codeString.push_back( (val/nbPointPerBit>0.5?'1':'0') );
      val = 0;
    }
  }
  return codeString;
}

/** Filter the code to avoid wrong bit segmentation
 * @param code is input/output cade value vector
 **/
void TargetDetector::filterCode(std::vector<int> &code) const
{
  int firstVal = code[0];

  //First filter, use (p)revious, (c)urrent and (n)ext value
  for(unsigned int c=0;c<code.size();c++){
    int p = (c-1+code.size())%code.size();
    int n = (c+1+code.size())%code.size();
    if(code[p]==firstVal && code[c]!=firstVal && code[n]==firstVal)
      code[c] = firstVal;
  }
  for(unsigned int c=0;c<code.size();c++){
    int p = (c-1+code.size())%code.size();
    int n = (c+1+code.size())%code.size();
    if(code[p]!=firstVal && code[c]==firstVal && code[n]!=firstVal)
      code[c] = !firstVal;
  }

  //Second filter, use (p)revious, (c)urrent, (n)ext  and (s)secondNext value
  for(unsigned int c=0;c<code.size();c++){
    int p = (c-1+code.size())%code.size();
    int n = (c+1+code.size())%code.size();
    int s = (c+2+code.size())%code.size();
    if(code[p]==firstVal && code[c]!=firstVal && code[n]!=firstVal && code[s]==firstVal){
      code[c] = firstVal;
      code[n] = firstVal;
    }
  }
  for(unsigned int c=0;c<code.size();c++){
    int p = (c-1+code.size())%code.size();
    int n = (c+1+code.size())%code.size();
    int s = (c+2+code.size())%code.size();
    if(code[p]!=firstVal && code[c]==firstVal && code[n]==firstVal && code[s]!=firstVal){
      code[c] = !firstVal;
      code[n] = !firstVal;
    }
  }
}

/** Align the code on the first point in the sequence
 * @param header is input header (code form)
 * @param code is input/ouput code value vector
 * @param points is input/output code point vector
 * @param minScore is the threshold to decide if the header match on code
 * @param shift is optional output represent the code shift to match the header
 * @return Return true if the header match on code
 **/
bool TargetDetector::alignPointsOnFirstBit(const std::vector<int> &header, std::vector<int> &code, std::vector<cv::Point2f> &points, const double &minScore, int *shift) const
{
  std::vector<int> codeDoubleLength = code;
  codeDoubleLength.insert(codeDoubleLength.end(),code.begin(),code.end());

  int bestPos = 0;
  double bestRatio = 0;
  int length = header.size();
  for(unsigned int i=0;i<code.size();i++){
    int nb = 0;
    for(unsigned int j=0;j<header.size();j++){
      if(codeDoubleLength[i+j]==header[j]){
	nb++;
      }
    }
    double ratio = (double)nb/(double)length;
    if(ratio>bestRatio){
      bestRatio = ratio;
      bestPos = i;
    }
  }

  if(bestRatio<minScore){
    if(shift)
      (*shift) = 0;
    return false;
  }

  std::vector<int> codeCopy;
  for(unsigned int i=bestPos;i<code.size();i++)
    codeCopy.push_back(code[i]);
  for(int i=0;i<bestPos;i++)
    codeCopy.push_back(code[i]);
  code = codeCopy;

  std::vector<cv::Point2f> pointsCopy;
  for(unsigned int i=bestPos;i<points.size();i++)
    pointsCopy.push_back(points[i]);
  for(int i=0;i<bestPos;i++)
    pointsCopy.push_back(points[i]);
  points = pointsCopy;

    if(shift)
      (*shift) = bestPos;
  return true;
}

/** Find the header bit sequence in the scanned code
 * @param str is the input code (without begin or end)
 * @param headerString is the header sequence
 * @return The position of header in code sequence, -1 if not found.
 **/
int TargetDetector::findHeader(const std::string &str, const std::string &headerString) const
{
  std::string code = str + str;
  std::size_t found = code.find(headerString);

  if(found==std::string::npos){
    return -1;
  }else{
    return found;
  }
}

/** Read code from code sequence nits
 * @param target input/output target
 * @param codeString code of the target under string format
 **/
bool TargetDetector::readCode(Target &target, const std::string &codeString) const
{
  //Search the header
  std::string header = target.headerString();
  int found = findHeader(codeString,header);
  if(found<0){
#ifdef DEBUG_MODE
    //std::cerr << "codeString: " << codeString << "   looking for: " << header << std::endl;
#endif
    return false;
  }

  std::string str = codeString + codeString;
  std::string message;

  message = str.substr(found+header.size(),target.nbBits()-header.size());

  //Check parity bit
  if(target.useParityBit()){
    int parityBit = message[message.size()-1]=='1'?1:0;
    message = message.substr(0,message.length()-1);

    if(target.inverseParityBit())
      parityBit = !parityBit;

    std::stringstream ss;
    ss << header << message;
    int codeParity = parity(ss.str());

    if(parityBit!=codeParity){
      //Wrong parity bits
#ifdef DEBUG_MODE
      std::cerr << "Wrong parity bit" << std::endl;
#endif
      return false;
    }
  }

  target.setMessage(message,target.useGrayCode());
  target.setFirstHeaderIndex(found*target.numberPointPerBit());

  return true;
}

/** Draw target on image
 * @param image input/output image
 * @param target target to draw
 **/
void TargetDetector::drawTarget(cv::Mat &image, const Target &target)
{
  std::vector<Target> targets;
  targets.push_back(target);
  drawTargets(image,targets);
}

/** Draw target on image
 * @param image input/output image
 * @param targets list of targets to tdraw
 **/
void TargetDetector::drawTargets(cv::Mat &image, const std::vector<Target> &targets)
{
  for(unsigned int i=0;i<targets.size();i++){
    //Ellipse center
    cv::circle(image,targets[i].center(),2,CV_RGB(255,0,0),-1);

    //Draw contour
    std::vector< std::vector<cv::Point> > contours;
    contours.push_back( targets[i].contours );
    cv::drawContours(image,contours,0,CV_RGB(255,0,0));

    if(!targets[i].readCode())
      continue;

    //Header points
    for(int j=0;j<targets[i].headerBits()*targets[i].numberPointPerBit();j++){
      int index = targets[i].firstHearderIndex() + j;
      cv::circle(image,targets[i].point( index%targets[i].nbPoints() ),4,CV_RGB(255,0,255),-1);
    }

    //Orientation
    cv::line(image,targets[i].center(),targets[i].point(targets[i].firstHearderIndex()),CV_RGB(0,0,255),3);

    //Points with value
    for(int j=0;j<targets[i].nbPoints();j++){
      if(targets[i].pointValue(j)==0)
	cv::circle(image,targets[i].point(j),2,CV_RGB(0,255,255),-1);
      else
	cv::circle(image,targets[i].point(j),2,CV_RGB(255,255,0),-1);
    }

    //Lines between bits
    if(targets[i].type()==Target::TwoRings){
      drawLineBetweenPoints(image,targets[i].point(0),targets[i].point(targets[i].headerBits()*targets[i].numberPointPerBit()-1));
      for(int j=1;j<targets[i].headerBits();j++){
	drawLineBetweenPoints(image,targets[i].point(targets[i].numberPointPerBit()*j-1),targets[i].point(targets[i].numberPointPerBit()*j));
      }
      drawLineBetweenPoints(image,targets[i].point(targets[i].headerBits()*targets[i].numberPointPerBit()),targets[i].point(targets[i].nbPoints()-1));
      for(int j=targets[i].headerBits()+1;j<targets[i].nbBits();j++){
	drawLineBetweenPoints(image,targets[i].point(targets[i].numberPointPerBit()*j-1),targets[i].point(targets[i].numberPointPerBit()*j));
      }
    }else{
      drawLineBetweenPoints(image,targets[i].point(0),targets[i].point(targets[i].nbPoints()-1));
      for(int j=1;j<targets[i].nbBits();j++){
	drawLineBetweenPoints(image,targets[i].point(targets[i].numberPointPerBit()*j-1),targets[i].point(targets[i].numberPointPerBit()*j));
      }
    }

    //Print message value
    if(targets[i].message()>=0){
      std::stringstream ss;
      ss << targets[i].message();
      cv::putText(image,ss.str(),targets[i].center(),cv::FONT_HERSHEY_SIMPLEX,1,CV_RGB(0,0,0),2);
      cv::putText(image,ss.str(),targets[i].center(),cv::FONT_HERSHEY_SIMPLEX,1,CV_RGB(255,255,255),1);
    }
  }
}

/** Draw red line to separate two points,
 * @param image is input/output image
 * @param p1,p2 are 2D points
 **/
void TargetDetector::drawLineBetweenPoints(cv::Mat &image, const cv::Point2f &p1, const cv::Point2f &p2)
{
  cv::Point2f n = p2-p1;
  std::swap(n.x,n.y);
  n.x *= -1;
  cv::Point2f p = 0.5 * (p1+p2);
  cv::line(image,p+n,p-n,CV_RGB(255,0,0),2);
}

/** Draw grid points on image
 * @param image is input/output image
 * @param size is size of the grid
 * @param centers are input centers
 * @param found set to true if grid was found
 */
void TargetDetector::drawTargetGrid(const cv::Mat &image, const cv::Size &size, const std::vector< cv::Point2f > &centers, const bool &found)
{
  cv::drawChessboardCorners(image,size,centers,found);
}




///////////////////////////////////// DEBUG FUNCTIONS ///////////////////////////////////////////



#ifdef DEBUG_MODE
bool isEllipse(const Blob &blob, const double &threshold)
{
  return (blob.area > threshold*blob.majorAxis*blob.minorAxis*M_PI);
}

void drawBlobRecursively(cv::Mat &image, const Blob &blob, int depth)
{
  cv::Scalar color1,color2;

  switch(depth%6){
    case 0:
      color1 = CV_RGB(255,0,0);
      color2 = CV_RGB(125,0,0);
      break;
    case 1:
      color1 = CV_RGB(0,255,0);
      color2 = CV_RGB(0,125,0);
      break;
    case 2:
      color1 = CV_RGB(0,0,255);
      color2 = CV_RGB(0,0,125);
      break;
    case 3:
      color1 = CV_RGB(255,255,0);
      color2 = CV_RGB(125,125,0);
      break;
    case 4:
      color1 = CV_RGB(0,255,255);
      color2 = CV_RGB(0,125,125);
      break;
    case 5:
      color1 = CV_RGB(255,0,255);
      color2 = CV_RGB(125,0,125);
      break;
    default:
      color1 = CV_RGB(255,255,255);
      color2 = CV_RGB(125,125,125);
  }


  std::vector< std::vector<cv::Point> > contours;
  contours.push_back( blob.contours );

  if(blob.isValid){
    if(isEllipse(blob,0.90))
      cv::drawContours(image,contours,0,color1,2);
    else
      cv::drawContours(image,contours,0,color2,2);
  }else{
    if(isEllipse(blob,0.90))
      cv::drawContours(image,contours,0,color1,-1);
    else
      cv::drawContours(image,contours,0,color2,-1);
  }

  for(unsigned int i=0;i<blob.children.size();i++){
    drawBlobRecursively(image,blob.children[i],depth+1);
  }
}

void drawEllipseRecursively(cv::Mat &image, const Blob &blob, int depth)
{
  cv::Scalar color;

  switch(depth%6){
    case 0:
      color = CV_RGB(255,0,0);
      break;
    case 1:
      color = CV_RGB(0,255,0);
      break;
    case 2:
      color = CV_RGB(0,0,255);
      break;
    case 3:
      color = CV_RGB(255,255,0);
      break;
    case 4:
      color = CV_RGB(0,255,255);
      break;
    case 5:
      color = CV_RGB(255,0,255);
      break;
    default:
      color = CV_RGB(255,255,255);
  }

  cv::ellipse(image, cv::RotatedRect(blob.center,cv::Size(blob.minorAxis*2.0,blob.majorAxis*2.0),blob.orientation), color, 2 );

  for(unsigned int i=0;i<blob.children.size();i++){
    drawEllipseRecursively(image,blob.children[i],depth+1);
  }
}

#endif
