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

#include "../include/targetdetector.h"

/** Create TargetDetector object
 * @param threshold define the binary threshold of gray scale image
 */
TargetDetector::TargetDetector(int threshold) : m_binaryThreshold(threshold), m_ellipseThreshold(0.95)
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
 */
bool TargetDetector::autoThreshold(cv::VideoCapture &capture, const std::vector<Target> &searches, const int &step, const int &nbIterationMax, const bool display)
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
 */
bool TargetDetector::autoThreshold(cv::VideoCapture &capture, Target search, const int &step, const int &nbIterationMax, const bool display)
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

#ifdef TARGET_DEBUG
  cv::Mat color = image.clone();
  for(unsigned int i=0;i<blobs.size();i++){
    std::vector< std::vector<cv::Point> > contours;
    contours.push_back( blobs[i].contours );
    if(isEllipse(blobs[i]))
      cv::drawContours(color,contours,0,CV_RGB(255,0,0),2);
    else
      cv::drawContours(color,contours,0,CV_RGB(255,255,255),2);
    for(unsigned int j=0;j<blobs[i].children.size();j++){
      contours.clear();
      contours.push_back( blobs[i].children[j].contours );
      if(isEllipse(blobs[i].children[j]))
	cv::drawContours(color,contours,0,CV_RGB(0,255,0),2);
      else
	cv::drawContours(color,contours,0,CV_RGB(255,255,0),2);
      for(unsigned int k=0;k<blobs[i].children[j].children.size();k++){
	contours.clear();
	contours.push_back( blobs[i].children[j].children[k].contours );
	if(isEllipse(blobs[i].children[j].children[k]))
	  cv::drawContours(color,contours,0,CV_RGB(0,0,255),2);
	else
	  cv::drawContours(color,contours,0,CV_RGB(0,255,255),2);
      }
    }
  }
  cv::imshow("color",color);
#endif

  //Call different method according to the type
  std::vector<Target> targets;
  for(unsigned int i=0;i<searches.size();i++){
    std::vector<Target> targetsTemp = checkBlobs(blobs,imageGray,searches[i]);
    if(targetsTemp.size()>0)
      targets.insert(targets.end(), targetsTemp.begin(), targetsTemp.end());
  }

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

#ifdef TARGET_DEBUG
  cv::imshow("Binary",segmented);
  cv::Mat segmentedClone = segmented.clone();
#endif

  std::vector<std::vector<cv::Point> > contours;
  std::vector<cv::Vec4i> hierarchy;

  cv::findContours(segmented, contours, hierarchy, CV_RETR_TREE, CV_CHAIN_APPROX_SIMPLE);

#ifdef TARGET_DEBUG
  cv::Mat contoursImage;
  cv::cvtColor(segmentedClone,contoursImage,CV_GRAY2BGR);
  
  int idx = 0;
  for( ; idx>=0; idx=hierarchy[idx][0]){
      cv::Scalar color( rand()&255, rand()&255, rand()&255 );
      cv::drawContours( contoursImage, contours, idx, color, CV_FILLED, 8, hierarchy );
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
  
#ifdef TARGET_DEBUG
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

  if(blob.area<5){
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

  if(blob.children.size()>0){
    for(unsigned int i=0;i<blob.children.size();i++){
      std::vector<Target> tmp = checkOneBlob(blob.children[i],frameGray,search);
      if(tmp.size()>0){
	targets.insert(targets.end(),tmp.begin(),tmp.end());
      }
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
    if(isEllipse(blob)){
      //Find internal circle
      for(unsigned int i=0;i<blob.children.size();i++){
	if(blob.children[i].index==blob.index || !blob.children[i].isValid)
	  continue;
	if(isConcentricEllispe(blob, blob.children[i], blob.minorAxis/10.0, 1.5/3.0, 2.5/3.0)){
	  //Find central blob
	  for(unsigned int j=0;j<blob.children[i].children.size();j++){
	    if(blob.children[i].children[j].index==blob.children[i].index || blob.children[i].children[j].index==blob.index || !blob.children[i].children[j].isValid)
	      continue;
	
	    if(isEllipseWithoutHole(blob.children[i].children[j],0.9) && isConcentricEllispe(blob, blob.children[i].children[j], blob.minorAxis/10.0, 0.5/3.0, 1.5/3.0)){
	      //Create target and set center
	      Target target = search;
	      target.setCenter(  0.5* (blob.children[i].children[j].center + blob.center) );
	      target.contours = blob.contours;
	
	      //Read code if needed
	      if(target.readCode()){
		if(extractCode(blob, frameGray, target)){
		  targets.push_back(target);
		}else{
#ifdef TARGET_DEBUG
		  std::cerr << "checkThreeBlobs: extractCode failed" << std::endl;
#endif
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
	      if(isEllipseWithoutHole(blob.children[i].children[j].children[k],0.9) && isConcentricEllispe(blob, blob.children[i].children[j].children[k], blob.minorAxis/10.0, 0.5/6.0, 1.5/6.0)){
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
	    if(isEllipseWithoutHole(blob.children[i].children[j],0.5) && isConcentricEllispe(blob, blob.children[i].children[j], blob.minorAxis/10.0, 0.5/6.0, 1.5/6.0)){
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

  std::vector<double> codeDouble;
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

      codeDouble.push_back(val);
      codePoints.push_back(cv::Point2f(xellipse,yellipse));

    }else{
      //Code outside of image
      return false;
    }
    t -= angle;
  }

  codeValues.resize(codeDouble.size());

  //Normalize values
  double threshold = 0.5;
  for(unsigned int i=0;i<codeDouble.size();i++){
    double val = ((codeDouble[i] - min)/(max - min));
    codeValues[i] = (val<threshold?1:0);
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

  //Extract code from image
  if(target.type()==Target::OneBlob || target.type()==Target::ThreeBlobs){
    std::vector<int> code;
    std::vector<cv::Point2f> points;

    if(!extractCodeFromImage(frameGray, center, a*target.radiusFactor(), b*target.radiusFactor(), orientation, codeAngle, code, points)){
      return false;
    }

    //Filter the code
    filterCode(code);
    alignPointsOnBit(code,points);

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

    return readCode(target,codeString);
  }
  if(target.type()==Target::TwoRings){
    std::vector<int> codeHeader, codeMessage;
    std::vector<cv::Point2f> pointsHeader, pointsMessage;

    int nbHeadetBits = target.headerBits();
    double angleHeader = 360.0/(nbHeadetBits*nbPointPerBit);
    if(!extractCodeFromImage(frameGray, center, a*target.radiusFactor(), b*target.radiusFactor(), orientation, angleHeader, codeHeader, pointsHeader))
      return false;

    //Filter the code
    filterCode(codeHeader);

    //Shift using angle not number of points
    double rotAlign = alignPointsOnBit(codeHeader,pointsHeader,angleHeader);
    std::string headerString = code2string(codeHeader,nbPointPerBit);

    int headerPos = findHeader(headerString,target.headerString());

    if(headerPos<0){
      return false;
    }

    std::vector<int> codeHeaderTemp;
    std::vector<cv::Point2f> pointsHeaderTemp;
    for(unsigned int i=headerPos*nbPointPerBit;i<codeHeader.size();i++){
      codeHeaderTemp.push_back(codeHeader[i]);
      pointsHeaderTemp.push_back(pointsHeader[i]);
    }
    for(int i=0;i<headerPos*nbPointPerBit;i++){
      codeHeaderTemp.push_back(codeHeader[i]);
      pointsHeaderTemp.push_back(pointsHeader[i]);
    }
    codeHeader = codeHeaderTemp;
    pointsHeader = pointsHeaderTemp;

    std::string str;
    for(unsigned int i=0;i<headerString.size();i++){
      int index = (i+headerPos+headerString.size())%headerString.size();
      str.push_back( headerString[index] );
    }
    headerString = str;

    int nbMessageBits = target.messageBits() + (target.useParityBit()?1:0);
    double angleMessage = 360.0/(nbMessageBits*nbPointPerBit);
    double shift = - (rotAlign + headerPos*angleHeader*nbPointPerBit);
    if(!extractCodeFromImage(frameGray, center, a*target.radiusFactor2(), b*target.radiusFactor2(), orientation, angleMessage, codeMessage, pointsMessage, shift))
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
 * @param code is input/ouput code value vector
 * @param points is input/output code point vector
 * @return Number of points shifted
 **/
int TargetDetector::alignPointsOnBit(std::vector<int> &code, std::vector<cv::Point2f> &points) const
{
  int firstVal = code[0];

  //Copy the first numbers at the end until value change
  int shift=0;
  while(code[shift]==firstVal)
    shift++;

  //Shift the code
  std::vector<int> codeCopy(code.size());
  for(unsigned int i=0;i<codeCopy.size();i++){
    int index = (i-shift+code.size())%code.size();
    codeCopy[index] = code[i];
  }
  code = codeCopy;

  //Shift points
  std::vector<cv::Point2f> pointsCopy;
  for(unsigned int i=shift;i<points.size();i++)
    pointsCopy.push_back(points[i]);
  for(int i=0;i<shift;i++)
    pointsCopy.push_back(points[i]);
  points = pointsCopy;

  return shift;
}

/** Align the code on the first point in the sequence
 * @param code is input/ouput code value vector
 * @param points is input/output code point vector
 * @param angle is angle between two points
 * @return Angle shifted
 **/
double TargetDetector::alignPointsOnBit(std::vector<int> &code, std::vector<cv::Point2f> &points, const double &angle) const
{
  int shift = alignPointsOnBit(code,points);
  return angle*shift;
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
  std::string str = codeString + codeString;

  std::string message;

  if(found>=0){
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
#ifdef TARGET_DEBUG
	std::cerr << "Wrong parity bit" << std::endl;
#endif
	return false;
      }
    }

    target.setMessage(message,target.useGrayCode());
    target.setFirstHeaderIndex(found*target.numberPointPerBit());

    return true;
  }else{
    return false;
  }
}

/** Draw target on image
 * @param image input/output image
 * @param target target to tdraw
 **/
void TargetDetector::drawTarget(cv::Mat &image, const Target &target) const
{
  std::vector<Target> targets;
  targets.push_back(target);
  drawTargets(image,targets);
}

/** Draw target on image
 * @param image input/output image
 * @param targets list of targets to tdraw
 **/
void TargetDetector::drawTargets(cv::Mat &image, const std::vector<Target> &targets) const
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
void TargetDetector::drawLineBetweenPoints(cv::Mat &image, const cv::Point2f &p1, const cv::Point2f &p2) const
{
  cv::Point2f n = p2-p1;
  std::swap(n.x,n.y);
  n.x *= -1;
  cv::Point2f p = 0.5 * (p1+p2);
  cv::line(image,p+n,p-n,CV_RGB(255,0,0),2);
}
