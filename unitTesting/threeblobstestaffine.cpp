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

#include "targetdetector.h"
#include "stdlib.h"

cv::Mat computeAffineTransform(const cv::Mat &input, cv::Mat &output, const double ratio, const double angle, const cv::Scalar color = cv::Scalar(125,125,125))
{
  cv::Point2f pt(input.cols/2., input.rows/2.);
  cv::RotatedRect rect(pt,input.size(),angle);
  cv::Mat transform = cv::getRotationMatrix2D(pt, angle, 1.0);
  
  cv::Rect bounding = rect.boundingRect();
  transform.at<double>(0,2) -= bounding.x;
  transform.at<double>(1,2) -= bounding.y;
  
  for(int i=0;i<3;i++)
    transform.at<double>(1,i) *= ratio;
  
  cv::Size size = bounding.size();
  size.height *= ratio;
  
  cv::warpAffine(input,output,transform,size,cv::INTER_LINEAR,cv::BORDER_CONSTANT,color);
  
  return transform;
}

cv::Point2f applyTransform(const cv::Point2f &in, const cv::Mat &transform)
{
  std::vector<cv::Point2f> pin,pout;
  pin.push_back(in);
  cv::transform(pin,pout,transform);
  return pout.at(0);
}

int main(int argc, char* argv[])
{
  
  cv::Mat I = cv::imread(argv[1]), T;
  cv::Point2f center(I.cols/2,I.rows/2);
  std::vector<Target> targets;
  cv::Mat gaussianNoise;
  
  //Define target
  TargetDetector detector;
  Target search(Target::ThreeBlobs,8,8,false,180);
  
  //Test with two value to display the result
  computeAffineTransform(I,T,0.5,45);
  gaussianNoise = T.clone();
  cv::randn(gaussianNoise,128,30);
  T = 0.5*T + 0.5*gaussianNoise; 
  targets = detector.track(T,search);
  detector.drawTargets(T,targets);
  cv::imwrite("ThreeBlobsAffine.png",T);
  
  //Scales to test
  std::vector<double> scales;
  for(double s=1.0;s>0.2;s-=0.1)
    scales.push_back(s);
  
  //Angles to test
  std::vector<double> angles;
  for(double a=0;a<90;a+=10)
    angles.push_back(a);
  
  //Test all possibilities
  for(unsigned int i=0;i<scales.size();i++){
    for(unsigned int j=0;j<angles.size();j++){
      
      //Apply affine transform on image
      cv::Mat transform = computeAffineTransform(I,T,scales[i],angles[j]);
      gaussianNoise = T.clone();
      cv::randn(gaussianNoise,128,30);
      T = 0.5* T + 0.5*gaussianNoise; 
      
      //Detect target
      targets = detector.track(T,search); 
      
      if(targets.size()!=1){
#ifdef DEBUG_MODE  
	detector.drawTargets(T,targets);
	cv::imshow("Image",T);
	cv::waitKey(0);
#endif
	return EXIT_FAILURE;
      }
      
      //Project center and compute distance
      cv::Point2f p = applyTransform(center,transform);
      double distance = cv::norm(p-targets[0].center());
      
#ifdef DEBUG_MODE  
      std::cout << "Scale: " << scales[i]  << std::endl;
      std::cout << "Angle: " << angles[j]  << std::endl;
      std::cout << "Distance: " << distance << std::endl << std::endl;
      
      detector.drawTargets(T,targets);
      cv::circle(T,p,3,CV_RGB(255,0,0),-1);
      cv::circle(T,targets[0].center(),3,CV_RGB(0,255,0),-1);
      cv::imshow("Image",T);
      cv::waitKey(0);
#endif
      
      if(distance>2)
	return EXIT_FAILURE;      
    }
  }

  return EXIT_SUCCESS;
}
