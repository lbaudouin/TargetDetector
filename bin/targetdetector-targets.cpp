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

#include <iostream>
#include <targetdetector.h>

#include "timer.h"

/**
* @author Léo Baudouin\n
* @em baudouin.leo@gmail.com
*/

void help(std::string exec);
bool isNumber(const std::string &str, int &nb);
std::string toStr(double d);

int main(int argc, char* argv[])
{
  cv::VideoCapture capture;
  bool isVideo = true;

  if(argc>1){
    //Try to read camera index
    int cameraIndex = 0;
    std::istringstream iss(argv[1]);
    iss >> std::ws >> cameraIndex >> std::ws;
    if(iss.eof()){
      //Open camera
      capture.open(cameraIndex);
    }else{
      //Open file
      capture.open(argv[1]);
      if(capture.get(CV_CAP_PROP_FRAME_COUNT)==1)
	isVideo = false;
    }
  }else{
    //Default
    capture.open(CV_CAP_ANY);
  }

  if(!capture.isOpened()){
    std::cerr << "Failed to open video capture" << std::endl;
    return 1;
  }


  //Create Detector
  TargetDetector targetDetector;
  cv::Mat image;

  if(!isVideo)
    capture >> image;

  //Set targets
  std::vector<Target> searches;
  searches.push_back( Target(Target::OneBlob,8,8,false,180) );
  searches.push_back( Target(Target::ThreeBlobs,8,8,false,180) );
  searches.push_back( Target(Target::ThreeBlobs,8,8,false,220) );
  searches.push_back( Target(Target::TwoRings,8,8,false,180) );

  //Read and set threshold
  int threshold = 125;
  if(argc>2){
    if(std::string(argv[2])=="auto"){
      if(isVideo)
	targetDetector.autoThreshold(capture,searches[0],10,250,true);
      else{
	targetDetector.autoThreshold(image,searches[0],10,250,true);
      }
    }else{
      if(isNumber(argv[2],threshold)){
	targetDetector.setThreshold(threshold);	
      }else{
	std::cerr << "'" << argv[2] << "' is not a number in [0-255] or 'auto'" << std::endl;
	help(argv[0]);
	return 1;
      }
    }
  }


  Timer fpsTimer,execTimer;
  char key = -1;

  int nbIter = 0;
  double execTimeMean = 0;
  double FPSmean = 0;

  do{
    //Grab image
    if(isVideo)
      capture >> image;

    //Check if the image is correct
    if(image.empty()){
      break;
    }

    execTimer.restart();

    //Detect targets
    std::vector<Target> targets = targetDetector.track(image,searches);

    double execTime = execTimer.elapsed();

    //Draw targets
    targetDetector.drawTargets(image,targets);

    //Compute and display FPS
    double fps = 1.0/fpsTimer.s_elapsed();

    fpsTimer.restart();
    cv::putText(image,toStr(fps),cv::Point2f(5,25),cv::FONT_HERSHEY_SIMPLEX,1,CV_RGB(255,255,255),3);
    cv::putText(image,toStr(fps),cv::Point2f(5,25),cv::FONT_HERSHEY_SIMPLEX,1,CV_RGB(0,0,0),2);
    cv::putText(image,toStr(execTime),cv::Point2f(5,50),cv::FONT_HERSHEY_SIMPLEX,1,CV_RGB(255,255,255),3);
    cv::putText(image,toStr(execTime),cv::Point2f(5,50),cv::FONT_HERSHEY_SIMPLEX,1,CV_RGB(0,0,0),2);
    cv::putText(image,toStr(targets.size()),cv::Point2f(5,75),cv::FONT_HERSHEY_SIMPLEX,1,CV_RGB(255,255,255),3);
    cv::putText(image,toStr(targets.size()),cv::Point2f(5,75),cv::FONT_HERSHEY_SIMPLEX,1,CV_RGB(0,0,0),2);

    //Display image
    cv::imshow("Detector",image);
    key = cv::waitKey(5);

    //Compute means
    execTimeMean = (nbIter*execTimeMean + execTime)/(nbIter+1);
    FPSmean = (nbIter*FPSmean + fps)/(nbIter+1);
    nbIter++;

#ifdef TARGET_DEBUG
    //Pause if more than 2
    if(targets.size()>2)
      cv::waitKey(0);
#endif
  }while(key<0 && isVideo);

  std::cout << "Mean FPS: " << FPSmean << "Hz" << std::endl;
  std::cout << "Mean exec: " << execTimeMean << "ms" << std::endl;

  if(!isVideo)
    cv::waitKey(0);

  return 0;
}


void help(std::string exec)
{

}

std::string toStr(double d)
{
    std::stringstream ss;
    ss << d;
    return ss.str();
}

bool isNumber(const std::string &number, int &nb)
{
  std::stringstream ss(number);
  ss >> nb;
  if(ss.fail())
    return false;
  if(nb<0 || nb>255)
    return false;
  return true;
}
