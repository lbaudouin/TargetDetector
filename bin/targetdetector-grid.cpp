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
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>

#include "timer.h"

/**
* @author Léo Baudouin\n
* @em baudouin.leo@gmail.com
*/

void help(std::string exec);
bool isNumber(const std::string &str, int &nb);
std::string toStr(double d);

#ifndef DISABLE_XML 
void createDefaultGridConfig(std::string filepath);
#endif

int main(int argc, char* argv[])
{  
  cv::VideoCapture capture;
  bool isVideo = false;

  std::string targetType = "3B";
  int headerBits = 8;
  int headerValue = 180;
  int messageBits = 8;
  bool useParityBit = false;
  std::string thresholdString; 
  cv::Size targetSize(4,5);
  std::vector<int> values;
  
  for(int i=1;i<argc;i++){
    if(!strcmp(argv[i],"-t") || !strcmp(argv[i],"--threshold")) {thresholdString = argv[++i]; continue;}
#ifndef DISABLE_XML 
    if(!strcmp(argv[i],"-c") || !strcmp(argv[i],"--config")) {
      cv::FileStorage fs(argv[++i], cv::FileStorage::READ);
      fs["targetType"] >> targetType;
      fs["headerBits"] >> headerBits;
      fs["headerValue"] >> headerValue;
      fs["messageBits"] >> messageBits;
      fs["useParityBit"] >> useParityBit;
      fs["targetSize"] >> targetSize;
      fs["values"] >> values;
      continue;
    }
    if(!strcmp(argv[i],"-g") || !strcmp(argv[i],"--generate")) {
      createDefaultGridConfig(argv[++i]);
      return 0;
    }
#endif
    if(!strcmp(argv[i],"-h") || !strcmp(argv[i],"--help") || !strcmp(argv[i],"-?")) {help(argv[0]); return 1;}
    
    //Try to read camera index
    int cameraIndex = 0;
    std::istringstream iss(argv[i]);
    iss >> cameraIndex;
    
    if(!iss.fail()){
      //Open camera
      int cameraIndex = atoi(argv[i]);
      capture.open(cameraIndex);
      isVideo = true;
    }else{
      //Open file
      capture.open(argv[i]);
      if(capture.get(CV_CAP_PROP_FRAME_COUNT)>1)
	isVideo = true;
    }
  }
  
  
  if(!capture.isOpened()){
      capture.open(CV_CAP_ANY);
	isVideo = true;
  }
  
  if(!capture.isOpened()){
    std::cerr << "Failed to open video capture" << std::endl;
    return 1;
  }
  
  std::cout << "Read from video: " << (isVideo?"true":"false") << std::endl;

  
  //Create Detector
  TargetDetector targetDetector;
  cv::Mat image;

  if(!isVideo)
    capture >> image;

  Target::Type type = Target::ThreeBlobs;
  if(targetType=="1B"){
    type = Target::OneBlob;
  }else if(targetType=="2R"){
    type = Target::TwoRings;
  }
    
  //Set target
  Target target = Target(type,headerBits,messageBits,useParityBit,headerValue);
  std::cout << "Target: " << target << std::endl;
  
  //Read and set threshold
  int threshold = 125;
  if(thresholdString!=""){
    if(thresholdString=="auto" || thresholdString=="Auto"){
      if(isVideo)
	targetDetector.autoThreshold(capture,target,10,250,true);
      else{
	targetDetector.autoThreshold(image,target,10);
      }
    }else{
      if(isNumber(thresholdString,threshold)){
	targetDetector.setThreshold(threshold);	
      }else{
	std::cerr << "'" << thresholdString << "' is not a number in [0-255] or 'auto'" << std::endl;
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
    std::vector<cv::Point2f> centers;

    bool found = targetDetector.findTargetsGrid(image,targetSize,centers,target,values);

    double execTime = execTimer.elapsed();

    //Draw grid
    cv::drawChessboardCorners(image,targetSize,centers,found);

    //Compute and display FPS
    double fps = 1.0/fpsTimer.s_elapsed();

    fpsTimer.restart();
    cv::putText(image,toStr(fps),cv::Point2f(5,25),cv::FONT_HERSHEY_SIMPLEX,1,CV_RGB(255,255,255),3);
    cv::putText(image,toStr(fps),cv::Point2f(5,25),cv::FONT_HERSHEY_SIMPLEX,1,CV_RGB(0,0,0),2);
    cv::putText(image,toStr(execTime),cv::Point2f(5,50),cv::FONT_HERSHEY_SIMPLEX,1,CV_RGB(255,255,255),3);
    cv::putText(image,toStr(execTime),cv::Point2f(5,50),cv::FONT_HERSHEY_SIMPLEX,1,CV_RGB(0,0,0),2);
    cv::putText(image,toStr(centers.size()),cv::Point2f(5,75),cv::FONT_HERSHEY_SIMPLEX,1,CV_RGB(255,255,255),3);
    cv::putText(image,toStr(centers.size()),cv::Point2f(5,75),cv::FONT_HERSHEY_SIMPLEX,1,CV_RGB(0,0,0),2);

    //Display image
    cv::imshow("Detector",image);
    key = cv::waitKey(5);

    //Compute means
    execTimeMean = (nbIter*execTimeMean + execTime)/(nbIter+1);
    FPSmean = (nbIter*FPSmean + fps)/(nbIter+1);
    nbIter++;

  }while(key<0 && isVideo);

  std::cout << "Mean FPS: " << FPSmean << "Hz" << std::endl;
  std::cout << "Mean exec: " << execTimeMean << "ms" << std::endl;

  if(!isVideo)
    cv::waitKey(0);

  return 0;
}

void help(std::string exec)
{
  std::cout << "Usage:\n\t" << exec << " [option] [source] " << std::endl;
  std::cout << "Source:\n\tCamera index (default: 0)\n\tVideo file\n\tImage file" << std::endl;
  std::cout << "Options:"<<std::endl;
  std::cout << "\t-t <value>\t\t\tBinary threshold in [0:255] or \"Auto\" (default: 125)" << std::endl;
#ifndef DISABLE_XML 
  std::cout << "\t-c <filename.xml>\t\tRead grid config" << std::endl;
  std::cout << "\t-g <filename.xml>\tGenerate default grid config" << std::endl;
#endif
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

#ifndef DISABLE_XML 
void createDefaultGridConfig(std::string filepath)
{
  cv::FileStorage file(filepath, cv::FileStorage::WRITE | cv::FileStorage::FORMAT_XML);
  file << "targetType" << "3B";
  file << "gridSize" << cv::Size(4,5);
  file << "values" << std::vector<int>();
  file << "headerBits" << 8;
  file << "headerValue" << 180;
  file << "messageBits" << 8;
  file << "useParityBit" << false;
}
#endif