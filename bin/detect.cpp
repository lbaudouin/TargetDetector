#include <iostream>
#include "../include/targetdetector.h"

#include "timer.h"

/**
* @author Léo Baudouin\n
* @em baudouin.leo@gmail.com
*/

std::string toStr(double d)
{
    std::stringstream ss;
    ss << d;
    return ss.str();
}

int main(int argc, char* argv[])
{
  //cv::VideoCapture capture("test.avi");
  cv::VideoCapture capture;
  
  if(argc>1){
    //Try to read camera index
    int cameraIndex;
    std::istringstream iss(argv[1]);
    iss >> std::ws >> cameraIndex >> std::ws;
    if(iss.eof()){
      //Open camera
      capture.open(cameraIndex); 
    }else{
      //Open file
      capture.open(argv[1]);          
    }
  }else{
    //Default
    capture.open(1);
  }
  
  if(!capture.isOpened()){
    std::cerr << "Failed to open video capture" << std::endl;
    return 1;
  }
  
  TargetDetector targetDetector;
  
  std::vector<Target> searches;
  searches.push_back( Target(Target::OneBlob,8,8,false,180) );
  searches.push_back( Target(Target::ThreeBlobs,8,8,false,180) );
  searches.push_back( Target(Target::ThreeBlobs,8,8,false,220) );
  searches.push_back( Target(Target::TwoRings,8,8,false,180) );
  
  //Run autoThreshold
  //if(!targetDetector.autoThreshold(capture,searches[0],10,250))
    targetDetector.setThreshold(125);
    
  
  Timer FPStimer,execTimer;
  char key;
  cv::Mat image;
  
  int nb = 0;
  double execTimeMean = 0;
  double FPSmean = 0;
  
  do{
    //Grab image
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
    double fps = 1.0/FPStimer.s_elapsed();

    FPStimer.restart();
    cv::putText(image,toStr(fps),cv::Point2f(5,25),cv::FONT_HERSHEY_SIMPLEX,1,CV_RGB(0,0,0),2);
    cv::putText(image,toStr(execTime),cv::Point2f(5,50),cv::FONT_HERSHEY_SIMPLEX,1,CV_RGB(0,0,0),2);
    cv::putText(image,toStr(targets.size()),cv::Point2f(5,75),cv::FONT_HERSHEY_SIMPLEX,1,CV_RGB(0,0,0),2);
    
    //Display image
    cv::imshow("Detector",image);
    key = cv::waitKey(5);
    
    execTimeMean = (nb*execTimeMean + execTime)/(nb+1);
    FPSmean = (nb*FPSmean + fps)/(nb+1);
    
    nb++;
  }while(key<0);
  
  std::cout << "Mean FPS: " << FPSmean << "Hz" << std::endl;
  std::cout << "Mean exec: " << execTimeMean << "ms" << std::endl;
  
  return 0;
}

