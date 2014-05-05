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

int main(int argc, char* argv[])
{
  TargetDetector detector;
  
  cv::Mat I1 = cv::imread(argv[1]);
  cv::Mat I2 = cv::imread(argv[2]);
  cv::Mat I3 = cv::imread(argv[3]);
  
  cv::Mat input(std::max(I1.cols+I2.cols,I3.cols),I1.rows+I3.rows, CV_MAKETYPE(I1.depth(), I1.channels()), CV_RGB(255,255,255) );
  
  cv::Mat mask1 = input( cv::Rect(0,0,I1.cols,I1.rows) );
  I1.copyTo(mask1);
  cv::Mat mask2 = input( cv::Rect(I1.cols,0,I2.cols,I2.rows) );
  I2.copyTo(mask2);
  cv::Mat mask3 = input( cv::Rect(0,I1.rows,I3.cols,I3.rows) );
  I3.copyTo(mask3);
    
  std::vector<Target> searches;
  searches.push_back( Target(Target::OneBlob,8,8,false,180) );
  searches.push_back( Target(Target::ThreeBlobs,8,8,false,180) );
  searches.push_back( Target(Target::TwoRings,8,8,false,180) );
  
  std::vector<Target> targets = detector.track(input,searches);
  detector.drawTargets(input,targets);
  
  cv::imwrite("Result.png",input);
  
  int nbTargets = targets.size();
  std::cout << nbTargets << std::endl;
  
#ifdef DEBUG_MODE
  cv::imshow("Result",input);
  cv::waitKey(0);
#endif
  
  return nbTargets!=3;
}
