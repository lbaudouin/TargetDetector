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

#include "targetdetector.h"

int main(int argc, char* argv[])
{
  TargetDetector detector;
  
  cv::Mat I1 = cv::imread(argv[1]);
  cv::Mat I2 = cv::imread(argv[2]);
  
  Target search(Target::OneBlob,8,8,false,180);
  
  cv::Size size5x5(5,5), size4x4(4,4);
  std::vector<cv::Point2f> centers;
  
  //Test ordered grid
  bool okOrdered = detector.findTargetsGrid(I1,size5x5,centers,search);
  cv::drawChessboardCorners(I1,size5x5,centers,okOrdered);
  cv::imwrite("OrderedGrid.png",I1);
  
#ifdef DEBUG_MODE
  cv::imshow("Result",I1);
  cv::waitKey(0);
#endif
  
  int val[16] = {5, 50, 55, 15, 10, 95, 125, 255, 25, 0, 65, 100, 105, 85, 30, 235};
  
  std::vector<int> values;
  for(int i=0;i<16;i++)
    values.push_back(val[i]);
  
  //Test unordered grid
  bool okUnordered = detector.findTargetsGrid(I2,size4x4,centers,search,values);
  cv::drawChessboardCorners(I2,size4x4,centers,okUnordered);
  cv::imwrite("UnorderedGrid.png",I2);
  
#ifdef DEBUG_MODE
  cv::imshow("Result",I2);
  cv::waitKey(0);
#endif
  
  if(!okUnordered && !okOrdered)
    return true;
  else
    return false;
}
