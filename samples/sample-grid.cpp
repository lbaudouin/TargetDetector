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

/**
* @author Léo Baudouin\n
* @em baudouin.leo@gmail.com
* 
* Sample showing 4 methods to find grid of targets
*/

int main(int argc, char* argv[])
{
  //Load image
  cv::Mat image = cv::imread(argv[1]);
  
  //Create TargetDetector
  int threshold = 125;
  TargetDetector targetDetector(threshold);

  //Create search target
  Target search(Target::OneBlob,8,8,false,180);

  //Set grid size
  cv::Size size(4,4);
  
  //Create vector of centers
  std::vector<cv::Point2f> centers;
  
  //Find grid using inscreasing value order
  bool found1 = targetDetector.findTargetsGrid(image,size,centers,search);
  
  //Find grid using inscreasing value order (static method)
  bool found2 = TargetDetector::findTargetsGrid(image,size,centers,search,threshold);
  
  //Create values vector
  static const int v[] = {6,12,10,24,30,20,18,48,54,60,40,46,36,34,96,102};
  std::vector<int> values (v, v + sizeof(v) / sizeof(v[0]) );
  
  //Find grid using input values
  bool found3 = targetDetector.findTargetsGrid(image,size,centers,search,values);

  //Find grid using input values (static method)
  bool found4 = TargetDetector::findTargetsGrid(image,size,centers,search,threshold,values);
    
  //Draw grid on image
  TargetDetector::drawTargetGrid(image,size,centers,found4);
  
  return !(found1 && found2 && found3 && found4);
}