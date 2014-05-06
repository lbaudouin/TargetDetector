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
* Sample showing 4 methods to find targets
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

  //Find targets
  std::vector<Target> targets1 = targetDetector.track(image,search);
  
  //Find tragets (static method)
  std::vector<Target> targets2 = TargetDetector::track(image,search,threshold);
  
  //Set vector of search targets
  std::vector<Target> searches;
  searches.push_back( Target(Target::OneBlob,8,8,false,180) );
  searches.push_back( Target(Target::ThreeBlobs,8,8,false,180) );
  
  //Find targets type
  std::vector<Target> targets3 = targetDetector.track(image,searches);
  
  //Find targets type (static method)
  std::vector<Target> targets4 = TargetDetector::track(image,search,threshold);
  
  //Draw targets on image
  TargetDetector::drawTargets(image,targets1);
  
  return 0;
}