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
#include "stdlib.h"

int main(int argc, char* argv[])
{
  if(argc<2) return EXIT_FAILURE;
  
  TargetDetector detector(0);
  
  cv::Mat I = cv::imread(argv[1]);
  
  Target search(Target::OneBlob,8,8,false,180);
  
  bool ok = detector.autoThreshold(I,search);
    
  if(!ok){
    return true;
  }
  
  int threshold = detector.threshold();
      
  return threshold!=127;
}
