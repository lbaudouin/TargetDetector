#include "targetdetector.h"

int main(int argc, char* argv[])
{
  TargetDetector detector;
  
  cv::Mat I = cv::imread(argv[1]);
  
  Target search(Target::OneBlob,8,8,false,180);
  
  std::vector<Target> targets = detector.track(I,search);
  detector.drawTargets(I,targets);
  
  cv::imwrite("OneBlob.png",I);
  
  int nbTargets = targets.size();
  std::cout << nbTargets << std::endl;
  
#if TARGET_DEBUG
  cv::imshow("Result",I);
  cv::waitKey(0);
#endif
  
  return nbTargets!=1;
}
