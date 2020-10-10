#include <iostream>
#include <string>

#include "Intersect.h"

int main() {
  std::string polygonAlpha = "0 0,0 1,2 0,0 0";
  std::string polygonBravo = "1 0,1 1,3 0,1 0";
  std::string polygonCharlie = "2 2,2 3,4 2,2 2";

  if (Intersect::overlap(polygonAlpha, polygonBravo)) {
    //std::cout << "OK" << std::endl;
  }

  // Alpha and Charlie are not overlapping, so if our function says there exists an overlap
  // then this is clearly not OK
  //if (Intersect::overlap(polygonAlpha, polygonCharlie)) {
   // std::cout << "NOK" << std::endl;
  //}

  return 0;
}
