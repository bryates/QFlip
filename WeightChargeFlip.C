#include <math.h>

Double_t pol1(Double_t x); 
Double_t pol2(Double_t x);
Double_t weight(Double_t x);

Double_t WeightChargeFlip(Double_t pt1 = 20, Double_t pt2 = 50) {
  Double_t W = (weight(pt1)/(1-weight(pt1)) + weight(pt2)/(1-weight(pt2)));
  //if(!isinf(W))
    return W;
  //else
    //return 0;
}


Double_t weight(Double_t pt) {
  if (1/pt > 0.01 && 1/pt < 0.015)
    return pol1(1/pt);
  else if(1/pt >= 0.015)
    return pol2(1/pt);
    return 0;
}

Double_t pol1(Double_t x) {
  return -0.279579*x + 0.00418904;
}

Double_t pol2(Double_t x) {
  return 0.00243305*x + -4.98711e-05;
}
/*
Double_t pol1(Double_t x) {
  return 41.3817*pow(x,2) - 1.38554*x + 0.0115067;
}

Double_t pol2(Double_t x) {
  return 0.246988*pow(x,2) - 0.0112697*x + 0.000138197;
}
*/
