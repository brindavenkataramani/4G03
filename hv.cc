#define BIT_SET(a,b) ((a) |= (1U<<(b)))
#define BIT_CLEAR(a,b) ((a) &= ~(1U<<(b))) 
#define BIT_FLIP(a,b) ((a) ^= (1U<<(b)))
#define BIT_CHECK(a,b) ((bool)((a) & (1U<<(b)))) // Set on the condition f else clear
//bool f; // conditional flag
//unsigned int m; // the bit mask
//unsigned int w; // the word to modify: if (f) w |= m; else w &= ~m;
#define COND_BIT_SET(a,b,f) ((a) = ((a) & ~(1U<<(b))) | ((-(unsigned int)f) & (1U<<(b)))) 

#include <vector>
#include <cmath>
#include <iostream>

using namespace std;

void hv(vector<double>& y, const vector<double>& x,int L) {

//for (vector<double>::iterator it = y.begin() ; it != y.end(); ++it) // *it=0.0;
  
  bool b ;
  unsigned int k;
  
  for (unsigned int i=0;i<x.size();i++){ 
    
    if (abs(x[i])>2.2e-16) {
      int jm = L-1;
      double xov2=x[i]/2.0; 
      
      for (int j=0;j<L;j++){
        k=i; 
        COND_BIT_SET(k,jm,BIT_CHECK(i,j)); 
        COND_BIT_SET(k,j,BIT_CHECK(i,jm)); 
        y[k%L] += xov2;
        jm = j;
    
        }
  
      }
  
    }

  for (unsigned int i=0;i<x.size();i++) {
    y[i]=y[i]-((double) L)/2.0*x[i]/2.0;
  
  }
  
}
