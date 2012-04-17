/*
 * transformation.cpp
 *
 *  Created on: Apr 17, 2012
 *      Author: thomas
 */

#include <math.h>
#include <iostream>
#include "transformation.h"
#include "data.h"



double xiof(double x){
  return x;

}
double etaof(double x, double y){
  return y/(x+1);

}
double xof(double xi){
  return xi;
}
double yof(double xi, double eta){
  return eta*(1+xi);

}


double dxi(sData* data, int i, int j, bool xdirection){

  if (xdirection){
      return (xiof(data->xi[i+1][j])-xiof(data->xi[i-1][j]))/(2*data->deltaXi);
  }
  else{
      return (xiof(data->eta[i][j+1])-xiof(data->eta[i][j-1]))/(2*data->deltaEta);
  }
}

double deta(sData* data, int i, int j, bool xdirection){

  if (xdirection){
      return (etaof(data->xi[i+1][j],data->eta[i][j])-etaof(data->xi[i-1][j],data->eta[i][j]))/(2*data->deltaXi);
  }
  else{
      return (etaof(data->xi[i][j],data->eta[i][j+1])-etaof(data->xi[i][j],data->eta[i][j-1]))/(2*data->deltaEta);
  }
}



double ddxi(sData* data, int i, int j, bool xdirection){

  if (xdirection){
      return (xiof(data->xi[i+1][j])-2*xiof(data->xi[i][j])+xiof(data->xi[i-1][j]))/(data->deltaXi*data->deltaXi);
  }
  else{
      return (xiof(data->eta[i][j+1])-2*xiof(data->eta[i][j])+xiof(data->eta[i][j-1]))/(data->deltaEta*data->deltaEta);
  }
}

double ddeta(sData* data, int i, int j, bool xdirection){

  if (xdirection){
      return (etaof(data->xi[i+1][j],data->eta[i][j])-
          2*etaof(data->xi[i][j],data->eta[i][j])+etaof(data->xi[i-1][j],data->eta[i][j]))/(data->deltaXi*data->deltaXi);
  }
  else{
      return (etaof(data->xi[i][j],data->eta[i][j+1])-2*etaof(data->xi[i][j],data->eta[i][j])
      +etaof(data->xi[i][j],data->eta[i][j-1]))/(data->deltaEta*data->deltaEta);
  }
}
