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
  double dx = data->finiteDiffDx;
  double dy = data->finiteDiffDy;

  if (xdirection){
      return (xiof(data->x[i][j]+dx)-xiof(data->x[i][j]-dx))/(2*dx);
  }
  else{
      return (xiof(data->y[i][j]+dy)-xiof(data->y[i][j]-dy))/(2*dy);
  }
}

double deta(sData* data, int i, int j, bool xdirection){
  double dx = data->finiteDiffDx;
  double dy = data->finiteDiffDy;

  if (xdirection){
      return (etaof(data->x[i][j]+dx,data->y[i][j])-etaof(data->x[i][j]-dx,data->y[i][j]))/(2*dx);
  }
  else{
      return (etaof(data->x[i][j],data->y[i][j]+dy)-etaof(data->x[i][j],data->y[i][j]-dy))/(2*dy);
  }
}



double ddxi(sData* data, int i, int j, bool xdirection){
  double dx = data->finiteDiffDx;
  double dy = data->finiteDiffDy;

  if (xdirection){
      return (xiof(data->x[i][j]+dx)-2*xiof(data->x[i][j])+xiof(data->x[i][j]-dx))/(dx*dx);
  }
  else{
      return (xiof(data->y[i][j]+dy)-2*xiof(data->y[i][j])+xiof(data->y[i][j]-dy))/(dy*dy);
  }
}

double ddeta(sData* data, int i, int j, bool xdirection){
  double dx = data->finiteDiffDx;
  double dy = data->finiteDiffDy;

  if (xdirection){
      return (etaof(data->x[i][j]+dx,data->y[i][j])
          -2*etaof(data->x[i][j],data->y[i][j])
      +etaof(data->x[i][j]-dx,data->y[i][j]))/(dx*dx);
  }
  else{
      return (etaof(data->x[i][j],data->y[i][j]+dy)
          -2*etaof(data->x[i][j],data->y[i][j])
      +etaof(data->x[i][j],data->y[i][j]-dy))/(dy*dy);
  }
}




void dxi(sData* data, double** dxidx, double** dxidy){
  double dx = data->finiteDiffDx;
  double dy = data->finiteDiffDy;
  double twicedx = 2*dx;
  double twicedy = 2*dy;
  double x;
  double y;
  for (int i=1; i<data->dimI-1;i++){
      for (int j =1; j<data->dimJ-1;j++){
          x = data->x[i][j];
          y = data->y[i][j];
          dxidx[i][j] = (xiof(x+dx)-xiof(x-dx)) /twicedx;
          dxidy[i][j] = (xiof(y+dy)-xiof(y-dy)) /twicedy;
      }
  }

}

void deta(sData* data, double** detadx, double** detady){
  double dx = data->finiteDiffDx;
  double dy = data->finiteDiffDy;
  double twicedx = 2*dx;
  double twicedy = 2*dy;
  double x;
  double y;
  for (int i=1; i<data->dimI-1;i++){
      for (int j =1; j<data->dimJ-1;j++){
          x = data->x[i][j];
          y = data->y[i][j];
          detadx[i][j] = (etaof(x+dx,y)-etaof(x-dx,y)) /twicedx;
          detady[i][j] = (etaof(x,y+dy)-etaof(x,y-dy)) /twicedy;
      }
  }

}

void ddxi(sData* data, double** ddxidx, double** ddxidy){
  double dx = data->finiteDiffDx;
  double dy = data->finiteDiffDy;
  double dx2 = dx*dx;
  double dy2 = dy*dy;
  double x;
  double y;
  for (int i=1; i<data->dimI-1;i++){
      for (int j =1; j<data->dimJ-1;j++){
          x = data->x[i][j];
          y = data->y[i][j];
          ddxidx[i][j] = (xiof(x+dx)-2*xiof(x)+xiof(x-dx)) /dx2;
          ddxidy[i][j] = (xiof(y+dy)-2*xiof(y)+xiof(y-dy)) /dy2;
      }
  }

}

void ddeta(sData* data, double** ddetadx, double** ddetady){
  double dx = data->finiteDiffDx;
  double dy = data->finiteDiffDy;
  double dx2 = dx*dx;
  double dy2 = dy*dy;
  double x;
  double y;
  for (int i=1; i<data->dimI-1;i++){
      for (int j =1; j<data->dimJ-1;j++){
          x = data->x[i][j];
          y = data->y[i][j];
          ddetadx[i][j] = (etaof(x+dx,y)-2*etaof(x,y)+etaof(x-dx,y)) /dx2;
          ddetady[i][j] = (etaof(x,y+dy)-2*etaof(x,y)+etaof(x,y-dy)) /dy2;
      }
  }

}
