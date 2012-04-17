/*
 * transformation.h
 *
 *  Created on: Apr 17, 2012
 *      Author: thomas
 */

#ifndef TRANSFORMATION_H_
#define TRANSFORMATION_H_


#include "data.h"

double xiof(double x );
double etaof(double x, double y);
double xof( double xi);
double yof( double eta, double xi);
double deta(sData* data, int i, int j, bool xdirection);
double dxi(sData* data, int i, int j, bool xdirection);
double ddeta(sData* data, int i, int j, bool xdirection);
double ddxi(sData* data, int i, int j, bool xdirection);


#endif /* TRANSFORMATION_H_ */
