/***************************************************************************
 *   Copyright (C) 2006-2011 by  Institute of Combustion Technology        *
 *   jens.henrik.goebbert@itv.rwth-aachen.de                               *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#include <math.h>
#include <iostream>

#include "setup.h"
#include "data.h"

//------------------------------------------------------
bool setup(sData* data)
{
        std::cout << "\nSetup:\n-------\n";

	///////////////////
	// SETUP XY-GRID //
	///////////////////
	for(int i=0; i<data->dimI; i++) {
		for(int j=0; j<data->dimJ; j++) {
			data->x[i][j] = data->deltaX*i;
		   	data->y[i][j] = data->deltaY*j;
		}
	}

	/////////////////////////////////
	// SETUP INITIAL SCALAR VALUES //
	/////////////////////////////////
	// set whole field 's1' to start value
	for(int i=0; i<data->dimI; i++) {
		for(int j=0; j<data->dimJ; j++) {
			data->s1[i][j] = 12;

		}
	}

	//////////////////////////////////
	// SETUP BOUNDARY SCALAR VALUES //
	//////////////////////////////////
	float dx;
	float dy;
	float E2PI = 1.0; 
	float UINF = 5.0;
	for(int i=0; i<data->dimI; i++) {
	/*	data->s1[i][0] = 5;
		data->s1[i][data->dimJ-1] = 5; */
		dx= 1.1 - data->x[i][0];
		dy = 0.5;
		data->s1[i][0] = UINF*i*data->deltaX + E2PI * log(pow(dx,2)+pow(dy,2))/2;
		dx = 1.1 -data->x[i][data->dimJ-1];
		data->s1[i][data->dimJ-1] = UINF*i*data->deltaX + E2PI * log(pow(dx,2)+pow(dy,2))/2;
	}
	for(int j=0; j<data->dimJ; j++) {
	/*	data->s1[0][j] = 5;
		data->s1[data->dimI-1][j] = 5; */

		dy = 0.5 - data->y[0][j];
		data->s1[0][j] = UINF*0 + E2PI * log(pow(1.1,2)+pow(dy,2))/2;
		dy =  data->y[data->dimI-1][j]-0.5;

		data->s1[data->dimI-1][j] = UINF*1.0 + E2PI * log(pow(0.1,2)+pow(dy,2))/2;
	}



	return true;
}
