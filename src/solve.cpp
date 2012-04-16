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

#include "solve.h"
#include "data.h"
#include "setup.h"

//------------------------------------------------------
bool solve(sData* data)
{
        std::cout << "\nSolve:\n-------\n";

	if(!gaussseidelMorphed(data,data->s1))	{ return false; }
	//if(!jacobi(data, data->s1))	{ return false; }
	//if(!thomas(data,data->s1))	{ return false; }

	return true;
}

//------------------------------------------------------
bool gaussseidel(sData* data, double** s)
{
	int curIter=0;
	float maxDiff = 0;
	float diff = s[0][0];
	float tmp;

	while(curIter<data->maxIter) {
                std::cout << "\r\tGauss-Seidel: Iteration " << ++curIter;
				maxDiff =0;
				for(int i = 1; i < data->dimI-1; i++)
				{
					for(int j = 1 ; j < data->dimJ-1; j++)
					{
					//Iterate over all values except border values
					tmp =(s[i-1][j]+s[i+1][j]+s[i][j+1]+s[i][j-1])/4;
					diff = diff-tmp;
					
					maxDiff += myAbs(diff);
					s[i][j] = tmp;
					
					}
				
				}

				
				if(maxDiff < data->residuum)
						break;
			
	}
	double ** tmp5 = new double*[data->dimI] ;
	for(int i=0; i < data->dimI; i++)
	{
	tmp5[i] = new double[data->dimJ];
	}
	for(int i = 0; i < data->dimI; i++)
	{
		tmp5[i][0] = s[i][0];
		tmp5[i][data->dimJ-1] = s[i][data->dimJ-1];
	}
	for(int j = 0 ; j < data->dimJ; j++)
	{
		tmp5[0][j] = s[0][j];
		tmp5[data->dimI-1][j] = s[data->dimJ-1][j];
	}

					for(int i = 1; i < data->dimI-1; i++)
					{
					for(int j = 0 ; j < data->dimJ; j++)
					{
						tmp5[i][j] = (s[i+1][j] - s[i-1][j])/(2);
						
					}
					}
					data->s1 = tmp5;
	return true;
}

bool gaussseidelMorphed(sData* data, double** s)
{
	int curIter=0;
	float maxDiff = 0;
	float diff = s[0][0];
	float tmp;
	while(curIter<data->maxIter) {
                std::cout << "\r\tGauss-Seidel: Iteration " << ++curIter;
				maxDiff =0;
				for(int i = 1; i < data->dimI-1; i++)
				{
					for(int j = 1 ; j < data->dimJ-1; j++)
					{
					//Iterate over all values except border values
						tmp = 1/(2*(a1+a2))*(a1 *(s[i+1][j]+s[i-1][j]) + a2* (s[i][j+1]+s[i][j-1]) +a3/4* (s[i+1][j+1]-s[i-1][j+1]-s[i+1][j-1]+s[i-1][j-1]) + a4/2 * (s[i+1][j]-s[i-1][j]) + a5/2 * (s[i][j+1]+s[i][j-1]));
					diff = diff-tmp;
					
					maxDiff += myAbs(diff);
					s[i][j] = tmp;
					
					}
				
				}

				
				if(maxDiff < data->residuum)
						break;
			
	}
	return true;
}

//------------------------------------------------------
bool jacobi(sData* data, double** s)
{
	int curIter=0;
	double** tmp_ptr;

	double** s_old = s;
	double** s_new = allocGrid1Mem(data,MAXDOUBLE);

	while(curIter<data->maxIter) {
            std::cout << "\r\tJakobi: Iteration " << ++curIter;

            // 1. loop over grid points
                  // TODO

            // 2. possible overrelaxation
                  // TODO

            // 3. new data-set becomes old data-set
		tmp_ptr = s_old;
		s_old = s_new;
		s_new = s_old;

            // 3. break on min diff
                  // TODO
	}

	// sync data-fields if nessessary
	if(s!=s_old) {
		for(int i=0; i<data->dimI; i++) {
			for(int j=0; j<data->dimJ; j++) {
				s_new[i][j] = s_old[i][j];
			}
		}
		tmp_ptr = s_old;
		s_old = s_new;
		s_new = s_old;
	}

	// free temp. memory (make sure you are not freeing memory s points to!)
	freeGrid1Mem(data,s_new);

	return true;
}

//------------------------------------------------------
bool thomas(sData* data, double** s)
{
	// optional
	return true;
}

float myAbs(float x)
{

    return (x > 0) ? x : -x;
}
