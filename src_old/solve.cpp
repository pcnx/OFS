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

if(!gaussseidel(data,data->s1))	{ return false; }
//if(!jacobi(data, data->s1))	{ return false; }
//if(!thomas(data,data->s1))	{ return false; }

return true;
}

//------------------------------------------------------
bool gaussseidel(sData* data, double** s)
{
int curIter=0;

while(curIter<data->maxIter) {
            std::cout << "\r\tGauss-Seidel: Iteration " << ++curIter ;
            float error = 0;
            float temp;


                for(int i = 1;i<data->dimI-1;i++){
                	for (int j=1;j<data->dimJ-1;j++){
                		temp =  (s[i][j+1]+s[i][j-1]+s[i+1][j]+s[i-1][j])/4;
                		error += abs(temp-s[i][j]);
                		s[i][j] = temp;

                	}
                		

                }




            if (error< data->residuum) return true;



} 

return true;
}

float abs(float x)
{
return x>0 ? x : -x;
}

//------------------------------------------------------
bool jacobi(sData* data, double** s)
{
int curIter=0;
double** tmp_ptr;

double** s_old = s;
double** s_new = allocGrid1Mem(data,MAXDOUBLE);
double error;
double temp;
for (int i= 0; i<data->dimI;i++){
    s_new[i][0] = s_old[i][0];
    s_new[i][data->dimJ-1]=s_old[i][data->dimJ-1];
}
for (int i=0;i<data->dimJ;i++){
    s_new[0][i] = s_old[0][i];
    s_new[data->dimI-1][i]=s_old[data->dimI-1][i];
}

while(curIter<data->maxIter ) {
        std::cout << "\r\tJakobi: Iteration " << ++curIter;
        error = 0;

        for (int i=1;i< data->dimI-1;i++){
            for (int j=1;j<data->dimJ-1;j++){
                temp = (s_old[i-1][j] +s_old[i+1][j] + s_old[i][j-1]+s_old[i][j+1])/4.0;
                error += abs( temp- s_old[i][j]);
                s_new[i][j] = temp;
            }
        }

	tmp_ptr = s_old;
	s_old = s_new;
	s_new  = tmp_ptr;

   // if (error < data->residuum) return true;

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
	s_new = tmp_ptr;
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
