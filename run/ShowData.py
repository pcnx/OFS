import os
from pylab import *


# set path
DATA_PATH ='../'

# load grid
MESHX=mlab.load(DATA_PATH+'phyGrid.meshX');
MESHY=mlab.load(DATA_PATH+'phyGrid.meshY');

subplot(1,1,1);
DATA_ZERO = zeros(MESHX.shape,float)
contourf(MESHX,MESHY,DATA_ZERO,colors = '0.75');
hold(True); 

# show potential lines
if os.path.isfile(DATA_PATH+'scalar.dat'):
    DATA_POT=mlab.load(DATA_PATH+'scalar.dat');

    subplot(1,1,1);

    contour(MESHX,MESHY,DATA_POT,20);
    title(' potential - isolines ');
    #colorbar(cpot)
    axis('equal')
    hold(True); 

# show stream lines
#if os.path.isfile(DATA_PATH+'scalar.dat'):
#    DATA_STR=numpy.load(DATA_PATH+'scalar.dat');
#    subplot(1,1,1);
#    cstr = contour(MESHX,MESHY,DATA_STR,20);
#    title(' stream - isolines ');
#    colorbar(cstr)
#    axis('equal')
#    hold(True);

title(' potential/stream - isolines ');
#grid(True)
show()
