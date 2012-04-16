% Plotting grid

clear all
shading interp

% load grid
MESHX=load('phyGrid.meshX');
MESHY=load('phyGrid.meshY');

% show data
if exist('str_sol.dat')
    DATA_STR=load('str_sol.dat');
    iso=subplot(2,2,1:2);
        [C_str,h_str]=contour(MESHX,MESHY,DATA_STR,20);
        title(' isolines ');
        set(iso,'DataAspectRatio',[1 1 1]);
        hold on;
    str=subplot(2,2,3);
        surf(MESHX,MESHY,DATA_STR);
        title(' stream ');
        set(str,'DataAspectRatio',[1 1 1]);   
    view(2);
    save 'str_cont.dat' -ASCII -DOUBLE C_str
end
if exist('pot_sol.dat')
    DATA_POT=load('pot_sol.dat');
    iso=subplot(2,2,1:2);
        [C_pot,h_pot]=contour(MESHX,MESHY,DATA_POT,20);
        title(' isolines ');
        set(iso,'DataAspectRatio',[1 1 1]);        
        hold on;
    pot=subplot(2,2,4);
        surf(MESHX,MESHY,DATA_POT);
        title(' potential ');
        set(pot,'DataAspectRatio',[1 1 1]);
    view(2);
    save 'pot_cont.dat' -ASCII -DOUBLE C_pot
end


plottools('on');


