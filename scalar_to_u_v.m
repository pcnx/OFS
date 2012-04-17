function  scalar_to_u_v
close all



load('scalar.dat');
scalar = scalar(1:10,1:10);

[U,V] = gradient(scalar);
% figure(1)
% mesh(U)
% figure(2)
% mesh(V)
% figure(3)
mesh(scalar)
figure(4)
quiver(U,V);
% 
% W = zeros(size(U));
% min(min(U))
% min(min(V))
% 
% tol = 1e-2;
% for i=1:length(U)
%     for j=1:length(U)
%         if abs(U(i,j))<tol && abs(V(i,j))<tol
%             W(i,j) = 1;
%         end
%     end
% end
% 
% surface(W)
% 
%     

% 
% POT = zeros(100,100);
% 
% for i=1:100
%     for j=1:100
%         POT(i,j) = i*0.01*1 + 1 * log((i*0.01-0.5)^2+(j*0.01-0.5)^2);
%     end
% end
% figure(4)
% mesh(POT)
% [W,Z] = gradient(POT);
% figure(5)
% mesh(W)
% figure(6)
% mesh(Z);






end