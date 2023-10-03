close all
 plo_ff
 cno=sqrt(mfxy);

% figure
% surf(X,Y,abs(Ef/max(max(Ef)))),
% shading('interp'), view(0,90),
% axis square, axis equal, grid,
% title(' Etot')
% colorbar
% axis([-1 1 -1 1]*axli/2),
% pausak

 figure
 maE=max(max(Efx));
 surf(X,Y,abs(Efx/maE)*abs(maE)/cno),
 shading('interp'), view(0,90),
 axis square, axis equal, grid,
 title(' Ex')
 colorbar
 axis([-1 1 -1 1]*axli/2),
 pausak

 figure
% surf(X,Y,abs(Efy)/cno),
 maE=max(max(Efy));
 surf(X,Y,abs(Efy/maE)*abs(maE)/cno),
 shading('interp'), view(0,90),
 axis square, axis equal, grid,
 title(' Ey')
 colorbar
 axis([-1 1 -1 1]*axli/2),
 pausak

 figure
% surf(X,Y,abs(Efz)/cno),
 maE=max(max(Efz));
 surf(X,Y,abs(Efz/maE)*abs(maE)/cno),
 shading('interp'), view(0,90),
 axis square, axis equal, grid,
 title(' Ez')
 colorbar
 axis([-1 1 -1 1]*axli/2),
 pausak

 nrff=55;
 teR=linspace(0,max(teR),nrff)';
 X1=teR*180/pi*cos(fian);
 Y1=teR*180/pi*sin(fian);
 Xpr=sin([teR])*cos(fian);
 Ypr=sin([teR])*sin(fian);
 Zpr=cos([teR])*ones(size(fian));

intgu_ff
npkf=nrff;
 po_pro
% ploff_nw
 cno=sqrt(mfxy);
% X1=sin([0; teR])*cos(fian);
% Y1=sin([0; teR])*sin(fian);
% X1=[0; teR]*180/pi*cos(fian);
% Y1=[0; teR]*180/pi*sin(fian);

% figure
% surf(X1,Y1,abs(Ef)),
% shading('interp'), view(0,90),
% axis square, axis equal, grid,
% title(' Etot')
% axis([-1 1 -1 1]*axli/2),
% colorbar
% pausak

 figure
 surf(X1,Y1,abs(Efx)/cno),
 shading('interp'), view(0,90),
 axis square, axis equal, grid,
 title(' Ex')
 axis([-1 1 -1 1]*axli/2),
 colorbar
 pausak

 figure
 surf(X1,Y1,abs(Efy)/cno),
 shading('interp'), view(0,90),
 axis square, axis equal, grid,
 title(' Ey')
 colorbar
 axis([-1 1 -1 1]*axli/2),
 pausak

 figure
 surf(X1,Y1,abs(Efz)/cno),
 shading('interp'), view(0,90),
 axis square, axis equal, grid,
 title(' Ez')
 colorbar
 axis([-1 1 -1 1]*axli/2),
 pausak
