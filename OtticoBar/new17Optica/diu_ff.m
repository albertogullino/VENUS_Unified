close all
 ploff_nw
 cno=sqrt(mfxy);

 figure
 surf(X,Y,real(Ef/max(max(Ef)))),
 shading('interp'), view(0,90),
 axis square, axis equal, grid,
 title(' Etot')
 colorbar
 axis([-1 1 -1 1]*axli/2),
 pausak

 figure
 maE=max(max(Efx));
 surf(X,Y,real(Efx/maE)*abs(maE)/cno),
 shading('interp'), view(0,90),
 axis square, axis equal, grid,
 title(' Ex')
 colorbar
 axis([-1 1 -1 1]*axli/2),
 pausak

 figure
% surf(X,Y,real(Efy)/cno),
 maE=max(max(Efy));
 surf(X,Y,real(Efy/maE)*abs(maE)/cno),
 shading('interp'), view(0,90),
 axis square, axis equal, grid,
 title(' Ey')
 colorbar
 axis([-1 1 -1 1]*axli/2),
 pausak

 figure
% surf(X,Y,abs(Efz)/cno),
 maE=max(max(Efz));
 surf(X,Y,real(Efz/maE)*abs(maE)/cno),
 shading('interp'), view(0,90),
 axis square, axis equal, grid,
 title(' Ez')
 colorbar
 axis([-1 1 -1 1]*axli/2),
 pausak

 teR=asin(KK(1:length(KK))*rr);

 nrff=55;
 teR=linspace(0,max(teR),nrff)';
 X1=teR*180/pi*cos(fian);
 Y1=teR*180/pi*sin(fian);
 teRg=repmat(teR,numodi,1);


 ploff_ul
 cno=sqrt(mfxy);

 figure
 surf(X1,Y1,real(Ef)),
 shading('interp'), view(0,90),
 axis square, axis equal, grid,
 title(' Etot')
 axis([-1 1 -1 1]*axli/2),
 colorbar
 pausak

 figure
 maE=max(max(Efx));
 surf(X1,Y1,real(Efx/maE)*abs(maE)/cno),
 shading('interp'), view(0,90),
 axis square, axis equal, grid,
 title(' Ex')
 axis([-1 1 -1 1]*axli/2),
 colorbar
 pausak

 figure
 maE=max(max(Efy));
 surf(X1,Y1,real(Efy/maE)*abs(maE)/cno),
 shading('interp'), view(0,90),
 axis square, axis equal, grid,
 title(' Ey')
 colorbar
 axis([-1 1 -1 1]*axli/2),
 pausak

 figure
 maE=max(max(Efz));
 surf(X1,Y1,real(Efz/maE)*abs(maE)/cno),
 shading('interp'), view(0,90),
 axis square, axis equal, grid,
 title(' Ez')
 colorbar
 axis([-1 1 -1 1]*axli/2),
 pausak
