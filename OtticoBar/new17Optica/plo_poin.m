
Poi=Pvectx+Pvecty+Pvectz;
Nor=max(max(Poi));

 figure;
            surf(X,Y,Poi/Nor),
            shading('interp'), view(0,90),
            colorbar
            title(' Ef_z ')
            axis square, axis equal, grid,
            axis([-1 1 -1 1]*axli/2),

 figure;
 pograp=[800   50   380   880];
 set(gcf,'Position',pograp)
 subplot(3,1,1)
            surf(X,Y,Pvectx/Nor),
            shading('interp'), view(0,90),
            colorbar
            title(' Pointing_x ')
            axis square, axis equal, grid,
            axis([-1 1 -1 1]*axli/2),
 subplot(3,1,2)
            surf(X,Y,Pvecty/Nor),
            shading('interp'), view(0,90),
            colorbar
            title(' Pointing_y ')
            axis square, axis equal, grid,
            axis([-1 1 -1 1]*axli/2),
 subplot(3,1,3)
            surf(X,Y,Pvectz/Nor),
            shading('interp'), view(0,90),
            colorbar
            title(' Pointing_z ')
            axis square, axis equal, grid,
            axis([-1 1 -1 1]*axli/2),
            print -dmeta
%            figure
%            surf(X,Y,fia), colorbar
%            shading('interp'), view(0,90),
