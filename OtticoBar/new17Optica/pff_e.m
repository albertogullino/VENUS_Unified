
 figure;
 pograp=[200   50   380   880];
 set(gcf,'Position',pograp)
 subplot(3,1,1)
            surf(X,Y,abs(Efz)),
            shading('interp'), view(0,90),
            colorbar
            title(' Ef_z ')
            axis square, axis equal, grid,
            axis([-1 1 -1 1]*axli/2),
 subplot(3,1,2)
            surf(X,Y,abs(Efx)),
            shading('interp'), view(0,90),
            colorbar
            title(' Ef_x ')
            axis square, axis equal, grid,
            axis([-1 1 -1 1]*axli/2),
 subplot(3,1,3)
            surf(X,Y,abs(Efy)),
            shading('interp'), view(0,90),
            colorbar
            title(' Ef_y ')
            axis square, axis equal, grid,
            axis([-1 1 -1 1]*axli/2),
