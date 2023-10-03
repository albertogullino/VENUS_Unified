
 figure;
 pograp=[600   50   380   880];
 set(gcf,'Position',pograp)
 subplot(3,1,2)
            surf(X,Y,abs(Hfx)),
            shading('interp'), view(0,90),
            colorbar
            title(' Ex ')
            axis square, axis equal, grid,
            axis([-1 1 -1 1]*axli/2),
 subplot(3,1,3)
            surf(X,Y,abs(Hfy)),
            shading('interp'), view(0,90),
            colorbar
            title(' Hy ')
            axis square, axis equal, grid,
            axis([-1 1 -1 1]*axli/2),
 subplot(3,1,1)
            surf(X,Y,abs(Hfz)),
            shading('interp'), view(0,90),
            colorbar
            title(' Hz ')
            axis square, axis equal, grid,
            axis([-1 1 -1 1]*axli/2),
