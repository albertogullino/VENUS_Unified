acli1=aax*2;
figure
 pograp=[300   0   380   680];
 set(gcf,'Position',pograp)
 cE=abs(max([max(max(Exdu)) max(max(Eydu)) max(max(Ezdu))]));
 cEx=max(max(Exdu));
 cEy=max(max(Eydu));
 cEz=max(max(Ezdu));
 fEx=abs(cEx/cE);
 fEy=abs(cEy/cE);
 fEz=abs(cEz/cE);

 subplot(3,1,1)
            surf(XP,YP,real(Exdu/cEx*fEx)),
            shading('interp'), view(0,90),
            colorbar
            title(' E_x ')
            axis square, axis equal, grid,
            axis([-1 1 -1 1]*acli1/2),
            axis off
 subplot(3,1,2)
            surf(XP,YP,real(Eydu/cEy*fEy)),
            shading('interp'), view(0,90),
            colorbar
            title(' E_y ')
            axis square, axis equal, grid,
            axis off
            axis([-1 1 -1 1]*acli1/2),
 subplot(3,1,3)
            surf(XP,YP,real(Ezdu/cEz*fEz)),
            shading('interp'), view(0,90),
            colorbar
            title(' E_z ')
            axis square, axis equal, grid,
            axis off
            axis([-1 1 -1 1]*acli1/2),

