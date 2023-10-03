
Poi=Pvectx+Pvecty+Pvectz;
Nor=max(max(Poi));
acli1=aax*1.7;
if iLP==1
 figure;
            surf(XP,YP,Poi/Nor),
            shading('interp'), view(0,90),
            colorbar
            title(' Ef_z ')
            axis square, axis equal, grid,
            axis([-1 1 -1 1]*acli1/2),
else

 figure;
 pograp=[300   0   680   680];
 set(gcf,'Position',pograp)

 subplot(3,3,1)
            surf(XP,YP,real(Pvectx/Nor)),
            shading('interp'), view(0,90),
            colorbar
            title(' Pointing_x ')
            axis square, axis equal, grid,
            axis off
            axis([-1 1 -1 1]*acli1/2),
 subplot(3,3,2)
            surf(XP,YP,real(Pvecty/Nor)),
            shading('interp'), view(0,90),
            colorbar
            title(' Pointing_y ')
            axis square, axis equal, grid,
            axis off
            axis([-1 1 -1 1]*acli1/2),
 subplot(3,3,3)
            surf(XP,YP,real(Pvectz/Nor)),
            shading('interp'), view(0,90),
            colorbar
            title(' Pointing_z ')
            axis square, axis equal, grid,
            axis off
            axis([-1 1 -1 1]*acli1/2),

 cE=real(max([max(max(Exdu)) max(max(Eydu)) max(max(Ezdu))]));
 cEx=max(max(Exdu));
 cEy=max(max(Eydu));
 cEz=max(max(Ezdu));
 fEx=abs(cEx/cE);
 fEy=abs(cEy/cE);
 fEz=abs(cEz/cE);

 subplot(3,3,4)
            surf(XP,YP,real(Exdu/cEx*fEx)),
            shading('interp'), view(0,90),
            colorbar
            title(' E_x ')
            axis square, axis equal, grid,
            axis([-1 1 -1 1]*acli1/2),
            axis off
 subplot(3,3,5)
            surf(XP,YP,real(Eydu/cEy*fEy)),
            shading('interp'), view(0,90),
            colorbar
            title(' E_y ')
            axis square, axis equal, grid,
            axis off
            axis([-1 1 -1 1]*acli1/2),
 subplot(3,3,6)
            surf(XP,YP,real(Ezdu/cEz*fEz)),
            shading('interp'), view(0,90),
            colorbar
            title(' E_z ')
            axis square, axis equal, grid,
            axis off
            axis([-1 1 -1 1]*acli1/2),


 cH=real(max([max(max(Hxdu)) max(max(Hydu)) max(max(Hzdu))]));
 cHx=max(max(Hxdu));
 cHy=max(max(Hydu));
 cHz=max(max(Hzdu));
 fHx=abs(cHx/cH);
 fHy=abs(cHy/cH);
 fHz=abs(cHz/cH);


 subplot(3,3,7)
            surf(XP,YP,real(Hxdu/cHx*fHx)),
            shading('interp'), view(0,90),
            colorbar
            title(' H_x ')
            axis square, axis equal, grid,
            axis([-1 1 -1 1]*acli1/2),
            axis off
 subplot(3,3,8)
            surf(XP,YP,real(Hydu/cHy*fHy)),
            shading('interp'), view(0,90),
            colorbar
            title(' H_y ')
            axis square, axis equal, grid,
            axis([-1 1 -1 1]*acli1/2),
            axis off
 subplot(3,3,9)
            surf(XP,YP,real(Hzdu/cHz*fHz)),
            shading('interp'), view(0,90),
            colorbar
            title(' H_z ')
            axis square, axis equal, grid,
            axis([-1 1 -1 1]*acli1/2),
            axis off


figure
 pograp=[100   300   880   300];
 set(gcf,'Position',pograp)
 cE=real(max([max(max(Exdu)) max(max(Eydu)) max(max(Ezdu))]));
 cEx=max(max(Exdu));
 cEy=max(max(Eydu));
 cEz=max(max(Ezdu));
 fEx=abs(cEx/cE);
 fEy=abs(cEy/cE);
 fEz=abs(cEz/cE);

 subplot(1,3,1)
 Ed=real(Exdu/cEx*fEx);
 Ed(1,1)=-max(max(Ed));
            surf(XP,YP,Ed),
            shading('interp'), view(0,90),
            colorbar
            title(' E_x ')
            axis square, axis equal, grid,
            axis([-1 1 -1 1]*acli1/2),
            axis off
 subplot(1,3,2)
 Ed=real(Eydu/cEy*fEy);
 Ed(1,1)=-max(max(Ed));

            surf(XP,YP,Ed),
            shading('interp'), view(0,90),
            colorbar
            title(' E_y ')
            axis square, axis equal, grid,
            axis off
            axis([-1 1 -1 1]*acli1/2),
 subplot(1,3,3)
 Ed=real(Ezdu/cEz*fEz);
 Ed(1,1)=-max(max(Ed));
            surf(XP,YP,Ed),
            shading('interp'), view(0,90),
            colorbar
            title(' E_z ')
            axis square, axis equal, grid,
            axis off
            axis([-1 1 -1 1]*acli1/2),

colormap gray
brighten(-.4)

figure
 pograp=[100   300   880   300];
 set(gcf,'Position',pograp)
 subplot(1,3,1)
 Ed=abs(Exdu/cEx*fEx);
            surf(XP,YP,Ed),
            shading('interp'), view(0,90),
            colorbar
            title(' E_x ')
            axis square, axis equal, grid,
            axis([-1 1 -1 1]*acli1/2),
            axis off
 subplot(1,3,2)
 Ed=abs(Eydu/cEy*fEy);

            surf(XP,YP,Ed),
            shading('interp'), view(0,90),
            colorbar
            title(' E_y ')
            axis square, axis equal, grid,
            axis off
            axis([-1 1 -1 1]*acli1/2),
 subplot(1,3,3)
 Ed=abs(Ezdu/cEz*fEz);
            surf(XP,YP,Ed),
            shading('interp'), view(0,90),
            colorbar
            title(' E_z ')
            axis square, axis equal, grid,
            axis off
            axis([-1 1 -1 1]*acli1/2),

 colormap gray
end
