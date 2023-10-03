%po_pro
%po_prou
po_pru

if iLP==0
 Poim=sqrt(Pvectx.^2+Pvecty.^2+Pvectz.^2);
 Poimn=sqrt(Pvectxn.^2+Pvectyn.^2+Pvectzn.^2);
 Pnorm=(Pvectx.*Xpr+Pvecty.*Ypr+Pvectz.*Zpr);
 Pnormn=(Pvectxn.*Xpr+Pvectyn.*Ypr+Pvectzn.*Zpr);
 fi0=find(Poim==0);
 Poim(fi0)=1e-50;
 fi0=find(Poimn==0);
 Poimn(fi0)=1e-50;
 Pver=(Pvectx.*Xpr+Pvecty.*Ypr+Pvectz.*Zpr)./Poim;
 Pvern=(Pvectxn.*Xpr+Pvectyn.*Ypr+Pvectzn.*Zpr)./Poimn;
 Potr=sqrt(Pvectx.^2+Pvecty.^2);
 Potrn=sqrt(Pvectxn.^2+Pvectyn.^2);
 te=acos(Pvectz./Poim);
 Poi=Pvectx+Pvecty+Pvectz;
 Poin=Pvectxn+Pvectyn+Pvectzn;
 Nor=max(max(Poim));
 Norn=max(max(Poimn));
 isalta=1;
else
 Pnorm=Pvectz;
 isalta=1;
end


if isalta==0
 figure;
 pograp=[400   0   580   680];
 set(gcf,'Position',pograp)
 subplot(3,2,1)
            surf(X1,Y1,Pvectx/Nor),
            shading('interp'), view(0,90),
            colorbar
            title(' Pointing_x ')
            axis square, axis equal, grid,
            axis([-1 1 -1 1]*axli/2),
 subplot(3,2,2)
            surf(X1,Y1,Pvectxn/Norn),
            shading('interp'), view(0,90),
            colorbar
            title(' Pointing_x ')
            axis square, axis equal, grid,
            axis([-1 1 -1 1]*axli/2),

 subplot(3,2,3)
            surf(X1,Y1,Pvecty/Nor),
            shading('interp'), view(0,90),
            colorbar
            title(' Pointing_y ')
            axis square, axis equal, grid,
            axis([-1 1 -1 1]*axli/2),
 subplot(3,2,4)
            surf(X1,Y1,Pvectyn/Norn),
            shading('interp'), view(0,90),
            colorbar
            title(' Pointing_y ')
            axis square, axis equal, grid,
            axis([-1 1 -1 1]*axli/2),

 subplot(3,2,5)
            surf(X1,Y1,Pvectz/Nor),
            shading('interp'), view(0,90),
            colorbar
            title(' Pointing_z ')
            axis square, axis equal, grid,
            axis([-1 1 -1 1]*axli/2),
 subplot(3,2,6)
            surf(X1,Y1,Pvectzn/Norn),
            shading('interp'), view(0,90),
            colorbar
            title(' Pointing_z ')
            axis square, axis equal, grid,
            axis([-1 1 -1 1]*axli/2),
            pausak


 cEx=max(max(Efx));
 cEy=max(max(Efy));
 cEz=max(max(Efz));
 cE=max([cEx cEy cEz]);
 fEx=abs(cEx/cE);
 fEy=abs(cEy/cE);
 fEz=abs(cEz/cE);

 subplot(1,3,1)
 Ed=real(Exdu/cEx*fEx);


 figure;
 pograp=[400   0   580   680];
 set(gcf,'Position',pograp)
 subplot(3,2,1)

            surf(X1,Y1,real(Efx/cEx*fEx)),
            shading('interp'), view(0,90),
            colorbar
            title(' Ef_x ')
            axis square, axis equal, grid,
            axis([-1 1 -1 1]*axli/2),
 subplot(3,2,2)
            surf(X1,Y1,real(Efxn/cEx*fEx)),
            shading('interp'), view(0,90),
            colorbar
            title(' Ef_x n ')
            axis square, axis equal, grid,
            axis([-1 1 -1 1]*axli/2),

 subplot(3,2,3)
            surf(X1,Y1,real(Efy/cEy*fEy)),
            shading('interp'), view(0,90),
            colorbar
            title(' Ef_y ')
            axis square, axis equal, grid,
            axis([-1 1 -1 1]*axli/2),

 subplot(3,2,4)
            surf(X1,Y1,real(Efyn/cEy*fEy)),
            shading('interp'), view(0,90),
            colorbar
            title(' Ef_y n ')
            axis square, axis equal, grid,
            axis([-1 1 -1 1]*axli/2),

 subplot(3,2,5)
            surf(X1,Y1,real(Efz/cEz*fEz)),
            shading('interp'), view(0,90),
            colorbar
            title(' Ef_z ')
            axis square, axis equal, grid,
            axis([-1 1 -1 1]*axli/2),
 subplot(3,2,6)
            surf(X1,Y1,real(Efzn/cEz*fEz)),
            shading('interp'), view(0,90),
            colorbar
            title(' Ef_zn ')
            axis square, axis equal, grid,
            axis([-1 1 -1 1]*axli/2),
            pausak

fi=find(Poim/Nor>1e-3);
lim=max(X1(fi));

 figure;
 pograp=[400   50   580   480];
 set(gcf,'Position',pograp)
 subplot(2,2,1)
            surf(X1,Y1,Pver),
            shading('interp'), view(0,90),
            colorbar
            title(' Verifica ')
            axis square, axis equal, grid,
            axis([-1 1 -1 1]*lim),
 subplot(2,2,2)
            surf(X1,Y1,Pvern),
            shading('interp'), view(0,90),
            colorbar
            title(' Verifica nuova ')
            axis square, axis equal, grid,
            axis([-1 1 -1 1]*lim),

 subplot(2,2,3)
            surf(X1,Y1,Poim/Nor),
            shading('interp'), view(0,90),
            colorbar
            title(' Pointing ')
            axis square, axis equal, grid,
            axis([-1 1 -1 1]*lim),
 subplot(2,2,4)
            surf(X1,Y1,Poimn/Norn),
            shading('interp'), view(0,90),
            colorbar
            title(' Pointing nuovo ')
            axis square, axis equal, grid,
            axis([-1 1 -1 1]*lim),




 figure;
 pograp=[400   50   580   280];
 set(gcf,'Position',pograp)
 subplot(1,2,1)
            surf(X1,Y1,Potr/Nor),
            shading('interp'), view(0,90),
            colorbar
            title(' Point_ trasv ')
            axis square, axis equal, grid,
            axis([-1 1 -1 1]*lim),
 subplot(1,2,2)
            surf(X1,Y1,Potrn/Norn),
            shading('interp'), view(0,90),
            colorbar
            title(' Point_ trasv nuovo ')
            axis square, axis equal, grid,
            axis([-1 1 -1 1]*lim),


else
 if iFFte==1
  X=X1;
  Y=Y1;
 else
  X=rox*fm0;
  Y=rox*gm0;
 end
% ' FFF'
% keyboard
 Ef=sqrt(Pnorm);

end

