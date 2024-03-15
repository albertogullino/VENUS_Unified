function [DTmax,iTz,iTr]=f_maxMatrix(DeltaT,r,z,IPLOT)

DTmax=max(max(DeltaT));
[iTz,iTr]=find(DeltaT==DTmax);

% iTr=iTr+1;

if IPLOT==1
    figure
    hold on
    surf(r,z,DeltaT)
    shading interp
    xlabel('\rho, \mum'),ylabel('z, \mum'),zlabel('\DeltaT, K')
    plot3(r(iTr),z(iTz),DeltaT(iTz,iTr),'k*')
end

return