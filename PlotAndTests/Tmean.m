
zbot=114.5;              % um
%     ztop=116.5;            % um

ztop=116;            % um (to verify the decrease of T after QWs)
%     ztop=y(end);            % um (to verify the decrease of T after QWs)

[~,izbot]=min(abs(zbot-y)); % extract the corresponding mesh index
[~,iztop]=min(abs(ztop-y)); % extract the corresponding mesh index
dz=y(iztop)-y(izbot);       % Length along z

% rho: r -> radial point for the area where the source is inserted
% Values used in "20230929_HeatSources_Stephan.pptx"
%     rSimpl=StrTT.Rox+1;          % um
rSimpl=StrTT.Rox;          % um
%     rSimpl=StrTT.ro_max;
%     rSimpl=StrTT.ro_mesa;

[~,iR]=min(abs(x-rSimpl));   % extract the corresponding mesh index

rho=x(1:iR);
drho=diff([0 rho]);
xdxN=rho.*drho;

z=y(izbot:iztop);
dzN=diff([z(1) z]);

DT=DeltaT(izbot:iztop,1:iR);

DTmean=2*pi*xdxN*DT'*dzN'/(pi*rSimpl^2*dz) % mW

plot(y(iTz),DTmean,'md','linewidth',1.5,'MarkerSize',10)