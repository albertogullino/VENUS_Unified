
x=mesh.xgrid;
y=mesh.ygrid;

JYn=reshape(mode.Jn_y,mesh.nny,mesh.nnx);
JYp=reshape(mode.Jp_y,mesh.nny,mesh.nnx);
JXn=reshape(mode.Jn_x,mesh.nny,mesh.nnx);
JXp=reshape(mode.Jp_x,mesh.nny,mesh.nnx);

% if length(MODEplot)>1
%     iT=input('iT?\n');
%     modep=MODEplot{iT};
% else
%     modep=MODEplot;
% end
% 
% JYn=squeeze(modep.JYn(iV,:,:));
% JYp=squeeze(modep.JYp(iV,:,:));
% JXn=squeeze(modep.JXn(iV,:,:));
% JXp=squeeze(modep.JXp(iV,:,:));

for indy=1:length(y)
%     Jn_y=squeeze(modep.JYn(iV,indy,:))';
%     Jp_y=squeeze(modep.JYp(iV,indy,:))';
    
    Jn_y=squeeze(JYn(indy,:));
    Jp_y=squeeze(JYp(indy,:));
    
    cury_n(indy)=1000*trapz(x,2.*pi.*x.*Jn_y);
    cury_p(indy)=1000*trapz(x,2.*pi.*x.*Jp_y);
end

for indx=1:length(x)
%     Jn_x=squeeze(modep.JXn(iV,:,indx))';
%     Jp_x=squeeze(modep.JXp(iV,:,indx))';
    
    Jn_x=squeeze(JXn(:,indx));
    Jp_x=squeeze(JXp(:,indx));
    
    curx_n(indx)=1000*trapz(y,Jn_x);
    curx_p(indx)=1000*trapz(y,Jp_x);
end

fP=abs((cury_n*mode.CarrierNorm)+(cury_p*mode.CarrierNorm));

if exist('StrTT')==0
    StrTT=modep.StrTT;
end

figure(322),clf
% 	subplot(211)

set(gcf,'Position',[1.3146e+03   6.6600e+01   5.6000e+02   4.2000e+02])
hold on,grid on,box on
plot(y*1e4,abs(cury_n)*mode.CarrierNorm,'b','linewidth',2)
plot(y*1e4,abs(cury_p)*mode.CarrierNorm,'r','linewidth',2)
plot(y*1e4,fP,'g--','linewidth',2)
%     axis([114.5 115.5 max(abs(cury_n*mode.CarrierNorm),[],'omitnan')/10 max(abs(cury_n*mode.CarrierNorm),[],'omitnan')])
xlim([StrTT.Tbuf+StrTT.Tdbr_inf-1 StrTT.Tbuf+StrTT.Tdbr_inf+StrTT.Tcav+1]),hold off
xlim([StrTT.Tbuf+StrTT.Tdbr_inf-1 StrTT.Tbuf+StrTT.Tdbr_inf+StrTT.Tcav+1]),hold off
% xlim([MODEplot{1}.StrTT.Tbuf+4.5 MODEplot{1}.StrTT.Tbuf+5.5])

% xlim([114 116])
xlabel('z, \mum'),ylabel('Current densities, A/cm^2')
set(gca,'yscale','log','fontsize',12)
drawnow

	
% 	subplot(212)
%     
%     plot(x,abs(JXn),'--','linewidth',1.5)
%     %chold
% 	hold on
%     plot(x,abs(JXp),'--','linewidth',1.5)
% 	
% 	chold
% 	plot(x,abs(JYn),'linewidth',1.5)
%     plot(x,abs(JYp),'linewidth',1.5)
% 		
% 	subplot(313)
% 	hold on
% 	plot(x,sqrt(abs(JXn).^2+abs(JYn).^2),'linewidth',1.5)
% 	plot(x,sqrt(abs(JXp).^2+abs(JYp).^2),'linewidth',1.5)	
% 	grid on

%     xcD=x*1e4;
%  figure
%  plot(xcD,2*pi*x.*curx_n,xcD,2*pi*x.*curx_p), 
%  hold on