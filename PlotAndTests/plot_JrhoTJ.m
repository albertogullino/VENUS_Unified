close all

if exist('JN_X')==0
    ITJ_finder

    mode=MODEplot{end};
    
    JN_X=mode.JXn;    % actually MODEplot{1} or modePlot
    JN_Y=(mode.JYn);
    JP_X=mode.JXp;
    JP_Y=(mode.JYp);

    % JN_X=reshape(mode.Jn_x,mesh.nny,mesh.nnx);
    % JN_Y=-reshape(mode.Jn_y,mesh.nny,mesh.nnx);
    % JP_X=reshape(mode.Jp_x,mesh.nny,mesh.nnx);
    % JP_Y=-reshape(mode.Jp_y,mesh.nny,mesh.nnx);

end

    x=mesh.xgrid*1e4;
    y=mesh.ygrid*1e4;

CORRENTI=[1 2 4];
pcor=CurrIndex(CORRENTI(1),mode.ii_dd*1e3);

% v=180:10:240;
v=198:2:210;    % inside the TJ
v=ITJ-5:3:ITJ(end)+1;
clear lgd
	
ND=reshape(mesh.dop_d,mesh.nny,mesh.nnx);
NA=reshape(mesh.dop_a,mesh.nny,mesh.nnx);

figure(v(1))
subplot(211)
hold on
plot(y,ND(:,1),'.-')
plot(y,NA(:,1),'.-')
plot(y(v),ND(v,1),'ko')
plot(y(v),NA(v,1),'ko')
xlim([y(v(1)) y(v(end))])
set(gca,'yscale','log')
subplot(212)
hold on
plot(ND(:,1),'.-')
plot(NA(:,1),'.-')
plot(v,ND(v,1),'ko')
plot(v,NA(v,1),'ko')
xlim([v(1) v(end)])
set(gca,'yscale','log')
% subplot(313)
% hold on

pausak

iiii=1;

for ilong=v
    figure(v(1))
    subplot(211)
    hold on
    plot(y(ilong),ND(ilong,1),'k*')
    plot(y(ilong),NA(ilong,1),'k*')
    subplot(212)
    hold on
    plot(ilong,ND(ilong,1),'k*')
    plot(ilong,NA(ilong,1),'k*')
%     subplot(313)
%     plot(x,ND(ilong,:),'.-')
%     hold on
%     plot(x,NA(ilong,:),'.-')
%     set(gca,'yscale','log')
%     hold off
    
    figure(iiii+500)
	set(gcf,'position',[iiii*100    -14    560    1000])
	
	subplot(411)
    plot(x,squeeze(JN_X(pcor,ilong,:)),'--','linewidth',1.5)
    %chold
	hold on,grid on
    plot(x,squeeze(JP_X(pcor,ilong,:)),'--','linewidth',1.5)
	    title(num2str(ilong))

	chold
	plot(x,abs(squeeze(JN_Y(pcor,ilong,:))),'linewidth',1.5)
    plot(x,abs(squeeze(JP_Y(pcor,ilong,:))),'linewidth',1.5)
	legend('J_n(\rho)','J_p(\rho)','J_n(z)','J_p(z)','location','best')

%     hold off
	
	subplot(412)
    plot(x,ND(ilong,:),'.-')
    hold on
    plot(x,NA(ilong,:),'.-')
    set(gca,'yscale','log')
    ylabel('Radial doping density, 1/cm^3'),legend('ND','NA','location','best')
%     hold off
    	
    subplot(413)     
    plot(x,abs(squeeze(JN_X(pcor,ilong,:))),'--','linewidth',1.5)
    %chold
	hold on
    plot(x,abs(squeeze(JP_X(pcor,ilong,:))),'--','linewidth',1.5)
	
	chold
	plot(x,abs(squeeze(JN_Y(pcor,ilong,:))),'linewidth',1.5)
    plot(x,abs(squeeze(JP_Y(pcor,ilong,:))),'linewidth',1.5)
	
% 	chold
% 	plot(x,ND(ilong,:)/max(ND(ilong,:))*max(abs(squeeze(JN_X(pcor,ilong,:)))),'--')
% 	plot(x,NA(ilong,:)/max(NA(ilong,:))*max(abs(squeeze(JP_X(pcor,ilong,:)))),'--')
	grid on%,legend('J_n(\rho)','J_p(\rho)','J_n(z)','J_p(z)','location','best')
%     hold off
    ylabel('Log. scale')
    set(gca,'yscale','log')
    
	
	subplot(414)
	hold on
	plot(x,sqrt(abs(squeeze(JN_X(pcor,ilong,:))).^2+abs(squeeze(JN_Y(pcor,ilong,:))).^2),'linewidth',1.5)
	plot(x,sqrt(abs(squeeze(JP_X(pcor,ilong,:))).^2+abs(squeeze(JP_Y(pcor,ilong,:))).^2),'linewidth',1.5)	
	grid on%,hold off
	xlabel('\rho, \mum')
	ylabel('SUM')
    lgd{iiii}=num2str(ilong);
%     set(gca,'FontSize',12)
	
    iiii=iiii+1;

    pausak
end

keyboard

%% RHO    
[~,nox]=min(abs(x-mode.rox)); % Rox index
vrho=[2 nox-5 nox-1 nox nox+1 nox+10];

clear lgd

iiir=1;

for irho=vrho
%     figure(vrho(1)+250)
%     subplot(211)
%     hold on
%     plot(y,ND(:,irho),'.-')
%     plot(y,NA(irho,:),'.-')
%     plot(y(irho),ND(irho,1),'k*')
%     plot(y(irho),NA(irho,1),'k*')
%     subplot(212)
%     hold on
%     plot(irho,ND(irho,1),'k*')
%     plot(irho,NA(irho,1),'k*')
%     subplot(313)
%     plot(x,ND(irho,:),'.-')
%     hold on
%     plot(x,NA(irho,:),'.-')
%     set(gca,'yscale','log')
%     hold off
    
    figure(iiir+1000)
	set(gcf,'position',[iiir*100    -14    560    1000])
	
	subplot(411)
    plot(y,squeeze(JN_X(pcor,:,irho)),'--','linewidth',1.5)
    %chold
	hold on,grid on
    plot(y,squeeze(JP_X(pcor,:,irho)),'--','linewidth',1.5)
    title(['irho = ',num2str(irho),' - \rho = ',num2str(x(irho)),'\mum'])

	chold
	plot(y,abs(squeeze(JN_Y(pcor,:,irho))),'linewidth',1.5)
    plot(y,abs(squeeze(JP_Y(pcor,:,irho))),'linewidth',1.5)
	legend('J_n(\rho)','J_p(\rho)','J_n(z)','J_p(z)','location','best')
    xlim([y(v(1))-.5 y(v(end))+.5])
%     xlim([y(v(1)) y(v(end))])    
%     hold off
	
	subplot(412)
    plot(y,ND(:,irho),'.-')
    hold on
    plot(y,NA(:,irho),'.-')
    set(gca,'yscale','log')
    ylabel('Radial doping density, 1/cm^3'),legend('ND','NA','location','best')
    xlim([y(v(1))-.5 y(v(end))+.5])
%     hold off
    	
    subplot(413)     
    plot(y,abs(squeeze(JN_X(pcor,:,irho))),'--','linewidth',1.5)
    %chold
	hold on
    plot(y,abs(squeeze(JP_X(pcor,:,irho))),'--','linewidth',1.5)
	
	chold
	plot(y,abs(squeeze(JN_Y(pcor,:,irho))),'linewidth',1.5)
    plot(y,abs(squeeze(JP_Y(pcor,:,irho))),'linewidth',1.5)
	xlim([y(v(1))-.5 y(v(end))+.5])
% 	chold
% 	plot(x,ND(irho,:)/max(ND(irho,:))*max(abs(squeeze(JN_X(pcor,irho,:)))),'--')
% 	plot(x,NA(irho,:)/max(NA(irho,:))*max(abs(squeeze(JP_X(pcor,irho,:)))),'--')
	grid on%,legend('J_n(\rho)','J_p(\rho)','J_n(z)','J_p(z)','location','best')
%     hold off
    ylabel('Log. scale')
    set(gca,'yscale','log')
    
	
	subplot(414)
	hold on
	plot(y,sqrt(abs(squeeze(JN_X(pcor,:,irho))).^2+abs(squeeze(JN_Y(pcor,:,irho))).^2),'linewidth',1.5)
	plot(y,sqrt(abs(squeeze(JP_X(pcor,:,irho))).^2+abs(squeeze(JP_Y(pcor,:,irho))).^2),'linewidth',1.5)
	grid on%,hold off
	zlabel('z, \mum')
	ylabel('SUM')
    lgd{iiir}=num2str(irho);
%     set(gca,'FontSize',12)
    xlim([y(v(1))-.5 y(v(end))+.5])
	
    iiir=iiir+1;

    pausak
end

