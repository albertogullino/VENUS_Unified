if length(MODEplot)>1
    iT=input('iT?\n');
    modep=MODEplot{iT};
else
    modep=MODEplot{1};
end
Pjoule=zeros(size(mode.ii_dd));
Prad=Pjoule;
Pcap=Pjoule;
Pfca=Pjoule;
Pnr=Pjoule;
Ptot=Pjoule;

% Computes and stores the 
for iV=1:length(mode.ii_dd)
    
    Joule=(squeeze(modep.HeatJoule(iV,:,:)));
    Rec_RAD=(squeeze(modep.HeatRec_RAD(iV,:,:)));
    Rec_Cap=(squeeze(modep.HeatRec_Cap(iV,:,:)));
    OptAbs=(squeeze(modep.HeatOptAbs(iV,:,:)));
    Rec_nr=(squeeze(modep.HeatRec_13(iV,:,:)));

    Pjoule(iV)=HeatIntegral(mesh,Joule);
    Prad(iV)=HeatIntegral(mesh,Rec_RAD);
    Pcap(iV)=HeatIntegral(mesh,Rec_Cap);
    Pfca(iV)=HeatIntegral(mesh,OptAbs);
    Pnr(iV)=HeatIntegral(mesh,Rec_nr);
    
    Ptot(iV)=HeatIntegral(mesh,Joule+Rec_RAD+Rec_Cap+OptAbs+Rec_nr);
end

I=mode.ii_dd*1e3;

figure
set(gcf,'position',[68 347 1166 601])
subplot(121)
hold on,grid on,box on
plot(I,Pjoule,'.-')
plot(I,Prad,'.-')
plot(I,Pcap,'.-')
plot(I,Pfca,'.-')
plot(I,Pnr,'.-')
plot(I,Ptot,'k.-')
xlabel('Current, mA'),ylabel('Ptherm'),axis([0 I(end) 0 max(Ptot)])
legend('Joule','RAD','Cap','FCA','NR','Total','location','best')

subplot(122)
hold on,box on,grid on
plot(I,mode.DeltaTmax_Joule,'.-')
plot(I,mode.DeltaTmax_RAD,'.-')
plot(I,mode.DeltaTmax_Ccap,'.-')
plot(I,mode.DeltaTmax_OptAbs,'.-')
plot(I,mode.DeltaTmax_srhAu,'.-')
plot(I,mode.DeltaTmax,'k.-')
ylabel('\DeltaT, K'),xlabel('Current, mA')
axis([0 max(I) 0 max(mode.DeltaTmax)])
legend('Joule','Rad','Ccap','FCA','NR','location','best')
