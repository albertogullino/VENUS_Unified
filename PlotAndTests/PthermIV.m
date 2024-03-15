if length(MODEplot)>1
    iT=input('iT?\n');
    modep=MODEplot{iT};
else
    modep=MODEplot{1};
end

I=mode.ii_dd'*1e3;

PTherm=zeros(length(I),1);
PJoule=PTherm; Prad=PTherm; Pfca=PTherm; Pcap=PTherm; Pnr=PTherm;

for iV=1:length(I)
    
    
    Joulep=(squeeze(modep.HeatJoule(iV,:,:)));
    Rec_RADp=(squeeze(modep.HeatRec_RAD(iV,:,:)));
    Rec_Capp=(squeeze(modep.HeatRec_Cap(iV,:,:)));
    OptAbsp=(squeeze(modep.HeatOptAbs(iV,:,:)));
    Rec_nrp=(squeeze(modep.HeatRec_13(iV,:,:)));
    
%     HeatSource=(squeeze(modep.HeatRec_RAD(iV,:,:)));
    HeatSource=Joulep+Rec_RADp+Rec_Capp+OptAbsp+Rec_nrp;
    
    PTherm(iV)=HeatIntegral(mesh,HeatSource);  % Extract Ptherm for a given thermal source
    
    PJoule(iV)=HeatIntegral(mesh,Joulep);  % Extract Ptherm for a given thermal source
    Prad(iV)=HeatIntegral(mesh,Rec_RADp);  % Extract Ptherm for a given thermal source
    Pnr(iV)=HeatIntegral(mesh,Rec_nrp);  % Extract Ptherm for a given thermal source
    Pfca(iV)=HeatIntegral(mesh,OptAbsp);  % Extract Ptherm for a given thermal source
    Pcap(iV)=HeatIntegral(mesh,Rec_Capp);  % Extract Ptherm for a given thermal source
end

figure,plot(I,PTherm,'k','linewidth',2)
grid on,chold
plot(I,PJoule,'c','linewidth',2)
plot(I,Pfca,'g','linewidth',2)
plot(I,Pcap,'r','linewidth',2)
plot(I,Prad,'b','linewidth',2)
plot(I,Pnr,'y','linewidth',2)
xlabel('Current, mA'),ylabel('Dissipated power, mW')
legend('Tot','Joule','FCA','Ccap','Rad','NR','location','northwest')
ylim([0 30])
set(gca,'FontSize',16,'FontName','Times','Box','on')