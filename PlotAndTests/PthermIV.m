if length(MODEplot)>1
    iT=input('iT?\n');
    modep=MODEplot{iT};
else
    modep=MODEplot{1};
end

for iV=1:length(mode.ii_dd)
    
    
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