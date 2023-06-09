
colordef white

%'entrp', keyboard
isperF='old';
%isperF='new';
irel=modePlot.irel;
indv=modePlot.indv;
corren=mode.ii_dd*1000;
PAc=4;
COR=CORRENTI;
if exist('COR')
    puC=[];
    for kc=1:length(COR)
        ci=COR(kc);
        [du, ic]=min(abs(corren-ci));
        puC=[puC ic];
    end
else
    puC=1:PAc:length(corren);
end

x=modePlot.x;
Te=mode.Temp;
Tqw=squeeze(Te(puC,150,:));
figure(30), plot(x,Tqw,'linewidth',2)
pausak

% len=modePlot.len;
vind=modePlot.vind;
NVbias=modePlot.NVbias;
%Imeas_res=modePlot.Imeas_res;
%Rmeas=modePlot.Rmeas;
if(modePlot.oflg)
    xfi=mode.x;
    Ca=mode.E2;
    
    COL='rgbcm';
    dox=modePlot.rox;
    Ca1=squeeze(Ca(puC,1,:));
    %       Ca2=squeeze(Ca(1:PAc:end,2,:));
    hc=figure;
    set(hc,'pos',[256   425   915   525])
    
    subplot(121)
    semilogy(xfi,Ca1,COL(1)), hold on,  ResetColor,
    for kc=2:size(Ca,2)
        Ca3=squeeze(Ca(puC,kc,:));
        semilogy(xfi,Ca3,COL(kc)), ResetColor,
    end
    xlim([1e-3 max([2*dox 5])])
    subplot(122)
    plot(xfi,Ca1,COL(1)), hold on,  ResetColor,
    for kc=2:size(Ca,2)
        Ca3=squeeze(Ca(puC,kc,:));
        plot(xfi,Ca3,COL(kc)), ResetColor,
    end
    xlim([1e-3 max([2*dox 5])])
    
    pausak
    Comp_Rec
    figure
    loglog(mode.ii_dd*1000,abs(Rnor),'LineWidth',1.5), hold on
    %         loglog(mode.ii_dd*1000,Nlea,'ko','LineWidth',2)
    legend('IntSRH','IntRad','IntAug','IntStim','IntLea','location','best')
    title('Recombinations/Ccap')
    xlabel('Corrente (mA)')
    xlim([1e-3 10])
    ylim([1e-4 1])
    grid
    pausak
    
    
    
    figure(1238),clf
    subplot(121)
    hold on
    grid on
    set(gcf,'pos',[1193         577         700         349])
    len=1:length(modePlot.ii_dd);
    plot(1000*modePlot.ii_dd,1e-12*modePlot.nMaxVet(len)/modePlot.NMQW,'o-',...
        1000*modePlot.ii_dd,1e-12*modePlot.pMaxVet(len)/modePlot.NMQW,'+-')
    hold on
    plot(1000*modePlot.ii_dd,1e-18*modePlot.n3MaxVet(len),'--',...
        1000*modePlot.ii_dd,1e-18*modePlot.p3MaxVet(len),'--')
    xlabel('Current (mA)')
    ylabel('N-P density (PER well) 1e12/cm^2)')
    
    subplot(122)
    hold on
    grid on
    
    %     nQW=0; pQW=0;
    %     for indQW=1:mesh.NMQW
    %         nQW=nQW+modePlot.nQW{end}{indQW};
    %         pQW=pQW+modePlot.pQW{end}{indQW};
    %     end
    leQW=length(modePlot.nQW{end}{2});
    x=modePlot.x;
    xQW=modePlot.x(1:leQW);
    plot(xQW,modePlot.nQW{end}{2},'b.-',xQW,modePlot.pQW{end}{2},'r.-')
    %hold on
    %plot(x,modePlot.Cn(end,:),'c.',x,modePlot.Cn(end,:),'m.')
    %    K>> figure, plot(x,mode.Cp*8e-7*1e-9/9)
    %K>> figure, plot(x,mode.Cn*8e-7*1e-9/9)
    
    xlim([xQW(1),xQW(end)])
    xlabel('\rho, nm')
    ylabel('N-P density (average) 1e12/cm^2)')
    
    %    if isopen(40)
    %    close(40)
    %    end
    figure(40),
    for kl=1:length(puC)
        plot(xQW,modePlot.nQW{puC(kl)}{2},'linewidth',2), hold on
    end
    ResetColor,
    for kl=1:length(puC)
        plot(xQW,modePlot.pQW{puC(kl)}{2},'o'),
    end
    xlabel('\rho, \mum')
    ylabel('Sheet carrier density, cm^{-2}')
    title('QW populations vs current, along \rho')
    grid on
    
    
    figure(1234),clf
    set(gcf,'Position',[75 552 600 376])
    subplot(1,2,1)
    hold on
    grid on
    box on
    plot(Vmeas,Imeas,'ro','LineWidth',2)
    plot(modePlot.vv_dd,modePlot.ii_dd*1000,'b-','LineWidth',2)
    % plot(modePlot.vv_dd,modePlot.ii_dd*1000,'b.','LineWidth',2)
    plot(modePlot.vv_dd(vind),modePlot.ii_dd(vind)*1000,'b.','markersize',20)
    %        axis([0 modePlot.vv_dd(vind(end))+.2 0 modePlot.ii_dd(vind(end))+1])
    axis([0 modePlot.vv_dd(end)+.2 0 1000*modePlot.ii_dd(end)+1])
    xlabel('Voltage, V')
    ylabel('Current, mA')
    %     title(['Vcurrent = ',num2str(v0_dd),' V'])
    
    
    subplot(1,2,2)
    hold on
    grid on
    box on
    faP=1;
    plot(Imeas,Lmeas,'ro','LineWidth',2)
    PPst=sum(modePlot.Pst_dd,1)+modePlot.Psp_dd;
    plot(modePlot.ii_dd*1000,PPst,'b-','LineWidth',2)
    plot(modePlot.ii_dd*1000,modePlot.Pst_dd,'-','LineWidth',2)
    plot(modePlot.ii_dd(vind)*1000,faP*PPst(vind),'b.','markersize',20)
    %        axis([0 modePlot.ii_dd(vind(end))+1 0 modePlot.Pst_dd(vind(end))+.5])
    axis([0 1000*modePlot.ii_dd(end)+1 0 max(sum(modePlot.Pst_dd,1))*1.1+.1])
    %        legend('Measurements','Simulation','Location','Best')
    xlabel('Current, mA')
    ylabel('Optical power, mW')
    title(['Voltage ',num2str(indv),' of ',num2str(NVbias)])
    drawnow
    nmod=mode.VelmOptions.num_azim;
    if nmod>1
        figure(1231),
        set(gcf,'pos',[83    98   529   367])
        %keyboard
        hold on
        grid on
        box on
        plot(Imeas,10*log10(Lmeas),'LineWidth',2)
        [du,icu]=min(abs(Imeas-Cur(fix(end/2))));
        
        %        'sia\.', keyboard
        P_lin=10.^(P_dB/10);
        Pls=P_lin(:,1);
        for kp=2:size(P_dB,2)
            fi=isnan(P_lin(:,kp))==0;
            if length(fi)>0
                Pls(fi)=Pls(fi)+P_lin(fi,kp);
            end
        end
        FatPow=Lmeas(icu)/Pls(fix(end/2));
        scLo=10*log10(FatPow);
        ax = gca;
        ax.ColorOrderIndex = 1;
        plot(Cur,P_dB+scLo,'o','LineWidth',2)
        
        sM=size(mode.lambda);
        %        'quo', keyboard
        %        PtotS(:,[1 3 4])=Ptot(:,[1 4 5]);
        %
        %        if sM(1)>2
        %        PtotS(:,2)=10*log10(sum(10.^(Ptot(:,[2:3])/10),2));
        %        else
        %         PtotS(:,2)=10*log10(sum(10.^(Ptot(:,[2:3])/10),2));
        %        end
        %        PtotS=PtotS+3;
        %        plot(Cur,PtotS,'-o','LineWidth',1)
        %        keyboard
        PPst=sum(modePlot.Pst_dd,1)+modePlot.Psp_dd;
        plot(modePlot.ii_dd*1000,10*log10(PPst),'b-','LineWidth',2),
        ax = gca;
        ax.ColorOrderIndex = 1;
        plot(modePlot.ii_dd*1000,10*log10(modePlot.Pst_dd),'-','LineWidth',2),
        %       plot(modePlot.ii_dd*1000,10*log10(modePlot.Pst_dd+ones(size(modePlot.Pst_dd,1),1)*modePlot.Psp_dd/size(modePlot.Pst_dd,1)),'-','LineWidth',2),
        plot(modePlot.ii_dd(vind)*1000,10*log10(PPst(vind)),'b.','markersize',15),
        %        axis([0 modePlot.ii_dd(vind(end))+1 0 modePlot.Pst_dd(vind(end))+.5])
        if 10*log10(max(sum(modePlot.Pst_dd,1))*2)>-40
            axis([0 1000*modePlot.ii_dd(end)+1 -40 10*log10(max(sum(modePlot.Pst_dd,1))*2)])
        end
        %        legend('Measurements','Simulation','Location','Best')
        
        
        xlabel('Current, mA')
        ylabel('Optical power, mW')
        title(['Voltage ',num2str(indv),' of ',num2str(NVbias)])
        drawnow
        
        %'VER:', keyboard
    end
    
    figure(1235),clf
    set(gcf,'Position',[  740    97   579   365])
    subplot(1,2,1)
    hold on
    grid on
    box on
    plot(modePlot.ii_dd*1000,modePlot.DeltaTmax,'k','LineWidth',2)
    plot(modePlot.ii_dd*1000,modePlot.DeltaTmax_srhAu,'r--','LineWidth',2)
    plot(modePlot.ii_dd*1000,modePlot.DeltaTmax_Ccap,'m--','LineWidth',2)
    plot(modePlot.ii_dd*1000,modePlot.DeltaTmax_RAD,'y--','LineWidth',2)
    plot(modePlot.ii_dd*1000,modePlot.DeltaTmax_Joule,'c--','LineWidth',2)
    plot(modePlot.ii_dd*1000,modePlot.DeltaTmax_OptAbs,'g--','LineWidth',2)
    %     suTer=modePlot.DeltaTmax_OptAbs+modePlot.DeltaTmax_Joule+modePlot.DeltaTmax_RAD+modePlot.DeltaTmax_Ccap+modePlot.DeltaTmax_srhAu;
    %     figure,     plot(modePlot.ii_dd*1000,modePlot.DeltaTmax,'k','LineWidth',2)
    %     hold on,     plot(modePlot.ii_dd*1000,suTer,'r','LineWidth',2)
    legend('Tota','SrhA','Ccap','RAD ','Joul','OpAb','Location','Best')
    plot(modePlot.ii_dd(vind)*1000,modePlot.DeltaTmax(vind),'ko','LineWidth',2)
    xlabel('Current, mA')
    ylabel('Temperature rise, K')
    
    subplot(1,2,2)
    
    %figure
    hold on
    grid on
    box on
    
    iloCur=1;
    if iloCur==0
        plot(modePlot.vv_dd,modePlot.Gmod,'k','LineWidth',2)
        plot(modePlot.vv_dd,modePlot.Lmod,'r--','LineWidth',2)
        if(modePlot.nmodes>1)
            plot(modePlot.vv_dd(vind),modePlot.Gmod(:,vind)','ko','LineWidth',2)
            plot(modePlot.vv_dd(vind),modePlot.Lmod(:,vind)','ro','LineWidth',2)
        else
            plot(modePlot.vv_dd(vind),modePlot.Gmod(vind),'ko','LineWidth',2)
            plot(modePlot.vv_dd(vind),modePlot.Lmod(vind),'ro','LineWidth',2)
        end
        xlabel('Voltage, V')
    else
        ccu=modePlot.ii_dd*1000;
        plot(ccu,modePlot.Gmod,'k','LineWidth',2)
        plot(ccu,modePlot.Lmod,'r--','LineWidth',2)
        ccu=modePlot.ii_dd(vind)*1000;
        if(modePlot.nmodes>1)
            plot(ccu,modePlot.Gmod(:,vind)','ko','LineWidth',2)
            plot(ccu,modePlot.Lmod(:,vind)','ro','LineWidth',2)
        else
            plot(ccu,modePlot.Gmod(vind),'ko','LineWidth',2)
            plot(ccu,modePlot.Lmod(vind),'ro','LineWidth',2)
        end
        xlabel('Current, mA')
    end
    ylabel('Modal gain & losses, cm^{-1}')
    legend('Gain','Losses','location','east')
    %Vv=modePlot.vv_dd(vind);
    %xlim=(ccu([1 end]));
    %keyboard
    
    %  figure, plot(modePlot.ii_dd*1000,modePlot.Gmod,'LineWidth',2)
    %   hold on
    %  ResetColor
    %  plot(1000*modePlot.ii_dd(vind),modePlot.Lmod(:,vind)','o')
    %  xlim([.2 14])
    %  ylim([30 110])
    %  grid
    %%set(gca,'yscale','lin')
    %%set(gca,'yscale','log')
    
    
    
    gipos=find(modePlot.Gmod(1,:)>0);
    if length(gipos)>1
        if iloCur==0
            xli=modePlot.vv_dd(gipos([ceil(end/2) end]) ) +[-.1,.1];
            % axis([xli 0 max(max(modePlot.Lmod))] ) % Debernardi
        else
            xli=1000*modePlot.ii_dd(gipos([ceil(end/8) end]) ) +[-1,1];
        end
        axis([xli (mean(modePlot.Gmod(1,fix(end/2):end)))*.2 modePlot.Lmod(end,end)*1.2] ) % Tibaldi
    end
    %        legend('Gain','Losses','Location','Best')
    %keyboard
    drawnow
    %'qui', keyboard
    
    if(indv>1)
        figure(1236),clf
        set(gcf,'Position',[636   560   585   365])
        subplot(1,2,1)
        
        grid on
        box on
        Resistence=diff(modePlot.vv_dd)./diff(modePlot.ii_dd);
        Resistence=[Resistence(1),Resistence];
        % plot(modePlot.ii_dd*1000,Resistence,Imeasinterp,Resmeasinterp)
        % plot(Imeasinterp,Resmeasinterp,'r-','LineWidth',2)
        semilogy(Imeas_res,Rmeas,'r.')
        hold on
        lres=1:length(Resistence);
        
        semilogy(modePlot.ii_dd(lres)*1000,Resistence(lres),'bo','LineWidth',2)
        ylim([70,500])
        %            xlim([0.5,1000*modePlot.ii_dd(end)+.1])
        xlabel('Current, mA')
        ylabel('Differential resistence, \Omega')
        %            legend('Measurements','Simulation','Location','Best')
        subplot(1,2,2)
        hold on
        grid on
        box on
        plot(CUR,LAM,'ro','LineWidth',2)
        fic=find(1000*mode.ii_dd>.1);
        plot(modePlot.ii_dd(fic)*1000,modePlot.lambda(:,fic),'b.','LineWidth',2)
        %            legend('Measurements','Simulation','Location','Best')
        if(modePlot.nmodes>1)
            plot(modePlot.ii_dd(vind)*1000,modePlot.lambda(:,vind),'b-','LineWidth',2)
        else
            plot(modePlot.ii_dd(vind)*1000,modePlot.lambda(vind),'b-','LineWidth',2)
        end
        %            xlim([0.,1000*modePlot.ii_dd(end)+.1])
        xlabel('Current, mA')
        ylabel('Modal wavelength, nm')
        drawnow
    end
    
else
    if(flagConv)
        figure(1234),clf
        set(gcf,'Position',[268 533 1096 420])
        subplot(1,2,1)
        hold on
        grid on
        box on
        plot(Vmeas,Imeas,'r','LineWidth',2)
        plot(modePlot.vv_dd,modePlot.ii_dd*1000,'b.','LineWidth',2)
        xlabel('Voltage, V')
        ylabel('Current, mA')
        title(['Vcurrent = ',num2str(v0_dd),' V'])
        %        legend('Measurements','Simulation','Location','Best')
        subplot(1,2,2)
        if(indv>1)
            hold on
            grid on
            box on
            Resistence=diff(modePlot.vv_dd)./diff(modePlot.ii_dd);
            Resistence=[Resistence(1),Resistence];
            plot(Imeasinterp,Resmeasinterp,'r-','LineWidth',2)
            plot(modePlot.ii_dd*1000,Resistence,'b.','LineWidth',2)
            ylim([75,120])
            xlabel('Current, mA')
            ylabel('Differential resistence, \Omega')
            %            legend('Measurements','Simulation','Location','Best')
            drawnow
        end
    end
end

PElec=modePlot.ii_dd*1000.*modePlot.vv_dd;
PDiss=PElec-PPst;
PDissmeas=Vmeas.*Imeas-Lmeas;
figure(1230),
if irel==0
    set(gcf,'pos',[1284         100         479         347])
else
    set(gcf,'pos',[84   102   479   347])
end


hold on
grid on
plot(Imeas,PDissmeas,'ro')
Cpl=modePlot.ii_dd.*1000;
plot(Cpl,PDiss,'b','LineWidth',2)
plot(Cpl,modePlot.PTherm,'k','LineWidth',2);
xlabel('Current, mA')
ylabel('Dissipated power, mW')
legend('Measurements','Electric - Optic','Thermic','Location','Best')

figure(1229),
set(1229,'pos',[ 1325          91         583         420])
subplot(221)
plot(Cpl,-modePlot.IntCcapN*1e3)
hold on
ax = gca;
ax.ColorOrderIndex = 1;
plot(Cpl,-modePlot.IntCcapP*1e3,'.')
title('IntCcap (mA)')
grid



subplot(222)



LeN=Cpl-sum(-modePlot.IntCcapN,2)'*1e3;
LeP=Cpl-sum(-modePlot.IntCcapP,2)'*1e3;
%  figure, subplot(121)
plot(Cpl,modePlot.Fleak.*modePlot.ii_dd.*1000,Cpl,LeN,Cpl,LeP,'--','linewidth',2)
legend('Jn/Jp','CcapN','CcapP','location','northwest')
title('Leak Cur (mA)')
%  xlim([.1 14]), subplot(122)
%    plot(Cpl,modePlot.Fleak.*modePlot.ii_dd.*100000./Cpl,Cpl,Le./Cpl*100,'linewidth',2)
%    title('Rel Leak Cur (mA)') , xlim([.1 14]),

% xlabel('Cur')

xT1=modePlot.DeltaTmax;
fiT=find(xT1>10);
xT=xT1(fiT);
lT=modePlot.lambda(1,fiT);
if length(fiT>3)
    subplot(223)
    plot(xT,lT)
    coT=polyfit(xT,lT,2);
    der=polyder(coT);
    plot(xT,polyval(der,xT))
    title('d\lambda/dT')
    xlabel('\DeltaT, K')
    ylim([.04 .07])
    grid
end
subplot(224)
Tma=modePlot.DeltaTmax;
Fat=exp(sqrt(Tma/modePlot.T_tauscat));
Fat=exp(1+Tma/modePlot.T_tauscat).^modePlot.Tcap_EXP;
plot(modePlot.ii_dd*1000,Fat,'k','LineWidth',2)
title('Fattore T_{cap}')
xlabel('Cur')

Cur=modePlot.ii_dd.*1000;
Ccap=-modePlot.IntCcapN*1e3;
Lea=modePlot.Fleak.*Cur;

CorLav=CORRENTI;
corrente=modePlot.ii_dd*1000;

pu=[];
for kC=1:length(CorLav)
    [du,puk]=min(abs(corrente-CorLav(kC)));
    pu=[pu puk];
end
fipst=pu;

pausak
close(1230:1239)
close(2:3)

Iniezione


Ifit=input(' Vuoi calcolare spettri ? [1, si;  Enter, NO]  ')
if length(Ifit)==0
    Ifit=0;
else
    Ifit=1;
end

if Ifit==1
    DD_FitSpettri
end