% clear
clear global
close all
colordef white
dbstop if error

addpathVENUS

INEW=1;

FLO=1;

IMAX=15;
VMAX=2.6;
% CORRENTI=[1 5 10];
% CORRENTI=[2 4];
CORRENTI=[1 4 8];

iDeLa=0;
TTT000=273;




IncFig=1;
F1002=1001;
F1003=2001;
F99=98;
F101=100;
F100=99;
Fcur=9;
Fall=19;

IPvet=11;

SETTA_loopsTemp


for IPAR=IPvet
    
    TMa=[];
    
    if IPAR==4
        Fpar=1e-18;
    else
        Fpar=1;
    end
    
    % keyboard
    colo='rgbcmrb';
    %col{1}='k'; col{2}='r'; col{3}='g'; col{4}='b'; col{5}='c'; col{6}='m'; col{7}='r.'; col{8}='b.';
    col{1}='r'; col{2}='g'; col{3}='b'; col{4}='c'; col{5}='m'; col{6}='k'; col{7}='y';
    if IPAR>0
        ppar=PMAT{IPAR}
    else
        ppar=[];
    end
    %keyboard
%     hh=figure;
%     set(hh,'position',[263          89        1517         891])
    
    
    %  legend(TITE{1},'location','best')
    %hh1=figure;
    %set(hh1,'position',[ 400 50 1400 400])
    %'ver', keyboard
    if isfield(mode,'nBTJ')
       nBTJ=mode.nBTJ; 
    end
        
    for kpar=1:length(MODEplot)
        kpar
        if exist('nBTJ')
            MODEplot{kpar}.nBTJ=nBTJ;
        end
        mode=MODEplot{kpar};
        if isfield(mode,'FLos')
            FLos=mode.FLos;
        else
            FLos=FLO;
        end
        modePlot=MODEplot{kpar};
        if isfield(modePlot,'Isize')
            Isize=modePlot.Isize;
%             T0=modePlot.T0-TTT000;
            T0=modePlot.Temp(1,1,1)-TTT000;
            %    'prima di load', keyboard
            %fibar=strfind(nomeSR,'\');
            %nomeSR=nomeSR(fibar(end)+1:end);
            if INEW==1
                eval(['load MarkusN_',num2str(Isize),'_T',num2str(T0),'.mat'])
            else
                eval(['load Markus_',num2str(Isize),'_T',num2str(T0),'.mat'])
            end
            PDB{kpar}=P_dB;
            CU{kpar}=Cur;
        else
            if irel==1
                eval(['load ',nomeSR,'Meas_ReliefVCSEL.mat'])
            else
                eval(['load ',nomeSR,'Meas_StandardVCSEL.mat'])
            end
        end
        CUR=Cur;
        Rmeas=diff(Vmeas)./diff(Imeas/1000);
        Rmeas=[Rmeas(1:end-2)];
        Imeas_res=Imeas(1:end-3);
        Imax=max(Imeas);
        Imax=IMAX;
        Pmax=max(Lmeas)*1.2;
        
        figure(1)
        set(gcf,'position',[263          89        1517         891])
        
        kcol=input('Color index?  \n');
        if isempty(kcol)
            kcol=kpar
        end            
        
        subplot(231)
        pa=1:5:length(Imeas);
        plot(Vmeas(pa),Imeas(pa),[colo(kpar),'.'],'markersize',10)
        hold on
        
        
        subplot(233)
        plot(Imeas(pa),Lmeas(pa),[colo(kpar),'.'],'markersize',10)
        hold on
        xlim([0 IMAX])
        %'tit', keyboard
        title(['settings-vari-',mode.settings(2:end)])
        
        if IPAR>0
            
            if kpar<length(ppar)
                tita=[tit{IPAR},'  Valore attuale= ',num2str(Fpar*ppar(kpar))];
            else
                tita=tit{IPAR};
            end
            % 'ver', keyboard
            % [fPES Gamma_z fPdif  E2 lambda TempVELM anti_gui]
            
            if IPAR==30
                titeff=['EFFETTI: ',TITE{PMAT{IPAR}(kpar)+2}]
            else
                titeff=[''];
            end
        else
            tita='Set 0';
            titeff='Set 0';
        end
        %mesh=MESH{kpar};
        I=mode.ii_dd*1000;
        %'pa', keyboard
        V=mode.vv_dd;
        P=sum(mode.Pst_dd,1)+mode.Psp_dd;
        Imax=max([IMAX max(I)]);
        Pmax=max([max(Lmeas) max(P)]) *1.1 ;
        Imod{kpar}=modePlot.ii_dd*1000;
        Pmod{kpar}=10*log10(modePlot.Pst_dd);
        
%         figure(hh)
        
        subplot(231)
        hold on
        grid on
        box on
        axis([1.4 VMAX 0 IMAX])
        plot(V,I,[col{kcol}],'linewidth',2)
%         plot(V,I,'k--','linewidth',2)
        ylabel('Current (mA)')
        xlabel('Voltage (V)')
        title(tita)
        % figure(11)
        subplot(233)
        hold on
        grid on
        box on
        if kpar==1
            axis([0 Imax 0 Pmax*1.2])
            Pold=Pmax;
        else
            if Pmax>Pold
                axis([0 Imax 0 Pmax*1.2])
            end
        end
        plot(I,P,[col{kcol},'.-'],'linewidth',2)
%         plot(I,P,'k--','linewidth',2)
        
        xlabel('Current, mA')
        ylabel('Optical power, mW')
        
        %legend(num2str([1:kpar]'),'location','best')
        %  legend([titeff0, TITE{PMAT{IPAR}(1:kpar)+1}],'location','best')
        
        subplot(232)
        hold on
        grid on
        box on
        %            xlim([0 20])
        len=1:length(mode.ii_dd);
        
        plot(1000*mode.ii_dd,1e-12*mode.nMaxVet(len)/mode.NMQW,[colo(kcol),'o-'],...
            1000*mode.ii_dd,1e-12*mode.pMaxVet(len)/mode.NMQW,[colo(kcol),'+--'])
        hold on
        plot(1000*mode.ii_dd,1e-18*mode.n3MaxVet(len),[colo(kcol),'-'],...
            1000*mode.ii_dd,1e-18*mode.p3MaxVet(len),[colo(kcol),'--'],'linewidth',2)
        
%         plot(1000*mode.ii_dd,1e-12*mode.nMaxVet(len)/mode.NMQW,'ko-',...
%             1000*mode.ii_dd,1e-12*mode.pMaxVet(len)/mode.NMQW,'k+--')
%         hold on
%         plot(1000*mode.ii_dd,1e-18*mode.n3MaxVet(len),'k-',...
%             1000*mode.ii_dd,1e-18*mode.p3MaxVet(len),'k--','linewidth',2)
        
        xlabel('Current (mA)')
        ylabel('N-P density (PER well) 1e12/cm^2)')
        legend('Electrons','Holes','location','best')
        if IPAR==30
            if kpar<=ppar(end)
                legend(TITE{[1 ppar(1:kpar)+2]},'location','best')
            else
                legend(TITE{[1 ppar(1):ppar(end)+2]},'location','best')
            end
        else
            fis= strfind(modePlot.structureName,'\');
%             fis= strfind(modePlot.structureName,'/');
			strName=modePlot.structureName(fis(end)+1:end);
            titeff=strName;
        end
        if length(titeff)>20
            title(titeff(15:end))
        else
            title(titeff)
        end
        % title(tit{IPAR})
        %figure(hh1)
        subplot(236)
        hold on
        grid on
        box on
        mI=mean(1000*mode.ii_dd);
        [du,iMii]=min(abs(1000*mode.ii_dd-mI));
        mI=1000*mode.ii_dd(iMii);
        [du,iMi]=min(abs(CUR-mI));
        LaC=mode.lambda(1,iMii);
        LaX=LAM(iMi,1);
        
        if iDeLa==1
            DeLam=LaX-LaC;
        else
            DeLam=0
            '  !!!!!!!!!!!!!!!!!!!!!! '
        end
        if isnan(LaX)
            DeLam=0;
        end
        TMa=[TMa mode.DeltaTmax(end)];
        fic=find(1000*mode.ii_dd>.1);
        plot(1000*mode.ii_dd(fic),mode.lambda(1,fic)+DeLam,col{kcol},'linewidth',2)
%         plot(1000*mode.ii_dd(fic),mode.lambda(1,fic)+DeLam,'k','linewidth',2)
        title(['DeLam = ',num2str(DeLam),'    TMa ',num2str(TMa,'%0.f- ')])
        plot(CUR,LAM(:,1),[colo(kpar),'o'],'Markersize',5)
        xlabel('Current, mA')
        ylabel('Wavelength, nm')
        subplot(234)
        %         hold on
        %          grid on
        %          box on
        
        Resistence=diff(mode.vv_dd)./diff(mode.ii_dd);
        Resistence=[Resistence(1),Resistence];
        % plot(mode.ii_dd*1000,Resistence,Imeasinterp,Resmeasinterp)
        %            plot(Imeasinterp,Resmeasinterp,'r-','LineWidth',2)
        pu=1:4:length(Rmeas);
        semilogy(Imeas_res(pu),Rmeas(pu),[colo(kpar),'o'],'Markersize',3)

        hold on
        lres=1:length(Resistence);
        
        semilogy(mode.ii_dd(lres)*1e3,Resistence(lres),[colo(kcol),'.-'],'linewidth',2)
%         semilogy(mode.ii_dd(lres)*1e3,Resistence(lres),'k--','linewidth',2)
        
%         ylim([30,500])
        ylim([30,5e3])
        grid on
        xlim([0 Imax])
        %            xlim([0.5,1000*mode.ii_dd(end)+.1])
        xlabel('Current, mA')
        ylabel('Differential resistence, \Omega')
        %            keyboard
        subplot(235)
        hold on
        grid on
        box on
        if kpar==1
            CurSav=Cur;
        end
        if FLos==1
            %figure(111)
            
            Cur=CurSav;
            
            [du,icu]=min(abs(Imeas-Cur(fix(end/2))));
            Smis='o+s';
            k=kpar;
            if IPAR==11 || IPAR==41
                
                P_dB=PDB{k};
                Cur=CU{k};
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
                plot(Cur,P_dB+scLo,[Smis(3),'--'])
%                 plot(Cur,P_dB(:,1)+scLo,[Smis(k),'--'])
                hold on
                ax = gca;
                ax.ColorOrderIndex = 1;
            else
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
                
%                 plot(Cur,P_dB+scLo,'o','LineWidth',2)
                plot(Cur,P_dB(:,1)+scLo,'o','LineWidth',2)
                hold on
            end
            
            
            ax = gca;
            ax.ColorOrderIndex = 1;
            lines{1}='-'; lines{2}='--'; lines{3}='-.';  lines{4}='.';  lines{5}='+';  lines{6}='o';
            lines{7}='p'; lines{8}='s';
            plot(Imod{k},Pmod{k},lines{mod(k,length(lines))},'linewidth',2)
            hold on
            ax = gca;
            ax.ColorOrderIndex = 1;
            %axis([0 Imod{end}(end)+1 -40 5])
            axis([0 15 -7 5])
            xlabel('Current, mA')
            ylabel('Optical power, mW')            
        
        end
        
        figure(2)
        plot(mode.ii_dd*1000,mode.DeltaTmax,col{kcol},'LineWidth',2)
        hold on
%         plot(mode.ii_dd*1000,mode.DeltaTmax,'k--','LineWidth',2)
            xlabel('Current, mA')
            ylabel('Temperature rise, K')
            
        pausak
        
%         PlotDyn
        
        
        iDDPlot=input(' Vuoi DD_PLOT? [Any key, NO; Enter, Si] ');
        
        if isempty(iDDPlot)
            Dyn_Comparison
            
            DD_PLOT_pris
%             DD_PLOT_Gullino
            pausak
        end
    end
    pausak
    
    Cur=CurSav;
    
    [du,icu]=min(abs(Imeas-Cur(fix(end/2))));
    
    figure(111)
    Smis='o+s';
    if IPvet==11 || IPvet==41
        
        for k=1:length(MODEplot)
            P_dB=PDB{k};
            Cur=CU{k};
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
            plot(Cur,P_dB+scLo,[Smis(3),'--'])
%             plot(Cur,P_dB+scLo,[Smis(k),'--'])
            hold on
            ax = gca;
            ax.ColorOrderIndex = 1;
        end
    else
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
        
        plot(Cur,P_dB+scLo,'o','LineWidth',2)
        hold on
    end
    
    
    ax = gca;
    ax.ColorOrderIndex = 1;
    lines{1}='-'; lines{2}='--'; lines{3}='-.';  lines{4}='.';  lines{5}='+';  lines{6}='o';
    for k=1:length(MODEplot)
        plot(Imod{k},Pmod{k},lines{mod(k,6)},'linewidth',2)
        hold on
        ax = gca;
        ax.ColorOrderIndex = 1;
    end
    %axis([0 Imod{end}(end)+1 -40 5])
    axis([0 15 -7 5])
    xlabel('Current, mA')
    ylabel('Optical power, mW')
    pausak
    % close all
end