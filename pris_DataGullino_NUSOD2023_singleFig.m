clear
clear global
close all
colordef white
dbstop if error

addpathVENUS


A=load('LW_MarkusN_FINALE_LgDBR_NUSOD2023.mat');
B=load('LW_MarkusN_TJ_oxBelow_LgDBR_NUSOD2023.mat');

colo='rgbcmrb';
%col{1}='k'; col{2}='r'; col{3}='g'; col{4}='b'; col{5}='c'; col{6}='m'; col{7}='r.'; col{8}='b.';
col{1}='r'; col{2}='g'; col{3}='b'; col{4}='c'; col{5}='m'; col{6}='k'; col{7}='y';


newcolors={'#FF8A00','#7840CC','#40CC8A','#D42424'};
% 3 colors plots
% newcolors = {'#ffb60a','m','c'};
clear col, col=newcolors;

IMAX=12;
VMAX=2.9;
iDeLa=0;

TTT000=273;

IPvet=11;

SETTA_loopsTemp

A.modeold=A.mode;
B.modeold=B.mode;

for IPAR=IPvet
    
    TMa=[]; TMb=[];
    
    if IPAR>0
        ppar=PMAT{IPAR}
    else
        ppar=[];
    end
    
    if isfield(B.mode,'nBTJ')
        nBTJ=B.mode.nBTJ;
    end
    
    % 20C is excluded in thesis, discussed separately
    for kpar=1:length(A.MODEplot)-1
        kpar
        if exist('nBTJ')
            B.MODEplot{kpar}.nBTJ=nBTJ;
        end
        A.mode=A.MODEplot{kpar};
        B.mode=B.MODEplot{kpar};
        if isfield(A.mode,'FLos')
            FLos=A.mode.FLos;
        else
            FLos=FLO;
        end
        A.modePlot=A.mode;
        B.modePlot=B.mode;
        if isfield(A.modePlot,'Isize')
            Isize=A.modePlot.Isize;
            %             T0=modePlot.T0-TTT000;
            T0=A.modePlot.Temp(1,1,1)-TTT000;
            %    'prima di load', keyboard
            %fibar=strfind(nomeSR,'\');
            %nomeSR=nomeSR(fibar(end)+1:end);
                eval(['load MarkusN_',num2str(Isize),'_T',num2str(T0),'.mat'])
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
        diffV=find(diff(Vmeas)==0,1);
        Vmeas=Vmeas(1:diffV);
        Imeas=Imeas(1:diffV);
        Lmeas=Lmeas(1:diffV);
        Rmeas=diff(Vmeas)./diff(Imeas/1000);
        Rmeas=[Rmeas(1:end-2)];
        Imeas_res=Imeas(1:end-3);
        Imax=max(Imeas);
        Imax=IMAX;
        Pmax=max(Lmeas);
        
        
        kcol=input('Color index?  \n');
        if isempty(kcol)
            kcol=kpar
        end
        
        
        
        %mesh=MESH{kpar};
        I_A=A.mode.ii_dd*1000;
        I_B=B.mode.ii_dd*1000;
        %'pa', keyboard
        V_A=A.mode.vv_dd;
        V_B=B.mode.vv_dd;
        P_A=sum(A.mode.Pst_dd,1)+A.mode.Psp_dd;
        P_B=sum(B.mode.Pst_dd,1)+B.mode.Psp_dd;
        Imax=max([IMAX max(I_B)]);
        Pmax=max([max(Lmeas) max(P_B)]);
        Imod_A{kpar}=A.modePlot.ii_dd*1000;
        Pmod_A{kpar}=10*log10(A.modePlot.Pst_dd);
        Imod_B{kpar}=B.modePlot.ii_dd*1000;
        Pmod_B{kpar}=10*log10(B.modePlot.Pst_dd);
        
        %% IV
        figure(70)
        hold on
        grid on
        box on
        pa=1:5:length(Imeas);
        plot(Vmeas(pa),Imeas(pa),'.','color',col{kcol},'markersize',10)
        hold on
        axis([1.4 VMAX 0 IMAX])
        plot(V_A,I_A,'color',col{kcol},'linewidth',2)
        plot(V_B,I_B,'--','color',col{kcol},'linewidth',2)
        %         plot(V,I,'k--','linewidth',2)
        ylabel('Current, mA')
        xlabel('Voltage, V')
        set(gcf,'pos',[220   131   560   522])
        set(gca,'FontSize',16,'FontName','Times new roman')

        %% LI
        figure(71)
        hold on
        grid on
        box on
        plot(Imeas(pa),Lmeas(pa),'.','color',col{kcol},'markersize',10)
        hold on

        if kpar==1
            axis([0 IMAX 0 Pmax*1.1])
            Pold=Pmax;
        else
            if Pmax>Pold
                axis([0 IMAX 0 Pmax*1.1])
            end
        end
        plot(I_A,P_A,'color',col{kcol},'linewidth',2)
        plot(I_B,P_B,'--','color',col{kcol},'linewidth',2)
        
        xlabel('Current, mA')
        ylabel('Optical power, mW')
        set(gcf,'pos',[220   131   560   522])
        set(gca,'FontSize',16,'FontName','Times new roman')
        
        %% T(I)
        figure(72)
        hold on
        grid on
        box on
        plot(A.mode.ii_dd*1000,A.mode.DeltaTmax,'color',col{kcol},'linewidth',2)
        plot(B.mode.ii_dd*1000,B.mode.DeltaTmax,'--','color',col{kcol},'linewidth',2)

        xlim([0 IMAX])
        xlabel('Current, mA')
        ylabel('Temperature rise, K')
        set(gcf,'pos',[220   131   560   522])
        set(gca,'FontSize',16,'FontName','Times new roman')
        
        %% lambda(I)
        figure(73)
        hold on
        grid on
        box on
        
        DeLam=0;
        
        TMa=[TMa A.mode.DeltaTmax(end)];
        fia=find(1000*A.mode.ii_dd>.1);
        plot(1000*A.mode.ii_dd(fia),A.mode.lambda(1,fia)+DeLam,'color',col{kcol},'linewidth',2)
        TMb=[TMb A.mode.DeltaTmax(end)];
        fib=find(1000*B.mode.ii_dd>.1);
        plot(1000*B.mode.ii_dd(fib),B.mode.lambda(1,fib)+DeLam+.4,'--','color',col{kcol},'linewidth',2)

        plot(CUR,LAM(:,1),'.','color',col{kcol},'Markersize',18)
        axis([1 IMAX 848 855])
        xlabel('Current, mA')
        ylabel('Wavelength, nm')
        set(gcf,'pos',[220   131   560   522])
        set(gca,'FontSize',16,'FontName','Times new roman')

        %% lambda(I)
        figure(74)
        hold on
        grid on
        box on
        
        ResistenceA=diff(A.mode.vv_dd)./diff(A.mode.ii_dd);
        ResistenceA=[ResistenceA(1),ResistenceA];
        ResistenceB=diff(B.mode.vv_dd)./diff(B.mode.ii_dd);
        ResistenceB=[ResistenceB(1),ResistenceB];

        pu=1:4:length(Rmeas);
        semilogy(Imeas_res(pu),Rmeas(pu),'.','color',col{kcol},'Markersize',12)
        
        lresA=1:length(ResistenceA);
        semilogy(A.mode.ii_dd(lresA)*1e3,ResistenceA(lresA),'color',col{kcol},'linewidth',2)
        lresB=1:length(ResistenceB);
        semilogy(B.mode.ii_dd(lresB)*1e3,ResistenceB(lresB),'--','color',col{kcol},'linewidth',2)
        
        axis([0 IMAX 20 300])
        
        xlabel('Current, mA')
        ylabel('Differential resistance, \Omega')
        set(gcf,'pos',[220   131   560   522])
        set(gca,'FontSize',16,'FontName','Times new roman')
        
        %% WPE
        figure(75)
        hold on
        grid on
        box on
        
        plot(Imeas(pa),Lmeas(pa)./(Vmeas(pa).*Imeas(pa))*100,'.','color',col{kcol},'markersize',12)
        plot(I_A,(P_A-A.mode.Psp_dd)./(I_A.*V_A)*100,'color',col{kcol},'linewidth',2)
        plot(I_B,(P_B-B.mode.Psp_dd)./(I_B.*V_B)*100,'--','color',col{kcol},'linewidth',2)
        ylabel('\eta_{WP}, %')
        xlabel('Current, mA')
        xlim([0 IMAX])
        legend('Exp - \it pin','VENUS - \it pin','VENUS - TJ','location','northeast')
        set(gcf,'pos',[220   131   560   522])
        set(gca,'FontSize',16,'FontName','Times new roman')
                
        pausak

    end
end

pausak
HeatSourcePlot(A.mesh,A.modeold,A.MODEplot{4})

HeatSourcePlot(B.mesh,B.modeold,B.MODEplot{4})