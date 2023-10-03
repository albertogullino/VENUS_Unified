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

IMAX=15;
VMAX=3.2;
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
    
    for kpar=1:length(A.MODEplot)
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
        
        figure(1)
        set(gcf,'position',[263          89        1517         891])
        
        kcol=input('Color index?  \n');
        if isempty(kcol)
            kcol=kpar
        end
        
        subplot(221)
        pa=1:5:length(Imeas);
        plot(Vmeas(pa),Imeas(pa),[colo(kpar),'.'],'markersize',10)
        hold on
        
        
        subplot(222)
        plot(Imeas(pa),Lmeas(pa),[colo(kpar),'.'],'markersize',10)
        hold on
        xlim([0 IMAX])
        %'tit', keyboard
        title(['settings\_vari\_',A.mode.settings(2:end)])
        
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
        %         figure(hh)
        
        subplot(221)
        hold on
        grid on
        box on
        axis([1.4 VMAX 0 IMAX])
        plot(V_A,I_A,[col{kcol}],'linewidth',2)
        plot(V_B,I_B,[col{kcol},'--'],'linewidth',2)
        %         plot(V,I,'k--','linewidth',2)
        ylabel('Current (mA)')
        xlabel('Voltage (V)')
        title(tita)
        % figure(11)
        subplot(222)
        hold on
        grid on
        box on
        if kpar==1
            axis([0 Imax 0 Pmax*1.1])
            Pold=Pmax;
        else
            if Pmax>Pold
                axis([0 Imax 0 Pmax*1.1])
            end
        end
        plot(I_A,P_A,[col{kcol},'-'],'linewidth',2)
        plot(I_B,P_B,[col{kcol},'--'],'linewidth',2)
        %         plot(I,P,'k--','linewidth',2)
        
        xlabel('Current, mA')
        ylabel('Optical power, mW')
        
        %legend(num2str([1:kpar]'),'location','best')
        %  legend([titeff0, TITE{PMAT{IPAR}(1:kpar)+1}],'location','best')
        
%         subplot(232)
%         hold on
%         grid on
%         box on
%         %            xlim([0 20])
%         plot(A.mode.ii_dd*1000,A.mode.DeltaTmax,[col{kcol}],'LineWidth',2)
%         plot(B.mode.ii_dd*1000,B.mode.DeltaTmax,[col{kcol},'--'],'LineWidth',2)
%         axis([0 Imax 0 max(B.mode.DeltaTmax)])
%         %         plot(A.mode.ii_dd*1000,A.mode.DeltaTmax,'k--','LineWidth',2)
%         xlabel('Current, mA')
%         ylabel('Temperature rise, K')
%         
%         if IPAR==30
%             if kpar<=ppar(end)
%                 legend(TITE{[1 ppar(1:kpar)+2]},'location','best')
%             else
%                 legend(TITE{[1 ppar(1):ppar(end)+2]},'location','best')
%             end
%         else
%             fis= strfind(A.modePlot.structureName,'\');
%             %             fis= strfind(A.modePlot.structureName,'/');
%             strName=A.modePlot.structureName(fis(end)+1:end);
%             titeff=strName;
%         end
%         if length(titeff)>20
%             title(titeff(15:end))
%         else
%             title(titeff)
%         end
        % title(tit{IPAR})
        %figure(hh1)
        
        subplot(224)
        hold on
        grid on
        box on
        
        DeLam=0;
        
        TMa=[TMa A.mode.DeltaTmax(end)];
        fia=find(1000*A.mode.ii_dd>.1);
        plot(1000*A.mode.ii_dd(fia),A.mode.lambda(1,fia)+DeLam,[col{kcol}],'linewidth',2)
        TMb=[TMb A.mode.DeltaTmax(end)];
        fib=find(1000*B.mode.ii_dd>.1);
        plot(1000*B.mode.ii_dd(fib),B.mode.lambda(1,fib)+DeLam+.4,[col{kcol},'--'],'linewidth',2)
%         title(['DeLam=',num2str(DeLam),' TMa ',num2str(TMa,'%0.f- '),' TMb ',num2str(TMb,'%0.f- ')])
        plot(CUR,LAM(:,1),[colo(kpar),'o'],'Markersize',5)
        axis([1 Imax-4 847 855])
        xlabel('Current, mA')
        ylabel('Wavelength, nm')
        
        subplot(223)
        %         hold on
        %          grid on
        %          box on
        
        ResistenceA=diff(A.mode.vv_dd)./diff(A.mode.ii_dd);
        ResistenceA=[ResistenceA(1),ResistenceA];
        ResistenceB=diff(B.mode.vv_dd)./diff(B.mode.ii_dd);
        ResistenceB=[ResistenceB(1),ResistenceB];
        % plot(A.mode.ii_dd*1000,Resistence,Imeasinterp,Resmeasinterp)
        %            plot(Imeasinterp,Resmeasinterp,'r-','LineWidth',2)
        pu=1:4:length(Rmeas);
        semilogy(Imeas_res(pu),Rmeas(pu),[colo(kpar),'o'],'Markersize',3)
        
        hold on
        lresA=1:length(ResistenceA);
        semilogy(A.mode.ii_dd(lresA)*1e3,ResistenceA(lresA),[colo(kcol),'-'],'linewidth',2)
        lresB=1:length(ResistenceB);
        semilogy(B.mode.ii_dd(lresB)*1e3,ResistenceB(lresB),[colo(kcol),'--'],'linewidth',2)
        
        axis([0 Imax 20,1e3])
        grid on
        
        xlabel('Current, mA')
        ylabel('Differential resistance, \Omega')
        %            keyboard
%         subplot(235)
%         hold on
%         grid on
%         box on
%         if kpar==1
%             CurSav=Cur;
%         end
%         if FLos==1
%             %figure(111)
%             
%             Cur=CurSav;
%             
%             [~,icu]=min(abs(Imeas-Cur(fix(end/2))));
%             Smis='o+s';
%             k=kpar;
%             if IPAR==11 || IPAR==41
%                 
%                 P_dB=PDB{k};
%                 Cur=CU{k};
%                 P_lin=10.^(P_dB/10);
%                 Pls=P_lin(:,1);
%                 for kp=2:size(P_dB,2)
%                     fi=isnan(P_lin(:,kp))==0;
%                     if length(fi)>0
%                         Pls(fi)=Pls(fi)+P_lin(fi,kp);
%                     end
%                 end
%                 FatPow=Lmeas(icu)/Pls(fix(end/2));
%                 scLo=10*log10(FatPow);
%                 plot(Cur,P_dB+scLo,[Smis(3),'--'])
%                 %                 plot(Cur,P_dB(:,1)+scLo,[Smis(k),'--'])
%                 hold on
%                 ax = gca;
%                 ax.ColorOrderIndex = 1;
%             else
%                 P_lin=10.^(P_dB/10);
%                 Pls=P_lin(:,1);
%                 for kp=2:size(P_dB,2)
%                     fi=isnan(P_lin(:,kp))==0;
%                     if length(fi)>0
%                         Pls(fi)=Pls(fi)+P_lin(fi,kp);
%                     end
%                 end
%                 FatPow=Lmeas(icu)/Pls(fix(end/2));
%                 scLo=10*log10(FatPow);
%                 
%                 %                 plot(Cur,P_dB+scLo,'o','LineWidth',2)
%                 plot(Cur,P_dB(:,1)+scLo,'o','LineWidth',2)
%                 hold on
%             end
%             
%             
%             ax = gca;
%             ax.ColorOrderIndex = 1;
%             lines{1}='-'; lines{2}='--'; lines{3}='-.';  lines{4}='.';  lines{5}='+';  lines{6}='o';
%             lines{7}='p'; lines{8}='s';
%             plot(Imod_A{k},Pmod_A{k},lines{mod(k,length(lines))},'linewidth',2)
%             chold
%             plot(Imod_B{k},Pmod_B{k},lines{mod(k,length(lines))},'linewidth',2)
%             ax = gca;
%             ax.ColorOrderIndex = 1;
%             %axis([0 Imod{end}(end)+1 -40 5])
%             axis([0 Imax -7 5])
%             xlabel('Current, mA')
%             ylabel('Optical power, mW')
%             
%         end
        
%         figure(2)
%         subplot(235)
%         hold on
%         grid on
%         box on
%         
%         lenA=1:length(A.mode.ii_dd);
%         
%         plot(1000*A.mode.ii_dd,1e-12*A.mode.nMaxVet(lenA)/A.mode.NMQW,[colo(kcol),'o-'],...
%             1000*A.mode.ii_dd,1e-12*A.mode.pMaxVet(lenA)/A.mode.NMQW,[colo(kcol),'+--'])
%         hold on
%         plot(1000*A.mode.ii_dd,1e-18*A.mode.n3MaxVet(lenA),[colo(kcol),'-'],...
%             1000*A.mode.ii_dd,1e-18*A.mode.p3MaxVet(lenA),[colo(kcol),'--'],'linewidth',2)
%         
%         lenB=1:length(B.mode.ii_dd);
%         
%         plot(1000*B.mode.ii_dd,1e-12*B.mode.nMaxVet(lenB)/B.mode.NMQW,[colo(kcol),'o-'],...
%             1000*B.mode.ii_dd,1e-12*B.mode.pMaxVet(lenB)/B.mode.NMQW,[colo(kcol),'+--'])
%         hold on
%         plot(1000*B.mode.ii_dd,1e-18*B.mode.n3MaxVet(lenB),[colo(kcol),'-'],...
%             1000*B.mode.ii_dd,1e-18*B.mode.p3MaxVet(lenB),[colo(kcol),'--'],'linewidth',2)
%         
% %                 plot(1000*A.mode.ii_dd,1e-12*A.mode.nMaxVet(len)/A.mode.NMQW,'ko-',...
% %                     1000*A.mode.ii_dd,1e-12*A.mode.pMaxVet(len)/A.mode.NMQW,'k+--')
% %                 hold on
% %                 plot(1000*A.mode.ii_dd,1e-18*A.mode.n3MaxVet(len),'k-',...
% %                     1000*A.mode.ii_dd,1e-18*A.mode.p3MaxVet(len),'k--','linewidth',2)
%         
%         xlabel('Current (mA)')
%         ylabel('N-P density (PER well) 1e12/cm^2)')
%         legend('Electrons','Holes','location','best')
        
        pausak

    end
end

return
HeatSourcePlot(A.mesh,A.modeold,A.MODEplot{4})

HeatSourcePlot(B.mesh,B.modeold,B.MODEplot{4})