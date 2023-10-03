    indP=[1,Nmodi+1];
    StP=[St_11(indP,indP,1),St_12(indP,indP,1);St_21(indP,indP,1),St_22(indP,indP,1)];
    % disp('errore unitarieta  matrice totale')
    % 
    % max(max(abs(StP'*StP-eye(4))))
    if strcmp(flagPlot,'Plot_SI')
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %  diagrammi campo alle varie interfacce
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        Npuntid1=101;
        Npuntid2=101;
        x1=linspace(0,d1,Npuntid1);
        x2=linspace(d1,d1+d2,Npuntid2);
        x=[x1,x2];
        
        cpiu=zeros(NmodiTETM,Nstratper+2);
        cmeno=zeros(NmodiTETM,Nstratper+2);
        
        if strcmp(polariz_inc,'TE') 
            cpiu(1,1)=1; % Incidenza TE
        elseif strcmp(polariz_inc,'TM')
            cpiu(Nmodi+1,1)=1; % Incidenza TM
        end
        
        cmeno(:,1)=St_11(:,:,1)*cpiu(:,1);
        Hiniz=50; %nm lunghezza del tratto di n3 prima della prima giunzione 
        Niniz=20;
        %   Calcolo delle funzioni modali dei modi di Floquet    
        [eFloqTExAn, eFloqTEyAn, hFloqTExAn, hFloqTEyAn] = FloqzAnComp(x1, d1, x2, d2,phi,csi, ky, k0, n3, 'TE',segem);
        [eFloqTMxAn, eFloqTMyAn, hFloqTMxAn, hFloqTMyAn] = FloqzAnComp(x1, d1, x2, d2,phi,csi, ky, k0, n3, 'TM',segem);
        %    Calcolo dei campi trasversali a z in funzione di x e z nel tratto Hiniz        
        z0=linspace(-Hiniz,0,Niniz);
        for indz=1:Niniz
            Propagz=exp(-j*kz(:,1)*z0(indz));
            cpiuz=cpiu(:,1).*Propagz;
            cmenoz=cmeno(:,1)./Propagz;
            Vz=sqrt(1./ymodale(1,:)).*(cpiuz + cmenoz).';
            Iz=sqrt(ymodale(1,:)).*(cpiuz - cmenoz).';
            Ex_z(indz,:)=Vz*[eFloqTExAn;eFloqTMxAn]; 
            Ey_z(indz,:)=Vz*[eFloqTEyAn;eFloqTMyAn]; 
            Hx_z(indz,:)=Iz*[hFloqTExAn;hFloqTMxAn];
            Hy_z(indz,:)=Iz*[hFloqTEyAn;hFloqTMyAn];
        end  %%% indz=1:Nz
        nor=sqrt(d);        
        Ex_ztot=Ex_z*nor;
        Ey_ztot=Ey_z*nor;
        Hx_ztot=Hx_z*nor;
        Hy_ztot=Hy_z*nor;
        assez=z0;
        clear Ex_z Ey_z Hx_z Hy_z 
        
        for ind=1:Nstratper
            cpiu(:,ind+1)=inv(eye(NmodiTETM)-S22(:,:,ind)*St_11(:,:,ind+1))*S21(:,:,ind)*cpiu(:,ind);
            cmeno(:,ind+1)=St_11(:,:,ind+1)*cpiu(:,ind+1);
            cpiuL(:,ind+1)=cpiu(:,ind+1)./Propag(:,ind);
            cmenoL(:,ind+1)=cmeno(:,ind+1).*Propag(:,ind);
        end
        cpiuL(:,Nstratper+2)=St_21(:,:,1)*cpiu(:,1);
        cmenoL(:,Nstratper+2)=zeros(NmodiTETM,1);
        
        ziniz=0;
        thicktot=sum(thick);
        Nztot=151;
        Nz=ceil(Nztot*thick/thicktot);
        for ind=1:Nstratper+1 % ciclo sulle interfacce per verificare la continuita'
            Vhat=sqrt(1./ymodale(ind,:)).*(cpiu(:,ind) + cmeno(:,ind)).';
            Ihat=sqrt(ymodale(ind,:)).*(cpiu(:,ind) - cmeno(:,ind)).';
            Vtilda=sqrt(1./ymodale(ind+1,:)).*(cpiuL(:,ind+1) + cmenoL(:,ind+1)).';
            Itilda=sqrt(ymodale(ind+1,:)).*(cpiuL(:,ind+1) - cmenoL(:,ind+1)).';
            if ind==1
                efunxL=[eFloqTExAn;eFloqTMxAn];
                efunyL=[eFloqTEyAn;eFloqTMyAn];
                hfunxL=[hFloqTExAn;hFloqTMxAn];
                hfunyL=[hFloqTEyAn;hFloqTMyAn];
                %       Valutazione delle funzioni modali PSW                
                polariz='TE';
                kzG_TE2=kz(1:Nmodi,ind+1).^2;
                kTG_TE2=kT2(1:Nmodi,ind+1);
                [hfunTEx, hfunTEy, efunTEx, efunTEy] = modefunComp(x1, d1, n1(ind), x2, d2, n2(ind), VTE1(:,:,ind), VTE2(:,:,ind),...
                    ITE1(:,:,ind), ITE2(:,:,ind),...
                    kzG_TE2, kTG_TE2, ky, k0, polariz, Clight);
                polariz='TM';
                kzG_TM2=kz(Nmodi+1:NmodiTETM,ind+1).^2;
                kTG_TM2=kT2(Nmodi+1:NmodiTETM,ind+1);
                [hfunTMx, hfunTMy, efunTMx, efunTMy] = modefunComp(x1, d1, n1(ind), x2, d2, n2(ind), ITM1(:,:,ind), ITM2(:,:,ind),...
                    VTM1(:,:,ind), VTM2(:,:,ind),...
                    kzG_TM2, kTG_TM2, ky, k0, polariz, Clight);
                
                efunxR=[efunTEx;efunTMx];
                efunyR=[efunTEy;efunTMy];
                hfunxR=[hfunTEx;hfunTMx];
                hfunyR=[hfunTEy;hfunTMy];
                
                zfin=ziniz+thick(ind);
                
                z=linspace(ziniz,zfin,Nz(ind));
                assez=[assez,z];
                for indz=1:Nz(ind)
                    Propagz=exp(-j*kz(:,ind+1)*(z(indz)-ziniz));
                    cpiuz=cpiuL(:,ind+1).*Propagz;
                    cmenoz=cmenoL(:,ind+1)./Propagz;
                    Vz=sqrt(1./ymodale(ind+1,:)).*(cpiuz + cmenoz).';
                    Iz=sqrt(ymodale(ind+1,:)).*(cpiuz - cmenoz).';
                    Ex_z(indz,:)=Vz*efunxR*nor; 
                    Ey_z(indz,:)=Vz*efunyR*nor; 
                    Hx_z(indz,:)=Iz*hfunxR*nor;
                    Hy_z(indz,:)=Iz*hfunyR*nor;
                end  %%% indz=1:Nz
                ziniz=zfin;
                Ex_ztot=[Ex_ztot;Ex_z];
                Ey_ztot=[Ey_ztot;Ey_z];
                Hx_ztot=[Hx_ztot;Hx_z];
                Hy_ztot=[Hy_ztot;Hy_z];
                clear Ex_z Ey_z Hx_z Hy_z
            elseif ind==Nstratper+1  % ind==1 ultima interfaccia
                efunxL=efunxR; % L = left, R = right
                efunyL=efunyR;
                hfunxL=hfunxR;
                hfunyL=hfunyR;
                
                efunxR=[eFloqTExAn;eFloqTMxAn];
                efunyR=[eFloqTEyAn;eFloqTMyAn];
                hfunxR=[hFloqTExAn;hFloqTMxAn];
                hfunyR=[hFloqTEyAn;hFloqTMyAn];
            else  %% ind==1 interfacce interne, diverse dalla prima e ultima
                efunxL=efunxR; % L = left, R = right
                efunyL=efunyR;
                hfunxL=hfunxR;
                hfunyL=hfunyR;
                
                polariz='TE';
                kzG_TE2=kz(1:Nmodi,ind+1).^2;
                kTG_TE2=kT2(1:Nmodi,ind+1);
                [hfunTEx, hfunTEy, efunTEx, efunTEy] = modefunComp(x1, d1, n1(ind), x2, d2, n2(ind), VTE1(:,:,ind), VTE2(:,:,ind),...
                    ITE1(:,:,ind), ITE2(:,:,ind),...
                    kzG_TE2, kTG_TE2, ky, k0, polariz, Clight);
                polariz='TM';
                kzG_TM2=kz(Nmodi+1:NmodiTETM,ind+1).^2;
                kTG_TM2=kT2(Nmodi+1:NmodiTETM,ind+1);
                [hfunTMx, hfunTMy, efunTMx, efunTMy] = modefunComp(x1, d1, n1(ind), x2, d2, n2(ind), ITM1(:,:,ind), ITM2(:,:,ind),...
                    VTM1(:,:,ind), VTM2(:,:,ind),...
                    kzG_TM2, kTG_TM2, ky, k0, polariz, Clight);
                
                efunxR=[efunTEx;efunTMx];
                efunyR=[efunTEy;efunTMy];
                hfunxR=[hfunTEx;hfunTMx];
                hfunyR=[hfunTEy;hfunTMy];
                
                zfin=ziniz+thick(ind);
                z=linspace(ziniz,zfin,Nz(ind));
                assez=[assez,z];
                for indz=1:Nz(ind)
                    Propagz=exp(-j*kz(:,ind+1)*(z(indz)-ziniz));
                    cpiuz=cpiuL(:,ind+1).*Propagz;
                    cmenoz=cmenoL(:,ind+1)./Propagz;
                    Vz=sqrt(1./ymodale(ind+1,:)).*(cpiuz + cmenoz).';
                    Iz=sqrt(ymodale(ind+1,:)).*(cpiuz - cmenoz).';
                    Ex_z(indz,:)=Vz*efunxR*nor; 
                    Ey_z(indz,:)=Vz*efunyR*nor; 
                    Hx_z(indz,:)=Iz*hfunxR*nor;
                    Hy_z(indz,:)=Iz*hfunyR*nor;
                end % indz=1:Nz
                ziniz=zfin;
                Ex_ztot=[Ex_ztot;Ex_z];
                Ey_ztot=[Ey_ztot;Ey_z];
                Hx_ztot=[Hx_ztot;Hx_z];
                Hy_ztot=[Hy_ztot;Hy_z];
                clear Ex_z Ey_z Hx_z Hy_z
                
                
            end % if ind==1
            Exhat=Vhat*efunxL;
            Eyhat=Vhat*efunyL;
            Hxhat=Ihat*hfunxL;
            Hyhat=Ihat*hfunyL;
            
            Extilda=Vtilda*efunxR;
            Eytilda=Vtilda*efunyR;
            Hxtilda=Itilda*hfunxR;
            Hytilda=Itilda*hfunyR;
            %       disegno dei campi trasversali a z alle varie giunzioni, per la verifica della continuita'        
            h=figure(2*ind-1);
            set(h,'pos',[201   112   947   568            ])
            subplot(2,2,1)
%            plot(x,abs(Exhat))
%            if strcmp(polariz,'TE')
            plot(x,real(Exhat))
            grid on
            title(['interfaccia  ',num2str(ind),'  campo Ex'])
            hold on
%            plot(x,abs(Extilda),'r')
            plot(x,real(Extilda),'r')
            subplot(2,2,2)
%            plot(x,abs(Eyhat))
            plot(x,real(Eyhat))
            grid on
            title(['interfaccia  ',num2str(ind),'  campo Ey'])
            hold on
%            plot(x,abs(Eytilda),'r')
            plot(x,real(Eytilda),'r')
%            figure(2*ind)
            subplot(2,2,3)
%            plot(x,abs(Hxhat))
            plot(x,real(Hxhat))
            grid on
            title(['interfaccia  ',num2str(ind),'  campo Hx'])
            hold on
%            plot(x,abs(Hxtilda),'r')
            plot(x,real(Hxtilda),'r')
            subplot(2,2,4)
%            plot(x,abs(Hyhat))
            plot(x,real(Hyhat))
            grid on
            title(['interfaccia  ',num2str(ind),'  campo Hy'])
            hold on
%            plot(x,abs(Hytilda),'r')
            plot(x,real(Hytilda),'r')
            
pausak                        
            
        end % for ind=1:Nstratper+1
    end %if strcmp(flagPlot,'Plot_SI')
% disegno dei parametri S in funzione di lambda
ifig=0;
if length(lambdavet)>1
 ifig=1;
end 

if ifig==1
figure
subplot(2,2,1)
plot(lambdavet,abs(S11TETEtotale),'k')
grid on
xlabel('lambda (nm)')
title('abs(S11TETE)')
sizefig=axis;
axis([sizefig(1),sizefig(2),0,1])

subplot(2,2,2)
plot(lambdavet,abs(S11TETMtotale),'k')
grid on
xlabel('lambda (nm)')
title('abs(S11TETM)')
sizefig=axis;
axis([sizefig(1),sizefig(2),0,1])

subplot(2,2,3)
plot(lambdavet,abs(S11TMTEtotale),'k')
grid on
xlabel('lambda (nm)')
title('abs(S11TMTE)')
sizefig=axis;
axis([sizefig(1),sizefig(2),0,1])

subplot(2,2,4)
plot(lambdavet,abs(S11TMTMtotale),'k')
grid on
xlabel('lambda (nm)')
title('abs(S11TMTM)')
sizefig=axis;
axis([sizefig(1),sizefig(2),0,1])
% plot di S21
figure
subplot(2,2,1)
plot(lambdavet,abs(S21TETEtotale),'k')
grid on
xlabel('lambda (nm)')
title('abs(S21TETE)')
sizefig=axis;
axis([sizefig(1),sizefig(2),0,1])

subplot(2,2,2)
plot(lambdavet,abs(S21TETMtotale),'k')
grid on
xlabel('lambda (nm)')
title('abs(S21TETM)')
sizefig=axis;
axis([sizefig(1),sizefig(2),0,1])

subplot(2,2,3)
plot(lambdavet,abs(S21TMTEtotale),'k')
grid on
xlabel('lambda (nm)')
title('abs(S21TMTE)')
sizefig=axis;
axis([sizefig(1),sizefig(2),0,1])

subplot(2,2,4)
plot(lambdavet,abs(S21TMTMtotale),'k')
grid on
xlabel('lambda (nm)')
title('abs(S21TMTM)')
sizefig=axis;
axis([sizefig(1),sizefig(2),0,1])
end  % ifig


% disegno dei diagrammi 3D: abs(E_x) in funzione di x e z
if strcmp(flagPlot,'Plot_SI')
 i3D=0;
 if i3D==1

    figure
    surf(x,assez,abs(Ex_ztot))
    title('abs(Extrasv)')
    xlabel('x')
    ylabel('z')
    
    % disegno dei diagrammi 3D: abs(E_y) in funzione di x e z
    figure
    surf(x,assez,abs(Ey_ztot))
    title('abs(Eytrasv)')
    xlabel('x')
    ylabel('z')
    
    % disegno dei diagrammi 3D: abs(H_x) in funzione di x e z
    figure
    surf(x,assez,abs(Hx_ztot))
    title('abs(Hxtrasv)')
    xlabel('x')
    ylabel('z')
    
    % disegno dei diagrammi 3D: abs(H_y) in funzione di x e z
    figure
    surf(x,assez,abs(Hy_ztot))
    title('abs(Hytrasv)')
    xlabel('x')
    ylabel('z')
    
    % contour plots di abs(E)^2 nel piano x z
    figure
    contour(x,assez,abs(Ex_ztot).^2+abs(Ey_ztot).^2,50)
    hold on
    zlimite=cumsum(thick);
    plot([0,d],[0,0],'k')
    for ind=1:Nstratper
        plot([0,d],[zlimite(ind),zlimite(ind)],'k')
    end
    plot([d1,d1],[0,zlimite(end)],'k')
    title('abs(Etrasv)^2')
    xlabel('x')
    ylabel('z')
    
    % contour plots di abs(H)^2 nel piano x z
    figure
    contour(x,assez,abs(Hx_ztot).^2+abs(Hy_ztot).^2,50)
    hold on
    zlimite=cumsum(thick);
    plot([0,d],[0,0],'k')
    
    for ind=1:Nstratper
        plot([0,d],[zlimite(ind),zlimite(ind)],'k')
    end
    plot([d1,d1],[0,zlimite(end)],'k')
    title('abs(Htrasv).^2')
    xlabel('x')
    ylabel('z')
 else
  h=figure, 
  set(h,'pos',[185          58        1055         730])
  subplot(2,2,1)

  plot(assez,abs(mean(Ex_ztot,2))), grid
  title(' Norm. Ex')
    subplot(2,2,2)
    plot(assez,abs(mean(Ey_ztot,2))), grid
  title(' Norm. Ey')
  subplot(2,2,3)
  plot(assez,abs(mean(Hx_ztot,2))), grid
  title(' Norm. Hx')
  subplot(2,2,4)
  plot(assez,abs(mean(Hy_ztot,2))), grid
  title(' Norm. Hy')
 end
 ar=assez;
 Ex=(mean(Ex_ztot,2));
 Ey=(mean(Ey_ztot,2));
 Hx=(mean(Hx_ztot,2));
 Hy=(mean(Hy_ztot,2));
 save caret ar Ex Ey Hx Hy Ex_ztot Ey_ztot Hx_ztot Hy_ztot
'fine plots', keyboard    
end % if strcmp(flagPlot,'Plot_SI')

