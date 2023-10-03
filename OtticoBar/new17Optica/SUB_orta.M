 
global ico
%ico=0;   %controllo Pierluigi

% Versione ORTA
% Torino, 14 maggio 2010
% ---------------------------------------------------------------
% ***************************************************************
% ---------------------------------------------------------------
%      BETA OTTENUTI CON UNO SVILUPPO SU POLINOMI DI LEGENDRE
% ---------------------------------------------------------------
% ***************************************************************
% ---------------------------------------------------------------
% Programma completo utile a ricavare le efficienze di trasmissio-
% ne e riflessione relative al caso di un'onda piana incidente su
% un reticolo a profilo rettangolare e con regione stratificata
% composta solo da due strati a indice diverso (reticolo rettango-
% lare classico).
% E` possibile considerare anche il caso con perdite (quindi indi-
% ci di rifrazione complessi).
% Il metodo usato e' quello del MODE-MATCHING in cui le funzioni
% in cui e' stato espanso il campo sull'apertura sono i modi della
% regione PERIODICA.
%
% I dati della struttura devono essere inseriti nel file testo
% datret2g.m, che deve quindi essere presente nella directory
% di lavoro. Inoltre si puo' scegliere la polarizzazione (TE o TM)
% dell'onda piana incidente e il tipo di formulazione per ricavare
% la matrice scattering (HFIE o EFIE).
%
% Le altre funzioni necessarie sono:
% msqrt.m               per il calcolo della radice quadra di un numero
%                               in modo tale che la parte immaginaria sia
%                               negativa
% stimbet.m     per il calcolo delle costanti di propagazione
%                               e dei modi nella regione stratificata periodica
% MatPr.m               per il calcolo della matrice di proiezione tra
%                               i modi della regione uniforme e quella strati-
%                               ficata
% GSMgiunzU_S.m per il calcolo della matrice scattering della
%                               giunzione reg.unif.-reg.stratif.
% GSMgiunzS_U.m per il calcolo della matrice scattering della
%                               giunzione reg.stratif.-reg.unif.
% andmodp.m             per calcolare l'andamento dei modi della
%                               regione stratificata periodica
% calcpotf.m    per il calcolo della potenza complessa agli
%                               estremi della matrice scattering completa
% EHplot.m              per plottare i campi elettrico e magnetico in
%                               modulo e fase su una giunzione
%
% Il programma permette di fare i calcoli al variare della
% lunghezza d'onda dell'onda piana incidente e della profondita'
% del reticolo.
% ---------------------------------------------------------------
% ***************************************************************
% ---------------------------------------------------------------
%

beta=zeros(Nmodi,Nstratper+2);
ymodale=zeros(Nstratper+2,Nmodi);

pl=1; % serve per matrice S



flagPlot='Plot_SI';
flagPlot='Plot_NO';




k0=2*pi/lambda;

KBd = 2*pi*d/lambda*real(n3)*sin(theta);   	% sfasamento per cella
k = 1:(Nmodi-1)/2;
Csi(2*k) = -k;
Csi(2*k+1) = k;
csi = KBd/d + 2*pi*Csi/d;

clear k

kz3 = msqrt((k0*n3)^2 - csi.^2); 
kz4 = msqrt((k0*n4)^2 - csi.^2); 
beta(:,1)=kz3.';
beta(:,Nstratper+2)=kz4.';
if strcmp(polariz,'TM')
    ymodale(1,:) = k0*n3^2./(kz3*c*mu0); 
    ymodale(Nstratper+2,:)= k0*n4^2./(kz4*c*mu0);
else % polariz='TE'
    ymodale(1,:)  = kz3./(k0*c*mu0);
    ymodale(Nstratper+2,:) = kz4./(k0*c*mu0);
end
zmodale(1,:)=1./ymodale(1,:);
zmodale(Nstratper+2,:)=1./ymodale(Nstratper+2,:);

for ind=1:Nstratper
    [betadiel2, Mdie1, Mdie2, Vdie1, Vdie2, Vdie1H, Vdie2H] = ...
        modiLeg(n1(ind),n2(ind),d1,d2,n3,theta,polariz,lambda,Nmodi,errore);
        
    beta(:,ind+1)=msqrt(betadiel2);
    V1(:,:,ind)=Vdie1;
    V2(:,:,ind)=Vdie2;
    V1H(:,:,ind)=Vdie1H;
    V2H(:,:,ind)=Vdie2H;
    
    if strcmp(polariz,'TM')
        ymodale(ind+1,:) = k0./(beta(:,ind+1).'*c*mu0);
    else % polariz='TE'
        ymodale(ind+1,:) = (beta(:,ind+1)).'./(k0*c*mu0);
    end
    
    zmodale(ind+1,:) = 1./ymodale(ind+1,:); 
    
    if strcmp(formulaz,'HFIE')
        if ind==1 % giunzione U-1^Periodo
            TeFG = MatProiezPerUnifAnalit(-csi,n1(1),d1,V1(:,:,1),n2(1),d2,V2(:,:,1),polariz,'HFIE');
            ThGF = MatProiezPerUnifAnalit(csi,n1(1),d1,V1H(:,:,1),n2(1),d2,V2H(:,:,1),polariz,'HFIE');
            [S11(:,:,ind), S12(:,:,ind), S21(:,:,ind), S22(:,:,ind)] = ...
                GSMgiunzU_Sorta(TeFG, ThGF, ymodale(1,:), ymodale(ind+1,:), formulaz);
        else % tutte le giunzioni eccetto la prima e l'ultima
        % matrici di proiezione L/R
            TeLR = MatProiezPerRL(n1(ind-1),d1,V1H(:,:,ind-1),n2(ind-1),d2,V2H(:,:,ind-1),...
            n1(ind),V1(:,:,ind),n2(ind),V2(:,:,ind),polariz);
		
            ThRL =  MatProiezPerRL(n1(ind-1),d1,V1(:,:,ind-1),n2(ind-1),d2,V2(:,:,ind-1),...
            n1(ind),V1H(:,:,ind),n2(ind),V2H(:,:,ind),polariz);

            [S11(:,:,ind), S12(:,:,ind), S21(:,:,ind), S22(:,:,ind)] = ...
                GSMgiunzorta(TeLR,ThRL, ymodale(ind,:), ymodale(ind+1,:), formulaz);
        end  %ind==1 % giunzione U-1^Periodo
        
    else % formulaz='EFIE'
        % da fare
        %[S111, S112, S121, S122] = GSMgiunzU_S(Teab, Thba, zFloq3, zmodale, formulaz);
        %[S2110, S2120, S2210, S2220] = GSMgiunzS_U(Teab, Thba, zFloq4, zmodale, formulaz);

    end  %if strcmp(formulaz,'HFIE')

end % for ind=1:Nstratper

% giunzione di destra, tra n1(end)/n2(end) e n4 uniforme
ThFG = MatProiezPerUnifAnalit(-csi,n1(Nstratper),d1,V1(:,:,Nstratper),n2(Nstratper),d2,V2(:,:,Nstratper),polariz,'EFIE');
TeGF = MatProiezPerUnifAnalit(csi,n1(Nstratper),d1,V1H(:,:,Nstratper),n2(Nstratper),d2,V2H(:,:,Nstratper),polariz,'EFIE');
[S11(:,:,Nstratper+1), S12(:,:,Nstratper+1), S21(:,:,Nstratper+1), S22(:,:,Nstratper+1)] = ...
    GSMgiunzS_Uorta(TeGF, ThFG, ymodale(Nstratper+2,:), ymodale(Nstratper+1,:), formulaz);

% spostiamo i piani della porta 2 di tutte le matrici S, eccetto l'ultima

for k=1:Nstratper
    Propag(:,k)=exp(-j*beta(:,k+1)*thick(k));
    matdiag = (Propag(:,k)*ones(1,Nmodi));
    S12(:,:,k)=S12(:,:,k).*matdiag.';
    S21(:,:,k)=matdiag.*S21(:,:,k);
    S22(:,:,k)=(Propag(:,k)*Propag(:,k).').*S22(:,:,k);
end
bet=beta(:,2:end-1)*1000;
%' propag', keyboard
%CASCATA PER DEFINIRE UN'UNICA MATRICE SCATTERING DI TUTTA LA STRUTTURA
St_11(:,:,Nstratper+1)=S11(:,:,Nstratper+1);
St_12(:,:,Nstratper+1)=S12(:,:,Nstratper+1);
St_21(:,:,Nstratper+1)=S21(:,:,Nstratper+1);
St_22(:,:,Nstratper+1)=S22(:,:,Nstratper+1);

for k=(Nstratper):-1:1
        St_11(:,:,k)=S11(:,:,k)+S12(:,:,k)*St_11(:,:,k+1)*(inv(eye(Nmodi)-S22(:,:,k)*St_11(:,:,k+1)))*S21(:,:,k);
        St_12(:,:,k)=S12(:,:,k)*(inv(eye(Nmodi)-St_11(:,:,k+1)*S22(:,:,k)))*St_12(:,:,k+1);
        St_21(:,:,k)=St_21(:,:,k+1)*(inv(eye(Nmodi)-S22(:,:,k)*St_11(:,:,k+1)))*S21(:,:,k);
        St_22(:,:,k)=St_22(:,:,k+1)+St_21(:,:,k+1)*(inv(eye(Nmodi)-S22(:,:,k)*St_11(:,:,k+1)))*S22(:,:,k)*St_12(:,:,k+1); 
end

% controllo unitarieta per i modi sopra taglio

% coefficienti di riflessione e trasmissione per il modo fondamentale
s11=St_11(1,1,1);
s21=St_21(1,1,1);
s12=St_12(1,1,1);
s22=St_22(1,1,1);



%figure, semilogy(csic,abs(S11p),'.',csic,abs(S11d),'.'), keyboard
        
        St(1,1)=s11;
        St(1,2)=s12;
        St(2,1)=s21;
        St(2,2)=s22;
   	icarico=0;
        if exist('Strut') 
         if length(Strut)>0
     	  icarico=Strut.icarico;       
  	 end
        end
%        'controllo sub_orta', keyboard
       if icarico==1
	L_i=Strut.t;
	n_i=Strut.n;
	ibast=Strut.ibast;

	Lp=flipud(L_i(2:ibast(1)-1));
	np=flipud(n_i(2:ibast(1)-1));
	np_out=n_i(1);
	np_in=n_i(ibast(1));

	Ld=(L_i(ibast(end)+1:end-1));
	nd=(n_i(ibast(end)+1:end-1));
	nd_out=n_i(end);
	nd_in=n_i(ibast(end));
	
	csic=csi(2:end);
	
	' in sub_ort, prima di Gamma ', keyboard
	[S11p]=Gam_ort(csic,np_in,np,np_out,Lp,lambda,polariz);
	[S11d]=Gam_ort(csic,nd_in,nd,nd_out,Ld,lambda,polariz);

	A11=St_11(:,:,1);
	A12=St_12(:,:,1);
	A21=St_21(:,:,1);
	A22=St_22(:,:,1);
	
	pu=2:length(A11);
	
	A11d=A11(pu,pu);
	A12d=A12(pu,pu);
	A21d=A21(pu,pu);
	A22d=A22(pu,pu);
	C22=[A11d A12d; A21d A22d];
	C11=[A11(1,1) A12(1,1); A21(1,1) A22(1,1)]; 
	C12=[A11(1,pu) A12(1,pu); A21(1,pu) A22(1,pu)]; 
	C21=[[A11(pu,1); A21(pu,1)] [A12(pu,1); A22(pu,1)]]; 
	
	pu=2:length(A11);

	A11d=A11(pu,pu);
	A12d=A12(pu,pu);
	A21d=A21(pu,pu);
	A22d=A22(pu,pu);
	C22=[A11d A12d; A21d A22d];
	C11=[A11(1,1) A12(1,1); A21(1,1) A22(1,1)]; 
	C12=[A11(1,pu) A12(1,pu); A21(1,pu) A22(1,pu)]; 
	C21=[[A11(pu,1); A21(pu,1)]  [A12(pu,1); A22(pu,1)]]; 
	
	Ga=diag([S11p S11d]);

	I12=inv(eye(length(Ga))-Ga*C22);

	Tcor=C12*I12*Ga*C21;
        Stc=C11+Tcor;       
        s11=Stc(1,1);
        s12=Stc(1,2);
        s21=Stc(2,1);
        s22=Stc(2,2);
        ' passo carico sub_orta', keyboard
       end

        sinv=inv(s12);        
	Tr(1,1)=s21-s22*sinv*s11;
	Tr(1,2)=s22*sinv;
	Tr(2,1)=-sinv*s11;
	Tr(2,2)=sinv;

%'fine Orta', keyboard
%ifp=-10
if ~exist('ifp')
 ifp=-4;
end 

if ifp==-10
  if strcmp(polariz,'TE')==1
   c1='r.';
   c2='c';
  else
   c1='b.';
   c2='g';
  end
%  'campi', keyboard
    %diagrammi campo alle varie interfacce
    Npuntid1=51;
    Npuntid2=51;
    x1=linspace(0,d1,Npuntid1);
    x2=linspace(d1,d1+d2,Npuntid2);
    x=[x1,x2];
    
    cpiu=zeros(Nmodi,Nstratper+2);
    cmeno=zeros(Nmodi,Nstratper+2);
    cpiu(:,1)=eye(Nmodi,1);
    cmeno(:,1)=St_11(:,:,1)*cpiu(:,1);
    % z0=linspace(-Hiniz,0,250);
    % cpiuz=[exp(-j*kzs(1)*z0)*cpiu(1,1);exp(-j*kzp(1)*z0)*cpiu(2,1)];
    % cmenoz=[exp(j*kzs(1)*z0)*cmeno(1,1);exp(j*kzp(1)*z0)*cmeno(2,1)];
    % psiz=M(:,:,1)*[cpiuz;cmenoz];
    % ziniz=0;
    % assez=0;
    %psitot=psiz;
    for ind=1:Nstratper
        cpiu(:,ind+1)=inv(eye(Nmodi)-S22(:,:,ind)*St_11(:,:,ind+1))*S21(:,:,ind)*cpiu(:,ind);
        cmeno(:,ind+1)=St_11(:,:,ind+1)*cpiu(:,ind+1);
        %     zfin=ziniz+thick(ind);
        %     z=linspace(ziniz,zfin,2);
        %     assez=[assez,z];
        cpiuL(:,ind+1)=cpiu(:,ind+1)./Propag(:,ind);
        cmenoL(:,ind+1)=cmeno(:,ind+1).*Propag(:,ind);
        %     psiz=M(:,:,ind+1)*[cpiuz;cmenoz];
        %     psitot=[psitot,psiz];
        %     ziniz=zfin;
    end
    cpiuL(:,Nstratper+2)=St_21(:,:,1)*cpiu(:,1);
    cmenoL(:,Nstratper+2)=zeros(Nmodi,1);
    
    % zfin=ziniz+H(N-1);
    % z=linspace(ziniz,zfin,50);
    %     assez=[assez,z];
    %     cpiuz=[exp(-j*kzs(N)*(z-ziniz))*cpiu(1,N);exp(-j*kzp(N)*(z-ziniz))*cpiu(2,N)];
    %     psiz=M(:,:,N)*[cpiuz;zeros(2,50)];
    %     psitot=[psitot,psiz];
    
    for ind=1:Nstratper+1
        Vhat=sqrt(zmodale(ind,:)).*(cpiu(:,ind) + cmeno(:,ind)).';
        Ihat=sqrt(ymodale(ind,:)).*(cpiu(:,ind) - cmeno(:,ind)).';
        Vtilda=sqrt(zmodale(ind+1,:)).*(cpiuL(:,ind+1) + cmenoL(:,ind+1)).';
        Itilda=sqrt(ymodale(ind+1,:)).*(cpiuL(:,ind+1) - cmenoL(:,ind+1)).';
        if ind==1
            % valutazione delle funzioni modali di Floquet (analitiche) 
            [hfunL, efunL] = FloqAnComp(x1, d1, x2, d2,csi, polariz);
            [hfunR, efunR] = modefunComp(x1, d1, n1(ind), x2, d2, n2(ind), V1(:,:,ind), V2(:,:,ind), polariz);
        elseif ind==Nstratper+1
            [hfunL, efunL] = modefunComp(x1, d1, n1(ind-1), x2, d2, n2(ind-1), V1(:,:,ind-1), V2(:,:,ind-1), polariz);
            [hfunR, efunR] = FloqAnComp(x1, d1, x2, d2,csi, polariz);
        else
            [hfunL, efunL] = modefunComp(x1, d1, n1(ind-1), x2, d2, n2(ind-1), V1(:,:,ind-1), V2(:,:,ind-1), polariz);
            [hfunR, efunR] = modefunComp(x1, d1, n1(ind), x2, d2, n2(ind), V1(:,:,ind), V2(:,:,ind), polariz);
        end
        Ehat=Vhat*efunL;
        Hhat=Ihat*hfunL;
        
        Etilda=Vtilda*efunR;
        Htilda=Itilda*hfunR;
        
        figure(2*ind-1)
        set(gcf,'pos',[63   560-(ind-1)*400   672   380])
        plot(x,abs(Ehat),c2)
        grid on
        title(['Interface  ',num2str(ind),'   |E|; dots in grating '])
        hold on
        plot(x,abs(Etilda),c1)
        xlabel(' Transvese coord (micron)')
        
        figure(2*ind)
        set(gcf,'pos',[ 824   560-(ind-1)*400    672   380])

        plot(x,abs(Hhat),c2)
        grid on
        title(['Interface  ',num2str(ind),'   |H|; dots in grating '])
        xlabel(' Transvese coord (micron)')

        hold on
        plot(x,abs(Htilda),c1)
        
    end % for ind=1:Nstratper+1

pausak

end  %ifp
        

