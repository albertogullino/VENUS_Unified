%'entro ML', keyboard
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
%
% Programma per l'analisi di un reticolo a rilievo di supeficie, per
% incidenza arbitraria. Versione base, con grafici dei campi
%
%                   R. Orta Dicembre 2010
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
%
%
%
%
%
% Il semispazio di sinistra ha indice n3 (reale), il
% semispazio di destra ha indice n4 (reale). 
% I denti hanno indice n1(ind), n2(ind), ind=1:Nstratper e possono essere complessi. 
% L'analisi e' effettuata con il mode matching su Nmodi modi TE e altrettanti TM, sempre
% dispari. I modi nella regione periodica sono trovati con il metodo degli
% elementi spettrali, basato su polinomi di Legendre. Il numero di polinomi
% e' determinato automaticamente per poter calcolare accuratamente il
% numero Nmodi di modi. 
% Il programma calcola la matrice scattering per
% tutti gli Nmodi TE (Ez=0) e TM (Hz=0), quindi ha dimensione 4*Nmodi. 
% Le porte sono ordinate nel modo seguente: da 1 a Nmodi TE0,TE-1, TE1,
% TE-2, TE2..... Poi da Nmodi+1 a 2*Nmodi TM0,TM-1,TM1,TM-2, TM2....
% Il flag 'PlotSI' abilita il disegno del campo elettrico e magnetico
% trasversale a ogni interfaccia, con la possibilita' di scegliere la polarizzazione del campo incidente.
% Poi disegna il grafico nel piano (x,z) delle
% componenti Ex, Ey, Hx, Hy del campo.
% Conviene usare questa opzione per una sola lambda. A questo scopo si fissa la lunghezza d'onda di interessa come lambda min, 
% Nlambda=2 e passo=4. 
% Il flag 'PlotNO' fa disegnare solo il grafico del moduli dei paramentri
% S11 e S21 per l'ordine 0.






S11TETEtotale=zeros(1,Nlambda); 
S11TETMtotale=zeros(1,Nlambda); 
S11TMTEtotale=zeros(1,Nlambda); 
S11TMTMtotale=zeros(1,Nlambda);

S21TETEtotale=zeros(1,Nlambda); 
S21TETMtotale=zeros(1,Nlambda);
S21TMTEtotale=zeros(1,Nlambda); 
S21TMTMtotale=zeros(1,Nlambda); 

S12TETEtotale=zeros(1,Nlambda); 
S12TETMtotale=zeros(1,Nlambda); 
S12TMTEtotale=zeros(1,Nlambda); 
S12TMTMtotale=zeros(1,Nlambda);

S22TETEtotale=zeros(1,Nlambda); 
S22TETMtotale=zeros(1,Nlambda);
S22TMTEtotale=zeros(1,Nlambda); 
S22TMTMtotale=zeros(1,Nlambda); 

kz=zeros(NmodiTETM,Nstratper+2); % costanti di propagazione nei vari modi
ymodale=zeros(Nstratper+2,NmodiTETM); % ammettenze modali
kT2=zeros(NmodiTETM,Nstratper+2); % k trasversale a z per i modi di Floquet (prima e ultima colonna)
%                                   k trasversale a (a x) per i modi della
%                                   PSW nelle colonne da 2 a Nstratper+1



n3=r_in(end);
n_in=r_in(1);
n4=r_out(1);
n_out=r_out(end);
n3ML=r_in;
n4ML=r_out;

phi=phi_in;
theta=asin(real(n_in/n3)*sin(theta_in));
%'teta', keyboard
%for indlambda=1:passo:Nlambda
    lambda=lambdavet;
    
    k0=2*pi/lambda;
    omega=k0*Clight;
    ky= k0*n3*sin(theta)*sin(phi);
    KBd = 2*pi*d/lambda*n3*sin(theta)*cos(phi);   	% sfasamento per cella
    
%    'cont teta', keyboard
    % Modi di Floquet nei semispazi
    k = 1:(Nmodi-1)/2;
    Csi(2*k) = -k;
    Csi(2*k+1) = k;
    csi = KBd/d + 2*pi*Csi/d;
    clear k
    ktF=psqrt(csi.^2+ky^2);
    kz3 = psqrt((k0*n3)^2 - ktF.^2); 
    kz30=k0*n3;
    fi3=length(find(kz3==0));
    %kz3(fi)=1e-6;
    kz4 = psqrt((k0*n4)^2 - ktF.^2); 
    fi4=length(find(kz4==0));
    %kz4(fi)=1e-6;
    if fi3>0 | fi4>0
     lambda=lambdavet+1e-3;
     k0=2*pi/lambda;
     omega=k0*Clight;
     ky= k0*n3*sin(theta)*sin(phi);
     KBd = 2*pi*d/lambda*n3*sin(theta)*cos(phi);   	% sfasamento per cella
     csi = KBd/d + 2*pi*Csi/d;
     clear k
     ktF=psqrt(csi.^2+ky^2);
     kz3 = psqrt((k0*n3)^2 - ktF.^2); 
    %kz3(fi)=1e-6;
     kz4 = psqrt((k0*n4)^2 - ktF.^2);         
    end
    

    %'qui', keyboard
    % caricamento dei parametri dei modi di Floquet
    kz(1:Nmodi,1)=kz3.';
    kz(Nmodi+1:NmodiTETM,1)=kz3.';
    kT2(1:Nmodi,1)=(ktF.^2).';
    kT2meta(1:Nmodi,1)=(ktF.^2).';
    kT2(Nmodi+1:NmodiTETM,1)=(ktF.^2).';
    kz(1:Nmodi,Nstratper+2)=kz4.';
    kz(Nmodi+1:NmodiTETM,Nstratper+2)=kz4.';
    kT2(1:Nmodi,Nstratper+2)=(ktF.^2).';
    kT2(Nmodi+1:NmodiTETM,Nstratper+2)=(ktF.^2).';
    ymodale(1,1:Nmodi)  = kz3./(k0*Clight*mu0); %yTE in n3
    ymodale(Nstratper+2,1:Nmodi) = kz4./(k0*Clight*mu0); %yTE in n4
    ymodale(1,Nmodi+1:NmodiTETM) = k0*n3^2./(kz3*Clight*mu0); % yTM in n3
    ymodale(Nstratper+2,Nmodi+1:NmodiTETM)= k0*n4^2./(kz4*Clight*mu0); %yTM in n4
    
    
    errorekz=1.e-6; % errore massimo sulle costanti di propagazione dei modi della guida PSW
    % loop sugli strati che compongono i denti 
    for ind=1:Nstratper
        % determinazione dei modi delle varie regioni PSW
        [kzG_TE2,kTG_TE2,ZG_TE,kzG_TM2,kTG_TM2,ZG_TM, M1, M2,VTE1loc, VTE2loc, VTE1Hloc, VTE2Hloc, ITE1loc, ITE2loc, ITE1Hloc, ITE2Hloc ,...
                VTM1loc, VTM2loc, VTM1Hloc, VTM2Hloc,ITM1loc, ITM2loc, ITM1Hloc, ITM2Hloc] = ...
            ModiPeriodLeg_s(n1(ind),n2(ind),d1,d2,n3,theta,phi,lambda,Nmodi,errorekz);
        kzG_TE=psqrt(kzG_TE2);
        kzG_TM=psqrt(kzG_TM2);
        % caricamento dei parametri dei modi della PSW    
        kz(1:Nmodi,ind+1)=kzG_TE;
        kz(Nmodi+1:NmodiTETM,ind+1)=kzG_TM;
        kT2(1:Nmodi,ind+1)=kTG_TE2;
        kT2(Nmodi+1:NmodiTETM,ind+1)=kTG_TM2;
        VTE1(:,:,ind)=VTE1loc;
        VTE2(:,:,ind)=VTE2loc;
        VTE1H(:,:,ind)=VTE1Hloc;
        VTE2H(:,:,ind)=VTE2Hloc;
        ITE1(:,:,ind)=ITE1loc;
        ITE2(:,:,ind)=ITE2loc;
        ITE1H(:,:,ind)=ITE1Hloc;
        ITE2H(:,:,ind)=ITE2Hloc;
        
        VTM1(:,:,ind)=VTM1loc;
        VTM2(:,:,ind)=VTM2loc;
        VTM1H(:,:,ind)=VTM1Hloc;
        VTM2H(:,:,ind)=VTM2Hloc;
        ITM1(:,:,ind)=ITM1loc;
        ITM2(:,:,ind)=ITM2loc;
        ITM1H(:,:,ind)=ITM1Hloc;
        ITM2H(:,:,ind)=ITM2Hloc;
        ymodale(ind+1,1:Nmodi) = (1./ZG_TE).'; %TE
        ymodale(ind+1,Nmodi+1:NmodiTETM) = (1./ZG_TM).'; %TM
        %   determinazione della matrice S delle varie giunzioni, partendo dalla
        %   piu' a sinistra. Il loop e' sulle interfacce
        if ind==1 % giunzione semispazio-1^Periodo
            ZG_TE=(1./ymodale(ind+1,1:Nmodi)).';
            ZG_TM=(1./ymodale(ind+1,Nmodi+1:NmodiTETM)).';
            %   Calcolo matrici di proiezione tra modi di Floquet e modi PSW alla prima giunzione       
            [TeLR,ThRL] = MatProiezPerUnifAnalitSkewLeft(phi,csi,ky,n1(1),d1,n2(1),d2,VTE1(:,:,ind), VTE2(:,:,ind), VTE1H(:,:,ind), VTE2H(:,:,ind),...
                ITE1(:,:,ind), ITE2(:,:,ind), ITE1H(:,:,ind), ITE2H(:,:,ind) ,kz(1:Nmodi,ind+1),ZG_TE,...
                VTM1(:,:,ind), VTM2(:,:,ind), VTM1H(:,:,ind), VTM2H(:,:,ind),ITM1(:,:,ind), ITM2(:,:,ind),...
                ITM1H(:,:,ind), ITM2H(:,:,ind),kz(Nmodi+1:NmodiTETM,ind+1),ZG_TM,segem);
        else % tutte le giunzioni eccetto la prima e l'ultima
            %   Calcolo matrici di proiezione tra modi PSW a sinistra e destra di ciascuna giunzione       
            ZL_TE=(1./ymodale(ind,1:Nmodi)).';
            ZR_TM=(1./ymodale(ind+1,Nmodi+1:NmodiTETM)).';
            [TeLR,ThRL]=MatProiezPerRLskew(ky,n1(ind),n2(ind),d1,d2,VTE1(:,:,ind-1),VTE2(:,:,ind-1),VTE1H(:,:,ind-1),VTE2H(:,:,ind-1),...
                ITE1(:,:,ind-1),ITE2(:,:,ind-1),ITE1H(:,:,ind-1),ITE2H(:,:,ind-1),ITM1(:,:,ind-1),ITM2(:,:,ind-1),ITM1H(:,:,ind-1),ITM2H(:,:,ind-1),...
                VTE1(:,:,ind),VTE2(:,:,ind),VTE1H(:,:,ind),VTE2H(:,:,ind),VTM1(:,:,ind),VTM2(:,:,ind),VTM1H(:,:,ind),VTM2H(:,:,ind),...
                ITM1(:,:,ind),ITM2(:,:,ind),ITM1H(:,:,ind),ITM2H(:,:,ind),kz(1:Nmodi,ind),ZL_TE,kz(Nmodi+1:NmodiTETM,ind+1),ZR_TM);
        end  %ind==1 % giunzione U-1^Periodo
        
        % Matrice scattering della giunzione ind        
        [S11(:,:,ind), S12(:,:,ind), S21(:,:,ind), S22(:,:,ind)] = ...
            GSMgiunzskew(TeLR,ThRL, ymodale(ind,:), ymodale(ind+1,:));
        
        %      ind1=[1 2 3 32 33 34];
        %      ind2=[1 2 32 33];
        %         ind1=find(abs(imag(kz(:,ind)))<=1e-10);
        %         ind2=find(abs(imag(kz(:,ind+1)))<=1e-10);
        %      S=[S11(ind1,ind1,ind),S12(ind1,ind2,ind);...
        %          S21(ind2,ind1,ind),S22(ind2,ind2,ind)];
        %   [m,n]=size(S);
        % disp(['errore unitarieta  giunzione ',num2str(ind)])
        %  max(max(abs(S'*S-eye(m))))
        
        % %         N1=3
        %         N2=2
        %        S12ref=S12([1:N1,Nmodi+1:Nmodi+N1],[1:N2,Nmodi+1:Nmodi+N2]);
        %        S22ref=S22([1:N2,Nmodi+1:Nmodi+N2],[1:N2,Nmodi+1:Nmodi+N2]);
        %        S21ref=S21([1:N2,Nmodi+1:Nmodi+N2],[1:N1,Nmodi+1:Nmodi+N1]);
        %        S11ref=S11([1:N1,Nmodi+1:Nmodi+N1],[1:N1,Nmodi+1:Nmodi+N1]);
        %     [S11(:,:,ind), S12(:,:,ind), S21(:,:,ind), S22(:,:,ind)]  = ...
        %         GSMgiunzskew(TeLR,ThRL,  ymodale(ind,:), ymodale(ind+1,:),NmodiAcc1,NmodiAcc2)
        %     NmodiAcc1=NmodiAcc2;
    end % for ind=1:Nstratper
    
    % giunzione di destra, tra n1(end)/n2(end) e n4 uniforme
    ZG_TE=(1./ymodale(Nstratper+1,1:Nmodi)).';
    ZG_TM=(1./ymodale(Nstratper+1,Nmodi+1:NmodiTETM)).';
    
    [TeGF, ThFG] = MatProiezPerUnifAnalitSkewRight(phi,csi,ky,n1(Nstratper),d1,n2(Nstratper),d2,VTE1(:,:,Nstratper), VTE2(:,:,Nstratper),...
        VTE1H(:,:,Nstratper), VTE2H(:,:,Nstratper),...
        ITE1(:,:,Nstratper), ITE2(:,:,Nstratper), ITE1H(:,:,Nstratper), ITE2H(:,:,Nstratper) ,kz(1:Nmodi,Nstratper+1),ZG_TE,...
        VTM1(:,:,Nstratper), VTM2(:,:,Nstratper), VTM1H(:,:,Nstratper), VTM2H,ITM1(:,:,Nstratper), ITM2(:,:,Nstratper),...
        ITM1H(:,:,Nstratper), ITM2H(:,:,Nstratper),kz(Nmodi+1:NmodiTETM,Nstratper+1),ZG_TM,segem);
    
    
    [S11(:,:,Nstratper+1), S12(:,:,Nstratper+1), S21(:,:,Nstratper+1), S22(:,:,Nstratper+1)] = ...
        GSMgiunzskew(TeGF,ThFG, ymodale(Nstratper+1,:), ymodale(Nstratper+2,:));
    
    %     ind1=find(abs(imag(kz(:,Nstratper+1)))<=1e-10);
    %         ind2=find(abs(imag(kz(:,Nstratper+2)))<=1e-10);
    %      S=[S11(ind1,ind1,Nstratper+1),S12(ind1,ind2,Nstratper+1);...
    %          S21(ind2,ind1,Nstratper+1),S22(ind2,ind2,Nstratper+1)];
    %   [m,n]=size(S);
    % disp(['errore unitarieta  giunzione ',num2str(Nstratper+1)])
    %  max(max(abs(S'*S-eye(m))))
    
%    'S11 prima', keyboard
    % Spostamento dei piani di riferimento della porta 2 di tutte le matrici S, eccetto l'ultima
    for k=1:Nstratper
        Propag(:,k)=exp(-j*kz(:,k+1)*thick(k));
        matdiag = (Propag(:,k)*ones(1,NmodiTETM));
        S12(:,:,k)=S12(:,:,k).*matdiag.';
        S21(:,:,k)=matdiag.*S21(:,:,k);
        S22(:,:,k)=(Propag(:,k)*Propag(:,k).').*S22(:,:,k);
    end
    
    %disp('prova alberto')
    %pausak

    
    % CASCATA PER DEFINIRE UN'UNICA MATRICE SCATTERING DI TUTTA LA STRUTTURA
    Str_11(:,:,Nstratper+1)=S11(:,:,Nstratper+1);
    Str_12(:,:,Nstratper+1)=S12(:,:,Nstratper+1);
    Str_21(:,:,Nstratper+1)=S21(:,:,Nstratper+1);
    Str_22(:,:,Nstratper+1)=S22(:,:,Nstratper+1);
    % St_ij(:,:,ind) e' la matrice scattering della struttura compresa tra
    % giunzione ind e Nstratper+1 (l'ultima)
    for k=(Nstratper):-1:1
        Str_11(:,:,k)=S11(:,:,k)+S12(:,:,k)*(inv(eye(NmodiTETM)-Str_11(:,:,k+1)*S22(:,:,k)))*Str_11(:,:,k+1)*S21(:,:,k);
        Str_12(:,:,k)=S12(:,:,k)*(inv(eye(NmodiTETM)-Str_11(:,:,k+1)*S22(:,:,k)))*Str_12(:,:,k+1);
        Str_21(:,:,k)=Str_21(:,:,k+1)*(inv(eye(NmodiTETM)-S22(:,:,k)*Str_11(:,:,k+1)))*S21(:,:,k);
        Str_22(:,:,k)=Str_22(:,:,k+1)+Str_21(:,:,k+1)*(inv(eye(NmodiTETM)-S22(:,:,k)*Str_11(:,:,k+1)))*S22(:,:,k)*Str_12(:,:,k+1); 
    end
   
     Stot0=[Str_11(1,1,1),Str_11(1,Nmodi+1,1),Str_12(1,1,1),Str_12(1,Nmodi+1,1);...
         Str_11(Nmodi+1,1,1),Str_11(Nmodi+1,Nmodi+1,1),Str_12(Nmodi+1,1,1),Str_12(Nmodi+1,Nmodi+1,1);...
         Str_21(1,1,1),Str_21(1,Nmodi+1,1),Str_22(1,1,1),Str_22(1,Nmodi+1,1);...
         Str_21(Nmodi+1,1,1),Str_21(Nmodi+1,Nmodi+1,1),Str_22(Nmodi+1,1,1),Str_22(Nmodi+1,Nmodi+1,1)];
         
         
%'dopo Str nuovo', keyboard
%%% n_in, theta_in, phi_in
%% n3ML(1:Nin), n4ML(1:Nout)
%% d_in,d_out, n_out
    SIN11=zeros(1:NmodiTETM,1:NmodiTETM);
    SIN12=SIN11;
    SIN21=SIN11;
    SIN22=SIN11;

    Nlayer_in=length(n3ML)-2;
    Nlayer_out=length(n4ML)-2;
    
    
    
    for kk=1:Nmodi
     beta2=kT2meta(kk);
     
     %multistrato ingresso
     %'fc_strat', keyboard
     [s11e,s12e,s21e,s22e]=fc_strat(beta2,Nlayer_in,n_in,n3ML(2:end-1),n3,d_in(1:end-1),lambda,'TE');
     [s11m,s12m,s21m,s22m]=fc_strat(beta2,Nlayer_in,n_in,n3ML(2:end-1),n3,d_in(1:end-1),lambda,'TM');
     if Nlayer_in>=0
      dexp=d_in(end);
     else
      dexp=0;
     end
     if kk==1
      Ds11=s11e;
      Ds22=s22e;
      Ds21=s21e;
     end
     SIN11(kk,kk)=s11e;
     SIN11(kk+Nmodi,kk+Nmodi)=s11m;
     SIN12(kk,kk)=s12e*exp(-j*kz3(kk)*dexp);
     SIN12(kk+Nmodi,kk+Nmodi)=s12m*exp(-j*kz3(kk)*dexp);
     SIN21(kk,kk)=s21e*exp(-j*kz3(kk)*dexp);
     SIN21(kk+Nmodi,kk+Nmodi)=s21m*exp(-j*kz3(kk)*dexp);
     SIN22(kk,kk)=s22e*exp(-j*2*kz3(kk)*dexp);
     SIN22(kk+Nmodi,kk+Nmodi)=s22m*exp(-j*2*kz3(kk)*dexp);    
  %   ' SIN keyboard', keyboard
     %multistrato uscita
     [s11e,s12e,s21e,s22e]=fc_strat(beta2,Nlayer_out,n4,n4ML(2:end-1),n_out,d_out(2:end),lambda,'TE');
     [s11m,s12m,s21m,s22m]=fc_strat(beta2,Nlayer_out,n4,n4ML(2:end-1),n_out,d_out(2:end),lambda,'TM');
     if Nlayer_out>0
      dexp=d_out(1);
     else
      dexp=0;
     end
     %'tetm', keyboard
     SOUT11(kk,kk)=s11e*exp(-j*2*kz4(kk)*dexp);
     SOUT11(kk+Nmodi,kk+Nmodi)=s11m*exp(-j*2*kz4(kk)*dexp);
     SOUT12(kk,kk)=s12e*exp(-j*kz4(kk)*dexp);
     SOUT12(kk+Nmodi,kk+Nmodi)=s12m*exp(-j*kz4(kk)*dexp);
     SOUT21(kk,kk)=s21e*exp(-j*kz4(kk)*dexp);
     SOUT21(kk+Nmodi,kk+Nmodi)=s21m*exp(-j*kz4(kk)*dexp);
     SOUT22(kk,kk)=s22e;
     SOUT22(kk+Nmodi,kk+Nmodi)=s22m;
     
    end

    St_11(:,:,3)=SOUT11;
    St_12(:,:,3)=SOUT12;
    St_21(:,:,3)=SOUT21;
    St_22(:,:,3)=SOUT22;

%        St_11(:,:,2)=Str_11(:,:,1)+Str_12(:,:,1)*St_11(:,:,3)*(inv(eye(NmodiTETM)-Str_22(:,:,1)*St_11(:,:,3)))*Str_21(:,:,1);
%        St_12(:,:,2)=Str_12(:,:,1)*(inv(eye(NmodiTETM)-St_11(:,:,3)*Str_22(:,:,1)))*St_12(:,:,3);
%        St_21(:,:,2)=St_21(:,:,3)*(inv(eye(NmodiTETM)-Str_22(:,:,1)*St_11(:,:,3)))*Str_21(:,:,1);
%        St_22(:,:,2)=St_22(:,:,3)+St_21(:,:,3)*(inv(eye(NmodiTETM)-Str_22(:,:,1)*St_11(:,:,3)))*Str_22(:,:,1)*St_12(:,:,3); 

%        St_11(:,:,1)=SIN11(:,:)+SIN12(:,:)*St_11(:,:,3)*(inv(eye(NmodiTETM)-SIN22(:,:)*St_11(:,:,3)))*SIN21(:,:);
%        St_12(:,:,1)=SIN12(:,:)*(inv(eye(NmodiTETM)-St_11(:,:,3)*SIN22(:,:)))*St_12(:,:,3);
%        St_21(:,:,1)=St_21(:,:,3)*(inv(eye(NmodiTETM)-SIN22(:,:)*St_11(:,:,3)))*SIN21(:,:);
%        St_22(:,:,1)=St_22(:,:,3)+St_21(:,:,3)*(inv(eye(NmodiTETM)-SIN22(:,:)*St_11(:,:,3)))*SIN22(:,:)*St_12(:,:,3);
    
        St_11(:,:,2)=Str_11(:,:,1)+Str_12(:,:,1)*inv(eye(NmodiTETM)-SOUT11*Str_22(:,:,1))*SOUT11*Str_21(:,:,1);
        St_12(:,:,2)=Str_12(:,:,1)*inv(eye(NmodiTETM)-SOUT11*Str_22(:,:,1))*SOUT12;
        St_21(:,:,2)=SOUT21*inv(eye(NmodiTETM)-Str_22(:,:,1)*SOUT11)*Str_21(:,:,1);
        St_22(:,:,2)=SOUT22+SOUT21*inv(eye(NmodiTETM)-Str_22(:,:,1)*SOUT11)*Str_22(:,:,1)*SOUT12; 

    
%        St_11(:,:,2)=SOUT11+SOUT12*St_11(:,:,3)*(inv(eye(NmodiTETM)-Str_22(:,:,1)*SOUT11))*Str_21(:,:,1);
%        St_12(:,:,2)=SOUT12*(inv(eye(NmodiTETM)-SOUT11*Str_22(:,:,1)))*SOUT12;
%        St_21(:,:,2)=SOUT21*(inv(eye(NmodiTETM)-Str_22(:,:,1)*St_11(:,:,3)))*Str_21(:,:,1);
%        St_22(:,:,2)=SOUT22+SOUT21*(inv(eye(NmodiTETM)-Str_22(:,:,1)*SOUT11))*Str_22(:,:,1)*SOUT12; 

        St_11(:,:,1)=SIN11+SIN12*inv(eye(NmodiTETM)-St_11(:,:,2)*SIN22)*St_11(:,:,2)*SIN21;
        St_12(:,:,1)=SIN12*inv(eye(NmodiTETM)-St_11(:,:,2)*SIN22)*St_12(:,:,2);
        St_21(:,:,1)=St_21(:,:,2)*inv(eye(NmodiTETM)-SIN22*St_11(:,:,2))*SIN21;
        St_22(:,:,1)=St_22(:,:,2)+St_21(:,:,2)*inv(eye(NmodiTETM)-SIN22*St_11(:,:,2))*SIN22*St_12(:,:,2);
    
    % controllo unitarieta per i modi sopra taglio supposti 2 in n3 e n4
     Stot=[St_11(1,1,1),St_11(1,Nmodi+1,1),St_12(1,1,1),St_12(1,Nmodi+1,1);...
         St_11(Nmodi+1,1,1),St_11(Nmodi+1,Nmodi+1,1),St_12(Nmodi+1,1,1),St_12(Nmodi+1,Nmodi+1,1);...
         St_21(1,1,1),St_21(1,Nmodi+1,1),St_22(1,1,1),St_22(1,Nmodi+1,1);...
         St_21(Nmodi+1,1,1),St_21(Nmodi+1,Nmodi+1,1),St_22(Nmodi+1,1,1),St_22(Nmodi+1,Nmodi+1,1)];
    U= Stot'*Stot;
    fi=find(abs(U)<1e-13);
    U(fi)=0;
%    '% Matrice S per lordine 0 TE e TM', keyboard
    
%    'ML', keyboard
    
  ilastm=1;
  
  if ilastm==0
    puor=[1:5];
    
    e1=exp(-j*kz3(puor)*1000).';
    e2=exp(-j*kzu(puor)*1000).';
    S11TETEtotale(puor)=St_11(puor,1,1).*e1;
    S11TETMtotale(puor)=St_11(puor,Nmodi+1,1).*e1;
    S11TMTEtotale(puor)=St_11(Nmodi+puor,1,1).*e1;
    S11TMTMtotale(puor)=St_11(Nmodi+puor,Nmodi+1,1).*e1;
    
    S21TETEtotale(puor)=St_21(puor,1,1).*e1;
    S21TETMtotale(puor)=St_21(puor,Nmodi+1,1).*e1;
    S21TMTEtotale(puor)=St_21(Nmodi+puor,1,1).*e1;
    S21TMTMtotale(puor)=St_21(Nmodi+puor,Nmodi+1,1).*e1;
    
    S12TETEtotale(puor)=St_12(puor,1,1).*e2;
    S12TETMtotale(puor)=St_12(puor,Nmodi+1,1).*e2;
    S12TMTEtotale(puor)=St_12(Nmodi+puor,1,1).*e2;
    S12TMTMtotale(puor)=St_12(Nmodi+puor,Nmodi+1,1).*e2;
    
    S22TETEtotale(puor)=St_22(puor,1,1).*e2;
    S22TETMtotale(puor)=St_22(puor,Nmodi+1,1).*e2;
    S22TMTEtotale(puor)=St_22(Nmodi+puor,1,1).*e2;
    S22TMTMtotale(puor)=St_22(Nmodi+puor,Nmodi+1,1).*e2;
    
 else
     e1=1;
     e2=1;
     pus=2:Nmodi;
     puor=1;
     S11TETEtotale(puor)=St_11(puor,1,1).*e1;
     S11TETMtotale(puor)=St_11(puor,Nmodi+1,1).*e1;
     S11TMTEtotale(puor)=St_11(Nmodi+puor,1,1).*e1;
     S11TMTMtotale(puor)=St_11(Nmodi+puor,Nmodi+1,1).*e1;
     
     S21TETEtotale(puor)=St_21(puor,1,1).*e1;
     S21TETMtotale(puor)=St_21(puor,Nmodi+1,1).*e1;
     S21TMTEtotale(puor)=St_21(Nmodi+puor,1,1).*e1;
     S21TMTMtotale(puor)=St_21(Nmodi+puor,Nmodi+1,1).*e1;
     
     S12TETEtotale(puor)=St_12(puor,1,1).*e2;
     S12TETMtotale(puor)=St_12(puor,Nmodi+1,1).*e2;
     S12TMTEtotale(puor)=St_12(Nmodi+puor,1,1).*e2;
     S12TMTMtotale(puor)=St_12(Nmodi+puor,Nmodi+1,1).*e2;
     
     S22TETEtotale(puor)=St_22(puor,1,1).*e2;
     S22TETMtotale(puor)=St_22(puor,Nmodi+1,1).*e2;
     S22TMTEtotale(puor)=St_22(Nmodi+puor,1,1).*e2;
     S22TMTMtotale(puor)=St_22(Nmodi+puor,Nmodi+1,1).*e2;
     
       puor=2;
           kzu = psqrt((k0*r_in(1))^2 - ktF.^2); 
           e1=exp(-j*kzu(pus)*1000).';
         S11TETEtotale(puor)=sqrt(sum(abs(St_11(pus,1,1).*e1).^2));
         S11TETMtotale(puor)=sqrt(sum(abs(St_11(pus,Nmodi+1,1).*e1).^2));
         S11TMTEtotale(puor)=sqrt(sum(abs(St_11(pus+Nmodi,1,1).*e1).^2));;
         S11TMTMtotale(puor)=sqrt(sum(abs(St_11(pus+Nmodi,1+Nmodi,1).*e1).^2));
         e1=1;
         
         S21TETEtotale(puor)=St_21(puor,1,1).*e1;
         S21TETMtotale(puor)=St_21(puor,Nmodi+1,1).*e1;
         S21TMTEtotale(puor)=St_21(Nmodi+puor,1,1).*e1;
         S21TMTMtotale(puor)=St_21(Nmodi+puor,Nmodi+1,1).*e1;
         
         S12TETEtotale(puor)=St_12(puor,1,1).*e2;
         S12TETMtotale(puor)=St_12(puor,Nmodi+1,1).*e2;
         S12TMTEtotale(puor)=St_12(Nmodi+puor,1,1).*e2;
         S12TMTMtotale(puor)=St_12(Nmodi+puor,Nmodi+1,1).*e2;
         
         S22TETEtotale(puor)=St_22(puor,1,1).*e2;
         S22TETMtotale(puor)=St_22(puor,Nmodi+1,1).*e2;
         S22TMTEtotale(puor)=St_22(Nmodi+puor,1,1).*e2;
         S22TMTMtotale(puor)=St_22(Nmodi+puor,Nmodi+1,1).*e2;
         
           puor=3;
           kzu = psqrt((k0*r_in(1))^2 - ktF.^2); 
           e1=exp(-j*kzu(pus)*1000).';
         S11TETEtotale(puor)=sqrt(mean(abs(diag(St_11(pus,pus,1)).*e1).^2));
         S11TETMtotale(puor)=sqrt(mean(abs(diag(St_11(pus,Nmodi+pus,1)).*e1).^2));
         S11TMTEtotale(puor)=sqrt(mean(abs(diag(St_11(pus+Nmodi,pus,1)).*e1).^2));
         S11TMTMtotale(puor)=sqrt(mean(abs(diag(St_11(pus+Nmodi,pus+Nmodi,1)).*e1).^2));
         
         
         e1=1;
         
         S21TETEtotale(puor)=St_21(puor,1,1).*e1;
         S21TETMtotale(puor)=St_21(puor,Nmodi+1,1).*e1;
         S21TMTEtotale(puor)=St_21(Nmodi+puor,1,1).*e1;
         S21TMTMtotale(puor)=St_21(Nmodi+puor,Nmodi+1,1).*e1;
         
         S12TETEtotale(puor)=St_12(puor,1,1).*e2;
         S12TETMtotale(puor)=St_12(puor,Nmodi+1,1).*e2;
         S12TMTEtotale(puor)=St_12(Nmodi+puor,1,1).*e2;
         S12TMTMtotale(puor)=St_12(Nmodi+puor,Nmodi+1,1).*e2;
         
         S22TETEtotale(puor)=St_22(puor,1,1).*e2;
         S22TETMtotale(puor)=St_22(puor,Nmodi+1,1).*e2;
         S22TMTEtotale(puor)=St_22(Nmodi+puor,1,1).*e2;
         S22TMTMtotale(puor)=St_22(Nmodi+puor,Nmodi+1,1).*e2;         
 end

if length(find(isnan(S11TETEtotale)==1))>0
%' qui ver', keyboard
end

if length(find(abs(S11TETEtotale)>1))>0
%' qui ver', keyboard
end

%' qui ver', keyboard
campi_reticolo_orta