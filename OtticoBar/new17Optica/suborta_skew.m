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

for indlambda=1:passo:Nlambda
    lambda=lambdavet(indlambda);
    
    k0=2*pi/lambda;
    omega=k0*Clight;
    ky= k0*n3*sin(theta)*sin(phi);
    KBd = 2*pi*d/lambda*n3*sin(theta)*cos(phi);   	% sfasamento per cella
    
    % Modi di Floquet nei semispazi
    k = 1:(Nmodi-1)/2;
    Csi(2*k) = -k;
    Csi(2*k+1) = k;
    csi = KBd/d + 2*pi*Csi/d;
    clear k
    ktF=psqrt(csi.^2+ky^2);
    kz3 = psqrt((k0*n3)^2 - ktF.^2); 
    kz4 = psqrt((k0*n4)^2 - ktF.^2); 
    % caricamento dei parametri dei modi di Floquet
    kz(1:Nmodi,1)=kz3.';
    kz(Nmodi+1:NmodiTETM,1)=kz3.';
    kT2(1:Nmodi,1)=(ktF.^2).';
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
    
    
    % Spostamento dei piani di riferimento della porta 2 di tutte le matrici S, eccetto l'ultima
    for k=1:Nstratper
        Propag(:,k)=exp(-j*kz(:,k+1)*thick(k));
        matdiag = (Propag(:,k)*ones(1,NmodiTETM));
        S12(:,:,k)=S12(:,:,k).*matdiag.';
        S21(:,:,k)=matdiag.*S21(:,:,k);
        S22(:,:,k)=(Propag(:,k)*Propag(:,k).').*S22(:,:,k);
    end
    
    
    % CASCATA PER DEFINIRE UN'UNICA MATRICE SCATTERING DI TUTTA LA STRUTTURA
    St_11(:,:,Nstratper+1)=S11(:,:,Nstratper+1);
    St_12(:,:,Nstratper+1)=S12(:,:,Nstratper+1);
    St_21(:,:,Nstratper+1)=S21(:,:,Nstratper+1);
    St_22(:,:,Nstratper+1)=S22(:,:,Nstratper+1);
    % St_ij(:,:,ind) e' la matrice scattering della struttura compresa tra
    % giunzione ind e Nstratper+1 (l'ultima)
    for k=(Nstratper):-1:1
        St_11(:,:,k)=S11(:,:,k)+S12(:,:,k)*St_11(:,:,k+1)*(inv(eye(NmodiTETM)-S22(:,:,k)*St_11(:,:,k+1)))*S21(:,:,k);
        St_12(:,:,k)=S12(:,:,k)*(inv(eye(NmodiTETM)-St_11(:,:,k+1)*S22(:,:,k)))*St_12(:,:,k+1);
        St_21(:,:,k)=St_21(:,:,k+1)*(inv(eye(NmodiTETM)-S22(:,:,k)*St_11(:,:,k+1)))*S21(:,:,k);
        St_22(:,:,k)=St_22(:,:,k+1)+St_21(:,:,k+1)*(inv(eye(NmodiTETM)-S22(:,:,k)*St_11(:,:,k+1)))*S22(:,:,k)*St_12(:,:,k+1); 
    end
    
    
    
    
    % controllo unitarieta per i modi sopra taglio supposti 2 in n3 e n4
    % Stot=[St_11(1,1,1),St_11(1,Nmodi+1,1),St_12(1,1,1),St_12(1,Nmodi+1,1);...
    %     St_11(Nmodi+1,1,1),St_11(Nmodi+1,Nmodi+1,1),St_12(Nmodi+1,1,1),St_12(Nmodi+1,Nmodi+1,1);...
    %     St_21(1,1,1),St_21(1,Nmodi+1,1),St_22(1,1,1),St_22(1,Nmodi+1,1);...
    %     St_21(Nmodi+1,1,1),St_21(Nmodi+1,Nmodi+1,1),St_22(Nmodi+1,1,1),St_22(Nmodi+1,Nmodi+1,1)];
    % Stot'*Stot
    % Matrice S per l'ordine 0 TE e TM
    S11TETEtotale(indlambda)=St_11(1,1,1);
    S11TETMtotale(indlambda)=St_11(1,Nmodi+1,1);
    S11TMTEtotale(indlambda)=St_11(Nmodi+1,1,1);
    S11TMTMtotale(indlambda)=St_11(Nmodi+1,Nmodi+1,1);
    
    S21TETEtotale(indlambda)=St_21(1,1,1);
    S21TETMtotale(indlambda)=St_21(1,Nmodi+1,1);
    S21TMTEtotale(indlambda)=St_21(Nmodi+1,1,1);
    S21TMTMtotale(indlambda)=St_21(Nmodi+1,Nmodi+1,1);
    
    S12TETEtotale(indlambda)=St_12(1,1,1);
    S12TETMtotale(indlambda)=St_12(1,Nmodi+1,1);
    S12TMTEtotale(indlambda)=St_12(Nmodi+1,1,1);
    S12TMTMtotale(indlambda)=St_12(Nmodi+1,Nmodi+1,1);
    
    S22TETEtotale(indlambda)=St_22(1,1,1);
    S22TETMtotale(indlambda)=St_22(1,Nmodi+1,1);
    S22TMTEtotale(indlambda)=St_22(Nmodi+1,1,1);
    S22TMTMtotale(indlambda)=St_22(Nmodi+1,Nmodi+1,1);




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
            figure(2*ind-1)
            subplot(2,1,1)
%            plot(x,abs(Exhat))
            plot(x,real(Exhat))
            grid on
            title(['interfaccia  ',num2str(ind),'  campo Ex'])
            hold on
%            plot(x,abs(Extilda),'r')
            plot(x,real(Extilda),'r')
            subplot(2,1,2)
%            plot(x,abs(Eyhat))
            plot(x,real(Eyhat))
            grid on
            title(['interfaccia  ',num2str(ind),'  campo Ey'])
            hold on
%            plot(x,abs(Eytilda),'r')
            plot(x,real(Eytilda),'r')
            figure(2*ind)
            subplot(2,1,1)
%            plot(x,abs(Hxhat))
            plot(x,real(Hxhat))
            grid on
            title(['interfaccia  ',num2str(ind),'  campo Hx'])
            hold on
%            plot(x,abs(Hxtilda),'r')
            plot(x,real(Hxtilda),'r')
            subplot(2,1,2)
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
end %%for indlambda=1:Nlambda
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

