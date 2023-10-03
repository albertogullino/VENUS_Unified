function [hfunx, hfuny, efunx, efuny] = modefunComp(x1, d1, n1, x2, d2, n2, P1C, P2C, Q1C, Q2C,...
    betaz2, kTG2, ky, k0, polariz, Clight )

%
%  Renato Orta Dicembre 2010
% 
% ------------------------------------------------------
% Funzione per il calcolo dell'andamento dei modi (sia
% campo elettrico che magnetico) di una regione stratifi-
% cata periodica calcolati come sviluppo in serie di 
% polinomi di Legendre e usando la formula di ricorrenza
% di Clenshaw.
% ------------------------------------------------------
% Variabili in ingresso:
%	x1 = vettore contenente i valori di x tra 0 e d1;
%	x2 = vettore contenente i valori di x tra d1 e d;
%	d1 = larghezza dello strato 1
%	d2 = larghezza dello strato 2
%	n1 = costante dielettrica strato 1
%	n2 = costante dielettrica strato 2
%	M1 = ordine massimo dei
%		polinomi di Legendre nello st due strati
%	P1C, P2C = matrici contenenti i coefficienti dello
%		sviluppo in serie nei due strati del 
% 		campo
%	polariz = polarizzazione dell'onda piana incidente
%				= 0 --> TE
%				= 1 --> TM
%
% Variabili in uscita:
%	hfun ecc = matrice contenente 
%		l'andamento delle autofunzioni h:
%  ind. riga --> ordine modo
%  ind. colonna --> valori di x
%	efun ecc= matrice contenente 
%		l'andamento delle autofunzioni e
%  ind. riga --> ordine modo
%  ind. colonna --> valori di x
%-------------------------------------------------------
mu0=4*pi*1.e-7;
omega=k0*Clight;
eps0=k0/(Clight*omega*mu0);
betaz=psqrt(betaz2);

d = d1+d2;
M1=size(P1C,1)-1;
M2=size(P2C,1)-1;
Nmodi =  size(P1C,2);
Nx1=length(x1);
Nx2=length(x2);

hfunx=zeros(Nmodi,Nx1+Nx2);
hfuny=zeros(Nmodi,Nx1+Nx2);
efunx=zeros(Nmodi,Nx1+Nx2);
efuny=zeros(Nmodi,Nx1+Nx2);


for ordine = 1:Nmodi
	cs1 = (2*x1 - d1)/d1;
	yx = zeros(M1+2,size(cs1,2));
  	yy = zeros(M1+2,size(cs1,2));
	N1 = 0:M1+1;
	gamma1 = ((2*N1+1)./(N1+1)).' * cs1;
	gamma2 = -N1./(N1+1);
   
  	for k = M1:-1:1
		yx(k,:) = gamma1(k+1,:).*yx(k+1,:)+gamma2(k+2).*yx(k+2,:)+P1C(k+1,ordine);
	end
	Xn1 = gamma2(2)*1*yx(2,:)+cs1.*yx(1,:)+1*P1C(1,ordine);

    for k = M1:-1:1
		yx(k,:) = gamma1(k+1,:).*yx(k+1,:)+gamma2(k+2).*yx(k+2,:)+Q1C(k+1,ordine);
	end
	Wn1 = gamma2(2)*1*yx(2,:)+cs1.*yx(1,:)+1*Q1C(1,ordine);
   
	cs2 = (2*x2 - d1-d)/d2;
	Yy = zeros(M2+2,size(cs2,2));
	Yx = zeros(M2+2,size(cs2,2));
	N2 = 0:M2+1;
	Gamma1 = ((2*N2+1)./(N2+1)).' * cs2;
	Gamma2 = -N2./(N2+1);
   
  	for k = M2:-1:1
		Yx(k,:) = Gamma1(k+1,:).*Yx(k+1,:)+Gamma2(k+2).*Yx(k+2,:)+P2C(k+1,ordine);
	end
	Xn2 = Gamma2(2)*1*Yx(2,:)+cs2.*Yx(1,:)+1*P2C(1,ordine);

    for k = M2:-1:1
		Yx(k,:) = Gamma1(k+1,:).*Yx(k+1,:)+Gamma2(k+2).*Yx(k+2,:)+Q2C(k+1,ordine);
	end
	Wn2 = Gamma2(2)*1*Yx(2,:)+cs2.*Yx(1,:)+1*Q2C(1,ordine);
    
	vettoreV = [Xn1,Xn2];
    if strcmp(polariz,'TE')
        efuny(ordine,:) = vettoreV;
        hfunx(ordine,:) = -vettoreV;
        ZinfG=omega*mu0*betaz(ordine)/kTG2(ordine);
        hfuny(ordine,:) = ZinfG*[Wn1,Wn2]*ky/betaz(ordine);
    elseif strcmp(polariz,'TM')
        efunx(ordine,:) = [Xn1/n1^2, Xn2/n2^2];
        hfuny(ordine,:) = vettoreV;
        ZinfG=kTG2(ordine)/betaz(ordine)/(omega*eps0);
        efuny(ordine,:) = -[Wn1,Wn2]*ky/betaz(ordine)/ZinfG;
    end
    
end  %%for ordine = 1:Nmodi
