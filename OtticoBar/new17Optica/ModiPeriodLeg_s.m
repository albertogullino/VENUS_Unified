function [kzG_TE2,kTG_TE2,ZG_TE,kzG_TM2,kTG_TM2,ZG_TM, M1, M2,VTE1, VTE2, VTE1H, VTE2H, ITE1, ITE2, ITE1H, ITE2H ,...
        VTM1, VTM2, VTM1H, VTM2H,ITM1, ITM2, ITM1H, ITM2H] =...
    ModiPeriodLeg(n1,n2,d1,d2,n3,theta,phi,lambda,Nmodi,errore)


%
% -------------------------------------------------------
% R. Orta Dicembre 2010
% -------------------------------------------------------
% Funzione per il calcolo delle costanti di propagazione
% in una regione stratificata periodica e delle relative
% autofunzioni. Il tutto e' ricavata da un'equazione agli
% autovalori.
% Si calcola inoltre il numero di polinomi di Legendre in 
% ogni strato in modo tale da avere un certo errore 
% massimo sulle costanti di propagazione tramite una
% formula empirica ricavata da Vito Lancellotti (vedi
% articolo).
% Le autofunzioni in uscita sono normalizzate in potenza.
% Gli indici di rifrazione possone essere complessi. Sono ricavati i modi
% della struttura originale e di quella aggiunta.
%
% -------------------------------------------------------
% Variabili d'ingresso:
% lambda	lunghezza d'onda a cui vengono fatti i calcoli
%
% Variabili d'uscita:
% kzG_TE2		vettore delle costanti di propagazione al
%			quadrato
% kTG_TE2   vettore dei k trasversali al quadrato
%
% ZG_TE     vettore delle impedenze modali
%
% M1, M2	ordine massimo dei polinomi di Legendre per
%			lo sviluppo nei due strati della regione period.
% VTE1, VTE2 ecc.	matrici degli autovettori relativi alle auto-
%			funzioni modali nei due strati della regione
%			stratificata.
% VTE1H,VTE2H ecc.  idem per il problema aggiunto
%  ind. riga --> ordine polinomio di Legendre
%  ind. colonna --> ordine del modo, quindi P1 e P2 hanno lo stesso numero
%  di colonne

% -------------------------------------------------------
%

Clight=2.99792458e8; %nm / ns, quindi frequenze in GHz
mu0=4*pi*1e-7;
d=d1+d2; %periodo
k0=2*pi/lambda;
omega=k0*Clight;
eps0=k0/(Clight*omega*mu0);

polariz='TE';
 
% determinazione dei modi TE della regione stratificata periodica
[kzG_TE2, kTG_TE2, M1, M2, VTE1, VTE2, VTE1H, VTE2H, ITE1, ITE2, ITE1H, ITE2H ] = ...
    modiLeg_s(n1,n2,d1,d2,n3,theta,phi,polariz,lambda,Nmodi,errore);
kzG_TE=psqrt(kzG_TE2);
 ZG_TE=omega*mu0*kzG_TE./kTG_TE2;
% determinazione dei modi TM della regione stratificata periodica
polariz='TM';
[kzG_TM2, kTG_TM2, M1, M2, ITM1, ITM2, ITM1H, ITM2H, VTM1, VTM2, VTM1H, VTM2H ] = ...
    modiLeg_s(n1,n2,d1,d2,n3,theta,phi,polariz,lambda,Nmodi,errore);

kzG_TM=psqrt(kzG_TM2);
 ZG_TM=kTG_TM2./kzG_TM/(omega*eps0);

