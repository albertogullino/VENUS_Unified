function [eFloqxAn, eFloqyAn, hFloqxAn, hFloqyAn] = FloqzAnComp(x1, d1, x2, d2,phi,csi, ky, k0, nindex, polariz,segem)

seg=segem;
'passo', keyboard
% segno relativo TE TM; =1, cambio segno relativo fra TE e TM rispetto convenzione fun. gen. scalare.

%
% Autore: Renato Orta
% Torino, 02 aprile 2010
% ------------------------------------------------------
% Riga di comando:
%
% [hFloqAn, eFloqAn] = FloqAnComp(x1, d1, n1, x2, d2, n2,Nmodi, polariz);
%
% ------------------------------------------------------
% Funzione per il calcolo dell'andamento dei modi di Floquet(sia
% campo elettrico che magnetico) TE e TM rispetto a z!!!!!!!!!!!!!!!
% ------------------------------------------------------
% Variabili in ingresso:
%	x1 = vettore contenente i valori di x tra 0 e d1;
%	x2 = vettore contenente i valori di x tra d1 e d;
%	d1 = larghezza dello strato 1
%	d2 = larghezza dello strato 2
%	n1 = costante dielettrica strato 1
%	n2 = costante dielettrica strato 2
%	polariz = polarizzazione dell'onda piana incidente
%				= 0 --> TE
%				= 1 --> TM
%
% Variabili in uscita:
%	hFloqAn = matrice contenente 
%		l'andamento delle autofunzioni h:
%  ind. riga --> ordine modo
%  ind. colonna --> valori di x
%	eFloqAn= matrice contenente 
%		l'andamento delle autofunzioni e
%  ind. riga --> ordine modo
%  ind. colonna --> valori di x
%-------------------------------------------------------

d = d1+d2;
x=[x1,x2];
kt=psqrt(csi.^2+ky^2);
Nmodi=length(csi);
kxsukt=cos(phi)*ones(1,Nmodi);
ind0=find(kt~=0);
kxsukt(ind0)=(csi(ind0)./kt(ind0));

kysukt=sin(phi)*ones(1,Nmodi);
ind0=find(kt~=0);
kysukt(ind0)=(ky./kt(ind0));


for ind=1:Nmodi
    Floq=exp(-j*csi(ind)*x)/sqrt(d);
    % indice di riga : modo di Floquet
    % indice di colonna: punto x
    
    if strcmp(polariz,'TE')==1
        eFloqxAn(ind,:)=-j*kysukt(ind)*Floq*seg;
        eFloqyAn(ind,:)=j*kxsukt(ind)*Floq*seg;
        hFloqxAn(ind,:)=-j*kxsukt(ind)*Floq*seg;
        hFloqyAn(ind,:)=-j*kysukt(ind)*Floq*seg;    
    else %TM
        eFloqxAn(ind,:)=j*kxsukt(ind)*Floq;
        eFloqyAn(ind,:)=j*kysukt(ind)*Floq;
        hFloqxAn(ind,:)=-j*kysukt(ind)*Floq;
        hFloqyAn(ind,:)=j*kxsukt(ind)*Floq;
        
    end
end