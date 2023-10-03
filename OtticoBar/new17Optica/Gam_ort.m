function [S11]=Gam_ort(kx,nindex_in,nindex,nindex_out,d,lambda0,pol)
%function [S11,S12,S21,S22,argS11,argS12,argS21,argS22]=Gam_ort(kx,nindex_in,nindex,nindex_out,d,lambda0,pol)

% La function calcola la matrice S di una statificazione dielettrica. I materiali possono avere perdite o guadagno.
% Vengono calcolate le fasi dei parametri S con unwrapping
% I piani di riferimento sono asimmetrici, come su figura di pagina 12
% appunti, ossia
%            nindex_in| nindex(1) |...... | nindex(Nlayer)| nindex_out
%                     |           |       |               | 
%                     |           |       |               | 
%                     |           |       |               | 
%                     |           |       |               | 
%                    -A                                   B+
% i piani di riferimento sono in A- e B+

% beta2 e' il quadrato di k0*nindex_in*sin(theta), cioe' la componente dei vettori
% d'onda nei vari strati parallela alle interfacce
%
% Nlayer e' il numero degli strati
% nindex(Nlayer) e d(Nlayer) sono vettori contenenti gli indici di rifrazione e gli spessori degli strati
% lambda0 e' la lunghezza d'onda nel vuoto e pol='TE' o pol='TM' specificano la polarizzazione

% %Struttura periodica con N dielettrici con perdite, ricava la matrice S

beta2=kx.^2;
Nlayer=length(nindex);
Nbeta2=length(beta2);

rFre=zeros(Nlayer+1,Nbeta2); % coeff. di Fresnel delle varie interfacce
Gam=zeros(Nlayer+1,Nbeta2); % coeffic. di riflessione nella sezione meno di ogni interfaccia
Gap=zeros(Nlayer,Nbeta2); % coeffic. di riflessione nella sezione piu' di ogni interfaccia
kx=zeros(Nlayer,Nbeta2); % kx nei vari strati


k00=2*pi/lambda0;

% Computation of S11 of the cell
if imag(nindex_out)<=0
    kx_out=psqrt(k00^2*nindex_out^2-beta2);
else % lo strato _out e' attivo
    kx_out=asqrt(k00^2*nindex_out^2-beta2);
end
if imag(nindex(Nlayer))<=0
    kx(Nlayer,:)=psqrt(k00^2*nindex(Nlayer)^2-beta2);
else % lo strato (Nlayer) e' attivo
    kx(Nlayer,:)=asqrt(k00^2*nindex(Nlayer)^2-beta2);
end

if strcmp(pol,'TE')
    rFre(Nlayer+1,:)=(kx(Nlayer,:)-kx_out)./(kx(Nlayer,:)+kx_out);
elseif strcmp(pol,'TM')
    rFre(Nlayer+1,:)=(kx_out/nindex_out^2-kx(Nlayer,:)/nindex(Nlayer)^2)...
        ./(kx_out/nindex_out^2+kx(Nlayer,:)/nindex(Nlayer)^2);
end
Gam(Nlayer+1,:)=rFre(Nlayer+1,:);
argS11=angle(Gam(Nlayer+1,:));

for ind=Nlayer:-1:2
    Gap(ind,:)=Gam(ind+1,:).*exp(-j*2*kx(ind,:)*d(ind));
    if imag(nindex(ind-1))<=0
        kx(ind-1,:)=psqrt(k00^2*nindex(ind-1)^2-beta2);
    else % lo strato Nlayer e' attivo
        kx(ind-1,:)=asqrt(k00^2*nindex(ind-1)^2-beta2);
    end

	if strcmp(pol,'TE')
        rFre(ind,:)=(kx(ind-1,:)-kx(ind,:))./(kx(ind-1,:)+kx(ind,:));
	elseif strcmp(pol,'TM')
        rFre(ind,:)=(kx(ind,:)/nindex(ind)^2-kx(ind-1,:)/nindex(ind-1)^2)...
            ./(kx(ind,:)/nindex(ind)^2+kx(ind-1,:)/nindex(ind-1)^2);
	end
    Gam(ind,:)=(rFre(ind,:)+Gap(ind,:))./(1+rFre(ind,:).*Gap(ind,:));
    argS11=argS11+angle(Gam(ind,:))-angle(Gap(ind,:))-2*real(kx(ind,:))*d(ind);
end % for ind=Nlayer:2
Gap(1,:)=Gam(2,:).*exp(-j*2*kx(1,:)*d(1));
if imag(nindex_in)<=0
    kx_in=psqrt(k00^2*nindex_in^2-beta2);
else % lo strato _in e' attivo
    kx_in=asqrt(k00^2*nindex_in^2-beta2);
end

if strcmp(pol,'TE')
    rFre(1,:)=(kx_in-kx(1,:))./(kx_in+kx(1,:));
elseif strcmp(pol,'TM')
    rFre(1,:)=(kx(1,:)/nindex(1)^2-kx_in/nindex_in^2)...
        ./(kx(1,:)/nindex(1)^2+kx_in/nindex_in^2);
end
Gam(1,:)=(rFre(1,:)+Gap(1,:))./(1+rFre(1,:).*Gap(1,:));
argS11=argS11+angle(Gam(1,:))-angle(Gap(1,:))-2*real(kx(1,:))*d(1);


% Gam(1) e' S11 e argS11 e' la sua fase unwrapped
S11=Gam(1,:);

return

% Computation of S21 of the cell
S21defas=ones(1,Nbeta2); % solo modulo
fasmen=zeros(1,Nbeta2); % accumula la fase di 1+Gam(ind,:)
faspiu=zeros(1,Nbeta2);  % accumula la fase di 1+Gap(ind,:)
fasexp=zeros(1,Nbeta2); % accumula la fase degli exp
for ind=1:Nlayer
    S21defas=S21defas.*abs((1+Gam(ind,:))./(1+Gap(ind,:)))...
        .*exp(imag(kx(ind,:))*d(ind));
    fasmen=fasmen+anglemod(1+Gam(ind,:),pi/8);
    faspiu=faspiu+anglemod(1+Gap(ind,:),pi/8);
    fasexp=fasexp-real(kx(ind,:))*d(ind);
end
S21defas=S21defas.*abs(1+Gam(Nlayer+1,:));
if strcmp(pol,'TE') % impedance factor
    impfact=sqrt(kx_out./kx_in);
elseif strcmp(pol,'TM')
    impfact=sqrt((kx_in./nindex_in^2)./(kx_out/nindex_out^2));
end
S21defas=S21defas.*abs(impfact);
fasmen=fasmen+angle(1+Gam(Nlayer+1,:));
argS21=fasexp+fasmen-faspiu+angle(impfact); %argS21 e' la fase unwrapped di S21
S21=S21defas.*exp(j*argS21);
S12=S21;
argS12=argS21;

% Computation of S22 of the cell

nindex_outrev=nindex_in;
nindex_inrev=nindex_out;

for ind=1:Nlayer
    nindexrev(ind)=nindex(Nlayer+1-ind);
    drev(ind)=d(Nlayer+1-ind);
    kxrev(ind,:)=kx(Nlayer+1-ind,:);
end
for ind=1:Nlayer+1
    rFrerev(ind,:)=-rFre(Nlayer+2-ind,:);
end
Gam(Nlayer+1,:)=rFrerev(Nlayer+1,:);
argS11=angle(Gam(Nlayer+1,:));

for ind=Nlayer:-1:2
    Gap(ind,:)=Gam(ind+1,:).*exp(-j*2*kx(ind,:)*d(ind));
    Gam(ind,:)=(rFrerev(ind,:)+Gap(ind,:))./(1+rFrerev(ind,:).*Gap(ind,:));
    argS11=argS11+angle(Gam(ind,:))-angle(Gap(ind,:))-2*real(kx(ind,:))*d(ind);
end % for ind=Nlayer:2
Gap(1,:)=Gam(2,:).*exp(-j*2*kx(1,:)*d(1));

Gam(1,:)=(rFrerev(1,:)+Gap(1,:))./(1+rFrerev(1,:).*Gap(1,:));
argS11=argS11+angle(Gam(1,:))-angle(Gap(1,:))-2*real(kx(1,:))*d(1);


% Gam(1) e' S22 e argS22 e' la sua fase unwrapped
S22=Gam(1,:);

