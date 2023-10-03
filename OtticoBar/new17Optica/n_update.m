fre=freq;
nv=nv0;
' nv0', keyboard
mass0=9.1e-31;
mi=pi*4e-7;
eps0=8.8541e-12;
c=1/sqrt(mi*eps0);
Z0=120*pi;
q=1.6e-19;

lambda_m=lambda*1e-6;

%si=size(nv);
%for ksi=si(2):-1:1
% xal=xm(:,ksi);
% imatw=xal;
% fial=find(xal>0);
% imatw(fial)=1;
% fisos=find(imatw~=0);
% nal=ref_ind(imatw,lambda_m,xal);
% nv(fisos,ksi)=nal(fisos);
%end

 si=size(nv);
 for ksi=1:si(1)
  xal=xm(ksi,:);
  imatw=xal;
  fial=find(xal>0);
  if length(fial)>0
   imatw(fial)=1;
   fisos=find(imatw~=0 & imatw~=-10);
   if ivfre==1
    nvv=nv0(ksi,fisos);
    nal=ref_ind(nvv,imatw(fisos),lambda0,xal(fisos));
    nv(ksi,fisos)=nal;
    nv0(ksi,fisos)=nal;
   end
  end  %length
 end




calfa=1e10*Z0*q^3*lambda_m^3/(mass0*2*pi*c)^2/(4*pi);
cenne=1e6*q^2*lambda_m^2/(2*eps0*mass0*(2*pi*c)^2);


fiA=find(dglo.Dop>0 & imatw>0);
palfae=calfa*abs(dglo.Dop(fiA))./...
(dglo.mue(fiA).*dglo.me(fiA).^2.*real(nv(fiA,1)))/(1+fre)^3;
pennee=cenne*abs(dglo.Dop(fiA))./(dglo.me(fiA).*real(nv(fiA,1)))/(1+fre)^3;
%figure, plot(palfae)
%figure, plot(pennee,'r')
nv(fiA,1)=nv(fiA,1)+pennee-j*palfae;

fiA=find(dglo.Dop<0);
palfah=calfa*abs(dglo.Dop(fiA))./...
(dglo.muh(fiA).*dglo.mh(fiA).^2.*real(nv(fiA,1)))/(1+fre)^3;
penneh=cenne*abs(dglo.Dop(fiA))./(dglo.mh(fiA).*real(nv(fiA,1)))/(1+fre)^3;
%figure, plot(palfah)
%figure, plot(penneh,'r')

nv(fiA,1)=nv(fiA,1)+penneh-j*palfah;
%disp('n_update')
%keyboard

iput=find(iauto(:,1)==1);
ipua=find(iauto(:,1)==2);
ipub=find(iauto(:,1)==3);
ipudd=find(iauto(:,2)==-4 | iauto(:,2)>0);

put=iput+1:ipua-1;
pub=ipua+1:ipub-1;
pua=ipua;

% ifimem=find(ifield~=0);
 ifimem=find(ifield==-2);
 istfiv=zeros(size(ifield));
% istfiv(ifimem)=ifimem-1-iput;
% ifip=find(ifield==-2)-1-iput;
 istfiv(ifimem)=ifimem-iput;
 ifip=find(ifield==-2)-iput;
 if length(ifip)>0
  if ifip>0
   istfiv(ifip)=-istfiv(ifip);
  end
 end
% disp(' verifica'), keyboard

put_pl=2:ipua-1;
pub_pl=ipua+1:length(iauto)-1;

% setto variabili per struttura completamente planare
%
% uscita

nt=nv(put_pl,1);
ft=fst(put_pl,:);
fdu=find(diff(ft(:,1))~=0);

if length(fdu)==2
 pu1=1:fdu(1);
 pu1=fliplr(pu1);
 npla.t.o=nt(pu1);

 pu1=fdu(1)+1:fdu(2);
 pu1=fliplr(pu1);
 npla.t.m=nt(pu1);

 pu1=fdu(2)+1:length(nt);
 pu1=fliplr(pu1);
 npla.t.i=nt(pu1);

else
 pu1=1:fdu(1);
 pu1=fliplr(pu1);
 npla.t.m=nt(pu1);

 pu1=fdu(1)+1:length(nt);
 pu1=fliplr(pu1);
 npla.t.i=nt(pu1);

 npla.t.o=[];
end


% sotto

nt=nv(pub_pl,1);
ft=fst(pub_pl,:);
fdu=find(diff(ft(:,1))~=0);

if length(fdu)==2
 pu1=1:fdu(1);
 npla.b.i=nt(pu1);

 pu1=fdu(1)+1:fdu(2);
 npla.b.m=nt(pu1);

 pu1=fdu(2)+1:length(nt);
 npla.b.o=nt(pu1);

elseif length(fdu)==1
 pu1=fdu(1)+1:length(nt);
 npla.b.m=nt(pu1);

 pu1=1:fdu(1);
 npla.b.i=nt(pu1);

 npla.b.o=[];

else

 pu1=fdu(1)+1:length(nt);
 npla.b.m=nt(pu1);

 pu1=1:fdu(1);
 npla.b.i=nt(pu1);

 npla.b.o=[];
end

nr.t=nv(put,:);
nr.b=nv(pub,:);
nr.a=nv(pua,:);




pudu=1:iput;
if length(pudu)>1
 pupt=find(abs(iauto(pudu,2))==3);
 pup.t.last=pudu(pupt(1));
 pup.t.o=fliplr(pudu(pupt(2:length(pupt))));
 pupt=find(abs(iauto(pudu,2))==2);
 pup.t.m=fliplr(pudu(pupt));
 pupt=find(abs(iauto(pudu,2))==1);
 pup.t.i=fliplr(pudu(pupt));
else
 pup.t.last=pudu(1);
 pup.t.o=[];
 pup.t.m=[];
 nmir.t=0;
 pup.t.i=[];
end

pudu=ipub:length(dv);
pupt=find(abs(iauto(pudu,2))==3);
if length(pupt)>0
 pup.b.o=pudu(pupt(1:length(pupt)-1));
 pup.b.last=pudu(pupt(length(pupt)));
else
 pup.b.o=[];
 pup.b.last=[];
end

pupt=find(abs(iauto(pudu,2))==2);
if length(pupt)>0
 pup.b.m=pudu(pupt);
 imir=find(abs(iauto(pudu,2))==2);
 nmir.b=fst(pudu(imir(1)),2);
else
 pup.b.m=[];
 nmir.b=0;
end
pupt=find(abs(iauto(pudu,2))==1);
if length(pupt)>0
 pup.b.i=pudu(pupt);
else
 pup.b.i=[];
end




npl.t.last=nv(pup.t.last,1);
npl.t.o=nv(pup.t.o,1);
npl.t.m=nv(pup.t.m,1);
npl.t.i=nv(pup.t.i,1);

npl.b.last=nv(pup.b.last,1);
npl.b.o=nv(pup.b.o,1);
npl.b.m=nv(pup.b.m,1);
npl.b.i=nv(pup.b.i,1);


  nvbr=npl.b.i;
  nuvb=npl.b.o;
  nbb =npl.b.m;
  rfd =npl.b.last;

  nvb=[nr.b(:,1); npl.b.i];

  nvtr=npl.t.i;
  nuv=npl.t.o;

  nbt =npl.t.m;
  rfu =npl.t.last;

% per struttura non planare e/o  anisotropa

  nitn=nr.t.';


  nib=nr.b.';

  niat=nr.a.';
  rqw=niat(1);

  nitot=[nitn.'; niat.'; nib.'];


% per struttura completamente planare (perdite e acc. uscita in TL_ef.m)
%
% sopra
%
nplit=npla.t.i;
nplt=npla.t.m;
nplut=npla.t.o;
%
% sotto
%
nplib=npla.b.i;
nplb=npla.b.m;
nplub=npla.b.o;



