%   Afzf=me_autu(Azvetf,pou,nso,sz,Fint,ze,if_only_out,icsmax,ifpsave,sgu);

function [A,sgu]=me_autn(Anu,pou,nso,sz,Fint,ze,if_only_out,icsmax,ifp,sgi)
if nargin<10
 ifsgi=zeros(1,length(Fint));
else
 ifsgi=sgi;
end
lep=size(pou);
 for kmo=1:icsmax
%  for iaf=1:lep(2)
  for iaf=1:size(pou,2)
   if length(if_only_out)==0
    Pun=pou(nso,iaf);
    if Pun>0
     Ada=reshape(Anu(:,Pun,kmo,iaf),sz(1),1);
    else
     Ada=zeros(sz(1),1);
    end
%    'ver ADA', keyboard
%    Ada=reshape(Anu(:,pou(nso,iaf),kmo,iaf),sz(1),1);
   else
    Pun=pou(nso,iaf);
    if Pun>0
     Ada=reshape(Anu(:,Pun,iaf),sz(1),1);
    else
     Ada=zeros(sz(1),1);
    end
%    Ada=reshape(Anu(:,pou(nso,iaf),iaf),sz(1),1);
   end
   if kmo==1
    As0(:,iaf)=Ada;
   end
%   if iaf==1
%    Ada=Ada/max(Ada)*abs(max(Ada));
%    [du,ima]=max(Ada);
%   else
%    Ada=Ada/Ada(ima)*abs(Ada(ima));
%   end
%' qui 0', keyboard
  iold=0;
  if iold==1
   if ifsgi(1,kmo)==0
    if iaf==1
     [du,ima]=max(Ada);
     sga0=sign(real(Ada(ima)));
     sgu(1,kmo)=sga0;
     sgu(2,kmo)=ima;
    else
     sgAi=sign(real(Ada(ima)));
     Ada=Ada*sgAi*sga0;
    end
   else
    if iaf==1
     ima=ifsgi(2,kmo);
     sga0=ifsgi(1,kmo);
     sgAi=sign(real(Ada(ima)));
     Ada=Ada*sgAi*sga0;
     sgu(1,kmo)=sga0;
     sgu(2,kmo)=ima;
    else
     sgAi=sign(real(Ada(ima)));
     Ada=Ada*sgAi*sga0;
    end
   end
  end  %iold

   As(:,iaf)=Ada;
%   'As', keyboard

  end

      Ares{kmo}=As;
 end %kmo
 Aresv=Ares;

if ifsgi(1)==0
  Ac=Ares{1};
  nf=size(pou,2);
  for kf=1:nf
   nor=sqrt(Ac(:,kf)'*Ac(:,kf));
   if nor~=0
    Ac(:,kf)=Ac(:,kf)/nor;
   else
    Ac(:,kf)=Ac(:,kf);
   end
   nvd(kf)=nor;
  end
%  'end primo', keyboard
  sA=sum(Ac);
%  [ma,ima]=max(Ac(:,1));
  [ma,ima1]=max(sA);
  [ma,ima]=max(Ac(:,ima1));  
  for kf=1:nf
   nor=Ac(ima,kf);
   if nor~=0
    Ac1(:,kf)=Ac(:,kf)/Ac(ima,kf);
    ac(kf)=Ac1(:,kf)'*Ac1(:,1);
    nv(kf)=1/nvd(kf)/Ac(ima,kf);
   else
    Ac1(:,kf)=Ac(:,kf);
    ac(kf)=Ac1(:,kf)'*Ac1(:,1);
    nv(kf)=0;   
   end
%   'cont ', keyboard
  end
%  figure, plot(Fint,ac/ac(1)), pausak
  Ares{1}=Ac1;
  for ka=2:kmo
   Ac=Ares{ka}*diag(nv);
   Ares{ka}=Ac;
  end
  sgu=nv;
else
 nv=ifsgi;
  for ka=1:kmo
   Ac=Ares{ka}*diag(nv);
   Ares{ka}=Ac;
  end
end
%Ares=Aresv;
%'dopo Ares', keyboard

 for kmo=1:icsmax
   As=Ares{kmo};
   [dude,isoz]=sort(abs(Fint-ze));
   DeFi=abs(diff(Fint(isoz(1:2))));
   if ze<Fint(end) & ze>Fint(1)
    A(:,kmo)=(1-dude(1)/DeFi)*As(:,isoz(1))+(1-dude(2)/DeFi)*As(:,isoz(2));
   else
    A(:,kmo)=As(:,isoz(1));
   end
 end
% for kmo=1:icsmax
%  if ifsgi(1,kmo)~=0
%   A(:,kmo)=A(:,kmo)*ifsgi(1,kmo)*sgu(1,kmo);
%  end
% end
% ' qui 1', keyboard
return

KK=1:length(A);
figure, plot(KK,real(Ares{1}),KK,real(A(:,1)),'w.-'),
title(' row eigenv. at output' )
pausak
[du,fimi]=sort(abs(Fint-ze));
figure, plot(KK,real(Ares{1}(:,fimi(1:5))),KK,real(A(:,1)),'w.-'),
%if length(Fint)>20
%hold  on
% plot(KK,real(Ares{1}(:,1:10:end)),'.'),
%end
pausak
%return
for kmo0=1:3
figure, plot(KK,abs(A(:,kmo0)),'w.',KK,abs(Ares{kmo0}*seg(kmo0)),'--'), pausak
figure, plot(KK,real(A(:,kmo0)),'w.',KK,real(Ares{kmo0}*seg(kmo0)),'--'), pausak
figure, plot(KK,imag(A(:,kmo0)),'w.',KK,imag(Ares{kmo0}*seg(kmo0)),'--'), keyboard
end

Ac=Ares{2};
nf=length(Fint);
for kf=1:nf
 nor=Ac(:,kf)'*Ac(:,kf);
 Ac(:,kf)=Ac(:,kf)/sqrt(nor);
end
[ma,ima]=max(Ac(:,1));
for kf=1:nf
 Ac1(:,kf)=Ac(:,kf)/Ac(ima,kf);
 ac(kf)=Ac1(:,kf)'*Ac1(:,1);
end
figure, plot(Fint,ac/ac(1))
keyboard