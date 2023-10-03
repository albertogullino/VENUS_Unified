function [A,sgu]=me_autn(Anu,pou,nso,sz,Fint,ze,if_only_out,icsmax,ifp,sgi)
if nargin<10
 ifsgi=zeros(2,icsmax);
else
 ifsgi=sgi;
end
 for kmo=1:icsmax
  for iaf=1:length(Fint)
   if length(if_only_out)==0
    Pun=pou(nso,iaf);
    if Pun>0
     Ada=reshape(Anu(:,Pun,kmo,iaf),sz(1),1);
    else
     Ada=zeros(sz(1),1);
    end
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

   As(:,iaf)=Ada;
%   [du,ima]=max(abs(Ada));
%   As(:,iaf)=Ada*sign(real(Ada(ima)));

%    ma=[real(Ada(ima)) imag(Ada(ima))];
%    [du,ima]=max(abs(ma));
%    As(:,iaf)=Ada*sign(ma(ima));
%    As(:,iaf)=Ada/Ada(ima);
%    As(:,iaf)=Ada;
  end

      Ares{kmo}=As;

   [dude,isoz]=sort(abs(Fint-ze));
   DeFi=abs(diff(Fint(isoz(1:2))));
   if ze<Fint(end) & ze>Fint(1)
    A(:,kmo)=(1-dude(1)/DeFi)*As(:,isoz(1))+(1-dude(2)/DeFi)*As(:,isoz(2));
   else
    A(:,kmo)=As(:,isoz(1));
   end
 end
 for kmo=1:icsmax
  if ifsgi(1,kmo)~=0
   A(:,kmo)=A(:,kmo)*ifsgi(1,kmo)*sgu(1,kmo);
  end
 end
return
KK=1:length(A);
figure, plot(KK,real(Ares{1}),KK,real(A(:,1)),'w'),
title(' row eigenv. at output' )
pausak
return
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
