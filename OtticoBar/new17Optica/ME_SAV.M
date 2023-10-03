function A=me_aut(Anu,pou,nso,sz,Fint,ze,if_only_out,icsmax,ifp)
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
   if iaf==1
    Ada=Ada/max(Ada)*abs(max(Ada));
    [du,ima]=max(Ada);
   else
    Ada=Ada/Ada(ima)*abs(Ada(ima));
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

return
KK=1:30;
for kmo0=1:3
figure, plot(KK,abs(A(:,kmo0)),'w.',KK,abs(Ares{kmo0}),'--'), pausak
figure, plot(KK,real(A(:,kmo0)),'w.',KK,real(Ares{kmo0}),'--'), pausak
figure, plot(KK,imag(A(:,kmo0)),'w.',KK,imag(Ares{kmo0}),'--'), keyboard
end
