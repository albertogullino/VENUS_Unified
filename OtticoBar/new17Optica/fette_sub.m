function nfette=fette_sub(R0v,Litot,aitot,nitot,ifp)
ll=length(Litot);
for ke=1:length(R0v)
R0=R0v(ke);
for kf=1:ll
 aid=aitot(kf,:);
 fiai=find(aid>0);
 ai=aid(fiai);
	 fiai=find(ai<R0);
 if length(fiai)>0
  fii=fiai(end)+1;
 else
  fii=1;
 end
nfette(kf,ke)=nitot(kf,fii);
% pausak
end
end

if ifp==-10
h=figure;
col='ymc';
for ke=1:length(R0v)
nto=nfette(:,ke);
nd=[nto nto];
ndplot=reshape(nd.',prod(size(nd)),1);
ddr=Litot;
dd=[zeros(size(ddr))  ddr];
ddplot=cumsum(reshape(dd.',prod(size(nd)),1));
%  plot(ddplot,real(ndplot),[col(ke),'.-']),
  plot(ddplot,real(ndplot)),
  title(' Struttura compattata ')
  if ke==1
   hold on
  end
  pausak

end  
% 'cont nfette', keyboard
end
