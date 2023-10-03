%fiai=find(aitot>4 & aitot<5);
%nired=nitot(:,1:end-1);

clear all
close all
load fette
ll=length(Litot);
R0v=[0 4];
%R0v=[4];
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
% ai
% fiai
% kf
nfette(kf,ke)=nitot(kf,fii);
% pausak
end
end

%figure, plot(cumsum(Litot),nfette(:,1),'.',cumsum(Litot),nfette(:,2)-.1,'.'), pausak
%figure, plot(aitot,cumsum(Litot),'.'), 

h=figure;
col='ymc';
for ke=1:length(R0v)
nto=nfette(:,ke)
nd=[nto nto];
ndplot=reshape(nd.',prod(size(nd)),1);
ddr=Litot;
dd=[zeros(size(ddr))  ddr];
ddplot=cumsum(reshape(dd.',prod(size(nd)),1));
  plot(ddplot,real(ndplot),[col(ke),'.-']),
  if ke==1
   hold on
  end
  pausak
end  
