G=gain;
F=delf';
 [Ftot,iF]=sort([F]);
 Gtot=G(iF);
 gmi=min(Gtot);
 iT=find(Gtot<30*gmi & Ftot<4.5e-3);
 Fto=Ftot(iT);
 Gto=Gtot(iT);
 figure, semilogy(Ftot,Gtot,'ro',Fto,Gto,'g.'), pausak
 D=diff(Fto);
 Ds0=sort(D);
 fact=3e5/lambda;
 iD=find(Ds0*fact>0);
 Ds=Ds0(iD);
% iT=iT(1:length(iT)-4);
% Ds=Ds(iT);
 maDs=max(Ds);

 nist=fix(length(Ds)/6);
b=figure;
%bn=figure;
for nist=5:40
 nist
 clear ps
 pist=linspace(0,maDs,nist);
 disp(' binsize GHz' )
 bz=diff(pist(1:2))*fact
 if bz<6, break, end
 for kp=1:nist-1
  ps(kp)=length(find(Ds>pist(kp) & Ds<=pist(kp+1)));
 end
 ps=[0 ps];
 s=pist/mean(Ds);
 area=sum(ps)*diff(s(1:2));
 figure(b)
 subplot(211)
 bar(s,ps/area)
 sw=s;
 W=pi/2*sw.*exp(-pi/4*sw.^2);
 We=exp(-sw);
 hold on
 plot(s,W,'r',s,We,'g')
   sym=['binsize = ' fcharc(bz,3) ' (GHz)'];
   ht=text(0.5*max(s),.9,sym); set(ht,'fontsize',12);
   hold off
 subplot(212)
 bar(s,ps)
%   sym=['binsize = ' fcharc(bz,3) ' (GHz)'];
%   ht=text(0.5*max(s),.9*max(ps),sym); set(ht,'fontsize',12);
 pausak
end
