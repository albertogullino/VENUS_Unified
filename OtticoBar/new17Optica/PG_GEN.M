%load pgam
%close all
%iLP=1;
% save pg  Nspla lambda nplt nplut nplit Lplut Lplit Lplt Lplut rfu ...
%           nstratiu rr rqw iLP1 freq
Kp=0;
clear GGev GGmv Tev Tmv

nstv=[0:5:25];
nstv=Nspla.t.m;

%Lplit=[];
%nplit=[];
%Lplut=[];
%Lplut(1:3)=0;
%nplut=[];
%Lplt=[0.069; 0.060];
%nplt=[3.044; 3.5155];
%lambda=2*real(sum(Lplt.*nplt))
%pausak
ireal=1;
if ireal==1
 nplt=real(nplt);
 nplut=real(nplut);
 nplit=real(nplit);
end
mdec=150;
dec=linspace(0,mdec,200);
Lplut0=Lplut;

s=0;
is=1;
lL=length(Lplut0);
while s<mdec
 s=sum(Lplut0(lL-is:lL))*1000
 is=is+1;
% pausak
end
%keyboard


for ist=1:length(nstv)
nstratiu=nstv(ist);

%Iad=length(Lplut)-4;
Iad=length(Lplut)-is;
Lplut=Lplut0;
 for it=1:length(dec)
  dp=dec(it)*1e-3;
  if is==4
    if dp<=Lplut0(4+Iad)
     Lplut(4+Iad)=Lplut0(4+Iad)-dp;
    elseif dp>Lplut0(4+Iad) & dp<=sum(Lplut0([3:4]+Iad))
     Lplut(4+Iad)=0;
     Lplut(3+Iad)=Lplut0(3+Iad)-dp+Lplut0(4+Iad);
    elseif dp>sum(Lplut0([3:4]+Iad))
     Lplut(4+Iad)=0;
     Lplut(3+Iad)=0;
     Lplut(2+Iad)=Lplut0(2+Iad)-dp+sum(Lplut0([3:4]+Iad));
    end
  elseif is==3
    if dp<=Lplut0(3+Iad)
     Lplut(3+Iad)=Lplut0(3+Iad)-dp;
    elseif dp>Lplut0(3+Iad) & dp<=sum(Lplut0([2:3]+Iad))
     Lplut(3+Iad)=0;
     Lplut(2+Iad)=Lplut0(2+Iad)-dp+Lplut0(3+Iad);
    end
  elseif is==2
    if dp<=Lplut0(2+Iad)
     Lplut(2+Iad)=Lplut0(2+Iad)-dp;
    elseif dp>Lplut0(2+Iad) & dp<=sum(Lplut0([1:2]+Iad))
     Lplut(2+Iad)=0;
     Lplut(1+Iad)=Lplut0(1+Iad)-dp+Lplut0(2+Iad);
    end
%    Lplut(lL-is:lL)
%    pausak
  elseif is==1

  end
%  Lplut
%  dp
%  pausak


   [GGe1,GGm1,Te1,Tm1]=gaemms_h(Kp,freq,lambda,Lplit,nplit,Lplt,...
              nplt,nstratiu,rfu,rr,iLP1,Lplut,nplut,real(rqw));

   GGev(it,ist)=GGe1;
%  GGev(:,ist)=GGe1;
%  GGmv(:,ist)=GGm1;
%  Tev(:,ist)=Te1;
%  Tmv(:,ist)=Tm1;
 end
end
figure,
plot(dec,abs(GGev)),
