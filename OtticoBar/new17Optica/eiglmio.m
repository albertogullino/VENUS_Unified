function [Ksi,lamres,Fi,uLong,uLong0,Fa]=eiglmio(lambda,uFunc,uF0,relPerm,hn)

relStor=length(relPerm);
h=1e-9;                                                                   % anger avståndet mellan två punkter (steglängden)
%h=1e-9*.885;                                                                   % anger avståndet mellan två punkter (steglängden)
%h=1e-9*.9;                                                                   % anger avståndet mellan två punkter (steglängden)
k0=2*pi/lambda;
nu=sqrt(relPerm(1));
A0=spalloc(relStor,relStor,relStor*3+3);  % Initierar egenvÄrdesmatrisen A ( Ax=lx)
lillaH=lambda^2/((2*pi)^2*h^2);                         % En konstant som ingår i fin. differens matrisen





%%%%    Bildar de tre diagonal elementen under-,over- & huvuddiagonal som fås vid finit differentiuering av
%%%%    ekv. 8 i artikeln

                A0(1,2)=lillaH/(relPerm(1))*(1+exp(+i*2*pi*sqrt(relPerm(1))/lambda*2*h));
                A0(1,1)=-2*lillaH/(relPerm(1))+1;

        for k=2:(relStor-1)

                A0(k,k+1)=lillaH/(relPerm(k));
                A0(k,k-1)=lillaH/(relPerm(k));
                A0(k,k)=-2*lillaH/(relPerm(k))+1;

        end

                A0(relStor,relStor-1)=lillaH/(relPerm(relStor))*(1+exp(+i*2*pi*sqrt(relPerm(relStor))/lambda*2*h));
                A0(relStor,relStor)=-2*lillaH/(relPerm(relStor))+1;



sigma.disp=0;
%disp(' longitudinal resonance evaluation ')
%keyboard
[eiVectors,Ksi,pppk]=eigs(A0,1,'SM',sigma);


totM=size(relPerm);
tot=max(totM);
stop=tot;
nu=sqrt(relPerm(1));
nl=sqrt(relPerm(end));
n=floor((tot-1)/2);
a=0;
b=tot;

upperInt=conj(eiVectors).*((relPerm).').*eiVectors;
lowInt=conj(eiVectors).*eiVectors;

taljare=h/2*(upperInt(1)+2*upperInt(2)+2*sum(upperInt(3:end-2))+...
   2*upperInt(end-1)+upperInt(end));

namnare=h/2*(lowInt(1)+2*lowInt(2)+2*sum(lowInt(3:end-2))+...
   2*lowInt(end-1)+lowInt(end));

eps=taljare/namnare;


VectKsi=nu*abs(eiVectors(1))^2/(k0*real(eps)*namnare);

%upperInt=conj(eiVectors).*(uFunc).'.*eiVectors;
%lowInt=conj(eiVectors).*eiVectors;

%upperInt=conj(eiVectors).*(uFunc).'.*eiVectors;
upperInt=conj(eiVectors).*eiVectors;
%lowInt=conj(eiVectors).*eiVectors.*relPerm.'/11.58;
lowInt=conj(eiVectors).*eiVectors;
WlowInt=conj(eiVectors).*eiVectors.*real(relPerm.');

fi=find(uFunc~=0);
fi0=find(uF0~=0);

%x=h*fi0;
%y=upperInt(fi0);
%xqw=h*linspace(fi0(1),fi0(end),50);
%yqw=spline(x,y,xqw);
%Ip=diff(xqw)*yqw(1:end-1)';
NQW=length(fi)/length(fi0);
dqw=h*(length(fi0)-1);
taljareU=NQW*dqw*mean(upperInt(fi));
taljareUqw=dqw*mean(upperInt(fi0));

WtaljareU=NQW*dqw*mean(WlowInt(fi));
WtaljareUqw=dqw*mean(WlowInt(fi0));

%taljareU=h/2*(upperInt(1)+2*upperInt(2)+2*sum(upperInt(3:end-2))+...
%   2*upperInt(end-1)+upperInt(end));

namnareU=h/2*(lowInt(1)+2*lowInt(2)+2*sum(lowInt(3:end-2))+...
   2*lowInt(end-1)+lowInt(end));

namnareW=h/2*(WlowInt(1)+2*WlowInt(2)+2*sum(WlowInt(3:end-2))+...
   2*WlowInt(end-1)+WlowInt(end));

%upperIntqw=conj(eiVectors).*(uF0).'.*eiVectors;

%taljareUqw=h/2*(upperIntqw(1)+2*upperIntqw(2)+2*sum(upperIntqw(3:end-2))+...
%   2*upperIntqw(end-1)+upperIntqw(end));

% vecchi
OuLong=taljareU/namnareU;
OuLong0=taljareUqw/namnareU;

uLong=WtaljareU/namnareW;
uLong0=WtaljareUqw/namnareW;



F=abs(eiVectors).^2;

Fi=3*F/max(F);
%' eigmio', keyboard
lamres=lambda/real(sqrt((1-(Ksi))));
%disp('eigmio')
%F=real(eiVectors);
F=abs(eiVectors);
Fa=3*F/max(F);
%keyboard
