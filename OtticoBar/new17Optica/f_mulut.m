function [fme,fmm]=f_mulut(z)


global fsto L_i n_i rr rfd rfu  iLP ifp iff Lam0 ifun i12 iorta ibast

s=size(z);

zu=1/rfd;
zi=1/rfu;
zr=1/rr;
rmi=zi/zr;
rmu=zu/zr;


%if i12==1
 ra=1/rmi;
 rs=rmu;
 rd=rfd;
%else
% ra=1/rmu;
% rs=rmi;
% rd=rfu;
%end

Vi=1;

Ii=-1*ra;
yi=[Vi; Ii];


Mva=0.5*sqrt(rd)*[1 rs; 1 -rs];

%' quio muly', keyboard
for k1=1:s(1)
for k2=1:s(2)

lambdai=Lam0+real(z(k1,k2));
nim=imag(z(k1,k2));
iel=1;
freq=0;
fiqw=find(fsto==-1);
ishow=0;
if iorta==0
 [Mte,Mtm]=gaz_mtu(fiqw,nim,L_i,n_i,rr,rfd,rfu,lambdai,freq,0,iLP,ifp,iff);
else 
% [Mte,Mtm]=gaz_mtor(fiqw,nim,L_i,n_i,rr,rfd,rfu,lambdai,freq,0,iLP,ifp,iff);
 [Mte,Mtm]=gaz_mtoru(fiqw,nim,L_i,n_i,rr,rfd,rfu,lambdai,freq,0,iLP,ifp,iff,ibast);
end
yu=Mte*yi;
fu=Mva*yu;
fme(k1,k2)=fu(2);

yu=Mtm*yi;
fu=Mva*yu;
fmm(k1,k2)=fu(2);
% 'dentro f_mulut', keyboard
end
end

%M=inv(Mt0);
%fmu1=zu*(M(1,1)+zi*M(2,1))-(M(1,2)+zi*M(2,2));


%fmu=zr*(M(1,1)+zi/zr*M(2,1))-zu*(M(1,2)-zi/zr*M(2,2));
if ifun==3
 'dentro f_mulu', keyboard
end

%'in fmu', keyboard
%
%'in f_mul1', keyboard