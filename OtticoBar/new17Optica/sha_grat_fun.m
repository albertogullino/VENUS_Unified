function [A,B,AB]=sha_grat_fun(PAR_sha,ifp)

tc1=clock;

nubesu=PAR_sha.nubesu;
r_add=PAR_sha.r_add;
Rma=PAR_sha.Rma;
rv_ret=PAR_sha.rv_ret;
rv=PAR_sha.rv;
shiret=PAR_sha.shiret;
cas=PAR_sha.cas;
d_stri=PAR_sha.d_stri;
D_stri=PAR_sha.D_stri;
X=PAR_sha.X;
Y=PAR_sha.Y;
Ylo=PAR_sha.Ylo;
orgra=PAR_sha.orgra;
kcav0=PAR_sha.kcav0;
pimu=PAR_sha.pimu;
meun=PAR_sha.meun;
mbv=PAR_sha.mbv;
iclo=PAR_sha.iclo;


muv=[0:2:2*(nubesu)];

    if r_add>0
     firv=find(rv_ret>r_add & rv_ret<=Rma);
     firv0=1:firv(1)-1;
    else
     firv=find(rv_ret>0 & rv_ret<=Rma);
     firv0=1;
    end

%AB=ones(length(rv),length(muv))*NaN;
AB=zeros(length(rv),length(muv));
AB(:,1)=pi;
im=0;
for mu=muv

im=im+1;
%mu, keyboard
 if mu==0

  fat=pi/2;
  if r_add==0
   if shiret==0
    fat1=fat;
   else
    fat1=0;
   end
  else
    fat1=fat;
  end



  F_sh(firv0)=fat1;
  for k=firv
%k
%tic
   r=rv(k);
   if cas==0
    a0=d_stri/(2*r);
   elseif cas==1
    H1=sqrt(X^2+(d_stri/2)^2);
    if r<=X
     a0=d_stri/(2*r);
    elseif r>X & r<=H1
     H=sqrt(r^2-Ylo^2);
     a0=H/r;
    elseif r>H1
     a0=0;
    end
   end
   nstr=fix((r+d_stri/2-shiret)/D_stri);
   som=0;
   som0=0;
   if  shiret==0
    if nstr==0 & shiret==0
     if a0<=1
      as0=asin(a0);
     else
      as0=fat;
     end
     som0=as0;
    else
     som0=asin(a0);
    end
   end
   for kj=1:nstr
    if cas==0
     a1=(kj*D_stri+d_stri/2+shiret)/r;
    elseif cas==1
     Y2=kj*D_stri+d_stri/2+shiret;
     H2=sqrt(X^2+(Y2)^2);
     if r<=H2
      a1=(kj*D_stri+d_stri/2+shiret)/r;
     elseif r>H1 & r<=H2
 %     H=sqrt(r^2-H1^2);
      H=sqrt(r^2-Ylo^2);
      a1=H/r;
     elseif r>H2
      a1=0;
     end
    end
    if a1<=1
     as1=asin(a1);
    else
     as1=fat;
    end

    if cas==0
     a2=(kj*D_stri-d_stri/2+shiret)/r;
    elseif cas==1
     Y1=kj*D_stri-d_stri/2+shiret;
     Y2=kj*D_stri+d_stri/2+shiret;
     H1=sqrt(X^2+(Y1)^2);
     H2=sqrt(X^2+(Y2)^2);
     if r<=H1
      a2=(kj*D_stri-d_stri/2+shiret)/r;
     elseif r>H1 & r<=H2
%      H=sqrt(r^2-H1^2);
      H=sqrt(r^2-Ylo^2);
      a2=H/r;
     elseif r>H2
      a2=0;
     end
    end
    as2=asin(a2);
    som=som+(as1-as2);
   end

   F_sh(k)=som+som0;
   if imag(F_sh(k))~=0
%    pausak
   end
%  toc
%  pausak
  end  %r
  if ifp>0
%   ver=4*(F_sh*rv')*diff(rv(1:2))/(pi*max(rv)^2)/(d_stri/D_stri)
%   pausak
  end
  F=F_sh;
 else

  fat=0;
  F_sh(1)=fat;
  if r_add>0
   F_sh(firv0)=pi/4; %ERRORE
   F_sh(firv0)=0;
  end
  for k=firv
   r=rv(k);
   if cas==0
    a0=d_stri/(2*r);
   elseif cas==1
    H1=sqrt(X^2+(d_stri/2)^2);
    if r<=X
     a0=d_stri/(2*r);
    elseif r>X & r<=H1
%     H=sqrt(r^2-(d_stri/2)^2);
     H=sqrt(r^2-Ylo^2);
     a0=H/r;
    elseif r>H1
     a0=0;
    end
   end
   nstr=fix((r+d_stri/2-shiret)/D_stri);
   som=0;
   som0=0;

   if shiret==0
    if nstr==0
     if a0<=1
      as0=asin(a0);
     else
      as0=fat;
     end
     som0=sin(mu*(as0+orgra))/mu;
    else
     as0=asin(a0);
     som0=sin(mu*(as0+orgra))/mu;
    end
   end

   for kj=1:nstr
    if cas==0
     a1=(kj*D_stri+d_stri/2+shiret)/r;
    elseif cas==1
     Y2=kj*D_stri+d_stri/2+shiret;
     H2=sqrt(X^2+(Y2)^2);
     if r<=H2
      a1=(kj*D_stri+d_stri/2+shiret)/r;
     elseif r>H1 & r<=H2
%      H=sqrt(r^2-H1^2);
      H=sqrt(r^2-Ylo^2);
      a1=H/r;
     elseif r>H2
      a1=0;
     end
    end

    if a1<=1
     as1=asin(a1);
    else
     as1=fat;
    end
    if cas==0
     a2=(kj*D_stri-d_stri/2+shiret)/r;
    elseif cas==1
     Y1=kj*D_stri-d_stri/2+shiret;
     Y2=kj*D_stri+d_stri/2+shiret;
     H1=sqrt(X^2+(Y1)^2);
     H2=sqrt(X^2+(Y2)^2);
     if r<=H1
      a2=(kj*D_stri-d_stri/2+shiret)/r;
     elseif r>H1 & r<=H2
%      H=sqrt(r^2-H1^2);
      H=sqrt(r^2-Ylo^2);
      a2=H/r;
     elseif r>H2
      a2=0;
     end
    end
    as2=asin(a2);
    som=som+(sin(mu*(as1+orgra))-sin(mu*(as2+orgra)))/mu;
   end
   F_sh(k)=som+som0;
  end  %r
%  F_sh=F_sh;
  F=F_sh;
 end
 AB(1:firv(end),im)=2*F_sh';
 clear F_sh
%etis2=etime(clock,tc1)
 if ifp>=2
  im
 end
%pausak

% figure, plot(rv,F_sh,'.g')
% pausak
 if ifp>=2
  if exist('F')
   figure, plot(rv,F,'.g')
   pausak
  end
 end
 clear F
end %im

r=rv;
if ifp>1
 if length(muv)>6
  figure, plot(r,AB(:,1:6)), hold on, plot(r,AB(:,7:length(muv)),'.-'),
 else
  figure, plot(r,AB),
 end
 title(' y0, m2, c4, r6, g8, b10, y12, m14, c16, r18, g20, b22')
 pausak
end
r=rv*kcav0;
rv=r;
%' dopo primo passo', keyboard
  if ifp>10
   figure
  end
 for imu=pimu
  jmu=imu-meun;
  mu=mbv(imu);
  for inu=pimu
   jnu=inu-meun;
   nu=mbv(inu);
   if (nu+mu)/2-fix((nu+mu)/2)==0
    dmn=abs(mu-nu);
    fim=find(dmn==muv);
    if length(fim)==1
     mfatd=AB(:,fim);
     A(:,jmu,jnu)=AB(:,fim);
    else
     disp('errore A mu in sha_grat ')
     pausak
    end
    dmn=abs(mu+nu);
    fim=find(dmn==muv);
    if length(fim)==1
     mfats=AB(:,fim);
     B(:,jmu,jnu)=AB(:,fim);
    else
     disp('errore B mu in sha_grat ')
     pausak
    end
   end  %if
%   [mu nu], pausak
  end
 end

 if ifp>=2
  ' fine ahgrat', keyboard
 end
%a=asav;
%disp(' fine sha_grat '), keyboard
%figure, plot(rv/a*avero,AB), keyboard
etis=etime(clock,tc1)
%if ifp~=-4
' fine sha_grat'
   if iclo==1
    close
   end
% keyboard
%end
