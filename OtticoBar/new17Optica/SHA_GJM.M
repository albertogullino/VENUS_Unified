sgimp=1;
igint=0;
ifalso=-1;

X=P.Rx;
Y=P.Ry;
Shap=P.shape;
D_stri=P.D;
d_stri=P.d;
r_add=P.Radd;
if P.orien~0
 orgra=pi/2;
 Xlo=X;
 Ylo=Y;
else
 orgra=0;
 Xlo=Y;
 Ylo=X;
end

 ns=round(2*Xlo/D_stri);
 if shiret==0
  if is_even(ns)==1
   ns=ns+1;
  end
 else
  if is_even(ns)==0
   ns=ns+1;
  end
 end
' parametri '
'[ ns Xlo Ylo D_stri d_stri Shap]'
[ ns Xlo Ylo D_stri d_stri]
if P.shif~0
 shiret=-D_stri/2;
else
 shiret=0;
end

   switch Shap
    case 'circle'
     Rma=Y;
     cas=0;
    case {'square','rectangle'}
     cas=1;
     Cmae=(ns-1)/2*D_stri+d_stri/2;

     if orgra==0
      Rma=sqrt(Cmae^2+X^2);
     else
      Rma=sqrt(Cmae^2+Y^2);
     end
   end

npr=101*fix(Rma/D_stri);
npr=500;
rv=linspace(0,Rma,npr);
pausak
keyboard

if ifp>-3 | ifp==-10 | ifp>=0

  figure
  yv0=d_stri/2*[-1 1];
  npf=50;
  rtot=[];
  rtoti=[];
  for k=1:ns
   Ym=(k-1-fix(ns/2))*D_stri;
   yvd=yv0+Ym;
   aas=yvd/Rma;
   fias=find(abs(aas)<=1);
   fiv=asin(aas(fias));

   if length(fiv)==1
    if fiv<0
     five=linspace(pi-fiv,2*pi+fiv,2*npf);
    else
     five=linspace(fiv,pi-fiv,2*npf);
    end
    rfiu=Rma*exp(j*five);
    fiex=Rma*cos(five([1 end]));
    rfiix=linspace(fiex(1),fiex(2),2*npf);
    rfiiy=ones(size(rfiix))*Rma*sin(five(1));
    rfii=rfiix+j*rfiiy;
    rtoti=[rfii rfiu];
    rtot=[rtot; rtoti];
   else
    five1=linspace(fiv(1),fiv(2),npf);
    five2=linspace(pi-fiv(2),pi-fiv(1),npf);
    rfiu1=Rma*exp(j*five1);
    rfiu2=Rma*exp(j*five2);
    fiex=Rma*cos([five1(end) five2(1)]);
    rfiix=linspace(fiex(1),fiex(2),npf);
    rfiiy=ones(size(rfiix))*Rma*sin(five1(end));
    rfii1=rfiix+j*rfiiy;

    fiex=Rma*cos([five2(end) five1(1)]);
    rfiix=linspace(fiex(1),fiex(2),npf);
    rfiiy=ones(size(rfiix))*Rma*sin(five1(1));
    rfii2=rfiix+j*rfiiy;
    rtoti=[rfiu1 rfii1 rfiu2 rfii2];
    rtot=[rtot; rtoti];
%    figure
%    rplo=rfiu1;
%    plot3(real(rplo),imag(rplo),ones(size(rplo)),'w'), view(2), hold on
%    keyboard

   end
%    plot3(real(rtot.'),imag(rtot.'),ones(size(rtot.')),'w'), view(2), hold on

%   pausak
  end
end

    xtot=real(rtot.');
    ytot=imag(rtot.');
    figure, plot(xtot,ytot)
    if ifp~=-4, pausak, end
muv=[0:2:2*nubesu];

%muv=2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if r_add>0
     firv=find(rv>r_add);
     firv0=1:firv(1)-1;
    else
     firv=2:npr;
     firv0=1;
    end

AB=ones(length(rv),length(muv))*NaN;
im=0;
for mu=muv

im=im+1;

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


   if cas==0
   elseif cas==1
   end

  F_sh(firv0)=fat1;
%  firv=100;
  for k=firv
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
    pausak
   end
  end  %r
  if ifp>0
   ver=4*(F_sh*rv')*diff(rv(1:2))/(pi*max(rv)^2)/(d_stri/D_stri)
   pausak
  end

 else

  fat=0;
  F_sh(1)=fat;
  if r_add>0
   F_sh(firv0)=pi/4;
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
  F_sh=F_sh;
 end
 AB(:,im)=2*F_sh';

% figure, plot(rv,F,'.g')
% pausak
end

ABv=AB;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a1=(d_stri/2)/Rma;
fi1=acos(a1);
ffi1=fi1-0.5*sin(2*fi1);
sv=Rma^2*(pi/2-ffi1);
for kj=1:fix(ns/2)
  a1=(kj*D_stri-d_stri/2)/Rma;
  a2=(kj*D_stri+d_stri/2)/Rma;
  fi1=acos(a1);
  fi2=acos(a2);
  ffi1=fi1-0.5*sin(2*fi1);
  ffi2=fi2-0.5*sin(2*fi2);
  svi=Rma^2*(ffi1-ffi2);
  pausak
  sv=sv+Rma^2*(ffi1-ffi2);
end



firv=find(rv>d_stri/2 & rv<=Rma);
firv0=1:firv(1)-1;
firv1=find(rv>Rma);

AB=ones(length(rv),length(muv))*NaN;
deb=figure;
im=0;
for mu=muv
mu
im=im+1;

 fat=0;
 fat1=pi/2;
 if mu==0
  fat01=pi/2;
 else
%  fat01=pi/4;
  fat01=0;
 end


  F_sh(firv0)=fat01;
  F_sh(firv1)=fat;
%  firv=100;
  for k=firv
   r=rv(k);

   a0=d_stri/(2*r);
   som=0;
   if mu==0
    as0=asin(a0);
    som0=as0;
   else
    as0=asin(a0);
    som0=sin(mu*as0)/mu;
   end

   nstr=fix((r+d_stri/2-shiret)/D_stri);
   for kj=1:nstr
    a1=(kj*D_stri+d_stri/2+shiret)/r;
    a2=(kj*D_stri-d_stri/2+shiret)/r;
    if a2<=1 & a1<=1
     as1=asin(a1);
     as2=asin(a2);
    elseif a2<=1 & a1>1
     as1=fat01;
     as2=asin(a2);
    elseif a2>1 & a1>1
     as1=fat;
     as2=fat;
    end

    if mu==0
     som=som+(as1-as2);
    else
     som=som+(sin(mu*as1)-sin(mu*as2))/mu;
    end

   end   %kj

   F_sh(k)=som+som0;
   if imag(F_sh(k))~=0
    pausak
   end
%   if r>1.13*10
   if r>0.5*10
   xp=linspace(0,2*pi,300);
   rp=r*exp(j*xp);
   figure(deb), plot(real(rp),imag(rp),'w',xtot,ytot,'r'), axis equal, pausak
   end

  end  %r


 AB(:,im)=2*F_sh';

 figure, plot(rv,F_sh,'.g')
 pausak
end

r=rv;
if ifp>1
 if length(muv)>6
  figure, plot(r,AB(:,1:6)), hold on, plot(r,AB(:,7:length(muv)),'.-'),
 else
  figure(fpro), hold on, plot(r,AB), plot(r,ABv,'.'),
 end
 title(' y0, m2, c4, r6, g8, b10, y12, m14, c16, r18, g20, b22')
 pausak
end

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
  end
 end

%a=asav;
%disp(' fine sha_grat '), keyboard
%figure, plot(rv/a*avero,AB)
