
%' inizio sha_grat'
%keyboard
tc1=clock;
global icambio_ref iNIexch

sgimp=1;
igint=0;
ifalso=-1;
clear A B AB F_sh

X=P.Rx;
Y=P.Ry;
Shap=P.shape;
D_stri=P.D;
d_stri=P.d;
r_add=P.Radd;
Rext=P.Rext;
if isempty(Rext)==1
 Rext=0;
end
if Rext<0
 Rext=0;
end
if Rext>0
 iNIexch=1;
end

if P.orien~0
 orgra=pi/2;
 Xlo=X;
 Ylo=Y;
else
 orgra=0;
 Xlo=Y;
 Ylo=X;
end
i_grap=P.grap;

global igr_app
if length(igr_app)==1
% if igr_app==1
  i_grap=igr_app;
% end
end


if iNIexch==1
 nitot(ind_vet-1,1:2)=nitot(ind_vet-1,[2 1]);
 nv(ind_vet,1:2)=nv(ind_vet,[2 1]);
 nitn=nv(put,:).';
 nib =nv(pub,:).';
 d_stri=P.D-P.d;
end

%' grat ', keyboard
ireturn=0;
 if ~exist('ifirstgrat')
  ifirstgrat=0;
 end
if i_grap>=1
 if iany==0 & igr_app==1
  ifirstgrat=1;
 end
 if ifp~=-4
%  ' SHAGRAT in', keyboard
 end
 t1=P.d;
 t2=P.D-P.d;
 tt=P.D;
 rad=radii.a;
 er1=nv(ind_vet,1)^2;
 er2=nv(ind_vet,2)^2;
 n_ve=tt*er1*er2/(t2*er1+t1*er2);
 n_pa=(t1*er1+t2*er2)/tt;
 n_me=((n_ve+n_pa)/2);
 n_di=((n_pa-n_ve)/2);
 if P.orien==0
  ey=n_ve;
  ex=n_pa;
 else
  ex=n_ve;
  ey=n_pa;
 end
 Delta_eps=(ex-ey)/2;


% ' igrapp', keyboard

 if igr_app==1
% ' igrapp', keyboard
  fi6=find(shailoop~=6);
  shailoopu=shailoop(fi6);
  sha=0;
  rad(ind_vet,:)=0;
  aitot(ind_vet,:)=0;
  nv(ind_vet,:)=0;
  if iLP==0
   nv(ind_vet,1)=sqrt(n_me);
%   anyf(ind_vet,:)=1;
%    nv(ind_vet,1)=sqrt(ex);
%    nv(ind_vet,1)=sqrt(ey);
  else
   if icambio_ref==1
    nv(ind_vet,1)=sqrt(ex);
   else
    nv(ind_vet,1)=sqrt(ey);
   end
  end
 elseif igr_app==2
  if iLP==0
   nv(ind_vet,1)=sqrt(n_me);
%'wuo', keyboard
  else
   if icambio_ref==1
    nv(ind_vet,1)=sqrt(ex);
   else
    nv(ind_vet,1)=sqrt(ey);
   end
  end
  sha=1;
  shag=sha;
  lxi=1;
  rgra=P.Ry;
  adis=rgra*kcav;
  rad(ind_vet,:)=rgra;
  aitot(ind_vet,:)=rgra;
%'wuo  las' , keyboard
%  nv(ind_vet,2)=sqrt(max([er1 er2]));
  nv(ind_vet,2)=sqrt(er2);

%    nv(ind_vet,1)=sqrt(ex);
%    nv(ind_vet,1)=sqrt(ey);

 end
 
 shavet(ind_vet)=sha;


 nitn=nv(put,:).';
 nib =nv(pub,:).';
 shav.t=shavet(put,:);
 shav.b=shavet(pub,:);
 if iLP==0
  anyret(ind_vet)=2*Delta_eps;
 end
 anyr.t=anyret(put,:);
 anyr.b=anyret(pub,:);
 aitn=rad(put,:).';
 aib=rad(pub,:).';

% bi=aral.y.t';
% pai=aral.p.t';
%'prima di ireturn =1 ', keyboard

 ireturn=1;
 if ifp~=-4
  ' SHAGRAT', keyboard
 end
 return
end


ilogra=0
if ilogra==1
 load retf
 return
 % save retf rv AB A B
end 

if Rext>Rma
 Rvero=Rext;
 shiff=1;
else
 shiff=P.shif;
 Rvero=Rma;
end

if shiff~0
 shiret=-D_stri/2;
else
 shiret=0;
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


%npr=201*ceil(Rma/D_stri);
%rv=linspace(0,Rma*1.01,npr);

npr=201*fix(Rvero/D_stri);
rv=linspace(0,Rvero,npr);
if Rext>0
 fireti=find(rv<=Rma);
 rv_ret=rv(fireti);
else
 rv_ret=rv;
end

%' rv ', keyboard

if ifp>-3 | ifp==-10 | ifp>=0

  figure

iold=0;
 if iold==1
%  yv0=[-d_stri -d_stri d_stri d_stri]/2;
  yv0=[-d_stri*ones(1,100) d_stri*ones(1,100)]/2;
  for k=1:ns
   Ym=(k-1-fix(ns/2))*D_stri-shiret;
   yvd=yv0+Ym;
    if cas==0
     dif=Xlo^2-(Ym)^2;
    elseif cas==1
     dif=Ylo^2;
    end
   if dif>0
    Xm=sqrt(dif);
    xpar=linspace(-Xm,Xm,100);
    xvd=[xpar fliplr(xpar)];
    if orgra==0
     yv=yvd;
     xv=xvd;
    else
     xv=yvd;
     yv=xvd;
    end
    if r_add>0
%    keyboard
     fiad=find(xv.^2+yv.^2<r_add^2);
     xv(fiad)=NaN;
     yv(fiad)=NaN;
    end
    plot3(xv,yv,ones(size(xv)),'w'), view(2), hold on
   end
  end

 else


%npr=201*fix(Rma/D_stri);
%rv=linspace(0,Rma,npr);


  yv0=d_stri/2*[-1 1];
  npf=50;
  rtot=[];
  rtoti=[];
  for k=1:ns
   Ym=(k-1-fix(ns/2))*D_stri-shiret;
   yvd=yv0+Ym;
   aas=yvd/Rma;
   fias=find(abs(aas)<=1);
   fiv=asin(aas(fias));
   if length(fiv)~=0
    if length(fiv)==1
     if fiv<0
      five=linspace(pi-fiv,2*pi+fiv,2*npf);
     else
      five=linspace(fiv,pi-fiv,2*npf);
     end
     five=five+orgra;
     rfiu=Rma*exp(j*(five));
     fiex=Rma*cos(five([1 end]));
     rfiix=linspace(fiex(1),fiex(2),2*npf);
     rfiiy=ones(size(rfiix))*Rma*sin(five(1));
     rfii=rfiix+j*rfiiy;
     rtoti=[rfii rfiu];
     rtot=[rtot; rtoti];
    else
     five1=linspace(fiv(1),fiv(2),npf)+orgra;
     five2=linspace(pi-fiv(2),pi-fiv(1),npf)+orgra;
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
    end
   end
  end
    if r_add>0
     fiad=find(abs(rtot)<r_add);
     rtot(fiad)=NaN;
    end

    plot3(real(rtot.'),imag(rtot.'),ones(size(rtot.')),'w'), view(2), hold on

 end

    if r_add>0
     fi=linspace(0,2*pi,155);
     cirv=r_add*exp(j*fi);
     plot3(real(cirv),imag(cirv),ones(size(cirv)),'w'), view(2)
    end

%  npr=1001;

  hold on,
   switch Shap
    case 'circle'
     polar(linspace(0,2*pi,100),max(rv)*ones(1,100),'r')
     if Rext>0
      polar(linspace(0,2*pi,100),Rma*ones(1,100),'r')
      title(' the areas within white and red lines is etched away ')
     end
    case {'square','rectangle'}
     if orgra==0
      xplo=[X X -X -X X];
      yplo=[-1 1 1 -1 -1]*Cmae;
     else
      xplo=[1 1 -1 -1 1]*Cmae;
      yplo=[-Y Y Y -Y -Y];
     end
     plot(xplo,yplo,'r')
   end
  axis equal
  drawnow
  if ifp>-3
   pausak
  end
%'cirv in sha_grat'  , keyboard
end

%if ifp~=-4
etis1=etime(clock,tc1)
%' verifica dove sono ', keyboard
%end
%muv=[0:2:2*(nubesu+1)];
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
%mu
%keyboard
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
%   if k==1000
%   etis3=etime(clock,tc1)
%   pausak
%   end
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
