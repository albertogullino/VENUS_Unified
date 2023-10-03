function [dilu,ailu,nilu,fstu]=lens_rel(th,h,r,N,nl,ifp,Hre,Npair,Rel,Rm_ring)

ife=0;
k0=5;
ideb=0;
iad0=0;

%   [[0; xu(pu,:)' ;0] real(nu(pu,:)')]
%   [[xu(pu,:)'] real(nu(pu,:)')]
%   [[0; xu(pu,:)' ;0] real(nu(pu,:)')]
%  puv=[450 550 700 800 900 1000 1080 1100 1120];
%  for pu=puv
%   [[ xu(pu,:)' ;0] real(nu(pu,:)')]
%   pu
%   figure(1)
%   xa=[1 1]*pu;
%   ya=[0 r+5];
%   hold on
%   plot(xa,ya,'w')
%   axis([350 pu+200 -1 32])
%   pausak
%  end
%

if r<=Rel
 ' lens_rel: r_max<Rel '
 keyboard
end
if r<Rm_ring
 ' lens_rel: r_max<Rm_ring '
 keyboard
end

R=(r^2+h^2)/(2*h);
xd=[0 r];
yd0=[1 1];


Hrei=Hre;
nin=nl;
rin=r;
rmax=r;
hin=h;
ha=h;
thin=th;
Nin=N;
Npa_in=Npair;
fracN=abs(Npair)-floor(abs(Npair));
if fracN~=0
 yfracN=th(1)/1000;
 duth=[th yfracN*1000];
else
 yfracN=0;
 duth=th;
end

  if Npair>0
   no_pair=1;
   no_rep=0;
  else
   per0=sum(th(1:2))/1000;
   y010=sqrt((R+ha).^2-Rel^2)-R;
   y020=sqrt((R+ha).^2-Rm_ring^2)-R;
   yrest=ha-y020;
   Npre=floor((Hrei-yrest)/per0);
   if Npre<2
    Npre=0;
    Had1=Hrei;
   else
    Had1=Hrei-(Npre-1)*per0-(ha-y010);
    Hre=Hre-(Npre-1)*per0;
   end
   Npev=floor(2*ha/per0)+1;
   Had2=2*ha;

    if Npre==0
     Npair=round(Had1/per0)+ round(Had2/per0/2);
%     Npair=ceil(Had1/per0)+ round(Had2/per0/2)+1;
    else
     Npair=round((Had1+Had2/2)/per0);
    end
    if Npair<1
     Npair=1;
    end
  Npair=Npair+fracN;
%  'cont', keyboard
  if floor(abs(Npa_in))*per0+yfracN<=Hrei
   ' Rilievo troppo spesso '
   keyboard
   return
  end
  if Npair+Npre<=abs(Npa_in) & Hrei>0
   Npair=Npair+1;
  end
%  'Npair', keyboard
%  Npair=4.5

    y010=sqrt((R+ha).^2-Rel^2)-R;
    y01tip=sqrt((R+ha).^2-Rel^2)-R+Hrei;
    CS=cumsum([0 th])/1000+ha;
    fii=find(CS-y01tip>0);
    y01=CS(fii);
    Npad=0;
    Npad1=0;
    if Hrei>0
     y01v=y01;
     iy=floor((y01-ha)/per0);
     sco=ha+iy*per0;
     if sco+th(1)/1000>y01
      y01=sco+th(1)/1000;
     else
      y01=sco+per0;
     end
     y01n=y01;
     Npad=1;
     Npad1=1+floor(Hrei/per0)-floor(Hrei/per0-1e-3);
     Npair=Npair+Npad;
    end
%
%    'ceil'
%    keyboard
    yrei=y01-ha;
    yreu=abs(Npa_in)*per0-yrei;
    Npu=floor(yreu/per0)-1;
    ym=h+floor(Npair)*per0+yfracN;
    if Npre==0
     if Hrei>0
      y01=sqrt((R+ha).^2-Rel^2)-R+Hrei;
      if y01<ha
       y01=ha;
      end
      Hs=Hrei
     else
      y01=ha;
      Hs=ha;
     end
%     NM=floor(((abs(Npa_in))*per0-Hs)/per0)-Npad;
%     NM=fix(abs(Npa_in))-Npad;
     NM=floor(((abs(Npa_in))*per0-(Hs-ha))/per0)-Npad1;
%     'NM1', keyboard
     Repet=[1 NM 1];
     yzo=[y01 y01+per0 ym];
     yzo=[y01 y01n ym];
      Ndif=NM+Npair-abs(Npa_in);
      'NM', keyboard
      if Ndif>0 & Hrei>0 & Npair>NM
       NM=NM-Ndif;
       Npair=abs(Npa_in);
       Hre=Hrei;
       ym=h+floor(Npair)*per0+yfracN;
       Repet=[1];
       yzo=[ym];
      end
    else
     NR=Npre;
%     NM=floor((abs(Npa_in)*per0-(Hrei+y010-ha)-ha)/per0);
     NM=round((floor(abs(Npa_in))*per0+yfracN-(Hrei+y010-ha)-ha)/per0);
     NM=floor((floor(abs(Npa_in))*per0+yfracN-Hrei-ha)/per0)-Npad;
%     'NM2', keyboard
     if NM>1
      NM=NM;
      Repet=[1 NR 1 NM 1];
      yzo=[ha ha+per0 y01 y01+per0 ym];
     else
      Repet=[1 NR 1];
      Ndif=NR+Npair-abs(Npa_in);
%      Npair=Npair+1;
%      ym=h+floor(Npair)*per0+yfracN;
      yzo=[ha ha+per0 ym];
%      'NR', keyboard
      if Ndif>0
       NR=NR-Ndif;
       if NR>1
        Npair=Npair+Ndif;
        Hre=Hrei-(NR-1)*per0;
        Repet=[1 NR 1];
        ym=h+floor(Npair)*per0+yfracN;
        yzo=[ha ha+per0 ym];
       else
        Npair=abs(Npa_in);
        Hre=Hrei;
        Repet=[1];
%        'yzo', keyboard
        ym=h+floor(Npair)*per0+yfracN;
        yzo=[ym];
       end
      end
     end
    end

  end
%  yzo=yzo-1e-3;
%  Npair=Npair+fracN;
      ym=h+floor(Npair)*per0+yfracN;
%  ' qui', keyboard

if Npair>0
 nps=2;
 duth=th(1:nps);
 dusa=duth;
 clear duth
 duth=repmat(dusa,1,fix(Npair));
 if Npair-floor(Npair)>0
  duth=[duth duth(1)];
 end
 if length(thin)>nps
  duth=[duth thin(nps+1)];
 end
 th=duth;
 dunr=[nl(1:nps+1) nl(end)];
 dusa0=dunr;
 dusa=dunr(2:end-1);
 clear dunr
 dunr=repmat(dusa,1,fix(Npair));
 if Npair-floor(Npair)>0
  dunr=[dunr dunr(1)];
 end
 if length(thin)>nps
  dunr=[dunr nin(nps+2)];
 end
 dunr=[dusa0(1) dunr dusa0(end)];
 nl=dunr;

end

st=per0/(N+1);
y=0:st:ym;
thm=th/1000;
st1=thm(1)/(N+1);
st2=thm(2)/(N+1);
y01=[st1:st1:thm(1)];
y02=[(st2:st2:thm(2))];
y0=[y01 y01(end)+y02];
st0=h/(N+1);
y=[0:st0:h];
for Np=1:fix(Npair)
 y=[y y(end)+y0];
end
if yfracN>0
 y=[y y(end)+y01];
end
%'y', keyboard
%if y(end)<ym
% y=[y ym];
%end

yi=cumsum([h*1000 duth])/1000;
hrel=h+Hre;

%yi=sort([yi hrel]);
duxr0=[Rel Rm_ring];

ic=0;
for k=1:length(y)
 fn01=0;
 fi01=[];
 dux=sqrt((R+yi).^2-(R+y(k)).^2);
% dux01=sqrt((rdux).^2+(R+y(k)).^2)-R-h;
 rdux=real(dux);
 idux=imag(dux);
 fi00=find(rdux<=rmax & idux==0 & rdux>=0);
 ifi0=0;
 if length(fi00)==0
  ifi0=1;
  fi00d=find(yi-y(k)>0);
  if length(find(yi-y(k)<0))==0 | length(find(yi-y(k)<0))>=length(yi)-1
   fi00=fi00d(1);
  else
   fi00=fi00d(1);
  end
 end
  fi01=fi00;
 if y(k)>=hrel
  dux1=[];
 else
  dux0=sqrt((duxr0).^2+(R+y(k)).^2)-R-h;
  fi0=find(real(dux0)<=Hre & real(dux0)>=0 );
  dux1=duxr0(fi0);
 end
 dux2=dux1;
 if length(dux1)==0
  fi=find(rdux<rmax & idux==0 );
 elseif length(dux1)==1
  if y(k)<h
   fi01=find(rdux<rmax & idux==0 & dux>dux1(1));
  else
   duxr=sqrt((R+hrel).^2-(R+y(k)).^2);
   dux2=sort([dux1 duxr]);
   fi01=find(rdux<rmax & idux==0 & (dux<dux2(1) | dux>dux2(2)));
  end
 else
  fi01=find(rdux<rmax & idux==0 & (dux<dux1(1) | dux>dux1(2)));
 end
 if ifi0==0
  [xad,isov]=sort([dux(fi01) dux2]);
 else
  [xad,isov]=sort([dux2]);
 end
   if ideb==1
    xad
   end
% pausak
% dunrea=real(dunr);
 fi01sa=fi01;
 dunrea=(dunr);
  if length(fi01)>0
   if ideb==1
    '  length(fi01)>0 '
   end
   ist=diff(fi01);
   istf=find(ist>1);
   if length(istf)==1
    if ist(istf)==2
     fir=istf;
     fn01=sort([fi01 fi01(fir)+1]);
     if dux(2)~=duxr0(2)
      fn01=[fn01(1:fir+1) 1 fn01(fir+2:end) fn01(end)+1];
     else
      fn01=[fn01(1:fir+1) 1 fn01(fir+2:end)];
     end
    else
     fir=istf;
     firu=istf+1;
     fn01=sort([fi01 fi01(fir)+1 ]);
     fn01=[fn01(1:fir+1) 1 fn01(fir+2:end) fn01(end)+1];
    end

   else
    if ideb==1
     ' esle 1'
    end
    if length(fi01)==1
     if length(fi00)==1
      fn01=([fi01 fi01(end)+1]);
      if length(dux2)==1
%       'uno'
       fn01=([1 fi01 fi01(end)+1 ]);
      elseif length(dux2)==2
        if ideb==1
         'due'
        end
       fia=find(xad==(dux2(2)));
       if fia==2
         pu=min([fia-1 length(fi01)]);
        if ifi0==0
         fn01=([fi01(1:pu) 1 fi01(pu:end) fi01(end)+1]);
        else
         fn01=([fi01(1:pu) 1 fi01(pu:end)]);
        end
       else
        fn01=([fi01 fi01+1 1 fi01+1]);
       end
%       pausak
      else
%       'tre'
       fn01=([fi01 fi01(end)+1]);
      end
     else
      fn01=([fi01 fi01(end)+1]);
      if length(dux2)==1
       fn01=([1 fi01 fi01(end)+1 ]);
      elseif length(dux2)==2
       fia=find(xad==(dux2(2)));
        if ideb==1
         ' sono qui?', %pausak
        end
        if ifi0==0
         pu=fia-1;
         fn01=([fi00(1:pu) 1 fi00(pu+1:end) fi01(end)+1]);
        else
         fn01=([fi01(1:pu) 1 fi01(pu:end)]);
        end
      else
       fn01=([fi01 fi01(end)+1]);
      end
     end
    else
     if ideb==1
      'dopo else fi01 ',
     end
     if length(dux2)==1
      fn01=([1 fi01 fi01(end)+1 ]);
     elseif length(dux2)==2
      fia=find(xad==(dux2(2)));
       if ideb==1
        ' dove pasa',
       end
        if ifi0==0
         pu=fia-1;
         fn01=([fi00(1:pu) 1 fi00(pu:end) fi01(end)+1]);
        else
         fn01=([fi01(1:pu) 1 fi01(pu:end)]);
        end
     else
      fn01=([fi01 fi01(end)+1]);
     end
    end
   end
   dunrea01=dunrea(fn01);

   if fn01(1)>1
    dudu=[dunrea01  dunrea(1)];
   else
    dudu=[dunrea01  dunrea(1)];
   end
   nad=dudu;
  else

    if ideb==1
     '  length(fi01)=0 '
    end
   fi01s=fi01;
   if length(fi00)>0
    fi01=fi00(1);
    if fi01(1)==1
     fn01=([fi01 fi01(end)+1]);
    else
     if length(dux2)==1
      fn01=([1 fi01 fi01(end)+1 ]);
     elseif length(dux2)==2
      fia=find(xad==(dux2(2)));
      pu=fia-1;
      fn01=([fi01(1:pu) 1  fi01(end)+1]);
        if ifi0==0
         pu=min([fia-1 length(fi01)]);
         fn01=([fi01(1:pu) 1 fi01(end)+1]);
        else
         fn01=([fi01(1:pu) 1 fi01(pu:end)]);
        end
     else
      fn01=([fi01 fi01(end)+1]);
     end
    end
    dunrea01=dunrea(fn01);

    if fn01(1)>1
     dudu=[dunrea01  dunrea(1)];
    else
     dudu=[dunrea01  dunrea(1)];
    end
    nad=dudu;
   else
    nad=[];
   end
  end
 if length(xad)>0
   if fn01(end)==length(dunr)
    ruad=[];
    nad=nad(1:end-1);
   else
    ruad=r;
   end
   if iad0>1
    ic=ic+1;
    yu(ic)=y(k-1);
%    ' qui ', keyboard
    xu(ic,1)=r ;
    du=(y(k-1)-ha);
    nper=floor(du/per0);
    drp=du-nper*per0;
    if drp*1000<th(1)
     nad0=dunr(2);
    else
     nad0=dunr(3);
    end
%    nad0
%    pausak
    nu(ic,1:2)=[nad0 nl(end)];
   end
   ic=ic+1;
   xad=[xad ruad];
   isa=0;
   if xad(1)==0
    if length(xad)>1
     xad=xad(2:end);
     nad=nad(2:end);
    else
     isa=1;
%     ' qui', keyboard
%     xad=xu(ic-1,:)/2;
%     nad=nu(ic-1,:);
%     yu(ic)=yu(ic-1)+diff(yu(ic-1:ic))*2/3;
%     yu(ic)=yu(ic-1)+diff(yu(ic-1:ic));
    end
   end
   if isa==0
   yu(ic)=y(k);
   xu(ic,1:length(xad))=xad ;
   nu(ic,1:length(nad))=nad;
   end
   iad0=0;
 else

   iad0=iad0+1;
%   ' qui ', keyboard
 end
  if ic-1==k0 & ife==1
   'k ferma',
   ic-1
   keyboard
  end
end
if ifp==-10
  figure, plot(xu,'.'), grid
% ' in sub_rel '
% pausak

 figure,
 plot(xu,yu,'.'), grid
 for k=1:length(yzo)
  hold on, plot(xd,yd0*yzo(k),'w')
 end
 pausak
end

yud=[yu];
dilu=diff(yud)';
dilu=[dilu; dilu(end)];
fstu=zeros(length(dilu),2);
yinf=0;

for k=1:length(yzo)
 ysup=yzo(k);
 fi=find(yu>=yinf & yu<ysup);
% fi=find(yu>yinf & yu<=ysup);
 if k==length(yzo)
  fi=find(yu>=yinf & yu<=ysup);
 end

 fstu(fi,2)=Repet(k);
 fstu(fi,1)=length(fi);
% [k ysup], pausak
 yinf=ysup;
end


ailu=(xu);
nilu=(nu);

sfst=sum(fstu,2);
itro=find(diff(sfst)~=0);
itro=[itro; length(fstu)];
ki=1;
xd=[];
yd=[];
nd=[];
dsu=cumsum(dilu);
sup=0;
for k=1:length(itro)
 ku=itro(k);
 pu=ki:ku;
 ydu=dsu(pu);
 xdu=ailu(pu,:);
 ndu=real(nilu(pu,1));
 RP=fstu(ki,2);
 for kk=1:RP
  if kk>1
   sup=sup+per0;
  end
  xd=[xd; xdu];
  yd=[yd; ydu+sup];
  nd=[nd; ndu];
 end
% if kk<RP
%  sup=sup+per01;
%  xd=[xd; xdu];
%  yd=[yd; ydu+sup];
% end
 ki=ku+1;
end

y0=[0 h];
thm=th/1000;
nre=real(nl);
n0=nre(1)*[1 1];
for k=1:fix(abs(Npa_in))
 y0=[y0 y0(end)+[0 thm(1)]];
 n0=[n0 nre(2)*[1 1]];
 y0=[y0 y0(end)+[0 thm(2)]];
 n0=[n0 nre(3)*[1 1]];
end
if yfracN~=0
 y0=[y0 y0(end)+[0 thm(1)]];
 n0=[n0 nre(2)*[1 1]];
end


if ifp==-10
 figure
 plot(xd,yd,'.'), grid, pausak
 figure
 plot(yd,nd,'.',y0,n0), grid
 pausak
 figure, plot(fstu(:,2),'.')
 pausak

end






dilu=[dilu(1); dilu];
nilu=[nilu(1,1)*ones(size(nilu(1,:))); nilu];
ailu=[ailu(1,:)*0; ailu];
fstu=[[0 1]; fstu];
ailu(end,:)=ailu(end,:)/2;

if ifp~=-4
 ' prima di return', keyboard
end

ailu=flipud(ailu);
nilu=flipud(nilu);
dilu=flipud(dilu)*1000;
fstu=flipud(fstu);



return




%keyboard
  pin=1
  pa=1
  puv=[pin:pa: length(yu)];
  for pu=puv
   [[ xu(pu,:)' ;0] real(nu(pu,:)')]
%   [[ xu(pu,:)'] real(nu(pu,:)')]
   pu
   figure(1)
   xa=[1 1]*pu;
   ya=[0 r+5];
   hold on
   plot(xa,ya,'w')
   axis([pu-10 pu+20 -1 32])
   pausak
  end

keyboard

return
for k=1:1:2000,
fi=find(nu(k,:)~=0);
plot(xu(k,fi(1:end-1)),nu(k,fi(1:end-1)),'.'), pausak,
end

