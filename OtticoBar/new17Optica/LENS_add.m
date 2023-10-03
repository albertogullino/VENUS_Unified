function [dilu,ailu,nilu,fstu]=lens_add(th_in,h,r,N,nl,ifp,Npair,Rflat)
%' entro add', keyboard
 if h==0
  Ni=abs(Npair);
  dN=Ni-floor(Ni);
  Nii=Ni-dN;
  dilu=th_in'; 
  nilu=nl(2:end-1)'; 
  nilu=[nilu nilu*0];
  ailu=zeros(size(dilu)); 
  fstu=zeros(length(dilu),2);
  fstu(:,1)=length(dilu);
  fstu(:,2)=Nii;
  if dN>0
   dilu=[dilu; dilu(1)];
   nilu=[nilu; nilu(1,:)];
   ailu=[ailu; ailu(1)];
   fstu=[fstu; [0 1]];
  end
  return
 end

th=th_in;

Rel=0;
Rm_ring=0;
Hre=0;

thmi=abs(th/1000);
ive=1;
ife=0;
k0=1;

ideb=0;
iad0=0;
k0=k0-1;

%h=R-sqrt(R^2-r^2)
R=(r^2+h^2)/(2*h);

'calcolato R', keyboard
xd=[0 r];
yd0=[1 1];


Hrei=Hre;

rin=r;
rmax=r;
hin=h;
ha=h;
fi=find(th>=0);
fimag=fi;
thin=th(fi);
thmi=abs(thin/1000);
nin=nl(fi+1);
fi=find(th<0);
fimin=fi;
thinu=abs(th(fi));
ninu=nl(fi+1);
suh=0;
nia=[];
tha=[];
if length(fimin)>0
hnm=h*1000;
while suh<hnm
 for kh=1:length(thinu)
  suh=suh+thinu(kh);
  tha=[tha thinu(kh)];
  nia=[nia ninu(kh)];
  if suh>hnm
   break
  end
 end
end
end %lenght
%'ninu', keyboard


Nin=N;
Npa_in=Npair;
fracN=abs(Npair)-floor(abs(Npair));
if fracN~=0
 yfracN=thmi(1);
 duth=[th yfracN*1000];
else
 yfracN=0;
 duth=th;
end
%'yfracN', yfracN
    thmit=repmat(thmi,1,floor(abs(Npa_in)));
    if yfracN>0
     thmit=[thmit yfracN];
    end
        thmit=[thmit tha];

CS0=cumsum([thmit]);



  if Npair>0
   no_pair=1;
   no_rep=0;
  else
   per0=sum(thmi(1:2));
   y010=sqrt((R+ha).^2-Rel^2)-R;
   y020=sqrt((R+ha).^2-Rm_ring^2)-R;
   yrest=ha-y020;
%   yrest=ha-y010;

    y01tip0=sqrt((R+ha).^2-Rel^2)-R+Hrei;
    y02tip0=sqrt((R+ha).^2-Rm_ring^2)-R+Hrei;
%   Npre=floor((Hrei-yrest)/per0)
   Npre=floor((y02tip0-ha)/per0);
%   Npre=floor((Hrei)/per0)
%   y020u=sqrt((R+ha+Hrei).^2-Rel^2)-R;
%   yrest=ha-y020;
%   Npre=floor((Hrei-yrest)/per0)
%   keyboard
%   if Npre<2
   if Npre<1
    Npre=0;
    Had1=Hrei;
   else
    Had1=Hrei-(Npre-1)*per0-(ha-y010);
    Hre=Hre-(Npre-1)*per0;
    if Hre<per0
     Hre=Hre+per0;
    end
   end
   Npev=floor(2*ha/per0)+1;
   Had2=2*ha;
%   ' yzo Npre', keyboard
    if Npre==0
     Npair=round(Had1/per0)+ round(Had2/per0/2);
%     Npair=ceil(Had1/per0)+ round(Had2/per0/2)+1;
    else
     Npair=round((Had1+Had2/2)/per0);
    end
    if Npair<1
     Npair=1;
    end
    if Npair>abs(Npa_in)
     Npair=abs(Npa_in);
    end    
  Npair=Npair+fracN;

  if Npair+Npre<=abs(Npa_in) & Hrei>0
   Npair=Npair+1;
  end


    y010=sqrt((R+ha).^2-Rel^2)-R;
%    'tip', keyboard
    y01tip=sqrt((R+ha).^2-Rel^2)-R+Hre;
    y02tip=sqrt((R+ha).^2-Rm_ring^2)-R+Hre;
    if Hrei==0
     y02tip=y01tip;
    end

    CS=cumsum([0 thmit])+ha;
%    fii=find(CS-y02tip0>0);
    fii=find(CS-y01tip0>=-1e-10);
    fii0=find(CS-y01tip0>0);
    fiia=find(CS0-ha>0);
    fiia=fii0;

%    Nsub=ceil((fiia(1)-1)/2)
    if ive==1
     Nsub=0;
    end

    NR=floor((fii(1)-2)/2);
%    NR=round((fii(1)-1)/2-.01);
    if NR<=0
     NR=1;
    end

    if NR<=1
     haz=CS(fii0(1))-(NR-1)*per0;
    else
     haz=ha;
    end
 %   'haz', keyboard
     if ive==1
      haz=ha;
     end


%    if NR>1
%     Nsub=Nsu0+floor((fiia(1)-1)/2)-NR-1
%    else
%     Nsub=Nsu0
%    end
%%    Nsub=ceil(y01tip/per0)-1
%    Nsub=0;

    y01=CS(fii(1))-(NR-1)*per0;
%    if NR<=2 & is_even(fii(1))==1
    if NR<=2
     y01=CS(fii(1));
    end
    ipai=0;
    if CS(fii(1))-y01tip0>0 & NR>2
     y01=CS(fii(1))-(NR-2)*per0;
     ipai=1;
    end
      if y01<ha
       y01=ha;
      end
%      'y01', keyboard
%      'y01', pausak
    if NR<floor(abs(Npa_in))
     Npair=ceil((y01-h)/per0)+fracN+1;
     Npair=floor((y01-h)/per0)+fracN+1;
     Npair=round((y01)/per0)+fracN+1;
    else
     Npair=ceil((y01-h)/per0)+fracN;
    end
if Npair<=1
 Npair=2;
end
    Npair0=Npair;
%    Npair=Npair+Nsub
%  Npair=4.5
%    pausak

    yrei=y01-ha;
    yreu=abs(Npa_in)*per0-yrei;
    Npu=floor(yreu/per0)-1;
    ym0=h+floor(Npair)*per0+yfracN;
    ym=h+floor(Npair)*per0+yfracN+suh/1000;
    if Npre==0
     if Hrei>0
      Hs=Hrei;
     else
      Hs=ha;
     end
     dtot_mir=(floor(abs(Npa_in))*per0+yfracN);
     dtot=(floor(abs(Npa_in))*per0+yfracN)+ha;
     NM=round((dtot-ym0)/per0);
     Npair=abs(Npa_in)-NM;
     
%     'cap', keyboard
     (y01-h+yfracN)/per0;
     Repet=[1 NM+1 1];
     if y01>haz
      haz=y01;
     end
     if ive==1
      haz=y01;
     end
     yzo=[haz haz+per0 ym0 ym];
     yzo0=[haz haz+per0 ym0];

%     'yzo', keyboard
%     'yzo', pausak
%     pausak

      Ndif=NM+Npair-abs(Npa_in);
%      'NM', keyboard
      if Ndif>0 & Hrei>0 & Npair>NM
       NM=NM-Ndif;
       Npair=abs(Npa_in);
       Hre=Hrei;
       ym=h+floor(Npair)*per0+yfracN;
       Repet=[1];
       yzo=[ym];
      end
%     ' fine Npre=0', keyboard

    else
     NR=Npre;
%     y01=CS(fii(1))-(NR-1)*per0;
%     NM=floor((abs(Npa_in)*per0-(Hrei+y010-ha)-ha)/per0);
     NM=round((floor(abs(Npa_in))*per0+yfracN-(Hrei+y010-ha)-ha)/per0);
     NM=floor((floor(abs(Npa_in))*per0+yfracN-Hrei-ha)/per0);

%     'NM2', keyboard
%     'NM2', pausak


     dtot=(floor(abs(Npa_in))*per0+yfracN)+ha;
     if NM>0
      NM=round((dtot-ha-Hrei-(ym-y01-per0))/per0);
     else
      NM=round((dtot-ha-Hrei-(ym-y01))/per0);
     end
%     Npair=abs(Npa_in)-NM;
     if NM+NR-2+Npair>abs(Npa_in)
      NM=NM-1;
     end
     if NM<0
      NM=0
     end
     if NM==0
      Npair=abs(Npa_in)-NR+1;
     end
%     Npair=4.5
%     'NM1', keyboard
%     'NM1', pausak



     if NM>=1
      NM=NM;
      Repet=[1 NR 1 NM 1];

      yzo=[haz haz+per0 y01 y01+per0 ym];
     else
      Repet=[1 NR 1];
      Ndif=NR+Npair-abs(Npa_in)-2;
%      Npair=Npair-1;
      ym=h+floor(Npair)*per0+yfracN;
      yzo=[haz haz+per0 ym];
%      'NR', keyboard
      if Ndif>0
%       NR=NR-Ndif;
%'qui', keyboard
       if NR>1
%        Npair=Npair+Ndif;
        Hre=Hrei-(NR-1)*per0;
        Repet=[1 NR 1];
        ym=h+floor(Npair)*per0+yfracN;
        yzo=[haz haz+per0 ym];
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
%       'yzo 2', keyboard
    end

  end
%  yzo=yzo-1e-3;
%  Npair=Npair+fracN;
%      ym=h+floor(Npair)*per0+yfracN;
%  ' qui ym', keyboard
if Npair>abs(Npa_in)
 Npair=abs(Npa_in);
 Repet=ones(size(Repet));
 ym=h+floor(Npair)*per0+yfracN;
 yzo=[haz haz+per0 ym];
end 

if Npair>0
 nps=2;
 nps=length(thmi);
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
 th=[duth tha];
 dunr=[nl(1:nps+1) nl(end)];
 dusa0=dunr;
 dusa=dunr(2:end-1);
%' problema', keyboard
 clear dunr
 dunr=repmat(dusa,1,fix(Npair));
 
 if Npair-floor(Npair)>0
  dunr=[dunr dunr(1)];
 end
 if length(thin)>nps
%  dunr=[dunr nin(nps+2)];
  dunr=[dunr nin(nps+1)];
 end
 dunr=[dusa0(1) dunr dusa0(end)];
 if length(nia)>0 
  dunr=[dunr(1:end-1) nia];
 end
 nl=dunr;

end
%Hre=y01tip-Hrei;
%Hre=Hrei;
%' st keyboard ', keyboard
if N<0
 stdu=abs(N)/1000;
% N=ceil(per0/stdu);
 N1=ceil(thmi(1)/stdu);
 N2=ceil(thmi(2)/stdu);
 N=N1+N2;
 end
 st=per0/(N+1);

y=0:st:ym;
thm=thmi;
%st1=thm(1)/(fix(N/2)+1);
%st2=thm(2)/(fix(N/2)+1);
st1=thm(1)/N1;
st2=thm(2)/N2;
y01d=[st1:st1:thm(1)];
y02d=[(st2:st2:thm(2))];
y0=[y01d y01d(end)+y02d];
%Npr=fix(mean([length(y01d) length(y02d)]));
%Npr=fix(2*h/per0*N);
Npr=fix(h/per0*N);
st0=h/(Npr+1);
%st0=h/(N+1);
%st0=st;
y=[0:st0:h];
for Np=1:fix(Npair)
 y=[y y(end)+y0];
end
if yfracN>0
 y=[y y(end)+y01d];
end

%'Nin', keyboard
if length(tha)>0 
 if Nin<0
  stin=abs(Nin)/1000;
  Npst=ceil((ym-ym0)/stin);
  stin=(ym-ym0)/(Npst+1);
 else
  stin=(ym-ym0)/(Nin+1);
 end
 y0add=ym0+stin:stin:ym;
else 
 y0add=[];
end
yold=y;
y=[y y0add];

%'y0add', keyboard
%if y(end)<ym
% y=[y ym];
%end

yi=cumsum([h*1000 duth tha])/1000;
hrel=h+Hre;

%yi=sort([yi hrel]);
duxr0=[Rel Rm_ring];
%' qui ic start', keyboard
ic=0;
for k=1:length(yold)
 duxr1=[];
 fn01=0;
 fi01=[];
 dux=sqrt((R+yi).^2-(R+y(k)).^2);
% dux=sqrt((yi).^2-(R+y(k)).^2);
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
    if y(k)>y02tip & y(k)<y01tip
     fi=[];
%     'uno'
     dur=rdux(fi01);
     fir=find(dur==0);
     fi01=fi01(fir);
     duxr=sqrt((R+Hre).^2-(R+y(k)-h).^2);
     fiv=find(rdux>duxr & rdux<rmax);
     duxr1=rdux(fiv);
     dux2=sort([dux1 duxr]);
    else
%     ' passo'
     fiv=find(rdux<=rmax & rdux>0);
%     fiv=find(rdux>duxr & rdux<rmax);
     duxr1=rdux(fiv);
     dux2=sort([dux1]);
%     fi01=[];
     dur=rdux(fi01);
     fir=find(dur==0);
     fi01=fi01(fir);
%     fi01=fiv;
%     k, pausak
    end
   elseif length(dux1)==1
    if y(k)<h
%    ' qui ', keyboard
     fi01=find(rdux<rmax & idux==0 & dux>dux1(1));
     if y(k)>y02tip & y(k)<y01tip
      duxr=sqrt((R+Hre).^2-(R+y(k)-h).^2);
      duxre=sqrt((R).^2-(R+y(k)-h).^2);
      fiti=find(duxre<duxr0(1));
      duxi=duxre(fiti);
      fiv=find(rdux>duxr & rdux<rmax);
      duxr1=rdux(fiv);
      dux2=sort([duxi duxr duxr0(1) ]);
%      fi01=[];
      dur=rdux(fi01);
      fir=find(dur==0);
      fi01=fi01(fir);
     end
    else
%     fi01=[];
     dur=rdux(fi01);
     fir=find(dur==0);
     fi01=fi01(fir);
     duxr=sqrt((R+Hre).^2-(R+y(k)-h).^2);
     fiv=find(rdux>duxr & rdux<rmax);
     duxr1=rdux(fiv);
     dux2=sort([dux1 duxr]);
     fi01=find(rdux<rmax & idux==0 & (dux<dux2(1) | dux>dux2(2)));
%     fi01=[];
    end
   else
    fi01=find(rdux<rmax & idux==0 & (dux<dux1(1) | dux>dux1(2)));
   end
   if ifi0==0
    [xad,isov]=sort([dux(fi01) dux2 duxr1]);
   else
    [xad,isov]=sort([dux2 duxr1]);
   end
   if length(find(diff(xad)==0))>0
     fiv=find(diff(xad)>0);
     fiv0=find(diff(xad)==0);
     fivu=sort([fiv fiv0]);
     xad=xad(fivu);
%    'xaa', keyboard
   end
  if ic==k0 & ife==1
   'k ferma xi',
   ic-1
   keyboard
  end

 fi01sa=fi01;
 if ifi0==0
  fi01=fi00;
 end
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
      else
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
    if y(k)>ym-ha
     ruad=[];
    else
     ruad=r;
    end
   end
   if iad0>1
    ic=ic+1;
    yu(ic)=y(k-1);
    xu(ic,1)=r ;
    if ife==1
    ' qui ', keyboard
    end

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
%   'prima volta', keyboard
  if ic-1==k0 & ife==1
   'k ferma ni',
   rnad=real(nad);
   ic-1
   keyboard
  end
   if length(nad)<=length(xad)
    nad=[nad(1:2) nad];
   end
   nu(ic,1:length(nad))=nad;
   end
   iad0=0;
 else

   iad0=iad0+1;
%   ' qui ', keyboard
 end
end
yus=yu;

fixa=find(xu>0);
amag=sort(xu(fixa));
damag=diff(amag);
NPra_vecchio=length(find(damag>0));
amin=min(amag);
st=.04;
st=.02;
st=.05;
%'st', keyboard
radig=sort([r:-st:amin-st]);
xsa=xu;

for kdi=1:length(radig)
 ai=radig(kdi)-st/2;
 af=radig(kdi)+st/2;
 fi=find(xu>=ai & xu<af);
 xu(fi)=radig(kdi);
end
fixa=find(xu>0);
amag=sort(xu(fixa));
damag=diff(amag);
NPra=length(find(damag>0));
[NPra NPra_vecchio]
sx=size(xu,2);
yusa=repmat(yu',1,sx);
dx=xu-xsa;
fi=find(abs(dx)>st);

lxu=length(nu);
nsav=nu;
for lpu=1:lxu
 fima=find(xu(lpu,:)==r);
 if length(fima)==1
  nu(lpu,fima+1)=nu(lpu,fima);
 end
end 

if ifp<=-10
figure, plot(xsa,yu,'.')
hold on, plot(xu(fi),yu(fi),'o')
%xu=xsa;
' qui dig fineeeeeeee', pausak
end



%yu=yu-yu(1)+diff(yu(1:2));

%xusa=xu;
% xun=xu;
%siz=size(xu);
%xud=reshape(xu,prod(siz),1);
%xud=sort(xud);
%firma=find(diff(xud)>.1);
%firmi=find(diff(xud)<=.1);
%for kk=1:length(firma)
% xui=xud(firma(kk));
% fii=find(abs(xu-xui)<=.1);
% xun(fii)=xui;
%end
%xu=xun;

%' xu ', keyboard



if ifp<=-11
  figure, plot(xu,'.'), grid
% ' in sub_rel '
% pausak

% figure,
% plot(xu,yu,'.'), grid
% for k=1:length(yzo)
%  hold on, plot(xd,yd0*yzo(k),'w')
% end
% pausak

psches

 for k=1:length(yzo)
  hold on, plot(xd,yd0*yzo(k),'w')
 end
 pausak
end

%yud=[yu];
%dilu=diff(yud)';
%dilu=[dilu; dilu(end)];
yud=[0 yu];
yud=[yu 0];
dilu=diff(yud)';
dilu(end)=dilu(end-1);
fstu=zeros(length(dilu),2);
ailu=(xu);
nilu=(nu);

yinf=0;
shi=0;
for k=1:length(yzo0)
 ysup=yzo(k);
 fi=find(yu>=yinf & yu<ysup);
% fi=find(yu>yinf & yu<=ysup);
 if k==length(yzo0)
  fi=find(yu>=yinf & yu<=ysup);
 end

 fstu(fi,2)=Repet(k);
 fstu(fi,1)=length(fi);
% [k yinf ysup], pausak
% cs=cumsum(dilu(fi));
% figure, plot(shi+cs,fstu(fi,2),cs+shi,nilu(fi,1),'r.'), pausak
% shi=shi+sum(dilu(fi));
 yinf=ysup;
end




%' quo fstu', keyboard
%cs=cumsum(dilu);
%figure, plot(cumsum(dilu),fstu(:,2),cs,nilu(:,1),'r.')
%keyboard


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


if ifp<=-11
 figure
 plot(xd,yd,'.'), grid, pausak
 figure
 plot(yd,nd,'.',y0,n0), grid
 pausak
% figure, plot(fstu(:,2),'.')
% pausak

end

isaq=1;

if isaq==0
dilu=[dilu(1); dilu];
nilu=[nilu(1,1)*ones(size(nilu(1,:))); nilu];
ailu=[ailu(1,:)*0; ailu];
fstu=[[0 1]; fstu];
ailu(end,:)=ailu(end,:)/2;

ailu=flipud(ailu);
nilu=flipud(nilu);
dilu=flipud(dilu)*1000;
fstu=flipud(fstu);
end  %isaq

sa=size(ailu);
if exist('Rflat')
for ks=1:sa(2)-1
 fi=find(ailu(:,ks)~=0);
 ailu(fi,ks)=ailu(fi,ks)+Rflat;
end
end
for ks=2:sa(2)
 fi=find(ailu(:,ks)==0);
 ailu(fi,ks)=r;
end

inew=1;
if inew==1
if ifp<=-11
%' modifica'
%keyboard
end
ailu_s=ailu;
nilu_s=nilu;
%fi=find(ailu==r);
%ailu(fi)=0;
%%ailu(:,:)=0;
%' cas', keyboard
clear ailu
clear nilu
for k=1:length(ailu_s)
 fi=find(ailu_s(k,:)==r);
 if length(fi)>0
  if fi>1
   ailu(k,1:fi-1)=ailu_s(k,1:fi-1);
   nilu(k,1:fi)=nilu_s(k,1:fi);
  else
   ailu(k,1)=0;
   nilu(k,1)=nilu_s(k,1);
  end
 else
  la=length(ailu_s(k,:));
  lan=length(nilu_s(k,:));
  ailu(k,1:la)=ailu_s(k,:);
%  nilu(k,1:la+1)=nilu_s(k,:);
  nilu(k,1:lan)=nilu_s(k,:);
 end
end
sa=size(ailu);
k=sa(2);
 fi=find(ailu(:,k)==0);
 if length(fi)==sa(1);
  ailu=ailu(:,1:end-1);
  nilu=nilu(:,1:end-1);
 end

if ifp<=-11
figure, plot(ailu,cumsum(dilu),'.')
a=axis;
a(1)=0;
axis(a);
' modifica'
keyboard
end
end
