function [gazk,enan,fak,Tue,Tb,Gue,Gum,Gbe,Gbm,Lf,Lcav,ztot,Ez,Hz,indz,nmean,Perd,Ge,Gm,lQW,Ezm,KKie,KKim]=...
         th_scattu3last(fiQWi,fiCavi,L_i,n_i,rr,nb,nu,lambda0,fr,kv,iLP,ifp,iff,ibasti,par_grat,NPXi)

ibast=abs(ibasti);
 if ifp==-10
  isto=1;
 else
  isto=0;
 end
ifer0=0;
%'debbtro orta', keyboard
global isempl icarico

if length(isempl)==0
isempl=1;
end
%isempl=1;
%isempl=0;

Ge=1;         
Gm=1;   
KKieB=1;         
KKimB=1;         
Gm=1;   
%ifp=2
%iff=1
%Lv=flipud(L_i);
%rv=flipud(n_i);
% prova plot
if iff==2
 figure, plot(cumsum(L_i(fiCavi)),n_i(fiCavi),'.'), pausak
 keyboard
end

ru=nu;
rb=nb;
Lu=[];
Lb=[];
%ru=[];
%rb=[];
fiQW=fiQWi+length(Lu);
fiCav=fiCavi+length(Lu);
Lv=[L_i/1000];
Lt=sum(Lv);
rv=[n_i];
z=[0; (cumsum(Lv))];
NQW=length(fiQW);
iqw=round(NQW/2);
fiQW_cen=fiQW(iqw);
dQW=L_i(fiQW_cen);
%if iff==1
nstu=[1:fiQW(iqw)];
%else
%nstu=[1:fiQW(iqw)-1];
%end
nstb=[1+fiQW(iqw):length(Lv)];

%'qui strati', keyboard

isap=0;
if isap==1
 nstu=[1:3];
 nstb=[1+3:length(Lv)];
end

%
Lvb=flipud(Lv(nstb));
rvb=flipud(rv(nstb));
Lvu=(Lv(nstu));
rvu=(rv(nstu));

ibastu=ibast;
ibastb=-100;
if ibast>0
 if ibast>length(Lvu) 
  ibastbd=ibast-length(Lvu);
  if iff==1
   ibastb=-ibastbd+length(Lvb)+1;
  else 
   ibastb=-ibastbd+length(Lvb)+1;
  end 

  ibastu=-100;
 end
end

  if iff==1
%     'qui strati scattu', keyboard
  end
%'qui strati', keyboard
if length(fiQW)>1
 fiQWd=[1:2:length(fiQW)*2-1];
 fiQWu=length(Lvu)+1-fiQWd(1:iqw);
else
 fiQWd=[];
 fiQWu=length(Lvu);
end

%rqw=rvu(fiQWu(1));

%'rqw', keyboard
rqw=rv(fiQW(iqw));

subfiqw=find(fiCav==fiQW(iqw));
ficavu=fiCav(1:subfiqw-1);

fiQWb=length(Lvb)-fiQWd(1:iqw-1);
subdu=length(fiCav)-subfiqw-1;
ficavb=(length(Lvb)-subdu):length(Lvb);



%if length(fiQWu)>length(fiQWd)
% fiQWb=length(Lvb)-fiQWd(1:iqw-1);
% subdu=length(fiCav)-subfiqw-1;
% ficavb=(length(Lvb)-subdu):length(Lvb);
%else
% ficavb=[];
% fiQWb=[];
%end
%' gaz7 ', keyboard

fiQW0=nstu(length(nstu));
%
%disp('pro_scatt entro'), keyboard
%
zu=[0; (cumsum(Lvu))];
zb=[0; (cumsum(Lvb))];
mu=rr/nu;
mub=rr/nb;

j=sqrt(-1.);
flp=1-iLP;
k0=2*pi/lambda0;


sginu=-1;

f1u=[1; sginu/mu];
f1um=f1u;
f1ue=f1u;
f1d=[1; sginu/mub];
f1dm=f1d;
f1de=f1d;
% coefficienti alle interfacce per riflessione

ic=1;
zi=0;
nz(1)=rv(1);

if exist('NPXi')
else
 NPXi=15;
end

F1u=[1; sginu/mu];
F1ue=F1u;
F1um=F1u;
zio=0;
nzo=[];
F1d=[1; sginu/mub];
F1de=F1d;
F1dm=F1d;
% campi in z nella struttura

zid=0;

nzd=[];

 nk=1;
 k=kv;
 Itue=0;
 Itum=0;
 Iau=0;
 Icau=0;
 Icauh=0;
 ic=1;

cont_grat=0;
ibastus=ibastu;
KKe=1;
KKm=1;
iferma=0; 
 for nst=1:length(Lvu)
  ic=ic+1;
  r=rvu(nst);
  dx=Lvu(nst); 
% if nst==ibastu
if length(find(ibastu-nst==0))==1
 if sum(ibastus)==sum(ibastu)
  KKe=1;
  KKm=1;
 end
 iferma=ifer0;
 if iferma==1
%  ' qui set ret ', keyboard
 end 
  if ibastus~=0
    cont_grat=cont_grat+1;
    ibadif=[ibastus(1)-1 ibastus];
    fico=find(diff(ibadif)~=1); 
    if length(fico)==0
     Nstrati=length(ibastus);
     nstu=ibastus(Nstrati)+1;
     ibastus=0;
    else 
     Nstrati=length(ibastus(1:fico-1));
     nstu=ibastus(Nstrati)+1;
     ibastus=ibastus(fico:end);
    end 
    
    nsti=nst-1;
    r_out=rvu(nstu);
    r_in=rvu(nsti);
    Mu=[1 1; r_out/rr -r_out/rr]/sqrt(r_out); 
    Mi=0.5*sqrt(r_in)*[1 rr/r_in; 1 -rr/r_in];     
%    Mi_Va=[1 1; rin -rin]/sqrt(rin); 
%    Mu_aV=0.5*sqrt(rout)*[1 1/rout; 1 -1/rout]; 
    puretloc=nsti+1:nstu-1;
    r=rvu(puretloc);
    n1=r;

    fimen=find(ibasti<0);
    if length(fimen)>0
     pun=1:fimen-1;
     ilay=puretloc(pun);
     n2=rvu(ilay+1)*ones(size(n1));      
     n2(pun)=r_in;
    else 
     n2=r_in*ones(size(n1));      
    end  
    

    dx=Lvu(puretloc); 


%  rout=rvu(nst-1);
%  rin=rvu(nst+1);
%   n2=rout;
%  n1=r;

  tetai=0.01;
  itetm=par_grat.itetm;
  px=par_grat.px;
  DC=par_grat.DC;
  NModi=par_grat.NModi;
  d1i=px*DC;
  d2i=px*(1-DC);


%' qui primo strato', nst
%' qui orta ', keyboard
Strut.t=L_i;
Strut.n=n_i;
Strut.ibast=ibast;
Strut.icarico=icarico;
 
 [Te,Tm,Ge,Gm]=Orta_tracar(tetai,r_in,r_out,n1,n2,d1i,d2i,dx,lambda0,NModi,itetm,Strut);

    KKie=(Mu*Te*Mi);
    KKim=(Mu*Tm*Mi);
if iferma==1
% 'dopo ret orta ', keyboard
end 

    iverif=0;
    if iverif==1
     k=0;
     r=nef(1);
     mi=rr/r;
     sq=sqrt(1-mi^2.*k^2)';
     sqz=sqrt(1-flp*mi^2.*k^2)';
     fi=1;
     be=k0*(1+fr)*r*sq*fi;
     dte=be*dx;
     KKi=[cos(dte),              j*sin(dte).*mi/sqz;...
       j*sin(dte)./mi*sqz,     cos(dte)];    
     'verifica KKi'
     KKi
     pausak
    end  %verif

   
   if iff==1
    dx=sum(dx);
    rm=mean([n1; n2]);
%   'dopo reticolo', keyboard
    if isempl==1
     [F1ue,zio,nzo]=fieldxt(dx,rm,KKie,F1ue,zio,nzo);
     [F1um]=fieldxt(dx,rm,KKim,F1um);
%     [F1de,zid,nzd]=fieldxt(dx,rm,KKie,F1de,zid,nzd);
%     [F1dm]=fieldxt(dx,rm,KKim,F1dm);
    else
     [dF1ue]=fieldxt(dx,rm,KKie,F1ue);
     [dF1um]=fieldxt(dx,rm,KKim,F1um);
     [F1ue,F1um,zio,nzo]=field_ret(NPXi,F1ue,F1um,zio,nzo,r_in,r_out,n1,n2,d1i,d2i,dx,lambda0,NModi,rr);
      F1ue(:,end)=dF1ue(:,end);
      F1um(:,end)=dF1um(:,end);
    end
   end  %iff
%   ' fine reticolo', keyboard
 else  %ibastus
   KKie=[1 0; 0 1];
   KKim=[1 0; 0 1];
%' salto strato', nst, pausak

%' qui orta salto ', keyboard

 end
 %  KKie=(Mi*Te*Mu);
 %  KKim=(Mi*Tm*Mu);
 % 'cont bast con conj per verifica sperimentale  ver_bast da capire (forse direzione inversa) ', keyboard 

else  %non reticolo
  mi=rr/r;
  sq=sqrt(1-mi^2.*k^2)';
  sqz=sqrt(1-flp*mi^2.*k^2)';
 % fi=r/real(r);
  fi=1;
  sqfi=sq*fi;
  if imag(fi)<0
   'fi ver',  pausak
  end
  if iff==1
   NPXin=NPXi;
%   if nst==fiQWu
   if length(find(nst==fiQWu))==1

    NPXin=-NPXi;
%' meno', keyboard
   end
%function [fz,hz,z,nz]=fieldz(L_i,n_i,NPXi,rr,nb,nu,lambda0,fr,k,iLP)
   [F1ue,zio,nzo]=fieldx(dx,r,NPXin,rr,sqfi,sqz,lambda0,fr,F1ue,zio,nzo);
   [F1um]=fieldx(dx,r,NPXin,rr,sqfi,sqz,lambda0,fr,F1um);
    if nst==length(Lvu)-1
     lpqw=length(zio);
    end 
  end
  be=k0*(1+fr)*r*sq*fi;
  dte=be*dx;
  KKi=[cos(dte),              -j*sin(dte).*mi/sqz;...
       -j*sin(dte)./mi*sqz,     cos(dte)];
%  if dx==100
%   [dx r]
%'dentro loop ',   keyboard
%  end

   nz(ic)=r;
%   f0(:,ic)=KKi*f0(:,ic-1);
   s2=sin(2*dte);
   Is=.5*[(dte+s2/2) (mi/sqz).^2*(dte-s2/2)];
   Inte=abs((Is*abs(f1ue(:,ic-1)).^2/mi^2)/be);
   Intm=abs((Is*abs(f1um(:,ic-1)).^2/mi^2)/be);
   Inteu(nst,1)=Inte;
   Intmu(nst,1)=Intm;
   Inteup(nst)=Inte*abs(mi^2)*imag(1/mi^2);
   Intmup(nst)=Intm*abs(mi^2)*imag(1/mi^2);
   Itue=Itue+Inte;
   Itum=Itum+Intm;
   Is=.5*[ (sqz/mi).^2*(dte-s2/2) (dte+s2/2)];
   Inthe=abs((Is*abs(f1ue(:,ic-1)).^2)/be);
   Inteu(nst,2)=Inthe;
   Inthm=abs((Is*abs(f1um(:,ic-1)).^2)/be);
   Intmu(nst,2)=Inthm;
   KKie=KKi;
   KKim=KKi;
  end 
   f1ue(:,ic)=KKie*f1ue(:,ic-1);
   f1um(:,ic)=KKim*f1um(:,ic-1);
   if iferma==1
     nst
    ' FERMA prima di KKm', pausak
   end 
   KKm=KKim*KKm;
   KKe=KKie*KKe;
   if iferma==1
     nst
    ' FERMA', pausak
     if nst==ibast(end)
     KKe
     KKm
       ' controlla K', pausak
     end
   end
 end  %L_v
 icu=ic;
 icumeno=ic-1;
 miu=mi;

 if iff==1
  lul=length(zio);
  lQW=[lpqw:lul];
% ' qui zio', keyboard
 else
   lQW=1;
 end

 
 miqw=rr/rqw;
 A0ufe=(f1ue(1,1)+f1ue(2,1)*mu)/sqrt(mu)/2;      % a_u
 A0ube=(f1ue(1,1)-f1ue(2,1)*mu)/sqrt(mu)/2;      % b_u
 A0ufh=(f1um(1,1)+f1um(2,1)*mu)/sqrt(mu)/2;      % a_u
 A0ubh=(f1um(1,1)-f1um(2,1)*mu)/sqrt(mu)/2;      % b_u

 Aufe=(f1ue(1,icumeno)+f1ue(2,icumeno)*miqw)/sqrt(miqw)/2; % a_i
 Aube=(f1ue(1,icumeno)-f1ue(2,icumeno)*miqw)/sqrt(miqw)/2; % b_i
 Aufh=(f1um(1,icumeno)+f1um(2,icumeno)*miqw)/sqrt(miqw)/2; % a_i
 Aubh=(f1um(1,icumeno)-f1um(2,icumeno)*miqw)/sqrt(miqw)/2; % b_i

 Tue(nk,1)=A0ufe/Aufe;
 Tum(nk,1)=A0ufh/Aufh;
 Gue(nk,1)=Aufe/Aube;
 Gum(nk,1)=Aufh/Aubh;


%' fine su ', keyboard

ibastbs=ibastb;

ic=1;
for nst=1:length(Lvb)
  ic=ic+1;
  r=rvb(nst);
  dx=Lvb(nst); 
 if length(find(ibastb-nst==0))==1
  if ibastbs~=0
    ibadif=[ibastbs(1)-1 ibastbs];
    fico=find(diff(ibadif)~=1); 
    if length(fico)==0
     Nstrati=length(ibastbs);
     ibastbs=0;
    else 
     Nstrati=length(ibastbs(1:fico));
     ibastbs=length(ibastbs(fico+1:end));
    end 
    nsti=nst-1;
    nstu=ibastb(Nstrati)+1;
    r_out=rvb(nstu);
    r_in=rvb(nsti);
    Mu=[1 1; r_out/rr -r_out/rr]/sqrt(r_out); 
    Mi=0.5*sqrt(r_in)*[1 rr/r_in; 1 -rr/r_in];     
    puretloc=nsti+1:nstu-1;
   
    r=rvb(puretloc);
    dx=Lvb(puretloc); 
    n1=r;
    fimen=find(ibasti<0);
    if length(fimen)>0
     pun=1:fimen-1;
     ilay=puretloc(pun);
     n2=rvu(ilay+1)*ones(size(n1));      
     n2(pun)=r_in;
    else 
     n2=r_in*ones(size(n1));      
    end  

    tetai=0.01;
    itetm=par_grat.itetm;
    px=par_grat.px;
    DC=par_grat.DC;
    NModi=par_grat.NModi;
    d1i=px*DC;
    d2i=px*(1-DC);

' qui orta sotto scattu3 ', keyboard
%' qui primo strato', nst
Strut.t=L_i;
Strut.n=n_i;
Strut.ibast=ibast;
Strut.icarico=icarico;
 
 [Te,Tm,Ge,Gm]=Orta_tracar(tetai,r_in,r_out,n1,n2,d1i,d2i,dx,lambda0,NModi,itetm,Strut);


%  [Te,Tm,Ge,Gm]=Orta_tras(tetai,r_in,r_out,n1,n2,d1i,d2i,dx,lambda0,NModi,itetm);
    KKie=(Mu*Te*Mi);
    KKim=(Mu*Tm*Mi);

   if iff==1
    dx=sum(dx);
    rm=mean([n1; n2]);
   'dopo reticolo', keyboard
    if isempl==1
     [F1de,zid,nzd]=fieldxt(dx,rm,KKie,F1de,zid,nzd);
     [F1dm]=fieldxt(dx,rm,KKim,F1dm);
    else 
     [dF1de]=fieldxt(dx,rm,KKie,F1de);
     [dF1dm]=fieldxt(dx,rm,KKim,F1dm);
     [F1de,F1dm,zid,nzd]=field_ret(NPXi,F1de,F1dm,zid,nzd,r_in,r_out,n1,n2,d1i,d2i,dx,lambda0,NModi,rr);
      F1de(:,end)=dF1de(:,end);
      F1dm(:,end)=dF1dm(:,end);     
    end
   end  %iff  
%   ' fine reticolo', keyboard
  else  %ibastbs
   KKie=[1 0; 0 1];
   KKim=[1 0; 0 1];
%' salto strato', nst, pausak

%' qui orta salto ', keyboard

 end

%' fine Orta post dopo reticolo scattu3 ', keyboard

%   ' fine Bast in th_scattu ', keyboard
 %  KKie=(Mi*Te*Mu);
 %  KKim=(Mi*Tm*Mu);
 % 'cont bast con conj per verifica sperimentale  ver_bast da capire (forse direzione inversa) ', keyboard 
 else
  mi=rr/r;
  sq=sqrt(1-mi^2.*k^2)';
  sqz=sqrt(1-flp*mi^2.*k^2)';
%  fi=r/real(r);
  fi=1;
   if iff==1
   NPXin=NPXi;
   if length(find(nst==fiQWb))==1
    NPXin=-NPXi;
%' meno', keyboard
   end   
%    [f1d,zid,nzd]=fieldx(dx,r,NPXi,rr,sqfi,sqz,lambda0,fr,f1d,zid,nzd);
    [F1de,zid,nzd]=fieldx(dx,r,NPXin,rr,sqfi,sqz,lambda0,fr,F1de,zid,nzd);
    [F1dm]=fieldx(dx,r,NPXin,rr,sqfi,sqz,lambda0,fr,F1dm);
%[length(F1de) length(F1dm)]
%    ' F1de F1dm', pausak
   end
   be=k0*(1+fr)*r*sq*fi;
   dte=be*dx;
   KKi=[cos(dte),              -j*sin(dte).*mi/sqz;...
       -j*sin(dte)./mi*sqz,     cos(dte)];
   nz(ic)=r;
%   f0(:,ic)=KKi*f0(:,ic-1);
   s2=sin(2*dte);
   Is=.5*[(dte+s2/2) (mi/sqz).^2*(dte-s2/2)];
   Int=abs((Is*abs(f1de(:,ic-1)).^2/mi^2)/be);
   Inteb(nst,1)=Int;
   Intebp(nst)=Int*abs(mi^2)*imag(1/mi^2);
   
   Int=abs((Is*abs(f1dm(:,ic-1)).^2/mi^2)/be);
   Intmb(nst,1)=Int;
   Intmbp(nst)=Int*abs(mi^2)*imag(1/mi^2);   
   
   Is=.5*[ (sqz/mi).^2*(dte-s2/2) (dte+s2/2)];
   Inth=abs((Is*abs(f1de(:,ic-1)).^2)/be);
   Inteb(nst,2)=Inth;
   Inth=abs((Is*abs(f1dm(:,ic-1)).^2)/be);
   Intmb(nst,2)=Inth;   
   
   KKie=KKi;
   KKim=KKi;
 end % ibast
%   f1ue(:,ic)=KKie*f1ue(:,ic-1);
%   f1um(:,ic)=KKim*f1um(:,ic-1);
%   KKp=KKie*KKp;
   f1de(:,ic)=KKie*f1de(:,ic-1);
   f1dm(:,ic)=KKim*f1dm(:,ic-1);
end  %L_v
 

if iff==1
 lfo=length(F1ue);
 lfd=length(F1de);
 fnce=F1ue(1,lfo)/F1de(1,lfd);
 fncm=F1um(1,lfo)/F1dm(1,lfd);
 Ez=[F1ue(1,:) fliplr(F1de(1,:))*fnce]/F1ue(1,lfo);
 Ezm=[F1um(1,:) fliplr(F1dm(1,:))*fncm]/F1um(1,lfo);
 Hz=[F1ue(2,:) fliplr(F1de(2,:))*fnce]/F1ue(1,lfo);
 Hzm=[F1um(2,:) fliplr(F1dm(2,:))*fncm]/F1um(1,lfo);

% lfo=length(f1oe);
% lfd=length(f1de);
% fnce=f1oe(1,lfo)/f1de(1,lfd);
% fncm=f1om(1,lfo)/f1dm(1,lfd);
% Ez=[f1oe(1,:) fliplr(f1de(1,:))*fnce]/f1oe(1,lfo);
% Ezm=[f1om(1,:) fliplr(f1dm(1,:))*fncm]/f1om(1,lfo);
% Hz=[f1oe(2,:) fliplr(f1de(2,:))*fnce]/f1oe(1,lfo);
% Hzm=[f1om(2,:) fliplr(f1dm(2,:))*fncm]/f1om(1,lfo);


 indz=[nzo fliplr(nzd)];
 zd=zid;
 pu1=[2:lfd 1];
 zd=zd(pu1);
 ztot=cumsum([zio; flipud(zd)]);
% ' dentro campi z', keyboard
%figure, plot(ztot,abs(Ez).^2*3,ztot,abs(Hz).^2*3,ztot,indz)
%keyboard
else
 ztot=[];
 indz=[];
 Ez=[];
 Ezm=[];
 Hz=[];
end

 icb=ic;
 mib=mi;
% fne=max(f1ue(1,:))/max(f1de(1,:));
 fne=(f1ue(1,icu))/(f1de(1,icb));
 fne2=abs(fne)^2;
 
% fnm=max(f1um(1,:))/max(f1dm(1,:));
 fnm=(f1um(1,icu))/(f1dm(1,icb));
 fnm2=abs(fnm)^2;

% pr=[Inteu' fliplr(Inteb')];
% figure, plot(abs(pr(1,:))), hold on, plot(Inteu,'.'),
%'qui scattu2', keyboard 

if isap==0

 Inteb=Inteb*fne2;
 Intebp=Intebp*fne2;

 Intmb=Intmb*fnm2;
 Intmbp=Intmbp*fnm2;

 Lca=sum(Lvu(ficavu))+sum(Lvb(ficavb));
 
 Icaee=sum(Inteu(ficavu,1))+sum(Inteb(ficavb,1));
 Icahe=sum(Inteu(ficavu,2))+sum(Inteb(ficavb,2));

 Icaem=sum(Intmu(ficavu,1))+sum(Intmb(ficavb,1));
 Icahm=sum(Intmu(ficavu,2))+sum(Intmb(ficavb,2));
 
 
if length(fiQWb)>0
 Iaee=sum(Inteu(fiQWu,1))+sum(Inteb(fiQWb,1));
 Iahe=sum(Inteu(fiQWu,2))+sum(Inteb(fiQWb,2));
 Iaem=sum(Intmu(fiQWu,1))+sum(Intmb(fiQWb,1));
 Iahm=sum(Intmu(fiQWu,2))+sum(Intmb(fiQWb,2));
else
 Iaee=sum(Inteu(fiQWu,1));
 Iahe=sum(Inteu(fiQWu,2));
 Iaem=sum(Intmu(fiQWu,1));
 Iahm=sum(Intmu(fiQWu,2));
end


 Icae=(Icaee+Icahe)/2;
 Icam=(Icaem+Icahm)/2;      %integrali su cavita
% Ica=Icae;

 Itee=sum(Inteu(:,1))+sum(Inteb(:,1));
 Itepe(1)=sum(Inteup);
 Itepe(2)=sum(Intebp);
 Ithe=sum(Inteu(:,2))+sum(Inteb(:,2));

 Item=sum(Intmu(:,1))+sum(Intmb(:,1));
 Itepm(1)=sum(Intmup);
 Itepm(2)=sum(Intmbp);
 Ithm=sum(Intmu(:,2))+sum(Intmb(:,2));


nmean=real((sum(Inteu(:,1).*rvu)+sum(Inteb(:,1).*rvb))/Itee);

%nmean=real((sum(Intzu(ficavu,1).*rvu(ficavu))+sum(Intzb(ficavb,1).*rvb(ficavb)))/Ica);
%' nmean ', keyboard

 ItE=(Itee+Ithe)/2;
 ItM=(Item+Ithm)/2;      % energia immagazzinata H totale
 
% It=Ite;
 Perd=Itepe/ItE;

% miqw2=(rr/rvu(fiQW0))^2;
% miqw2=1;
% miqw2=(rr/n_i(fiQW_cen))^2;
% keyboard

 IaE0=Inteu(fiQW0,1);
 IaM0=Intmu(fiQW0,1);
% Iau0=Intzu(fiQW0,1);
% Ia=(Iae+Iah)/2;
 IaE=Iaee;
 IaM=Iaem;
 gazk0(1,1)=IaE0/ItE;
 gazk0(2,1)=IaM0/ItM;

 gazk(1,1)=IaE/ItE;
 gazk(2,1)=IaM/ItM;
 fak(1,1)=IaE/(IaE0*NQW);

 fak(2,1)=IaM/(IaM0*NQW);
 Lf(1,1)=ItE/Icae*Lca; 
 Lf(2,1)=ItM/Icam*Lca; 
if iff==1
%'nuovo fak ', keyboard
end

else
gazk=1;
fak=1;
Iku=1;
Ikb=1;
Lf=1;
nmean=1;
Perd=1;
end
%' fak ', keyboard

% A0ufe=(f1ue(1,1)+f1ue(2,1)*mu)/sqrt(mu)/2;      % a_u
% A0ube=(f1ue(1,1)-f1ue(2,1)*mu)/sqrt(mu)/2;      % b_u
%
% Aufe=(f1ue(1,icu-1)+f1ue(2,icu-1)*miqw)/sqrt(miqw)/2; % a_i
% Aube=(f1ue(1,icu-1)-f1ue(2,icu-1)*miqw)/sqrt(miqw)/2; % b_i
%
% Tue(nk,1)=A0ufe/Aufe;
% Gue(nk,1)=Aube/Aufe;


 A0bfe=(f1de(1,1)+f1de(2,1)*mub)/sqrt(mub)/2;
 A0bbe=(f1de(1,1)-f1de(2,1)*mub)/sqrt(mub)/2;
 A0bfm=(f1dm(1,1)+f1dm(2,1)*mub)/sqrt(mub)/2;
 A0bbm=(f1dm(1,1)-f1dm(2,1)*mub)/sqrt(mub)/2;

 Abfe=(f1de(1,icb)+f1de(2,icb)*miqw)/sqrt(miqw)/2;
 Abbe=(f1de(1,icb)-f1de(2,icb)*miqw)/sqrt(miqw)/2;
 Abfm=(f1dm(1,icb)+f1dm(2,icb)*miqw)/sqrt(miqw)/2;
 Abbm=(f1dm(1,icb)-f1dm(2,icb)*miqw)/sqrt(miqw)/2;

 Tbe(nk,1)=A0bfe/Abfe;
 Tbm(nk,1)=A0bfm/Abfm;
 Gbe(nk,1)=Abfe/Abbe;
 Gbm(nk,1)=Abfm/Abbm;

 
 Tb(nk,1)=Tbe(nk,1);

%return
%Lcav=sum(L_i(fiCavi));
Lcav=1;

% IaE=Iaee*miqw2;
% IaM=Iaem*miqw2;
%
% gazk0(1,1)=IaE0/ItE;
% gazk0(2,1)=IaM0/ItM;
%
% gazk(1,1)=IaE/ItE;
% gazk(2,1)=IaM/ItM;
%
% fak(1,1)=IaE/(IaE0*NQW);
% fak(2,1)=IaM/(IaM0*NQW);
%
% Lf(1,1)=ItE/Icae*Lca; 
% Lf(2,1)=ItM/Icam*Lca; 


%fat=eat/eaca;
enan=[1 1];
return
fat=IaE0/Icae;
rd=1e-3*dQW/Lca;
enan(1)=fat/rd;
fat=IaM0/Icam;
enan(2)=fat/rd;


global istopk
if istopk==1
' fak nuovo', keyboard
end




