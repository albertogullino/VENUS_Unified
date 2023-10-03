vgconv=3e10/rr;
xri=xvero;
 xdx3=xri.^3*diff(xri(1:2));
 xdx1=xri.*diff(xri(1:2));
%iWidth=1;
kMO=0;  %setta quanti modi hanno iWidth
iNObre=0;

iUNO_long=0;  %trova tutti i modi longitudinali, 1 solo il dominante
%global uL Ppol
nAris=2;   % 1: uscita, 2: QW
%nAris=1;   % 1: uscita, 2: QW
%fPmax=fix(length(KK)/3);
%fPmax=fix(length(KK)/2);
%fPmax=fix(length(KK));
fPmax=10;
fPmax=0;
iexpe=1;      %=1 usa interp su alvet e gvet
inor=0;
if ~exist('ipost')
 ipost=0;
end 
   isav_Az=0;
   if length(Ps)>0
    if isfield(Ps,'isav_Az')==1
      isav_Az=Ps.isav_Az;
    end
   end
      isav_Post=0;
      if length(Ps)>0
       if isfield(Ps,'isav_Post')==1 & ipost==0
         isav_Post=Ps.isav_Post;
       end
      end
   if isav_Post==1
      nosamix=Ps.fi_sa;
      nosamix=[nosamix,'-',num2str(lp2),'-',num2str(lp1)];
      eval(['save ',nosamix])
   end 
icam_fr=0            %fa vedere i campi in frequenza
%' diss_nst', pausak
global GMA

fig=find(Gvet~=0);
gmi=min(min(Gvet(fig)));

if length(GMA)==0
 GMA=1e6;
end 
GMA=1e4/2;
GMA=gmi*20/vgconv;
%GMA=1e12;
clear ya aou gou fou figdisp
ikstop=13;
ikstop=0;

% 'qui', keyboard
 faGam=uL(2)/uL(1);
 vg=3e10/rr*faGam; 
 vgconv=vg;
isadiss=0;  %salvo per controllo trova soluzioni
%so=0;
%so=1e6;
%so=10;
% IMPORTANTE!!! fissa la soglia per analisi autovettori. Serve per eliminare l'effetto di forti picchi
%so=3;

autexp=1000;
ifiz=[];
if ~exist('ICON')==1
 ICON=0
end 
    if ICON==1 | ikstop>0
      hf1=figure; 
%            set(hf1,'pos',[  265   427   359   504])
            set(hf1,'pos',[ 363   380  321   339])
     h1=figure;
%      set(h1,'pos',[  12   424   246   504])
      set(h1,'pos',[  117   380   232   340])
      hf=figure; 
%            set(hf,'pos',[ 600   425   583   504])      
            set(hf,'pos',[ 684   380   463   351])      
    end 

al_shi=-2
al_shi=0
%pausak
%ifp=1
fso=[];
global del_Frea
    del_Frea=0;
if ~exist('icam_fr')==1
 icam_fr=0;
end
if ifp==-10
 icontrc=1;
else
 icontrc=0;
end
if ifp>=-1
 icontr=1;
else
 icontr=0;
end
icontrc=1;
icontr=0;
icontrc=ICON;
icontr=ICON;
%icontr=2;         % mostra ricerca zeri
if ifp==-4
 icontr=0;
end
alMAX=1e15;
%alMAX=200;
global isoga
isofg=0;
if length(isoga)==0
 isoga=0;
% isoga=1;
end 
isoga=1
% disoga=1
if ifp==-10
%disp(' attenzione !!!!!!!!!!!: isoga = 1 in diss_nst1.m ')
%keyboard
end
ierr=0;
dus=.7;
%dus=.95;
%dus=.5;
dus=.0;
%dus=.9;
%dus=.7;


% icontr=1;
% icontrc=1;
dpeak0=1.5;
%dpeak0=30;
dpeak0=10;
dpeak0=100;
dpeak0=1000;
inorm=0;
kex=1;
%' GMA ', keyboard
%GMA=1e4;
%if inorm==1
% vg=3e10/rqw;
%else
% vg=1e9;
%end

global alpha_th
if length(alpha_th)==0
 alpha_th=0;
end

pos=20;
vg1=1e13;
if exist('iconfina')==0
 iconfina=0;
end
if exist('iins')==0
 iins=1;
iins=0;
end
iins=0;
if exist('icampi')==0
 icampi=0;
end
if exist('isolut')==0
 isolut=1;
end
clear he
remi=.05;
rema=.95;
clear nsv ipuv fsov gsov asov tsov Amv Am1 Am Am1d Amvca Amvme Amr Amsort
clear gamsov gamosov XE CE CEv XEv CEu CEuv
clear rtetmsov mazisov mradsov polasov polratv

%keyboard

s=size(Gvet);

if length(s)==2
 siz=s(1:2);
  sz=size(Azvet);

%  sz(1)=25;
  Anu=reshape(Azvet(:,:,nAris,:),sz([1 2 4]));
  
%  Anu=reshape(Azvet(1:25,:,nAris,:),sz([1 2 4]));
  sA=size(Anu);
  a1=sA(1);
  a2=sA(2);
  a3=sA(3);
else
 disp(' errore size(Gvet) in diss_new '), keyboard
end

% figure, for k=1:length(Gvet), plot(abs(Azvet(:,k,2,1))),
% text(50,max(abs(Azvet(:,k,2,1))),num2str(Gvet(k,1),'%0.5g')), pausak, end

its=0;
fi0=find(Fint==0);
if length(fi0)~0
 Fint(fi0)=1e-10;
end

Fints=Fint;
sF=size(Fints);
if exist('pia0')==1
  pia=pia0;
end
if exist('istopi')==0
 istopi=0;
end
icontrsa=icontr;

pia=1;

for it=pia
% if istopi==1 & it==itsd
 if istopi==1 & it==its
  icontr=1;
 else
  icontr=icontrsa;
 end
 if exist('ISOm')
  ISO=ISOm(it,kex);
 end


if exist('nk1mat')==1
 if exist('alimati')==1
  alimi=alimati(it,kex);
 end
 if exist('alimatu')==1
  alim=alimatu(it,kex);
 end
 npk=nk1mat(it,kex);
end
 its=its+1;
 clear fov gov alv Am tetm pus
% pausak
Gvv=zeros(siz)*NaN;
avv=Gvv;
s=size(Gvet);
if sF(1)~=1
% npF=npFmat(it,kex);
% Fint=Fints(it,npF);
 Fint=Fints(it,:);
end
  if iLP==1
   dnum=numodi;
  else
   dnum=2*numodi;
  end

  Gve=Gvet;
%  alve=alvet;
  alve=alvet+al_shi;
  if exist('Kvet')
   Kve=Kvet;
  end
  npk=a1;

 Gvep=Gve;
 if iLP==0
  npk=npk/2;
 end
% npk
% keyboard

if exist('ikiautv')==1
 if length(ikiautv)>0
  ikiaut=ikiautv(it,kex);
 end
end

if exist('ikiaut')==0
 isopik=0;
else
 if length(ikiaut)==0
  isopik=0;
 else
  isopik=ikiaut;
 end
end
if iLP==1
 puA=[1:npk];
else
 puA=[1:2*npk];
end
if exist('ISO')==1
 if length(ISO)>0
   isopik=ISO;
 end
end
Gves=Gve;
alves=alve;
fi=find(Gve/vgconv>GMA);
Gve(fi)=0;
alve(fi)=0;
%keyboard
if isopik==0

isal=1;
if isal==0
Anus=Anu;
pak=1;
for k=1:pak:length(Fint)

  p=Gve(:,k);
  [du,iso1]=sort(p);
  fima=find(du>0);
  iso=iso1(fima);
  lis=1:length(iso);
  Gve(lis,k)=Gves(iso,k);
  alve(lis,k)=alve(iso,k);
  Anu(:,lis,k)=Anus(:,iso,k);
   Adisp=Anu(:,:,k);
      figure
        semilogy(abs(Adisp)), 
        a=axis;
        a=[1 length(Adisp) 1e-2 max(max(abs(Adisp)))];
        axis(a)
        pausak    
        IS=input(' modo =')
        isv(k)=IS;
      save ISsav isv
end 
'fine 1', keyboard 
end

gsogv=[1e14 1e15 1e16 1e17 1e18 1e19];
%psog0=3*gsogv(1);
psog0=gsogv(4);
%psog=gsogv(4);
%gsog1=gsogv(3);
%gsog2=gsogv(3);
isol=0;
fiz=find(Fint==0);
if length(fiz)>0
 Fint(fiz)=1e-10;
end
pak=1;
ks=0;
fiP=find(diff(Pusas)>1);    
if length(fiP)>0
 fP=1:fiP(1);
 fP=fiP(1)+1:length(Pusas);
else
 fP=1:length(Pusas);
end
if length(fP)>length(Ar)
 fP=1:length(Ar);
end
 fP=1:length(Ar);
 if fPmax>0
  fP=1:fPmax;
 end 
% 'q1', keyboard
for k=1:pak:length(Fint)
% p=alve(:,k);
if k<=kMO
iWidth=1;
else
iWidth=0;
end
%'iWidth', keyboard
if ifp~=-4
k
end
ks=ks+1;
 %p=Gve(:,k).*alve(:,k);
 %[du,iso]=sort(-p);
 %g=Gve(iso,k);
 %al1=alve(iso,k);
% p=-du;
% fi1=find(g>0 & g/vg<GMA &  p>-psog0);
  p=Gve(:,k);
  [du,iso]=sort(p);
  g=Gve(iso,k);
 al1=alve(iso,k);
 fi1=find(g>0 & g/vg<GMA);
% keyboard
% fi1=find(g>0 & p>-psog0);
 gv=g(fi1);
 al=al1(fi1);
 Ami=Anu(:,:,k);
 Am1(:,1:length(fi1))=Ami(:,iso(fi1));
 lfi1=length(fi1);
%  'q2', keyboard

 if lfi1>0 & k<length(Fint)
 ifi=1;
 while ifi<=lfi1
   Amr=Am1(puA,ifi);
   if sum(Amr)~=0
%    Amr=Amr/sqrt(abs(Amr'*Amr));
       if inor==1
        cn=sqrt(abs(Amr(fP)'*Amr(fP)));
        cn=median(abs(Amr(fP)));
       else 
        cn=sqrt(abs(Amr'*Amr));
        cn=median(abs(Amr));
       end 
       Amr=Amr/cn;    
   else
    Amr=Amr*0;
   end
%      'Amr', pausak

   Amrtm=Am1(:,ifi);
%   size(Amr)
%pausak
   if iLP==0
    l=fix(length(Amrtm)/2);
    le=1:l;
    lm=l+1:l*2;
    A=abs(Amrtm);
    Ate=A(le);
    Atm=A(lm);
    se=sum(Ate);
    sm=sum(Atm);
    if se~=0
     rem=se/(se+sm);
    else
     rem=0;
    end
   else
    rem=0.5;
   end
      nonval=0;
      if iLP==0 & nube==0
       if (pola==1 & rem>rema) | (pola==-1 & rem<remi)
        nonval=1;
%         disp(' dissu3 nonval ')
%         keyboard
       end
      end
%           disp(' dissu3 nonval '), keyboard 
  if nonval==0
   isol=isol+1;
   tetm(isol,ks)=rem;
   pus(isol,ks)=iso(fi1(ifi));
   fov(isol,ks)=Fint(k);
   gov(isol,ks)=gv(ifi);
%   [isol ks Fint(k)]
%   pausak
   alv(isol,ks)=al(ifi);
%%   Am(:,isol,ks)=Amr;
   ks1=ks;
   [du,idu]=max(abs(Amr));
   if idu>npk
    ipsci=[1:idu-npk-fix(npk/20)];
    ipscu=[idu+1-npk+fix(npk/20):npk];
   else
    ipsci=[1:idu-fix(npk/20)];
    ipscu=[idu+1+fix(npk/20):npk];
   end
   ipsc=[ipsci ipscu];
   if iLP==0
    ipsc=[ipsc ipsc+npk];
   end
%   'pa',pausak
%fP=find(diff(Pusc)>1);
%PUac=1:fP(1);
%   rppr=peAu(Amr);  
 if iWidth==1
   rppr=sumFsemp(Amr,Mvefm0,Mvefp0,Mvegm0,Mvegp0,besp,besm,xdx1,xdx3);      
 end   
   for ki=k+pak:pak:length(Fint)
    ks1=ks1+1;
    p=Gve(:,ki).*alve(:,ki);
    gd=Gve(:,ki);
    ald=alve(:,ki);
    fi2=find(gd>0);
    if length(fi2)>0
     gv2=gd(fi2);
     alv2=ald(fi2);
     Amir=Anu(:,:,ki);   
     clear Am1d
     Am1d(:,1:length(fi2))=Amir(:,fi2);
     lfi2=length(fi2);
     if lfi2>0
      clear autc autc1 autcp autca dpea
      lF=1:lfi2;
%    'prima', keyboard
  
     for ifii=lF
     
      Amix=Am1d(puA,ifii);
     if iWidth==1
      rppx=sumFsemp(Amix,Mvefm0,Mvefp0,Mvegm0,Mvegp0,besp,besm,xdx1,xdx3);      
     end
%      rppx=peAu(Amix);

      if sum(Amix)~=0
       if inor==1
        cn=sqrt(abs(Amix(fP)'*Amix(fP)));
       else 
        cn=sqrt(abs(Amix'*Amix));
       end 
       cn=median(abs(Amix));
       Amix=Amix/cn;
      else
       Amix=Amix*0;
      end
%      'Amix', pausak
      autca(ifii)=abs((filtA(Amix,so)))'*abs((filtA(Amr,so)));
   %   autc(ifii)=log10(abs(filtA(Amix,so)))'*log10(abs(filtA(Amr,so)));
      autc1(ifii)=(filtA(Amix,so))'*(filtA(Amr,so)); 

%      autdiff=sqrt(sum((abs(Amix(fP))-abs(Amr(fP))).^2)/length(fP));
%      autc(ifii)=1./autdiff;
      DIF=(abs(Amix(fP))-abs(Amr(fP)));
      DIF=(log10(abs(Amix(fP)))-log10(abs(Amr(fP))));
      autdiff=sqrt(sum(DIF.^2)/length(fP));
      autc(ifii)=1./autdiff;
     if iWidth==1
      dpea(ifii)=abs(rppx/rppr);
     end
%      nuo=1./sqrt(sum(abs(1-Amix./Amr).^2)/length(Amix))
% figure, plot(xla,abs(Amix),xla,abs(Amr),'w'), pausak
%      autca(ifii)=abs(filtA(Amix,so))'*abs(filtA(Amr,so));
%      autc(ifii)=(filtA(Amix,so))'*(filtA(Amr,so));
     if   length(find(ki-ikstop==0))==1
           xla=[1:length(Amix)];
       if ~exist('figdisp')
         figdisp=figure;
         set(figdisp,'pos',[12     9   625   366])
       else
        figure(figdisp)
       end
       semilogy(xla,abs(Amix),xla,abs(Amr),'w.-',xla,abs(abs(Amr)-abs(Amix)),'r.-'), 
       expla=' White: reference, Red: difference ';
       title([expla,' ifii ',num2str(ifii),';          difference ',num2str(autdiff)]), pausak     
     end


%      autc(ifii)=Amix(fP)'*Amr(fP);

%      autc(ifii)=log10(Amix(fP)')*log10(Amr(fP));
  %    autc(ifii)=(Amix)'*(Amr);
%      autc1(ifii)=Amix(ipsc)'*Amr(ipsc);
     end
%     autcp=autc.*autc1;
%     [du,fiz]=max(abs(autcp));
    
%     [du,fiz]=max(abs(autc));
%     [du,fiz0]=sort(1-abs(autc));
%     [du,fiz0]=sort(1-abs(autca));
%     du=1-du;
%     fiz=fiz0(1:2);
%     fiz=fiz0(1);
Tagg=0;
       Aexp=ones(size(Ar))*NaN;     
inew=iexpe;
if inew==1
      kg0=2;
      if ki>3
     
       if ki<kg0+2
        kg=ki-2;
        puf=1:ki-1;
        pufu=1:ki;
       else
        kg=kg0;
        puf=1:ki-1;
        puf=puf(end-kg0+1:end);
        pufu=[puf puf(end)+1];
       end      
      fprec=fov(isol,puf)*lambda*1000;
      aprec=alv(isol,puf);
      gprec=gov(isol,puf);
      fpi=Fint(pufu)*lambda*1000;
      fpid=Fint(pufu);
%      'cont fpo', keyboard

       coa=polyfit(fprec,aprec,kg);
       cog=polyfit(fprec,gprec,kg);
       aexp=polyval(coa,fpi);
       gexp=polyval(cog,fpi);
       Ae=aexp(end);
       Ge=gexp(end);
       Tagg=1./[(abs(alv2/Ae))+ 10*abs(log10(gv2/Ge))]';
%      Tagg=0;
 %     Tagg=1./[(abs(alv2-Ae))]';
%       Tagg=1./[(abs(alv2-Ae))*100+ abs(log10(gv2/Ge))]';
 %      Tagg=1./[abs(log10(gv2/Ge))]';
%        Aexp=fitAut(Am(:,puf(end-1:end)),fprec(end-1:end),fpi(end),1);
        
      else
       Ae=0;
       Ge=0;      
       fpid=[];
       aexp=[];
       gexp=[];      
       Aexp=ones(size(Ar))*NaN;      
       Tagg=0;

      end

end
     if iWidth==1
       fi=find(dpea==0);
       dpea(fi)=1e-2;
       Tagg=1./dpea;
     end 
%      'fine Tagg', keyboard       

  autc=autc/mean(abs(autc));
%  if Tagg~0
%   Tagg=Tagg/mean(Tagg);
%  end
 % autu=autc/max(autc)+Tagg/max(Tagg);
 if iWidth==1
  autu=Tagg;
 else
  autu=autc+Tagg;
 end
%      'fine Tagg', keyboard       
 
     [du1,ias]=sort(1-abs(autu));
%     [du1,ias]=sort(1./autu);
%     'fine', keyboard
%     [du1,ias]=sort(1-abs(autc));
%     [du1,ias]=sort(1-abs(autca));
     ipma=min([length(ias) 5]);
     iasd=ias(1:ipma);
     icok=1;
     clear Adis
%     for kau=iasd
     for kau=1:length(iasd)
      Adum=Am1d(:,kau);
%      Adis(:,icok)=Adum/sqrt(abs(Adum'*Adum));
      Adis(:,icok)=Adum/median(abs(Adum));
      icok=icok+1;
     end
     [du1,fiz1]=max(abs(autca));
%     [du,fiz]=max(abs(autc));
     [du,fiz]=max(abs(autu));
   %  [du,fiz]=max(abs(autca));
      gu1=gv2(fiz1);
      gu=gv2(fiz);
%      [gu gu1], pausak
%presto
%      if gu1<gu
%       fiz=fiz1;
%       du=du1;
%      end

     Amix0=Am1d(:,fiz);
     for kso=1:length(fiz)
%       Amix0(:,kso)=Amix0(:,kso)/sqrt(abs(Amix0(:,kso)'*Amix0(:,kso)));
       Amix0(:,kso)=Amix0(:,kso)/median(abs(Amix0(:,kso)));
     end  
%'Amiz', keyboard
%     Amix0=Am1d(:,fiz);
     
     if  icontrc==1 | ki==ikstop
      icontr=2;
      disp('[k isol ki du]')
      [k isol ki du]
     figure(h1);
%     h1=figure;
%      set(h1,'pos',[  12   424   246   504])
%      plot(lF,abs(autc),lF,abs(autc1),lF,abs(autcp),lF(fiz),abs(autc(fiz)),'*'), pausak
      plot(lF,abs(autc),lF(fiz),abs(autc(fiz)),'*',lF,Tagg,lF,abs(autu),'w'), 
      xla=[1:length(Amix)];
      figure(hf1)
      clg
      subplot(211),
       if isol>1
        fovp=fov(1:isol-1,:)'*lambda*1000;
        govp=gov(1:isol-1,:)'/vgconv;
        alvp=alv(1:isol-1,:)';
        plot(fovp,alvp,'--'), hold on
       end
      plot(fov(isol,1:ki-1)*lambda*1000,alv(isol,1:ki-1),'.-',Fint(ki)*lambda*1000,alv2(fiz),'wo','linewidth',2),
      hold on, 
      if inew==1, plot(fpid*lambda*1000,aexp,'w--','linewidth',2), end
      plot(Fint(ki)*lambda*1000,alv2,'x'), grid
      subplot(212)
       if isol>1
        semilogy(fovp,govp,'--'), hold on
       end
      semilogy(fov(isol,1:ki-1)*lambda*1000,gov(isol,1:ki-1)/vgconv,'.-',Fint(ki)*lambda*1000,gv2(fiz)/vgconv,'wo','linewidth',2), 
      hold on, 
      if inew==1, semilogy(fpid*lambda*1000,gexp/vgconv,'w--','linewidth',2),       end
      semilogy(Fint(ki)*lambda*1000,gv2/vgconv,'x'), grid, 
      figure(hf)
%      semilogy(xla,abs(Adis),xla,abs(Amr),'wo',xla,abs(Amix0),'wx') %,xla,abs(Adis(:,fiz)),'rx'), 
      semilogy(xla,abs(Adis),xla,abs(Amr),'wo',xla,abs(Amix0),'wx',xla,abs(Aexp),'ro'), 
      title([' o: previous,  x: next,     point #  ',num2str(ki)])
      a=axis;
%      a=[fP(1) fP(end) 1e-3 1];
      a=[fP(1) length(Adis) 1e-2 max(max(abs(Adis)))];
      axis(a)
%      pausak
      if isnan(Aexp(1))==0
            autc(fiz)
            ' confronto valore stimato con vero '
            
             autexp=1./sqrt(sum(abs(Amix0(fP)-Aexp(fP)).^2)/length(fP))
             
      end
%      if autc(fiz)>autexp
pausak
       ifiz=input(' ifz (enter se va bene) ');
%      else
%       pausak
%      end
      if length(ifiz)~=0 
       fiz=ifiz;
      end
%      close(h1:hf)
     end
       if ks1==2
        Amsort(:,1,isol)=Amr*median(abs(Amr));
        Amsort(:,1,isol)=Amr;
       end 
        Amsort(:,ks1,isol)=Amix0*median(abs(Amix0));
        Amsort(:,ks1,isol)=Amix0;
     difacc=max(diff(abs(Amr)));
%     [du difacc], pausak
if ks1<3
 aprec=alv2(fiz)/2; 
end
if iNObre==0
 condibre=du>dus & difacc<dpeak0 & alv2(fiz)>mean(aprec);
else
 condibre=du>dus & difacc<dpeak0;
end
     if condibre
%     if du>dus & difacc<dpeak0 
%     if difacc<dpeak0
%     [du difacc], pausak

      Amr=Am1d(puA,fiz);
%       cn=sqrt(abs(Amr(fP)'*Amr(fP)));
       if inor==1
       cn=median(abs(Amr(fP)));
       else
       cn=median(abs(Amr));
%       cn=sqrt(abs(Amr'*Amr));
       end
       Amr=Amr/cn;
%      Amr=Amr/sqrt(abs(Amr'*Amr));

   [du,idu]=max(abs(Amr));
   if idu>npk
    ipsci=[1:idu-npk-fix(npk/20)];
    ipscu=[idu+1-npk+fix(npk/20):npk];
   else
    ipsci=[1:idu-fix(npk/20)];
    ipscu=[idu+1+fix(npk/20):npk];
   end
   ipsc=[ipsci ipscu];
   if iLP==0
    ipsc=[ipsc ipsc+npk];
   end
      A=abs(Amr);
      if iLP==0
       Ate=A(le);
       Atm=A(lm);
       se=sum(Ate);
       sm=sum(Atm);
       rem=se/(se+sm);
      else
       rem=0.5;
      end
       tetm(isol,ks1)=rem;
       pus(isol,ks1)=fi2(fiz);
       gov(isol,ks1)=gv2(fiz);
       alv(isol,ks1)=alv2(fiz);
     
       fov(isol,ks1)=Fint(ki);

       Gve(fi2(fiz),ks1)=0;
%       disp(' Gve '), pausak
%       pausak
     else
%'     prima di brak', keyboard
      break
     end
    end %if lfi2
   end %if fi2
   end %ki
%   disp(' ki')
%   pausak
  end  %nonval
 ifi=ifi+1;
 end  %ifi
 end   %if
%   pausak
% [gv gs]
 
end
%' fine ordine', keyboard
% figure, fcl=1000*lambda; subplot(121), plot(Fint*fcl,alvet'), grid,
% vg=3e10/rqw;        subplot(122), semilogy(Fint*fcl,Gvet'/vg)
%



s=size(alv);
% if isolut==1
%  nummo=s(1);
% elseif isolut==0
%  nummo=1;
% else
%  nummo=isolut;
% end
  nummo=s(1);
else
 psou
 if isolut==1
  nummo=2;
 else
  nummo=1;
 end
end %sol picchi

ics=0;
clear fou gou aou Amu tetmu pou
for is=1:nummo
 fi=find(alv(is,:)~=0 );
 fi1=find(alv(is,:)==0);
% fi=find(fov(is,:)~=0 );
% fi1=find(fov(is,:)==0);

% if icontr==1
%  figure, plot(fov(is,fi),alv(is,fi)), pausak
% end
%pausak
%fi, pausak
 if length(fi)>2
% if length(fi)==length(Fint)
  ics=ics+1;
   fou(ics,:)=fov(is,:);
   pou(ics,:)=pus(is,:);
   gou(ics,:)=gov(is,:);
   aou(ics,:)=alv(is,:);
   tetmu(ics,:)=tetm(is,:);
%   Amu(:,ics,:)=Am(:,is,:);
   fou(ics,fi1)=fov(is,fi1)*NaN;
   gou(ics,fi1)=gov(is,fi1)*NaN;
   aou(ics,fi1)=alv(is,fi1)*NaN;
   tetmu(ics,fi1)=tetm(is,fi1)*NaN;
 end
end
if exist('isoga')==0
 isoga=0;
end
%' qui', keyboard
if ~exist('gou')
 figure, fcl=1000*lambda; subplot(121), plot(Fint*fcl,alvet'), grid,
         subplot(122), semilogy(Fint*fcl,Gvet')
 'manca gou !!!!!!!'
 keyboard
 keyboard
end

sga=size(gou);
%' pausak'
if min(sga)>1
clear meg
 if isoga>=1
  if isoga==2 
   prga=abs(gou.*aou); 
   fina=find(isnan(prga(:,1))==0); 
  else 
   fina=find(isnan(gou(:,1))==0); 
  end 
  for kca=1:sga(1)
   if isoga==2  
    fina=find(isnan(prga(kca,:))==0); 
    meg(kca)=mean(prga(kca,fina));
   else 
    fina=find(isnan(gou(kca,:))==0); 
    meg(kca)=mean(gou(kca,fina));
   end 
  end 
%  [du,ics]=sort(mean(gou,2));
  [du,icsd]=sort(gou(:,1));
%'wui meg', keyboard
%  [du,icsd]=sort(meg);
%   ics=fina(icsd);
   ics=icsd;
   fou=fou(ics,:);
   pou=pou(ics,:);
   gou=gou(ics,:);
   aou=aou(ics,:);
   tetmu=tetmu(ics,:);
%      fou=fou(:,ics);
%      pou=pou(:,ics);
%      gou=gou(:,ics);
%      aou=aou(:,ics);
%      tetmu=tetmu(:,ics);
 end
  fia=find((-aou(:,1))>0);
  sa=size(aou);
  aoud=aou;
  foud=fou;
  goud=gou;
  poud=pou;
  tetmud=tetmu;
  aou=[];
  fou=[];
  gou=[];
  pou=[];
  tetmu=[];

  iaoud=0;
  if sa(1)>1
   for isa=1:sa(1)
    fitum=find((aoud(isa,:))>alpha_th);
    fitup=find((aoud(isa,:))<alpha_th);
%    fitum=find((aoud(isa,:))>0);
%    fitup=find((aoud(isa,:))<0);
%     if length(fitup)*length(fitum)>0 | ( min(abs(aoud(isa,:)))<1 & length(find(goud(isa,:)/vg<GMA/2))==length(Fint) )
     if length(fitup)*length(fitum)>0
      iaoud=iaoud+1;
      aou(iaoud,:)=aoud(isa,:);
      gou(iaoud,:)=goud(isa,:);
      fou(iaoud,:)=foud(isa,:);
      pou(iaoud,:)=poud(isa,:);
      tetmu(iaoud,:)=tetmud(isa,:);
     end
   end
  end
%  disp('key'),  keyboard
  saou=prod(size(aou));

%    ' algag  prewa' , keyboard
%  ex_min=
  del_Frea=0;
  if saou==0 & length(fso)==0
    maal=max(max(alv(1,:)));
    mial=min(min(alv(1,:)));
    mave=[mial maal];
%    [du,iaoup]=min(abs(mave));
    al_shi=-mean(mave);
    al_shi=-mean(mial*1.2);
%    ' algag' , keyboard
%    if iaoup==1
%     al_shi=-daal;
%    else
%     al_shi=daal;
%    end
    Fi1=Fint(1:size(alv,2));
    coau=polyfit(Fi1,alv(1,:),1);
    del_Frea=al_shi/coau(1);
    if ifp==-10
     ' algag' , keyboard
    end
    diss_las
    return
  end
  if saou==0 & length(fso)==0
   iiplot=0;
%   if ifp~=-4
   iiplot=1;
%   end
   if iiplot==1
   figure,
   subplot(121), plot(fov',alv',fov(1,:)',alv(1,:)','wo')
   subplot(122), semilogy(fov',gov',fov(1,:)',gov(1,:)','wo')

   disp(' ')
   disp(' ')
   disp(' ')
   disp(' ******************  Wrong frequency region to find a mode ')
   disp(' ')
   disp(' --->  --->    --->  Restart with appropriate frequency range ')
   disp(' ')
   disp(' ')
   disp(' ')
   end %iiplot
   fb=fov(1,:);
   ab=alv(1,:);
   gb=gov(1,:);
%   figure, plot(fb,ab)
  if ireset_int>=0
   proble=100;
    meao=mean(ab);
    if meao>0
     sg_add=-1;
    else
     sg_add=1;
    end
    if ireset_int==0
     sg_add=sg_add/2;
    end
    st_add=sg_add*abs(diff(Dlam_mod(1:2)));
    if is_even(mmint)==1
      Dlao(4)=Dlao(4)+st_add;
      Frisi=(Dlao(1)+Dlao(4))/fala;
      Frisu=(Dlao(2)+Dlao(4))/fala;
    else
      Dlao(1:2)=Dlao(1:2)+st_add;
      Frisi=(Dlao(1))/fala;
      Frisu=(Dlao(2))/fala;
    end
   end %ireset_int
   if ifp~=-4
    keyboard
   end
   'ret', keyboard
   return
  end
%   'ret 1', keyboard  
  if saou==0, return, end
  fia=find((-aou(:,1))>0);
  if iaoud>1 & isoga==0
   [du,icsa]=sort(-aou(fia,1));
    ics=fia(icsa);
    fou=fou(ics,:);
    pou=pou(ics,:);
    gou=gou(ics,:);
    aou=aou(ics,:);
    tetmu=tetmu(ics,:);
   end
  ics=length(fia);

end


if iins>=4
 fov=fou;
 pus=pou;
 gov=gou;
 alv=aou;
 tetm=tetmu;
end
%'qui', keyboard
s=size(aou);
nma0=s(1);
 if isolut==1
  nma=s(1);
 elseif isolut==0
  nma=1;
 else
  nma=abs(isolut);
 end
 if nma0<nma
  nma=nma0;
 end

 if ~exist('nmasce')
  nmascel=1:nma;
 else
  if nma>nmasce
   nmascel=1:nmasce;
  else
   nmascel=1:nma;
  end
 end
 nmascel=1:nma;
 fcov=1000*lambda;
 if ifp==-10
  figure, subplot(211), plot(fcov*fou',aou','.-'), grid
  subplot(212), semilogy(fcov*fou',2*gou'/vgconv,'.-'), grid
  title(' Curves to be seached ')
  pausak
 end
if icontr>=1 & iins>=1
% p=aou;
% figure, plot(fou',p')
% axis([min(Fint) max(Fint) -100 100]), pausak
% figure, semilogy(fou',gou'/vg)
% pausak

 p=aou;
 figure, subplot(211), plot(fou',p','.-'), grid
% axis([min(Fint) max(Fint) -50 50]),
 subplot(212), semilogy(fou',gou'/vg,'.-')
 pausak
 if istopi==1
  icambia=input(' icambia = [0/1] ');
  if icambia==1
   pud=pou(1,6);
   pou(1,6)=pou(2,6);
   pou(2,6)=pud;
   pud=gou(1,6);
   gou(1,6)=gou(2,6);
   gou(2,6)=pud;
   pud=aou(1,6);
   aou(1,6)=aou(2,6);
   aou(2,6)=pud;
   alv=aou;
   gov=gou;
   p=aou;
   figure, subplot(211), plot(fou',p','.-'), grid
   subplot(212), semilogy(fou',gou'/vg,'.-')
   pausak
  end
 end
% figure
% sA=size(Amu);
% lF=sA(3)-1;
%' passo ?? ', keyboard
 for nm=1:nma
% for nm=2:nma
  nm
  clear V
  fip=find(pou(nm,:)~=0);
  for ifip=1:length(fip)
   V(:,ifip)=abs(Anu(:,pou(nm,fip(ifip)),fip(ifip)));
  end
  if icontrc>=1
   figure, plot(abs(V)), pausak
  end
%  disp(' Ami')
%  pausak
%  V=abs(reshape(Amu(:,nm,:),sA(1),sA(3)));
  lF=length(fip)-1;
  pd=p(nm,:);
  lF1=fip(1:lF);
  lF2=fip(2:lF+1);
  lFt=fip(1:lF+1);
  fiz=find(pd(lF1).*pd(lF2)<0 & diff(pd(lFt))>0 );
  fiz1=find(abs(pd(lFt))<alMAX);
%  fiz1=find(abs(pd(lFt))<2e7);

  if (length(fiz)>0) | (length(fiz1)>0 & iins>=1)
   if length(fiz)==0
    [du,is]=sort(abs(pd(lFt)));
    if length(is)>5
     ps=sort(is(1:5));
    else
     ps=sort(is);
    end
    psV=ps;
    ps=lFt(ps);

   else
    if length(lFt)>2
     ps=fiz(1)+[-1 0 1];
    else
     ps=fiz(1)+[0 1];
    end
     if ps(1)<1
      ps=ps+1;
     end
     if ps(length(ps))>lF+1 & ps(1)>1
      ps=ps-1;
     end
     psV=ps;
     ps=lFt(ps);
   end
%   if mean(gou(nm,ps))/vg<GMA

   subplot(311),
   plot(fou(nm,:)',p(nm,:)',fou(nm,ps)',p(nm,ps)','r*'), grid
   subplot(312),
   semilogy(fou(nm,:)',gou(nm,:)'/vg,'m',fou(nm,ps)',gou(nm,ps)'/vg,'r*'), grid
%   l1=min(p(nm,ps)); l2= max(p(nm,ps));
%   dl=l2-l1;
%   axis([min(Fint) max(Fint) l1-dl l2+dl]),
   subplot(313), plot(V(:,psV)),
%   pausak
   if iins>=4
   iacc=input(' accetto soluzione? [1]  ')
   if isempty(iacc)==1
    iacc=0;
   end
    if iacc==0
     fov(nm,:)=fov(nm,:)*0;
    end
   end
%   end  %<GMA
  else
    fov(nm,:)=fov(nm,:)*0;
  end
%  pausak
 end
end

if iins>=4
ics=0;
clear fou gou aou Amu tetmu pou
for is=1:nma
% for is=2:nma
% fi=find(fov(is,:)~=0 & gov(is,:));
 fi=find(fov(is,:)~=0 );
 fi1=find(fov(is,:)==0);
% if icontr==1
%  figure, plot(fov(is,fi),alv(is,fi)), pausak
% end
 if length(fi)>1
  ics=ics+1;
  fou(ics,:)=fov(is,:);
  pou(ics,:)=pus(is,:);
  gou(ics,:)=gov(is,:);
  aou(ics,:)=alv(is,:);
  tetmu(ics,:)=tetm(is,:);
%  Amu(:,ics,:)=Am(:,is,:);
  fou(ics,fi1)=fov(is,fi1)*NaN;
  gou(ics,fi1)=gov(is,fi1)*NaN;
  aou(ics,fi1)=alv(is,fi1)*NaN;
  tetmu(ics,fi1)=tetm(is,fi1)*NaN;
 end
end
end %iins

%s=size(aou);

nma=ics;
tso=[];
gso=[];
M2l=[];
mradv=[];
maziv=[];
rtetmv=[];
polcav=[];
polrat=[];
aso=[];
gamso=[];
gamoso=[];
XE=[];
CE=[];
CEu=[];
fso=[];
ipu=[];
ns=[];
pak1=1;
Anso=[];
Acso=[];
Amso=[];
sMem=size(aou);
%keyboard

%if s(1)<max(nmascel)
% nmascel=1:s(1);
%end

%disp('dissu3'), keyboard
nsc=0;
Fint1=Fint(1:pak1:length(Fint));
%for nso=1:nma
%' nmasce ', keyboard
save sacont

for nso=nmascel
clear Fi1du
%disp(' acktung !!!!!!!!!!!! dissus '), keyboard
% for nso=2:nma

   ps=1:pak1:sMem(2);
   ya1=aou(nso,ps)-alpha_th;

%&&&&&&&&&&&&&
   fi=find(isnan(ya1)==0);
%   ya1=aou(nso,fi);
   ya1=ya1(fi);
%   ta=tetmu(nso,fi);
%   ga=gou(nso,fi);
%   Fi1=Fint1(fi);
   ta=tetmu(nso,ps);
   ta=ta(fi);
   ga=gou(nso,ps);
   ga=ga(fi);
   Fi1=Fint1(fi);
%   'ga', keyboard
   fisa=fi;
  if icontr==2
   h1=figure;
   subplot(121), plot(Fi1,ya1), grid,
   subplot(122), semilogy(Fi1,ga), grid,
  end
 dya1=[1 diff(ya1)./diff(Fi1)];
 dya1p=[diff(ya1)./diff(Fi1) 1];
 ly=[1:length(dya1)-1];
 dya2=[1 ya1(ly).*ya1(ly+1)];
% f10v=find(dya1>0 & dya1p>0 & dya2<0)
 f10v=find(dya1>0 & dya2<0);
% f10v=find(dya2<0);
  if length(f10v)>0
     fpu=f10v(1);
     fpA=fpu+[-2:2];
     fiA=find(fpA>=1 & fpA<=length(Fint));
     fpA=fpA(fiA);
     Adispla=abs(Amsort(:,fpA,nso));   
%     'A1', keyboard
  end

  if icontr==2
   figure(h1),
   subplot(121), hold on, plot(Fi1(f10v),ya1(f10v),'wo'), ,
   subplot(122), hold on, semilogy(Fi1(f10v),ga(f10v),'wo'), ,
   title(' controllo 5')
   disp(' controllo ')
   pausak
   close(h1)
  end


  ifnz=0;
  fsovi=[];
  if length(f10v)>1

   lumax=15;
%  lumax=270;
   if length(ya1)>lumax

    Fi1s=Fi1;
    ya1s=ya1;
    gasa=ga;

    zev=Fi1(f10v);

    gra=ceil(length(ya1)-1);
    pfo=linspace(Fi1(1),Fi1(end),100);

    co=polyfit(Fi1,log10(ga),gra);
    gf=10.^(polyval(co,pfo));
    g0v=10.^(polyval(co,zev));
    [du,igmi]=min(g0v);


    if icontr==2
     h1=figure;
     subplot(121), plot(Fi1s,ya1s,'.-',Fi1s(f10v),ya1s(f10v),'wo',zev,0,'r+'),
     subplot(122), semilogy(Fi1s,gasa,'.-',Fi1s(f10v),gasa(f10v),'wo'),
     disp(' controllo 10')
     pausak
    end
    ksi=0;
    if icontr==2
    hdu=figure;
    end
    clear zeu g0u
    for ks=1:length(zev)
     zel=zev(ks);
     [du,fi]=sort(abs(Fi1-zel));
     pu=[1:7];
     grl=length(pu)-1;
     fi=sort(fi(pu));
     fz=Fi1(fi);
     az=ya1(fi);

     co=polyfit(fz,az,grl);
     zed=roots(co);
     fia=find(abs(imag(zed))<=abs(real(zed)) & ...
         (real(zed)<fz(end) & real(zed)>fz(1)));
     if length(fia)>0
      ksi=ksi+1;
      zedu=zed(fia);
      [du,izeus]=min(abs(zedu-zel));
      zeus=zedu(izeus);
      zeu(ksi)=zeus;
      co=polyfit(fz,log10(ga(fi)),grl);
      g0k=10.^(polyval(co,zeus));
      g0u(ksi)=g0k;
       if icontr==2
        figure(hdu)
        subplot(121), plot(Fi1s,ya1s,'.-',Fi1s(f10v),ya1s(f10v),'wo',fz,az,'c.',zeus,0,'r*'),
        subplot(122), semilogy(Fi1s,gasa,'.-',Fi1s(f10v),gasa(f10v),'wo',fz,ga(fi),'c.',zeus,g0k,'r*'),
        disp(' controllo ')
        pausak
       end

     end
    end

    [g0s,ifs]=min(g0u);
    zes=zeu(ifs);

    ze=zes;
    gg0=g0s;
    [du,ivero]=min(abs(Fi1-ze));
    fsovi=ivero;
     
    if icontr==2
     close(hdu)
     h1=figure;
     subplot(121), plot(Fi1s,ya1s,'.-',Fi1s(f10v),ya1s(f10v),'wo',zeu,0,'r+',zes,0,'go'),
     subplot(122), semilogy(Fi1s,gasa,'.-',Fi1s(f10v),gasa(f10v),'wo',zeu,g0u,'r+',zes,g0s,'go'),
     disp(' controllo ULTIMO')
     pausak

     [du,igmi]=min(gasa(f10v));
     fpu=f10v(igmi);
     fpA=fpu+[-2:2];
     fiA=find(fpA>=1 & fpA<=length(Fint));
     fpA=fpA(fiA);     
     Adispla=abs(Amsort(:,fpA,nso));
%          'A5', keyboard

%     figure, semilogy(Adispla), pausak

     close(h1)
    end

    fr_pun=[];
    for kze=1:length(zeu)
     zelo=zeu(kze);
     [du,ize]=min(abs(Fint-zelo));
     fr_pun=[fr_pun ize];
    end
%    'passo 1', keyboard
     clear As
     for iaf=1:length(Fint)
      Pun=pou(nso,iaf);
      if Pun>0
       Ada=reshape(Anu(:,Pun,iaf),sA(1),1);
      else
       Ada=zeros(sA(1),1);
      end
      [du,ima]=max(abs(Ada));
      As(:,iaf)=Ada*sign(real(Ada(ima)));
     end
    lep=size(pou);     
 for iaf=1:lep(2)
  Pun=pou(nso,iaf);
  if Pun>0
   Ada=reshape(Anu(:,Pun,iaf),sA(1),1);
  else
   Ada=zeros(sA(1),1);
  end
%  [du,ima]=max(abs(Ada));
%  As(:,iaf)=Ada*sign(real(Ada(ima)))/median(abs(Ada));
  As(:,iaf)=Ada/median(abs(Ada));
 end
 
 Apre=As(:,1);
 corp=Apre'*Apre;
for ksor=2:lep(2)
  Apo=As(:,ksor);
  corpi=Apre'*Apo;
  Apre=Apo;
  corp=corpi;
  if real(corpi)<0
   As(:,ksor)=-As(:,ksor);
   Apre=-Apo;
   corp=-corpi; 
  end
 end
 
%    field_fr

   else

    Fi1s=Fi1;
    ya1s=ya1;
    gasa=ga;

    if f10v(1)>2
     pumi=f10v(1)-2;
    else
     pumi=f10v(1)-1;
    end
    if f10v(end)<length(Fi1)
     puma=f10v(end)+1;
    else
     puma=f10v(end);
    end
    puf=pumi:puma;
    Fi1=Fi1(puf);
    ya1=ya1(puf);
    ga=ga(puf);
%    'ga2', keyboard


    gra=ceil(length(ya1)-1);
    co=polyfit(Fi1,ya1,gra);
    pfo=linspace(Fi1(1),Fi1(end),100);
    af=polyval(co,pfo);
    zed=roots(co);
    fia=find(imag(zed)==0 & (real(zed)<Fi1(end) & real(zed)>Fi1(1)));

    if icontr==2
     figure;
     plot(Fi1s,ya1s,Fi1,ya1,'w.',pfo,af,zed(fia),0,'ro'),
%     plot(Fi1s,ya1s,Fi1,ya1,'w.',pfo,af,Fi1s(f10v),0,'ro'),
    end

    df=diff(Fi1(1:2));
    if length(fia)==0
     fia=find(imag(zed)==0 & (real(zed)<Fi1(end)+df & real(zed)>Fi1(1)-df ));
    end
    zev=zed(fia);

%     figure, plot(Fi1,ya1,'g',pfo,af,'r',zev,zev*0,'w*'), pausak

    co=polyfit(Fi1,log10(ga),gra);
    gf=10.^(polyval(co,pfo));
    g0v=10.^(polyval(co,zev));
    [du,igmi]=min(g0v);
     fpu=zev(igmi);
    [du,fpu]=min(abs(fpu-Fint));
     fpA=fpu+[-2:2];
     fiA=find(fpA>=1 & fpA<=length(Fint));
     fpA=fpA(fiA);     
     Adispla=abs(Amsort(:,fpA,nso));
     if ifp==-10
     'A6', keyboard
     end
%     figure, plot(Fi1,ga,'c',ze,gg0,'wo',pfo,gf,'m',zev,g0v,'w*'), pausak
    if icontr==2
     h1=figure;
     subplot(121), plot(Fi1s,ya1s,Fi1s(f10v),ya1s(f10v),'wo',zev,0,'r+'),
     subplot(122), semilogy(Fi1s,gasa,Fi1s(f10v),gasa(f10v),'wo'),
     disp(' controllo questo ')
     pausak
    end
    ksi=0;
    for ks=1:length(zev)
     zel=zev(ks);
     [du,fi]=sort(abs(Fi1-zel));
     pu=[1:3];
     grl=length(pu)-1;
     fi=sort(fi(pu));
     fz=Fi1(fi);
     az=ya1(fi);
     co=polyfit(fz,az,grl);
     zed=roots(co);
%     fia=find(imag(zed)==0 & (real(zed)<fz(end) & real(zed)>fz(1)));
     fia=find(abs(imag(zed))<=abs(real(zed)) & ...
         (real(zed)<fz(end) & real(zed)>fz(1)));
     if length(fia)>0
      ksi=ksi+1;
      zedu=zed(fia);
      [du,izeus]=min(abs(zedu-zel));
      zeus=zedu(izeus);
      zeu(ksi)=zeus;
      co=polyfit(fz,log10(ga(fi)),grl);
      g0k=10.^(polyval(co,zeus));
      g0u(ksi)=g0k;
     end
    end
    if icontr==2
     figure(h1)
     subplot(122), hold on, plot(zeu,g0u,'r+'),
     disp(' controllo ')
     pausak
     close(h1)
    end

    [g0s,ifs]=min(g0u);
    zes=zeu(ifs);

    fr_pun=[];
    for kze=1:length(zeu)
     zelo=zeu(kze);
     [du,ize]=min(abs(Fint-zelo));
     fr_pun=[fr_pun ize];
    end
%    'passo 2', keyboard

%    field_fr


    ze=zes;
    gg0=g0s;
    [du,ivero]=min(abs(Fi1-ze));
    fsovi=ivero;
%   else
%    [du,imi]=min(ga(f10v));
%    fsovi=f10v(imi);
   end

  elseif length(f10v)==1
   fsovi=f10v;
  elseif length(f10v)==0
   fsovi0=find(abs(ya1)<alMAX);
%   fsovi0=find(abs(ya1)<2e7);
   npma=2;
   if length(fsovi0)>=npma
    fsovi=[];
     [du,idu]=min(abs(ya1));
     fsovi00=ya1(idu);
     if fsovi00>0
      ifnz=-1;
     else
      ifnz=1;
     end
   end
  end

sau=size(aou);
%'sau', keyboard
if sau(1)==1
 icas=length(find(aou>0))==length(Fint) | length(find(aou<0))==length(Fint);
 ya1=aou;
 Fi1du=Fint(1:length(ya1));
% 'Fi1 3', keyboard
 dya1=[1 diff(ya1)./diff(Fi1du)];
 dya1p=[diff(ya1)./diff(Fi1du) 1];
 ly=[1:length(dya1)-1];
 dya2=[1 ya1(ly).*ya1(ly+1)];
 f10v=find(dya1>0 & dya2<0);

 if icas==1 | length(f10v)==0
   ya=aou;
   yg=gou;
  if icas==1
   dya=diff(ya)./diff(Fint);
   fim=find(dya>0);
   if length(fim)==0
    ipu0=1;
    puas=[1 2];
   else
    if fim(end)+1==length(Fint)
      puas=[1:length(Fint)];
      if length(find(ya>0))==0
       ipu0=length(Fint);
      else
       ipu0=1;
      end
    else
     if fim(1)==1
      ipu0=1;
      puas=[fim fim(end)+1];
     else
      ipu0=length(Fint);
      puas=[fim(1)-1 fim];
     end
    end
   end
  else
   ipu0=1;
   puas=[1 2];
  end

   fiAz=ipu0;
   yas=ya(puas);
   ygs=yg(puas);
   xgn=Fint(puas);
   coz=polyfit(xgn,yas,1);
   ze=roots(coz);
   [mi,imi]=min(abs(ze-Fint));
   if mi>abs(diff(Fint(1:2))*2)
    ze=Fint(imi);
   end

   cog=polyfit(xgn,log10(ygs),length(xgn)-1);
   gg0=10^(polyval(cog,ze));
   ta=tetmu(nso,puas);
   cot=polyfit(xgn,ta,length(xgn)-1);
   tt0=(polyval(cot,ze));
   if gg0/vg<GMA
    nsc=nsc+1;
    ns=[ ns nso];
    ipu=[ipu ipu0];
    fso=[fso ze];
    gso=[gso gg0];
    aso=[aso 0];
    tso=[tso tt0];
   end
   sA=length(puA);
   An0=reshape(Anu(:,pou(nso,fiAz),fiAz),sA(1),1);
   Anso=[Anso An0];
   
' capeoi', keyboard   

   if icampi>=1

%    fieval
 %  disp('camdu in diss_nst 1'), keyboard
%ifp=-10
    lep=size(pou);
    Fint=Fint(1:lep(2));
    iLP=iLP1;
    fie_new
    iLP=iLPr;

%    'ferma', keyboard
    if iLP==1
     rtetm=0.5;
     nuazi=0;
     polca=0;
     polratio=0;
     mrad=0;
    end
    M2l=[M2l M2];
    rtetmv=[rtetmv rtetm];
    maziv=[maziv nuazi];
    mradv=[mradv mrad];
    polcav=[polcav polca];
    polrat=[polrat polratio];
%        ' cont diss', keyboard
 %  disp('polcav'), pausak
   end





%   Dla_new=Dlam_mod;
%   dF=diff(Dlam_mod(1:2))/(Dlam_mod(3)-1);
%   dFv=[-dF dF];
%   Dla_new(1:2)=ze*1000*lambda+dFv;


 end

end  %sau

if iUNO_long==0
%'fsovi', keyboard
 fsovi=f10v;
end

if ~exist('Fi1du')==1 
 Fi1du=Fi1;
end
for f10=fsovi
 kp=f10;

 ip=1;
 ys=[];
 fs=[];
 while kp<=min([length(Fi1du) length(ya1)])-1
% while kp<=length(ya1)-1
  de=ya1(kp+1)-ya1(kp);
  if isnan(de)==0
  if de>0
   ys(ip)=ya1(kp);
   ys(ip+1)=ya1(kp+1);
   fs(ip)=Fi1du(kp);
   fs(ip+1)=Fi1du(kp+1);
  else
   break
  end
  end  
  ip=ip+1;
  kp=kp+1;
 end
%keyboard
 kp=f10;
 if kp>length(ya1)
  kp=length(ya1);
 end 
 fd=[];
 yd=[];
 ip=1;
% 'KP', keyboard
% kp=min([length(ya1) length(Fi1du)]);
% fina=find(ya1>0 & isnan(ya1)==0);
% kp=fina(end);
 while kp>1
  de=ya1(kp)-ya1(kp-1);
%  pausak
  if de>0
   yd(ip+1)=ya1(kp-1);
   yd(ip)=ya1(kp);
   fd(ip+1)=Fi1du(kp-1);
   fd(ip)=Fi1du(kp);
  else
   break
  end
  kp=kp-1;
  ip=ip+1;
 end
 if length(fd)>0 & length(fs)>0
  Fi=[fliplr(fd) fs(2:length(fs))];
  ya=[fliplr(yd) ys(2:length(fs))];
 end
 if length(fd)>0 & length(fs)==0
  Fi=[fliplr(fd)];
  ya=[fliplr(yd)];
 end
 if length(fd)==0 & length(fs)>0
  Fi=fs;
  ya=ys;
 end

 if exist('ya')==0
  ya=[];
 end

 if length(ya)>=1
  dya=[1 diff(ya)./diff(Fi)];
  f2mv=(find(dya>0 & ya<0));
  f2pv=(find(dya>0 & ya>0));
  f2m1=length(f2mv);
  f2p1=length(f2pv);
  ie=0;
%  f1=find(dya>0 & abs(ya)<300);
  f1=find(dya>0 & abs(ya)<alMAX);
%  f1=find(dya>0 & abs(ya)<2e7);
  st=diff(Fint1(1:2));
  if f2m1*f2p1==0
   f2=find(dya>0 & abs(ya)<5);
%   f2=find(dya>0 & abs(ya)<100);
   ie=1;
  else
%   f2=find(dya>0 & abs(ya)<30);
   f2=find(dya>0);
  end
%  if icontr>=1
%   if exist('he'), close(he), clear he, end
%   he=figure; plot(Fi,ya,'*',Fi1,ya1), grid, pausak
%  end
  ya2=ya(f1);
  fa2=Fi(f1);
  dya3=[1 diff(ya2)./diff(fa2)];
  f2mv=(find(dya3>0 & ya2<0));
  f2pv=(find(dya3>0 & ya2>0));
  f2m=length(f2mv);
  f2p=length(f2pv);
  clear ya
  nosol=0;
 else
  nosol=1;
  f1=0;
  f2=0;
 end

if length(f1)>=2 & length(f2)>0 & nosol==0
%  [du,izp]=min(abs(ya2));
%  zes=fa2(izp);
  ly=length(ya2);
  fiy=find(ya2(1:ly-1).*ya2(2:ly)<0);
  puy=[fiy fiy+1];
  cozs=polyfit(fa2(puy),ya2(puy),1);
  zes=roots(cozs);
 if f2m>=2 & f2p>=2
  [du,iso]=sort(abs(ya2(f2mv)));
  yas1=ya2(f2mv(iso(1:2)));
  xas1=fa2(f2mv(iso(1:2)));
  [du,iso]=sort(abs(ya2(f2pv)));
  yas2=ya2(f2pv(iso(1:2)));
  xas2=fa2(f2pv(iso(1:2)));
  yf=[yas1 yas2];
  xf=[xas1 xas2];
  [du,isce]=sort(xf);
  xfu=xf(isce(1:3));
  yfu=yf(isce(1:3));
  coz=polyfit(xfu,yfu,length(xfu)-1);
  zep=roots(coz);
  ifz=find(imag(zep)==0);
  zep1=zep(ifz);
  if length(zep1)>1
   [du,izp1]=min(abs(zep1-zes));
   ze=zep1(izp1);
   if ze>fa2(puy(1)) & ze<fa2(puy(2))
   else
    [du,izp1]=max(abs(zep1-zes));
    ze=zep1(izp1);
   end
  else
   ze=zep1;
  end
 elseif f2m>=2 & f2p<2
  [du,iso]=sort(abs(ya2(f2mv)));
  yas1=ya2(f2mv(iso(1:2)));
  xas1=fa2(f2mv(iso(1:2)));
  yf=[yas1 ya2(f2pv)];
  xf=[xas1 fa2(f2pv)];
  coz=polyfit(xf,yf,length(xf)-1);
  zep=roots(coz);
  ifz=find(imag(zep)==0);
  zep1=zep(ifz);
  if length(zep1)>1
   [du,izp1]=min(abs(zep1-zes));
   ze=zep1(izp1);
   if ze>fa2(puy(1)) & ze<fa2(puy(2))
   else
    [du,izp1]=max(abs(zep1-zes));
    ze=zep1(izp1);
   end
  else
   ze=zep1;
  end
 elseif f2p>=2 & f2m<2
  [du,iso]=sort(abs(ya2(f2pv)));
  yas1=ya2(f2pv(iso(1:2)));
  xas1=fa2(f2pv(iso(1:2)));
  yf=[yas1 ya2(f2mv)];
  xf=[xas1 fa2(f2mv)];
  coz=polyfit(xf,yf,length(xf)-1);
  zep=roots(coz);
  ifz=find(imag(zep)==0);
  zep1=zep(ifz);
  if length(zep1)>1
   [du,izp1]=min(abs(zep1-zes));
   ze=zep1(izp1);
   if ze>fa2(puy(1)) & ze<fa2(puy(2))
   else
    [du,izp1]=max(abs(zep1-zes));
    ze=zep1(izp1);
   end
  else
   ze=zep1;
  end
 else
  [du,iso]=sort(abs(ya2));
  yas=ya2(iso);
  xas=fa2(iso);
  yf=yas(1:2);
  xf=xas(1:2);
  coz=polyfit(xf,yf,1);
  ze=roots(coz);
 end
  if icontr>=1
%   he=figure; subplot(211), plot(xf,yf,'*',Fi1,ya1,ze,0,'wo'), grid,
%   ax=axis; ax(3:4)=[-5 2]; axis(ax);
%   pausak
   if iins==4
    ichg=input(' cambio punti ? [0/1] ');
    if isempty(ichg)==1
     ichg=0;
    end
    if ichg==1
     pu=input(' punti = ');
     yf=Fi1(pu);
     xf=ya1(pu);
     coz=polyfit(xf,yf,1);
     ze=roots(coz);
     plot(xf,yf,'*',Fi1,ya1,ze,0,'wo'),
    end
   end
%   close(he), clear he
  end

if (ie==1 & min(abs(xf-ze))<st) | ie==0
%  if icontr==1
%   he=figure; plot(Fi,ya,'*',Fi1,ya1), grid, pausak
%   close(he), clear he
%  end

 [du,fd]=sort(abs(Fi-ze));
 fd0=fd(1);
 [du,fd]=sort(abs(Fint1-ze));
 ffA=(fd(1:2));
 FiA=Fint1(ffA);

 [du,fd]=sort(abs(xf-ze));
 clear ff
 for kf=1:length(xf)
  ff(kf)=find(xf(kf)==Fi1);
 end
 ff=sort(ff);
 xg=Fi1(ff);
 yg=ga(ff);
 ygm=max(yg);
 yt=ta(ff);

%disp('key')
%keyboard
% if ygm>20*min(yg) & length(yg)>2
%  f=find(yg~=ygm);
%  yg=yg(f);
%  xg=xg(f);
%  yt=yt(f);
% end

% cog=polyfit(xg,yg,length(xg)-1);
% gg0=polyval(cog,ze);
  if length(xg)>2
   xgma=max(abs(xg));
   xgn=xg/xgma;
   cog=polyfit(xgn,log10(yg),length(xg)-1);
   gg0=10^(polyval(cog,ze/xgma));
  else
   cog=polyfit(xg,yg,length(xg)-1);
   gg0=polyval(cog,ze);
  end
 cot=polyfit(xg,yt,1);
 tt0=polyval(cot,ze);
  if icontr>=1
   falam=lambda*1000;
   he=figure; subplot(211), plot(xf*falam,yf,'*',Fi1du*falam,ya1,ze*falam,0,'wo'), grid,
   figure(he); subplot(212)
   semilogy(xg*falam,yg,ze*falam,gg0,'wo'), grid,
   title(' control 35')
   pausak
%   close(he), clear he
  end
% disp(' verifica '), pausak

% pausak


 if gg0/vg<GMA
  nsc=nsc+1;
  ns=[ ns nso*fd0./fd0];
%  ipu=[ipu fi(fd0)];
  ipu=[ipu fd0];
  fso=[fso ze];
  gso=[gso gg0];
%  disp('gso'), pausak
  aso=[aso 0];
  tso=[tso tt0];
%  fip=find(pou(nso,:)~=0);
%  for ifip=1:length(fip)
%   V(:,ifip)=(Anu(:,pou(nso,fip(ifip)),fip(ifip)));
%  end
%  sA=size(Amu);
  if icontr>=1
%   he=figure; subplot(211), plot(xf,yf,'*',Fi1,ya1,ze,0,'wo'), grid,
%   ax=axis; ax(3:4)=[-5 2]; axis(ax);
%   pausak
%   close(he), clear he
  end
  lF=length(Fint1);

  fiAz=find(FiA(1)==Fint);
  if exist('Kvet')
   KKdu=Kve(:,fiAz);
   if KKdu(1)<1e-10
    KKdu(1)=1e-10;
   end
   fiK=find(KKdu~=0);
   KK=KKdu(fiK);
   npk=length(KK);
  end
  pbA=[1:dnum*length(KK)];
  sA=length(puA);

%  An0=reshape(Anu(:,pou(nso,fiAz),fiAz),sA(1),1);
clear As
lep=size(pou);

% for iaf=1:length(Fint)
 for iaf=1:lep(2)
  Pun=pou(nso,iaf);
  if Pun>0
   Ada=reshape(Anu(:,Pun,iaf),sA(1),1);
  else
   Ada=zeros(sA(1),1);
  end
%  [du,ima]=max(abs(Ada));
%  As(:,iaf)=Ada*sign(real(Ada(ima)))/median(abs(Ada));
  As(:,iaf)=Ada/median(abs(Ada));
 end
 
 Apre=As(:,1);
 corp=Apre'*Apre;
for ksor=2:lep(2)
  Apo=As(:,ksor);
  corpi=Apre'*Apo;
  Apre=Apo;
  corp=corpi;
  if real(corpi)<0
   As(:,ksor)=-As(:,ksor);
   Apre=-Apo;
   corp=-corpi; 
  end
 end
 
  duded=((Fint-ze));
  fimaz=find(duded>0);
  isoz=fimaz(1)+[-1 0];
  DeFi=(diff(Fint(isoz)));
  dude=abs(Fint(isoz)-ze);
%  An0=(dude(1)/DeFi)*As(:,isoz(1))+(1-dude(2)/DeFi)*As(:,isoz(2));
  An0=(dude(2)/DeFi)*As(:,isoz(1))+(dude(1)/DeFi)*As(:,isoz(2));
 if ifp==-1
  KKd=1:length(An0);
  figure, semilogy(KKd,abs(As),KKd,abs(An0),'w.-'), pausak
 end
  if icam_fr==1
   KKd=1:length(An0);
   haut=figure;
%   plot(KKd,abs(As),KKd,abs(An0),'w.-'), pausak
%   ' An0 disply '
   fr_pun =1:length(Fint);
   field_fr
   keyboard
  end

  if icontr>=1
   figure, semilogy(abs(An0)/median(abs(An0)),'ro'), hold on,
   semilogy(Adispla),
   semilogy(abs(As(:,isoz(1:2))),'.')
%   figure, semilogy(abs(An0),'ro'), hold on,
%   semilogy(Adispla),   
   expla=' Coeff to be interpolated,  Red dots: solution ';
   title(expla)
   
   
   pausak
  end
%  'prima Anso'
 % keyboard
fiAv=length(find(isnan(An0))); 

if fiAv==0
  Anso=[Anso An0];
else 
 gso=gso(1:end-1);
 fso=fso(1:end-1);
end  
%  'dopo Anso'
%  keyboard

  if ifp~=-4
  disp('campi in diss_new')
  end
%' capeoi questa', keyboard   

  if icampi>=1 & fiAv==0
%%%%%%%%% routine plot campi
%disp(' prima di fieval in dissus '), keyboard

%   fieval
%   disp('camdu in diss_nst 2'), keyboard
   fa1d=fou(nso,ps);
   fa1d=fa1d(fi);
%   figure, plot(fa1d,ya1), pausak, plot(fa1d,ga), pausak
%   disp('fie_new in diss_nst 1'), keyboard
%ifp=-10
    Fint=Fint(1:lep(2));
    iLP=iLP1;
    fie_new
    iLP=iLPr;   

%    'ferma', keyboard
   if iLP==1
    rtetm=0.5;
    nuazi=0;
    polca=0;
    polratio=0;
    mrad=0;
   end
   M2l=[M2l M2];
   rtetmv=[rtetmv rtetm];
   maziv=[maziv nuazi];
   mradv=[mradv mrad];
   polcav=[polcav polca];
   polrat=[polrat polratio];
%  disp('polcav'), keyboard
  end


 end %if gg0
end %ie
end %length(f1)
end %for f10

if iins>1 & nsc==0
 fst_d=5;
else
 fst_d=1;
end

if abs(ifnz)==1 & length(fsovi)==0 & iins>=1

  ifi=length(find(ya1<0));
  if length(ya1)==ifi
   iso=[ifi-1 ifi];
   isop=iso;
  elseif ifi==0
   iso=[1 2];
   isop=iso;
  else
   [du,fid]=sort(abs(ya1));
   isop=sort(fid(1:2));
%   isop=fid;
   if isop(1)<length(isop)
    iso=[isop(1) isop(1)+1];
   else
    iso=[isop(1)-1 isop(1)];
   end
  end
  isoA=isop;
  yt=ta(iso);
  yg=ga(iso);
  yf=ya1(iso);
  xf=Fi1(iso);
  xg=xf;
  coz=polyfit(xf,yf,npma-1);
  zevt=sort(roots(coz));
%  pausak
  iz=find(imag(zevt)==0);
  if length(iz)==0
   isacc=0
  elseif length(iz)==1
   ze=zevt;
   if min(abs(ze-xf))<2*diff(sort(xf))*fst_d
    isacc=1;
   else
    isacc=0;
   end
  elseif length(iz)==2
   zev=sort(zevt(iz));
   xfs=(sort(xf));
   stv=diff(sort(xf));
   st2=stv(1);
   st1=stv(1);
   if zev(2)<Fint(1)
    zev=zev(2);
    st1=stv(1)*5;
   elseif zev(2)>Fint(length(Fint))
    zev=zev(1);
    if iins>1
     st1=stv(1)*fst_d;
    end
   end
   izea=find(zev>xfs(1)-st1 & zev<xfs(length(xfs))+st2);
   zevs=zev(izea);
   isacc=1;
   if length(zevs)>1
    if ifnz==-1
     ze=zevs(1);
    else
     ze=zevs(2);
    end
   elseif length(zevs)==1
     ze=zevs;
   elseif length(zevs)==0
    isacc=0;
   end
  end
 if nsc==0
  isacc=1;
 end
 if iins>2
  if icontr>=1
   he=figure; plot(xf,yf,'*',Fi1,ya1,ze,0,'wo'), grid,
%   ax=axis; ax(3:4)=[-5 2]; axis(ax);
   pausak
   isacc=input(' accetti soluzione ? [0/1]');
   if isempty(isacc)==1
    isacc=1;
   end
   close(he), clear he
  end
 end

 if isacc==1

  if iins<=2
   if icontr>=1
    he=figure; plot(xf,yf,'*',Fi1,ya1,ze,0,'wo'), grid,
    ax=axis; ax(3:4)=[-5 2]; axis(ax);
    pausak
    close(he), clear he
   end
  end
  ygm=max(yg);
  if ygm>10*min(yg) & length(yg>2)
   f=find(yg~=ygm);
   yg=yg(f);
   xg=xg(f);
   yt=yt(f);
  end
  if length(xg)>2
   cog=polyfit(xg,log10(yg),length(xg)-1);
   gg0=10^(polyval(cog,ze));
  else
   cog=polyfit(xg,yg,length(xg)-1);
   gg0=polyval(cog,ze);
  end
  cot=polyfit(xg,yt,1);
  tt0=polyval(cot,ze);
  if icontr>=1
   he=figure; semilogy(xg,yg,ze,gg0,'wo'), pausak
   close(he), clear he
  end
   if gg0/vg<GMA
    nsc=nsc+1;
    ns=[ ns nso];
    [du,fiz]=min(abs(Fint-ze));
    ipu=[ipu fiz];
    fso=[fso ze];
    gso=[gso gg0];
 %  disp('gso'), pausak
    aso=[aso 0];
    tso=[tso tt0];
 %   sA=size(Amu);
    sA=length(puA);

 %   An1=reshape(Amu(:,nso,fisa(isoA)),sA(1),length(isoA));
 %   An0=An0d(:,1);
 %   if icontr==1
 %    figure, plot(abs(An0d)), hold on, plot(abs(An0),'w*'), pausak
 %   end
   fib=find(pou(nso,:)~=0);
   Fz=Fint./Fint*1e10;
   Fz(fib)=Fint(fib);
   [du,fiAz]=min(abs(ze-Fz));

'alsp metodo', keyboard
    An0=reshape(Anu(:,pou(nso,fiAz),fiAz),sA(1),1);
    Anso=[Anso An0];


   disp('campi in diss_new altro')

 %  imod=1;
 %  disp('camdu in dissu3')
 %  keyboard
   if icampi>=1

%    fieval
%   disp('camdu in diss_nst 3'), keyboard
%ifp=-10
    Fint=Fint(1:lep(2));
    iLP=iLP1;
    fie_new
    iLP=iLPr;    
    'ferma', keyboard

    if iLP==1
     rtetm=0.5;
     nuazi=0;
     polca=0;
     polratio=0;
     mrad=0;
    end

    rtetmv=[rtetmv rtetm];
    maziv=[maziv nuazi];
    mradv=[mradv mrad];
    polcav=[polcav polca];
    polrat=[polrat polratio];
    M2l=[M2l M2];    
%    ' cont diss', keyboard
 %  disp('polcav'), pausak
   end

%%%%%%%%%%%%%%%%%%%%%%%%%

   if exist('Kvet')
    KKdu=Kve(:,fiAz);
    KKdu(1)=1e-10;
    fiK=find(KKdu~=0);
    KK=KKdu(fiK);
    npk=length(KK);
   end
%   pbA=[1:dnum*length(KK)];
%   An0c=An0(pbA);

   if icontr>=1
    figure, plot(abs(An0),'w.'),
    pausak
   end
  end  %gg0
 end % isacc
end  % if ifnz

end % nso

%'qui cont', keyboard
avv=aou;
Gvv=gou;
F=fou;
if exist('isofg')==0
 isofg=1;
end
gsoM=gso.*M2l;
gsoM=M2l;
' guadagno pesato con M2'
%gsoM=gso;
if isofg==1
[du,iso]=sort(fso);
else
[du,iso]=sort(gsoM);
end
%iso=iso(1:nmasce);

%'qui nmaze', keyboard
gso=gso(iso);
if length(gso)>nmasce
 gso(nmasce+1:end)=-gso(nmasce+1:end);
end 

if exist('Gsov')
%'Gsov', keyboard
%if length(Gsov)>length(gso)
for ksov=1:length(Gsov)
 if length(find(Gsov(ksov)/2==-gso))==1 
  Gsov(ksov)=-1;
 end
end
end
%end
fso=fso(iso);
mso=M2l(iso);

%fso=fso(1:end-1);
%gso=gso(1:end-1);

%fou=fou(iso,:);
%aou=aou(iso,:);
%gou=gou(iso,:);

ns=ns(iso);
ipu=ipu(iso);
aso=aso(iso);


if icampi==0 | length(icampi)==0
 polcav=zeros(size(gso));
 polrat=polcav;
 rtetmv=polcav;
 maziv=polcav;
 mradv=polcav;
end
polcav=polcav(iso);
polrat=polrat(iso);
rtetmv=rtetmv(iso);
maziv=maziv(iso);
mradv=mradv(iso);
%gamso=gamso(iso);
%gamoso=gamoso(iso);
tso=tso(iso);
Am=Anso(:,iso);

nsv(1:length(fso),its)=ns';
ipuv(1:length(fso),its)=ipu';
fsov(1:length(fso),its)=fso';
gsov(1:length(fso),its)=gso';
rtetmsov(1:length(fso),its)=rtetmv';
mazisov(1:length(fso),its)=maziv';
mradsov(1:length(fso),its)=mradv';
polasov(1:length(fso),its)=polcav';
polratv(1:length(fso),its)=polrat';
%gamsov(1:length(fso),its)=gamso';
%gamosov(1:length(fso),its)=gamoso';
asov(1:length(fso),its)=aso';
tsov(1:length(fso),its)=tso';
sm=size(Am);
Amv(its,1:sm(1),1:sm(2))=Am;

 %pv0=avv.*Gvv/vg;
 %z0=aso.*gso/vg;
 pv0=avv;
 fie1=find(tso>=.5);
 z0e=aso(fie1);
 fim1=find(tso<.5);
 z0m=aso(fim1);
 sG=size(Gvv);
end

if iLP==1
 Fsov=fsov;
 Gsov=gsov;
 Gsosa=gsov;
end
La=lambda*1000./(1+Fsov);
RESU=[Gsov'/vgconv La'];

 


if isadiss==1
 save diss
end 

if length(gso)==0
%' fine diss_nst', keyboard
 'fine diss1, no modi', keyboard
end