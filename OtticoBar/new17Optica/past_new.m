function [fatqw,nref,aref,Lr,Lpl,nr,npl,xr,xpl,nmir,ar,shav,frp,lambda,Fi,any,...
          fmu,z,uL,fiQ,fiC,Lpla,npla,Nspla,aral,...
          nv0,dv,xm,radii,fst,perd,anyf,iauto,shavet,ipar,ifield,istfi,...
          icaplo,icsfit,icsfib,icsfi,istmet,dglo,iFFsect,L_i,n_i,confzv,lambda0]=...
          past_new(nomeF,ifp,iany,nomeFe,Bi,iloa,ilo,imet,saELM);

fatqw=[];nref=[]; aref=[]; Lr=[]; 
Lpl=[]; nr=[]; npl=[]; xr=[]; xpl=[]; nmir=[]; ar=[]; shav=[]; 
frp=[]; lambda=[]; Fi=[]; any=[]; 
fmu=[]; z=[]; uL=[]; fiQ=[]; fiC=[]; Lpla=[]; npla=[]; Nspla=[]; aral=[]; 
nv0=[]; dv=[]; xm=[]; radii=[]; fst=[]; perd=[]; anyf=[]; iauto=[];
icaplo=[]; icsfit=[]; icsfib=[]; icsfi=[]; istmet=[]; dglo=[]; iFFsect=[]; 
L_i=[]; n_i=[]; confzv=[]; lambda0=[];  shavet=[]; ipar=[]; ifield=[]; istfi=[]; 

iNUO_bir=1; 
% ' dentro past_new ', keyboard
% ' dentro past_new ', keyboard
% ' dentro past_new ', keyboard
if ilo==1
 eval(['load ' saELM]);
else

icalcFi=1;
% icalcFi=0;                                             %
%disp(' ')                                               %
%disp(' ')                                               %
%disp(' ')                                               %
%disp(' !!!!!!!!!!!!!!!!!! attenzione,    attenzione ')  %
%disp(' ')                                               %
%disp('          non si calcola il campo ')              %

% read structure data

%[nref,aref,nv0,dv,xm,Dov,radii,fst,perd,anyf,iauto,dvw,xmw,fsw,dov,shavet,...
%   ipar,ifield,lambda0,igrad_tr]=/igrad_trLay_new(nomeF);
%save lay
%'prima di layt', keyboard
%'prima di layt', keyboard
%'prima di layt', keyboard

[nref,aref,nv0,dv,xm,Dov,radii,fst,perd,anyf,iauto,dvw,xmw,fsw,dov,shavet,...
   ipar,ifield,lambda0,igrad_tr]=Lay_tran(nomeF);

global FCA
FCA.on=0;
FCA.e=3;
FCA.h=7;
if FCA.on==1
    a_find=find(Dov>0);
    d_find=find(Dov<0);
    
    alpha=zeros(size(Dov));
    alpha(a_find)=FCA.h.*abs(Dov(a_find));
    alpha(d_find)=FCA.e.*abs(Dov(d_find));
    alpha=FCA.f*alpha/1e18;
    kNoll=2*pi/lambda0;
    
%     fullVectorConversion

    perd_new=alpha*100/(2*kNoll);
    nv0_new=nv0(:,1)-1i*alpha*100/(2*kNoll);
%     nv0_new(1:4)=nv0(1:4,1);
    
    % Overwrite hard coded vector
%     nv0(:,1)=nv0_new;
%     perd(5:end,1)=alpha(5:end);
    
    nv0(:,1)=nv0_new;  
    perd(:,1)=perd(:,1)+alpha;
    
    figure
    hold on
    plot(cumsum(dv),perd(:,1))
    plot(cumsum(dv),alpha)
    ylabel('log(FCA)')
    set(gca,'yscale','log')
    set(gca,'FontName','Times New Roman','Fontsize',13)
    keyboard
end


%disp(' contr dopo lay_new '), keyboard

if ifp~=-4
%disp(' contr dopo lay_new '), pausak
end


if length(find(ifield==-2))>1
 disp(' past_new.m: error: only one print section ')
 keyboard
end
%'dopo lay', keyboard


for kstra=1:length(nv0)
 idu=find(xm(kstra,:)~=0);
 [du,ip]=max(xm(kstra,idu));
 n1(kstra,1)=nv0(kstra,idu(ip));
 c1(kstra,1)=xm(kstra,idu(ip));
end
%'fst', keyboard
fima=max(fst(:,2));
fim0=find(fst(:,2)==fima);
du=find(diff(fim0)>1);
if length(du)>0
 fim=fim0(1:du(1));
else
 fim=fim0;
end
lastim=dv(fim)'*real(n1(fim))*2e-9;
if abs(1-lastim/lambda0)>.01
% disp(' lam_stim =  '), lastim
% disp(' restart with better lam0 guess ')
% disp(' restart with better lam0 guess ')
% keyboard
% disp(' restart with better lam0 guess ')
else
 %lambda0=lastim;
end

fiff=find(abs(ifield)>=20);
if length(fiff)>0
 if fiff<3
  iFFsect=1;
 else
  iFFsect=3;
 end
else
 iFFsect=1;
end
%' iFFsect', keyboard

T=300;
% trovo lo strato del contatto
sn=size(nv0);
if sn(2)>1
icont=find(nv0(:,2)==1);  % struttura con rilievo
if length(icont)==0
 icont=find(imag(nv0(:,2))<.1); % struttura senza rilievo
end
if length(icont)==0
 icont=find(nv0(:,1)==1); % struttura senza contatto incluso nel file
end
end

global ivfre


 nv=nv0;
 si=size(nv0);
 for ksi=1:si(1)
  xal=xm(ksi,:);
  imatw=xal;
  fial=find(xal>0);
  if length(fial)>0
   imatw(fial)=1;
   fisos=find(imatw~=0 & imatw~=-10);
   if ivfre==1
    nvv=nv0(ksi,fisos);
    nal=ref_ind(nvv,imatw(fisos),lambda0,xal(fisos));
    nv(ksi,fisos)=nal;
    nv0(ksi,fisos)=nal;
   end
  end  %length
 end

% nv=nv0;
% si=size(nv0);
% for ksi=si(2):-1:1
%  xal=xm(:,ksi);
%  imatw=xal;
%  fial=find(xal>0);
%   imatw(fial)=1;
%   fisos=find(imatw~=0);
%   if ivfre==1
%    nal=ref_ind(imatw,lambda0,xal);
%    nv(fisos,ksi)=nal(fisos);
%    nv0=nv;
%   end
% end

Dop=Dov(fial,1);

%yal=[];
%[mue,muh,me,mh,tc]=momatc(imatw,T,xal,yal);
%
%dglo.mue=mue;
%dglo.muh=muh;
%dglo.me=me;
%dglo.mh=mh;
%dglo.tc=tc;
dglo.Dop=Dov;


%T=300;
%x=linspace(0,1,101)';
%[mue,muh,me,mh,tc]=momatc(x,T);
%figure
%plot(x,mue,x,muh), title(' mobility')
%
%figure
%plot(x,me,x,mh), title(' masses')
%
%figure
%plot(x,tc,'r'), title(' ther. cond ')
%disp(' momatc in past_new ')
%keyboard

if length(find(ifield~=0))==0
 ifi=find(shavet~=0);
 ifield(ifi(1))=-2;
end

fatqw=0;
iput=find(iauto(:,1)==1);
ipua=find(iauto(:,1)==2);
ipub=find(iauto(:,1)==3);

ipudd=find(iauto(:,2)==-4 | iauto(:,2)>0);
%fiC=[ipudd(1):ipudd(length(ipudd))]-1;
%fiQ=find(fst(:,2)==-1)-1;
fi1=find(iauto(:,1)==1)-1;
fiCdu=find(iauto(:,2)==-4)-fi1;
fiC=fiCdu(1):fiCdu(end);
%fiQ=find(fst(:,2)==-1);
%fim4=find(iauto(:,2)==-4)-fi1;
%ficavit=fim4-fi1;
%pcavit=ficavit(1):ficavit(2);
%'fiC',keyboard

put=iput+1:ipua-1;
putm=iput:fiC(1)-1;

pub=ipua+1:ipub-1;
pubm=fiC(end)+1:ipub-1;

pua=ipua;
if iany>=2
 anyv=anyf;
else 
 anyv=zeros(size(dv));
end 


 any.t=anyv(put,:);
 any.b=anyv(pub,:);
 any.a=anyv(pua,:);
 global EO
%EO.Dov=Dov; 
% 'prima di return any', keyboard
% 'prima di return', keyboard
if iNUO_bir==1
 return
end 

 fip=find(nomeF=='.');
 nomeP=[nomeF(1:fip) 'par'];
 nomePsa=[nomeF(1:fip-1),'par'];
 
    DirNam=[cd, '/D1ANA/GenerateGeom/'];
   DR=dir(DirNam);
   lD=length(DR);
   fl_find=0;
   for kD=3:lD
    dun=getfield(DR(kD),'name');
     if strcmpi(dun,nomeP)==1
      fl_find=1;
      break
     end
   end
   if fl_find==1
    if ifp~=-4
%     isovr=input(' Vuoi sovrascrivere file ? [1]');
    isovr=1;
    else
    isovr=[];
    end
    if length(isovr)==0
     isovr=0;
    end
    if isovr~=1
     nomeP='fileDUM.par';
    end
   end

 fidfile = fopen(nomeP,'w');

if length(ipar)~=1
  sip=size(ipar);
  ip1=reshape(ipar(:,1,:),sip(1)*sip(3),1);
  if sip(2)>=2
   ip2=reshape(ipar(:,2,:),sip(1)*sip(3),1);
  else
   ip2=ip1;
  end
  if sip(2)>=3
   ip3=reshape(ipar(:,3,:),sip(1)*sip(3),1);
  else
   ip3=ip1;
  end

  riga=[1:sip(1)]';
  righe=riga*ones(1,sip(3));
  ip4=reshape(righe,sip(1)*sip(3),1);
  fi1=find(ip1~=0);
  [dup,fi1p]=sort(ip1(fi1));
  fi1=fi1(fi1p);

  fid=find(diff([0; dup])~=0);
  fiud=fi1(fid);
  [ip11,fiudd]=sort(ip1(fiud));
  fiu=fiud(fiudd);

  npar=length(find(diff([0; dup])~=0));

if ifp~=-4
 disp(' ')
 disp(' ')
 disp(' ')
 disp(' ')

 screen('         # of parameters = ',npar)
 disp(' ')
 disp(' ')
 disp(' ')
end
 type(1,:)='thickness    ';
 type(2,:)='y-axis       ';
 type(3,:)='x-axis       ';
 type(4,:)='shape param  ';
 type(5,:)='ref. index   ';
 type(6,:)='array centers';
 type(7,:)='shape        ';
 type(8,:)='pair number  ';
 type(9,:)='rad. cir. gr.';
 type(10,:)='Absorption   ';

 typem(1,:)='thi';
 typem(2,:)='yax';
 typem(3,:)='xax';
 typem(4,:)='spa';
 typem(5,:)='ind';
 typem(6,:)='cen';
 typem(7,:)='sha';
 typem(8,:)='#pa';
 typem(9,:)='rcg';
 typem(10,:)='Abs';

 typeg(2,:)='grat. orien';
 typeg(3,:)='grat. shift';
 typeg(4,:)='cent. r_cir';
 typeg(5,:)='grat. peach';
 typeg(6,:)='grat. width';
 typeg(7,:)='grat. est.y';
 typeg(8,:)='grat. est.x';
 typeg(9,:)='grat. est.s';
typeg(11,:)='grat.  Rext';

 typea(2,:)='number x pixels';
 typea(3,:)='number y pixels';
 typea(4,:)='pixel shape    ';
 typea(5,:)='pixel Ry       ';
 typea(6,:)='pixel Rx       ';
 typea(7,:)='pixel Delta    ';
 typea(8,:)='pixel spacing  ';

 typegm(2,:)='ori';
 typegm(3,:)='shi';
 typegm(4,:)='ccr';
 typegm(5,:)='per';
 typegm(6,:)='wid';
 typegm(7,:)='esy';
 typegm(8,:)='esx';
 typegm(9,:)='esh';
 typegm(11,:)='Rex';

 typeam(2,:)='#px';
 typeam(3,:)='#py';
 typeam(4,:)='sha';
 typeam(5,:)='pRy';
 typeam(6,:)='pRx';
 typeam(7,:)='pDe';
 typeam(8,:)='spa';

 shad=shavet;

 fis=find(shavet==0);
 shad(fis)=20;
 shad=[shad ones(size(shad(:,1)))*20];

 types(1,:)='  circle  ';
 types(2,:)=' rectangle';
 types(3,:)='  ellipse ';
 types(4,:)='  rhombus ';
 types(5,:)='   array  ';
 types(6,:)='  grating ';
 types(7,:)='  lens    ';
 types(8,:)='  DOE     ';
 types(9,:)=  '  particle';
 types(11,:)='   vortex ';
 types(20,:)='  planar  ';

 typesm(1,:)='ci';
 typesm(2,:)='re';
 typesm(3,:)='el';
 typesm(4,:)='rh';
 typesm(5,:)='ar';
 typesm(6,:)='gr';
 typesm(7,:)='le';
 typesm(8,:)='DO';
 typesm(9,:)='pa';
 typesm(11,:)='vx';
 typesm(20,:)='pl';




 nul='';
if ifp~=-4
if npar>0

 disp(' ')
 disp(' ')
 eval(['fprintf(fidfile,''' nul '\r\n'');']);
 prp=(['      Parameters to be supplied in par_in ']);
 disp(prp);
 eval(['fprintf(fidfile,''' prp '\r\n'');']);
 eval(['fprintf(fidfile,''' nul '\r\n'');']);
 disp(' ')
 prp=(['      Par,  Layer,   Tran_sect,   Par. type,       shape,       short_name']);
 eval(['fprintf(fidfile,''' prp '\r\n'');']);
 eval(['fprintf(fidfile,''' nul '\r\n'');']);
 disp(prp);
 disp(' ')
end
end
 tabm='      ';
 tab='       ';
 tab2='  ';
 tab1='   ';
 tab0='';
% fiu=fi1;
 fiu_sav=fiu;
 [du,puparam]=sort(ip4(fiu));
  icsa=1;
% for k=1:length(fiu)

 for k=puparam'
 TYPE=[];
 TYPEm=[];
% if  k==16,
%  ' stop ', keyboard
% end
% if length( find( ip2(fiu(k))<0 & ip2(fiu(k))>-11))>0
 if length( find( ip2(fiu(k))<0))>0
  fis=find(shad(ip4(fiu(k)),:)==6);
  if length(fis)>0
%  'quo', keyboard
%   TYPE=typeg(abs(ip2(fiu(k))),:);
%   TYPEm=typegm(abs(ip2(fiu(k))),:);
   if abs(ip2(fiu(k)))<=11
   TYPE=typeg(abs(ip2(fiu(k))),:);
   TYPEm=typegm(abs(ip2(fiu(k))),:);
   end   
  end
  fis=find(shad(ip4(fiu(k)),:)==5);
  if length(fis)>0
   TYPE=typea(abs(ip2(fiu(k))),:);
   TYPEm=typeam(abs(ip2(fiu(k))),:);
%   TYPEm=TYPE;
  end
 else
  TYPE=type(ip2(fiu(k)),:);
  TYPEm=typem(ip2(fiu(k)),:);
 end
 %type(ip2(fiu(k)),:)
 if ip1(fiu(k))>9
  ta1=tabm;
 else
  ta1=tab;
 end
 if ip4(fiu(k))>9
  ta2=tabm;
 else
  ta2=tab;
 end

%' qui ', keyboard

      if ip4(fiu(k))<=9
       ch1=['0',num2str(ip4(fiu(k)))];
      elseif ip4(fiu(k))>=100
       ch1=[num2str(ip4(fiu(k))-50)];
      else
       ch1=[num2str(ip4(fiu(k)))];
      end
      ch2=num2str(ip3(fiu(k)));
      shada=abs(shad(ip4(fiu(k)),ip3(fiu(k))));
  sac=[TYPEm,ch2,typesm(shada,:),ch1];
%  pausak
  du=sac(1:min([10 end]));
  sacm(icsa,1:length(du))=du;
  sapm(icsa,:)=k;
  icsa=icsa+1;

  chio=[tab1,'@ '];
  prp=([ta1 num2str(ip1(fiu(k))) ta2 tab1  num2str(ip4(fiu(k))) tab ...
       num2str(ip3(fiu(k)))  tab TYPE ...
       tab2 types(shada,:),chio,sac]);
 if ifp~=-4
  disp(prp);
 end
  eval(['fprintf(fidfile,''' prp '\r\n'');']);
 end
 fiu=fi1;

else
 prp=(['     NO Parameters to be supplied in par_in ']);
 if ifp~=-4
 disp(' ')
 disp(' ')
  disp(prp);
 end
  eval(['fprintf(fidfile,''' prp '\r\n'');']);
 sacm=[];
 sapm=[];

end

 eval(['fprintf(fidfile,''' nul '\r\n'');']);
 fclose(fidfile);


% eval(['save ',nomePsa,' sacm sapm']);

 if ifp~=-4
 disp(' ')
 disp(' ')
  disp(' ')
  disp(' ')
  end



%fis=find(shavet~=0);
%shd=shavet(fis);
%if length(find(diff(shd)~=0))>0
% istrumix=1;
%else
% istrumix=0;
%end
%istrumix=1;

%' prima di reassign ', keyboard

reassign
fiQ=find(fst(:,2)==-1)-1;

%'fiQ ', keyboard


xpl.b.last=xm(pup.b.last,1);
xpl.b.o=xm(pup.b.o,1);
xpl.b.m=xm(pup.b.m,1);
xpl.b.i=xm(pup.b.i,1);

xpl.t.o=xm(pup.t.o,1);
xpl.t.m=xm(pup.t.m,1);
xpl.t.i=xm(pup.t.i,1);


for kstra=1:length(nv0)
% idu=find(xm(kstra,:)~=0);
 idu=1;
 [du,ip]=max(xm(kstra,idu));
 n1(kstra,1)=nv0(kstra,idu(ip));
 c1(kstra,1)=xm(kstra,idu(ip));
end
%'dopo'
%keyboard

%n1=(nv0(:,1));
%c1=xm(:,1);

fima=find(dv>1000);
if length(fima)>0
% 'attenzione riduco dv per campo numerico' , keyboard
% dv(fima)=1000;
end

xd=[zeros(size(dv)) dv];
fl=abs(fst(:,2));
fim=find(fl>1);
%fim=find(fl>=1);
if length(fim)>0
 if length(find(diff(fim)>1))==0 & sum(diff(fim))>length(fim)-1
  fim=[fim; length(dv)-1];
 end
else
 fim=length(dv);
end
fw=length(fim);
%fw=length(fl)-1;

%'fim', keyboard

miup=fst(fim(1),2);
midw=fst(fim(end),2);
poup=miup-4;
podw=4;


dto=[];
nto=[];
ato=[];
cto=[];
puo=[];

%fista=find(abs(Dov)>0);
%piv=fista(1);
piv=1;
ico=1;
%'fim ', keyboard
fsto=[];
while fw>0
 pd=fim(ico);
 nl=fst(pd,1);
 nlay=abs(fst(pd,2));
% [pd piv nl nlay], pausak
 if piv-pd~=0
  pup0=piv:pd-1;
  ddu=dv(pup0);
  imz=find(ddu>0);
  pup0=pup0(imz);
  if length(find(pup0==pua))==1
   fi0=find(pup0==pua);
   leqw=length(fsto)+fi0;
  end
  ay=anyf(pup0);
  nto=[nto; n1(pup0)];
  cto=[cto; c1(pup0)];
  fsto=[fsto; fst(pup0,2)];
  dto=[dto; dv(pup0)];
  ato=[ato; dv(pup0).*ay];
  puo=[puo; pup0'.*ay];
 end
%'ICI 0 bir', keyboard

 pupe=pd:pd+nl-1;
%'ICI 0 bir', keyboard
%if nlay>1
% [pd piv nl nlay fw], pausak
%end 

 ay=anyf(pupe);
 if length(pupe)==0
  nlay=0;
 end
 for inp=1:nlay
  nto=[nto; n1(pupe)];
  cto=[cto; c1(pupe)];
  fsto=[fsto; fst(pupe,2)];
  dto=[dto; dv(pupe)];
  if nlay==miup | nlay==midw
   if (inp==poup & nlay==miup) | (inp==podw & nlay==midw)
    ato=[ato; dv(pupe).*ay];
    puo=[puo; pupe'.*ay];
   else
    ato=[ato; -dv(pupe).*ay];
    puo=[puo; pupe'*0];
   end
  else
   if inp==fix(nlay/2)
    ato=[ato; dv(pupe).*ay];
    puo=[puo; pupe'.*ay];
   else
    ato=[ato; -dv(pupe).*ay];
    puo=[puo; pupe'*0];
   end
  end
 end
 pupe
 
% 'pupe', pausak
 if length(pupe)>0
 piv=pupe(length(pupe))+1;
 fw=fw-nl;
 ico=ico+nl;
 else
  fw=0;
  piv=length(dv);
 end
%  'pupe dopo', pausak
end

%'ICI 1 bir', keyboard

if piv<length(dv)+1
 pup0=piv:length(dv);
 nto=[nto; n1(pup0)];
 cto=[cto; c1(pup0)];
 fsto=[fsto; fst(pup0,2)];
 dto=[dto; dv(pup0)];
 ay=anyf(pup0);
 ato=[ato; dv(pup0).*ay];
 puo=[puo; pup0'.*ay];
end

puac=puo;
puae=find(ato>0);
%puae=find(ato>=0);
%' puae ', keyboard


nd=[nto nto];
ccd=[cto cto];
%pudto=
xd=[zeros(size(dto)) dto];
xt=cumsum(reshape(xd',2*length(dto),1));
nt=real(reshape(nd',2*length(dto),1));
ct=(reshape(ccd',2*length(dto),1));
L_i=dto;
n_i=nd(:,1);
nto(leqw)=real(nto(leqw));

lmaxim=1e6;
%lmaxim=1e3;
fitro=find(dto>lmaxim);
if length(fitro)>0
 dto(fitro)=lmaxim;
end
dJ=cumsum([0; dto]);
fJ=([0; fsto]);
nJ=([nto(1); nto]);
xJ=([cto(1); cto]);
lJ=length(dJ);

%' hz set', keyboard
sthz=1;
if length(find(dto>2e4)>0)
sthz=3;
end

hz=[dJ(1):sthz:dJ(lJ)];
uFunc=hz*0;
uF0=hz*0;
perm=[];
xfit=[];
for il=1:lJ-1
 fiz=find(hz>=dJ(il) & hz<dJ(il+1));
 if length(fiz)>0
  perm=[perm nJ(il+1)*ones(size(fiz))];
  xfit=[xfit xJ(il+1)*ones(size(fiz))];
  if fJ(il+1)==-1
   uFunc([fiz(1)-1 fiz])=1;
  end
  if fJ(il+1)==-1 & il==leqw
   uF0([fiz(1)-1 fiz])=1;
  end
 end  %if 
end
xfit(length(hz))=xfit(length(hz)-1);
perm(length(hz))=perm(length(hz)-1);
uFunc(length(hz))=0;
uF0(length(hz))=0;

%' primapo fipe', keyboard
%' primapo fipe', keyboard

fipe=find(real(perm)<0);
perm(fipe)=-perm(fipe);

relPerm=conj(perm.^2);
%relPerm=real(relPerm);

xdu=[dJ(1:end-1) dJ(2:end)];
ydu=[nJ(2:end) nJ(2:end)];

x=reshape(xdu',length(xdu)*2,1);
y=reshape(ydu',length(xdu)*2,1);

fipe=find(real(y)<0);
y(fipe)=-y(fipe);

%'pastnew ICI 3 bir', keyboard
%lambda0=974.5/1e9

if icalcFi==1
 if ifp~=-4
 disp(' prima di risonanza ')
 figure, plot(hz,real(sqrt(relPerm)),x,y)
 keyboard
 end
% ' mio prima', keyboard
% lambda0p=lambda0; uFuncp=uFunc; uF0p=uF0; relPermp=relPerm;
% save pri lambda0p uFuncp uF0p relPermp

 [Ksi,lambdas,Fi,uLong,uLong0,Fa]=eiglmio(lambda0,uFunc,uF0,relPerm);
 uL(1)=uLong0;
 uL(2)=uLong;
 if exist('rr')
  rg=rr;
 else
  fia=find(iauto(:,1)==2);
  rg=nv0(fia,1);
 end 
 fatqw=uLong/(uLong0*nmir.a);
 lambda=lambdas*1e6;
 gpla=2e4*pi*rg/lambda*imag(Ksi)/uLong0/fatqw;
 gpQW=gpla/nmir.a;


 if ifp~=-4
 disp(' fine Ksi '),
 end
 disp(' fine Ksi '), pausak


else
 uL(1)=0;
 uL(2)=0;
 fatqw=0;
 lambda=lambda0*1e6;
 Fi=0;
end

%'lambda', pausak
 confzv=uL;
%if ifp>=-1 | ifp==-10
man=max(uF0.*Fi.')/3.5;
 figure, plot(hz,real(perm),'w',hz,Fi/man,'r',hz,-imag(perm)*6000,'c'),
 ax=axis;
 ax(4)=5;
 axis(ax)
  title([' lambda_{res} = ',num2str(lambda),'     Gth = ',num2str(gpQW)]), pausak
vedis=[real(perm); Fi'/man];
faperd=4e4*pi/lambda;
perd=-imag(perm)*faperd;
fi=find(perd>10000);
perd(fi)=NaN;
figure, [ax,h1,h2]=plotyy(hz,vedis,hz,perd),
  title([' lambda_{res} = ',num2str(lambda),'     Gth = ',num2str(gpQW)]), 
  xlabel(' Long coord (um)')
  ylabel(' Standing Wave & Ref. index profiles')  
%set(h1(1),'linewidth',2)
%set(h1(2),'linewidth',2)
%set(h2(1),'linewidth',2,'linestyle','--')
ylabel(ax(2),' Abs (1/cm)')  
  pausak  
 figure, plot(hz,real(perm),'w',hz,Fa*1.3,'g'), pausak
% figure, plot(xt,nt,'r.-',hz,real(perm),'w',hz,uFunc*4,hz,Fi)
% if ifp>1, pausak, end
%end
%fatqw=1;
%Fi=1;
%lambda=lambda0;
 subperd
 'subperd ', keyboard


 iloas=iloa;
 ianysa=iany;
 ilos=ilo;
 clear ilo iloa iany
 vers=version
 versi=str2num(vers(1));
 %keyboard
 if versi<7
  eval(['save ' saELM]);
  'vecchio'
 else
  eval(['save ' saELM, ' -v6' ]);
  'nuovo'
 end

 iloa=iloas;
 ilo=ilos;
 iany=ianysa;


end

disp(' dentro past_new '), keyboard

if iloa==0

if iany==1
 axorien=1;  % 1: x=[110],  -1: x=[1-10],
 r411=-1.6e-6;  % in micron/V
 r412=-1.1e-6;
 core=polyfit([0 1],[r411 r412],1);
 r41vet=polyval(core,xfit);
 n2r41=r41vet.*real(perm).^2;
 epsr41=axorien*2*r41vet.*real(perm).^4;
 fi0=find(xfit<0);
 n2r41(fi0)=0;
 epsr41(fi0)=0;


igrad_trans=igrad_tr;   %interfaccie graded


%if igrad_trans==0
%  wriwins
%else
ive=0

if ive==1
% ' sono qui ', keyboard
  wriwingr
%end

if igrad_trans==1
  wriwimir
end

 if iloa==0
  ico=0;
  for ifil=Bi
   ico=ico+1;
   if ifil==0
    chBi='00';
   else
    chBi=num2str(ifil*10);
   end
   fileName=[nomeFe,'_',chBi,'.dat'];
   disp([' loading file at V_bias = ',num2str(ifil)])
   [edu]=readsim(fileName);
   f(:,ico)=edu(:,2);
  end
  z=edu(:,1);
  parfe=' z f';
  if versi<7
   eval(['save ' nomeFe parfe]);
  else
   eval(['save ' nomeFe parfe , ' -v4']);
  end
 else
  eval(['load ' nomeFe]);
 end
 fmu=f*1e-4;
global uL
else
 global EO
 z=EO.zT';
 fmu(:,1)=mean(EO.Es,2);
 
end
'qui EO vecchio', keyboard

 fista=find(abs(Dov)>0);
 finiz=1:fista(1)-1;

 dorigin=sum(dv(finiz));
 dorigin1=sum(dv(finiz(2:end)));
 x=(hz'-dorigin)*1e-3;
 Fiany=Fi.*n2r41';
 Fianyeps=Fi.*epsr41';

 Fiwaeps=spline(x,Fianyeps,z);
 Fiwa=spline(x,Fiany,z);
 Fiw=spline(x,Fi,z);
% pewa=spline(x,imag(relPerm),z);

clear amedis ame amz amzd ameps

 dz=[diff(z); 0];
 smu=size(fmu);
 F_dz=((Fiwa.*dz)*ones(1,smu(2))).*fmu;
 F_dzeps=((Fiwaeps.*dz)*ones(1,smu(2))).*fmu;
 F_dis=((Fiw.*dz)*ones(1,smu(2))).*fmu;
 Fnorm=Fiw.*dz;
 Fnorms=Fiw'*dz;
 dmic=(dto)/1000;
 fid=find(dmic==0);
 if length(fid)>0
  dmic(fic)=1.2e-7;
 end

 fista=find(abs(Dov)>0);
 kin=fista(1);

 for k=1:length(puae)
%  li=sum(dmic(1:puae(k)-1))-dto(1)/1000;
  li=sum(dmic(1:puae(k)-1))-dorigin/1000;
  lu=li+dmic(puae(k));
  fi=find(z>=li & z<lu);
  amedis(k,:)=sum(F_dis(fi,:))/sum(Fnorm(fi));
  du=sum(F_dzeps(fi,:))/sum(Fnorm(fi));
  if isnan(du)
  du=0;
%  ' du = 0 ', keyboard
  end
  ameps(k,:)=du;
  ame(k,:)=sum(F_dz(fi,:))/sum(Fnorm(fi));
  amz(k,:)=li;
  amzd(k,:)=dmic(puae(k));
 end

 iaa=1;
 amedis=amedis(:,iaa)';
 ad=[amedis amedis];
 xd=[amz' (amz+amzd)'];
 xat=(reshape(xd',2*length(amz),1));
 %xd=[zeros(size(amz')) amz'];
 %xat=cumsum(reshape(xd',2*length(amz),1));
 at=(reshape(ad',2*length(amz),1));

 xts=(xt-dorigin)/1000;
if ifp>=-1 | ifp==-10

 figure, plot(z,fmu,x,real(sqrt(relPerm)),'r',z,Fiw*5,'w',...
          xts,nt,'c',xts,nt,'m.',xat,at,'r.'),
          pausak
 figure, plot(z,fmu(:,iaa),'r',z,Fiw*50,'b',xts,nt.^2,'g'),
 pausak
% figure, plot(z,fmu/100,x,sqrt(real(relPerm)),'y',z,Fiw*5,'w')
%          xts,nt*10,'c',xts,nt,'m.',xat,at,'r.'),
% if ifp>1, pausak, end
end


 anyv=zeros(length(dv),smu(2));

 puc=puac(find(puac>0));

 anyv(puc,:)=ameps;

% figure, plot(ameps), keyboard
 puaf=[kin:length(dto)];
 for k=1:length(puaf)
%  li=sum(dmic(1:puaf(k)-1))-dto(1)/1000;
  li=sum(dmic(1:puaf(k)-1))-dorigin/1000;
  lu=li+dmic(puaf(k));
  fi=find(z>=li & z<lu);
  if length(fi)==0
   Amedis(k)=0;
   Ame(k)=0;
  else
   sFn=sum(Fnorm(fi));
   if sFn~=0
    Amedis(k)=sum(F_dis(fi))/sFn;
    Ame(k)=sum(F_dz(fi))/sFn;
%    plot(z(fi),F_dis(fi),'.-',z(fi),Amedis(k)*ones(size(fi)),'w'), pausak
   else
    Amedis(k)=0;
    Ame(k)=0;
   end
  end
   Amz(k)=li;
   Amzd(k)=dmic(puaf(k));
 end
 fizd=find(Amzd==0);
 Amzd(fizd)=1e-7;

 Ad=[Amedis' Amedis'];
 Ade=[Ame' Ame'];
 Xd=[Amz' (Amz+Amzd)'-1e-6];
 Xat=(reshape(Xd',2*length(Amz),1));
 At=(reshape(Ad',2*length(Amz),1));
 Ate=(reshape(Ade',2*length(Amz),1));

 Atf=interp1(Xat,At,z,'linear');

 fiz=find(z<Xat(1));
 z(fiz)=Xat(1);
 Atfe=interp1(Xat,Ate,z,'linear');


if ifp>=-1 | ifp==-10
%if ifp>=-1
 figure, plot(z,fmu,x,3*real(sqrt(relPerm)),'c',z,Fiw*5,'w',...
         xts,3*nt,'b.',xat,at,'wo',Xat,At,'r.',z,Atf),
pausak
% figure, plot(z,fmu,xat,at,'wo',Xat,At,'r.',z,Atf),


  fip=find(fmu>0);
  fim=find(fmu<0);
  Fp=fmu(fip)/3;
  Fm=-fmu(fim)/3;
  figure, plot(z(fip),Fp,'y',z(fim),Fm,'r',z,Fiw*3,'w',xts,nt*2,'c')
         pausak
 if ifp>1, pausak, end
end


% return

 F_dz=Fiwa.*dz;
 ne=min(size(fmu));
 Fm=ones(ne,1)*F_dz';
 cSup=(Fm'.*fmu)/Fnorms;
 Sup=(F_dz'*fmu)/Fnorms;
 cobir=3e5/lambda;
 birt=Sup*cobir;

% fi0=find(pewa~=0);
 zcav_inf=sum(dv(1:fiC(1)-1).*abs(fst(1:fiC(1)-1,2)))-dorigin;
 zcav_sup=sum(dv(1:fiC(end)).*abs(fst(1:fiC(end),2)))-dorigin;
 z1000=z*1000;
 fi0=find(z1000<zcav_inf | z1000>zcav_sup);
 F_dz0=Fiwa.*dz;
 F_dz0(fi0)=0;
 Sup0=(F_dz0'*fmu)/Fnorms;
 bircav=Sup0*cobir

% fi0=find(pewa==0);
 fi0=find(z1000>=zcav_inf & z1000<=zcav_sup);
 F_dz0=Fiwa.*dz;
 F_dz0(fi0)=0;
 Sup0=(F_dz0'*fmu)/Fnorms;
 birmir=Sup0*cobir;


% verifica

 F_dz=Atfe.*Fiw.*dz;
 ne=min(size(Atfe));
 Fm=ones(ne,1)*F_dz';
 inonan=find(isnan(F_dz)==0);
 Supv=sum(F_dz(inonan))/Fnorms;
 cSupv=(F_dz)/Fnorms;
 cobir=3e5/lambda;
 birt_verifica=Supv*cobir;
 disp(' ');
 disp(' verifica Bir. ');
 [birt       birt_verifica]
 disp(' ');
 keyboard

 if ifp>=-1 | ifp==-10
  figure, plot(z,cumsum(cSup)*cobir,z,cumsum(cSupv)*cobir,'-'), title(' fy-fx')
  pausak
  if ifp>1, pausak, end
 end

 any.t=anyv(put,:);
 any.b=anyv(pub,:);
 any.a=anyv(pua,:);

end  %iany

else
 z=0;
 any=0;
 fmu=0;
end

if iany>=2
 anyv=anyf;
 any.t=anyv(put,:);
 any.b=anyv(pub,:);
 any.a=anyv(pua,:);
end 



disp(' fine past_new  e birif: ferma per verif')
pausak
keyboard
keyboard
%keyboard
