
ipim=1;         %=1 uses pimra as QW losses for N=0, =0, compute losses at N=0
%pimra=-0.01;
pimra=0;
if iany==1 | idyn==1
 icalcFi=1;
end
icalcFi=1;

istrumix=1;
icampi=1;   % 0: no,  1: solo uscita,  2:tutti

fatd=1;
imet=0;    % dice se metallo uscita e` rilevante
Temp=Tvet(1);

  disp([' pausa verifica input ', nvar])
%  ilo
%  keyboard

  nomeFs='STEIN.STR';

  rad=nomeFs(1:end-4);
  nomeFe=rad;
  nomeloELM=[rad,'_ELM.mat'];

  lambda0=847e-9;
  Bi=2;


das=0;
if length(ilo)==0
   ibd=0;
   ibs=0;
   DirNam=[cd, '/D1ANA/GenerateGeom'];
   DR=dir(DirNam);
   lD=length(DR);
   for kD=3:lD
    dun=getfield(DR(kD),'name');
    if length(dun)==length(nomeFs)
     if dun==nomeFs
      dad=datenum(getfield(DR(kD),'date'));
      ibd=1;
     end
    end
    if length(dun)==length(nomeloELM)
     if dun==nomeloELM
      das=datenum(getfield(DR(kD),'date'));
      ibs=1;
     end
     if ibs*ibd==1, break, end
    end
   end
   if das-dad<0
    ilo=0
   else
    ilo=1
   end
else
 ilo=1;
end


  iloa=ilo;  % carica anysotropia memorizzata
%  iloa=0;  % carica anysotropia memorizzata


[du,iaa]=min(abs(Bi-Bias));

 if ilo==0
  [fatqw,rr,a0ref,Lr,Lpl,nr,npl,xr,xpl,nmir,ar,shav,frp,lambdame,Fi,any,...
   Fist,zsta,anys,gam_v,fiQ0,fiC0,Lpla,npla,Nspla,aral,...
   nv,dv,xm,radii,fst,perd,anyf,iauto,shavet,ipar,ifield,istfi,...
   icaplo,icsfit,icsfib,icsfi,istmet]=...
   parst_pa(nomeFs,lambda0,ifp,iany,ianys,nomeFe,Bi,iloa,icalcFi,imet);

%  [fatqw,rr,a0ref,Lr,Lpl,nr,npl,xr,xpl,nmir,ar,shav,frp,lambdame,Fi,any,Fist,...
%  zsta,anys,gam_v,fiQ0,fiC0,Lpla,npla,Nspla,shav,aral,shavet]=...
%  parst_pa(nomeFs,lambda0,ifp,iany,ianys,nomeFe,Bi,iloa,icalcFi);
  disp('dopo parst_n ')
  keyboard
  sva=' fatqw rr Lr Lpl nr npl xr xpl nmir ar frp lambdame ';
  sva1=' Fi any Fist zsta anys gam_v fiQ0 fiC0 Lpla npla Nspla';
  sva2=' shav aral shavet radii a0ref ';
  sva3=' nv dv xm radii fst perd anyf iauto shavet ipar ifield istfi';
  sva4=' icaplo icsfit icsfib icsfi istmet ';
  eval(['save ' nomeloELM  sva sva1 sva2 sva3 sva4 ]);
 else
  eval(['load ' nomeloELM]);
 end

 if icalcFi==0
  lambda=lambda0*1e6;
 else
  lambda=lambdame*1e6;
 end

% ishamem=find(ifield<0);
% sad=[1 0 -1];
% for k=1:3
%  fiir(k)=find(iauto(:,1)==k)+sad(k);
% end

% cfi=0;
% if length(ishamem)>0
%  for k=1:length(ishamem)
%   fi=find(fiir==ishamem(k));
%   if length(fi)==0
%    cfi=cfi+1;
%    fime(cfi)=ishamem(k);
%   end
%  end
%  if cfi>=1
%   istfiv=fime-1;
%  end
% end

% ifimem=find(ifield~=0);
% fime=zeros(size(ifield));
% fime(ifimem)=ifimem-1;








 if length(par_in)>0

  sip=size(ipar);
  riga=[1:sip(1)]';
  righe=riga*ones(1,sip(3));
  ip1=reshape(ipar(:,1,:),sip(1)*sip(3),1);   % par. #
  ip2=reshape(ipar(:,2,:),sip(1)*sip(3),1);   % par. type
  ip3=reshape(ipar(:,3,:),sip(1)*sip(3),1);   % trans. sect.
  ip4=reshape(righe,sip(1)*sip(3),1);         % layer
  fi1=find(ip1~=0);
  [dup,fi1p]=sort(ip1(fi1));
  fi1=fi1(fi1p);

%  fid=find(diff([0; dup])~=0);
%  fiud=fi1(fid);
%  [ip11,fiudd]=sort(ip1(fiud));
%  fiu=fiud(fiudd);
  fiu=fi1;

  npar=length(find(diff([0; dup])~=0));


%  dup=ip1(fi1);
%  fid=find(diff([0; dup])~=0);
%  fiud=fi1(fid);
%  [ip11,fiudd]=sort(ip1(fiud));
%  fiu=fiud(fiudd);
%
%  npar=length(find(diff([0; dup])~=0));

  if length(par_in)~=npar

   disp(' error: par_in does not match the number of parameters ')
   keyboard
   return
  else

   for ip=1:length(par_in)
    fip=find(ip1(fiu)==ip);
%    disp(' loo 1 '), pausak
    if length(fip)>0
     pup=fiu(fip);
     pur=ip4(pup);
     puc=ip3(pup);
     ipad=ip2(pup);
     ipa=ipad(1);
%     disp(' loo 2 '), pausak
     if ipa==1
      dv(pur)=par_in(ip);
     elseif ipa==2
      radii.a(pur,puc)=par_in(ip);
     elseif ipa==3
      radii.b(pur,puc)=par_in(ip);
     elseif ipa==4
      radii.c(pur,puc)=par_in(ip);
     elseif ipa==5
      nv(pur,puc)=par_in(ip);
     end
    end
   end

   reassign

  end

 end  %length(par_in)>0


  if iany>0
   aniatv=any.a;
   deltanynv=any.t;
   deltanybv=any.b;

   aniat    =aniatv(:,iaa);
   deltanyn=deltanynv(:,iaa);
   deltanyb=deltanybv(:,iaa);
   Xiu=zsta;
   Fiu=Fist(:,iaa);
  end

  if ianys>0
   aniatsv=anys.a;
   deltanynsv=anys.t;
   deltanybsv=anys.b;
   aniats    =aniatsv*Strain;
   deltanyns=deltanynsv*Strain;
   deltanybs=deltanybsv*Strain;
  end

% Set parametri per autoconsistenza (rela_new.m)
%
% per struttura planare - autoconsistenza

  nstratid=nmir.b;
  Lvbr=Lpl.b.i;
  Luvb=Lpl.b.o;
  Lbb =Lpl.b.m;
  nvbr=npl.b.i;
  nuvb=npl.b.o;
  nbb =npl.b.m;
  rfd =npl.b.last;

  Lvb=[Lr.b; Lpl.b.i];
  nvb=[nr.b(:,1); npl.b.i];

  nstratiu=nmir.t;
  Lvtr=Lpl.t.i;
  Luv=Lpl.t.o;
  Lbt =Lpl.t.m;
  nvtr=npl.t.i;
  nuv=npl.t.o;

  nbt =npl.t.m;
  rfu =npl.t.last;

% per struttura non planare e/o  anisotropa

  Litn=Lr.t;
  nitn=nr.t.';


  aitn=ar.t.';

  fmlsp=abs(frp.t);
  fmlst=abs(frp.t);
  fi1=find(fmlsp(:,1)==0);
  fmlst(fi1,1)=1;
  fi2=find(fmlsp(:,2)==1);
  fmlst(fi1,2)=0;



  Lib=Lr.b;
  nib=nr.b.';
  aib=ar.b.';

  fmlsp=abs(frp.b);
  fmlsb=abs(frp.b);
  fi1=find(fmlsp(:,1)==0);
  fmlsb(fi1,1)=1;
  fi2=find(fmlsp(:,2)==1);
  fmlsb(fi1,2)=0;

  d=Lr.a*1e-6;
  L_QW=Lr.a;
  NQW=nmir.a;
  if NQW>1
   dbar=Lr.t(length(Lr.t))*1e-6;
  else
   dbar=0;
  end
  niat=nr.a.';
  aiat=ar.a.';
  rqw=niat(1);
  ipilat=0;
  if length(aiat)>1
   if aiat(2)>0
    ipilat=1;
   end
  end


% per g0

    if ipim==0
     lax=(lambda+Dla)/(1+iafr*frequer);
     if ifit==0
      [coeffgv,coeffnv]=coef(lax);
      fi=find(coeffgv-coeffnv~=0);
      coeffgv(fi)=fatqw*coeffgv(fi);
      [g0]=NQW*ppvald(coeffgv,N00);
     else
          if isolu==1
           Gf=interpn(port,lav,Tv,Gtot,N00,lax,Temp,'spline');
          else
           Gf=interp2(port',lav,Gtot.',N00,lax);
          end
      g0=iperat*NQW*fatqw*real(Gf);
     end

     c0=3e5; %in micron/ns
     fconv=1/(c0*k0);
     ra=rqw*(1+j*g0*fconv);
    else
     c0=3e5; %in micron/ns
     ra=rqw+j*pimra;
     cconv=k0*c0/rr;
     g0=pimra*cconv;
    end


   pucavi=1;

   Litot=[Litn; d*1e6; Lib];
   aitot=[aitn.'; aiat'; aib.'];
   nitot=[nitn.'; niat.'; nib.'];
   fmlstot=[fmlst; [1 0] ; fmlsb];


% per struttura completamente planare (perdite e acc. uscita in TL_ef.m)
%
% sopra
%
Lplit=Lpla.t.i;
nplit=npla.t.i;
%
Lplt=Lpla.t.m;
nplt=npla.t.m;
Nsplt=Nspla.t.m;
%
Lplut=Lpla.t.o;
nplut=npla.t.o;
%
%%
%
% sotto
%
Lplib=Lpla.b.i;
nplib=npla.b.i;
%
Lplb=Lpla.b.m;
nplb=npla.b.m;
Nsplb=Nspla.b.m;
%
Lplub=Lpla.b.o;
nplub=npla.b.o;


% per struttura completamente planare (perdite e acc. uscita in gaz7.m)

  L_ud=Lvtr;
  n_ud=nvtr;
  for isr=1:nstratiu
   L_ud=[L_ud; Lbt];
   n_ud=[n_ud; nbt];
  end
  L_ud=[L_ud; Luv];
  n_ud=[n_ud; nuv];
  L_u=flipud(L_ud);
  n_u=flipud(n_ud);
  a_u=zeros(size(n_u));

  L_b=Lvbr;
  n_b=nvbr;
  for isr=1:nstratid
   L_b=[L_b; Lbb];
   n_b=[n_b; nbb];
  end
  L_b=[L_b; Luvb];
  n_b=[n_b; nuvb];
  a_b=zeros(size(n_b));

  fitustra=find(fmlstot(:,1)>1);

  if length(fitustra)>0
    lLi=length(Litot);
    Litottu=[];
    nitottu=[];
    aitottu=[];
    kL=1;
    istadd=0;
    while kL<lLi
     fdu=fmlstot(kL,1);
     if fdu==1
      Litottu=[Litottu; Litot(kL)];
      nitottu=[nitottu; nitot(kL,1)];
      aitottu=[aitottu; aitot(kL)];
      kL=kL+1;
     else
      pus=kL:kL+fdu-1;
      pit=fmlstot(kL,2);
      for kp=1:pit
       Litottu=[Litottu; Litot(pus)];
       nitottu=[nitottu; nitot(pus,1)];
       aitottu=[aitottu; aitot(pus,1)];
       if kL<fiQ0(1)
        istadd=istadd+length(pus);
       end
      end
       if kL<fiQ0(1)
        istadd=istadd-length(pus);
       end
      kL=kL+fdu;
     end  % fdu
    end   % while

  else
    Litottu=Litot;
    nitottu=nitot(:,1);
    aitottu=aitot;
  end    % if fitustra



  L_i=[L_u; Litottu; L_b];
  n_i=[n_u; nitottu; n_b];
  a_i=[a_u; aitottu; a_b];

  fiQ=fiQ0+istadd;
  fiCav=fiC0+istadd;
  if ifp>-1
    Ldu=[L_i*0 L_i];
    Ld=reshape(Ldu',prod(size(Ldu)),1);
    ndu=[n_i n_i];
    nd=reshape(ndu.',prod(size(Ldu)),1);
    disp(' fiQ ')
    xplo=cumsum(Ld);
    puQve=[2*fiQ-1 2*fiQ];
    puCve=[2*fiCav-1 2*fiCav];
    figure, plot(xplo,real(nd),xplo(puQve),real(nd(puQve)),'ro',...
                 xplo(puCve),real(nd(puCve)),'g.')
    keyboard
   end


   if ifp>-4, pausak, end




iprorel=1;
%keyboard
if iprorel==1 & ilo==0
 disp(' qui per relief '), keyboard
 pro_rel
end

