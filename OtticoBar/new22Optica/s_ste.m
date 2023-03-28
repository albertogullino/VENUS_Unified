
ipim=1;         %=1 uses pimra as QW losses for N=0, =0, compute losses at N=0
%pimra=-0.01;
pimra=0;
icalcFi=1;
iany0=1;
if iany==1 | idyn==1
 icalcFi=1;
end

istrumix=1;
icampi=1;   % 0: no,  1: solo uscita,  2:tutti

fatd=1;
imet=0;    % dice se metallo uscita e` rilevante
Temp=Tvet(1);

  disp([' pausa verifica input ', nvar])
%  ilo
%  keyboard

  nomeFs='STEIN.STR';
  lambda0=847e-9;

%  nomeFs='PROA.STR';
%  lambda0=9.809298482353390e-007;
%  nomeFs='ARRAY.STR';    %josip
%  lambda0=9.70e-007;

  rad=nomeFs(1:end-4);
  nomeFe=rad;
  nomeloELM=[rad,'_ELM.mat'];

  Bi=2;


das=0;
if length(ilo)==0
   ibd=0;
   ibs=0;
   DR=dir;
   lD=length(DR);
   for kD=3:lD
    dun=getfield(DR(kD),'name');
    if length(dun)==length(nomeFs)
     if strcmpi(dun,nomeFs)==1
      dad=datenum(getfield(DR(kD),'date'));
      ibd=1;
     end
    end
    if length(dun)==length(nomeloELM)
     if strcmpi(dun,nomeloELM)==1
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


%[du,iaa]=min(abs(Bi-Bias));
iaa=1;

 if ilo==0
  [fatqw,rr,a0ref,Lr,Lpl,nr,npl,xr,xpl,nmir,ar,shav,frp,lambdame,Fi,any,...
   Fist,zsta,anys,gam_v,fiQ0,fiC0,Lpla,npla,Nspla,aral,...
   nv,dv,xm,radii,fst,perd,anyf,iauto,shavet,ipar,ifield,istfi,...
   icaplo,icsfit,icsfib,icsfi,istmet]=...
   parst_pa(nomeFs,lambda0,ifp,iany0,ianys,nomeFe,Bi,iloa,icalcFi,imet);

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

  lambda0=9.686737242433540e-007;
 if icalcFi==0
  lambda=lambda0*1e6;
 else
  lambda=lambdame*1e6;
 end
 k0=2*pi/lambda*1e6;
 confztot=gam_v(2);

 ass_par


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
   if length(Lr.t)>0
    dbar=Lr.t(length(Lr.t))*1e-6;
   else
    dbar=Lr.b(1)*1e-6;
   end
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

%  L_ud=Lvtr;
%  n_ud=nvtr;
%  for isr=1:nstratiu
%   L_ud=[L_ud; Lbt];
%   n_ud=[n_ud; nbt];
%  end
%  L_ud=[L_ud; Luv];
%  n_ud=[n_ud; nuv];
%  L_u=flipud(L_ud);
%  n_u=flipud(n_ud);
%  a_u=zeros(size(n_u));
%
%  L_b=Lvbr;
%  n_b=nvbr;
%  for isr=1:nstratid
%   L_b=[L_b; Lbb];
%   n_b=[n_b; nbb];
%  end
%  L_b=[L_b; Luvb];
%  n_b=[n_b; nuvb];
%  a_b=zeros(size(n_b));
%
%  fitustra=find(fmlstot(:,1)>1);
%
%  if length(fitustra)>0
%    lLi=length(Litot);
%    Litottu=[];
%    nitottu=[];
%    aitottu=[];
%    kL=1;
%    istadd=0;
%    while kL<lLi
%     fdu=fmlstot(kL,1);
%     if fdu==1
%      Litottu=[Litottu; Litot(kL)];
%      nitottu=[nitottu; nitot(kL,1)];
%      aitottu=[aitottu; aitot(kL)];
%      kL=kL+1;
%     else
%      pus=kL:kL+fdu-1;
%      pit=fmlstot(kL,2);
%      for kp=1:pit
%       Litottu=[Litottu; Litot(pus)];
%       nitottu=[nitottu; nitot(pus,1)];
%       aitottu=[aitottu; aitot(pus,1)];
%       if kL<fiQ0(1)
%        istadd=istadd+length(pus);
%       end
%      end
%       if kL<fiQ0(1)
%        istadd=istadd-length(pus);
%       end
%      kL=kL+fdu;
%     end  % fdu
%    end   % while
%
%  else
%    Litottu=Litot;
%    nitottu=nitot(:,1);
%    aitottu=aitot;
%  end    % if fitustra
%
%
%
%  L_i=[L_u; Litottu; L_b];
%  n_i=[n_u; nitottu; n_b];
%  a_i=[a_u; aitottu; a_b];
%
%  fiQ=fiQ0+istadd;
%  fiCav=fiC0+istadd;
%
%  if ifp>-1
%    Ldu=[L_i*0 L_i];
%    Ld=reshape(Ldu',prod(size(Ldu)),1);
%    ndu=[n_i n_i];
%    nd=reshape(ndu.',prod(size(Ldu)),1);
%    disp(' fiQ ')
%    xplo=cumsum(Ld);
%    puQve=[2*fiQ-1 2*fiQ];
%    puCve=[2*fiCav-1 2*fiCav];
%    figure, plot(xplo,real(nd),xplo(puQve),real(nd(puQve)),'ro',...
%                 xplo(puCve),real(nd(puCve)),'g.')
%    keyboard
%   end
%



isk=1;
if isk==0
disp(' nuovo '), keyboard

  fstr=fst(:,1);
  fi=find(fstr==0);
  fstr(fi)=1;
  cstr=abs(fst(:,2));

  L_i=[];
  n_i=[];
  a_i=[];
  ics=0;
  icsc=0;
  istac=0;

  for iv=1:length(fstr)
  iv
  pausak
   for ist=1:cstr(iv)
  ist
  pausak
    for ist1=1:fstr(iv)
     iv1=iv-1+ist1;
     ics=ics+1;
     L_i=[L_i; dv(iv1)];
     n_i=[n_i; nv(iv1,1)];
     a_i=[a_i; radii.a(iv1,1)];
     if (iauto(iv1,2)==-4 | iauto(iv1,2)>0)
      if istac==0
       istac=1;
      else
       icsc=icsc+1;
       fiC0n(icsc)=ics;
       istac=0;
      end
     end
     if iauto(iv1,1)==2
      istadd=ics;
     end
     if istac==1
      icsc=icsc+1;
      fiC0n(icsc)=ics;
     end
    end
   end
  end

%  fiQ=fiQ0+istadd;
%  fiCav=fiC0+istadd;
  fiQ=istadd;
  fiCav=fiC0n;
  L_i=L_i/1000;

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




iprorel=0;
%keyboard
if iprorel==1 & ilo==0
 disp(' qui per relief '), keyboard
 pro_rel
end
 disp(' s_ste')
 keyboard
end
