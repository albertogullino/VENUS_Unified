%ipim=1;         %=1 uses pimra as QW losses for N=0, =0, compute losses at N=0
%pimra=-0.01;
ifit=1;
iperat=1;


global pimra ipim fattany
global alpha_th g0

%'pimra', keyboard
%'pimra', keyboard
%'pimra', keyboard
%'pimra', keyboard

if length(pimra)==0
 pimra=0;
end
%' pimra ', keyboard
if length(ipim)==0
 ipim=1;
end

if length(fattany)==0 & iany==1
 ' aggiungere nel main: '
 ' global fattany; fattany=1 '
 ' fattany=1: (x -> 011, y -> 01-1);  fattany=-1 (y -> 011, x -> 01-1) '
 keyboard
end
%fattany=1;


Strain=0;

%icalcFi=0;
%if iany==1 | idyn==1
% icalcFi=1;
%end

if ianys>0
% disp(' ')
% disp(' predisporre per strain !!! ')
% keyboard
 disp(' ')

 global Ps
 Strain=Ps.Strain
  disp(' Fissato strain in Ps.Strain !!! ')
  if ifp==-10 
   pausak
  end
% return
end

istrumix=1;
if ~exist('icampi')
 icampi=1;   % 0: no,  1: solo uscita,  2:tutti
end
fatd=1;
imet=1;    % dice se metallo uscita e` rilevante
%'imet ', keyboard

Temp=Tvet(1);


  fileName=nvar;
  
     fis= strfind(fileName,'\');
     filNam=fileName(fis(end)+1:end);
     DirNam=fileName(1:fis(end));
     
       nomeFs=filNam;
       rad=nomeFs(1:end-4);
       nomeloFEM=[rad,'_FEM.mat'];
     
   DR=dir(DirNam);


  rad=nomeFs(1:end-4);
  nomeFe=rad;
  saELM=[rad,'_SA.mat'];
  nomereadELM=[rad,'_ELM.mat'];
  nomereadANY=[rad,'_ANY.mat'];
  nomeloELM=[DirNam,rad,'_ELM.mat'];
  nomeloANY=[DirNam,rad,'_ANY.mat'];  
%'ELM', keyboard

  %disp([' pausa verifica input ', nvar]),  keyboard


  Bi=2;
% Bi=[0 1 2 3];


das=0;
das1=0;

%if length(ilo)==0

   ibd=0;
   ibs=0;
   ibs1=0;
   lD=length(DR);
   for kD=3:lD
    dun=getfield(DR(kD),'name');
%    ' dad', pausak
     if strcmpi(dun,nomeFs)==1
      dad=datenum(getfield(DR(kD),'date'));
      ibd=1;
     end
     if strcmpi(dun,nomereadELM)==1
      das=datenum(getfield(DR(kD),'date'));
      ibs=1;
     end
     if strcmpi(dun,nomereadANY)==1
      das1=datenum(getfield(DR(kD),'date'));
      ibs1=1;
     end
     if ibs*ibd*ibs1==1, break, end
   end
   
   if dad-das>0
    ilo=0;
   else
    ilo=1;
   end
   if das1-dad<0
    iload=0;
   else
    iload=1;
   end
%else
% ilo=1;
% iload=0;
%end

if iany==1 & iload==0
 iloa=0;
elseif iany==1 & iload==1
 iloa=1;
else
 iloa=-1;
end

%disp(' contr set_stru'), keyboard
%disp(' contr set_stru'), keyboard

 if ilo==0
% disp(' contr set_stru dentro ilo'), keyboard
 nomeFs=fileName;
  [fatqw,rr,a0ref,Lr,Lpl,nr,npl,xr,xpl,nmir,ar,shav,frp,lambda,Fi,any,...
   Fist,zsta,gam_v,fiQ0,fiCav0,Lpla,npla,Nspla,aral,...
   nv0,dv,xm,radii,fst,perd,anyf,iauto,shavet,ipar,ifield,istfi,...
   icaplo,icsfit,icsfib,icsfi,istmet,dglo,iFFsect,du,du,confzv,lambda0]=...
   past_new(nomeFs,ifp,iany,nomeFe,Bi,iloa,ilo,imet,saELM);

%   if idyn==1
    per_dyn
%   end

 if ifp~=-4
  disp('dopo past_new '), pausak
%  keyboard
 end
  sva=' fatqw rr Lr Lpl nr npl xr xpl nmir ar frp lambda confzv lambda0 ';
  sva1=' Fi gam_v fiQ fiCav Lpla npla Nspla';
  sva2=' shav aral shavet radii a0ref L_i n_i';
  sva3=' nv0 dv xm radii fst perd anyf iauto shavet ipar ifield istfi';
  sva4=' icaplo icsfit icsfib icsfi istmet dglo iFFsect any';
  vers=version
%  keyboard
  versi=str2num(vers(1));
  %'quao save\', keyboard
  if versi<7
%  if vers(1)=='6'
   eval(['save ' nomeloELM  sva sva1 sva2 sva3 sva4 ]);
  else
   eval(['save ' nomeloELM  sva sva1 sva2 sva3 sva4, ' -v6'  ]);
  end

%  eval(['save ' nomeloELM  sva sva1 sva2 sva3 sva4 ]);
  if iloa==0
   sva0=' any Fist zsta ';
   if vers(1)=='6'
    eval(['save ' nomeloANY sva0 ]);
   else
%    eval(['save ' nomeloELM  sva sva1 sva2 sva3 sva4, ' -v6'  ]);
    eval(['save ' nomeloANY sva0 , '-v6']);
   end
   iloa=-2;
  end
  ilo=-2;
 
 else
%  'qui carico', keyboard
 
  eval(['load ' nomeloELM]);
  
   
   global nrLosses
   if length(nrLosses)>0
    nv0=nrLosses;
   end 
  
 
 end
% 'TOGWELERUI!!!!!!!!'
%ilo=0

 if iloa==0
% disp(' contr set_stru'), keyboard
  [fatqw,rr,a0ref,Lr,Lpl,nr,npl,xr,xpl,nmir,ar,shav,frp,lambda,Fi,any,...
   Fist,zsta,gam_v,fiQ0,fiCav0,Lpla,npla,Nspla,aral,...
   nv0,dv,xm,radii,fst,perd,anyf,iauto,shavet,ipar,ifield,istfi,...
   icaplo,icsfit,icsfib,icsfi,istmet,dglo,iFFsect,du,du,confzv]=...
   past_new(nomeFs,ifp,iany,nomeFe,Bi,iloa,ilo,imet,saELM);

%   if idyn==1
    per_dyn
%   end

 if ifp~=-4
  disp('dopo past_new ')
  keyboard
 end

   sva0=' any Fist zsta  ';
   eval(['save ' nomeloANY sva0 ]);
   iloa=-2;
   if ilo==0
    sva=' fatqw rr Lr Lpl nr npl xr xpl nmir ar frp lambda confzv ';
    sva1=' Fi gam_v fiQ fiCav Lpla npla Nspla';
    sva2=' shav aral shavet radii a0ref L_i n_i';
    sva3=' nv0 dv xm radii fst perd anyf iauto shavet ipar ifield istfi';
    sva4=' icaplo icsfit icsfib icsfi istmet dglo iFFsect ';
    eval(['save ' nomeloELM  sva sva1 sva2 sva3 sva4 ]);
    ilo=-2;
   end
 elseif iloa==1
  eval(['load ' nomeloANY]);
 end

%'set_stru', keyboard
 global Ps
 if isfield(Ps,'loss_ch')==0
  Ps.loss_ch=0;
 end

 if Ps.loss_ch==1
  fipe=find(abs(imag(nv0)<0.1) & real(nv0)~=0 & real(nv0)>2);
  nv0(fipe)=real(nv0(fipe));
 end
 if Ps.loss_ch<0
  fipe=find(abs(imag(nv0)<0.1) & real(nv0)~=0 & real(nv0)>2);
  k0=2*pi/lambda*1e6;
  nv0(fipe)=real(nv0(fipe))+j*Ps.loss_ch/(2*k0);
 end
iold=0;
if iold==1
 if isfield(Ps,'radd')==1
  radd=Ps.radd;
  ra0=radd(1);
  imesaK=radd(4);
  if ra0<0
   imesaK=1;
  end
  if ra0<0
   fira=find(radii.a>=abs(ra0));
   if length(fira)>0
    radii.a(fira)=0;
   end
  else
   rad=radd(2);
   stra=radd(3);
   fira=find(radii.a>=ra0);
   if length(fira)>0
    asho=radii.a(fira);
    firep=length(find(asho==asho(1)));
    npra=length(asho)/firep;
    stv=linspace(0,stra,npra);
    stvt=repmat(stv,firep,1);
    stvtu=reshape(stvt,prod(size(stvt)),1);
    radii.a(fira)=asho(1)+rad+stvtu;
    if ifp==-10
    'radd', keyboard
    end
   end
  end  %ra0
 else
  ra0=1000;
  imesaK=0;
 end
end

  ra0=1000;
  imesaK=0;

%disp(' set_stru'), keyboard
 if ~exist('iFFsect')
  iFFsect=1;
 end

 fief=find(abs(ifield)>10);
 if length(fief)>0
  iFF=1;
  z0c=(ifield(fief)/1e4);
  if z0c<0
   iFFte=1;
   z0c=abs(z0c);
  else
   iFFte=0;
  end


  ifield(fief)=0;
 else
  iFF=0;
 end

%'setstru prima asspar', keyboard
 k0=2*pi/lambda0*1e6;
 k0cm=2*pi/lambda0*1e-2;
 confztot=1;


 ass_par

if i1D==1
  return
end
 
%'setstru dopo asspar', keyboard


%[du,iaa]=min(abs(Bi-Bias));
iaa=1;

  if iany>0
   aniatv=any.a;
   deltanynv=any.t;
   deltanybv=any.b;

   aniat    =aniatv(:,iaa);
   deltanyn=deltanynv(:,iaa);
   deltanyb=deltanybv(:,iaa);
   if iany==100
    Xiu=zsta;
    Fiu=Fist(:,iaa);
   end 
  end

  if ianys>0
   p44=-0.072;
   anysv=-p44*ones(size(dv)).*real(nv0(:,1)).^4;
   anys.t=anysv(put,:);
   anys.b=anysv(pub,:);
   anys.a=anysv(pua,:);

   aniatsv=anys.a;
   deltanynsv=anys.t;
   deltanybsv=anys.b;
   aniats    =aniatsv*Strain;
   deltanyns=deltanynsv*Strain;
   deltanybs=deltanybsv*Strain;
  else
   anys=0;
  end



% '% per g0', keyboard

   if ipim==0
    if length(g0)==0
      eval(['load ',fileeps]);
      imatold=0;
      if size(Gtot)>2
       isolu=1;
      else
       isolu=0;
      end
      lax=(lambda+Dla);
      if imatold==0
       [du,iT]=min(abs(Tv-Temp));
       [du,iL]=min(abs(lav-lax));
       Gs=reshape(Gtot(:,iL,iT),length(port),1);
       [Gs0,iG0]=min(real(Gs));
       g0=iperat*NQW*fatqw*Gs0/2;
       if ifp>=0
        figure, plot(port,real(Gs),port(iG0),Gs0,'wo'), pausak
        Gn=imag(Gs);
        Gn=Gn-Gn(1);
        figure, plot(port,Gn)
        pausak
        vg0=3e8;
        convn=iperat*NQW*fatqw*1e9/(vg0*k0);
        figure, plot(port,Gn*convn)
        pausak
       end

      else


       N00=min(port);
      'N00', keyboard
% Dla=-.01;
% lax=linspace(.82,.87,50);
% N00=[.01 .5 5 10];
%           Gf=interpn(port,lav,Tv,Gtot,N00,lax+Dla,Temp,'spline');
% figure, plot(lax,real(Gf))
% Dla=-.01;
% N00=linspace(0,15,50);
% lax=.85;
%           Gf=interpn(port,lav,Tv,Gtot,N00,lax+Dla,Temp,'spline');
% figure, plot(N00,real(Gf))

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
       g0=g0/2;

      end  %imatold

    end

     c0=3e-1;   % m/ns
     fconv=100/(2*k0*rqw);
     ra=rqw*(1+j*g0*fconv);
     fconv=rr/(c0*k0);
     ra=real(rqw)+j*g0*fconv;


   else

%'qui cont', keyboard
%     pimra=imag(niat(1));
     NQW_eq=uL(2)/uL(1);
     ra=real(rqw)+j*pimra*NQW_eq;
     cconv=k0cm;
     g0=pimra*cconv*NQW_eq;
   end
   niat(1)=ra;

% move 2
%g0=0
%' ipim', keyboard
%keyboard
%keyboard

% per struttura completamente planare (perdite e acc. uscita in gaz7.m)


isk=1;
if isk==0
 disp(' nuovo '), keyboard
 per_dyn

 disp(' fine '), keyboard



 iprorel=1;
 %if iprorel==1 & ilo==0
 if iprorel==1
  disp(' qui per relief '), keyboard
  pro_rel
 end
end  %isk

if ~exist('iPERD_vol')
 iPERD_vol=1;
end
% 'qui rfd',keyboard

if iPERD_vol==0
 fi=find(imag(nitn(1,:))>-1e-1);
 nitn(1,fi)=real(nitn(1,fi));
 fi=find(imag(nib(1,:))>-1e-1);
 nib(1,fi)=real(nib(1,fi));

 % metto a 0 le parti immaginarie
 if exist('n_i')
  n_i=ima0(n_i);
 end

 if exist('nbb')
  nbb=ima0(nbb);
 end
 if exist('nuvb')
  nuvb=ima0(nuvb);
 end
 if exist('rfd')
  rfd=ima0(rfd);
 end
 if exist('nvb')
  nvb=ima0(nvb);
 end
 if exist('nvbr')
  nvbr=ima0(nvbr);
 end

 if exist('nbt')
  nbt=ima0(nbt);
 end
 if exist('nuv')
  nuv=ima0(nuv);
 end
 if exist('rfu')
  rfu=ima0(rfu);
 end
 if exist('nvt')
  nvt=ima0(nvt);
 end
 if exist('nvtr')
  nvtr=ima0(nvtr);
 end
end


% disp(' set_stru'), keyboard
