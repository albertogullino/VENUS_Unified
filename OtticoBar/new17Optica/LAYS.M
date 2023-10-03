function [z,ro,th,Sig,NQW,lambdaNot,d,radii,shav,flsv,iauto]=Lay_ct(fileName);

  nomeFs=fileName;
  rad=nomeFs(1:end-4);
  nomeloFEM=[rad,'_FEM.mat'];


  das=0;
  das1=0;


   ibd=0;
   ibs=0;
   ibs1=0;
   DR=dir;
   lD=length(DR);
   for kD=3:lD
    dun=getfield(DR(kD),'name');
     if strcmpi(dun,nomeFs)==1
      dad=datenum(getfield(DR(kD),'date'));
      ibd=1;
     end
     if strcmpi(dun,nomeloFEM)==1
      das=datenum(getfield(DR(kD),'date'));
      ibs=1;
     end
     if ibs*ibd*ibs1==1, break, end
   end
   if das-dad<0
    ilo=0;
   else
    ilo=1;
   end
   if das1-dad<0
    iload=0;
   else
    iload=1;
   end

 if ilo==0

  fid=fopen(fileName,'r');
  Nline=fgetl(fid);

  ch0='% DATAFORFEM_START';
  ch1='% DATAFORFEM_END';
  if strcmp(Nline,ch0)~=1
   'lay_ct.m: file non adatto', keyboard
  end


   while feof(fid)==0

      Nline=fgetl(fid);
      if strcmp(Nline,ch1)==1, break, end



      igg=0;
      if length(findstr(Nline,'mesa'))==1
        ro.mesa=extrac_n(Nline,'radius=');
        th.mesa=extrac_n(Nline,'depth=');
      elseif length(findstr(Nline,'ring'))==1
        ro.met=extrac_n(Nline,'radius=');
        if length(findstr(Nline,'temp_include=yes'))==0
         ro.met=-ro.met;
        end
      elseif length(findstr(Nline,'buffer'))==1
       if findstr(Nline,'thick=')
        th.buf=extrac_n(Nline,'thick=');
       else
        th.Btran=extrac_n(Nline,'trans=');
       end
      elseif length(findstr(Nline,'conduction'))==1
        ro.cond=extrac_n(Nline,'radius=');
      elseif length(findstr(Nline,'confin'))==1
         ipu=findstr(Nline,'GG=');
         [Tem,remainder]=strtok(Nline(ipu+3));
         if Tem=='B'
          igg=-1;
         elseif Tem=='T'
          igg=1;
         else
          ' errore file .str  da lay_ct ', keyboard
         end
      end

   end


   fclose(fid);


   [nref,aref,n,D,xm,Dov,radii,flsv,perd,anyf,iauto,...
         dw,xw,fsw,dov,shav,ipar,ifield,lambdaNot]= Lay_new(fileName);
%   'lay_new in lay_ct', keyboard

   d=D*1e-3;


   sn=size(n);
   dv=repmat(d,sn(2),1);
   sa=size(radii.a);
   av=repmat(radii.a,sa(2),1);
  sn=size(n)
%  if sn(2)>1
  if igg==0
   fi=find(n(:,2)>1.2 & n(:,2)<1.8);
   if length(fi)>0
    doxiv=dv(fi);
    aoxiv=av(fi);
    [ro.ox,ipu]=min(aoxiv);
    ip=fi(ipu);
    th.ox=doxiv(ipu);
    fi0=fi-length(d);
    fim=find(n(:,1)>1);
    ipox=fim(1):ip-1;
    ds=d(ipox).*flsv(ipox,2);
    z.ox=sum(ds);
   else
    th.ox=0;
    z.ox=[];
   end
   fia=find(flsv(:,2)==-1);
   if igg==1
    fii=find(imag(n(:,1))<-.5);
    fiu=fia(1)-2;
   elseif igg==-1
    fii=fia(end)+2;
    fiu=find(imag(n(:,1))<-.5)-1;
   end
  else
   fia=find(flsv(:,2)==-1);
   if igg==1
    fii=find(imag(n(:,1))<-.5);
    fiu=fia(1)-2;
%    puf=[fiu+1:length(n)];
%    ddu=d(puf).*abs(flsv(puf,2));
%    z.ox=sum(ddu);
    z.ox=0;
   elseif igg==-1
    fii=fia(end)+2;
    fiu=find(imag(n(:,1))<-.5)-1;
%    z.ox=0;
    puf=[2:fii-1];
    ddu=d(puf).*abs(flsv(puf,2));
    z.ox=sum(ddu);
   end
    puf=[fii:fiu];
    ddu=d(puf).*flsv(puf,2);
    th.ox=sum(ddu);
  end

' pus' ,keyboard

   pud=2:length(d)-1;
   ds=d(pud).*abs(flsv(pud,2));
   fiD=find(Dov>.1e18);
   dsD=d(fiD).*abs(flsv(fiD,2));

   z.DOP=sum(dsD);
   z.tot=sum(ds);
   z.igg=igg;
%   ' ma! ossido'
%   keyboard

   fiqw=find(iauto(:,1)==2);
   th.qw=d(fiqw);
   NQW=length(find(flsv(:,2)==-1));

   iprqw=find(flsv(:,2)==-1);
   sn=size(n);
%   if sn(2)>1
%    imet=find(imag(n(:,2))<-.5);
%    if length(imet)==0
%     imet=1;
%    end
%   else
    imet=1;
%   end

   ipqw=imet+1:iprqw(1)-1;
   ds=d(ipqw).*abs(flsv(ipqw,2));
   z.qw=sum(ds);

%   ficavp=find(iauto(:,2)==-4);
%   ficav=ficavp(1):ficavp(2);
%   th.cav=sum(d(ficav));

%   iDov=find(Dov(2:end-1)==0)+1;
%   if length(iDov)>2
%    ficav=iDov;
%   else
    ficavp=find(iauto(:,2)==-4);
    ficav=ficavp(1):ficavp(end);
%   end
   th.cav=sum(d(ficav));

   ficavp=find(flsv(:,2)==-1);
%   ficav=ficavp(1)-1:ficavp(end)+1;
   ficav=ficavp(1):ficavp(end);
   th.qwb=sum(d(ficav));


   iDov=find(Dov<0);
%   if length(iDov)>2
%    ipqw=1:iDov(1)-1;
%   else
    ficavp=find(iauto(:,2)==-4);
    ipqw=1:ficavp(end)+1;
 %
 %   ipqw=1:ficavp(end)+1;
 %   ipqw=imet+1:iprqw(1)-1;
%   end
   if igg~=-1
    ipqw=imet+1:iprqw(1)-1;
   else
    ipqw=imet+1:fiu;
   end
   ds=d(ipqw).*abs(flsv(ipqw,2));
   th.elct=sum(ds);

   fiDH=find(Dov>9e18 & d>0);
   if length(fiDH)>0
    ds=d(fiDH).*abs(flsv(fiDH,2));
   else
    ds=0
    fiDH=1;
   end
   th.DH=sum(ds);

   fiD=find(Dov>0 & d>0);
   ds=d(fiD).*abs(flsv(fiD,2));
   th.D=sum(ds);


   fiA=find(Dov<0 & d>0);
   ds=d(fiA).*abs(flsv(fiA,2));
   th.A=sum(ds);
   fiI=find(abs(Dov)<=1e17 & n(:,1)>1.2 & d>0);
   ds=d(fiI).*abs(flsv(fiI,2));
   th.I=sum(ds);



   mass0=9.1e-31;
   mi=pi*4e-7;
   eps0=8.8541e-12;
   c=1/sqrt(mi*eps0);
   Z0=120*pi;
   q=1.6e-19;
   lambda_m=lambdaNot;
   calfa=1e10*Z0*q^3*lambda_m^3/(mass0*2*pi*c)^2/(4*pi);
   cenne=1e6*q^2*lambda_m^2/(2*eps0*mass0*(2*pi*c)^2);

   fip=find(xw(:,1)>0);
   xcp=xw(fip);
   [mued,muhd,med,mhd,tcd]=mmtalga(xcp,300);
   mue=zeros(size(d));
   muh=mue;
   me=mue;
   mh=mue;
   tc=mue;
   mue(fip)=mued;
   muh(fip)=muhd;
   me(fip)=med;
   mh(fip)=mhd;
   tc(fip)=tcd;


   palfae=calfa*abs(Dov(fiD))./(mue(fiD).*me(fiD).^2.*real(n(fiD,1)));
   pennee=cenne*abs(Dov(fiD))./(me(fiD).*real(n(fiD,1)));
   sige=100*q*mue.*abs(Dov)./(1+sqrt(abs(Dov*1e-17)));
   sigh=100*q*muh.*abs(Dov)./(1+sqrt(abs(Dov*1e-17)));
   ds=d.*abs(flsv(:,2));
   dsh=d.*abs(flsv(:,2))./sigh;
   dse=d.*abs(flsv(:,2))./sige;
   dshs=d.*abs(flsv(:,2)).*sigh;
   dses=d.*abs(flsv(:,2)).*sige;


   sighs=sum(ds(fiD))./sum(dsh(fiD)) ;  %1/(Ohm*m)
   siges=sum(ds(fiA))./sum(dse(fiA)) ;  %1/(Ohm*m)

   sighp=sum(dshs(fiD))./sum(ds(fiD));  %1/(Ohm*m)
   sigep=sum(dses(fiA))./sum(ds(fiA));  %1/(Ohm*m)

   Sig.Mup=sighp;
   Sig.Mus=sighs;

   Sig.Mdp=sigep;
   Sig.Mds=siges;
   Sig.buf=sige(end);
   Sig.DH=mean(sige(fiDH));

' LAY', keyboard

   dse=d.*abs(flsv(:,2))./sige;
   dshs=d.*abs(flsv(:,2)).*sigh;
   dses=d.*abs(flsv(:,2)).*sige;

   sighs=sum(ds(fiD))./sum(dsh(fiD)) ;  %1/(Ohm*m)
   siges=sum(ds(fiA))./sum(dse(fiA)) ;  %1/(Ohm*m)

   sighp=sum(dshs(fiD))./sum(ds(fiD));  %1/(Ohm*m)
   sigep=sum(dses(fiA))./sum(ds(fiA));  %1/(Ohm*m)

   Sig.Mup=sighp;
   Sig.Mus=sighs;

   Sig.Mdp=sigep;
   Sig.Mds=siges;
   Sig.buf=sige(end);
   Sig.DH=mean(sige(fiDH));

%fiA=find(dglo.Dop>0 & imatw>0);
%palfae=calfa*abs(dglo.Dop(fiA))./...
%(dglo.mue(fiA).*dglo.me(fiA).^2.*real(nv(fiA,1)))/(1+fre)^3;
%pennee=cenne*abs(dglo.Dop(fiA))./(dglo.me(fiA).*real(nv(fiA,1)))/(1+fre)^3;
%nv(fiA,1)=nv(fiA,1)+pennee-j*palfae;
%
%fiA=find(dglo.Dop<0);
%palfah=calfa*abs(dglo.Dop(fiA))./...
%(dglo.muh(fiA).*dglo.mh(fiA).^2.*real(nv(fiA,1)))/(1+fre)^3;
%penneh=cenne*abs(dglo.Dop(fiA))./(dglo.mh(fiA).*real(nv(fiA,1)))/(1+fre)^3;
%
%nv(fiA,1)=nv(fiA,1)+penneh-j*palfah;


%   N=logspace(15,19,100);
%   mu=1./(1+sqrt(N*1e-17));
%   figure, semilogx(N,mu)


  ' th.elct=sum(ds); '
  keyboard
   th.elct=z.DOP;



   if sn(2)>=2
    fiox=find(n(:,2)>1 & n(:,2)<2);
   else
    fiox=[];
   end
   if length(fiox)==0
    fiox=ficavp(1)
   end
   ipdif=fiox:iprqw(1)-1;

   ico=ipdif(1);
   nrif=n(ico,1);

   while nrif>n(ico-1,1)
    nrif=n(ico-1,1);
    ico=ico-1;
   end
   ipdif=[ico:ipdif(end)];


   ds=d(ipdif).*abs(flsv(ipdif,2));
   th.diff=sum(ds);

   if th.mesa<0
    th.mesa=z.ox;
   end
%   'lay_ct elect', keyboard


  if sn(2)>1
   fimed=find(imag(n(:,2))<-.5);
   if length(fimed)>=1
    fimet=fimed(1);
    th.met=d(fimet);
   else
    th.met=0.1;
   end
  else
    th.met=0.1;
  end

   eval(['save ' nomeloFEM]);
   ilo=-2;
 else
  eval(['load ' nomeloFEM]);
 end

%' dopo', keyboard

% xma=ro.mesa             +
% xox=ro.ox               +
% xcui=ro.met             +

% doxi=th.ox              +
% dmesa=th.mesa           +
% dqw=th.qw               +
% thick_met=th.met
% thick_cav=th.cav        +
% thick_bu=th.bu          +

% zox=z.ox                +
% zqw=z.qw                +


%' Lay_ct ', keyboard



