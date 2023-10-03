function [Kte,Ktm,f1o,zio,nzo]=gaz_mtor(fiQ,ni,L_i,n_i,rr,nb,nu,lambda0,fr,kv,iLP,ifp,iff,ibasti)


global isempl  par_grat itetmt icomp_grat icrescita icarico
global KKies KKims

ibast=abs(ibasti);
mu=rr/nu;
mub=rr/nb;
sginu=-1;

f1u=[1; sginu/mu];
f1um=f1u;
f1ue=f1u;

ifer0=0;

if length(isempl)==0
isempl=1;
end
isto=0;

j=sqrt(-1.);
flp=1-iLP;
k0=2*pi/lambda0;

ibasts=ibast;


ic=1;
zi=0;

if iff==1
sginu=-1;
mu=rr/nu;
mub=rr/nb;
nz(1)=n_i(1);
NPXi=20;
f1o=[1; sginu/mu];
zio=0;
nzo=[];
puL=1:length(L_i);
else
f1o=0;
zio=0;
nzo=0;
puL=2:length(L_i)-1;
end

ncom=n_i;
ncom(fiQ)=ncom(fiQ)+j*ni;
%fiQ
%ni
if iff==1
%' ni dentro', keyboard
end
%' ni dentro Orta', pausak

cont_grat=0;
 k=0;

 ic=1;
 Kte=1;
 Ktm=1;
 iferma=0;
for nst=puL
  ic=ic+1;
%  r=n_i(nst);
  r=ncom(nst);
  mi=rr/r;
  dx=L_i(nst)/1000;
%  'qoi', keyboard
if length(find(ibast-nst==0))==1
  iferma=ifer0;
  if iferma==1
   ' qui set ret ', pausak
  end

  if ibasts~=0
    cont_grat=cont_grat+1;
    ibadif=[ibasts(1)-1 ibasts];
    fico=find(diff(ibadif)~=1); 
    if length(fico)==0
     Nstrati=length(ibasts);
     nstu=ibasts(Nstrati)+1;
     ibasts=0;
    else 
     Nstrati=length(ibasts(1:fico-1));
     nstu=ibasts(Nstrati)+1;
     ibasts=ibasts(fico:end);
    end   
  
   if icomp_grat==1
    nsti=nst-1;
    r_out=ncom(nstu);
    r_in=ncom(nsti);
    Mu=[1 1; r_out/rr -r_out/rr]/sqrt(r_out); 
    Mi=0.5*sqrt(r_in)*[1 rr/r_in; 1 -rr/r_in];     
    puretloc=nsti+1:nstu-1;
   
    r=ncom(puretloc);
    dx=L_i(puretloc)/1000; 
    n1=r;
    if strcmp(icrescita,'bottom') 
     n2d=r_out;  
    else 
     n2d=r_in;  
    end     
%    n2=n2d*ones(size(n1)); 

    fimen=find(ibasti<0);
    if length(fimen)>0
     pun=1:fimen-1;
     ilay=puretloc(pun);
     n2=ncom(ilay+1)*ones(size(n1));      
     n2(pun)=n2d;
    else 
     n2=n2d*ones(size(n1));      
    end  

    tetai=0.01;
    itetm=par_grat.itetm;
    px=par_grat.px;
    DC=par_grat.DC;
    NModi=par_grat.NModi;
    d1i=px*DC;
    d2i=px*(1-DC);


%' qui orta f_mul ', keyboard
%' qui primo strato', nst
Strut.t=L_i;
Strut.n=n_i;
Strut.ibast=ibast;
Strut.icarico=icarico;
 
 [Te,Tm,Ge,Gm]=Orta_tracar(tetai,r_in,r_out,n1,n2,d1i,d2i,dx,lambda0,NModi,itetm,Strut);

%  [Te,Tm,Ge,Gm]=Orta_tras(tetai,r_in,r_out,n1,n2,d1i,d2i,dx,lambda0,NModi,itetm);
%  [Te,Tm,Ge,Gm,du,du,du,du,nef]=bast_tras(tetai,r_in,r_out,n1,n2,d1i,d2i,dx,lambda0,NModi,isto,itetm);

    KKie=(Mu*Te*Mi);
    KKim=(Mu*Tm*Mi);
    KKies{cont_grat}=KKie;
    KKims{cont_grat}=KKim;
     if iferma==1
      ' qui trasm controllo Orta', pausak  
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
      KKi=[cos(dte),              -j*sin(dte).*mi/sqz;...
        -j*sin(dte)./mi*sqz,     cos(dte)];    
      'verifica KKi'
      KKi
      'dopo reticolo', keyboard
     end

   
   if iff==1
     dx=sum(dx);
     rm=mean([n1; n2]);
     if isempl==1
%     [F1ue,zio,nzo]=fieldxt(dx,rm,KKie,F1ue,zio,nzo);
%     [F1um]=fieldxt(dx,rm,KKim,F1um);
       if itetmt==1
        KKf=KKie;
       else
        KKf=KKim;
       end
%      ' qui campo controllo Orta', keyboard
      [f1o,zio,nzo]=fieldxt(dx,rm,KKf,f1o,zio,nzo);
%     [F1um]=fieldxt(dx,rm,KKim,F1um);
     else  %isempl
      [dF1ue]=fieldxt(dx,rm,KKie,F1ue);
      [dF1um]=fieldxt(dx,rm,KKim,F1um);
      [F1ue,F1um,zio,nzo]=field_ret(NPXi,F1ue,F1um,zio,nzo,r_in,r_out,n1,n2,d1i,d2i,dx,lambda0,NModi,rr);
      F1ue(:,end)=dF1ue(:,end);
      F1um(:,end)=dF1um(:,end);
     end  %isempl
   end  %iff
  else  %icomp_grat
    KKie=KKies{cont_grat};
    KKim=KKims{cont_grat};
   if iferma==1
    ' KK salvati', nst, pausak
   end
  end %icomp_grat 
  
 else  %ibasts
   KKie=[1 0; 0 1];
   KKim=[1 0; 0 1];

   if iferma==1
    ' KK = Ide', nst, pausak
   end
%' qui orta salto ', keyboard
 end   %ibasts

else %ibast
%'ver',  pausak
  if iff==1
%function [fz,hz,z,nz]=fieldz(L_i,n_i,NPXi,rr,nb,nu,lambda0,fr,k,iLP)
   [f1o,zio,nzo]=fieldxu(dx,r,NPXi,rr,1,1,lambda0,fr,f1o,zio,nzo);
  end
  be=k0*(1+fr)*r;
  dte=be*dx;
  C=cos(dte);
  S=-j*sin(dte);
  KKi=[C,             S.*mi;...
      S./mi,     C];
  KKie=KKi;    
  KKim=KKi;    
 end  %ibast
  Kte=KKie*Kte;
  Ktm=KKim*Ktm;
  f1ue(:,ic)=Kte*f1ue(:,1);
  f1um(:,ic)=Ktm*f1um(:,1);  
   if iferma==1
    nst
    ' FERMA', pausak
    if nst==ibast(end)+2
     iferma=0;
    end
   end 
end  %L_v
 

%disp('gaz_mtor fine'), keyboard

if ifp>-3
% keyboard
end