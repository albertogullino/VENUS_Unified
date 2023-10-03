%close all
%clear

% inanom = 0, thick in nm, da trasformare in micron

function [n_pa,n_ve,n_pag,n_veg,Ga_in,SW,neff,T,nmo]=giro_eq(lambda,thicki,period,DC,n1,n2,ifp,nin,nout,Gac,inanom,itetm,Nmodi)

%'entro cont', keyboard

if thicki==0
 n_pa=0;
 n_ve=0;
 n_pag=0;
 n_veg=0;
 return
end 
if nargin<=12
 Nmodi=31;
end
if nargin<=10
 itetm=3;
end
%' entro giro_eq', keyboard
%' entro giro_eq', keyboard
if nargin<10
 inanon=0;
end
if inanom==0
 thick=thicki/1000;
else
 thick=thicki;
end

if nargin<8
 ifp=-4;
end

if nargin>=8
  r_in=nin; 
else
  r_in=n1;
end
if nargin>=9
  r_out=nout; 
else
  r_out=n2;
end


if nargin>=10
 iGac=1;
else
 iGac=0;
 Gac=0;
end
iga=iGac;

if Gac==0
 iGac=0;
 iga=0;
 nu=n2;
end

%path(path,'e:\mvcsel\ret_inf')

if ifp==-10
 iplo=1;
% iplo=10;
else
 iplo=0;
end
%iplo=1

r_inint=r_out;
icomp=1;

lami=lambda;
  dv0=period;
  d1=dv0*DC;
  d2=dv0*(1-DC);
  
% ' in giroeq parametro ', keyboard


 teta=0;
 it=1;



 tt=d1+d2;
 t1=d1;
 t2=d2;
 er1=n1.^2;
 er2=n2.^2;
 e_ve=tt*er1*er2./(t2*er1+t1*er2);
 e_pa=(t1*er1+t2*er2)/tt;
 n_ve=sqrt(e_ve);
 n_pa=sqrt(e_pa);
 if length(n_ve)>1
  n_ve=n_ve(2);
  n_pa=n_pa(1);
 end
 nvBV=n_ve;
 npBV=n_pa;
 if Gac==-1
  n_pag=0;
  n_veg=0;
  Ga_in=0;
  SW=0;
%  'quit giro', keyboard
  return
 end
 
Lpm=[];
npm=[];
npair=0;
lu=[];
nu=[];
rr=r_in;



    if ifp==-10
%      ' in giro_eq: verifica parametri', keyboard 
    end  
%      ' in giro_eq: verifica parametri', keyboard 

 if iGac==0
   [Glae,Glam,se,sm,Te,Tm,neff,nmo]=bast_gen(teta,real(r_in),real(r_inint),n1,n2,d1,d2,thick,lami,itetm,Nmodi);
   SW=[se sm];
   T=[Te Tm];
% ' prima di cerca_n0 ',   keyboard
    cerca_n0


 else
%   [Glae,Glam,se,sm,neff]=bast_genc(teta,r_in,r_inint,n1,n2,d1,d2,thick,lami,Gac,itetm,Nmodi);
   [Glae,Glam,se,sm,neff]=bast_genc(teta,real(r_in),real(r_inint),n1,n2,d1,d2,thick,lami,Gac,itetm,Nmodi);
   SW=[se sm];
%   ' qio cpmt', keyboard
    cerca_n0gau
 end
 
Ga_in=[Glae Glam]; 
if length(Gac)>1
 Gace=Gac(1);
 Gacm=Gac(2);
else
 Gace=Gac(1);
 Gacm=Gac(1);
end

    Nef=[nf1 nf2];  
    n_pag=nf1;
    n_veg=nf2;
%  [gd,du,tp,du,gdu]=gaemt(0,0,lami,thick,n_pag,Lpm,npm,npair,r_out,rr,0,lu,nu,r_inint);
  [gd,du,tp,du,gdu]=gaemt(0,0,lami,thick,n_pag,Lpm,npm,npair,r_inint,rr,0,lu,nu,r_in);   
  if iga==0
   Ga_p= gd;
  else 
   Ga_p=gd+tp^2*Gace/(1-gdu*Gace);
  end   
   Gae_ap=Ga_p;
   
%    [gd,du,tp]=gaemms(0,0,lami,thick,n_veg,Lpm,npm,npair,r_out,rr,0,lu,nu,r_inint);
    [gd,du,tp,du,gdu]=gaemt(0,0,lami,thick,n_veg,Lpm,npm,npair,r_inint,rr,0,lu,nu,r_in);   

     if iga==0
      Ga_p= gd;
     else 
      Ga_p=gd+tp^2*Gacm/(1-gdu*Gacm);
     end   
    Gam_ap=Ga_p;
    if ifp==-10
%     ' quop fim', keyboard
    end 
if itetm==1
 Glam=0;
 Gam_ap=0;
end
if itetm==2
 Glae=0;
 Gae_ap=0;
end
%' cont Glae', keyboard
if abs(Glae)-abs(Gae_ap)>1e-1
  'Problemi determinazione N_reticolo in giro_ep TE',
  %keyboard
end
if  abs(Glam)-abs(Gam_ap)>1e-1
  'Problemi determinazione N_reticolo in giro_ep TM', 
  %keyboard
end

if ifp==-11
 'Confronto reticolo equivalente'
'differenza relativa valori veri/app TE ', abs(1-Glae/Gae_ap), pausak
'n_ve ', n_veg, 
'n_ve B&W ', nvBV, 
pausak
'differenza relativa valori veri/app TM ', abs(1-Glam/Gam_ap)
'n_pa ', n_pag, 
'n_ap B&W ', npBV, 
pausak

end
