%ifp=-10;
%ifp=-4;


if dret0==0 | imet>0
 return
end
Gac=0;
%iGac=0;
%iGac=1;
inanom=1;
itetm=3;
%'GIRO', keyboard

 lami=lambda_cen;
 
%if ite==1
% lv=[1 ];
%elseif ite==2
% lv=[1 2];
%else
% lv=[2];
%end 




%n1=3.16;
%n2=1.34;
%nu=n1;

rfu=n2;
r_in=n1;
nr=r_in;
rr=nr;
teta=0;
Nmodi=31;


[a,b]=nau(lambda_cen*1.e-6);
ngold=a-j*b;


 r_in=n1;
 r_out=n2;

% carico
 nh=2.28;
 nl=1.34;
 tsep=.0;
 nsep=nl;
 npairc=3;
 r_outc=ngold;
 luc=.264;
 nuc=nl;


rrg=nh;
npmc=[nl nh];
Lpmc=lami/4./npmc;



r_inc=nsep;

if iGac==1
%'cont', keyboard
 if nc(1)==nc(2)
  fiF=find(iparv==41);
  if length(fiF)==1
   fiF1=find(Fiev==41);
   dc(1)=Pa_v(fiF1);
   if exist('vparei')
   if iPe==41
     dc(1)=vparei;
   end 
   end
   if exist('vparii')
   if iPi==41
     dc(1)=vparii;
   end    
   end
   if exist('vpare1i')
   if iPe1==41
     dc(1)=vpare1i;
   end    
   end
  end
  
 end
   [Gac]=gaemms(0,0,lambda_cen,dc,nc,[],[],0,r_outc,rr,iLP,[],[],r_incc); 
%[Gac]=gaemms(0,0,lami,tsep,nsep,Lpmc,npmc,npairc,r_outc,rrg,0,luc,nuc,r_inc);
%'cont dopo gac', keyboard
end
if ~exist('thick')
 thick=dret/1000;
end

%'prima di giro', keyboard
%[n_paBW,n_veBW,n_pa1,n_ve1]=giro_eq(lambda_cen,thick,period,DC,n1,n2,ifp,r_incc,Gac,inanom)
if iGac==-1
 Gac=-1;
end

[n_paBW,n_veBW,n_pa1,n_ve1,du,SW]=giro_eq(lambda_cen,thick,period,DC,n1,n2,ifp,r_incc,npost,Gac,inanom,itetm,Nmodi)

if iGac>=0
n_pa=real(n_pa1);
n_ve=real(n_ve1);
else
n_pa=real(n_paBW);
n_ve=real(n_veBW);
end
%' GIRO dopo', keyboard

close all