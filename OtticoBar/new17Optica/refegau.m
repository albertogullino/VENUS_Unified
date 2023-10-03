clear Gaf
inot=1;
if ~exist('iplo0')
 iplo0=1+inot;
end 
stres=stre;
nf1s=nf1;
if lref>1
 enner=linspace(-stre/2,stre/2,NPF)+real(nf1);
else
%'qui enner', keyboard
 enner=linspace(-stre*2,stre*2,NPF)+real(nf1);
end 
%enner=linspace(-2*stre,2*stre,NPF)+real(nf1);
stre=diff(enner(1:2));
if icomp==1
 if imag(nf1)==0
  ennei=0;
 else
%  ennei=linspace(-stie,stie,NPF)+imag(nf1);
%  ennei=linspace(-2*stie,2*stie,NPF)+imag(nf1);
  if lref>1
   ennei=linspace(-stie/2,stie/2,NPF)+imag(nf1);
  else
   ennei=linspace(-stie*2,stie*2,NPF)+imag(nf1);
  end  
  stie=diff(ennei(1:2));
 end
 
else 
 ennei=0;
 stie=0;
end


 for kl1=1:length(enner)
 for kl2=1:length(ennei)
 nif=enner(kl1)+i*ennei(kl2);
 [gd,du,tp,du,gdu]=gaemt(0,0,lami,thick,nif,Lpm,npm,npair,r_out,rr,0,lu,nu,r_in);
% [gd,du,tp,du,gdu]=gaemt(0,0,lami,thick,nif,Lpm,npm,npair,r_inint,rr,0,lu,nu,r_in);
  if iga==0
   Ga_p= gd;
  else 
%   fiu=pi-angle(gd)+2*angle(tp);
%   gdu=abs(gd)*exp(j*fiu);
%  [gdu,du,tpi]=gaemms(0,0,lami,thick,nif,Lpm,npm,npair,r_in,rr,0,lu,nu,r_inint);
   Ga_p=gd+tp^2*Gace/(1-gdu*Gace);  
  end   
   Gaf(kl1,kl2)=[Ga_p];
 end
 end
if icomp==1
 Dpe=Gaf-Glae;
else
 Dpe=Gaf-Glae;
end

[du,i1v]=(min(abs(Dpe)));
[du,i2]=(min(du));
ire=i1v(i2);
iie=i2;
nf1=enner(i1v(i2))+i*ennei(i2);
Gaef=Gaf(ire,iie);
Dpee=Dpe;
ennere=enner;
iprec=abs(Dpee(ire,iie));
if iplo==iplo0
 figure, plot(enner,abs(Dpe),enner(ire),abs(Dpee(ire,iie)),'wo'), 
 title(' E zoom ')
 pausak
end