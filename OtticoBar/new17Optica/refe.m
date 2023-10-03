clear Gaf
%lref, pausak
if lref>1
 enner=linspace(-stre/2,stre/2,NPF)+real(nf1);
else
 enner=linspace(-stre*2,stre*2,NPF)+real(nf1);
end 
stre=diff(enner(1:2));
if icomp==1
 if imag(nf1)==0
  ennei=0;
 else
%  ennei=linspace(-stie/2,stie/2,NPF)+imag(nf1);
%  ennei=linspace(-stie,stie,NPF)+imag(nf1);
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
%'refe', keyboard
clear Gaf Taf
 for kl1=1:length(enner)
 for kl2=1:length(ennei)
 nif=enner(kl1)+i*ennei(kl2);
  [Ga_p,du,Tr]=gaemms(0,0,lami,thick,nif,Lpm,npm,npair,r_out,rr,0,lu,nu,r_in);
   Gaf(kl1,kl2)=[Ga_p];
   Taf(kl1,kl2)=Tr;
 end
 end
Dpe=Gaf-Glae;
 if iGT==2
  Dpe=Taf-Te;
 end 
%Dpe=abs(Gaf)-abs(Glae);

[du1,i1v]=(min(abs(Dpe)));
[du,i2]=(min(du1));
ire=i1v(i2);
iie=i2;
nf1=enner(i1v(i2))+i*ennei(i2);
Gaef=Gaf(ire,iie);
Dpee=Dpe;
ennere=enner;
iprec=abs(Dpee(ire,iie));
if iplo==10
 figure, plot(enner,abs(Dpe),enner(ire),abs(Dpee(ire,iie)),'wo'), 
 title(' E zoom ')
 pausak
end
%'E', keyboard