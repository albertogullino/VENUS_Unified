clear Gaf
inot=1;
%inot=0;
%'inot', keyboard
%if ~exist('iplo0')
 iplo0=1+inot;
%end 
%enner=linspace(-2*strm,2*strm,NPF)+real(nf2);
if lref>1
 enner=linspace(-strm/2,strm/2,NPF)+real(nf2);
else
 enner=linspace(-strm*2,strm*2,NPF)+real(nf2);
end 
strm=diff(enner(1:2));
if icomp==1
 if imag(nf2)==0
  ennei=0;
 else
  if lref>1
   ennei=linspace(-stim/2,stim/2,NPF)+imag(nf2);
  else
   ennei=linspace(-stim*2,stim*2,NPF)+imag(nf2);
  end  
%  ennei=linspace(-2*stim,2*stim,NPF)+imag(nf2);
  stim=diff(ennei(1:2));
 end
else 
 ennei=0;
 stim=0;
end

%'qui TM', keyboard

 for kl1=1:length(enner)
 for kl2=1:length(ennei)
 nif=enner(kl1)+i*ennei(kl2);
 [gd,du,tp,du,gdu]=gaemt(0,0,lami,thick,nif,Lpm,npm,npair,r_inint,rr,0,lu,nu,r_in); 
  if iga==0
   Ga_p= gd;
  else 
%   fiu=pi-angle(gd)+2*angle(tp);
%   gdu=abs(gd)*exp(j*fiu);
%  [gdu,du,tpi]=gaemms(0,0,lami,thick,nif,Lpm,npm,npair,r_in,rr,0,lu,nu,r_inint);
   Ga_p=gd+tp^2*Gacm/(1-gdu*Gacm);
  end   
   Gaf(kl1,kl2)=[Ga_p]; 
 end
 end
Dpe=Gaf-Glam;
%Dpe=abs(Gaf)-abs(Glam);

[du,i1v]=(min(abs(Dpe)));
[du,i2]=(min(du));
irm=i1v(i2);
iim=i2;
Dpem=Dpe;
iprec=abs(Dpem(irm,iim));
ennerm=enner;
nf2=enner(i1v(i2))+i*ennei(i2);
Gaem=Gaf(irm,iim);
if iplo==iplo0
 figure, plot(enner,abs(Dpe),enner(irm),abs(Dpe(irm,iim)),'wo'), 
 title(' H zoom ')
 pausak
%'H', keyboard
end