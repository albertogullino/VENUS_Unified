if lref>1
 enner=linspace(-strm/2,strm/2,NPF)+real(nf2);
else
 enner=linspace(-strm*4,strm*4,NPF)+real(nf2);
end 
%enner=linspace(-stre/2,stre/2,NPF)+real(nf1);
%enner=linspace(-strm,strm,NPF)+real(nf2);
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
  stim=diff(ennei(1:2));
 end
else 
 ennei=0;
 stim=0;
end

clear Gaf Taf
 for kl1=1:length(enner)
 for kl2=1:length(ennei)
 nif=enner(kl1)+i*ennei(kl2);
  [Ga_p,du,Tr]=gaemms(0,0,lami,thick,nif,Lpm,npm,npair,r_out,rr,0,lu,nu,r_in);
   Gaf(kl1,kl2)=[Ga_p];
   Taf(kl1,kl2)=Tr;
 end
 end
Dpe=Gaf-Glam;
 if iGT==2
  Dpe=Taf-Tm;
 end 
%Dpe=abs(Gaf)-abs(Glam);

[du,i1v]=(min(abs(Dpe)));
[du,i2]=(min(du));
irm=i1v(i2);
iim=i2;
Dpem=Dpe;
ennerm=enner;
nf2=enner(i1v(i2))+i*ennei(i2);
Gaem=Gaf(irm,iim);
iprec=abs(Dpem(irm,iim));

if iplo==10
 figure, plot(enner,abs(Dpe),enner(irm),abs(Dpem(irm,iim)),'wo'), 
 title(' H zoom ')
 pausak
end
%'H', keyboard