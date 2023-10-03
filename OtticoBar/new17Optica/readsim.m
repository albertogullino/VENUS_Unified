%Grid Pos (microns),Ec (eV),Efn (eV),Ev (eV),Efp (eV)
%0.0000e+00,-4.0700e+00,-5.4292e+00,-5.4940e+00,-5.4292e+00
%2.6500e-03,-4.0700e+00,-5.4292e+00,-5.4940e+00,-5.4292e+00
%5.3000e-03,-4.0699e+00,-5.4292e+00,-5.4939e+00,-5.4292e+00
%7.9500e-03,-4.0698e+00,-5.4292e+00,-5.4938e+00,-5.4292e+00
%1.0600e-02,-4.0697e+00,-5.4292e+00,-5.4937e+00,-5.4292e+00

function [out]=readsim(fileName)

count=1;
firstWord=[];

fid=fopen(fileName,'r');
NL=0;
nllo=0;
lv=0;
while feof(fid)==0
      Nline=fgetl(fid);
      NL=NL+1;
      ln=length(Nline);
      if ln<5
       break
      end
      if NL>1
       if abs(lv-ln)>0
        clear pu
        fiv=find(Nline==',');
        lv=ln;
        l=length(fiv);
        piv=1;
        for k=1:l
         pudu=[piv:fiv(k)-1];
         pu(k,1:length(pudu))=pudu;
         piv=fiv(k)+1;
        end
        pudu=[piv:length(Nline)];
        pu(l+1,1:length(pudu))=pudu;
       end
        for k=1:l+1
         pul=pu(k,find(pu(k,:)>0));
         out(NL-1,k)=str2num(Nline(pul));
        end
      end
end

fclose(fid);
