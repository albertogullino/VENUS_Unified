load sa

 eii=(diag(ei));
 reii=real(eii);
 fipu=find(abs(reii)<1e-10);
 reii(fipu)=0;
 imii=imag(eii);
 fipu=find(abs(imii)<1e-10);
 imii(fipu)=0; 
 fig=find(abs(reii)>1e-5);
 fip=find(abs(reii)<1e-5);
 fiat=find(reii<0);
 eiis=eii;
 eii=reii+j*imii;
% fip=find(imag(eii)>0 & real(eii)>0);
 fipd=find(imii(fip)<0);
 fipr=fip(fipd);
 fip=[fiat; fipr];
 eip=j*eii(fip)/k0/dos;
 
  Adupd1=(Adu(:,fip));
  Adupd=abs(real(Adupd1)+imag(Adupd1));
  puK1=1:length(KK)*2;
  puK2=puK1+length(KK)*2;
  Adup=(Adupd(puK1,:)+Adupd(puK2,:));
  
  [du,mad]=max(abs(Adup),[],1);
  %[Asod,isod]=sort(du);
  %mad=madi(isod);
  fiTE=(find(mad<=length(KK)));
  fiTM=(find(mad>length(KK)));
  
  nefTE=flipud(sort(eip(fiTE)));
 nefTM=flipud(sort(eip(fiTM)));