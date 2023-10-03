clear 
close all
load grat
Np0=Npair;
Np0=-1;

thi=[th(1)/5 th(2)];
nli=[nla(1) -1 nla(3:4)];
thg=0;

'ingresso', keyboard
   [dilua,ailua,nilua,fstua]=lens_add(thi,h,r,N,nli,ifps,Np0);
   
   fia=find(fstua(:,2)==1);
   
   csu=cumsum(dilua(fia));
   figure, plot(ailua(fia),csu,'r.'), pausak
   keyboard

Np0=-5;

'ingresso', keyboard
   [dilub,ailub,nilub,fstub]=lens_add(thg,th,h,r,N,nla,ifps,Np0);   
      csu=cumsum(dilub);
   figure, plot(ailub,csu,'g.'), pausak
   
 fidu=find(diff(fstub(:,2))~=0);
 fi=fidu(1)+1:length(fstub);   
 dilu=[dilua(fia); dilub(fi)];
 ailu=[ailua(fia,:); ailub(fi,:)];
 nilu=[nilua(fia,:); nilub(fi,:)];
 fstu=[fstua(fia,:); fstub(fi,:)];

   csu=cumsum(dilu);
   figure, plot(ailu,csu,'c.'), pausak 
   
  