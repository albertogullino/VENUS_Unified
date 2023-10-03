load sa
h=50e-2;
 tha=th;
 nla=nl;
 tha=[thgv -th(1:2)];
 nla=[nl(1) -1 nl(2:3) nl(2:3)];
 Np0=-1.5
   [dilua,ailua,nilua,fstua]=lens_add(tha,h,r,N,nla,ifps,Np0);
   figure, plot(ailua,cumsum(dilua),'r.'), pausak