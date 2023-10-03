function [dilu,ailu,nilu,flu]=particle(R,N,cs,dlay,npar,nlay,ifp)

if R==0
 dilu=dlay*1000;
 ailu=[0];
 nilu=[nlay 0];
 flu=[1 0];
 return
end

ds=dlay;

D=R/(N+1);

nD=0:D:R;

nDv=[0 nD(2:end)-D/2];

h=sqrt(R^2-nD.^2);



hh=diff(-h);

dsu=cs+R;
if dsu>ds
  'error particle too high',
  ['Rp Thick']
  [R ds]
  keyboard
end
dup=ds-dsu;

dgu=cs-R;
if dgu<0
  'error particle too low', keyboard
end
dlo=dgu;
 
dilup=[dup hh(1:end-1) 2*hh(end) fliplr(hh(1:end-1))  dlo]';
ailup=[nDv(2:end) fliplr(nDv)  ]';

dilu=1000*[dup hh(1:end-1) 2*hh(end) fliplr(hh(1:end-1))  dlo]';
ailu=[nDv(1:end-1) fliplr(nDv)  ]';
one=ones(size(ailu));
nilu=[one*npar one*nlay];
nilu(1,:)=nlay;
nilu(end,:)=nlay;
flu=[one one*0];

p=sqrt(R^2-h.^2);

if ifp==-10
 x=linspace(0,R,300);
 y=sqrt(R^2-x.^2);

 xu=[x fliplr(x) -x -fliplr(x)];
 yu=[y -fliplr(y) -y fliplr(y)]+cs;
 pp=[p(1:end); [p(2:end) R]];
 hp=[h(1:end); [h(1:end) ]];
 hr=reshape(hp,prod(size(hp)),1);
 pr=reshape(pp,prod(size(hp)),1)-D/2;
 pr(1)=0;

 Xu=[pr' fliplr(pr') -pr' -fliplr(pr')];
 Yu=[hr' -fliplr(hr') -hr' fliplr(hr')]+cs;

 Xu1=[pr' fliplr(pr'); -pr' -fliplr(pr')];
 Yu1=[hr' -fliplr(hr'); hr' -fliplr(hr')]+cs;

shic=0;
 figure, plot(Xu+shic,Yu,'r','linewidth',1.5)
 hold on
 plot(xu+shic,yu,'w--',Xu1+shic,Yu1,'w'), 
 axis equal
 pausak
 figure, plot(p,h,'r.'), pausak
 figure, semilogy(hh,'r.'), pausak

 figure, plot(ailup,ds-cumsum(dilup),'.',xu,yu), 
 axis equal
% axis([0 R 0 ds])
 pausak
end


