%load GRA
if fish>0
 Ngra=fish-1;
 if isfield(par_grat,'th')==1
  Lg=par_grat.th'/1000;
 else
  Lg=Litot(Ngra);
 end
 perio=par_grat.per;
 dg1=perio*par_grat.DC;
 dg2=perio*(1-par_grat.DC);
 ng1=par_grat.n1.';
 ng2=par_grat.n2.';
 O=ones(size(Lg));
 ag0=[dg1*O (dg1+dg2)*O];
 ng0=[ng1 ng2];

 ag=[dg1/2*O (dg1/2+dg2)*O];
 ng=ng0;
 rlat=sum(ag(1,:));
 if Amax/perio>500
  AmaxP=100*perio;
 else
  if max(aitot(fish,:))>0
   AmaxP=max(aitot(fish,:));
  else
   AmaxP=Amax;
  end 
 end
 while ag(1,end)<=AmaxP
  ag=[ag ag(1,end)+ag0];
  ng=[ng ng0];
  rlat=rlat+perio;
 end
 nidopo=nitot(fish+1,:);
 finz=find(abs(nidopo)>0);
 nlat=nidopo(finz(end));
 ng=[ng nlat];
 ncol=size(ng,2); 
 nrig=size(Lto,1)-1+length(Lg); 
  LiP=zeros(nrig,1);
  aiP=zeros(nrig,ncol-1);
  niP=zeros(nrig,ncol);
 if Ngra==1
  pug=[1:length(Lg)];
  puv=length(Lg)+1:nrig;
  LiP(pug,:)=Lg;
  aiP(pug,:)=ag;
  niP(pug,:)=ng; 
  LiP(puv,1)=Lto(2:end);
  aiP(puv,1:size(aitot,2))=au(2:end,:);
  niP(puv,1:size(nitot,2))=nu(2:end,:);   
 else
  pug=[1:length(Lg)]+Ngra-1;
  puv1=1:Ngra-1;
  puv2=Ngra+length(Lg):nrig;
  puv1o=1:Ngra-1;
  puv2o=(Ngra+1):length(Lto);  
  LiP(pug,:)=Lg;
  aiP(pug,:)=ag;
  niP(pug,:)=ng; 
  LiP(puv1,1)=Lto(puv1o);
  aiP(puv1,1:size(au,2))=au(puv1o,:);
  niP(puv1,1:size(nitot,2))=nu(puv1o,:); 
  LiP(puv2,1)=Lto(puv2o);
  aiP(puv2,1:size(aitot,2))=au(puv2o,:);
  niP(puv2,1:size(nitot,2))=nu(puv2o,:);    
 end 
else
% n_distr(Lto,au,abs(nu),n_type,coord_num,'v');
 Ngra=0;
 par_grat=0;
 Litotp=Litot;
 aitotp=aitot;
 nitotp=nitot;
 ratot=[radii.a radii.b];
 ratotc=ratot(2:end-1,:);
 shatotc=shavet(2:end-1);
 fi=find(shatotc==11); 
 if length(fi)>0
  aitotp=[aitotp; zeros(size(aitotp))];
  aitotp(fi,2)=ratotc(fi,2);
  nitotp(fi,3)=nitotp(fi,1);
 end 
 
end

ames=max(max(aitot));
fi=find(aiP>ames);
aiP(fi)=ames;
fi=find(imag(niP)<-2);
niP(fi)=j*4;
shif=sum(LiP(1:2));
%n_distrG(Litotp,aitotp,abs(nitotp),n_type,coord_num,'v',Amax);
n_distrG1(LiP,aiP,abs(niP),n_type,coord_num,'v',Amax);
shading flat
 title('VCSEL transverse ref. index ')
 xlabel('radius (um)')
 ylabel('depth (um)')
pausak