%kiniz=.25;
%kfin=.4;
if icalcola==1
%'acc01', keyboard
  iLP1=iLP;
  idis=0;
  if isfield(Ps,'alim_ridotto')==1
   if length(Ps.alim_ridotto)==2
    idis=100;
   end
  end
  if ~exist('ipeak')
   if isfield(Ps,'KKadv')==1
    KKadv=Ps.KKadv;
    idis=20;
   end
  else
   if ipeak==1
    if isfield(Ps,'Nkf')==1
     Nkf=Ps.Nkf;
     DK=Ps.DK;
     KKadv=Ps.KKadv;
     KKo=KK;
     if Nkf>0
          idis=20;
     else
          idis=21;
     end
    end 
   end 
  end 

  

%'qui', keyboard  
  
  
  if idis==20
   k=linspace(kiniz,kfin,npuntik+1);
   dk=diff(k(1:2));
   NK=ceil(DK/dk);
   if NK==0
    NK=1;
   end

   kdu=[];
%   'qui', keyboard
%   keyboard
   kip=1;
   for ik=1:length(KKadv)
    ki=KKadv(ik);
    [du,fik]=min(abs(k-ki));
%    pausak

%    kfid=linspace(k(fik-1),k(fik+1),Nkf+2);
%    kfi=kfid(2:end-1);
    kfi=linspace(k(fik)-DK,k(fik)+DK,Nkf);
    kdu=[kdu k(kip:fik-NK) kfi];
    kip=fik+NK;
%    pausak
   end
   kdu=[kdu k(kip:end)];
   
   kdu=sort(kdu);
   
   k=kdu';  
   
   dkt=[diff(k); 1];
   
   fiva=find(dkt>DK/Nkf);
   k=k(fiva);
   dkt=diff(k);
   
   KK=k(2:end);
   pes00=dkt;
   
   

   nk1max=length(KK);
   if ifp==-10
    Asp=spline(KKo,Ape,KK),
    figure, plot(KKo,Ape,KKo(fiFI),Ape(fiFI),'ro',KK,Asp,'g.'), %pausak    
%    figure, plot(KK,'.'), pausak
%    figure, plot(KK(1:end-1),diff(KK),'.'), 
%     'ver KK', keyboard
%     'ver KK', keyboard
     'ver KK', keyboard
   end

%    'ver KK', keyboard
   if iLP==1
    nni=numodi;
   else
    nni=2*numodi;
   end
   pes0=[];
   for ipe=1:nni
    pes0=[pes0; pes00];
   end
   pes=pes0;
   if ipolar==0
    pes=[pes; pes];
   end

  end

  if idis==21
   k=linspace(kiniz,kfin,npuntik+1);
   dk=diff(k(1:2));
   NK=ceil(DK/dk);
   if NK==0
    NK=1;
   end

   kdu=[];
%   'qui', keyboard
%   keyboard
   kip=1;
   for ik=1:length(KKadv)
    ki=KKadv(ik);
    [du,fik]=min(abs(k-ki));
%    pausak

%    kfid=linspace(k(fik-1),k(fik+1),Nkf+2);
%    kfi=kfid(2:end-1);
    kfi=[];
    kdu=[kdu k(kip:fik-NK) kfi];
    kip=fik+NK;
%    pausak
   end
   kdu=[kdu k(kip:end)];
   
   k=kdu';  
   
   
   KK=k(2:end);
   pes00=dk*ones(size(KK));
   
      if ifp==-10
       figure, plot(KK,'.'), pausak
       figure, plot(pes00,'.'), 
        'ver KK', keyboard
      end
   

   nk1max=length(KK);


%    'ver KK', keyboard
   if iLP==1
    nni=numodi;
   else
    nni=2*numodi;
   end
   pes0=[];
   for ipe=1:nni
    pes0=[pes0; pes00];
   end
   pes=pes0;
   if ipolar==0
    pes=[pes; pes];
   end

  end

  if idis==100
   kval=Ps.alim_ridotto;
   k=linspace(kval(1),kval(2),npuntik+1)';
%   k=linspace(kiniz,sqrt(kfin),npuntik+1)';
%   k1=linspace(kiniz,kfin*.7,fix((npuntik+1)/3))';
%   k2=linspace(kfin*.7,kfin,fix((npuntik+1)/3*2))';
%   k=[k1; k2(2:end)];
%   'qui KK', keyboard  
%  k=k([1:35 37:end]);
   dkt=diff(k);
   KK=k(2:end);

   nk1max=length(KK);
   pes00=dkt(1)*ones(size(KK));

   if iLP==1
    nni=numodi;
   else
    nni=2*numodi;
   end
   pes0=[];
   for ipe=1:nni
    pes0=[pes0; pes00];
   end
   pes=pes0;
   if ipolar==0
    pes=[pes; pes];
   end

  end 


  if idis==0
   k=linspace(kiniz,kfin,npuntik+1)';
%   k=linspace(kiniz,sqrt(kfin),npuntik+1)';
%   k1=linspace(kiniz,kfin*.7,fix((npuntik+1)/3))';
%   k2=linspace(kfin*.7,kfin,fix((npuntik+1)/3*2))';
%   k=[k1; k2(2:end)];
%   'qui KK', keyboard  
%  k=k([1:35 37:end]);
  dkt=diff(k);
   KK=k(2:end);

   %fie=[1:50 60:100];
   if npuntik==1
    KK=.000001;
    %KK=.01;
   %'qui KK', keyboard  
   end 
   %'qui KK', keyboard     
%   KK=k(2:end).^2;
%   KK=k(1:end);
%   KK=k;
%   KK=k(1:end-1);
%   KK(1)=1e-4;
   if exist('alimu')
    k=linspace(kfin,alimu,nk2+1)';
    KK=[KK; k(2:nk2+1)];
    nk1max=length(KK);
   end
%   'kk',keyboard


%   pes00=[dkt(1:2); dkt(1:end-1)];
   %pes00=[dkt(1); dkt];
   pes00=dkt(1)*ones(size(KK));
%   pes00=[dkt; dkt(end)];

   if iLP==1
    nni=numodi;
   else
    nni=2*numodi;
   end
   pes0=[];
   for ipe=1:nni
    pes0=[pes0; pes00];
   end
   pes=pes0;
   if ipolar==0
    pes=[pes; pes];
   end

   igauleg=0;
   if igauleg==1
    [KK,pes]=gauleg(0,eint/a,npuntik);
    KK=KK';
    pes=pes';
   end
  elseif idis==1
   km=eint/a;
   kfit=kmau-kmad;
   nk12=fix(1/3*nk1);
   nk11=fix(2/3*nk1);
   st1=(kfit)/(nk11);
   st2=(km-kfit)/(nk12);
   k=0;
   ic=0;
   while k<=km
   ic=ic+1;
    if k<=kmad | k>kmau
     k=k+st2;
    else
     k=k+st1;
    end
    KK(ic)=k;
   end
   nk1s=length(KK);
   dkt=diff(KK);
   dkt=[dkt; dkt(nk1s-1)];
   pes00=dkt;


   nk1s=nk11+nk12-2;
   k1=linspace(0,1/3*eint/a,nk11)';
   k2=linspace(1/3*eint/a,eint/a,nk12)';
   KK=[k1;k2(2:nk12)];
   dkt=diff(KK);
   dkt=[dkt; dkt(nk1s-1)];
   pes00=dkt;

    if iLP==1
     nni=numodi+1;
    else
     nni=2*(numodi+1);
    end
    pes0=[];
    for ipe=1:nni
     pes0=[pes0; pes00];
    end
    pes=pes0;
  elseif idis==2
   dkt=diff(KK);
   pes00=[dkt; dkt(length(dkt))];
    if iLP==1
     nni=numodi;
    else
     nni=2*numodi;
    end
    pes0=[];
    for ipe=1:nni
     pes0=[pes0; pes00];
    end
    pes=pes0;
  end

  ifaccsemat=0;

  global iKexact
  if length(iKexact)==0
   iKexact=0;
  end


  if iKexact==0
   icalcola=0;
  end
  if ifr>1
   icalcola=0;
  end  
   if ifr==1 & Pf.ipolar~=3
     icalcola=0;
   end  
%  ' iKex ',  
%  ifr
%  Pf.ipolar
%  ipolar
%  keyboard  

else

  ifaccsemat=1;

end


%'prima acc', keyboard
      acc

