if iplan==1
               iama=[];
 return
end
pu=1:length(KK);
pa=length(pu);

KKK=[KAp KAp; KAp KAp];
%npia=length(KKK)/pa;
npia=numodi*4;
if ipolar==0 
 npia=npia*2;
end

lP=length(Pus0);
kei=0;
for ke=1:npia/2
 pue=pu+pa*(ke-1);
 Fi=[];
 for kk=pu 
  fi=find(pue(kk)==Pus0);
  Fi=[Fi fi];
 end
 %if length(Fi)>0
  kei=kei+1;
  PuI{kei}=Fi;
  PuI{kei+npia/2}=Fi+lP;
 %end 
 %PuI{ke+npia}=Fi+2*lP;
 %PuI{ke+npia*3/2}=Fi+3*lP;
end

npia=kei*2;
clear mema
for ke=1:npia
 pue=PuI{ke};
 for ki=1:npia
  pui=PuI{ki};
  memd=mean(mean(abs(KKK(pui,pue)))); 
  if memd>1e-10
   mema(ki,ke)=1;
  else
   mema(ki,ke)=0;   
  end 
 end
end

memas=mema;

kei=0;
for ke=1:npia
 fi=find(memas(ke,:)==1)
 if length(fi)>0
  kei=kei+1;
   psce=[];
   for kk=fi
    psce=[psce PuI{kk}];
   end
  FiMATo{kei}=psce;
 end 
 if length(fi)>0
 memas(fi,fi)=0;
 end
end

%'Fune scemp', keyboard