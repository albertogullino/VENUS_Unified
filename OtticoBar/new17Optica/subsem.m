pu=1:length(KK);
pa=length(pu);

KKK=[KTe KTe; KTe KTe];
npia=length(KKK)/pa;

for ke=1:npia
 pue=pu+pa*(ke-1);
 for ki=1:npia
  pui=pu+pa*(ki-1);
  memd=mean(mean(abs(KKK(pui,pue)))); 
  if memd>1e-10
   mema(ki,ke)=1; 
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
    psce=[psce pu+pa*(kk-1)];
   end
  FiMAT{kei}=psce;
 end 
 if length(fi)>0
 memas(fi,fi)=0;
 end
end