load e
   eM=IdeOo;
   mapab(eM), pausak
   sia=size(iacc);
   for ki=1:sia(1)
    puu=[];
    fi=find(iacc(ki,:)~=0);
    iaccv=iacc(ki,fi);
    if length(iaccv)>0
     for ipu=iaccv
      for ipui=iaccv
       add=((ipu-1)*2*lKA+ipui-1)*nk1max
       puu=[puu pMc+add];
      end
     end
     puu=sort(puu);
     sim=sim0*length(iaccv);
     xx=reshape(x(puu),sim,sim);
     eM(puu)=expm(xx);
   mapab(eM), pausak
    end
   end
   eMb=expm(x);
   mapab(eM-eMb)
