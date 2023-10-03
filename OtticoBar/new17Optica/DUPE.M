
   N=N/max(N);

   fir=find(rdis<=r_pil);
   %  for ir=1:length(fir)-1
    for ir=1:length(rdis)-1
     fi=find(ro_in>=rdis(ir) & ro_in<rdis(ir+1));

     rodisN(ir)=mean(rdis(ir:ir+1));
     if length(Nref)>1
      if ianti_gui==1
       Ndis_ag(ir)=mean(Nref(fi));
      end
     end
     if length(N)~=1
      Ndis(ir)=mean(N(fi));
     end
    end
    rodisN(ir+1)=rdis(ir+1);
    if length(N)~=1
     Ndis(ir+1)=0;
    end
     if length(Nref)>1
      if ianti_gui==1
       Ndis_ag(ir+1)=0;
      end
     end

         yiN=[];
         if length(Ndis)>1
          xiN=rodisN(2:length(rodisN));
          NN=Ndis/max(Ndis);
          yiN=-diff(NN);
         end
         if ianti_gui==1
          if length(Ndis_ag)>1
           xiN=rodisN(2:length(rodisN));
            NN_ag=Ndis_ag/max(Ndis_ag);
            yiN_ag=-diff(NN_ag);
          end
         end
if ifp==-10
 figure, plot(xiN,1-cumsum(yiN),'.',ro_in,N,'-'), pausak
 if exist('yiN_ag')
  figure, plot(xiN,Nref(1)-cumsum(yiN_ag),'.',ro_in,Nref,'-'), pausak
 end
end

