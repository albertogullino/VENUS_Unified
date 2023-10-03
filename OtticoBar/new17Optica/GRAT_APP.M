%' qui gratalpp', keyboard
     if igr_app<=1
      return
     end

     if inume==0
      if iLP==1
       Kmat_LP
      else
      tic
       Kmat_ve
      toc
      end
     else
      if iLP==1
       Kmat_LP2
      else
      tic
       Kmat_ve2
      toc
      end
     end

if ifp==-10
'fine qui gratalpp', keyboard
end

if length(find(aytot(:,shag)==rgra))==0   %aggingo parte circolare solo se non c'e` gia`

sM=size(MKosp);
if sM(2)>=1
 sM=size(MKosp{1});
 if length(sM)>2
  pu=sM(3)+1;
 else
  pu=2;
 end
 du=MKosp{shag};
 du(:,:,pu)=Kosp;
 MKosp{shag}=du;
 du=MKosm{shag};
 du(:,:,pu)=Kosm;
 MKosm{shag}=du;
 if iztm==1
  du=MKoszp{shag};
  du(:,:,pu)=Koszp;
  MKoszp{shag}=du;
  du=MKoszm{shag};
  du(:,:,pu)=Koszm;
  MKoszm{shag}=du;
 end
else
 MKosp{shag}=Kosp;
 MKosm{shag}=Kosm;
 if iztm==1
  MKoszp{shag}=Koszp;
  MKoszm{shag}=Koszm;
 end
end
aytot(pu,shag)=rgra;
%axtot(pu,shag)=0;
end