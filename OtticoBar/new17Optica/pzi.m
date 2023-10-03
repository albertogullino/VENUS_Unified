load sai
 if izer==0
 nzi=[];
 fia=find(diff(azi)~=0);
 azi=azi(fia);
 nzii=nzii(fia);
 
   for ks=1:length(azi)-1
    fi=find(alati>=azi(ks) & alati<azi(ks+1))
    if length(fi)>0
      nad=nzii(ks)*ones(size(fi));
      nzi=[nzi nad];
      ks
      pausak
    end
   end
   keyboard
  if length(nzi)<length(alati)-1
      fiu=find(diff(azid)<0)+1;
   nzi=[nzi nziis(fiu)*ones(1,(length(alati)-1-length(nzi)))];
  end 
 else 
   nzi=ones(1,sl)*nzii(1);
 end
 