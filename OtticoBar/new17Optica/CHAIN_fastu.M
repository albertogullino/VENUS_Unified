'prima', keyboard
ilaymem=ones(size(Li));
ilayfast=ones(size(Li));
if ifr==1 & pol==pvet(1)
 ilayfast=ilayfast*0;
end
layskip=0*ones(size(Li));
ifplatot=0;
ifplaco=0;
ifpla=0;
ncop=1;
izi=1;
icoustor=0;
 for istr=iloo
 
    DelT_z0=0;
    if igau==4
     KTe=0;
     KTez=0;
     DelT_z0=Tdis(end,izi);
      KTe=reshape(KTemp(:,:,izi),si2);
      if iztm==1
       KTez=reshape(KTempz(:,:,izi),si2);
      else
       KTez=0;
      end
    end
    
     eval(emme)
     
     Tstof{istr}=Oo;
    icmem=0; 
    if length(istfie)>0
     if length(find(istr==istfie))==1
      icmem=icmem+1;
      Tmeduf(:,:,icmem)=Oo;
%'memorizzo metallo', pausak
     end
    end   
     
%     'fine emme', 
%     Tstor(:,:,istr)=Oo;
size(Tstor)
pausak
end
