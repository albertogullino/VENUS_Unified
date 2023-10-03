 h1=map(flipud(log10(abs(Es)).'),xvero_tot,-fliplr(ztot));
set(h1,'pos',[200  50   950   650]);
 %ylabel([' long. coord. (micron), Gain = ',num2str(gamm{ix,iy}/vg,4),' 1/cm'])
 ylabel([' long. coord. (\mum)'])
 xlabel('  transv. coord. (\mum)')
 title([tit,'  (log units)'])
%keyboard
     oldlimits=caxis;
     caxis([-4 0])
     aitot=aitot(:,1);

axis equal
 axis([-ama,ama,-ztot((end)),-ztot((1))])
%a(1:2)=[-ama ama];
%     caxis([-4 0]-1);
%     caxis([(oldlimits(2)-10) oldlimits(2)-5]);
%     axis([ztot((1)),ztot((end)),0,max(xvero)])
    colorbar;

 
% map(log10(abs(Es(:,fiz))),ztot(fiz),xvero)
% axis([ztot((1)),ztot((end)),0,max(xvero)])
 hold on,
axis normal
if ischema==1
add_schema_Mayer
end
 pausak
