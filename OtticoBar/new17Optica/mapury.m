ha1=figure;
set(ha1,'pos',[593         612        1313         367]);
etot=sum(effintv(2:end,nsol));
s3=size(effintv,1);
for imazver=1:s3
             figure(ha1)
%             subplot(1,3,imazver)
             subplot('position',[(imazver-1)*.3 0.1 .3 .75])                   
             surf(X,Y,EFFy{imazver,nsol}/effma), 
             axis off

%             colorbar('SouthOutside')

             
             
%             surf(X,Y,real(EF(:,:,imazver))/sqrt(effma)), colorbar
                         shading('interp'), view(0,90),
	                 axis square, axis equal, grid,
	                 axis([-1 1 -1 1]*axlim*3/4),
             hc=colorbar;
                          set(gca,'fontsize',14)
             phc=get(hc,'pos');
%      keyboard
             phc=get(hc,'pos');
             
             phc(1)=phc(1)+.02;
             set(hc,'pos',phc)
	                 if imazver==1
	                  tit=[' Total far-field '];
	                 elseif imazver==2
	                  tit=[' Purity = ',num2str(effintv(imazver,nsol)/etot*100,3),' %'];
	                 elseif imazver==3
	                  tit=['	         Other comp. = ',num2str(effintv(imazver,nsol)/etot*100,3),' %'];	                  
	                 end
	                 title(tit,'fontsize',14)
end            
