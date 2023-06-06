SEGNO_elettroni=1;  % -1 particelle
Cur_2D=1;  % 0 se non la devo plottare
 Cn=mode.Cn;
 Cn(:,1)=Cn(:,2);
 Cp=mode.Cp;
 Cp(:,1)=Cp(:,2);
 %puc=[8 19];
 %puc=[4 8];
 puc=pu;
 Fcur=Fcur+1;
 hm=Fcur;
 figure(Fcur)
 set(Fcur,'pos',[349         398        1525         525 ])
subplot(121)
 plot(x,1e-12*Cn(puc,:),'o')
         ax = gca;
       ax.ColorOrderIndex = 1
  hold on, plot(x,1e-12*Cp(puc,:),'.')
  legend(num2str(Cpl(puc)',2))
  xlabel('\rho, \mum')
  ylabel('2D Injection, (cm^2 ns)^{-1}')
    
    hold off

J_XN=SEGNO_elettroni*mode.JXn;
J_YN=SEGNO_elettroni*(mode.JYn);
J_XP=mode.JXp;
J_YP=(mode.JYp);

%J_X=mode.JX;
%J_Y=mode.JY;
y=modePlot.y;
[X,Y]=meshgrid(x,y);

SZ=.02

 X0=x;
 Y0=y;
 dx=diff(X0);
 sogx=.3;
 su=0;
 sX=[];
 for k=1:length(dx)
  su=su+dx(k);
  if su>sogx
   su=0;
   sX=[sX k];
  end
 end 
 sX0=sX;
  ifiver=0;
  
  if ifiver==1
   figure, plot(X0,'.')
   hold on
   plot(sX,X0(sX),'ro'), pausak
  end 
 
  dx=diff(Y0);
  sogx=.3;
  su=0;
  sY=[];
  for k=1:length(dx)
   su=su+dx(k);
   if su>sogx
    su=0;
    sY=[sY k];
   end
  end 
  sY0=sY;
  
   X0=x;
   Y0=y;
   dx=diff(X0);
   sogx=SZ;
   su=0;
   sX=[];
   for k=1:length(dx)
    su=su+dx(k);
    if su>sogx
     su=0;
     sX=[sX k];
    end
   end 
   sX1=sX;
    ifiver=0;
    
    if ifiver==1
     figure, plot(X0,'.')
     hold on
     plot(sX,X0(sX),'ro'), pausak
    end 
   
    dx=diff(Y0);
    sogx=SZ;
    su=0;
    sY=[];
    for k=1:length(dx)
     su=su+dx(k);
     if su>sogx
      su=0;
      sY=[sY k];
     end
    end 
  sY1=sY;
  
  
  if ifiver==1
   figure, plot(Y0,'.')
   hold on
   plot(sY,Y0(sY),'ro'), pausak
  end 
 
 
 SX=X(sY0,sX0);
 SY=Y(sY0,sX0);
 SX1=X(sY1,sX1);
 SY1=Y(sY1,sX1);
 
 geom=mode.geom;
 mesh.node=mode.node;
 mesh.triangle=mode.triangle;
 mesh.nn=mode.nn;
 mesh.nt=mode.nt;
 
 xMax=10;
  node='off'; triangle='off'; color='on'; vpath='off'; arrows='off';
 sd=1:geom.nd; scale=0.000002; cmap=[]; grido='off'; cbar='off';
 
 iar=0
 if iar==1
 Jt_x=pdeintrp(node,triangle(1:4,:),J_XN.'); % T on triangles
 Jt_y=pdeintrp(node,triangle(1:4,:),J_YN.'); % T on triangles
 uvetN=[Jt_x;Jt_y];
 Jt_x=pdeintrp(node,triangle(1:4,:),J_XP.'); % T on triangles
 Jt_y=pdeintrp(node,triangle(1:4,:),J_YP.'); % T on triangles
 uvetP=[Jt_x;Jt_y]; 
else
 uvet=[];
end  


 pcor=8;
%for pcor=[8 15 25]
Vcor=pu;
ipcor=0;
for pcor=Vcor
 ipcor=ipcor+1; 
figure(hm)
subplot(121)

 plot(x,1e-12*Cn(pcor,:),'o-',x,1e-12*Cp(pcor,:),'.-')
  legend('Electrons','Holes')
  xlabel('\rho, \mum')
  ylabel('2D Injection, (cm^2 ns)^{-1}') 
  hold off
subplot(122)
 PlotStru
 quiver(SX,SY,squeeze(J_XN(pcor,sY0,sX0)),squeeze(J_YN(pcor,sY0,sX0)),.2),
 ylim([yo-6 yc]), xlim([0 xMax]),  
 pausak
 quiver(SX,SY,squeeze(J_XP(pcor,sY0,sX0)),squeeze(J_YP(pcor,sY0,sX0)),.2),
 hold on

 if SEGNO_elettroni==-1
     title(['Particle Currents @ ',num2str(Cpl(Vcor(ipcor)),2),', mA'])
 else
     title(['Currents @ ',num2str(Cpl(Vcor(ipcor)),2),', mA'])
 end
 


    xcm=x*1e-4;
    xDiv=(xcm+diff(xcm(1:2))/2)*2*pi;
    jnQ_x=sum(squeeze(mode.JnQW(pcor,:,:)),2)*1000;
    jpQ_x=sum(squeeze(mode.JpQW(pcor,:,:)),2)*1000;
    jQn=SEGNO_elettroni*[jnQ_x' zeros(1,length(xcm)-length(xQW))]./xDiv;
    jQp=[jpQ_x' zeros(1,length(xcm)-length(xQW))]./xDiv;
    
   if Cur_2D==1
    yQWmedioE=mode.yMQW{3}*ones(size(sX0))*1e4;
    yQWmedioH=mode.yMQW{1}*ones(size(sX0))*1e4;
    quiver(x(sX0),yQWmedioE,jQn(sX0)*WQW,zeros(size(sX0)),.3,'b')
    quiver(x(sX0),yQWmedioH,jQp(sX0)*WQW,zeros(size(sX0)),.3,'r')
   end
   
    hold off, 
    pausak
% xlim([0 6])
% ylim([114,116])
% axis normal
 
 iZoom=input(' iZoom ');
 if length(iZoom)==0
  iZoom=0;
 else
  iZoom=1;
 end
 if iZoom==1
 figure
 plot_tri_meshPLOT(geom,mesh,uvet,sd,color,triangle,node,grido,cbar,cmap,scale,vpath,arrows)
  hold on
  ze=y(end);
  thvis=.1;

  if isfield(mesh,'yBTJ')
	  ITJ_finder
	  
      yBTJ=ze-mesh.ygrid(iTJ(1))*1e4;
      yo=ze-yBTJ;
      
      thTJ=.03;
      Yc=[yo yo  yo+thTJ yo+thTJ];
  else
      yo=ze-mode.zox;
      Yc=[yc yc  yc+thvis yc+thvis];

  end
  
  xo=mode.rox;
  
  xc=mode.Contact_i;
  yc=ze;
  wc=mode.Contact_e-xc;
  wo=mode.Contact_e-xo;
  RecCont=[xc yc wc thvis];
  Xc=[xc xc+wc xc+wc xc];
  RecOx=[xo yo wo thvis];
  Xo=[xo xo+wo xo+wo xo]; 
  Yo=[yo yo  yo+thvis yo+thvis];
  hhco=patch(Xc,Yc,'y'); set(hhco,'EdgeColor','none');
  hhco=patch(Xqw,Yqw,'c'); set(hhco,'EdgeColor','none'); 
  ylim(yo+[-.8 1]), xlim(xo+[-2 2])
  quiver(SX1,SY1,squeeze(J_XN(pcor,sY1,sX1)),squeeze(J_YN(pcor,sY1,sX1)),.15),
  quiver(SX1,SY1,squeeze(J_XP(pcor,sY1,sX1)),squeeze(J_YP(pcor,sY1,sX1)),.15),
  axis normal
   if Cur_2D==1
    yQWmedioE=mode.yMQW{3}*ones(size(sX1))*1e4;
    yQWmedioH=mode.yMQW{1}*ones(size(sX1))*1e4;
    quiver(x(sX1),yQWmedioE,jQn(sX1)*WQW,zeros(size(sX1)),.3,'g')
    quiver(x(sX1),yQWmedioH,jQp(sX1)*WQW,zeros(size(sX1)),.3,'m')
   end 
    hold off, pausak
  pausak
 end
end 


