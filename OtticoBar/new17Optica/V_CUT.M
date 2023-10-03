close all
clear
load mes
ipver=1;
pume=1:length(dv);

     xa=[];
     na=[];
     ra=[];
     sa=[];
     dvu=[];
     rvu=[];
     ic=1;
     k=pume(ic);
     if ipver==1
      hc=figure;
      set(hc,'pos',[452   271   560   420])
     end
     while ic<=length(pume)
      k=pume(ic);

      its=fst(k,1);
      if its<=1
       ic=ic+1;
       xd=dv(k);
       nd=nv0(k,:);
       shd=shavet(k,:);
       rd=radii.a(k,:);
       if length(xa)==0
        xs=0;
       else
        xs=xa(end);
       end
       xa=[xa; xs+xd];
       na=[na; nd];
       ra=[ra; rd];
       sa=[sa; shd];
      else
       kv=k:k+its-1;
       ic=ic+its;
       xd=dv(kv);
       nd=nv0(kv,:);
       shd=shavet(kv,:);
       rd=radii.a(kv,:);
       for kkm=1:fst(k,2)
        xs=xa(end);
        xs=xs+cumsum(xd);
        xa=[xa; xs];
        ra=[ra; rd];
       end
     end
     end

 %    dvp=[dv(1:end-1) dv(2:end)];
 %    dvp=reshape(dvp,prod(size(dvp)),1);
 %    figure, plot(dvp,'r.')
     fi=find(ra(:,1)==0);
     ra(fi,1)=20;
     fi=find(shavet(:,1)>4);
     ra(fi,1)=NaN;
     fi=find(ra(:,2)==0);
     ra(fi,2)=NaN;
%     fi=find(shavet(:,1)==6 & dv~=0);
     fi=find(shavet(:,1)==6 );
     if length(fi)>0
     rar=radii.array{fi};
     dgr=dv(fi);
     dgr=400;
     p1=rar{5};
     p2=rar{6};
     p0=p1-p2; p3=rar{7}; x0=0; nvg=nv(fi,:);
     if nvg(1)<nvg(2)
      gr0_1=[0 0 1];
      gr0_2=[1 1 0];
     else
      gr0_1=[1 1 0];
      gr0_2=[0 0 1];
     end
     gr_sha=gr0_1;
     xr_sha=[0 p2/2*[1 1]];
     xr0_1=[0 p0*[1 1]];
     xr0_2=[0 p2*[1 1]];
     x0=p2/2;
     while x0<p3
      gr_sha=[gr_sha gr0_2];
      xr_sha=[xr_sha x0+xr0_2];
      x0=x0+sum(xr0_2(1:end-1));
      if x0>p3
       break
      end
      gr_sha=[gr_sha gr0_1];
      xr_sha=[xr_sha x0+xr0_1];
      x0=x0+sum(xr0_1(1:end-1));
     end
     if x0>p3
      xr_sha(end+[-1 0])=p3;
     end
     figure, plot(xr_sha,gr_sha,'r.-')
     ra(fi,1:length(gr_sha))=xr_sha;
     xa(fi,1:length(gr_sha))=gr_sha*dgr;
     end

     figure, plot(ra,xa,'w')
