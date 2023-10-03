function u=angwr(f);

ud=unwrap(angle(f));
u=ud;
%return

d=diff(ud);
lf=length(f);
fi=find(d<0);
if length(fi)>0
 for k=fi
  pu=k+1:lf;
  ud(pu)=ud(pu)+2*pi;
 end
end
u=ud;
if mean(u)<0
 u=u+2*pi;
end
