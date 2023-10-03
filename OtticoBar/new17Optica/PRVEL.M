v=ones(100,1);
v2=ones(200,1);
M=diag(v);
M2=[M M; M M];

tic
for k=1:100000
 r=M2*v2;
end
toc

tic
for k=1:100000
 r=v.*v+v.*v;
 r=v.*v+v.*v;
end
toc
