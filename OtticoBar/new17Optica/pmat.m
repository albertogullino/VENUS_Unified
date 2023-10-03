tic
e0=expm(Mo);
toc
tic
e3=expm3(Mo);
toc

mapab(e3-e0)
p=1;

tic
for k=1:10
 p=e0*p;
end
toc

