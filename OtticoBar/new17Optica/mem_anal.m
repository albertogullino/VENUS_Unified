maz=input(' min size [Mb] = ');
maz=maz*1e6;
clear var By mvn
var=whos;

by={var.bytes};
mv={var.name};
for k=1:length(by), By(k)=by{k}; end
[By,iso]=sort(By);
for k=1:length(by), mvn{k}=mv{iso(k)}; end
M=sum(By)*1.e-6;

disp([' Mem. occupation = ', num2str(M),' [Mb]']),
fi=find(By>maz);
ma=0;
for k=1:length(fi)
 mvn{fi(k)}
 By(fi(k))
 bites=By(fi(k));
 ma=ma+bites;

% pausak
end
disp([' Mem. occupation parziale = ', num2str(ma*1e-6),' [Mb]']),
disp([' Mem. occupation = ', num2str(M),' [Mb]']),
