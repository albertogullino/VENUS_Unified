var=whos;
by={var.bytes};
for k=1:length(by), By(k)=by{k}; end
M=sum(By)*1.e-6;
disp([' Mem. occupation = ', num2str(M),' [Mb]']),
