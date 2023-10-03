function m=matmul(M,ind)
m=1;
for k=ind
 m=M(:,:,k)*m;
end
