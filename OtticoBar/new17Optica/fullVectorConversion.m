% function [Dov_full,dv_full,perd_full,alpha_full]=fullVectorConversion(Dov,dv,perd,alpha,ABS,fst)
perd1=perd(:,1);
Dov_full=zeros(sum(abs(fst(:,2))),1);
dv_full=Dov_full;perd_full=Dov_full;nv0_full=Dov_full;

ii=1;jj=1;
while ii<length(fst)+1
    if fst(ii,2)>1
        rep=fst(ii,1);
        repS=fst(ii,2);
        
        Dov_full(jj:jj+(rep*repS)-1)=repmat(Dov(ii:ii+rep-1),repS,1);
        dv_full(jj:jj+(rep*repS)-1)=repmat(dv(ii:ii+rep-1),repS,1);
        perd_full(jj:jj+(rep*repS)-1)=repmat(perd1(ii:ii+rep-1),repS,1);
        nv0_full(jj:jj+(rep*repS)-1)=repmat(nv0(ii:ii+rep-1),repS,1);
%         x_full(jj:jj+(rep*repS)-1)=repmat(x_short(ii:ii+rep-1),repS,1);
        
        jj=jj+rep*repS;
        ii=ii+rep;
    else
        Dov_full(jj)=Dov(ii);
        dv_full(jj)=dv(ii);
        perd_full(jj)=perd1(ii);
        nv0_full(jj)=nv0(ii);
        
        ii=ii+1;
        jj=jj+1;
    end
end

a_find=find(Dov_full>0);
d_find=find(Dov_full<0);

alpha_full=zeros(size(Dov));
alpha_full(a_find)=FCA.h.*abs(Dov_full(a_find));
alpha_full(d_find)=FCA.e.*abs(Dov_full(d_find));
alpha_full=ABS.f*alpha_full/1e18;
kNoll=2*pi/lambda0;

perd_new=alpha_full*100/(2*kNoll);
nv0_newFull=nv0_full(:,1)-1i*alpha_full*100/(2*kNoll);

% Overwrite hard coded vector
nv0_full(:,1)=nv0_newFull;
perd_full(:,1)=perd_full(:,1)+alpha_full;

figure
hold on
% plot(cumsum(dv_full),perd_full-alpha_full)
plot(cumsum(dv_full),perd_full)
ylabel('log(FCA)')
set(gca,'yscale','log')
set(gca,'FontName','Times New Roman','Fontsize',13)
legend('str','FCA+str')
keyboard