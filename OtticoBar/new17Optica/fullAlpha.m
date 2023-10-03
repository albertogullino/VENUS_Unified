clear 
load dVELM

dfull=zeros(sum(abs(fst(:,2))),1);
ii=1;jj=1;
while ii<length(fst)+1
    if fst(ii,2)>1
        rep=fst(ii,1);
        repS=fst(ii,2);
        dfull(jj:jj+(rep*repS)-1)=repmat(dv(ii:ii+rep-1),repS,1);
        jj=jj+rep*repS;
        ii=ii+rep;
    else
        dfull(jj)=dv(ii);
        ii=ii+1;
        jj=jj+1;
    end
end
