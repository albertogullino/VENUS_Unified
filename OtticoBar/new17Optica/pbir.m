load pbir

dto=[];
nto=[];
ato=[];
cto=[];
puo=[];

fista=find(abs(Dov)>0);
piv=fista(1);
%piv=1;
ico=1;
'ICI bir', keyboard

fsto=[];
while fw>0
 pd=fim(ico);
 nl=fst(pd,1);
 nlay=abs(fst(pd,2));
% [pd piv nl nlay], pausak
 if piv-pd~=0
  pup0=piv:pd-1;
  if length(find(pup0==pua))==1
   fi0=find(pup0==pua);
   leqw=length(fsto)+fi0;
  end
  ay=anyf(pup0);
  nto=[nto; n1(pup0)];
  cto=[cto; c1(pup0)];
  fsto=[fsto; fst(pup0,2)];
  dto=[dto; dv(pup0)];
  ato=[ato; dv(pup0).*ay];
  puo=[puo; pup0'.*ay];
 end
 pupe=pd:pd+nl-1;
 ay=anyf(pupe);
 for inp=1:nlay
  nto=[nto; n1(pupe)];
  cto=[cto; c1(pupe)];
  fsto=[fsto; fst(pupe,2)];
  dto=[dto; dv(pupe)];
  if (inp==poup & nlay==miup) | (inp==podw & nlay==midw)
   ato=[ato; dv(pupe).*ay];
   puo=[puo; pupe'.*ay];
  else
   ato=[ato; -dv(pupe).*ay];
   puo=[puo; pupe'*0];
  end
 end
 piv=pupe(length(pupe))+1;
 fw=fw-nl;
 ico=ico+nl;
end

'ICI 1 bir', keyboard

if piv<length(dv)+1
 pup0=piv:length(dv);
 nto=[nto; n1(pup0)];
 cto=[cto; c1(pup0)];
 fsto=[fsto; fst(pup0,2)];
 dto=[dto; dv(pup0)];
 ay=anyf(pup0);
 ato=[ato; dv(pup0).*ay];
 puo=[puo; pup0'.*ay];
end

puac=puo;
puae=find(ato>0);


nd=[nto nto];
ccd=[cto cto];
xd=[zeros(size(dto)) dto];
xt=cumsum(reshape(xd',2*length(dto),1));
nt=real(reshape(nd',2*length(dto),1));
ct=(reshape(ccd',2*length(dto),1));

dJ=cumsum([0; dto]);
fJ=([0; fsto]);
nJ=([nto(1); nto]);
xJ=([cto(1); cto]);
lJ=length(dJ);
hz=[dJ(1):dJ(lJ)];
uFunc=hz*0;
uF0=hz*0;
perm=[];
xfit=[];
for il=1:lJ-1
 fiz=find(hz>=dJ(il) & hz<dJ(il+1));
 perm=[perm nJ(il+1)*ones(size(fiz))];
 xfit=[xfit xJ(il+1)*ones(size(fiz))];
 if fJ(il+1)==-1
  uFunc([fiz(1)-1 fiz])=1;
 end
 if fJ(il+1)==-1 & il==leqw
  uF0([fiz(1)-1 fiz])=1;
 end
end
xfit(length(hz))=xfit(length(hz)-1);
perm(length(hz))=perm(length(hz)-1);
uFunc(length(hz))=0;
uF0(length(hz))=0;
relPerm=conj(perm.^2);

'ICI 3 bir', keyboard
