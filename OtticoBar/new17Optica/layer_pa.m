function [nref,aref,n,d,xm,radii,flsv,perd,anyf,iauto,...
         dw,xw,fsw,dov,shav,ipar,ifield]= Layer_pa(lambdaNot,fileName)

kNoll=2*pi/lambdaNot;
radii.array=0;
count=1;
countw=1;
firstWord=[];

fid=fopen(fileName,'r');
NL=0;
nllo=0;
iautoci=0;
setarray=0;
while feof(fid)==0
dovdu=0;
clear  radir radirb radirc ram
   typelay=[];
   while isempty(typelay)

      Nline=fgetl(fid);
      [typelay,remainder]=strtok(Nline);
      NL=NL+1;
      if nllo~=0 & isempty(typelay)
%       [NL nllo]
%       pausak
       puv=ipmem(find(ipmem~=0));
       fls(puv,2)=nllo;
       clear ipmem
       nllo=0;
      end

   end

    fan=findstr(typelay,'%');
 if length(typelay)<8 & length(fan)==0
%   typelay
%   pausak
   iany=0;
   iautoc=0;
   itip=0;
  switch typelay

  case 'Centers'

   [tycen,remainder]=strtok(remainder);
   labl=str2num(tycen);
   setarray=0;
   ipcen=1;
   clear cemem
   while ~isempty(remainder)
    [tycen,remainder]=strtok(remainder);
    [xy]=strtok(tycen,'()');
    [bc,cc]=strtok(xy,',');
    cc=cc(2:end);
    x=str2num(bc);
    y=str2num(cc);
    if x==0 & y==0
     cen=0;
    elseif y==0
     cen=x;
    elseif x==0
     cen=j*y;
    else
     if x<0
      suf=pi;
     else
      suf=0;
     end
     cen=sqrt(x^2+y^2)*exp(j*(atan(y/x)+suf));
    end
    cemem(ipcen,1)=cen;
    ipcen=ipcen+1;
   end

    radii.array(counces,1:length(cemem)+2)=[labl; ipcen-1; cemem];
  otherwise

   if setarray==1
      disp(' forgotten array centers ')
      keyboard
   end

%   if length(typelay)>2
%    fan=typelay(3);
%    if fan=='a'
%     iany=1;
%     if length(typelay)==4
%      fan=typelay(4);
%      if fan=='r'
%       iautoc=1;
%      elseif fan=='C'
%       itip=-4;
%      end
%     end
%    elseif fan=='C'
%       itip=-4;
%    elseif fan=='P'
%     if length(typelay)>=4
%      fan=typelay(4);
%      if fan=='o'
%       itip=-3;
%      elseif fan=='m'
%       itip=-2;
%      elseif fan=='i'
%       itip=-1;
%      end
%      if length(typelay)==5
%       iautoc=1;
%      end
%     end
%    elseif fan=='r'
%       iautoc=1;
%    end
%
%    typelay=typelay(1:2);
%   end

   if length(typelay)>2
    fan=findstr(typelay,'a');
    if length(fan)==1
      iany=1;
    end
    fan=findstr(typelay,'r');
    if length(fan)==1
      iautoc=1;
    end
    fan=findstr(typelay,'Po');
    if length(fan)==1
      itip=-3;
    end
    fan=findstr(typelay,'Pm');
    if length(fan)==1
      itip=-2;
    end
    fan=findstr(typelay,'Pi');
    if length(fan)==1
      itip=-1;
    end
    fan=findstr(typelay,'C');
    if length(fan)==1
      if itip<0
       itip=-itip;
      else
       itip=-4;
      end
    end

    fan=findstr(typelay,'Lc');
    if length(fan)==1
     typelay='Lc';
    else
     typelay='Lg';
    end

   end

   if iautoc>0
    iautoci=iautoci+iautoc;
   end

  end  %switch

 end


 switch typelay

  case 'Reference'
   [type,remainder]=strtok(remainder);

     switch type

      case 'AlGaAs'

        [secondWord,remainder]=strtok(remainder);
        AlCont=str2num(secondWord);
        if AlCont==0
         AlCont=1e-8;
        end

        nref=nAlGaAs(lambdaNot,AlCont);


      case 'Ref'

        [secondWord,remainder]=strtok(remainder);
        nref=str2num(secondWord);

     end

     [Word,remainder]=strtok(remainder);
     aref=str2num(Word);

  case 'Lc'

   ictP=0;


   [word2,remainder]=strtok(remainder);

    if isequal(word2(1),'F')
     if length(word2)==2
      if isequal(word2(2),'s')
       ifield(count)=-1;
      else
       ifield(count)=-2;
      end
     else
       ifield(count)=1;
     end
     [word2,remainder]=strtok(remainder);
    else
     ifield(count)=0;
    end


   if length(word2)>5
    dop=fix(str2num(word2(4:end)));
    if isequal(word2(1:2),'Na')
      dov(countw)=dop;
    elseif isequal(word2(1:2),'Nd')
      dov(countw)=-dop;
    end
%   word2, pausak
    [word2,remainder]=strtok(remainder);
   end

   nlong=fix(str2num(word2));
   if abs(nlong)>1
    nllo=nllo+1;
%    nllo
    ipmem(count)=count;
   end

   fls(count,1:2)=[nlong 0];
   flw(countw)=[nlong];

   [word2,remainder]=strtok(remainder);

   pubar=find(word2=='|');
   if length(pubar)==1
    ipar(count,1,1)=str2num(word2(pubar+2:end));
    ipar(count,2,1)=1;
    ictP=ictP+1;
    word2=word2(1:pubar-1);
   end

   thick=(str2num(word2));
   d(count,1)=thick;
   dw(countw,1)=thick;
   anyf(count,1)=iany;
   if iautoc>0
    iauto(count,1)=iautoci;
    iauto(count,2)=itip;
   else
    iauto(count,2)=itip;
   end

   [type,remainder]=strtok(remainder);


   ictr=0;
   radii.a(count,1)=0;
   radii.b(count,1)=0;
   radii.c(count,1)=0;
%   shav(count,1:4)='plan';
   shav(count,1)=0;

   while isempty(type)==0

   if length(type)==1
    ifty=type=='%';
    if ifty==1
     break
    end
   end

   ictr=ictr+1;

     switch type

     case 'AlGaAs'

        [secondWord,remainder]=strtok(remainder);
        AlCont=str2num(secondWord);
        if AlCont==0
         AlCont=1e-8;
        end

        [fourthWord,remainder]=strtok(remainder);

        Absorp=str2num(fourthWord);

        n(count,ictr)=nAlGaAs(lambdaNot,AlCont)-i*Absorp*100/(2*kNoll);
        xm(count,ictr)=AlCont;
        xw(countw,ictr)=AlCont;


     case 'InGaAs'

        [secondWord,remainder]=strtok(remainder);
        AlCont=str2num(secondWord);
        if AlCont==0
         AlCont=1e-8;
        end

        [fourthWord,remainder]=strtok(remainder);

        Absorp=str2num(fourthWord);

        n(count,ictr)=nInGaAs980(lambdaNot,AlCont)-i*Absorp*100/(2*kNoll);
        xm(count,ictr)=AlCont;
        xw(countw,ictr)=AlCont;


     case 'Ref'

        [secondWord,remainder]=strtok(remainder);
        word2=secondWord;
        pubar=find(word2=='|');
        if length(pubar)==1
         ictP=ictP+1;
         ipar(count,1,ictP)=str2num(word2(pubar+2:end));
         ipar(count,2,ictP)=5;
         ipar(count,3,ictP)=ictr;
         secondWord=word2(1:pubar-1);
        end

        ni=str2num(secondWord);
        [secondWord,remainder]=strtok(remainder);
        Absorp=str2num(secondWord);
        n(count,ictr)=ni-i*Absorp*100/(2*kNoll);
        xm(count,ictr)=-1;
        xw(countw,ictr)=-1;

     case 'Au'

        [nIRe,keIm]=nAu(lambdaNot);
        n(count,ictr)=nIRe-i*keIm;
        xm(count,ictr)=-1;
        xw(countw,ictr)=-1;

     case 'Cr'

        [nIRe,keIm]=nCr(lambdaNot);
        n(count,ictr)=nIRe-i*keIm;
        xm(count,ictr)=-1;
        xw(countw,ictr)=-1;

     case 'Ti'

        [nIRe,keIm]=nTi(lambdaNot);
        n(count,ictr)=nIRe-i*keIm;
        xm(count,ictr)=-1;
        xw(countw,ictr)=-1;

     case 'Pt'

        [nIRe,keIm]=nPt(lambdaNot);
        n(count,ictr)=nIRe-i*keIm;
        xm(count,ictr)=-1;
        xw(countw,ictr)=-1;

     otherwise
%       keyboard
       ictr=ictr-1;
       word2=type;
       pubar=find(word2=='|');
       if length(pubar)==1
        ictP=ictP+1;
        ipar(count,1,ictP)=str2num(word2(pubar+2:end));
        ipar(count,2,ictP)=2;
        ipar(count,3,ictP)=ictr;
        type=word2(1:pubar-1);
       end
       radii.a(count,ictr)=str2num(type);
       [Word,remainder]=strtok(remainder);
%       shav(count,1:length(Word))=Word;
%Word
%keyboard

       fan=findstr(Word,'|P');
       if length(fan)>0
        ictP=ictP+1;
        pubar=find(Word=='|');
        ipar(count,1,ictP)=str2num(Word(pubar+2:end));
        ipar(count,2,ictP)=7;
        ipar(count,3,ictP)=ictr;
        Word=Word(1:fan-1);
       end

       switch Word
       case {'ellipse','rectangle','square'}
        [Word1,remainder]=strtok(remainder);

         word2=Word1;
         pubar=find(word2=='|');
         if length(pubar)==1
          ictP=ictP+1;
          ipar(count,1,ictP)=str2num(word2(pubar+2:end));
          ipar(count,2,ictP)=3;
          ipar(count,3,ictP)=ictr;
          Word1=word2(1:pubar-1);
         end

        radii.b(count,ictr)=str2num(Word1);
        switch Word
         case 'ellipse'
         shav(count,ictr)=3;
         case {'rectangle','square'}
         shav(count,ictr)=2;
        end
       case 'rhombus'
        shav(count,ictr)=4;
        [Word1,remainder]=strtok(remainder);
         word2=Word1;
         pubar=find(word2=='|');
         if length(pubar)==1
          ictP=ictP+1;
          ipar(count,1,ictP)=str2num(word2(pubar+2:end));
          ipar(count,2,ictP)=3;
          ipar(count,3,ictP)=ictr;
          Word1=word2(1:pubar-1);
         end

        radii.b(count,ictr)=str2num(Word1);
        [Word2,remainder]=strtok(remainder);

         word2=Word2;
         pubar=find(word2=='|');
         if length(pubar)==1
          ictP=ictP+1;
          ipar(count,1,ictP)=str2num(word2(pubar+2:end));
          ipar(count,2,ictP)=4;
          ipar(count,3,ictP)=ictr;
          Word2=word2(1:pubar-1);
         end

        radii.c(count,ictr)=str2num(Word2);
       case {'array'}
        setarray=1;
        [Word1,remainder]=strtok(remainder);

         word2=Word1;
         pubar=find(word2=='|');
         if length(pubar)==1
          ictP=ictP+1;
          ipar(count,1,ictP)=str2num(word2(pubar+2:end));
          ipar(count,2,ictP)=6;
          ipar(count,3,ictP)=ictr;
          Word1=word2(1:pubar-1);
         end
         counces=count;
         shav(count,ictr)=5;
       case 'circle'
        shav(count,ictr)=1;
        radii.b(count,ictr)=0;
        radii.c(count,ictr)=0;
       otherwise
        disp(' wrong type of transverse shape ')
        keyboard
       end
     end
   [type,remainder]=strtok(remainder);

   end  % while trasv.
  count=count+1;
  countw=countw+1;

  case 'Lg'       % graded layer
   ictP=0;

   [word2,remainder]=strtok(remainder);

    if isequal(word2(1),'F')
     if length(word2)==2
      if isequal(word2(2),'s')
       ifield(count)=-1;
      else
       ifield(count)=-2;
      end
     else
       ifield(count)=1;
     end
     [word2,remainder]=strtok(remainder);
    else
     ifield(count)=0;
    end

   if length(word2)>5

    dop=fix(str2num(word2(4:end)));
    if isequal(word2(1:2),'Na')
      dovdu=dop;
    elseif isequal(word2(1:2),'Nd')
      dovdu=-dop;
    end
%   word2, pausak
    [word2,remainder]=strtok(remainder);
   end
   nlong=fix(str2num(word2));

%   fls(count,1)=nlong;

   [word2,remainder]=strtok(remainder);

   pubar=find(word2=='|');
   if length(pubar)==1
    ipar(count,1,1)=str2num(word2(pubar+2:end));
    ipar(count,2,1)=1;
    ictP=ictP+1;
    word2=word2(1:pubar-1);
   end


   thick=(str2num(word2));
   radir=0;
   nir=0;
%   d(count,1)=thick;

   [type,remainder]=strtok(remainder);

   ictr=0;


   while isempty(type)==0

   if length(type)==1
    ifty=type=='%';
    if ifty==1
     break
    end
   end

   ictr=ictr+1;

     switch type

     case 'AlGaAs'

        [secondWord,remainder]=strtok(remainder);
        AlCont1=str2num(secondWord);
        if AlCont1==0
         AlCont1=1e-8;
        end

        [secondWord,remainder]=strtok(remainder);
        AlCont2=str2num(secondWord);
        if AlCont2==0
         AlCont2=1e-8;
        end

        [fourthWord,remainder]=strtok(remainder);
        Absorp1=str2num(fourthWord);

        [fourthWord,remainder]=strtok(remainder);
        Absorp2=str2num(fourthWord);

        [fourthWord,remainder]=strtok(remainder);
        ndis=fix(str2num(fourthWord));

%        n(count,ictr)=nAlGaAs(lambdaNot,AlCont)-i*Absorp*100/(2*kNoll);
%        xm(count,ictr)=AlCont;

     case 'Ref'
        [secondWord,remainder]=strtok(remainder);

        word2=secondWord;
        pubar=find(word2=='|');
        if length(pubar)==1
         ictP=ictP+1;
         pucount=count:count+ndis-1;
         ipar(pucount,1,ictP)=str2num(word2(pubar+2:end));
         ipar(pucount,2,ictP)=5;
         ipar(pucount,3,ictP)=ictr;
         secondWord=word2(1:pubar-1);
        end


        ni=str2num(secondWord);
        [secondWord,remainder]=strtok(remainder);
        Absorp=str2num(secondWord);
        nir(ictr)=ni-i*Absorp*100/(2*kNoll);
%        n(count,ictr)=ni;
%        xm(count,ictr)=-1;

     case 'Au'

        [nIRe,keIm]=nAu(lambdaNot);
        ni=nIRe-i*keIm;
        nir(ictr)=ni;
%        n(count,ictr)=nIRe-i*keIm;
%        xm(count,ictr)=-1;

     case 'Cr'

        [nIRe,keIm]=nCr(lambdaNot);
        ni=nIRe-i*keIm;
        nir(ictr)=ni;
%        n(count,ictr)=nIRe-i*keIm;
%        xm(count,ictr)=-1;

     case 'Ti'

        [nIRe,keIm]=nTi(lambdaNot);
        ni=nIRe-i*keIm;
        nir(ictr)=ni;
%        n(count,ictr)=nIRe-i*keIm;
%        xm(count,ictr)=-1;

     case 'Pt'

        [nIRe,keIm]=nPt(lambdaNot);
        ni=nIRe-i*keIm;
        nir(ictr)=ni;
%        n(count,ictr)=nIRe-i*keIm;
%        xm(count,ictr)=-1;

     otherwise

       ictr=ictr-1;
%       radii(count,ictr)=str2num(type);


       word2=type;
       pubar=find(word2=='|');
       if length(pubar)==1
        ictP=ictP+1;
        pucount=count:count+ndis-1;
        ipar(pucount,1,ictP)=str2num(word2(pubar+2:end));
        ipar(pucount,2,ictP)=2;
        ipar(pucount,3,ictP)=ictr;
        type=word2(1:pubar-1);
       end

       radir(ictr)=str2num(type);
       [Word,remainder]=strtok(remainder);

       fan=findstr(Word,'|P');
       if length(fan)>0
        ictP=ictP+1;
        pubar=find(Word=='|');
        pucount=count:count+ndis-1;
        ipar(pucount,1,ictP)=str2num(Word(pubar+2:end));
        ipar(pucount,2,ictP)=7;
        ipar(pucount,3,ictP)=ictr;
        Word=Word(1:fan-1);
       end

       switch Word
       case {'ellipse','rectangle','square'}
        [Word1,remainder]=strtok(remainder);
         word2=Word1;
         pubar=find(word2=='|');
         if length(pubar)==1
          ictP=ictP+1;
          pucount=count:count+ndis-1;
          ipar(pucount,1,ictP)=str2num(word2(pubar+2:end));
          ipar(pucount,2,ictP)=2;
          ipar(pucount,3,ictP)=ictr;
          Word1=word2(1:pubar-1);
         end
        radirb(ictr)=str2num(Word1);
        radirc(ictr)=0;
        switch Word
         case 'ellipse'
         shavr(ictr)=3;
         case {'rectangle','square'}
         shavr(ictr)=2;
        end
       case 'rhombus'
        [Word1,remainder]=strtok(remainder);
         word2=Word1;
         pubar=find(word2=='|');
         if length(pubar)==1
          ictP=ictP+1;
          pucount=count:count+ndis-1;
          ipar(pucount,1,ictP)=str2num(word2(pubar+2:end));
          ipar(pucount,2,ictP)=3;
          ipar(pucount,3,ictP)=ictr;
          Word1=word2(1:pubar-1);
         end
        radirb(ictr)=str2num(Word1);
        [Word2,remainder]=strtok(remainder);

         word2=Word2;
         pubar=find(word2=='|');
         if length(pubar)==1
          ictP=ictP+1;
          pucount=count:count+ndis-1;
          ipar(pucount,1,ictP)=str2num(word2(pubar+2:end));
          ipar(pucount,2,ictP)=4;
          ipar(pucount,3,ictP)=ictr;
          Word2=word2(1:pubar-1);
         end
        radirc(ictr)=str2num(Word2);
        shavr(ictr)=4;

       case {'array'}
        setarray=1;
        [Word1,remainder]=strtok(remainder);

         word2=Word1;
         pubar=find(word2=='|');
         if length(pubar)==1
          ictP=ictP+1;
          pucount=count:count+ndis-1;
          ipar(pucount,1,ictP)=str2num(word2(pubar+2:end));
          ipar(pucount,2,ictP)=6;
          ipar(pucount,3,ictP)=ictr;
          Word1=word2(1:pubar-1);
         end
         counces=pucount;
         shavr(ictr)=5;

       case 'circle'
        radirb(ictr)=0;
        radirc(ictr)=0;
        shavr(ictr)=1;
       otherwise
        disp(' wrong type of transverse shape ')
        keyboard
       end
     end
   [type,remainder]=strtok(remainder);

   end  % while trasv.


   x=[AlCont1 AlCont2];
   P=-100/(2*kNoll)*[Absorp1 Absorp2];
   if P(1)==P(2)
    cxp=[0 P(1)];
   else
    cxp=polyfit(x,P,1);
   end
   [Ld ,nd, xc]=gradpol(thick,AlCont1,AlCont2,lambdaNot,ndis,cxp);


   icloc=0;
    fiqw=find(fls(:,1)<0);
%    disp( ' layer '), keyboard
    for counti=count:count+ndis-1

     if iautoc>0
      if length(fiqw)>0
       ipuco=count+ndis-1;
      else
       ipuco=count;
      end
     else
      ipuco=counti;
     end

     icloc=icloc+1;
     fls(counti,1:2)=[nlong 0];
     if abs(nlong)>1
      ipmem(counti)=counti;
     end
     d(counti,1)=thick/ndis;
     anyf(counti,1)=iany;
     if iautoc>0
      iauto(ipuco,1)=iautoci;
      iauto(ipuco,2)=itip;
     else
      iauto(ipuco,2)=itip;
     end
    end


   if abs(nlong)>1
    nllo=nllo+ndis;
   end

    if ictr>1
     ram=ones(ndis,1)*radir;
     sram=ones(ndis,1)*shavr;
     nim=ones(ndis,1)*nir(2:length(nir));
    else
%     sn=size(n);
     nco=1;
     ram=zeros(ndis,nco);
     nim=zeros(ndis,nco);
     radirb=ram;
     radirc=ram;
     sram=ram;
    end
   nit=[nd nim];
   xct=[xc nim*0];

   si=size(nit);
   co=1:si(2);
   if si(2)>1
    cor=1:si(2)-1;
   else
    cor=co;
   end
   n(count:counti,co)=nit;
   xm(count:counti,co)=xct;
   if AlCont1==AlCont2
    xw(countw,co)=AlCont1;
    dw(countw,1)=thick;
    dov(countw)=dovdu;
    flw(countw)=nlong;
    countw=countw+1;
   else
    xw(countw:countw+ndis-1,co)=xct;
    dw(countw:countw+ndis-1,1)=thick/ndis;
    dov(countw:countw+ndis-1)=dovdu;
    flw(countw:countw+ndis-1)=nlong;
    countw=countw+ndis;
   end
   radii.a(count:counti,cor)=ram;
   radii.b(count:counti,cor)=radirb;
   radii.c(count:counti,cor)=radirc;
   shav(count:counti,cor)=sram;
%   keyboard
%   n=[n; nit];
%   xm=[xm; xct];
%   radii=[radii; ram];
   count=counti+1;
 otherwise
%   if nllo~=0
%    keyboard
%    puv=ipmem(find(ipmem~=0));
%    fls(puv,2)=nllo;
%    clear ipmem
%   end
%
%
%   [NL nllo]
%   pausak
%   nllo=0;

 end % switch typelay

end
flsv=fliplr(fls);
fsw=flw';
if exist('dov')
 dov=dov';
else
 dov=zeros(size(d));
end
lf=length(ifield);
if lf~=count-1
 ifield(lf+1:count-1)=0;
end
ifield=ifield';

si=size(radii.array);
lf=si(2);
if lf~=count-1
 radii.array(lf+1:count-1,:)=0;
end


if exist('ipar')
spa=size(ipar);
lf=spa(1);
if lf~=count-1
 for k=1:spa(3)
  ipar(lf+1:count-1,:,k)=0;
 end
end
else
 ipar=zeros(count-1,3,3);
end
perd=-imag(n)*2*kNoll/100;

fclose(fid);
%disp('layer_pa')
%keyboard
