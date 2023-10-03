function [lambda,rr]=Lay_par(fileName);

  nomeFs=fileName;
  rad=nomeFs(1:end-4);
  nomeloFEM=[rad,'_ELM.mat'];


  das=0;
  das1=0;


   ibd=0;
   ibs=0;
   ibs1=0;
   DR=dir;
   lD=length(DR);
   for kD=3:lD
    dun=getfield(DR(kD),'name');
     if strcmpi(dun,nomeFs)==1
      dad=datenum(getfield(DR(kD),'date'));
      ibd=1;
     end
     if strcmpi(dun,nomeloFEM)==1
      das=datenum(getfield(DR(kD),'date'));
      ibs=1;
     end
     if ibs*ibd*ibs1==1, break, end
   end
   if das-dad<0
    ilo=0;
   else
    ilo=1;
   end
   if das1-dad<0
    iload=0;
   else
    iload=1;
   end

 if ilo==0

  fid=fopen(fileName,'r');
  Nline=fgetl(fid);
  %pausak



   while feof(fid)==0

      Nline=fgetl(fid);
    if length(findstr(Nline,'Wavelength'))==1 | length(findstr(Nline,'Reference'))==1

     [typelay,remainder]=strtok(Nline);
     switch typelay
      case 'Wavelength'
       [type,remainder]=strtok(remainder);
       lam=str2num(type);
       lambdaNot=lam*1e-9;
%       'qi', keyboard
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
         lambda=lam/1000;
         rr=nref;
        break
      end
    end
   end
 else
  eval(['load ',nomeloFEM]);
 end

