function [ictP,ipar,PAR]=fipar0(remainder,ictP,count,ipar,sha,par,ictr);

      switch par
       case 'Isha_me='
       nupar=-2;
       case 'Dr_me='
       nupar=-3;
       case 'StDmin='
       nupar=-4;
       case 'Include_me='
       nupar=-5;
       case 'Ra_me='
       nupar=-6;
       case 'n_ext='
       nupar=-7;
      end

     fstr=findstr(remainder,par);
     if length(fstr)>0
      du=strtok(remainder(fstr+length(par):end));
         pubar=find(du=='|');
         if length(pubar)==1
          ictP=ictP+1;
          duco=count-1;
          ipar(duco,1,ictP)=str2num(du(pubar+2:end));
          ipar(duco,2,ictP)=nupar;
%          if duco==3
%           'fipar'
%           keyboard
%          end
%          [duco ictr ictP]
          ipar(duco,3,ictP)=ictr;
          du=du(1:pubar-1);
         end

      if isequal(par,'shape=')
       if sha==5
        PAR=shaf(du);
       else
        PAR=du;
       end
      elseif isequal(par,'Shar=')
       'dentro', pausak
        PAR=du;
      else
       PAR=str2num(du);
      end
     else
      PAR=[];
     end

%' fipar0, count ', count,  pausak

return
     if sha==6
      switch par
       case 'D='
       nupar=-5;
       case 'd='
       nupar=-6;
       case 'Ry='
       nupar=-7;
       case 'Rx='
       nupar=-8;
       case 'shape='
       nupar=-9;
       case 'orientation='
       nupar=-2;
       case 'shift='
       nupar=-3;
       case 'circle='
       nupar=-4;
       case 'Rext='
       nupar=-11;
      end
     elseif sha<0
      switch par
       case 'Isha='
       nupar=-2;
       case 'Dr='
       nupar=-3;
       case 'Nrdis='
       nupar=-4;
       case 'Iappg='
       nupar=-5;
      end
     elseif sha==7
      switch par
       case 'H='
       nupar=-5;
       case 'D='
       nupar=-6;
       case 'Ndis='
       nupar=-7;
       case 'Nlay='
       nupar=-8;
       case 'ud='
       nupar=-9;
       case 'Rel='
       nupar=-10;
      end
     elseif sha==8
      switch par
       case 'D='
       nupar=-5;
       case 'd='
       nupar=-6;
       case 'Ry='
       nupar=-7;
      end
     elseif sha==5
      switch par
       case 'mx='
       nupar=-2;
       case 'my='
       nupar=-3;
       case 'shape='
       nupar=-5;
       case 'Ry='
       nupar=-5;
       case 'Rx='
       nupar=-6;
       case 'Delta='
       nupar=-7;
       case 'spac='
       nupar=-8;
       case 'bchess='
       nupar=-10;
       case 'Shar='
       nupar=-11;
      end
     end

     fstr=findstr(remainder,par);
     if length(fstr)>0
      du=strtok(remainder(fstr+length(par):end));
         pubar=find(du=='|');
         if length(pubar)==1
          ictP=ictP+1;
          duco=count-1;
          ipar(duco,1,ictP)=str2num(du(pubar+2:end));
          ipar(duco,2,ictP)=nupar;
%          if duco==3
%           'fipar'
%           keyboard
%          end
%          [duco ictr ictP]
          ipar(duco,3,ictP)=ictr;
          du=du(1:pubar-1);
         end

      if isequal(par,'shape=')
       if sha==5
        PAR=shaf(du);
       else
        PAR=du;
       end
      elseif isequal(par,'Shar=')
       'dentro', pausak
        PAR=du;
      else
       PAR=str2num(du);
      end
     else
      PAR=[];
     end

%' fipar, count ', count,  pausak
