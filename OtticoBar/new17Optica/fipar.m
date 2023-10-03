function [ictP,ipar,PAR]=fipar(remainder,ictP,count,ipar,sha,par,ictr);

%'sha'
%sha
%keyboard
     if sha==6
      switch par
       case 'D='
       nupar=-5;
       case 'PE='
       nupar=-5;       
       case 'd='
       nupar=-6;
       case 'Ry='
       nupar=-7;
       case 'Rx='
       nupar=-8;
       case 'shape='
       nupar=-9;
       case 'orien.='
       nupar=-2;
       case 'shift='
       nupar=-3;
%       'shift', keyboard
       case 'circle='
       nupar=-4;
       case 'Rext='
       nupar=-11;
       case 'nlay='
       nupar=0;       
       case 'Tilt='
       nupar=-15; 
       case 'NPET='
       nupar=-12;  
       case 'NVER='
       nupar=-13;         
      end
%      'GRATING fiopar', keyboard
     elseif sha<0 & sha~=-8
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
       case 'Npair='
       nupar=-11;
       case 'Rflat='
       nupar=-12;
       case 'Nrel='
       nupar=-13;
       case 'Rm_rel='
       nupar=-14;
       case 'Rax='
       nupar=-15;
       case 'Mis='
       nupar=-16;
       case 'pos='
       nupar=-17;
       case 'thg='
       nupar=-18;
       case 'thga='
       nupar=-19;
       case 'LA='
       nupar=-20;
       case 'd='
       nupar=-21;       
       case 'orien='
       nupar=-22;   
       case 'Base='
       nupar=-23;  
       case 'Alt='
       nupar=-24;  
       case 'Sth='
       nupar=-25;  
      end
     elseif sha==10
%Tilt _pa  ud=0 An=.05|P1 D=17|P8 n_ex=1 Ndis=-5|P3 Nlay=2  Npair=-6|P5
      switch par
       case 'An='
       nupar=-5;
       case 'D='
       nupar=-6;
       case 'Ndis='
       nupar=-7;
       case 'Nlay='
       nupar=-8;
       case 'ud='
       nupar=-9;
       case 'n_ex='
       nupar=-13;
       case 'Npair='
       nupar=-11;   
       case 'Method='
       nupar=-12;   
       case 'Lbuf='
       nupar=-14;          
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
     elseif sha==-8

      switch par
       case 'RoC='
       nupar=-5;
       case 'R='
       nupar=-6;
       case 'DC='
       nupar=-7;
       case 'Pe='
       nupar=-8;   
       case 'Ns='
       nulpar=-9;      
      end 
%      'sha=-8', keyboard
  
     elseif sha==9        %particella
%       radii.arrayd{1}=Dp;
%       radii.arrayd{2}=Cp;
%       radii.arrayd{3}=Ndispa; 
%       radii.arrayd{4}=Refp;
%       radii.arrayd{5}=Th;
%       radii.arrayd{6}=Refl;
%       radii.arrayd{7}=0;
      switch par
       case 'Dp='
       nupar=-1;
       case 'Cp='
       nupar=-2;
       case 'Ndis='
       nupar=-3;
       case 'Refp='
       nupar=-4;
       case 'Thick='
       nupar=-5;
       case 'Ref='
       nupar=-6;
       case 'sha='
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

%' fipar: count ', count,   
%' fipar: par ', par,  pausak 
