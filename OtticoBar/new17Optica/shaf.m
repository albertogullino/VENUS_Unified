function sh=shaf(s)

       switch s
        case {'planar'}
         sh=0;
        case {'circle'}
         sh=1;
        case {'ellipse'}
         sh=3;
        case {'rectangle','square'}
         sh=2;
        case {'rhombus'}
         sh=4;
        case {'array'}
         sh=5;
        case {'grating'}
         sh=6;
%        case {'off_axis'}
%         sh=7;
%        case {'ring'}
%         sh=-1;
        otherwise
         disp(' shaf: this shape is not implemented ')
         keyboard
       end
