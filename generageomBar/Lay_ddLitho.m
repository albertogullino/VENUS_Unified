function [ParDD,ParOpt,ParMore,paTJ]= Lay_ddLitho(fileName)



%' lay_tran dd BAR' , keyboard
count_lay=0;
paTJ=[];
iDD=[];
identroNEFG=0;
nlongcnt=0;
iETCHdd=[];
iNEGF=[];
Imesa=0;
labl=-1;
igrad=0;
fattz=1;
icrit=1;
ipar=zeros(1,3,5);
%ipar=0;
count=1;
isavcount=0;
countw=1;
Dov=0;
firstWord=[];
radii.array{1,1}{1}=0;
radii.a=0;
RaDD=0;
radii.b=0;
radii.c=0;
shav=0;

fid=fopen(fileName,'r');
if fid==-1
    'Non existing file  ', keyboard
end
NL=0;
nllo=0;
iautoci=0;
setarray=0;

while feof(fid)==0
    dovdu=0;
    igrad=0;
    clear  radir radirb radirc ram radirar shavr
    typelay=[];
    while isempty(typelay)
        
        Nline=fgetl(fid);
%		pausak
        NL=NL+1;
        ch0='$ DATAFORDD_START';
        ch1='$ DATAFORDD_END';
        if (strcmp(Nline,ch0)~=1 & NL==1)
            'lay_ct.m: file non adatto', keyboard
        end
        
        if (strcmp(Nline,ch0)==1)
            ienab=1;
            %   'set', keyboard
        end
        
        if (strcmp(Nline,ch1)==1)
            ienab=0;
            %   'reset', keyboard
        end
        fidollaro=length(findstr(Nline,'$'));
        %'qui', keyboard
        if ienab==1 & fidollaro >0
            if length(findstr(Nline,'mesa'))==1
                mesa.radius=extrac_n(Nline,'radius=');
                mesa.Dpass=extrac_n(Nline,'Dpass=');
                %      elseif length(findstr(Nline,'Metalring'))==1
                %        contact.radius_i=extrac_n(Nline,'radius_i=');
                %        contact.width=extrac_n(Nline,'width=');
            elseif length(findstr(Nline,'PassivationSection'))==1
                ParMore.r_passiv=extrac_n(Nline,'passiv_size=');
                ParMore.depth_passiv=extrac_n(Nline,'passiv_depth=');
                ParMore.mesh_passiv=extrac_n(Nline,'passiv_mesh=');
                %StrDD.r_passiv=5;
                %StrDD.mesh_passiv=10;
                %$ PassivationSection passiv_size=23 passiv_mesh=5
                %        'qui', keyboard
                
            elseif length(findstr(Nline,'DD_radial_sections'))==1
                
                fimes=findstr(Nline,'Radial_mesh=');
                contact.Nrmesh=str2num(Nline(fimes+12:end));
                
                %         'qui str', keyboard
                %         fimer=findstr(Nline,'radii=');
                ragAddi=extrac_n(Nline,'radii=');
                ragAdde=extrac_n(Nline,'radie=');
                %ParDD.ragAdd=ragAdd;
                ragAdd=[-ragAddi ragAdde];
                
                %        'qui', keyboard
            elseif length(findstr(Nline,'buffer'))==1
                Buffer.th=extrac_n(Nline,'thick=');
                Buffer.thdd=extrac_n(Nline,'thick_dd=');
                Buffer.size=extrac_n(Nline,'size=');
            elseif length(findstr(Nline,'Thermal'))==1
                Bcond0.CondZc=extrac_n(Nline,'cavity=');
                Bcond0.Cond_air=extrac_n(Nline,'air=');
                Bcond0.CondZm=extrac_n(Nline,'dbr=');
                Bcond0.Condpas=extrac_n(Nline,'pass=');
                Bcond0.CondZmet=extrac_n(Nline,'metal=');
                Bcond0.CondZb=extrac_n(Nline,'substrate=');
                
                %	'Bcond', keyboard
                % Added for double n-doped DBRs
            elseif length(findstr(Nline,'mirroreq_n2'))==1
                mireq.n2.P_i=extrac_n(Nline,'P_i=');
                mireq.n2.P_f=extrac_n(Nline,'P_f=');
                mireq.n2.x_eq=extrac_n(Nline,'x_eq=');
                mireq.n2.Dop_eq=extrac_n(Nline,'Dop_eq=');
                mireq.n2.Nmesh=extrac_n(Nline,'Nmesh=');
                mireq.n2.f_any=extrac_n(Nline,'f_any=');
                mireq.n2.Mateq=extrac_n(Nline,'Mateq=');
                % Standard VCSEL
            elseif length(findstr(Nline,'mirroreq_p'))==1
                mireq.p.P_i=extrac_n(Nline,'P_i=');
                mireq.p.P_f=extrac_n(Nline,'P_f=');
                mireq.p.x_eq=extrac_n(Nline,'x_eq=');
                mireq.p.Dop_eq=extrac_n(Nline,'Dop_eq=');
                mireq.p.Nmesh=extrac_n(Nline,'Nmesh=');
                mireq.p.f_any=extrac_n(Nline,'f_any=');
                mireq.p.Mateq=extrac_n(Nline,'Mateq=');
            elseif length(findstr(Nline,'mirroreq_n'))==1
                mireq.n.P_i=extrac_n(Nline,'P_i=');
                mireq.n.P_f=extrac_n(Nline,'P_f=');
                mireq.n.x_eq=extrac_n(Nline,'x_eq=');
                mireq.n.Dop_eq=extrac_n(Nline,'Dop_eq=');
                mireq.n.Nmesh=extrac_n(Nline,'Nmesh=');
                mireq.n.f_any=extrac_n(Nline,'f_any=');
                mireq.n.Mateq=extrac_n(Nline,'Mateq=');
            end
        end  % ienab
        if NL==10
            %       'fine', keyboard
        end
        
        
        % Nline, pausak
        [typelay,remainder]=strtok(Nline);
        Rem0=remainder;
        
        %       [NL nllo]
        %       pausak
        if nllo~=0 & isempty(typelay)
            %       [NL nllo]
            %       pausak
            %       'dentro'
            %       keyboard
            puv=ipmem(find(ipmem~=0));
            fls(puv,2)=nllo;
%            'fls1',       pausak
            clear ipmem
            nllo=0;
        end
        
    end
    
    % typelay
    if isequal(typelay(1),'%')==0 & length(findstr(typelay,'FATTORE_Z'))>0;
        fattz=str2num(remainder);
        %       lambdaNot=lambdaNot*fattz;
        rigaf=count;
        % ' fattore z', keyboard
    end
    if isequal(typelay(1),'%')==0 & length(findstr(typelay,'CRIT'))>0;
        riga_crit(icrit)=count;
        icrit=icrit+1;
        %       remainder='1';
        % ' CRIT', keyboard
    end
    
    if isequal(typelay(1),'%')==0 & length(findstr(typelay,'ARRAY_GENERIC'))>0;
        pubar=find(typelay=='|');
        if length(pubar)==1
            ictP=ictP+1;
            if setarray<=1
                duco=count;
            elseif setarray==2
                duco=count:count+ndis-1;
            end
            duco=duco-1;
            ipar(duco,1,ictP)=str2num(typelay(pubar+2:end));
            ipar(duco,2,ictP)=6;
            ipar(duco,3,ictP)=ictr;
            typelay=typelay(1:pubar-1);
        end
    end
    switch typelay
        % case {'MESA_PARAM'}
        case 'TJ_PARAM'	
         [star,remainder]=strtok(remainder);
         if strcmp(star,'Start')==1		
          Mline=[];
          while isempty(Mline) | isequal(Mline(1),'%')==1
%         ' leggo '
           Mline=fgetl(fid);
%         pausak
          end
%       ' dopo ', keyboard
          remainder=Mline;		


%       ' dopo ', keyboard
      remainder=Mline;

      ictr0=1;
      ictP0=0;
      ipar0=[];
      count0=2;
      sha=1;
 [ictP0,ipar0,Isham]=fipar0(remainder,ictP0,count0,ipar0,sha,'Isha_me=',ictr0);
 [ictP0,ipar0,Drm]=fipar0(remainder,ictP0,count0,ipar0,sha,'Dr_me=',ictr0);
% [ictP0,ipar0,StDm]=fipar0(remainder,ictP0,count0,ipar0,sha,'StDmin=',ictr0);
% [ictP0,ipar0,Inclm]=fipar0(remainder,ictP0,count0,ipar0,sha,'Include_me=',ictr0);
 [ictP0,ipar0,Ram]=fipar0(remainder,ictP0,count0,ipar0,sha,'Ra_me=',ictr0);
 [ictP0,ipar0,next]=fipar0(remainder,ictP0,count0,ipar0,sha,'n_ext=',ictr0);
 [ictP0,ipar0,OxVero]=fipar0(remainder,ictP0,count0,ipar0,sha,'OxVero=',ictr0);
 [ictP0,ipar0,Nlay_ethc]=fipar0(remainder,ictP0,count0,ipar0,sha,'Nlay_ethc=',ictr0);

       paTJ.shape=Isham;
       paTJ.Dm=Drm;
       paTJ.Ram=Ram;
       paTJ.nex=next;
       paTJ.OxVero=OxVero;
       paTJ.Nlay_ethc=Nlay_ethc;
       paTJ.ipar0=ipar0;


	
%          [~,~,OxVero]=fipar0(remainder,0,0,0,0,'OxVero=',0);	

%'qio OxVero', keyboard		
         end % Start	
         if strcmp(star,'End')==1	
%		  'sono END', keyboard
		  identroNEFG=0;
		  
		 end
        case 'MESA_PARAM'
            [star,remainder]=strtok(remainder);
            if strcmp(star,'Start')==1
                lay_mes(1)=count;
                Imesa=1;
                %      ' qui |MESA', keyboard
                
                Mline=[];
                while isempty(Mline) | isequal(Mline(1),'%')==1
                    %         ' leggo '
                    Mline=fgetl(fid);
                    %         pausak
                end
                %       ' dopo ', keyboard
                remainder=Mline;
                
                ictr0=1;
                ictP0=0;
                ipar0=[];
                count0=2;
                sha=1;
                [ictP0,ipar0,Isham]=fipar0(remainder,ictP0,count0,ipar0,sha,'Isha_me=',ictr0);
                [ictP0,ipar0,Drm]=fipar0(remainder,ictP0,count0,ipar0,sha,'Dr_me=',ictr0);
                [ictP0,ipar0,StDm]=fipar0(remainder,ictP0,count0,ipar0,sha,'StDmin=',ictr0);
                [ictP0,ipar0,Inclm]=fipar0(remainder,ictP0,count0,ipar0,sha,'Include_me=',ictr0);
                [ictP0,ipar0,Ram]=fipar0(remainder,ictP0,count0,ipar0,sha,'Ra_me=',ictr0);
                [ictP0,ipar0,next]=fipar0(remainder,ictP0,count0,ipar0,sha,'n_ext=',ictr0);
                
                pames{1}=Isham;
                pames{2}=Drm;
                pames{3}=StDm;
                pames{4}=Inclm;
                pames{5}=Ram;
                pames{6}=next;
                pames{7}=ipar0;
                
                
                %      'ICI MESA lay_tran', keyboard
                
            else
                lay_mes(2)=count;
                pames{8}=lay_mes;
            end
            
            
            %      'ICI ultimo MESA lay_tran', keyboard
            
            
        case {'ARRAY_GENERIC','ARRAY_MATRIX','GRATING_PARAM','HCG_PARAM','LENS_PARAM','DOE_PARAM','GRAD_LAY','PARTICLE'}
            
            if isavcount==count
                ictrar=ictrar+1;
            else
                ictrar=1;
            end
            % ictrar
            % pausak
            
            isavcount=count;
            
            [tyar,remainder]=strtok(remainder);
            fiar=findstr(tyar,'Lable=');
            if length(fiar)==0
                'Lay_tran: manca Lable ', keyboard
            else
                labl=str2num(strtok(tyar(fiar+6:end)));
                ralabl(count-1,ictrar)=labl;
            end
            %   ralabl(count-1,1)=labl;
            %   count-1
            %   '1'
            %   keyboard
            
            typemat=[];
            
            while isequal(typemat,'END_data')~=1
                
                typemat=[];
                while isempty(typemat)
                    Nline=fgetl(fid);
                    [typemat,remainder]=strtok(Nline);
                end
                if isequal(typemat(1),'%')==0
                    pubar=find(typemat=='|');
                    if length(pubar)==1
                        ictP=ictP+1;
                        if setarray<=1
                            duco=count;
                        elseif setarray==2
                            duco=count:count+ndis-1;
                        end
                        duco=duco-1;
                        ipar(duco,2,ictP)=6;
                        ipar(duco,3,ictP)=ictr;
                        typemat=typemat(1:pubar-1);
                    end
                end
                setarray=0;
                iseta=0;
                ipcen=1;
                clear cemem
                %     typemat,     pausak
                
                switch typemat
                    case 'Centers'
                        
                        if length(pubar)==1
                            ipar(duco,2,ictP)=6;
                        end
                        
                        %      [tycen,remainder]=strtok(remainder);
                        %      labl=str2num(tycen);
                        lablar=1;
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
                        
                        for pcounces=counces;
                            %       radii.array{pcounces,9,ictrar}=lablar;
                            %       radii.array{pcounces,1,ictrar}=labl;
                            %       radii.array{pcounces,2,ictrar}=ipcen-1;
                            %       radii.array{pcounces,3,ictrar}=cemem;
                            radii.arrayd{9}=lablar;
                            radii.arrayd{1}=labl;
                            radii.arrayd{2}=ipcen-1;
                            radii.arrayd{3}=cemem;
                        end
                        
                    case 'Param'
                        
                        lablar=2;
                        
                        sha=5;
                        [ictP,ipar,mx]=fipar(remainder,ictP,count,ipar,sha,'mx=',ictrar);
                        
                        [ictP,ipar,my]=fipar(remainder,ictP,count,ipar,sha,'my=',ictrar);
                        
                        [ictP,ipar,spac]=fipar(remainder,ictP,count,ipar,sha,'spac=',ictrar);
                        
                        [ictP,ipar,shape]=fipar(remainder,ictP,count,ipar,sha,'shape=',ictrar);
                        
                        [ictP,ipar,Ry]=fipar(remainder,ictP,count,ipar,sha,'Ry=',ictrar);
                        
                        [ictP,ipar,Rx]=fipar(remainder,ictP,count,ipar,sha,'Rx=',ictrar);
                        
                        [ictP,ipar,Delta]=fipar(remainder,ictP,count,ipar,sha,'Delta=',ictrar);
                        
                        [ictP,ipar,bchess]=fipar(remainder,ictP,count,ipar,sha,'bchess=',ictrar);
                        
                        
                        [ictP,ipar,Shar]=fipar(remainder,ictP,count,ipar,sha,'Shar=',ictrar);
                        
                        
                        %      fstr=findstr(remainder,'mx=');
                        %      mx=str2num(strtok(remainder(fstr+3:end)));
                        %      fstr=findstr(remainder,'my=');
                        %      my=str2num(strtok(remainder(fstr+3:end)));
                        %      fstr=findstr(remainder,'spac=');
                        %      spac=str2num(strtok(remainder(fstr+5:end)));
                        %      fstr=findstr(remainder,'shape=');
                        %      shapec=strtok(remainder(fstr+6:end));
                        %      shape=shaf(shapec);
                        %      fstr=findstr(remainder,'Ry=');
                        %      Ry=str2num(strtok(remainder(fstr+3:end)));
                        %      fstr=findstr(remainder,'Rx=');
                        %      Rx=str2num(strtok(remainder(fstr+3:end)));
                        %      fstr=findstr(remainder,'Delta=');
                        %      Delta=str2num(strtok(remainder(fstr+6:end)));
                        
                        for pcounces=counces;
                            %       radii.array{pcounces,1,ictrar}=labl;
                            %       radii.array{pcounces,2,ictrar}=mx;
                            %       radii.array{pcounces,3,ictrar}=my;
                            %       radii.array{pcounces,4,ictrar}=shape;
                            %       radii.array{pcounces,5,ictrar}=Ry;
                            %       radii.array{pcounces,6,ictrar}=Rx;
                            %       radii.array{pcounces,7,ictrar}=Delta;
                            %       radii.array{pcounces,8,ictrar}=spac;
                            %       radii.array{pcounces,9,ictrar}=lablar;
                            radii.arrayd{1}=labl;
                            radii.arrayd{2}=mx;
                            radii.arrayd{3}=my;
                            radii.arrayd{4}=shape;
                            radii.arrayd{5}=Ry;
                            radii.arrayd{6}=Rx;
                            radii.arrayd{7}=Delta;
                            radii.arrayd{8}=spac;
                            radii.arrayd{9}=lablar;
                            radii.arrayd{10}=bchess;
                            if length(Shar)>0
                                radii.arrayd{11}=Shar;
                            end
                        end
                        
                        
                        %       radii.arrayd{1}=labl;
                        %       radii.arrayd{2}=ori;
                        %       radii.arrayd{3}=shif;
                        %       radii.arrayd{4}=circle;
                        %       radii.arrayd{5}=D;
                        %       radii.arrayd{6}=drea;
                        %       radii.arrayd{7}=Ry;
                        %       radii.arrayd{9}=lablar;
                        
                    case 'Lens'
                        %    'typemat0',  keyboard
                        
                        ictrs=ictr;
                        ictr=labl;
                        sha=7;
                        
                        [typemat1,remainder]=strtok(remainder);
                        %      'typemat1', typemat1, keyboard
                        
                        
                        switch typemat1
                            
                            case '_grat'
                                ipardu=ipar;
                                ictPdu=ictP;
                                coup=count+1;
                                
                                [ictPdu,ipardu,L_pos]=fipar(remainder,ictPdu,coup,ipardu,sha,'pos=',ictr);
                                [ictPdu,ipardu,L_thg]=fipar(remainder,ictPdu,coup,ipardu,sha,'thg=',ictr);
                                [ictPdu,ipardu,L_thga]=fipar(remainder,ictPdu,coup,ipardu,sha,'thga=',ictr);
                                [ictPdu,ipardu,L_LA]=fipar(remainder,ictPdu,coup,ipardu,sha,'LA=',ictr);
                                [ictPdu,ipardu,L_d]=fipar(remainder,ictPdu,coup,ipardu,sha,'d=',ictr);
                                [ictPdu,ipardu,L_orien]=fipar(remainder,ictPdu,coup,ipardu,sha,'orien=',ictr);
                                PAlens.pos=L_pos;
                                PAlens.thg=L_thg;
                                PAlens.thga=L_thga;
                                PAlens.LA=L_LA;
                                PAlens.d=L_d;
                                PAlens.orien=L_orien;
                                gra.thg=L_thg;
                                gra.thga=L_thga;
                                gra.pos=L_pos;
                                gra.d=L_d;
                                gra.LA=L_LA;
                                ipar=ipardu;
                                ictP=ictPdu;
                                %          ' palens laytr', keyboard
                                
                            case '_pa'
                                
                                ipardu=ipar;
                                ictPdu=ictP;
                                coup=count+1;
                                %      'ICI lens _pa prima', pausak
                                [ictPdu,ipardu,UD]=fipar(remainder,ictPdu,coup,ipardu,sha,'ud=',ictr);
                                
                                %      [ictPdudu,ipardu,Ra]=fipar(remainder,ictPdu,coup,ipar,sha,'R=',ictr);
                                [ictPdu,ipardu,ha]=fipar(remainder,ictPdu,coup,ipardu,sha,'H=',ictr);
                                
                                [ictPdu,ipardu,raD]=fipar(remainder,ictPdu,coup,ipardu,sha,'D=',ictr);
                                ra=raD;
                                %      'ra ',keyboard
                                
                                [ictPdu,ipardu,Ndisc]=fipar(remainder,ictPdu,coup,ipardu,sha,'Ndis=',ictr);
                                
                                [ictPdu,ipardu,NlAR]=fipar(remainder,ictPdu,coup,ipardu,sha,'Nlay=',ictr);
                                
                                [ictPdu,ipardu,Rel]=fipar(remainder,ictPdu,coup,ipardu,sha,'Rel=',ictr);
                                if length(Rel)==0
                                    Rel=0;
                                end
                                
                                [ictPdu,ipardu,Nrel]=fipar(remainder,ictPdu,coup,ipardu,sha,'Nrel=',ictr);
                                if length(Nrel)==0
                                    Nrel=1;
                                end
                                
                                [ictPdu,ipardu,Npair]=fipar(remainder,ictPdu,coup,ipardu,sha,'Npair=',ictr);
                                if length(Npair)==0
                                    Npair=0;
                                end
                                
                                [ictPdu,ipardu,Rflat]=fipar(remainder,ictPdu,coup,ipardu,sha,'Rflat=',ictr);
                                if length(Rflat)==0
                                    Rflat=0;
                                end
                                
                                [ictPdu,ipardu,Rm_rel]=fipar(remainder,ictPdu,coup,ipardu,sha,'Rm_rel=',ictr);
                                if length(Rm_rel)==0
                                    Rm_rel=0;
                                end
                                
                                [ictPdu,ipardu,Rax]=fipar(remainder,ictPdu,coup,ipardu,sha,'Rax=',ictr);
                                if length(Rax)==0
                                    Rax=1;
                                end
                                
                                [ictPdu,ipardu,Mis]=fipar(remainder,ictPdu,coup,ipardu,sha,'Mis=',ictr);
                                if length(Mis)==0
                                    Mis=0;
                                end
                                if Mis~=0
                                    Rax=1.001;
                                end
                                
                                [ictPdu,ipardu,Base]=fipar(remainder,ictPdu,coup,ipardu,sha,'Base=',ictr);
                                if length(Base)==0
                                    Base=0;
                                end
                                
                                [ictPdu,ipardu,Alt]=fipar(remainder,ictPdu,coup,ipardu,sha,'Alt=',ictr);
                                if length(Alt)==0
                                    Alt=0;
                                end
                                
                                [ictPdu,ipardu,Sth]=fipar(remainder,ictPdu,coup,ipardu,sha,'Sth=',ictr);
                                if length(Sth)==0
                                    Sth=0;
                                end
                                %      ' laytran Mis', keyboard
                                
                                %      'ICI lens _pa lay_new 1', pausak
                                %      keyboard
                                
                            case '_th'
                                
                                pubar=find(remainder=='|');
                                if length(pubar)==1
                                    ictP=ictP+1;
                                    ipar(count,1,ictP)=str2num(remainder(pubar+2:end));
                                    ipar(count,2,ictP)=-3;
                                    ipar(count,3,ictP)=ictr;
                                    %          ' ictr ', keyboard
                                end
                                clear duth
                                for kpar=1:NlAR
                                    [pach,remainder]=strtok(remainder);
                                    duth(kpar)=str2num(pach);
                                end
                                icunt=1;
                                dutha=[];
                                while length(remainder)>1
                                    [pach,remainder]=strtok(remainder);
                                    dutha(icunt)=str2num(pach);
                                    icunt=icunt+1;
                                end
                                
                                duth
                                %      'ICI lens _th lay_new', pausak
                                %      keyboard
                                sNpa=sign(Npair);
                                Npair=abs(Npair);
                                %      keyboard
                                if sNpa==0
                                    dusa=duth;
                                    %        duth=repmat(dusa,1,fix(Npair));
                                    duth=repmat(dusa,1,1);
                                    if Npair-floor(Npair)>0
                                        duth=[duth duth(1)];
                                    end
                                end
                                if length(dutha)>0
                                    duth=[duth dutha];
                                end
                                radii.array{count,1}{12}=duth;
                                
                                if length(find(duth<0))>0
                                    Nline=fgetl(fid);
                                    [du,Np_du]=strtok(Nline);
                                    Np_ad=-str2num(Np_du);
                                else
                                    Np_ad=[];
                                end
                                ipar=ipardu;
                                ictP=ictPdu;
                                %           'ICI lens _pa lay_new', pausak
                                PAlens.UD=UD;
                                PAlens.H=ha;
                                PAlens.D=raD;
                                PAlens.Ndis=Ndisc;
                                PAlens.NlAR=NlAR;
                                PAlens.Rel=Rel;
                                PAlens.Nrel=Nrel;
                                PAlens.Npair=sNpa*Npair;
                                PAlens.Np_ad=Np_ad;
                                PAlens.Rflat=Rflat;
                                PAlens.Rm_rel=Rm_rel;
                                PAlens.Rax=Rax;
                                PAlens.Base=Base;
                                PAlens.Alt=Alt;
                                PAlens.Sth=Sth;
                                
                                radii.array{count,1}{11}=PAlens;
                                
                            case '_ab'
                                pubar=find(remainder=='|');
                                if length(pubar)==1
                                    ictP=ictP+1;
                                    ipar(count,1,ictP)=str2num(remainder(pubar+2:end));
                                    ipar(count,2,ictP)=-40;
                                    ipar(count,3,ictP)=ictr;
                                end
                                
                                clear duab
                                for kpar=1:NlAR+2
                                    [pach,remainder]=strtok(remainder);
                                    duab(kpar)=str2num(pach);
                                end
                                icunt=1;
                                duaba=[];
                                while length(remainder)>0
                                    [pach,remainder]=strtok(remainder);
                                    duaba(icunt)=str2num(pach);
                                    icunt=icunt+1;
                                end
                                if length(duaba)>0
                                    duab=[duab duaba];
                                end
                                
                            case '_nr'
                                pubar=find(remainder=='|');
                                if length(pubar)==1
                                    ictP=ictP+1;
                                    ipar(count,1,ictP)=str2num(remainder(pubar+2:end));
                                    ipar(count,2,ictP)=-4;
                                    ipar(count,3,ictP)=ictr;
                                end
                                
                                clear dunr
                                for kpar=1:NlAR+2
                                    [pach,remainder]=strtok(remainder);
                                    dunr(kpar)=str2num(pach);
                                end
                                icunt=1;
                                dunra=[];
                                while length(remainder)>0
                                    [pach,remainder]=strtok(remainder);
                                    dunra(icunt)=str2num(pach);
                                    icunt=icunt+1;
                                end
                                if length(dunra)>0
                                    dunr=[dunr dunra];
                                end
                                if exist('duab')
                                    dunre=dunr;
                                    clear dunr
                                    for kpar=1:length(dunre)
                                        dunr(kpar)=dunre(kpar)-j*100*duab(kpar)/(2*kNoll);
                                    end
                                end
                                
                                %        ' qui npari', keyboard
                                %       if Npair>0
                                if sNpa>0
                                    dusa0=dunr;
                                    dusa=dunr(2:1+NlAR);
                                    dunr=repmat(dusa,1,1);
                                    %        dunr=repmat(dusa,1,fix(Npair));
                                    if Npair-floor(Npair)>0
                                        dunr=[dunr dunr(1)];
                                    end
                                    dunr=[dusa0(1) dunr dusa0(end)];
                                end
                                %       if length(dutha)>0
                                %        dunr=[dunr dunra];
                                %       end
                                
                                
                                radii.array{count,1}{13}=dunr;
                                
                                %       dunr
                                %    'ICI lens _nr lay_new', pausak
                                %      [dilu,ailu,nilu]=lens_sub(duth,Ra,ra,Ndisc,dunr,UD);
                                if ~exist('duth')
                                    duth=[];
                                end
                                
                                
                                if length(Np_ad)>0
                                    Npair(2)=sNpa*Np_ad;
                                end
                                if exist('gra')~=1
                                    gra.pos=0;
                                    gra.thg=0;
                                    gra.thga=0;
                                    gra.d=0;
                                    gra.LA=0;
                                end
                                
                                ' prima di lens_sub in lay_tran  OOOOOOO ', keyboard
                                [dilu,ailu,nilu,filu]=lens_sub(duth,ha,ra,Ndisc,dunr,UD,-10,Rel,Nrel,sNpa*Npair,Rflat,Rm_rel,gra);
                                %      ' ipar ',keyboard
                                %      ' dopo di lens_sub ', keyboard
                                %      ' dopo di lens_sub ', keyboard
                                %     ' dopo di lens_sub ', keyboard
                                d=[d; dilu];
                                sd=size(ipar)
                                dua=zeros(length(dilu),3,sd(3));
                                ipar=[ipar; dua];
                                Dov=[Dov zeros(size(dilu'))];
                                anyf=[anyf; zeros(size(dilu))];
                                %      fls=[fls; [ones(size(dilu)) zeros(size(dilu))] ];
                                %' fls', keyboard
                                %' fls', keyboard
                                %      fls=[fls; [abs(Npair)*ones(size(dilu)) zeros(size(dilu))] ];
                                fls=[fls; fliplr(filu) ];
%								'fls2', keyboard
                                iauto=[iauto; [zeros(size(dilu)) zeros(size(dilu))] ];
                                ifield=[ifield zeros(size(dilu'))];
                                %      shav=[shav; ones(size(dilu))];
                                sra=size(ailu);
                                srad=size(radii.a);
                                srn=size(n);
                                if sra(2)>srad(2)
                                    radii.a=[[radii.a zeros(srad(1),sra(2)-srad(2))]; ailu];
                                    radii.b=[[radii.b zeros(srad(1),sra(2)-srad(2))]; ailu*0];
                                    radii.c=[[radii.c zeros(srad(1),sra(2)-srad(2))]; ailu*0];
                                    shav=[[shav zeros(srad(1),sra(2)-srad(2))]; ones(size(ailu))*sha];
                                else
                                    radii.a=[radii.a ; [ailu zeros(sra(1),srad(2)-sra(2))]];
                                    radii.b=[radii.b; zeros(sra(1),srad(2))];
                                    radii.c=[radii.c; zeros(sra(1),srad(2))];
                                    shav=[shav;  ones(sra(1),srad(2))*sha];
                                end
                                sen=size(n);
                                senu=size(nilu);
                                if senu(2)>sen(2)
                                    n=[[n zeros(srn(1),sra(2)-srn(2)+1)]; nilu];
                                    xm=[[xm -10*ones(srn(1),sra(2)-srn(2)+1)]; -10*ones(size(nilu))];
                                else
                                    n=[n; [nilu zeros(senu(1),sen(2)-senu(2))]];
                                    xm=[xm; -10*ones(senu(1),sen(2))];
                                end
                                count_lay=count_lay+1;
                                
                                'lay1', keyboard
                                
                                count=count+length(dilu);
                                countw=countw+length(dilu);
                                %      'ICI lens fine dopo in lay_tran', pausak
                                %      ' ipar dopo ',keyboard
                        end
                        
                        % [nref,aref,n,d,xm,Dov,radii,flsv,perd,anyf,iauto,...
                        %         dw,xw,fsw,dov,shav,ipar,ifield,lambdaNot]= Lay_new(fileName)
                        % *Dov *flsv *anyf *iauto ipar *ifield
                        %  dw xw fsw dov : questi no!
                        
                        
                        ictr=ictrs;
                        %**************************************
                    case 'DOE'
                        
                        
                        
                        ictrs=ictr;
                        ictr=ictrar;
                        sha=8;
                        
                        [typemat1,remainder]=strtok(remainder);
                        switch typemat1
                            
                            case '_period'
                                
                                [ictP,ipar,Ddu]=fipar(remainder,ictP,count,ipar,sha,'D=',ictr);
                                [ictP,ipar,ddu]=fipar(remainder,ictP,count,ipar,sha,'d=',ictr);
                                [ictP,ipar,Ry]=fipar(remainder,ictP,count,ipar,sha,'Ry=',ictr);
                                
                                pe=Ddu+ddu;
                                NS=ceil(Ry/pe);
                                pe0=[Ddu Ddu+ddu];
                                duth=pe0;
                                for k=1:NS-1
                                    duth=[duth pe0+k*pe];
                                end
                                
                                %       dup=[[0; duth'] [duth'; duth(end)+10]];
                                dup=[[0; duth(1:end-1)'] [duth(1:end-1)'; duth(end)+10]];
                                val=zeros(size(dup));
                                val(1:2:end,:)=1;
                                ru=reshape(dup',1,prod(size(dup)));
                                pu=reshape(val',1,prod(size(dup)));
                                
                                figure, plot(ru,pu,'.-r'), axis([0 ru(end) 0 1.1]); pausak
                                
                                radii.array{count-1,1}{5}=Ddu;
                                radii.array{count-1,1}{6}=ddu;
                                radii.array{count-1,1}{7}=Ry;
                                radii.array{count-1,1}{10}=duth;
                                radii.a(count-1,1)=duth(1);
								'qui;', keyboard
                                
                            case '_ri'
                                
                                pubar=find(remainder=='|');
                                if length(pubar)==1
                                    ictP=ictP+1;
                                    ipar(count,1,ictP)=str2num(remainder(pubar+2:end));
                                    ipar(count,2,ictP)=-4;
                                    ipar(count,3,ictP)=ictr;
                                end
                                
                                
                                clear duth
                                kpar=0;
                                while length(remainder)>0
                                    kpar=kpar+1;
                                    [pach,remainder]=strtok(remainder);
                                    duth(kpar)=str2num(pach);
                                end
                                %       duth
                                radii.array{count,1}{12}=duth;
                                
                                dup=[[0; duth'] [duth'; duth(end)+10]];
                                val=zeros(size(dup));
                                val(1:2:end)=1;
                                ru=reshape(dup',1,prod(size(dup)));
                                pu=reshape(val',1,prod(size(dup)));
                                
                                radii.array{count-1,1}{10}=duth;
                                radii.a(count-1,1)=duth(1);
								'qui;', keyboard								
                                %      radii.b(count-1,1:length(duth))=duth*0;
                                %      radii.c(count-1,1:length(duth))=duth*0;
                                
                                %      figure, plot(ru,pu,'.-r'), axis([0 ru(end) 0 1.1]); pausak
                                
                            case '_lenssw'
                                sha=-8;
                                %    'DOEET', keyboard
                                %DOE _lenssw  RoC=0 R=20|P8 Ndis=0 DC=0.5 Pe=0.5
                                
                                [ictP,ipar,RoC]=fipar(remainder,ictP,count,ipar,sha,'RoC=',ictrar);
                                
                                [ictP,ipar,R]=fipar(remainder,ictP,count,ipar,sha,'R=',ictrar);
                                
                                [ictP,ipar,Ns]=fipar(remainder,ictP,count,ipar,sha,'Ndis=',ictrar);
                                
                                [ictP,ipar,DC]=fipar(remainder,ictP,count,ipar,sha,'DC=',ictrar);
                                
                                [ictP,ipar,Pe]=fipar(remainder,ictP,count,ipar,sha,'Pe=',ictrar);
                                
                                
                                pgrat.igra=1;
                                pgrat.DC=DC;
                                pgrat.n1=n(end,1);
                                pgrat.n2=n(end,2);
                                pgrat.period=Pe;
                                %       radii.array
                                
                                %       [h,R]=doeRoC(RoC,ra,Ns);
                                if Ns>0
                                    [h,R]=doeRoCg(RoC,ra,Ns,pgrat);
                                    d(end)=h*1000;
                                end
                                duth=R;
                                
                                
                                %       'DOEET', keyboard
                                %       duth
                                radii.array{count-1,1}{12}=duth;
                                radii.array{count,1}{12}=duth;
                                radii.array{count-1,1}{13}=pgrat;
                                
                                dup=[[0; duth'] [duth'; duth(end)+10]];
                                val=zeros(size(dup));
                                val(1:2:end)=1;
                                ru=reshape(dup',1,prod(size(dup)));
                                pu=reshape(val',1,prod(size(dup)));
                                
                                radii.a(count-1,1)=duth(1);
								'qui;', keyboard								
                                %      radii.b(count-1,1:length(duth))=duth*0;
                                %      radii.c(count-1,1:length(duth))=duth*0;
                                
                                %      figure, plot(ru,pu,'.-r'), axis([0 ru(end) 0 1.1]); pausak
                                fib=find(ipar(end,2,:)==2 | ipar(end,2,:)==-6 );
                                %              if length(fib)==1
                                %              fib
                                %              'fib', keyboard
%%%%                                IPar(count_lay,fib)=squeeze(ipar(count,1,fib))';
                                
                                %              IPar(count_lay,:)=ipar(end,1:2,fib);
                                %              IPar(count_lay,2)=2;
                                %              end
                                
                                %'DEO', keyboard
                                
                            case '_nr'
                                clear dunr
                                for kpar=1:NlAR+2
                                    [pach,remainder]=strtok(remainder);
                                    dunr(kpar)=str2num(pach);
                                end
                                %       dunr
                                %      'ICI lens _nr lay_new', pausak
                                %      'ICI lens fine', pausak
                                %      [dilu,ailu,nilu]=lens_sub(duth,Ra,ra,Ndisc,dunr,UD);
                                [dilu,ailu,nilu]=lens_sub(duth,ha,ra,Ndisc,dunr,UD,-10,Rel,Nrel);
                                ' dopo lenst', keyboard
                                %      ' ipar ',keyboard
                                d=[d; dilu];
                                dua=zeros(length(dilu),3,5);
                                ipar=[ipar; dua];
                                Dov=[Dov zeros(size(dilu'))];
                                anyf=[anyf; zeros(size(dilu))];
                                fls=[fls; [zeros(size(dilu)) ones(size(dilu))] ];
%									'fls3', keyboard
                                iauto=[iauto; [zeros(size(dilu)) zeros(size(dilu))] ];
                                ifield=[ifield zeros(size(dilu'))];
                                %      shav=[shav; ones(size(dilu))];
                                sra=size(ailu);
                                srad=size(radii.a);
                                srn=size(n);
                                radii.a=[[radii.a zeros(srad(1),sra(2)-srad(2))]; ailu];
                                radii.b=[[radii.b zeros(srad(1),sra(2)-srad(2))]; ailu*0];
                                radii.c=[[radii.c zeros(srad(1),sra(2)-srad(2))]; ailu*0];
                                shav=[[shav zeros(srad(1),sra(2)-srad(2))]; ones(size(ailu))];
                                sen=size(n);
                                senu=size(nilu);
                                if senu(2)>sen(2)
                                    n=[[n zeros(srn(1),sra(2)-srn(2)+1)]; nilu];
                                    xm=[[xm -10*ones(srn(1),sra(2)-srn(2)+1)]; -10*ones(size(nilu))];
                                else
                                    n=[n; [nilu zeros(senu(1),sen(2)-senu(2))]];
                                    xm=[xm; -10*ones(senu(1),sen(2))];
                                end
                                count_lay=count_lay+1;
                                'lay2', keyboard
                                
                                count=count+length(dilu);
                                countw=countw+length(dilu);
                        end
                        
                        ictr=ictrs;
                        
                        %**************************************
                        
                        
                    case 'Grating'
                        
                        ictrs=ictr;
                        ictr=ictrar;
                        sha=6;
                        [ictP,ipar,D]=fipar(remainder,ictP,count,ipar,sha,'D=',ictr);
                        
                        [ictP,ipar,drea]=fipar(remainder,ictP,count,ipar,sha,'d=',ictr);
                        
                        [ictP,ipar,Shap]=fipar(remainder,ictP,count,ipar,sha,'shape=',ictr);
                        
                        [ictP,ipar,Ry]=fipar(remainder,ictP,count,ipar,sha,'Ry=',ictr);
                        
                        [ictP,ipar,Rx]=fipar(remainder,ictP,count,ipar,sha,'Rx=',ictr);
                        
                        [ictP,ipar,ori]=fipar(remainder,ictP,count,ipar,sha,'orien.=',ictr);
                        
                        [ictP,ipar,shif]=fipar(remainder,ictP,count,ipar,sha,'shift=',ictr);
                        
                        [ictP,ipar,circle]=fipar(remainder,ictP,count,ipar,sha,'circle=',ictr);
                        
                        [ictP,ipar,Rext]=fipar(remainder,ictP,count,ipar,sha,'Rext=',ictr);
                        if ori~=0
                            keyboard
                        end
                        
                        
                        par='i_grap=';
                        fstr=findstr(remainder,par);
                        if length(fstr)>0
                            du=strtok(remainder(fstr+length(par):end));
                            i_grap=str2num(du);
                        else
                            i_grap=0;
                        end
                        
                        ictr=ictrs;
                        for pcounces=counces
                            radii.arrayd{1}=labl;
                            radii.arrayd{2}=ori;
                            radii.arrayd{3}=shif;
                            radii.arrayd{4}=circle;
                            radii.arrayd{5}=D;
                            radii.arrayd{6}=drea;
                            radii.arrayd{7}=Ry;
                            radii.arrayd{8}=Rx;
                            radii.arrayd{9}=Shap;
                            radii.arrayd{10}=i_grap;
                            radii.arrayd{11}=Rext;
                        end
                        %   'ICI grat lay_new'
                        %   keyboard
                        
                    case 'hcg'
                        
                        %'ICI HCG inizio'
                        %keyboard
                        ictrs=ictr;
                        ictr=ictrar;
                        sha=6;
                        
                        
                        
                        [typemat1,remainder]=strtok(remainder);
                        %      'typemat1', typemat1, keyboard
                        
                        
                        %'inizio typemat1', keyboard
                        
                        switch typemat1
                            
                            case '_pe'
                                
                                %'inizio fipar', keyboard
                                
                                [ictP,ipar,Pe]=fipar(remainder,ictP,count,ipar,sha,'D=',ictr);
                                
                                [ictP,ipar,d1]=fipar(remainder,ictP,count,ipar,sha,'d=',ictr);
                                
                                [~,~,NlAR]=fipar(remainder,ictP,count,ipar,sha,'nlay=',ictr);
                                
                            case '_par'
                                
                                %'inizio fipar', keyboard
                                
                                [ictP,ipar,Pe]=fipar(remainder,ictP,count,ipar,sha,'PE=',ictr);
                                
                                [~,~,NlAR]=fipar(remainder,ictP,count,ipar,sha,'nlay=',ictr);
                                
                            case '_th'
                                
                                
                                clear duth
                                for kpar=1:NlAR
                                    
                                    [pach,remainder]=strtok(remainder);
                                    pubar=find(pach=='|');
                                    
                                    rema=pach;
                                    if length(pubar)>0
                                        rema=pach(1:pubar-1);
                                        param=pach(pubar+2:end);
                                    end
                                    %         'pubar', keyboard
                                    if length(pubar)==1
                                        ictP=ictP+1;
                                        ipar(count-1,1,ictP)=str2num(param);
                                        ipar(count-1,2,ictP)=-3;
                                        ipar(count-1,3,ictP)=kpar;
                                        %          ' ictr ', keyboard
                                    end
                                    
                                    duth(kpar)=str2num(rema);
                                end
                                
                                
                            case '_ab'
                                
                                clear duab
                                for kpar=1:NlAR
                                    
                                    [pach,remainder]=strtok(remainder);
                                    pubar=find(pach=='|');
                                    
                                    rema=pach;
                                    if length(pubar)>0
                                        rema=pach(1:pubar-1);
                                        param=pach(pubar+2:end);
                                    end
                                    %         'pubar', keyboard
                                    if length(pubar)==1
                                        ictP=ictP+1;
                                        ipar(count-1,1,ictP)=str2num(param);
                                        ipar(count-1,2,ictP)=-40;
                                        ipar(count-1,3,ictP)=kpar;
                                        %          ' ictr ', keyboard
                                    end
                                    
                                    duab(kpar)=str2num(rema);
                                end
                                
                                
                            case '_nr'
                                
                                
                                clear dunr
                                for kpar=1:NlAR
                                    
                                    [pach,remainder]=strtok(remainder);
                                    pubar=find(pach=='|');
                                    
                                    rema=pach;
                                    if length(pubar)>0
                                        rema=pach(1:pubar-1);
                                        param=pach(pubar+2:end);
                                    end
                                    %         'pubar', keyboard
                                    if length(pubar)==1
                                        ictP=ictP+1;
                                        ipar(count-1,1,ictP)=str2num(param);
                                        ipar(count-1,2,ictP)=-4;
                                        ipar(count-1,3,ictP)=kpar;
                                        %          ' ictr ', keyboard
                                    end
                                    
                                    dunrdu=str2num(rema);
                                    
                                    
                                    if exist('duab')
                                        dunrdu=dunrdu-j*100*duab(kpar)/(2*kNoll);
                                    end
                                    dunr(kpar)=dunrdu;
                                    
                                end
                                
                                
                                
                                
                                %' qui npari', keyboard
                                
                                par_grat.nv=dunr;
                                par_grat.tv=duth;
                                par_grat.d1=d1;
                                par_grat.d2=Pe-d1;
                                
                                ictr=ictrs;
                                Ry=0;
                                for pcounces=counces
                                    radii.arrayd{5}=Pe;
                                    radii.arrayd{6}=d1;
                                    radii.arrayd{7}=Ry;
                                    radii.arrayd{12}=duth;
                                    radii.arrayd{13}=dunr;
                                end
                        end
                        
                        %'ICI grat HCG fine'
                        %keyboard
                        
                        
                        
                    case 'Particle'
                        ipardu=ipar;
                        ictPdu=ictP;
                        coup=count+1;
                        
                        %          [ictPdu,ipardu,L_pos]=fipar(remainder,ictPdu,coup,ipardu,sha,'pos=',ictr);
                        %          [ictPdu,ipardu,L_thg]=fipar(remainder,ictPdu,coup,ipardu,sha,'thg=',ictr);
                        %          [ictPdu,ipardu,
                        
                        %      ictrs=ictr;
                        %      ictr=ictrar;
                        sha=9;
                        
                        %Thick=5000|P6 Ref=1.5 Dp=10|P21 Cp=2000|P22 Ndis=20|P23 Refp=1.7
                        
                        [ictPdu,ipardu,Th]=fipar(remainder,ictPdu,coup,ipardu,sha,'Thick=',ictr);
                        [ictPdu,ipardu,Refl]=fipar(remainder,ictPdu,coup,ipardu,sha,'Ref=',ictr);
                        %      [ictP,ipar,Shap]=fipar(remainder,ictP,count,ipar,sha,'shape=',ictr);
                        [ictPdu,ipardu,Dp]=fipar(remainder,ictPdu,coup,ipardu,sha,'Dp=',ictr);
                        [ictPdu,ipardu,Cp]=fipar(remainder,ictPdu,coup,ipardu,sha,'Cp=',ictr);
                        [ictPdu,ipardu,Ndispa]=fipar(remainder,ictPdu,coup,ipardu,sha,'Ndis=',ictr);
                        [ictPdu,ipardu,Refp]=fipar(remainder,ictPdu,coup,ipardu,sha,'Refp=',ictr);
                        [ictPdu,ipardu,sha_p]=fipar(remainder,ictPdu,coup,ipardu,sha,'sha=',ictr);
                        if length(sha_p)==0
                            sha_p=7;
                        end
                        ipar=ipardu;
                        ictP=ictPdu;
                        radii.arrayd{1}=Dp;
                        radii.arrayd{2}=Cp;
                        radii.arrayd{3}=Ndispa;
                        radii.arrayd{4}=Refp;
                        radii.arrayd{5}=Th;
                        radii.arrayd{6}=Refl;
                        radii.arrayd{7}=sha_p;
                        Rp=Dp/2;
                        
                        [dilu,ailu,nilu,filu]=particle(Rp,Ndispa,Cp,Th,Refp,Refl,-10);
                        
                        ' dopo particle ', keyboard
                        %             sha=1;
                        
                        d=[d; dilu];
                        sd=size(ipar)
                        dua=zeros(length(dilu),3,sd(3));
                        ipar=[ipar; dua];
                        Dov=[Dov zeros(size(dilu'))];
                        anyf=[anyf; zeros(size(dilu))];
                        fls=[fls; fliplr(filu) ];
%							'fls4', keyboard
                        iauto=[iauto; [zeros(size(dilu)) zeros(size(dilu))] ];
                        ifield=[ifield zeros(size(dilu'))];
                        sra=size(ailu);
                        srad=size(radii.a);
                        srn=size(n);
                        if sra(2)>srad(2)
                            radii.a=[[radii.a zeros(srad(1),sra(2)-srad(2))]; ailu];
                            radii.b=[[radii.b zeros(srad(1),sra(2)-srad(2))]; ailu*0];
                            radii.c=[[radii.c zeros(srad(1),sra(2)-srad(2))]; ailu*0];
                            shadd=ones(size(ailu))*sha;
                            %        shadd(1,:)=0;
                            %        shadd(end,:)=0;
                            shav=[[shav zeros(srad(1),sra(2)-srad(2))]; shadd];
                        else
                            radii.a=[radii.a ; [ailu zeros(sra(1),srad(2)-sra(2))]];
                            radii.b=[radii.b; zeros(sra(1),srad(2))];
                            radii.c=[radii.c; zeros(sra(1),srad(2))];
                            shadd=ones(sra(1),srad(2))*sha;
                            %        shadd(1,:)=0;
                            %        shadd(end,:)=0;
                            shav=[shav;  shadd];
                        end
                        sen=size(n);
                        senu=size(nilu);
                        if senu(2)>sen(2)
                            n=[[n zeros(srn(1),sra(2)-srn(2)+1)]; nilu];
                            xm=[[xm -10*ones(srn(1),sra(2)-srn(2)+1)]; -10*ones(size(nilu))];
                        else
                            n=[n; [nilu zeros(senu(1),sen(2)-senu(2))]];
                            xm=[xm; -10*ones(senu(1),sen(2))];
                        end
                        count_lay=count_lay+1;
                        'lay3', keyboard
                        
                        counces=count;
                        count=count+length(dilu);
                        countw=countw+length(dilu);
                        'fine partocle', keyboard
                        sha=1;
                        
                    case 'PARAB'
                        
                        %      lablar=3;
                        
                        ictrs=ictr;
                        ictr=ictrar;
                        sha=-shav(end);
                        shav(end)=sha;
                        igrad=1;
                        %      ' qui |GRAD', keyboard
                        
                        [ictP,ipar,Ishag]=fipar(remainder,ictP,count,ipar,sha,'Isha=',ictr);
                        
                        [ictP,ipar,Drg]=fipar(remainder,ictP,count,ipar,sha,'Dr=',ictr);
                        
                        [ictP,ipar,Nrg]=fipar(remainder,ictP,count,ipar,sha,'Nrdis=',ictr);
                        
                        [ictP,ipar,Iappg]=fipar(remainder,ictP,count,ipar,sha,'Iappg=',ictr);
                        
                        ictr=ictrs;
                        
                        %      'ICI grad 0 lay_new'
                        %      keyboard
                        
                        parray{1}=Ishag;
                        parray{2}=Drg;
                        parray{3}=Nrg;
                        parray{4}=Iappg;
                        radii.array{count-1}=parray;
                        
                        %      'ICI grat lay_new'
                        %      keyboard
                        
                    case 'Shapes'
                        
                        while ~isempty(remainder)
                            [tycen,remainder]=strtok(remainder);
                            cen=shaf(tycen);
                            cemem(ipcen,1)=cen;
                            ipcen=ipcen+1;
                        end
                        
                        for pcounces=counces
                            %       radii.array{pcounces,4,ictrar}=cemem;
                            radii.arrayd{4}=cemem;
                        end
                        
                    case 'Ry'
                        
                        while ~isempty(remainder)
                            [tycen,remainder]=strtok(remainder);
                            cen=str2num(tycen);
                            cemem(ipcen,1)=cen;
                            ipcen=ipcen+1;
                        end
                        
                        for pcounces=counces
                            %       radii.array{pcounces,5,ictrar}=cemem;
                            radii.arrayd{5}=cemem;
                        end
                        
                    case 'Rx'
                        while ~isempty(remainder)
                            [tycen,remainder]=strtok(remainder);
                            cen=str2num(tycen);
                            cemem(ipcen,1)=cen;
                            ipcen=ipcen+1;
                        end
                        
                        for pcounces=counces
                            %       radii.array{pcounces,6,ictrar}=cemem;
                            radii.arrayd{6}=cemem;
                        end
                        
                    case 'Delta'
                        while ~isempty(remainder)
                            [tycen,remainder]=strtok(remainder);
                            cen=str2num(tycen);
                            cemem(ipcen,1)=cen;
                            ipcen=ipcen+1;
                        end
                        
                        for pcounces=counces
                            %       radii.array{pcounces,7,ictrar}=cemem;
                            radii.arrayd{7}=cemem;
                        end
                        
                    otherwise
                        
                        if setarray>0
                            disp(' forgotten array centers ')
                            keyboard
                        end
                        
                end
                
            end
            if isfield(radii,'arrayd')==1
											'qui;', keyboard
                for pcounces=counces
                    %    radii.a(pcounces,ictrar)=radii.array{pcounces,5,ictrar}(1);
                    %    radii.b(pcounces,ictrar)=radii.array{pcounces,6,ictrar}(1);
                    %    radii.c(pcounces,ictrar)=radii.array{pcounces,7,ictrar}(1);
                    %'wui',keyboard
                    radii.array{pcounces,ictrar}=radii.arrayd;
                    radii.a(pcounces,ictrar)=radii.arrayd{5}(1);
                    radii.b(pcounces,ictrar)=radii.arrayd{6}(1);
                    radii.c(pcounces,ictrar)=radii.arrayd{7}(1);
                end
            end
            %  disp(' fine array')
            %  pausak
            
        otherwise
            
            fan=findstr(typelay,'%');
            if (length(typelay)<8 | length(findstr(typelay,'Centers'))>0) & length(fan)==0
                iany=0;
                iautoc=0;
                itip=0;
                
                if length(typelay)>1
                    fan=findstr(typelay,'a');
                    if length(fan)==1
                        iany=1;
                    end
                    fan=findstr(typelay,'r');
                    if length(fan)==1
                        iautoc=1;
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
%                            'C',keyboard
                        if itip<0
                            itip=-itip;
                        else
                            itip=-4;
                        end
                    end
                    fan=findstr(typelay,'DD');
                    if length(fan)==1
                        iDD=[iDD count_lay+1];
                        itip=-100;
                        %     'qui DD' ,  keyboard
                    end
                    fan=findstr(typelay,'NEGF');
					if itip==-500
					itip=0;
					end
                    if length(fan)==1 | identroNEFG==1
					 if length(fan)==1
					  identroNEFG=1;
                        itip=-500;
					 end

%                         iETCHdd=[iETCHdd count_lay+2-iDD];
                         iETCHdd=[iETCHdd count_lay+2];
%                            'qui DD' ,  keyboard
                    end                    
					
                    if length(fan)==1 
                        iNEGF=[iNEGF nlongcnt];
%                             'qui DD' ,  keyboard
                    end  					
                    
                    fan=findstr(typelay,'Po');
                    if length(fan)==1
                        itip=-3;
                    end
                    
                    fan=findstr(typelay,'Lc');
                    ifif=0;
                    if length(fan)==1
                        typelay='Lc';
                        ifif=ifif+1;
                    end
                    fan=findstr(typelay,'Lg');
                    if length(fan)==1
                        typelay='Lg';
                        ifif=ifif+1;
                    end
                    if ifif==0
                        typelay='';
                    end
                    
                end  % if length(typelay)>2
                
                if iautoc>0
                    iautoci=iautoci+iautoc;
                end
                
                
            end  % if lenght
            
            %  typelay,
            %  pausak
            switch typelay
                
                case 'Wavelength'
                    [type,remainder]=strtok(remainder);
                    lam=str2num(type);
                    
                    
                    lambdaNot=lam*1e-9;
                    kNoll=2*pi/lambdaNot;
                    
                case 'Igrad_transport'
                    [type,remainder]=strtok(remainder);
                    igrad=str2num(type);
                    
                case 'Reference'
                    [type,remainder]=strtok(remainder);
                    
                    switch type
                        
                        case 'AlGaAs'
                            
                            [secondWord,remainder]=strtok(remainder);
                            AlCont=str2num(secondWord);
                            if AlCont==0
                                AlCont=1e-8;
                            end
                            
                            nref=real(nAlGaAs(lambdaNot,AlCont));
                            
                        case 'InGaAs'
                            
                            [secondWord,remainder]=strtok(remainder);
                            AlCont=str2num(secondWord);
                            if AlCont==0
                                AlCont=1e-8;
                            end
                            
                            nref=real(nInGaAs(lambdaNot,AlCont));
                            
                        case 'Ref'
                            
                            [secondWord,remainder]=strtok(remainder);
                            nref=str2num(secondWord);
                            Alcont=-1;
                            
                    end
                    
                    [Word,remainder]=strtok(remainder);
                    aref=str2num(Word);
                    
                case 'Lc'
                    
                    ictP=0;
                    ictr=1;
                    
                    [word2,remainder]=strtok(remainder);
                    %pausak
                    
                    
                    if isequal(word2(1),'F')
                        %     word2
                        %     pausak
                        
                        if length(word2)==2
                            if isequal(word2(2),'s')
                                ifield(count)=-1;
                            else
                                ifield(count)=-2;
                            end
                        elseif length(word2)==1
                            ifield(count)=1;
                        elseif length(word2)>2
                            fist=findstr(word2,'@');
                            zff=(str2num(word2(fist+1:end-2)));
                            if length(findstr(word2,'ang'))==1
                                ifield(count)=-zff*1e4;
                            else
                                ifield(count)=zff*1e4;
                            end
                            
                        end
                        [word2,remainder]=strtok(remainder);
                    end
                    
                    
                    if length(word2)>5
                        dop=fix(str2num(word2(4:end)));
                        %'NA', keyboard
                        if isequal(word2(1:2),'Na')
                            dov(countw)=dop;
                            Dov(count)=dop;
                        elseif isequal(word2(1:2),'Nd')
                            dov(countw)=-dop;
                            Dov(count)=-dop;
                        end
                        %   word2, pausak
                        [word2,remainder]=strtok(remainder);
                    end
                    
                    pubar=find(word2=='|');
                    if length(pubar)==1
                        ictP=ictP+1;
                        ipar(count,1,ictP)=str2num(word2(pubar+2:end));
                        ipar(count,2,ictP)=8;
                        ipar(count,3,ictP)=ictr;
                        word2=word2(1:pubar-1);
                    end
                    
                    nlong=fix(str2num(word2));
                    % nlongcnt=nlongcnt+nlong; % ATTENZIONE! viene usato un indice negativo per flaggare le QWs! Se voglio contare giusto, devo usare abs
                    nlongcnt=nlongcnt+abs(nlong);
					
%					'nlong2', keyboard
                    					
                    if abs(nlong)>1
                        nllo=nllo+1;
                        %    nllo
                        ipmem(count)=count;
                    end
                    
                    fls(count,1:2)=[nlong 0];
                    if itip==-100
%                       'fls', keyboard
                        flsKazzo=fls;
                        flsKazzo(count,2)=-1;
                    end
                    
%                       'fls 9', keyboard
                    if length(nlong)==0
                        nlong=0;
                    end
                    flw(countw)=[nlong];
                    
                    [word2,remainder]=strtok(remainder);
                    
                    pubar=find(word2=='|');
                    %   count
                    %   'lay', pausak
                    if length(pubar)==1
                        ictP=ictP+1;
                        ipar(count,1,ictP)=str2num(word2(pubar+2:end));
                        ipar(count,2,ictP)=1;
                        
                        ipar(count,3,ictP)=ictr;
                        word2=word2(1:pubar-1);
                    end
                    
                    %   count
                    thick=(str2num(word2));
                    d(count,1)=thick;
                    %   d
                    %   [count thick]
                    %   pausak
                    dw(countw,1)=thick;
                    anyf(count,1)=iany;
                    if iautoc>0
                        iauto(count,1)=iautoci;
                        iauto(count,2)=itip;
                    else
                        iauto(count,2)=itip;
                    end
                    
                    [type,remainder]=strtok(remainder);
                    %   pausak
                    
                    ictr=0;
                    ictrar=0;
                    
                    %******aggiunto


					count_lay=count_lay+1;

		            RaDD(count_lay,1)=0;					
                    radii.a(count,1)=0;
                    radii.b(count,1)=0;
                    radii.c(count,1)=0;
                    shav(count,1)=0;
%                    '******aggiunto 1673', keyboard                    
                    while isempty(type)==0
                        
                        if length(type)==1
                            ifty=type=='%';
                            if ifty==1
                                break
                            end
                        end
                        
                        ictr=ictr+1;
                        
                        %type
                        %pausak
                        
                        switch type
                            
                            %mater{count_lay}=type;
                            %'vero', keyboard
                            
                            case 'AlGaAs'
                                materi{ictr}=type;
                                
                                [secondWord,remainder]=strtok(remainder);
                                AlCont=str2num(secondWord);
                                if AlCont==0
                                    AlCont=1e-8;
                                end
                                %%
                                [secondWord,remainder]=strtok(remainder);
                                word2=secondWord;
                                pubar=find(word2=='|');
                                if length(pubar)==1
                                    ictP=ictP+1;
                                    ipar(count,1,ictP)=str2num(word2(pubar+2:end));
                                    ipar(count,2,ictP)=10;
                                    ipar(count,3,ictP)=ictr;
                                    %         ipar(count,3,ictP)=ictr-1;
                                    secondWord=word2(1:pubar-1);
                                end
                                
                                Absorp=str2num(secondWord);
                                
                                %'nAlGa', keyboard
                                n(count,ictr)=nAlGaAs(lambdaNot,AlCont)-i*Absorp*100/(2*kNoll);
                                %         n(count,ictr)
                                %        'n', pausak
                                xm(count,ictr)=AlCont;
                                xw(countw,ictr)=AlCont;
                                
                                
                            case 'InGaAs980'
                                materi{ictr}=type;
                                
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
                                
                            case 'InGaAs'
                                materi{ictr}=type;
                                [secondWord,remainder]=strtok(remainder);
                                AlCont=str2num(secondWord);
                                if AlCont==0
                                    AlCont=1e-8;
                                end
                                %%
                                [secondWord,remainder]=strtok(remainder);
                                word2=secondWord;
                                pubar=find(word2=='|');
                                if length(pubar)==1
                                    ictP=ictP+1;
                                    ipar(count,1,ictP)=str2num(word2(pubar+2:end));
                                    ipar(count,2,ictP)=10;
                                    ipar(count,3,ictP)=ictr;
                                    %         ipar(count,3,ictP)=ictr-1;
                                    secondWord=word2(1:pubar-1);
                                end
                                
                                Absorp=str2num(secondWord);
                                
                                
                                %        [fourthWord,remainder]=strtok(remainder);
                                %        Absorp=str2num(fourthWord);
                                
                                n(count,ictr)=real(nInGaAs(lambdaNot,AlCont))-i*Absorp*100/(2*kNoll);
                                xm(count,ictr)=AlCont;
                                xw(countw,ictr)=AlCont;
                                
                            case 'Ref'
                                materi{ictr}=type;
                                %      'materi', keyboard
                                [secondWord,remainder]=strtok(remainder);
                                word2=secondWord;
                                pubar=find(word2=='|');
                                if length(pubar)==1
                                    ictP=ictP+1;
                                    ipar(count,1,ictP)=str2num(word2(pubar+2:end));
                                    ipar(count,2,ictP)=5;
                                    ipar(count,3,ictP)=ictr;
                                    %         ipar(count,3,ictP)=ictr-1;
                                    secondWord=word2(1:pubar-1);
                                end
                                
                                ni=str2num(secondWord);
                                materi{ictr}=[type,'_',num2str(ni)];
                                %      'materi', keyboard
                                [secondWord,remainder]=strtok(remainder);
                                word2=secondWord;
                                pubar=find(word2=='|');
                                if length(pubar)==1
                                    ictP=ictP+1;
                                    ipar(count,1,ictP)=str2num(word2(pubar+2:end));
                                    ipar(count,2,ictP)=10;
                                    ipar(count,3,ictP)=ictr;
                                    %         ipar(count,3,ictP)=ictr-1;
                                    secondWord=word2(1:pubar-1);
                                end
                                
                                Absorp=str2num(secondWord);
                                n(count,ictr)=ni-i*Absorp*100/(2*kNoll);
                                xm(count,ictr)=-10;
                                xw(countw,ictr)=-10;
                                %'laytran'
                                %Absorp
                                %count
                                %pausak
                                
                            case 'Au'
                                materi{ictr}=type;
                                [nIRe,keIm]=nAu(lambdaNot);
                                n(count,ictr)=nIRe-i*keIm;
                                xm(count,ictr)=-1;
                                xw(countw,ictr)=-1;
                                
                            case 'Aur'
                                materi{ictr}=type;
                                
                                [nIRe,keIm]=nAur(lambdaNot);
                                n(count,ictr)=nIRe-i*keIm;
                                xm(count,ictr)=-1;
                                xw(countw,ictr)=-1;
                                
                            case 'Cr'
                                materi{ictr}=type;
                                
                                [nIRe,keIm]=nCr(lambdaNot);
                                n(count,ictr)=nIRe-i*keIm;
                                xm(count,ictr)=-4;
                                xw(countw,ictr)=-4;
                                
                            case 'Ti'
                                materi{ictr}=type;
                                [nIRe,keIm]=nTi(lambdaNot);
                                n(count,ictr)=nIRe-i*keIm;
                                xm(count,ictr)=-2;
                                xw(countw,ictr)=-2;
                                
                            case 'Pt'
                                materi{ictr}=type;
                                [nIRe,keIm]=nPt(lambdaNot);
                                n(count,ictr)=nIRe-i*keIm;
                                xm(count,ictr)=-3;
                                xw(countw,ictr)=-3;
                                
                            otherwise
                                %       keyboard
                                ictr=ictr-1;
                                Word=type;
                                
                                fan=findstr(Word,'|P');
                                if length(fan)>0
                                    ictP=ictP+1;
                                    pubar=find(Word=='|');
                                    ipar(count,1,ictP)=str2num(Word(pubar+2:end));
                                    ipar(count,2,ictP)=7;
                                    ipar(count,3,ictP)=ictr;
                                    Word=Word(1:fan-1);
                                end
                                
                                %   '2'
                                %        keyboard
                                fan=findstr(Word,':Lable=');
                                iseta=0;
                                if length(fan)>0
                                    pubar=find(Word==':');
                                    %        iseta=1;
                                    iseta=0;
                                    setarray=0;
                                    du=str2num(Word(pubar+7:end));
                                    ralabl(count,ictr)=-du;
                                    Word=Word(1:fan-1);
                                end
                                
                                switch Word
                                    
                                    case {'circle','rhombus','ellipse','rectangle','square','vortex'}
                                        
                                        [Word1,remainder]=strtok(remainder);
                                        
                                        word2=Word1;
                                        pubar=find(word2=='|');
                                        if length(pubar)==1
                                            ictP=ictP+1;
                                            ipar(count,1,ictP)=str2num(word2(pubar+2:end));
                                            ipar(count,2,ictP)=2;
                                            ipar(count,3,ictP)=ictr;
                                            Word1=word2(1:pubar-1);
                                        end
                                        radii.a(count,ictr)=str2num(Word1);

					                    if ictr==1
%										count_lay=count_lay+1
										end
		                                RaDD(count_lay,ictr)=str2num(Word1);
%										'qui radiale', keyboard
                                end
                                
                                switch Word
                                    
                                    case {'ellipse','rectangle','square','vortex'}
                                        
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
                                            case {'vortex'}
                                                shav(count,ictr)=11;
                                        end
                                    case {'rhombus'}
                                        
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
                                        
                                        
                                        %        if isequal(Word,'grating')
                                        %         [Word2,remainder]=strtok(remainder);
                                        %          word2=Word2;
                                        %          pubar=find(word2=='|');
                                        %          if length(pubar)==1
                                        %           ictP=ictP+1;
                                        %           ipar(count,1,ictP)=str2num(word2(pubar+2:end));
                                        %           ipar(count,2,ictP)=9;
                                        %           ipar(count,3,ictP)=ictr;
                                        %           Word2=word2(1:pubar-1);
                                        %         end
                                        %        end
                                        
                                        
                                    case {'array'}
                                        if iseta==0
                                            setarray=1;
                                        end
                                        counces=count;
                                        shav(count,ictr)=5;
                                    case {'grating'}
                                        shav(count,ictr)=6;
                                        counces=count;
                                    case {'HCG'}
                                        shav(count,ictr)=6;
                                        counces=count;
                                    case 'circle'
                                        shav(count,ictr)=1;
                                        radii.b(count,ictr)=0;
                                        radii.c(count,ictr)=0;
                                    case 'doe'
                                        shav(count,ictr)=-8;
                                        counces=count;
                                    otherwise
                                        fi=strncmpi(Word,'Nmesh=',6);
                                        if fi==1
                                            if ~exist('count_lay')
                                                cl=1;
                                            else
                                                cl=count_lay;
                                            end
                                            DDmeshDU=str2num(Word(7:end));
                                            
                                            if fls(end,1)==-1 & is_even(DDmeshDU)==0
                                                'NUMBER of z points in QW must be even: correct the str file'
                                                %keyboard
                                                DDmeshDU= DDmeshDU+1;
                                            end
                                            DDmesh(cl,1)=DDmeshDU;
%											'DDmesh DU;', keyboard
                                            
                                            clear Word
                                        else
                                            disp(' wrong type of transverse shape 1')
                                            keyboard
                                        end
                                end
                        end
                        [type,remainder]=strtok(remainder);
                        
                    end  % while trasv.
                    %'count',keyboard
                    if count==1
                        count_lay=0;
                    end
                    if count>=2
                        %    if count>=0
%                        count_lay=count_lay+1;
                        DDd(count_lay,1)=d(count);
                        if itip==-4
                            DDd(count_lay,1)=-d(count);
                        end
                        DDmat{count_lay,1}=materi;
                        % count
                        if count==8
                            %          '8;\;', keyboard
                        end
%                                         'IPar 2014', keyboard    
                        if size(ipar,1)>=count
                            fib=find(ipar(count,2,:)==2 | ipar(count,2,:)==-6 );
                            %              if length(fib)==1
                            %              fib
                            %              'fib', keyboard
                            IPar(count_lay,fib)=squeeze(ipar(count,1,fib))';
                            %              end
                            %             'IPar', keyboard
                            %        'mat', keyboard
						else
						 IPar(count_lay,:)=0;
						
                        end
                       Pun_cont(count_lay)=count;
                    
                        clear materi
                        if ~exist('AlCont')
                            Alx=-1;
                        else
                            Alx=AlCont;
                        end
                        
                        if Alx==1e-8;
                            Alx=0;
                        end
                        DDx(count_lay,1:2)=Alx*[1 1];
                        if exist('dop')
                            if length(Dov)==count
                                Dddov=Dov(count);
                            else
                                Dddov=0;
                            end
                            if Dddov>0
                                dddra=Dddov *[1 1];
                                dddrn=[0 0];
                            else
                                dddrn=Dddov *[1 1];
                                dddra=[0 0];
                            end
                            DDna(count_lay,:)=dddra;
                            DDnd(count_lay,:)=dddrn;
                        end
                        
                                 %'lay4', pausak
                        clear AlCont AlCont1 AlCont2
 
                    end
                    count=count+1;
                    %pausak
                    countw=countw+1;
                    
                case 'Lg'       % graded layer
                    
                    ictP=0;
                    
                    fiN=findstr(remainder,'Ndis=');
                    %'fiN', keyboard
                    if length(fiN)~=1
                        disp( ' error in "file.str": Ndis is missing')
                        Nline
                        NL
                        keyboard
                    end
                    [wordN,remainder1]=strtok(remainder(fiN+5:end));
                    fiM=findstr(remainder1,'=');
                    if length(fiM)~=1
                        disp( ' error in "file.str": Nmesh is missing')
                        Nline
                        NL
                        keyboard
                    end
                    Nmesh=str2num(remainder1(fiM+1:end));
                    Ndis=str2num(wordN);
                    ndis=Ndis;
                    pucount=count:count+Ndis-1;
                    remainder=(remainder(1:fiN-1));
                    %    disp(' set Ndis ')
                    %   keyboard
                    
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
                    end
                    
                    if length(word2)>5
                        
                        dop=fix(str2num(word2(4:end)));
                        %'DOP 2', keyboard
                        if isequal(word2(1:2),'Na')
                            dovdu=dop;
                        elseif isequal(word2(1:2),'Nd')
                            dovdu=-dop;
                        end
                        [do2c,remainder]=strtok(remainder);
                        do2=str2num(do2c(4:end));
                        %   word2, pausak
                        [word2,remainder]=strtok(remainder);
                    end
                    if exist('do2')
                        if length(do2)==0
                            clear do2
                        end
                    end
                    
                    pubar=find(word2=='|');
                    if length(pubar)==1
                        ictP=ictP+1;
                        ipar(pucount,1,ictP)=str2num(word2(pubar+2:end));
                        ipar(pucount,2,ictP)=8;
                        ipar(pucount,3,ictP)=ictr;
                        word2=word2(1:pubar-1);
                    end
                    
                    nlong=fix(str2num(word2));
                    nlongcnt=nlongcnt+nlong;
                    %   fls(count,1)=nlong;
%					'nlong1', keyboard
                    
                    [word2,remainder]=strtok(remainder);
                    
                    pubar=find(word2=='|');
                    if length(pubar)==1
                        ictP=ictP+1;
                        ipar(pucount,1,ictP)=str2num(word2(pubar+2:end));
                        ipar(pucount,2,ictP)=1;
                        ipar(pucount,3,ictP)=ictr;
                        word2=word2(1:pubar-1);
                    end
                    
                    
                    thick=(str2num(word2));
                    radir=0;
                    %   nir=0;
                    clear nit xct
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
                        %   type
                        % 'lay'
                        % keyboard
                        
                        
                        switch type
                            
                            case 'AlGaAs'
                                materi{ictr}=type;
                                
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
                                
                                %        [fourthWord,remainder]=strtok(remainder);
                                %        ndis=fix(str2num(fourthWord));
                                %         ndis=Ndis;
                                
                                %        n(count,ictr)=nAlGaAs(lambdaNot,AlCont)-i*Absorp*100/(2*kNoll);
                                %        xm(count,ictr)=AlCont;
                                
                                x=[AlCont1 AlCont2];
                                
                                
                                P=-100/(2*kNoll)*[Absorp1 Absorp2];
                                if P(1)==P(2)
                                    cxp=[0 P(1)];
                                else
                                    cxp=polyfit(x,P,1);
                                end
                                
                                %   ' gradpol', pausak
                                %   keyboard
                                [Ld ,nd, xc]=gradpol(thick,AlCont1,AlCont2,lambdaNot,ndis,cxp);
                                nit(1:ndis,ictr)=nd;
                                xct(1:ndis,ictr)=xc;
                                
                                
                            case 'InGaAs'
                                materi{ictr}=type;
                                
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
                                
                                %        [fourthWord,remainder]=strtok(remainder);
                                %        ndis=fix(str2num(fourthWord));
                                %         ndis=Ndis;
                                
                                %        n(count,ictr)=nAlGaAs(lambdaNot,AlCont)-i*Absorp*100/(2*kNoll);
                                %        xm(count,ictr)=AlCont;
                                
                                x=[AlCont1 AlCont2];
                                P=-100/(2*kNoll)*[Absorp1 Absorp2];
                                if P(1)==P(2)
                                    cxp=[0 P(1)];
                                else
                                    cxp=polyfit(x,P,1);
                                end
                                
                                %   ' gradpol', pausak
                                [Ld ,nd, xc]=gradpol(thick,AlCont1,AlCont2,lambdaNot,ndis,cxp,1);
                                nit(1:ndis,ictr)=nd;
                                xct(1:ndis,ictr)=xc;
                                
                                
                                
                                
                            case 'Ref'
                                materi{ictr}=type;
                                %'ref', keyboard
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
                                
                                
                                ni1=str2num(secondWord);
                                materi{ictr}=[type,'_',num2str(ni1)];
                                
                                [secondWord,remainder]=strtok(remainder);
                                ni2=str2num(secondWord);
                                
                                [secondWord,remainder]=strtok(remainder);
                                Absorp1=str2num(secondWord);
                                
                                [secondWord,remainder]=strtok(remainder);
                                Absorp2=str2num(secondWord);
                                nirdu1=ni1-i*Absorp1*100/(2*kNoll);
                                nirdu2=ni2-i*Absorp2*100/(2*kNoll);
                                
                                [Ld ,nd, xc]=gradref(thick,nirdu1,nirdu2,ndis);
                                nit(1:ndis,ictr)=nd;
                                xct(1:ndis,ictr)=-10;
                                %        n(count,ictr)=ni;
                                %        xm(count,ictr)=-1;
                                %        ictr=2;
                                AlCont1=-1;
                                AlCont2=-1;
                                %         ' fine qui', keyboard
                                
                            case 'Au'
                                materi{ictr}=type;
                                [nIRe,keIm]=nAu(lambdaNot);
                                ni=nIRe-i*keIm;
                                nit(1:ndis,ictr)=ni;
                                xct(1:ndis,ictr)=-1;
                                
                                
                            case 'Aur'
                                materi{ictr}=type;
                                [nIRe,keIm]=nAur(lambdaNot);
                                ni=nIRe-i*keIm;
                                nit(1:ndis,ictr)=ni;
                                xct(1:ndis,ictr)=-1;
                                
                            case 'Cr'
                                materi{ictr}=type;
                                [nIRe,keIm]=nCr(lambdaNot);
                                ni=nIRe-i*keIm;
                                %        nir(ictr)=ni;
                                nit(1:ndis,ictr)=ni;
                                xct(1:ndis,ictr)=-4;
                                %        n(count,ictr)=nIRe-i*keIm;
                                %        xm(count,ictr)=-1;
                                
                            case 'Ti'
                                materi{ictr}=type;
                                [nIRe,keIm]=nTi(lambdaNot);
                                ni=nIRe-i*keIm;
                                %        nir(ictr)=ni;
                                nit(1:ndis,ictr)=ni;
                                xct(1:ndis,ictr)=-2;
                                %        n(count,ictr)=nIRe-i*keIm;
                                %        xm(count,ictr)=-1;
                                
                            case 'Pt'
                                materi{ictr}=type;
                                [nIRe,keIm]=nPt(lambdaNot);
                                ni=nIRe-i*keIm;
                                %        nir(ictr)=ni;
                                nit(1:ndis,ictr)=ni;
                                xct(1:ndis,ictr)=-3;
                                %        n(count,ictr)=nIRe-i*keIm;
                                %        xm(count,ictr)=-1;
                                
                            otherwise
                                
                                ictr=ictr-1;
                                Word=type;
                                
                                fan=findstr(Word,':Lable=');
                                iseta=0;
                                if length(fan)>0
                                    pubar=find(Word==':');
                                    iseta=0;
                                    %        iseta=1;
                                    setarray=0;
                                    du=str2num(Word(pubar+7:end));
                                    pucount=count:count+ndis-1;
                                    ralabl(pucount,ictr)=-du;
                                    Word=Word(1:fan-1);
                                end
                                
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
                                
                                %       Word
                                %       'Word'
                                %       pausak
                                switch Word
                                    
                                    case {'circle','rhombus','ellipse','rectangle','square'}
                                        
                                        [word2,remainder]=strtok(remainder);
                                        pubar=find(word2=='|');
                                        if length(pubar)==1
                                            ictP=ictP+1;
                                            pucount=count:count+ndis-1;
                                            ipar(pucount,1,ictP)=str2num(word2(pubar+2:end));
                                            ipar(pucount,2,ictP)=2;
                                            ipar(pucount,3,ictP)=ictr;
                                            apert=word2(1:pubar-1);
                                        end
                                        
                                        radir(ictr)=str2num(apert);
                                        %'set radir', keyboard
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
                                            ipar(pucount,2,ictP)=3;
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
                                    case {'rhombus'}
                                        
                                        shavr(ictr)=4;
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
                                        
                                    case {'array'}
                                        if iseta==1
                                            setarray=2;
                                        end
                                        counces=pucount;
                                        shavr(ictr)=5;
                                        
                                    case {'grating'}
                                        
                                        counces=pucount;
                                        shavr(ictr)=6;
                                        
                                    case {'HCG'}
                                        
                                        counces=pucount;
                                        shavr(ictr)=6;
                                        
                                    case 'circle'
                                        radirb(ictr)=0;
                                        radirc(ictr)=0;
                                        shavr(ictr)=1;
                                        %        'circle', keyboard
                                    otherwise
                                        fi=strncmpi(Word,'Nmesh=',6)
                                        if fi==1
                                            DDmesh(count_lay,1)=str2num(Word(7:end));
                                            'DD Wor', keyboard
                                        else
                                            disp(' wrong type of transverse shape 1')
                                            keyboard
                                        end
                                        disp(' wrong type of transverse shape 2')
                                        keyboard
                                end
                        end
                        [type,remainder]=strtok(remainder);
                        
                    end  % while trasv.
                    
                    
                    
                    
                    icloc=0;
                    fiqw=find(fls(:,1)<0);
                    %    disp( ' layer '), keyboard
                    for counti=count:count+ndis-1
                        %    for counti=count:count+ndis
                        
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
%							'fls5', keyboard
                        if abs(nlong)>1
                            ipmem(counti)=counti;
                        end

                        d(counti,1)=thick/ndis;
                        %     d
                        %   [counti thick]
                        %   pausak
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
%											nlong=1;
                    %'ictr'
                    %keyboard
                    
                    if ictr>1
                        %'ictr'
                        %keyboard
                        ram=ones(ndis,ictr-1)*radir;
                        %     ram=ones(ndis,ictr)*radir;
                        radirb=ram;
                        radirc=ram;
                        %     sram=ones(ndis,1)*shavr;
                        %     nim=ones(ndis,1)*nir;
                        %     sram=ones(ndis,ictr)*shavr;
                        sram=ones(ndis,ictr-1)*shavr;
                        %     nim=ones(ndis,1)*nir;
                        %     nit=[nim nd];
                        %     nit=nir;
                        %     xct=[nim*0 xc];
                    else
                        %     sn=size(n);
                        nco=1;
                        ram=zeros(ndis,nco);
                        nim=zeros(ndis,nco);
                        radirb=ram;
                        radirc=ram;
                        radirar=ram;
                        sram=ram;
                        nit=[nd];
                        xct=[xc];
                    end
                    
                    %   'nit'
                    %   keyboard
                    
                    si=size(nit);
                    co=1:si(2);
                    if si(2)>1
                        cor=1:si(2)-1;
                    else
                        cor=co;
                    end
                    
                    n(count:counti,co)=nit;
                    xm(count:counti,co)=xct;
                    Dov(count:counti)=dovdu;
                    
                    if AlCont1==AlCont2
                        xw(countw,co)=AlCont1;
                        dw(countw,1)=thick;
                        dov(countw)=dovdu;
                        flw(countw)=nlong;
                        countw=countw+1;
                    else
                        if igrad==1
                            xw(countw,co)=AlCont1+j*AlCont2;
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
                    end
%					'qui ram', keyboard

                    radii.a(count:counti,cor)=ram;
                    radii.b(count:counti,cor)=radirb;
                    radii.c(count:counti,cor)=radirc;
                    %   radii.array(count:counti,cor)=radirar;
                    shav(count:counti,cor)=sram;
                    %   keyboard
                    %   n=[n; nit];
                    %   xm=[xm; xct];
                    %   radii=[radii; ram];
                    
                    %'set Lg', keyboard
					if count_lay==0
					count_lay=count_lay+1;
		            RaDD(count_lay,1)=ram(1);
					else
                    count_lay=count_lay+1;
					end	

                    DDd(count_lay,1)=d(count)*Ndis;
					
					RaDD(count_lay,1)=ram(1);
                    if itip==-4
                        DDd(count_lay,1)=-DDd(count_lay,1);
                    end
                    
                    %         if count==8
                    %          '8;\;', keyboard
                    %         end
%                                         'IPar 2600', keyboard     
                    if size(ipar,1)>=count
                        %              IPar(count_lay,:)=ipar(count,1:2,1);
                        fib=find(ipar(count,2,:)==2 | ipar(count,2,:)==-6 );
                        %              'fib', keyboard
                        IPar(count_lay,fib)=squeeze(ipar(count,1,fib))';
                        
                        %              if length(fib)==1
                        %              IPar(count_lay,:)=ipar(count,1:2,fib);
                        %              end
				    else
						 IPar(count_lay,:)=0;						
                    end
                   Pun_cont(count_lay)=count;
               
                    DDmat{count_lay,1}=materi;
                    
                    clear materi
                    DDx(count_lay,1)=AlCont1;
                    DDx(count_lay,2)=AlCont2;
                    if exist('dop')
                        if ~exist('do2')
                            do2=Dov(count);
                        end
                        if Dov(count)>0
                            dddra=[Dov(count) do2];
                            dddrn=[0 0];
                        else
                            dddrn=[Dov(count) do2];
                            dddra=[0 0];
                        end
%   'lay5', pausak						
                        clear do2
                        DDna(count_lay,:)=dddra;
                        DDnd(count_lay,:)=dddrn;
%count_lay
%                              'lay5', pausak
                    end
 
                    if length(fiM)==1
                        
                        DDmesh(count_lay,1)=Nmesh;
%                                                 'DD que', keyboard
                        fiM=[];
                    end
                    

                    clear AlCont AlCont1 AlCont2
                    
                    if setarray==0
                        count=counti+1;
                    end
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
                    
            end % switch typelay  (Lc o Lg)
            
    end % switch typelay (Array o layer)
    % if count>50
    % if count==56
    %  count, keyboard
    % end
    
end  %fid

%' qui ', keyboard

flsv=fliplr(fls);
fsw=flw';
if exist('dov')
    dov=dov';
else
    dov=zeros(size(d));
end
lf=length(ifield);


%' ifield ', keyboard

if lf~=count-1
    ifield(lf+1:count-1)=0;
end
ifield=ifield';



lf=length(Dov);
if lf~=count-1
    Dov(lf+1:count-1)=0;
end
Dov=Dov';

%if exist('radii')
% si=size(radii.array);
% lf=si(1);
% if length(si)==3
%  ltc=si(3);
% else
%  ltc=1;
% end
%else
% lf=0;
%end
%if lf~=count-1
% for ka=lf+1:count-1
%  for ktc=1:ltc
%   radii.array{ka,1,ltc}=0;
%  end
% end
%end

if exist('radii')
    si=size(radii.array);
    lf=si(1);
    if length(si)==2
        ltc=si(2);
    else
        ltc=1;
    end
else
    lf=0;
end
zecel{1}=0;
if lf~=count-1
    for ka=lf+1:count-1
        for ktc=1:ltc
            radii.array{ka,ktc}=zecel;
        end
    end
    for ka=1:lf
        for ktc=1:ltc
            if length(radii.array{ka,ktc})==0
                radii.array{ka,ktc}=zecel;
            end
        end
    end
end

if exist('radii')
    si=size(radii.a);
    lf=si(1);
else
    lf=0;
end
if lf~=count-1
'qi', keyboard
%    RaDD(lf)=0;
    for ka=lf+1:count-1
        radii.a(ka,1)=0;
        radii.b(ka,1)=0;
        radii.c(ka,1)=0;
        shav(ka,1)=0;
        dov(ka,1)=0;
    end
end

if exist('ralabl')
    si=size(ralabl);
    lf=si(1);
else
    lf=0;
end
if lf~=count-1
    for ka=lf+1:count-1
        ralabl(ka,1)=0;
    end
end
%' prima ipar', keyboard

if exist('ipar')
    spa=size(ipar);
    % if length(spa)>2
    %  lf=spa(1);
    %  if lf~=count-1
    %   for k=1:spa(3)
    %    ipar(lf+1:count-1,:,k)=0;
    %   end
    %  else
    %   ipar=zeros(count-1,3,3);
    %  end
    % end
    
    if length(spa)<=2
        spa(3)=2;
    end
    lf=spa(1);
    if lf~=count-1
        for k=1:spa(3)
            ipar(lf+1:count-1,:,k)=0;
        end
        %  else
        %   ipar=zeros(count-1,3,3);
    end
    
else
    ipar=zeros(count-1,3,3);
end

perd=-imag(n)*2*kNoll/100;

fclose(fid);

%disp('lay_new prima')
%keyboard


if exist('ralabl')
    filab=find(ralabl(:,1)>0);
    for ila=filab'
        dur=ipar(ila,:,:);
        sd=size(dur);
        filabn=find(ralabl==-ralabl(ila));
        % pausak
        radii.a(filabn)=radii.a(ila);
        radii.b(filabn)=radii.b(ila);
        radii.c(filabn)=radii.c(ila);
        for ilan=filabn'
            radii.array{ilan}=radii.array{ila};
            du=ipar(ilan,:,:);
            for kas=1:sd(3)
                du=ipar(ilan,:,kas);
                ik=find(du==0);
                if length(ik)==3
                    ipar(ilan,:,kas)=dur(1,:,kas);
                end
            end
        end
    end
end

if exist('rigaf')
    % 'controllo fattz', keyboard
    ds=d;
    dws=dw;
    d(rigaf:end)=d(rigaf:end)*fattz;
    dw(rigaf:end)=dw(rigaf:end)*fattz;
    lambdaNot=lambdaNot*fattz;
end

if exist('pames')
    radii.mesa=pames;
end

% [nref,aref,n,d,xm,Dov,radii,flsv,perd,anyf,iauto,...
%         dw,xw,fsw,dov,shav,ipar,ifield,lambdaNot]= Lay_new(fileName)
% Dov flsv anyf iauto ipar ifield
%  dw xw fsw dov : questi no!
%  ip1=reshape(ipar(:,1,:),sip(1)*sip(3),1);   % par. #
%  ip2=reshape(ipar(:,2,:),sip(1)*sip(3),1);   % par. type
%  ip3=reshape(ipar(:,3,:),sip(1)*sip(3),1);   % trans. sect.

fiau=find(iauto(:,2)==-4);
if length(fiau)>4
    iauto(fiau(2:end-1),2)=0;
end
cocrit=3;
if icrit>1
    fiatti=find(iauto(:,1)==2);
    if max(riga_crit)>fiatti
        %  'ICRIT solo parte top; sotto non ancora predisposto'
        %  keyboard
        %  riga_crit=fiatti;
    end
    for k=riga_crit
        iauto(k-1,cocrit)=-10;
    end
end


%keyboard
if length(find(iauto(:,1)==1))==0
    ' upper investigation limit is missing '
    ' correct the .str file !!!!! '
    Nline
    NL
    keyboard
end

ParOpt.nref=nref;
ParOpt.aref=aref;
ParOpt.n=n;
ParOpt.d=d;
ParOpt.xm=xm;
ParOpt.Dov=Dov;
ParOpt.radii=radii;
ParOpt.flsv=flsv;
ParOpt.perd=perd;
ParOpt.anyf=anyf;
ParOpt.iauto=iauto;
ParOpt.dw=dw;
ParOpt.xw=xw;
ParOpt.fsw=fsw;
ParOpt.dov=dov;
ParOpt.shav=shav;
ParOpt.ipar=ipar;
ParOpt.ifield=ifield;
ParOpt.lambdaNot=lambdaNot;
ParOpt.igrad=igrad;

puDD=1:length(DDd);
%fiDD=find(Pun_cont==iDD);
%puDD=Pun_cont(fiDD:end);


if length(iDD)>0
     if length(iDD)==1
%         if iDD<=length(DDd)/2
           puDD=puDD(iDD:end);
%         else
%             puDD=puDD(1:iDD);
%         end
     else
         puDD=iDD(1):iDD(2);
     end
end

if length(RaDD)<length(DDna)
 RaDD=[RaDD; zeros(length(DDna)-length(RaDD),1)];
end
%'qui iDD xm', keyboard

ParDD.d=DDd(puDD)*fattz;
IPar=IPar(puDD,:);
ParDD.Na=DDna(puDD,:);
DDnd=abs(DDnd(puDD,:));
ParDD.Nd=abs(DDnd);
ParDD.x=DDx(puDD,:);

rag=radii.a;
purag=Pun_cont(puDD); 
%'DDrag', keyboard

% ParDD.xm=xm(purag,2:end);
ParDD.xm=xm(purag,:);
%ragTJ=rag;
if RaDD(1)~0
%RaDD=[zeros(size(RaDD(1,:))); RaDD];
end
ragTJ=RaDD;

ragTJ(iETCHdd(1:end+paTJ.Nlay_ethc)-1,1)=paTJ.Ram;

%'qui MAT DD rep', keyboard
%DDrag=ragTJ(purag,:);
DDrag=ragTJ(puDD,:);
ParDD.rag=DDrag;

 
%rafake=[];
%if purag(1)>1
% pupri=1:purag(1)-1;
% rafaked=rag(pupri);
% fi=find(rafaked>0);
% rafake=rafaked(fi);
%end
%ParDD.ragAdd=rafake';
ParDD.ragAdd=ragAdd;

DDrep=fls(purag,:);
DDrep0=DDrep;

dd=DDrep(:,2);
fi=find(dd>1);
dd=DDrep(:,2);
fit=find(diff(fi)>1 | diff(dd(fi))~=0);
if length(fi)>0
    fin=fi(1);
    for k=1:length(fit)
        fu=fi(fit(k));
        pu=[fin:fu];
        DDrep(pu,2)=length(pu);
        fin=fi(fit(k)+1);
    end
    fu=fi(end);
    pu=[fin:fu];
    DDrep(pu,2)=length(pu);
end

fia500=find(iauto(purag,2)==-500);
DDrep(fia500,1)=-500;

ParDD.rep=DDrep;
ParDD.mesh=DDmesh(puDD) ;



DDdu=DDmat(puDD,:) ;


for k=1:length(DDdu)
    dch=DDdu{k};
    dchn=dch;
    l=length(dch);
    for ki=1:l
        if strcmp(char(dch(ki)),'Ref_1')
            dchn{ki}='vacuum';
        elseif strcmp(char(dch(ki)),'Ref_1.6')
            dchn{ki}='AlOx';
        elseif strcmp(char(dch(ki)),'Au')
            dchn{ki}='Line';
        elseif strcmp(char(dch(ki)),'Ti')
            dchn{ki}='Line';
        elseif strcmp(char(dch(ki)),'Pt')
            dchn{ki}='Line';
        elseif strcmp(char(dch(ki)),'AlGaAs')
            %   dchn{ki}='AlOx';
        elseif strcmp(char(dch(ki)),'InGaAs')
            %   dchn{ki}='AlOx';
        else
            dchn{ki}='vacuum';
            %   'quiii', keyboard
        end
    end
    DDdu{k}=dchn;
end

ParDD.mat=DDdu;
ParDD.iDD=iDD;

%IPard=IPar;

%IPard=squeeze(ipar(:,1,1:2));


%IPar=zeros(length(puDD),2);
%if length(iDD)==0
%    iDD=1;
%end


fisiz=find(ipar0(1,2,:)==1);
IPdu=IPar;
fiMa=find(IPar(:,1)>0);
%'Ipar; PRIMAAAAAAAA', keyboard
if length(fisiz)==1
IPar(1:iETCHdd,1)=1;
end
IPar(fiMa,1)=IPdu(fiMa,1);
%'Ipar;', keyboard
%IPar(1:length(IPard(iDD+1:end,:)),1:size(IPard,2))=IPard(iDD+1:end,:);
ParDD.IPar=IPar;
if ~exist('ParMore')
    ParMore.r_passiv=[];
    ParMore.mesh_passiv=[];
end

%'par;', keyboard

ParMore.Bcond=Bcond0;
ParMore.mesa=mesa;
ParMore.Buffer=Buffer;
ParMore.contact=contact;
if ~exist('mireq')
    mireq.p.P_i=0;
    mireq.p.P_f=0;
    mireq.p.x_eq=[];
    mireq.p.Dop_eq=[];
    mireq.p.Nmesh=[];
    mireq.p.Mateq='niente';
    
    mireq.n.P_i=0;
    mireq.n.P_f=0;
    mireq.n.x_eq=[];
    mireq.n.Dop_eq=[];
    mireq.n.Nmesh=[];
    mireq.n.Mateq='niente';
    
    mireq.n2.P_i=0;
    mireq.n2.P_f=0;
    mireq.n2.x_eq=[];
    mireq.n2.Dop_eq=[];
    mireq.n2.Nmesh=[];
    mireq.n2.Mateq='niente';
    
    mireq.exist=0;
else
    mireq.exist=1;
end
if abs(max(fls(:,1))>1)
    mireq.exist=1;
end
ParMore.mireq=mireq;
% disp('fine lay_tran stop')
% pausak
%  'Laydd', keyboard

if exist('flsKazzo')
    fiDD=find(flsKazzo(:,2)==-1);
    if length(fiDD)>2
        fid=[2:fiDD(1)-1];
        dz_Term=sum(flsv(fid,2).*d(fid));
    else
        dz_Term=0;
    end
    ParMore.dz_Term=dz_Term;
else
    ParMore.dz_Term=0;
end
%iparDD=ipar(iDD:end,:,:);

if length(iNEGF)>0
%     if length(find(diff(iNEGF)~=1))>0
%         'error NEGF not contiguous', keyboard
%     end
%     
%     pu=1:iNEGF(1)-1;
%     
%     %  flsN=fls(pu,1);
%     %  fiM=find(flsN>1);
%     %  if length(fiM)>0
%     %   nsAdNegf=flsN(fiM)*length(fiM)-length(fiM);
%     %  else
%     %   nsAdNegf=0;
%     %  end
%     
%     ParDD.iNEGF=iNEGF-iDD(1)+1;
%     
%     %  indNEGF = length(iNEGF)-1;
%     %  iNEGF = length(ParDD.d)-1-iNEGF(1);
%     %  ParDD.iNEGF = [iNEGF:(iNEGF+indNEGF)];
%     
%     fiCav=find(iauto(:,2)==-4);
%     
%     if length(fiCav)==0
%         le=length(iauto);
%         pu=fix(le/3*[1 2]);
%         iauto(pu,2)=-4;
%     end
    if length(find(diff(iNEGF)~=1))>0
        
        % Modifications by Torrelli Valerio, 02/11/2021 -> reading of
        % multiple tunnel junctions
        
        % Within iNEGF I have all of the information I need to address the
        % problem
        
        % 1) I need to count the number of BTJs
        nBTJ=length(find(diff(iNEGF)~=1))+1;
        
        % 2) Now, I need to put the different indices in iNEGF in different
        % vectors associated to each BTJ. To do so, a cell formalism seems
        % to be convenient.
        InitialPosBTJ=[1,find(diff(iNEGF)~=1)+1];
        
        for ii = 1:nBTJ
            if ii == nBTJ
                iNEGF_multiple{ii}=iNEGF(InitialPosBTJ(ii):end)-iDD(1)+1;
            else
                iNEGF_multiple{ii}=iNEGF(InitialPosBTJ(ii):(InitialPosBTJ(ii+1)-1))-iDD(1)+1;
            end
        end
        
        % Now, the cell iNEGF_multiple contains many vectors, corresponding
        % to the BTJ layers in the VCSEL
        
        % I'm saving the most important quantities I derived
        ParDD.iNEGF_multiple = iNEGF_multiple;
        ParDD.nBTJ = nBTJ;      
        
    else
        nBTJ = 1;
        iNEGF_multiple = {iNEGF-iDD(1)+1};
        
        ParDD.iNEGF_multiple = iNEGF_multiple;
        ParDD.nBTJ=nBTJ;
        
    end
	if exist('OxVero')
    	ParDD.OxVero=OxVero;
	end
    
    pu=1:iNEGF(1)-1;
    
    %  flsN=fls(pu,1);
    %  fiM=find(flsN>1);
    %  if length(fiM)>0
    %   nsAdNegf=flsN(fiM)*length(fiM)-length(fiM);
    %  else
    %   nsAdNegf=0;
    %  end
    
    
%     ParDD.iNEGF=iNEGF-iDD(1)+1;  % what is the point of this -1+1? Keep in mind that iDD = 1.
    
    
    %  indNEGF = length(iNEGF)-1;
    %  iNEGF = length(ParDD.d)-1-iNEGF(1);
    %  ParDD.iNEGF = [iNEGF:(iNEGF+indNEGF)];
    
    fiCav=find(iauto(:,2)==-4);
    
    if length(fiCav)==0
        le=length(iauto);
        pu=fix(le/3*[1 2]);
        iauto(pu,2)=-4;
    end
end

ParDD.iNEGFdd=iETCHdd-iDD;
ParOpt.iauto=iauto;
  %'Laydd Litho end', keyboard