function npar=gest_opt(nomeF,ifp)
fip=find(nomeF=='.');
nomeP=[nomeF(1:fip) 'par'];
nomePsa=[nomeF(1:fip-1),'par'];

DR=dir;
lD=length(DR);
fl_find=0;
for kD=3:lD
    dun=getfield(DR(kD),'name');
    if strcmpi(dun,nomeP)==1
        fl_find=1;
        break
    end
end
if fl_find==1
    if ifp~=-4
        isovr=input(' Vuoi sovrascrivere file ? [1]');
        if length(isovr)==0
            isovr=0;
        end
    else
        isovr=0;
    end
    
    if isovr~=1
        nomeV=nomeP;
        nomeP='fileDUM.par';
    end
else
    isovr=1;
end

% iold=1;
iold=0;

if isovr==1
    if iold==1
        [nref,aref,nv0,dv,xm,Dov,radii,fst,perd,anyf,iauto,dvw,xmw,fsw,dov,shavet,...
            ipar,ifield,lambda0]=Lay_new(nomeF);
    else
        [nref,aref,nv0,dv,xm,Dov,radii,fst,perd,anyf,iauto,dvw,xmw,fsw,dov,shavet,...
            ipar,ifield,lambda0]=Lay_tran(nomeF);
    end
    
    fidfile = fopen(nomeP,'w');
    
    if length(ipar)~=1
        sip=size(ipar);
        ip1=reshape(ipar(:,1,:),sip(1)*sip(3),1);
        if sip(2)>=2
            ip2=reshape(ipar(:,2,:),sip(1)*sip(3),1);
        else
            ip2=ip1;
        end
        if sip(2)>=3
            ip3=reshape(ipar(:,3,:),sip(1)*sip(3),1);
        else
            ip3=ip1;
        end
        
        riga=[1:sip(1)]';
        righe=riga*ones(1,sip(3));
        ip4=reshape(righe,sip(1)*sip(3),1);
        fi1=find(ip1~=0);
        [dup,fi1p]=sort(ip1(fi1));
        fi1=fi1(fi1p);
        
        fid=find(diff([0; dup])~=0);
        fiud=fi1(fid);
        [ip11,fiudd]=sort(ip1(fiud));
        fiu=fiud(fiudd);
        
        npar=length(find(diff([0; dup])~=0));
        
        type(1,:)='thickness    ';
        type(2,:)='y-axis       ';
        type(3,:)='x-axis       ';
        type(4,:)='shape param  ';
        type(5,:)='ref. index   ';
        type(6,:)='array centers';
        type(7,:)='shape        ';
        type(8,:)='pair number  ';
        type(9,:)='rad. cir. gr.';
        
        typem(1,:)='thi';
        typem(2,:)='yax';
        typem(3,:)='xax';
        typem(4,:)='spa';
        typem(5,:)='ind';
        typem(6,:)='cen';
        typem(7,:)='sha';
        typem(8,:)='#pa';
        typem(9,:)='rcg';
        
        typeg(2,:)='grat orient  ';
        typeg(3,:)='grat shift   ';
        typeg(4,:)='cent cir rad.';
        typeg(5,:)='grat peach   ';
        typeg(6,:)='grat width   ';
        typeg(7,:)='grat este. y ';
        typeg(8,:)='grat este. x ';
        typeg(9,:)='grat est. sha';
        
        typea(2,:)='number x pixels';
        typea(3,:)='number y pixels';
        typea(4,:)='pixel shape    ';
        typea(5,:)='pixel Ry       ';
        typea(6,:)='pixel Rx       ';
        typea(7,:)='pixel Delta    ';
        typea(8,:)='pixel spacing  ';
        
        typegm(2,:)='ori';
        typegm(3,:)='shi';
        typegm(4,:)='ccr';
        typegm(5,:)='per';
        typegm(6,:)='wid';
        typegm(7,:)='esy';
        typegm(8,:)='esx';
        typegm(9,:)='esh';
        
        typeam(2,:)='#px';
        typeam(3,:)='#py';
        typeam(4,:)='sha';
        typeam(5,:)='pRy';
        typeam(6,:)='pRx';
        typeam(7,:)='pDe';
        typeam(8,:)='spa';
        
        shad=shavet;
        
        fis=find(shavet==0);
        shad(fis)=20;
        shad=[shad ones(size(shad(:,1)))*20];
        
        types(1,:)='  circle  ';
        types(2,:)=' rectangle';
        types(3,:)='  ellipse ';
        types(4,:)='  rhombus ';
        types(5,:)='   array  ';
        types(6,:)='  grating ';
        types(8,:)='    doe   ';
        types(20,:)='  planar  ';
        
        typesm(1,:)='ci';
        typesm(2,:)='re';
        typesm(3,:)='el';
        typesm(4,:)='rh';
        typesm(5,:)='ar';
        typesm(6,:)='gr';
        typesm(8,:)='do';
        typesm(20,:)='pl';
        
        
        
        
        nul='';
        
        % eval(['fprintf(fidfile,''' nul '\r\n'');']);
        prp=(['      Parameters=',num2str(npar),' in par_in: PARV (active parameter marked by at) must have the same order as the listed parameters ']);
        eval(['fprintf(fidfile,''' prp '\r\n'');']);
        if npar>0
            eval(['fprintf(fidfile,''' nul '\r\n'');']);
            prp=(['      Par,  Layer,   Tran_sect,   Par. type,       shape,       short_name']);
            eval(['fprintf(fidfile,''' prp '\r\n'');']);
            eval(['fprintf(fidfile,''' nul '\r\n'');']);
        end
        tabm='      ';
        tab='       ';
        tab2='  ';
        tab1='   ';
        tab0='';
        
        fiu_sav=fiu;
        [du,puparam]=sort(ip4(fiu));
        icsa=1;
        % for k=1:length(fiu)
        % ' puparam', keyboard
        
        for k=puparam'
            
            if length( find( ip2(fiu(k))<0))>0
                fis=find(shad(ip4(fiu(k)),:)==6);
                if length(fis)>0
                    TYPE=typeg(abs(ip2(fiu(k))),:);
                    %   TYPEm=TYPE;
                    TYPEm=typegm(abs(ip2(fiu(k))),:);
                end
                fis=find(shad(ip4(fiu(k)),:)==5);
                if length(fis)>0
                    TYPE=typea(abs(ip2(fiu(k))),:);
                    TYPEm=typeam(abs(ip2(fiu(k))),:);
                    %   TYPEm=TYPE;
                end
            else
                TYPE=type(ip2(fiu(k)),:);
                TYPEm=typem(ip2(fiu(k)),:);
            end
            %type(ip2(fiu(k)),:)
            if ip1(fiu(k))>9
                ta1=tabm;
            else
                ta1=tab;
            end
            if ip4(fiu(k))>9
                ta2=tabm;
            else
                ta2=tab;
            end
            if ip4(fiu(k))<=9
                ch1=['0',num2str(ip4(fiu(k)))];
            else
                ch1=[num2str(ip4(fiu(k)))];
            end
            ch2=num2str(ip3(fiu(k)));
            if iold==1
                sac=[TYPEm,ch2,typesm(shad(ip4(fiu(k)),ip3(fiu(k))),:),ch1];
            else
                shada=abs(shad(ip4(fiu(k)),ip3(fiu(k))));
                sac=[TYPEm,ch2,typesm(shada,:),ch1];
            end
            sacm(icsa,:)=sac;
            sapm(icsa,:)=k;
            icsa=icsa+1;
            
            chio=[tab1,'@ '];
            if iold==1
                prp=([ta1 num2str(ip1(fiu(k))) ta2 tab1  num2str(ip4(fiu(k))) tab ...
                    num2str(ip3(fiu(k)))  tab TYPE ...
                    tab2 types( shad( ip4(fiu(k)),ip3(fiu(k)) ),:),chio,sac]);
            else
                prp=([ta1 num2str(ip1(fiu(k))) ta2 tab1  num2str(ip4(fiu(k))) tab ...
                    num2str(ip3(fiu(k)))  tab TYPE ...
                    tab2 types( shada,:),chio,sac]);
            end
            
            %  k
            %  ' prima di scrivere ', pausak
            
            eval(['fprintf(fidfile,''' prp '\r\n'');']);
            
        end
        fiu=fi1;
        
    else
        prp=(['     NO Parameters to be supplied in par_in ']);
        eval(['fprintf(fidfile,''' prp '\r\n'');']);
        sacm=[];
        sapm=[];
        
    end
    
    eval(['fprintf(fidfile,''' nul '\r\n'');']);
    
    fclose(fidfile);
    
else
    
    fidfile = fopen(nomeV,'r');
    
    Nline=fgetl(fidfile);
    [pa,remainder]=strtok(Nline);
    npar=str2num(pa(12:end));
    fclose(fidfile);
end  %isovr
%'fine gest_opt', keyboard
