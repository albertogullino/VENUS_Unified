% The function n_distr plots layers position.
% #####################################################################
% Each layer can be a set of rings with different radii and refractive inexes.
% The function requires the presence of the function bdfCallback(), 
% which provides the possibility of representation of the data for each
% segment on click.
% #####################################################################
% 
% Input data:
% Litot - vector of layer widths
% 
% aitot(i,j) - matrix of radii: i - number of layer, j - number of ring
% 
% nitot(i,j) - matrix of refractive indexes: i - number of layer, j -number of ring
% 
% n_type - flag for representation of different parts of the refractive
% index:
% 'real'  - will be shown real part of the refractive index;
% 'imag' - will be shown real part of the refractive index;
% 'abs' - will be calculated abs for  the refractive index
% 
% coord_num - is a key, that allow to show numbers of the layers on the axis instead of the coordinate.
% if coord_num defined as 'num' it will be layer number on the axis, and just coordinate otherwise.
% 
% vert - flag for rotated representation of the structure, if vert=='v' ,
% the structure will be plotted vertically, and otherwise - horizontally

function n_distr(Litot, aitot, nitot,n_type,coord_num, vert,Ngrat,par_grat)

global currPoint position click_counter t1

% loop to remove central rings from the consideration, if their radius=0
for ii=1:(length(aitot(1,:))-1)
    skip=find(((aitot(:,ii)==0))&(sum(aitot(:,ii:end),2)~=0));
    aitot(skip,ii)=aitot(skip,ii+1);
    aitot(skip,ii+1)=0;
    nitot(skip,ii:end-1)=nitot(skip,ii+1:end);
end

click_counter=0;
clear t1;
position=vert;
if strcmp(n_type,'real')
    nitot=real(nitot);
elseif strcmp(n_type,'imag')
    nitot=imag(nitot);
elseif strcmp(n_type,'abs')
    nitot=abs(nitot);
end
R_max=max(max(aitot))+2; % radius of the structure
aitot=[aitot zeros(2,length(aitot))'];
X=[0;cumsum(Litot)];
figure
hold on
for ii=1:length(Litot) %loop for each layer
    if Litot(ii)~=0
        for jj=1:length(nitot(1,:)) % loop for each ring within layer
            n=[nitot(ii,jj) nitot(ii,jj); nitot(ii,jj) nitot(ii,jj)];
            if jj==1 %central part
                x=[X(ii) X(ii);X(ii+1) X(ii+1)];
                if (aitot(ii,jj)==0) %just one ring in the layer
                    y=[-R_max R_max;-R_max R_max];
                else % not the last ring
                    y=[-aitot(ii,jj) aitot(ii,jj);-aitot(ii,jj) aitot(ii,jj)];
                end
                if strcmp(vert,'v')
                    pcolor(y,-x,n);
                else
                    pcolor(x,y,n);
                end
            elseif aitot(ii,jj-1)~=0 %side parts if not zero
                if sum(aitot(ii,jj+1:end))==0 %last ring
                    y1=[(aitot(ii,jj-1)) R_max; (aitot(ii,jj-1)) R_max];
                    y2=[-R_max -(aitot(ii,jj-1));-R_max -(aitot(ii,jj-1))];
               else %not the last ring
                    y1=[(aitot(ii,jj-1)) (aitot(ii,jj)); (aitot(ii,jj-1)) (aitot(ii,jj))];
                    y2=[-(aitot(ii,jj)) -(aitot(ii,jj-1));-(aitot(ii,jj)) -(aitot(ii,jj-1))];
               end
               if strcmp(vert,'v')
                   pcolor(y1,-x,n);
                   pcolor(y2,-x,n);
               else
                   pcolor(x,y1,n);
                   pcolor(x,y2,n);
               end

           end
        end
    end
end

non_zero=find(Litot>0);
if strcmp(coord_num,'num')
    for ii=1:length(Litot)

        if strcmp(vert,'v')
           tick_label(ii)={num2str(length(Litot)-ii+1)};
        else
           tick_label(ii)={num2str(ii)};
        end
    end
    
    tick_coord=cumsum([0;Litot(non_zero)]+[Litot(non_zero);0])/2;
    if strcmp(vert,'v')
        set(gca,'YTick',-tick_coord(end-1:-1:1));
        set(gca,'YTickLabel',{tick_label{length(Litot)+1-non_zero(end:-1:1)}});  
    else
        set(gca,'XTick',tick_coord(1:end-1));
        set(gca,'XTickLabel',tick_label(non_zero));
    end
end
colorbar
set(gcf,'WindowButtonDownFcn',{@bdfCallback},'Units','normalized');
axis([-R_max R_max -sum(Litot) 0])
caxis([1 max(max(nitot))])