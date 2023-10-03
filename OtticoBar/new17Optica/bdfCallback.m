function bdfCallback(hFig,event,src)
global currPoint t1 filename position click_counter
load(filename)
if click_counter~=0
 delete(t1)
else 
    click_counter=1;
end

aitot(find(aitot==0))=max(max(aitot))+2;

currPoint = get(gca,'CurrentPoint');
currPoint(:,3) = [];
currPoint(2,:) = [];


 
% find the data for the layer
if strcmp(position,'v')
    curr_x=currPoint(1); %along radius
    curr_y=abs(currPoint(2));%along the structure
else
    curr_x=currPoint(2);%along radius
    curr_y=abs(currPoint(1));%along the structure
end
coord_x=cumsum([0;Litot]);
layer_num=min(find(curr_y<coord_x))-1;
clear coord_x
if max(aitot(layer_num,:))==0
    ring_num=1;
else
    ring_num=min(find(abs(curr_x)<aitot(layer_num,:)));
end
%ring_num
n=nitot(layer_num,ring_num);
w=Litot(layer_num);
line1=['{\bfLayer #',num2str(layer_num),'}'];
line2=['{\bfn=',num2str(n),'}'];
line3=['{\bfw=',num2str(w),'}'];
t1=text(currPoint(1),(currPoint(2)),{line1;line2;line3},'Color','black','EdgeColor','black',...
'LineWidth',3);
curr_x;
curr_y;
end

