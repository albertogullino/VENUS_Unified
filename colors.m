% Change the colors and their order in the plots

% Rainbow 
newcolors = {'#F00','#F80','#FF0','#0B0','#00F','#50F','#A0F'};

% rgbcmykw
newcolors = {'r','g','b','c','m','y','k','w'};

% Tango Project
newcolors = {'#edd400','#f57900','#c17d11','#73d216','#3465a4','#75507b','#cc0000','#d3d7cf','#555753'};

% Tango Project (light)
newcolors = {'#fce94f','#fcaf3e','#e9b96e','#8ae234','#729fcf','#ad7fa8','#ef2929','#eeeeec','#888a85'};

% colororder help (very good for 4-curve plots)
newcolors = [0.83 0.14 0.14
             1.00 0.54 0.00
             0.47 0.25 0.80
             0.25 0.80 0.54];

% 3 colors plots
newcolors = {'#ffb60a','m','c'};

% Default Matlab
newcolors = {'#0072BD','#D95319','#EDB120','#7E2F8E','#77AC30','#4DBEEE','#A2142F','#eeeeec','#888a85'};
         
colororder(newcolors)