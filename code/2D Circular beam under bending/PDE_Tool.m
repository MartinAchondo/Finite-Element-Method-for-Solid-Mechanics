% This script is written and read by pdetool and should NOT be edited.
% There are two recommended alternatives:
% 1) Export the required variables from pdetool and create a MATLAB script
%    to perform operations on these.
% 2) Define the problem completely using a MATLAB script. See
%    https://www.mathworks.com/help/pde/examples.html for examples
%    of this approach.
function pdemodel
[pde_fig,ax]=pdeinit;
pdetool('appl_cb',3);
pdetool('snapon','on');
set(ax,'DataAspectRatio',[1 1.5 1]);
set(ax,'PlotBoxAspectRatio',[6.5 4.333333333333333 1]);
set(ax,'XLim',[-1 12]);
set(ax,'YLim',[-1 12]);
set(ax,'XTick',[ -12,...
 -11,...
 -10,...
 -9,...
 -8,...
 -7,...
 -6,...
 -5,...
 -4,...
 -3,...
 -2,...
 -1,...
 0,...
 1,...
 2,...
 3,...
 4,...
 5,...
 6,...
 7,...
 8,...
 9,...
 10,...
 11,...
 12,...
]);
set(ax,'YTick',[ -12,...
 -11,...
 -10,...
 -9,...
 -8,...
 -7,...
 -6,...
 -5,...
 -4,...
 -3,...
 -2,...
 -1,...
 0,...
 1,...
 2,...
 3,...
 4,...
 5,...
 6,...
 7,...
 8,...
 9,...
 10,...
 11,...
 12,...
]);
pdetool('gridon','on');

% Geometry description:
pdeellip(0,0,5,5,...
0,'E1');
pdeellip(0,0,10,10,...
0,'E2');
pderect([0 -10 10 -10],'R1');
pderect([10 0 0 -10],'SQ1');
set(findobj(get(pde_fig,'Children'),'Tag','PDEEval'),'String','E2-E1-R1-SQ1')

% Boundary conditions:
pdetool('changemode',0)
pdesetbd(4,...
'dir',...
2,...
char('0','0','0','0'),...
char('0','0'))
pdesetbd(3,...
'dir',...
2,...
char('0','0','0','0'),...
char('0','0'))
pdesetbd(2,...
'dir',...
2,...
char('1','0','0','1'),...
char('0','0'))
pdesetbd(1,...
'dir',...
2,...
char('1','0','0','0'),...
char('-0.01','0'))

% Mesh generation:
setappdata(pde_fig,'Hgrad',1.1000000000000001);
setappdata(pde_fig,'refinemethod','regular');
setappdata(pde_fig,'jiggle',char('on','mean',''));
setappdata(pde_fig,'MesherVersion','preR2013a');
pdetool('initmesh')

% PDE coefficients:
pdeseteq(1,...
char('2*((1E4)./(2*(1+(0.25))))+(2*((1E4)./(2*(1+(0.25)))).*(0.25)./(1-(0.25)))','0','(1E4)./(2*(1+(0.25)))','0','(1E4)./(2*(1+(0.25)))','2*((1E4)./(2*(1+(0.25)))).*(0.25)./(1-(0.25))','0','(1E4)./(2*(1+(0.25)))','0','2*((1E4)./(2*(1+(0.25))))+(2*((1E4)./(2*(1+(0.25)))).*(0.25)./(1-(0.25)))'),...
char('0.0','0.0','0.0','0.0'),...
char('0.0','0.0'),...
char('1.0','0','0','1.0'),...
'0:10',...
'0.0',...
'0.0',...
'[0 100]')
setappdata(pde_fig,'currparam',...
['1E4 ';...
'0.25';...
'0.0 ';...
'0.0 ';...
'1.0 '])

% Solve parameters:
setappdata(pde_fig,'solveparam',...
char('0','1000','10','pdeadworst',...
'0.5','longest','0','1E-4','','fixed','Inf'))

% Plotflags and user data strings:
setappdata(pde_fig,'plotflags',[18 1 1 1 1 1 1 1 0 0 0 1 1 0 0 0 1 1]);
setappdata(pde_fig,'colstring','');
setappdata(pde_fig,'arrowstring','');
setappdata(pde_fig,'deformstring','');
setappdata(pde_fig,'heightstring','');