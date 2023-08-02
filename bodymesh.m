%% bodymesh
function [xbody,xbodyL,xbodyB,xbodyR,ybody,ybodyL,ybodyB,ybodyR] = bodymesh(xfl,etaL,xfr,etaR,xbodyL,xbodyB,xbodyR,ybodyL,ybodyB,ybodyR)
startpointY = etaL(end);
startpointX = xfl(end);

finishpointY = etaR(1);
finishpointX = xfr(1);

hold on
xleft = [startpointX xbodyL];
yleft = [startpointY ybodyL];

%   z = linspace(0,10,length(x));
  N = length(xleft);
  pl = [xleft',yleft'];
  ql = curvspace(pl,N);
%   figure(1);
%   plot(pl(:,1),pl(:,2),'*b',ql(:,1),ql(:,2),'.r');
%   axis equal;
%   legend('Original Points','Interpolated Points');
  
xbodyB; ybodyB;
hold on
xbot = xbodyB;
ybot = ybodyB;
%   z = linspace(0,10,length(x));
  N = length(xbot);
  pbot = [xbot',ybot'];
  qbot = curvspace(pbot,N);
%   figure(1);
%   plot(pbot(:,1),pbot(:,2),'*b',qbot(:,1),qbot(:,2),'.r');
%   axis equal;
%   legend('Original Points','Interpolated Points');
  



xright = [ xbodyR finishpointX] ;
yright = [ ybodyR finishpointY] ;

%   z = linspace(0,10,length(x));
  N = length(xright);
  pr = [xright',yright'];
  qr = curvspace(pr,N);
%   figure(1);
% %   plot(pr(:,1),pr(:,2),'*b',qr(:,1),qr(:,2),'.r');
%   axis equal;
%   legend('Original Points','Interpolated Points');
% hold on

xbodyL = ql(2:end,1)'; ybodyL =  ql(2:end,2)';
xbodyB = qbot(:,1)'; ybodyB =  qbot(:,2)';
xbodyR = qr(1:end-1,1)'; ybodyR =  qr(1:end-1,2)';
% 
% %%
xbody = [xbodyL xbodyB  xbodyR];
ybody = [ybodyL ybodyB  ybodyR];
% 
