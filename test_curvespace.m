 %% Test Curvespace 3D
  x = -2*pi:0.5:2*pi;
  y = 10*sin(x);
  z = linspace(0,10,length(x));
  N = 50;
  p = [x',y',z'];
  q = curvspace(p,N);
  figure;
  plot3(p(:,1),p(:,2),p(:,3),'*b',q(:,1),q(:,2),q(:,3),'.r');
  axis equal;
  legend('Original Points','Interpolated Points');

%   See also LINSPACE.


 %% Test Curvespace 2D
%   x = -2*pi:0.5:2*pi;
%   y = 10*sin(x);
tic
hold on
xleft = xbodyL
yleft = ybodyL
%   z = linspace(0,10,length(x));
  N = length(xleft);
  pl = [xleft',yleft'];
  ql = curvspace(pl,N);
  figure(1);
  plot(pl(:,1),pl(:,2),'*b',ql(:,1),ql(:,2),'.r');
%   axis equal;
  legend('Original Points','Interpolated Points');
  
xbodyB; ybodyB;
hold on
xbot = xbodyB
ybot = ybodyB
%   z = linspace(0,10,length(x));
  N = length(xbot);
  pbot = [xbot',ybot'];
  qbot = curvspace(pbot,N);
  figure(1);
  plot(pbot(:,1),pbot(:,2),'*b',qbot(:,1),qbot(:,2),'.r');
%   axis equal;
  legend('Original Points','Interpolated Points');
  


xright = xbodyR
yright = ybodyR
%   z = linspace(0,10,length(x));
  N = length(xright);
  pr = [xright',yright'];
  qr = curvspace(pr,N);
  figure(1);
  plot(pr(:,1),pr(:,2),'*b',qr(:,1),qr(:,2),'.r');
%   axis equal;
  legend('Original Points','Interpolated Points');
hold on
toc
%% Next iter

ybodyL(1) = Amp;
ybodyR(end) = -Amp;
 hold on
xleft = xbodyL
yleft = ybodyL
%   z = linspace(0,10,length(x));
  N = length(xleft);
  pl = [xleft',yleft'];
  ql = curvspace(pl,N);
  figure(1);
  plot(pl(:,1),pl(:,2),'*b',ql(:,1),ql(:,2),'.r');
%   axis equal;
  legend('Original Points','Interpolated Points');
  
xbodyB; ybodyB;
hold on
xbot = xbodyB
ybot = ybodyB
%   z = linspace(0,10,length(x));
  N = length(xbot);
  pbot = [xbot',ybot'];
  qbot = curvspace(pbot,N);
  figure(1);
  plot(pbot(:,1),pbot(:,2),'*b',qbot(:,1),qbot(:,2),'.r');
%   axis equal;
  legend('Original Points','Interpolated Points');
  


xright = xbodyR
yright = ybodyR
%   z = linspace(0,10,length(x));
  N = length(xright);
  pr = [xright',yright'];
  qr = curvspace(pr,N);
  figure(1);
  plot(pr(:,1),pr(:,2),'*b',qr(:,1),qr(:,2),'.r');
%   axis equal;
  legend('Original Points','Interpolated Points');
hold on