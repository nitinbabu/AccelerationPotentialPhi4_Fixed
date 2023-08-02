
clc
close all
clear
tic

% Works for these combinations|| not sure why                                      
% especially when fx is 20/40/60/80/100
% fx = 40; L = 100; h = 20;
% fx = 60; L = 300; h = 30;
% fx = 40; L = 300; h = 60;


fx =60*1; L = 6 ;h = 1;
time = 0%:0.1:10;
%%  eta = 0.5.*sin(0.09.*xf-1*t); Lambda = 70; T = 6;

%% Code breaks if there is some mismatch between designated elements and
%% actual elements due to positions
for m = 1
    t = time(m);
[x,y,xf,xd,xb,xp,yf,uin,uout,pot] = boundary(fx,h,L,t); %dx changed

pot = 5; uin = 0; uout = 0;
[ktvec,knvec]=bcs(xf, xd, xb, xp,pot,uout,uin);
[C1, C2, C3, C4] = corners(xf,xd,xb,xp);
nnode = length(x);  nelem  = ((length(x))/2);
corner4 = [C1 C2 C3 C4];

figure(1)
plot(x,y,'.-'); grid on; title('2D Domain'); xlabel('x'); ylabel('y');

%%
% Meshing
qnode =zeros(nelem,3);
for n = 1:nelem
     if n ~= nelem
            for i = 1:3
                qnode(n,1) = 2*(n-1)+1;
                qnode(n,2) = 2*(n-1)+2;
                qnode(n,3) = 2*(n-1)+3;
            end
     else
            for i = 1:3
                qnode(n,1) = 2*(n-1)+1;
                qnode(n,2) = 2*(n-1)+2;
                qnode(n,3) = 1;%2*(n-1)+3; % Closing final node for geometry
            end
     end
end

[p,dpdn,MatA,MatB,Ivec,Jmid]=HOBEM(x,y,nnode,nelem,qnode,ktvec,knvec,corner4,0);
toc


if m <2 
figure(2);
subplot(1,2,1)
Handle1 = plot(1:length(x),p,'.-'); hold on

title('Potential ')
xlabel('nodes')
% pa = .2.*y+5;
% pa = [5 5 5 4 3 3 3 4];
hold on
% plot(pa,'o r')

% legend('mf-HOBEM','Analytical')
grid on
subplot(1,2,2)
Handle2 = plot(1:length(x),dpdn,'.-'); hold on
grid on
title('Gradient ')
hold on

% plot(qa,'o r')
xlabel('nodes')
% legend('mf-HOBEM','Analytical')

else
    set(Handle1,'XData',1:length(x),'YData', p);
                set(Handle2,'XData',1:length(x),'YData', dpdn);
                drawnow;
end
end
%%
[U,V,NX,NY] = sderivatives(xf,yf,x,y,p,dpdn,qnode);
