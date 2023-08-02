function [x,y,xfl, nodexbody,xfr,xd,xb,xp,yfl,nodeybody,yfr,yd,yb,yp,cx,cy,xbody,ybody,wpoints] = INITnwtbody(fx,h,L,w,k,Amp,wpoints)
        

        ; theta = 0; sway = 0; heave = 0;
        [xbody,ybody,wetx,wety,cx,cy] = bodyXY(L,wpoints,theta,sway,heave);     % Body discretization

        nodePortX = wetx(1);    nodePortY = wety(1);
        nodeStarX = wetx(end);  nodeStarY = wety(end);

        dy = h/8;
        xfl =  linspace(0,nodePortX,fx);        yfl = zeros(1,length(xfl));
        xfr = linspace(nodeStarX ,L,fx);        yfr = zeros(1,length(xfr));

        nodexbody = wetx(2:end-1)';  nodeybody = wety(2:end-1)'; 
        yd = linspace(-dy,-h+dy,17);         xd = L+zeros(1,length(yd)); 
        xb = linspace(L,0,121);              yb = -h+zeros(1,length(xb));              
        yp = linspace(-h+dy,-dy,11);         xp = zeros(1,length(yp));    
        
        x = [xfl nodexbody xfr xd xb xp];           y = [yfl nodeybody yfr yd yb yp];
       
end