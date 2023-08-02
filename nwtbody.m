function [x,y,xfl, xbody,xbodyL,xbodyB,xbodyR, xfr,xd,xb,xp,yfl,ybody,ybodyL,ybodyB,ybodyR,yfr,yd,yb,yp] = nwtbody(fx,h,L,w,k,Amp)
        
        R = 0.064;
        dx = L/6;
        dy = h/8;
%         draugt = 0.25

        BL = L/2-0.25; BR = L/2+0.25; % Keeping D/L ~> 0.2
        
        xfl =  linspace(0,BL,fx);
%         xfl = exp(xfl)
        yfl = zeros(1,length(xfl));

        
        ybodyL =linspace(-0.025,-0.15, 9);    xbodyL =  BL+zeros(1,length(ybodyL));
        

        
        xbodyB = linspace(BL,BR, 9 );        ybodyB = -0.25+zeros(1,length(xbodyB));
        ybodyR = linspace(-0.15,-0.025,9); xbodyR =  BR+zeros(1,length(ybodyR));
        
        
        xbody = [xbodyL xbodyB xbodyR]; ybody = [ybodyL ybodyB ybodyR];
        
  
     
        xfr = linspace(BR,L,fx);
        yfr = zeros(1,length(xfr));
        
        yd = linspace(-dy,-h+dy,11);         xd = L+zeros(1,length(yd)); 
        xb = linspace(L,0,121);        yb = -h+zeros(1,length(xb));              
        yp = linspace(-h+dy,-dy,11);      xp = zeros(1,length(yp));    
        
        x = [xfl xbody xfr xd xb xp];           y = [yfl ybody yfr yd yb yp];
        
end