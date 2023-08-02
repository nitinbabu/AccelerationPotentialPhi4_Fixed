function [xcr,ycr,pcr] = CSregrid(xo,yo,poto)


nodes = length(xo); nnode = zeros(1,nodes); chordL = zeros(1,length(nnode));
segment = length(xo)-1;
for k = 2:length(nnode)
    chordL(k) = sqrt((xo(k)-xo(k-1))^2+(yo(k)-yo(k-1))^2); % This is the chord length parameter
end

t = chordL; % t2 = tk+1 or t_current and tk = 0;
[Bx,By,Bp] = CsplineInterpXYP(nodes,t,xo,yo,poto);  %% Gets Spline functions based on chord length t


Px = 11;
%% Creates the Spline points based on polygonal arclength or chordlength
for m = 1:segment
    if m == 1 
        t_end = t(m+1);
        tvec = 0:t_end/Px:t_end;
        for i =  1:length(tvec)
            tv = tvec(i);
            x_t(i,m) = [1 tv tv^2 tv^3]*Bx(:,m);
            y_t(i,m) = [1 tv tv^2 tv^3]*By(:,m);
            p_t(i,m) = [1 tv tv^2 tv^3]*Bp(:,m);
            
            xp_t(i,m) = [1 2*tv 3*tv^2]*Bx(2:4,m);
            yp_t(i,m) = [1 2*tv 3*tv^2]*By(2:4,m);
            
            
            
        end
    
    elseif m == segment
        t_end = t(m+1);
        tvec = 0:t_end/Px:t_end;
        for i =  1:length(tvec)
            tv = tvec(i);
            x_t(i,m) = [1 tv tv^2 tv^3]*Bx(:,m);
            y_t(i,m) = [1 tv tv^2 tv^3]*By(:,m);
            p_t(i,m) = [1 tv tv^2 tv^3]*Bp(:,m);

            xp_t(i,m) = [1 2*tv 3*tv^2]*Bx(2:4,m);
            yp_t(i,m) = [1 2*tv 3*tv^2]*By(2:4,m);            

        end
    
    else
        t_end = t(m+1);
        tvec = 0:t_end/Px:t_end;
        for i =  1:length(tvec)
            tv = tvec(i);
            x_t(i,m) = [1 tv tv^2 tv^3]*Bx(:,m);
            y_t(i,m) = [1 tv tv^2 tv^3]*By(:,m);
            p_t(i,m) = [1 tv tv^2 tv^3]*Bp(:,m);
            
            xp_t(i,m) = [1 2*tv 3*tv^2]*Bx(2:4,m);
            yp_t(i,m) = [1 2*tv 3*tv^2]*By(2:4,m);
            
        end
    
    end
    
    
    
end

%% Arc Length calculation from derivative
[wg, xig,~,~] = GuassianQuad8points();
%% t = (tbar).(tend/2) + tend/2
for m = 1:segment
    if m == 1 
        t_end = t(m+1);
        tvec = xig;
        arcm = 0;
        for i =  1:length(tvec)
            tvbar = tvec(i)*t_end/2 +t_end/2;  
            
            %% xbar = x - (a+b)/2  / (b-a)/2       here a = 0 b = chord length t(i) each segment
            %% -1 to 1
            xp_tbar(i,m) = [1 2*tvbar 3*tvbar^2]*Bx(2:4,m);
            yp_tbar(i,m) = [1 2*tvbar 3*tvbar^2]*By(2:4,m);   
            arcf(i,m) = (t_end/2)*sqrt(xp_tbar(i,m)^2 + yp_tbar(i,m)^2);
            arcm = arcm+wg(i)*arcf(i,m);
        end
        s(m) = arcm; 
    
    elseif m == segment
        t_end = t(m+1);
        tvec = xig;
        arcm = 0;
        for i =  1:length(tvec)
            tvbar = tvec(i)*t_end/2 +t_end/2;   
            xp_tbar(i,m) = [1 2*tvbar 3*tvbar^2]*Bx(2:4,m);
            yp_tbar(i,m) = [1 2*tvbar 3*tvbar^2]*By(2:4,m);  
            arcf(i,m) = (t_end/2)*sqrt(xp_tbar(i,m)^2 + yp_tbar(i,m)^2);
            arcm = arcm+wg(i)*arcf(i,m);
        end
        s(m) = arcm; 
    else
        t_end = t(m+1);
        tvec = xig;
        arcm = 0;
        for i =  1:length(tvec)
            tvbar = tvec(i)*t_end/2 +t_end/2;                           
            xp_tbar(i,m) = [1 2*tvbar 3*tvbar^2]*Bx(2:4,m);
            yp_tbar(i,m) = [1 2*tvbar 3*tvbar^2]*By(2:4,m);
            arcf(i,m) = (t_end/2)* sqrt(xp_tbar(i,m)^2 + yp_tbar(i,m)^2);
            arcm = arcm+wg(i)*arcf(i,m);
        end  
        s(m) = arcm; 
    end   
end
s = [0 s];
arclength = sum(s);
polygonal_arclength  = sum(t);

%% ^ Uptill here all correct


%% refils all points in a single curve % only needed to plot
k =1;i = 1; m = 1;
while k < segment*Px+1
  
   
    if i == Px+1
        i = 1;
        if m<segment
        m = m+1; 
        end
    else
        
    end
    xs(k) = x_t(i,m); 
    ys(k) = y_t(i,m); 
    ps(k) = p_t(i,m);
    
    xsp(k) = xp_t(i,m);
    ysp(k) = yp_t(i,m);
    B(k)= m; % not needed
    i = i+1;
    k = k+1;
    
end

 xs(k) = x_t(end,m); 
 ys(k) = y_t(end,m); 
 ps(k) = p_t(end,m);
 
 xsp(k) = xp_t(end,m);
 ysp(k) = yp_t(end,m);
 B(k) = m;
 

 
 
 %% 2nd Spline
 [Cx,Cy,Cp] = CsplineInterpXYP(nodes,s,xo,yo,poto); % gets us spline functions based on arc length 's'
 Px = 20;
%% Creates the Spline based on arc length parameter 's' points
for m = 1:segment
    if m == 1 
        s_end = s(m+1);
        svec = 0:s_end/Px:s_end;
        for i =  1:length(svec)
            sv = svec(i);
            x_c(i,m) = [1 sv sv^2 sv^3]*Cx(:,m);
            y_c(i,m) = [1 sv sv^2 sv^3]*Cy(:,m);
            p_c(i,m) = [1 sv sv^2 sv^3]*Cp(:,m);            
            
        end
    
    elseif m == segment
        s_end = s(m+1);
        svec = 0:s_end/Px:s_end;
        for i =  1:length(svec)
            sv = svec(i);
            x_c(i,m) = [1 sv sv^2 sv^3]*Cx(:,m);
            y_c(i,m) = [1 sv sv^2 sv^3]*Cy(:,m);
            p_c(i,m) = [1 sv sv^2 sv^3]*Cp(:,m);

        end
    
    else
        s_end = s(m+1);
        svec = 0:s_end/Px:s_end;
        for i =  1:length(svec)
            sv = svec(i);
            x_c(i,m) = [1 sv sv^2 sv^3]*Cx(:,m);
            y_c(i,m) = [1 sv sv^2 sv^3]*Cy(:,m);
            p_c(i,m) = [1 sv sv^2 sv^3]*Cp(:,m);
                    
            
        end
    
    end
    
    
    
end
 



%% Gets all segmented values in a single vector to plot
k =1;i = 1; m = 1;
while k < segment*Px+1
  
   
    if i == Px+1
        i = 1;
        if m<segment
        m = m+1; 
        end
    else
        
    end
    xsc(k) = x_c(i,m); 
    ysc(k) = y_c(i,m);
    psc(k) = p_c(i,m);


    i = i+1;
    k = k+1;
    
end

 xsc(k) = x_c(end,m); 
 ysc(k) = y_c(end,m); 
 psc(k) = p_c(end,m);



 
 
 %% Regrid Part
%  Divide arc length into fx equal segments to obtain segment arc length
 eqarcs = arclength/segment;
 sq = [0 eqarcs+zeros(1,segment)];


%% Get Spline based on arc length parameter 's' points overall
% Clue: We are only interested in new grids

mainarc =  cumsum(s);
desgarc = cumsum(sq);
;
for i = 2:segment

    if  (desgarc(i)>mainarc(i))
        if i ==2
        sqnew(i) = desgarc(i)-s(i);
        locs(i) = i+1;
        else
        sqnew(i) = s(i+1)-(mainarc(i+1)-desgarc(i));
        locs(i) = i+1;
        end
    elseif desgarc(i)<mainarc(i)
        sqnew(i) = s(i) - (mainarc(i)-desgarc(i));
        locs(i) = i;
    elseif desgarc(i)==mainarc(i)
        sqnew(i) = s(i);
        locs(i)  = i;
    end 


end
s;
sqnew;
locs;



xcr(1) = xo(1);
ycr(1) = yo(1);
pcr(1) = poto(1);
% Get equal arc length in terms of main arc segment
for i = 2:segment
    sv =sqnew(i) ;
    m = locs(i)-1;
    x_cr = [1 sv sv^2 sv^3]*Cx(:,m);
    y_cr = [1 sv sv^2 sv^3]*Cy(:,m);
    p_cr = [1 sv sv^2 sv^3]*Cp(:,m);
    
    

    xcr(i) = x_cr;
    ycr(i) = y_cr;
    pcr(i) = p_cr;
    
    

end

xcr(nodes) = xo(end);
ycr(nodes) = yo(end);
pcr(nodes) = poto(end);



end