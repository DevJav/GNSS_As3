function r2 = R2(alpha)

r2 = [cos(2*pi*(alpha/360)) 0 -sin(2*pi*(alpha/360));
      0 1 0; 
      sin(2*pi*(alpha/360)) 0 cos(2*pi*(alpha/360))];  

end