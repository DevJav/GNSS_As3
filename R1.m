function r1 = R1(alpha)

r1 = [1 0 0;
      0 cos(2*pi*(alpha/360)) sin(2*pi*(alpha/360)); 
      0 -sin(2*pi*(alpha/360)) cos(2*pi*(alpha/360))];  

end