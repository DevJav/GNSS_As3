function r3 = R3(alpha)

r3 = [cos(2*pi*(alpha/360)) sin(2*pi*(alpha/360)) 0;
      -sin(2*pi*(alpha/360)) cos(2*pi*(alpha/360)) 0; 
      0 0 1];  

end