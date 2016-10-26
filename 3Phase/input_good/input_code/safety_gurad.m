function [valiation,x,length ] = safety_gurad(x,x_o,max,min,length,theta )
valiation=0;
 if x<min
     length=(x_o-min)/theta;
     valiation=1;
     x=min;
 end
 if x>max
     length=(max-x_o)/theta;
     valiation=1;
     x=max;
 end
 
 
end

