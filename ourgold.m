function [xmin,ymin] = ourgold(f,a,b,tolx,toly)
% This MATLAB functions finds the minimum point
% of f in [a,b].
r=(3-sqrt(5))/2;
while abs(b-a)>tolx && abs(f(b)-f(a))>toly
    x1=a+r*(b-a);
    x2=a+(1-r)*(b-a);
    if f(x1)<f(x2)
        b=x2;
    else
        a=x1;
    end
end
xmin=(a+b)/2;
ymin=f(xmin);
end

