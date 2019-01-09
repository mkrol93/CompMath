function shoot(f,a,b,alpha,beta,g1,g2,n,tol)
% This MATLAB function solves the boundary-value
% problem y''=f(x,y,y'), y(a)=alpha, y(b)=beta
% by using the shooting method with initial 
% guess-directions g1 and g2 and n subintervals.
h=(b-a)/n; % step size
x(1)=a;
x1(1)=alpha;
x2(1)=g1;
for i=1:n
    x(i+1)=x(i)+h;
    x1(i+1)=x1(i)+h*x2(i);
    x2(i+1)=x2(i)+h*f(x(i),x1(i),x2(i));
end
y1b=x1(n+1);
x2(1)=g2;
for i=1:n
    x(i+1)=x(i)+h;
    x1(i+1)=x1(i)+h*x2(i);
    x2(i+1)=x2(i)+h*f(x(i),x1(i),x2(i));
    fprintf('%6.0f, %6.2f, %12.2f \n',i,x(i+1),x1(i+1))
end
y2b=x1(n+1);
while abs(y2b-beta)>=tol
    k=g2;
    g2=g1+(g2-g1)*(beta-y1b)/(y2b-y1b);
    g1=k;
    x2(1)=g2;
    for i=1:n
    x(i+1)=x(i)+h;
    x1(i+1)=x1(i)+h*x2(i);
    x2(i+1)=x2(i)+h*f(x(i),x1(i),x2(i));
    end
    y1b=y2b;
    y2b=x1(n+1);
end
plot (x,x1)
end

