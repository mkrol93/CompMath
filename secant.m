function secant(f,x1,x2,tol,n)
% This MATLAB function solves f(x)=0 by the secant method.
u=f(x1);
v=f(x2);
err=abs(x2-x1);
iter=0;
disp('________________________________')
disp('iter     xn      f(xn)    err')
disp('________________________________')
fprintf('%2.0f %10.6f %10.6f %10.6f\n',iter,x1,u,err)
fprintf('%2.0f %10.6f %10.6f %10.6f\n',iter,x2,v,err)
while (err>tol)&&(iter<=n)
    if u==v
        disp('Division by 0')
        break
    else
    x=x2-v*(x2-x1)/(v-u);
    x1=x2;
    x2=x;
    u=v;
    v=f(x2);
    err=abs(x2-x1);
    iter=iter+1;
    fprintf('%2.0f %10.6f %10.6f %10.6f\n',iter,x2,v,err)
    end
end
if iter>n
    disp('Method failed to converge. Increase n.')
end
end

