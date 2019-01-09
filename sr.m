function S = sr(f,a,b,n)
%this matlab function aproximates the integral of f from a to b by Simpsons rule. 

if mod(n,2)==1
  disp('n is odd')
  return
end

m=n/2; 
h=(b-a)/n; 
S=f(a)+f(b); 
x=linspace(a,b,n+1); 
S=S+2*sum(f(x(1:2:n-1)))+4*sum(f(x(2:2:n)));
%matlab starts couting at x1 not x0
S=S*h/3

end

