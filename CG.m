function [x] = CG(A,b,x0,eps,itmax)
r0 = b - A*x0;
p0 = r0;
for i = 1:itmax
    if abs(p0) < eps
        break;
    end
    a0 = r0'*r0/(p0'*A*p0);
    x1 = x0 + a0*p0;
    r1 = r0 -a0*A*p0;
    b0 = r1'*r1/(r0'*r0);
    p1 = r1 + b0*p0;
    x0 = x1;
    r0 = r1;
    p0 = p1;
end
x = x0;
end