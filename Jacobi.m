% Das Jacobi Verfahren 
% A ist die Systemmatrix
% b ist die rechte Seite
% x_0 ist ein Anfangswert
% Itmax ist die maximale Anzahl der Iterationen
function [x1,resmax] = Jacobi(A,b,x0,Itmax)
    D = diag(diag(A));
    L = -tril(A, -1);
    U = -triu(A, 1);
    G = D\(L+U);
    f = D\b;
    resmax = zeros(Itmax,1);
    bnorm = norm(b);
    resmax(1) = norm(b - A*x0)/bnorm;
    for i = 1:Itmax
        x1 = G*x0 + f;
%         if i == 1
%             r1 = b - A*x1;
%             %r1 = x1;
%         elseif i ==10
%             r10 = b - A*x1;
%             %r10 = x1;
%         elseif i == 35
%             %r35 = x1;
%             r1 = b - A*x1;
%         end
        resmax(i+1) = norm(b - A*x1)/bnorm;
        if resmax(i+1) < 1e-6
            resmax = resmax(resmax > 0);
            break;
        end
        x0 = x1;
    end
end