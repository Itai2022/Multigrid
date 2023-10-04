function [x1, resmax] = GaussSeidel(A,b,x0,Itmax)
    D = diag(diag(A));
    L = -tril(A, -1);
    U = -triu(A, 1);
    G = (D - L)\U;
    f = (D-L)\b;
    resmax = zeros(Itmax, 1);
    bnorm = norm(b); 
    resmax(1) = norm(b - A*x0)/bnorm;
    for i = 1:Itmax
        x1 = G*x0 + f;
%         if i == 1
%             r1 = x1;
%         elseif i == 10
%             r10 = x1;
%         elseif i == 35
%             r35 = x1;
%         end
        resmax(i+1) = norm(b - A*x1)/bnorm;
        if resmax(i+1) < 1e-6
            resmax = resmax(resmax > 0);
            break;
        end
        x0 = x1;
    end
end