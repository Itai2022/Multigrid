%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mehrgitterverfahren in 2 Dimension
%% A die Systemmatrix, b die rechte Seite
%% u_0 Anfangswert tol die Abbruchskriterium 
%% maxIter die maximale Anzahl von Iterationen
%% V_Cycle: mu = 1; W_Cycle: mu = 2
%% A ist eine Menge von Systemmatrix aller Knoten
%% min ist die minimale Gitterebene
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [u_0, resmax] = MGM2d(A, u_0, b, dof, P, tol, maxIte,mu,min)
    resmax = zeros(maxIte,1);
    resmax(1) = norm(b(dof{end}) - A{end}(dof{end},dof{end})*u_0(dof{end}))/norm(b(dof{end})); 
    %resmax(1) = max(abs(b(dof{end}) - A{end}(dof{end},dof{end})*u_0(dof{end})));
    l = size(A,1);
    for ite = 1:maxIte
        u_0 = W_V_Cycle(l, A, dof, P,u_0, b,mu,min);
        res = b(dof{end}) - A{end}(dof{end},dof{end})*u_0(dof{end});
        %resmax(ite+1) = max(abs(res)); 
        resmax(ite + 1) = norm(res)/norm(b(dof{end})); 
        if resmax(ite + 1) < tol
            resmax = resmax(resmax>0);
            break;
        end
    end
    if ite == maxIte
        fprintf("Der Algorithmus konvergiert vielleicht nicht.\n");
    end
end

function u0 = W_V_Cycle(l, A, dof, P, u0, b,mu,min)
    % 2-mal Vorglättung mit GaussSeidel
    u0(dof{l}) = gaussSeidel(A{l}(dof{l},dof{l}),b(dof{l}),u0(dof{l}),2);

    %Grob-Gitter_korrektur

    % Berechnung der Defekt.
    d1 = zeros(size(A{l},1),1);
    d1(dof{l}) = b(dof{l}) - A{l}(dof{l},dof{l})...
                 *u0(dof{l});           

    % Restriktion des Defekts.
    d0 = P{l}'*d1;              
    e0 = zeros(length(d0),1);
    if l == min
        e0(dof{min - 1}) =   A{min - 1}(dof{min - 1},dof{min - 1}) \...
                             d0(dof{min - 1});
    else
        for i = 1:mu
            e0 = W_V_Cycle(l-1, A, dof, P, e0, d0,mu,min);
        end
    end
    v0 = P{l}*e0;
    u0(dof{l}) = u0(dof{l}) + v0(dof{l});  % Prolongation der Fehler.

    %A-Posteriori-Glättung mit GaussSeidel
    u0(dof{l}) = gaussSeidel(A{l}(dof{l},dof{l}),b(dof{l}),u0(dof{l}),2);
end

function x1 = gaussSeidel(A,b,x0,Itmax)
    D = diag(diag(A));
    L = -tril(A, -1);
    U = -triu(A, 1);
    G = (D - L)\U;
    f = (D-L)\b;
    for i = 1:Itmax
        x1 = G*x0 + f;
        x0 = x1;
    end
end
