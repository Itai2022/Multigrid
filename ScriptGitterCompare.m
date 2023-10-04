function ScriptGitterCompare

% Einheitsquadrat
c4nQ = [0 0; 1 0; 0 1; 1 1];  
n4eQ = [2 3 1; 3 2 4];        
n4DbQ = [1 2; 2 4; 4 3; 3 1]; 

% L-förmig 
c4nL = [-1 -1; 0 -1; -1 0; 0 0; 1 0; -1 1; 0 1; 1 1];
n4eL = [2 3 1; 3 2 4; 4 6 3; 6 4 7; 5 7 4; 7 5 8];
n4DbL = [1 2; 2 4; 4 5; 5 8; 8 7; 7 6; 6 3; 3 1];

for level = 7:10
    gebiet = "Lshape";
    %gebiet = "Dreieck";
    
    if gebiet == "Quadrat"
        [A, P, dof, b, ~, ~] = Generate(level, c4nQ, n4eQ, n4DbQ, @f, @boundary);
    else
        [A, P, dof, b, ~, ~] = Generate(level, c4nL, n4eL, n4DbL, @f, @boundary);
    end
    
    um_0 = boundary(c4nL);
    um_0(dof{end}) = zeros(length(dof{end}),1);
    umw_0 = um_0;
    ug_0 = um_0;
    uj_0 = um_0;
    uc_0 = um_0;
    tol = 1e-8;
    
    [~,resm] = MGM2d(A, um_0, b, dof, P, tol, 100, 1, 2);
    [~,resmw] = MGM2d(A, umw_0, b, dof, P, tol, 100, 2, 2);
    [~, resg] = GaussSeidel(A{end}(dof{end},dof{end}), b(dof{end}), ug_0(dof{end}), 1000);
    [~, resj] = Jacobi(A{end}(dof{end},dof{end}), b(dof{end}), uj_0(dof{end}), 1000);
    [~,~, ~, ~, resvec] = pcg(A{end}(dof{end}, dof{end}), b(dof{end}), tol, 1000, [], [],uc_0(dof{end}));
    resc = resvec / norm(b(dof{end}));
    
    disp(["Level = ", num2str(level)]);
    disp(['Anzahl der Iterationen von Multigv: ', num2str(size(resm,1) - 1)]);
    disp(['normv:',num2str(resm(end))])
    disp(['Anzahl der Iterationen von Multigw: ', num2str(size(resmw,1) - 1)]);
    disp(['normw:',num2str(resmw(end))])
    disp(['Anzahl der Iterationen von GS: ', num2str(size(resg,1) - 1)]);
    disp(['normgs:',num2str(resg(end))])
    disp(['Anzahl der Iterationen von Ja: ', num2str(size(resj,1) - 1)]);
    disp(['normja:',num2str(resj(end))])
    disp(['Anzahl der Iterationen von cg: ', num2str(size(resc,1) - 1)]);
    disp(['normcg:',num2str(resc(end))])
    disp(['No. dofs: ', num2str(size(dof{end},2))]);
end
end

function val = f(x)
    val = 0;
    % Für Einheisquadrat
   %val = -pi^2 * cos(pi*x(1)) - 4*pi^2 * sin(2 * pi * x(2));
end

% Randbedingung 
function u = boundary(x)
    [phi,r] = cart2pol(x(:,1),x(:,2));
    phi( phi<0 ) = phi( phi<0 ) + 2*pi;
    u = r.^(2/3).*sin(2/3*phi); 
    % Für Einheisquadrat
    %u = zeros(size(x,1),1); 
end



