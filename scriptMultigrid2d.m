%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Skript für Mehrgitterverfahren in 2D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function scriptMultigrid2d

% Quadrat
c4nQ = [0 0; 1 0; 0 1; 1 1];  
%c4nQ = [-1 -1; 1 -1; -1 1; 1 1];  
n4eQ = [2 3 1; 3 2 4];        
n4DbQ = [1 2; 2 4; 4 3; 3 1]; 

% L-förmig 
c4nL = [-1 -1; 0 -1; -1 0; 0 0; 1 0; -1 1; 0 1; 1 1];
n4eL = [2 3 1; 3 2 4; 4 6 3; 6 4 7; 5 7 4; 7 5 8];
n4DbL = [1 2; 2 4; 4 5; 5 8; 8 7; 7 6; 6 3; 3 1];



level = 4;
gebiet = "Lshape";
%gebiet = "Dreieck";
%gebiet = "Quadrat";

if gebiet == "Quadrat"
    [A, P, dof, b, c4n, n4e] = Generate(level, c4nQ, n4eQ, n4DbQ, @f, @boundary);
else
    [A, P, dof, b, c4n, n4e] = Generate(level, c4nL, n4eL, n4DbL, @f,@boundary);
end

um_0 = boundary(c4n);
%u = real(c4n); für Einheitsquadrat
u = um_0;
um_0(dof{end}) = zeros(length(dof{end}),1);
ug_0 = um_0;
uj_0 = um_0;
uc_0 = um_0;
tol = 1e-8;

[um_0,resm] = MGM2d(A, um_0, b, dof, P, tol, 100, 1, 2);
[ug_0(dof{end}), resg] = GaussSeidel(A{end}(dof{end},dof{end}), b(dof{end}), ug_0(dof{end}), 100);
[uj_0(dof{end}), resj] = Jacobi(A{end}(dof{end},dof{end}), b(dof{end}), uj_0(dof{end}), 100);
[~,~, ~, ~, resvec] = pcg(A{end}(dof{end}, dof{end}), b(dof{end}), tol, 100, [], []);
[uc_0(dof{end})] = CG(A{end}(dof{end}, dof{end}), b(dof{end}), uc_0(dof{end}),1e-8,400);
resc = resvec / norm(b(dof{end}));


figure(1);
X = c4n(:,1); Y = c4n(:,2); Z1 = um_0;
patch(X(n4e)',Y(n4e)',Z1(n4e)',Z1(n4e)');
%shading interp;

figure(2);
Z2 = u - ug_0;
patch(X(n4e)',Y(n4e)',Z2(n4e)',Z2(n4e)');
%shading interp;

figure(3);
Z3 = u - uj_0;
patch(X(n4e)',Y(n4e)',Z3(n4e)',Z3(n4e)');
%shading interp;

figure(4);
Z4 = u - um_0;
patch(X(n4e)',Y(n4e)',Z4(n4e)',Z4(n4e)');
%shading interp;

figure(7)
u_diff = u;
patch(X(n4e)',Y(n4e)',u_diff(n4e)',u_diff(n4e)')
%shading interp;


figure(5);
semilogy((0:length(resm)-1), resm, ".-", (0:length(resg)-1), resg, ".-", (0:length(resj)-1), resj, ".-",(0:length(resc) - 1), resc,".-",  LineWidth=2);
xlabel("Anzahl der Iterationen");
yline(tol,"r--",LineWidth=2);
legend("Multigrid", "Gauss-Seidel", "Jacobi", "CG-Verfahren","Toleranz");
ylabel("Relatives Residuum");

figure(6)
semilogy((0:length(resg)-1), resg, ".-", (0:length(resj)-1), resj, ".-");
legend("Gauss-Seidel", "Jacobi")
xlabel("Anzahl der Iterationen");
ylabel("Relatives Residuum");

disp(['Anzahl der Iterationen: ', num2str(size(resm,1) - 1)]);
disp(['No. dofs: ', num2str(size(dof{end},2))]);
end

function val = f(x)
%   val = -(2 -2 * pi^2) * exp(x(1)-x(2)) * sin(pi*x(1)) * sin(pi*x(2))...
%         -2 * pi * exp(x(1)-x(2)) * sin(pi * (x(2) - x(1)));
    val = 0;
end

% Randbedingung 
function u = boundary(x)
    [phi,r] = cart2pol(x(:,1),x(:,2));
    phi( phi<0 ) = phi( phi<0 ) + 2*pi;
    u = r.^(2/3).*sin(2/3*phi); 
    %u = zeros(size(x,1),1);
end


