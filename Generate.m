function [A4lvl, P4lvl, dof4lvl, b, c4n, n4e] = Generate(level, c4n, n4e, n4Db,f,bound)

[n4ed,ed4e] = SIDES(n4e,n4Db); 
A4lvl = cell(level,1);
P4lvl = cell(level,1);
dof4lvl = cell(level,1);

% generieren die Systematrix für jede Gitterebene
for k = 1 : level
    % uniforme Verfeinerung.
    [c4n,n4e,n4Db,P] = BISECTP(c4n,n4e,n4Db,n4ed,ed4e,true(size(n4e,1),1));
    [n4ed,ed4e] = SIDES(n4e,n4Db);
   
    [A,dof,b] = CFEM(c4n,n4e,n4Db,f,bound);
    A4lvl{k} = A;
    dof4lvl{k} = dof;
    P4lvl{k} = P;
end

end