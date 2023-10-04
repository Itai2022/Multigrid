%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Modifiziert aus Code von Prof. Dr. Joscha Gedicke
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [A,dof,b] = CFEM(c4n,n4e,n4Db,f, boundary)
n = size(c4n,1); m = size(n4e,1); d = size(c4n,2);
A = zeros(d+1,d+1,m); b = zeros(m,d+1); %x = zeros(n,1);
for j = 1 : m
    area = abs(det([ones(1,d+1);c4n(n4e(j,:),:)'])/prod(2:d));
    grads = [ones(1,d+1);c4n(n4e(j,:),:)']\[zeros(1,d);eye(d)];
    A(:,:,j) = area*(grads*grads');
    mid = 1/3*sum(c4n(n4e(j,:),:),1);
    b(j,:) = f(mid).*area/(d+1);
end
J = n4e(:,cumsum(ones(d+1,d+1)))'; K = n4e(:,cumsum(ones(d+1,d+1))')';
A = sparse(J(:),K(:),A(:),n,n);
b = accumarray(n4e(:),b(:),[n 1]);
dof = setdiff(1:n,n4Db(:));
boundary_index = 1:n;
boundary_index(dof) = [];
b = b - A(:, boundary_index) * boundary(c4n(boundary_index,:));
end 