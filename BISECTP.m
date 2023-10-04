function [c4n,n4e,n4Db,P] = BISECTP(c4n,n4e,n4Db,n4ed,ed4e,r4e)
n = size(c4n,1);
ref4ed = false(size(n4ed,1),1);
ref4ed(ed4e(r4e,:)) = true;
while any(ref4ed(ed4e(:,3))==0 & any(ref4ed(ed4e(:,[1,2])),2))
    ref4ed(ed4e(ref4ed(ed4e(:,3))==0 & any(ref4ed(ed4e(:,[1,2])),2),3))= true;
end
m4ed = zeros(size(n4ed,1),1);
m4ed(ref4ed==1) = n + (1:nnz(ref4ed));
c4n = [c4n; (c4n(n4ed(ref4ed,1),:) + c4n(n4ed(ref4ed,2),:))/2];
m4n = sparse(n4ed(:,1),n4ed(:,2),m4ed,size(c4n,1),size(c4n,1));
m4n = m4n + m4n';
while any(r4e)
    m = full(m4n(uint64(n4e(r4e,1)) + size(m4n,1)*uint64(n4e(r4e,2)-1)));
    n4e = [n4e(~r4e,:); n4e(r4e,3),n4e(r4e,1),m; n4e(r4e,2),n4e(r4e,3),m];
    r4e = full(m4n(uint64(n4e(:,1)) + size(m4n,1)*uint64(n4e(:,2)-1))) > 0;
end
r4Db = full(m4n(uint64(n4Db(:,1)) + size(m4n,1)*uint64(n4Db(:,2)-1))) > 0;
m = full(m4n(uint64(n4Db(r4Db,1)) + size(m4n,1)*uint64(n4Db(r4Db,2)-1)));
n4Db = [n4Db(~r4Db,:);n4Db(r4Db,1),m; m,n4Db(r4Db,2)];
I = [(1:n)',(1:n)';m4ed,m4ed];
J = [(1:n)',(1:n)';n4ed(ref4ed==1,:)];
P = sparse(I(:),J(:),.5*ones(2*size(I,1),1));
end