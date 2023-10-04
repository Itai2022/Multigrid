function [n4s,s4e,e4s,s4Db] = SIDES(n4e,n4Db)
d = size(n4e,2)-1;
sides = n4e(:,fliplr(nchoosek(1:d+1,d)'));
sides = sort([reshape(sides',d,[]),n4Db'],1)';
[n4s, e4s1, s4e] = unique(sides,'rows','stable');
[~, perm] = unique(sides,'rows','first');
[~, perm] = sort(perm);
[~, e4s2] = unique(sides,'rows','last');
e4s2 = e4s2(perm);
e4s2(e4s2>(d+1)*size(n4e,1)) = 0;
e4s = ceil([e4s1,e4s2]/(d+1));
s4Db = s4e((d+1)*size(n4e,1)+(1:size(n4Db,1)));
s4e = reshape(s4e(1:(d+1)*size(n4e,1)),d+1,size(n4e,1))';
end