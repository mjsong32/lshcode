%% Make LSH
n = 1e5;
d = 2;
X = -1+2*rand(n,d);
k = 20;
L = 20;
T = lsh('e2lsh',L,k,d,X','range',10,'w',-100);
%% Make query with LSH
q = -1+2*rand(1,2);
tic;
[nnlsh,numcand]=lshlookup(q',X',T,'k',1,'sel','best');
r_lsh=sum((X(nnlsh,:)-q).^2);
fprintf('LSH l2 distance: %d\n',r_lsh);
toc;
%% Nearest neighbor 
tic;
[r_l1,min_idx] = min(sum(abs(X-q),2)); % note that this is l1 distance
r_l2=sum((X(min_idx,:)-q).^2);
fprintf('L1-NN l2 distance: %d\n',r_l2);
toc;