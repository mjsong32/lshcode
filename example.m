%% Make LSH
n = 1e5;
d = 10;
X = -1+2*rand(n,d);
k = 10;
L = 10;
T = lsh('e2lsh',L,k,d,X','range',10,'w',-50);
%% Make queries
n_test = 1e2;
q = -1+2*rand(n_test,d);
%% LSH lookup
tic;
for i=1:n_test
    [nnlsh,numcand]=lshlookup(q(i,:)',X',T,'k',1,'sel','best');
end
toc;
rnd_idx=randi(n_test);
r_lsh=sum((X(nnlsh,:)-q(rnd_idx,:)).^2);
fprintf('LSH l2 distance: %d\n',r_lsh);
%% Faster LSH lookup (return only the top-1)
tic;
for i=1:n_test
    [nnlsh,numcand]=fastlshlookup(q(i,:)',X',T,'k',1,'sel','best');
end
toc;
r_lsh=sum((X(nnlsh,:)-q(rnd_idx,:)).^2);
fprintf('fast LSH l2 distance: %d\n',r_lsh);
%% Nearest neighbor 
tic;
for i=1:n_test
    [r_l1,min_idx] = min(sum(abs(X-q(i,:)),2)); % note that this is l1 distance
end
toc;
r_l2=sum((X(min_idx,:)-q(rnd_idx,:)).^2);
fprintf('L1-NN l2 distance: %d\n',r_l2);