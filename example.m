%% Make LSH
n = 5e5;
d = 10;
X = -1+2*rand(n,d);
k = 10;
L = 10;
T = lsh('e2lsh',L,k,d,X','range',10,'w',-50);
%% Make queries
n_test = 100;
q = -1+2*rand(n_test,d);

timeStamp=zeros(5);
%% LSH lookup
tStart = tic;
dist=0;
for i=1:n_test
    [nnlsh,numcand]=lshlookup(q(i,:)',X',T,'k',1,'sel','best');
    dist = dist + sum((X(nnlsh,:)-q(i,:)).^2);
end
timeStamp(1)=toc(tStart);
fprintf('average LSH time: %.4f\n',timeStamp(1)/n_test);

rnd_idx=randi(n_test);
r_lsh=sum((X(nnlsh,:)-q(rnd_idx,:)).^2);
fprintf('LSH l2 distance: %d\n',r_lsh);
fprintf('average LSH l2 distance: %f\n',dist/n_test);

%% Faster LSH lookup (return only the top-1)
timeStamp(2)=toc(tStart);
dist=0;
for i=1:n_test
    [nnlsh,numcand]=fastlshlookup(q(i,:)',X',T,'k',1,'sel','best');
    dist = dist + sum((X(nnlsh,:)-q(i,:)).^2);
end
timeStamp(3)=toc(tStart);
fprintf('average fastLSH time: %.4f\n',(timeStamp(3)-timeStamp(2))/n_test);

r_lsh=sum((X(nnlsh,:)-q(rnd_idx,:)).^2);
fprintf('fast LSH l2 distance: %d\n',r_lsh);
fprintf('average fast LSH l2 distance: %f\n',dist/n_test);

%% Nearest neighbor
timeStamp(4)=toc(tStart);
dist=0;
for i=1:n_test
    [r_l1,min_idx] = min(sum(abs(X-q(i,:)),2)); % note that this is l1 distance
    dist = dist + sum((X(min_idx,:)-q(i,:)).^2);
end

timeStamp(5)=toc(tStart);
fprintf('average l1 exact time: %.4f\n',(timeStamp(5)-timeStamp(4))/n_test);

r_l2=sum((X(min_idx,:)-q(rnd_idx,:)).^2);
fprintf('L1-NN l2 distance: %d\n',r_l2);
fprintf('average L1-NN l2 distance: %f\n',dist/n_test);