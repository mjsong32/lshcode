function [iNN,cand] = fastlshlookup(x0,x,T,varargin)
% [iNN,cand] = fastlshlookup(x0,x,T)
% A simple faster version of lsh lookup applicable if only interested
% in the closest (top-1) point. Instead of sorting the distances of
% candidates in O(logn), simply return the min in O(n).

% measureTime=zeros(1,4);
% tStart=tic;

distfun='lpnorm';
switch T(1).type
 case 'lsh', distargs={1};
 case 'e2lsh', distargs={2};
 case 'hyperplane', distargs={2};
end
k=1;
r=inf;
r=5000;
sel='best';
f=[];
fargs=[];
verb=T(1).verbose;

% parse args.
for a=1:2:length(varargin)
  eval(sprintf('%s = varargin{a+1};',varargin{a}));
end

l = length(T);
iNN=[];

% measure time
% measureTime(1)=toc(tStart);

% find the union of buckets in all tables that match query
for j=1:l
  % look up T_j
  % buck is the # of bucket in T{j}
  buck = findbucket(T(j).type,x0,T(j).I);
  % find the bucket in j-th table
  key = lshhash(buck);
  ihash = T(j).bhash{key}; % possible matching buckets
  
  if (~isempty(ihash)) % nothing matches -> MJ: multiple matches?
    b = ihash(find(all(bsxfun(@eq,buck,T(j).buckets(ihash,:)),2)));
    if (~isempty(b))
      iNN = [iNN T(j).Index{b}];
    end
  end
end

% delete duplicates
[iNN,iu]=unique(iNN);
cand = length(iNN);

% measure time
% measureTime(2)=toc(tStart);

% now iNN has the collection of candidate indices 
% we can start examining them

if (verb > 0)
  fprintf('Examining %d candidates\n',cand);
end

if (~isempty(iNN))
  
  if (strcmp(sel,'best'))
    D=feval(distfun,x0,Xsel(x,iNN),distargs{:});
    [~,min_idx]=min(D);
    iNN=iNN(min_idx);
    % measure time
    % measureTime(3)=toc(tStart);
    
  else % random
    
    rp=randperm(cand);
    choose=[];
    for i=1:length(rp)
      d = feval(distfun,x0,Xsel(x,iNN(rp(i))),distargs{:});
      if (d <= r)
	choose = [choose iNN(rp(i))];
	if (length(choose) == k)
	  break;
	end
      end
    end
    iNN = choose;
  end
end

% measure time
% measureTime(4)=toc(tStart);
% disp('Measure time');
% disp(measureTime);

%%%%%%%%%%%%%%%%%%%%%%%%55 
function x=Xsel(X,ind)
% x=Xsel(X,ind)
% selects (i.e. collects) columns of cell array X
% (automatically determining the class, and looking for each column in
% the right cell.)

if (~iscell(X))
  x=X(:,ind);
  return;
end

d=size(X{1},1);

if (strcmp(class(X{1}),'logical'))
  x=false(d,length(ind));
else
  x=zeros(d,length(ind),class(X{1}));
end
sz=0; % offset of the i-th cell in X
collected=0; % offset within x
for i=1:length(X)
  thisCell=find(ind > sz & ind <= sz+size(X{i},2));
  if (~isempty(thisCell))
    x(:,thisCell)=X{i}(:,ind(thisCell)-sz);
  end
  collected=collected+length(thisCell);
  sz=sz+size(X{i},2);      
end
