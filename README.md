This directory contains an simple implementation of e2lsh (exact euclidean lsh). It is a minor modification of Greg Shakhnarovich's LSH codes.

## Main files

1. `lsh.m` - creates the LSH data structure for a given data set.
2. `lshlookup.m` - the approximate similarity search using the LSH. 
3. `lshstats.m` - examine statistics of LSH data structure

### Helper files (not-so-important)

1. `lshprep.m`, `lshfunc.m` - these functions, respectively, set up the LSH hashing functions and the hash tables based on these functions. You may not need to call these directly (they are called inside lsh.m)
2. `lshins.m` - insert data into the LSH data structure. Again, usually you may not need this function since it, too, is called inside lsh.m, but you will need it if you need to insert additional data points after the hash tables are created.

## LSH parameters

Should be determined experimentally (although recommendations for worst-case data are known)

1. r: closeness threshold.
2. c: multiplicative gap for (r, cr)-ANN. note that c > 1.
3. w: discretization chunk size. Recall that the e2lsh hash function is defined as h_{a,b}(p) = \lfloor (a\cdot p + b)/w \rfloor​, where a is a Gaussian random vector and b is a uniformly random number taking values in [0,w]​.
4. k: number of bits for each hash function (table). larger k amplifies the collision probability **gap** between r-close points and cr-far points. However, large k requires large L since it brings down the collision probability (for both close and far points), a.k.a many empty buckets.
5. L: number of hash tables.

### Theoretical recommendation

Assume r=1​ without loss of generality (can scale the entire dataset by r). For any point p,q \in R^d s.t. dist(p,q) \le 1, there exists at least one collision between p​ and q among L hash tables {g_j : R^d -> Z^k\} with probability 1-\delta if

L > \frac{\ln(1/\delta)}{-\ln(1-P_1^k)}.

**Caveat**: P_1 = P_1(w,c)​ is a bit complicated to determine. It's simpler to just run grid search over w, c, k, L.

## Examples

Let X \in R^{d \times N}​ be the dataset. To construct L=30 hash tables, each with k=25 bits,

```matlab
load data;
X = data;
k=25;
L=30;
T=lsh('e2lsh',L,k,size(X,1),X,'range',255,'w',-4);
```

To perform a lookup on query q on LSH data structure T, which consists of L hash tables,

```matlab
[nnlsh,numcand]=lshlookup(q,X,T,'k',10,'sel','best');
```

`nnlsh` is an array of k​-nearest neighbors (k=10) sorted by distance from query q.

To check the health of the created LSH hash table T, use `lshstats.m`

```matlab
m=2;
Xtest=X(:,1:1000);
lshstats(T,'test',X,Xtest,m);
```

Returns median, max, and expected occupancy of the buckets in each hash table. It also returns the mean and max number of candidates returned for query points in `Xtest`. Lastly, it reports the number of failures for retrieving the m-nearest neighbor for each query point in `Xtest`.

## Function explanations

1. `T=lsh(TYPE, L, k, d, X, varargin)`: only TYPE available is 'e2lsh'. May extend to multiprobe LSH, if necessary.

   **input**

   - `L`: number of hash tables, `k`: number of coordinates (bits) for each hash table.
   - `d`: dimension of each data point, `X`: dataset.

   **output**

   - `T`: LSH data structure.

   **varargs**

   - `range`: specify the range of values in each data dimension. Can either be scalar, d-dimensional vector, or 2xd-dimensional matrix. 2xd-dimensional matrix is interpreted as [r_min, r_max] for each coordinate in [d].
   - `w`: specify chunk size. if given as a negative value, the absolute value of w is interpreted as the desired number of intervals for each random projection.

2. `[nnlsh, num_cand]=lshlookup(query, X, T, varargin)`

   **input**

   - `query`: query point, `X`: dataset, `T`: LSH data structure.

   **output**

   - `nnlsh`: array of retrieved nearest neighbor candidates, `num_cand`: number of candidates.

   **varargs**

   - `k`: returns k​ nearest neighbors (BAD NOTATION, different from k used in LSH hash table).
   - `sel`: selects closest neighbors if set to "best". Otherwise, selects random neighbors. random neighbors makes query time faster, but the quality of retrieval can be pretty bad.
   - `r`: distance cutoff threshold. will not return any neighbors exceeding distance r (even for random selection).

3. `[mi,ma,me]=lshstats(T, B, xref, xtst, minNN)`

   **input**

   - `xref`: reference dataset
   - `xtst`: test query points
   - `minNN`: integer which specifies desired nearest neighbor rank m-NN. If specified, lshstat returns the number of failed query points for which the `minNN`-th nearest neighbor was not retrieved.

   **output**

   - `mi`:median bucket occupancy, `ma`:max bucket occupancy, `me`:expected bucket occupancy.