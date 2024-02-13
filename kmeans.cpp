/*
Author : Hocky Yudhiono
Tue Feb 13 19:19:58 2024

1. You can sort the query if offline!
2. Don't bring the dp remaining state when dfsing on DP on Tree.
3. Try to reverse (Think from the back) if you stuck.
4. Be careful when submitting Div. 2 D-F, dont waste it on stupid WAs.
5. Try to reduce a problem, think of it when you're making subtasks
   like when problemsetting.
6. Matching problems can be solved with DP and vice versa.
   Counting and optimizing problems can be solved with DP.
   Try bitmasking when N is small. When big, consider greedy-ing.
7. map<string,int> is faster than you think

*/

#include <algorithm>
#include <iostream>
#include <numeric>
#include <cstdlib>
#include <cassert>
#include <cstring>
#include <iomanip>
#include <cstdio>
#include <limits>
#include <string>
#include <vector>
#include <cmath>
#include <deque>
#include <queue>
#include <stack>
#include <map>
#include <set>
using namespace std;

// Suffix Literal:
// U: Unsigned
// L: Long double
// LL: Long long
typedef long long ll;
typedef long long LL;
typedef vector<int> vi;
typedef vector<ll> vl;
typedef vector<vi> vvi;
typedef vector<vl> vvl;
typedef pair<int, int> PII;
typedef pair<int, int> pii;
typedef pair<ll, ll> PLL;
typedef pair<ll, ll> pll;
typedef long double ld;

#define rep(i, a, b) for(int i = a; i < (b); ++i)
#define rrep(i, a, b) for(int i = a; i > (b); --i)
#define trav(a, x) for(auto& a : x)
#define all(x) begin(x), end(x)
#define sz(x) (int)(x).size()
#define popf pop_front
#define pf push_front
#define popb pop_back
#define mp make_pair
#define pb push_back
#define fi first
#define se second

// If the time limit is strict, try not to use long double

const double EPS = 1e-9;
const int INFMEM = 63;
// inf constant close to, but not 2^30 - 1
// the value is 0x3f3f... and aligns with INFMEM
// Can replace with INT_MAX or LLONG_MAX for 2^31-1 and 2^63-1
const int INF = 1061109567;
const ll LINF = 4557430888798830399LL;
const double DINF = numeric_limits<double>::infinity();
const int MOD = 1000000007;
// const int MOD = 998244353;

// Do dir^1 to get reverse direction
const int dx[8] = {0, 0, 1, -1, 1, -1, 1, -1};
const int dy[8] = {1, -1, 0, 0, 1, -1, -1, 1};
const char dch[4] = {'R', 'L', 'D', 'U'};

// Do (dir + 2)%4 to get reverse direction
// const int dx[8] = {-1,0,1,0,-1,1,1,-1};
// const int dy[8] = {0,1,0,-1,1,1,-1,-1};
// const char dch[4] = {'U','R','D','L'};
const double PI = 3.141592653589793;

inline void fasterios() {
  cin.tie(0)->sync_with_stdio(0);
  cin.exceptions(cin.failbit);
}

// Remember to undefine if the problem is interactive
#define endl '\n'

// In multi dimensional vector, avoid using
// vector.resize() to reset between each testcases
// use assign instead
typedef double LD;
struct Point {
  LD x, y, z;
  Point() : x(0), y(0), z(0) {};
  LD dist2(const Point &other) const {
    return (x - other.x) * (x - other.x) + (y - other.y) * (y - other.y)
           + (z - other.z) * (z - other.z);
  }
  void add(const Point &other) {
    x += other.x;
    y += other.y;
    z += other.z;
  }
  void div(const long long divisor) {
    x /= divisor;
    y /= divisor;
    z /= divisor;
  }
};


#include <random>
#include <chrono>
// mt19937 rng(chrono::steady_clock::now().time_since_epoch().count()); //For int
mt19937_64 rng(chrono::steady_clock::now().time_since_epoch().count()); //For LL
// cout << rng() << endl;
// shuffle(isi.begin(),isi.end(),rng);

LD fullBf(const vector <Point> &points) {
  LD ans = 99999999;
  for (int i = 0, e = (1 << sz(points)); i < e; i++) {
    vector <Point> centroids(2);
    vector <int> divisors(2);
    rep(j, 0, sz(points)) {
      centroids[(i >> j) & 1].add(points[j]);
      divisors[(i >> j) & 1]++;
    }
    if(divisors[0]) centroids[0].div(divisors[0]);
    if(divisors[1]) centroids[1].div(divisors[1]);
    LD curans = 0;
    rep(j,0,sz(points)){
      curans += centroids[(i >> j) & 1].dist2(points[j]);
    }
    ans = min(ans, curans);
  }
  return ans;
}

LD currentCost(const vector <Point> &points, int pivot, int pivot2) {
  vector <Point> centroids(2);
  // do here
  PLL tmPair = {pivot, pivot2};
  if (tmPair.se >= tmPair.fi) tmPair.se++;
  centroids[0] = points[tmPair.fi];
  centroids[1] = points[tmPair.se];
  rep(tmp, 0, 100) {
    vector <Point> newCentroids(2);
    vector <int> divisors(2);
    trav(cur, points) {
      if (centroids[0].dist2(cur) < centroids[1].dist2(cur)) {
        newCentroids[0].add(cur);
        divisors[0]++;
      } else {
        newCentroids[1].add(cur);
        divisors[1]++;
      }
    }
    newCentroids[0].div(divisors[0]);
    newCentroids[1].div(divisors[1]);
    centroids = move(newCentroids);
    // cout << newCentroids[0].x << " " << newCentroids[0].y << " " << newCentroids[0].z << endl;
  }
  LD tmpans = 0;
  trav(cur, points) {
    if (centroids[0].dist2(cur) < centroids[1].dist2(cur)) {
      tmpans += centroids[0].dist2(cur);
    } else {
      tmpans += centroids[1].dist2(cur);
    }
  }
  return tmpans;
}

int main() {
  cout << fixed << setprecision(6);
  int n, k; cin >> n >> k;
  vector <Point> points(n);
  trav(cur, points) cin >> cur.x >> cur.y >> cur.z;
  if (n <= k) {
    cout << 0.0 << endl;
    return 0;
  }
  Point avg;
  if (k == 1) {
    trav(cur, points) {
      avg.x += cur.x;
      avg.y += cur.y;
      avg.z += cur.z;
    }
    avg.x /= n;
    avg.y /= n;
    avg.z /= n;
    LD ans = 0;
    trav(cur, points) {
      ans += (avg.dist2(cur));
    }
    cout << ans << endl;
    return 0;
  }
  LD ans = currentCost(points, 0, 1);
  rep(i,0,n) rep(j,0,n) {
    if(i == j) continue;
    LD curans = (currentCost(points, i, j));
    ans = min(ans, curans);
    // cout << ans << " " << curans << endl;
  }
  cout << ans << endl;
}