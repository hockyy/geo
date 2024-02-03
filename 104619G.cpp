/*
Author : Hocky Yudhiono
Rab 18 Okt 2023 07:22:01
*/
#pragma GCC optimize("no-stack-protector,O3,unroll-loops")
#pragma GCC target("sse,sse2,sse3,ssse3,sse4,popcnt,abm,mmx,avx,avx2,tune=native")
#include <bits/stdc++.h>
using namespace std;

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
#define trav(a, x) for(auto& a : x)
#define all(x) begin(x), end(x)
#define sz(x) (int)(x).size()
#define popf pop_front
#define pf push_front
#define popb pop_back
#define pb push_back
#define fi first
#define se second

const double EPS = 1e-9;
const int INFMEM = 63;

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
// #define endl '\n'
const int MOD = 1000000007;
// const int MOD = 998244353;

template<class T> inline bool eq(T x, T y) { return fabs(x - y) < EPS; }
template<class T> inline bool le(T x, T y) { return x < y + EPS; }
template<class T> inline bool lt(T x, T y) { return x + EPS < y; }
template<class T> inline T doubleMax(T x, T y) { return lt(x, y) ? y : x; }
template<class T> inline T doubleMin(T x, T y) { return lt(x, y) ? x : y; }
template <class T> int sgn(T x) { return (x > 0) - (x < 0); }
template<class T>
struct Point {
  typedef Point P;
  T x, y;
  explicit Point(T X = 0, T Y = 0) : x(X), y(Y) {}
  bool operator<(P p) const { return tie(x, y) < tie(p.x, p.y); }
  bool operator==(P p) const { return tie(x, y) == tie(p.x, p.y); }
  P operator+(P p) const { return P(x + p.x, y + p.y); }
  P operator-(P p) const { return P(x - p.x, y - p.y); }
  P operator*(T d) const { return P(x * d, y * d); }
  P operator/(T d) const { return P(x / d, y / d); }
  T dot(P p) const { return x * p.x + y * p.y; }
  T cross(P p) const { return x * p.y - y * p.x; }
  T cross(P a, P b) const { return (a - *this).cross(b - *this); }
  // return the orientation of (*this,a,b), 1 if ccw
  int ccw(P a, P b) const { return sgn((a - *this).cross(b - *this)); }
  T dist2() const { return x * x + y * y; }
  double dist() const { return sqrt((double)dist2()); }
  // angle to x-axis in interval [-pi, pi]
  double angle() const { return atan2(y, x); }
  P unit() const { return *this / dist(); } // makes dist()=1
  P perp() const { return P(-y, x); } // rotates +90 degrees
  P normal() const { return perp().unit(); }
  // returns point rotated 'a' radians ccw around the origin
  P rotate(double a) const {
    return P(x * cos(a) - y * sin(a), x * sin(a) + y * cos(a));
  }
  friend ostream& operator<<(ostream& os, P p) {
    return os << "(" << p.x << "," << p.y << ")";
  }
};

typedef Point<LL> P;
typedef pair<P, bool> Gadget;
typedef pair<P, PII> Slope;
typedef array<int, 3> Container;
const int LIM = 500;

LL cntTop[LIM + 5][LIM + 5];
LL memo[2][LIM + 5][LIM + 5];
int cntIn(int i, int j, int k) {
  if (j > k) swap(j, k);
  if (i > j) swap(i, j);
  assert(i <= j && j <= k);
  int cnt = cntTop[i][j] + cntTop[j][k] - cntTop[i][k];
  return cnt >= 0 ? cnt : -cnt - 1;
}

int main() {
  cin.tie(0)->sync_with_stdio(0);
  cout.tie(0);
  int n; cin >> n;
  vector<Gadget> isi(n);
  trav(cur, isi) {
    cin >> cur.fi.x >> cur.fi.y;
    char ch; cin >> ch;
    cur.se = ch == 'b';
  }
  sort(all(isi));
  vector<Slope> slopes;
  rep(i, 0, n) rep(j, i + 1, n) {
    if (isi[i].se != isi[j].se) {
      slopes.pb({isi[j].fi - isi[i].fi, {i, j}});
    }
  }
  // Sort ccw
  sort(all(slopes), [&](const Slope & a, const Slope & b) {
    return a.fi.cross(b.fi) > 0;
  });

  vector <LL> pow2(LIM + 5, 1);
  rep(i, 1, LIM + 1) {
    pow2[i] = pow2[i - 1] << 1;
    if (pow2[i] >= MOD) pow2[i] -= MOD;
  }
  rep(i, 0, n) rep(j, i + 1, n) rep(k, i + 1, j) {
    if (isi[i].fi.ccw(isi[j].fi, isi[k].fi) > 0)
      cntTop[i][j]++;
  }
  rep(i, 0, n) memo[0][i][i] = memo[1][i][i] = 1;
  trav(slope, slopes) {
    int i = slope.se.fi, j = slope.se.se;
    rep(k, 0, n) {
      // cout << "Computing " << i << " " << j << " " << k << endl;
      // cout << isi[i].fi << " " << isi[j].fi << " " << isi[k].fi << " " << memo[0][i][k] << " += " << memo[0][j][k] << " * 2^" << cntIn(i, j, k) << endl;
      int cost = pow2[cntIn(i, j, k)];
      memo[0][i][k] = (memo[0][i][k] + memo[0][j][k] * cost) % MOD;
      memo[1][j][k] = (memo[1][j][k] + memo[1][i][k] * cost) % MOD;
    }
  }
  LL ans = 0;
  rep(i, 0, n) rep(j, 0, n) {
    ans += memo[0][i][j] * memo[1][j][i] % MOD;
    // cout << i << " " << j << " " << memo[0][i][j] * memo[1][j][i] << endl;
  }
  cout << (ans - n - sz(slopes) + MOD) % MOD << endl;
}