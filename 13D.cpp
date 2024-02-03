/*
Author : Hocky Yudhiono
Rab 18 Okt 2023 04:38:42
*/

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
#define endl '\n'
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

LL solve(vector <P> &red, vector <P> &blue) {
  sort(all(red)); sort(all(blue));
  int n = sz(red), m = sz(blue);
  vvi dp(n, vi(n));
  rep(i, 0, n) rep(j, i + 1, n) rep(k, 0, m) {
    if (!(red[i] < blue[k] && blue[k] < red[j])) continue;
    if (red[i].ccw(red[j], blue[k]) > 0) dp[i][j]++;
  }
  LL ans = 0;
  rep(i, 0, n) rep(j, i + 1, n) rep(k, i + 1, j) {
    ans += (dp[i][k] + dp[k][j] == dp[i][j]);
  }
  return ans;
}

int main() {
  fasterios();
  int n; cin >> n;
  int m; cin >> m;
  vector <P> red(n), blue(m);
  trav(cur, red) cin >> cur.x >> cur.y;
  trav(cur, blue) cin >> cur.x >> cur.y;
  cout << solve(red, blue) << endl;
}