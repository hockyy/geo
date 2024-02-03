#include <bits/stdc++.h>
using namespace std;

typedef long long ll;
typedef long double ld;
typedef long long LL;
typedef pair<ll, ll> pll;
typedef pair<ll, ll> PLL;
#define fi first
#define se second
#define pb push_back
#define rep(i,a,b) for(int i = a;i < b;i++)
#define all(a) begin(a), end(a)
#define sz(a) (int) a.size()
#define trav(nx, v) for(auto &nx : v)

template <class T> int sgn(T x) {
  return (x > 0) - (x < 0);
}

const double eps = 1e-9;

template <class T>
struct Point {
  typedef Point P;
  T x, y;
  explicit Point(T _x = 0, T _y = 0) : x(_x), y(_y) {

  }
  bool operator<(const P &other)const {
    return tie(x, y) < tie(other.x, other.y);
  }
  bool operator==(const P &other)const {
    return tie(x, y) == tie(other.x, other.y);
  }
  P operator+(P p) const {
    return P(x + p.x, y + p.y);
  }
  P operator-(P p) const {
    return P(x - p.x, y - p.y);
  }
  P operator*(T d) const {
    return P(x * d, y * d);
  }
  P operator/(T d) const {
    return P(x / d, y / d);
  }
  T cross(P p) const {
    return x * p.y - y * p.x;
  }
  T dot(P p) const {
    return x * p.x + y * p.y;
  }
  T cross(P a, P b) const {
    return (a - *this).cross(b - (*this));
  }
  T dist2() const {
    return x * x + y * y;
  }
  double dist() const {
    return sqrt(double(dist2()));
  }
  P unit() const {
    return *this / dist();
  }
};

template <class P> bool segInter(P a, P b, P c, P d) {
  auto oa = c.cross(d, a), ob = c.cross(d, b);
  auto oc = a.cross(b, c), od = a.cross(b, d);
  if (sgn(oa) * sgn(ob) < 0 && sgn(oc) * sgn(od) < 0)
    return 1;
  return 0;
}

typedef Point<double> PD;
double segDist(PD& s, PD & e, PD &p) {
  if (s == e) return (p - s).dist();
  auto d = (e - s).dist2(), t = min(d, max(.0, (p - s).dot(e - s )));
  return ((p - s) * d - (e - s) * t).dist() / d;
}

template<class P> bool onSegment(P s, P e, P p) {
  return segDist(s, e, p) < eps;
}

template<class P>
bool inPolygon(const vector <P> &p, P a, bool strict = true) {
  int cnt = 0, n = sz(p);
  rep(i, 0, n) {
    P q = p[(i + 1) % n];
    if (onSegment(p[i], q, a)) return !strict;
    cnt ^= ((a.y < p[i].y) - (a.y < q.y)) * a.cross(p[i], q) > 0;
  }
  return cnt;
}

int isBad[205][205];

template<class P>
int solve(int u, int v, const vector <P> &isi) {
  if (u == v) return 0;
  if (u > v) return solve(v, u, isi);
  int &ret = isBad[u][v];
  if (~ret) return ret;
  int n = sz(isi);
  rep(k, 0, n) {
    if (k == u || k == v || (k + 1) % n == u || (k + 1) % n == v) continue;
    auto &a = isi[k], &b = isi[(k + 1) % n];
    auto &c = isi[u], &d = isi[v];
    if (onSegment(c, d, a)) return ret = solve(u, k, isi) || solve(k, v, isi);
    if (onSegment(c, d, b)) return ret = solve(u, (k + 1) % n, isi) || solve((k + 1) % n, v, isi);
    auto res = segInter(a, b, c, d);
    if (res) return ret = 1;
  }

  return ret = !inPolygon(isi, (isi[u] + isi[v]) / 2, false);
}

const ll MOD = 998244353;
vector <int> edge[205];

inline void add(int &a, int &b) {
  a = a + b >= MOD ? a + b - MOD : a + b;
}

typedef pair<int, int> PII;
typedef pair<PD, PII> Slope;
typedef array<int, 3> Container;
const int LIM = 200;
LL memo[2][LIM + 5][LIM + 5];
int main() {
  memset(isBad, -1, sizeof(isBad));
  cin.tie(0)->sync_with_stdio(0);
  cout.tie(0);
  int n;
  cin >> n;
  vector <PD> isi(n);
  trav(cur, isi) cin >> cur.x >> cur.y;
  rep(i, 0, n) rep(j, 0, n) isBad[i][j] = isBad[j][i] = solve(i, j, isi);
  // rep(i,0,n) rep(j,0,n) cout << i << " " << j << " " << isBad[i][j] << endl;
  vector <int> perm(n);
  iota(all(perm), 0);
  sort(all(perm), [&](int a, int b) {
    return isi[a] < isi[b];
  });
  vector <int> pow2(LIM * 2 + 5);
  pow2[0] = 1;
  rep(i, 1, LIM * 2 + 5) {
    pow2[i] = pow2[i - 1] << 1;
    if (pow2[i] >= MOD) pow2[i] -= MOD;
  }
  LL ans = 0;
  vector<Slope> slopes;
  for (int i = n - 1; i >= 0; i--) rep(j, i + 1, n) {
    if (!isBad[perm[i]][perm[j]]) {
      slopes.pb({isi[perm[j]] - isi[perm[i]], {perm[i], perm[j]}});
      int sameGradient = 0;
      rep(k, i + 1, j) {
        if (isi[perm[k]].cross(isi[perm[i]], isi[perm[j]]) == 0) {
          sameGradient++;
        }
      }
      ans = ans + MOD - pow2[sameGradient * 2] + pow2[sameGradient];
      while (ans >= MOD) ans -= MOD;
    }
  }
  rep(i, 0, n) memo[0][i][i] = memo[1][i][i] = 1;
  rep(iter, 0, 2) {
    // Sort ccw
    stable_sort(all(slopes), [&](const Slope & a, const Slope & b) {
      return a.fi.cross(b.fi) > 0;
    });
    trav(slope, slopes) {
      int i = slope.se.fi, j = slope.se.se;
      if(iter) swap(i, j);
      rep(k, 0, n) {
        memo[iter][i][k] = (memo[iter][i][k] + memo[iter][j][k]) % MOD;
      }
    }
    reverse(all(slopes));
  }
  rep(i, 0, n) rep(j, 0, n) {
    ans += memo[0][i][j] * memo[1][j][i] % MOD;
  }
  cout << (ans - n + MOD) % MOD << endl;
  return 0;
}