// F. Islands from the Sky
/*
Author : Hocky Yudhiono
Sel 17 Okt 2023 11:29:42
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
    return os << p.x << " " << p.y << endl;
  }
};

typedef Point<double> PD;
vector <PD> pulau[105];
typedef pair<PD, double> P3D;
struct Pesawat {
  P3D st;
  P3D ed;
  PD perp;
};
bool bisa[105];
Pesawat pesawat[105];

typedef long double LD;

int main() {
  int n, m; cin >> n >> m;
  cout << fixed << setprecision(12);
  for (int i = 1; i <= n; i++) {
    int k; cin >> k;
    pulau[i].resize(k);
    trav(cur, pulau[i]) cin >> cur.x >> cur.y;
  }
  for (int i = 1; i <= m; i++) {
    PD &stt = pesawat[i].st.fi;
    PD &edd = pesawat[i].ed.fi;
    cin >> stt.x >> stt.y >> pesawat[i].st.se;
    cin >> edd.x >> edd.y >> pesawat[i].ed.se;
    pesawat[i].perp = (pesawat[i].ed.fi - pesawat[i].st.fi).perp().unit();
  }
  LD kiri = 0;
  const LD LIM = 1e8;
  LD kanan = LIM;
  rep(iter, 0, 100) {
    // cout << kiri << " " << kanan << endl;
    LD mid = (kiri + kanan) / 2;
    // cout << iter << " " << mid << "--binser" << endl;
    // Binser tannya
    memset(bisa, 0, sizeof(bisa));
    rep(i, 1, m + 1) {
      vector <PD> aerial;
      LD lebarSt = mid * pesawat[i].st.se;
      LD lebarEd = mid * pesawat[i].ed.se;
      auto perpSt = pesawat[i].perp * lebarSt;
      auto perpEd = pesawat[i].perp * lebarEd;
      aerial.pb({pesawat[i].st.fi + perpSt});
      aerial.pb({pesawat[i].st.fi - perpSt});
      aerial.pb({pesawat[i].ed.fi - perpEd});
      aerial.pb({pesawat[i].ed.fi + perpEd});
      auto pivot = aerial[0];
      // sort(begin(aerial) + 1, begin(aerial) + 4, [&](const PD & a, const PD & b) {
      //   return pivot.ccw(a, b) > 1;
      // });
      // cout << "Computing pesawat " << i << endl;
      // trav(cur, aerial) cout << cur << " ";
      // cout << endl;
      rep(j, 1, n + 1) {
        // cout << i << " " << j << endl;
        if (bisa[j]) continue;
        set <int> res;
        rep(dir, 0, 4) {
          trav(titik, pulau[j]) {
            // if (i == 2 && j == 3) {
            //   cout << "Here-- " << endl;
            //   cout << aerial[dir] << " " << aerial[(dir + 1) % 4] << " " << titik << endl;
            // }
            int val = aerial[dir].ccw((aerial[(dir + 1) % 4]), titik);
            res.insert(val);
          }
        }
        if (res.count(1) && res.count(-1)) {
          continue;
        }
        // cout << "here " << j << endl;
        bisa[j] = 1;
        continue;
      }
    }
    bool can = 1;
    rep(i, 1, n + 1) if (!bisa[i]) {
      // cout << i << " " << bisa[i] << endl;
      can = 0;
      break;
    }
    // cout << kiri << " " << kanan << " " << can << endl;
    if (can) kanan = mid;
    else kiri = mid;
  }
  if (kiri == LIM) cout << "impossible" << endl;
  else cout << atan(kiri) / PI * 180 << endl;
  return 0;
}