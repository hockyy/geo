/*
Author : Hocky Yudhiono
9/14/2021 1:06:46 PM

1. You can sort the query if offline!
2. Don't bring the dp remaining state when dfsing on DP on Tree.
3. Try to reverse (Think from the back) if you stuck.
4. Be careful when submitting Div. 2 D-F, dont waste it on stupid WAs.
5. Try to reduce a problem, think of it when you're making subtasks
   like when problemsetting.
6. Matching problems can be solved with DP and vice versa.
   Counting and optimizing problems can be solved with DP.
   Try bitmasking when N is small. When big, consider greedy-ing.

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
#define remove erase
#define fi first
#define se second

// If the time limit is strict, try not to use long double

// Remember to undefine if the problem is interactive
#define endl '\n'
#define DEBUG(X) cout << ">>> DEBUG(" << __LINE__ << ") " << #X << " = " << (X) << endl

const double EPS = 1e-9;
const int INFMEM = 63;
const int INF = 1061109567;
const ll LINF = 4557430888798830399LL;
const double DINF = numeric_limits<double>::infinity();
const ll MOD = 1000000007;
const int dx[8] = {0,0,1,-1,1,-1,1,-1};
const int dy[8] = {1,-1,0,0,1,-1,-1,1};
// Do dir^1 to get reverse direction
const char dch[4] = {'R','L','D','U'};
// const string ds[8] = {"E","W","S","N","SE","NW","SW","NE"};
const double PI = 3.141592653589793;

inline void open(string a){
  freopen((a+".in").c_str(),"r",stdin);
  freopen((a+".out").c_str(),"w",stdout);
}

inline void fasterios(){
  // Do not use if interactive
  ios_base::sync_with_stdio(0);
  cin.tie(0); cout.tie(0);
  // cout << fixed << setprecision(10);
}


template<class T> inline bool eq(T x, T y) { return fabs(x-y) < EPS; }
template<class T> inline bool le(T x, T y) { return x < y + EPS; }
template<class T> inline bool lt(T x, T y) { return x + EPS < y; }
template<class T> inline T doubleMax(T x, T y) { return lt(x,y) ? y : x; }
template<class T> inline T doubleMin(T x, T y) { return lt(x,y) ? x : y; }
template <class T> int sgn(T x) { return (x > 0) - (x < 0); }
template<class T>
struct Point {
    typedef Point P;
    T x, y;
    explicit Point(T X=0, T Y=0) : x(X), y(Y) {}
    bool operator<(P p) const { return tie(x,y) < tie(p.x,p.y); }
    bool operator==(P p) const { return tie(x,y)==tie(p.x,p.y); }
    P operator+(P p) const { return P(x+p.x, y+p.y); }
    P operator-(P p) const { return P(x-p.x, y-p.y); }
    P operator*(T d) const { return P(x*d, y*d); }
    P operator/(T d) const { return P(x/d, y/d); }
    T dot(P p) const { return x*p.x + y*p.y; }
    T cross(P p) const { return x*p.y - y*p.x; }
    T cross(P a, P b) const { return (a-*this).cross(b-*this); }
    // return the orientation of (*this,a,b), 1 if ccw
    int ccw(P a, P b) const { return sgn((a-*this).cross(b-*this)); }
    T dist2() const { return x*x + y*y; }
    double dist() const { return sqrt((double)dist2()); }
    // angle to x-axis in interval [-pi, pi]
    double angle() const { return atan2(y, x); }
    P unit() const { return *this/dist(); } // makes dist()=1
    P perp() const { return P(-y, x); } // rotates +90 degrees
    P normal() const { return perp().unit(); }
    // returns point rotated 'a' radians ccw around the origin
    P rotate(double a) const {
        return P(x*cos(a)-y*sin(a),x*sin(a)+y*cos(a)); }
    friend ostream& operator<<(ostream& os, P p) {
        return os << "(" << p.x << "," << p.y << ")"; }
};

template<class P>
P lineProj(P a, P b, P p, bool refl=false) {
  P v = b - a;
  return p - v.perp()*(1+refl)*v.cross(p-a)/v.dist2();
}

template<class P>
int sideOf(P s, P e, P p) { return sgn(s.cross(e, p)); }

template<class P>
int sideOf(const P& s, const P& e, const P& p, double eps) {
  auto a = (e-s).cross(p-s);
  double l = (e-s).dist()*eps;
  return (a > l) - (a < -l);
}

typedef Point<double> P;
double rat(P a, P b) { return sgn(b.x) ? a.x/b.x : a.y/b.y; }
double polyUnion(vector<vector<P>>& poly) {
  double ret = 0;
  rep(i,0,sz(poly)) rep(v,0,sz(poly[i])) {
    P A = poly[i][v], B = poly[i][(v + 1) % sz(poly[i])];
    vector<pair<double, int>> segs = {{0, 0}, {1, 0}};
    rep(j,0,sz(poly)) if (i != j) {
      rep(u,0,sz(poly[j])) {
        P C = poly[j][u], D = poly[j][(u + 1) % sz(poly[j])];
        int sc = sideOf(A, B, C), sd = sideOf(A, B, D);
        if (sc != sd) {
          double sa = C.cross(D, A), sb = C.cross(D, B);
          if (min(sc, sd) < 0)
            segs.emplace_back(sa / (sa - sb), sgn(sc - sd));
        } else if (!sc && !sd && j<i && sgn((B-A).dot(D-C))>0){
          segs.emplace_back(rat(C - A, B - A), 1);
          segs.emplace_back(rat(D - A, B - A), -1);
        }
      }
    }
    sort(all(segs));
    for (auto& s : segs) s.first = min(max(s.first, 0.0), 1.0);
    double sum = 0;
    int cnt = segs[0].second;
    rep(j,1,sz(segs)) {
      if (!cnt) sum += segs[j].first - segs[j - 1].first;
      cnt += segs[j].second;
    }
    ret += A.cross(B) * sum;
  }
  return fabs(ret) / 2;
}

template<class T>
T polygonArea2(vector<Point<T>>& v) {
  T a = v.back().cross(v[0]);
  rep(i,0,sz(v)-1) a += v[i].cross(v[i+1]);
  return a;
}

int main(){
  fasterios();
  int n; cin >> n;
  vector <P> poly1(n);
  trav(cur, poly1) cin >> cur.x >> cur.y;
  P st; P ed; cin >> st.x >> st.y >> ed.x >> ed.y;
  vector <P> poly2 = poly1;
  reverse(all(poly2));
  trav(cur, poly2) cur = lineProj(st, ed, cur, true);
  vector <vector<P>> allpoly = {poly1, poly2};
  double ans = fabs(polygonArea2(poly2)); 
  // cout << ans << " " << polyUnion(allpoly) << endl;
  cout << fixed << setprecision(10) << (ans - polyUnion(allpoly)) / 2.0 << endl;
}