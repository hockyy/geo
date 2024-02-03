// Spaceship exploration
#pragma GCC optimize("O3,unroll-loops")
#pragma GCC target("avx2,popcnt,lzcnt,abm,bmi,bmi2,fma,tune=native")
#include "bits/stdc++.h"
using namespace std;

#define rep(i,a,b) for(int i = a;i < b;i++)
#define sz(a) (int) a.size()
#define all(a) begin(a), end(a)
#define trav(nx, v) for(auto &nx : v)
#define pb push_back
#define fi first
#define se second
#define endl '\n'

typedef long double LD;
const LD EPS = 1e-9;

typedef long long LL;
template<class T> inline bool eq(T x, T y) { return fabs((long double)(x - y)) < EPS; }
template<class T> inline bool le(T x, T y) { return x < y + EPS; }
template<class T> inline bool lt(T x, T y) { return x + EPS < y; }
template<class T> inline bool gt(T x, T y) { return lt(y, x); }
template<class T> inline bool ge(T x, T y) { return le(y, x); }
template<class T> inline T doubleMax(T x, T y) { return lt(x, y) ? y : x; }
template<class T> inline T doubleMin(T x, T y) { return lt(x, y) ? x : y; }
template <class T> int sgn(T x) { return gt(x, T(0)) - lt(x, T(0)); }
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
  LD dist() const { return sqrtl((LD)dist2()); }
  // angle to x-axis in interval [-pi, pi]
  double angle() const { return atan2(y, x); }
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
#define cmp(i,j) sgn(dir.perp().cross(poly[(i)%n]-poly[(j)%n]))
#define extr(i) cmp(i + 1, i) >= 0 && cmp(i, i - 1 + n) < 0
template <class P> int extrVertex(vector<P>& poly, P dir) {
  int n = sz(poly), lo = 0, hi = n;
  if (extr(0)) return 0;
  while (lo + 1 < hi) {
    int m = (lo + hi) / 2;
    if (extr(m)) return m;
    int ls = cmp(lo + 1, lo), ms = cmp(m + 1, m);
    (ls < ms || (ls == ms && ls == cmp(lo, m)) ? hi : lo) = m;
  }
  return lo;
}

#define cmpL(i) sgn(a.cross(poly[i], b))
template <class P>
array<int, 2> lineHull(P a, P b, vector<P>& poly) {
  int endA = extrVertex(poly, (a - b).perp());
  int endB = extrVertex(poly, (b - a).perp());
  if (cmpL(endA) < 0 || cmpL(endB) > 0)
    return { -1, -1};
  array<int, 2> res;
  rep(i, 0, 2) {
    int lo = endB, hi = endA, n = sz(poly);
    while ((lo + 1) % n != hi) {
      int m = ((lo + hi + (lo < hi ? 0 : n)) / 2) % n;
      (cmpL(m) == cmpL(endB) ? lo : hi) = m;
    }
    res[i] = (lo + !cmpL(hi)) % n;
    swap(endA, endB);
  }
  if (res[0] == res[1]) return {res[0], -1};
  if (!cmpL(res[0]) && !cmpL(res[1]))
    switch ((res[0] - res[1] + sz(poly) + 1) % sz(poly)) {
    case 0: return {res[0], res[0]};
    case 2: return {res[1], res[1]};
    }
  return res;
}


typedef Point<LL> PL;
// typedef long double LD;
template <class P>
double segDist(P& s, P& e, P& p) {
  if (s == e) return (p - s).dist();
  auto d = (e - s).dist2();
  auto t = min((LD) d, max((LD).0L, (LD) (p - s).dot(e - s)));
  return ((p - s) * d - (e - s) * t).dist() / d;
}

template<class P> bool onSegment(P s, P e, P p) {
  return fabs(segDist(s, e, p)) < EPS;
}

typedef Point<LD> PD;
template<class P>
pair<int, PD> lineInter(P s1, P e1, P s2, P e2) {
  auto d = (e1 - s1).cross(e2 - s2);
  if (eq(d, (LD)0.0)) // if parallel
    return { eq(-s1.cross(e1, s2), (LD)0), PD(0, 0)};
  auto p = s2.cross(e1, e2), q = s2.cross(e2, s1);
  return {1, (s1 * p + e1 * q) / d};
}

template<class P> vector<P> segInter(P a, P b, P c, P d) {
  auto oa = c.cross(d, a), ob = c.cross(d, b),
       oc = a.cross(b, c), od = a.cross(b, d);
  // Checks if intersection is single non-endpoint point.
  if (sgn(oa) * sgn(ob) < 0 && sgn(oc) * sgn(od) < 0)
    return {(a * ob - b * oa) / (ob - oa)};
  set<P> s;
  if (onSegment(c, d, a)) s.insert(a);
  if (onSegment(c, d, b)) s.insert(b);
  if (onSegment(a, b, c)) s.insert(c);
  if (onSegment(a, b, d)) s.insert(d);
  return {all(s)};
}

PD toPD(PL a) {
  return PD(a.x, a.y);
}
#define endl '\n'
const int stRep = 0;
const int edRep = 1;
int main() {
  cin.tie(0)->sync_with_stdio(0);
  cout.tie(0);
  cout << fixed << setprecision(15);
  int n; cin >> n;
  // if(n != 5 && n != 4 && n != 8) assert(0);
  vector <PD> poly(n);
  trav(cur, poly) {
    PL tmp; cin >> tmp.x >> tmp.y;
    cur = toPD(tmp);
  }
  // reverse(all(poly));
  int q; cin >> q;

  // auto isTwoIntersection = [&](PD & a, PD & b) {
  //   auto res = lineHull(a, b, poly);
  //   // cout << res[0] << " " << res[1] << endl;
  //   if (res[1] != -1 && res[0] != res[1]) return 1;
  //   return 0;
  // };
  auto isTwoIntersection = [&](PD & a, PD & b) {
    auto res = lineHull(a, b, poly);
    // cout << res[0] << " " << res[1] << endl;
    if (res[1] != -1 && res[0] != res[1]) return 1;
    return 0;
  };

  auto binserGo = [&](PD & titik, int kiri, int kanan) {
    while (kiri < kanan) {
      int mid = (kiri + kanan) >> 1;
      auto a = mid % n;
      auto b = (mid + 1) % n;
      if (titik.ccw(poly[a], poly[b]) == 1) kanan = mid;
      else kiri = mid + 1;
    }
    return kiri;
  };
  auto binserCome = [&](PD & titik, int kiri, int kanan) {
    while (kiri < kanan) {
      int mid = (kiri + kanan + 1) >> 1;
      auto a = mid % n;
      auto b = (mid - 1 + n) % n;
      // cout << mid%n << " " << titik.ccw(poly[a], poly[b]) << endl;
      if (titik.ccw(poly[a], poly[b]) == -1) kiri = mid;
      else kanan = mid - 1;
    }
    return kiri;
  };

  auto binserGo2 = [&](PD & titik, int kiri, int kanan) {
    while (kiri < kanan) {
      int mid = (kiri + kanan) >> 1;
      auto a = mid % n;
      auto b = (mid + 1) % n;
      if (titik.ccw(poly[a], poly[b]) == 1) kanan = mid;
      else kiri = mid + 1;
    }
    return kiri;
  };
  auto binserCome2 = [&](PD & titik, int kiri, int kanan) {
    while (kiri < kanan) {
      int mid = (kiri + kanan + 1) >> 1;
      auto a = mid % n;
      auto b = (mid - 1 + n) % n;
      // cout << mid%n << " " << titik.ccw(poly[a], poly[b]) << endl;
      if (titik.ccw(poly[a], poly[b]) == -1) kiri = mid;
      else kanan = mid - 1;
    }
    return kiri;
  };

  while (q--) {
    PL st, ed; cin >> st.x >> st.y >> ed.x >> ed.y;
    PD std, edd;
    std = toPD(st);
    edd = toPD(ed);
    auto res = lineHull(std, edd, poly);
    // cout << "Here intersection " << res[0] << " " << res[1] << endl;
    if (res[1] == -1 || res[0] == res[1]) {
      cout << (long double)(st - ed).dist() << endl;
    } else {
      auto &aa = poly[res[0]];
      auto &bb = poly[(res[0] + 1) % n];
      auto &cc = poly[res[1]];
      auto &dd = poly[(res[1] + 1) % n];
      auto inter1 = segInter(aa, bb, std, edd);
      auto inter2 = segInter(cc, dd, std, edd);
      inter1.insert(end(inter1), all(inter2));
      sort(all(inter1));
      inter1.erase(unique(all(inter1)), inter1.end());
      if (sz(inter1) <= 1) {
        cout << (long double)(st - ed).dist() << endl;
        continue;
      }
      const LD DINF = HUGE_VALL;
      LD ans = DINF;
      // cout << "Yay " << st << " " << ed << endl;
      vector <int> stIdx;
      vector <int> edIdx;
      {
        auto lo = res[0] + 1 + n;
        auto hi = res[1] + n;
        if (lo > hi) hi += n;
        // if (poly[hi % n] == edd) hi--;
        // cout << "Binsering " << lo % n << " " << hi % n << endl;
        int kiri = binserGo(std, lo, hi);
        int kanan = binserCome(edd, lo, hi);
        kiri %= n, kanan %= n;
        // cout << kiri << " " << kanan << endl;

        rep(i, stRep, edRep) {
          if (!isTwoIntersection(std, poly[(kiri + i + n) % n])) {
            stIdx.pb((kiri + i + n) % n);
          }
        }
        rep(i, stRep, edRep) {
          if (!isTwoIntersection(edd, poly[(kanan + i + n) % n])) {
            edIdx.pb((kanan + i + n) % n);
          }
        }
        // cout << kiri << " " << kanan << endl;
      }
      {
        auto lo = res[1] + 1 + n;
        auto hi = res[0] + n;
        if (lo > hi) hi += n;
        // if (poly[(lo - 1) % n] == std) lo--;
        // if (poly[hi % n] == edd) hi--;
        // cout << "Binsering " << lo % n << " " << hi % n << endl;
        int kiri = binserCome2(std, lo, hi);
        int kanan = binserGo2(edd, lo, hi);

        kiri %= n, kanan %= n;
        // cout << "Binser come 2 " << std << " " << lo << " " << hi << " " << kiri << " " << kanan << endl;


        rep(i, stRep, edRep) {
          if (!isTwoIntersection(std, poly[(kiri + i + n) % n])) {
            stIdx.pb((kiri + i + n) % n);
          }
        }
        rep(i, stRep, edRep) {
          if (!isTwoIntersection(edd, poly[(kanan + i + n) % n])) {
            edIdx.pb((kanan + i + n) % n);
          }
        }
      }
      sort(all(stIdx));
      stIdx.erase(unique(all(stIdx)), stIdx.end());
      sort(all(edIdx));
      edIdx.erase(unique(all(edIdx)), edIdx.end());
      trav(tangent1, stIdx) {
        trav(tangent2, edIdx) {
          // cout << tangent1 << " " << tangent2 << endl;
          auto ret = lineInter(std, poly[tangent1], edd, poly[tangent2]);
          if (ret.fi > 0) {
            LD curans =  ((std - ret.se).dist() + (edd - ret.se).dist());
            ans = min(ans, curans);
          }
        }
      }
      // cout << log10(ans) << endl;
      // Get intersection of kiri and kanan
      if (ans == DINF) {
        cout << -1 << endl;
      } else {
        cout << (long double) ans << endl;
      }
    }
  }
}