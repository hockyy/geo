// D. Guardians of the Gallery

#include "bits/stdc++.h"

using namespace std;
#define rep(i,a,b) for(int i = a;i < b;i++)
#define trav(nx, v) for(auto &nx : v)
#define all(a) begin(a), end(a)
#define sz(a) (int)a.size()
#define fi first
#define se second
#define pb push_back
typedef long double LD;
typedef long long LL;
typedef pair<LL, LL> PLL;
typedef pair<int, int> PII;
const LD eps = 1e-13;
LD babs(LD x){return max(x,-x);}
template <class T> inline bool eq(T x, T y) {
  return babs(x - y) < eps;
}
template <class T> inline bool le(T x, T y){
  return x < y + eps;
}
template <class T> inline bool lt(T x, T y){
  return x + eps < y;
}
template <class T> int sgn(T x){
  return (lt(LD(0.0), x)) - (lt(x, (LD)(0.0)));
}
template <class T>
struct Point {
  typedef Point P;
  T x, y;
  bool operator <(P p) const {
    return tie(x, y) < tie(p.x, p.y);
  }
  bool operator==(P p){
    return tie(x, y) == tie(p.x, p.y);
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
  T dist2() const {
    return x * x + y * y;
  }
  LD dist() const {
    return sqrtl((LD)dist2());
  }
  T cross(P p) const {
    return x * p.y - y * p.x;
  }
  T cross(P a, P b) const {
    return(a - *this).cross(b - *this);
  }
  T ccw(P a, P b) const {
    return sgn((a - *this).cross(b - *this));
  }
  T dot(P p) const{
    return x * p.x + y*p.y;
  }
  P perp() const {
    return P(-y, x);
  }
};


typedef Point<LD> P;

pair<int, P> lineInter(P s1, P e1, P s2, P e2){
  auto d = (e1 - s1).cross(e2 - s2);
  if(eq(d, (LD)0.0)){
    return {-eq(s1.cross(e1, s2),(LD)0.0), P(0,0)};
  }
  auto p = s2.cross(e1, e2), q = s2.cross(e2, s1);
  return {1, (s1 * p + e1 * q) / d};
}



P lineProj(P a, P b, P p){
  P v = b - a;
  return p - v.perp()*v.cross(p - a) / v.dist2();
}

LD segDist(P s, P e, P p){
  if(s == e) return (p - s).dist();
  auto d = (e - s).dist2(), t = min(d, max((LD)0.0, (p-s).dot(e-s)));
  return ((p - s) * d - (e - s) * t).dist() / d;
}

LD segDist2(P s, P e, P p){
  if(s == e) return (p - s).dist2();
  auto d = (e - s).dist2(), t = min(d, max((LD)0.0, (p-s).dot(e-s)));
  return ((p - s) * d - (e - s) * t).dist2();
}

typedef pair<P, PII> Intersect;
pair<LD, P> isSafe(P a, P b, const vector <P> &isi){
  vector <Intersect> inter;
  int n = sz(isi);
  rep(i,0,n){
    auto &nx = isi[i+1 == n ? 0 : i + 1];
    auto res = lineInter(a, b, isi[i], nx);
    if(res.fi == 1 && eq(segDist2(isi[i], nx, res.se), (LD)0.0)) {
      //~ cout << a.x << " " << a.y << " got " << res.se.x << " " << res.se.y  << endl;
      inter.pb({res.se, {a.ccw(b, isi[i]), a.ccw(b, nx)}});
    }
  }
  auto vec = b - a;
  
  // Add this to stop at b
  inter.pb({b, {0,0}});
  
  sort(all(inter), [&](const Intersect&aa, const Intersect&bb){
    return lt(vec.dot(aa.fi - a), vec.dot(bb.fi - a));
  });
  int furthest = 0;
  int pre = -1; // Set to -1 if on a segment
  int inside = -1;
  rep(i,0,sz(inter)){
    // If a is NOT on segment, you can enable this
    //~ if(vec.dot(inter[i].fi - a) < 0) continue;
    //~ cout << inter[i].fi.x << " " << inter[i].fi.y << endl;
    if(inside >= 0) { furthest = i - 1;/* inter[i].fi, isi[i - 1].fi is an active segment */}
    if(inter[i].se.fi * inter[i].se.se == -1) inside = -inside;
    else if(inside == 0) inside = (inter[i].se.fi+inter[i].se.se) * pre;
    else {
      pre = (inter[i].se.fi + inter[i].se.se) * inside;
      inside = 0;
    }
    // do operations here, inside == 0 means parallel to segment
    //~ if(inside == -1) {}
    
    // add this to stop at b;
    if(inter[i].fi == b) return {(inter[i].fi - inter[furthest].fi).dist2(), b};
    if(inside == -1) {
      furthest = i + 1;
    }
  }
  return {0, b};
}

typedef pair<P, PII> Intersect;
P stab(P a, P b, const vector <P> &isi){
  vector <Intersect> inter;
  int n = sz(isi);
  rep(i,0,n){
    auto &nx = isi[(i + 1)%n];
    auto res = lineInter(a, b, isi[i], nx);
    if(res.fi == 1 && eq(segDist(isi[i], nx, res.se), (LD)0.0))
      inter.pb({res.se, {a.ccw(b, isi[i]), a.ccw(b, nx)}});
  }
  auto vec = b - a;
  
  sort(all(inter), [&](const Intersect&aa, const Intersect&bb){
    return lt(vec.dot(aa.fi - a), vec.dot(bb.fi - a));
  });
  
  set <int> existingView;
  int inside = 1; int pre;
  rep(i,0,sz(inter)){
    // If a is NOT on segment, you can enable this
    //~ cout << "Got " << inter[i].fi.x << " " << inter[i].fi.y << endl;
    if(lt(vec.dot(inter[i].fi - a), (LD)0.0)) continue;
    //~ cout << "OK " << inter[i].fi.x << " " << inter[i].fi.y << endl;
    if(inside >= 0) {/* isi[i], isi[i - 1] is an active segment */}
    if(inter[i].se.fi * inter[i].se.se == -1) inside = -inside;
    else if(inside == 0) inside = (inter[i].se.fi+inter[i].se.se) * pre;
    else {
      existingView.insert(inter[i].se.fi + inter[i].se.se);
      pre = (inter[i].se.fi + inter[i].se.se) * inside;
      inside = 0;
    }
    if(inside == -1) return inter[i].fi;
    if(existingView.size() == 2) return inter[i].fi;
  }
  return b;
}

int main(){
  int n; cin >> n;
  vector <P> isi(n);
  trav(cur, isi) {
    LL a, b; cin >> a >> b;
    cur.x = a; cur.y = b;
  }
  P guard; P diamond;
  LL tmpa, tmpb; cin >> tmpa >> tmpb;
  guard.x = tmpa; guard.y = tmpb;
  cin >> tmpa >> tmpb;
  diamond.x = tmpa; diamond.y = tmpb;
  vector <vector <LD>> dist(n + 5, vector <LD>(n + 5, 1e19));
  rep(i,0,n){
    dist[i][i] = 0;
    // tembak ini ke setiap sisi
    rep(j,i + 1,n){
      auto ret = isSafe(isi[i], isi[j], isi);
      //~ cout << i << " " << j << endl;
      //~ cout << ret.se.x << " " << ret.se.y << " " << ret.fi << endl;
      if(ret.se == isi[j] && le((isi[i] - isi[j]).dist2(), ret.fi)) {
        // ok safe
        //~ cout << i << " " << j << " safe " << endl;
        dist[i][j] = dist[j][i] = (isi[i] - isi[j]).dist();
      }
    }
    auto ret = isSafe(guard, isi[i], isi);
    if(ret.se == isi[i] && le((guard - isi[i]).dist2(), ret.fi)) {
      //~ cout << "guard" << " " << i << endl;
      dist[n][i] = dist[i][n] = (isi[i] - guard).dist();
    }
    ret = isSafe(diamond, isi[i], isi);
    if(ret.se == isi[i] && le((diamond - isi[i]).dist2(), ret.fi)) {
      //~ cout << "diamond " << " " << i << endl;
      dist[n + 1][i] = dist[i][n + 1] = (isi[i] - diamond).dist();
    }
  }
  dist[n][n] = 0;
  dist[n + 1][n + 1] = 0;
  rep(k,0,n + 2) rep(i,0, n + 2) rep(j,0,n + 2) {
    dist[i][j] = min(dist[i][j], dist[i][k] + dist[k][j]);
  }
  vector <P> viewPoints;
  //~ cout << "stabbing " << endl;
  auto tmpView = stab(diamond, guard, isi);
  viewPoints.pb(tmpView);
  //~ cout << tmpView.x << " " << tmpView.y << endl;
  //~ cout << "End " << endl;
  rep(i,0,n){
    // tembak dari diamond bikin segment garis
    auto view = stab(diamond, isi[i], isi);
    //~ cout << i << " " << view.x << " " << view.y << endl;
    viewPoints.pb(view);
  }
  LD ans = dist[n][n + 1];
  //~ rep(i,0,n + 2) {
    //~ rep(j,0,n + 2) cout << dist[i][j] << " ";
    //~ cout << endl;
  //~ }
  rep(i,0,n + 1){
    auto &cureval = i == n ? guard : isi[i];
      //~ cout << "Computing " << i << endl;
    //~ cout << cureval.x << " " << cureval.y << endl;
    trav(view, viewPoints){
      //~ cout << view.x << " " << view.y << endl;
      {
        P jato = lineProj(diamond, view, cureval);
        LD cost2 = (jato - cureval).dist2();
        //~ cout << "Using " << view.x << " " << view.y << " " << cost << " " << jato.x << " " << jato.y << endl;
        if(!lt(LD(0.0), segDist(diamond, view, jato))) {
          auto ret = isSafe(cureval, jato, isi);
          //~ cout << "OK use " << dist[n][i] << " " << cost << endl;
          if(ret.se == jato && le(cost2, ret.fi)) {
            ans = min(ans, dist[n][i] + sqrtl(cost2));
          }
        }
      }
      {
        LD cost2 = (view - cureval).dist2();
        auto ret = isSafe(cureval, view, isi);
        //~ cout << "OK Ga " << ret.se.x << " " << ret.se.y << " " << ret.fi << endl;
        //~ cout << ret.se.x << " " << ret.se.y << " " << cost << " " << ret.fi << endl;
        if(ret.se == view && le(cost2, ret.fi)) {
          //~ cout << "???? " << n << " " << i << " " << dist[n][i] << " " << cost << endl;
          ans = min(ans, dist[n][i] + sqrtl(cost2));
        }
      }
    }
  }
  cout << fixed << setprecision(12);
  long double tmp = ans;
  cout << tmp << endl;
  
}