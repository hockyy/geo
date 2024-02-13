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
template <class T> int sgn(T x) { return (x > 0) - (x < 0); }
template<class T>
struct Point {
    typedef Point P;
    T x, y;
    bool operator<(P p) const { return tie(x,y) < tie(p.x,p.y); }
    bool operator==(P p) const { return tie(x,y)==tie(p.x,p.y); }
    P operator+(P p) const { return P(x+p.x, y+p.y); }
    P operator-(P p) const { return P(x-p.x, y-p.y); }
    P operator*(T d) const { return P(x*d, y*d); }
    P operator/(T d) const { return P(x/d, y/d); }
    // Same direction means a.b > 0, when perp, a.b = 0
    T dot(P p) const { return x*p.x + y*p.y; } // a.b = |a|.|b|cos(t)
    // b left of a, a^b > 0, ccw sort comparator. When perp a.b = 1
    T cross(P p) const { return x*p.y - y*p.x; } // a^b = |a|.|b|sin(t)
    T cross(P a, P b) const { return (a-*this).cross(b-*this); }
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
pair<int, P> lineInter(P s1, P e1, P s2, P e2) {
    auto d = (e1 - s1).cross(e2 - s2);
    if (d == 0) // if parallel
        return {-(s1.cross(e1, s2) == 0), P(0, 0)};
    auto p = s2.cross(e1, e2), q = s2.cross(e2, s1);
    return {1, (s1 * p + e1 * q) / d};
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
typedef array<P, 2> Line;
#define sp(a) a[0], a[1]
#define ang(a) (a[1] - a[0]).angle()

int angDiff(Line a, Line b) { return sgn(ang(a) - ang(b)); }
bool cmp(Line a, Line b) {
    int s = angDiff(a, b);
    return (s ? s : sideOf(sp(a), b[0])) < 0;
}
vector<P> halfPlaneIntersection(vector<Line> vs) {
    const double EPS = sqrt(2) * 1e-8;
    sort(all(vs), cmp);
    vector<Line> deq(sz(vs) + 5);
    vector<P> ans(sz(vs) + 5);
    deq[0] = vs[0];
    int ah = 0, at = 0, n = sz(vs);
    rep(i,1,n+1) {
        //~ cout << i << " " << ah << " " << deq.size() << endl;
        if (i == n) vs.push_back(deq[ah]);
        if (angDiff(vs[i], vs[i - 1]) == 0) continue;
        while (ah<at && sideOf(sp(vs[i]), ans[at-1], EPS) < 0)
            at--;
        while (i!=n && ah<at && sideOf(sp(vs[i]),ans[ah],EPS)<0)
            ah++;
        auto res = lineInter(sp(vs[i]), sp(deq[at]));
        if (res.first != 1) continue;
        ans[at++] = res.second, deq[at] = vs[i];
    }
    //~ cout << "Heree " << endl;
    if (at - ah <= 2) return {};
    //~ cout << ah << " " << at << " " << sz(ans) << endl;
    return {ans.begin() + ah, ans.begin() + at};
}
template<class T>
T polygonArea2(vector<Point<T>>& v) {
    T a = v.back().cross(v[0]);
    rep(i,0,sz(v)-1) a += v[i].cross(v[i+1]);
    return a;
}

int main() {
    P ed; cin >> ed.x >> ed.y;
    LL n; cin >> n;
    vector <P> poly(n);
    trav(cur, poly) cin >> cur.x >> cur.y;
    vector <LL> order(n);
    vector <Line> halfplane;
    halfplane.pb({P(), P(ed.x, 0)});
    halfplane.pb({P(ed.x, 0), P(ed.x, ed.y)});
    halfplane.pb({P(ed.x, ed.y), P(0, ed.y)});
    halfplane.pb({P(0, ed.y), P()});
    rep(i,0,n){
        cin >> order[i];
        order[i]--;
        rep(j,0,i){
            halfplane.pb({poly[order[i]], poly[order[j]]});
        }
    }
    //~ cout << "hp " << endl;
    //~ trav(cur, halfplane) cout << cur[0].x << " " << cur[0].y  << " to " << cur[1].x << " " << cur[1].y << endl;
    auto res = halfPlaneIntersection(halfplane);
    //~ cout << "OK " << endl;
    cout << fixed << setprecision(10);
    if(!sz(res)) {
        cout << 0 << endl;
        return 0;
    }
    //~ trav(cur, res) cout << cur.x << " " << cur.y << endl;
    auto ans = polygonArea2(res);
    cout << ans / 2 << endl;
    return 0;
    
}
