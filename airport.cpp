/*
Author : Hocky Yudhiono
Min 08 Mar 2020 03:55:47  WIB

1. You can sort the query if offline!
2. Don't bring the dp remaining state when dfsing on DP on Tree.
3. Try to reverse (Think from the back) if you stuck.
4. Be careful when submitting Div. 2 D-F, dont waste it on stupid WAs.
5. Try to reduce a problem
*/

#include <algorithm>
#include <iostream>
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

typedef long long LL;
typedef long double LD;
typedef unsigned long long ULL;
typedef pair<int,int> PII;
typedef pair<LL,LL> PLL;
//If the time limit is strict, try not to use long double

#define popf pop_front
#define pf push_front
#define popb pop_back
#define mp make_pair
#define pb push_back
#define remove erase
#define fi first
#define se second

//Remember to undefine if the problem is interactive
#define endl '\n'
#define DEBUG(X) cout << ">>> DEBUG(" << __LINE__ << ") " << #X << " = " << (X) << endl

const double EPS = 1e-9;
const int INFMEM = 63;
const int INF = 1061109567;
const LL LINF = 4557430888798830399LL;
const double DINF = numeric_limits<double>::infinity();
const LL MOD = 1000000007;
const int dx[8] = {1,0,-1,0,1,1,-1,-1};
const int dy[8] = {0,1,0,-1,1,-1,1,-1};
const double PI = 3.141592653589793;

inline void open(string a){
    freopen((a+".in").c_str(),"r",stdin);
    freopen((a+".out").c_str(),"w",stdout);
}

inline void fasterios(){
    //Do not use if interactive
    ios_base::sync_with_stdio(0);
    cin.tie(0); cout.tie(0);
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

typedef Point<LD> P;
typedef pair<P,pair<int,int> > intersectPoint;

template<class P>
pair<int, P> lineInter(P s1, P e1, P s2, P e2) {
    auto d = (e1 - s1).cross(e2 - s2);
    if (d == 0) // if parallel
        return {-(s1.cross(e1, s2) == 0), P(0, 0)};
    auto p = s2.cross(e1, e2), q = s2.cross(e2, s1);
    return {1, (s1 * p + e1 * q) / d};
}

LL n;
P isi[205];
LD segDist(P& s, P& e, P& p) {
    if (s==e) return (p-s).dist();
    auto d = (e-s).dist2(), t = min(d,max(LD(.0),(p-s).dot(e-s)));
    return ((p-s)*d-(e-s)*t).dist()/d;
}
template<class P> bool onSegment(P s, P e, P p) {
    return eq(segDist(s,e,p),(LD)(0));
}

int main(){
    fasterios();
    cin >> n;
    for(int i = 1;i <= n;i++){
        cin >> isi[i].x >> isi[i].y;
    }
    LD ans = 0;
    for(int i = 1;i <= n;i++){
        for(int j = i+1;j <= n;j++){
            vector<intersectPoint> all;
            for(int k = 1;k <= n;k++){
                auto res = lineInter(isi[i],isi[j],isi[k],isi[k+1 <= n ? k+1 : 1]);
                // if(res.fi == 1) cout << res.se << endl;
                if(res.fi == 1 && onSegment(isi[k],isi[k+1 <= n ? k+1 : 1],res.se)){
                    all.pb({res.se,{isi[i].ccw(isi[j],isi[k]),isi[i].ccw(isi[j],isi[k+1 <= n ? k+1 : 1])}});
                }
            }
            // cout << "Currently computing " << isi[i] << " " << isi[j] << endl;
            sort(all.begin(),all.end());
            // for(auto poi : all){
            //     cout << poi.fi << " " << poi.se.fi << " " << poi.se.se << endl;
            // }
            int isInside = -1,pre;
            LD curAns = 0; P lst;
            for(int k = 0;k < all.size();k++){
                if(isInside >= 0) curAns += (all[k].fi-all[k-1].fi).dist();
                if(all[k].se.fi*all[k].se.se == -1) isInside = -isInside;
                else if(isInside == 0) isInside = (all[k].se.fi+all[k].se.se)*pre;
                else{
                    pre = (all[k].se.fi+all[k].se.se)*isInside;
                    isInside = 0;
                }
                // cout << curAns << endl;
                ans = doubleMax(ans,curAns);
                if(isInside == -1) curAns = 0,lst = all[k+1 <= (int)(all.size())-1 ? k+1 : 0].fi;
            }
        }
    }
    cout << fixed << setprecision(10) << ans << endl;
    return 0;
}