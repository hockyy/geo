/*
Author : Hocky Yudhiono
Sel 03 Okt 2023 07:59:08
*/

#include <bits/stdc++.h>
using namespace std;
typedef long long LL;
#define sz(a) (int) a.size()
#define rep(i,a,b) for(int i = a;i < b;i++)
#define pb push_back

template <class T>
struct Point {
  typedef Point P;
  T x, y;
  explicit Point(T _x = 0, T _y = 0): x(_x), y(_y) {
  }
  bool operator<(const Point &other) const {
    return tie(x, y) < tie(other.x, other.y);
  }
  P operator-(const Point other)const {
    return P(x - other.x, y - other.y);
  }
  P operator+(const Point other)const {
    return P(x + other.x, y + other.y);
  }
  P operator*(const T other) const {
    return P(x * other, y * other);
  }
  P operator/(const T other) const {
    return P(x / other, y / other);
  }
  T cross(const Point other) const {
    return x * other.y - y * other.x;
  }
  T cross(const Point b, const Point c) const {
    return (b - *this).cross(c - *this);
  }
};

typedef Point<LL> PL;

struct Segment {
  PL v;
  int a, b;
  bool operator<(const Segment &other)const {
    return v.cross(other.v) > 0;
  }
};
LL target;
LL getArea2(const PL &a, const PL &b, const PL &c) {
  LL curArea = abs(a.cross(b) + b.cross(c) + c.cross(a));
  // cout << "Here " << curArea << endl;
  if (curArea == target) {
    cout << "Yes" << endl;
    cout << a.x << " " << a.y << endl;
    cout << b.x << " " << b.y << endl;
    cout << c.x << " " << c.y << endl;
    exit(0);
  }
  return curArea;
}

int main() {
  LL n; cin >> n >> target;
  target *= 2;
  vector <PL> isi(n);
  for(auto &cur : isi) cin >> cur.x >> cur.y;
  sort(begin(isi), end(isi));
  vector <Segment> events;
  vector <int> pointOrder(n), positionOfIndex(n);
  rep(i,0,n) rep(j,i + 1, n) events.pb({isi[j] - isi[i], i, j});
  rep(i,0,n) pointOrder[i] = positionOfIndex[i] = i;
  sort(begin(events), end(events));

  rep(i,0,sz(events)) {
    int &a = positionOfIndex[events[i].a];
    int &b = positionOfIndex[events[i].b];
    assert(a < b);
    {
      int kiri = 0, kanan = min(a, b) - 1;
      while (kiri <= kanan) {
        LL mid = (kiri + kanan) >> 1;
        LL curArea = getArea2(isi[events[i].a], isi[events[i].b], isi[pointOrder[mid]]);
        if (curArea < target) kanan = mid - 1;
        else kiri = mid + 1;
      }
    }

    {
      int kiri = max(a, b) + 1, kanan = n - 1;
      while (kiri <= kanan) {
        LL mid = (kiri + kanan) >> 1;
        LL curArea = getArea2(isi[events[i].a], isi[events[i].b], isi[pointOrder[mid]]);
        if (curArea < target) kiri = mid + 1;
        else kanan = mid - 1;
      }
    }
    swap(pointOrder[a], pointOrder[b]);
    swap(a, b);
  }
  cout << "No" << endl;
}