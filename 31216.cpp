#pragma region template
#pragma clang diagnostic ignored "-Wno-deprecated"
#pragma clang diagnostic ignored "-Werror"
#pragma warning(disable : 4996)
/*
* Countable. Loop known at run time. Must stay constant. Cannot exit by variable. (break then can't use)
* switch statements cannot be used. else if cannot be used
* no function calls
#pragma GCC target ("avx2")
#pragma GCC optimization ("O3")
#pragma GCC optimization ("unroll-loops")
*/
// c#.cpp : This file contains the 'main' function. Program execution begins and ends there.
#include <iomanip>
#include <iostream>
#include <cctype>
#include <string>
#include <math.h> 
#include <cmath>
#include <sstream>
#include <vector>
#include <algorithm>
#include <fstream>
#include <exception>
#include <limits>
#include <new>
#include <typeinfo>
#include <bitset>
#include <deque>
#include <iterator>
#include <numeric>
#include <iosfwd>
#include <ios>
#include <istream>
#include <streambuf>
#include <strstream>
#include <array>
#include <regex>
#include <cassert>
#include <cerrno>
#include <climits>
#include <cfenv>
#include <type_traits>
#include <chrono>
#include <functional>
#include <memory>
#include <set>
#include <map>
#include <stack>
#include <queue>
#include <utility>
#include <cwctype>
#include <cwchar>
#include <unordered_map>
#include <unordered_set>
#include <complex>
#include <bit>
#include <cstdio>
#include <cstring>
#include <stdio.h>
using namespace std; //名前空間の宣言
using std::cout; using std::cin;
using std::endl; using std::string;
using std::vector; using std::istringstream;
using std::ios;
using std::stringstream;
using std::chrono::duration_cast;
using namespace std::chrono;
using ll = long long; using ld = long double;
using ull = unsigned long long; using str = string;
using bl = bool; using ch = char;
using cd = complex<ld>;
typedef vector<ll> vl;
typedef vector<ld> vd;
typedef set<ll> sl;
typedef unordered_set<ll> usl;
typedef unordered_set<str> uss;
typedef vector<vector<ll>> vl2;
typedef vector<str> vs;
typedef vector<ch> vc;
typedef vector<pair<ll, ll>> vp;
typedef vector<bl> vb;
typedef map<ll, str> mls;
typedef map<str, str> mss;
typedef map<ll, ll> mll;
typedef map<str, ll> msl;
typedef map<ch, ll> mcl;
typedef map<ld, ld> mdd;
typedef queue<ll> ql;
typedef deque<ll> dql;
typedef priority_queue<ll> pqlg;
typedef priority_queue<ll, vl, greater<ll>> pqls;
#define _USE_MATH_DEFINES
#define M_PI 3.14159265358979323846
#define LL_MAX 9223372036854775807
#define ScanfNowTryToBeatMyCompileSpeed std::ios::sync_with_stdio(EXIT_SUCCESS); std::cin.tie(EXIT_SUCCESS); std::cout.tie(EXIT_SUCCESS);
#define up(initial, n, step) for(ll i = (ll)(initial);i < (ll)(n);i+=(ll)(step))
#define up2(initial, n, step) for(ll j = (ll)(initial);j < (ll)(n);j+=(ll)(step))
#define up3(initial, n, step) for(ll k = (ll)(initial);k < (ll)(n);k+=(ll)(step))
#define down(initial, n, step) for(ll i = (ll)(initial) - 1;i >= (ll)(n);i-=(ll)(step))
#define all(x) (x).begin(), (x).end()
#define vget(v) for(auto& element : v) element = get();
#define vcin(v) for(auto& element : v) cin >> (element);
#define YES(a) ((a)?"YES":"NO")
#define Yes(a) ((a)?"Yes":"No")
#define yes(a) ((a)?"yes":"no")
#define stop return EXIT_SUCCESS;
#define rev(s) reverse((s).begin(), (s).end());
#define sortCol(v, v2, col) sort(v, v[(col)] < v2[(col)])
#define vsort(v) sort((v).begin(), (v).end())
#define vsum(v) accumulate((v).begin(), (v).end(), 0)
#define vprint(v, spacing) foreach(x, v) cout << (x) << (spacing)
#define toStr(s) to_string((s))
#define throwErr(s) throw invalid_argument(s)
#define setdecimal(n) cout << fixed << setprecision((n))
#define fillAsc(v, start) iota((v).begin(), (v).end(), (start))
#define nextPerm next_permutation
#define prevPerm prev_permutation
#define varConcat(a, b) a##b
#define foreach(a, v) for(auto& a : v)
#define skip continue
#define stc static
#define _m_ int main(void)
#define tc \
  ll testcase = get();  \
  while (testcase--)
#define tcin \
  ll testcase;\
  cin >> testcase;\
  while(testcase--)
auto _ = NULL;//Shame that C++ does not have Discards like C#.
const str alphabet = "abcdefghijklmnopqrstuvwxyz";
const str upAlphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
const ld PHI = (1 + sqrt(5)) / 2;
const ll Mod = 998244353;
const ll Mod2 = pow(10, 9) + 7;
const ld EPS = 1e-9;
const ld PI = 2 * acos(0.0);
str toLower(str s){str result = "";for (auto& c : s)result += tolower(c);return result;}
str toUpper(str s) {str result = "";for (auto& c : s)result += toupper(c);return result;}
ld TAN(ld degree) { return tan(degree * PI / 180.0); }; 
ld SIN(ld degree) { return sin(degree * PI / 180.0); }; 
ld COS(ld degree) { return cos(degree * PI / 180.0); };
ld ATAN(ld len) { return atan(len) * 180 / PI; }
ld ASIN(ld len) { return asin(len) * 180 / PI; }
ld ACOS(ld len) { return acos(len) * 180 / PI; }
ll lexOrder(str s, str s1) { if (s == s1) stop; if (s.length() <= s1.length()) { up(0, s.length(), 1) { if (s[i] < s1[i]) return -1; else if (s[i] > s1[i]) return 1; } return -1; } else { for (int i = 0; i < s1.length(); i++) { if (s[i] > s1[i]) return -1; else if (s[i] < s1[i]) return 1; } return 1; } }
ll revBinarySearch(vector<int> v, int X) { int start = 0, end = v.size() - 1; while (start <= end) { int mid = start + (end - start) / 2; if (X == v[mid]) return mid; else if (X < v[mid]) start = mid + 1; else end = mid - 1; } return -1; }
bl isNumber(const str& st) { for (ch const& c : st) if (isdigit(c) == 0) return false; return true; }
ll getMonth(str m) { m = toLower(m); if (m == "january") return 1; else if (m == "february") return 2; else if (m == "march") return 3; else if (m == "april") return 4; else if (m == "may") return 5; else if (m == "june") return 6; else if (m == "july") return 7; else if (m == "august") return 8; else if (m == "sepember") return 9; else if (m == "october") return 10; else if (m == "november") return 11; return 12; }
bl isLeapYear(int n) { return (n % 4 == 0 ? n % 100 == 0 && n % 400 == 0 ? true : n % 100 != 0 ? true : false : false); }
bl sortcol(const vector<int>& v, const vector<int>& v2, int col) { return v[col] < v2[col]; }
bl cmp(pair<str, ll>& a, pair<str, ll>& b) { if (a.second == b.second) { for (ll i = 0; i < min(a.first.length(), b.first.length()); i++) { if (a.first[i] < b.first[i])return 0; else if (a.first[i] > b.first[i])return 1; } } return a.second < b.second; }
ll gcd(ll a, ll b) { while (b) b ^= a ^= b ^= a %= b; return a;}
ll lcm(ll a, ll b) { return (abs(a * b)/ gcd(a, b)); }
bl islower(str s) { for (const auto& c : s) { if (c > 'z' || c < 'a') return false; } return true; }
bl isupper(str s) { for (const auto& c : s) { if (c > 'Z' || c < 'A') return false; } return true; }
ll ascSum(ll n) { return n * (n + 1) / 2; }
ll fib(ll num) { const ld f = sqrt(5); return (pow(1 + f, num) - pow(1 - f, num)) / (pow(2, num) * f); } //1,1,2,3,5,8,13,21,34,55,...
ll vColSum(vl2 v, ll col) { ll sum = 0; foreach(x, v) sum += x[col]; return sum; }
ll mapMaxElement(mll x) { mll::iterator best = max_element(x.begin(), x.end(), [](const pair<ll, ll>& a, const pair<ll, ll>& b)->bool { return a.second < b.second; }); return best->second; }
ll mapMaxElement(msl x) { msl::iterator best = max_element(x.begin(), x.end(), [](const pair<str, ll>& a, const pair<str, ll>& b)->bool { return a.second < b.second; }); return best->second; }
ll mapMinElement(mll x) { mll::iterator best = min_element(x.begin(), x.end(), [](const pair<ll, ll>& a, const pair<ll, ll>& b)->bool { return a.second < b.second; }); return best->second; }
bl isPrime(ll n) { if (n == 1) return 0; for (ll i = 2; i * i <= n; i++) if (!(n % i))return 0; return 1; }
vl eratosthenesSieve(ll limit) { vl prime(limit + 1); for (ll i = 2; i <= limit; i++)prime[i] = 1; for (ll p = 2; p * p <= limit; p++)if (prime[p])for (ll i = p * p; i <= limit; i += p)prime[i] = 0; return prime; }
vl linearSieve(ll limit) {vl lp(limit + 1), pr; up(2, limit + 1, 1) { if (!lp[i]) lp[i] = i, pr.push_back(i); for (ll j = 0; i * pr[j] <= limit; ++j) { lp[i * pr[j]] = pr[j]; if (pr[j] == lp[i]) break; } } return pr;}
bl isPow2(ll i) { return i && (i & -i) == i; }
bl isPal(str s) { str s1 = s; rev(s1); return s1 == s; }
bl isParenthesis(str s) { ll p = 0; foreach(t, s) { if (t == '(')p++; else if (t == ')') { p--; if (p < 0)  return false; } }return p == 0; }
vl factor(ll n) { vl v; up(1, sqrt(n) + 1, 1) if (!(n % i)) { v.push_back(i); if(n != i * i) v.push_back(n / i); } return v; }
ull fa(ull n) {ull i, fa = 1;for (i = n; i > 1; i--) fa *= i; return fa;}
ull nCr(ull n, ull r) { ull nume = 1, i; for (i = n; i > r; i--) nume *= i; return ull(nume / fa(n - r)); }
ll bigmod(ll a, ll b, ll m) { ll res = 1 % m; while (b) { if (b & 1) res = (res * a) % m; a = (a * a) % m, b >>= 1; }return res; }
auto vmin = [](vl v) {return *min_element(all(v)); };
auto vmax = [](vl v) {return *max_element(all(v)); };
auto multiply = [](ll n, ll n2, ll m = 1) {ll result = 0; while (n2 > 0) { if (n2 & 1)  result += n; n = n << 1; n2 = n2 >> 1; result %= m; } return result; };
ull faM(ull n, ull m) { ull i, fa = 1; for (i = n; i > 1; i--) fa = multiply(fa, i, m); return fa; }
ll lis(vl v) { vl s2(v.size(), 0); ll L = 1; s2[0] = v[0]; up(1, v.size(), 1) { auto it = lower_bound(s2.begin(), s2.begin() + L, v[i]); if (it == s2.begin() + L) s2[L++] = v[i]; else *it = v[i]; }return L; }
vl primeFactors(ll n) { vl res; while (!(n % 2)) res.push_back(2), n /= 2; for (ll i = 3; i * i <= n; i += 2) while (!(n % i)) res.push_back(i), n /= i; if (n > 2) res.push_back(n); return res; }
ll longest_common_subsequence(string X, string Y, ll m, ll n){vl2 L(m, vl(n));for (int i = 0; i <= m; i++) {for (int j = 0; j <= n; j++) {if (i == 0 || j == 0) L[i][j] = 0; else if (X[i - 1] == Y[j - 1]) L[i][j] = L[i - 1][j - 1] + 1; else L[i][j] = max(L[i - 1][j], L[i][j - 1]);}} return L[m][n];}
ll MEX(vl& A) { sl b(A.begin(), A.end()); ll result = 0; while (b.count(result)) ++result; return result; }
//Unity Engine based Vector2
class Vector2 {
public:
    ld x, y;
    friend Vector2 operator-(const Vector2&, const Vector2&);
};
Vector2 operator+(const Vector2& a, const Vector2& b) {
    Vector2 v;
    v.x = a.x + b.x, v.y = a.y + b.y;
    return v;
}
Vector2 operator-(const Vector2& a, const Vector2& b) {
    Vector2 v;
    v.x = a.x - b.x, v.y = a.y + b.y;
    return v;
}
Vector2 operator*(const Vector2& a, const Vector2& b) {
    Vector2 v;
    v.x = a.x * b.x, v.y = a.y * b.y;
    return v;
}
Vector2 operator/(const Vector2& a, const Vector2& b) {
    Vector2 v;
    v.x = a.x / b.x, v.y = a.y / b.y;
    return v;
}
ld Distance(Vector2& a, Vector2& b) {
    return sqrtl(powl(a.x - b.x, 2) + powl(a.y - b.y, 2));
}
ld Magnitude(Vector2* a) {
    return sqrtl(powl(a->x, 2) + powl(a->y, 2));
}
ld ManhattanDistance(Vector2& a, Vector2& b) {
    return fabsl(a.x - b.x) + fabsl(a.y - b.y);
}
bl Collinear(Vector2 a, Vector2 b, Vector2 c) {
    return (b.y - a.y) * (c.x - b.x) == (c.y - b.y) * (b.x - a.x);
}

//fastInput, fastOutput and BigInt provided by GeeksForGeeks.
//fastInput and fastOutput cannot be used interchangeably with std::cin and std::cout
inline ll get(void) { ch t = getchar(); ll x = 0, neg = 0; while ((t < 48 || t>57) && t != '-') t = getchar(); if (t == '-') { neg = 1; t = getchar(); } while (t >= 48 && t <= 57) { x = (x << 3) + (x << 1) + t - 48; t = getchar(); } if (neg) x = -x; return x; }
inline void out(ll x, ll mode = 1) { ch a[20]; a[0] = '0'; ll i = 0, j; if (x < 0) { putchar('-'); x = -x; }if (x == 0) putchar('0'); while (x) { a[i++] = x % 10 + 48; x /= 10; } for (j = i - 1; j >= 0; j--) putchar(a[j]); putchar(mode ? '\n' : ' '); }
//fastPow provided by rookiesLab
ull fastPow(ull b, ull power) {
    ull result = 1;
    while (power > 0) {
        if (power % 2 == 1) result = (result * b);
        b = (b * b);
        power /= 2;
    }
    return result;
}
ull fastPowMod(ull b, ull power, ull mod) {
    b %= mod;
    ull result = 1;
    while (power > 0) {
        if (power % 2 == 1) result = (result * b) % mod;
        b = (b * b) % mod;
        power /= 2;
    }
    return result;
}
bl isInt(ld n) {return n - (ll)n <= EPS;}
template <typename T>
class fenwick {
public:
    vector<T> fenw;
    int n;

    fenwick(int _n) : n(_n) {
        fenw.resize(n);
    }

    void modify(int x, T v) {
        while (x < n) {
            fenw[x] += v;
            x |= (x + 1);
        }
    }

    T get(int x) {
        T v{};
        while (x >= 0) {
            v += fenw[x];
            x = (x & (x + 1)) - 1;
        }
        return v;
    }
};
struct FenwickTree2D {
    vector<vector<int>> bit;
    vector<int>g;
    int n, m;
    void init(int n, int m) {
        this->n = n;
        this->m = m;
        for (int i = 0; i < n; i++)
        {
            g.assign(m, 0);
            bit.push_back(g);
        }
    }
    FenwickTree2D(int n, int m, vector<vector<int>>a) {
        init(n, m);
        for (size_t i = 0; i < a.size(); i++)
            for (size_t j = 0; j < a[i].size(); j++)
                add(i, j, a[i][j]);
    }
    int sum(int x, int y) {
        int ret = 0;
        for (int i = x; i >= 0; i = (i & (i + 1)) - 1)
            for (int j = y; j >= 0; j = (j & (j + 1)) - 1)
                ret += bit[i][j];
        return ret;
    }

    void add(int x, int y, int delta) {
        for (int i = x; i < n; i = i | (i + 1))
            for (int j = y; j < m; j = j | (j + 1))
                bit[i][j] += delta;
    }
};
ll reverse(ll num, ll lg_n) {
    ll res = 0;
    up(0, lg_n, 1) if (num & (1 << i)) res |= 1 << (lg_n - 1 - i);
    return res;
}
void FFT(vector<cd>& a, bl invert) {
    ll n = a.size(), lg_n = 0;
    while ((1 << lg_n) < n) lg_n++;
    for (ll i = 0; i < n; i++) if (i < reverse(i, lg_n)) swap(a[i], a[reverse(i, lg_n)]);
    for (ll len = 2; len <= n; len <<= 1) {
        ld ang = 2 * PI / len * (invert ? -1 : 1);
        cd wlen(cos(ang), sin(ang));
        for (int i = 0; i < n; i += len) {
            cd w(1);
            for (ll j = 0; j < len / 2; j++) {
                cd u = a[i + j], v = a[i + j + len / 2] * w;
                a[i + j] = u + v;
                a[i + j + len / 2] = u - v;
                w *= wlen;
            }
        }
    }
    if (invert) for (cd& x : a) x /= n;
}
vector<int> multiplyFFT(vector<int> const& a, vector<int> const& b) {
    vector<cd> fa(a.begin(), a.end()), fb(b.begin(), b.end());
    int n = 1;
    while (n < a.size() + b.size())
        n <<= 1;
    fa.resize(n);
    fb.resize(n);

    FFT(fa, false);
    FFT(fb, false);
    for (int i = 0; i < n; i++)
        fa[i] *= fb[i];
    FFT(fa, true);

    vector<int> result(n);
    for (int i = 0; i < n; i++)
        result[i] = round(fa[i].real());
    int carry = 0;
    for (int i = 0; i < n; i++){
        result[i] += carry;
        carry = result[i] / 10;
        result[i] %= 10;
    }
    return result;
}
void setIO(string s) {
    freopen((s + ".in").c_str(), "r", stdin);
    freopen((s + ".out").c_str(), "w", stdout);
}
#pragma endregion
void p1(){
    /*
Problem Statement (Draft):
Given an array of n integers,  a[1], a[2] ... a[n],
answer q queries.
For each query, multiply all integers between the Lth element to the Rth element inclusively, and print the parity of the result.

Solution:
Prefix Sum. Even = 1, Odd = 0. If prefix sum >= 1 then there exist even between the integers, thus output even.

Problem Statement 2(Draft):
Given an array of n integers,  a[1], a[2] ... a[n],
answer q queries.
For each query, XOR all integers between the Lth element to the Rth element inclusively, and print the parity of the result.

Solution:
2 Prefix Sum for even and odd counts separately. If prefix sum of odd count mod 2 = 0, output even.
*/
    ll n = get();
vl v(n); vget(v);
vl prefixSum(n + 1, 0);
up(1, n + 1, 1) prefixSum[i] += prefixSum[i - 1] + (v[i - 1] + 1) % 2;
vprint(prefixSum, ' ');

tc{
    ll l = get(), r = get();
cout << (prefixSum[r] - prefixSum[l - 1] ? "even" : "odd") << '\n';
}
}


ll LP_binary_search(ld slope, ld yintercept, ld x, ld y) {
    /*
    -1: y < mx + c (Below)
    0: y = mx + c
    1: y > mx + c (Above)
    */
    ld p = slope * x + yintercept - y;
    if (EPS < p) {
        return -1;
    }
    if (-EPS > p) {
        return 1;
    }
    return 0;
}
void p2(){
    /*
Problem Statement (Draft):
Given a n-sided closed convex polygon, and a point P(x, y) lying inside the polygon, find a pair of vertices such that these 3 points are collinear.
No other pairs exists.
n >= 4.
All points' x-coordinate and y-coordinate are integers.
Vertices are given in clockwise order (Starting from an arbitary vertex)
For simplicity, no 2 points with the same y-coordinate is given.

Please output the answer in the order of input (Clockwise).
    */
    ll n = get();
vl2 v(n);
up(0, n, 1) {
    ll x = get(), y = get();
    v[i] = { x, y };
}
ld p1 = get(), p2 = get();
up(0, n, 1) {
    ld x1 = v[i][0], y1 = v[i][1];
    ll L = 0, R = n - 1;
    while (L + 1 < R) {
        ll l = (L + i), r = (R + i);
        ll m = (l + r) / 2 % n;
        ld x2 = v[m][0], y2 = v[m][1], slope = (y2-y1)/(x2-x1);
        /*
        y - y1 = m(x - x1)
        y = mx - m * x1 + y1
        */
        ll res = LP_binary_search(slope, -slope * x1 + y1, p1, p2);
        if (res == 0) {
            cout << x1 << ' ' << y1 << '\n' << x2 << ' ' << y2 << '\n';
            return; //return 0.
        }
        else if (res == -1) {
            m = (L + R) / 2;
            R = m;
        }
        else {
            m = (L + R) / 2;
            L = m;
        }
        l = (L + i), r = (R + i);
        m = (l + r) / 2 % n;
        if (m != i) {
            x2 = v[m][0], y2 = v[m][1], slope = (y2 - y1) / (x2 - x1);
            res = LP_binary_search(slope, -slope * x1 + y1, p1, p2);
            if (res == 0) {
                cout << x1 << ' ' << y1 << '\n' << x2 << ' ' << y2 << '\n';
                return; //change to return 0.
            }
        }
    }
}
cout << -1 << '\n';
}
void p3(){
    /*
Theorem:
Let f(x) be the number of factors of x.
f(n^2) % f(n) == 0 if k is integer, then it must be 3, and n must be a perfect square.

Solution:
Let numbers of prime be p.
Iterate through all primes (denote by pr). For each iterate through all primes again, (pr * pq^2)^2 is a solution. O(p^2)

Problem Statement (Draft):
Denote f(x) be the number of factors of x.
A number n is called Marvelous if f(n^2) mod f(n) == 0 and 1 <= f(n) <= .

q queries
l r
Find the number of Marvelous integers between l and r inclusively.

1 <= l <= r <= 2^64
1 <= q <= 10^7

144 324 400 784 1936 2025
12 = 2 * 2 * 3
18 = 2 * 3 * 3 * 5

*/
vl sol;
vl pr = linearSieve(1100);
vsort(pr);
foreach(x, pr) foreach(y, pr) if (x != y) sol.push_back(x * x * x * x * y * y);
vsort(sol);
    
ll q = get();
while (q--) {
    ll l = get(), r = get();
    auto u = lower_bound(all(sol), r), d = lower_bound(all(sol), l);
    cout << distance(d, u) + 1 << '\n';
}
}
_m_{ //main函數的開始
ScanfNowTryToBeatMyCompileSpeed
    vl v = linearSieve(10000000);
    vsort(v);
    vl res;
    for(int i = 2;i <= v.size();i++) if(isPrime(i)) res.push_back(v[i-1]);
    tc{
    ll n = get();
    cout << res[n - 1] << '\n';
}
stop  //編碼程序結束
}
// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu
// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file   #pragma region template
