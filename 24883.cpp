#pragma region template
#pragma clang diagnostic ignored "-Wno-deprecated"
#pragma clang diagnostic ignored "-Werror"
// cpp.cpp : This file contains the 'main' function. Program execution begins and ends there.
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
#include <memory.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <set>
#include <map>
#include <stack>
#include <queue>
#include <utility>
#include <cwctype>
#include <cwchar>
#include <unordered_set>
using namespace std;
using std::cout; using std::cin;
using std::endl; using std::string;
using std::vector; using std::istringstream;
using std::ios;
using std::stringstream;
using ll = long long; using ld = long double;
using ull = unsigned long long; using str = string;
using bl = bool; using db = double;
using ch = char; using sh = short;
typedef vector<ll> vl;
typedef vector<ld> vd;
typedef set<ll> sl;
typedef unordered_set<ll> usl;
typedef vector<vector<ll>> vl2;
typedef vector<str> vs;
typedef vector<ch> vc;
typedef map<ll, str> mls;
typedef map<str, str> mss;
typedef map<ll, ll> mll;
typedef map<str, ll> msl;
typedef map<ch, ll> mcl;
typedef stack<ll> stl;
typedef queue<ll> ql;
typedef deque<ll> dq;
typedef priority_queue<ll> pq;
#define _USE_MATH_DEFINES
#define M_PI 3.14159265358979323846
#define INT_MAX 2147483647
#define LLONG_MAX 9223372036854775807
#define ULLONG_MAX 18446744073709551615
#define LL_MAX LLONG_MAX
#define up(initial, n, step) for(ld i = (initial);i < n;i+=(step))
#define down(initial, n, step) for(ld o = (initial) - 1;o >= n;o-=(step))
#define forCond(initial, cond, step) for(ld i = (initial);cond;i+=(step))
#define floor(a) (ll)a
#define vget(v, n) for(ll i = 0;i < n;i++) cin>>(v[i]);
#define YES(a) ((a)?"YES":"NO")
#define Yes(a) ((a)?"Yes":"No")
#define yes(a) ((a)?"yes":"no")
#define stop return 0;
#define rev(s) reverse((s).begin(), (s).end());
#define sortCol(v, v2, col) sort(v, v[(col)] < v2[(col)])
#define vsort(v) sort(v.begin(), v.end())
#define pb(a) push_back((a));
#define vpushf(v, a) (v).insert((v).begin(), (a))
#define vsum(v) accumulate(v.begin(), v.end(), 0)
#define vavg(v) accumulate(v.begin(), v.end(), 0) / v.size()
#define vremoveDupe(v) (v).erase(unique((v).begin(), (v).end()), (v).end());
#define vrand(v) (v).random_shuffle(v.begin(), v.end());
#define vfind(v, val) find((v).begin(), v.end(), val)
#define spush(s, a) (s).insert((a));
#define stTop(s) (s).top();
#define toStr(s) to_string((s))
#define throwErr(s) throw invalid_argument(s);
#define setdecimal(n) cout << fixed << setprecision((n))
#define fillAsc(v, start) iota((v).begin(), (v).end(), (start));
#define all(x) (x).begin(), (x).end()
#define isPerm is_permutation
#define nextPerm next_permutation
#define prevPerm prev_permutation
#define testcase \
  ll t = get();  \
  while (t--)
#define varConcat(a, b) a##b
#define setup ios::sync_with_stdio(false); cin.tie(0); cout.tie(0);
#pragma endregion
int main(void) {
    setup;
    ch c; cin >> c; c = tolower(c);
    cout << (c == 'n' ? "Naver D2\n" : "Naver Whale\n");
    stop;
}
// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu
// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
