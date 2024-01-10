#include <bits/stdc++.h>
using namespace std;
typedef long long ll;
int main(void) {
    ll a, b, c, d; cin >> a >> b >> c >> d;
    cout << min((b - a) + (d - c), max(b, d) - min(a, c)) << '\n';
    return 0;
}
