#include <bits/stdc++.h>
using namespace std;
typedef long long ll;
typedef vector<ll> vl;
int main(void) {
ll n, m, start = 0, res = 0; cin >> n >> m;
	vector<int> limit(101);
	for (int i = 0; i < n; i++) {
		ll l, s; cin >> l >> s;
		for (int j = start; j < start + l; j++) limit[j] = s;
		start += l;
	}
	start = 0;
	vl bessie(101);
	for (int i = 0; i < m; i++) {
		ll l, s; cin >> l >> s;
		for (int j = start; j < start + l; j++) bessie[j] = s;
		start += l;
	}
	for (int i = 0; i < 100; i++) res = max(res, bessie[i] - limit[i]);
	cout << res << '\n';
    return 0;
}
