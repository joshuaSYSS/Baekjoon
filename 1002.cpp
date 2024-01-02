#include <iostream>
#include <cmath>

using namespace std;

int main() {
    int T;
    cin >> T;

    while (T--) {
        int x1, y1, r1, x2, y2, r2;
        cin >> x1 >> y1 >> r1 >> x2 >> y2 >> r2;

        // 두 점 사이의 거리 계산
        double d = sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));

        // 무한대인 경우
        if (d == 0 && r1 == r2) {
            cout << -1 << endl;
        }
        // 두 원이 외접하거나 내접하는 경우
        else if (d == r1 + r2 || d == abs(r1 - r2)) {
            cout << 1 << endl;
        }
        // 두 원이 만날 수 있는 경우
        else if (d < r1 + r2 && d > abs(r1 - r2)) {
            cout << 2 << endl;
        }
        // 두 원이 만나지 않는 경우
        else {
            cout << 0 << endl;
        }
    }

    return 0;
}
