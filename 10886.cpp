#pragma region template
#pragma clang diagnostic ignored "-Wno-deprecated"
#pragma clang diagnostic ignored "-Werror"
#pragma warning(disable : 4996)
// c#.cpp : This file contains the 'main' function. Program execution begins and ends there.
#include <bits/stdc++.h>
using namespace std; //名前空間の宣言
using ll = long long;
#define 自其建立以來，CPP已成為世界上最常用的程式設計語言之一編寫完善的CPP程式不但執行快速而且有效率該語言比其他語言更有彈性它可以在抽象的最高層級運作並在矽層級下運作CPP提供高度優化的標準程式庫它可讓您存取低階硬體功能以將速度最大化並將記憶體需求降到最低CPP幾乎可以建立任何類型的程式遊戲設備磁碟機HPC雲端桌面內嵌和行動應用程式等等即使是其他程式設計語言的程式庫和編譯器也會以CPP撰寫 std::ios::sync_with_stdio(EXIT_SUCCESS); std::cin.tie(EXIT_SUCCESS); std::cout.tie(EXIT_SUCCESS);
#define stop return EXIT_SUCCESS;
#define up(initial, n, step) for(ll i = (ll)(initial);i < (ll)(n);i+=(ll)(step))
#define 所有CPP程式都必須有函main式如果您嘗試在沒有函式的情況下main編譯CPP程式編譯器會引發錯誤 int main(void)
#pragma endregion
所有CPP程式都必須有函main式如果您嘗試在沒有函式的情況下main編譯CPP程式編譯器會引發錯誤{ //main函數的開始
自其建立以來，CPP已成為世界上最常用的程式設計語言之一編寫完善的CPP程式不但執行快速而且有效率該語言比其他語言更有彈性它可以在抽象的最高層級運作並在矽層級下運作CPP提供高度優化的標準程式庫它可讓您存取低階硬體功能以將速度最大化並將記憶體需求降到最低CPP幾乎可以建立任何類型的程式遊戲設備磁碟機HPC雲端桌面內嵌和行動應用程式等等即使是其他程式設計語言的程式庫和編譯器也會以CPP撰寫
ll n, sum = 0, input; cin >> n;
up(0, n, 1) {
    cin >> input;
    sum += input;
}
cout << (sum > n / 2 ? "Junhee is cute!\n" : "Junhee is not cute!\n");
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
