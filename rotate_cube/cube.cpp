#include <bits/stdc++.h>
#include <windows.h>

using namespace std;

#define ll long long
#define N 45
#define N_2 (1 << 18)
#define pi 3.14159265
#define oo 0x3f3f3f3f
#define ool 0x3f3f3f3f3f3f3f3f
#define P (ll)(998244353)
#define F first
#define S second
#define pii pair<int, int>
#define pll pair<ll, ll>
#define lowbit(x) (x&-x)
// #define max(a, b) (a > b ? a : b)
// #define min(a, b) (a < b ? a : b)
#define lc (id << 1)
#define rc (id << 1 | 1)
// #define m ((l + r) / 2)

int H = 25, W = 80;
double a = 6;
double A = 0, B = 0, C = 0;
double dA = 0.09, dB = 0.1, dC = 0.1;
int delay = 10;

struct point{
    int x, y;
    double z;
};

struct matrix{
    double a, b, c, d, e, f, g, h, i;
}M;

point tran(double x, double y, double z){
    point p;

    p.x = W/2 + (M.a*x + M.b*y + M.c*z)*2;
    p.y = H/2 + M.d*x + M.e*y + M.f*z;
    p.z = 100 + M.g*x + M.h*y + M.i*z;
    return p;
}

int main() {
    ios::sync_with_stdio(0);
    cin.tie(0);

    char b[H][W];
    double d[H][W];
    int type[2] = {-1, 1};
    printf("\033[?25l");
    printf("\x1b[2J");

    while(1){
        memset(b, ' ', sizeof(b));
        memset(d, 0, sizeof(d));

        M.a = cos(C)*cos(B);
        M.b = -sin(C)*cos(A) + cos(C)*sin(B)*sin(A);
        M.c = sin(C)*sin(A) + cos(C)*sin(B)*cos(A);
        M.d = sin(C)*cos(B);
        M.e = cos(C)*cos(A) + sin(C)*sin(B)*sin(A);
        M.f = -cos(C)*sin(A) + sin(C)*sin(B)*cos(A);
        M.g = -sin(B);
        M.h = cos(B)*sin(A);
        M.i = cos(B)*cos(A);

        // top, bottom, left, right, front, back
        for(int t=0; t<6; t++){
            for(double i=-a; i<a; i+=0.1){
                for(double j=-a; j<a; j+=0.1){
                    double x, y, z;
                    if(t == 0 || t == 1)
                        x=i, y=a*type[t&1], z=j;
                    else if(t == 2 || t == 3)
                        x=a*type[t&1], y=i, z=j;
                    else
                        x=i, y=j, z=a*type[t&1];

                    point p = tran(x, y, z);
                    if(0 < p.x && p.x < W && 0 <= p.y && p.y < H && p.z > d[p.y][p.x]){
                        d[p.y][p.x] = p.z;
                        b[p.y][p.x] = ":;-~*^"[t];
                    }
                }
            }
        }

        printf("\033[?25l");
        printf("\x1b[H");
        for(int i=0; i<H; i++){
            for(int j=0; j<W; j++)
                printf("%c", b[i][j]);
            printf("\n");
        }
        A += dA;
        B += dB;
        C += dC;
        printf("%f %f %f", A, B, C);

        Sleep(delay);
    }

    return 0;
}