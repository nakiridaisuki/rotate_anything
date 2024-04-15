#include <bits/stdc++.h>
#include <windows.h>

using namespace std;

#define max(a, b) (a > b ? a : b)
#define min(a, b) (a < b ? a : b)

// the coordinate of every point or the vector of the normal
struct point{
    double x, y, z;
};

// the 3x3 rotation matrix
struct matrix{
    double a, b, c, d, e, f, g, h, i;
}M;

// the parameters of this heart
int H = 30, W = 80; // length and width of the canvas
double times = 10; // how large of your heart
double A = 0, B = 0, C = 0; // the initialize angle of the rotation matrix
double dA = 0, dB = 0.12, dC = 0.06; // the change of angle in every loop
int delay = 100; // the time sleep of every loop
point light1 = {1, 1, 1}, light2 = {-1, -1, 1};

// Newton's method to solve the root of y_axis coordinate
// there are two root of every x, but I can't solve it in every x
// still fixing

// the function
double f(double x, double y){
    return pow(y,6) + (3*pow(x,2)-3)*pow(y,4) - pow(x,2)*pow(y,3) + (3*pow(x,4)-6*pow(x,2)+3)*pow(y,2) + pow(x,6)-3*pow(x,4)+3*pow(x,2)-1;
}

// the differential of the function
double df(double x, double y){
    return 6*pow(y,5) + 4*(3*pow(x,2)-3)*pow(y,3) - 3*pow(x,2)*pow(y,2) + 2*(3*pow(x,4)-6*pow(x,2)+3)*y;
}

// to get the root from the function with Newton's method
double get_y(double x, double y){
    double now;
    do{
        now = y;
        y = y - f(x, y)/df(x, y);
    }while(abs(now - y) > 1e-3);
    return y;
}

// the matrix multiplication function
// here is how I rotate the heart

// for every point
point point_rotate(double x, double y, double z){
    point p;

    p.x = W/2 + times * (M.a*x + M.b*y + M.c*z)*2;
    p.y = H/2 + times * (M.d*x + M.e*y + M.f*z);
    p.z = 100 + times * (M.g*x + M.h*y + M.i*z);
    return p;
}

// for every vector
point vector_rotate(point p){
    double x = p.x, y = p.y, z = p.z;
    p.x = (M.a*x + M.b*y + M.c*z)*2;
    p.y = (M.d*x + M.e*y + M.f*z);
    p.z = (M.g*x + M.h*y + M.i*z);
    return p;
}

// the dot and cross product of vectors for calculation of lumination
// I use cross product to get the surface normal
// and use dot product to get the cosine of the angle between the light direction and the surface normal
// the higher the value, the more light falls on the surface

point cross(point u, point v){
    return {u.y*v.z-u.z*v.y, u.z*v.x-u.x*v.z, u.x*v.y-u.y*v.x};
}

double dot(point u, point v){
    return u.x*v.x + u.y*v.y + u.z*v.z;
}

// the Fast Inverse Square Root 
float Q_Inv_sqrt(float x){
    float x2 = 0.5f * x;
    int i = *(int*)&x;       // evil floating point bit level hacking
    i = 0x5f3759df - (i>>1); // what the fuck?
    x = *(float*)&i;
    x = x * (1.5f - x2*x*x); // 1st iteration
    return x;
}

// get the length of the vector
double value(point p){
    return Q_Inv_sqrt(p.x*p.x + p.y*p.y + p.z*p.z);
}

// precanculate the length of light vector
double light1_val = value(light1);
double light2_val = value(light2);

// calculate cosine 
int get_luminance(point v){
    double L1 = max(0, dot(v, light1) * value(v) * light1_val) * 12;
    double L2 = max(0, dot(v, light2) * value(v) * light2_val) * 12;
    return int(max(L1, L2));
}

// chack if the coordinates is in the canvas
bool chack_x(int x){
    return 0 <= x && x < W;
}
bool chack_y(int y){
    return 0 <= y && y < H;
}

int main() {
    ios::sync_with_stdio(0);
    cin.tie(0);

    char canvas[H][W];
    double deepth[H][W];
    printf("\033[?25l"); // hide the cursor
    printf("\x1b[2J");   // clear the terminal

    LARGE_INTEGER t1,t2,tc;
    QueryPerformanceFrequency(&tc);

    int loop = 0;
    while(1){
        QueryPerformanceCounter(&t1);

        // initialize the canvas and deepth
        memset(canvas, ' ', sizeof(canvas));
        memset(deepth, 0, sizeof(deepth));

        // precalculate sine and cosine
        double sA = sin(A), sB = sin(B), sC = sin(C);
        double cA = cos(A), cB = cos(B), cC = cos(C);

        // calculate the rotation matrix
        M.a = cC*cB;
        M.b = -sC*cA + cC*sB*sA;
        M.c = sC*sA + cC*sB*cA;
        M.d = sC*cB;
        M.e = cC*cA + sC*sB*sA;
        M.f = -cC*sA + sC*sB*cA;
        M.g = -sB;
        M.h = cB*sA;
        M.i = cB*cA;

        // creat the heart
        vector<double> px(90), py(90), pz(90); // record the coordinate to calculate the surface normal
        for(double i=0; i<1.139; i+=0.001){
            double y1 = get_y(i, 1.5);
            double y2 = get_y(i, -1.5);
            double y_mean = (y1 + y2) / 2, r = abs(y_mean - y1);
            if(abs(r) < 1e-3) continue;

            int idx = 0;
            for(double j=0; j<6.28; j+=0.07, idx++){
                double sj = sin(j), cj = cos(j);

                // the coordinate
                double x, y, z;
                x = i;
                y = y_mean + r*cj;
                z = r*sj*0.8;

                // get the right half surface normal
                point u = {0, -r*sj, r*cj*0.8}, v, w;
                if(i == 0) v = {1, 0, 0};
                else v = {x-px[idx], y-py[idx], z-pz[idx]};
                w = cross(u, v);

                // get the right half heart
                point p1 = point_rotate(x, y, z), n1 = vector_rotate(w);
                int L = get_luminance(n1);
                if(chack_x(int(p1.x)) && chack_y(int(p1.y)) && deepth[int(p1.y)][int(p1.x)] < p1.z){
                    deepth[int(p1.y)][int(p1.x)] = p1.z;
                    canvas[int(p1.y)][int(p1.x)] = ".,-~:;=!*#$@"[L];
                }

                // get the left half surface normal
                v.x = -v.x;
                w = cross(v, u);

                // get the left half heart
                point p2 = point_rotate(-x, y, z), n2 = vector_rotate(w);
                L = get_luminance(n2);
                if(chack_x(int(p2.x)) && chack_y(int(p2.y)) && deepth[int(p2.y)][int(p2.x)] < p2.z){
                    deepth[int(p2.y)][int(p2.x)] = p2.z;
                    canvas[int(p2.y)][int(p2.x)] = ".,-~:;=!*#$@"[L];
                }

                // record the coordinate
                px[idx] = x;
                py[idx] = y;
                pz[idx] = z;
            }

            i += r * 0.003;
        }

        printf("\033[?25l"); // hide the cursor
        printf("\x1b[H");    // move the cursor to home position
        printf("\x1b[;3%dm", loop%6+1);
        loop = (loop+1) % 6;
        // print the heart
        for(int i=H-1; i>=0; i--){
            for(int j=0; j<W; j++){
                // char tmp = canvas[i][j];
                // if(tmp == '.' || tmp == ',')
                //     printf("\x1b[;31m%c\x1b[0;m", canvas[i][j]);
                // else if(tmp == '-' || tmp == '~')
                //     printf("\x1b[;36m%c\x1b[0;m", canvas[i][j]);
                // else if(tmp == ';' || tmp == ':')
                //     printf("\x1b[;34m%c\x1b[0;m", canvas[i][j]);
                // else if(tmp == '=' || tmp == '!')
                //     printf("\x1b[;36m%c\x1b[0;m", canvas[i][j]);
                // else if(tmp == '*' || tmp == '#')
                //     printf("\x1b[;32m%c\x1b[0;m", canvas[i][j]);
                // else if(tmp == '$' || tmp == '@')
                //     printf("\x1b[;33m%c\x1b[0;m", canvas[i][j]);
                // else
                //     printf("%c", tmp);
                printf("%c", canvas[i][j]);
            }
            printf("\n");
        }
        A += dA;
        B += dB;
        C += dC;
        printf("%f %f %f\n", A, B, C);

        QueryPerformanceCounter(&t2);
        int time = (double)(t2.QuadPart - t1.QuadPart) / (double)tc.QuadPart * 1000;
        Sleep(max(0, delay - time));
    }

    return 0;
}