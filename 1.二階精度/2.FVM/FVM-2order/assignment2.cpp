//利用有限體積法求解不可壓縮理想氣體的二維穩態熱擴散對流方程式 
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>
#include <iomanip>
using namespace std;

//會用到的參數
const int NX = 80;
const int NY = 80;
const int n = NX * NY;
double a[n+1][n+1], b[n+1]; //a二維矩陣為等式左側係數矩陣;b一維矩陣為等式右側之外源項，方程組之非齊性項
double x[n+1], x_old[n+1], T[NX][NY];
const double  Pelect = 0.1;//熱對流係數
const double  Gamma = (sqrt(2.0)/Pelect) ;//pelect number
const double dx = 1.0/double(NX) ;
const double dy = 1.0/double(NY) ;
const double a_W = (Gamma + dy/2.0);
const double a_E = (Gamma - dy/2.0);
const double a_S = (Gamma + dx/2.0);
const double a_N = (Gamma - dx/2.0);
int G, max_G = 10000; //限制最大迭代次數
double maxerror; 
const double tolerance = 1e-6;
bool steadystate;

//初始化矩陣
void initial(double a[][n+1], double b[], int n) {
    // 初始化所有元素為0
    for(int i = 1; i <= n; i++) {
        for(int j = 1; j <= n; j++) {
            a[i][j] = 0.0;
        }
        b[i] = 0.0;
        x[i] = 0.0;
        x_old[i] = 0.0;
    }
    
    // 設定邊界條件和係數矩陣
    // 四個角點
    a[1][1] = (6.0*Gamma) + (dx/2.0) + (dy/2.0) ; 
    a[1][2] = -1.0;
    a[1][NX+1] = -1.0;
    b[1] = 2*a_W; // 邊界溫度
    
    a[NX][NX] = (6.0*Gamma) - (dx/2.0) + (dy/2.0);
    a[NX][NX-1] = -1.0;
    a[NX][2*NX] = -1.0;
    b[NX] = 0.0;// 邊界溫度
    
    a[n-NX+1][n-NX+1] = (6.0*Gamma) + (dx/2.0) - (dy/2.0);
    a[n-NX+1][n-NX+2] = -1.0;
    a[n-NX+1][n-2*NX+1] = -1.0;
    b[n-NX+1] = 2*(a_W+a_N); //邊界效應引起之外源項
    
    a[n][n] =(6.0*Gamma) - (dx/2.0) - (dy/2.0);
    a[n][n-1] = -1.0;
    a[n][n-NX] = -1.0;
    b[n] = 2.0*a_N;// 邊界溫度
    
    // 下邊界 (除角點外)
    for(int i = 2; i <= NX-1; i++) {
        a[i][i] = (5*Gamma)+(dx/2.0);
        a[i][i+1] = -1.0;
        a[i][i-1] = -1.0;
        a[i][i+NX] = -1.0;
        b[i] = 0.0; // 邊界溫度
    }
    
    // 上邊界 (除角點外)
    for(int i = n-NX+2; i < n; i++) {
        a[i][i] = (5*Gamma)-(dx/2.0);
        a[i][i+1] = -1.0;
        a[i][i-1] = -1.0;
        a[i][i-NX] = -1.0;
        b[i] = 2.0*a_N;
    }
    
    // 左邊界 (除角點外)
    for(int i = 1; i <= NY-2; i++) {
        int idx = NX*i+1;
        a[idx][idx] = (5*Gamma)+(dy/2.0);
        a[idx][idx+1] = -1.0;
        a[idx][idx+NX] = -1.0;
        a[idx][idx-NX] = -1.0;
        b[idx] = 2*a_W; //邊界溫度
    }
    
    // 右邊界 (除角點外)
    for(int i = 2; i <= NY-1; i++) {
        int idx = NX*i;
        a[idx][idx] = (5*Gamma)-(dy/2.0);
        a[idx][idx-1] = -1.0;
        a[idx][idx-NX] = -1.0;
        a[idx][idx+NX] = -1.0;
        b[idx] = 0.0;
    }
    
    // 內點
    for(int i = 2; i <= NX-1; i++) {
        for(int j = 1; j <= NY-2; j++) {
            int idx = j*NX+i;
            a[idx][idx] = (4*Gamma);
            a[idx][idx+1] = -1.0;
            a[idx][idx-1] = -1.0;
            a[idx][idx+NX] = -1.0;
            a[idx][idx-NX] = -1.0;
            b[idx] = 0.0; // 內點無熱源
        }
    }
}

void Jacobi(double a[][n+1], double b[], double x[], int n) {
    // 先複製當前解到x_old
    for(int k = 1; k <= n; k++) {
        x_old[k] = x[k];
    }
    
    // 計算新的解
    for(int k = 1; k <= n; k++) {
        double sum = 0;
        for(int p = 1; p <= n; p++) {
            if(p != k) {
                sum += a[k][p] * x_old[p];
            }
        }
        x[k] = (b[k] - sum) / a[k][k];
    }
    
    // 計算最大誤差
    maxerror = 0;
    for(int k = 1; k <= n; k++) {
        double error = fabs(x[k] - x_old[k]);
        if(maxerror < error) {
            maxerror = error;
        }
    }
}

void output(int m) {
    // 將一維解轉換為二維溫度場
    for(int j = 0; j < NY; j++) {
        for(int i = 0; i < NX; i++) {
            T[i][j] = x[j*NX + i + 1];
        }
    }
    
    ostringstream name;
    name << "steady_diffusion_2D_" << setfill('0') << setw(6) << m << ".vtk";
    ofstream out(name.str().c_str());
    
    // VTK 文件頭
    out << "# vtk DataFile Version 3.0\n";
    out << "steady_diffusion_2D\n";
    out << "ASCII\n";
    out << "DATASET STRUCTURED_POINTS\n";
    out << "DIMENSIONS " << NX << " " << NY << " 1\n";
    out << "ORIGIN 0 0 0\n";
    out << "SPACING 1 1 1\n";
    out << "POINT_DATA " << NX * NY << "\n";
    
    // 輸出溫度場
    out << "SCALARS Temperature double 1\n";
    out << "LOOKUP_TABLE default\n";
    for(int j = 0; j < NY; j++) {
        for(int i = 0; i < NX; i++) {
            out << scientific << T[i][j] << "\n";
        }
    }
    
    out.close();
    cout << "VTK 文件已輸出: " << name.str() << endl;
}

int main() {
    cout << "程式碼開始執行...." << endl;
    steadystate = false;
    initial(a, b, n);
    
    for(G = 0; G < max_G; G++) {
        if(G > 0) {
            Jacobi(a, b, x, n);
        }
        
        if(G % 1000 == 0) {
            cout << "迭代次數 = " << G << endl;
            if(G > 0) {
                cout << "最大變化量max_change = " << maxerror << endl;
            }
            output(G);
        }
        
        if(G > 100 && maxerror < tolerance) {
            steadystate = true;
            cout << "已經達到迭代終點，溫度場收斂!!" << endl;
            break;
        }
    }
    
    if(!steadystate) {
        cout << "達到最大迭代次數，但未達到穩態!" << endl;
    }
    
    output(G);
    cout << "格子大小 " << NX << "x" << NY << " 計算完成\n" << endl;
    return 0;
}