//採用標準四階九點方程求解Laplace方程 
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>
#include <iomanip>
#define pi 3.14159265358979323846

using namespace std;

int nx_data[] = {11, 21, 41, 81};
int ny_data[] = {11, 21, 41, 81};

vector<int> Nx(nx_data, nx_data + sizeof(nx_data)/sizeof(nx_data[0]));
vector<int> Ny(ny_data, ny_data + sizeof(ny_data)/sizeof(ny_data[0]));
double FinalL1error[4];
int n, NX, NY;
double dx, dy; 
vector<vector<double> > a;
vector<double> b;
vector<double> x, x_old;
vector<vector<double> > T;
int G, max_G = 100000;
double T_left = 10.0;
double T_right = 10.0;
double T_bottom = 10.0;
double L1sum;
double maxerror;
const double tolerance = 1e-10; // 迭代收斂判據
bool steadystate;

double cfd[6] = {(10.0/12.0) , (-15.0/12.0) , (-4.0 / 12.0) , (14.0 / 12.0) ,  (-6.0 / 12.0) , (1.0/12.0) } ;
//四階精度差分係數可依序排列做前向差分與後向差分 

// 解析解
double T_analytical_fixed(double x_pos, double y_pos){
    return sin(pi * x_pos) * (sinh(pi * y_pos) / sinh(pi)) + 10.0;
}

double T_analytical(int k){
    int i = ((k-1) % NX) + 1;
    int j = ((k-1) / NX) + 1;
    
    double x_pos = (i-1) * dx;
    double y_pos = (j-1) * dy;

    return T_analytical_fixed(x_pos, y_pos);
}

// 上邊界邊界條件
double T_up(int i){
    double x_pos = (i-1) * dx;
    return 10.0 + sin(pi * x_pos);
}


//初始化迭代矩陣
void initial(vector<vector<double> >& a, vector<double>& b, int n) {
    
    double H = 1.0 / (6*dx*dx); //用於九點方程 
    double H2 = 1.0/ (dx*dx) ; //四階前向差分或四階精度後項差分 
    //1/(6*dx*dx) = 1/(6*dy*dy)
    double a1 = -20 ; //本點係數採用a1*H
	double a2 = 4.0 ; //直角方向臨計算點a2*H
	double a3 = 1.0 ; //斜角方向臨計算點a3*H
    
    

    // 初始化矩陣
    for(int i = 0; i <= n+1; i++) {
        for(int j = 0; j <= n+1; j++) {
            a[i][j] = 0.0;
        }
        b[i] = 0.0;
        x[i] = 0.0;
        x_old[i] = 0.0;
    }
    
    cout << "Setting boundary conditions..." << endl;

    // 左邊界條件 - 直接賦值
    // 左邊界 (i=1)
    for(int j = 1; j <= NY; j++){
        int idx = (j-1)*NX + 1;
        a[idx][idx] = 1.0; 
        b[idx] = T_left;
    }

    // 右邊界 (i=NX)
    for(int j = 1; j <= NY; j++){
        int idx = (j-1)*NX + NX;
        a[idx][idx] = 1.0;
        b[idx] = T_right;
    }

    // 下邊界 (j=1)
    for(int i = 1; i <= NX; i++){
        int idx = i;
        a[idx][idx] = 1.0;
        b[idx] = T_bottom;
    }

    // 上邊界 (j=NY)
    for(int i = 1; i <= NX; i++){
        int idx = (NY-1)*NX + i;
        a[idx][idx] = 1.0 ;
        b[idx] = T_up(i);
    }
    
    cout << "Setting interior points with consistent 6th-order 9 points difference shceme ..." << endl;
    for(int j = 3; j <= NY-2; j++) {
        for(int i = 3; i <= NX-2; i++) {
            int idx = (j-1)*NX + i;
            a[idx][idx] = a1*H ;  
            a[idx][idx+1] = a2*H ;//(1)右 
            a[idx][idx+NX] = a2*H ;//(2)上
			a[idx][idx-1] = a2*H ;//(3)左 
			a[idx][idx-NX] = a2*H ;//(4)下
			a[idx][idx+1+NX] = a3*H ;//(5)右上
			a[idx][idx-1+NX] = a3*H ;//(6)左上
			a[idx][idx-1-NX] = a3*H ;//(7)左下
		    a[idx][idx+1-NX] = a3*H ;//(8)右下 
            b[idx] = 0.0;
        }
    }
    double A = -2.0 * ((1.0/(dx*dx))+(1.0/(dy*dy))) ;
    double B = (1.0/(dx*dx)) ;
    double C = (1.0/(dy*dy)) ;
	//左邊界計算點 
	//處理方式:東西方向採用單邊插分，南北方向採用中心差分
	for(int j = 3 ; j <= NY-2 ; j++){
		int idx = (j-1)*NX + 2 ; //i = 2 ;
		b[idx] = 0.0 ;
		int idx_plus = idx-1 ;
		for(int ac = 0 ; ac <= 5 ; ac++){
			a[idx][idx_plus+ac] = H2*cfd[ac] ; 
		}
		a[idx][idx] += -2.0*C ;
		a[idx][idx-NX] = C ;
		a[idx][idx+NX] = C ;
		
	} 
	//右邊界計算點 
	//處理方式:東西方向採用單邊插分，南北方向採用中心差分
	for(int j = 3 ; j <= NY-2 ; j++){
		int idx = (j-1)*NX + (NX-1) ; // i = (NX-1)
		b[idx] = 0.0 ;
		int idx_plus = idx+1 ;
		for(int ac = 0 ; ac <= 5 ; ac++){
			a[idx][idx_plus-ac] = H2*cfd[ac] ;
		}
		a[idx][idx] += -2.0*C ;
		a[idx][idx-NX] = C ;
		a[idx][idx+NX] = C ;
	} 
    //下邊界計算點
	//處理方式:南北方向採用單邊差分，東西方向採用中心差分
	for(int i = 3 ; i <= NX-2 ; i++){
		int idx = (2-1)*NX + i ; // j = 2
		b[idx] = 0.0 ;
		int idx_plus = idx - NX ;
		for(int ac = 0 ; ac <= 5 ; ac++){
			a[idx][idx_plus+NX*ac] = H2*cfd[ac] ;
		}
		a[idx][idx] += -2.0*B ;
		a[idx][idx-1] = B ;
		a[idx][idx+1] = B ; 
	} 
	//上邊界計算點
	//處理方式:南北方向採用單邊差分，東西方向採用中心差分
	for(int i = 3 ; i <= NX-2 ; i++){
		int idx = (NY-2)*NX + i ; // j = (NY-2)
		b[idx] = 0.0 ;
		int idx_plus = idx + NX ;
		for(int ac = 0 ; ac <= 5 ; ac++){
			a[idx][idx_plus-NX*ac] = H2*cfd[ac] ; 
		}
		a[idx][idx-1] = B ;
		a[idx][idx] += -2.0*B ;
		a[idx][idx-1] = B ;
		a[idx][idx+1] = B ; 
	}
	
    //左下角點
	int p1 = (2-1)*NX + 2 ;
    b[p1] = 0.0 ;
    a[p1][p1] = A ;
	a[p1][p1-1] = B ;
	a[p1][p1+1] = B ;
	a[p1][p1-NX] = C ;
	a[p1][p1+NX] = C ;
    //右下角點
	int p2 = (2-1)*NX + (NX-1) ; 
	b[p2] = 0.0 ;
	a[p2][p2] = A ;
	a[p2][p2-1] = B ;
	a[p2][p2+1] = B ;
	a[p2][p2-NX] = C ;
	a[p2][p2+NX] = C ;
	//左上角點
	int p3 = (NY-2)*NX + 2;
    b[p3] = 0.0;
    a[p3][p3] = A;        // 修正：應該是 a[p3][p3] 而不是 a[p3][p3-1]
    a[p3][p3+1] = B;
    a[p3][p3-1] = B;      // 這行是正確的
    a[p3][p3-NX] = C;
    a[p3][p3+NX] = C;
	int p4 = (NY-2)*NX + (NX-1);
    b[p4] = 0.0;
    a[p4][p4] = A;
    a[p4][p4-1] = B;
    a[p4][p4+1] = B;
    a[p4][p4-NX] = C;
    a[p4][p4+NX] = C;
    cout << "Matrix initialization completed." << endl;
    cout << "Total equations: " << n << endl;
    cout << " for 9 points difference shceme(4-order bound) ,Interior  points: " << (NX-4)*(NY-4) << endl;
}

//超鬆弛迭代法_自適應性鬆弛因子
void SOR(vector<vector<double> >& a, vector<double>& b, vector<double>& x, int n) {
    // 自適應松弛因子
    double omega;
    if(NX <= 21) omega = 0.5;
    else if(NX <= 41) omega = 0.8;
    else omega = 1.2;
    
    //保存舊的解
    for(int k = 1; k <= n; k++) {
        x_old[k] = x[k];
    }

    // SOR迭代
    for(int k = 1; k <= n; k++) {
        if(fabs(a[k][k]) < 1e-15) continue; // 跳過奇異矩陣
        
        double sum = 0;
        for(int p = 1; p <= n; p++) {
            if(p != k) {
                sum += a[k][p] * x[p];
            }
        }
        double x_new = (b[k] - sum) / a[k][k];
        x[k] = x_old[k] + omega * (x_new - x_old[k]);
    }

    // 計算最大收斂誤差
    maxerror = 0;
    for(int k = 1; k <= n; k++) {
        double error = fabs(x[k] - x_old[k]);
        if(maxerror < error) {
            maxerror = error;
        }
    }

    // 計算L1誤差
    double sum = 0;
    for(int k = 1; k <= n; k++) {
        sum += fabs(x[k] - T_analytical(k));
    }
    L1sum = sum / double(n);
}


void output(int m) {
    for(int j = 1; j <= NY; j++){
        for(int i = 1; i <= NX; i++){
            T[i-1][j-1] = x[(j-1)*NX + i];
        }
    }
    
    ostringstream name;
    name << "9 points difference shceme(4-order bound)_" << NX << "x" << NY << "_" << setfill('0') << setw(6) << m << ".vtk";
    ofstream out(name.str().c_str());
    
    out << "# vtk DataFile Version 3.0\n";
    out << "9 points difference shceme(4-order bound)\n";
    out << "ASCII\n";
    out << "DATASET STRUCTURED_POINTS\n";
    out << "DIMENSIONS " << NX << " " << NY << " 1\n";
    out << "ORIGIN 0.0 0.0 0.0\n";
    out << "SPACING " << dx << " " << dy << " 1.0\n";
    out << "POINT_DATA " << NX * NY << "\n";
    
    out << "SCALARS Temperature double 1\n";
    out << "LOOKUP_TABLE default\n";
    for(int j = 0; j < NY; j++) {
        for(int i = 0; i < NX; i++) {
            out << scientific << setprecision(6) << T[i][j] << "\n";
        }
    }
    
    out.close();
    cout << "VTK document output: " << name.str() << endl;
}

void output_gnuplot_data() {
    bool valid_data = true;
    for(int grid_idx = 0; grid_idx < 4; grid_idx++) {
        if(FinalL1error[grid_idx] <= 0 || isnan(FinalL1error[grid_idx]) || isinf(FinalL1error[grid_idx])) {
            cout << "Warning: Invalid L1 error for grid " << grid_idx << ": " << FinalL1error[grid_idx] << endl;
            valid_data = false;
        }
    }
    
    if(!valid_data) {
        cout << "Cannot generate convergence analysis due to invalid data." << endl;
        return;
    }

    ofstream data_file("grid_convergence_9 points difference shceme(4-order bound).dat");
    data_file << "# Grid_Size dx log(dx) L1_Error log(L1_Error)" << endl;
    
    for(int grid_idx = 0; grid_idx < 4; grid_idx++) {
        double dx_value = 1.0 / (Nx[grid_idx]-1);
        double log_dx = log(dx_value);
        double log_error = log(FinalL1error[grid_idx]);
        
        data_file << Nx[grid_idx] << "\t" 
                 << scientific << setprecision(6) << dx_value << "\t"
                 << scientific << setprecision(6) << log_dx << "\t"
                 << scientific << setprecision(6) << FinalL1error[grid_idx] << "\t"
                 << scientific << setprecision(6) << log_error << endl;
    }
    data_file.close();
    cout << "Data file output: grid_convergence_9 points difference shceme(4-order bound).dat" << endl;

    double sum_x = 0, sum_y = 0, sum_xy = 0, sum_x2 = 0;
    int n_points = 4;
    
    for(int grid_idx = 0; grid_idx < 4; grid_idx++) {
        double x = log(1.0 / (Nx[grid_idx]-1));
        double y = log(FinalL1error[grid_idx]);
        
        sum_x += x;
        sum_y += y;
        sum_xy += x * y;
        sum_x2 += x * x;
    }
    
    double slope = (n_points * sum_xy - sum_x * sum_y) / (n_points * sum_x2 - sum_x * sum_x);
    double intercept = (sum_y - slope * sum_x) / n_points;

    ofstream gnuplot_file("plot_convergence_9 points difference shceme(4-order bound).plt");
    gnuplot_file << "set terminal png enhanced size 800,600" << endl;
    gnuplot_file << "set output 'grid_convergence_9 points difference shceme(4-order bound).png'" << endl;
    gnuplot_file << "set title 'Improved Grid Convergence Analysis: L1 Error vs Grid Spacing'" << endl;
    gnuplot_file << "set xlabel 'log(dx)'" << endl;
    gnuplot_file << "set ylabel 'log(L1 Error)'" << endl;
    gnuplot_file << "set grid" << endl;
    gnuplot_file << "set key left top" << endl;
    
    double x_min = log(1.0 / (Nx[3]-1));
    double x_max = log(1.0 / (Nx[0]-1));
    double y_ref = log(FinalL1error[1]);
    double x_ref = log(1.0 / (Nx[1]-1));
    
    gnuplot_file << "f(x) = " << slope << " * x + " << intercept << endl;
    gnuplot_file << "g(x) = 4.0 * (x - " << x_ref << ") + " << y_ref << endl;
    
    gnuplot_file << "plot 'grid_convergence_9 points difference shceme(4-order bound).dat' using 3:5 with linespoints pt 7 ps 1.5 lw 2 title sprintf('Improved (slope = %.2f)', " << slope << "), \\" << endl;
    gnuplot_file << "     f(x) with lines lw 2 lc rgb 'red' title sprintf('Linear Fit (slope = %.2f)', " << slope << "), \\" << endl;
    gnuplot_file << "     g(x) with lines lw 2 lc rgb 'green' dashtype 2 title '9 points difference shceme(4-order bound)(slope = 4.0)'" << endl;
    
    gnuplot_file.close();
    cout << "Gnuplot script output: plot_convergence_9 points difference shceme(4-order bound).plt" << endl;
    
    cout << "\n=== Improved Grid Convergence Analysis ===" << endl;
    cout << "Linear regression results:" << endl;
    cout << "Slope = " << fixed << setprecision(3) << slope << " (theoretical 4.0)" << endl;
    cout << "Order of accuracy = " << fixed << setprecision(3) << slope << endl;
}

int main() {
    for(int grid_idx = 0; grid_idx < Nx.size(); grid_idx++) {
        cout << "\n========================================" << endl;
        cout << "Grid size: " << Nx[grid_idx] << "x" << Ny[grid_idx] << endl;

        NX = Nx[grid_idx];
        NY = Ny[grid_idx];
        n = NX * NY;
        dx = 1.0 / (NX-1);
        dy = 1.0 / (NY-1);
        
        cout << "Grid spacing: dx = " << dx << ", dy = " << dy << endl;
        
        a.assign(n+2, vector<double>(n+2, 0.0));
        b.assign(n+2, 0.0);
        x.assign(n+2, 0.0);
        x_old.assign(n+2, 0.0);
        T.assign(NX, vector<double>(NY, 0.0));

        cout << "Program execution started with  9 points difference shceme(4-order bound) ...." << endl;
        steadystate = false;
        initial(a, b, n);  //初始化

        for(G = 0; G < max_G; G++) {
            SOR(a, b, x, n);  // 使用改進的SOR

            if(G % 1000 == 0) {
                cout << "Iteration = " << G;
                cout << ", Convergence error = " << scientific << setprecision(3) << maxerror;
                cout << ", L1 error = " << scientific << setprecision(3) << L1sum << endl;
                
                if(G % 5000 == 0) {
                    output(G);
                }
            }
            
            if(G > 100 && maxerror < tolerance) {
                steadystate = true;
                cout << "Steady state reached!" << endl;
                cout << "Final iteration: " << G << ", Final convergence error : " << maxerror << endl;
                cout << "Final L1 error: " << L1sum << endl;
                FinalL1error[grid_idx] = L1sum;
                break;
            }
        }
        
        if(!steadystate) {
            cout << "Maximum iteration reached!" << endl;
            cout << "Final convergence error: " << maxerror << endl;
            cout << "Final L1 error: " << L1sum << endl;
            FinalL1error[grid_idx] = L1sum;
        }
        
        output(G);
        cout << "Grid size " << NX << "x" << NY << " computation completed" << endl;
        cout << "======================================" << endl;
    }
    
    output_gnuplot_data();
    cout << "\nAll computations completed!" << endl;
    
    return 0;
}