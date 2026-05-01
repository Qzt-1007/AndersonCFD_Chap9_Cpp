// 不可压库埃特流的压力修正法求解
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

// 3D 数组类型别名 (x, y, component)
using Array3D = std::vector<std::vector<std::vector<double>>>;

// 2D 数组类型别名
using Array2D = std::vector<std::vector<double>>;

// 计算两个二维数组的欧几里得范数（用于判断收敛）
double norm(const Array2D& a, const Array2D& b) {
    double sum = 0.0;
    int nx = a.size();
    int ny = a[0].size();
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            double diff = a[i][j] - b[i][j];
            sum += diff * diff;
        }
    }
    return std::sqrt(sum);
}

// 创建并初始化数组
Array3D createArray3D(int nx, int ny, int ncomp, double initVal = 0.0) {
    return Array3D(nx, std::vector<std::vector<double>>(ny, std::vector<double>(ncomp, initVal)));
}

Array2D createArray2D(int nx, int ny, double initVal = 0.0) { return Array2D(nx, std::vector<double>(ny, initVal)); }

// 复制二维数组
void copyArray2D(const Array2D& src, Array2D& dst) {
    int nx = src.size();
    int ny = src[0].size();
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            dst[i][j] = src[i][j];
        }
    }
}

void outputResult(const Array3D& pset, int i, const std::vector<double>& y) {
    // 构造文件名
    std::ostringstream oss;
    oss << "../output/C2_yuv_i" << i << ".csv";
    std::string filename = oss.str();

    // 打开文件
    std::ofstream outFile(filename);
    if (!outFile.is_open()) {
        std::cerr << "无法打开文件: " << filename << std::endl;
        return;
    }

    // 写入数据，保留 3 位小数
    outFile << std::fixed << std::setprecision(3);
    outFile << "y\tu\tv\n";
    int ynum = y.size();
    for (int j = 0; j < ynum; ++j) {
        double yVal = y[j];
        double uVal = pset[i][j][1];
        double vVal = pset[i][j][2];
        outFile << yVal << "\t" << uVal << "\t" << vVal << std::endl;
    }

    outFile.close();
    std::cout << "断面数据已保存至: " << filename << std::endl;
}

int main() {
    // 几何条件与物性
    constexpr double L = 0.5;         // ft
    constexpr double D = 0.01;        // ft
    constexpr double rho = 0.002377;  // slug/ft^3
    constexpr double mu = 3.74e-7;    // slug/(ft·s)

    // 边界条件
    constexpr double ue = 1.0;        // ft/s, 上边界速度
    constexpr double ve = 0.0;        // ft/s, 下边界速度
    constexpr double pe_prime = 0.0;  // 上边界压力修正量
    constexpr double u0 = 0.0;
    constexpr double v0 = 0.0;
    constexpr double p0_prime = 0.0;
    constexpr double p_inlet_prime = 0.0;
    constexpr double v_inlet = 0.0;
    constexpr double p_outlet_prime = 0.0;

    // 网格参数
    constexpr int xnum = 21;  // x方向网格数（以压力点为基准）
    constexpr int ynum = 11;  // y方向网格数
    constexpr double dx = L / (xnum - 1);
    constexpr double dy = D / (ynum - 1);
    constexpr double dt = 0.001;

    std::vector<double> y(ynum);
    for (int j = 0; j < ynum; ++j) y[j] = j * dy;

    // 主变量 pset (x, y, 3): 0->p, 1->u, 2->v
    Array3D pset = createArray3D(xnum, ynum, 3, 0.0);

    // 计算参数与临时变量
    constexpr int max_iteration = 800;
    constexpr int max_sub_iteration = 400;

    Array2D p_prime = createArray2D(xnum, ynum, 0.0);
    Array2D p_prime_old = createArray2D(xnum, ynum, 0.0);

    // 存储历史数据
    std::vector<Array3D> pset_history(max_iteration);
    std::vector<double> d_history(max_iteration, 0.0);

    // 初始化边界条件
    // 上下边界速度
    for (int i = 0; i < xnum; ++i) {
        pset[i][ynum - 1][1] = ue;  // u 上边界
        pset[i][ynum - 1][2] = ve;  // v 上边界
        pset[i][0][1] = u0;         // u 下边界
        pset[i][0][2] = v0;         // v 下边界
    }
    // 入口边界 v 速度
    for (int j = 0; j < ynum; ++j) {
        pset[0][j][2] = v_inlet;
    }
    // 压力修正量边界
    for (int j = 0; j < ynum; ++j) {
        p_prime[0][j] = p_inlet_prime;
        p_prime[xnum - 1][j] = p_outlet_prime;
    }
    for (int i = 0; i < xnum; ++i) {
        p_prime[i][0] = p0_prime;
        p_prime[i][ynum - 1] = pe_prime;
    }
    // 初始扰动：点 (15,5) 的 v 分量
    pset[15][5][2] = ue / 2.0;

    // 复制初始 p_prime 到 old 数组
    copyArray2D(p_prime, p_prime_old);

    // 主迭代（SIMPLE 算法）
    for (int iter = 0; iter < max_iteration; ++iter) {
        // 保存当前迭代步的场
        pset_history[iter] = pset;

        //  第二步：计算 u*, v*（动量预测）
        //  临时存储 A_star 和 B_star（内部点）
        Array2D A_star = createArray2D(xnum - 2, ynum - 2, 0.0);
        Array2D B_star = createArray2D(xnum - 2, ynum - 2, 0.0);

        // 计算 A_star 和 B_star（内点索引从 1 到 nx-2, 1 到 ny-2）
        for (int i = 1; i <= xnum - 2; ++i) {
            for (int j = 1; j <= ynum - 2; ++j) {
                // 用于 u 方程的对流项系数
                double v_bar = 0.5 * (pset[i][j][2] + pset[i + 1][j][2]);
                double v = 0.5 * (pset[i][j - 1][2] + pset[i + 1][j - 1][2]);
                double conv_u = -(
                    rho * (pset[i + 1][j][1] * pset[i + 1][j][1] - pset[i - 1][j][1] * pset[i - 1][j][1]) / (2.0 * dx) -
                    rho * (pset[i][j + 1][1] * v_bar - pset[i][j - 1][1] * v) / (2.0 * dy));

                double diff_u = mu * ((pset[i + 1][j][1] - 2.0 * pset[i][j][1] + pset[i - 1][j][1]) / (dx * dx) +
                                      (pset[i][j + 1][1] - 2.0 * pset[i][j][1] + pset[i][j - 1][1]) / (dy * dy));

                A_star[i - 1][j - 1] = conv_u + diff_u;

                // 用于 v 方程的对流项系数
                double u_bar = 0.5 * (pset[i][j][1] + pset[i][j + 1][1]);
                double u = 0.5 * (pset[i - 1][j][1] + pset[i - 1][j + 1][1]);
                double conv_v = -(
                    rho * (pset[i][j + 1][2] * pset[i][j + 1][2] - pset[i][j - 1][2] * pset[i][j - 1][2]) / (2.0 * dy) -
                    rho * (pset[i + 1][j][2] * u_bar - pset[i - 1][j][2] * u) / (2.0 * dx));

                double diff_v = mu * ((pset[i + 1][j][2] - 2.0 * pset[i][j][2] + pset[i - 1][j][2]) / (dx * dx) +
                                      (pset[i][j + 1][2] - 2.0 * pset[i][j][2] + pset[i][j - 1][2]) / (dy * dy));

                B_star[i - 1][j - 1] = conv_v + diff_v;
            }
        }

        // 更新内部点的 u, v（预测步）
        for (int i = 1; i <= xnum - 2; ++i) {
            for (int j = 1; j <= ynum - 2; ++j) {
                double dpdx = (pset[i + 1][j][0] - pset[i][j][0]) / dx;
                double dpdy = (pset[i][j + 1][0] - pset[i][j][0]) / dy;
                pset[i][j][1] = (rho * pset[i][j][1] + A_star[i - 1][j - 1] * dt - dt * dpdx) / rho;
                pset[i][j][2] = (rho * pset[i][j][2] + B_star[i - 1][j - 1] * dt - dt * dpdy) / rho;
            }
        }

        // 零阶外插：入口与出口边界（上下边界不变）
        for (int j = 1; j <= ynum - 2; ++j) {
            pset[0][j][1] = pset[1][j][1];
            pset[xnum - 1][j][1] = pset[xnum - 2][j][1];
            pset[xnum - 1][j][2] = pset[xnum - 2][j][2];
        }

        // 第三步：计算压力修正量 p'
        //  系数 a, b, c（常数）
        double a = 2.0 * (dt / (dx * dx) + dt / (dy * dy));
        double b = -dt / (dx * dx);
        double c = -dt / (dy * dy);

        // 计算质量源项 d（内部点大小）
        Array2D d = createArray2D(xnum - 2, ynum - 2, 0.0);
        for (int i = 1; i <= xnum - 2; ++i) {
            for (int j = 1; j <= ynum - 2; ++j) {
                double du = (rho * pset[i][j][1] - rho * pset[i - 1][j][1]) / dx;
                double dv = (rho * pset[i][j][2] - rho * pset[i][j - 1][2]) / dy;
                d[i - 1][j - 1] = du + dv;
            }
        }
        // 保存某点 (15,5) 的质量源项
        d_history[iter] = d[14][4];

        // 子迭代求解 p'（Gauss-Seidel）
        for (int sub = 0; sub < max_sub_iteration; ++sub) {
            // 更新内部点（边界保持不变）
            for (int i = 1; i <= xnum - 2; ++i) {
                for (int j = 1; j <= ynum - 2; ++j) {
                    double rhs = b * p_prime[i + 1][j] + b * p_prime[i - 1][j] + c * p_prime[i][j + 1] +
                                 c * p_prime[i][j - 1] + d[i - 1][j - 1];
                    p_prime[i][j] = -rhs / a;
                }
            }

            // 检查收敛
            if (norm(p_prime, p_prime_old) < 1e-9) {
                break;
            }
            copyArray2D(p_prime, p_prime_old);
        }

        // 第四步：更新压力（低松弛因子 0.1）
        for (int i = 0; i < xnum; ++i) {
            for (int j = 0; j < ynum; ++j) {
                pset[i][j][0] += 0.1 * p_prime[i][j];
            }
        }

        // 子迭代结束后重置 p_prime_old 为当前 p_prime（下一主循环会重新复制）
        copyArray2D(p_prime, p_prime_old);
    }
    // 输出断面 i=15 的数据
    outputResult(pset, 15, y);
    return 0;
}