#include "Output.h"
#include "Thomas.h"
using namespace std;

constexpr double ReD = 5000.00;
constexpr double E = 1.0;
constexpr double NY = 20;  // 共有 NY+1 个网格点
constexpr double deltaY = 1.00 / NY;
constexpr int timesteps = 1000;

int main() {
    double deltaT = E * ReD * deltaY * deltaY;
    double A0 = -0.5 * E;
    double B0 = 1.0 + E;
    // u[0] = u1,...,u[NY]=u_(NY+1)
    vector<double> u(NY + 1, 0.0);
    // K[0] = K2,K[1] = K3,...
    vector<double> K(NY - 1, 0.0);
    // A,B,D分别为下对角元,对角元,右侧元素
    vector<double> A(NY - 1, A0), B(NY - 1, B0), D(NY - 1, 0.0);
    // 托马斯算法的解向量
    vector<double> S(NY - 1, 0.0);
    // 需要输出的时间步列表
    vector<int> out_steps = {1,12,36,60,240};

    // 初值条件
    u[NY] = 1.0;

    for (int tt = 0; tt < timesteps; tt++) {
        for (int i = 0; i < NY - 1; i++) {
            K[i] = (1.0 - E) * u[i + 1] + 0.50 * E * (u[i] + u[i + 2]);
            if (i < NY - 2)
                D[i] = K[i];
            else
                D[i] = K[i] - A0;
        }

        S = ThomasAlgorithm(A, B, A, D);
        for (int i = 0; i < NY - 1; i++) {
            u[i + 1] = S[i];
        }

        // --- 输出控制：检查当前步数是否需要保存结果 ---
        int current_step = tt + 1;  // 当前已完成的时间步编号
        if (std::find(out_steps.begin(), out_steps.end(), current_step) != out_steps.end()) {
            OutputU(u, current_step, deltaY, E);
        }
    }

    return 0;
}
