#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <vector>

void OutputU(const std::vector<double> &u, int step, double dy, double E)
{
    // 将 E 格式化为一位小数的字符串（用于文件名）
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(1) << E;
    std::string e_str = oss.str();

    std::string filename = "../output/u_step_" + std::to_string(step) + "_E_" + e_str + ".csv";

    std::ofstream fout(filename);
    if (!fout)
    {
        std::cerr << "无法创建文件: " << filename << std::endl;
        return;
    }

    // 数据列统一保留三位小数
    fout << std::fixed << std::setprecision(3);
    // 恢复数据列的精度为三位小数
    fout << std::fixed << std::setprecision(3);

    fout << "j\ty\tu\n";   // 表头

    for (size_t i = 0; i < u.size(); ++i)
    {
        int j = i + 1;
        double y = i * dy;
        fout << j << "\t" << y << "\t" << u[i] << "\n";
    }

    fout.close();
    std::cout << "已保存 " << filename << std::endl;
}