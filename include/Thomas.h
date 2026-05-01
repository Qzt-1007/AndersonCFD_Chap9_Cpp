#pragma once
#include <iostream>
#include <vector>
#include <stdexcept>

/**
 * 托马斯算法（追赶法）求解三对角线性方程组
 *
 * @param a 下对角线元素，长度 n，a[0] 未使用
 * @param b 主对角线元素，长度 n
 * @param c 上对角线元素，长度 n，c[n-1] 未使用
 * @param d 右端项，长度 n
 * @return  解向量 x，长度 n
 * @throws std::invalid_argument 若输入向量长度不一致或主对角线元素为零
 */
std::vector<double> ThomasAlgorithm(const std::vector<double> &a,
                                    const std::vector<double> &b,
                                    const std::vector<double> &c,
                                    const std::vector<double> &d)
{
    size_t n = b.size();
    // 检查所有向量长度是否匹配
    if (a.size() != n || c.size() != n || d.size() != n)
    {
        throw std::invalid_argument("All input vectors must have the same length.");
    }
    if (n == 0)
        return {}; // 空方程组

    // 拷贝，避免修改原始数据
    std::vector<double> bb = b;
    std::vector<double> dd = d;
    std::vector<double> x(n);

    // 向前消元,将 A 化为上三角
    for (size_t i = 1; i < n; ++i)
    {
        // 计算消去因子 w
        double w = a[i] / bb[i - 1];
        // 更新主对角线和右端项
        bb[i] = bb[i] - w * c[i - 1];
        dd[i] = dd[i] - w * dd[i - 1];
        // 检查 bb[i] 是否接近零，若为零则矩阵奇异
        if (std::abs(bb[i]) < 1e-12)
        {
            throw std::runtime_error("Singular matrix or zero pivot encountered.");
        }
    }

    // 回代求解
    x[n - 1] = dd[n - 1] / bb[n - 1];
    for (size_t i = n - 1; i-- > 0;)
    { // i = n-2, n-3, ..., 0
        x[i] = (dd[i] - c[i] * x[i + 1]) / bb[i];
    }

    return x;
}