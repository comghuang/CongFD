#pragma once

#include "macro.hpp"
#include <concepts>

const real EPSILON = 1e-10;

template <typename F, size_t ncomp>
concept PositiveQuantityFunc = requires(F f, const std::array<real, ncomp>& W) {
    { f(W) } -> std::same_as<real>; // 要求 f(W) 返回 real 类型
};

template <size_t ncomp, PositiveQuantityFunc<ncomp> PositiveQuantityFuncType>
 void variable_positive_limiter(
    std::array<real,ncomp>const& W, // 原始变量值
    std::array<real, ncomp>& W_interpolated, // 高阶插值后的变量值
    PositiveQuantityFuncType computePositiveQuantity, // 计算正值的函数
    real EPSILON = 1e-10 // 正值的阈值
) {
    // 计算插值后变量中的正值
    real positiveQuantity_W_interpolated = computePositiveQuantity(W_interpolated);

    // 如果插值后的正值大于或等于阈值，则无需调整
    if (positiveQuantity_W_interpolated >= EPSILON) {
        return; // 直接返回，无需调整
    }

    // 计算原始变量中的正值
    real positiveQuantity_W = computePositiveQuantity(W);

    // 如果原始变量中的正值小于阈值，则完全使用原始值
    if (positiveQuantity_W < EPSILON) {
        for (size_t i = 0; i < ncomp; ++i) {
            W_interpolated[i] = W[i];
        }
        return; // 直接返回，无需进一步计算
    }

    // 如果插值后的正值小于阈值，则计算调整系数 theta
    real theta = (EPSILON - positiveQuantity_W) / (positiveQuantity_W_interpolated - positiveQuantity_W);

    // 调整插值后的变量值
    for (size_t i = 0; i < ncomp; ++i) {
        W_interpolated[i] = (1 - theta) * W[i] + theta * W_interpolated[i];
    }
}

