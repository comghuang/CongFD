#pragma once

#include <array>
#include <vector>

typedef std::array<real,2> arr2;
std::vector<real> roeFlux1D(arr2 r, arr2 u, arr2 p, arr2 H, arr2 RT);
std::vector<real> roeFlux1D2(arr2 r, arr2 u, arr2 p, arr2 H, arr2 RT);
std::vector<real> HLLCFlux1D(arr2 r, arr2 u, arr2 p, arr2 H, arr2 RT);