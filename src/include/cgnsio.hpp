#pragma once
#include "block.hpp"
#include "cgnslib.h"
#include "data.hpp"
#include "info.hpp"
#include "macro.hpp"

class CgnsIO {
public:
  void BlockCgnsOutput(Block *block, Info *info);
  void solCgnsOutput(Data *data, Info *info);
  void oneDsolOutput(Info *info);

private:
};