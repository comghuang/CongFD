#pragma once
#include "macro.hpp"
#include "cgnslib.h"
#include "block.hpp"
#include "data.hpp"
#include "info.hpp"


class CgnsIO
{
    public:
    void BlockCgnsOutput(std::shared_ptr<Block> block,std::shared_ptr<Info> info);

    private:
};