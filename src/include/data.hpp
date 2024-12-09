#pragma once
#include "macro.hpp"

class Data {
public:
    Data() {};
    Data(int, int);
    Data(const Data&);

    void solInit(int, int);
    void init(int, int);
    void setValue(std::vector<real>);
    void setvarName(std::vector<std::string>&&);
    // void setDim(int,std::vector<int>);
    real& operator()(int, int);
    real& operator[](int);

    // void operator= (Data&);
    void operator+=(std::vector<real>);
    void setZeros();
    int size();
    int getNVar();
    int getN();

    std::vector<real> getIVar(int);
    std::vector<real>::iterator begin();
    std::vector<real>::iterator random(int i);
    std::vector<real>::iterator end();

    // for global LF flux in burgers equation
    real maxElement(int);

    real getL2(int ivar);
    real getLinf(int ivar);
    std::vector<std::string> varName;

private:
    std::vector<real> data;
    int n = 0;
    int nVar = 0;
};

class DataReader {
public:
    DataReader() { }
    DataReader(int n_, int i0_, int offset_, int idim_, Data* data_)
        : n(n_)
        , offset(offset_)
        , data(data_)
        , idim(idim_)
    {
        i0 = i0_;
    }

    std::vector<real>::iterator operator()(int i)
    {
        return data->random(i0 + offset * i);
    }
    int getN() { return n; }
    int getNVar() { return data->getNVar(); }
    int getOffset() { return offset * data->getNVar(); }
    int idim;

protected:
    int n, i0, offset;

    Data* data;
};
