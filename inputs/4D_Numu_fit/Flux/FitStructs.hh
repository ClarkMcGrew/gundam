#ifndef FITSTRUCTS_HH
#define FITSTRUCTS_HH

namespace xsllh
{

struct FitBin
{
    double D1low, D1high;
    double D2low, D2high;
    double D3low, D3high;
    double D4low, D4high;

    FitBin() : D1low(0), D1high(0), D2low(0), D2high(0), D3low(0), D3high(0), D4low(0), D4high(0) {}
    FitBin(const double D1_L, const double D1_H,
           const double D2_L, const double D2_H,
           const double D3_L, const double D3_H,
           const double D4_L, const double D4_H)
          : D1low(D1_L), D1high(D1_H),
            D2low(D2_L), D2high(D2_H),
            D3low(D3_L), D3high(D3_H),
            D4low(D4_L), D4high(D4_H)
          {}
};

enum Flavor
{
    kNumu = 0,
    kNumubar = 1,
    kNue = 2,
    kNuebar = 3,
};

}

#endif
