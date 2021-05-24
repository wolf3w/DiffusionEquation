#include "headers/solver.h"

Solver::Solver (int N, int Nzero, int tMax)
    : N_(N), Nzero_(Nzero), tMax_(tMax)
{
    h_   = 1. / N_;
    tau_ = tMax_ / Nzero_;
}
