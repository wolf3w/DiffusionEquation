#ifndef SOLVER_H
#define SOLVER_H

#include <vector>

class Solver
{
private:
    int N_, Nzero_, tMax_;
    double h_, tau_;
public:

    Solver(int N, int Nzero, int tMax);

    Solver() = delete;
    Solver(Solver const& arg) = delete;
    Solver(Solver&& arg) = delete;
    ~Solver() {}
};

#endif // SOLVER_H
