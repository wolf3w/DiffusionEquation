#ifndef SOLVER_H
#define SOLVER_H

#include <vector>

class Solver
{
private:
    /*! Переменные для сетки */
    int N_, tN_;
    double L_, tMax_;
    double h_, tau_;

    /*! Физические переменные */
    double beta_, mu_, pressure_, a_, b_, c_,
           aMax_, kZero_, D_, mZero_;

    /*! Трёхдиагональная матрица */
    std::vector<double> aTMA_, bTMA_, cTMA_, fTMA_;
    /*! Прогоночные коэффициенты */
    std::vector<double> ksi_, phi_;
    /*! Скорость фильтрации, концетрация вещества и пористость */
    std::vector<std::vector<double>> Vel_, Conc_, mP_, aC_;
public:

    Solver(int N, int tN, double L, double tMax, double beta, double mu,
           double press, double a, double b, double c, double aMax,
           double kZero, double D, double mZero);

    Solver() = delete;
    Solver(Solver const& arg) = delete;
    Solver(Solver&& arg) = delete;
    ~Solver() {}

    /*! @brief Найти концетрацию методом прогонки
     *  @param tInd Индекс по времени (обычно j) */
    void processTMA (int tInd);
};

#endif // SOLVER_H
