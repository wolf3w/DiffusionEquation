#include "headers/solver.h"

Solver::Solver (int N, int tN, double L, double tMax, double beta, double mu,
                double press, double a, double b, double c, double aMax,
                double kZero, double D, double mZero)
    : N_(N), tN_(tN), L_(L), tMax_(tMax), beta_(beta), mu_(mu), pressure_(press),
      a_(a), b_(b), c_(c), aMax_(aMax), kZero_(kZero), D_(D), mZero_(mZero)
{
    h_   = L_ / (N_ - 1);
    tau_ = tMax_ / tN_;

    /* Инициализация векторов */
    aTMA_.resize(N_ + 1);
    bTMA_.resize(N_ + 1);
    cTMA_.resize(N_ + 1);
    fTMA_.resize(N_ + 1);
    ksi_.resize(N_ + 1);
    phi_.resize(N_ + 1);

    Vel_.resize(tN_ + 1);
    Conc_.resize(tN_ + 1);
    mP_.resize(tN_ + 1);
    aC_.resize(tN_ + 1);

    for (int i = 0; i <= N_; i++)
    {
        Vel_[i].resize(N_ + 1);
        Conc_[i].resize(N_ + 1);
        mP_[i].resize(N_ + 1);
        aC_[i].resize(N_ + 1);
    }

    /* Временное заполнение скоростей */
    for (int i = 1; i <= N_; i++)
    {
        for (int j = 0; j < tN_; j++)
        {
            Vel_[j][i] = 1.;
        }
    }
}

void Solver::processTMA(int tInd)
{
    /* Заполнение диагональной матрицы и свободных членов */
    for (int i = 1; i <= N_; i++)
    {
        double Ad = beta_ * L_ * mZero_ / Vel_[tInd][i],
               Pe = Vel_[tInd][i] * L_ / D_;

        aTMA_[i] = tau_ * Vel_[tInd + 1][i - 1] / (Ad * h_) + tau_ / (Ad * Pe * h_ * h_);
        bTMA_[i] = tau_ / (Ad * Pe * h_ * h_);
        cTMA_[i] = mP_[tInd][i] * (1 + tau_) + tau_ * Vel_[tInd + 1][i - 1] / (Ad * h_)
                   + 2 * tau_ / (Ad * Pe * h_ * h_);
        fTMA_[i] = mP_[tInd][i] * Conc_[tInd][i] + tau_ * a_ / ( b_ * (aMax_ - a_) );
    }

    /* Прямой ход прогонки */
    ksi_[1] = 0;
    phi_[1] = 1.;
    for (int i = 2; i <= N_; i++)
    {
        double den = cTMA_[i - 1] - aTMA_[i - 1] * ksi_[i - 1];
        ksi_[i] = phi_[i - 1] / den;
        phi_[i] = (aTMA_[i - 1] * phi_[i - 1] + fTMA_[i - 1]) / den;
    }

    /* Обратный ход прогонки */
    Conc_[tInd][N_] = phi_[N_] / (1 - ksi_[N_]);
    for (int i = N_ - 1; i != 0; i--)
    {
        Conc_[tInd][i] = ksi_[i + 1] * Conc_[tInd][i + 1] + phi_[i + 1];
    }
}
