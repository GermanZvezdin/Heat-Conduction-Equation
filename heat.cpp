//
// Created by German Zvezdin on 30/10/2020.
//

#include "z-arrays.hpp"
#include <iostream>
#include <iomanip>
ArrD<double, 1> real_solurion(int N = 10000, double x0 = 0.525) {
    ArrD<double, 1> x;
    x.init(ind(N + 1));
    double h = 1.0 / N;
    double a1 = x0 * x0 + 1;
    double b1 = exp(-x0);
    double c1 = 1;
    double d1 = c1 / b1;

    double a2 = x0;
    double b2 = exp(-x0);
    double c2 = x0 * x0 * x0;
    double d2 = c2 / b2;

    double k1 = sqrt(b1 / a1);
    double k2 = sqrt(b2 / a2);

    double A1 = exp(-k1 * x0) - exp(k1 * x0);
    double B1 = exp(-2 * k2 + k2 * x0) - exp(-k2 * x0);
    double D1 = d1 * (exp(k1 * x0) - 1) + exp(k2 * x0 - k2) * (1 - d2) + d2;

    double A2 = -a1 * (k1 * exp(k1 * x0) + k1 * exp(-k1 * x0));
    double B2 = a2 * (k2 * exp(-2 * k2 + k2 * x0) + k2 * exp(-k2 * x0));
    double D2 = a1 * d1 * k1 * exp(k1 * x0) + a2 * k2 * exp(-k2 + k2 * x0) * (1 - d2);

    double mainopred = A1 * B2 - B1 * A2;

    double const2 = (D1 * B2 - B1 * D2) / mainopred;
    double const4 = (A1 * D2 - A2 * D1) / mainopred;
    double const1 = -d1 - const2;
    double const3 = exp(-k2) * (1 - d2 - const4 * exp(-k2));

    for(int i = 1; i < N + 1; i++) {
        x[ind(i)] = i * h;
    }

    ArrD<double, 1> u;
    u.init(ind(N + 1));
    u[ind(0)] = 0;
    u[ind(N)] = 1;
    for (int i = 1; i < N; i++) {
        if (x[ind(i)] < x0)
            u[ind(i)] = const1 * exp(k1 * x[ind(i)]) + const2 * exp(-k1 * x[ind(i)]) + d1;
        else
            u[ind(i)] = const3 * exp(k2 * x[ind(i)]) + const4 * exp(-k2 * x[ind(i)]) + d2;
    }
    std::cout << const1 << " " << const2 << " " << const3 << " " << const4 << std::endl;
    std::cout << k1 << " " << k2 <<" " << d1<< " " << d2 << std::endl;
    FILE*real_sol;
    real_sol = fopen("real_sol.txt", "w");

    for(int i = 0; i < N + 1; i++){
        fprintf(real_sol, "%lf %lf\n", i*h, u[ind(i)]);
    }
    fclose(real_sol);
    return u;

}
ArrD<double, 1> counter_sweeps(int N = 1000, double x0 = 0.525) {
    ArrD<double, 1> q, k, y, u, x, fi;
    ArrD<double, 1> alf, bet, a, b, c;
    double h = 1.0 / N;

    //Инициализация массивов с индексацией по Z-кривой
    q.init(ind(N + 1));
    k.init(ind(N + 1));
    y.init(ind(N + 1));
    u.init(ind(N + 1));
    x.init(ind(N + 1));
    fi.init(ind(N + 1));

    alf.init(ind(N + 2));
    bet.init(ind(N + 2));
    a.init(ind(N + 2));
    b.init(ind(N + 2));
    c.init(ind(N + 2));

    for (int i = 0; i < N + 1; i++) {
        x[ind(i)] = k[ind(i)] = q[ind(i)] = fi[ind(i)] = 0.0;
    }


    for (int i = 1; i < N + 1; i++) {
        x[ind(i)] = i * h;

        if (x[ind(i)] < x0) {
            k[ind(i)] = x0*x0 + 1;
            fi[ind(i)] = 1.0;

        } else {
            k[ind(i)] = x0;
            fi[ind(i)] = x0*x0*x0;
        }

        q[ind(i)] = exp(-x0);
    }
    for (int i = 1; i < N; i++) {
        b[ind(i)] = (k[ind(i + 1)]) / (h * h);
        c[ind(i)] = ((k[ind(i + 1)] + k[ind(i)]) / (h * h)) + q[ind(i)];
        a[ind(i)] = (k[ind(i)]) / (h * h);

    }
    alf[ind(1)] = 0.0;
    bet[ind(1)] = 0.0;

    for (int i = 1; i < N; i++) {
        alf[ind(i + 1)] = b[ind(i)] / (c[ind(i)] - a[ind(i)] * alf[ind(i)]);
        bet[ind(i + 1)] = (fi[ind(i)] + a[ind(i)] * bet[ind(i)]) / (c[ind(i)] - a[ind(i)] * alf[ind(i)]);
    }

    y[ind(N)] = 1.0;
    for (int i = N - 1; i >= 0; i--) {
        y[ind(i)] = alf[ind(i + 1)] * y[ind(i + 1)] + bet[ind(i + 1)];

    }

    FILE *res;
    char path[40];
    sprintf(path, "%d_res.txt", N);
    res = fopen(path, "w");

    for(int i = 0; i < N + 1; i++){
        fprintf(res, "%lf %lf\n", i*h, y[ind(i)]);
    }
    fclose(res);
    return y;

}
int heat(){
    ArrD<double, 1> u, y;
    real_solurion(1000);
    counter_sweeps(1000);
    return 0;
}
int main(){
    heat();
    return 0;
}


