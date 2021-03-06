//
// Created by German Zvezdin on 21/10/2020.
//

#ifndef HEAT_CONDUCTION_EQUATION_Z_ARRAYS_HPP
#define HEAT_CONDUCTION_EQUATION_Z_ARRAYS_HPP
#include <vector>
#include <cmath>
#include <cstdlib>
#include <cstdio>
template <int D> class Ind {
    int p[D];

    void set_x() {}

    template<typename T2, typename ... Args>
    void set_x(const T2 &x, const Args &... xxx) {
        p[D - 1 - sizeof...(Args)] = x;
        set_x(xxx...);
    }

public:
    template<typename ... Args>
    explicit Ind(const Args &... xxx) {
        set_x(xxx...);
    }

    int &operator[](int i) { return p[i]; }

    const int &operator[](int i) const { return p[i]; }

    int operator*(const Ind &b) const {
        int res = 0;
        for (int i = 0; i < D; i++) res += p[i] * b[i];
        return res;
    }

    int prod() const {
        int res = p[0];
        for (int i = 1; i < D; i++) res *= p[i];
        return res;
    }
};

template<typename ... Args> inline Ind<sizeof...(Args)> ind(Args ... args) { return Ind<sizeof...(Args)>(args...); }


template <typename T, int D> class ArrD {
    std::vector<T> data;
    Ind<D> N, mul;
public:
    const Ind<D> &size() const { return N; }

    void init(const Ind<D> &N_) {
        N = N_;
        mul[0] = 1;
        for (int i = 1; i < D; i++) mul[i] = mul[i - 1] * N[i - 1];
        data.resize(N.prod());
    }

    T &operator[](const Ind<D> &pos) { return data[pos * mul]; }

    const T &operator[](const Ind<D> &pos) const { return data[pos * mul]; }
};
#endif //HEAT_CONDUCTION_EQUATION_Z_ARRAYS_HPP
