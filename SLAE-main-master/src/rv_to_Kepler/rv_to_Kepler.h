//
// Created by perseverance on 31.03.2023.
//

#ifndef SATELLITE_SYSTEM_KEPPLER_TO__R_V_H
#define SATELLITE_SYSTEM_KEPPLER_TO__R_V_H
#define pi 3.14159265358979323846


#include <iostream>
#include <cmath>
#include <limits>
#include "../src/Vector.h"

template<typename T>
struct basis {
    std::vector<T> I;// = {1, 0, 0};
    std::vector<T> J;// = {0, 1, 0};
    std::vector<T> K;// = {0, 0, 1};
};

template<typename T>
std::vector<T> operator-(const std::vector<T> &a, const std::vector<T> &b) {
    std::vector<T> c(a.size());
    for (auto i = 0; i < a.size(); ++i) {
        c[i] = a[i] - b[i];
    }
    return c;
}


//Works
template<typename T>
std::vector<T> operator+(const std::vector<T> &a, const std::vector<T> &b) {
    std::vector<T> c(a.size());
    for (auto i = 0; i < c.size(); ++i) {
        c[i] = a[i] + b[i];
    }
    return c;
}

//Works
template<typename T>
std::vector<T> operator*(const T &a, const std::vector<T> v) noexcept {
    std::vector<T> result(v.size());
    for (auto i = 0; i < v.size(); ++i) {
        result[i] = a * v[i];
    }
    return result;
}

//Works
template<typename T>
T scal_dot(const std::vector<T> &a, const std::vector<T> &b) noexcept {
    T scal = 0;
    if (a.size() == b.size()) {
        for (auto i = 0; i < a.size(); ++i) {
            scal += a[i] * b[i];
        }
    } else { std::cout << "Different size of vectors in scal dot" << std::endl; }
    return scal; //вместо std::cout ASSERT загугли
}

// Works
template<typename T>
std::vector<T> vect_dot(const std::vector<T> &a, const std::vector<T> &b, const basis<T> &basis) noexcept {
    std::vector<T> vect(a.size());
    vect = (a[1] * b[2] - a[2] * b[1]) * basis.I + (a[2] * b[0] - a[0] * b[2]) * basis.J +
           (a[0] * b[1] - b[0] * a[1]) * basis.K;
    return vect;
}

template<typename T>
std::vector<T> vect_dot(const std::vector<T> &a, const std::vector<T> &b) noexcept {
    std::vector<T> vect(a.size());
    vect[0] = (a[1] * b[2] - a[2] * b[1]);
    vect[1] = (a[2] * b[0] - a[0] * b[2]);
    vect[2] = (a[0] * b[1] - b[0] * a[1]);
    return vect;
} //array (знаем размер массива)

//Works
template<typename T>
T norm(const std::vector<T> &a) {
    T norm = 0;
    for (auto i: a) {
        norm += pow(i, 2); // sqrt и обычное возведение в степень в 20 раз быстрее
    }
    return pow(norm, 0.5);
}

struct KeplerOrbit{
        double inclination;
        double semimajorAxis;
    };

/*struct CartesianOrbit{
    Vector3d position;
    Vector3d velocity;
};*/

/*
KeplerOrbit cartToKep(const CartesianOrbit& cart, double mu);
CartesianOrbit kepToCart(const KeplerOrbit& kep, double mu);

*/

template<typename T>
class Kepler_elements {
private:
    const double mu = 398600.4418; // Earth's mass dot G // гравитационный параметр принимать на вход в конструктор
    T p, a, ex, i, W, w, nu, ksi;
    //double cos_w;
    std::vector<T> e;
    std::vector<T> h;
    std::vector<T> N;
    T norm_N;
    T norm_h;
    std::vector<T> orbital_elements;

public:
    Kepler_elements(const basis<T> &basis, const std::vector<T> &r, const std::vector<T> &velocity) noexcept {
        h = vect_dot(r, velocity);
        norm_h = norm(h);
        N = vect_dot(basis.K, h);
        norm_N = norm(N);
        T norm_r = norm(r);
        T mu_r = mu / norm_r;
        T v_square = scal_dot(velocity, velocity);
        ksi = v_square / 2 - mu_r;

        // Вектор эксцентриситета + сам эксцентриситет
        e = 1 / mu * ((v_square - mu_r) * r - scal_dot(r, velocity) * velocity);
        ex = norm(e);

        // Полуоси
        if (ex != 1.0) {
            a = -mu / (2 * ksi);
            //p = a * (1 - pow(ex, 2));
            p = scal_dot(h, h) / mu;
        } else {
            p = scal_dot(h, h) / mu;
            a = std::numeric_limits<int>::max();
        }

        // Наклонение
/*        T cos_i = h[2] / norm(h);
        T sin_i = pow((pow(h[0], 2) + pow(h[1], 2)), 0.5);
        i = atan2(sin_i, cos_i);*/
        //double cos_i = h[2]/ norm(h);
        i = acos(h[2] / norm(h));

        // Долгота восходящего узла
        //double cos_W = N[0] / norm(N);
        W = acos(N[0] / norm_N);
        if (N[1] < 0) {
            W = 2 * pi - W;
        }

        //Аргумент перицентра
        w = acos(scal_dot(N, e) / (norm(N) * ex));
        if (e[2] < 0) {
            w = 2 * pi - 2;
        }

        //Истинная аномалия
        nu = acos(scal_dot(e, r) / (norm_r * ex));
        if (scal_dot(r, velocity) < 0) {
            nu = 2 * pi - nu;
        }

        //Вектор элементов
        orbital_elements = {i, W, w, ex, a, nu, p};
    }

    std::vector<T> &get_Excentryvect() {
        return e;
    }

    std::vector<T> &get_h() {
        return h;
    }

    T get_h_norm() {
        return norm_h;
    }

    std::vector<T> &get_N() {
        return N;
    }

    T get_N_norm() {
        return norm_N;
    }

    T get_Excentry() {
        return ex;
    }

    T get_ksi() {
        return ksi;
    }

    T get_a() {
        return a;
    }

    T get_p() {
        return p;
    }

    T get_i() {
        return i;
    }

    double get_W() {
        return W;
    }

    double get_w() {
        return w;
    }

    double get_nu() {
        return nu;
    }

    std::vector<T> &get_orbital_elements() {
        return orbital_elements;
    }

    void cout_orbital_elements() {
        std::vector<std::string> elem = {"Inclination, rad", "Ascending node longitude, rad",
                                         "Argument of perigee, rad", "Eccentricity", "Semimajor axis, km",
                                         "True anomaly, rad", "Semiparametr, km"};
        for(auto i = 0; i < orbital_elements.size(); ++i){
            std::cout << elem[i] << " = " << orbital_elements[i] << std::endl;
        }
    }
};

#endif //SATELLITE_SYSTEM_KEPPLER_TO__R_V_H