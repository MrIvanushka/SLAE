//
// Created by petrov on 19.02.2022.
//

#ifndef SLAE_NORM_HPP
#define SLAE_NORM_HPP

#include <vector>
#include <cmath>

class Norm {
public:
    template<typename T>
    static const T &scoreInfiniteNorm(const std::vector<T> &vector) {
        T norm = static_cast<T>(0);
        T elementAbs;

        for(const auto& element : vector)
        {
            elementAbs = std::abs(element);

            if(elementAbs > norm)
                norm = elementAbs;
        }
        return norm;
    }

    template<typename T>
    static const T &scoreFirstNorm(const std::vector<T> &vector) {
        T norm = static_cast<T>(0);

        for(const auto& element : vector)
        {
            norm += element;
        }
        return norm;
    }

    template<typename T>
    static const T &scoreSecondNorm(const std::vector<T> &vector) {
        T norm = static_cast<T>(0);

        for(const auto& element : vector)
        {
            norm += element * element;
        }
        return std::sqrt(norm);
    }
};

#endif//SLAE_NORM_HPP
