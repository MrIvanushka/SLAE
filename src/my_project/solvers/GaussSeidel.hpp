//
// Created by petrov on 19.02.2022.
//

#ifndef SLAE_GAUSSSEIDEL_HPP
#define SLAE_GAUSSSEIDEL_HPP
#include "../sparse/CSR.hpp"
#include "../utility/Norm.hpp"

template<typename T>
std::vector<T> GaussSeidel(const CSR<T> &A, const std::vector<T> &b, const std::vector<T>& initialState, const T& tolerance)
{
    T sum;
    std::vector<T> currentState;

    for(auto r = A * initialState - b; Norm::scoreSecondNorm(r) > tolerance; r = A * currentState - b)
    {
        for(auto i = 0u; i < A.sizeH(); i++)
        {
            sum = static_cast<T>(0);
            auto skip = A._rows[i];

            for(auto k = skip; k < A._rows[i+1]; ++k)
            {
                if(A._cols[k] != i)
                    sum += A._values[k] * b[i];
            }
            currentState[i] = (b[i] - sum) / A(i, i);
        }
    }
    return currentState;
}

#endif//SLAE_GAUSSSEIDEL_HPP
