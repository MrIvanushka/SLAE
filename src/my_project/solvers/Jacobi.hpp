//
// Created by petrov on 19.02.2022.
//

#ifndef SLAE_JACOBI_HPP
#define SLAE_JACOBI_HPP
#include "../sparse/CSR.hpp"
#include "../utility/Norm.hpp"

template<typename T>
std::vector<T> Jacobi(const CSR<T> &A, const std::vector<T> &b, const std::vector<T>& initialState, const T& tolerance)
{
    std::vector<T> tempState(initialState.size());
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
            tempState[i] = (b[i] - sum) / A(i, i);
            currentState = tempState;
        }
    }
    return currentState;
}

#endif//SLAE_JACOBI_HPP
