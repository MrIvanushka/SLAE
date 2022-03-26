//
// Created by ivanb on 26.03.2022.
//

#ifndef MY_PROJECT_CONJUGATEGRADIENT_HPP
#define MY_PROJECT_CONJUGATEGRADIENT_HPP
#include "../sparse/CSR.hpp"
#include "../utility/Norm.hpp"
#include "../utility/overloads.hpp"

class ConjugateGradient {
public:
    template<typename T>
    static std::vector<T> solve(const CSR<T> &A, const std::vector<T> &b, const std::vector<T> &initialState,
                                     const T &tolerance) {
        std::vector<T> currentState;
        std::vector<T> previousResult;
        std::vector<T> result = A * initialState - b;
        std::vector<T> d = result;
        bool go = true;

        while(go)
        {
            currentState -= (result * result / (d * (A * d))) * d;
            previousResult = result;
            result = A * currentState - b;

            if(Norm::scoreSecondNorm(result) > tolerance)
            {
                d = result + (result*result)/ (previousResult * previousResult) * d;
            }
            else
            {
                go = false;
            }
        }

        return result;
    }
};

#endif //MY_PROJECT_CONJUGATEGRADIENT_HPP
