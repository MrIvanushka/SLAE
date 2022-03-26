//
// Created by petrov on 12.02.2022.
//

#ifndef SLAE_DENSEMATRIX_HPP
#define SLAE_DENSEMATRIX_HPP
#include <vector>
#include <set>
#include "../utility/triplet.hpp"
#include "../Exceptions/UndeclaredElementCallException.hpp"

template<typename T>
class DenseMatrix {
public:
    using elm_t = T;
    using idx_t = std::size_t;

private:

    std::vector<T> matrix;
    idx_t H, W;

public:
    DenseMatrix(const idx_t &h, const idx_t &w) : matrix(h * w), H(h), W(w) {}

    DenseMatrix(const idx_t &h, const idx_t &w, const std::set<Triplet<T>> &in) : matrix(h * w), H(h), W(w) {
        for (auto &element: in) {
            matrix[element.j * w + element.i] = element.value;
        }
    }

    elm_t &operator()(const idx_t &i, const idx_t &j) {

#ifndef NDEBUG
        if(i >= W || j >= H)
            throw Slae::UndeclaredElementCallException("Cannot find element in dense matrix");
#endif

        return matrix[j * W + i];
    }

    const elm_t &operator()(const idx_t &i, const idx_t &j) const
    {
        return this(i, j);
    }

    [[nodiscard]] const idx_t &sizeH() const {
        return H;
    }

    [[nodiscard]] const idx_t &sizeW() const
    {
        return W;
    }

    void swap(const idx_t& first, const idx_t& second) {
        for (auto i = 0u; i < W; i++)
        {
            elm_t temp = matrix[first * W + i];
            matrix[first * W + i] = matrix[second * W + i];
            matrix[second * W + i] = temp;
        }
    }

    void deleteLastRow()
    {
        H -= 1;
        matrix.resize(W * H);
    }

};
#endif//SLAE_DENSEMATRIX_HPP
