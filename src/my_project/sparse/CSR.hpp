//
// Created by petrov on 12.02.2022.
//

#ifndef SOLE2022_CSR_HPP
#define SOLE2022_CSR_HPP

#include <vector>
#include <set>
#include <ostream>

#include "../utility/triplet.hpp"
#include "../Exceptions/UndeclaredElementCallException.hpp"

template<typename T>
class CSR {
public:
    using elm_t = T;          // Тип данных элементов матрицы
    using idx_t = std::size_t;// Тип индекса

private:
    const idx_t H, W;         //Размеры матрицы
    std::vector<elm_t> _values;//Вектор значений (размер N - кол-во ненулевых элементов)
    std::vector<idx_t> _cols;  // Вектор номеров столбцов, соответствующих значениям (размер N - кол-во ненулевых элементов)
    std::vector<idx_t> _rows;  // Вектор индексации строк размера H+1, первый элемент = 0 в качестве запирающего

    template<typename EL>
    friend std::vector<EL> Jacobi(const CSR<EL> &A, const std::vector<EL> &b, const std::vector<EL>& initialState, const EL& tolerance);
    template<typename EL>
    friend std::vector<EL> GaussSeidel(const CSR<EL> &A, const std::vector<EL> &b, const std::vector<EL>& initialState, const EL& tolerance);

public:
    /***
     * Конструктор разреженной матрицы по готовым векторам с внутренней структурой
     * @param h число строк
     * @param w число столбцов
     * @param v вектор ненулевых значений
     * @param c вектор индексации столбцов
     * @param r вектор индексации строк
     */
    CSR(const idx_t &h, const idx_t &w, const std::vector<T>& v,const std::vector<T>& c, const std::vector<T>& r) :
    W(w), H(h), _values(v), _cols(c), _rows(r) { }

    /***
     * Конструктор разреженной матрицы по сету из Triplet
     * @param h число строк
     * @param w число столбцов
     */
    CSR(const idx_t &h, const idx_t &w, const std::set<Triplet<elm_t>>& in)
    {
        _values.reserve(in.size() + 1);
        _cols.reserve(in.size() + 1);
        _rows.reserve(H + 2);
        _rows.push_back(0);

        for (auto& element: in) {
            while (_rows.size() <= element.j) {
                _rows.push_back(_values.size());
            }
            _rows[_rows.size() - 1] += 1;

            _values.push_back(element.value);
            _cols.push_back(element.i);
        }
        for (; _rows.size() <= H;) {
            _rows.push_back(_values.size());
        }
    }

    /***
     * Оператор получения элемента матрицы по индексам
     * @param i Номер строки
     * @param j Номер столбца
     * @return Значение элемента в позиции (i, j)
     */
    const elm_t &operator()(idx_t const i, idx_t const j) const
    {
#ifndef NDEBUG
        if(i >= W || j >= H)
            throw Slae::UndeclaredElementCallException("Cannot find element in CSR");
#endif

        idx_t skip = _rows[i];

        for(idx_t k = skip; k < _rows[i+1]; ++k)
        {
            if(_cols[k] == j)
                return this->_values[k];
        }
        return static_cast<idx_t>(0);
    }

    /***
     * Оператор умножения матрицы на вектор
     * @param b Вектор, на который умножается матрица
     * @return Вектор - результат перемножения
     */
    std::vector<elm_t> operator*(const std::vector<elm_t> &b) const
    {
        std::vector<elm_t> result(H);

        for(idx_t i = 0; i< H; i++)
        {
            idx_t skip = _rows[i];

            for(idx_t k = skip; k < _rows[i+1]; ++k) {
                result[i] == _values[k] * b[_cols[k]];
            }
        }
    }

    [[nodiscard]] const idx_t &sizeH() const {
        return H;
    }

    [[nodiscard]] const idx_t &sizeW() const
    {
        return W;
    }
};

template<typename T>
std::ostream &operator<<(std::ostream &os, const CSR<T> &A);

#endif//SOLE2022_CSR_HPP
