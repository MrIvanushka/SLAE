//
// Created by ivanb on 11.02.2022.
//

#ifndef SLAE_UNDECLAREDELEMENTCZLLEXCEPTION_HPP
#define SLAE_UNDECLAREDELEMENTCZLLEXCEPTION_HPP

#include <string>

#include "SlaeBaseException.hpp"

namespace Slae {
    class UndeclaredElementCallException : public SlaeBaseExceptionCpp {
    public:
        UndeclaredElementCallException(std::string message) : SlaeBaseExceptionCpp(message)
        { }
    };
}

#endif //SLAE_UNDECLAREDELEMENTCZLLEXCEPTION_HPP
