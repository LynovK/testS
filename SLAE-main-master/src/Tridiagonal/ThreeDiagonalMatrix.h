#ifndef SLAE_THREEDIAGONALMATRIX_H
#define SLAE_THREEDIAGONALMATRIX_H

#include <vector>

template<typename T>
struct elements {
    T a;
    T b;
    T c;
};

template<typename T>
class ThreeDiagonalMatrix {
private:
    std::vector <elements<T>> triple_;
public:
    ThreeDiagonalMatrix(const std::vector <elements<T>> &triple) noexcept: triple_(triple) {}

    const elements<T>& getString(unsigned int i) const  noexcept {
        return triple_[i];
    }
    unsigned int size() const{
        return triple_.size();
    }
};

#endif //SLAE_THREEDIAGONALMATRIX_H
