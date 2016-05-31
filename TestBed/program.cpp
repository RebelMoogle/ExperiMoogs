#ifndef NDEBUG
#include <vld.h>
#endif

#include <iostream>
#include <assert.h>
#include <vector>

template <typename T>
inline T mySum(const std::vector<T>& elementsToSum)
{
    T result = 0;
    for (T element : elementsToSum) {
        result += element;
    }
    return result;
}

template <typename T>
inline void DivideAll(const T factor, std::vector<T>& elementsToSum)
{
    if (factor == 0) {
       throw std::exception("Division by zero!!");
    }
    for (T& element : elementsToSum) {
        element /= factor;
    }
}

void main()
{
    std::vector<double> testContainer = { 1.0, 0.0, 2.0, 3.5 };
    assert(mySum(testContainer) == 6.5 && "Sum not correctly calculated.");

    try {
        DivideAll<double>(0.0, testContainer);
        std::cout << "Fail: No exception thrown!" << "\n";
    }
    catch (std::exception& e) {
        std::cout << "Successful: Exception thrown: " << e.what() << "\n";
    }
    catch (...) {
        std::cout << "Fail: unexpected exception thrown! " << "\n";
    }



}