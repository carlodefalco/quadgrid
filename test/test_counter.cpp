#include "counter.h"
#include <iostream>
#include <iterator>
#include <algorithm>

int main() {
    // Test range<int> from 0 to 5
    range<int> r(0, 5);
    std::cout << "Testing range<int> iterator with std::next:\n";
    //auto it = r.begin();
    for (auto it = r.begin(); it != r.end();  ++it) {
        std::cout << "Value: " << *it << ", std::next: ";
        if (std::next(it) != r.end())
            std::cout << *std::next(it);
        else
            std::cout << "(end)";
        std::cout << "\n";
    }

    // Check usability in std::for_each
    std::cout << "\nTesting usability in std::for_each:\n";
    std::for_each(r.begin(), r.end(), [](int v) {
        std::cout << v << " ";
    });
    std::cout << "\n";
    return 0;
}
