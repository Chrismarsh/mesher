#include "mesh.h"
#include <iostream>
#include <stdexcept>
#include <vector>

void require(bool condition, const char* message)
{
    if (!condition) throw std::runtime_error(message);
}

int main()
{
    try {
        std::vector<Criteria::Quality> values;
        for (double size : {.01, .5, 1., 2., 100.}) {
            for (double sine : {.001, .125, .5, .75}) {
                Criteria::Quality q;
                q._area = 6250000 * size;
                q._size = size;
                q._sine = sine;
                values.push_back(q);
                q._tolerance.emplace_back(50, 30);
                q._category_tol.emplace_back(.4, .6);
                values.push_back(q);
            }
        }
        auto equivalent = [](const Criteria::Quality& a, const Criteria::Quality& b) {
            return !(a < b) && !(b < a);
        };
        for (auto a : values) {
            require(!(a < a), "Queue order must be irreflexive");
            for (auto b : values) {
                require(!(a < b && b < a), "Queue order must be asymmetric");
                if (a._size > 1 && b._size <= 1)
                    require(a < b, "Oversized faces must be refined first");
                if (a._size > 1 && b._size > 1 && a._size > b._size)
                    require(a < b, "Largest area violation must be refined first");
                if (a._size <= 1 && b._size <= 1 && a._sine < b._sine)
                    require(a < b, "Smallest angle must be refined first below the area limit");
                for (auto c : values) {
                    require(!(a < b && b < c) || a < c, "Queue order must be transitive");
                    require(!(equivalent(a, b) && equivalent(b, c)) || equivalent(a, c),
                            "Queue equivalence must be transitive");
                }
            }
        }
        std::cout << "Mesh criteria ordering regressions passed" << std::endl;
    } catch (const std::exception& error) {
        std::cerr << error.what() << std::endl;
        return 1;
    }
}
