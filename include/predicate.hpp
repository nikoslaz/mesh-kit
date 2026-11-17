#ifndef PREDICATE_HPP
#define PREDICATE_HPP

extern "C" {
    void   exactinit();
    double orient2d(const double* pa, const double* pb, const double* pc);
    double incircle(const double* pa, const double* pb, const double* pc, const double* pd);
}

#endif