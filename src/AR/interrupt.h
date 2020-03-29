#pragma once

//! binary interrupt status
enum class InterruptStatus{change=1, merge=2, none=0};

//! Binary interrupt recoreder
template <class Tparticle>
struct BinaryInterrupt {
    COMM::BinaryTree<Tparticle>* adr;
    Float time_now;
    Float time_end;
    InterruptStatus status;
};
