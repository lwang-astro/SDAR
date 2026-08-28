#pragma once

// g-function (time transformation for kick) method selection.
// Exactly one (or none) of these compile-time macros selects the method; each
// build compiles a single method:
//   AR_G_FUNC_BLOGH  : BLogH  - product of the innermost binary pair potentials
//   AR_G_FUNC_BTLOGH : BTLogH - tree-level product (inner pairs x outer node
//                       potentials); no auto switch (hyperbolic outer orbits
//                       are handled by the peri-center cap in processOuterNode)
// Runtime option (--g-func / g_func): 0 = standard LogH, 1 = the method of this
// build, 2 = auto switch between 1 and 0 (not available for BTLogH).
//
// NOTE: this header MUST be included (directly or transitively) before any use
// of the derived macros below; both AR/information.h and
// AR/symplectic_integrator.h include it first for this reason.
// History: NORM_BLOGH / MUL_ALL_POT / MAX_POT / ADD_INNER_POT variants were
// removed in the 2026-08-28 simplification; the last version implementing them
// is preserved at the git tag gfunc-archive.

#if (defined AR_G_FUNC_BLOGH) + (defined AR_G_FUNC_BTLOGH) > 1
#error "at most one AR_G_FUNC_* method macro can be defined for one build"
#endif

#if (defined AR_G_FUNC_BLOGH) || (defined AR_G_FUNC_BTLOGH)
#define AR_G_FUNC
#endif

