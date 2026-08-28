#pragma once

// g-function (time transformation for kick) method selection.
// Exactly one (or none) of these compile-time macros selects the method; each
// build compiles a single method, matching the MAX_POT convention:
//   AR_G_FUNC_BLOGH         : BLogH      - product of the innermost binary pair potentials
//   AR_G_FUNC_NORM_BLOGH    : normalized BLogH - (product of innermost pair potentials)^(1/n_bin)
//   AR_G_FUNC_MUL_ALL_POT   : product over ALL pair potentials (auto ds unsupported, needs --s)
//   AR_G_FUNC_BTLOGH        : tree-level product (inner pairs x outer node potentials); no auto switch
//   AR_G_FUNC_MAX_POT       : maximum of the innermost binary pair potentials (smoothed for TTL)
//   AR_G_FUNC_ADD_INNER_POT : sum of the innermost binary pair potentials
//     (implemented indirectly: base sum accumulation + skipping cross pairs via calc_gt_cross)
// Runtime option (--g-func / g_func): 0 = standard LogH, 1 = the method of this
// build, 2 = auto switch between 1 and 0 (not available for BTLogH).
//
// NOTE: this header MUST be included (directly or transitively) before any use
// of the derived macros below; both AR/information.h and
// AR/symplectic_integrator.h include it first for this reason.

#if (defined AR_G_FUNC_BLOGH) + (defined AR_G_FUNC_NORM_BLOGH) + (defined AR_G_FUNC_MUL_ALL_POT) \
    + (defined AR_G_FUNC_BTLOGH) + (defined AR_G_FUNC_MAX_POT) + (defined AR_G_FUNC_ADD_INNER_POT) > 1
#error "at most one AR_G_FUNC_* method macro can be defined for one build"
#endif

#if (defined AR_G_FUNC_BLOGH) || (defined AR_G_FUNC_NORM_BLOGH) || (defined AR_G_FUNC_MUL_ALL_POT) \
    || (defined AR_G_FUNC_BTLOGH) || (defined AR_G_FUNC_MAX_POT) || (defined AR_G_FUNC_ADD_INNER_POT)
#define AR_G_FUNC
#endif

// family macro: the four product-form methods share the GtKickInv product
// structure (nbin / mul_pot_no_pow) and the pair-accumulation branch.
// MAX_POT has its own structure; ADD_INNER_POT uses the base sum structure.
#if (defined AR_G_FUNC_BLOGH) || (defined AR_G_FUNC_NORM_BLOGH) || (defined AR_G_FUNC_MUL_ALL_POT) || (defined AR_G_FUNC_BTLOGH)
#define AR_G_FUNC_MUL_POT_FAMILY
#endif

#ifdef AR_G_FUNC_MUL_POT
#error "AR_G_FUNC_MUL_POT is split into AR_G_FUNC_BLOGH / AR_G_FUNC_NORM_BLOGH / AR_G_FUNC_MUL_ALL_POT / AR_G_FUNC_BTLOGH (one per build)"
#endif
#ifdef AR_G_FUNC_ADD_POT
#error "AR_G_FUNC_ADD_POT is renamed to AR_G_FUNC_ADD_INNER_POT"
#endif
