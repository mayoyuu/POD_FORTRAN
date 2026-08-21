# ADS split via DA composition

## Goal

Make the Fortran ADS split follow the reference C++ implementation and stop
calling DACE `translateVariable`. No DACE core source or library is modified.
The public Fortran interface remains
`patch_split(parent, direction, left, right)`.

## Split algorithm

For each split, construct an identity DA map with one entry per active DACE
variable. For direction `dir`, replace only that entry:

```text
left_map(dir)  = 0.5 * DA(dir) - 0.5
right_map(dir) = 0.5 * DA(dir) + 0.5
left.da_vec    = parent.da_vec.eval(left_map)
right.da_vec   = parent.da_vec.eval(right_map)
```

This is the same polynomial composition used by C++ `Patch::split`. The Patch
output-component count comes from `parent%da_vec`; the identity-map length
comes from `dace_max_variables()`. They remain independent.

## Fortran ownership

Add a no-temporary DA-vector composition subroutine in
`pod_dace_classes`. It passes the map handles to the existing
`fdace_eval_da_vec` wrapper and writes directly into preallocated output DA
handles. `patch_split` uses this subroutine for both children, avoiding
finalizable function results and raw handle replacement.

Identity variables are created without overloaded-expression temporaries. A
separate DA temporary forms `0.5*DA(dir) +/- 0.5`, so input and output dummy
arguments never alias.

`patch_init` transfers the allocated DA element array with `move_alloc`
instead of allocating destination handles and overwriting them. This preserves
the existing move semantics while removing the leaked handles.

The obsolete ADS-facing `da_translate_variable` wrapper may remain available
for unrelated callers, but ADS no longer imports or calls it.

## Error handling

`patch_split` rejects directions outside `1:dace_max_variables()`. It also
requires an allocated parent DA vector and preserves the parent's output
component count in both children.

## Tests

Extend `test_ads_core_dynamic_dimensions` as the regression test:

- initialize DACE with order 2 and seven variables;
- split a three-component Patch whose first component depends on variable 7;
- verify split direction 7 is selected;
- verify both children retain three output components;
- evaluate left and right children at representative global-domain points and
  compare with the unsplit polynomial;
- verify Patch cleanup does not leave additional active DA handles.

The existing test currently fails inside `translateVariable`, so it supplies
the RED state. Per the user's instruction, Codex will edit but will not compile
or run; the user performs RED/GREEN verification locally.

## Non-goals

- No modification or rebuild of DACE core.
- No change to ADS normalization or physical uncertainty scaling.
- No change to HFEM force-model behavior in this step.
- No full-history replay during a single split; each child is composed from
  its current parent, matching C++ `Patch::split`.
