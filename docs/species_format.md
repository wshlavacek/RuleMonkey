# `.species` output format and species enumeration

RuleMonkey can enumerate the **live species** of a paused simulation —
group the molecules in the agent pool into connected complexes, collapse
graph-isomorphic complexes into one species, and report each as a
BNG-format pattern line with its instance count. The result is written
to a `.species` file readable by BioNetGen's `BNG2.pl` (its
`readNFspecies` reader). This document is the authoritative reference
for the format and the enumeration semantics.

This is RuleMonkey's answer to [issue #9 §2][i9]. It is the analogue of
NFsim's `-ss` species dump, but RuleMonkey deduplicates by *true* graph
isomorphism (see [Canonical species identity](#canonical-species-identity)).

## File shape

One or more `#` comment lines, then one data line per species:

```
# RuleMonkey generated species list
# <N> species, <M> complexes
<canonical pattern><SPACE><SPACE><count>
<canonical pattern><SPACE><SPACE><count>
...
```

- Every line beginning with `#` is a comment. `BNG2.pl readNFspecies`
  strips `#` to end-of-line, so the header is informational only.
- Each data line is a **canonical BNGL species pattern** followed by an
  **integer count** — the number of complex instances of that species
  currently in the pool. The separator written is two spaces; any run of
  whitespace is accepted on read.
- Data lines are sorted lexicographically by the species pattern, so the
  file is deterministic for a given pool state.
- The count is always a non-negative integer (network-free particles
  have unit population). Only **live** species appear — there are no
  zero-count rows for extinct seed species, matching NFsim `-ss`.

Example (`A(a) + A(a) <-> A(a!1).A(a!1)`, mid-run):

```
# RuleMonkey generated species list
# 2 species, 754 complexes
A(a!1).A(a!1)  246
A(a)  508
```

Molecules are conserved: `508·1 + 246·2 = 1000`, the seed `A_tot`.

## Canonical species identity

Two complexes are the **same species** iff their molecule graphs are
isomorphic (same molecule types, component names and states, and bond
topology — independent of the order molecules and components happen to
occupy in the pool). RuleMonkey assigns each complex a *canonical label*
— a normalized BNGL string that is byte-identical for all isomorphic
complexes and distinct for non-isomorphic ones — and uses it as both the
deduplication key and the emitted pattern text. Symmetric complexes
(rings, homo-oligomers) are resolved by individualization-refinement, so
the count for, say, a symmetric homodimer is exact rather than inflated
by molecule-ordering accidents.

The canonical string is a valid BNG-format pattern. RuleMonkey's
canonical format need not byte-match NFsim's internal label: `BNG2.pl`
re-canonicalizes every species it reads, and `readNFspecies` sums the
counts of any lines that turn out to be isomorphic. A parity check
against NFsim `-ss` output must therefore compare the two files as a
*multiset of canonicalized species*, not as text.

## How to produce a `.species` file

### `rm_driver` CLI

```
rm_driver <model.xml> <t_end> <n_steps> [seed] --species <path>
```

`--species` writes the final-state census to `<path>` after the run.
Passing it forces session mode (the pool must outlive the run), so it
composes with `--save-state` / `--load-state` / `--t-start`.

### Embedding API

`RuleMonkeySimulator` exposes the census to in-process callers on an
active session:

- `std::vector<SpeciesRow> enumerate_species() const` — the sorted
  `(species, count)` rows. Each `SpeciesRow` (see `types.hpp`) is a
  canonical pattern string and an instance count.
- `void write_species_file(const std::string& path) const` — writes the
  rows to `path` in the format above.
- `long species_count(const std::string& canonical_species) const` —
  the count for a single species, looked up by its canonical string.
- `long total_complex_count() const` — the total number of live
  complexes in the pool (the sum of every row's count).

All four require a live session (`initialize()` / `simulate()` first)
and throw `std::runtime_error` otherwise. Enumeration is a one-shot pool
walk, intended to be called while the simulation is paused (between
`simulate()` segments or after a run) — not per SSA event.

### Looking up a single species — `species_count`

`species_count` answers "how many of *this* species are live" without
the caller scanning the whole `enumerate_species()` vector. Its key is a
**canonical string RuleMonkey itself emitted** — a `SpeciesRow::species`
value or a data-line pattern from a written `.species` file:

```cpp
for (const auto& row : sim.enumerate_species())
  assert(sim.species_count(row.species) == row.count);
```

RuleMonkey has **no runtime BNGL-pattern parser** (that is [issue #9
§1][i9], separate work). `species_count` therefore does *not*
canonicalize its argument: a hand-written pattern that is not
byte-identical to a label RuleMonkey would emit returns `0`, even when
it denotes the same species. An embedder that needs pattern-keyed lookup
(`get_species_count("A(b!1).B(a!1)")`-style, NFsim-parity) must
canonicalize the pattern on its own side and pass the canonical form.

`species_count` is a *batch query*: each call is internally a full pool
walk, the same cost as `enumerate_species()`. To read many species,
call `enumerate_species()` once and index its rows rather than calling
`species_count` per species. `total_complex_count()` is cheaper than
either — it counts live complexes without canonicalizing them.

[i9]: https://github.com/richardposner/RuleMonkey/issues/9
