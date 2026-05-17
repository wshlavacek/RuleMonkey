// Runtime pattern-parser unit test (issue #9 §1 step 2).
//
// parse_species_pattern(text, model) is pure — string + Model -> Pattern
// — so this test hand-builds a small Model (no XML, no engine) and
// drives the parser directly: the accepted grammar, model resolution
// (molecule types, components, states, symmetric components, bonds),
// and the full error surface that exact-species scope (design
// decision A) rejects.
//
// Molecule types in the fixture model:
//   A(b)        — one stateless binding component
//   B(a)        — one stateless binding component
//   L(r,r)      — two same-named (symmetric) stateless components
//   T(s)        — one stateful component, states {u,p}

#include "pattern_parser.hpp"

#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

int g_failures = 0;

void check(bool ok, const std::string& msg) {
  if (!ok) {
    std::fprintf(stderr, "FAIL: %s\n", msg.c_str());
    ++g_failures;
  }
}

rulemonkey::MoleculeType make_type(const std::string& name,
                                   const std::vector<rulemonkey::MoleculeTypeComponent>& comps) {
  rulemonkey::MoleculeType mt;
  mt.id = name;
  mt.name = name;
  mt.components = comps;
  return mt;
}

rulemonkey::MoleculeTypeComponent comp(const std::string& name,
                                       const std::vector<std::string>& states = {}) {
  rulemonkey::MoleculeTypeComponent mc;
  mc.name = name;
  mc.allowed_states = states;
  return mc;
}

rulemonkey::Model make_model() {
  rulemonkey::Model m;
  m.molecule_types.push_back(make_type("A", {comp("b")}));
  m.molecule_types.push_back(make_type("B", {comp("a")}));
  m.molecule_types.push_back(make_type("L", {comp("r"), comp("r")}));
  m.molecule_types.push_back(make_type("T", {comp("s", {"u", "p"})}));
  for (int i = 0; i < static_cast<int>(m.molecule_types.size()); ++i)
    m.molecule_type_index[m.molecule_types[i].name] = i;
  return m;
}

// Assert that `text` parses cleanly; return the Pattern for further checks.
rulemonkey::Pattern accept(const rulemonkey::Model& model, const std::string& text) {
  try {
    return rulemonkey::parse_species_pattern(text, model);
  } catch (const std::exception& e) {
    std::fprintf(stderr, "FAIL: expected '%s' to parse, but it threw: %s\n", text.c_str(),
                 e.what());
    ++g_failures;
    return rulemonkey::Pattern{};
  }
}

// Assert that `text` is rejected with a std::runtime_error.
void reject(const rulemonkey::Model& model, const std::string& text, const std::string& why) {
  try {
    rulemonkey::parse_species_pattern(text, model);
  } catch (const std::runtime_error&) {
    return; // expected: the parser rejected the malformed pattern
  } catch (const std::exception& e) {
    std::fprintf(stderr, "FAIL: '%s' threw a non-runtime_error: %s\n", text.c_str(), e.what());
    ++g_failures;
    return;
  }
  std::fprintf(stderr, "FAIL: expected '%s' to be rejected (%s), but it parsed\n", text.c_str(),
               why.c_str());
  ++g_failures;
}

void test_single_molecule(const rulemonkey::Model& model) {
  const auto pat = accept(model, "A(b)");
  check(pat.molecules.size() == 1, "A(b): one molecule");
  if (pat.molecules.size() == 1) {
    const auto& mol = pat.molecules[0];
    check(mol.type_index == 0, "A(b): A resolves to type index 0");
    check(mol.components.size() == 1, "A(b): one component");
    if (mol.components.size() == 1) {
      check(mol.components[0].comp_type_index == 0, "A(b): component b is type index 0");
      check(mol.components[0].bond_constraint == rulemonkey::BondConstraint::Free,
            "A(b): unbonded component is Free");
    }
  }
  check(pat.bonds.empty(), "A(b): no bonds");
}

void test_two_molecule_bond(const rulemonkey::Model& model) {
  const auto pat = accept(model, "A(b!1).B(a!1)");
  check(pat.molecules.size() == 2, "A(b!1).B(a!1): two molecules");
  check(pat.bonds.size() == 1, "A(b!1).B(a!1): exactly one bond");
  if (pat.bonds.size() == 1) {
    // Flat indices: A.b == 0, B.a == 1.
    const auto& bnd = pat.bonds[0];
    const bool wired = (bnd.comp_flat_a == 0 && bnd.comp_flat_b == 1) ||
                       (bnd.comp_flat_a == 1 && bnd.comp_flat_b == 0);
    check(wired, "A(b!1).B(a!1): bond connects the two flat components");
  }
  if (pat.molecules.size() == 2) {
    check(pat.molecules[0].components[0].bond_constraint == rulemonkey::BondConstraint::BoundTo,
          "A(b!1).B(a!1): bonded component is BoundTo");
    check(pat.molecules[0].components[0].bond_label == 0,
          "A(b!1).B(a!1): bonded component carries bond label 0");
    check(pat.molecules[1].components[0].bond_label == 0,
          "A(b!1).B(a!1): both endpoints share bond label 0");
  }
}

void test_symmetric_components(const rulemonkey::Model& model) {
  const auto pat = accept(model, "L(r,r)");
  check(pat.molecules.size() == 1, "L(r,r): one molecule");
  if (!pat.molecules.empty()) {
    const auto& mol = pat.molecules[0];
    check(mol.components.size() == 2, "L(r,r): two components");
    if (mol.components.size() == 2) {
      // The two same-named components must bind distinct type indices.
      check(mol.components[0].comp_type_index == 0 && mol.components[1].comp_type_index == 1,
            "L(r,r): symmetric components map to distinct type indices by occurrence");
    }
  }

  // An intramolecular bond between the two symmetric components.
  const auto ring = accept(model, "L(r!1,r!1)");
  check(ring.bonds.size() == 1, "L(r!1,r!1): one intramolecular bond");
}

void test_stateful_component(const rulemonkey::Model& model) {
  const auto pat = accept(model, "T(s~p)");
  check(pat.molecules.size() == 1, "T(s~p): one molecule");
  if (!pat.molecules.empty() && !pat.molecules[0].components.empty()) {
    const auto& c = pat.molecules[0].components[0];
    check(c.required_state == "p", "T(s~p): required_state captured");
    check(c.required_state_index == 1, "T(s~p): state 'p' resolves to index 1");
  }
}

void test_multi_molecule_chain(const rulemonkey::Model& model) {
  // L bridges A and B through its two r components — a connected complex.
  const auto pat = accept(model, "L(r!1,r!2).A(b!1).B(a!2)");
  check(pat.molecules.size() == 3, "L(r!1,r!2).A(b!1).B(a!2): three molecules");
  check(pat.bonds.size() == 2, "L(r!1,r!2).A(b!1).B(a!2): two bonds");
}

void test_whitespace_tolerance(const rulemonkey::Model& model) {
  const auto pat = accept(model, "  A( b ! 1 ) . B( a ! 1 )  ");
  check(pat.molecules.size() == 2 && pat.bonds.size() == 1,
        "whitespace around tokens is tolerated");
}

void test_error_surface(const rulemonkey::Model& model) {
  reject(model, "", "empty string");
  reject(model, "Z(x)", "unknown molecule type");
  reject(model, "A(q)", "unknown component");
  reject(model, "A()", "component b not specified (under-specified)");
  reject(model, "A", "bare molecule name with components unspecified");
  reject(model, "A(b,b)", "component b listed more times than A has it");
  reject(model, "L(r)", "only one of L's two r components specified");
  reject(model, "T(s)", "stateful component needs a ~state");
  reject(model, "T(s~z)", "z is not an allowed state of s");
  reject(model, "A(b~u)", "stateless component given a ~state");
  reject(model, "A(b!+)", "bond wildcard !+ rejected");
  reject(model, "A(b!?)", "bond wildcard !? rejected");
  reject(model, "A(b!1)", "dangling bond label (appears once)");
  reject(model, "L(r!1,r!1).A(b!1)", "bond label on three components");
  reject(model, "A(b!1).B(a!1).L(r,r)", "disconnected complex (L not bonded)");
  reject(model, "A(b!x)", "non-numeric bond label");
  reject(model, "A(b(", "unbalanced parenthesis");
}

} // namespace

int main() {
  const rulemonkey::Model model = make_model();

  try {
    test_single_molecule(model);
    test_two_molecule_bond(model);
    test_symmetric_components(model);
    test_stateful_component(model);
    test_multi_molecule_chain(model);
    test_whitespace_tolerance(model);
    test_error_surface(model);
  } catch (const std::exception& e) {
    std::fprintf(stderr, "EXCEPTION: %s\n", e.what());
    return 2;
  }

  if (g_failures > 0) {
    std::fprintf(stderr, "\n%d assertion(s) failed\n", g_failures);
    return 1;
  }
  std::fprintf(stderr, "OK: pattern-parser assertions all passed\n");
  return 0;
}
