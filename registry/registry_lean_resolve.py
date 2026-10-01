#!/usr/bin/env python3
"""Resolve every registry node's `lean:` field against the Lean sources.

WHY THIS EXISTS.  `registry_validate.py` never reads the `lean` field -- grep it: the only
occurrence outside the argument parser is line 318, inside report formatting.  So a node can
carry `trust: lean-verified` and `lean: <name that does not exist>` and the validator prints
`OK`, exit 0.  Measured 2026-10-01: promoting `lemma-T-interior-concavity` to `lean-verified`
with `lean: TworowD4Kernel.LemmaT.Setup.thm_main_does_not_exist` passes `registry_validate.py`
cleanly.  A validator's name is a claim, and nothing grades the name.

WHAT THIS GRADES.  For each node with a `lean` field: does that fully-qualified name occur as a
declaration in some .lean file under the given roots?  Declaration names are reconstructed by
tracking `namespace`/`end` nesting, and a declaration written `Foo.bar` inside `namespace Baz`
resolves as `Baz.Foo.bar`.

WHAT THIS DOES NOT GRADE.  (a) Whether the declaration is sorry-free -- grep the build for
`sorryAx`.  (b) Whether its STATEMENT is the paper's lemma -- no tool can do that; read it.
(c) Whether `trust: lean-verified` is deserved; it only checks the pointer resolves.

Usage:  python3 code/registry_lean_resolve.py [registry.json ...] [--lean-root DIR ...]
Exit 0 if every `lean` pointer resolves, 1 otherwise.
"""

import argparse
import json
import os
import re
import sys

DECL = re.compile(
    r"^\s*(?:@\[[^\]]*\]\s*)?(?:private\s+|protected\s+|noncomputable\s+|partial\s+|unsafe\s+)*"
    r"(?:theorem|lemma|def|abbrev|instance|structure|inductive|class|opaque|axiom)\s+"
    r"(\{[^}]*\}\s*)?([^\s({\[:]+)"
)
NS = re.compile(r"^\s*namespace\s+(\S+)")
END = re.compile(r"^\s*end\b\s*(\S*)")


def declarations(roots):
    """Set of fully-qualified declaration names found under roots."""
    names = set()
    for root in roots:
        for dirpath, dirnames, filenames in os.walk(root):
            dirnames[:] = [d for d in dirnames if d != ".lake"]
            for fn in filenames:
                if not fn.endswith(".lean"):
                    continue
                stack = []
                with open(os.path.join(dirpath, fn), errors="replace") as fh:
                    for line in fh:
                        m = NS.match(line)
                        if m:
                            stack.append(m.group(1))
                            continue
                        m = END.match(line)
                        if m:
                            # `end Foo.Bar` closes exactly the namespace frames it names.  A
                            # named `end` that does NOT match an open namespace suffix is a
                            # SECTION end and must pop nothing -- getting this wrong silently
                            # drops the namespace prefix from every later declaration in the
                            # file, which is how this checker first reported 35 phantom
                            # dangling pointers (all five hand-spot-checks "confirmed" it,
                            # because the confirming greps were broken the same way).
                            tgt = m.group(1)
                            if tgt and stack:
                                for k in range(1, len(stack) + 1):
                                    if ".".join(stack[-k:]) == tgt:
                                        del stack[-k:]
                                        break
                            continue
                        m = DECL.match(line)
                        if m:
                            short = m.group(2)
                            names.add(".".join(stack + [short]) if stack else short)
    return names


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("registries", nargs="*", default=None)
    ap.add_argument("--lean-root", action="append", default=None)
    args = ap.parse_args()

    regs = args.registries or sorted(
        os.path.join("proofs/registry", f)
        for f in os.listdir("proofs/registry") if f.endswith(".json")
    )
    roots = args.lean_root or ["lean"]
    names = declarations(roots)
    print(f"{len(names)} declarations indexed under {', '.join(roots)}")

    problems = []
    checked = 0
    for reg in regs:
        with open(reg) as fh:
            tree = json.load(fh).get("tree")
        if not isinstance(tree, dict):
            continue
        stack = [(tree, reg + ":" + str(tree.get("id")))]
        while stack:
            node, path = stack.pop()
            decl = node.get("lean")
            if decl:
                # Several existing nodes hold a comma-separated LIST in one `lean` field.
                for one in [d.strip() for d in decl.split(",") if d.strip()]:
                    checked += 1
                    if one not in names:
                        problems.append(f"{path}: lean '{one}' does not resolve")
            if node.get("trust") == "lean-verified" and not decl:
                problems.append(f"{path}: trust lean-verified with no `lean` field")
            for c in node.get("children", []):
                stack.append((c, path + "/" + str(c.get("id"))))

    print(f"{checked} `lean` pointer(s) checked across {len(regs)} registry file(s)")
    if problems:
        print(f"\n{len(problems)} problem(s):")
        for p in problems:
            print("  - " + p)
        return 1
    print("OK: every `lean` pointer resolves")
    return 0


if __name__ == "__main__":
    sys.exit(main())
