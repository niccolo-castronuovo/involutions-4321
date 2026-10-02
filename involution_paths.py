#!/usr/bin/env python3
"""
involution_paths.py
===================

Companion code for the paper
  "Involutions avoiding 4321 and another pattern of length four"
  (M. Barnabei, N. Castronuovo, M. Silimbani).

For a set T of patterns (avoided together with 4321) and every n <= N, the
script

  1. generates all connected involutions in I_n(4321, T);
  2. computes the associated (unlabelled) Motzkin paths via Biane's bijection;
  3. finds the minimal words over {U, H, D} of length <= L that occur as a
     (not necessarily contiguous) subsequence of none of these paths
     (the procedure guess[T] of the paper);
  4. checks, for every n <= N_check, whether the connected paths avoiding
     these words are exactly the paths of the connected involutions in the
     class (same sets, hence same numbers);
  5. optionally checks whether the set of connected paths is closed under
     taking connected subsequences, a necessary and sufficient condition for
     the existence of ANY characterization by forbidden subsequences;
  6. optionally prints the joint distribution of length and descents over the
     whole class I_n(4321, T), to be compared with the generating functions.

Main design choice
------------------
Involutions avoiding 4321 are in bijection with UNLABELLED Motzkin paths
(every down step closes the earliest open transposition). So instead of
generating all involutions and discarding those containing 4321, we generate
Motzkin paths and invert the bijection: this is cheaper.

Usage examples
--------------
    python involution_paths.py 2143
    python involution_paths.py 2341 4123 -N 12 -L 4 --check-N 14 --closure
    python involution_paths.py 1324 --distribution 8
"""

from __future__ import annotations

import argparse
import sys
from collections import deque
from functools import lru_cache
from itertools import product
from typing import Dict, Iterable, List, Sequence, Set, Tuple

Perm = Tuple[int, ...]


# ----------------------------------------------------------------------
# 1. Motzkin paths and Biane's bijection (unlabelled case)
# ----------------------------------------------------------------------

@lru_cache(maxsize=None)
def motzkin_paths(n: int) -> Tuple[str, ...]:
    """All Motzkin paths of length n, as words over {U, H, D}."""
    out: List[str] = []

    def rec(prefix: List[str], height: int, remaining: int) -> None:
        if remaining == 0:
            if height == 0:
                out.append("".join(prefix))
            return
        if height > remaining:          # cannot come back to 0
            return
        prefix.append("U"); rec(prefix, height + 1, remaining - 1); prefix.pop()
        prefix.append("H"); rec(prefix, height, remaining - 1); prefix.pop()
        if height > 0:
            prefix.append("D"); rec(prefix, height - 1, remaining - 1); prefix.pop()

    rec([], 0, n)
    return tuple(out)


def connected_motzkin_paths(n: int) -> Tuple[str, ...]:
    """Connected Motzkin paths of length n: H (n = 1) or U m D."""
    if n == 1:
        return ("H",)
    if n < 2:
        return ()
    return tuple("U" + m + "D" for m in motzkin_paths(n - 2))


def biane_inverse(path: str) -> Perm:
    """
    Involution avoiding 4321 associated with an unlabelled Motzkin path:
    H = fixed point, U = opening, D = closing of the EARLIEST open
    transposition (all labels equal to 1).
    """
    n = len(path)
    perm = [0] * n
    open_positions: deque = deque()
    for i, step in enumerate(path):
        if step == "H":
            perm[i] = i + 1
        elif step == "U":
            open_positions.append(i)
        else:
            j = open_positions.popleft()
            perm[i], perm[j] = j + 1, i + 1
    return tuple(perm)


def fmt_perm(perm: Sequence[int]) -> str:
    """One-line notation; entries separated by spaces when some entry exceeds 9."""
    sep = "" if len(perm) < 10 else " "
    return sep.join(map(str, perm))


def descents(perm: Sequence[int]) -> int:
    return sum(perm[i] > perm[i + 1] for i in range(len(perm) - 1))


# ----------------------------------------------------------------------
# 2. Pattern containment (backtracking with value windows)
# ----------------------------------------------------------------------

def _pattern_plan(pattern: Sequence[int]) -> List[Tuple[int, int]]:
    """
    For each position t of the pattern, the indices s < t of the earlier
    entries that are the closest smaller and the closest larger values
    (or -1). They give the window in which the t-th chosen value must lie.
    """
    plan = []
    for t, v in enumerate(pattern):
        lower = [s for s in range(t) if pattern[s] < v]
        upper = [s for s in range(t) if pattern[s] > v]
        lo = max(lower, key=lambda s: pattern[s]) if lower else -1
        hi = min(upper, key=lambda s: pattern[s]) if upper else -1
        plan.append((lo, hi))
    return plan


def make_containment_test(pattern: Sequence[int]):
    """Return a fast function perm -> bool telling whether perm contains pattern."""
    k = len(pattern)
    plan = _pattern_plan(pattern)

    def contains(perm: Sequence[int]) -> bool:
        n = len(perm)
        if k > n:
            return False
        chosen = [0] * k

        def rec(t: int, start: int) -> bool:
            lo, hi = plan[t]
            low_val = chosen[lo] if lo >= 0 else 0
            high_val = chosen[hi] if hi >= 0 else n + 1
            for p in range(start, n - (k - t) + 1):
                v = perm[p]
                if low_val < v < high_val:
                    chosen[t] = v
                    if t == k - 1 or rec(t + 1, p + 1):
                        return True
            return False

        return rec(0, 0)

    return contains


def make_class_test(patterns: Iterable[Sequence[int]]):
    tests = [make_containment_test(p) for p in patterns]
    return lambda perm: not any(test(perm) for test in tests)


# ----------------------------------------------------------------------
# 3. Data: connected paths of the class
# ----------------------------------------------------------------------

def class_connected_paths(n: int, in_class) -> Set[str]:
    """Paths beta(pi) of the connected involutions pi in I_n(4321, T)."""
    return {w for w in connected_motzkin_paths(n) if in_class(biane_inverse(w))}


# ----------------------------------------------------------------------
# 4. Minimal forbidden subsequences (procedure guess[T])
# ----------------------------------------------------------------------

def subsequences_up_to(word: str, max_len: int) -> Set[str]:
    """All distinct subsequences of `word` of length <= max_len."""
    subs = {""}
    for c in word:
        subs |= {s + c for s in subs if len(s) < max_len}
    return subs


def is_subsequence(small: str, big: str) -> bool:
    it = iter(big)
    return all(c in it for c in small)


def minimal_forbidden_words(paths: Iterable[str], max_len: int,
                            alphabet: str = "UHD") -> List[str]:
    """
    Minimal words of length <= max_len that are a subsequence of no path.
    Instead of testing every word against every path, we collect once all
    subsequences that DO occur, which is much faster.
    """
    occurring: Set[str] = set()
    for p in set(paths):
        occurring |= subsequences_up_to(p, max_len)
    minimal: List[str] = []
    for length in range(1, max_len + 1):
        for letters in product(alphabet, repeat=length):
            w = "".join(letters)
            if w in occurring:
                continue
            if any(is_subsequence(m, w) for m in minimal):
                continue
            minimal.append(w)
    return minimal


# ----------------------------------------------------------------------
# 5. Checks
# ----------------------------------------------------------------------

def check_characterization(words: List[str], in_class, n_max: int,
                           data: Dict[int, Set[str]]) -> bool:
    """
    For every n <= n_max compare the connected paths of the class with the
    connected paths avoiding all `words`. Prints a table; returns True iff
    the two sets coincide for every n.
    """
    print(f"\n--- Check: class paths  vs  paths avoiding the words (n <= {n_max}) ---")
    print(f"{'n':>3} {'#class':>9} {'#avoiders':>10}  result")
    all_ok = True
    for n in range(1, n_max + 1):
        A = data[n] if n in data else class_connected_paths(n, in_class)
        B = {w for w in connected_motzkin_paths(n)
             if not any(is_subsequence(x, w) for x in words)}
        ok = A == B
        all_ok &= ok
        msg = "OK" if ok else "MISMATCH"
        print(f"{n:>3} {len(A):>9} {len(B):>10}  {msg}")
        if not ok:
            extra = sorted(B - A)[:3]
            missing = sorted(A - B)[:3]
            if extra:
                print("      avoid the words but not in the class, e.g.:",
                      ", ".join(f"{w} -> {fmt_perm(biane_inverse(w))}"
                                for w in extra))
            if missing:
                print("      in the class but containing a word, e.g.:", ", ".join(missing))
    return all_ok


def check_closure(data: Dict[int, Set[str]], n_max: int):
    """
    The connected paths of the class can be characterized by forbidden
    subsequences iff they are closed under taking connected subsequences.
    Look for a connected path q NOT in the class which is a subsequence of
    a path p IN the class (n <= n_max). Returns (q, p) or None.
    """
    for n in range(1, n_max):
        for q in connected_motzkin_paths(n):
            if q in data[n]:
                continue
            for m in range(n + 1, n_max + 1):
                for p in data[m]:
                    if is_subsequence(q, p):
                        return q, p
    return None


def descent_distribution(n: int, in_class) -> Dict[int, int]:
    """Number of involutions in I_n(4321, T) with k descents, for every k."""
    dist: Dict[int, int] = {}
    for w in motzkin_paths(n):
        perm = biane_inverse(w)
        if in_class(perm):
            d = descents(perm)
            dist[d] = dist.get(d, 0) + 1
    return dist


def format_poly(dist: Dict[int, int]) -> str:
    terms = []
    for k in sorted(dist):
        c = dist[k]
        if k == 0:
            terms.append(f"{c}")
        elif k == 1:
            terms.append(f"{c}*y" if c != 1 else "y")
        else:
            terms.append(f"{c}*y^{k}" if c != 1 else f"y^{k}")
    return " + ".join(terms) if terms else "0"


# ----------------------------------------------------------------------
# 6. Driver
# ----------------------------------------------------------------------

def parse_pattern(s: str) -> Perm:
    p = tuple(int(c) for c in s)
    if sorted(p) != list(range(1, len(p) + 1)):
        raise argparse.ArgumentTypeError(f"{s} is not a permutation pattern")
    return p


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("patterns", nargs="+", type=parse_pattern,
                    help="patterns avoided together with 4321, e.g. 2341 4123")
    ap.add_argument("-N", type=int, default=12,
                    help="max length of the generated involutions (default 12)")
    ap.add_argument("-L", type=int, default=4,
                    help="max length of the forbidden words (default 4)")
    ap.add_argument("--check-N", type=int, default=None,
                    help="max length for the final check (default: N)")
    ap.add_argument("--closure", action="store_true",
                    help="test closure under connected subsequences")
    ap.add_argument("--distribution", type=int, default=0, metavar="M",
                    help="print the descent distribution over I_n(4321,T), n <= M")
    args = ap.parse_args(argv)

    patterns = [(4, 3, 2, 1)] + list(args.patterns)
    in_class = make_class_test(patterns)
    label = ", ".join("".join(map(str, p)) for p in args.patterns)
    check_n = args.check_N or args.N

    print("=" * 70)
    print(f"Connected involutions in I_n(4321, {label}),  n = 1..{args.N}")
    print("=" * 70)

    data = {n: class_connected_paths(n, in_class) for n in range(1, max(args.N, check_n) + 1)}
    print("\nNumber of connected involutions:")
    print("  n :", list(range(1, args.N + 1)))
    print("  #  :", [len(data[n]) for n in range(1, args.N + 1)])

    all_paths = [w for n in range(1, args.N + 1) for w in data[n]]
    words = minimal_forbidden_words(all_paths, args.L)
    print(f"\nMinimal forbidden subsequences of length <= {args.L} (guess[T]):")
    print("  ", ", ".join(sorted(words, key=lambda w: (len(w), w))) or "none")

    ok = check_characterization(words, in_class, check_n, data)
    print("\nConclusion:",
          "the words characterize the connected paths for every n <=" if ok
          else "the words do NOT characterize the connected paths (n <=",
          f"{check_n}" + ("." if ok else ")."))

    if args.closure:
        ex = check_closure(data, check_n)
        print("\n--- Closure under connected subsequences ---")
        if ex is None:
            print(f"  closed for n <= {check_n}: a characterization by forbidden "
                  "subsequences may exist.")
        else:
            q, p = ex
            print(f"  NOT closed: {q} (not in the class, involution "
                  f"{fmt_perm(biane_inverse(q))}) is a subsequence of "
                  f"{p} (in the class, involution {fmt_perm(biane_inverse(p))}).")
            print("  Hence NO set of forbidden subsequences can characterize the class.")

    if args.distribution:
        print(f"\n--- Descent distribution over I_n(4321, {label}) ---")
        for n in range(0, args.distribution + 1):
            dist = {0: 1} if n == 0 else descent_distribution(n, in_class)
            print(f"  n={n:>2}: {format_poly(dist)}")

    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
