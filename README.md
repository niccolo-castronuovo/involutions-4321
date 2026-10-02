# Involutions avoiding 4321 and another pattern of length four

Companion code for the paper

> M. Barnabei, N. Castronuovo, M. Silimbani,
> *Involutions avoiding 4321 and another pattern of length four*.

The paper enumerates, according to length and number of descents, the
involutions avoiding 4321 and another pattern of length four. Its main tool is
Biane's bijection: involutions avoiding 4321 correspond to (unlabeled) Motzkin
paths, and the avoidance of a second pattern translates into conditions on the
paths. The code in this repository was used to guess these conditions, to test
them, and to check the generating functions obtained in the paper.

## Contents

| File | Language | Purpose |
|---|---|---|
| `involution_paths.py` | Python 3 | Procedure `guess[T]`: forbidden subsequences of the Motzkin paths of a class, and checks of the resulting characterization |
| `mathematica/check_identities.wl` (and `.nb`) | Mathematica | Verification of the identities used in the proof of the continued fraction for I_n(4321,4231) |
| `mathematica/continued_fraction.wl` | Mathematica | Taylor expansion of the generating function F(x,y) of I_n(4321,4231) from the continued fraction |
| `mathematica/recurrence_search.wl` | Mathematica | Search for linear recurrences with polynomial coefficients for F(x,1) (evidence of non-D-finiteness) |

## The Python script

### Requirements

Python 3.8 or later. No external packages are needed.

### What it does

For a set `T` of patterns, avoided together with 4321, and for every length
`n <= N`, the script

1. generates the connected involutions in I_n(4321, T). Instead of generating
   all involutions, it generates Motzkin paths and inverts Biane's bijection,
   since involutions avoiding 4321 correspond exactly to unlabeled Motzkin paths;
2. computes the minimal words over {U, H, D} of length at most `L` that occur as
   a (not necessarily contiguous) subsequence of none of the associated paths:
   this is the procedure `guess[T]` of the paper;
3. checks, for every `n <= N_check`, whether the connected paths avoiding these
   words are **exactly** the paths of the connected involutions of the class.
   The two *sets* of paths are compared, hence also their sizes; in case of a
   mismatch, examples are printed.

Optionally, it also

* tests whether the set of connected paths is closed under taking connected
  subsequences (`--closure`). This condition is necessary and sufficient for the
  existence of *any* characterization by forbidden subsequences;
* prints the joint distribution of length and descents over the whole class
  (`--distribution M`), to be compared with the generating functions of the paper.

### Usage

```
python involution_paths.py PATTERNS [-N N] [-L L] [--check-N N_CHECK] [--closure] [--distribution M]
```

(Depending on your system, the command may be `python3` or `py`.)

| Argument | Meaning | Default |
|---|---|---|
| `PATTERNS` | patterns avoided together with 4321, e.g. `2143` or `2341 4123` | required |
| `-N` | maximal length of the involutions used for the guess | 12 |
| `-L` | maximal length of the forbidden words | 4 |
| `--check-N` | maximal length for the final check | same as `-N` |
| `--closure` | test closure under connected subsequences | off |
| `--distribution M` | print the descent distribution over I_n(4321, T), n <= M | off |

For involutions, avoiding a pattern is equivalent to avoiding its inverse, so
`2341` and `2341 4123` define the same class.

The script returns exit code 0 if the words characterize the connected paths
up to `N_check`, and 1 otherwise.

### Examples

**A class characterized by forbidden subsequences**

```
python involution_paths.py 2143 --check-N 14 --closure
```

The forbidden words are `DU, HUDH`; the check succeeds for every n <= 14 and
the class is closed under connected subsequences. Using `--check-N` larger
than `-N` tests the guessed words on lengths that were not used to find them.

**A class that is not characterized by forbidden subsequences**

```
python involution_paths.py 2341 4123 --closure
```

The check fails at n = 5: the path `UHHHD` avoids all the words, but its
involution 52341 contains 2341. The closure test explains why no set of
forbidden subsequences can work: `UHHHD` is a subsequence of `UHUHDHD`,
whose involution 5274163 belongs to the class.

**The descent distribution**

```
python involution_paths.py 4231 --distribution 7
```

prints, for each n <= 7, the polynomial in y counting the involutions in
I_n(4321,4231) by number of descents, e.g. `n= 6: 1 + 9*y + 19*y^2 + 8*y^3 + y^4`.

### Reproducing the table of the paper

Each command below prints the forbidden words listed in the summary table of
the paper (with `-N 12 --check-N 13`, running time: a few seconds each).

| Class | Command | Forbidden words | Characterization |
|---|---|---|---|
| {1234} | `python involution_paths.py 1234` | DDU, DUU, DDDD, DDDH, DDHH, DHHH, DHHU, HHHH, HHHU, HHUU, HUUU, UUUU | yes |
| {1243} | `python involution_paths.py 1243` | DDH, DDU, DHU, HHU | yes |
| {1342, 1423} | `python involution_paths.py 1342 1423` | DU, DHH, HUU, DHDD, HUHH | yes |
| {1432} | `python involution_paths.py 1432` | DUH, HUH | yes |
| {2143} | `python involution_paths.py 2143` | DU, HUDH | yes |
| {2413, 3142} | `python involution_paths.py 2413 3142` | DH, DU, HU | yes |
| {3241, 4213} | `python involution_paths.py 3241 4213` | HU, HDH | yes |
| {3412} | `python involution_paths.py 3412` | DD, DH, DU, HU, UU | yes |
| {1324} | `python involution_paths.py 1324 --closure` | DHU, DHDH, DUDU, HUHU | no |
| {2341, 4123} | `python involution_paths.py 2341 4123 --closure` | DU, DDD, DDH, DHH, HHU, HUU, UUU, HHHH | no |
| {3421, 4312} | `python involution_paths.py 3421 4312 --closure` | none | no |
| {4231} | `python involution_paths.py 4231 --closure` | none | no |

As always for experimental results, a successful check for n <= N_check is
evidence, not a proof: the characterizations are proved in the paper.

## The Mathematica scripts

Open a `.wl` file in Mathematica and evaluate it (or run
`wolframscript -file mathematica/<name>.wl` from a terminal).

* **`check_identities.wl`** verifies the identities used in the proof of the
  continued fraction: the behavior of the kernel K, of lambda and of omega
  under the change of variable u = (u0+t)/(1+u0 t) (omega becomes the dilation
  t -> q t), the identities giving the explicit coefficients, the identity
  between the two expressions of xi (in the statement and in the proof), and
  the resulting formula for Q(0). Each check returns `0`; the last one is a
  deliberately false identity and returns a nonzero expression. The method
  is exact: the numerator of each expression is reduced modulo K(u0).
  `check_identities.nb` contains the same code as a notebook.
* **`continued_fraction.wl`** computes the coefficients of F(x,y) up to
  x^`ord` (default 15). For larger `ord` with symbolic y the computation
  becomes slow; set `y = 1` at the beginning to get only the totals.
* **`recurrence_search.wl`** computes 421 coefficients of F(x,1), checks them
  against a direct enumeration, and shows that no linear recurrence with
  polynomial coefficients of order r <= 12 and (r+1)(d+1) <= 380 fits them.
  The computation is done modulo a large prime: a nullity 0 modulo p proves
  that no such recurrence exists over the rationals. The same procedure finds
  the recurrence of the Motzkin numbers.

## How to cite

If you use this code, please cite the paper:

```bibtex
@article{BarnabeiCastronuovoSilimbani2026,
  author  = {Barnabei, Marilena and Castronuovo, Niccol{\`o} and Silimbani, Matteo},
  title   = {Involutions avoiding 4321 and another pattern of length four},
  year    = {2026},
  note    = {Submitted}
}
```


## License

MIT License, see `LICENSE`.

## Contact

Niccolò Castronuovo, Università di Bologna — niccolo.castronuovo2@unibo.it
