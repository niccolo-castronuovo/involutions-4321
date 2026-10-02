(* ::Package:: *)

(* recurrence_search.wl
   Companion code for "Involutions avoiding 4321 and another pattern of length four"
   (M. Barnabei, N. Castronuovo, M. Silimbani).

   Experimental evidence that F(x,1), the generating function of
   I_n(4321,4231), is not D-finite:
     1. compute the first 421 coefficients from the continued fraction
        (corollary for y = 1), and check them against a direct enumeration;
     2. search for linear recurrences  sum_{i=0..r} P_i(n) a(n+i) = 0  with
        deg P_i <= d, for every order r <= 12 and every degree d with
        (r+1)(d+1) <= 380, by computing the null space of the associated
        linear system modulo a large prime.
   A nullity 0 modulo p PROVES that no such recurrence fits the computed
   terms (a rational solution would reduce to a nonzero solution mod p).
   As a sanity check, the same procedure detects the Motzkin recurrence.  *)

ClearAll[B, dpCounts, nullity];

(* ---------------- 1. Coefficients of F(x,1) ---------------- *)
L   = 421;                 (* coefficients of x^0 .. x^420 *)
ord = L - 1;

s = Sum[CatalanNumber[k] x^(2 k + 1), {k, 0, Ceiling[ord/2]}] + O[x]^(ord + 2);
q = s^2;
B[n_] := (1 - q)^2 - s q^n (1 + q);

depth = 14;    (* 2N^2+8N+6 > 422 already for N = 13 *)
Psi = B[depth];
Do[Psi = B[n] - s^2 q^(2 n + 2)/Psi, {n, depth - 1, 1, -1}];
F = (1 - q^2)/(1 - q - s - s^2 q^2/Psi);

a = CoefficientList[Normal[F], x][[1 ;; L]];
Print["First terms: ", Take[a, 12]];

(* independent check: direct enumeration of admissible Motzkin paths;
   state (h, r) = (height, down steps still needed before the next H) *)
dpCounts[nmax_] := Module[{cur, nxt, res = {}},
  cur = ConstantArray[0, {nmax + 2, nmax + 2}]; cur[[1, 1]] = 1;
  Do[
   AppendTo[res, Total[cur[[1]]]];
   nxt = ConstantArray[0, {nmax + 2, nmax + 2}];
   Do[
    If[cur[[h + 1, r + 1]] != 0,
     With[{v = cur[[h + 1, r + 1]]},
      If[h + 1 <= nmax, nxt[[h + 2, r + 1]] += v];          (* U *)
      If[h > 0, nxt[[h, Max[r - 1, 0] + 1]] += v];          (* D *)
      If[r == 0, nxt[[h + 1, h + 1]] += v]]],               (* H *)
    {h, 0, nmax}, {r, 0, h}];
   cur = nxt, {n, 0, nmax}];
  res];
Print["Agreement with direct enumeration (n <= 60): ", a[[1 ;; 61]] == dpCounts[60]];

(* ---------------- 2. Recurrence search modulo p ---------------- *)
p = NextPrime[2^31];

nullity[seq_, r_, d_, p_] := Module[{Lq = Length[seq], sp = Mod[seq, p], M},
  M = Table[
    Flatten@Table[
      Mod[sp[[n + i + 1]] If[j == 0, 1, PowerMod[n, j, p]], p],
      {i, 0, r}, {j, 0, d}],
    {n, 0, Lq - r - 1}];
  Length[NullSpace[M, Modulus -> p]]];

results = Table[
   With[{d = Floor[380/(r + 1)] - 1},
    {r, d, (r + 1) (d + 1), L - r, nullity[a, r, d, p]}],
   {r, 1, 12}];
Print[TableForm[results,
  TableHeadings -> {None, {"order r", "max degree d", "unknowns",
     "equations", "nullity mod p"}}]];
(* expected: nullity 0 in every row *)

(* ---------------- 3. Sanity check: Motzkin numbers ---------------- *)
mot = Table[Sum[Binomial[n, 2 k] CatalanNumber[k], {k, 0, Floor[n/2]}], {n, 0, ord}];
Print["Motzkin, order 2, degree 1 (expected 1): ", nullity[mot, 2, 1, p]];
Print["Motzkin, order 1, degree 189 (expected 0): ", nullity[mot, 1, 189, p]];

Mq = Table[Flatten@Table[mot[[n + i + 1]] If[j == 0, 1, n^j], {i, 0, 2}, {j, 0, 1}],
   {n, 0, 60}];
Print["Motzkin recurrence (coefficients of P_0, P_1, P_2): ", NullSpace[Mq]];
(* expected {{-3, -3, -5, -2, 4, 1}}, i.e.
   (n+4) M(n+2) = (2n+5) M(n+1) + 3(n+1) M(n) *)
