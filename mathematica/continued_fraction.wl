(* ::Package:: *)

(* continued_fraction.wl
   Companion code for "Involutions avoiding 4321 and another pattern of length four"
   (M. Barnabei, N. Castronuovo, M. Silimbani).

   Taylor expansion of F(x,y), the generating function of the involutions in
   I_n(4321,4231) counted by length (x) and descents (y), computed from the
   continued fraction of the theorem in the section on I_n(4321,4231).

   The coefficient of x^n is the descent polynomial of I_n(4321,4231): the
   coefficient of y^k is the number of such involutions with k descents.  *)

ClearAll[A];

ord = 15;   (* highest power of x to compute *)

(* truncating the continued fraction after N levels gives F exactly up to
   x^(2N^2+8N+5), so a small depth is enough *)
depth = Module[{n = 1}, While[2 n^2 + 8 n + 6 <= ord + 2, n++]; n];

u0  = Series[(1 + (1 - y) x^2 - Sqrt[(1 + (1 - y) x^2)^2 - 4 x^2])/(2 x),
             {x, 0, ord + 2}];
sig = 1 - (1 - y) x u0;
q   = y u0^2/sig^2;
A[n_] := (1 - u0^2)^2 - y u0 (1 + u0^2) q^n;

(* Phi = A_1 - y^2 u0^4 q^3/(A_2 - y^2 u0^4 q^5/(A_3 - ...)), built bottom-up *)
Phi = A[depth];
Do[Phi = A[n] - y^2 u0^4 q^(2 n + 1)/Phi, {n, depth - 1, 1, -1}];

xi = y u0 (1 - u0^2)/(sig (1 - u0^2) - y u0 - y^2 u0^4 q/Phi);
F  = ((1 + xi)/(1 - (1 - y) x (1 + xi)) - 1)/x;

coeffs = Table[Expand[SeriesCoefficient[F, {x, 0, i}]], {i, 0, ord}];

Print[TableForm[Table[{i, coeffs[[i + 1]]}, {i, 0, ord}],
  TableHeadings -> {None, {"n", "descent polynomial"}}]];

(* total numbers of involutions (y = 1):
   expected 1, 1, 2, 4, 8, 18, 38, 87, 192, 440, 998, 2290, ... *)
Print[coeffs /. y -> 1];
