(* ::Package:: *)

(* check_identities.wl
   Companion code for "Involutions avoiding 4321 and another pattern of length four"
   (M. Barnabei, N. Castronuovo, M. Silimbani).

   Verifies the identities used in the proof of the continued fraction
   theorem (Section on I_n(4321,4231)). All of them hold because u0 is the
   power-series root of the kernel K.

   Method: every identity is written as "expr = 0"; we put expr over a common
   denominator and take the remainder of its numerator modulo K(u0), viewed
   as a polynomial in u0. The identity holds (when K(u0) = 0) iff the
   remainder is 0. Every check below must return 0, except the last one,
   which is a deliberately false identity used as a sanity check.

   Here K(u) = u - x - x u^2 + (1-y) x^2 u = -x u^2 + (1+(1-y)x^2) u - x,
   exactly as in the paper. *)

ClearAll[K, omega, lambda, uhat, check, q, zeta, sigma, delta, r, gamma,
  xiProof, xiThm, Qu0, Q0, Phi];

(* ---------------- Definitions ---------------- *)
K[u_] := u - x - x u^2 + (1 - y) x^2 u;
omega[u_] := x (1 - (1 - y) x u)/(1 - x u);
lambda[u_] := x (1 - (1 - y) x u)/(1 - x u0);
uhat[t_] := (u0 + t)/(1 + u0 t);
q = x^2 y/(1 - x u0)^2;
zeta = x (u0 - (1 - y) x)/(1 - x u0);
sigma = 1 - (1 - y) x u0;

(* output 0  <-->  expr = 0 when K(u0) = 0 *)
check[expr_] := Simplify[
   PolynomialRemainder[Expand[Numerator[Together[expr]]], K[u0], u0]];

(* ---------------- Identities used in the proof ---------------- *)

(* 1. K(uhat(t)) = x t (1-u0^2)^2 / (u0 (1+u0 t)^2) *)
Print["1. K(uhat(t)):          ",
  check[K[uhat[t]] - x t (1 - u0^2)^2/(u0 (1 + u0 t)^2)]];

(* 2. (1 + u0 t) lambda(uhat(t)) = u0 + zeta t *)
Print["2. lambda(uhat(t)):     ",
  check[(1 + u0 t) lambda[uhat[t]] - (u0 + zeta t)]];

(* 3. omega(uhat(t)) = uhat(q t): after the change of variable,
      omega acts as a dilation *)
Print["3. omega(uhat(t)):      ",
  check[omega[uhat[t]] - uhat[q t]]];

(* ---------------- Identities for the explicit coefficients ---------------- *)
Print["4a. x/(1-x u0) = u0/sigma:          ", check[x/(1 - x u0) - u0/sigma]];
Print["4b. q = y u0^2/sigma^2:             ", check[q - y u0^2/sigma^2]];
Print["4c. 1-zeta = (1-u0^2)/sigma:        ", check[(1 - zeta) - (1 - u0^2)/sigma]];

(* ---------------- The quantity xi ---------------- *)
(* 5. In the proof, substituting d1 = r d0 into the first equation of the
      recurrence gives Q(u0) = gamma E, with E = 1 + x(1-y) Q(0) and
         gamma = (1-u0^2)^2 / ((1-u0^2)^2 - y u0 (1-zeta) - y u0^2 q r),
      and one sets xi = x y gamma/(1 - x u0).
      Here r = y u0^2 (1-zeta)/Phi, where Phi = A_1 - y^2 u0^4 q^3/(A_2 - ...)
      is the continued fraction of the statement, kept as a symbol.
      We check that this xi coincides with the expression of xi given in the
      statement of the theorem. *)
delta = (1 - u0^2)^2;
r = y u0^2 (1 - zeta)/Phi;
gamma = delta/(delta - y u0 (1 - zeta) - y u0^2 q r);
xiProof = x y gamma/(1 - x u0);
xiThm = y u0 (1 - u0^2)/(sigma (1 - u0^2) - y u0 - y^2 u0^4 q/Phi);
Print["5. xi (proof) = xi (statement):     ", check[xiProof - xiThm]];

(* 6. The formula for Q(0). The two relations
         (1 - (1-y) x) Q(0) = 1 + x y Q(u0)/(1 - x u0)     (from u = 0)
         Q(u0) = gamma (1 + x (1-y) Q(0))                  (from the recurrence)
      are satisfied by Q(0) = (1 + xi)/(1 - (1-y) x (1 + xi)). *)
Qu0 = gamma (1 + x (1 - y) Q0);
Print["6. formula for Q(0):                ",
  check[((1 - (1 - y) x) Q0 - 1 - x y Qu0/(1 - x u0)) /.
     Q0 -> (1 + xiProof)/(1 - (1 - y) x (1 + xiProof))]];

(* ---------------- Sanity check ---------------- *)
(* A FALSE identity (the factor (1 + u0 t) is missing):
   the output must be a NONZERO expression. *)
Print["Sanity check (must be nonzero):     ",
  check[lambda[uhat[t]] - (u0 + zeta t)]];
