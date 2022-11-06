/************************************/
/* COMPARING IN K_2(E) FOR E = 24a1 */
/************************************/

/*

Author: François Brunault
Date: November 2022

This file contains the PARI/GP scripts used in the appendix (by F. Brunault) of the preprint

M. Asakura, M. Chida, "A numerical approach toward the p-adic Beilinson conjecture for elliptic curves over Q", arXiv:2003.08888 [math.NT].

In this appendix, we compare two elements in K_2 of the elliptic curve X_0(24): the element \xi in (6.10) and the Beilinson–Kato element. This comparison is used in the proof of Theorem 6.2.

To use these programs, you should start PARI/GP in the directory containing this file and the file "relationsB3E.gp". Then type the command

\r compareK2E.gp

This will run all the computations needed to prove the relation mentioned above.

*/


\\ Load the scripts related to B_3(E)

read("relationsB3E.gp");


/********************************/
/* MODULAR PARAMETRISATION OF E */
/********************************/


\\ Torsion subgroup of E

E = ellinit("24a1");

EQ_tors = elltors(E);

print("Torsion subgroup of E = 24a1: ", EQ_tors);
\\ [8, [4, 2], [[0, 2], [1, 0]]]

[P1, P2] = EQ_tors[3];
print("p1 = ", P1);
print("p2 = ", P2);
print();

\\ The modular parametrisation X_0(24) -> E

mf = mfinit([24, 2]);
f = mfeigenbasis(mf)[1];
s = mfsymbol(mf, f);

\\ The command
\\ mfsymboleval(s, [a, b])*2*Pi*I
\\ returns the integral of
\\ omega_f = 2*Pi*I*f(z)*dz
\\ over the modular symbol {a,b}.
\\ Note that mfsymboleval returns
\\ a polynomial, so we use polcoef.

{
phiE(c) =
  ellztopoint(E, polcoef(mfsymboleval(s, [oo, c])*2*Pi*I, 0));
}

cusps = [oo, 0, 1/2, 1/3, 1/4, 1/6, 1/8, 1/12];

print("Cusps of X_0(24): ", cusps);
print();

phiE_cusps = apply(phiE, cusps);

print("Images of the cusps:");
print();

{
for(i = 1, #cusps,
  print(cusps[i], ": ", phiE_cusps[i]));
print();
}

\\ List of torsion points of E

{
Etors = [];
for(a = 0, 3,
  my(P = ellmul(E, P1, a));
  Etors = concat(Etors, [[P, [a, 0]], [elladd(E, P, P2), [a, 1]]]));
Etors = vecsort(Etors, 2);
}

\\ Recognise the images of the cusps
\\ as rational points on E

print("In what follows [a,b] means the point a * p1 + b * p2");
print();

{
for(i = 1, #cusps,
  my(p = phiE_cusps[i]);
  p = apply(bestappr, p);
  for(j = 1, #Etors,
    my(P = Etors[j]);
    if(p == P[1],
      print(cusps[i], ": ", P[2]))));
print();
}


/************************************/
/* THE MODULAR UNIT u_24 ON X_0(24) */
/************************************/


\\ Order of vanishing of u_24 at the cusps of X_0(24)

B2(x) = x^2-x+1/6; \\ Bernoulli polynomial B_2

\\ Order of vanishing of
\\ u_N = \prod_{(a,N)=1} \prod_b g_{a,b}
\\ at the cusp 1/d with d dividing N

{
ord(N, d) =
  d * eulerphi(N) / (2 * gcd(d, N/d) * eulerphi(d)) * sum(a = 1, d, (gcd(a, d) == 1) * B2(a/d));
}

print("Order of vanishing of u_24 at the cusps:");
print();

{
my(c, d);
for(i = 1, #cusps,
  c = cusps[i];
  if(c == oo,
    d = 24,
    d = denominator(c));
  print(c, ": ", ord(24, d)));
print();
}


/************************************/
/* COMPARING THE ELEMENTS IN K_2(E) */
/************************************/


print("Comparing the elements in K_2(E):");
print();

\\ Computation of \beta(z_E) where z_E
\\ is the Beilinson-Kato element from
\\ Section 8.2.3.

div_u24 = [[1, [0,0]], [4, [3,1]], [-2, [1,1]], [-4, [3,0]], [2, [1,0]], [-1, [2,1]], [-1, [0,1]], [1, [2,0]]];

div_W24_u24 = [[1, [3,1]], [4, [0,0]], [-2, [2,0]], [-4, [0,1]], [2, [2,1]], [-1, [1,0]], [-1, [3,0]], [1, [1,1]]];

b24 = beta(div_u24, div_W24_u24, [4, 2]);
b24s = beta(div_u24, div_W24_u24, [4, 2], , 1);

print("36 * beta(u_24, W_24 u_24) = ", b24);
print("In B_3(E), simplifies to: ", b24s);
print();

\\ Computation \beta(\xi) where \xi is
\\ the symbol {f, g} from Section 8.2.4.

div_f = [[-1, [3,0]], [1, [1,0]], [-1, [3,1]], [1, [1,1]]];

div_g = [[4, [2,1]], [-4, [0,1]]];

bE = beta(div_f, div_g, [4, 2]);
bEs = beta(div_f, div_g, [4, 2], , 1);

print("beta(f, g) = ", bE);
print("In B_3(E), simplifies to: ", bEs);
print();

\\ Bloch-Wigner dilogarithm D

D(z) = polylog(2, z, 2); 

\\ Elliptic dilogarithm D_E(p)

{
DE(e, p) =
  my(q, x);
  \\ Write E(C) = C^*/q^Z
  q = exp(2*Pi*I * -e.omega[2]/e.omega[1]);
  \\ Write p = [x] with x in C^*
  x = exp(2*Pi*I * ellpointtoz(e, p)/e.omega[1]);
  suminf(n = 0, D(x*q^n)) - suminf(n = 1, D(q^n/x));
}

DE_P1 = DE(E, P1);
DE_P1P2 = DE(E, elladd(E, P1, P2));

\\ Find numerical relation between
\\ D_E(p1) and D_E(p1+p2)

print("D_E(p1) = ", DE_P1);
print("D_E(p1+p2) = ", DE_P1P2);
print("lindep: ", lindep([DE_P1, DE_P1P2]));
print();

print("Finding relations in B_3(E)...");
print();

\\ E = ellinit("24a1");

relationsB3E(E, 1);
