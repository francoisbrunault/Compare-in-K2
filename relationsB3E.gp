/*******************************************
*                                          *
* FINDING RELATIONS IN THE GONCHAROV-LEVIN *
* GROUP B_3 OF ELLIPTIC CURVES USING       *
* MELLIT'S METHOD OF INCIDENT LINES        *
*                                          *
*******************************************/

/*

Author: François Brunault
Date: November 2022

This file contains PARI/GP code for computing relations in the Goncharov-Levin group B_3(E), where E is an elliptic curve, using Mellit's technique of incident lines in

A. Mellit, Elliptic dilogarithms and parallel lines, J. Number Theory 204 (2019), 1-24.

Such relations have the following applications:

1) Comparing elements in K_2 of elliptic curves. See my appendix for the preprint:

M. Asakura, M. Chida, "A numerical approach toward the p-adic Beilinson conjecture for elliptic curves over Q", arXiv:2003.08888 [math.NT].

2) Proving exotic relations between values of the elliptic dilogarithm.

Note: to use these programs, you should start PARI/GP in the directory containing this file and type the command

\r relationsB3E.gp

Example of use:

E = ellinit("11a3");
relationsB3E(E)

returns

[[[-3, [1]], [2, [2]]]]

meaning that -3[p] + 2[2p] = 0 in B_3(E) \otimes Q, where p = (0,0) is the generator of E(Q)_tors given by elltors(E)[3].

This implies the exotic relation -3 D_E(p) + 2 D_E(2p) = 0, where D_E is the elliptic dilogarithm on E.

*/


\\ The function beta below computes the
\\ convolution \beta(D1, D2) of two
\\ divisors D1 and D2 on the group
\\ G = Z/N or Z/N1 * Z/N2, encoded as
\\ [N] or [N1, N2] respectively.
\\ The argument G_list is the vector,
\\ sorted lexicographically, of the
\\ elements of G, represented as [a] or
\\ [a, b] respectively (if G_list is not
\\ given, the function creates it first).
\\ A divisor D = \sum_{p \in G} n_p (p)
\\ is encoded as the vector of 2-component
\\ vectors [n_p, p] with p in G_list.
\\ The output is a divisor on G represented
\\ in the same way but sorted
\\ lexicographically with respect to p in
\\ G_list.
\\ If the optional argument simple is
\\ nonzero, the function simplifies the
\\ image of the divisor in the group B_3(E).

\\ The function G_to_list(G) creates the list of
\\ the elements of G, sorted lexicographically.

{
G_to_list(G) =
  my(G_list = []);
  forvec(X = vector(#G, i, [0, G[i]-1]),
    G_list = concat(G_list, [X]));
  G_list;
}

\\ Reduce a point p mod G

{
redG(p, G) =
  vector(#G, i, Mod(p[i], G[i]));
}

\\ Simplify in B_3(E) a divisor on the
\\ group G = Z/N or Z/N1 * Z/N2. The divisor
\\ is given as a vector L of integers
\\ of length #G, representing the divisor
\\ \sum_{i=1}^{#G} L[i] (G_list[i]).
\\ The arguments G and G_list are as
\\ explained above for the function beta.
\\ The function returns a divisor encoded
\\ as a vector of 2-component vectors
\\ [n_p, p] with p in G_list.

{
simplifyB3E(L, G, G_list) =
  my(p, k, pts = [], Lsimple = [], Lfinal = []);
  \\ Create G_list if not given
  if(G_list === 0,
    G_list = G_to_list(G));
  for(i = 1, #G_list,
    p = redG(G_list[i], G);
    \\ Skip the 2-torsion points
    if(2*p == 0,
      next());
    \\ Check if -p already appears in Lsimple
    k = vecsearch(pts, lift(-p));
    if(k == 0,
      \\ Incorporate p into Lsimple
      pts = concat(pts, [lift(p)]);
      Lsimple = concat(Lsimple, L[i]),
      \\ Otherwise update Lsimple
      Lsimple[k] -= L[i]));
  \\ Remove the zero terms in the divisor
  for(i = 1, #pts,
    if(Lsimple[i] != 0,
      Lfinal = concat(Lfinal, [[Lsimple[i], pts[i]]])));
  Lfinal;
}

{
beta(D1, D2, G, G_list, simple = 0) =
  my(L, n1, p1, n2, p2, p, k);
  \\ Create G_list if not given
  if(G_list === 0,
    G_list = G_to_list(G));
  \\ Reduce D1 and D2 mod G
  for(i = 1, #D1,
      D1[i][2] = redG(D1[i][2], G));
  for(i = 1, #D2,
      D2[i][2] = redG(D2[i][2], G));
  L = vector(#G_list);
  for(i = 1, #D1,
    n1 = D1[i][1];
    p1 = D1[i][2];
    for(j = 1, #D2,
      n2 = D2[j][1];
      p2 = D2[j][2];
      p = lift(p1-p2);
      k = vecsearch(G_list, p);
      L[k] += n1*n2));
  \\ If simple = 0, return the divisor L
  if(simple == 0,
    return(vector(#G_list, i, [L[i], G_list[i]])));
  \\ Otherwise, simplify the divisor L
  simplifyB3E(L, G, G_list);
}

\\ Equation of the line passing through two
\\ given points p and q.
\\ If p = q, returns the tangent of E at p.
\\ The function returns [a, b, c], where
\\ aX + bY + cZ = 0 is a homogeneous equation
\\ of the line.

{
line(e, p, q) =
  my(xp, yp, xq, yq, Fx, Fy);
  if((p == [0]) && (q == [0]),
    return([0, 0, 1]));
  if(p == [0],
    return([1, 0, -q[1]]));
  if(q == [0],
    return([1, 0, -p[1]]));
  xp = p[1];
  yp = p[2];
  if(p == q,
    Fx = e.a1*yp - 3*xp^2 - 2*e.a2*xp - e.a4;
    Fy = 2*yp + e.a1*xp + e.a3;
    return([Fx, Fy, - xp*Fx - yp*Fy]));
  xq = q[1];
  yq = q[2];
  [yp - yq, xq - xp, xp*yq - yp*xq];
}

\\ Given an elliptic curve E and a set
\\ of points p_list, the function lines
\\ returns the list of lines l such
\\ that l \cap E is contained in p_list.
\\ The set p_list should form a subgroup
\\ G of E. The arguments G and G_list
\\ are as in the function beta above,
\\ and the elements of p_list are the
\\ actual points of E corresponding to
\\ the elements of G_list.
\\ The function returns a vector
\\ of vectors [eq, [p, q, r]], where
\\ eq is an equation of the line, and
\\ p, q, r are the points (encoded as
\\ elements of G_list) through which
\\ the line passes.

{
lines(e, p_list, G, G_list) =
  my(l_list = [], n, p, q, r, k, eq);
  \\ Create G_list if not given
  if(G_list === 0,
    G_list = G_to_list(G));
  \\ Loop over triples of points in G_list
  for(i = 1, #G_list,
    p = G_list[i];
    for(j = i, #G_list,
      q = G_list[j];
      r = lift(redG(-(p+q), G));
      k = vecsearch(G_list, r);
      if(k >= j,
        eq = line(e, p_list[i], p_list[j]);
        l_list = concat(l_list, [[eq, [p, q, r]]]))));
  l_list;
}

\\ Finding the relations in B_3(E)

\\ Given an elliptic curve E over Q,
\\ the command relationsB3E(E) finds the
\\ incident lines supported in E(Q)_tors,
\\ and compute the associated relations
\\ in the group B_3(E) \otimes Q.
\\
\\ The function returns a basis of these
\\ relations, in the form of a vector of
\\ divisors D. Each divisor D is encoded
\\ as a vector of 2-component vectors
\\ [n_p, p] with n_p in Z and p in G_list,
\\ where G_list represents the group
\\ E(Q)_tors (see the description of the
\\ function 'beta' above).
\\
\\ If the optional argument verb is
\\ nonzero, the function prints each
\\ triple of incident lines whose
\\ intersection point does not lie on E
\\ (with the same format as in the
\\ function 'lines' above), as well as
\\ the associated relation in B_3(E),
\\ given in non-simplified form as a
\\ vector of integers of size #G, and in
\\ simplified form as a divisor D as
\\ above.

{
relationsB3E(e, verb = 0) =
  my(e_tors, G, G_list, p_list = [], p, l_list, r_list = [], l1, l2, l3, S1, S2, S3, M, p0, D1, D2, D3, b12, b23, b31, rel, rel_simple, v, pos, rels);
  \\ Compute E(Q)_tors
  e_tors = elltors(e);
  G = e_tors[2];
  G_list = G_to_list(G);
  \\ Create p_list = E(Q)_tors
  for(i = 1, #G_list,
    p = [0];
    for(j = 1, #G,
      p = elladd(e, p, ellmul(e, e_tors[3][j], G_list[i][j])));
    p_list = concat(p_list, [p]));
  l_list = lines(e, p_list, G, G_list);
  \\ Loop over triples of lines in l_list
  for(i = 1, #l_list,
    l1 = l_list[i];
    for(j = i+1, #l_list,
      l2 = l_list[j];
      \\ Test if l1 and l2 intersect on E
      S1 = Set(l1[2]);
      S2 = Set(l2[2]);
      if(#setintersect(S1, S2),
        next());
      \\ Compute the intersection point
      \\ p0 = l1 \cap l2
      M = Mat([l1[1], l2[1]]~);
      p0 = matker(M)[, 1];
      for(k = j+1, #l_list,
        l3 = l_list[k];
        \\ Test if l3 intersects l1 or l2 on E
        S3 = Set(l3[2]);
        if(#setintersect(S3, setunion(S1, S2)),
          next());
        \\ Test if l1, l2, l3 are incident
        if(l3[1] * p0 == 0,
          \\ Print the incident lines
          if(verb,
            print("Incident lines: l1 = ", l1, ", l2 = ", l2, ", l3 = ", l3));
          \\ Compute the associated relation
          D1 = vector(3, m, [1, l1[2][m]]);
          D2 = vector(3, m, [1, l2[2][m]]);
          D3 = vector(3, m, [1, l3[2][m]]);
          rel = vector(#G_list);
          b12 = beta(D1, D2, G, G_list);
          b23 = beta(D2, D3, G, G_list);
          b31 = beta(D3, D1, G, G_list);
          rel = vector(#G_list, j, b12[j][1] + b23[j][1] + b31[j][1]);
          rel_simple = simplifyB3E(rel, G, G_list);
          \\ Print the relation
          if(verb,
            print("Relation: ", rel);
            print("In B_3(E), simplifies to: ", rel_simple);
            print());
          \\ Write rel_simple in vector form
          v = vectorv(#G_list);
          for(m = 1, #rel_simple,
            pos = vecsearch(G_list, rel_simple[m][2]);
            v[pos] = rel_simple[m][1]);
          r_list = concat(r_list, [v])))));
  \\ Span of the simplified relations
  rels = Vec(matimage(Mat(r_list)));
  rels = apply(v -> simplifyB3E(v, G, G_list), rels);
  if(verb,
    print("Basis of the relations in B_3(E): ");
    for(i = 1, #rels,
      print("  ", rels[i])));
  rels;
}
