# Compare-in-K2

*Author:* François Brunault

*Date:* November 2022

This repository contains PARI/GP scripts to prove relations in $K_2$ of elliptic curves. They are used in the appendix (by F. Brunault) of the preprint

M. Asakura, M. Chida, [A numerical approach toward the p-adic Beilinson conjecture for elliptic curves over Q](https://arxiv.org/abs/2003.08888).

In this appendix, we compare two elements in $K_2$ of the elliptic curve $X_0(24)$: the element $\xi$ in (6.10) and the Beilinson–Kato element. This comparison is used in the proof of Theorem 6.2.

The file `compareK2E.gp` contains all the computations needed to establish this relation. To use it, you should start PARI/GP in the directory containing this file and the file `relationsB3E.gp`. Then type the command

```
\r compareK2E.gp
```

These computations rely on the PARI/GP scripts in the file `relationsB3E.gp`. These scripts allow to find relations in the Goncharov-Levin group $B_3(E)$, where $E$ is an elliptic curve, using Mellit's technique of incident lines in

A. Mellit, Elliptic dilogarithms and parallel lines, J. Number Theory 204 (2019), 1-24.

Relations in $B_3(E)$ imply exotic relations between values of the elliptic dilogarithm on $E$.

Note: to use these programs, you should start PARI/GP in the directory containing this file, and then type the command

```
\r relationsB3E.gp
```

**Example of use:**

```
E = ellinit("11a3");
relationsB3E(E)
```

returns

```
[[[-3, [1]], [2, [2]]]]
```

meaning that $-3[p] + 2[2p] = 0$ in $B_3(E) \otimes \mathbb{Q}$, where $p = (0,0)$ is the generator of $E(\mathbb{Q})_{\textrm{tors}}$ given by `elltors(E)[3]`. This implies the exotic relation $-3 D_E(p) + 2 D_E(2p) = 0$, where $D_E$ is the elliptic dilogarithm on $E$. More information can be obtained by typing instead

```
relationsB3E(E, 1)
```



Shield: [![CC BY 4.0][cc-by-shield]][cc-by]

This work is licensed under a
[Creative Commons Attribution 4.0 International License][cc-by].

[![CC BY 4.0][cc-by-image]][cc-by]

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-image]: https://i.creativecommons.org/l/by/4.0/88x31.png
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg
