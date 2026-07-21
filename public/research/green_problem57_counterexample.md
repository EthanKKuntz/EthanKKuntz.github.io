# A counterexample to Green's Problem 57

**Working note, 20 July 2026.** This note gives an explicit negative answer to
Problem 57 in Ben Green's *100 Open Problems*. The example is over
\(G=\mathbb Z/3\mathbb Z\). The numerical search that found the example is not
part of the proof: the witness and the separating upper bound below have exact,
finite certificates.

The result has not yet been independently checked. A targeted web search on 20
July 2026 found no prior solution beyond Green's statement of the problem.

## 1. The two spaces

For functions bounded pointwise in modulus by one, put

\[
 \phi(g)=\mathbb E_{x_1+x_2+x_3=g}
 f_1(x_2,x_3)f_2(x_1,x_3)f_3(x_1,x_2).
\]

Let \(\Phi(G)\) be the convex hull of these functions. Let \(\Phi'(G)\) be
defined in the same way, except that \(f_3(x_1,x_2)\) must be a function of
\(x_1+x_2\). Clearly \(\Phi'(G)\subseteq\Phi(G)\). We prove that the inclusion
is strict when \(G=\mathbb Z/3\mathbb Z\).

All indices below are taken modulo 3. Let

\[
 h(0)=20,\qquad h(1)=-14,\qquad h(2)=5.
\]

We use the real linear functional

\[
 L(\phi)=\operatorname{Re}\sum_{g\in G}h(g)\phi(g).
\]

The normalization in the conditional expectation contributes the same factor
\(1/9\) on both sides. It is therefore enough to compare the corresponding
unnormalized trilinear score

\[
 T=\operatorname{Re}\sum_{x_1,x_2,x_3\in G}
 h(x_1+x_2+x_3)f_1(x_2,x_3)f_2(x_1,x_3)f_3(x_1,x_2).
\]

## 2. An unrestricted witness with score greater than 309

Let \(\zeta=e^{2\pi i/48}\). Define \(f_r=\zeta^{E_r}\) entrywise, with
\(f_1\) indexed by \((x_2,x_3)\), \(f_2\) by \((x_1,x_3)\), and \(f_3\) by
\((x_1,x_2)\), where

\[
E_1=\begin{pmatrix}
43&35&38\\
1&12&9\\
47&20&4
\end{pmatrix},\quad
E_2=\begin{pmatrix}
13&40&30\\
6&5&36\\
7&16&35
\end{pmatrix},\quad
E_3=\begin{pmatrix}
41&8&37\\
21&31&46\\
46&42&12
\end{pmatrix}.
\]

Every entry has modulus one. Direct collection of the 27 terms gives

\[
\begin{aligned}
T={}&30+147\cos\frac{\pi}{24}+62\cos\frac{\pi}{12}
 +40\cos\frac{\pi}{8}+15\cos\frac{5\pi}{12}\\
&-10\cos\frac{11\pi}{24}+14\cos\frac{\pi}{4}
 +28\cos\frac{\pi}{6}.
\end{aligned}
\]

The elementary rational bounds

\[
\begin{gathered}
\cos(\pi/24)>.991,\quad \cos(\pi/12)>.965,\quad
\cos(\pi/8)>.923,\quad \cos(5\pi/12)>.258,\\
\cos(11\pi/24)<.131,\quad \cos(\pi/4)>.707,\quad
\cos(\pi/6)>.866
\end{gathered}
\]

give

\[
 T>\frac{309133}{1000}>309.
\]

For example, the first and fifth inequalities follow from
\(\cos x>1-x^2/2\), \(\sin x<x\), and \(\pi<22/7\); the remaining bounds
follow immediately from the standard radical values at multiples of
\(\pi/24\). Direct numerical evaluation is \(T\approx309.3102071492\).

## 3. Every restricted witness has score less than 309

Suppose now that \(f_3(x_1,x_2)=c(x_1+x_2)\), and write \(c_s=c(s)\). For
fixed \(x_3=z\), the relevant \(3\times3\) matrix is

\[
 H_z(x,y)=h(x+y+z)c_{x+y}.
\]

The two vectors supplied by \(f_1(\cdot,z)\) and \(f_2(\cdot,z)\) have
Euclidean norm at most \(\sqrt3\). Hence their contribution has modulus at
most \(3\lVert H_z\rVert_{\mathrm{op}}\). Since \(H_z\) is a reverse
circulant matrix, if \(\omega=e^{2\pi i/3}\), then

\[
 \lVert H_z\rVert_{\mathrm{op}}
 =\max_{j\in\mathbb Z/3\mathbb Z}
 \left|\sum_{s\in G}h(s+z)c_s\omega^{js}\right|.
\]

Choose a maximizing \(j_z\) for each \(z\), and define

\[
 A_{\mathbf j}(z,s)=h(s+z)\omega^{j_zs}.
\]

It remains to bound \(\sum_z|(A_{\mathbf j}c)_z|\). Introduce phases
\(d_z\) so that this sum is \(\operatorname{Re}(d^*A_{\mathbf j}c)\), and
put

\[
 C_{\mathbf j}=\frac12
 \begin{pmatrix}0&A_{\mathbf j}^*\\A_{\mathbf j}&0\end{pmatrix}.
\]

If \(y_1,\dots,y_6>0\) satisfy

\[
 \operatorname{diag}(y_1,\dots,y_6)-C_{\mathbf j}\succeq0,
\]

then, because \(|c_s|,|d_z|\leq1\),

\[
 \sum_z|(A_{\mathbf j}c)_z|\leq\sum_{i=1}^6y_i.
\]

Up to adding a constant to every \(j_z\), negating every \(j_z\), and
cyclically shifting the \(z\)-coordinates, there are only three patterns:
all equal, exactly two equal, and all distinct. Representatives and exact dual
certificates are

| pattern \(\mathbf j\) | \(100(y_1,\ldots,y_6)\) | \(\sum_i y_i\) |
|---|---:|---:|
| \((0,0,0)\) | \((1478,1478,1478,1478,1478,1478)\) | \(88.68\) |
| \((0,0,1)\) | \((1772,1926,1448,1542,1786,1817)\) | \(102.91\) |
| \((0,1,2)\) | \((1713,1713,1713,1713,1713,1713)\) | \(102.78\) |

Each claimed matrix is positive definite. This can be checked exactly over
\(\mathbb Q(\omega)\) by Sylvester's criterion. The six leading principal
minors in the only nonuniform case \((0,0,1)\) are

\[
\frac{443}{25},\quad
\frac{426609}{1250},\quad
\frac{77216229}{15625},\quad
\frac{105027914611}{3125000},\quad
\frac{8430332160187}{39062500},\quad
\frac{241353436796433}{7812500000},
\]

all positive. The accompanying exact checker verifies this and the two other
patterns using rational arithmetic in \(\mathbb Q(\omega)\).

It follows uniformly in \(\mathbf j\) that

\[
 \sum_z|(A_{\mathbf j}c)_z|<103,
\]

and therefore every restricted elementary function has \(T<3\cdot103=309\).
The same is true of every convex combination of restricted elementary
functions.

## 4. Conclusion

The unrestricted witness has \(T>309\), whereas every element of
\(\Phi'(\mathbb Z/3\mathbb Z)\) has \(T<309\). Consequently

\[
 \boxed{\Phi'(\mathbb Z/3\mathbb Z)\subsetneq
 \Phi(\mathbb Z/3\mathbb Z).}
\]

Thus Green's Problem 57 has a negative answer.

## 5. Reproducibility and sources

The file `problem57_counterexample_check.py` checks the collection of the 27
witness terms and all exact positive-definiteness certificates. It uses only the
Python standard library.

- Ben Green, [*100 Open Problems*](https://people.maths.ox.ac.uk/greenbj/papers/open-problems.pdf), Problem 57.
- Kevin Buzzard, [*Human mathematicians are being outcounterexampled*](https://xenaproject.wordpress.com/2026/07/20/human-mathematicians-are-being-outcounterexampled/), for the methodological context motivating the search.
