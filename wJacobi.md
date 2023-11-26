---
title: "Partial and Weighted Jacobi"
author: 
- Jan de Leeuw - University of California Los Angeles
date: 'Started November 07 2023, Version of November 16, 2023'
output:
  bookdown::html_document2:
    keep_md: yes
    css: preamble.css
    toc: true
    toc_depth: 3
    number_sections: yes
  bookdown::pdf_document2:
    latex_engine: lualatex 
    includes:
      in_header: preamble.tex
    keep_tex: yes
    toc: true
    toc_depth: 3
    number_sections: yes
graphics: yes
mainfont: Times New Roman
fontsize: 12pt
bibliography: ["mypubs.bib","total.bib"]
abstract: TBD
---





**Note:** This is a working paper which will be expanded/updated frequently. All suggestions for improvement are welcome. 

# Introduction

@deleeuw_ferrari_R_08a
@deleeuw_E_17o
@deleeuw_E_18g

# Just the formulas

Suppose $A$ is a symmetric matrix of order $n$, not necessarily positive semi-definite. We want to find a square orthonormal $K$ (a *rotation matrix*) such that
\begin{equation}
\sigma(K)=\sum_{i=1}^n\sum_{j=1}^n w_{ij}\{K'AK\}_{ij}^2
\end{equation}
is minimized. Here the $W=\{w_{ij}\}$ are symmetric non-negative weights. In most cases that interest us the matrix $W$ is *hollow* (i.e. has zero diagonal) and *binary* (i.e. all elements are either zero or one), but we shall first look at the general case. The minimum of $\sigma$ exists, because the set of rotation matrices is compact and $\sigma$ is continuous and bounded below by zero.

Following @jacobi_46 we build up $K$ using elementary rotations, constructed by using the unit vectors
$e_i$, which have element $i$ equal to one and all other elements equal to zero.
Suppose $X$ is a matrix
with column $x_i$ equal to $\alpha e_i+\beta e_j$ and column $x_j$ equal to $\alpha e_j-\beta e_i$,
where $e_i$ and $e_j$ are units vectors, where $\alpha^2+\beta^2=1$, and where $i<j$. Column $k$ for $k\not\in\{i, j\}$ is equal to $e_k$. More explicitly we could write $X$ as $X_{ij}(\alpha,\beta)$. Clearly $X$ is square orthonormal. 

The general idea of the Jacobi method is that we have an infinite sequence
$(i(\nu),j(\nu))$ of *pivots*, leading to the infinite sequence
$X^{\nu}_{i(\nu),j(\nu)}(\alpha^{\nu},\beta^{\nu})$ where &\alpha^{\nu}$ and $\beta^{\nu}$
are chosen to minimize $\sigma$. We then replace $A^{(\nu)}$ by 
${A^{(\nu+1)}=$ and $X^{(\nu)}$ by $X^{(\nu+1)}=$
$\overline{X}$ and go to the next pivot in the sequence.

Let's look at the problem of optimizing , i.e. at a singler pivot. Then $\sigma$ is a function of $(\alpha,\beta)$ on the unit cicle. Define the symmetric matrix $\overline{A}=X'AX$. Then $\overline{a}_{kl}=x_k'Ax_l$,
which means that for $k\not= i$ and $k\not= j$ as well as $l\not= i$ and $l\not= j$ we have
$\overline{a}_{kl}=a_{kl}$. For $k\not\in\{i,j\}$ we have
\begin{align}
\overline{a}_{ik}&=x_i'Ae_k=x_i'a_k=\alpha a_{ik}+\beta a_{jk},\\
\overline{a}_{jk}&=x_j'Ae_k=x_j'a_k=-\beta a_{ik}+\alpha a_{jk}.
\end{align}
Moreover
\begin{align}
\overline{a}_{ij}&=x_i'Ax_j=(\alpha e_i+\beta e_j)'A(\alpha e_j-\beta e_i)=
(\alpha^2-\beta^2)a_{ij}+\alpha\beta(a_{jj}-a_{ii}),\\
\overline{a}_{ii}&=x_i'Ax_i=(\alpha e_i+\beta e_j)'A(\alpha e_i+\beta e_j)=
\alpha^2a_{ii}+\beta^2a_{jj}+2\alpha\beta a_{ij},\\
\overline{a}_{jj}&=x_j'Ax_j=(-\beta e_i+\alpha e_j)'A(-\beta e_i+\alpha e_j)=
\beta^2a_{ii}+\alpha^2a_{jj}-2\alpha\beta a_{ij}.
\end{align}
In summary
\begin{equation}
\overline{a}_{kl}=\begin{cases}
a_{kl}&\text{ if }k\not=\{i,j\}\text{ and }l\not=\{i,j\},\\
&\text{ if }k=i\text{ and }l\not=\{i,j\}\text{ or }k\not=\{i,j\}\text{ and }l=i,\\
&\text{ if }k=j\text{ and }l\not=\{i,j\},\\
&\text{ if }k=i\text{ and }l=j,\\
&\text{ if }k=i\text{ and }l=i,\\
&\text{ if }k=j\text{ and }l=j.
\end{cases}
\end{equation}

Thus
\begin{align}
\sigma(X)&=\sum_{_{k\not\in\{i,j\}}}^n\sum_{_{l\not\in\{i,j\}}}^nw_{kl}^{\ }a_{kl}^2+\\
&+2\sum_{k\not\in\{i,j\}}^nw_{ik}(\alpha a_{ik}+\beta a_{jk})^2+\\
&+2\sum_{k\not\in\{i,j\}}^nw_{jk}(-\beta a_{ik}+\alpha a_{jk})^2+\\
&+2w_{ij}\{(\alpha^2-\beta^2)a_{ij}+\alpha\beta(a_{jj}-a_{ii})\}^2+\\
&+w_{ii}(\alpha^2a_{ii}+\beta^2a_{jj}+2\alpha\beta a_{ij})^2+\\
&+w_{jj}(\beta^2a_{ii}+\alpha^2a_{jj}-2\alpha\beta a_{ij})^2
\end{align}

Trigonometry

$\alpha=\sin(\theta)$ and $\beta=\cos(\theta)$

# The Sequence

# Special Cases

Hollow
Binary

# Applications

## Symmetric Matrices

All eigenvalues
Some eigenvalues

## Pairs of Matrices

## Rectangular Matrices

## Simultaneous Diagonalization

## DMCA

## GCCA

# Appendix: Code

## template.R


```r
mPrint <- function(x,
                   digits = 6,
                   width = 8,
                   format = "f",
                   flag = "+") {
  print(noquote(
    formatC(
      x,
      digits = digits,
      width = width,
      format = format,
      flag = flag
    )
  ))
}

butLast <- function(x, m = 1) {
  return(rev(rev(x)[-(1:m)]))
}

butFirst <- function(x, m = 1) {
  return(x[-(1:m)])
}
```

# References
