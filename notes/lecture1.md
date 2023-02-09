# 18.656 Lecture 1 Reading Notes: Introduction
We note the first lecture covers 1.2.2 and 2.1.1. 
We will briefly summarize the core ideas of other sections for completion.

## 1.1 Classical vs high-dimensional theory
* ***Classical theory*** provides statements that apply to a fixed class of models, parameterized by an index $n$ that is allowed to increase.
    * For example, the ***law of large numbers*** guarantees the sample mean $\hat{\mu}_n := \frac{1}{n} \sum_{i=1}^{n} X_i$ converges in probability to $\mu$.
    * Consequently, $\hat{\mu}_n$ is a consistent estimator of the unknown population mean.
    * The more refined statement is provided by the ***central limit theorem***, guaranteeing the rescaled devitation $\sqrt{n}(\hat{\mu}_n-\mu)$ converges in distribution to a centered Gaussian with covariance matrix $\Sigma = cov(X_1)$. 
> The two theoretical statements underlie the analysis of various statistical estimators, in particular to ensure their consistency and aymptotic normality, respectively.

* Some essential facts that motivate the discussion are
    1. Datasets nowadays frequently have a "high-dimensional flavor", with $d$ on the same order or larger than the sample size $n$.
    2. Classical "large $n$, fixed $d$" theory fails to provide useful prediction frequently.
    3. Classical methods can break down dramatically in high-dimensional regimes.

## 1.2 What can go wrong in high dimensions?

### 1.2.2 Covariance estimation
* **Given:** a collection of random vectors $\{x_1, \cdots x_n\}$, where each $x_i$ is drawn in an i.i.d. manner from some zero-mean distribution in $R^d$.
* **Goal:** estimate the unknown covariance matrix $\Sigma = cov(X)$.
* A natrual (unbiased, $E[\hat{\Sigma}] = \Sigma$) estimator is the *sample covariance matrix* $$\hat{\Sigma}:=\frac{1}{n} \sum_{i=1}^{n} x_ix_i^T $$ a $d \times d$ random matrix corresponding to the sample average of the outer products $x_ix_i^T \in R^{d \times d}$.
* **[Classical analysis]** considers the behavior of $\hat{\Sigma}$ as $n \uparrow$ while $d$ stays fixed.
    * ie. $l_2$-operator norm is a matrix norm used to measure the distance between $\hat{\Sigma}$ and $\Sigma$ $$\||\hat{\Sigma}-\Sigma\||_2 := \sup_{u \neq 0} \frac{\|(\hat{\Sigma}-\Sigma)u\|_2}{\|u\|_2}$$
    * Under the classical law of large numbers, $\||\hat{\Sigma}-\Sigma\||_2$ converges to zero almost surely as $n \rightarrow \infty$.
* **[High-dimensional asymtotics] We ask:** is this consistency preserved when $d \rightarrow \infty$?
    * Suppose both $n$ and $d$ increase with their ratio $d/n=\alpha \in (0,1)$ remaining fixed. We run simulations for a random ensemble $\Sigma = I_d$ with each $X_i \sim N(0,I_d)$ for $i = 1, \cdots, n$.
    * If the sample covariance converges to the identity matrix, then the vector of eigenvalues $\gamma (\hat{\Sigma})$ should converge to the all-ones vector. 
    * However, as shown below, the histograms are highly dispersed around 1.


<div style="text-align:center;">
    <img src="https://i.imgur.com/BcgfKwH.png" width="600">
</div>


* This is characterized by the **Marcenko-Pastur Law**, an asymptotic distribution.
    * Under some mild moment conditions, this theory predicts convergence to a strictly positive density supported on the interval $[t_{min}(\alpha), t_{max}(\alpha)]$ where $$t_{min}(\alpha) := (1-\sqrt{\alpha})^2 \text{ and  } t_{max}(\alpha) := (1+\sqrt{\alpha})^2$$

* **[Non-asymptotic]** We seek results that hold for *all* choices of the pair $(n,d)$, ie. more in Chapter 6, the setting above can be shown that the max eigenvalue $\gamma_{max} (\hat{\Sigma})$ satisfies the upper deviation inequality $$P[\gamma_{max} (\hat{\Sigma}) \geq (1+\sqrt{d/n} + \delta)^2] \leq e^{-n\delta^2/2} \, \text{ for all } \delta \geq 0,$$ with an analogous lower deviation inequality for $\gamma_{min} (\hat{\Sigma})$ in the regime $n \geq d$.
* This gives more refined information about the max eigenvalue, showing 
    1. the probability that it deviates above $(1+\sqrt{d/n} + \delta)^2$ is exponentially small in the sample size $n$.
    2. $\hat{\Sigma}$ is an operator-norm-consistent estimate of the population matrix $\Sigma$ as long as $d/n \rightarrow 0$


## 2.1 Classical bounds
### 2.1.1 From *Markov* to *Chernoff*

### 2.1.2 *Sub-Gaussian variables* and *Hoeffding bounds*