# 18.656 Lecture 1 Reading Notes: Introduction
###### tags: `18.656`
```
Summary of lecture:

1. Discussion of course content, requirements, grading, and structure. 
* See the Syllabus page on Canvas for full details.

2. Motivation for non-asymptotic and/or high-dimensional analysis in statistics.  
* Relevant reading:  Chapter 1 of HighDimStat.  
* See Section 1.4 for more details on "three modes" of analysis.

3. More involved discussion using `covariance estimation` as a motivating example.  
* Section 1.2.2 for relevant details, 
* Figure 1.2 for Marcenko-Pastur law.  
* Section 1.3.2 for the possibility of structure in covariance matrices.

4. End of lecture: moved onto basic tail bounds in Chapter 2, 
* Partial Section 2.1.  
* Please read ahead of Chapter 2; 
we will cover large parts (but not all) of it, 
but it is all extremely relevant and useful.
```
[toc]

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
> <u>Overview:</u> Here, we anlyze the behavior of the engenvalues of  a sample covariance matrix $\hat{\Sigma}$ based on $n$ samples of a $d$-dimensional random vector with the identity matrix as its covariance.

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

## 1.3 What can help us in high dimension?
### 1.3.2 Structure in covariance matrices
> <u>Overview:</u> What types of ***low-dimensional structure*** might be appropriate for modeling covariance matrices? And how can they can be exploited to construct better estimators?

> The low-dimensional structure of matrices refers to the ways in which matrices can be characterized or represented using a smaller number of parameters or features. In other words, it refers to the finding of a low-dimensional representation of a matrix that captures its important characteristics. ie. SVD, matrix factorization.

* **Example 1:** if the goal is to estimate a *diagonal* $\Sigma$, $\hat{\Sigma}$ can be improved by zeroing out the non-diagonal entries, leading to a diagonal covariance estimate $\hat{D}$.
* **Example 2:** if $\Sigma$ were assumed to be *sparse* with unknown positions, we can have a hard-thresholded version of $\hat{\Sigma}$, $$\tilde{\Sigma}:=T_{\lambda_n}(\hat{\Sigma}), \text{ say with } \lambda_n = \sqrt{\frac{2 \, log \,d}{n}}$$
    * Figure 1.5(a)
        * We see the eigenspectrum $\gamma(\tilde{\Sigma})$ is sharply concentrated around the point mass at 1.
        * Tail bounds and theory can be used to show that, with high probability,  $$\||\hat{\Sigma}-\Sigma\||_2 \lesssim \sqrt{\frac{log \, d}{n}} $$ 

<div style="text-align:center;">
    <img src="https://i.imgur.com/bkGJahq.png" width="700">
</div>



* **Example 3:** symmetric matrices with fast decay in their eigenspectra. Suppose sequence of covariance matrices have a bounded trace -- that is, $trace(\Sigma) \leq R$, independent of $d$. Then the ordered $\gamma_j(\Sigma)$ must decay quicker than $j^{-1}$.
    * Figure 1.5(b)
        * We plot $\||\hat{\Sigma}-\Sigma\||_2$ over a range of pairs $(n,d)$ with fixed ratio $d/n = 0.2$, for a sequence of covariance matrices that satisfy the constraint $trace( \Sigma ) \leq 20$.
        * Theoretical results in Ch6 predict for such sequence of covariance matrices with a bounded trace -- that is, $trace( \Sigma ) \leq R$, independent of the dimension $d$ -- **the error  $\||\hat{\Sigma}-\Sigma\||_2$ should decay as $n^{-1/2}$ even if $d$ grows in proportion to the sample size $n$**.


## 2.1 Classical bounds
One way to control a tail probability $P[X \geq t]$ is by controlling the moments of R.V. $X$. This ranges from Markov's inequality (require $1^{st}$ moment) to the Chernoff bound (requirement moment generating function).

### 2.1.1 From *Markov* to *Chernoff*
#### A) Markov's inequality (upper tail bound)
* Given a non-negative random variable $X$ with *finite mean*, we have $$P[X \geq t] \leq \frac{E[X]}{t} \, \text{ for all } t> 0.$$

#### B) Chebyshev's inequality (concentration inequality)
* Given a random variable $X$ that also has a finite variance, we have $$P[|X-\mu |\geq t] \leq \frac{var(X)}{t^2} \, \text{ for all } t> 0.$$
* This guarantees that $X$ is close to its mean whenever the variance is small.

#### C) Markov's inequality to R.V. with higher-order moments
* Whenever $X$ has a central moment of order $k$, an application of Markov's inequality to the R.V. $|X - \mu |^k$ yields that $$P[|X-\mu |\geq t] \leq \frac{E[|X - \mu|^k]}{t^k} \, \text{ for all } t> 0.$$

#### D) Chernoff bound
$$P[(X - \mu) \geq t] = P[e^{\lambda (X - \mu)} \geq e^{\lambda t}] \leq \frac{E[e^{\lambda (X - \mu)}]}{e^{\lambda t}} $$ 

* Optimizing the choice of $\lambda$ to obtain the tightest result yields the *Chernoff bound* -- namely,

$$log \, P[(X - \mu) \geq t] \leq \inf_{\lambda \in [0,b]} \{log \, E[e^{\lambda (X - \mu)}] - \lambda t\}$$


<!-- ### 2.1.2 *Sub-Gaussian variables* and *Hoeffding bounds* -->
