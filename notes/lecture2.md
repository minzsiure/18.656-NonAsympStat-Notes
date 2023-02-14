# 18.656 Lecture 2 Reading Notes: Sub-Gaussian random variables and tail bounds
###### tags: `18.656`
```
Section 2.1, Basic definitions and results: 
* Sub-Gaussian variables.  
* Hoeffding bounds.  
* Symmetrization.  
* Sub-exponential tail bounds.

Motivating uses:  
* testing multiple hypotheses (discussed in class) at length.
```
## 1. Sub-Gaussian tail bounds *(Chaper 2.1.2)*

### A. Gaussian tail bounds
Here we illustrate the use of the Chernoff bound in deriving tail bounds for a Gaussian variable.

Let $X \sim N(\mu, \sigma^2)$ be a Gaussian R.V. with mean $\mu$ and variance $\sigma^2$. We find that $X$ has the moment generating function $$E[e^{\lambda X}] = e^{\mu \lambda + \frac{\sigma^2 \lambda^2}{2}}, \text{valid for all $\lambda \in R$}$$

Substitute this into the optimization problem definining the optimized Chernoff bound, we have $$\inf_{\lambda \geq 0}\{log E[e^{\lambda (X - \mu)}] - \lambda t\} = \inf_{\lambda \geq 0}\{\frac{\lambda^2 \sigma^2}{2} - \lambda t\} = -\frac{t^2}{2 \sigma^2}$$

where we have taken derivatives to find the optimum of this quadratic function. Returning to the Chernoff bound, we conclude any $N(\mu, \sigma^2)$ R.V. satisfied the *upper deviation inequality* $$P[X \geq \mu + t] \leq e^{-\frac{t^2}{2 \sigma^2}} \text{ for all $t \geq 0$.} \, \text{ (2.7)}$$

![](https://i.imgur.com/jAU7cgT.png)
* The constant $\sigma$ is the *sub-Gaussian parameter*. We say $X$ is *sub-Gaussian* with parameter $\sigma$ when (2.8) holds. 
* Any Gaussian variable with variance $\sigma^2$ is *sub-Gaussian* with parameter $\sigma$. 
    * A large number of non-Gaussian R.V. also satisfy the condition (2.8).
* Condition (2.8) with Chernoff bound shows that if $X$ is sub-Gaussian  with parameter $\sigma$, it satisfies the *upper deviation inequality* (2.7).
* By symmetry of the definition, $-X$ is sub-Gaussian iif $X$ is sub-Gaussian, so we also have *lower deviation inequality* $$P[X \leq \mu + t] \leq e^{-\frac{t^2}{2 \sigma^2}} \text{ for all $t \geq 0$.}$$

**Concentation Inequality**
* Combined, any sub-Gaussian variables satisfies the *concentation inequality* $$P[|X - \mu| \geq t] \leq 2e^{-\frac{t^2}{2 \sigma^2}} \text{ for all $t \in R$.} \, \text{  (2.9)}$$

> Let's consider sub-Gaussian variables that are non-Gaussian.
### B) Rademacher variables
> **Defn:** a Rademacher R.V. $\epsilon$ takes the value $\{-1, +1\}$ equiprobably. 
> It is sub-Gaussian with parameter $\sigma = 1$.

A Rademacher variable $\epsilon \in \{-1, +1\}$ with equal probability such that 
$$P(\epsilon = +1) = P(\epsilon = -1) = \frac 1 2$$

Then, by taking expectations and play with Taylor expansions, we have that $\forall \lambda \in R$

$$\begin{align*}
    \mathbb{E}[e^{\lambda \epsilon}] 
    &= \frac 1 2 (e^\lambda + e^{-\lambda}) \\
    &= \frac 1 2 \left\{\sum_{k=0}^\infty \frac{(-\lambda)^k}{k!} +\sum_{k=0}^\infty \frac{\lambda^k}{k!} \right\} \\
    &= \sum_{k=0}^\infty \frac {\lambda^{2k}}{(2k)!} 
    \text{(Only the even powers remain)}\\
    &\leq 1+\sum_{k=1}^\infty \frac{\lambda^{2k}}{2^k k!}\\
    &= e^{\frac {\lambda^2} 2}
\end{align*}$$ 

Therefore, $\epsilon$ is sub-Gaussian with parameter $\sigma=1$.

> Generalized, any bounded R.V. is also sub-Gaussian
### C) Bounded random variables
Let $X$ be zero-mean on some bounded interval $[a, b]$.

> We use the technique of **Symmetrization** in here.

Letting $X'$ be an independent copy (_Symmetry_), for any $\lambda \in \mathbb{R}$, we have

$$\mathbb{E}_X[e^{\lambda X}] = \mathbb{E}_X[e^{\lambda(X-E_X'[X'])}] \leq \mathbb{E}_{X, X'} [e^{\lambda(X-X')}]$$ (by _Jensen's inequality_)
> Reminder: **Jensen's Inequality**
    $$\phi(\mathbb{E}[X]) \leq \mathbb{E}[\phi(X)]$$ for convex function $\phi$.
    
Here the inequality holds by convexity of exponential functions.

Let $\epsilon$ be an independent Rademacher variable, note that $\epsilon(X-X')$ has the same distribution as $X-X'$. 
> This is because $X$ and $X'$ are from the same distribution, and are independent, so the distribution of $X - X'$ is the same as that of $\epsilon(X' - X$). 

Then, we have
$$ \mathbb{E}_{X, X'} [e^{\lambda(X-X')}] = \mathbb{E}_{X, X'} [E_{\epsilon}[e^{\lambda\epsilon(X-X')}]]$$ 

If condition on $X$, $X'$, we could evaluate the inner expectation. By the previous conclusion, 
$$E_{\epsilon}[e^{\lambda\epsilon(X-X')}] \leq e^{\frac{\lambda^2(X-X')^2}{2}}$$

and since $|X-X'| < b-a$ we have that 

$$\mathbb{E}_{X, X'} [E_{\epsilon}[e^{\lambda\epsilon(X-X')}]] \leq E_{X, X'} [e^{\frac{\lambda^2(X-X')^2}{2}}] \leq e^{\frac {\lambda^2}{2}(b-a)^2}$$

$X$ is $(b-a)$-sub-Gaussian.


> **Is this the "best" bound we can get?**
No! We can actually improve it to $\sigma = \frac{(b-a)}{2}$. For details, see exercise 2.4.

## 2. Hoeffding bounds *(Chaper 2.1.2)*
> The property of Gaussianity is preserved by linear operations, so is the property of sub-Gaussianity! 
If $X_1$ and $X_2$ are indepedent sub-Gaussian variables with parameters $\sigma_1$ and $\sigma_2$, then $X_1 + X_2$ is sub-Gaussian with parameter $\sqrt{\sigma_1^2 + \sigma_2^2}.$

> *Hoeffding bounds* is applicable to sums of independent sub-Gaussian R.V.

![](https://i.imgur.com/fAMf4E4.png)

The general technique we use is to bound the MGF and apply Chernoff bound.

***Proof.***
> Note: This proof is not stated in the textbook.

$$\mathbb{E}[e^{\lambda\sum_{i=1}^n(X_i-\mu_i)}] = \mathbb{E}[\Pi_{i=1}^n e^{\lambda(X_i-\mu_i)}] = \Pi_{i=1}^n \mathbb{E}[e^{\lambda(X_i-\mu_i)}]$$  
>The first equality is by the property of exponentials, and the second equality is by the independence of $X_i$'s.

Each $X_i$ is $\sigma_i$-sub-Gaussian so we have

$$\Pi_{i=1}^n \mathbb{E}[e^{\lambda(X_i-\mu_i)}] \leq \Pi_{i=1}^n e^{\frac{\lambda^2\sigma_i^2} 2} = e^{\frac {\lambda^2}{2}(\sum_{i=1}^n \sigma_i^2)}$$

The key take-away is that for $\sigma_i$-sub-Gaussian variables $X_i$, 
$\sum_{i=1}^n X_i$ is another sub-Gaussian variable with parameter $\sqrt{\sum_{i=1}^n \sigma_i^2}$.

> We could use independence to manipulate MGFs like this.

After this, complete the proof by applying the Chernoff techinque from $e^{\frac {\lambda^2}{2}(\sum_{i=1}^n \sigma_i^2)}$. 
He skipped as an exercise.

By Chernoff's bound, $P((X-\mu) \geq t) = P(e^{\lambda(X-\mu)} \geq e^{\lambda t}) \leq \frac{E[e^{\lambda(X-\mu)}]}{e^{\lambda t}}$
$$\mathbb{P}[\sum_{i=1}^n (X_i-\mu_i) \geq t] \leq \frac{E[e^{\lambda\sum_{i=1}^n (X_i-\mu_i)}]}{e^{\lambda t}} \leq e^{\frac {\lambda^2}{2}(\sum_{i=1}^n \sigma_i^2) - \lambda t} \leq e^\frac{-t^2}{2 \sum_{i=1}^n \sigma_i^2}$$

The last inequality comes from bounding the value $\frac {\lambda^2}{2}(\sum_{i=1}^n \sigma_i^2) - \lambda t$ *with respect to $\lambda$*. 
By taking derivative, the upper bound is optimized when $\lambda=\frac{t}{\sum_{i=1}^n \sigma_i^2}$. 
We then plug it back into the original expression and obtain the inequality we desired.


## 3. Application/Motivating Example: Multiple Testing *(Not in textbook)*
Suppose we have random variables $Y_i = \theta_i + W_i$. 
> $\theta_i$: unknown scalars
$W_i$: independent, zero-mean, and $\sigma$-sub-Gaussian.


The problem is we would like to differentiate these two scenarios:
- $H_0$: $\theta_i=0$ for all i.
- $H_1$: At least one such $j \in \{1, 2, ..., n\}, \theta_j \geq t > 0$.

> - Useful for drug companies and treatment, experiments and etc. They see outcome and try to judge if the drugs are trash or t-effective.
> - Another application is ads from large corps like Facebook/Google. They figure out products and click through rates. Martin doesn't like that 👎

**Class Question:** *What if we want to figure out which $j$ is actually effective and how effective?*
**Answer:** *Localization (find $j$) is harder. Estimation (how effective) is even harder.*


> **Idea**: Make the threshold high enough so that the noises won't be above.
To accomplish this, we want to estimate the *probability* that the noise is greater than the threshold ($t$). 
Essentially, <u>we want a high enough $t$ so that this probability is small enough</u>.

**Here is a reasonable testing procedure:** 

Declare $H_0$ if $\max_{i=1, .., n} Y_i \leq t$, $H_1$ otherwise. 


$P_{H_0}[\text{failing}] = P[\max_{i=1, ..., n} W_i \geq t] \leq \sum_{i=1}^n P(W_i \geq t) \leq \sum_{i=1}^n e^{\frac{-t^2}{2\sigma^2}} = ne^{\frac{-t^2}{2\sigma^2}}$

The conclusion is
$$P_{H_0}[\text{failing}] \leq ne^{\frac{-t^2}{2\sigma^2}}$$

In practice, if the failing rate is some $\delta \in (0, 1)$, we could set $ne^{\frac{-t^2}{2\sigma^2}} \leq \delta$ and solve for $t$.

By algebra, we see that $t=\sqrt{2\sigma^2 \log(\frac n \delta)}$ is the threshold we want to fail within rate $\delta$.


## 4. What is not sub-Gaussian? Sub-exponentials! *(Chaper 2.1.3)*
![](https://i.imgur.com/FlhaJ9I.png)
> It follows immediately that any sub-Gaussian variable is also sub-exponential (with $v = \sigma$ and $\alpha = 0$), but the converse is not true.

![](https://i.imgur.com/UTVVoS0.png)