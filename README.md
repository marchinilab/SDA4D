# SDA4D

An R package providing an implementation of a variational Bayes algorithm for 
four dimensional Parallel Factor Analysis.

The model is as follows for a data tensor $Y = (y_{nlmt})_{n,l,m,t}$:

### The Likelihood

$$
\begin{aligned}
\mathbb{P}(Y\vert\Theta) &= \prod_{n,l,m,t}N\left(y_{nlmt} \lvert \sum_{c=1}^C a_{nc}b_{tc}d_{mc}w_{cl}s_{cl},\lambda_{lt}^{-1}\right). 
\end{aligned}
$$

### The Prior

$$
\begin{aligned}
\mathbb{P}(\lambda_{lt}) &= \mbox{Gamma}(\lambda_{lt}\vert u,v);\\
\mathbb{P}(a_{nc}) &= N(a_{nc}\vert 0,1);\\
\mathbb{P}(b_{tc}) &= N(b_{tc}\vert 0,1);\\
\mathbb{P}(d_{mc}) &= N(d_{mc}\vert 0,1);\\
\mathbb{P}(w_{cl}\vert \beta_c) &= N(w_{cl}\vert 0,\beta_c^{-1});\\
\mathbb{P}(s_{cl}\vert \phi_{cl},\psi_{cl}) &= \mbox{Beta}(s_{cl}\vert \psi_{cl}\phi_{cl});\\
\mathbb{P}(\beta_c) &= \mbox{Gamma}(\beta_c\vert e,f);\\
\mathbb{P}(\phi_{cl}\vert\rho_c) &= \mbox{Bernoulli}(\phi_{cl}\vert \rho_c);\\
\mathbb{P}(\psi_{cl}) &= \mbox{Beta}(\psi_{cl}\vert g,h);\\
\mathbb{P}(\rho_c) &= \mbox{Beta}(\rho_c\vert r,z). 
\end{aligned}
$$

We apply a structured mean field approach and partition the variables as follows:

 $$\{a_{nc}\},\{b_{tc}\},\{d_{mc}\},\{\phi_{cl}\},\{\psi_{cl}\},\{\rho_c\},\{\beta_c\},\{w_{cl},s_{cl}\}$$