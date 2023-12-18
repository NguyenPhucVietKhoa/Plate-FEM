# Elément de plaque de Kirchhoff-Love {#KC}

Le vecteur des déformations est décrit par :
$$\boldsymbol{\epsilon}^h = \begin{bmatrix}
        \epsilon_{xx} \\
        \epsilon_{yy} \\
        \epsilon_{zz} \\
        \gamma_{xz} \\
        \gamma_{yz} \\
        \gamma_{xy} \\
    \end{bmatrix} = \begin{bmatrix}
        \frac{\partial u}{\partial x} \\
        \frac{\partial v}{\partial y} \\
        \frac{\partial w}{\partial z} \\
        \frac{\partial u}{\partial z} + \frac{\partial w}{\partial x} \\
        \frac{\partial v}{\partial z} + \frac{\partial w}{\partial y} \\
        \frac{\partial u}{\partial y} + \frac{\partial v}{\partial x} \\
    \end{bmatrix} = \begin{bmatrix}
        -z \, \frac{\partial^2 w^h}{\partial x^2} \\
        -z \, \frac{\partial^2 w^h}{\partial y^2} \\
         0 \\
         0 \\
         0 \\
         -2z\,\frac{\partial^2 w^h}{\partial x \partial y}\\
    \end{bmatrix}$$

On peut réécrire ce vecteur: $$\boldsymbol{\epsilon}^h = \begin{bmatrix}
        \epsilon_{xx} \\
        \epsilon_{yy} \\
        \gamma_{xy} \\ 
    \end{bmatrix} = z \begin{bmatrix}
        -\frac{\partial^2}{\partial x^2} \\
        -\frac{\partial^2}{\partial y^2} \\
        -2\frac{\partial^2}{\partial x \partial y}\\
    \end{bmatrix}w^h = z {\bf D}w^h$$

Le vecteur des contraintes est donc :
$\boldsymbol{\sigma}^h = z{\bf H}_f^* {\bf D}w^h$. La formulation par
éléments finis de $w^h$ est :
$$w^h = \sum\limits_{i=1}^{n^e} N_i(\xi, \eta) u_i^e = \begin{bmatrix}
        N_1 & N_2 & \ldots & N_{n^e} \\
    \end{bmatrix}\begin{bmatrix}
         w_1^e \\ {\theta_x}_1^e \\ {\theta_y}_1^e \\ \vdots \\ w_{n^e}^e \\ {\theta_x}_{n^e}^e \\ {\theta_y}_{n^e}^e \\
    \end{bmatrix} =  {\bf N}{\bf u^e}$$

La matrice de rigidité est donc: $$\begin{split}
        \int\limits_{\Omega} (\delta{\bf w}^h)^T {\bf D}^T {\bf H} {\bf D} {\bf w}^h \, dV & = \int\limits_{\Omega} (\delta{\bf u}^e)^T({\bf N}^T{\bf D}^T){\bf H}_f^*({\bf D} {\bf N}){\bf u}^e z^2\, dz\,dS \\
        & = (\delta{\bf u}^e)^T \int\limits_{S}{\bf P}^T{\bf H}_f^*{\bf P}\,dS \int\limits_{-h/2}^{h/2}z^2\,dz\, {\bf u}^e \\
        & = (\delta{\bf u}^e)^T \, \frac{h^3}{12}\int\limits_{S}{\bf P}^T{\bf H}_f^*{\bf P}\,dS\,{\bf u}^e \\
    \end{split}$$

Donc on déduit:
$${\bf K^e} = h\int\limits_{-1}^{1}\int\limits_{-1}^{1}{\bf P}^T{\bf H}_f^e{\bf P}J^e \,d\xi\,d\eta$$

Avec $${\bf P} = \begin{bmatrix}
        -\frac{\partial^2 {\bf N}}{\partial x^2} \\
        -\frac{\partial^2 {\bf N}}{\partial y^2} \\
        -2\frac{\partial^2 {\bf N}}{\partial x \partial y}\\
    \end{bmatrix}$$

Le calcul des dérivées du second ordre s'effectue de la manière
suivante: $$\begin{aligned}
    \frac{\partial^2}{\partial x^2} &= \left(\frac{\partial \xi}{\partial x}\right)^2 \frac{\partial^2}{\partial \xi^2} + \left(\frac{\partial \eta}{\partial x}\right)^2 \frac{\partial^2}{\partial \eta^2} + 2\left(\frac{\partial \xi}{\partial x}\frac{\partial \eta}{\partial x}\right)\frac{\partial^2}{\partial \xi \partial \eta} + \left(\frac{\partial^2 \xi}{\partial x^2}\right)\frac{\partial}{\partial \xi} + \left(\frac{\partial^2 \eta}{\partial x^2}\right)\frac{\partial}{\partial \eta} \\
     \frac{\partial^2}{\partial y^2} &= \left(\frac{\partial \xi}{\partial y}\right)^2 \frac{\partial^2}{\partial \xi^2} + \left(\frac{\partial \eta}{\partial y}\right)^2 \frac{\partial^2}{\partial \eta^2} + 2\left(\frac{\partial \xi}{\partial y}\frac{\partial \eta}{\partial y}\right)\frac{\partial^2}{\partial \xi \partial \eta} + \left(\frac{\partial^2 \xi}{\partial y^2}\right)\frac{\partial}{\partial \xi} + \left(\frac{\partial^2 \eta}{\partial y^2}\right)\frac{\partial}{\partial \eta} \\
      \frac{\partial^2}{\partial x \partial y} &= \left(\frac{\partial \xi}{\partial x}\frac{\partial \xi}{\partial y}\right) \frac{\partial^2}{\partial \xi^2} + \left(\frac{\partial \eta}{\partial x}\frac{\partial \eta}{\partial y}\right)\frac{\partial^2}{\partial \eta^2} + \left(\frac{\partial \xi}{\partial x}\frac{\partial \eta}{\partial y} + \frac{\partial \xi}{\partial y}\frac{\partial \eta}{\partial x}\right)\frac{\partial^2}{\partial \xi \partial \eta} + \left(\frac{\partial^2 \xi}{\partial x \partial y}\right)\frac{\partial}{\partial \xi} + \left(\frac{\partial^2 \eta}{\partial x \partial y}\right)\frac{\partial}{\partial \eta}
\end{aligned}$$

Dans nos calculs, on maille avec des rectangles de côtés $l^e_x$ et
$l^e_y$ 4 noeuds. Les fonctions d'interpolation sont suivante [@cours]
(ces fonctions sont utilisés pour interpoler les coordonnées physiques
seulement. Pour les déplacements, les fonctions d'interpolation sont
décrits plus tard): $$\begin{aligned}
    N_1^e &= \frac{1}{4}(1-\xi)(1-\eta) \\
    N_2^e &= \frac{1}{4}(1+\xi)(1-\eta) \\
    N_3^e &= \frac{1}{4}(1+\xi)(1+\eta) \\
    N_4^e &= \frac{1}{4}(1-\xi)(1+\eta)
\end{aligned}$$

La matrice Jacobienne se simplifie comme suivant:
$${\bf J}^e = \begin{bmatrix}
        \frac{\partial x}{\partial \xi} & \frac{\partial x}{\partial \eta} \\
        \frac{\partial y}{\partial \xi} & \frac{\partial y}{\partial \eta} \\
    \end{bmatrix} = \begin{bmatrix}
        \frac{l^e_x}{2} & 0 \\
        0 & \frac{l^e_y}{2} \\
    \end{bmatrix} \Rightarrow {{\bf J}^e}^{-T} = \begin{bmatrix}
        \frac{\partial \xi}{\partial x} & \frac{\partial \eta}{\partial x} \\
        \frac{\partial \xi}{\partial y} & \frac{\partial \eta}{\partial y} \\
    \end{bmatrix} = \begin{bmatrix}
        \frac{2}{l^e_x} & 0 \\
        0 & \frac{2}{l^e_y} \\
    \end{bmatrix}$$

Donc le calcul des dérivées du second ordre se simplifie :
$$\begin{aligned}
    \frac{\partial^2}{\partial x^2} &= \left(\frac{2}{l^e_x}\right)^2 \frac{\partial^2}{\partial \xi^2} \\
    \frac{\partial^2}{\partial y^2} &= \left(\frac{2}{l^e_y}\right)^2 \frac{\partial^2}{\partial \eta^2} \\
    \frac{\partial^2}{\partial x \partial y} &= \left(\frac{4}{l^e_x\,l^e_y} \right) \frac{\partial^2}{\partial \xi \partial \eta} \\
\end{aligned}$$

L'expression du vecteur de force thermique: $$\begin{split}
        {\bf F}_{th}^e &=  \int\limits_{\Omega} {\bf N}^T{D}^T {\bf H} \boldsymbol{\epsilon}_{th} \, dV \\
    &= \int\limits_{\Omega} z{\bf P}^T {\bf H}_f^* \, z \alpha \frac{\Delta T_z}{h} \begin{bmatrix}
        1 & 1 & 0 \\
    \end{bmatrix}^T \, dz\,dS \\
    &= \int\limits_{S} {\bf P}^T {\bf H}_f^* \, \alpha \frac{\Delta T_z}{h} \begin{bmatrix}
        1 & 1 & 0 \\
    \end{bmatrix}^T\, dS \int\limits_{-h/2}^{h/2}z^2 \, dz \\
    & = \int\limits_{-1}^{1}\int\limits_{-1}^{1}{\bf P}^T{\bf H}_f^e \,  \alpha \Delta T_z \begin{bmatrix}
        1 & 1 & 0 \\
    \end{bmatrix}^T\, J^e\,d\xi \, d\eta \\
    \end{split}$$

Les fonctions d'interpolation pour un élément de plaque Kirchhoff-Love
(4 noeuds et 3 DDLs/noeud) [@ferreira2009matlab]: $$\begin{aligned}
    N_1^e &= -\frac{1}{8}(\xi-1)(\eta-1)(\xi^2 + \xi + \eta^2 + \eta -2) \\
    N_2^e &= -\frac{l^e_y}{16}(\xi-1)(\eta-1)^2(\eta+1) \\
    N_3^e &= \frac{l^e_x}{16}(\xi-1)^2(\xi+1)(\eta-1) \\
    N_4^e &= \frac{1}{8}(\xi+1)(\eta-1)(\xi^2 - \xi + \eta^2 + \eta -2) \\
    N_5^e &= \frac{l^e_y}{16}(\xi+1)(\eta-1)^2(\eta+1) \\
    N_6^e &= \frac{l^e_x}{16}(\xi-1)(\xi+1)^2(\eta-1) \\
    N_7^e &= \frac{1}{8}(\xi+1)(\eta+1)(-\xi^2 + \xi - \eta^2 + \eta +2) \\
    N_8^e &= \frac{l^e_y}{16}(\xi+1)(\eta-1)(\eta+1)^2 \\
    N_9^e &= -\frac{l^e_x}{16}(\xi-1)(\xi+1)^2(\eta-1) \\
    N_{10}^e &= \frac{1}{8}(\xi-1)(\eta+1)(\xi^2 + \xi + \eta^2 - \eta -2) \\
    N_{11}^e &= -\frac{l^e_y}{16}(\xi-1)(\eta-1)(\eta+1)^2 \\
    N_{12}^e &= -\frac{l^e_x}{16}(\xi-1)^2(\xi+1)(\eta+1) \\
    \Rightarrow {\bf N} &= \left[N_1^e \; N_2^e \; N_3^e \; N_4^e \; N_5^e \; N_6^e \; N_7^e \; N_8^e \; N_9^e \; N_{10}^e \; N_{11}^e \; N_{12}^e\right] 
\end{aligned}$$
