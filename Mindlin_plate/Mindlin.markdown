# Elément de plaque de Reissner-Mindlin

Le champ de déplacements interpolé est :
$$w^h = \sum\limits_{i=1}^{n^e} N_i(\xi, \eta) w_i^e \; , \; \theta_x^h = \sum\limits_{i=1}^{n^e} N_i(\xi, \eta) {\theta_x}_i^e \; , \; \theta_y^h = \sum\limits_{i=1}^{n^e} N_i(\xi, \eta) {\theta_y}_i^e$$

Le vecteur des déplacement est donc: $${\bf u}^h = \begin{bmatrix}
        N_1 & 0 & 0 & N_2 & 0 & 0 & \ldots & N_{n^e} & 0 & 0 \\
        0 & N_1 & 0 & 0 & N_2 & 0 & \ldots & 0 & N_{n^e} & 0 \\
        0 & 0 & N_1 & 0 & 0 & N_2 & \ldots & 0 & 0 & N_{n^e} \\
    \end{bmatrix} \begin{bmatrix}
        w_1^e \\ {\theta_x}_1^e \\ {\theta_y}_1^e \\ \vdots \\ w_{n^e}^e \\ {\theta_x}_{n^e}^e \\ {\theta_y}_{n^e}^e \\
    \end{bmatrix} = {\bf N}{\bf u}^e$$

Donc, on a le vecteur des déformations de flexion:
$$\boldsymbol{\epsilon}_f = \begin{bmatrix}
        \epsilon_{xx} \\
        \epsilon_{yy} \\
        \gamma_{xy} \\
    \end{bmatrix} = z \begin{bmatrix}
        \frac{\partial \theta_y}{\partial x} \\
        -\frac{\partial \theta_x}{\partial y} \\
        -\frac{\partial \theta_x}{\partial x} + \frac{\partial \theta_y}{\partial y}\\
    \end{bmatrix} = z \begin{bmatrix}
        0 & 0 & \frac{\partial}{\partial_x} \\
        0 & -\frac{\partial}{\partial_y} & 0\\
        0 & -\frac{\partial}{\partial_x} & -\frac{\partial}{\partial_y} \\
    \end{bmatrix}\begin{bmatrix}
        w \\ \theta_x \\ \theta_y \\
    \end{bmatrix} = z {\bf D}_f{\bf u}^h$$

Le vecteur des déformations en cisaillement:
$$\boldsymbol{\epsilon}_s = \begin{bmatrix}
        \gamma_{xz} \\
        \gamma_{yz} \\
    \end{bmatrix} = \begin{bmatrix}
        \frac{\partial w}{\partial x} + \theta_y \\
        \frac{\partial w}{\partial y} - \theta_x \\
    \end{bmatrix} = \begin{bmatrix}
        \frac{\partial}{\partial_x} & 0 & 1\\
        \frac{\partial}{\partial_y} & -1 & 0 \\
    \end{bmatrix}\begin{bmatrix}
         w \\ \theta_x \\ \theta_y \\
    \end{bmatrix}={\bf D}_s {\bf u}^h$$

Donc on trouve les vecteur contraintes de flexion et en cisaillement
correspondants: $$\begin{aligned}
    \boldsymbol{\sigma}_f &= {\bf H}_f^* \boldsymbol{\epsilon}_f = z {\bf H}_f^* {\bf D}_f{\bf u}^h \\
    \boldsymbol{\sigma}_s &= {\bf H}_s^* \boldsymbol{\epsilon}_s = {\bf H}_s^* {\bf D}_s{\bf u}^h
\end{aligned}$$

Avec

$${\bf H}_f^* = \frac{E}{1-\nu^2}\begin{bmatrix}
        1 & \nu & 0 \\
        \nu & 1 & 0 \\
        0 & 0 & \frac{1-\nu}{2} \\
    \end{bmatrix} \quad \hbox{et} \quad {\bf H}_s^* = \frac{E \kappa}{2(1+\nu)}\begin{bmatrix}
        1 & 0 \\
        0 & 1 \\
    \end{bmatrix}$$

On peut déduire que : $$\begin{gathered}
    \int\limits_{\Omega} (\delta{\bf u}^h)^T {\bf D}^T {\bf H} {\bf D} {\bf u}^h \, dV =  \int\limits_{\Omega} (\delta{\bf u}^e)^T({\bf N}^T{\bf D}_f^T){\bf H}_f^*({\bf D}_f {\bf N}){\bf u}^e z^2\, dz\,dS \\
    + \int\limits_{\Omega} (\delta{\bf u}^e)^T({\bf N}^T{\bf D}_s^T){\bf H}_s^*({\bf D}_s {\bf N}){\bf u}^e\, dz\,dS \\
    =  (\delta{\bf u}^e)^T \left(\int\limits_{\Omega} ({\bf N}^T{\bf D}_f^T){\bf H}_f^*({\bf D}_f {\bf N}) z^2\, dz\,dS + \int\limits_{\Omega}({\bf N}^T{\bf D}_s^T){\bf H}_s^*({\bf D}_s {\bf N})\, dz\,dS\right){\bf u}^e
\end{gathered}$$

On défini les matrices ${\bf P}_f$ et ${\bf P}_s$ par :
$$\begin{aligned}
    {\bf P}_f &= {\bf D}_f{\bf N} = \begin{bmatrix}
        0 & 0 & \frac{\partial N_1}{\partial x} & 0 & 0 & \frac{\partial N_2}{\partial x} & \ldots & 0 & 0 & \frac{\partial N_{n^e}}{\partial x} \\
        0 & -\frac{\partial N_1}{\partial y} & 0 & 0 & -\frac{\partial N_2}{\partial y} & 0 & \ldots & 0 & -\frac{\partial N_{n^e}}{\partial y} & 0 \\
        0 & -\frac{\partial N_1}{\partial x} & \frac{\partial N_1}{\partial y} & 0 & -\frac{\partial N_2}{\partial x} & \frac{\partial N_2}{\partial y} & \ldots & 0 & -\frac{\partial N_{n^e}}{\partial x} & \frac{\partial N_{n^e}}{\partial y} \\
    \end{bmatrix} \\
    {\bf P}_s &= {\bf D}_s{\bf N} = \begin{bmatrix}
        \frac{\partial N_1}{\partial x} & 0 & N_1 & \frac{\partial N_2}{\partial x} & 0 & N_2 & \ldots & \frac{\partial N_{n^e}}{\partial x} & 0 & N_{n^e} \\
        \frac{\partial N_1}{\partial y} & -N_1 & 0 & \frac{\partial N_2}{\partial y} & -N_2 & 0 & \ldots & \frac{\partial N_{n^e}}{\partial y} & -N_{n^e} & 0 \\
    \end{bmatrix}
\end{aligned}$$

Dans le cas 2D, on a l'expression de la matrice Jacobienne pour passer
entre les coordonnées paramétriques $(\xi, \eta)$ et les coordonnées
physiques $(x,y)$ est: $${\bf J}^e (\xi, \eta) = \begin{bmatrix}
        \frac{\partial x}{\partial \xi} & \frac{\partial x}{\partial \eta} \\
        \frac{\partial y}{\partial \xi} & \frac{\partial y}{\partial \eta} \\
    \end{bmatrix} = \begin{bmatrix}
        \sum\limits_{i=1}^{n^e}\frac{\partial N_i(\xi, \eta)}{\partial \xi}x^e_i & \sum\limits_{i=1}^{n^e}\frac{\partial N_i(\xi, \eta)}{\partial \eta}x^e_i \\
        \sum\limits_{i=1}^{n^e}\frac{\partial N_i(\xi, \eta)}{\partial \xi}y^e_i & \sum\limits_{i=1}^{n^e}\frac{\partial N_i(\xi, \eta)}{\partial \eta}y^e_i \\
    \end{bmatrix}$$

Les dérivées $\frac{\partial}{\partial x}$ et
$\frac{\partial}{\partial y}$ sont calculés par les expression suivants:
$$\begin{bmatrix}
        \frac{\partial}{\partial x} \\
        \frac{\partial}{\partial y} \\
    \end{bmatrix}={\bf J}^e (\xi, \eta)^{-T} \begin{bmatrix}
        \frac{\partial}{\partial \xi} \\
        \frac{\partial}{\partial \eta} \\
    \end{bmatrix}$$

On a ${\bf K}^e = {\bf K}_f^e + {\bf K}_s^e$, avec: $$\begin{aligned}
    {\bf K}_f^e &= \int\limits_{S}{\bf P}_f^T{\bf H}_f^*{\bf P}_f\,dS \int\limits_{-h/2}^{h/2}z^2\,dz =  \frac{h^3}{12}\int\limits_{S}{\bf P}_f^T{\bf H}_f^*{\bf B}_f\,dS \\
    {\bf K}_s^e &= \int\limits_{S}{\bf P}_s^T{\bf H}_s^*{\bf P}_s\,dS \int\limits_{-h/2}^{h/2}\,dz = h \int\limits_{S}{\bf P}_s^T{\bf H}_s^*{\bf P}_s\,dS
\end{aligned}$$

Aux coordonnées paramétriques $(\xi, \eta)$: $$\begin{aligned}
    {\bf K}_f^e &= h \int\limits_{-1}^{1}\int\limits_{-1}^{1}{\bf P}_f^T{\bf H}_f^e{\bf P}_f\,J^e \,d\xi \, d\eta \\
    {\bf K}_s^e &= h \int\limits_{-1}^{1}\int\limits_{-1}^{1}{\bf P}_s^T{\bf H}_s^e{\bf P}_s\, J^e \,d\xi \, d\eta \\
\end{aligned}$$ (car on a $dS = J^e\,d\xi\,d\eta$)

où $${\bf H}_f^e = \frac{Eh^2}{12(1-\nu^2)}\begin{bmatrix}
        1 & \nu & 0 \\
        \nu & 1 & 0 \\
        0 & 0 & \frac{1-\nu}{2} \\
    \end{bmatrix} \; et \; {\bf H}_s^e = \frac{E \kappa}{2(1+\nu)}\begin{bmatrix}
        1 & 0 \\
        0 & 1 \\
    \end{bmatrix}$$

Similaire, on peut trouver l'expression pour le chargement thermique. On
a le vecteur des déformations thermiques défini par:

$$\boldsymbol{\epsilon}_{th} = z \alpha \frac{\Delta T_z}{h}\begin{bmatrix}
        1 & 1 & 0 \\
    \end{bmatrix}^T
    \label{defor thermique}$$

Donc, on a: $$\begin{split}
        {\bf F}_{th}^e &=  \int\limits_{\Omega} {\bf N}^T{D}^T {\bf H} \boldsymbol{\epsilon}_{th} \, dV \\
    &= \int\limits_{\Omega} z{\bf P_f}^T {\bf H}_f^* \, z \alpha \frac{\Delta T_z}{h} \begin{bmatrix}
        1 & 1 & 0 \\
    \end{bmatrix}^T \, dz\,dS \\
    &= \int\limits_{S} {\bf P_f}^T {\bf H}_f^* \, \alpha \frac{\Delta T_z}{h} \begin{bmatrix}
        1 & 1 & 0 \\
    \end{bmatrix}^T\, dS \int\limits_{-h/2}^{h/2}z^2 \, dz \\
    & = \int\limits_{-1}^{1}\int\limits_{-1}^{1}{\bf P_f}^T{\bf H}_f^e \,  \alpha \Delta T_z \begin{bmatrix}
        1 & 1 & 0 \\
    \end{bmatrix}^T\, J^e\,d\xi \, d\eta \\
    \end{split}$$

Pour un élément de plaque Mindlin Serendip (un quadrilatère à 8 noeuds
et 3 DDLs par noeud), les fonctions d'interpolation son données par
[@cours]: $$\begin{aligned}
    N_1^e &= -\frac{1}{4}(1-\xi)(1-\eta)(1+\xi+\eta) \\
    N_2^e &= -\frac{1}{4}(1+\xi)(1-\eta)(1-\xi+\eta) \\
    N_3^e &= -\frac{1}{4}(1+\xi)(1+\eta)(1-\xi-\eta) \\
    N_4^e &= -\frac{1}{4}(1-\xi)(1+\eta)(1+\xi-\eta) \\
    N_5^e &= \frac{1}{2}(1-\xi)(1+\xi)(1-\eta) \\
    N_6^e &= \frac{1}{2}(1+\xi)(1-\eta)(1+\eta) \\
    N_7^e &= \frac{1}{2}(1-\xi)(1+\xi)(1+\eta) \\
    N_8^e &= \frac{1}{2}(1-\xi)(1-\eta)(1+\eta) \\
\end{aligned}$$

Pour économiser le charge de calcul, on va utiliser ces fonctions
d'interpolation pour interpoler le champ des déplacements ainsi que la
coordonnée physique.
