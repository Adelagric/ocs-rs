# Sélection à contribution optimale exacte à l'échelle génomique : un solveur « support-first », sans matrice, validé contre optiSel et AlphaMate sur panels de marqueurs réels

**Adel Kaleche**

*Chercheur indépendant, France.*

*Version française ; le manuscrit de référence est en anglais (`MANUSCRIPT.md`).*


## Résumé

La sélection à contribution optimale (OCS, *optimum contribution selection*) — maximiser le gain génétique sous une borne sur la coancestrie de la génération suivante — est l'outil central pour arbitrer entre gain et perte de diversité dans les programmes de sélection. Son goulot d'étranglement computationnel est la matrice de parenté qu'elle contraint : à l'échelle génomique, cette matrice est dense et n×n, coûtant O(n²m) à construire et O(n²) à stocker, de sorte que les outils établis — exact (optiSel) et heuristique (AlphaMate) — heurtent un mur mémoire précisément quand le nombre de candidats croît. Nous présentons **support-first**, un solveur exact qui ne la forme jamais : sur un ordinateur portable, il résout une instance (n = 40 000) dont la matrice dense (11,9 Gio) dépasse ce que les outils établis peuvent allouer, le support optimal se maintenant autour de quinze individus — accomplissant ce qui est autrement infaisable, pas seulement plus vite. Il exploite deux faits structurels de l'OCS : le vecteur de contributions optimal a un support restreint — quelques dizaines d'individus aux bornes qu'emploient les programmes de sélection (19 à 28 sur n pour les panels réels ici), ne croissant qu'à mesure que la borne approche la coancestrie minimale — et la seule difficulté est combinatoire — quels candidats sont libres versus fixés à une borne (zéro, ou un plafond par candidat). Support-first fait croître le support par une règle d'ensemble actif / génération de colonnes et résout chaque sous-problème à support fixé — maximiser un objectif linéaire sur un ellipsoïde intersecté avec les contraintes affines de somme (et de sexe) — en forme fermée, sans solveur itératif interne ; les produits de parenté G·c sont formés sans matrice à partir de la matrice de génotypes Z, ce qui permet de ne jamais stocker la matrice dense n×n G (l'accélération d'un ordre de grandeur vient du petit ensemble actif, pas du produit sans matrice). Sur données réelles — un panel de blé CIMMYT (n = 599), un panel de porc PIC (n = 3534, 52k SNP) et un panel de souris en population hétérogène (n = 1814, avec sexe réel) — il atteint l'optimum exact sur la frontière de coancestrie, en accord à 10⁻⁸ près avec un solveur conique à point intérieur qui, lui, s'arrête juste en deçà de cette frontière. Il est aussi rapide : 13× à 182× plus vite qu'optiSel de bout en bout sans former aucun **G** (et, recevant **G** comme l'exigent les outils à point intérieur, jusqu'à ~2500×), et ~76000× plus vite à n = 10000 qu'un solveur conique généraliste, dans une comparaison Rust contre Rust exempte de tout biais de langage. Face à AlphaMate — une heuristique pour le problème distinct de l'allocation discrète des accouplements — l'optimum exact n'est pas moins bon, à coancestrie appariée, sur la relaxation continue commune aux deux méthodes ; les plafonds de contribution par candidat (0 ≤ c ≤ u) sont gérés. Mais la vitesse est le moindre des résultats. Parce qu'une résolution coûte le support et le nombre de marqueurs plutôt que la matrice n×n, support-first reste peu coûteux précisément là où les outils établis s'arrêtent : il rend l'OCS exacte et reproductible praticable à l'échelle génomique, sur un ordinateur portable.

## 1. Introduction

L'amélioration génétique du bétail et des plantes doit tenir deux objectifs en tension : maximiser le mérite génétique de la génération suivante, et conserver la diversité génétique dont l'érosion entraîne la dépression de consanguinité et ferme la porte à la réponse future à la sélection. La sélection à contribution optimale (OCS ; Meuwissen 1997) explicite ce compromis. Étant donné la valeur génétique estimée de chaque candidat et la matrice de parenté entre candidats, l'OCS choisit les contributions génétiques proportionnelles **c** qui maximisent le gain espéré **bᵀc** sous une borne sur la coancestrie moyenne de la descendance, **cᵀGc ≤ k**, les contributions étant non négatives et de somme un — et, dans tout schéma d'accouplement réel, réparties de sorte que pères et mères fournissent chacun la moitié. L'OCS est le cadre de fait pour gérer la diversité, tant en sélection commerciale qu'en programmes de conservation, et la matrice de parenté génomique (VanRaden 2008) a largement remplacé son prédécesseur fondé sur le pedigree comme **G** qu'elle contraint.

Avec des données génomiques, il s'agit d'un programme quadratique convexe à contrainte quadratique, équivalent à un programme conique du second ordre, et deux coûts en viennent à dominer à mesure que le nombre de candidats croît. Le premier est la matrice de parenté elle-même : dense et n×n, elle coûte O(n²) à stocker et O(n²m) à construire à partir de m marqueurs, et à des dizaines de milliers de candidats elle ne tient plus dans la mémoire d'une station de travail. Le second est la résolution conique. Or le vecteur de contributions optimal est presque toujours porté par une poignée d'individus — une parcimonie qui, comme nous le montrons, persiste à mesure que n croît — et les solveurs généralistes n'exploitent ni cela ni la structure peu coûteuse de la contrainte de parenté.

Les outils en usage couvrent l'exact et l'heuristique, mais partagent un même schéma : assembler la matrice de parenté complète et confier le problème entier à un optimiseur générique portant sur tous les candidats. La méthode originale de Meuwissen (1997), et son successeur Gencont2 (Dagnachew & Meuwissen 2016), itèrent les conditions lagrangiennes sur l'ensemble complet des candidats. Le paquet optiSel (Wellmann 2019), l'outil exact de référence, formule l'OCS comme un programme conique pour un solveur primal-dual à point intérieur. Les formulations semi-définies (Pong-Wong & Woolliams 2007) et de récents solveurs ADMM/JuMP (Waldmann 2025) suivent la même voie, ce dernier formant explicitement la matrice dense **G** et ne tronquant les petites contributions qu'a posteriori. AlphaMate (Gorjanc & Hickey 2018) troque l'exactitude contre une heuristique à évolution différentielle qui alloue aussi les accouplements. La seule tentative d'exploiter une parcimonie, Yamashita, Mullin & Safarina (2018), exploite la parcimonie de l'*inverse du pedigree* **A⁻¹** — une propriété de la matrice de données — au sein d'une résolution à point intérieur sur tous les candidats ; ce n'est ni la parcimonie de la *solution*, ni une méthode sans matrice fondée sur les génotypes. Aucun de ces outils n'exploite les deux faits qui rendent l'OCS génomique peu coûteuse : que le support optimal est minuscule et borné en n, et qu'à support fixé le sous-problème contraint en parenté admet une solution en forme fermée. Aucun n'évite de matérialiser **G**, qui devient la contrainte limitante — en mémoire et en temps de construction — précisément aux tailles de population où l'OCS importe.

Nous présentons **support-first** (« le support d'abord »), un solveur OCS exact bâti sur ces deux faits. Il maintient un petit support de travail, initialisé avec le meilleur candidat de chaque sexe et étendu par une règle d'ensemble actif fondée sur le coût réduit : les candidats sont évalués contre les multiplicateurs courants, le plus avantageux est introduit, ceux qui deviennent négatifs sont retirés, le tout vers une unique borne de coancestrie fixée k. Chaque sous-problème à support fixé — maximiser une forme linéaire sur un ellipsoïde intersecté avec les contraintes affines de somme (et de sexe) — est résolu en forme fermée, en éliminant les contraintes d'égalité et en se ramenant à une équation quadratique scalaire dans le multiplicateur actif, sans solveur itératif interne. Les produits de parenté **G·c** sont formés sans matrice à partir de la matrice de génotypes centrée **Z** sous la forme ε**c** + **Z**(**Zᵀc**)/s, de sorte que le solveur ne construit ni ne stocke jamais la matrice dense n×n **G**.

Les ingrédients pris isolément sont classiques, et nous n'en revendiquons aucun. Les méthodes d'ensemble actif qui suivent un petit ensemble de travail d'actifs détenus descendent de l'algorithme de la ligne critique de Markowitz (1956) pour la sélection moyenne–variance long-only, dont l'OCS est l'analogue contraint en parenté ; la forme fermée pour un objectif linéaire sur un ellipsoïde intersecté avec des contraintes linéaires est un résultat de valeur propre contrainte / équation séculaire (Gander, Golub & von Matt 1989) ; et le produit sans matrice **Z**(**Zᵀc**) est standard en prédiction génomique (VanRaden 2008 ; Legarra & Misztal 2008). Notre contribution est leur **synthèse, spécialisée à l'OCS génomique** : une génération de colonnes fondée sur le coût réduit qui fait croître un support minuscule vers une *unique* borne de coancestrie — par opposition au balayage paramétrique de *toute* la frontière efficace de l'algorithme de la ligne critique — avec le produit de parenté évalué sans matrice, de sorte que le coût d'une résolution suit la taille du support et le nombre de marqueurs, plutôt que la matrice dense n×n. À notre connaissance, aucune méthode OCS antérieure n'exploite la parcimonie du support de la solution, et aucune n'est sans matrice fondée sur les génotypes ; le bénéfice, le produit sans matrice étant ce qui rend une résolution à support borné rentable à n génomique, tient à la combinaison plutôt qu'à l'une quelconque de ses parties.

Nous l'évaluons sur trois panels de marqueurs réels : un panel de blé CIMMYT (n = 599), un panel de porc PIC (n = 3534, 52k SNP) et un panel de souris en population hétérogène (n = 1814, avec sexe réel). Il atteint l'optimum exact sur chacun, en accord à 10⁻⁸ près avec l'optimum conique. À coancestrie réalisée appariée, il coïncide avec les méthodes à point intérieur ; là où celles-ci s'arrêtent juste à l'intérieur de la contrainte, il atteint la frontière, de sorte que son léger avantage y est le budget de diversité qu'elles laissent inemployé, et non un optimum différent.

Il est aussi rapide. Le solveur sans matrice livré, qui ne forme aucun **G**, s'exécute 13× à 182× plus vite qu'optiSel de bout en bout ; recevant **G** comme l'exigent les outils à point intérieur, l'ensemble actif seul atteint le même optimum jusqu'à ~2500× plus vite. Face à un solveur conique généraliste à n = 10000, le facteur est de ~76000×, dans une comparaison Rust contre Rust en un seul processus, donc exempte de tout biais de langage. Face à AlphaMate, une heuristique pour le problème distinct de l'allocation discrète des accouplements, l'optimum exact n'est pas moins bon, à coancestrie appariée, sur la relaxation continue commune aux deux méthodes.

Mais la vitesse est le moindre des résultats. Parce qu'une résolution coûte le support et le nombre de marqueurs plutôt que la matrice n×n, l'avantage s'élargit précisément là où les outils établis s'arrêtent. Sur des populations synthétiques, le support optimal reste entre 14 et 19 quand n croît de 1000 à 40000, et support-first résout encore en moins de 0,1 s — tandis que la matrice dense **G** que les alternatives doivent matérialiser atteint 11,9 Gio, une empreinte 40× plus grande que **Z**, au-delà de la mémoire vive d'un portable de 8–16 Go. À cette taille, la question n'est plus la vitesse des outils établis, mais s'ils s'exécutent seulement. Support-first rend la sélection à contribution optimale exacte et reproductible praticable à l'échelle génomique, sur un ordinateur portable.

## 2. Méthodes

### Sélection à contribution optimale

Pour n candidats à la sélection de valeurs génétiques estimées **b** ∈ ℝⁿ et de matrice de parenté génomique **G** ∈ ℝⁿˣⁿ, la sélection à contribution optimale choisit les contributions génétiques proportionnelles **c** qui résolvent

  maximiser **bᵀc**  sous  **Ac = d**,  **0 ≤ c ≤ u**,  **cᵀGc ≤ k**.

Les contraintes affines **Ac = d** encodent le budget de contribution. Dans la forme simplexe, une seule ligne **A = 𝟙ᵀ**, **d = 1**, impose Σcᵢ = 1. Dans la forme *sexuée* — la vraie OCS, puisque chaque accouplement tire un parent de chaque sexe — **A** est la matrice d'incidence de sexe 2×n (ligne 1 l'indicatrice mâle, ligne 2 l'indicatrice femelle) et **d** = (½, ½)ᵀ, imposant Σ_{mâles}cᵢ = Σ_{femelles}cᵢ = ½. Le plafond par candidat **u** borne la contribution d'un individu donné — un sélectionneur laisse rarement un seul parent dominer une génération — et est inactif à **u = 1**. La borne de parenté **cᵀGc ≤ k** plafonne la coancestrie moyenne de la descendance. Nous prenons **G** comme la matrice de parenté génomique de VanRaden, **G = ZZᵀ/s + εI**, où **Z** est la matrice n×m des génotypes centrés par deux fois les fréquences alléliques, s = 2Σⱼpⱼ(1−pⱼ), et ε une faible régularisation ridge assurant la définie-positivité (la même matrice régularisée que contraignent les références coniques, de sorte que tous les solveurs affrontent un problème identique). Le programme est un programme quadratique convexe à contrainte quadratique, équivalent à un programme conique du second ordre.

### Deux faits structurels

Support-first repose sur deux propriétés de l'optimum **c\***. Premièrement, à une borne de parenté active, **c\*** a pour support un petit ensemble S = {i : 0 < c\*ᵢ < uᵢ} de contributions *libres*, et S est empiriquement borné quand n croît (Résultats). Deuxièmement, sur la face où les contraintes de borne inactives sont supprimées (candidats fixés à 0 ou à leur plafond **u**), le problème restreint à S — maximiser une forme linéaire sur un ellipsoïde intersecté avec les contraintes affines — admet une solution en forme fermée. Toute la difficulté est donc d'identifier S et de savoir quels candidats sont à leur plafond ; tout le reste est une petite résolution directe.

### Forme fermée sur un support fixé

Fixons un support S de taille |S|, et notons **G_S**, **b_S**, **A_S** les restrictions à S (**A_S** est q×|S| avec q = 1 ou 2 lignes d'égalité). À l'optimum du problème restreint, la contrainte de parenté est active, et les conditions de stationnarité s'écrivent

  **b_S = A_Sᵀ μ + 2λ G_S c_S**,   **A_S c_S = d**,   **c_Sᵀ G_S c_S = k**,

avec multiplicateurs **μ** ∈ ℝ^q et λ ≥ 0. Résoudre la stationnarité pour **c_S** donne **c_S = (1/2λ) G_S⁻¹(b_S − A_Sᵀμ)**. Imposer les contraintes affines élimine **μ** : avec la matrice q×q **P = A_S G_S⁻¹ A_Sᵀ** et le vecteur **q_v = A_S G_S⁻¹ b_S**,

  **μ = P⁻¹(q_v − 2λ d)**,   d'où   **c_S = g/(2λ) + h**,

où **g = G_S⁻¹(b_S − A_Sᵀ P⁻¹ q_v)** et **h = G_S⁻¹ A_Sᵀ P⁻¹ d** vérifient **A_S g = 0** et **A_S h = d**, de sorte que **A_S c_S = d** pour tout λ. En substituant dans l'ellipsoïde actif **c_Sᵀ G_S c_S = k** et en posant α = **gᵀG_S g**, β = **gᵀG_S h**, γ = **hᵀG_S h**, on obtient une unique équation quadratique scalaire dans le multiplicateur :

  **4(γ − k) λ² + 4β λ + α = 0.**

On prend la racine positive maximisant **b_Sᵀ c_S**. Comme **G_S g = b_S − A_SᵀP⁻¹q_v** et **G_S h = A_SᵀP⁻¹d** sont déjà disponibles, α, β, γ ne coûtent que des produits scalaires — pas de seconde factorisation. Une résolution restreinte est donc : une factorisation de Cholesky de **G_S**, des remontées pour les q+1 seconds membres **[A_Sᵀ | b_S]**, une inversion q×q (q ≤ 2) et une équation quadratique scalaire — sans itération interne. Le cas simplexe (q = 1, **A_S = 𝟙ᵀ**) est identique et se réduit à une équation quadratique scalaire en μ. Quand les égalités déterminent entièrement **c_S** (|S| = q + 1), les contributions sont forcées — (½, ½) dans le cas sexué — et l'on pose λ = 0 (le sous-cas actif est de mesure nulle, comme pour un support singleton). La résolution renvoie « infaisable sur S » quand S manque d'un sexe ou que l'ellipsoïde ne rencontre pas l'enveloppe affine, signal pour agrandir S.

### L'algorithme support-first

Le support est trouvé par génération de colonnes en ensemble actif. Il est initialisé avec le meilleur candidat de chaque sexe, S = {argmax_{mâles} b, argmax_{femelles} b}. Chaque itération résout la forme fermée sur le S courant et se ramifie :

- **Infaisable sur S** (le support est trop apparenté pour atteindre la borne k) : ajouter les candidats les moins apparentés, ceux de plus petit **(Gc)ⱼ**, par paquet doublant |S|. Cette phase de faisabilité, lorsqu'elle se faisait un candidat à la fois, était le coût dominant ; le traitement par paquets a réduit le nombre total d'itérations d'environ un facteur trente, pour un optimum identique.
- **Une contribution est négative** : retirer tout i tel que c_Sᵢ ≤ 0, en les marquant tabous à la réintroduction jusqu'à la prochaine amélioration véritable, puis re-résoudre ; si un retrait vide un sexe, son meilleur candidat est réintroduit.
- **Sinon** : évaluer chaque candidat j ∉ S par son coût réduit **rⱼ = bⱼ − μ_{sexe(j)} − 2λ (Gc)ⱼ**, et ajouter le maximiseur si rⱼ > tol ; si aucun n'est positif, arrêter. Un ajout véritable vide l'ensemble tabou.

Le problème étant convexe, un point réalisable où aucun candidat n'a de coût réduit positif satisfait les conditions de Karush–Kuhn–Tucker et constitue l'optimum global ; la correction n'est donc pas heuristique. L'ensemble tabou, conjugué à un seuil de dégénérescence calibré sur l'ordre de grandeur des coefficients, garantit une terminaison finie. L'évaluation par coût réduit vers une *unique* borne fixée k distingue la mise à jour du support de l'algorithme de la ligne critique, qui modifie l'ensemble actif de ±1 le long d'un multiplicateur balayé pour tracer toute la frontière efficace.

Un plafond par candidat **c ≤ u** est absorbé par le même ensemble actif sans perte de la forme fermée. Une contribution libre qui dépasse son plafond est fixée à uⱼ et déplacée dans un ensemble *supérieur* U ; un candidat déjà dans U est relâché vers l'ensemble libre quand son coût réduit devient négatif. À U fixé, la résolution restreinte conserve sa structure — les contributions fixées entrent dans le second membre affine sous la forme **d − A_U u_U** et dans l'ellipsoïde actif comme un décalage sans matrice **G_{·U} u_U** plus une constante, laissant la même factorisation de Cholesky et la même équation quadratique scalaire. La résolution bornée se réduit à la résolution non bornée quand aucun plafond n'est actif, et atteint le même optimum qu'un solveur conique à point intérieur sur le programme conique **c ≤ u** correspondant, sous plafonds actifs.

### Produits de parenté sans matrice

Tout produit de parenté complet est formé sans matérialiser **G** :

  **Gc = ε c + Z(Zᵀc) / s**,

deux produits matrice–vecteur contre la matrice de génotypes n×m, en O(n·m) en temps et O(n·m) en mémoire pour **Z** contre O(n²) pour une **G** dense. La **G_S** restreinte n'est assemblée que sur le support. C'est l'élément habilitant à grand n, où la matrice dense n×n **G** est infaisable à stocker ou coûteuse (O(n²m)) à construire ; ce n'est pas une accélération de boucle interne quand m > n, et l'avantage d'un ordre de grandeur sur les références coniques est algorithmique — le petit ensemble actif — indépendant de ce choix.

### Implémentation et reproductibilité

Le solveur est implémenté en Rust sur une pile purement Rust : algèbre linéaire dense (la factorisation de Cholesky et les remontées sur le support) via `faer`, sans dépendance système BLAS/LAPACK, sans `unsafe`, et avec un générateur aléatoire reproductible amorcé pour les générateurs synthétiques. Un solveur conique à point intérieur (`clarabel`) n'est inclus que comme oracle indépendant de contre-vérification. La correction est attestée par des tests unitaires à certificat KKT sur une plage de bornes k — un résultat `Solved` est toujours réalisable, sur le budget et réparti par sexe — et par un accord inter-langage : sur une instance active, le `solve_sexed` en Rust reproduit un optimum de référence NumPy à une différence de gain de 1,5×10⁻¹⁴ près. La compilation est conditionnée à `cargo fmt`, `cargo clippy -D warnings` et la suite de tests.

### Données et références de comparaison

La méthode est évaluée sur trois panels de marqueurs publics : un panel de blé CIMMYT (`BGLR`, n = 599), un panel de porc PIC (n = 3534, 52k SNP, avec valeurs génétiques estimées réelles) et un panel de souris en population hétérogène (`BGLR`, n = 1814, avec sexe enregistré, 934 mâles et 880 femelles, l'indice de masse corporelle servant de critère de sélection). Pour chacun, **G** est la matrice de VanRaden avec régularisation ridge ε = 10⁻⁵.

Une quatrième instance est dérivée du panel PIC et est la seule à porter simultanément une vraie valeur génétique et un sexe réel. Les données PIC sont livrées avec un pedigree, et un pedigree assigne le sexe sans ambiguïté : un animal figurant en colonne père est mâle, en colonne mère femelle. Le croisement avec les animaux génotypés restitue un sexe enregistré pour 1194 des 3534 (390 pères, 804 mères, aucun n'apparaissant comme les deux), tous porteurs d'une valeur génétique estimée réelle issue de l'évaluation associée (précision moyenne 0,811 sur le caractère 3). Les fréquences alléliques, et donc **G**, sont recalculées sur les animaux retenus. Le sexe n'étant récupérable que pour les animaux devenus parents, ce sous-panel est un sous-ensemble post-sélection plutôt qu'un échantillon aléatoire de candidats — une propriété de la récupération, non un choix, qui porte sur l'interprétation du vecteur de contributions et non sur ce que mesure la comparaison. Les valeurs fournies avec ce panel provenant de l'évaluation propre au programme de sélection et n'étant pas des prédictions génomiques, le sous-panel est également traité avec un critère calculé comme le calcule la sélection génomique : un GBLUP ajusté par `BGLR` sur les phénotypes du caractère 3 du panel (3141 des 3534 animaux) contre sa propre matrice de parenté génomique, le critère étant lu comme la moyenne a posteriori de l'effet génétique. L'héritabilité génomique ajustée vaut 0,247, et 119 des 1194 animaux du sous-panel ne portent aucun phénotype propre : leur valeur provient donc uniquement des parentés génomiques — le cas même pour lequel la sélection génomique existe. Les références de comparaison sont optiSel (R, le solveur à point intérieur Nesterov–Todd `cccp` ; la référence exacte du domaine), Clarabel (la contre-vérification conique) et AlphaMate (évolution différentielle en Fortran ; exécuté depuis son binaire Linux sous émulation, faute de compilation macOS et la source étant verrouillée à une chaîne d'outils Intel). Réserves honnêtes reportées à la Discussion : le sexe enregistré est réel pour le panel souris et pour le sous-panel porcin sexué, tandis que le blé (autogame) et le panel porcin complet reçoivent une répartition équilibrée arbitraire ; le critère de sélection **b** est une valeur génétique estimée réelle sur les panels porcins et un phénotype enregistré en tenant lieu sur le blé et la souris ; et le Tableau 2 rapporte deux chronométrages de support-first face à optiSel — l'ensemble actif avec **G** dense fournie (le régime d'optiSel, un prototype NumPy) et le solveur sans matrice en Rust livré, de bout en bout — séparant ainsi l'avantage algorithmique de ce que l'outil livré paie réellement, au lieu de les confondre.

## 3. Résultats

### Exactitude

Support-first est exact par construction : pour un programme convexe, un point réalisable où aucun candidat n'a de coût réduit positif satisfait les conditions KKT, et c'est ce que l'ensemble actif certifie à la terminaison (et ce que les tests unitaires vérifient sur une plage de bornes de parenté). Empiriquement, sur l'exemple `Cattle` d'optiSel (bovins Angler, n = 268) avec le sexe réel enregistré fourni par le paquet, support-first et optiSel sélectionnent le *même* support de 36 individus et s'accordent sur les contributions à 3×10⁻⁴ près (contribution maximale 0,0664 contre 0,0661). Ils ne diffèrent qu'à la frontière de parenté : l'optimum OCS se situe sur la borne (le gain y est monotone), et support-first l'atteint (cᵀGc = k), tandis que le solveur à point intérieur d'optiSel s'arrête à sa tolérance de convergence juste à l'intérieur de la région réalisable (parenté de groupe 0,0576 contre une borne de 0,0578). À coancestrie réalisée appariée, les deux coïncident — sur le panel souris, tous deux atteignent un gain de −0,29754 à une coancestrie de 0,0346 — de sorte que le léger avantage de support-first à son propre point de fonctionnement est le budget de diversité qu'optiSel laisse inemployé sur la frontière, et non un optimum différent : la différence est la saturation de la frontière versus un arrêt intérieur. Face à un solveur conique à point intérieur indépendant (Clarabel) sur données synthétiques, les deux s'accordent à une différence de gain inférieure à 10⁻⁸ sur toute la plage de bornes de parenté, et l'implémentation Rust du solveur sexué reproduit un optimum de référence NumPy à 1,5×10⁻¹⁴ près sur une instance active. Support-first est donc au moins aussi exact que l'outil du domaine, et exact là où les méthodes à point intérieur ne font qu'approcher.

### Vitesse face à un solveur conique généraliste

Cette comparaison tranche en outre une question que la confrontation à optiSel ne peut pas trancher : les deux solveurs sont ici en Rust, appelés dans un même processus sur les mêmes données générées, de sorte qu'aucune différence de langage ne peut être prise pour une différence algorithmique. Sur des instances génomiques synthétiques, support-first est plus rapide que le solveur conique Clarabel d'un facteur qui croît fortement avec n (Tableau 1) : de ~290× à n = 1000 à ~76000× à n = 10000, où la voie conique prend une demi-heure et support-first 24 ms. Ce passage à l'échelle est structurel — Clarabel factorise un système KKT dense en O(n³) à chaque itération de point intérieur, tandis que support-first effectue un nombre quasi constant de produits matrice–vecteur (3 à 6 ici), chacun en O(n·m).

**Tableau 1.** support-first contre Clarabel, optimum identique (accord en gain de 1,2×10⁻⁹ à 4,8×10⁻⁹), un balayage en série sur machine au repos. Les deux solveurs sont en Rust dans un même processus sur la même instance synthétique (m = 20000, borne k = 0,6·moyenne diag **G**) : la comparaison ne comporte aucun biais de langage.

| n | m | Clarabel¹ | support-first | accélération |
|---|---|---|---|---|
| 1000 | 20000 | 1,70 s | 0,006 s | 290× |
| 2000 | 20000 | 12,66 s | 0,012 s | 1093× |
| 5000 | 20000 | 182,9 s | 0,018 s | 9917× |
| 10000 | 20000 | 1792 s | 0,024 s | **76179×** |

¹ toute la voie conique : construction de **G** (7,0 s à n = 10000), sa factorisation de Cholesky, l'assemblage du cône et la résolution à point intérieur. La résolution seule fait 1,60 / 12,31 / 180,9 / 1784,2 s, soit 274× à 75831×. La colonne support-first est son exécution entière, puisqu'il ne construit rien.

Deux réserves honnêtes. Les temps de support-first sont déterministes à cette précision (trois exécutions à n = 1000 et à n = 2000 reproduisent 0,006 s et 0,012 s à la milliseconde) tandis que Clarabel varie de quelques pour cent ; les cellules à grand n sont des exécutions uniques, Clarabel demandant une demi-heure à n = 10000. Et la borne est ici lâche, de sorte que le support optimal vaut 2 à 5 — le bout facile du spectre, ce qui explique le très petit nombre de produits de parenté nécessaires. Le régime à borne serrée, où le support monte à 59–168, est la moitié synthétique du Tableau 2, et le facteur y diminue en conséquence ; nous rapportons les deux régimes plutôt que le seul favorable.

### Vitesse face à l'outil du domaine optiSel

Un solveur générique est une cible facile ; la comparaison informative se fait contre optiSel, l'outil OCS exact de référence, sur sa propre formulation. Après extension de support-first aux deux contraintes d'égalité de sexe, il renvoie le même optimum qu'optiSel. Le solveur sans matrice livré — qui ne forme aucun **G** — s'exécute 13× à 182× plus vite sur des panels réels et synthétiques (Tableau 2), le facteur s'élargissant avec n à mesure que le coût du point intérieur d'optiSel augmente. Avec **G** précalculée, comme l'exigent les outils à point intérieur, l'ensemble actif seul atteint le même optimum ~90× à 2500× plus vite : l'avantage algorithmique que cette section dissèque, dont le solveur sans matrice échange une partie pour ne jamais construire la matrice.

**Tableau 2.** support-first contre optiSel, même optimum partout, mesuré en série sur une machine au repos (Apple M4 Max) : le solveur sans matrice et optiSel dans une même session R, le prototype algorithmique en Python sur les mêmes instances. Chaque temps est une médiane (sans matrice et algorithme 5×, optiSel 3×). optiSel contraint cᵀ·sKin·c ≤ ub sur sKin = **G**/2 + 10⁻⁵**I**, de sorte que support-first reçoit le problème identique sous la forme k = 2·ub sur **G** + 2×10⁻⁵**I**. Deux colonnes support-first : l'ensemble actif avec **G** dense fournie (le régime d'optiSel), et le solveur sans matrice livré, qui ne forme aucun **G**, chronométré de bout en bout depuis les génotypes, copie du *binding* R comprise.

| jeu de données | n | m | algorithme¹ | sans matrice² | optiSel | accél.³ |
|---|---|---|---|---|---|---|
| synthétique (structuré) | 1000 | 500 | 0,023 s | 0,117 s | 2,06 s | 18× |
| synthétique (structuré) | 2000 | 500 | 0,102 s | 0,280 s | 13,16 s | 47× |
| synthétique (structuré) | 5000 | 500 | 0,506 s | 0,902 s | 163,8 s | 182× |
| blé CIMMYT (GRM réelle) | 599 | 1279 | 0,003 s | 0,008 s | 0,595 s | 74× |
| porc PIC (GRM réelle, 52k SNP) | 3534 | 52843 | 0,020 s | 0,622 s | 49,9 s | 80× |
| souris HS (GRM réelle, sexe réel) | 1814 | 10346 | 0,008 s | 0,057 s | 6,93 s | **122×** |
| porc PIC, sous-panel sexué (VGE réelle **et** sexe réel) | 1194 | 52843 | 0,004 s | 0,227 s | 2,85 s | 13× |

¹ ensemble actif avec **G** dense précalculée, solve seul — le régime d'optiSel (prototype NumPy) · ² sans matrice, ne forme aucun **G**, bout en bout, copie des génotypes du *binding* R comprise (Rust livré) · ³ optiSel / sans matrice

La colonne accélération est celle du solveur livré, qui ne forme rien, face à optiSel : 13× à 182×. La colonne algorithme — l'ensemble actif recevant **G** comme les outils à point intérieur — atteint le même optimum ~90× à 2500× plus vite (le ~2500× du porc le plus grand, et l'origine du chiffre antérieur de ce travail), mais paie le O(n²m) de construction et le O(n²) de stockage de **G** que le solveur sans matrice évite ; recevant **G**, le solve évite de streamer **Z**, il est donc le plus rapide des deux sur chaque panel, et la colonne sans matrice est le prix pour ne rien garder de dense en mémoire. Sur le porc, une grande part des 0,622 s livrés est la copie par le *binding* R de la matrice de génotypes de 1,5 Go, de sorte qu'à m = 52 843 c'est le mouvement des données, non le solve, qui fixe le temps bout en bout. Les deux colonnes renvoient le même optimum sur chaque ligne, atteignant la frontière où optiSel s'arrête juste à l'intérieur.

La dernière ligne est la seule instance portant à la fois une vraie valeur génétique et un sexe réel, donc la plus proche de ce qu'exécute réellement un programme de sélection. Le sexe n'y est pas supposé mais lu dans le pedigree fourni avec les données PIC — père ⇒ mâle, mère ⇒ femelle —, ce qui résout 1194 des 3534 animaux génotypés (390 pères, 804 mères, aucun n'apparaissant comme les deux), tous porteurs d'une VGE réelle issue de l'évaluation associée, de précision moyenne 0,811. Support-first y atteint un gain de 1,844614 contre 1,842975 pour optiSel, sur la frontière de contrainte où optiSel s'arrête à 99,53 % de celle-ci, le budget par sexe étant réparti exactement par moitié et 6 mâles et 21 femelles étant retenus sur 1194 candidats. C'est aussi la ligne la moins favorable à la voie sans matrice : à m/n ≈ 44 la matrice de marqueurs domine, d'où une accélération (13×) qui est la plus faible du tableau alors que la même résolution avec **G** fournie prend 4 ms — l'arbitrage énoncé dans les Méthodes, visible là où il mord le plus. La ligne souris (934 mâles, 880 femelles) reste l'autre instance véritablement sexuée. Réserves reportées à la Discussion : sur le blé et la souris, le critère de sélection est un phénotype enregistré plutôt qu'une valeur génétique, et dans aucun panel il ne s'agit d'une valeur génétique *génomique* ; et le sous-panel porcin sexué, constitué des animaux devenus parents, est un sous-ensemble post-sélection.

Une dernière vérification porte sur le critère lui-même plutôt que sur le solveur. Reprendre le sous-panel sexué avec une vraie GEBV — un GBLUP ajusté sur les phénotypes et la matrice de parenté du panel, corrélé à 0,64 avec la VGE de l'évaluation — ne change rien de matériel : support-first atteint un gain de 0,836775 contre 0,836057 pour optiSel, de nouveau sur la frontière de contrainte où optiSel s'arrête juste à l'intérieur (100 % contre 99,54 % de la borne), avec un support de 28, le budget par sexe réparti exactement, en 0,233 s contre 2,850 s (12×). Le comportement rapporté ici n'est donc pas un artefact du choix de la valeur génétique servant de critère.

### Comparaison avec l'heuristique AlphaMate

AlphaMate optimise un problème apparenté mais distinct — l'allocation discrète des accouplements par un algorithme évolutionnaire stochastique — de sorte que ce n'est pas une confrontation à armes égales, et nous l'interprétons en conséquence. Nous évaluons *son* vecteur de contributions dans notre métrique et comparons, à coancestrie de groupe appariée, uniquement sur la relaxation à contributions continues commune aux deux méthodes. Sur cette relaxation, l'optimum exact de support-first a un gain supérieur à celui du vecteur d'AlphaMate en chaque point de sa frontière (Tableau 3) : une faible marge au compromis à 45° (Δgain +0,004) et des marges plus grandes aux coins — comme il se doit d'un optimum convexe exact face à une heuristique stochastique portant en outre des contraintes d'accouplement entières. Nous lisons cela comme un contrôle de cohérence (support-first est exact sur la relaxation commune), et non comme une victoire sur AlphaMate à sa propre tâche, l'allocation discrète des accouplements, que support-first n'effectue pas et qui est hors champ ici. AlphaMate s'est en outre révélé nettement fragile sur données génomiques réelles : une exécution réussie a exigé six configurations et trois contournements distincts — plafonner les accouplements sous n ; restaurer l'ensemble complet des parents, pour éviter une faute de segmentation à l'initialisation déclenchée par le nombre réduit de parents ; et décaler positivement le critère de sélection, pour défaire une inversion de signe valeur-sur-maximum qui faisait maximiser l'heuristique dans le mauvais sens sur les VGE centrées et négatives — alors que support-first et optiSel ont tourné sans modification. AlphaMate a calculé sa frontière entière en 882 s de temps CPU (un binaire x86 émulé, faute de compilation native) ; support-first trace la frontière exacte à ≤ 1,1 s par point.

**Tableau 3.** Gain génétique à coancestrie de groupe appariée, panel souris (évalué dans la même métrique à partir du vecteur de contributions de chaque méthode).

| coancestrie de groupe cᵀKc | AlphaMate | support-first | Δgain |
|---|---|---|---|
| 0,000272 (coancestrie min. AlphaMate) | −0,45885 | **−0,37843** | +0,080 |
| 0,001317 (optimum 45° AlphaMate) | −0,35748 | **−0,35318** | +0,004 |
| 0,007574 (critère max. AlphaMate) | −0,34066 | **−0,32240** | +0,018 |

### Passage à l'échelle et avantage sans matrice (Figure 1)

L'avantage en coût se maintient — et croît — à l'échelle. En balayant le nombre de candidats de 1000 à 40000 à panel de marqueurs fixé (m = 1000) sous une borne de parenté active, le support optimal reste entre 14 et 19 (Figure 1A) et la résolution sans matrice reste sous 0,1 s. La matrice de parenté dense que tout autre solveur doit matérialiser raconte l'histoire inverse. Sa seule *construction* coûte O(n²m) : 3,6 s à n = 30000, déjà 63× la résolution support-first entière à cette taille, et croissant quadratiquement (Figure 1A). Son *stockage* coûte O(n²) : 11,9 Gio à n = 40000, où l'empreinte sans matrice de Z est de 0,30 Gio et où la matrice dense ne tient plus dans la mémoire vive d'un portable (Figure 1B). La matrice dense devient la contrainte limitante, en temps d'initialisation et en mémoire, précisément dans le régime où support-first reste peu coûteux — parce que le coût du solveur suit le support et le nombre de marqueurs, jamais la matrice n×n.

### Comportement du support

L'avantage repose sur le support, dont la taille est une propriété du point de fonctionnement plutôt qu'une constante universelle — aussi le caractérisons-nous le long de toute la frontière, et non à une seule borne. Sur le panel souris, le support croît régulièrement et de façon monotone à mesure que la borne de coancestrie se resserre : 19 candidats sur 1814 à la coancestrie de travail (0,0346), 61 à 0,010, 133 à 0,003, 189 à 0,0015, 473 à 0,0002, et environ 1163 à l'approche de zéro — où minimiser l'apparentement force la solution à s'étaler sur une grande partie de la population — et il s'effondre à deux pour une borne lâche. L'affirmation est donc double, et nous gardons les axes séparés : à une borne de fonctionnement *fixée*, le support est borné quand n croît (14–19 jusqu'à n = 40000, Figure 1A), tandis que *le long* de la frontière, il croît à mesure qu'on pousse la diversité. Le coût par résolution suit le support, qui est petit dans le régime que les sélectionneurs utilisent réellement, et que nous rapportons sur toute la frontière plutôt qu'en un seul point. La question de savoir si cela est démontrable — indépendant de n à borne fixée, croissant quand la borne se resserre — est reprise dans la Discussion.

## 4. Discussion

Support-first rend la sélection à contribution optimale exacte peu coûteuse à l'échelle génomique. Sur les panels testés, il atteint le même optimum que l'outil exact du domaine tout en s'exécutant un à deux ordres de grandeur plus vite de bout en bout, sans former aucune matrice de parenté — et, recevant cette matrice comme l'exigent les outils à point intérieur, son ensemble actif seul atteint le même optimum jusqu'à trois ordres de grandeur plus vite. Il reste peu coûteux précisément là où la matrice de parenté dense que tout autre solveur forme devient infaisable — le régime à grand nombre de candidats qui motive l'OCS génomique en premier lieu. La méthode étant exacte et déterministe — un ensemble actif certifié Karush–Kuhn–Tucker plutôt qu'une recherche stochastique — sa sortie est reproductible jusqu'au dernier chiffre, contrairement aux outils heuristiques de sélection d'accouplements auxquels elle est comparée.

Les ingrédients de la méthode sont individuellement classiques, et nous n'en revendiquons que la synthèse. Le suivi par ensemble actif d'un petit ensemble de travail descend de l'algorithme de la ligne critique pour la sélection moyenne–variance long-only ; la forme fermée par support est un résultat de valeur propre contrainte ; le produit sans matrice est standard en prédiction génomique. La contribution est leur combinaison, spécialisée à l'OCS : une génération de colonnes par coût réduit qui fait croître un support minuscule vers une unique borne de coancestrie, avec le produit de parenté évalué sans matrice de sorte que le coût suit le support et le nombre de marqueurs plutôt que la matrice n×n. Située face aux travaux antérieurs les plus proches, ce n'est ni la parcimonie de la matrice de pedigree d'un solveur OCS à point intérieur, ni les résolutions denses sur toute la population des outils de référence.

Plusieurs limites bornent ces résultats, et nous les énonçons clairement. Premièrement, le critère de sélection est une valeur génétique estimée réelle sur les panels porcins et un phénotype enregistré en tenant lieu sur le blé et la souris. Sur le sous-panel porcin sexué, nous avons de plus utilisé une vraie GEBV, ajustée par GBLUP sur les phénotypes et la matrice de parenté de ce panel, sans que le résultat change : les conclusions ne dépendent donc pas de la valeur génétique retenue comme critère. Sur le blé et la souris, en revanche, le critère demeure un phénotype, et les vecteurs de contributions restent dans tous les cas illustratifs plutôt que des recommandations de sélection. Deuxièmement, un sexe réellement enregistré est disponible pour le panel souris et pour le sous-panel porcin sexué — la seule instance combinant une vraie valeur génétique et un sexe réel, et donc ce qui se rapproche le plus ici d'un problème OCS opérationnel —, tandis que le blé et le panel porcin complet reçoivent une répartition équilibrée arbitraire. Ce sous-panel porte sa propre restriction : le sexe se déduit du pedigree précisément pour les animaux devenus parents, c'est donc un sous-ensemble post-sélection des candidats et non un échantillon aléatoire de ceux-ci. Troisièmement, les chronométrages en confrontation directe comparent un prototype NumPy à R/optiSel : l'écart d'un ordre de grandeur est algorithmique — tous deux réalisent le même ensemble actif — mais une comparaison mono-langage le placerait hors de doute, et nos propres mesures sont explicites sur le fait que le produit sans matrice n'est *pas* une accélération de boucle interne quand les marqueurs dépassent les candidats (m > n), où le parcours en flux de la matrice de génotypes coûte plus qu'un produit dense résident. La voie sans matrice est l'élément habilitant en mémoire et à grand n ; l'avantage en vitesse sur les solveurs coniques est le petit ensemble actif. Quatrièmement, le caractère borné du support optimal est, dans cet article, une observation empirique sur un balayage synthétique et les panels réels, et pas encore un théorème, et la voie est plus subtile qu'un comptage de contraintes. Sans ridge (ε = 0, G = ZZᵀ/s de rang r), un argument propre le borne : l'optimum maximise un lagrangien qui ne dépend de c qu'à travers (Zᵀc, bᵀc), de sorte que la tranche optimale {Ac = d, c ≥ 0, Zᵀc = Zᵀc\*, bᵀc = bᵀc\*} est un polytope LP dont les sommets portent ≤ q + r + 1 composantes non nulles — une borne de point extrême / Carathéodory, indépendante de n, sur exactement le sommet que renvoie un solveur à ensemble actif comme support-first (les solveurs à point intérieur et ADMM renvoient au contraire des points intérieurs non parcimonieux qu'ils seuillent a posteriori). Le ridge opératoire ε > 0 rend toutefois G de plein rang, et nous constatons numériquement que le support réalisé n'est alors pas gouverné par rang(G₀) : il reste petit (quelques dizaines au plus) et stable en n sur les panels ici, mais ne se réduit pas à une unique quantité spectrale propre — le rang effectif est suggestif mais, sur un balayage de spectres, n'est pas un prédicteur fiable — de sorte qu'une borne serrée pour le problème ridgé reste ouverte. La génétique rend compte de la croissance : ΔF ∝ Σcᵢ² (Wray & Thompson 1990 ; Woolliams & Bijma 2000) avec Ne ≈ 1/(2ΔF), de sorte que resserrer la borne étale les contributions et fait croître le support, comme observé. Nous développons cette caractérisation dans des travaux ultérieurs. Enfin, le solveur gère une unique contrainte quadratique de parenté, des plafonds de contribution par candidat (0 ≤ c ≤ u) et les égalités de sexe, mais pas encore plusieurs contraintes quadratiques à la fois — matrices de parenté multiples, limites de coancestrie spécifiques à des groupes — ni l'allocation entière des accouplements que fournit AlphaMate ; il renvoie des contributions continues.

Chaque limite indique une extension. Le cœur en forme fermée par support admet déjà plus d'une contrainte d'égalité via la même élimination P = A G⁻¹ Aᵀ ; des plafonds quadratiques actifs supplémentaires transforment la recherche de racine scalaire en un petit système de faible degré plutôt qu'ils n'en changent la structure, de sorte que des contraintes de parenté multiples sont à portée. Une couche d'arrondi ou de séparation et évaluation (branch-and-bound) au-dessus de l'optimum continu ajouterait l'allocation des accouplements tout en gardant la relaxation exacte comme borne serrée — combinant l'exactitude de support-first avec le plan discret que les heuristiques visent directement. Le basculement sans matrice en m versus n mérite une caractérisation explicite, tout comme la borne de support elle-même. Et la méthode devrait être validée sur des jeux de données portant de vraies valeurs génomiques et sur des populations un ordre de grandeur plus grandes, où l'écart mémoire de la Figure 1B passe de visible à décisif.

Support-first fait une chose — la sélection à contribution optimale exacte, à contrainte unique, continue — à une échelle et une vitesse qui mettent l'OCS génomique à portée d'un ordinateur portable et d'un script reproductible. Cette affirmation étroite et vérifiée, et les extensions ouvertes qu'elle invite, constituent la contribution.

## Disponibilité des données et du code

Le solveur support-first ainsi que tous les scripts de comparaison et de reproduction sont disponibles à https://github.com/Adelagric/ocs-rs (licence MIT) et archivés sur Zenodo (https://doi.org/10.5281/zenodo.20746987). Chaque tableau et figure est régénéré par `research/repro/repro.sh`, dont les étapes portent le nom du tableau ou de la figure qu'elles produisent ; une étape dont la chaîne d'outils ou le jeu de données manque est ignorée avec un message. Deux exceptions sont délibérées : le point n = 10000 de Clarabel (Tableau 1) ne s'exécute que sous `REPRO_FULL=1` (Clarabel y demande ~26 minutes), et la comparaison AlphaMate se lance manuellement depuis un binaire émulé. Les panels de marqueurs sont publics : les panels blé et souris sont fournis avec le paquet R BGLR (Pérez & de los Campos 2014) ; le panel de porc PIC est le jeu de données commun de Cleveland, Hickey & Forni (2012) (DOI 10.1534/g3.111.001453).

## Références

### Méthodes et outils OCS

1. **Meuwissen, T.H.E. (1997).** Maximizing the response of selection with a
   predefined rate of inbreeding. *Journal of Animal Science* 75(4):934–940.
   DOI 10.2527/1997.754934x.

2. **Dagnachew, B.S. & Meuwissen, T.H.E. (2016).** A fast Newton–Raphson based
   iterative algorithm for large scale optimal contribution selection. *Genetics
   Selection Evolution* 48(1):70. DOI 10.1186/s12711-016-0249-2.

3. **Pong-Wong, R. & Woolliams, J.A. (2007).** Optimisation of contribution of
   candidate parents to maximise genetic gain and restricting inbreeding using
   semidefinite programming. *Genetics Selection Evolution* 39(1):3–25.
   DOI 10.1186/1297-9686-39-1-3.

4. **Wellmann, R. (2019).** Optimum contribution selection for animal breeding and
   conservation: the R package optiSel. *BMC Bioinformatics* 20(1):25.
   DOI 10.1186/s12859-018-2450-5.

5. **Gorjanc, G. & Hickey, J.M. (2018).** AlphaMate: a program for optimizing
   selection, maintenance of diversity and mate allocation in breeding programs.
   *Bioinformatics* 34(19):3408–3411. DOI 10.1093/bioinformatics/bty375.

6. **Waldmann, P. (2025).** Genomic optimum contribution selection and mate
   allocation using JuMP. *Bioinformatics Advances* 5(1):vbaf259.
   DOI 10.1093/bioadv/vbaf259.

7. **Yamashita, M., Mullin, T.J. & Safarina, S. (2018).** An efficient
   second-order cone programming approach for optimal selection in tree breeding.
   *Optimization Letters* 12(7):1683–1697. DOI 10.1007/s11590-018-1229-y ;
   arXiv:1506.04487. — Exploite la parcimonie de l'*inverse du pedigree*, pas le support de la solution.

### Parenté génomique et produits sans matrice

8. **VanRaden, P.M. (2008).** Efficient methods to compute genomic predictions.
   *Journal of Dairy Science* 91(11):4414–4423. DOI 10.3168/jds.2007-0980.

9. **Legarra, A. & Misztal, I. (2008).** Technical note: Computing strategies in
   genome-wide selection. *Journal of Dairy Science* 91(1):360–366.
   DOI 10.3168/jds.2007-0403.

### Filiation ensemble actif et forme fermée (antériorité reconnue)

10. **Markowitz, H. (1956).** The optimization of a quadratic function subject to
    linear constraints. *Naval Research Logistics Quarterly* 3(1–2):111–133.
    DOI 10.1002/nav.3800030110. — L'algorithme de la ligne critique.

11. **Gander, W., Golub, G.H. & von Matt, U. (1989).** A constrained eigenvalue
    problem. *Linear Algebra and its Applications* 114–115:815–839.
    DOI 10.1016/0024-3795(89)90494-1. — La forme fermée par support (équation séculaire).

### Logiciels et données

12. **Goulart, P.J. & Chen, Y. (2024).** Clarabel: an interior-point solver for
    conic programs with quadratic objectives. arXiv:2405.12762. (Version évaluée par
    les pairs : *Mathematical Programming Computation*, 2026,
    DOI 10.1007/s12532-026-00320-7.) — Utilisé uniquement comme oracle de contre-vérification indépendant.

13. **Quiñones, S.** faer: a linear-algebra library for the Rust programming
    language. Logiciel, version 0.24.0 (épinglée dans `Cargo.toml`) ;
    <https://github.com/sarah-quinones/faer-rs>. (Soumission JOSS en cours d'examen ; pas
    de DOI publié à ce jour.)

14. **Pérez, P. & de los Campos, G. (2014).** Genome-wide regression and
    prediction with the BGLR statistical package. *Genetics* 198(2):483–495.
    DOI 10.1534/genetics.114.164442. — Source des panels blé (CIMMYT) et souris.

15. **Cleveland, M.A., Hickey, J.M. & Forni, S. (2012).** A common dataset for
    genomic analysis of livestock populations. *G3: Genes|Genomes|Genetics*
    2(4):429–435. DOI 10.1534/g3.111.001453. — Le panel de porc PIC.

## Annexe : la taille du support — un encadrement

Le coût de support-first suit la taille du support |S| ; on l'encadre (preuves complètes, identité KKT de liaison et carte empirique dans la note compagne [`support_bound_fr.md`](support_bound_fr.md), scripts reproductibles `research/bound_*.py`).

**Théorème 1 (borne sans ridge).** *Pour ε = 0 (G = ZZᵀ/s de rang r ≤ m), l'OCS atteint son optimum en un vecteur de contributions de support |S| ≤ q + r + 1, indépendant de n* (q ∈ {1, 2} lignes de budget). À ε = 0 le lagrangien ne dépend de c qu'à travers (Zᵀc, bᵀc) ∈ ℝ^{r+1} ; fixer ces quantités linéarise la face optimale en un polytope LP dont les sommets ont ≤ q + r + 1 composantes non nulles — le sommet que renvoie un solveur à ensemble actif. (Carathéodory / Barvinok–Pataki sur la face linéarisée, pas sur la frontière courbe de l'ellipsoïde, où tout point est extrême.)

**Théorème 2 (aucune borne universelle avec le ridge).** *Pour G = ZZᵀ/s + εI (ε > 0), aucune borne f(q, r) indépendante de n ne tient — le support peut valoir n.* Témoin : G = εI (aucun marqueur, tous non apparentés) réduit le cap de parenté à ε‖c‖² ≤ k, où ‖c‖² est le proxy du taux de consanguinité ; juste au-dessus du minimum uniforme la boule réalisable force le support plein (vérifié à n = 300). C'est la limite sans structure.

Entre les deux, le support réalisé est fixé par la géométrie conjointe de (spectre, b, cap k) : petit et stable en n précisément quand le spectre décroît et que b évite les directions propres dominantes — le régime des panels réels (|S| ≈ 15–30 ici) — mais sans loi scalaire prédictive unique. Le borner sous cette structure réaliste est le problème ouvert.

## Figure

**Figure 1.** Passage à l'échelle sans matrice contre G dense (`research/fig_scaling.pdf`) ; décrit dans les Résultats, section *Passage à l'échelle et avantage sans matrice*.
