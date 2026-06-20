# Sur la taille du support de l'optimum OCS génomique

> Note autoportante — une preuve, un contre-exemple et le régime empirique — pour le nombre
> d'individus qu'active un optimum de sélection à contribution optimale. Note compagne du
> solveur support-first (`MANUSCRIPT_fr.md`). Version française ; note de référence en anglais
> (`support_bound_sketch.md`). Toute affirmation empirique est reproductible (voir *Évidence
> numérique*).

## 1. La question

La sélection à contribution optimale (OCS) génomique résout, pour n candidats de valeurs
génétiques **b** ∈ ℝⁿ et de matrice de parenté **G = ZZᵀ/s + εI** (**Z** la matrice n×m des
génotypes centrés, s une échelle, ε un petit ridge pour la définie-positivité),

  maximiser **bᵀc**  s.c.  **Ac = d**,  **0 ≤ c ≤ u**,  **cᵀGc ≤ k**,    (OCS)

où **A** est la matrice de budget q×n (q = 1 pour le simplexe 𝟙ᵀc = 1 ; q = 2 pour la forme
sexuée Σ_mâles = Σ_femelles = ½). Empiriquement l'optimum **c\*** active très peu de candidats :
son support S = {i : c\*ᵢ > 0} est ≈ 15–30 sur les panels réels et reste borné quand n croît à
cap fixé. Qu'est-ce qui borne |S| ?

La réponse est un **encadrement**. Dans la limite sans ridge, une borne propre indépendante de n
existe (Théorème 1). Dès que le ridge est présent, il n'y a prouvablement aucune borne
universelle — le support peut valoir toute la population (Théorème 2). Le régime pratique entre
les deux est fixé par la géométrie conjointe de (spectre, **b**, cap k) ; il n'admet aucune loi
scalaire unique, et on le cartographie empiriquement (§5). Une identité KKT (§4) relie les deux
théorèmes.

## 2. Théorème 1 — la borne ε = 0

**Théorème 1.** *Pour ε = 0 (donc* **G** = **G₀** = **ZZᵀ**/s *de rang r ≤ m), (OCS) admet un
optimum* **c\*** *avec* **|S| ≤ q + r + 1**, *indépendant de n.*

*Démonstration.* Par dualité forte (Slater tient pour toute instance intérieurement réalisable),
il existe un multiplicateur λ\* ≥ 0 tel que **c\*** maximise le lagrangien
L(**c**) = **bᵀc** − λ\* **cᵀG₀c** sur le polytope {**Ac = d**, **c ≥ 0**}. Comme
**cᵀG₀c** = ‖**Zᵀc**‖²/s, la valeur L(**c**) ne dépend de **c** qu'à travers le couple
(**Zᵀc**, **bᵀc**) ∈ ℝ^{r+1}. Donc la tranche

  P = { **c ≥ 0** : **Ac = d**, **Zᵀc = Zᵀc\***, **bᵀc = bᵀc\*** }

est **entièrement optimale** : tout **c** ∈ P est réalisable (**cᵀG₀c** = ‖**Zᵀc\***‖²/s ≤ k) et
atteint la valeur optimale **bᵀc\***. P est non vide (il contient **c\***), borné (le budget avec
**c ≥ 0**), et découpé par q + r + 1 lignes d'égalité ; il possède donc un sommet, et un sommet de
{**c ≥ 0** : **Mc = e**} à q + r + 1 lignes a au plus q + r + 1 composantes non nulles. Ce sommet
est un optimum de support annoncé. ∎

C'est l'idée de Carathéodory / Barvinok–Pataki appliquée correctement : on ne compte *pas* sur la
frontière courbe de l'ellipsoïde (où, **G** étant définie positive, tout point est déjà extrême et
le support n'est pas borné) — on *fige* l'image de bas rang **Zᵀc** et l'objectif, ce qui rend la
tranche optimale affine, et on y compte les sommets LP. Les solveurs à ensemble actif
(support-first ; ligne critique) renvoient un tel sommet ; les solveurs à point intérieur et ADMM
renvoient des points intérieurs non parcimonieux seuillés a posteriori — la borne porte donc
précisément sur ce que calcule un solveur à ensemble actif.

## 3. Théorème 2 — aucune borne universelle pour ε > 0

**Théorème 2.** *Pour ε > 0, aucune borne sur |S| de la forme f(q, r) indépendante de n n'existe :
le support peut valoir n.*

*Démonstration (contre-exemple).* Prenons **G** = εI — le cas dégénéré m = 0 (aucun marqueur,
toute paire également non apparentée ; rang(**G₀**) = 0). Alors (OCS), forme simplexe, devient

  maximiser **bᵀc**  s.c.  𝟙ᵀ**c** = 1, **c ≥ 0**, ε‖**c**‖² ≤ k,

et ‖**c**‖² = Σ c²ᵢ est exactement le proxy du taux de consanguinité. Le minimum de ‖**c**‖² sur
le simplexe est 1/n, atteint uniquement au plan uniforme **c** = 𝟙/n. Pour k légèrement supérieur
à ε/n, l'ensemble réalisable est le simplexe intersecté avec une boule qui se contracte sur 𝟙/n,
forçant |S| = n. Numériquement (`bound_validation.py`, bloc 7, n = 300) : |S| = 300, 295, 256,
151, 52, 12 quand k se desserre de 1,05× à 50× le minimum. Donc |S| parcourt [1, n] sans plafond
indépendant de n. ∎

Génétiquement, **G** = εI est la limite *sans structure* ; sans schéma d'apparentement, le cap de
diversité ne peut être satisfait qu'en étalant les contributions sur toute la population. Le petit
support observé sur données réelles est donc une propriété de la **structure** de **G** (et de
**b**), pas une garantie de pire cas.

## 4. Le pont — une identité KKT

**Proposition.** *À un optimum où la contrainte de parenté est active, les contributions sur le
support sont une fonction affine du feature augmenté* (bᵢ, **zᵢ**) ∈ ℝ^{m+1} *de chaque candidat*
(**zᵢ** *la i-ème ligne de* **Z**) :  **c\*ᵢ = α bᵢ + wᵀzᵢ + β_{sexe(i)}**  *pour i ∈ S.*

*Démonstration.* Les conditions KKT donnent μ ∈ ℝ^q, λ ≥ 0, **s ≥ 0** avec
**b** = **Aᵀμ** + 2λ**Gc\*** − **s** et sᵢc\*ᵢ = 0. Sur S (sᵢ = 0) : bᵢ = (**Aᵀμ**)ᵢ + 2λ(**Gc\***)ᵢ.
En posant **y** := **Zᵀc\*** et **Gc\*** = **Z**y/s + ε**c\***, on obtient, pour i ∈ S,
ε c\*ᵢ = (bᵢ − (**Aᵀμ**)ᵢ)/(2λ) − (**zᵢ**ᵀ**y**)/s, soit **c\*ᵢ = α bᵢ + wᵀzᵢ + β_{sexe(i)}** avec
α = 1/(2λε), **w** = −**y**/(sε), β_{sexe(i)} = −(**Aᵀμ**)ᵢ/(2λε). ∎

Les contributions du support vivent donc sur une famille à (m + 2) paramètres. Quand ε → 0 le terme
εc\*ᵢ s'annule et l'identité force **b** − **Aᵀμ** dans l'espace-ligne de rang r de **Z**, ce qui
retrouve le Théorème 1 et fait apparaître m (le nombre de marqueurs) comme la dimension pertinente.
Pour ε > 0 la même identité ne borne **pas** |S|, en accord avec le Théorème 2.

## 5. Le régime empirique entre les théorèmes

Entre les deux énoncés propres se trouve le régime que les sélectionneurs utilisent, cartographié
sur des spectres synthétiques, des alignements de **b** et des caps, et sur les panels réels :

- **L'alignement de b avec le spectre gouverne |S|, inversement.** À **G** et k fixés, placer **b**
  sur la direction propre dominante (chère en coancestrie) force un grand support (moyenne |S| ≈ 146
  sur n = 800 pour **b** le vecteur propre de tête) ; étaler **b** sur de nombreuses directions bon
  marché l'effondre à ≈ 4. Le gain cherché dans une direction chère doit être dilué sur beaucoup de
  candidats pour respecter le cap (`bound_balign.py`).
- **Aucun scalaire unique ne prédit |S| à travers les régimes.** Ni rang(**G₀**), ni le rang
  effectif / ratio de participation, ni le coût directionnel **bᵀGb/bᵀb**, ni le support de la
  direction non contrainte (**G⁻¹b**)₊ ne suit |S| sur les instances synthétiques et réelles ; la
  meilleure combinaison sans dimension, |S| ∼ (**bᵀGb/bᵀb** · 1/k)^{0,8}, n'atteint que R² ≈ 0,58 et
  rate le support relatif des panels réels (`bound_predictor.py`, `bound_lawfit.py`).
- **La courbe |S|(k) est une loi de puissance par instance, d'exposant non universel.** |S| ∼ k^{−α}
  ajuste bien chaque instance (R² ≈ 0,9–0,98 là où le support est bien résolu), avec α ≈ 1 sur les
  deux panels réels mais variant de 0,5 à 1,7 selon les spectres et alignements synthétiques
  (`bound_curve.py`).
- **Panels réels.** Blé (n = 599) et souris (n = 1814) : |S| = 25/26 et 17/26 (simplexe / optiSel
  sexé) au cap de travail, de l'ordre de la dizaine — petits et stables en n, comme le régime à
  spectre décroissant le prédit qualitativement, sans qu'aucune formule n'en fixe la valeur
  (`bound_real.py`).

Le Théorème 2 explique *pourquoi* aucune loi scalaire n'existe : il n'y a pas de borne universelle
à prédire, donc la valeur pratique est fonction de la géométrie complète (spectre, **b**, k).

## 6. Problème ouvert

La question universelle est close, négativement, par le Théorème 2 ; la question utile est
*conditionnelle* : sous des hypothèses structurelles que les vraies matrices de parenté génomique
satisfont — un spectre décroissant sans grand plancher dégénéré, et un vecteur **b** non concentré
sur les directions propres de tête — prouver |S| ≤ (quelque chose de petit, indépendant de n). Une
voie naturelle couple le gap spectral à la projection de **b** sur le sous-espace dominant ; la
moitié croissante (|S| augmentant quand k → 0) est gouvernée par la filiation classique
contributions ↔ ΔF ↔ Nₑ (Wray & Thompson 1990 ; Woolliams & Bijma 2000). C'est le point de rencontre
de l'optimisation (la dimension de face du cône perturbé) et de la génétique quantitative (le nombre
effectif de lignées contributrices).

## Évidence numérique (reproductible)

`bound_validation.py` (la borne ε=0, la n-indépendance, les balayages ridge et décroissance
spectrale, et le contre-exemple du Théorème 2) ; `bound_real.py` (blé et souris, après
`research/repro/*_export.R`) ; `bound_balign.py` (le balayage d'alignement de **b**) ;
`bound_predictor.py` et `bound_lawfit.py` (recherche de prédicteur / de loi) ; `bound_curve.py`
(l'ajustement de |S|(k) par instance). Journal complet : `research/support_bound_sketch.md`.

## Références

- Pataki, G. (1998). *Math. Oper. Res.* 23(2):339–358. DOI 10.1287/moor.23.2.339. — comptage de rang des points extrêmes.
- Markowitz, H. (1956). *Naval Res. Logist. Q.* 3:111–133. DOI 10.1002/nav.3800030110. — algorithme de la ligne critique.
- Wray, N.R. & Thompson, R. (1990). *Genet. Res.* 55(1):41–54. DOI 10.1017/S0016672300025180. — ΔF ∝ Σ(contribution²).
- Woolliams, J.A. & Bijma, P. (2000). *Genetics* 154(4):1851–1864. DOI 10.1093/genetics/154.4.1851. — contributions ↔ ΔF ↔ Nₑ.
- Yamashita, M., Mullin, T.J. & Safarina, S. (2018). *Optim. Lett.* 12(7):1683–1697. DOI 10.1007/s11590-018-1229-y. — OCS comme programme conique du second ordre.
- Waldmann, P. (2025). *Bioinform. Adv.* 5(1):vbaf259. DOI 10.1093/bioadv/vbaf259. — comptages de support empiriques indépendants.
