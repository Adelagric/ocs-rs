# Sur la taille du support de l'optimum OCS génomique

> Note autoportante — une preuve, un contre-exemple et le régime empirique — pour le nombre
> d'individus qu'active un optimum de sélection à contribution optimale. Note compagne du
> solveur support-first (`MANUSCRIPT_fr.tex`). Version française ; note de référence en anglais
> (`support_bound_sketch.md`). Les résultats historiques sont conservés avec leur contexte ; scripts et prérequis figurent sous *Évidence numérique*.

## 1. La question

La sélection à contribution optimale (OCS) génomique résout, pour n candidats de valeurs
génétiques **b** ∈ ℝⁿ et de matrice de parenté **G = ZZᵀ/s + εI** (**Z** la matrice n×m des
génotypes centrés, s une échelle, ε un petit ridge pour la définie-positivité),

  maximiser **bᵀc**  s.c.  **Ac = d**,  **0 ≤ c ≤ u**,  **cᵀGc ≤ k**,    (OCS)

où **A** est la matrice de budget q×n (q = 1 pour le simplexe 𝟙ᵀc = 1 ; q = 2 pour la forme
sexuée Σ_mâles = Σ_femelles = ½). Empiriquement l'optimum **c\*** active très peu de candidats :
son support total S = {i : c\*ᵢ > 0} était ≈ 15–30 dans les essais historiques à cap lâche. Aux points de proxy opérationnels du manuscrit, il va de quelques dizaines à quelques centaines (souris : 47/89/182 à ΔF = 2/1/0,5 %). Le balayage structuré en n à proxy fixé donne environ 80–103 ; ces observations ne prouvent pas une borne universelle. Qu'est-ce qui borne |S| ?

La réponse est un **encadrement**. Sans plafonds effectifs et sans ridge, une borne existentielle indépendante de n
existe (Théorème 1). Dès que le ridge est présent, il n'y a prouvablement aucune borne
universelle — le support peut valoir toute la population (Théorème 2). Le régime pratique entre
les deux est fixé par la géométrie conjointe de (spectre, **b**, cap k) ; les prédicteurs scalaires testés ne sont pas fiables entre régimes et on le cartographie empiriquement (§5). Une identité KKT conditionnelle (§4) décrit les contributions libres.

## 2. Théorème 1 — la borne ε = 0

**Théorème 1 (contributions non plafonnées).** *Supposons ε = 0, s > 0, un ensemble réalisable non vide et les budgets simplex ou sexués, de sorte que {Ac = d, c ≥ 0} est borné. Sans plafonds supérieurs effectifs, il existe un optimum vérifiant*
**|S| ≤ rang([A; Zᵀ; bᵀ]) ≤ q + r + 1**, *où r = rang(Z) ≤ m.*

*Démonstration.* Partons d'un optimum atteint **c\*** et figeons son image de bas rang et sa valeur objective dans

  P = { **c ≥ 0** : **Ac = d**, **Zᵀc = Zᵀc\***, **bᵀc = bᵀc\*** }.

Tout point de P a la même parenté ‖**Zᵀc\***‖²/s et le même objectif : il est donc optimal. P est non vide et borné. La matrice M = [A; Zᵀ; bᵀ] est de rang au plus q + r + 1. Un sommet ayant plus de rang(M) coordonnées positives admettrait une direction non nulle du noyau portée par ces coordonnées ; deux petites perturbations opposées conserveraient les égalités et la positivité, contredisant son caractère extrême. Un sommet de P fournit la borne. Aucune dualité ni qualification de contraintes n'est requise. ∎

**Extension aux plafonds finis.** En intersectant P avec **c ≤ u**, un sommet vérifie seulement **|F| ≤ rang(M) ≤ q + r + 1** pour F = {i : 0 < cᵢ < uᵢ}, par le même argument de perturbation. Les coordonnées U₊ = {i : cᵢ = uᵢ > 0} comptent également dans le support, donc **|S| = |F| + |U₊| ≤ rang(M) + |U₊|**. Des plafonds finis arbitraires interdisent une borne du support total indépendante de n : q = 1, r = 0, G = 0, k = 0 et uᵢ = 2/n imposent au moins ceil(n/2) contributions positives à tout plan réalisable, alors que q + r + 1 = 2. Un b strictement décroissant atteint ce minimum en un sommet LP optimal plafonné.

Ce résultat porte sur l'existence d'un sommet d'une *tranche optimale*. Il ne garantit pas que l'ensemble actif implémenté renvoie ce sommet, ni qu'un solveur à point intérieur ou ADMM ne puisse rendre une solution parcimonieuse. Une extraction séparée de sommet serait nécessaire pour imposer cette garantie à une sortie arbitraire.

## 3. Théorème 2 — aucune borne universelle pour ε > 0

**Théorème 2.** *Pour ε > 0, aucune borne du support total de la forme f(q, r) indépendante de n n'existe, même sans plafonds supérieurs effectifs : le support peut valoir n.*

*Démonstration (contre-exemple).* Prenons l'instance OCS abstraite q = 1, G = εI, r = 0 et uᵢ = 1. Ce contre-exemple algébrique ne représente pas un panel de marqueurs utilisable dans une interface publique. Sur le simplexe, la contrainte est ε‖c‖² ≤ k. Cauchy–Schwarz donne ‖c‖² ≥ 1/n ; si une coordonnée est nulle, ‖c‖² ≥ 1/(n−1). Ainsi, pour

  **ε/n < k < ε/(n−1)**,

le plan uniforme est strictement réalisable et tout plan réalisable, donc tout optimum pour n'importe quel b, a un support plein. Cela prouve l'énoncé pour tout n ≥ 2. ∎

L'ancien balayage numérique (`bound_validation.py`, bloc 7, n = 300) donnait |S| = 300, 295, 256, 151, 52, 12 quand k se desserrait de 1,05× à 50× le minimum. Ce sont des observations pour un b ; 1,05× est hors de l'intervalle analytique de support plein et ne constitue pas un témoin universel. Le bloc 8 vérifie directement cet intervalle et l'obstruction des plafonds finis.

La somme Σcᵢ² mesure la concentration des contributions. Avec Nₑ = 1/Σcᵢ², ΔF = 1/(2Nₑ) est un proxy fondé sur les contributions, pas un taux mesuré de consanguinité génomique réalisée. Le petit support observé dépend de la géométrie de l'instance et n'est pas une garantie de pire cas.

## 4. Une identité KKT conditionnelle

**Proposition.** *Supposons ε > 0 et un optimum disposant de multiplicateurs KKT avec λ > 0 pour la parenté. Sur F = {i : 0 < c\*ᵢ < uᵢ}, les contributions sont affines en (bᵢ, zᵢ, aᵢ) :*

  **c\*ᵢ = α bᵢ + wᵀzᵢ + vᵀaᵢ**, *i ∈ F,*

*aᵢ désignant la colonne i de A. Pour les budgets d'incidence de sexe, le dernier terme est β_{sexe(i)}.*

*Démonstration.* Avec les multiplicateurs inférieur et supérieur ℓ,h ≥ 0, la stationnarité KKT est **b = Aᵀμ + 2λGc\* − ℓ + h**, avec ℓᵢc\*ᵢ = 0, hᵢ(uᵢ−c\*ᵢ) = 0 et λ(c\*ᵀGc\*−k) = 0. Sur F, les deux multiplicateurs de borne s'annulent. Posons y = Zᵀc\*. Comme Gc\* = Zy/s + εc\*, la division par 2λε donne α = 1/(2λε), w = −y/(sε) et v = −μ/(2λε). ∎

L'identité comporte m + q + 1 coefficients (m + 2 en simplex, m + 3 pour deux budgets sexués). Une parenté active ne suffit pas à imposer λ > 0. L'existence des multiplicateurs exige une qualification appropriée ou un certificat fourni ; elle est supposée ici. L'identité ne s'applique pas aux contributions au plafond et ne borne aucune cardinalité. Le passage ε → 0 ne retrouve pas le Théorème 1 : les multiplicateurs et supports peuvent changer ou diverger et la stationnarité ne vaut que sur les coordonnées actuellement libres. La preuve par tranche optimale est indépendante.

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
- **Panels réels.** Les expériences historiques à borne lâche sur le blé (n = 599) et la souris
  (n = 1814) donnent |S| = 25/26 et 17/26 (simplexe / forme sexuée ; `bound_real.py`).
  Deux panels de taille fixée ne démontrent pas la stabilité en n. Les supports aux cibles
  opérationnelles du manuscrit sont plus grands ; ces nombres historiques ne les remplacent pas.

Le Théorème 2 exclut une borne universelle en q et r seuls. Les balayages empiriques montrent séparément que les prédicteurs scalaires testés ne prédisent pas le support de façon fiable entre régimes.

## 6. Une borne conditionnelle : l'obstruction et la conjecture

Une borne conditionnelle reste ouverte. Si la contrainte de parenté est *redondante sur tout le polytope budgétaire*, le LP non plafonné admet un sommet optimal de support |S| ≤ rang(A) ≤ q ; avec plafonds finis, seule la borne |F| ≤ rang(A) reste valable. Une parenté non saturée en un optimum ne suffit pas : pour b constant, G = I et 1/n < k < 1/(n−1), l'optimum uniforme est strictement réalisable mais tout optimum réalisable a un support plein. Au rang exact r et à ε = 0, le Théorème 1 donne la borne existentielle non plafonnée et son extension borne les coordonnées libres.

Une borne en *rang effectif* — l'objectif réaliste, puisque les vrais spectres décroissent — ne
découle **pas** du même argument, et l'obstruction est précise. Le Théorème 1 linéarise la tranche
optimale en épinglant Zᵀc, soit *toute* l'image de G₀ ; pour un G de plein rang il faudrait épingler
les n projections uᵢᵀc, ce qui ne donne que le trivial |S| ≤ q + n + 1. On ne peut pas n'épingler que
les d directions de tête et négliger la queue, car celle-ci contribue cᵀEc = Σ_{i>d} λᵢ(uᵢᵀc)² — une
*vraie quadratique*, pas un terme de bas rang — donc la fixer n'est pas une contrainte linéaire et la
tranche n'est plus un polytope. La courbure de la queue est exactement ce qui laisse le support
s'étaler ; le rang effectif seul ne remplace pas le rang exact dans ce comptage de sommets. Une borne conditionnelle
exige donc un autre argument (perturbation spectrale).

Deux conditions qualitatives motivent un éventuel théorème conditionnel ; elles ne sont ni quantitativement spécifiées ni prouvées suffisantes :
- **(A1) Pas de grand plateau dégénéré bon marché.** Un bloc de nombreuses petites valeurs propres
  quasi-égales est un sous-espace bon marché où l'optimum s'étale ; G = εI est le cas extrême et donne
  |S| = n (Théorème 2). Un spectre régulièrement décroissant n'a pas un tel plateau.
- **(A2) b évite les directions dominantes.** Même avec un gap, b sur le vecteur propre de tête (cher)
  force la dilution et un grand support (moyenne |S| ≈ 146 sur 800 dans le balayage d'alignement).

Le contre-exemple du ridge pur et le balayage d'alignement motivent (A1) et (A2), sans prouver la nécessité de ces conditions informelles. Sous des versions quantitatives à point de fonctionnement spécifié, on conjecture un petit support ; les 15–30 historiques concernent les caps lâches, pas tout panel réel. Il serait de
l'ordre des directions propres dominantes que b excite réellement. La croissance observée quand k approche le minimum réalisable est indexée ici par le proxy de contributions ΔF = 1/(2Nₑ), distinct de la consanguinité génomique réalisée. La littérature génétique fournit le contexte (Wray & Thompson 1990 ; Woolliams & Bijma 2000).
Prouver la borne — coupler le gap spectral à la projection de b sur le sous-espace dominant, en
résistant à la courbure de la queue — est le problème ouvert, point de rencontre de l'optimisation
(la face du cône perturbé) et de la génétique quantitative (lignées effectives).

**Une prise exacte sur le support effectif.** Écrivons le nombre *effectif* de contributeurs comme le
ratio de participation PR(c\*) = 1/‖c\*‖² (égal à |S| pour des contributions uniformes, moindre quand
elles sont inégales). Au cap actif positif k > 0 avec Σc\* = 1, une identité tient : **PR(c\*) = R(c\*)/k**, où
R(c\*) = c\*ᵀGc\*/c\*ᵀc\* est le quotient de Rayleigh propre de l'optimum (car R·‖c\*‖² = c\*ᵀGc\* = k).
Donc **PR(c\*) ≤ λ₁/k** sous ces conditions de cap actif positif, et une borne du support *effectif* se réduit à borner R(c\*) — le
coût directionnel de coancestrie que l'optimum est forcé d'encaisser. C'est exactement là qu'agit
(A2) : dans les cas mesurés où b évite les directions dominantes, R(c\*) est petit (observé : R ≈ 0,02–0,4 sur les cas structurés contre λ₁ ≈ 45–143, donnant PR ≈ 2–12), tandis
que G = εI force R = ε = λ₁. Au cap actif, PR = ε/k appartient à (n−1, n) dans tout
l'intervalle de support plein ε/n < k < ε/(n−1), et tend vers n lorsque k décroît vers ε/n.
Le support total vaut déjà n dans tout cet intervalle ; PR = n exige des contributions
uniformes. Par exemple, n = 3 et k = 5ε/12 donnent |S| = 3 mais PR = 12/5.
La pièce ouverte résiduelle est désormais nette et scalaire —
borner R(c\*) a priori sous (A1)+(A2). (La cardinalité brute |S| valait historiquement ≈ 1,5–2,5× PR dans ces cas. En général PR ≤ |S| : une borne supérieure de PR ne borne pas |S| sans hypothèse supplémentaire sur les contributions.) Vérifié dans `bound_kernel.py`.

**Attaquer ce résidu.** Pour G définie positive, la relaxation simplex sans bornes inférieures ni supérieures a, dans le cas à parenté active avec λ > 0, l'expression du quotient de Rayleigh
R_relax = (βᵀG⁻¹β)/(βᵀG⁻²β) avec β = b − μ𝟙 — une moyenne des λ pondérée par βᵢ²/λᵢ², donc petite
sous (A2). Mais elle ne prédit R(c\*) que lorsque c ≥ 0 est quasi-inactif : sur les panels réels
l'optimum sans contrainte porte de grandes composantes négatives (min c_relax ≈ −1 à −2,6) et
R_relax tombe 60–300× sous R(c\*). Donc **la non-négativité est le cœur irréductible** — de fait,
|S| est la cardinalité de l'optimum d'un QP contraint en positivité, dur en général — et une preuve
doit être *genuinely c ≥ 0-aware*, pas un court-circuit par relaxation. La prise a priori (la formule
de relaxation) et l'obstacle irréductible (c ≥ 0) sont désormais tous deux identifiés.

## Évidence numérique (reproductible)

`bound_validation.py` (observations de bas rang avec ridge, qui ne testent pas le théorème de sommet à ε=0, le balayage empirique en n, les balayages ridge et décroissance
spectrale, le balayage historique du Théorème 2 et les témoins analytiques avec plafonds / support plein du bloc 8) ; `bound_real.py` (blé et souris, après
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
