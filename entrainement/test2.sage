import random
from sage.all import Permutation, matrix

"""
KEM de Niederretier
"""

def keygen(g, L):
    # Génération de la matrice de parité H (r x n)
    C = codes.GoppaCode(g, L)
    H = C.parity_check_matrix()
    r, n = H.dimensions()

    # Tirage de la permutation jusqu'à obtenir une forme systématique [I_r | M]
    while True:
        l_perm = [(i + 1) for i in range(n)]
        random.shuffle(l_perm)
        Perm = Permutation(l_perm).to_matrix().change_ring(GF(2))
        
        H1 = H * Perm
        H_pub = H1.echelon_form()
        
        if H_pub.matrix_from_columns(range(r)) == identity_matrix(GF(2), r):
            break

    # Extraction de S^-1 (les r premières colonnes de H * Perm)
    S_inv = H1.matrix_from_columns(range(r))

    # On inclut H_pub dans sk pour la vérification finale de decaps
    return H_pub, (H_pub, S_inv, H, Perm, g, L)


def encaps(H_pub, t):
    r, n = H_pub.dimensions()

    # Génération du vecteur d'erreur secret e de poids exact t
    l_err = list(range(n))
    random.shuffle(l_err)
    e = vector(GF(2), [0] * n)
    for indice in l_err[:t]:
        e[indice] = 1

    # Calcul de la capsule (syndrome) c = H_pub * e
    c = H_pub * e

    return c, e


def compute_synd(y, g2_poly, L):
    res = P(0)
    for i in range(len(L)):
        if y[i] != 0:
            alpha_i = L[i]
            res += (F(y[i]) / g2_poly(alpha_i)) * ((g2_poly(z) - g2_poly(alpha_i)) // (z - alpha_i))
    return res


def rebuild_err(loc, L):
    n = len(L)
    err = vector(GF(2), [0] * n)
    for i in range(n):
        if loc(L[i]) == 0:
            err[i] = 1
    return err


def solve_eq(g2_poly, s, t):
    if s == 0:
        return P(1)
    
    u1, v1 = P(1), P(0)
    u2, v2 = P(0), P(1)
    a, b = s, g2_poly

    while b.degree() >= t:
        q, r = a.quo_rem(b)
        u1, u2 = u2, u1 - q * u2
        v1, v2 = v2, v1 - q * v2
        a, b = b, r

    if u2(0) != 0:
        u2 = u2 / u2(0)
    return u2


def decaps(c, sk):
    H_pub, S_inv, H, Perm, g, L = sk
    r, n = H.dimensions()
    t = g.degree()

    # 1. Retrouver le syndrome dans le code de Goppa d'origine : s_goppa = S^-1 * c
    s_goppa = S_inv * c

    # 2. Trouver un antécédent y de taille n tel que H * y = s_goppa (sans crash m4ri)
    pivots_H = H.pivots()
    H_sq = H.matrix_from_columns(pivots_H)
    y_pivots = H_sq.inverse() * s_goppa

    y = vector(GF(2), [0] * n)
    for idx, val in zip(pivots_H, y_pivots):
        y[idx] = val

    # 3. Décoder y pour retrouver l'erreur permutée e' = Perm * e
    g2_poly = g * g
    s = compute_synd(y, g2_poly, L)
    loc = solve_eq(g2_poly, s, t)
    e_prime = rebuild_err(loc, L)

    # 4. Retirer la permutation : e = Perm^-1 * e'
    e_rec = Perm.inverse() * e_prime

    # 5. Vérification de sécurité KEM (rejet explicite/implicite si invalide)
    if e_rec.hamming_weight() != t or H_pub * e_rec != c:
        return None

    return e_rec


# --- TEST DE VALIDATION (20 encapsulations) ---
F.<a> = GF(2^6)
P.<z> = PolynomialRing(F)
t = 4
g = P.irreducible_element(t)
L = [x for x in F.list() if g(x) != 0]

H_pub, sk = keygen(g, L)

for i in range(20):
    c, e = encaps(H_pub, t)
    e_rec = decaps(c, sk)
    assert e_rec == e, f"Erreur d'encapsulation au test {i}"

print("SUCCÈS : Les 20 cycles KEM Niederreiter ont réussi à 100 % !")