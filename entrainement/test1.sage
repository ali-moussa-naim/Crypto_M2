import random

"""
Cryptosystème McEliece complet
"""

def keygen(g, L):
    C = codes.GoppaCode(g, L)
    # Mettre G sous forme échelonnée garantit une base canonique stable
    G = C.generator_matrix().echelon_form()
    k, n = G.dimensions()

    S = random_matrix(GF(2), k, k)
    while not S.is_invertible():
        S = random_matrix(GF(2), k, k)

    l_perm = [i + 1 for i in range(n)]
    random.shuffle(l_perm)
    Perm = Permutation(l_perm).to_matrix().change_ring(GF(2))

    G_prime = S * G * Perm

    # On passe G directement dans sk !
    return G_prime, (S, G, Perm, g, L)


def encrypt(m, G_pub, t):
    k, n = G_pub.dimensions()

    # Génération du vecteur erreur de poids exact t
    l_err = list(range(n))
    random.shuffle(l_err)
    e = vector(GF(2), [0] * n)
    for indice in l_err[:t]:
        e[indice] = 1

    return m * G_pub + e


def compute_synd(c, g2_poly, L):
    res = P(0)
    for i in range(len(L)):
        if c[i] != 0:
            alpha_i = L[i]
            res += (F(c[i]) / g2_poly(alpha_i)) * ((g2_poly(z) - g2_poly(alpha_i)) // (z - alpha_i))
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


def decrypt(c_prime, sk):
    S, G, Perm, g, L = sk
    g2_poly = g * g

    # 1. Retirer la permutation
    c_1_prime = c_prime * Perm.inverse()

    # 2. Décoder pour retrouver l'erreur
    s = compute_synd(c_1_prime, g2_poly, L)
    poly_loc_erreur = solve_eq(g2_poly, s, g.degree())
    erreur = rebuild_err(poly_loc_erreur, L)

    # 3. Corriger l'erreur -> c1 = (m * S) * G
    c_1 = c_1_prime + erreur

    # 4. Comme G est sous forme échelonnée réduite (echelon_form),
    # les coordonnées de m_S sont directement les bits de c_1 aux positions pivots !
    pivots = G.pivots()
    m_S = vector(GF(2), [c_1[i] for i in pivots])

    # 5. Retirer S
    m_dechiffre = m_S * S.inverse()

    return m_dechiffre


# --- TEST DE VALIDATION (20 messages) ---
F.<a> = GF(2^6)
P.<z> = PolynomialRing(F)
t = 4
g = P.irreducible_element(t)
L = [x for x in F.list() if g(x) != 0]

G_pub, sk = keygen(g, L)
k, n = G_pub.dimensions()

for i in range(20):
    m = random_vector(GF(2), k)
    c_prime = encrypt(m, G_pub, t)
    m_rec = decrypt(c_prime, sk)
    assert m == m_rec, f"Erreur au test {i}"

print("SUCCÈS : Les 20 messages ont été chiffrés et déchiffrés avec succès !")