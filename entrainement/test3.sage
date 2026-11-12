import random
from sage.all import Permutation

"""
Attaque par ensemble d'informations, algorithme de prange
"""

def keygen(g, L):
    C = codes.GoppaCode(g, L)
    G = C.generator_matrix().echelon_form()
    k, n = G.dimensions()

    S = random_matrix(GF(2), k, k)
    while not S.is_invertible():
        S = random_matrix(GF(2), k, k)

    l_perm = [i + 1 for i in range(n)]
    random.shuffle(l_perm)
    Perm = Permutation(l_perm).to_matrix().change_ring(GF(2))

    G_prime = S * G * Perm
    return G_prime, (S, G, Perm, g, L)


def encrypt(m, G_pub, t):
    k, n = G_pub.dimensions()
    l_err = random.sample(range(n), t)
    e = vector(GF(2), [0] * n)
    for indice in l_err:
        e[indice] = 1
    return m * G_pub + e


def prange_attack(G_pub, c_prime, t, max_iter=10000):
    k, n = G_pub.dimensions()

    for nb_iter in range(1, max_iter + 1):
        # 1. Tirer aléatoirement k positions (un ensemble d'information candidat)
        I = random.sample(range(n), k)

        # 2. Extraire la sous-matrice carrée k x k correspondante
        G_sub = G_pub.matrix_from_columns(I)

        # 3. Vérifier l'inversibilité sur GF(2)
        if not G_sub.is_invertible():
            continue

        # 4. Extraire le sous-vecteur c'_I et résoudre m_cand = c'_I * (G'_I)^-1
        c_sub = vector(GF(2), [c_prime[i] for i in I])
        m_cand = c_sub * G_sub.inverse()

        # 5. Calculer l'erreur associée et vérifier son poids de Hamming
        e_cand = c_prime + m_cand * G_pub
        if e_cand.hamming_weight() == t:
            return m_cand, e_cand, nb_iter

    return None, None, max_iter


# --- 2. VALIDATION SUR UN PETIT CODE DE GOPPA (m=5, t=2) ---
F.<a> = GF(2^5)
P.<z> = PolynomialRing(F)
t = 2
g = P.irreducible_element(t)
L = [x for x in F.list() if g(x) != 0]

G_pub, sk = keygen(g, L)
k, n = G_pub.dimensions()
print(f"Paramètres du code : n = {n}, k = {k}, t = {t}")

m = random_vector(GF(2), k)
c_prime = encrypt(m, G_pub, t)

m_rec, e_rec, iters = prange_attack(G_pub, c_prime, t)
assert m_rec == m, "Échec de l'attaque de Prange !"
print(f"Attaque réussie sur un essai en {iters} itérations !\n")


# --- 3. ÉTUDE EXPÉRIMENTALE VS COMPLEXITÉ THÉORIQUE (50 essais) ---
nb_tests = 50
total_iters = 0

for _ in range(nb_tests):
    m_test = random_vector(GF(2), k)
    c_test = encrypt(m_test, G_pub, t)
    m_found, _, nb_iter = prange_attack(G_pub, c_test, t)
    assert m_found == m_test
    total_iters += nb_iter

moyenne_exp = float(total_iters / nb_tests)
n_theorique = float(binomial(n, t) / binomial(n - k, t))

print(f"Moyenne expérimentale ({nb_tests} essais) : {moyenne_exp:.2f} itérations")
print(f"Espérance théorique combinatoire   : {n_theorique:.2f} itérations")