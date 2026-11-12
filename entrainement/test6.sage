import random
from sage.all import Permutation

"""
Cryptanalyse : Attaque par re-chiffrement / Messages liés 
(Brier-Clavier ou attaque sur McEliece sans bourrage)
"""

def brier_attack(G_pub, c1_prime, c2_prime, t, max_iter=1000):
    k, n = G_pub.dimensions()
    
    # 1. Somme des deux chiffrés : c1' + c2' = e1 + e2
    diff = c1_prime + c2_prime

    # 2. Garder uniquement les indices où diff[i] == 0 (quasi sans erreur)
    U = [i for i in range(n) if diff[i] == 0]

    # 3. Chercher un ensemble d'information inversible dans U
    for nb_iter in range(1, max_iter + 1):
        I = random.sample(U, k)
        G_sub = G_pub.matrix_from_columns(I)

        if not G_sub.is_invertible():
            continue

        c_sub = vector(GF(2), [c1_prime[i] for i in I])
        m_cand = c_sub * G_sub.inverse()

        # Vérification : l'erreur reconstituée sur c1' doit être de poids exact t
        e1_cand = c1_prime + m_cand * G_pub
        if e1_cand.hamming_weight() == t:
            return m_cand, nb_iter

    return None, max_iter


# --- MAIN DE TEST (SUR UN CODE DE GOPPA RÉEL m=6, t=4) ---
# Remarque : Prange prendrait plusieurs milliers d'itérations sur ces paramètres,
# alors que Brier va le casser instantanément !
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
    return S * G * Perm

def encrypt(m, G_pub, t):
    k, n = G_pub.dimensions()
    e = vector(GF(2), [0] * n)
    for idx in random.sample(range(n), t):
        e[idx] = 1
    return m * G_pub + e

F.<a> = GF(2^6)
P.<z> = PolynomialRing(F)
t = 4
g = P.irreducible_element(t)
L = [x for x in F.list() if g(x) != 0]

G_pub = keygen(g, L)
k, n = G_pub.dimensions()

# Alice chiffre DEUX FOIS le même message m avec deux erreurs aléatoires distinctes
m = random_vector(GF(2), k)
c1_prime = encrypt(m, G_pub, t)
c2_prime = encrypt(m, G_pub, t)

m_rec, iters = brier_attack(G_pub, c1_prime, c2_prime, t)
assert m_rec == m, "Échec de l'attaque de Brier !"
print(f"SUCCÈS : Message retrouvé sur un code [64, 40, 9] en seulement {iters} itération(s) !")