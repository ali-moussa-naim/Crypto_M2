import random

def rs_encode(m, L):
    F = m.base_ring()
    R.<X> = PolynomialRing(F)
    
    # Construction du polynôme P(X) de degré < k et évaluation sur L
    P = R(list(m))
    return vector(F, [P(alpha) for alpha in L])


def rs_decode_bw(r, L, k, t):
    F = r.base_ring()
    R.<X> = PolynomialRing(F)
    n = len(L)

    # 1. Construction du système linéaire M * x = b
    # Inconnues : [q_0, ..., q_{k+t-1}, e_0, ..., e_{t-1}] (k + 2t colonnes)
    lignes_M = []
    b_list = []

    for i in range(n):
        alpha = L[i]
        ri = r[i]
        
        # Partie Q(alpha) : [1, alpha, ..., alpha^{k+t-1}]
        ligne_Q = [alpha^j for j in range(k + t)]
        
        # Partie -r_i * E_tronque(alpha) : [-r_i, -r_i*alpha, ..., -r_i*alpha^{t-1}]
        ligne_E = [-ri * (alpha^j) for j in range(t)]
        
        lignes_M.append(ligne_Q + ligne_E)
        b_list.append(ri * (alpha^t))

    M = matrix(F, lignes_M)
    b = vector(F, b_list)

    # 2. Résolution du système linéaire sur F
    try:
        sol = M.solve_right(b)
    except ValueError:
        # Système incompatible (trop d'erreurs)
        return None

    # 3. Reconstruction de Q(X) et E(X) (unitaire de degré t)
    Q = R(list(sol[:k + t]))
    E = R(list(sol[k + t:]) + [1])

    # 4. Division euclidienne Q(X) = P(X) * E(X)
    P_rec, reste = Q.quo_rem(E)

    if reste != 0 or P_rec.degree() >= k:
        return None

    # padded_list(k) garantit un vecteur de taille exactement k
    return vector(F, P_rec.padded_list(k))


# --- MAIN DE TEST ---
# Paramètres : F_17, code RS [12, 6, 7] corrigeant jusqu'à t = 3 erreurs
F = GF(17)
n = 12
k = 6
t = (n - k) // 2
L = [F(i) for i in range(n)]

print(f"Code Reed-Solomon sur GF(17) : n = {n}, k = {k}, t = {t}")

# Test sur 20 messages aléatoires avec 0, 1, 2 ou 3 erreurs
for test_idx in range(20):
    m = random_vector(F, k)
    c = rs_encode(m, L)
    
    # Injection d'un nombre aléatoire d'erreurs w <= t
    nb_erreurs = random.randint(0, t)
    r = copy(c)
    positions = random.sample(range(n), nb_erreurs)
    
    for pos in positions:
        # Ajout d'une erreur non nulle dans F_17*
        err_val = F.random_element()
        while err_val == 0:
            err_val = F.random_element()
        r[pos] += err_val

    m_rec = rs_decode_bw(r, L, k, t)
    assert m_rec == m, f"Échec au test {test_idx} avec {nb_erreurs} erreurs !"

print("SUCCÈS : Les 20 messages bruités ont été parfaitement décodés par Berlekamp-Welch !")