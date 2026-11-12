from itertools import permutations

def build_goppa_matrices(g, L):
    """
    Construit manuellement les matrices binaires de parité H et génératrice G
    d'un code de Goppa Gamma(L, g) sans utiliser codes.GoppaCode.
    """
    F = g.base_ring()
    m = F.degree()
    t = g.degree()
    n = len(L)

    # 1. Matrice de parité HH sur le corps étendu F = GF(2^m) (taille t x n)
    HH = matrix(F, [[alpha^i / g(alpha) for alpha in L] for i in range(t)])

    # 2. Projection binaire sur GF(2) (taille (m*t) x n)
    H = matrix(GF(2), m * t, n)
    for i in range(t):
        for j in range(n):
            rep_bin = list(HH[i, j])
            for aux in range(m):
                H[i * m + aux, j] = rep_bin[aux]

    # 3. Mise sous forme échelonnée réduite et suppression des éventuelles lignes nulles
    H = H.echelon_form()
    r = H.rank()
    H = H[:r, :]

    # 4. Déduction de la matrice génératrice G (taille (n - r) x n)
    if H.pivots() == tuple(range(r)):
        A = H.matrix_from_columns(range(r, n))
        G = block_matrix([[A.transpose(), matrix.identity(GF(2), n - r)]])
    else:
        G = H.right_kernel().basis_matrix()

    return H, G

"""
Décodage par syndrome d'un code quelconque
"""
def minimum_distance(G):
    k, n = G.dimensions()
    min_w = n
    V = VectorSpace(GF(2), k)
    
    for m in V:
        if not m.is_zero():
            w = (m * G).hamming_weight()
            if w < min_w:
                min_w = w
    return min_w

def build_syndrome_table(H, t):
    dic_erreur = {}
    r_dim, n = H.dimensions()
    
    # Pour chaque poids w de 0 à t inclus
    for w in range(t + 1):
        # On choisit w positions parmi n pour placer les '1'
        for positions in itertools.combinations(range(n), w):
            e = vector(GF(2), [0] * n)
            for pos in positions:
                e[pos] = 1
            s = H * e
            dic_erreur[tuple(s)] = e
            
    return dic_erreur

def syndrome_decode(r, H, syndrome_table):
    s = tuple(H * r)
    if s in syndrome_table:
        e = syndrome_table[s]
        return r + e  # Renvoie le mot de code corrigé (ou return e si tu veux juste l'erreur)
    return None



# --- TEST DE LA FONCTION ---
m = 4
t = 3
F.<a> = GF(2^m)
R.<z> = PolynomialRing(F)

g = R.irreducible_element(t)
L = [b for b in F.list() if g(b) != 0]

H, G = build_goppa_matrices(g, L)

"""print("H :\n", H)
print("\nG :\n", G)

# Vérification d'orthogonalité (G * H^T = 0) et comparaison avec SageMath natif
assert G * H.transpose() == 0, "Erreur : G et H ne sont pas orthogonaux !"
assert G.row_space() == codes.GoppaCode(g, L).generator_matrix().row_space(), "Code différent !"
print("\nSUCCÈS : H et G sont exactes et orthogonales !")"""

print(H)
print("\n")
err = build_syndrome_table(H, t)
print(err)
