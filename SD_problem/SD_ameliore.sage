from sage.all_cmdline import *   # import sage library

import random

"""
Fonction qui génère une matrice de taille (n - k) * n à partir de H
H : matrice de taille (n - k) * n
n : entier, nombre de colonnes
k : entier, dimension du code
"""
def gen_M(H, n, k):
    # On permute aléatoirement une liste de n nombres
    p = list(range(n))
    random.shuffle(p)

    # On garde les (n - k) premières colonnes
    H_perm = H.matrix_from_columns(p)
    limite = n - k
    M = H_perm[:, :limite]
    Q = H_perm[:, limite:]

    return M, Q, p

"""
Fonction qui ajoute s à droite de M
M : matrice taille (n - k) * (n - k)
s : vecteur de taille (n - k)
Q : matrice taille (n - k) * k
"""
def echelon_M_S(M, Q, s):
    # On vérifie si M inversible
    if not M.is_invertible():
        print("La matrice M n'est pas inversible")
        return None

    # On applique les opérations de Gauss à M, Q et s
    M_inv = M.inverse()
    Q_tilde = M_inv * Q
    s_tilde = M_inv * s

    return Q_tilde, s_tilde


"""
Fonction qui réitère l'expérience jusqu'à ce que M soit inversible
H : matrice, taille (n - k) * n
n : entier
k : entier
s : vecteur de taille (n - k)
"""
def verif_M_S(H, n, k, s):
    # On récupère M, Q et la permutation
    M, Q, p = gen_M(H, n, k)
    M_S = echelon_M_S(M, Q, s)

    # On teste jusqu'à ce qu'on tombe sur un M inversible
    while M_S is None:
        M, Q, p = gen_M(H, n, k)
        M_S = echelon_M_S(M, Q, s)

    return M_S, p

"""
Fonction qui caclcule le poids d'un vecteur
s : liste, le syndrome
"""
def poids_hamming(s):
    return s.hamming_weight()

"""
Fonction qui énumère tous les e2 de poids p
n : entier
k : entier
p : entier, poids hamming
"""
def gen_e2_p(n, k, p):
    # Dimension de la matrice res de toutes les possibilités
    limite = binomial(n - k, p)
    res = matrix(GF(2), limite, k, 0)

    # Vecteur de base pour les possibilités
    vec = vector(GF(2), k, [1] * p + [0] * (k - p))

    # Ajout dans la matrice res
    res[0] = vec

    # Boucle pour faire toutes les possibilités
    for i in range(1, limite):

        random.shuffle(vec)
        while(vec in res.rows()):
            random.shuffle(vec)

        res[i] = vec

    return res


"""
Fonction qui retourne e2
H : matrice (n - k) * n
n : entier
k : entier
s : liste, le syndrome
p : entier, poids hamming de e2
"""
def rech_e(H, n, k, s, p, t):
    # On prend la matrice de toutes les possibilités pour e2
    mat_e2 = gen_e2_p(n, k, p)

    # On récupère la matrice M_S, liste permutation, Q et s_tilde
    M_S, l_perm = verif_M_S(H, n, k, s)
    Q, s_tilde = M_S

    e1 = 0

    # On fait une recherche exhaustive pour e1
    for i in range(mat_e2.nrows()):
        e2 = mat_e2[i]
        tmp = s_tilde - Q * e2

        # On vérifie si on a trouvé un bon e1
        if(poids_hamming(tmp) <= (t - p)):
            e1 = tmp
            break

    # On vérifie si on a pu trouvé e1
    if(e1 == 0):
        print("On n'a pas trouvé e1")
        return None
    
    e_perm = vector(GF(2), list(e1) + list(e2))
    e = vector(GF(2), n)

    # On remet les coordonnées dans l'ordre original
    for i in range(n):
        e[l_perm[i]] = e_perm[i]

    # On vérifie les résultats
    if(H * e == s and poids_hamming(e) <= t):
        return e
    
    return None


n, k, t = (32, 16, 3)
syndrome = [0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0]
H = [[0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1], [0, 1, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1], [0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1], [0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1], [1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0], [1, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1], [1, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 1], [0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 1], [0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0], [1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0], [1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0], [0, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 0, 1], [0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0], [0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1], [0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0], [0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0]]
p = 2

H = matrix(GF(2), H)
s = vector(GF(2), syndrome)

# On répète l'algo tant que l'on a pas trouvé e
e = None
while(e == None):
    e = rech_e(H, n, k, s, p, t)

print(e)
print("\n")
if(e is not None):
    print("Vérification H * e = s:", H * e == s)





