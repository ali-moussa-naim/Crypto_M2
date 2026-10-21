from sage.all_cmdline import *   # import sage library

import random

"""
Fonction qui génère une matrice de taille n * k à partir de H
H : matrice de taille n * k
n : entier, nombre de colonnes
k : entier, nmobre de ligne
"""
def gen_M(H, n, k):
    # On permute aléatoirement une liste de n nombres
    p = list(range(n))
    random.shuffle(p)

    # On garde seulement les (n - k) premières colonnes
    H_perm = H.matrix_from_columns(p)
    limite = n - k
    M = H_perm[:, :limite]

    return M, p

"""
Fonction qui ajoute s à droite de M
M : matrice taille (n - k) * (n - k)
s : vecteur de taille (n - k)
"""
def echelon_M_S(M, s):
    # On vérifie si M inversible
    if(not M.is_invertible()):
        print("La matrice M n'est pas inversible")
        return None

    # On ajoute le syndrome comme la colonne de droite M | s
    S = vector(GF(2), s)  
    M_S = M.augment(S.column())

    # Réduction gaussienne
    M_S = M_S.echelon_form()

    return M_S

"""
Fonction qui réitère l'expérience jusqu'à ce que M soit inversible
H : matrice, taille (n - k) * n
n : entier
k : entier
s : liste
"""
def verif_M_S(H, n, k, s):
    # On récupère M | S
    M, p = gen_M(H, n, k)
    M_S = echelon_M_S(M, s)

    # On teste jusqu'à ce qu'on tombe sur un M inversible
    while(M_S == None):
        M, p = gen_M(H, n, k)
        M_S = echelon_M_S(M, s)
    
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
    limite = factorial(n - k) / (factorial(p) * factorial(n - k - p))
    res = matrix(GF(2), limite, n-k, 0)

    # Vecteur de base pour les possibilités
    vec = vector(GF(2), n - k, [1] * p + [0] * (n - k - p))

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
def rech_e2(H, n, k, s, p, t):
    # On prend la matrice de toutes les possibilités pour e2
    mat_e2 = gen_e2_p(n, k, p)

    # On récupère la matrice M_S et permutation
    M_S, l_perm = verif_M_S(H, n, k, s)

    # On reconstruit la matrice Q 
    limite_Q = n - k
    Q = matrix(GF(2), n - k, n - k)
    l_perm = l_perm[limite_Q:]
    Q = H.matrix_from_columns(l_perm)
    
    # On récupère s_tilde
    s_tilde = M_S.column(n - k)

    continuer = True
    i = 0
    while(continuer and i < (n - k)):
        e2 = mat_e2[i]
        vecteur = vector(GF(2), n - k)
        vecteur = s_tilde - Q * e2
        
        limite_poids = t - p
        if(poids_hamming(vecteur) <= limite_poids):
            continuer = false
        
        i += 1
    
    print(vecteur)
    vecteur_e = vector(GF(2), list(s_tilde) + list(e2))
    print(vecteur_e)
    

    return Q



n, k, t = (32, 16, 3)
syndrome = [0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0]
H = [[0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1], [0, 1, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1], [0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1], [0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1], [1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0], [1, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1], [1, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 1], [0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 1], [0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0], [1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0], [1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0], [0, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 0, 1], [0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0], [0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1], [0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0], [0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0]]
p = 2


H = matrix(GF(2), H)

permutation = rech_e2(H, n, k, syndrome, p, t)
print(permutation)




