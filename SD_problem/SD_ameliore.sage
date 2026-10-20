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
def rech_e2(H, n, k, s, p):
    # On prend la matrice de toutes les possibilités pour e2
    mat_e2 = gen_e2_p(n, k, p)

    
    


n, k, t = (8, 4, 3)
p = 2

vec = gen_e2_p(n, k, p)
print(vec)





