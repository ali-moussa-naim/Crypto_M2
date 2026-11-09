import random

# renvoie le syndrome du mot c
def compute_synd(c):
    # On vérifie que la taille de c et de L soit la même
    if(len(c) != len(L)):
        print("Le mot c et la Liste L doivent avoir la même taille")
        return None

    res = P(0)
    for i in range(len(L)):
        if(c[i] != 0):
            alpha_i = L[i]
            res += (c[i] / g2(alpha_i)) * ((g2(z) - g2(alpha_i)) // (z - alpha_i))

    return res

# renvoie un mot de longueur n et de poids t
def gen_err(n,t):
    liste_poids = list(range(n))
    random.shuffle(liste_poids)
    mot = [0] * n

    for i in range(t):
        indice = liste_poids[i]
        mot[indice] = 1
    
    return mot

def rebuild_err(loc,n):
# renvoie le mot erreur de longueur n
# a partir du polynome localisateur

    err=CS([0]*n)
    i = 0
    while i < len(L):
        if loc(L[i]) == 0:
            err[i] = 1
        i += 1
    return err
    
# resolution l'equation cle des codes de goppa
# et renvoie le polynome localisateur
def solve_eq(g,s,t):
    # Si syndrome nul
    if(s == 0):
        return P(1)
    
    u1, v1 = P(1), P(0)
    u2, v2 = P(0), P(1)
    a, b = s, g

    while(b.degree() >= t):
        q, r = a.quo_rem(b)
        u1, u2 = u2, u1 - q * u2
        v1, v2 = v2, v1 - q * v2

        a, b = b, r

    return u2

# renvoie le mot de code C
# associe à m
def encode(m):
    c = m * G
    return c


F.<a>=GF(2^10)
P.<z>=PolynomialRing(F)
g = z^50+z^4+z^3+z^2+1
g2 = g*g
L = [a for a in F.list() if g(a) != 0]
C=codes.GoppaCode(g,L)
# espace des messages GF(2)^k
MS = GF(2)^C.dimension()
# espace des mots de code GF(2)^n
CS = GF(2)^C.length()
# Matrice génératrice de CS
G = C.generator_matrix()
m = MS.random_element()

for cpt in range(10):
    # on tire au hasard un message m
    c = encode(m)
    # on genere l'erreur
    e = gen_err(C.length(),g.degree())
    # on simule un envoi avec erreur
    e = vector(e)
    r = c + e 
    synd = compute_synd(r)
    loc = solve_eq(g2,synd,g.degree())
    erreur = rebuild_err(loc,C.length())
    print(cpt,':',erreur == e)


