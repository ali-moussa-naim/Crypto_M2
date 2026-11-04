import random

# renvoie le syndrome du mot c
def compute_synd(c):
    pass

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
    

def solve_eq(g,s,t):
    pass
# resolution l'equation cle des codes de goppa
# et renvoie le polynome localisateur

def encode(m):
    pass
# renvoie le mot de code C
# associe à m

F.<a>=GF(2^10)
P.<z>=PolynomialRing(F)
g = z^50+z^4+z^3+z^2+1
g2 = g*g
print("g2:",g2)
L = [a for a in F.list() if g(a) != 0]
C=codes.GoppaCode(g,L)
print("C:", C)
# espace des messages GF(2)^k
MS = GF(2)^C.dimension()
# espace des mots de code GF(2)^n
CS = GF(2)^C.length()
for cpt in range(10):
    # on tire au hasard un message m
    c = encode(m)
    # on genere l'erreur
    e = gen_err(C.length(),g.degree())
    # on simule un envoi avec erreur
    r = c + e 
    synd = compute_synd(r)
    loc = solve_eq(g2,synd,g.degree())
    erreur = rebuild_err(loc,C.length())
    print(cpt,':',erreur == e)
"""

# TEST
c = [1, 0, 0, 1, 1, 0]
M = compute_synd(c)
print(M.parity_check_matrix())
"""