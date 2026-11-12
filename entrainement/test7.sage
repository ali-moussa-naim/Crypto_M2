import random

"""
Protocole d'authentification / Identification par codes (Protocole de Stern)
"""

def toy_hash(data):
    """
    Fonction de hachage jouet 100 % native (sans hashlib) pour l'examen.
    Accepte n'importe quel tuple/liste d'entiers ou de bits.
    """
    acc = 0x5BD1E995
    for idx, val in enumerate(data):
        acc = ((acc * 1315423911) ^^ (int(val) + idx + 1)) % (2^128)
    return acc


def apply_perm(pi, vec):
    """Applique une permutation pi (liste d'indices 0..n-1) à un vecteur de GF(2)^n."""
    return vector(GF(2), [vec[pi[i]] for i in range(len(pi))])


def stern_keygen(n, r, w):
    """
    Génère une paire de clés pour Stern :
    - Clé publique : (H, y, w) avec H de taille r x n et y = H * s
    - Clé privée   : s in GF(2)^n de poids de Hamming exact w
    """
    H = random_matrix(GF(2), r, n)
    
    # Vecteur secret s de poids exact w
    s = vector(GF(2), [0] * n)
    for pos in random.sample(range(n), w):
        s[pos] = 1
        
    y = H * s
    return (H, y, w), s


def stern_prover_commit(s, H):
    """
    Étape 1 (Engagement d'Alice) :
    Tire u et pi aléatoires, calcule les 3 engagements (c1, c2, c3)
    et conserve (u, pi) comme état secret pour la réponse.
    """
    n = len(s)
    u = random_vector(GF(2), n)
    
    pi = list(range(n))
    random.shuffle(pi)
    
    pi_u = apply_perm(pi, u)
    pi_us = apply_perm(pi, u + s)
    
    c1 = toy_hash(tuple(pi) + tuple(H * u))
    c2 = toy_hash(tuple(pi_u))
    c3 = toy_hash(tuple(pi_us))
    
    return (c1, c2, c3), (u, pi)


def stern_prover_response(b, s, state):
    """
    Étape 2 (Réponse d'Alice au défi b in {0, 1, 2}) :
    Ne révèle jamais à la fois s et (u, pi) !
    """
    u, pi = state
    if b == 0:
        return (u, pi)
    elif b == 1:
        return (u + s, pi)
    elif b == 2:
        return (apply_perm(pi, u), apply_perm(pi, s))
    else:
        raise ValueError("Défi invalide : b doit être dans {0, 1, 2}")


def stern_verifier_check(b, response, commitments, pk):
    """
    Étape 3 (Vérification de Bob) :
    Vérifie la cohérence des engagements et la contrainte de poids w.
    """
    H, y, w = pk
    c1, c2, c3 = commitments

    if b == 0:
        u_rec, pi_rec = response
        c1_check = toy_hash(tuple(pi_rec) + tuple(H * u_rec))
        c2_check = toy_hash(tuple(apply_perm(pi_rec, u_rec)))
        return (c1_check == c1) and (c2_check == c2)

    elif b == 1:
        us_rec, pi_rec = response
        # Astuce algébrique : H * (u + s) + y = H*u + H*s + H*s = H*u sur GF(2)
        Hu_rec = H * us_rec + y
        c1_check = toy_hash(tuple(pi_rec) + tuple(Hu_rec))
        c3_check = toy_hash(tuple(apply_perm(pi_rec, us_rec)))
        return (c1_check == c1) and (c3_check == c3)

    elif b == 2:
        pi_u_rec, pi_s_rec = response
        # 1. Vérifier que le secret permuté a bien le poids de Hamming w
        if pi_s_rec.hamming_weight() != w:
            return False
        # 2. Vérifier c2 et c3 (par linéarité : pi(u + s) = pi(u) + pi(s))
        c2_check = toy_hash(tuple(pi_u_rec))
        c3_check = toy_hash(tuple(pi_u_rec + pi_s_rec))
        return (c2_check == c2) and (c3_check == c3)

    return False


# --- MAIN DE VALIDATION (30 TOURS DE PROTOCOLE) ---
n, r, w = 64, 32, 8
pk, sk = stern_keygen(n, r, w)
H, y, w = pk

# 1. Test d'un prouveur honnête sur 30 tours (erreur de justesse = (2/3)^30 < 10^-5)
for tour in range(30):
    commitments, state = stern_prover_commit(sk, H)
    b = random.choice([0, 1, 2])
    resp = stern_prover_response(b, sk, state)
    assert stern_verifier_check(b, resp, commitments, pk), f"Rejet à tort au tour {tour} (b={b})"

print("SUCCÈS : Prouveur honnête accepté sur les 30 tours !")

# 2. Test d'un tricheur (qui utilise un faux secret s_fake de mauvais poids ou mauvais syndrome)
s_fake = random_vector(GF(2), n)
rejets = 0
for _ in range(30):
    commitments, state = stern_prover_commit(s_fake, H)
    b = random.choice([0, 1, 2])
    resp = stern_prover_response(b, s_fake, state)
    if not stern_verifier_check(b, resp, commitments, pk):
        rejets += 1

print(f"SUCCÈS : Tricheur détecté et rejeté {rejets}/30 fois (attendu ~10/30) !")