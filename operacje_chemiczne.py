def wzor_lancucha_aminokwasow(x):
    Smiles={
        "S":"N[C@H](C(O)=O)CO",
        "F":"NC(C(O)=O)CC1=CC=C(C=C1)",
        "L":"N[C@H](C(O)=O)CC(C)C",
        "Y":"N[C@H](C(O)=O)CC1=CC=C(O)C=C1",
        "C":"NC(C(O)=O)CS",
        "W":"N[C@@H](C(O)=O)CC1C(=CC=CC=2)C2NC=1",
        "P":"NC(C(O)=O)C",
        "H":"NC(C(O)=O)CC1NC=NC=1",
        "Q":"NC(C(O)=O)CCC(=O)N",
        "R":"NC(C(O)=O)CCCN=C(N)N",
        "I":"NC(C(O)=O)CC(C)C",
        "T":"N[C@H](C(O)=O)[C@H](O)C",
        "N":"NC(C(O)=O)CC(=O)N",
        "K":"N[C@H](C(O)=O)CCCCN",
        "V":"N[C@H](C(O)=O)C(C)C",
        "A":"NC(C(O)=O)C",
        "D":"NC(C(O)=O)CC(=O)O",
        "E":"NC(C(O)=O)CCC(=O)O",
        "G":"NC(C(O)=O)",
        "M":"N[C@H](C(O)=O)CCSC"
    } #słownik który każdemu kodonowi ich wzór chemiczny
    PozFirC={
        "S":9,
        "F":5,
        "L":9,
        "Y":9,
        "C":5,
        "W":10,
        "P": 5,
        "H": 5,
        "Q": 5,
        "R": 5,
        "I": 5,
        "T": 9,
        "N": 5,
        "K": 9,
        "V": 9,
        "A": 5,
        "D": 5,
        "E": 5,
        "G": 5,
        "M":9
    } # słownik który przechowuje pozycje grupy OH każdego kodonu we wzorach smiles
    slowo=""
    liczba=0
    for znak in x:
        if znak=="*": #jeśli kodon stop, to przestań kodować
            return slowo+"*"
        slowo=slowo[0:liczba]+Smiles[znak]+slowo[liczba+1:] # zastępuje grupę OH następnym związkiem
        liczba += PozFirC[znak]
    return slowo

def rozklad_na_bialka(lanc): #funkcja przyjmująca łańcuch kodonów, a zwracająca kandydatów na białka w tablicy
    bialka=[]

    j=0
    czy_bialko=False
    for i in lanc: #for który każdy kodon podłacza do danego białka
        if (i == "*"): #kodony stop oznaczone są w biopythonie *
            czy_bialko = False
            j += 1
        if (czy_bialko):
            bialka[j] = bialka[j] + i
        if(i=="M"):
            bialka.append("M")
            czy_bialko=True
    return bialka