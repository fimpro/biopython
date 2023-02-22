from rdkit import Chem
def wzor_lancucha_aminokwasow(x):
    Smiles={
        "S":"N[C@H](C(O)=O)CO",
        "F":"NC(C(O)=O)CC1=CC=C(C=C1)",
        "L":"N[C@H](C(O)=O)CC(C)C",
        "Y":"N[C@H](C(O)=O)CC1=CC=C(O)C=C1",
        "C":"NC(C(O)=O)CS",
        "W":"N[C@@H](C(O)=O)CC1C(=CC=CC=2)C2NC=1",
        "P":"N1CCCC1C(O)=O",
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
        "M":"N[C@H](C(O)=O)CCSC",
        "Q":"I"
    } #słownik który każdemu kodonowi ich wzór chemiczny
    PozFirC={
        "S":9,
        "F":5,
        "L":9,
        "Y":9,
        "C":5,
        "W":10,
        "P": 9,
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
        slowo=slowo[0:liczba]+Smiles[znak]+slowo[liczba+1:]# zastępuje grupę OH następnym związkiem
        liczba += PozFirC[znak]
    lanc=slowo[0:liczba+1]
    dlugosc=0
    for x in lanc:
        if (x=='N' or x=='C'):
            dlugosc+=1
    return slowo, dlugosc

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
        elif(i=="M"):
            bialka.append("M")
            czy_bialko=True
    return bialka
def indeks_hydrofobowy(lanc_Kodonow):  # funkcja ze słownikiem zwraca gotowe tablice do wykorzystwania w wykresie
    kd = {'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
          'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
          'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
          'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2}
    dane_y = []
    for x in lanc_Kodonow:
        dane_y.append(kd[x])
    return dane_y
def listaAtomow(wzor):
    mol=Chem.MolFromSmiles(wzor)
    for atom in mol.GetAtoms():
        print(atom.GetAtomicNum())