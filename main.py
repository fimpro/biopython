from Bio.Seq import Seq
import time
import sys
import customtkinter as ctk
from PIL import ImageTk,Image
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
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
    }
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
    }
    slowo=""
    liczba=0
    for znak in x:
        if znak=="*":
            return slowo+"*"
        slowo=slowo[0:liczba]+Smiles[znak]+slowo[liczba+1:]
        liczba += PozFirC[znak]
    return slowo
def interfejs(oknoR, okno): #funkcja obsługująca główny interfejs programu
    okno.destroy()
    okno = ctk.CTkFrame(oknoR)
    okno.pack()
    seq=Seq("AUGGAACGCGAACCCUAC")
    translated=seq.translate()
    mol = Chem.MolFromSmiles(wzor_lancucha_aminokwasow(translated))
    napis=ctk.CTkLabel(okno, text="Wybierz Opcje")
    plikButton = ctk.CTkButton(okno, text="Dane z Pliku", command=wczytaj_z_pliku)
    recznieButton = ctk.CTkButton(okno, text="Dane Ręcznie", command=lambda:wczytaj_recznie(oknoR,okno,wzorpoczatkowy="CACUAC"))
    instrukcjeButton = ctk.CTkButton(okno, text="Instrukcja",command= instrukcja)
    opcjeButton = ctk.CTkButton(okno, text="Opcje", command=interfaceOpcje)
    img1=ImageTk.PhotoImage(Draw.MolToImage(mol, size=(500, 500), kekulize=True, wedgeBonds=True, fitImage=True))
    obraz=ctk.CTkLabel(okno,image=img1)
    obraz.grid(row=0, column=1, rowspan=10)
    plikButton.grid(row=1, column=0)
    recznieButton.grid(row=2, column=0)
    instrukcjeButton.grid(row=3, column=0)
    opcjeButton.grid(row=4, column=0)
    napis.grid(row=0, column=0)
    oknoR.mainloop()
def wczytaj_z_pliku(): #funkcja wczytująca dane z pliku
    print("wczytaj_z_pliku")
def wczytaj_recznie(oknoR,okno,wzorpoczatkowy=""): #funkcja wczytująca dane z "palca"
    okno.destroy()
    okno = ctk.CTkFrame(oknoR)
    okno.pack()
    genom=Seq(wzorpoczatkowy)
    mol = Chem.MolFromSmiles(wzor_lancucha_aminokwasow(genom.translate()))
    img = ImageTk.PhotoImage(Draw.MolToImage(mol, size=(500, 500), kekulize=True, wedgeBonds=True, fitImage=True))
    obraz = ctk.CTkLabel(okno, image=img)
    obraz.grid(row=1, column=0,rowspan=5, columnspan=5)
    opis=ctk.CTkLabel(okno, text="Wpisz kod nici Rna",width=100)
    wejscie = ctk.CTkEntry(okno, width=398)
    wejscie.insert(0, wzorpoczatkowy)
    generujButton = ctk.CTkButton(okno, text="Generuj", width=140,command=lambda:generuj(okno,wejscie))
    wykresyButton = ctk.CTkButton(okno, text="Wykresy", width=140,command=lambda:wykresy(okno,oknoR,wejscie))
    wrocButton = ctk.CTkButton(okno, text="Wroc", width=140,command=lambda:interfejs(oknoR, okno))

    opis.grid(row=0,column=0)
    wejscie.grid(row=0,column=1, columnspan=4 )
    generujButton.grid(row=0,column=5)
    wykresyButton.grid(row=1,column=5)
    wrocButton.grid(row=2,column=5)
    okno.mainloop()
def generuj(okno,wejscie):
    t0 = time.time()
    genom=Seq(wejscie.get())
    kodon=genom.translate()
    Smiles=wzor_lancucha_aminokwasow(kodon)
    mol = Chem.MolFromSmiles(Smiles)
    img = ImageTk.PhotoImage(Draw.MolToImage(mol, size=(500, 500), kekulize=True, wedgeBonds=True, fitImage=True))
    obraz = ctk.CTkLabel(okno,image=img)
    obraz.grid(row=1, column=0, columnspan=5)
    print (time.time() - t0)
    okno.mainloop()
def wykresy(okno, oknoR,wejscie):
    wzor=wejscie.get()
    okno.destroy()
    okno = ctk.CTkFrame(oknoR)
    okno.pack()
    napis=ctk.CTkLabel(okno,text=wzor, width=100)
    returnButton=ctk.CTkButton(okno, text="Wróć",command=lambda:wczytaj_recznie(oknoR,okno,wzorpoczatkowy=wzor))
    napis.grid(row=0, column=0)
    returnButton.grid(row=1, column=0)
    okno.mainloop()
def interfaceOpcje(): #funkcja wczytująca interface opcji
    print("interfaceOpcje")
def instrukcja(): #funkcja wczytująca instruckje użytkowania
    print("instrukcja")
if __name__ == '__main__':
    oknoR = ctk.CTk()
    okno = ctk.CTkFrame(oknoR)
    okno.pack()
    oknoR.minsize(height=500, width=650)
    oknoR.title("Katolgenom")
    ctk.set_appearance_mode("System")
    ctk.set_default_color_theme("blue")
    interfejs(oknoR, okno)



