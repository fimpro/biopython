import math
from customtkinter import *
import operacje_chemiczne
import main
import wykresy
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import rdCoordGen
from rdkit.Chem import Descriptors
def rysuj_dziwnie_interface(ramex, lanc_Kodonow,czy_podswietlaj):
    okno= CTkFrame(ramex)
    size = (ramex.winfo_width(), ramex.winfo_height())
    okno.grid(row=0, column=0, rowspan=10, columnspan=10, sticky=N + W + S + E)
    for widget in okno.winfo_children():  # czycimy okno
        widget.destroy()
    #ustawiam okno
    for x in range(10):
        okno.columnconfigure(x,weight=0)
    for x in range(10):
        okno.rowconfigure(x,weight=0)
    okno.columnconfigure(0,weight=25)
    okno.columnconfigure(1, weight=6)
    okno.rowconfigure(0,weight=1)
    okno.rowconfigure(1, weight=1)

    Dzielimy_na = math.ceil(len(lanc_Kodonow) / 10) #licze na ile podobrazkow mam podzielić nasz łańcuch
    bialeczka = []#lista podzielonych łańcuchów
    poczatek = 0
    dokiedy = math.ceil(len(lanc_Kodonow) / Dzielimy_na)#dlugosc pojedynczego lancucha(na jednym obrazku)
    for x in range(Dzielimy_na): #dzielimy łancuch białek na części, odpowiedniej długości
        bialeczka.append(lanc_Kodonow[poczatek:dokiedy])
        poczatek = dokiedy
        dokiedy += math.floor(len(lanc_Kodonow) / Dzielimy_na)
    i = 0
    k = 0
    #tworzymy ramke gdzie wrzucimy obrazy
    obrazy=CTkScrollableFrame(okno,width=(int((size[0]/31)*25)),height=size[1])
    obrazy.grid(row=0,column=0,rowspan=2,sticky=W+E+S+N)
    obrazy.columnconfigure(0,weight=1)
    #standardowo rysujemy każde białko, a następnie wklejami do ramki
    for x in bialeczka:
        Smiles, ostatni = operacje_chemiczne.wzor_lancucha_aminokwasow(x)
        mol = Chem.MolFromSmiles(Smiles)
        rdCoordGen.AddCoords(mol)
        if (czy_podswietlaj):
            hit_bonds = []
            atomy = (0, ostatni)
            for x in range(ostatni):
                hit_bonds.append(x)
            obraz = Draw.MolToImage(mol, highlightBonds=hit_bonds, highlightAtoms=atomy, size=(int((size[0]/31)*25), 200))
        else:
            obraz = Draw.MolToImage(mol, size=(int((size[0]/31)*25), 200))
        wzor_strukturalny = CTkImage(light_image=obraz, dark_image=obraz, size=(int((size[0]/31)*25), 200))
        obraz = CTkLabel(obrazy, text="", image=wzor_strukturalny)
        obraz.grid(row=i, column=0)
        i += 1

    #wrzucamy pozostałe widgety
    wrocButton = CTkButton(okno, text="Wróć", width=(size[0] / 31) * 6,command=lambda: zniszcz(okno))
    wrocButton.grid(row=0, column=1, sticky=W + E)
    WykresyButton = CTkButton(okno, text="Wykresy", width=int((size[0] / 31) * 6),
                              command=lambda: wykresy.wykresy(okno, Smiles, lanc_Kodonow))

    WykresyButton.grid(row=1, column=1, sticky=W + E)
def zniszcz(okno):
    okno.destroy()
def rysuj_szybko_interface(ramex, lanc_Kodonow,czy_podswietlaj):
    okno = CTkFrame(ramex)
    size = (ramex.winfo_width(), ramex.winfo_height())
    okno.grid(row=0, column=0, rowspan=10, columnspan=10, sticky=N + W + S + E)
    for widget in okno.winfo_children():  # czycimy okno
        widget.destroy()
    # to jest lista bialek pelnych typu bialka = ('MIIIIIIII','MIIF')
    # bialka to lista jkbc
    # resetujemy formatowanie okna i tworzymy nowe
    for x in range(10):
        okno.columnconfigure(x, weight=0)
    for x in range(10):
        okno.rowconfigure(x, weight=0)
    okno.columnconfigure(0, weight=25)
    okno.columnconfigure(1, weight=6)
    okno.rowconfigure(0, weight=1)
    okno.rowconfigure(1, weight=1)
    Dzielimy_na = math.ceil(len(lanc_Kodonow) / 10)
    bialeczka = []
    poczatek = 0
    jakdlugo = math.ceil(len(lanc_Kodonow) / Dzielimy_na)
    if(Dzielimy_na<5): #jest to równoważne że długość jest mniejsza niż 50, jeśli jest mniejsza to dzielimy tak samo jak w funkcji rysuj_dziwnie_interface
        for x in range(Dzielimy_na):
            bialeczka.append(lanc_Kodonow[poczatek:jakdlugo])
            poczatek = jakdlugo
            jakdlugo += math.floor(len(lanc_Kodonow) / Dzielimy_na)
        i = 0
        k = 0
        obrazy=CTkScrollableFrame(okno,width=int((size[0]/31)*25),height=500)
        obrazy.grid(row=0,column=0,rowspan=2,sticky=N+W+E+S)
        for x in bialeczka:
            Smiles, ostatni = operacje_chemiczne.wzor_lancucha_aminokwasow(x)
            mol = Chem.MolFromSmiles(Smiles)
            rdCoordGen.AddCoords(mol)
            if (czy_podswietlaj):
                hit_bonds = []
                atomy = (0, ostatni)
                for x in range(ostatni):
                    hit_bonds.append(x)
                obraz = Draw.MolToImage(mol, highlightBonds=hit_bonds, highlightAtoms=atomy, size=(int((size[0]/31)*25), 200))
            else:
                obraz = Draw.MolToImage(mol, size=(int((size[0]/31)*25), 200))
            wzor_strukturalny = CTkImage(light_image=obraz, dark_image=obraz, size=(int((size[0]/31)*25), 200))
            obraz = CTkLabel(obrazy, text="", image=wzor_strukturalny)
            obraz.grid(row=i, column=0)
            i += 1
    else: #tworzymy 2 obrazki po 10 aminokwasów, wstawiamy 3 kropki i kolejne 2 po 10.
        bialeczka=[lanc_Kodonow[0:10],lanc_Kodonow[10:20],lanc_Kodonow[len(lanc_Kodonow)-21:len(lanc_Kodonow)-11],lanc_Kodonow[len(lanc_Kodonow)-11:len(lanc_Kodonow)-1]]#tworzymy podzieloną liste białek
        i = 0
        k = 0
        obrazy=CTkScrollableFrame(okno,width=int((size[0]/31)*25),height=500)
        obrazy.grid(row=0,column=0,rowspan=2,sticky=N+W+E+S)
        #robimy to samo co przy krótszym łancuchu, tylko z przerwą na "..." i na innej liście bialeczka
        for x in [bialeczka[0],bialeczka[1]]:
            Smiles, ostatni = operacje_chemiczne.wzor_lancucha_aminokwasow(x)
            mol = Chem.MolFromSmiles(Smiles)
            rdCoordGen.AddCoords(mol)
            if (czy_podswietlaj):
                hit_bonds = []
                atomy = (0, ostatni)
                for x in range(ostatni):
                    hit_bonds.append(x)
                obraz = Draw.MolToImage(mol, highlightBonds=hit_bonds, highlightAtoms=atomy, size=(int((size[0]/31)*25), 200))
            else:
                obraz = Draw.MolToImage(mol, size=(int((size[0]/31)*25), 200))
            wzor_strukturalny = CTkImage(light_image=obraz, dark_image=obraz, size=(int((size[0]/31)*25), 200))
            obraz = CTkLabel(obrazy, text="", image=wzor_strukturalny)
            obraz.grid(row=i, column=0)
            i += 1
        tekst=CTkLabel(obrazy,text="...", font=("Arial", 50))
        tekst.grid(row=i,column=0)
        i+=1
        for x in [bialeczka[len(bialeczka)-2],bialeczka[len(bialeczka)-1]]:
            Smiles, ostatni = operacje_chemiczne.wzor_lancucha_aminokwasow(x)
            mol = Chem.MolFromSmiles(Smiles)
            rdCoordGen.AddCoords(mol)
            if (czy_podswietlaj):
                hit_bonds = []
                atomy = (0, ostatni)
                for x in range(ostatni):
                    hit_bonds.append(x)
                obraz = Draw.MolToImage(mol, highlightBonds=hit_bonds, highlightAtoms=atomy, size=(int((size[0]/31)*25), 200))
            else:
                obraz = Draw.MolToImage(mol, size=(int((size[0]/31)*25), 200))
            wzor_strukturalny = CTkImage(light_image=obraz, dark_image=obraz, size=(int((size[0]/31)*25), 200))
            obraz = CTkLabel(obrazy, text="", image=wzor_strukturalny)
            obraz.grid(row=i, column=0)
            i += 1
    #reszta widgetów
    wrocButton = CTkButton(okno, text="Wróć", width=int((size[0]/31)*6), command=lambda: zniszcz(okno))
    WykresyButton = CTkButton(okno, text="Wykresy",width=int((size[0]/31)*6), command=lambda: wykresy.wykresy(okno, Smiles, lanc_Kodonow))

    WykresyButton.grid(row=1, column=1, sticky=W + E)
    wrocButton.grid(row=0, column=1, sticky=W + E)
def rysuj_interface(ramex, lanc_Kodonow,czy_podswietlaj): #funkcja rysująca cały łancuch polipeptydowy w jednym obrazku
    okno = CTkFrame(ramex)
    size = (ramex.winfo_width(), ramex.winfo_height())
    okno.grid(row=0, column=0, rowspan=10, columnspan=10, sticky=N + W + S + E)
    for widget in okno.winfo_children():  # czycimy okno
        widget.destroy()
    for widget in okno.winfo_children():  # czyscimy okno
        widget.destroy()
    for x in range(10):
        okno.columnconfigure(x, weight=0)
    for x in range(10):
        okno.rowconfigure(x, weight=0)
    okno.columnconfigure(0, weight=25)
    okno.columnconfigure(1, weight=6)
    okno.rowconfigure(0, weight=1)
    okno.rowconfigure(1, weight=1)
    Smiles, ostatni = operacje_chemiczne.wzor_lancucha_aminokwasow(lanc_Kodonow) #uzyskujemy wzór smiles
    mol = Chem.MolFromSmiles(Smiles) #tworzymy obiekt mol
    if(len(lanc_Kodonow)<40):
        rdCoordGen.AddCoords(mol) #uwzględniamy coordynaty 2D
    #rysujemy obrazek (z opcją podświetlania lub bez)
    if (czy_podswietlaj):
        hit_bonds = []
        atomy = (0, ostatni)
        for x in range(ostatni):
            hit_bonds.append(x)
        obraz = Draw.MolToImage(mol, highlightBonds=hit_bonds, highlightAtoms=atomy, size=(int(size[0]*10/13),int(size[1])))
    else:
        obraz = Draw.MolToImage(mol, size=(int(size[0]*10/13),int(size[1])))
    #wrzucamy widgety na okno, tworzymy przyciski
    wzor_strukturalny = CTkImage(light_image=obraz, dark_image=obraz, size=(int(size[0]*10/13),int(size[1])))
    WrocButton = CTkButton(okno, text="Wróć", command=lambda: zniszcz(okno),width=size[0]*3/13)
    WykresyButton = CTkButton(okno, text="Wykresy", command=lambda: wykresy.wykresy(okno, Smiles,lanc_Kodonow))
    obraz = CTkLabel(okno, image=wzor_strukturalny, text="")
    obraz.grid(row=0, column=0, rowspan=2)
    WrocButton.grid(row=0, column=1,sticky=E+W)
    WykresyButton.grid(row=1, column=1,sticky=E+W)