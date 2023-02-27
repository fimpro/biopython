
from Bio.Seq import Seq
import time
import sys
import os
from customtkinter import *
from PIL import ImageTk, Image
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from rdkit.Chem import rdCoordGen
import operacje_chemiczne
from tkinter import filedialog
import tkinter
import math
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import time
mode = "dark"  # zmienna przechowująca tryb aplikacji
base_color = "dark-blue"
czy_podswietlaj = True

class przycisk_dane():  # klasa, która tworzy przycisk ctk z zapamiętaniem danych. Używam do zapisu danych z fora w przycisku (rozwiązuje errora)
    def __init__(self,okno, oknofunkcja, zwiazek, lanc):
        self.przycisk = CTkButton(okno, text="Generuj", command=lambda: rysuj_interface(oknofunkcja, zwiazek, lanc))

    def gridbutton(self, x, y):
        self.przycisk.grid(row=x, column=y,sticky=E+W)


class przycisk_dziwne_dane():  # ten sam przycisk co na gorze tylko z inna funkcja
    def __init__(self, okno, oknofunkcja, zwiazek, lanc):
        self.przycisk = CTkButton(okno, text="Generuj_dziwnie",
                                  command=lambda: rysuj_dziwnie_interface(oknofunkcja, zwiazek, lanc))

    def gridbutton(self, x, y):
        self.przycisk.grid(row=x, column=y,sticky=E+W)


def interfejs(okno):  # funkcja obsługująca główny interfejs programu
    for widget in okno.winfo_children():
        widget.destroy()
    global size
    size=(okno.winfo_width(),okno.winfo_height())
    if(size==(200,200)):
        size=(650,500)
    okno.columnconfigure(0, weight=3)
    okno.columnconfigure(1, weight=10)
    for x in range(10):
        okno.rowconfigure(x, weight=1)
    napis = CTkLabel(okno, text="Wybierz Opcje")
    plikButton = CTkButton(okno, text="Dane z Pliku", command=lambda: wczytaj_z_pliku(okno),width=size[0]*3/13)
    recznieButton = CTkButton(okno, text="Dane Ręcznie", command=lambda: wczytaj_recznie(okno, lancpoczatkowy="AUGUAA"))
    instrukcjeButton = CTkButton(okno, text="Instrukcja", command=lambda: instrukcja(okno))
    opcjeButton = CTkButton(okno, text="Opcje", command=lambda: interfaceOpcje(okno))
    plikButton.grid(row=1, column=0, sticky=W + E)
    recznieButton.grid(row=2, column=0, sticky=W + E)
    instrukcjeButton.grid(row=3, column=0, sticky=W + E)
    opcjeButton.grid(row=4, column=0, sticky=W + E)
    napis.grid(row=0, column=0, sticky=W + E)

    seq = Seq("AUGGAACGCGAACCCUAC")
    translated = seq.translate()
    wzor, ostatni = operacje_chemiczne.wzor_lancucha_aminokwasow(translated)
    mol = Chem.MolFromSmiles(wzor)
    rdCoordGen.AddCoords(mol)
    hit_bonds = []
    for x in range(ostatni):
        hit_bonds.append(x)
    atomy = (0, ostatni)
    if (czy_podswietlaj):
        obrazStartowy = Draw.MolToImage(mol, highlightBonds=hit_bonds, highlightAtoms=atomy, highlightColor=((0, 1, 0)),
                                        size=(int(size[0]*10/13),int(size[1]))) #wyrażenie dostosowywuje rozmiar obrazka do okna, niestety nie dynamicznie :(
    else:
        obrazStartowy = Draw.MolToImage(mol, size=(int(size[0]*10/13),int(size[1])))
    img1 = CTkImage(light_image=obrazStartowy, dark_image=obrazStartowy, size=(int(size[0]*10/13),int(size[1])))
    obraz = CTkLabel(okno, text="", image=img1)
    obraz.grid(row=0, column=1, rowspan=10,sticky=NW)




def wczytaj_z_pliku(okno):  # funkcja wczytująca dane z pliku
    for widget in okno.winfo_children():
        widget.destroy()

    napisTytul = CTkLabel(okno, text="wpisz lokalizacje pliku:")
    sciezkaWejscie = CTkEntry(okno, width=500)
    przyciskPrzegladaj = CTkButton(okno, text="przegladaj", command=lambda: przegladaj(sciezkaWejscie))
    przyciskSzukaj = CTkButton(okno, text="Szukaj", command=lambda: szukaj_z_pliku(okno))
    przyciskWroc = CTkButton(okno, text="Wróc", command=lambda: interfejs(okno))
    napisTytul.grid(row=0, column=0)
    sciezkaWejscie.grid(row=1, column=0, columnspan=5)
    przyciskPrzegladaj.grid(row=1, column=6)
    przyciskSzukaj.grid(row=2, column=6)
    przyciskWroc.grid(row=3, column=6)


def przegladaj(wejscie):
    filename = filedialog.askopenfilename(initialdir=os.getcwd(), title="Select a File",
                                          filetypes=(("Text files", "*.txt*"), ("all files", "*.*")))
    wejscie.insert(0, filename)


def szukaj_z_pliku():
    print("hendorzyc disa psa syna diabla")


def wczytaj_recznie(okno, lancpoczatkowy=""):  # funkcja wczytująca dane z "palca"
    for widget in okno.winfo_children():
        widget.destroy()
    global size
    size = (okno.winfo_width(), okno.winfo_height())
    for x in range(5):
        okno.columnconfigure(x, weight=2)
    okno.columnconfigure(5, weight=3)
    genom = Seq(lancpoczatkowy)
    wejscie = CTkEntry(okno, width=((size[0]*8)/15))
    wejscie.insert(0, lancpoczatkowy)
    szukaj_interface(okno, wejscie)
    opis = CTkLabel(okno, text="Wpisz kod nici Rna", width=((size[0]*2)/15))
    szukajButton = CTkButton(okno, text="Szukaj Białek", width=((size[0]*3)/15), command=lambda: szukaj_interface(okno, wejscie))
    wrocButton = CTkButton(okno, text="Wroc", width=15, command=lambda: interfejs(okno))

    opis.grid(row=0, column=0,sticky=W+E)
    wejscie.grid(row=0, column=1, columnspan=4,sticky=W+E)
    szukajButton.grid(row=0, column=5,sticky=W+E)
    wrocButton.grid(row=1, column=5,sticky=W+E)


def szukaj_interface(okno, wejscie):  # funkcja wypisująca wszystkie łańcuchy w interfejsie
    global size
    size = (okno.winfo_width(), okno.winfo_height())
    for x in range(10):
        okno.rowconfigure(x,weight=0)
    okno.rowconfigure(1, weight=1)
    lanc=Seq(wejscie.get()) #czysty lancuch wpisany
    lanc1,lanc2,lanc3 = operacje_chemiczne.translacjaBezBugow(wejscie.get())#obrobione lancuchy
    bialka1 = operacje_chemiczne.rozklad_na_bialka(lanc1) #3 łancuchy białek dla 3 przesunięć
    bialka2 = operacje_chemiczne.rozklad_na_bialka(lanc2)
    bialka3=operacje_chemiczne.rozklad_na_bialka(lanc3)
    bialka_napisy = []
    bialka_przyciski = []
    bialka_dziwne_przyciski = []
    # to jest lista bialek pelnych typu bialka = ('MIIIIIIII','MIIF') a x oznacza ktore z tych
    BialkaRamka=CTkScrollableFrame(okno,width=(size[0]/15)*12,height=(size[1]-50))
    Bwidth=(size[0]/15)*12#zmienna przechowująca szerokość naszej ramki
    BialkaRamka.grid(row=1,column=0,columnspan=4,sticky=E+W+N+S)
    for x in range(5):
        BialkaRamka.columnconfigure(x, weight=1)
    for x in bialka1:
        if(len(x)>10): #jeśli dłuższe niż 10, to go skracamy
            bialka_napisy.append(CTkLabel(BialkaRamka, text=(x[:5]+"..."+x[len(x)-5:]+" P:0"),width=(Bwidth/5)*3))    #napisy do białek
        else:
            bialka_napisy.append(CTkLabel(BialkaRamka, text=(x +" P:0"),width=(Bwidth/5)*3))  # napisy do białek
        bialka_przyciski.append(przycisk_dane(BialkaRamka,okno, x, lanc)) #przyciski do normalnej generacji
        bialka_dziwne_przyciski.append(przycisk_dziwne_dane(BialkaRamka,okno, x, lanc)) #przyciski do "dziwnej" generacji
        #przyciski do szybkiej generacji
    for x in bialka2:
        if(len(x) > 10):  # jeśli dłuższe niż 10, to go skracamy
            bialka_napisy.append(CTkLabel(BialkaRamka, text=(x[:5] + "..." + x[len(x) - 5:] + " P:1"),width=(Bwidth/5)*3))  # napisy do białek
        else:
            bialka_napisy.append(CTkLabel(BialkaRamka, text=(x + " P:1"),width=(Bwidth/5)*3))  # napisy do białek
        bialka_przyciski.append(przycisk_dane(BialkaRamka, okno, x, lanc))  # przyciski do normalnej generacji
        bialka_dziwne_przyciski.append(
        przycisk_dziwne_dane(BialkaRamka, okno, x, lanc))  # przyciski do "dziwnej" generacji
        # przyciski do szybkiej generacji
    for x in bialka3:
        if(len(x) > 10):  # jeśli dłuższe niż 10, to go skracamy
            bialka_napisy.append(CTkLabel(BialkaRamka, text=(x[:5] + "..." + x[len(x) - 5:] + " P:2"),width=(Bwidth/5)*3))  # napisy do białek
        else:
            bialka_napisy.append(CTkLabel(BialkaRamka, text=(x + " P:2"),width=(Bwidth/5)*3))  # napisy do białek
        bialka_przyciski.append(przycisk_dane(BialkaRamka, okno, x, lanc))  # przyciski do normalnej generacji
        bialka_dziwne_przyciski.append(
        przycisk_dziwne_dane(BialkaRamka, okno, x, lanc))  # przyciski do "dziwnej" generacji
        # przyciski do szybkiej generacji
    i = 2
    for x in bialka_napisy:
        x.grid(row=i, column=0, columnspan=3,sticky=E+W)
        i += 1
    i = 2
    for x in bialka_dziwne_przyciski:
        x.gridbutton(i, 3) #funkcja klasy, działa jak zwykły grid
        i += 1
    i = 2
    for x in bialka_przyciski:
        x.gridbutton(i, 4) #funkcja klasy, działa jak zwykły grid
        i += 1


def rysuj_dziwnie_interface(okno, lanc_Kodonow, lancpowrotny):
    for widget in okno.winfo_children():  # czycimy okno
        widget.destroy()
    # to jest lista bialek pelnych typu bialka = ('MIIIIIIII','MIIF')
    # bialka to lista jkbc
    Dzielimy_na = math.ceil(len(lanc_Kodonow) / 10)
    bialeczka = []
    poczatek = 0
    jakdlugo = math.ceil(len(lanc_Kodonow) / Dzielimy_na)
    for x in range(Dzielimy_na):
        bialeczka.append(lanc_Kodonow[poczatek:jakdlugo])
        poczatek = jakdlugo
        jakdlugo += math.floor(len(lanc_Kodonow) / Dzielimy_na)
    i = 0
    k = 0
    obrazy=CTkScrollableFrame(okno,width=500,height=500)
    obrazy.grid(row=0,column=0)
    for x in bialeczka:
        Smiles, ostatni = operacje_chemiczne.wzor_lancucha_aminokwasow(x)
        mol = Chem.MolFromSmiles(Smiles)
        rdCoordGen.AddCoords(mol)
        if (czy_podswietlaj):
            hit_bonds = []
            atomy = (0, ostatni)
            for x in range(ostatni):
                hit_bonds.append(x)
            obraz = Draw.MolToImage(mol, highlightBonds=hit_bonds, highlightAtoms=atomy, size=(500, 200))
        else:
            obraz = Draw.MolToImage(mol, size=(500, 200))
        wzor_strukturalny = CTkImage(light_image=obraz, dark_image=obraz, size=(500, 200))
        obraz = CTkLabel(obrazy, text="", image=wzor_strukturalny)
        obraz.grid(row=i, column=0)
        i += 1


    wrocButton = CTkButton(okno, text="Wroc", width=120, command=lambda: wczytaj_recznie(okno,lancpowrotny))
    wrocButton.grid(row=0, column=1)
    okno.mainloop()


def rysuj_interface(okno, lanc_Kodonow, lancpowrotny):
    global size
    size = (okno.winfo_width(), okno.winfo_height())
    for widget in okno.winfo_children():  # czyscimy okno
        widget.destroy()
    okno.columnconfigure(0, weight=10)
    okno.columnconfigure(1, weight=3)

    Smiles, ostatni = operacje_chemiczne.wzor_lancucha_aminokwasow(lanc_Kodonow)
    mol = Chem.MolFromSmiles(Smiles)
    rdCoordGen.AddCoords(mol)
    if (czy_podswietlaj):
        hit_bonds = []
        atomy = (0, ostatni)
        for x in range(ostatni):
            hit_bonds.append(x)
        obraz = Draw.MolToImage(mol, highlightBonds=hit_bonds, highlightAtoms=atomy, size=(int(size[0]*10/13),int(size[1])))
    else:
        obraz = Draw.MolToImage(mol, size=(int(size[0]*10/13),int(size[1])))
    wzor_strukturalny = CTkImage(light_image=obraz, dark_image=obraz, size=(int(size[0]*10/13),int(size[1])))
    WrocButton = CTkButton(okno, text="wroc", command=lambda: wczytaj_recznie(okno, lancpowrotny),width=size[0]*3/13)
    WykresyButton = CTkButton(okno, text="wykresy", command=lambda: wykresy(okno, Smiles, lancpowrotny,lanc_Kodonow))
    obraz = CTkLabel(okno, image=wzor_strukturalny, text="")
    obraz.grid(row=0, column=0, rowspan=6)
    WrocButton.grid(row=0, column=1,sticky=E+W)
    WykresyButton.grid(row=1, column=1,sticky=E+W)
    okno.mainloop()


def wykresy(okno, wzor, lancpowrotny, lanc_Kodonow):  # robocza funkcja do wykresów

    for widget in okno.winfo_children():
        widget.destroy()

    analizuj = ProteinAnalysis(lanc_Kodonow)  # tablica do analizowania bialek
    dw_y = []  # dane wykres w osi y; tablica 2d
    dw_y.append(analizuj.instability_index())
    dw_y.append(analizuj.isoelectric_point())
    dw_y.append(analizuj.get_amino_acids_percent())
    indeksh = operacje_chemiczne.indeks_hydrofobowy(lanc_Kodonow)
    dw_y.append(indeksh)
    print(lancpowrotny)
    print(dw_y)

    napis = CTkLabel(okno, text=wzor, width=100)
    returnButton = CTkButton(okno, text="Wróć", command=lambda: rysuj_interface(okno, lanc_Kodonow, lancpowrotny))
    AtomyButton = CTkButton(okno, text="Ilość Atomów", command=lambda: WypiszAtomy(okno,wzor, lanc_Kodonow, lancpowrotny))
    napis.grid(row=0, column=0)
    returnButton.grid(row=1, column=0)
    AtomyButton.grid(row=2, column=0)

def WypiszAtomy(okno, wzor,lancpowrotny, lanc_Kodonow):
    for widget in okno.winfo_children():
        widget.destroy()
    operacje_chemiczne.listaAtomow(wzor)
    returnButton = CTkButton(okno, text="Wróć", command=lambda: wykresy(okno,wzor, lanc_Kodonow, lancpowrotny))
    returnButton.grid(row=1, column=0)
def interfaceOpcje(okno):  # funkcja wczytująca interface opcji
    for widget in okno.winfo_children():
        widget.destroy()
    napisMode = CTkLabel(okno, text="zmien kolor tla:")
    napisKolor = CTkLabel(okno, text="zmien kolor przycisków:")
    napisWroc = CTkLabel(okno, text="wróc do menu:")
    napisPodswietlaj = CTkLabel(okno, text="Podświetlać główny łańcuch:")
    przyciskMode = CTkButton(okno, text=mode, command=lambda: zmien_tryb(okno))
    przyciskKolor = CTkButton(okno, text=base_color, command=lambda: zmien_kolor(okno))
    przyciskWroc = CTkButton(okno, text="Wróc", command=lambda: zapis(okno))
    if (czy_podswietlaj == True):
        przyciskPodswietlaj = CTkButton(okno, text="Tak", command=lambda: podswietlenie(okno))
    else:
        przyciskPodswietlaj = CTkButton(okno, text="Nie", command=lambda: podswietlenie(okno))
    przyciskMode.grid(row=1, column=0)
    przyciskKolor.grid(row=3, column=0)
    przyciskPodswietlaj.grid(row=5, column=0)
    przyciskWroc.grid(row=7, column=0)
    napisMode.grid(row=0, column=0)
    napisKolor.grid(row=2, column=0)
    napisPodswietlaj.grid(row=4, column=0)
    napisWroc.grid(row=6, column=0)


def podswietlenie(okno):
    global czy_podswietlaj
    if (czy_podswietlaj == True):
        czy_podswietlaj = False
        interfaceOpcje(okno)  # ta funkcja aktualizuje napisy na przyciskach
        return 0
    if (czy_podswietlaj == False):
        czy_podswietlaj = True
        interfaceOpcje(okno)  # ta funkcja aktualizuje napisy na przyciskach
        return 0


def zapis(okno):
    plik = open("opcje.txt", 'w')
    plik.write(mode + '\n')
    plik.write(base_color + '\n')
    plik.write(str(czy_podswietlaj) + '\n')
    plik.close()
    interfejs(okno)


def zmien_tryb(okno):
    global mode

    if (mode == "dark"):
        set_appearance_mode("light")
        mode = "light"
        interfaceOpcje(okno)  # ta funkcja aktualizuje napisy na przyciskach
        return 0
    if (mode == "light"):
        set_appearance_mode("dark")
        mode = "dark"
        interfaceOpcje(okno)  # ta funkcja aktualizuje napisy na przyciskach
        return 0


def zmien_kolor(okno):
    global base_color
    if (base_color == "dark-blue"):
        set_default_color_theme("blue")
        base_color = "blue"
        print("XD")
        interfaceOpcje(okno)  # ta funkcja aktualizuje napisy na przyciskach
        return 0
    if (base_color == "blue"):
        set_default_color_theme("green")
        base_color = "green"
        interfaceOpcje(okno)  # ta funkcja aktualizuje napisy na przyciskach
        return 0
    if (base_color == "green"):
        set_default_color_theme("dark-blue")
        base_color = "dark-blue"
        interfaceOpcje(okno)  # ta funkcja aktualizuje napisy na przyciskach


def instrukcja(okno):  # funkcja wczytująca instruckje użytkowania
   print("HWDP")

if __name__ == '__main__':
    oknoR = CTk()
    okno = CTkFrame(oknoR)
    oknoR.geometry("650x500")
    oknoR.columnconfigure(0,weight=1)
    oknoR.rowconfigure(0, weight=1)
    okno.grid(row=0,column=0,sticky=W+E+N+S)

    oknoR.minsize(height=500, width=650)
    oknoR.title("Katolgenom")
    plik = open("opcje.txt")
    opcje = plik.readlines()
    plik.close()
    mode = opcje[0].rstrip()
    base_color = opcje[1].rstrip()
    if (opcje[2] == "True\n"):
        czy_podswietlaj = True
    else:
        czy_podswietlaj = False
    set_appearance_mode(mode)
    set_default_color_theme(base_color)
    interfejs(okno)
    oknoR.mainloop()