
from Bio.Seq import Seq
from customtkinter import *
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import rdCoordGen
import operacje_chemiczne
import rysowania
from tkinter import filedialog
import math
from Bio.SeqUtils.ProtParam import ProteinAnalysis
mode = "dark"  # zmienna przechowująca tryb aplikacji
base_color = "dark-blue"
czy_podswietlaj = True

class przycisk_dane():  # klasa, która tworzy przycisk ctk z zapamiętaniem danych. Używam do zapisu danych z fora w przycisku (rozwiązuje errora)
    def __init__(self,okno, oknofunkcja, zwiazek, lanc,czy):
        self.przycisk = CTkButton(okno, text="Generuj", command=lambda: rysowania.rysuj_interface(oknofunkcja, zwiazek, lanc,czy))

    def gridbutton(self, x, y):
        self.przycisk.grid(row=x, column=y,sticky=E+W)
class przycisk_szczegol_dane():  # ten sam przycisk co na gorze tylko z inna funkcja
    def __init__(self, okno, oknofunkcja, zwiazek, lanc,czy):
        self.przycisk = CTkButton(okno, text="G. Szczegółowo",
                                  command=lambda: rysowania.rysuj_dziwnie_interface(oknofunkcja, zwiazek, lanc,czy))

    def gridbutton(self, x, y):
        self.przycisk.grid(row=x, column=y,sticky=E+W)
class przycisk_szybkie_dane():  # ten sam przycisk co na gorze tylko z inna funkcja
    def __init__(self, okno, oknofunkcja, zwiazek, lanc,czy):
        self.przycisk = CTkButton(okno, text="G. Szybko",
                                  command=lambda: rysowania.rysuj_szybko_interface(oknofunkcja, zwiazek, lanc,czy))

    def gridbutton(self, x, y):
        self.przycisk.grid(row=x, column=y,sticky=E+W)

def interfejs(okno):  # funkcja obsługująca główny interfejs programu(menu widoczne po włączeniu
    for widget in okno.winfo_children():
        widget.destroy()
    global size #potrzebuje wielkosci, żeby widgety miały odpowiednie wymiary
    size=(okno.winfo_width(),okno.winfo_height())
    if(size==(200,200)): #przy pierwszym włączeniu jest bug funkcji winfo, wynik to zawsze (200,200)
        size=(650,500)
    #2 fory czyszczą ustawienie wagi kolumn i wierszy
    for x in range(10): #
        okno.columnconfigure(x, weight=0)
    okno.columnconfigure(0, weight=3)
    okno.columnconfigure(1, weight=10)
    for x in range(10):
        okno.rowconfigure(x, weight=1)
    #tworze odpowiednie przyciski
    napis = CTkLabel(okno, text="Wybierz Opcje")
    plikButton = CTkButton(okno, text="Dane z Pliku", command=lambda: wczytaj_z_pliku(okno),width=int(size[0]*3/13))
    recznieButton = CTkButton(okno, text="Dane Ręcznie", command=lambda: wczytaj_recznie(okno, lancpoczatkowy="AUGAAUGCAUGUAGAUAGAUAGAUGUGA"))
    instrukcjeButton = CTkButton(okno, text="Instrukcja", command=lambda: instrukcja(okno))
    opcjeButton = CTkButton(okno, text="Opcje", command=lambda: interfaceOpcje(okno))
    #rysuje przyciski
    plikButton.grid(row=1, column=0, sticky=W + E)
    recznieButton.grid(row=2, column=0, sticky=W + E)
    instrukcjeButton.grid(row=3, column=0, sticky=W + E)
    opcjeButton.grid(row=4, column=0, sticky=W + E)
    napis.grid(row=0, column=0, sticky=W + E)
    #tworze obraz który widać na starcie:
    #tworze wzór smiles tego związku
    seq = Seq("AUGGAACGCGAACCCUAC")
    translated = seq.translate()
    wzor, ostatni = operacje_chemiczne.wzor_lancucha_aminokwasow(translated)
    #tworze obiekt molecule który później rysuje
    mol = Chem.MolFromSmiles(wzor)
    rdCoordGen.AddCoords(mol)
    #podświetlam łańcuch główny
    hit_bonds = []
    for x in range(ostatni):
        hit_bonds.append(x)
    atomy = (0, ostatni)
    #ostatecznie go podświetlam
    if (czy_podswietlaj):
        obrazStartowy = Draw.MolToImage(mol, highlightBonds=hit_bonds, highlightAtoms=atomy, highlightColor=((0, 1, 0)),
                                        size=(int(size[0]*10/13),int(size[1]))) #wyrażenie dostosowywuje rozmiar obrazka do okna, niestety nie dynamicznie :(
    else:
        obrazStartowy = Draw.MolToImage(mol, size=(int(size[0]*10/13),int(size[1])))
    img1 = CTkImage(light_image=obrazStartowy, dark_image=obrazStartowy, size=(int(size[0]*10/13),int(size[1])))
    #wrzucam rysunek w okno
    obraz = CTkLabel(okno, text="", image=img1)
    obraz.grid(row=0, column=1, rowspan=10,sticky=NW)




def wczytaj_z_pliku(okno):  # funkcja wczytująca dane z pliku
    for widget in okno.winfo_children():
        widget.destroy()
    global size
    size=(okno.winfo_width(), okno.winfo_height())
    #czyszcze i ustawiam nowe wagi kolumn
    for x in range(10):
        okno.columnconfigure(x, weight=0)
    for x in range(10):
        okno.rowconfigure(x, weight=0)
    for x in range(1,4):
        okno.rowconfigure(x, weight=1)
    for x in range(6):
        okno.columnconfigure(x, weight=1)
    #tworze przyciski, i pole wejsciowe:
    napisTytul = CTkLabel(okno, text="wpisz lokalizacje pliku:")
    sciezkaWejscie = CTkEntry(okno, width=(size[0]/6)*5)
    przyciskPrzegladaj = CTkButton(okno, text="przegladaj",width=(size[0]/6), command=lambda: przegladaj(sciezkaWejscie)) #otwiera możliwość wyboru pliku z windowsowego explorera
    przyciskSzukaj = CTkButton(okno, text="Szukaj", command=lambda: szukaj_z_pliku())
    przyciskWroc = CTkButton(okno, text="Wróc", command=lambda: interfejs(okno)) #wraca do menu głównego
    #wrzucam widgety na okno
    napisTytul.grid(row=0, column=0,sticky=W+E)
    sciezkaWejscie.grid(row=1, column=0, columnspan=5,sticky=W+E+N)
    przyciskPrzegladaj.grid(row=1, column=5,sticky=W+E+N)
    przyciskSzukaj.grid(row=2, column=5,sticky=W+E+N)
    przyciskWroc.grid(row=3, column=5,sticky=W+E+N)


def przegladaj(wejscie):
    filename = filedialog.askopenfilename(initialdir=os.getcwd(), title="Select a File",
                                          filetypes=(("Text files", "*.txt*"), ("all files", "*.*"))) #funkcja otwierająca wyszukiwarke plików
    wejscie.insert(0, filename)


def szukaj_z_pliku():
    print("hendorzyc disa psa syna diabla")


def wczytaj_recznie(okno, lancpoczatkowy=""):  # funkcja wczytująca dane z "palca"
    #czyszcze okno
    for widget in okno.winfo_children():
        widget.destroy()
    global size
    size = (okno.winfo_width(), okno.winfo_height())
    #ustawiam wagi kolumn
    for x in range(10):
        okno.columnconfigure(x, weight=0)
    for x in range(10):
        okno.rowconfigure(x, weight=0)
    for x in range(5):
        okno.columnconfigure(x, weight=2)
    okno.columnconfigure(5, weight=3)
    okno.rowconfigure(1,weight=1)

    #tworze elementy okna
    wejscie = CTkEntry(okno, width=((size[0]*8)/15))
    wejscie.insert(0, lancpoczatkowy)
    szukaj_interface(okno, wejscie)
    opis = CTkLabel(okno, text="Wpisz kod nici Rna", width=((size[0]*2)/15))
    szukajButton = CTkButton(okno, text="Szukaj Białek", width=((size[0]*3)/15), command=lambda: szukaj_interface(okno, wejscie)) #wyszukuje bialka w wpisanym lancuchu, i rysuje je w ramce
    wrocButton = CTkButton(okno, text="Wroc", width=15, command=lambda: interfejs(okno))

    opis.grid(row=0, column=0,sticky=W+E)
    wejscie.grid(row=0, column=1, columnspan=4,sticky=W+E)
    szukajButton.grid(row=0, column=5,sticky=W+E)
    wrocButton.grid(row=1, column=5,sticky=W+E)


def szukaj_interface(okno, wejscie):  # funkcja wypisująca wszystkie łańcuchy w interfejsie
    global size
    global czy_podswietlaj
    size = (okno.winfo_width(), okno.winfo_height())
    Bwidth = (size[0] / 15) * 12  # zmienna przechowująca szerokość naszej ramki

    lanc=Seq(wejscie.get()) #czysty lancuch wpisany
    lanc1,lanc2,lanc3 = operacje_chemiczne.translacjaBezBugow(wejscie.get())#obrobione lancuchy(podzielne przez 3, zdebugowane, rodzielone na 3 podlinie)
    bialka1 = operacje_chemiczne.rozklad_na_bialka(lanc1) #3 łancuchy białek dla 3 przesunięć
    bialka2 = operacje_chemiczne.rozklad_na_bialka(lanc2)
    bialka3 =   operacje_chemiczne.rozklad_na_bialka(lanc3)
    #tworze listy w których przechowywane będą WIDGETY
    bialka_napisy = []
    bialka_przyciski = []
    bialka_szczegol_przyciski = []
    bialka_szybkie_przyciski=[]
    #tworze ramke w której będą wyświetlane nasze białeczka
    BialkaRamka=CTkScrollableFrame(okno,width=(size[0]/15)*12,height=(size[1]-50))
    BialkaRamka.grid(row=1,column=0,columnspan=4,sticky=E+W+N+S)
    for x in range(5):
        BialkaRamka.columnconfigure(x, weight=1)
    #poniższe 3 fory tworzą przyciski i napisy do białek
    for x in bialka1:
        if(len(x)>10): #jeśli dłuższe niż 10, to go skracamy
            bialka_napisy.append(CTkLabel(BialkaRamka, text=(x[:5]+"..."+x[len(x)-5:]+" P:0"),width=(Bwidth/5)*2))    #napisy do białek
        else:
            bialka_napisy.append(CTkLabel(BialkaRamka, text=(x +" P:0"),width=(Bwidth/5)*2))  # napisy do białek
        bialka_przyciski.append(przycisk_dane(BialkaRamka,okno, x, lanc,czy_podswietlaj)) #przyciski do normalnej generacji
        bialka_szczegol_przyciski.append(przycisk_szczegol_dane(BialkaRamka,okno, x, lanc,czy_podswietlaj)) #przyciski do "szczegółowej" generacji
        bialka_szybkie_przyciski.append(przycisk_szybkie_dane(BialkaRamka,okno, x, lanc,czy_podswietlaj))#przyciski do szybkiej generacji
    for x in bialka2:
        if(len(x) > 10):  # jeśli dłuższe niż 10, to go skracamy
            bialka_napisy.append(CTkLabel(BialkaRamka, text=(x[:5] + "..." + x[len(x) - 5:] + " P:1"),width=(Bwidth/5)*2))  # napisy do białek
        else:
            bialka_napisy.append(CTkLabel(BialkaRamka, text=(x + " P:1"),width=(Bwidth/5)*2))  # napisy do białek
        bialka_przyciski.append(przycisk_dane(BialkaRamka, okno, x, lanc,czy_podswietlaj))  # przyciski do normalnej generacji
        bialka_szczegol_przyciski.append(
        przycisk_szczegol_dane(BialkaRamka, okno, x, lanc,czy_podswietlaj))  # przyciski do "szczegolj" generacji
        bialka_szybkie_przyciski.append(przycisk_szybkie_dane(BialkaRamka,okno, x, lanc,czy_podswietlaj))# przyciski do szybkiej generacji
    for x in bialka3:
        if(len(x) > 10):  # jeśli dłuższe niż 10, to go skracamy
            bialka_napisy.append(CTkLabel(BialkaRamka, text=(x[:5] + "..." + x[len(x) - 5:] + " P:2"),width=(Bwidth/5)*2))  # napisy do białek
        else:
            bialka_napisy.append(CTkLabel(BialkaRamka, text=(x + " P:2"),width=(Bwidth/5)*2))  # napisy do białek
        bialka_przyciski.append(przycisk_dane(BialkaRamka, okno, x, lanc,czy_podswietlaj))  # przyciski do normalnej generacji
        bialka_szczegol_przyciski.append(
        przycisk_szczegol_dane(BialkaRamka, okno, x, lanc,czy_podswietlaj))  # przyciski do "szczegolj" generacji
        bialka_szybkie_przyciski.append(przycisk_szybkie_dane(BialkaRamka,okno, x, lanc,czy_podswietlaj))# przyciski do szybkiej generacji
    i = 0
    #poniższe fory rysują nasze przyciski i napisy w ramce
    for x in bialka_napisy:
        x.grid(row=i, column=0, columnspan=2,sticky=E+W)
        i += 1
    i = 0
    for x in bialka_szczegol_przyciski:
        x.gridbutton(i, 2) #funkcja klasy, działa jak zwykły grid
        i += 1
    i = 0
    for x in bialka_przyciski:
        x.gridbutton(i, 3) #funkcja klasy, działa jak zwykły grid
        i += 1
    i = 0
    for x in bialka_szybkie_przyciski:
        x.gridbutton(i, 4) #funkcja klasy, działa jak zwykły grid
        i += 1

def interfaceOpcje(okno):  # funkcja wczytująca interface opcji
    for widget in okno.winfo_children():
        widget.destroy()
    #ustawiamy wagi
    for x in range(10):
        okno.columnconfigure(x, weight=0)
    for x in range(10):
        okno.rowconfigure(x, weight=0)
    okno.columnconfigure(0,weight=1)
    #tworze widgety
    napisMode = CTkLabel(okno, text="zmien kolor tla:")
    napisKolor = CTkLabel(okno, text="zmien kolor przycisków:")
    napisWroc = CTkLabel(okno, text="wróc do menu:")
    napisPodswietlaj = CTkLabel(okno, text="Podświetlać główny łańcuch:")
    przyciskMode = CTkButton(okno, text=mode, command=lambda: zmien_tryb(okno))#zmienia tryb jasny na ciemny i viceversa
    przyciskKolor = CTkButton(okno, text=base_color, command=lambda: zmien_kolor(okno))#zmienia bazowy kolor
    przyciskWroc = CTkButton(okno, text="Wróc", command=lambda: zapis(okno))
    if (czy_podswietlaj == True): #zmienia czy podswietla się główny łańcuch
        przyciskPodswietlaj = CTkButton(okno, text="Tak", command=lambda: podswietlenie(okno))
    else:
        przyciskPodswietlaj = CTkButton(okno, text="Nie", command=lambda: podswietlenie(okno))
    #wrzucam widgety w okno
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