from Bio import SeqIO
from Bio.Seq import Seq
from customtkinter import *
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import rdCoordGen
import operacje_chemiczne
import rysowania
from tkinter import filedialog
from PIL import Image
from tkinter import ttk
from tkinter.messagebox import showinfo
import math
import time
mode = "dark"  # zmienna przechowująca tryb aplikacji
base_color = "dark-blue"
czy_podswietlaj = True


def interfejs(okno):  # funkcja obsługująca główny interfejs programu(menu widoczne po włączeniu
    for widget in okno.winfo_children():
        widget.destroy()
    global size  # potrzebuje wielkosci, żeby widgety miały odpowiednie wymiary
    size = (okno.winfo_width(), okno.winfo_height())
    if (size[0] < 650 or size[1] < 500):  # przy pierwszym włączeniu jest bug funkcji winfo, wynik to zawsze (200,200)
        size = (750, 500)
    # 2 fory czyszczą ustawienie wagi kolumn i wierszy
    for x in range(10):  #
        okno.columnconfigure(x, weight=0)
    okno.columnconfigure(0, weight=3)
    okno.columnconfigure(1, weight=10)
    for x in range(10):
        okno.rowconfigure(x, weight=1)
    # tworze odpowiednie przyciski
    napis = CTkLabel(okno, text="Wybierz Opcje")
    plikButton = CTkButton(okno, text="Dane z Pliku", command=lambda: wczytaj_z_pliku(okno),
                           width=int(size[0] * 3 / 13))
    recznieButton = CTkButton(okno, text="Dane Ręcznie",
                              command=lambda: wczytaj_recznie(okno, lancpoczatkowy="AUGAAUGCAUGUAGAUAGAUAGAUGUGA"))
    instrukcjeButton = CTkButton(okno, text="Instrukcja", command=lambda: instrukcja(okno))
    opcjeButton = CTkButton(okno, text="Opcje", command=lambda: interfaceOpcje(okno))
    # rysuje przyciski
    plikButton.grid(row=1, column=0, sticky=W + E)
    recznieButton.grid(row=2, column=0, sticky=W + E)
    instrukcjeButton.grid(row=3, column=0, sticky=W + E)
    opcjeButton.grid(row=4, column=0, sticky=W + E)
    napis.grid(row=0, column=0, sticky=W + E)
    # tworze obraz który widać na starcie:
    # tworze wzór smiles tego związku
    seq = Seq("AUGGAACGCGAACCCUAC")
    translated = seq.translate()
    wzor, ostatni = operacje_chemiczne.wzor_lancucha_aminokwasow(translated)
    # tworze obiekt molecule który później rysuje
    mol = Chem.MolFromSmiles(wzor)
    rdCoordGen.AddCoords(mol)
    # podświetlam łańcuch główny
    hit_bonds = []
    for x in range(ostatni):
        hit_bonds.append(x)
    atomy = (0, ostatni)
    # ostatecznie go podświetlam
    if (czy_podswietlaj):
        obrazStartowy = Draw.MolToImage(mol, highlightBonds=hit_bonds, highlightAtoms=atomy, highlightColor=((0, 1, 0)),
                                        size=(int(size[0] * 10 / 13), int(size[
                                                                              1])))  # wyrażenie dostosowywuje rozmiar obrazka do okna, niestety nie dynamicznie :(
    else:
        obrazStartowy = Draw.MolToImage(mol, size=(int(size[0] * 10 / 13), int(size[1])))
    img1 = CTkImage(light_image=obrazStartowy, dark_image=obrazStartowy, size=(int(size[0] * 10 / 13), int(size[1])))
    # wrzucam rysunek w okno
    obraz = CTkLabel(okno, text="", image=img1)
    obraz.grid(row=0, column=1, rowspan=10, sticky=NW)


def wczytaj_recznie(okno, lancpoczatkowy=""):  # funkcja wczytująca dane z "palca"
    file_path = "x"
    # czyszcze okno
    for widget in okno.winfo_children():
        widget.destroy()
    global size
    size = (okno.winfo_width(), okno.winfo_height())
    # ustawiam wagi kolumn
    for x in range(10):
        okno.columnconfigure(x, weight=0)
    for x in range(10):
        okno.rowconfigure(x, weight=0)
    for x in range(4):
        okno.columnconfigure(x, weight=2)
    okno.columnconfigure(4, weight=3)
    okno.rowconfigure(0, weight=1)
    okno.rowconfigure(1, weight=1)
    okno.rowconfigure(2, weight=5)

    # tworze elementy okna
    wejscie = CTkEntry(okno, width=((size[0] * 8) / 15))
    wejscie.insert(0, lancpoczatkowy)

    opis = CTkLabel(okno, text="Wpisz kod nici RNA lub DNA:", width=((size[0] * 2) / 15))
    szukajButton = CTkButton(okno, text="Szukaj Białek", width=((size[0] * 3) / 15),
                             command=lambda: otwieranie_lancucha(okno, file_path,
                                                                wejscie.get()))  # wyszukuje bialka w wpisanym lancuchu, i rysuje je w ramce
    wrocButton = CTkButton(okno, text="Wróć", width=15, command=lambda: interfejs(okno))

    wzor, ostatni = operacje_chemiczne.wzor_lancucha_aminokwasow("MPPPPPPPNACR")
    # tworze obiekt molecule który później rysuje
    mol = Chem.MolFromSmiles(wzor)
    rdCoordGen.AddCoords(mol)
    # podświetlam łańcuch główny
    hit_bonds = []
    for x in range(ostatni):
        hit_bonds.append(x)
    atomy = (0, ostatni)
    # ostatecznie go podświetlam
    if (czy_podswietlaj):
        obrazStartowy = Draw.MolToImage(mol, highlightBonds=hit_bonds, highlightAtoms=atomy, highlightColor=((0, 1, 0)),
                                        size=(int(size[0] * 10 / 13), int(size[
                                                                              1]) - 75))  # wyrażenie dostosowywuje rozmiar obrazka do okna, niestety nie dynamicznie :(
    else:
        obrazStartowy = Draw.MolToImage(mol, size=(int(size[0] * 10 / 13), int(size[1]) - 75))
    img1 = CTkImage(light_image=obrazStartowy, dark_image=obrazStartowy,
                    size=(int(size[0] * 10 / 13), int(size[1]) - 75))
    # wrzucam rysunek w okno
    obraz = CTkLabel(okno, text="", image=img1)
    obraz.grid(row=2, column=0, columnspan=4, sticky=N)
    opis.grid(row=0, column=0, sticky=W + E)
    wejscie.grid(row=1, column=0, columnspan=4, sticky=W + E)
    szukajButton.grid(row=1, column=4, sticky=W + E)
    wrocButton.grid(row=2, column=4, sticky=W + E)


def wczytaj_z_pliku(okno):  # funkcja wybierająca plik
    for widget in okno.winfo_children():
        widget.destroy()
    global size  # potrzebuje wielkosci, żeby widgety miały odpowiednie wymiary
    size = (okno.winfo_width(), okno.winfo_height())
    for x in range(10):
        okno.columnconfigure(x, weight=0)
    for x in range(10):
        okno.rowconfigure(x, weight=0)
    for x in range(3):
        okno.columnconfigure(x, weight=1)
    for x in range(4):
        okno.rowconfigure(x, weight=1)
    okno.rowconfigure(3, weight=5)
    okno.columnconfigure(4, weight=1)
    # elementy okna:
    napisTytul = CTkLabel(okno, text="wpisz lokalizacje pliku:")
    sciezkaWejscie = CTkEntry(okno, width=int(size[0] * 10 / 13))
    przyciskPrzegladaj = CTkButton(okno, text="Przegladaj", command=lambda: przegladaj(okno))
    przyciskSzukaj = CTkButton(okno, text="Szukaj", command=lambda: otwieranie_pliku(okno, sciezkaWejscie.get()))
    przyciskWroc = CTkButton(okno, text="Wstecz", command=lambda: interfejs(okno))
    wzor, ostatni = operacje_chemiczne.wzor_lancucha_aminokwasow("MHFVDR")
    # tworze obiekt molecule który później rysuje
    mol = Chem.MolFromSmiles(wzor)
    rdCoordGen.AddCoords(mol)
    # podświetlam łańcuch główny
    hit_bonds = []
    for x in range(ostatni):
        hit_bonds.append(x)
    atomy = (0, ostatni)
    # ostatecznie go podświetlam
    if (czy_podswietlaj):
        obrazStartowy = Draw.MolToImage(mol, highlightBonds=hit_bonds, highlightAtoms=atomy, highlightColor=((0, 1, 0)),
                                        size=(int(size[0] * 10 / 13), int(size[
                                                                              1]) - 100))  # wyrażenie dostosowywuje rozmiar obrazka do okna, niestety nie dynamicznie :(
    else:
        obrazStartowy = Draw.MolToImage(mol, size=(int(size[0] * 10 / 13), int(size[1]) - 100))
    img1 = CTkImage(light_image=obrazStartowy, dark_image=obrazStartowy,
                    size=(int(size[0] * 10 / 13), int(size[1]) - 100))
    # wrzucam rysunek w okno
    obraz = CTkLabel(okno, text="", image=img1)
    obraz.grid(row=3, column=0, columnspan=4, sticky=N)
    napisTytul.grid(row=0, column=0, sticky=W + E)
    sciezkaWejscie.grid(row=1, column=0, rowspan=2, columnspan=4, sticky=W + E)
    przyciskPrzegladaj.grid(row=1, column=4, pady=3, sticky=W + E)
    przyciskSzukaj.grid(row=2, column=4, pady=3, sticky=W + E)
    przyciskWroc.grid(row=3, column=4, sticky=W + E)


def przegladaj(okno):  # funkcja pomocnicza do otwierania menu wyboru pliku

    file_path = filedialog.askopenfilename(filetypes=[("All Files", "*.*")])
    if (file_path == ""):
        wczytaj_z_pliku(okno)
    else:
        otwieranie_pliku(okno, file_path)


def otwieranie_pliku(okno, file_path):  # funkcja wybierająca konkretny ciag z pliku
    for widget in okno.winfo_children():
        widget.destroy()
    for x in range(10):
        okno.columnconfigure(x, weight=0)
    for x in range(10):
        okno.rowconfigure(x, weight=0)

    okno.columnconfigure(0, weight=4)
    okno.columnconfigure(1, weight=1)
    okno.rowconfigure(0, weight=1)
    size = (okno.winfo_width(), okno.winfo_height())
    ramka = CTkScrollableFrame(okno, width=(size[0] / 15) * 12, height=(size[1]))
    Bwidth = (size[0] / 15) * 12
    ramka.grid(row=0, column=0, columnspan=4, sticky=E + W + N + S)
    przyciskWroc = CTkButton(okno, text="Wstecz", command=lambda: wczytaj_z_pliku(okno))
    przyciskWroc.grid(row=0, column=4)
    ramka.columnconfigure(0, weight=1)
    if file_path.endswith(".txt"):
        with open(file_path, 'r') as plik:
            zawartosc = plik.read()
        if zawartosc.startswith(">"):  # jeśli plik tekstowy zaczyna się od > (układ pliku tytuł + zawartość)
            headers = []
            sequences = []
            with open(file_path, 'r') as file:
                seq = ""
                for line in file:
                    line = line.strip()
                    if line.startswith(">"):
                        line = line[1:]
                        headers.append(line)
                        if seq:
                            sequences.append(seq)
                            seq = ""
                    else:
                        seq += line
                sequences.append(seq)

            i = 0
            for string in headers:
                x = CTkButton(ramka, text=headers[i], width=(size[0] / 15) * 12 - 50,
                              command=lambda numer=i: otwieranie_lancucha(okno, file_path, sequences[numer]))
                x.grid(row=i, column=0, sticky=W + E)
                i += 1
        else:  # jeśli plik nie zaczyna się od > (same geny bez opisu)
            with open(file_path, 'r') as plik2:
                i = 0
                for linia in plik2:
                    i_str = str(i + 1)
                    linia_stripped = linia.strip()
                    b = CTkButton(ramka, text=("genom " + i_str),
                                  command=lambda k=linia_stripped: otwieranie_lancucha(okno, file_path, k))
                    b.grid(row=i, column=0, sticky=W + E)
                    i += 1

    elif (file_path.endswith(".fasta") or file_path.endswith(".fast") or file_path.endswith(
            ".seq") or file_path.endswith(".fa")
          or file_path.endswith(".fsa") or file_path.endswith(".nt") or file_path.endswith(".aa")):  # plik fasta
        i = 0
        with open(file_path, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                x = CTkButton(ramka, text=record.id, width=(size[0] / 15) * 12 - 50,
                              command=lambda sekwencja=record.seq: otwieranie_lancucha(okno, file_path, sekwencja))
                x.grid(row=i, column=0, sticky=W + E)
                i += 1
    elif (file_path.endswith(".gb") or file_path.endswith(".gbk")):
        i = 0
        with open(file_path, "r") as handle:
            for record in SeqIO.parse(handle, "genbank"):
                x = CTkButton(ramka, text=record.id, width=(size[0] / 15) * 12 - 50,
                              command=lambda sekwencja=record.seq: otwieranie_lancucha(okno, file_path, sekwencja))
                x.grid(row=i, column=0, sticky=W + E)
                i += 1
    else:
        wczytaj_z_pliku(okno)


def rysuj_przyciski(ramka, i, bialko, przesuniecie):
    if (len(bialko) > 10):
        Label = CTkLabel(ramka, text=(
                    bialko[:5] + "..." + bialko[len(bialko) - 5:] + " P:" + przesuniecie + " Len: " + str(len(bialko))))
        Label.grid(row=i, padx=2, column=0)
    else:
        Label = CTkLabel(ramka, text=(bialko + " P:" + przesuniecie))
        Label.grid(row=i, padx=2, column=0)
    Button = CTkButton(ramka, text="G.szybko", width=(size[0] / 15) * 2,
                       command=lambda y=bialko: rysowania.rysuj_szybko_interface(okno, y, czy_podswietlaj))
    Button.grid(row=i, column=1, padx=2, sticky=W + E)
    Button = CTkButton(ramka, text="Generuj", width=(size[0] / 15) * 2,
                       command=lambda y=bialko: rysowania.rysuj_interface(okno, y, czy_podswietlaj))
    Button.grid(row=i, column=2, padx=2, sticky=W + E)
    Button = CTkButton(ramka, text="G.dziwnie", width=(size[0] / 15) * 2,
                       command=lambda y=bialko: rysowania.rysuj_dziwnie_interface(okno, y, czy_podswietlaj))
    Button.grid(row=i, column=3, padx=2, sticky=W + E)


def otwieranie_lancucha(okno, file_path, sequence):
    global size
    size = (okno.winfo_width(), okno.winfo_height())
    for x in range(10):
        okno.rowconfigure(x, weight=0)
    for x in range(10):
        okno.columnconfigure(x, weight=0)
    for x in range(5):
        okno.columnconfigure(x, weight=1)
    okno.rowconfigure(1, weight=1)
    okno.rowconfigure(0, weight=0)
    for widget in okno.winfo_children():
        widget.destroy()
    # tworzenie ramki:
    global aktualne_bialko
    aktualne_bialko = 0
    ramka = CTkFrame(okno, width=(size[0] / 15) * 12, height=(size[1] - 50))
    przyciskNext = CTkButton(ramka, text="Następna", command=lambda: rysuj_bialka(ramka, bialka, aktualne_bialko + 10))
    przyciskBack = CTkButton(ramka, text="Poprzednia",
                             command=lambda: rysuj_bialka(ramka, bialka, aktualne_bialko - 10))
    Bwidth = (size[0] / 15) * 12
    ramka.grid(row=1, column=0, columnspan=4, sticky=E + W + N + S)
    for x in range(12):
        ramka.rowconfigure(x, weight=1)
    for x in range(4):
        ramka.columnconfigure(x, weight=1)
    lanc1, lanc2, lanc3 = operacje_chemiczne.translacja_bez_bugow(sequence)  # obrobione lancuchy
    bialka1 = operacje_chemiczne.rozklad_na_bialka(lanc1)  # 3 łancuchy białek dla 3 przesunięć
    bialka2 = operacje_chemiczne.rozklad_na_bialka(lanc2)
    bialka3 = operacje_chemiczne.rozklad_na_bialka(lanc3)
    dl1=len(bialka1)
    dl2=len(bialka2)
    dl3 = len(bialka3)
    if(dl1==1):
        text1 = " białko dla łańcucha zaczynającego \n się od 0 aminokwasu"
    else:
        text1=" białek dla łańcucha zaczynającego \n się od 0 aminokwasu"
    Napis="W łańcuchu jest:"+"\n"+"- "+str(dl1)+text1+"\n"\
          +"- "+str(dl2)+" od 1 aminokwasu"+"\n"+"- "+str(dl3)+" od 2 aminokwasu"
    bialka = []
    liczbaBialek=CTkLabel(okno,text=Napis)
    liczbaBialek.grid(row=1, column=4, sticky=W+E)
    for x in bialka1:
        bialka.append([x, "0"])
    for x in bialka2:
        bialka.append([x, "1"])
    for x in bialka3:
        bialka.append([x, "2"])
    rysuj_bialka(ramka, bialka, 0)
    opis = CTkLabel(okno, text="Wybierz co chcesz przeanalizować:")
    if (file_path == "x"):
        przycisk_wroc = CTkButton(okno, text="Wstecz", command=lambda: wczytaj_recznie(okno, sequence))
        przycisk_wroc.grid(row=0, column=4, sticky=W + E, padx=3, pady=25)
    else:
        przycisk_wroc = CTkButton(okno, text="Wstecz", command=lambda: otwieranie_pliku(okno, file_path))
        przycisk_wroc.grid(row=0, column=4, sticky=W + E, padx=3, pady=25)
    opis.grid(row=0, column=0, sticky=W + E)


def rysuj_bialka(ramka, bialka, i):
    global aktualne_bialko
    for widget in ramka.winfo_children():
        widget.destroy()
    przyciskNext = CTkButton(ramka, text="Następne", command=lambda: rysuj_bialka(ramka, bialka, aktualne_bialko + 10))
    przyciskBack = CTkButton(ramka, text="Poprzednie",
                             command=lambda: rysuj_bialka(ramka, bialka, aktualne_bialko - 10))
    wpiszStrone = CTkEntry(ramka, width=5, height=1)
    przyciskIdzDo = CTkButton(ramka, text="Idź do:",
                              command=lambda: rysuj_bialka(ramka, bialka, int(wpiszStrone.get()) * 10 - 10))
    przyciskIdzDo.grid(row=12, column=0, pady=10, sticky=E + W)
    wpiszStrone.grid(row=12, column=1, pady=10, sticky=E + W)
    przyciskNext.grid(row=11, column=2, pady=10, sticky=E + W)
    przyciskBack.grid(row=11, column=0, pady=10, sticky=E + W)
    if (len(bialka) >= 10):
        if (i >= 0 and i < len(bialka)):
            aktualne_bialko = i
        strona = CTkLabel(ramka,
                          text="Strona: " + str(int(aktualne_bialko / 10) + 1) + "/" + str(math.ceil(len(bialka) / 10)))
        strona.grid(row=11, column=1)
        bialka_w_ramce = bialka[aktualne_bialko:aktualne_bialko + 10]
        if ((len(bialka) - aktualne_bialko) < 10):
            for x in range(len(bialka) - aktualne_bialko):
                rysuj_przyciski(ramka, x, bialka_w_ramce[x][0], bialka_w_ramce[x][1])
        else:
            for x in range(10):
                rysuj_przyciski(ramka, x, bialka_w_ramce[x][0], bialka_w_ramce[x][1])
    if (len(bialka) < 10):
        aktualne_bialko = 0
        strona = CTkLabel(ramka, text="Strona: 1/1")
        strona.grid(row=11, column=1)
        for x in range(len(bialka)):
            rysuj_przyciski(ramka, x, bialka[x][0], bialka[x][1])


def interfaceOpcje(okno):  # funkcja wczytująca interface opcji
    for widget in okno.winfo_children():
        widget.destroy()
    # ustawiamy wagi
    for x in range(10):
        okno.columnconfigure(x, weight=0)
    for x in range(10):
        okno.rowconfigure(x, weight=0)
    okno.columnconfigure(0, weight=1)
    for x in range(8):
        okno.rowconfigure(x, weight=1)
    # tworze widgety
    napisMode = CTkLabel(okno, text="Zmień kolor tła:")
    napisKolor = CTkLabel(okno, text="Zmień kolor przycisków:")
    napisWroc = CTkLabel(okno, text="Wróć do menu:")
    napisPodswietlaj = CTkLabel(okno, text="Czy podświetlać główny łańcuch:")
    przyciskMode = CTkButton(okno, text=mode,
                             command=lambda: zmien_tryb(okno))  # zmienia tryb jasny na ciemny i viceversa
    przyciskKolor = CTkButton(okno, text=base_color, command=lambda: zmien_kolor(okno))  # zmienia bazowy kolor
    przyciskWroc = CTkButton(okno, text="Wstecz", command=lambda: zapis(okno))
    if (czy_podswietlaj == True):  # zmienia czy podswietla się główny łańcuch
        przyciskPodswietlaj = CTkButton(okno, text="Tak", command=lambda: podswietlenie(okno))
    else:
        przyciskPodswietlaj = CTkButton(okno, text="Nie", command=lambda: podswietlenie(okno))
    # wrzucam widgety w okno
    przyciskMode.grid(row=1, column=0, sticky=W + E, padx=100)
    przyciskKolor.grid(row=3, column=0, sticky=W + E, padx=100)
    przyciskPodswietlaj.grid(row=5, column=0, sticky=W + E, padx=100)
    przyciskWroc.grid(row=7, column=0, sticky=W + E, padx=100)
    napisMode.grid(row=0, column=0, sticky=W + E, padx=100)
    napisKolor.grid(row=2, column=0, sticky=W + E, padx=100)
    napisPodswietlaj.grid(row=4, column=0, sticky=W + E, padx=100)
    napisWroc.grid(row=6, column=0, sticky=W + E, padx=100)


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
    for widget in okno.winfo_children():
        widget.destroy()

    global size
    size = (okno.winfo_width(), okno.winfo_height())
    # ustawiam wagi kolumn
    for x in range(10):
        okno.columnconfigure(x, weight=0)
    for x in range(10):
        okno.rowconfigure(x, weight=0)
    for x in range(5):
        okno.columnconfigure(x, weight=2)
    okno.columnconfigure(5, weight=3)
    okno.rowconfigure(0, weight=1)
    okno.rowconfigure(1, weight=1)
    okno.rowconfigure(2, weight=5)

    ramka = CTkScrollableFrame(okno, width=500, height=1000)
    ramka.grid(row=1, column=0, columnspan=4, sticky=E + W + N + S)
    for x in range(5):
        ramka.columnconfigure(x, weight=1)

    # tworze combobox
    selected_month = StringVar()
    month_cb = ttk.Combobox(ramka, textvariable=selected_month)
    month_cb['values'] = ["Wstep", "Opcje", "Dane Ręcznie", "Dane z pliku", "Generowanie", "Wykresy"]
    month_cb['state'] = 'readonly'
    month_cb.grid(row=1, column=3, sticky=W + E)
    def month_changed(event):
        if(selected_month.get()=="Wstep"):
            img1 = CTkImage(light_image=Image.open("podswietl.PNG"),
                            dark_image=Image.open("podswietl.PNG"),
                            size=(600, 800))
            obraz = CTkLabel(ramka, text="", image=img1)
            obraz.grid(row=2, column=2, rowspan=3, sticky=NW)
        if (selected_month.get() == "Opcje"):
            img1 = CTkImage(light_image=Image.open("opcje.PNG"),
                        dark_image=Image.open("opcje.PNG"),
                        size=(600, 800))
            obraz = CTkLabel(ramka, text="", image=img1)
            obraz.grid(row=2, column=2, rowspan=3, sticky=NW)
        if (selected_month.get() == "Dane Ręcznie"):
            img1 = CTkImage(light_image=Image.open("twojastara.PNG"),
                        dark_image=Image.open("twojastara.PNG"),
                        size=(600, 800))
            obraz = CTkLabel(ramka, text="", image=img1)
            obraz.grid(row=2, column=2, rowspan=3, sticky=NW)
        if (selected_month.get() == "Dane z pliku"):
            img1 = CTkImage(light_image=Image.open("opcje.PNG"),
                            dark_image=Image.open("opcje.PNG"),
                            size=(600, 800))
            obraz = CTkLabel(ramka, text="", image=img1)
            obraz.grid(row=2, column=2, rowspan=3, sticky=NW)
        if (selected_month.get() == "Generowanie"):
            img1 = CTkImage(light_image=Image.open("opcje.PNG"),
                            dark_image=Image.open("opcje.PNG"),
                            size=(600, 800))
            obraz = CTkLabel(ramka, text="", image=img1)
            obraz.grid(row=2, column=2, rowspan=3, sticky=NW)
        if (selected_month.get() == "Wykresy"):
            img1 = CTkImage(light_image=Image.open("opcje.PNG"),
                            dark_image=Image.open("opcje.PNG"),
                            size=(600, 800))
            obraz = CTkLabel(ramka, text="", image=img1)
            obraz.grid(row=2, column=2, rowspan=3, sticky=NW)

    month_cb.bind('<<ComboboxSelected>>', month_changed)

    napis = CTkLabel(ramka, text="Instrukcja", fg_color=("white","black"))
    napis.grid(row=0, column=2, sticky=W + E)
    img1 = CTkImage(light_image=Image.open("podswietl.png"),
                    dark_image=Image.open("podswietl.png"),
                    size=(600, 800))
    obraz = CTkLabel(ramka, text="", image=img1)

    przerwa = CTkLabel(ramka, text="", width=15)
    przerwa.grid(row=0, column=0, sticky=W + E)

    obraz.grid(row=2, column=2, rowspan=3, sticky=NW)
    wrocButton = CTkButton(ramka, text="Wroc", width=15, command=lambda: interfejs(okno))
    wrocButton.grid(row=0, column=3, sticky=W + E)

    ramka.mainloop()
if __name__ == '__main__':
    oknoR = CTk()
    okno = CTkFrame(oknoR)
    oknoR.geometry("650x500")
    oknoR.columnconfigure(0, weight=1)
    oknoR.rowconfigure(0, weight=1)
    okno.grid(row=0, column=0, sticky=W + E + N + S)

    oknoR.minsize(height=500, width=750)
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
