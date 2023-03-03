import customtkinter
import tkinter
from customtkinter import *
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqUtils.ProtParam import ProtParamData
import operacje_chemiczne
import main
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.ticker import (MultipleLocator)
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors as _rdMolDescriptors

MolWt = lambda *x, **y: _rdMolDescriptors._CalcMolWt(*x, **y)


def dane_wykres(lanc, suwak):
    dane = []
    window = int(suwak.get())  # suwak
    analizuj = ProteinAnalysis(lanc)  # analiza proteiny(białka)
    dane.append(analizuj.protein_scale(ProtParamData.kd, window))  # indeks hydrofobowy
    dane.append(analizuj.protein_scale(ProtParamData.em, window))  # surface accessibility
    dane.append(analizuj.protein_scale(ProtParamData.Flex, window))  # Normalized flexibility parameters
    dane.append(analizuj.protein_scale(ProtParamData.ja, window))  # Janin Interior to surface transfer energy scale
    dane.append(analizuj.protein_scale(operacje_chemiczne.gra, window))  # instability index by Grantham R.
    dane.append(len(dane[0]))  # pomaga stworzyc oś X w wykresie
    return dane


def checkbox_event1(okno, lanc_Kodonow):  # funkcja działająca po kliknięciu checkboxa
    global c1
    if (c1 == True):
        c1 = False
    else:
        c1 = True
    rysowanie_wykresu(okno, lanc_Kodonow)


def checkbox_event2(okno, lanc_Kodonow):
    global c2
    if (c2 == True):
        c2 = False
    else:
        c2 = True
    rysowanie_wykresu(okno, lanc_Kodonow)


def checkbox_event3(okno, lanc_Kodonow):
    global c3
    if (c3 == True):
        c3 = False
    else:
        c3 = True
    rysowanie_wykresu(okno, lanc_Kodonow)


def checkbox_event4(okno, lanc_Kodonow):
    global c4
    if (c4 == True):
        c4 = False
    else:
        c4 = True
    rysowanie_wykresu(okno, lanc_Kodonow)


def checkbox_event5(okno, lanc_Kodonow):
    global c5
    if (c5 == True):
        c5 = False
    else:
        c5 = True
    rysowanie_wykresu(okno, lanc_Kodonow)


def rysowanie_wykresu(okno, lanc_Kodonow):  # funkcja rysująca wykres ze zmiennymi danymi
    global c1, c2, c3, c4, c5
    dane_do_wykresu = dane_wykres(lanc_Kodonow, suwak)
    wartości_hydrofobia = (dane_do_wykresu[0])
    wartości_dostępność = (dane_do_wykresu[1])  # załadowanie danych o białku
    wartości_parametry = (dane_do_wykresu[2])
    wartości_skala = (dane_do_wykresu[3])
    wartości_niestabilność = (dane_do_wykresu[4])
    dlugosc = len(wartości_hydrofobia)  # funkcja licząca długość, którą będzie miał wykres
    plt.rcParams.update({'font.size': 18})
    aminokwas = []
    for x in range(dlugosc):
        aminokwas.append(x + 1)  # każdej sekwencji przypisujemy numer jednostkowy
    dane_hydrofobia = {'Aminokwasy': aminokwas, 'Indeks hydrofobowy': wartości_hydrofobia}
    dane_dostępność = {'Aminokwasy': aminokwas,
                       'Dostępność powierzchniowa': wartości_dostępność}  # zamiana danych na przystępne do załadowania do wykresu
    dane_parametry = {'Aminokwasy': aminokwas, 'Zogólnione parametry elastyczności': wartości_parametry}
    dane_skala = {'Aminokwasy': aminokwas, 'Skala Janin transferu energii': wartości_skala}
    dane_niestabilność = {'Aminokwasy': aminokwas, 'indeks niestabilności Granthama R.': wartości_niestabilność}
    Obraz_cech = plt.Figure(figsize=(10, 5), dpi=50)  # stworzenie figury, w której będzie wykres
    Wypisz_wykres_cech = Obraz_cech.add_subplot(1, 1, 1)  # dodanie wykresu do figury
    Wykres_szczegółowy = FigureCanvasTkAgg(Obraz_cech, okno)  # przełożenie wykresu na format TkIntera
    Wykres_szczegółowy.get_tk_widget().grid(row=2, column=0, sticky=N + W + E + S)  # ułożenie wykresu w oknie

    Dane_indeks_hydrofobowy = pd.DataFrame(dane_hydrofobia)  # załadowanie danych
    Dane_indeks_hydrofobowy = Dane_indeks_hydrofobowy[['Aminokwasy', 'Indeks hydrofobowy']].groupby(
        'Aminokwasy').sum()  # pogrupowanie danych
    if (c1):
        Dane_indeks_hydrofobowy.plot(kind='line', legend=False, ax=Wypisz_wykres_cech,
                                     color='#581266')  # narysowanie wykresu od danych, jeśli odpowiadający checkbox jest zaznaczony
    Dane_dostępność_powierzchniowa = pd.DataFrame(dane_dostępność)
    Dane_dostępność_powierzchniowa = Dane_dostępność_powierzchniowa[
        ['Aminokwasy', 'Dostępność powierzchniowa']].groupby('Aminokwasy').sum()
    if (c2):
        Dane_dostępność_powierzchniowa.plot(kind='line', legend=False, ax=Wypisz_wykres_cech, color='#3fc2d1')
    Dane_zogólnione_parametry_elastyczności = pd.DataFrame(dane_parametry)
    Dane_zogólnione_parametry_elastyczności = Dane_zogólnione_parametry_elastyczności[
        ['Aminokwasy', 'Zogólnione parametry elastyczności']].groupby('Aminokwasy').sum()
    if (c3):
        Dane_zogólnione_parametry_elastyczności.plot(kind='line', legend=False, ax=Wypisz_wykres_cech, color='#057a13')
    Dane_skala_transferu_energii_Janin = pd.DataFrame(dane_skala)
    Dane_skala_transferu_energii_Janin = Dane_skala_transferu_energii_Janin[
        ['Aminokwasy', 'Skala Janin transferu energii']].groupby('Aminokwasy').sum()
    if (c4):
        Dane_skala_transferu_energii_Janin.plot(kind='line', legend=False, ax=Wypisz_wykres_cech, color='#db1a1a')
    Dane_indeks_niestabilności = pd.DataFrame(dane_niestabilność)
    Dane_indeks_niestabilności = Dane_indeks_niestabilności[
        ['Aminokwasy', 'indeks niestabilności Granthama R.']].groupby('Aminokwasy').sum()
    if (c5):
        Dane_indeks_niestabilności.plot(kind='line', legend=False, ax=Wypisz_wykres_cech, color='#d0e831')
    Wypisz_wykres_cech.xaxis.set_major_locator(MultipleLocator(1))  # wypisanie wszystkich indeksów aminokwasów na osi x
    Wypisz_wykres_cech.yaxis.set_minor_locator(
        MultipleLocator(1))  # ustawienie pomniejszych 'ticków' na osi y w odległości 1 od siebie

    Wypisz_wykres_cech.set_title('Cechy fizyczno-chemiczne aminokwasów')  # tytuł

    analizuj = ProteinAnalysis(lanc_Kodonow)  # tablica do analizowania bialek
    dane_kwasy = analizuj.get_amino_acids_percent()  # ilosc kazdego z aminokwasow
    lista = []  # częstotliwość występowania aminokwasów w tabeli
    for i in dane_kwasy:   dane_kwasy[i] = dane_kwasy[i] * 100
    for i in dane_kwasy:
        lista.append(dane_kwasy.get(i))
    dane = {'Aminokwas': ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W',
                          'Y', ],
            'Częstotliwość występowania': lista}

    Obraz_aminokwasów = plt.Figure(figsize=(10, 4), dpi=50)  # utworzenie figury
    Wypisz_wykres_aminokwasów = Obraz_aminokwasów.add_subplot(1, 1, 1)  # dodanie wykresu do figury
    Wykres_aminokwasów = FigureCanvasTkAgg(Obraz_aminokwasów, okno)  # zamienienie wykresu na format TkIntera
    Wykres_aminokwasów.get_tk_widget().grid(row=1, column=0, sticky=N + W + E + S)  # ułożenie wykresu w oknie
    Dane_aminokwasów = pd.DataFrame(dane)  # dodanie danych
    Dane_aminokwasów = Dane_aminokwasów[['Aminokwas', 'Częstotliwość występowania']].groupby(
        'Aminokwas').sum()  # grupowanie danych
    Dane_aminokwasów.plot(kind='bar', legend=False, ax=Wypisz_wykres_aminokwasów, rot=0)  # wypisanie danych na wykresie
    Wypisz_wykres_aminokwasów.set_title('Częstotliwość występowania aminokwasów w białku w procentach')  # tytuł


def wykresy(ramex, wzor, lanc_Kodonow):  # robocza funkcja do wykresów
    global c1, c2, c3, c4, c5, suwak
    c1 = 1  # zmienna oznaczająca stan checkboxa
    c2 = 1
    c3 = 1
    c4 = 1
    c5 = 1
    okno = CTkFrame(ramex)
    okno.grid(row=0, column=0, rowspan=10, columnspan=10, sticky=N + W + S + E)
    size = (okno.winfo_width(), okno.winfo_height())
    if (size[0] < 650 or size[1] < 500):  # przy pierwszym włączeniu jest bug funkcji winfo
        size = (650, 500)
    for widget in okno.winfo_children():  # czyscimy okno
        widget.destroy()
    for x in range(10):
        okno.rowconfigure(x, weight=0)
    for x in range(10):
        okno.columnconfigure(x, weight=0)
    okno.rowconfigure(1, weight=1)
    okno.rowconfigure(2, weight=1)
    okno.columnconfigure(0, weight=3)
    okno.columnconfigure(1, weight=1)

    ramka1 = customtkinter.CTkFrame(master=okno)  # ramka na przyciski
    ramka_wykres = customtkinter.CTkFrame(master=okno)  # ramka na checkboxy
    napis = CTkLabel(okno, text="Właściwości Białka")  # tytuł
    napiss1 = CTkLabel(ramka_wykres, text="Ustal rozmiar okna (długość interwału")
    napiss2 = CTkLabel(ramka_wykres, text="używanego do obliczania profilu wykresu)")
    suwak = CTkSlider(master=ramka_wykres, width=200, from_=2, to=len(lanc_Kodonow) - 1, orientation="horizontal")
    suwak.set(2)  # wartosc domyslna
    wykresButton = CTkButton(master=ramka_wykres, text="Rysuj", command=lambda: rysowanie_wykresu(okno, lanc_Kodonow))
    returnButton = CTkButton(master=ramka1, text="Wróć", command=lambda: zniszcz(okno), width=100)
    DaneButton = CTkButton(master=ramka1, text="Dane", command=lambda: dane_interfejs(okno, wzor, lanc_Kodonow),
                           width=100)
    napis.grid(row=0, column=0, columnspan=2, sticky=E + W)
    ramka1.grid(row=1, column=1, sticky=N + W + E + S)
    ramka_wykres.grid(row=2, column=1, sticky=N + W + E + S)
    ramka_wykres.columnconfigure(0, weight=1)
    for x in range(15):
        ramka_wykres.rowconfigure(x, weight=1)  # ustawienie wartości dla wszystich rzędów i kolumn
    ramka1.columnconfigure(0, weight=1)
    for x in range(5):
        ramka1.rowconfigure(x, weight=1)
    returnButton.grid(row=1, column=0, sticky=E + W)
    DaneButton.grid(row=3, column=0, sticky=E + W)

    check_var1 = tkinter.IntVar(value=1)  # zmienna bool w tkinterze
    check_var2 = tkinter.IntVar(value=1)
    check_var3 = tkinter.IntVar(value=1)
    check_var4 = tkinter.IntVar(value=1)
    check_var5 = tkinter.IntVar(value=1)
    checkbox_indeks_hydrofobowy = customtkinter.CTkCheckBox(master=ramka_wykres, text="Indeks hydrofobowy",
                                                            text_color='#581266',
                                                            command=lambda: checkbox_event1(okno, lanc_Kodonow),
                                                            variable=check_var1, onvalue=1, offvalue=0)  # checkboxy :)
    checkbox_dostępność_powierzchniowa = customtkinter.CTkCheckBox(master=ramka_wykres,
                                                                   text="Dostępność powierzchniowa",
                                                                   text_color='#3fc2d1',
                                                                   command=lambda: checkbox_event2(okno, lanc_Kodonow),
                                                                   variable=check_var2, onvalue=1, offvalue=0)
    checkbox_zogólnione_parametry_elastyczności = customtkinter.CTkCheckBox(master=ramka_wykres,
                                                                            text="Zogólnione parametry elastyczności",
                                                                            text_color='#057a13',
                                                                            command=lambda: checkbox_event3(okno,
                                                                                                            lanc_Kodonow),
                                                                            variable=check_var3, onvalue=1, offvalue=0)
    checkbox_skala_transferu_energii_Janin = customtkinter.CTkCheckBox(master=ramka_wykres,
                                                                       text="Skala Janin transferu energii",
                                                                       text_color='#db1a1a',
                                                                       command=lambda: checkbox_event4(okno,
                                                                                                       lanc_Kodonow),
                                                                       variable=check_var4, onvalue=1, offvalue=0)
    checkbox_indeks_niestabilności = customtkinter.CTkCheckBox(master=ramka_wykres,
                                                               text="Indeks niestabilności Granthama R.",
                                                               text_color='#d0e831',
                                                               command=lambda: checkbox_event5(okno, lanc_Kodonow),
                                                               variable=check_var5, onvalue=1, offvalue=0)
    rysowanie_wykresu(okno, lanc_Kodonow)
    checkbox_indeks_hydrofobowy.grid(row=0, column=0, rowspan=2, sticky=W + E, pady=5, padx=5)
    checkbox_dostępność_powierzchniowa.grid(row=2, column=0, rowspan=2, sticky=W + E, pady=5, padx=5)
    checkbox_zogólnione_parametry_elastyczności.grid(row=4, column=0, rowspan=2, sticky=W + E, pady=5, padx=5)
    checkbox_skala_transferu_energii_Janin.grid(row=6, column=0, rowspan=2, sticky=W + E, pady=5, padx=5)
    checkbox_indeks_niestabilności.grid(row=8, column=0, rowspan=2, sticky=W + E, pady=5, padx=5)
    napiss1.grid(row=10, column=0, rowspan=2, sticky=W + E, pady=5, padx=5)
    napiss2.grid(row=11, column=0, rowspan=2, sticky=W + E, pady=5, padx=5)
    suwak.grid(row=12, column=0, rowspan=2, sticky=W + E, pady=5, padx=5)
    wykresButton.grid(row=14, column=0, rowspan=2, sticky=W + E, pady=5, padx=5)


def dane_interfejs(ramex, wzor, lanc_Kodonow):
    okno = CTkFrame(ramex)
    okno.grid(row=0, column=0, rowspan=10, columnspan=10, sticky=N + W + S + E)
    size = (okno.winfo_width(), okno.winfo_height())
    if (size[0] < 650 or size[1] < 500):  # przy pierwszym włączeniu jest bug funkcji winfo, wynik to zawsze (200,200)
        size = (650, 500)
    for widget in okno.winfo_children():  # czyscimy okno
        widget.destroy()
    for widget in okno.winfo_children():
        widget.destroy()
    for x in range(10):
        okno.rowconfigure(x, weight=0)
    for x in range(10):
        okno.columnconfigure(x, weight=0)
    okno.columnconfigure(0, weight=3)
    analizuj = ProteinAnalysis(lanc_Kodonow)
    dw = []
    dw.append(analizuj.secondary_structure_fraction())
    dw.append(analizuj.isoelectric_point())
    dw.append(operacje_chemiczne.pH_bialka(lanc_Kodonow))  # obliczanie wszystkich potrzebnych danych
    dw.append(analizuj.instability_index())
    if (analizuj.instability_index() <= 40):
        dw.append("białko stabilne")
    elif (analizuj.instability_index() > 40):
        dw.append("białko niestabilne")
    dw.append(analizuj.molecular_weight())
    struktura = dw[0]
    linijka_1 = CTkLabel(okno, text="Na strukturę drugorzędową tego białka składa się:")
    linijka_2 = CTkLabel(okno, text=str(
        "- " + str(round(struktura[0] * 100, 3)) + " % harmonijek beta (pofałdowanej płaszczyzny)"))
    linijka_3 = CTkLabel(okno, text=str("- " + str(round(struktura[1] * 100, 3)) + " % helis alfa (helis pi)"))
    linijka_4 = CTkLabel(okno, text=str("- " + str(
        round(struktura[2] * 100, 3)) + " % beta zakrętów (pętli omega)"))  # wypisanie wszystkich potrzebnych danych
    linijka_5 = CTkLabel(okno, text="Punkt izoelektryczny tego białka wynosi " + str(round(dw[1], 3)))
    linijka_6 = CTkLabel(okno, text="W punkcie izoelektrycznym to białko jest " + dw[2])
    linijka_7 = CTkLabel(okno, text="Indeks niestabilności Guruprasada wynosi " + str(round(dw[3], 3)))
    linijka_8 = CTkLabel(okno, text="Jest to " + dw[4])
    linijka_9 = CTkLabel(okno, text="Masa tego białka wynosi " + str(round(dw[5], 3)) + " u")
    linijka_1.grid(row=0, column=0, sticky=W + E)
    linijka_2.grid(row=1, column=0, sticky=W + E)
    linijka_3.grid(row=2, column=0, sticky=W + E)
    linijka_4.grid(row=3, column=0, sticky=W + E)
    linijka_5.grid(row=4, column=0, sticky=W + E)  # ułożenie danych w oknie
    linijka_6.grid(row=5, column=0, sticky=W + E)
    linijka_7.grid(row=6, column=0, sticky=W + E)
    linijka_8.grid(row=7, column=0, sticky=W + E)
    linijka_9.grid(row=8, column=0, sticky=W + E)
    returnButton = CTkButton(okno, text="Wróć", command=lambda: zniszcz(okno),
                             width=200)
    returnButton.grid(row=9, column=0, pady=50)


def zniszcz(okno):
    okno.destroy()