import customtkinter
from customtkinter import *
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqUtils.ProtParam import ProtParamData
import operacje_chemiczne
import main
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg


def pH_bialka(lanc):
    analizuj = ProteinAnalysis(lanc)
    pH = (analizuj.isoelectric_point())
    if pH < 5.0:
        return "kwasowe"
    elif 5.0 <= pH <= 6.5:
        return "obojetne"
    elif pH > 6.5:
        return "zasadowe"


def dane_wykres(lanc):
    dane = []
    # suwak = slider()    w przyszłosci funkcja do ustawiania drugiej zmiennej tego gowna na dole
    analizuj = ProteinAnalysis(lanc)  # nie ma za uj pojecia co to sa te skale ale brzmi madrze
    dane.append(analizuj.protein_scale(ProtParamData.kd, 3))  # indeks hydrofobowy
    dane.append(analizuj.protein_scale(ProtParamData.em, 3))  # surface accessibility
    dane.append(analizuj.protein_scale(ProtParamData.Flex, 3))  # Normalized flexibility parameters
    dane.append(analizuj.protein_scale(ProtParamData.ja, 3))  # Janin Interior to surface transfer energy scale
    dane.append(analizuj.protein_scale(operacje_chemiczne.gra, 3))  # instability index by Grantham R.
    dane.append(len(dane[0]))  # pomaga stworzyc oś X w wykresie
    return dane


def wykresy(okno, wzor, lancpowrotny, lanc_Kodonow):  # robocza funkcja do wykresów

    for widget in okno.winfo_children():
        widget.destroy()

    analizuj = ProteinAnalysis(lancpowrotny)  # tablica do analizowania bialek
    # pojedyncze wartości
    dw = []
    dw.append(analizuj.secondary_structure_fraction())
    # zwraca tablice z 3 wartościami, ktore zawieraja %: sheets, helixes, turns cokolwiek by to nie było XD
    dw.append(analizuj.isoelectric_point())
    dw.append(analizuj.instability_index())
    dw.append(pH_bialka(lancpowrotny))
    if (analizuj.instability_index() <= 40):
        dw.append("Białko stabilne")
    elif (analizuj.instability_index() > 40):
        dw.append("Białko niestabilne")
    dane_kwasy = analizuj.get_amino_acids_percent()  # ilosc kazdego z kwasow
    lista = []
    for i in dane_kwasy:   dane_kwasy[i] = dane_kwasy[i] * 100
    for i in dane_kwasy:
        lista.append(dane_kwasy.get(i))
    dane = {'Aminokwas': ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W',
                          'Y', ],
            'Częstotliwość występowania': lista}

    Obraz_aminokwasów = plt.Figure(figsize=(10, 4), dpi=50)
    Wypisz_wykres_aminokwasów = Obraz_aminokwasów.add_subplot(1, 1, 1)
    Wykres_aminokwasów = FigureCanvasTkAgg(Obraz_aminokwasów, okno)
    Wykres_aminokwasów.get_tk_widget().grid(row=1, column=0)
    Dane_aminokwasów = pd.DataFrame(dane)
    Dane_aminokwasów = Dane_aminokwasów[['Aminokwas', 'Częstotliwość występowania']].groupby('Aminokwas').sum()
    Dane_aminokwasów.plot(kind='bar', legend=False, ax=Wypisz_wykres_aminokwasów)
    Wypisz_wykres_aminokwasów.set_title('Częstotliwość występowania aminokwasów w białku w procentach')

    dane_do_wykresu = dane_wykres(lanc_Kodonow)
    wartości_hydrofobia = (dane_do_wykresu[0])
    wartości_dostępność = (dane_do_wykresu[1])
    wartości_parametry = (dane_do_wykresu[2])
    wartości_skala = (dane_do_wykresu[3])
    wartości_niestabilność = (dane_do_wykresu[4])
    dlugosc = len(wartości_hydrofobia)
    aminokwas = []
    for x in range(dlugosc):
        aminokwas.append(x)
    dane_hydrofobia = {'Aminokwasy': aminokwas, 'Indeks hydrofobowy': wartości_hydrofobia}
    dane_dostępność = {'Aminokwasy': aminokwas, 'Dostępność powierzchniowa': wartości_dostępność}
    dane_parametry = {'Aminokwasy': aminokwas, 'Zogólnione parametry elastyczności': wartości_parametry}
    dane_skala = {'Aminokwasy': aminokwas, 'Skala Janin transferu energii': wartości_skala}
    dane_niestabilność = {'Aminokwasy': aminokwas, 'indeks niestabilności Granthama R.': wartości_niestabilność}
    Obraz_cech = plt.Figure(figsize=(10, 5), dpi=50)
    Wypisz_wykres_cech = Obraz_cech.add_subplot(1, 1, 1)
    Wykres_szczegółowy = FigureCanvasTkAgg(Obraz_cech, okno)
    Wykres_szczegółowy.get_tk_widget().grid(row=2, column=0)
    Dane_indeks_hydrofobowy = pd.DataFrame(dane_hydrofobia)
    Dane_indeks_hydrofobowy = Dane_indeks_hydrofobowy[['Aminokwasy', 'Indeks hydrofobowy']].groupby('Aminokwasy').sum()
    Dane_indeks_hydrofobowy.plot(kind='line', legend=True, ax=Wypisz_wykres_cech)
    Dane_dostępność_powierzchniowa = pd.DataFrame(dane_dostępność)
    Dane_dostępność_powierzchniowa = Dane_dostępność_powierzchniowa[
        ['Aminokwasy', 'Dostępność powierzchniowa']].groupby('Aminokwasy').sum()
    Dane_dostępność_powierzchniowa.plot(kind='line', legend=True, ax=Wypisz_wykres_cech)
    Dane_zogólnione_parametry_elastyczności = pd.DataFrame(dane_parametry)
    Dane_zogólnione_parametry_elastyczności = Dane_zogólnione_parametry_elastyczności[
        ['Aminokwasy', 'Zogólnione parametry elastyczności']].groupby('Aminokwasy').sum()
    Dane_zogólnione_parametry_elastyczności.plot(kind='line', legend=True, ax=Wypisz_wykres_cech)
    Dane_skala_transferu_energii_Janin = pd.DataFrame(dane_skala)
    Dane_skala_transferu_energii_Janin = Dane_skala_transferu_energii_Janin[
        ['Aminokwasy', 'Skala Janin transferu energii']].groupby('Aminokwasy').sum()
    Dane_skala_transferu_energii_Janin.plot(kind='line', legend=True, ax=Wypisz_wykres_cech)
    Dane_indeks_niestabilności = pd.DataFrame(dane_niestabilność)
    Dane_indeks_niestabilności = Dane_indeks_niestabilności[
        ['Aminokwasy', 'indeks niestabilności Granthama R.']].groupby('Aminokwasy').sum()
    Dane_indeks_niestabilności.plot(kind='line', legend=True, ax=Wypisz_wykres_cech)
    Wypisz_wykres_cech.set_title('Cechy fizyczno-chemiczne aminokwasów')
    ramka1 = customtkinter.CTkFrame(master=okno, width=100, height=200)
    napis = CTkLabel(okno, text=wzor)
    returnButton = CTkButton(master=ramka1, text="Wróć", command=lambda: main.wczytaj_recznie(okno, lancpowrotny), width=100)
    AtomyButton = CTkButton(master=ramka1, text="Ilość Atomów",
                            command=lambda: WypiszAtomy(okno, wzor, lanc_Kodonow, lancpowrotny), width=100)
    napis.grid(row=0, column=0, columnspan=2, sticky=E+W)
    ramka1.grid(row=1,column=1)
    returnButton.grid(row=0, column=0, sticky=W + E, pady=20)
    AtomyButton.grid(row=2, column=0, sticky=W + E)


def WypiszAtomy(okno, wzor, lancpowrotny, lanc_Kodonow):
    for widget in okno.winfo_children():
        widget.destroy()
    operacje_chemiczne.listaAtomow(wzor)
    returnButton = CTkButton(okno, text="Wróć", command=lambda: wykresy(okno, wzor, lanc_Kodonow, lancpowrotny))
    returnButton.grid(row=1, column=0)
