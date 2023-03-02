import customtkinter
import tkinter
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
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors as _rdMolDescriptors
MolWt = lambda *x,**y:_rdMolDescriptors._CalcMolWt(*x,**y)


def pH_bialka(lanc):
    analizuj = ProteinAnalysis(lanc)
    pH = (analizuj.isoelectric_point())
    if pH < 5.0:
        return "kwasowe"
    elif 5.0 <= pH <= 6.5:
        return "obojętne"
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

def checkbox_event1(okno,lanc_Kodonow):
    global c1
    if(c1==True):
        c1=False
    else:
        c1=True
    rysowanie_wykresu(okno,lanc_Kodonow)

def checkbox_event2(okno,lanc_Kodonow):
    global c2
    if (c2 == True):
        c2 = False
    else:
        c2 = True
    rysowanie_wykresu(okno,lanc_Kodonow)

def checkbox_event3(okno,lanc_Kodonow):
    global c3
    if (c3 == True):
        c3 = False
    else:
        c3 = True
    rysowanie_wykresu(okno,lanc_Kodonow)

def checkbox_event4(okno,lanc_Kodonow):
    global c4
    if (c4 == True):
        c4 = False
    else:
        c4 = True
    rysowanie_wykresu(okno,lanc_Kodonow)

def checkbox_event5(okno,lanc_Kodonow):
    global c5
    if (c5 == True):
        c5 = False
    else:
        c5 = True
    rysowanie_wykresu(okno,lanc_Kodonow)

def rysowanie_wykresu(okno, lanc_Kodonow):
    global c1,c2,c3,c4,c5
    dane_do_wykresu = dane_wykres(lanc_Kodonow)
    wartości_hydrofobia = (dane_do_wykresu[0])
    wartości_dostępność = (dane_do_wykresu[1])
    wartości_parametry = (dane_do_wykresu[2])
    wartości_skala = (dane_do_wykresu[3])
    wartości_niestabilność = (dane_do_wykresu[4])
    dlugosc = len(wartości_hydrofobia)
    aminokwas = []
    for x in range(dlugosc):
        aminokwas.append(x+1)
    dane_hydrofobia = {'Aminokwasy': aminokwas, 'Indeks hydrofobowy': wartości_hydrofobia}
    dane_dostępność = {'Aminokwasy': aminokwas, 'Dostępność powierzchniowa': wartości_dostępność}
    dane_parametry = {'Aminokwasy': aminokwas, 'Zogólnione parametry elastyczności': wartości_parametry}
    dane_skala = {'Aminokwasy': aminokwas, 'Skala Janin transferu energii': wartości_skala}
    dane_niestabilność = {'Aminokwasy': aminokwas, 'indeks niestabilności Granthama R.': wartości_niestabilność}
    Obraz_cech = plt.Figure(figsize=(10, 5), dpi=50)
    Wypisz_wykres_cech = Obraz_cech.add_subplot(1, 1, 1)
    Wykres_szczegółowy = FigureCanvasTkAgg(Obraz_cech, okno)
    Wykres_szczegółowy.get_tk_widget().grid(row=2, column=0,sticky=N+W+E+S)

    Dane_indeks_hydrofobowy = pd.DataFrame(dane_hydrofobia)
    Dane_indeks_hydrofobowy = Dane_indeks_hydrofobowy[['Aminokwasy', 'Indeks hydrofobowy']].groupby('Aminokwasy').sum()
    if (c1):
        Dane_indeks_hydrofobowy.plot(kind='line', legend=False, ax=Wypisz_wykres_cech, color='#581266')
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
    Wypisz_wykres_cech.xaxis.set_major_locator(MultipleLocator(1))
    Wypisz_wykres_cech.yaxis.set_minor_locator(MultipleLocator(1))

    Wypisz_wykres_cech.set_title('Cechy fizyczno-chemiczne aminokwasów')

def wykresy(okno, wzor, lanc_Kodonow, lancpowrotny):  # robocza funkcja do wykresów
    global c1,c2,c3,c4,c5
    c1 = 1
    c2 = 1
    c3 = 1
    c4 = 1
    c5 = 1
    for widget in okno.winfo_children():
        widget.destroy()
    for x in range(10):
        okno.rowconfigure(x,weight=0)
    for x in range(10):
        okno.columnconfigure(x,weight=0)
    okno.rowconfigure(1,weight=1)
    okno.rowconfigure(2, weight=1)
    okno.columnconfigure(0, weight=1)
    okno.columnconfigure(1, weight=1)
    analizuj = ProteinAnalysis(lanc_Kodonow)  # tablica do analizowania bialek
    # pojedyncze wartości

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
    Wykres_aminokwasów.get_tk_widget().grid(row=1, column=0,sticky=N+W+E+S)
    Dane_aminokwasów = pd.DataFrame(dane)
    Dane_aminokwasów = Dane_aminokwasów[['Aminokwas', 'Częstotliwość występowania']].groupby('Aminokwas').sum()
    Dane_aminokwasów.plot(kind='bar', legend=False, ax=Wypisz_wykres_aminokwasów)
    Wypisz_wykres_aminokwasów.set_title('Częstotliwość występowania aminokwasów w białku w procentach')


    ramka1 = customtkinter.CTkFrame(master=okno)
    ramka_wykres = customtkinter.CTkFrame(master=okno)
    napis = CTkLabel(okno, text="Właściwości Białka")
    returnButton = CTkButton(master=ramka1, text="Wróć", command=lambda: main.wczytaj_recznie(okno, lancpowrotny), width=100)
    DaneButton = CTkButton(master=ramka1, text="Dane", command=lambda: dane_interfejs(okno, wzor, lanc_Kodonow, lancpowrotny), width=100)
    napis.grid(row=0, column=0, columnspan=2, sticky=E+W)
    ramka1.grid(row=1,column=1,sticky=N+W+E+S)
    ramka_wykres.grid(row=2,column=1,sticky=N+W+E+S)
    returnButton.grid(row=0, column=0, pady=10)
    DaneButton.grid(row=2, column=0, pady=10)
    check_var1 = tkinter.IntVar(value=1)
    check_var2 = tkinter.IntVar(value=1)
    check_var3 = tkinter.IntVar(value=1)  #
    check_var4 = tkinter.IntVar(value=1)
    check_var5 = tkinter.IntVar(value=1)
    checkbox_indeks_hydrofobowy = customtkinter.CTkCheckBox(master=ramka_wykres, text="Indeks hydrofobowy",text_color='#581266', command=lambda: checkbox_event1(okno,lanc_Kodonow),
                                         variable=check_var1, onvalue=1, offvalue=0)
    checkbox_dostępność_powierzchniowa = customtkinter.CTkCheckBox(master=ramka_wykres, text="Dostępność powierzchniowa",text_color='#3fc2d1', command=lambda: checkbox_event2(okno,lanc_Kodonow),
                                         variable=check_var2, onvalue=1, offvalue=0)
    checkbox_zogólnione_parametry_elastyczności = customtkinter.CTkCheckBox(master=ramka_wykres, text="Zogólnione parametry elastyczności",text_color='#057a13', command=lambda: checkbox_event3(okno,lanc_Kodonow),
                                         variable=check_var3, onvalue=1, offvalue=0)
    checkbox_skala_transferu_energii_Janin = customtkinter.CTkCheckBox(master=ramka_wykres, text="Skala Janin transferu energii",text_color='#db1a1a', command=lambda: checkbox_event4(okno,lanc_Kodonow),
                                         variable=check_var4, onvalue=1, offvalue=0)
    checkbox_indeks_niestabilności = customtkinter.CTkCheckBox(master=ramka_wykres, text="Indeks niestabilności Granthama R.",text_color='#d0e831', command=lambda: checkbox_event5(okno,lanc_Kodonow),
                                         variable=check_var5, onvalue=1, offvalue=0)
    rysowanie_wykresu(okno,lanc_Kodonow)
    checkbox_indeks_hydrofobowy.grid(row=0,column=0,rowspan=2,sticky=W,pady=5,padx=5 )
    checkbox_dostępność_powierzchniowa.grid(row=2,column=0,rowspan=2,sticky=W,pady=5,padx=5 )
    checkbox_zogólnione_parametry_elastyczności.grid(row=4,column=0, rowspan=2,sticky=W,pady=5,padx=5 )
    checkbox_skala_transferu_energii_Janin.grid(row=6,column=0, rowspan=2,sticky=W,pady=5,padx=5 )
    checkbox_indeks_niestabilności.grid(row=8,column=0, rowspan=2,sticky=W,pady=5,padx=5)
def dane_interfejs(okno, wzor, lanc_Kodonow, lancpowrotny):
    for widget in okno.winfo_children():
        widget.destroy()
    for x in range(10):
        okno.rowconfigure(x,weight=0)
    for x in range(10):
        okno.columnconfigure(x,weight=0)
    analizuj = ProteinAnalysis(lanc_Kodonow)
    mol = Chem.MolFromSmiles(wzor)
    dw = []
    dw.append(analizuj.secondary_structure_fraction())
    # zwraca tablice z 3 wartościami, ktore zawieraja %: sheets, helixes, turns cokolwiek by to nie było XD
    dw.append(analizuj.isoelectric_point())
    dw.append(pH_bialka(lanc_Kodonow))
    dw.append(analizuj.instability_index())
    if (analizuj.instability_index() <= 40):
        dw.append("białko stabilne")
    elif (analizuj.instability_index() > 40):
        dw.append("białko niestabilne")
    dw.append(MolWt(Chem.MolFromSmiles(wzor)))
    struktura=dw[0]
    linijka_1 = CTkLabel(okno, text="Na strukturę drugorzędową tego białka składa się:")
    linijka_2 = CTkLabel(okno, text= str(struktura[0]*100)+"% harmonijek beta (pofałdowanej płaszczyzny)")
    linijka_3 = CTkLabel(okno, text= str(struktura[1]*100)+"% helis alfa (helis pi)")
    linijka_4 = CTkLabel(okno, text= str(struktura[2]*100)+"% beta zakrętów (pętli omega)")
    linijka_5 = CTkLabel(okno, text="Punkt izoelektryczny tego białka wynosi "+str(dw[1]))
    linijka_6 = CTkLabel(okno, text="W punkcie izoelektrycznym to białko jest "+dw[2])
    linijka_7 = CTkLabel(okno, text="Indeks niestabilności Guruprasada wynosi "+str(dw[3]))
    linijka_8 = CTkLabel(okno, text="Jest to "+dw[4])
    linijka_9 = CTkLabel(okno, text="Masa tego białka wynosi "+str(dw[5]))
    linijka_1.grid(row=0,column=0, sticky=W)
    linijka_2.grid(row=1, column=0, sticky=W)
    linijka_3.grid(row=2, column=0, sticky=W)
    linijka_4.grid(row=3, column=0, sticky=W)
    linijka_5.grid(row=4, column=0, sticky=W)
    linijka_6.grid(row=5, column=0, sticky=W)
    linijka_7.grid(row=6, column=0, sticky=W)
    linijka_8.grid(row=7, column=0, sticky=W)
    linijka_9.grid(row=8, column=0, sticky=W)


    returnButton = CTkButton(okno, text="Wróć", command=lambda: main.wczytaj_recznie(okno, lancpowrotny),
                             width=100)
    returnButton.grid(row=0, column=1)
