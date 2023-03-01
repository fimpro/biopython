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

def WypiszAtomy(okno, wzor,lanc_powrotny, lanc_kodonow):
    for widget in okno.winfo_children():
        widget.destroy()
    operacje_chemiczne.listaAtomow(wzor)
    returnButton = CTkButton(okno, text="Wróć", command=lambda: wykresy(okno,wzor, lanc_powrotny, lanc_kodonow))
    returnButton.grid(row=1, column=0)

def dane_wykres(lanc):
    dane = []
    #suwak = slider()    w przyszłosci funkcja do ustawiania drugiej zmiennej tego gowna na dole
    analizuj = ProteinAnalysis(lanc)        #nie ma za uj pojecia co to sa te skale ale brzmi madrze
    dane.append(analizuj.protein_scale(ProtParamData.kd, 3))         #indeks hydrofobowy
    dane.append(analizuj.protein_scale(ProtParamData.em, 3))         #surface accessibility
    dane.append(analizuj.protein_scale(ProtParamData.Flex, 3))       #Normalized flexibility parameters
    dane.append(analizuj.protein_scale(ProtParamData.ja, 3))         #Janin Interior to surface transfer energy scale
    dane.append(analizuj.protein_scale(operacje_chemiczne.gra, 3))   #instability index by Grantham R.
    dane.append(len(dane[0]))  #pomaga stworzyc oś X w wykresie
    return dane




def wykresy(okno, wzor, lanc_kodonow, lanc_powrotny):  # robocza funkcja do wykresów
    global size
    size = (okno.winfo_width(), okno.winfo_height())
    for widget in okno.winfo_children():  # czyscimy okno
        widget.destroy()
    okno.columnconfigure(0, weight=100)
    okno.columnconfigure(1, weight=2)


    analizuj = ProteinAnalysis(lanc_kodonow)  # tablica do analizowania bialekdane_kwasy = []
    dw = [] #pojedyncze wartości
    dw.append(analizuj.secondary_structure_fraction())
    # zwraca tablice z 3 wartościami, ktore zawieraja %: sheets, helixes, turns cokolwiek by to nie było XD
    dw.append(analizuj.instability_index())
    dw.append(analizuj.isoelectric_point())
    dw.append(pH_bialka(lanc_kodonow))
    if (analizuj.instability_index() <= 40):
        dw.append("Białko stabilne")
    elif (analizuj.instability_index() > 40):
        dw.append("Białko niestabilne")
    print(lanc_powrotny)
    print(dw)
    print(dane_wykres(lanc_kodonow))
    dane_kwasy1 = analizuj.get_amino_acids_percent()  # ilosc kazdego z kwasow
    lista = []
    for i in dane_kwasy1:   dane_kwasy1[i] = dane_kwasy1[i] * 100
    for i in dane_kwasy1:
        lista.append(dane_kwasy1.get(i))
    dane = {'Aminokwas':['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y',],
            'Częstotliwość występowania':lista}

    figure1 = plt.Figure(figsize=(12, 4), dpi=50)
    ax1 = figure1.add_subplot(1,1,1)
    bar1 = FigureCanvasTkAgg(figure1, okno)
    bar1.get_tk_widget().grid(row=1, column=0)
    df1 = pd.DataFrame(dane)
    df1 = df1[['Aminokwas', 'Częstotliwość występowania']].groupby('Aminokwas').sum()
    df1.plot(kind='bar', legend=True, ax=ax1)
    ax1.set_title('Częstotliwość występowania aminokwasów w białku w procentach')

    dane_do_wykresu = dane_wykres(lanc_kodonow)
    danedo1=(dane_do_wykresu[0])
    danedo2=(dane_do_wykresu[1])
    danedo3=(dane_do_wykresu[2])
    danedo4=(dane_do_wykresu[3])
    danedo5=(dane_do_wykresu[4])
    dlugosc=len(danedo1)
    LISTA=[]
    for x in range(dlugosc):
        LISTA.append(x)
    dane1 = {'Aminokwasy': LISTA, 'Indeks hydrofobowy': danedo1}
    dane2 = {'Aminokwasy': LISTA, 'Surface accesibility': danedo2}
    dane3 = {'Aminokwasy': LISTA, 'Normalized flexibility parameters': danedo3}
    dane4 = {'Aminokwasy': LISTA, 'Janin Interior to surface transfer energy scale': danedo4}
    dane5 = {'Aminokwasy': LISTA, 'instability index by Grantham R.': danedo5}
    figure2 = plt.Figure(figsize=(12, 4), dpi=50)
    ax2 = figure2.add_subplot(2, 2, 1)
    print(dane1)
    print(dane2)
    print(dane3)
    print(dane4)
    print(dane5)
    bar2 = FigureCanvasTkAgg(figure2, okno)
    bar2.get_tk_widget().grid(row=2, column=0)
    df2 = pd.DataFrame(dane1)
    df2 = df2[['Aminokwasy', 'Indeks hydrofobowy']].groupby('Aminokwasy').sum()
    df2.plot(kind='line', legend=True, ax=ax2)
    df3 = pd.DataFrame(dane2)
    df3 = df3[['Aminokwasy', 'Surface accesibility']].groupby('Aminokwasy').sum()
    df3.plot(kind='line', legend=True, ax=ax2)
    df4 = pd.DataFrame(dane3)
    df4 = df4[['Aminokwasy', 'Normalized flexibility parameters']].groupby('Aminokwasy').sum()
    df4.plot(kind='line', legend=True, ax=ax2)
    df5 = pd.DataFrame(dane4)
    df5 = df5[['Aminokwasy', 'Janin Interior to surface transfer energy scale']].groupby('Aminokwasy').sum()
    df5.plot(kind='line', legend=True, ax=ax2)
    df6 = pd.DataFrame(dane5)
    df6 = df6[['Aminokwasy', 'instability index by Grantham R.']].groupby('Aminokwasy').sum()
    df6.plot(kind='line', legend=True, ax=ax2)
    ax2.set_title('bekawchuj')

    napis = CTkLabel(okno, text=wzor, width=550)
    returnButton = CTkButton(okno, text="Wróć", command=lambda: main.wczytaj_recznie(okno, lanc_powrotny), width=100)
    AtomyButton = CTkButton(okno, text="Ilość Atomów", command=lambda: WypiszAtomy(okno,wzor, lanc_kodonow, lanc_powrotny), width=100)
    napis.grid(row=0, column=0, sticky=W+E)
    returnButton.grid(row=0, column=1, sticky=W+E)
    AtomyButton.grid(row=1, column=1, sticky=W+E)
