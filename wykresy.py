from customtkinter import *
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import operacje_chemiczne
import main
def WypiszAtomy(okno, wzor,lanc_powrotny, lanc_kodonow):
    for widget in okno.winfo_children():
        widget.destroy()
    operacje_chemiczne.listaAtomow(wzor)
    returnButton = CTkButton(okno, text="Wróć", command=lambda: wykresy(okno,wzor, lanc_kodonow, lanc_powrotny))
    returnButton.grid(row=1, column=0)
def wykresy(okno, wzor, lanc_kodonow, lanc_powrotny):  # robocza funkcja do wykresów

    for widget in okno.winfo_children():
        widget.destroy()

    analizuj = ProteinAnalysis(lanc_kodonow)  # tablica do analizowania bialek
    dw_y = []  # dane wykres w osi y; tablica 2d
    dw_y.append(analizuj.instability_index())
    dw_y.append(analizuj.isoelectric_point())
    dw_y.append(analizuj.get_amino_acids_percent())
    indeksh = operacje_chemiczne.indeks_hydrofobowy(lanc_kodonow)
    dw_y.append(indeksh)
    print(lanc_powrotny)
    print(dw_y)

    napis = CTkLabel(okno, text=wzor, width=100)
    returnButton = CTkButton(okno, text="Wróć", command=lambda: main.wczytaj_recznie(okno, lanc_powrotny))
    AtomyButton = CTkButton(okno, text="Ilość Atomów", command=lambda: WypiszAtomy(okno,wzor, lanc_kodonow, lanc_powrotny))
    napis.grid(row=0, column=0)
    returnButton.grid(row=1, column=0)
    AtomyButton.grid(row=2, column=0)