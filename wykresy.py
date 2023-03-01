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
    #suwak = slider()    w przyszłosci funkcja do ustawiania drugiej zmiennej tego gowna na dole
    analizuj = ProteinAnalysis(lanc)        #nie ma za uj pojecia co to sa te skale ale brzmi madrze
    dane.append(analizuj.protein_scale(ProtParamData.kd, 7, 1))         #indeks hydrofobowy
    dane.append(analizuj.protein_scale(ProtParamData.em, 7, 1))         #surface accessibility
    dane.append(analizuj.protein_scale(ProtParamData.Flex, 7, 1))       #Normalized flexibility parameters
    dane.append(analizuj.protein_scale(ProtParamData.ja, 7, 1))         #Janin Interior to surface transfer energy scale
    dane.append(analizuj.protein_scale(operacje_chemiczne.gra, 7, 1))   #instability index by Grantham R.
    dane.append(len(dane[0]))   #pomaga stworzyc oś X w wykresie
    return dane

def wykresy(okno, wzor, lancpowrotny, lanc_Kodonow):  # robocza funkcja do wykresów

    for widget in okno.winfo_children():
        widget.destroy()

    analizuj = ProteinAnalysis(lanc_Kodonow)  # tablica do analizowania bialek
    dane_kwasy = []
    dane_kwasy.append(analizuj.count_amino_acids())   #ilosc kazdego z kwasow
    dane_kwasy.append(analizuj.get_amino_acids_percent())   #procent kazdego z kwasow
    # pojedyncze wartości
    dw = []
    dw.append(analizuj.secondary_structure_fraction())
    # zwraca tablice z 3 wartościami, ktore zawieraja %: sheets, helixes, turns cokolwiek by to nie było XD
    dw.append(analizuj.isoelectric_point())
    dw.append(pH_bialka(lancpowrotny))
    if (analizuj.instability_index() <= 40):
        dw.append("Białko stabilne")
    elif (analizuj.instability_index() > 40):
        dw.append("Białko niestabilne")
    print(dane_wykres(lancpowrotny))

    napis = CTkLabel(okno, text="", width=100)    #dopisać wzór tam gdzie puste!!!
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
