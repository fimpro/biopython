from Bio.Seq import Seq
import time
import sys
import os
from customtkinter import *
from PIL import ImageTk,Image
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
import operacje_chemiczne
from tkinter import filedialog
mode="dark" #zmienna przechowująca tryb aplikacji
base_color="dark-blue"
class przycisk_dane():                #klasa, która tworzy przycisk ctk z zapamiętaniem danych. Używam do zapisu danych z fora w przycisku (rozwiązuje errora)
    def __init__(self,okno,zwiazek,lanc):
        self.przycisk = CTkButton(okno,text="Generuj",width=100,command=lambda:rysuj_interface(okno,zwiazek,lanc))
    def gridbutton(self,x,y):
        self.przycisk.grid(row=x, column=y)

def interfejs(okno): #funkcja obsługująca główny interfejs programu
    for widget in okno.winfo_children():
        widget.destroy()
    seq=Seq("AUGGAACGCGAACCCUAC")
    translated=seq.translate()
    mol = Chem.MolFromSmiles(operacje_chemiczne.wzor_lancucha_aminokwasow(translated))
    napis=CTkLabel(okno, text="Wybierz Opcje")
    plikButton = CTkButton(okno, text="Dane z Pliku", command=lambda:wczytaj_z_pliku(okno))
    recznieButton = CTkButton(okno, text="Dane Ręcznie", command=lambda:wczytaj_recznie(okno,lancpoczatkowy="AUGUAA"))
    instrukcjeButton = CTkButton(okno, text="Instrukcja",command= instrukcja)
    opcjeButton = CTkButton(okno, text="Opcje", command=lambda:interfaceOpcje(okno))
    img1=ImageTk.PhotoImage(Draw.MolToImage(mol, size=(500, 500), kekulize=True, wedgeBonds=True, fitImage=True))
    obraz=CTkLabel(okno,image=img1)
    obraz.grid(row=0, column=1, rowspan=10)
    plikButton.grid(row=1, column=0)
    recznieButton.grid(row=2, column=0)
    instrukcjeButton.grid(row=3, column=0)
    opcjeButton.grid(row=4, column=0)
    napis.grid(row=0, column=0)
    oknoR.mainloop()
def wczytaj_z_pliku(okno): #funkcja wczytująca dane z pliku
    for widget in okno.winfo_children():
        widget.destroy()
    napisTytul = CTkLabel(okno, text="wpisz lokalizacje pliku:")
    sciezkaWejscie= CTkEntry(okno, width=500)
    przyciskPrzegladaj = CTkButton(okno, text="przegladaj", command=lambda:przegladaj(sciezkaWejscie))
    przyciskSzukaj = CTkButton(okno, text="Szukaj", command=lambda: szukaj_z_pliku(okno))
    przyciskWroc = CTkButton(okno, text="Wróc", command=lambda: interfejs(okno))
    napisTytul.grid(row=0,column=0)
    sciezkaWejscie.grid(row=1,column=0,columnspan=5)
    przyciskPrzegladaj.grid(row=1,column=6)
    przyciskSzukaj.grid(row=2,column=6)
    przyciskWroc.grid(row=3,column=6)
def przegladaj(wejscie):
    filename = filedialog.askopenfilename(initialdir=os.getcwd(), title="Select a File", filetypes=(("Text files","*.txt*"), ("all files","*.*")))
    wejscie.insert(0,filename)
def szukaj_z_pliku():
    print("hendorzyc disa psa syna diabla")

def wczytaj_recznie(okno,lancpoczatkowy=""): #funkcja wczytująca dane z "palca"
    for widget in okno.winfo_children():
        widget.destroy()
    genom=Seq(lancpoczatkowy)
    wejscie = CTkEntry(okno, width=398)
    wejscie.insert(0, lancpoczatkowy)
    szukaj_interface(okno,wejscie)
    opis=CTkLabel(okno, text="Wpisz kod nici Rna",width=100)
    szukajButton = CTkButton(okno, text="Szukaj Białek", width=140,command=lambda:szukaj_interface(okno,wejscie))
    wrocButton = CTkButton(okno, text="Wroc", width=140,command=lambda:interfejs(okno))

    opis.grid(row=0,column=0)
    wejscie.grid(row=0,column=1, columnspan=4 )
    szukajButton.grid(row=0,column=5)
    wrocButton.grid(row=1,column=5)
    okno.mainloop()
def szukaj_interface(okno,wejscie): #funkcja wypisująca wszystkie łańcuchy w interfejsie
    lanc=Seq(wejscie.get())

    bialka=operacje_chemiczne.rozklad_na_bialka(lanc.translate())
    bialka_napisy=[]
    bialka_przyciski=[]
    for x in bialka:
        bialka_napisy.append(CTkLabel(okno,text=x,width=398))
        bialka_przyciski.append(przycisk_dane(okno,x,lanc))
    i=1
    for x in bialka_napisy:
        x.grid(row=i,column=0,columnspan=4)
        i+=1
    i=1
    for x in bialka_przyciski:
        x.gridbutton(i,4)
        i+=1
def rysuj_interface(okno,lanc_Kodonow,lancpowrotny):
    for widget in okno.winfo_children(): #czycimy okno
        widget.destroy()
    Smiles=operacje_chemiczne.wzor_lancucha_aminokwasow(lanc_Kodonow)
    mol = Chem.MolFromSmiles(Smiles)
    wzor_strukturalny = ImageTk.PhotoImage(Draw.MolToImage(mol, size=(500, 500), kekulize=True, wedgeBonds=True, fitImage=True))
    WrocButton = CTkButton(okno,text="wroc",command=lambda:wczytaj_recznie(okno,lancpowrotny))
    WykresyButton = CTkButton(okno,text="wykresy",command=lambda:wykresy(okno,lanc_Kodonow,lancpowrotny))
    obraz=CTkLabel(okno, image=wzor_strukturalny)
    obraz.grid(row=0,column=0,rowspan=6)
    WrocButton.grid(row=0,column=1)
    WykresyButton.grid(row=1,column=1)
    okno.mainloop()
def wykresy(okno,wzor,lancpowrotny): #robocza funkcja do wykresów

    for widget in okno.winfo_children():
        widget.destroy()
    napis=CTkLabel(okno,text=wzor, width=100)
    returnButton=CTkButton(okno, text="Wróć",command=lambda:rysuj_interface(okno,wzor,lancpowrotny))
    napis.grid(row=0, column=0)
    returnButton.grid(row=1, column=0)
    okno.mainloop()
def interfaceOpcje(okno): #funkcja wczytująca interface opcji
    for widget in okno.winfo_children():
        widget.destroy()
    napisMode=CTkLabel(okno,text="zmien kolor tla:")
    napisKolor = CTkLabel(okno, text="zmien kolor przycisków:")
    napisWroc = CTkLabel(okno, text="wróc do menu:")
    przyciskMode= CTkButton(okno,text=mode,command=lambda:zmien_tryb(okno))
    przyciskKolor = CTkButton(okno,text=base_color,command=lambda:zmien_kolor(okno))
    przyciskWroc= CTkButton(okno,text="Wróc",command=lambda:zapis(okno),border_color=okno.fg_color)
    przyciskMode.grid(row=1,column=0)
    przyciskKolor.grid(row=3,column=0)
    przyciskWroc.grid(row=5,column=0)
    napisMode.grid(row=0, column=0)
    napisKolor.grid(row=2, column=0)
    napisWroc.grid(row=4, column=0)
def zapis(okno):
    plik = open("opcje.txt", 'w')
    plik.write(mode+'\n')
    plik.write(base_color+'\n')
    plik.close()
    interfejs(okno)
def zmien_tryb(okno):
    global mode

    if(mode=="dark"):
        set_appearance_mode("light")
        mode="light"
        interfaceOpcje(okno)   #ta funkcja aktualizuje napisy na przyciskach
        return 0
    if(mode=="light"):
        set_appearance_mode("dark")
        mode = "dark"
        interfaceOpcje(okno) #ta funkcja aktualizuje napisy na przyciskach
        return 0
def zmien_kolor(okno):
    global base_color
    if (base_color == "dark-blue"):
        set_default_color_theme("blue")
        base_color="blue"
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
def instrukcja(): #funkcja wczytująca instruckje użytkowania
    print("instrukcja")
if __name__ == '__main__':
    oknoR = CTk()
    okno = CTkFrame(oknoR,fg_color=oknoR.fg_color,relief=RAISED)
    okno.pack()
    oknoR.minsize(height=500, width=650)
    oknoR.title("Katolgenom")
    plik = open("opcje.txt")
    opcje=plik.readlines()
    plik.close()
    mode=opcje[0].rstrip()
    base_color = opcje[1].rstrip()
    set_appearance_mode(mode)
    set_default_color_theme(base_color)
    interfejs(okno)
    oknoR.mainloop()

