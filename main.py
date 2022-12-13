from Bio.Seq import Seq
import time
import sys
from customtkinter import *
from PIL import ImageTk,Image
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
import operacje_chemiczne
#zmienna przechowująca tryb aplikacji
mode="dark"
def interfejs(okno): #funkcja obsługująca główny interfejs programu
    for widget in okno.winfo_children():
        widget.destroy()
    seq=Seq("AUGGAACGCGAACCCUAC")
    translated=seq.translate()
    mol = Chem.MolFromSmiles(operacje_chemiczne.wzor_lancucha_aminokwasow(translated))
    napis=CTkLabel(okno, text="Wybierz Opcje")
    plikButton = CTkButton(okno, text="Dane z Pliku", command=wczytaj_z_pliku)
    recznieButton = CTkButton(okno, text="Dane Ręcznie", command=lambda:wczytaj_recznie(okno,wzorpoczatkowy="CACUAC"))
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
def wczytaj_z_pliku(): #funkcja wczytująca dane z pliku
    print("wczytaj_z_pliku")
def wczytaj_recznie(okno,wzorpoczatkowy=""): #funkcja wczytująca dane z "palca"
    for widget in okno.winfo_children():
        widget.destroy()
    genom=Seq(wzorpoczatkowy)
    mol = Chem.MolFromSmiles(operacje_chemiczne.wzor_lancucha_aminokwasow(genom.translate()))
    img = ImageTk.PhotoImage(Draw.MolToImage(mol, size=(500, 500), kekulize=True, wedgeBonds=True, fitImage=True))
    obraz = CTkLabel(okno, image=img)
    obraz.grid(row=1, column=0,rowspan=5, columnspan=5)
    opis=CTkLabel(okno, text="Wpisz kod nici Rna",width=100)
    wejscie = CTkEntry(okno, width=398)
    wejscie.insert(0, wzorpoczatkowy)
    generujButton = CTkButton(okno, text="Generuj", width=140,command=lambda:generuj(okno,wejscie))
    wykresyButton = CTkButton(okno, text="Wykresy", width=140,command=lambda:wykresy(okno,wejscie))
    wrocButton = CTkButton(okno, text="Wroc", width=140,command=lambda:interfejs(okno))

    opis.grid(row=0,column=0)
    wejscie.grid(row=0,column=1, columnspan=4 )
    generujButton.grid(row=0,column=5)
    wykresyButton.grid(row=1,column=5)
    wrocButton.grid(row=2,column=5)
    okno.mainloop()
def generuj(okno,wejscie): #funkcja generująca obrazek danego związku z ciągu aminokwasów
    t0 = time.time()
    genom=Seq(wejscie.get())
    kodon=genom.translate()
    Smiles=operacje_chemiczne.wzor_lancucha_aminokwasow(kodon)
    mol = Chem.MolFromSmiles(Smiles)
    img = ImageTk.PhotoImage(Draw.MolToImage(mol, size=(500, 500), kekulize=True, wedgeBonds=True, fitImage=True))
    obraz = CTkLabel(okno,image=img)
    obraz.grid(row=1, column=0, columnspan=5, rowspan=5)
    print (time.time() - t0)
    okno.mainloop()
def wykresy(okno,wejscie): #robocza funkcja do wykresów
    wzor=wejscie.get()
    for widget in okno.winfo_children():
        widget.destroy()
    napis=CTkLabel(okno,text=wzor, width=100)
    returnButton=CTkButton(okno, text="Wróć",command=lambda:wczytaj_recznie(okno,wzorpoczatkowy=wzor))
    napis.grid(row=0, column=0)
    returnButton.grid(row=1, column=0)
    okno.mainloop()
def interfaceOpcje(okno): #funkcja wczytująca interface opcji
    for widget in okno.winfo_children():
        widget.destroy()
    przyciskMode= CTkButton(okno,text=mode,command=lambda:zmien_tryb(okno))
    przyciskWroc= CTkButton(okno,text="Wróc",command=lambda:interfejs(okno))
    przyciskMode.grid(row=0,column=0)
    przyciskWroc.grid(row=1,column=0)

def zmien_tryb(okno):
    global mode
    plik = open("opcje.txt")
    opcje = plik.readlines()
    plik.close()

    if(mode=="dark"):
        set_appearance_mode("light")
        mode="light"
        opcje[0] ="light\n"
        plik = open("opcje.txt",'w')
        plik.write(opcje[0])
        plik.write(opcje[1])
        plik.close()
        return 0
    if(mode=="light"):
        set_appearance_mode("dark")
        mode = "dark"
        opcje[0] = "dark\n"
        plik = open("opcje.txt", 'w')
        plik.write(opcje[0])
        plik.write(opcje[1])
        plik.close()
        return 0


def instrukcja(): #funkcja wczytująca instruckje użytkowania
    print("instrukcja")
if __name__ == '__main__':
    oknoR = CTk()
    okno = CTkFrame(oknoR)
    okno.pack()
    oknoR.minsize(height=500, width=650)
    oknoR.title("Katolgenom")

    plik = open("opcje.txt")
    opcje=plik.readlines()
    set_appearance_mode(opcje[0].rstrip())
    set_default_color_theme(opcje[1].rstrip())
    interfejs(okno)



