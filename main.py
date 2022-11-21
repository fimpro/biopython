from Bio.Seq import Seq
import customtkinter as ctk
from PIL import ImageTk,Image
def rozklad_na_bialka(lanc): #funkcja przyjmująca łańcuch kodonów, a zwracająca kandydatów na białka w tablicy
    bialka=[]

    j=0
    czy_bialko=False
    for i in lanc: #for który każdy kodon podłacza do danego białka
        if (i == "*"): #kodony stop oznaczone są w biopythonie *
            czy_bialko = False
            j += 1
        if (czy_bialko):
            bialka[j] = bialka[j] + i
        if(i=="M"):
            bialka.append("M")
            czy_bialko=True
    return bialka
def interfejs(): #funkcja obsługująca główny interfejs programu
    okno = ctk.CTk()
    okno.minsize(height=500, width = 650)
    okno.title("Katolgenom")
    ctk.set_appearance_mode("System")
    ctk.set_default_color_theme("blue")
    napis=ctk.CTkLabel(okno, text="Wybierz Opcje")
    plikButton = ctk.CTkButton(okno, text="Dane z Pliku", command=wczytaj_z_pliku)
    recznieButton = ctk.CTkButton(okno, text="Dane Ręcznie", command=wczytaj_recznie)
    instrukcjeButton = ctk.CTkButton(okno, text="Instrukcja",command= instrukcja)
    opcjeButton = ctk.CTkButton(okno, text="Opcje", command=interfaceOpcje)
    imgHelisa=ImageTk.PhotoImage((Image.open("helisa.png")).resize((500,500)))
    obraz=ctk.CTkLabel(image=imgHelisa)
    obraz.grid(row=0, column=1, rowspan=10)
    plikButton.grid(row=1, column=0)
    recznieButton.grid(row=2, column=0)
    instrukcjeButton.grid(row=3, column=0)
    opcjeButton.grid(row=4, column=0)
    napis.grid(row=0, column=0)
    okno.mainloop()
def wczytaj_z_pliku(): #funkcja wczytująca dane z pliku
    print("wczytaj_z_pliku")
def wczytaj_recznie(): #funkcja wczytująca dane z "palca"
    print("wczytaj_recznie")
def interfaceOpcje(): #funkcja wczytująca interface opcji
    print("interfaceOpcje")
def instrukcja(): #funkcja wczytująca instruckje użytkowania
    print("instrukcja")
if __name__ == '__main__':
    interfejs()



