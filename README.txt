Please use the "Download ZIP" button of github page
---------------------------------------------------

to compile:
make

to run:
depms FSP_INSTANCE_FILE CSV_OUTPUT_FILE --heu NEH_FILE
deptft FSP_INSTANCE_FILE CSV_OUTPUT_FILE --heu LR_FILE

FSP_INSTANCE_FILE, NEH_FILE, LR_FILE are contained in the directory "test_problems"

----------------------------------------------------

the directories "lr" and "neh" contain the software to generate LR_FILE and NEH_FILE given a FSP_INSTANCE_FILE

----------------------------------------------------

DEFINES PER COMPILAZIONE:

1) ONLINEPRINT: attiva stampa durante l'esecuzione
2) MYDEBUG: attiva stampa/stop per debug
3) GFC: attiva l'uso della struttura dati GFC per il calcolo fitness veloce in local search

NOTA: Si attivano aggiungendo "-DONLINEPRINT" (e simili) alla fine della riga "CPPFLAGS" nel "Makefile"

-------------------------------------------------------------------------------------------------------

MUTAZIONI DIFFERENTI:

1) dep.cpp.ORIGINAL
contiene il codice per la differential mutation originale dell'articolo

2) dep.cpp.PATHRELINKING 
contiene il codice per sostituire la diff. mutation con il path relinking

3) dep.cpp.RANDOMINDIVIDUAL 
contiene il codice per sostituire la diff. mutation con un individuo random dalla popolazione corrente

4) dep.cpp.COMPLETELYRANDOM
contiene il codice per sostituire la diff. mutation con un nuovo individuo completamente random

NOTA: Ognuna delle 4 diverse mutazioni va attivata copiando il relativo file su dep.cpp prima della compilazione
