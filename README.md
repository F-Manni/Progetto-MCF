# Progetto-MCF
Progetto per l'esame di metodi computazionali per la fisica.\
Tutte le analisi sono eseguite sia per i dati raccolti su base mensile che settimanale.\
Lo script python analizza le curve di luce delle quattro sorgenti:\
-Calcola la trasformata di Fourier ai dati del flusso di fotoni cui è stata prima sottratta la media per annullare la componente continua.\
-Calcola le frequenze.\
-Calcola lo spettro di potenza.\
-Elimina il trend dello spettro di potenza delle sorgenti: esegue il fit dello spettro in scala loglog con una retta; ritrasforma la retta in scala lineare; calcola i valori sulla curva per le frequenze trovate; divide il valore dello spettro originale per i valori calcolati sulla curva ottimizzata.\
-Trova media e deviazione standard dell'insieme calcolato al punto sopra.\
-Individua gli indici dei picchi nello spettro di potenza, trattato come descritto, che sono maggiori della somma tra la maedia dello spettro corretto e due deviazioni standard dello stesso.\
-Calcola per ognuni picco individuato la probabilità che generando curve di luce sintetiche, in corrispondenza delle loro frequenze, si crei un picco maggiore di quello presente nei dati originali.\

Il programma può mandare a schermo vari risultati sotto forma di grafici o DataFrame e per altri creare file .csv che li contengono.\
La gestione di queste funzionalità avviene tramite riga di comando da terminale.\
Tutte le funzionalità possono essere visualizzate tramite comando  "--help".\
Di seguito un elenco.\
Manda a schermo i seguenti grafici:\
	-Flusso di fotoni\
	-Spettro di potenza (Scala lineare o logaritmica)\
	-Spettro di potenza loglog con retta di fit sovrapposta\
Manda a schermo i DataFrame contenenti le informazioni sui picchi selezionati per l'analisi, il valore dello spettro originale con cui sono confrontati, la probabilità trovata.\
Crea un file .csv per ogni sorgente contenente i valori del suo spettro di potenza e le relative frequenze.
