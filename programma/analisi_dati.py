import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse 
from scipy import fft
from scipy import optimize



    #######################################################
    #                                                     #
    #                  SIGLE USATE                        #
    #                                                     #
    # # W -> campionamenti settimanali                    #
    # # M -> campionamenti mensili                        #
    #                                                     # 
    #                                                     #
    # # Jd -> Julian Date                                 # 
    # # pf -> photon flux                                 #
    # # Ft -> Fourier tranform                            #
    # # ps -> power spectrum                              #
    # # frq -> frequencies                                #
    # # log_fit -> valori alle freq. sulla curva fittata #
    # # wps -> spettro di potenza sbiancato               #
    # # idx_fda -> indici frequenze da analizzare         #
    # # p0 -> parametri iniziali per il fit               #
    # # color -> colore per i grafici                     # 
    #                                                     #
    #######################################################


def corr_pf(file_path):

    """
    
    Esegue: Crea il df e lo pulisce dai simboli ' < ', il valore assegnato è il limite superiore
    Restituisce: Il df pulito dai simboli '<'
    
    Input:
    
    #file_path -> percorso relativo del file     

    """
    pf = pd.read_csv(str(file_path), usecols = ['Photon Flux [0.1-100 GeV](photons cm-2 s-1)'], dtype = str).squeeze() #il metodo str.lstrip() funziona solo sulle padnas series, mentre quello creato è un df.
    pf = pf.str.lstrip('<')                                                                                            #alternativa: pf['nome colonna'], estrae una series

    return pf


def corr_Ft(sorg):
    
    """
    Esegue: calcola la Ft dei dati e sottrae la media. Il valore per frequenza 0 è posto uguale a 0
    Restituisce: array contenente la Ft corretta

    Input:
    
    #sorg -> sorgente di cui effettuare l'analisi
    
    """
    
    FtW = fft.rfft(sorg.pfW - np.mean(sorg.pfW))
    FtW[0] = 0.0

    FtM = fft.rfft(sorg.pfM - np.mean(sorg.pfM))
    FtM[0] = 0.0

    return FtW, FtM


def F_fit(x, m, q):
 
#    Funzione per il fit dello spettro di potenza 
    
    y = m * x + q

    return y


def fit(sorg, F_fit, p0):

    """
    Esegue: fit in sacal loglog dello spettro di potenza
    Restituisce: punti calcolati in corrispondenza delle frequenze sulla curva fittata
    
    Input: 
        
        # sorg -> sorgente di cui effettuare l'analisi   #  
        # f -> funzione per il fit                       #
        # p0 -> parametri iniziali per il fit            #
    
    """
    xW = np.log10(sorg.frqW[1:])
    yW = np.log10(sorg.psW[1:])
    paramsW, pcovW = optimize.curve_fit(F_fit, xW, yW, p0)
    y_fitW = np.power(sorg.frqW[1:], paramsW[0]) * ( 10 ** (paramsW[1]) )

    xM = np.log10(sorg.frqM[1:])
    yM = np.log10(sorg.psM[1:])
    paramsM, pcovM = optimize.curve_fit(F_fit, xM, yM, p0)
    y_fitM = np.power(sorg.frqM[1:], paramsM[0]) * ( 10 ** (paramsM[1]) )

    return y_fitW, paramsW, y_fitM, paramsM



def pre_wht(sorg):

    """
    Esegue: Sbiancamento dello spettro di potenza
    Restituisce: 
    
    # yw: rapporto tra lo potenza della sorgente diviso per e i punti sulla curva di fit #
    
    Input: 
    
    # sorg -> sorgente di cui effettuare l'analisi  
   
    
    """ 
    ywW = sorg.psW[1:] / sorg.log_fitW 
    ywM = sorg.psM[1:] / sorg.log_fitM

    return ywW, ywM



def find_idx(sorg):

    """
    Esegue: cerca gli indici dei picchi da analizzare
    Restituisce: array contenente gli indici

    Input:
    
    #sorg -> sorgente di cui effettuare l'analisi

    """
    maskW = sorg.wpsW > sorg.wpsW_avg + 2 * sorg.wpsW_std
    idx_fdaW = np.array(np.nonzero(maskW)).flatten()

    maskM = sorg.wpsM > sorg.wpsM_avg + 2 * sorg.wpsM_std
    idx_fdaM = np.array(np.nonzero(maskM)).flatten()

    return idx_fdaW, idx_fdaM



#############################################
# Creazione classe per la gestione dei dati #
#############################################


class Sorgente:

   
    p0 = [-1, -14]

    
    def  __init__(self, file_path_W, file_path_M, color ):
        
       
       self.JdW = pd.read_csv(str(file_path_W), usecols = ['Julian Date']).to_numpy(dtype = int).flatten()
       self.JdM = pd.read_csv(str(file_path_M), usecols = ['Julian Date']).to_numpy(dtype = int).flatten()

       self.pfW = corr_pf(file_path_W).to_numpy(dtype = float)
       self.pfM = corr_pf(file_path_M).to_numpy(dtype = float)

       self.FtW, self.FtM = corr_Ft(self)

       self.psW = np.abs(self.FtW)**2 
       self.psM = np.abs(self.FtM)**2

       self.frqW = fft.rfftfreq(len(self.pfW), d = 7)
       self.frqM = fft.rfftfreq(len(self.pfM), d = 30)

       self.log_fitW,self.paramsW, self.log_fitM, self.paramsM = fit(self, F_fit, self.p0)

       self.wpsW, self.wpsM = pre_wht(self)

       self.wpsW_std = np.std(self.wpsW)
       self.wpsM_std = np.std(self.wpsM)

       self.wpsW_avg = np.average(self.wpsW)
       self.wpsM_avg = np.average(self.wpsM)

       self.idx_fdaW, self.idx_fdaM = find_idx(self)

       self.color = str(color)




def argument_parser():

    parser = argparse.ArgumentParser(description='Analisi e plot curve di luce. Sorgenti: J1104, J1256, J1555, J2253', usage ='python3 light_curves.py  -opzione')

    parser.add_argument('-pscsvW', '--PsToCsvW', action='store_true', help='Copia i dati di ogni sorgente in un file .csv e manda a schermo un riassunto; campionamento settimanale')
    parser.add_argument('-pscsvM', '--PsToCsvM', action='store_true', help='Copia i dati di ogni sorgente in un file .csv e manda a schermo un riassunto; campionamento mensile')
    parser.add_argument('-dfsW', '--DataFrameSignificativitàW', action='store_true',help='mostra i DataFrame con i valori delle frequenze, il picco e la significatività dello stesso; campionamento settimanale')
    parser.add_argument('-dfsM', '--DataFrameSignificativitàM', action='store_true',help='mostra i DataFrame con i valori delle frequenze, il picco e la significatività dello stesso; campionamento mensile')
    parser.add_argument('-dfsWM', '--DataFrameSignificativitàWM',action='store_true',help='mostra i DataFrame con i valori delle frequenze, il picco e la significatività dello stesso; campionamento mensile e settimanale')
    parser.add_argument('-pfW', '--PhotonFluxW', action='store_true', help='mostra i grafici del flusso; campionamento settimanale')
    parser.add_argument('-pfM', '--PhotonFluxM', action='store_true', help='mostra i grafici del flusso; campionamento mensile')
    parser.add_argument('-psW', '--PowerSpectrumW', action='store_true',help='mostra i grafici dello spettro di potenza; campionamento settimanale')
    parser.add_argument('-psM', '--PowerSpectrumM', action='store_true',help='mostra i grafici dello spettro di potenza; campionamento mensile')
    parser.add_argument('-pfWM', '--PhotonFluxWM', action='store_true', help='mostra i grafici del flusso campionato settimanalmente e mensilmente sovrapposti per ogni sorgente')
    parser.add_argument('-psWM', '--PowerSpectrumWM', action='store_true', help='mostra i grafici dello spettro di potenza per campionamenti mensili e settimanali sovrapposti per ogni sorgente')
    parser.add_argument('-logpsW', '--LogPowerSpectrumW', action='store_true',help='mostra i grafici dello spettro di potenza in scala logaritmica; campionamento settimanale')
    parser.add_argument('-logpsM', '--LogPowerSpectrumM', action='store_true',help='mostra i grafici dello spettro di potenza in scala logaritmica; campionamento mensile')
    parser.add_argument('-flogpsW', '--FitANDLogPowerSpectrumW', action='store_true',help='mostra i grafici dello spettro di potenza in scala logaritmica e la retta di fit; campionamento settimanale')
    parser.add_argument('-flogpsM', '--FitANDLogPowerSpectrumM', action='store_true',help='mostra i grafici dello spettro di potenza in scala logaritmica e la retta di fit; campionamento mensile')
  

    return parser.parse_args()





S11 = Sorgente('sorgenti/4FGL_J1104.4+3812_weekly_2_20_2025.csv', 'sorgenti/4FGL_J1104.4+3812_monthly_2_20_2025.csv', 'blue' )
S12 = Sorgente('sorgenti/4FGL_J1256.1-0547_weekly_2_20_2025.csv', 'sorgenti/4FGL_J1256.1-0547_monthly_2_20_2025.csv', 'orange')
S15 = Sorgente('sorgenti/4FGL_J1555.7+1111_weekly_2_20_2025.csv', 'sorgenti/4FGL_J1555.7+1111_monthly_2_20_2025.csv', 'green')
S22 = Sorgente('sorgenti/4FGL_J2253.9+1609_weekly_2_20_2025.csv', 'sorgenti/4FGL_J2253.9+1609_monthly_2_20_2025.csv','purple')

Sorgenti = {'J1104':S11, 'J1256':S12, 'J1555':S15, 'J2253':S22 }

args = argument_parser()


#######################
# Calcolo probabilità #
#######################

"""
Ad ogni sorgente, in un dizionario (dict_probsM, dict_probsW, key1) inzializzato con valori None, viene associato a sua volta un dizionario (dict_prob) contenente l'indice delle frequenze emerse dalla
 scrematura come chiave (key2) e un numero contenente il valore della probabilità di ottenere un picco maggiore di quello originale. Per ogni sorgente si esegue lo shuffle e si valuta per ogni frequenza se il picco 
dello spettro di potenza della curva di luce sintetica (slc) è maggiore o minore di quello della curva di luce originale. In caso affermativo il counter viene aumentato di uno. Alla fine del ciclo il 
dizionario che contiene le sorgenti viene aggioranto e come valore assegnato il dizionario 
contenente le frequenze e i conteggi divisi per il numero di shuffle  
"""

np.random.seed(67)


dict_probsM = { name : None for name in Sorgenti.keys()}

for key1 in dict_probsM:
    dict_prob = dict.fromkeys(Sorgenti[key1].idx_fdaM, 0)
    slc = Sorgenti[key1].pfM.copy() #python crea un riferimento in questo caso, con copy si impedsice di modificare i dati originali
    for i in range(1000):
        np.random.shuffle(slc)
        slcfft = fft.rfft(slc - np.mean(slc))
        slcps = np.abs(slcfft)**2
        
        for key2 in dict_prob:
            if Sorgenti[key1].psM[key2] < slcps[key2]:
                dict_prob[key2] += 1
        
    dict_probsM[key1] = {a: b / 1000 for a,b in dict_prob.items()} 


dict_probsW = {name : None for name in Sorgenti.keys()}
for key1 in dict_probsW:
    dict_prob = dict.fromkeys(Sorgenti[key1].idx_fdaW, 0)
    slc = Sorgenti[key1].pfW.copy()
    for i in range(1000):
        np.random.shuffle(slc)
        slcfft = fft.rfft(slc - np.mean(slc))
        slcps = np.abs(slcfft)**2
        
        for key2 in dict_prob:
            if Sorgenti[key1].psW[key2] < slcps[key2]:
                dict_prob[key2] += 1
        
    dict_probsW[key1] = {a: b / 1000 for a,b in dict_prob.items()} 



###############################
# Creazione grafici e tabelle #
###############################


if args.PhotonFluxW:
    fig1, pfgW = plt.subplots(2,2)
    
    fig1.suptitle('Curva di luce, campionamento settimanale')
    
    fig1.supxlabel('time (Julian date)')
    fig1.supylabel('Photon Flux [0.1-100 GeV](photons cm-2 s-1)')
    for key,j,k in zip(Sorgenti, (0,0,1,1),(0,1,0,1)):
        pfgW[j][k].plot(Sorgenti[key].JdW, Sorgenti[key].pfW, c = Sorgenti[key].color, label = key)
        pfgW[j][k].grid(True)
        pfgW[j][k].legend()

        dictpf = {'time (Julian Date)': Sorgenti[key].JdW, 'Photon Flux [0.1-100 GeV](photons cm-2 s-1)': Sorgenti[key].pfW}
        df = pd.DataFrame(dictpf)
        print(df.describe())


if args.PhotonFluxM:
    fig2, pfgM = plt.subplots(2,2)
    
    fig2.suptitle('Curve di luce, campionamento mensile')
    
    fig2.supxlabel('time (Julian date)')
    fig2.supylabel('Photon Flux [0.1-100 GeV](photons cm-2 s-1)')
    for key,j,k in zip(Sorgenti, (0,0,1,1),(0,1,0,1)):
        pfgM[j][k].plot(Sorgenti[key].JdM, Sorgenti[key].pfM, c = Sorgenti[key].color, label =key)
        pfgM[j][k].grid(True)
        pfgM[j][k].legend()

        dictpf = {'time (Julian Date)': Sorgenti[key].JdM, 'Photon Flux [0.1-100 GeV](photons cm-2 s-1)': Sorgenti[key].pfM}
        df = pd.DataFrame(dictpf)
        print(df.describe())


if args.PowerSpectrumW:
    fig3, psgW = plt.subplots(2,2)
    
    fig3.suptitle("Spettro di potenza, campionamento settimanale." + "\n" + " In rosso le frequenze evidenziate dall'analisi")
    
    fig3.supylabel('Power spectrum (photons cm-2 s-1)^2')
    fig3.supxlabel('frequencies (cycles/day)')
    for key,j,k in zip(Sorgenti, (0,0,1,1),(0,1,0,1)):
        psgW[j][k].bar(Sorgenti[key].frqW, Sorgenti[key].psW, color = Sorgenti[key].color, label = key, width = 0.0001)
        psgW[j][k].bar(Sorgenti[key].frqW[Sorgenti[key].idx_fdaW], Sorgenti[key].psW[Sorgenti[key].idx_fdaW], color = 'red', width = 0.0001)
        psgW[j][k].legend()


if args.PowerSpectrumM:
    fig4, psgM = plt.subplots(2,2)
    
    fig4.suptitle("Spettro di potenza, campionamento mensile." + "\n" + " In rosso le frequenze evidenziate dall'analisi")
    
    fig4.supxlabel('frequencies (cycles/day)')
    fig4.supylabel('Power spectrum (photons cm-2 s-1)^2')
    for key,j,k in zip(Sorgenti, (0,0,1,1),(0,1,0,1)):
        psgM[j][k].bar(Sorgenti[key].frqM, Sorgenti[key].psM, color = Sorgenti[key].color, label = key, width = 0.0001)
        psgM[j][k].bar(Sorgenti[key].frqM[Sorgenti[key].idx_fdaM], Sorgenti[key].psM[Sorgenti[key].idx_fdaM], color = 'red', width = 0.0001)
        psgM[j][k].legend()


if args.PhotonFluxWM:
    
    fig5, pfgWM = plt.subplots(2,2)
    fig5.suptitle('Curve di luce, campionamento mensile e settimanale')
    
    fig5.supxlabel('time (Julian date)')
    fig5.supylabel('Photon Flux [0.1-100 GeV](photons cm-2 s-1)')
    for key,j,k in zip(Sorgenti, (0,0,1,1),(0,1,0,1)):
        pfgWM[j][k].plot(Sorgenti[key].JdW, Sorgenti[key].pfW, label = key + 'W')
        pfgWM[j][k].plot(Sorgenti[key].JdM, Sorgenti[key].pfM, label = key + 'M')
        pfgWM[j][k].grid(True)
        pfgWM[j][k].legend()


if args.PowerSpectrumWM:
    fig6, pfgWM = plt.subplots(2,2)
    
    fig6.suptitle = ("Spettro di potenza; campionamento mensile e sttimanale." + "\n" + "In rosso le frequenze evidenziate dall'analisi")
    
    fig6.supxlabel('frequencies (cycles/day)')
    fig6.supylabel('Power spectrum ((photons cm-2 s-1)^2)')
    for key,j,k in zip(Sorgenti, (0,0,1,1),(0,1,0,1)):
        pfgWM[j][k].bar(Sorgenti[key].frqW, Sorgenti[key].psW, label = key + 'W', width = 0.0003)
        pfgWM[j][k].bar(Sorgenti[key].frqM, Sorgenti[key].psM, label = key + 'M', width = 0.0003)
        pfgWM[j][k].grid(True)
        pfgWM[j][k].legend()


if args.LogPowerSpectrumW:
    fig7, logpsgW = plt.subplots(2,2)

    fig7.suptitle('Spettro di potenza (loglog), campionamento settimanale. ')
    
    fig7.supxlabel('frequencies(cycles/day)')
    fig7.supylabel('Power spectrum  ((photons cm-2 s-1)^2 ))')
    for key,j,k in zip(Sorgenti, (0,0,1,1),(0,1,0,1)):
        logpsgW[j][k].loglog(Sorgenti[key].frqW[1:], Sorgenti[key].psW[1:], color = Sorgenti[key].color, label = i)
        logpsgW[j][k].legend()


if args.LogPowerSpectrumM:
    fig8, logpsgM = plt.subplots(2,2)
    
    fig8.suptitle('Spettro di potenza (loglog), campionamento mensile')
    
    fig8.supxlabel('frequencies (cycles/day) ')
    fig8.supylabel('Power spectrum  ((photons cm-2 s-1)^2)')
    for key,j,k in zip(Sorgenti, (0,0,1,1),(0,1,0,1)):
        logpsgM[j][k].loglog(Sorgenti[key].frqM[1:], Sorgenti[key].psM[1:], color = Sorgenti[key].color, label = i)
        logpsgM[j][k].legend()

    
if args.FitANDLogPowerSpectrumW:
    fig9, flogpsgW = plt.subplots(2,2)
    
    fig9.suptitle('Spettro di potenza (loglog) & fit, campionamento settimanale')
    
    fig9.supxlabel('frequencies (cycles/day)')
    fig9.supylabel('Power spectrum  ((photons cm-2 s-1)^2)')
    for key,j,k in zip(Sorgenti, (0,0,1,1),(0,1,0,1)):
        flogpsgW[j][k].loglog(Sorgenti[key].frqW, Sorgenti[key].psW, color = Sorgenti[key].color, label = key + '\n' + 'm: ' + str(Sorgenti[key].paramsW[0]) + '\n' + 'q: ' + str(Sorgenti[key].paramsW[1]))
        flogpsgW[j][k].loglog(Sorgenti[key].frqW[1:], Sorgenti[key].log_fitW, color = 'red')
        flogpsgW[j][k].set_label('pendenza: ' + str(Sorgenti[key].paramsW[0]))
        flogpsgW[j][k].legend()
        
        
if args.FitANDLogPowerSpectrumM:
    fig10, flogpsgM = plt.subplots(2,2)
    
    fig10.suptitle('Spettro di potenza (loglog) & fit, campionamento mensile')
    
    fig10.supxlabel('frequencies (cycles/day) ')
    fig10.supylabel('Power spectrum  ((photons cm-2 s-1)^2 )')
    for key,j,k in zip(Sorgenti, (0,0,1,1),(0,1,0,1)):
        flogpsgM[j][k].loglog(Sorgenti[key].frqM, Sorgenti[key].psM, color = Sorgenti[key].color, label = key + '\n' + 'm: ' + str(Sorgenti[key].paramsM[0]) + '\n' + 'q: ' + str(Sorgenti[key].paramsM[1]))
        flogpsgM[j][k].loglog(Sorgenti[key].frqM[1:], Sorgenti[key].log_fitM, color = 'red')
        flogpsgM[j][k].legend()


if args.DataFrameSignificativitàW:
    for key in Sorgenti:
        indx = Sorgenti[key].idx_fdaW
        d = {'frequenze [cicli/giorni]' : Sorgenti[key].frqW[indx], 'Picchi spettro di potenza ' : Sorgenti[key].psW[indx], 'Significatività':dict_probsW[key].values()}
        df = pd.DataFrame(d)
        print('Sorgente: ', key, ' analisi settimanale')
        print(df, '\n', '\n', '\n')


if args.DataFrameSignificativitàM:
    for key in Sorgenti:
        indx = Sorgenti[key].idx_fdaM
        d = {'frequenze [cicli/giorno]' : Sorgenti[key].frqM[indx], 'Picchi spettro di potenza' : Sorgenti[key].psM[indx], 'Significatività':dict_probsM[key].values()}
        df = pd.DataFrame(d)
        print('Sorgente: ', key, ' analisi mensile')
        print(df, '\n', '\n', '\n')

if args.DataFrameSignificativitàWM:
    for key in Sorgenti:
        indxW = Sorgenti[key].idx_fdaW
        indxM = Sorgenti[key].idx_fdaM
        dW = {'frequenze [cicli/giorno]' : Sorgenti[key].frqW[indxW], 'Picchi spettro di potenza ' : Sorgenti[key].psW[indxW], 'Significatività':dict_probsW[key].values()}
        dM = {'frequenze [cicli/giorno]' : Sorgenti[key].frqM[indxM], 'Picchi spettro di potenza' : Sorgenti[key].psM[indxM], 'Significatività':dict_probsM[key].values()}
        dfW = pd.DataFrame(dW)
        dfM = pd.DataFrame(dM)
        print('Sorgente: ', key, ' analisi mensile')
        print(dfM, '\n', '\n','analisi settimanale: ', '\n')
        print(dfW, '\n', '\n', '\n')
    

if args.PsToCsvW:
    for key in Sorgenti:
        dictW = {'Frequencies (cycles/day)' : Sorgenti[key].frqW, 'Power spectrum (photons cm-2 s-1)^2' : Sorgenti[key].psW}
        df = pd.DataFrame(dictW)
        df.to_csv('Ps' + key + 'W' + '.csv')
        print(key, '  mensile.' , '\n',df.describe(),'\n',)

if args.PsToCsvM:
    for key in Sorgenti:
        dictM = {'Frequencies (cycles/day)' : Sorgenti[key].frqM, 'Power spectrum (photons cm-2 s-1)^2' : Sorgenti[key].psM}
        df = pd.DataFrame(dictM)
        df.to_csv('Ps' + key + 'M' '.csv')
        print(key, '  mensile.' , '\n',df.describe, '\n')
        
plt.show()




    
