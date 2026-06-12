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
    # # log_fit -> valori alle freq. sulla curva fittata  #
    # # logps -> spettro di potenza in scala  log10       #
    # # dmax ->  max. dist. dalla retta di fit del ps     #
    # # idmax -> indice nell'array di dmax                #
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
    pf = pd.read_csv(str(file_path), usecols = ['Photon Flux [0.1-100 GeV](photons cm-2 s-1)'], dtype = str).squeeze() #str.lstrip() funziona solo sulle padnas series, mentre quello creato è un df.
    pf = pf.str.lstrip('<')                                                                                            #alternativa: pf['nome colonna'], estrae una series

    return pf


def Ft(pf):
    
    """
    Esegue: calcola la Ft dei dati cui è sottratta la media. Il valore per frequenza 0, quindi la media dei dati, è posto uguale a 0.0 per evitare errori di arrotondamento invece di essere sottratto
    Restituisce: array contenente la Ft corretta

    
    """
    
    Ft = fft.rfft(pf)
    Ft[0] = 0.0

    
    return Ft


def F_fit(x, m, q):
 
#    Funzione per il fit dello spettro di potenza 
    
    y = m * x + q

    return y

def fit(logps, frq, F_fit, p0):

    """
    Esegue: fit con una retta in scala loglog dello spettro di potenza
    Restituisce: dizionari con punti calcolati in corrispondenza delle frequenze sulla curva fittata e i parametri
    
    Input: 
        
        # F_fit -> funzione per il fit                   #
        # p0 -> parametri iniziali per il fit:           #
        #        [0]: pendenza; [1]: intercetta          #
    
    """
    x = np.log10(frq[1:]) #il valore per [0] è zero
    params, pcov = optimize.curve_fit(F_fit, x, logps, p0)
    y_fit = F_fit(x, params[0], params[1]) #una retta in scala loglog sarebbe del tipo log(y) = mlog(x) + q. 
       
    return y_fit, params



def iddmax(logps, fitted_data):

    """
    Esegue: trova il valore massimo della distanza dei punti del ps in scala logaritmica dalla retta di fit e il relativo indice
    Restituisce: dizionari con valore della distanza massima e l'indice sull'array del logaritmo dei valori. Per trovare l'indice originale serve aggiungere 1, poichè log(0) non è definito e il primo valore della Ft è 0 il log si calcola su [1:].

    """
    arr = logps - fitted_data
    dmax = np.max(arr)
    idmax = np.argmax(arr)
    
    return dmax, idmax







#############################################
# Creazione classe per la gestione dei dati #
#############################################


class Sorgente:

   
    p0 = [-1, -14]

    
    def  __init__(self, file_path_W, file_path_M, color):
        
       
       self.JdW = pd.read_csv(str(file_path_W), usecols = ['Julian Date']).to_numpy(dtype = int).flatten()
       self.JdM = pd.read_csv(str(file_path_M), usecols = ['Julian Date']).to_numpy(dtype = int).flatten()

       self.pf = {'W' : corr_pf(file_path_W).to_numpy(dtype = float), 'M' : corr_pf(file_path_M).to_numpy(dtype = float)}
       
       self.Ft = {'W' : Ft(self.pf['W']), 'M': Ft(self.pf['M'])}

       self.ps = {'W': np.abs(self.Ft['W'])**2, 'M': np.abs(self.Ft['M'])**2}

       self.logps = {'W': np.log10(self.ps['W'][1:]), 'M': np.log10(self.ps['M'][1:])}

       self.frq ={'W': fft.rfftfreq(len(self.pf['W']), d = 7), 'M' : fft.rfftfreq(len(self.pf['M']), d = 30)}

       self.log_fitW, self.paramsW = fit(self.logps['W'], self.frq['W'], F_fit,  self.p0)
       self.log_fitM, self.paramsM = fit(self.logps['M'], self.frq['M'], F_fit,  self.p0)
       self.log_fit = { 'W' : self.log_fitW, 'M' : self.log_fitM }
       self.params = { 'W' : self.paramsW, 'M' : self.paramsM }

      
       self.dmaxW, self.idmaxW = iddmax(self.logps['W'], self.log_fit['W'])
       self.dmaxM, self.idmaxM = iddmax(self.logps['M'], self.log_fit['M'])
       self.dmax = {'W':self.dmaxW, 'M':self.dmaxM}
       self.idmax = { 'W' :self.idmaxW, 'M': self.idmaxM}

       self.sign = {'W' : None, 'M' : None }

       self.hist = {'W': [], 'M' : [] }

       self.color = str(color)




def argument_parser():

    parser = argparse.ArgumentParser(description='Analisi e plot curve di luce. Sorgenti: J1104, J1256, J1555, J2253', usage ='python3 light_curves.py  -opzione')

    parser.add_argument('-dfs', '--DataFrameSignificatività',action='store_true',help='mostra i DataFrame con i valori delle frequenze, il picco e la significatività dello stesso; campionamento mensile e settimanale')
    parser.add_argument('-pfW', '--PhotonFluxW', action='store_true', help='mostra i grafici del flusso; campionamento settimanale')
    parser.add_argument('-pfM', '--PhotonFluxM', action='store_true', help='mostra i grafici del flusso; campionamento mensile')
    parser.add_argument('-psW', '--PowerSpectrumW', action='store_true',help='mostra i grafici dello spettro di potenza; campionamento settimanale')
    parser.add_argument('-psM', '--PowerSpectrumM', action='store_true',help='mostra i grafici dello spettro di potenza; campionamento mensile')
    parser.add_argument('-pfWM', '--PhotonFluxWM', action='store_true', help='mostra i grafici del flusso campionato settimanalmente e mensilmente sovrapposti per ogni sorgente')
    parser.add_argument('-logpsW', '--LogPowerSpectrumW', action='store_true',help='mostra i grafici dello spettro di potenza in scala logaritmica; campionamento settimanale')
    parser.add_argument('-logpsM', '--LogPowerSpectrumM', action='store_true',help='mostra i grafici dello spettro di potenza in scala logaritmica; campionamento mensile')
    parser.add_argument('-flogpsW', '--FitANDLogPowerSpectrumW', action='store_true',help='mostra i grafici dello spettro di potenza in scala logaritmica e la retta di fit; campionamento settimanale')
    parser.add_argument('-flogpsM', '--FitANDLogPowerSpectrumM', action='store_true',help='mostra i grafici dello spettro di potenza in scala logaritmica e la retta di fit; campionamento mensile')
    parser.add_argument('-histW', '--histogramW', action = 'store_true', help= "mostra istogramma delle distanze massime dal fit delle curve di luce sintetiche; campionamento settimanale") 
    parser.add_argument('-histM', '--histogramM', action = 'store_true', help= "mostra istogramma delle distanze massime dal fit delle curve di luce sintetiche; campionamento mensile") 

    return parser.parse_args()



S11 = Sorgente('sorgenti/4FGL_J1104.4+3812_weekly_2_20_2025.csv', 'sorgenti/4FGL_J1104.4+3812_monthly_2_20_2025.csv', 'blue' )
S12 = Sorgente('sorgenti/4FGL_J1256.1-0547_weekly_2_20_2025.csv', 'sorgenti/4FGL_J1256.1-0547_monthly_2_20_2025.csv', 'orange')
S15 = Sorgente('sorgenti/4FGL_J1555.7+1111_weekly_2_20_2025.csv', 'sorgenti/4FGL_J1555.7+1111_monthly_2_20_2025.csv', 'green')
S22 = Sorgente('sorgenti/4FGL_J2253.9+1609_weekly_2_20_2025.csv', 'sorgenti/4FGL_J2253.9+1609_monthly_2_20_2025.csv','purple')

Sorgenti = {'4FGL_J1104.4+3812':S11, '4FGL_J1256.1-0547':S12, '4FGL_J1555.7+1111':S15, '4FGL_J2253.9+1609':S22 }

args = argument_parser()


###########################
# Calcolo significatività #
###########################

#np.random.seed(67) #Fissandolo si dovrebbero ottenere sempre gli stessi risultati che altrimenti varierebbero leggermente   

W = 'W'
M = 'M'
for sorg in Sorgenti.values():
    it = 0
    counterW = 0
    counterM = 0
    while it < 1000: # dovendo saltare alcuni cicli per via degli zeri nel ps il numero di iterazioni necessarie potrebbe essere maggiore di 1000, quindi uso while poichè il ciclo for ne ha uno fisso
        
        slcW = np.copy(sorg.pf[W])
        slcM = np.copy(sorg.pf[M])
        np.random.shuffle(slcW)
        np.random.shuffle(slcM)
        
        slcFtW = Ft(slcW)
        slcFtM = Ft(slcM)
        
        maskM = slcFtM[1:] == 0.0
        maskW = slcFtW[1:] == 0.0
        if np.any(maskW) or np.any(maskM): #potrebbe capitare che un valore sia esattamente 0.0 e il log restituirebbe -inf. Va controllato a monte. La funzione any è utile perché restituisce un booleano
            continue
        else:
            slcpsW = np.abs(slcFtW)**2
            slcpsM = np.abs(slcFtM)**2

            slclogpsW = np.log10(slcpsW[1:])
            slclogpsM = np.log10(slcpsM[1:])
        
            slclog_fitW, slcparamsW  = fit(slclogpsW, sorg.frq[W], F_fit, sorg.p0)
            slclog_fitM, slcparamsM = fit(slclogpsM, sorg.frq[M], F_fit, sorg.p0)
            slcdmaxW, slcidmaxW = iddmax(slclogpsW, slclog_fitW)
            slcdmaxM, slcidmaxM = iddmax(slclogpsM, slclog_fitM)
            if slcdmaxW > sorg.dmaxW:
                counterW += 1
            if slcdmaxM > sorg.dmax[M]:
                counterM += 1
            sorg.hist[W].append(slcdmaxW)
            sorg.hist[M].append(slcdmaxM)
            it += 1
            sorg.sign[W] = counterW / 1000
            sorg.sign[M] = counterM / 1000

    
###############################
# Creazione grafici e tabelle #
###############################


if args.PhotonFluxW:
    fig1, pfgW = plt.subplots(2,2)
    
    fig1.suptitle('Curva di luce, campionamento settimanale')
    
    fig1.supxlabel('time (Julian date)')
    fig1.supylabel('Photon Flux [0.1-100 GeV](photons cm-2 s-1)')
    for key,j,k in zip(Sorgenti, (0,0,1,1),(0,1,0,1)):
        pfgW[j][k].plot(Sorgenti[key].JdW, Sorgenti[key].pf['W'], c = Sorgenti[key].color, label = key)
        pfgW[j][k].grid(True)
        pfgW[j][k].legend()




if args.PhotonFluxM:
    fig2, pfgM = plt.subplots(2,2)
    
    fig2.suptitle('Curve di luce, campionamento mensile')
    
    fig2.supxlabel('time (Julian date)')
    fig2.supylabel('Photon Flux [0.1-100 GeV](photons cm-2 s-1)')
    for key,j,k in zip(Sorgenti, (0,0,1,1),(0,1,0,1)):
        pfgM[j][k].plot(Sorgenti[key].JdM, Sorgenti[key].pf['M'], c = Sorgenti[key].color, label = key)
        pfgM[j][k].grid(True)
        pfgM[j][k].legend()



if args.PowerSpectrumW:
    fig3, psgW = plt.subplots(2,2)
    
    fig3.suptitle("Spettro di potenza, campionamento settimanale." + "\n" + " In rosso le frequenze evidenziate dall'analisi")
    
    fig3.supylabel('Power spectrum (photons cm-2 s-1)^2')
    fig3.supxlabel('frequencies (cycles/day)')
    for key,j,k in zip(Sorgenti, (0,0,1,1),(0,1,0,1)):
        psgW[j][k].bar(Sorgenti[key].frq['W'][Sorgenti[key].idmax['W'] + 1], Sorgenti[key].ps['W'][Sorgenti[key].idmax['W'] + 1], color = 'red', width = 0.0004)
        psgW[j][k].plot(Sorgenti[key].frq['W'], Sorgenti[key].ps['W'], color = Sorgenti[key].color, label = key)
        psgW[j][k].grid(True)
        psgW[j][k].legend()


if args.PowerSpectrumM:
    fig4, psgM = plt.subplots(2,2)
    
    fig4.suptitle("Spettro di potenza, campionamento mensile." + "\n" + " In rosso le frequenze evidenziate dall'analisi")
    
    fig4.supxlabel('frequencies (cycles/day)')
    fig4.supylabel('Power spectrum (photons cm-2 s-1)^2')
    for key,j,k in zip(Sorgenti, (0,0,1,1),(0,1,0,1)):
        psgM[j][k].bar(Sorgenti[key].frq['M'][Sorgenti[key].idmax['M'] + 1], Sorgenti[key].ps['M'][Sorgenti[key].idmax['M'] + 1], color = 'red', width = 0.0004)
        psgM[j][k].plot(Sorgenti[key].frq['M'], Sorgenti[key].ps['M'], color = Sorgenti[key].color, label = key)
        psgM[j][k].legend()
        psgM[j][k].grid(True)


if args.PhotonFluxWM:
    
    fig5, pfgWM = plt.subplots(2,2)
    fig5.suptitle('Curve di luce, campionamento mensile e settimanale')
    
    fig5.supxlabel('time (Julian date)')
    fig5.supylabel('Photon Flux [0.1-100 GeV](photons cm-2 s-1)')
    for key,j,k in zip(Sorgenti, (0,0,1,1),(0,1,0,1)):
        pfgWM[j][k].plot(Sorgenti[key].JdW, Sorgenti[key].pf['W'], label = key + ' W', color = Sorgenti[key].color)
        pfgWM[j][k].plot(Sorgenti[key].JdM, Sorgenti[key].pf['M'], label = key + ' M', color = 'red')
        pfgWM[j][k].grid(True)
        pfgWM[j][k].legend()


if args.LogPowerSpectrumW:
    fig7, logpsgW = plt.subplots(2,2)

    fig7.suptitle('Spettro di potenza, campionamento settimanale. (LogLog)')
    
    fig7.supxlabel('frequencies(cycles/day)')
    fig7.supylabel('Power spectrum (log)  ((photons cm-2 s-1)^2 ))')
    for key,j,k in zip(Sorgenti, (0,0,1,1),(0,1,0,1)):
        logpsgW[j][k].semilogx(Sorgenti[key].frq['W'][1:], Sorgenti[key].logps['W'], color = Sorgenti[key].color, label = key)
        logpsgW[j][k].legend()
        logpsgW[j][k].grid(True)
        

if args.LogPowerSpectrumM:
    fig8, logpsgM = plt.subplots(2,2)
    
    fig8.suptitle('Spettro di potenza, campionamento mensile. (LogLog)')
    
    fig8.supxlabel('frequencies (cycles/day) ')
    fig8.supylabel('Power spectrum (log)  ((photons cm-2 s-1)^2)')
    for key,j,k in zip(Sorgenti, (0,0,1,1),(0,1,0,1)):
        logpsgM[j][k].semilogx(Sorgenti[key].frq['M'][1:], Sorgenti[key].logps['M'], color = Sorgenti[key].color, label = key)
        logpsgM[j][k].legend()
        logpsgM[j][k].grid(True)
    
if args.FitANDLogPowerSpectrumW:
    fig9, flogpsgW = plt.subplots(2,2)
    
    fig9.suptitle('Spettro di potenza  & fit, campionamento settimanale. (LogLog)')
    
    fig9.supxlabel('frequencies (cycles/day)')
    fig9.supylabel('Power spectrum (log) ((photons cm-2 s-1)^2)')
    for key,j,k in zip(Sorgenti, (0,0,1,1),(0,1,0,1)):
        flogpsgW[j][k].semilogx(Sorgenti[key].frq['W'][1:], Sorgenti[key].logps['W'], color = Sorgenti[key].color, label = key + '\n' + 'm: ' + str(Sorgenti[key].params['W'][0]) + '\n' + 'q: ' +                                str(Sorgenti[key].params['W'][1]))
        flogpsgW[j][k].semilogx(Sorgenti[key].frq['W'][1:], Sorgenti[key].log_fit['W'], color = 'red')
        flogpsgW[j][k].semilogx(Sorgenti[key].frq['W'][Sorgenti[key].idmax['W'] + 1], Sorgenti[key].logps['W'][Sorgenti[key].idmax['W']], marker = '.', color = 'red') #deve essere traslato di 1 sulle frequenze poichè quelle
        #partono dallo [0] mentre gli array logaritmici partono da [1] perché non possono essere definiti per zero 
        flogpsgW[j][k].legend()
        flogpsgW[j][k].grid(True)
        
if args.FitANDLogPowerSpectrumM:
    fig10, flogpsgM = plt.subplots(2,2)
    
    fig10.suptitle('Spettro di potenza & fit, campionamento mensile')
    
    fig10.supxlabel('frequencies (cycles/day) ')
    fig10.supylabel('Power spectrum (log)  ((photons cm-2 s-1)^2 )')
    for key,j,k in zip(Sorgenti, (0,0,1,1),(0,1,0,1)):
        flogpsgM[j][k].semilogx(Sorgenti[key].frq['M'][1:], Sorgenti[key].logps['M'], color = Sorgenti[key].color, label = key + '\n' + 'm: ' + str(Sorgenti[key].params['M'][0]) + '\n' + 'q: ' +                                str(Sorgenti[key].params['M'][1]))
        flogpsgM[j][k].semilogx(Sorgenti[key].frq['M'][1:], Sorgenti[key].log_fit['M'], color = 'red')
        flogpsgM[j][k].semilogx(Sorgenti[key].frq['M'][Sorgenti[key].idmax['M']+1], Sorgenti[key].logps['M'][Sorgenti[key].idmax['M']], marker = '.', color = 'red')
        flogpsgM[j][k].legend()
        flogpsgM[j][k].grid(True)

if args.histogramW:
    fig11, histW = plt.subplots(2,2)

    fig11.suptitle('Isogrammi distanze massime dalla retta di fit SLC; campionamento settimanale')
    fig11.supxlabel('distanza dalla retta di fit')
    for key,j,k in zip(Sorgenti, (0,0,1,1),(0,1,0,1)):
        histW[j][k].hist(Sorgenti[key].hist['W'], bins = 100, color = Sorgenti[key].color)
        histW[j][k].bar(Sorgenti[key].dmax['W'], 7 , color = 'red', width = 0.003, label = 'distanza curva di luce originale')
        histW[j][k].legend()
        histW[j][k].grid(True)

if args.histogramM:
    fig12, histM = plt.subplots(2,2)

    fig12.suptitle('Isogrammi distanze massime dalla retta di fit SLC; campionamento mensile')
    fig12.supxlabel('distanza dalla retta di fit')
    for key,j,k in zip(Sorgenti, (0,0,1,1),(0,1,0,1)):
        histM[j][k].hist(Sorgenti[key].hist['M'], bins = 70, color = Sorgenti[key].color)
        histM[j][k].bar(Sorgenti[key].dmax['M'], 7 , color = 'red', width = 0.003, label = 'distanza curva di luce originale')
        histM[j][k].legend()
        histM[j][k].grid(True)

if args.DataFrameSignificatività:
    keys = ('W', 'M')
    datiW = { 'sorgenti' : [], 'frequenze[cicli/giorno]' : [], 'Picco di potenza((photons cm-2 s-1)^2 )': [], 'distanza dal fit': [], 'Indice': [], 'Significatività' : []}
    datiM = { 'sorgenti' : [], 'frequenze[cicli/giorno]' : [], 'Picco di potenza((photons cm-2 s-1)^2 )': [], 'distanza dal fit': [], 'Indice': [], 'Significatività' : []}
    dati = {'W' : datiW, 'M' : datiM}
    for key in Sorgenti:
        for j in keys:
            dati[j]['sorgenti'].append(key)
            dati[j]['frequenze[cicli/giorno]'].append(Sorgenti[key].frq[j][Sorgenti[key].idmax[j] + 1])
            dati[j]['Picco di potenza((photons cm-2 s-1)^2 )'].append(Sorgenti[key].ps[j][Sorgenti[key].idmax[j] + 1])
            dati[j]['distanza dal fit'].append(Sorgenti[key].dmax[j])
            dati[j]['Indice'].append(Sorgenti[key].idmax[j] + 1)
            dati[j]['Significatività'].append(1 - Sorgenti[key].sign[j])
        
    dfW = pd.DataFrame(dati['W'])
    dfM = pd.DataFrame(dati['M'])
    print('Analisi Settimanale:  ', '\n')
    print(dfW)
    print( '\n', '\n', 'Analisi mensile:   ', '\n')
    print(dfM)




plt.show()

