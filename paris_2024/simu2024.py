import time
import matplotlib.pyplot as plt
import os
import numpy as np
from pyexcel_ods import get_data
import pandas as pd


ficLignes2 = "C:\\Users\\guezg\\Desktop\\tipe\\code_tipe\\paris_2024\\parametre\\lignes2024.ods"
ficReseau2 = "C:\\Users\\guezg\\Desktop\\tipe\\code_tipe\\paris_2024\\parametre\\reseau2024.ods"
ficLignes3 = "D:\\TIPE\\code_tipe\\paris_2024\\parametre\\lignes2024.ods"
ficReseau3 = "D:\\TIPE\\code_tipe\\paris_2024\\parametre\\reseau2024.ods"

def start_simulation(id_evenement):
        #données d'entrée:
    reseau = get_data(ficReseau3)
    lignes = get_data(ficLignes3)
    dt = reseau.get("duree")[1][1]
    NbrIt = reseau.get("duree")[2][1]
    name_stations,entree_annuelle, type_stations = reseau.get("stations")[1:4]
    rames = reseau.get("rames")[2]
    aire_quai_ligne = reseau.get("rames")[0:2]
    repartition_horaire = reseau.get("repartition_horaire")[1:5]
    NbrS = len(entree_annuelle)
    NbrL = len(rames)
    ens_lignes = [0 for i in range(NbrL)]
    ens_lignes[0] = lignes.get("1")
    ens_lignes[1] = lignes.get("2")
    ens_lignes[2] = lignes.get("3")
    ens_lignes[3] = lignes.get("4")
    ens_lignes[4] = lignes.get("5")
    ens_lignes[5] = lignes.get("6")
    ens_lignes[6] = lignes.get("7")
    ens_lignes[7] = lignes.get("8")
    ens_lignes[8] = lignes.get("9")
    ens_lignes[9] = lignes.get("10")
    ens_lignes[10] = lignes.get("11")
    ens_lignes[11] = lignes.get("12")
    ens_lignes[12] = lignes.get("13")
    ens_lignes[13] = lignes.get("14")
    ens_lignes[14] = lignes.get("15")
    ens_lignes[15] = lignes.get("16")
    ens_lignes[16] = lignes.get("17")
    ens_lignes[17] = lignes.get("3bis")
    ens_lignes[18] = lignes.get("7bis")
    correspondances = reseau.get("correspondances")
    S_eve = reseau.get("evenement")[0][id_evenement]
    H_eve = reseau.get("evenement")[1][id_evenement]
    entree_eve = reseau.get("evenement")[2][id_evenement]
        #applications sur les données:
    type_stations = conversion(type_stations)
    M_adja, temps_trajet = creer_matrice_adjacente_et_temps_trajet(ens_lignes,NbrS,NbrL)
    evenement = creer_evenement(S_eve,H_eve,entree_eve,M_adja,NbrS,NbrIt,dt)
    entree = creer_entree(entree_annuelle,repartition_horaire,type_stations,evenement,NbrS,NbrIt,dt)
    arcs = creer_arcs(ens_lignes,rames,type_stations,M_adja, temps_trajet,evenement,S_eve,H_eve,NbrS,NbrL,NbrIt,dt)
    aire_stations = creer_aire_stations(correspondances,aire_quai_ligne,NbrS,NbrL)
        #initialisation:
    donnees = []
    stations = np.zeros((NbrIt,NbrS),int)
    sortie = np.zeros((NbrIt,NbrS),int)
    densite = np.zeros((NbrIt,NbrS))
    couleur = np.zeros((NbrIt,NbrS),str)
    stations[0],sortie[0] = iteration(np.zeros(NbrS),entree,arcs,NbrS,NbrIt,0)
    densite[0] = creer_densite_it(stations[0],aire_stations,NbrS)
    couleur[0] = codage_couleur_it(densite[0],NbrS)
    for u in range(NbrS):
        donnees += [[0,u,stations[0][u],densite[0][u],entree[0][u],sortie[0][u]]]
        #tour:
    for it in range(1,NbrIt):
        stations[it],sortie[it] = iteration(stations[it-1],entree,arcs,NbrS,NbrIt,it)
        densite[it] = creer_densite_it(stations[it],aire_stations,NbrS)
        couleur[it] = codage_couleur_it(densite[it],NbrS)
        for u in range(NbrS):
            donnees += [[it,u,stations[it][u],densite[it][u],entree[it][u],sortie[it][u]]]
    #print(stations,entree,sortie)
        #resultats
    return  np.array(donnees),couleur,name_stations


def conversion(type_stations):
    for u in range(len(type_stations)):
        if type_stations[u] == "C":
            type_stations[u] = 3
        elif type_stations[u] == "G":
            type_stations[u] = 2
        elif type_stations[u] == "B":
            type_stations[u] = 1
        else :
            type_stations[u] = 0
    return type_stations

def creer_loi_entree(repartition_horaire,NbrIt,dt):
    loi_entree = np.zeros((NbrIt,4),float) #4 types de stations
    #On parcours les 4 types de stations
    for it in range(NbrIt):
        for ts in range(4):
            hc = (it*dt)//60 #l'heure qui correspond à l'IT
            if ((it*dt)%60)%5 == 0 :
                loi_entree[it][ts] = repartition_horaire[ts][hc]/(100*12)
            else :
                loi_entree[it][ts] = (repartition_horaire[ts][hc]/100)*dt/60
    #print(np.sum(loi_entree_temp))
    return loi_entree

def creer_entree(entree_annuelle,repartition_horaire,type_stations,evenement,NbrS,NbrIt,dt):
    print("debut entree")
    entree = np.zeros((NbrIt,NbrS),int)
    loi_entree = creer_loi_entree(repartition_horaire,NbrIt,dt)
    NbrT = 2000000/(330*24*int(1*60/dt)) #repartition uniforme sur l'ensemble de paris du nombre de touristes attendues pour les JO par iteration
    for it in range(NbrIt):
        for u in range(NbrS):
            #entree pour station 0 à 303 sur l'annee; par jour pour les autres
            if u < 304 :
                entree[it][u] = int(loi_entree[it][type_stations[u]]*entree_annuelle[u]/228) + evenement[it][u][1]
            else :
                entree[it][u] = int(loi_entree[it][type_stations[u]]*entree_annuelle[u]) + evenement[it][u][1]
            if type_stations[u] != 3 :
                entree[it][u] += NbrT
    return entree

def creer_matrice_adjacente_et_temps_trajet(ens_lignes,NbrS,NbrL):
    print("debut M adja")
    M = np.zeros((NbrS,NbrS),bool)
    temps_trajet = np.zeros((NbrS,NbrS),int)
    for li in range(1,NbrL+1):
        #disjonction pour les lignes 10 et 7bis
        if li == 10 :
            for x in ens_lignes[li-1]:
                M[x[0]][x[1]] = True
                temps_trajet[x[0]][x[1]] = int(x[2])
        elif li == NbrL :
            for x in ens_lignes[li-1]:
                M[x[0]][x[1]] = True
                temps_trajet[x[0]][x[1]] = int(x[2])
        else :
            for x in ens_lignes[li-1]:
                M[x[0]][x[1]] = True
                temps_trajet[x[0]][x[1]] = int(x[2])
                M[x[1]][x[0]] = True
                temps_trajet[x[1]][x[0]] = int(x[2])
    return M,temps_trajet

def creer_sphere_influence(S_eve,M_adja,NbrS):
    verif = np.zeros(NbrS,bool)
    verif[S_eve] = True
    sphere_influence = [[] for i in range(11)] #sphere rayon 10 (arbitraire)
    sphere_influence[0] = [S_eve]
    for r in range(1,11):
        #on cherche dans les voisins des stations d'un rayon ceux qui n'ont pas encore étaient comptés et on les ajoute au rayon suivant
        for u in sphere_influence[r-1]:
            for v in range(NbrS):
                if M_adja[u][v] and not verif[v]:
                    sphere_influence[r] += [v]
                    verif[v] = True
    return sphere_influence

def creer_evenement(S_eve,H_eve,entree_eve,M_adja,NbrS,NbrIt,dt):
    print("debut evenement")
    evenement = [[[0,0] for u in range(NbrS)] for it in range(NbrIt)]
    sphere_influence = creer_sphere_influence(S_eve,M_adja,NbrS)
    it_eve = int((H_eve-5)*60/dt) #nombre d'iterations avant le debut de l'évenement
    entree_rayon = entree_eve/10 #repartition uniforme des entrees entre les rayons de la sphere d'influence
    for r in range(len(sphere_influence)):
        N = len(sphere_influence[r]) #nombre de stations sur rayon r
        for u in sphere_influence[r]:
            entree = entree_rayon/N
            fact = 2
            for it in range(NbrIt):
                if r == 0:
                        evenement[it][u] = [10,0]
                else :
                    if it_eve-int(1*60/dt)+(10-r)*3 <= it < it_eve: #entre H_eve-1+h(r) et H_eve
                        entree = int(entree/fact)
                        fact = 2*fact
                        evenement[it][u] = [2*(10-r),entree] #poids de la station,entre de la station par iteration
    return evenement

def creer_capacite(ens_lignes,rames,NbrS,NbrL):
    capacite = np.zeros((NbrS,NbrS),int)
    for li in range(1,NbrL+1):
        #disjonction pour les lignes 10 et 7bis
        #possible de sommer la capacité de plusieurs rames, ex: M8 et M9
        if li == 10 :
            for x in ens_lignes[li-1]:
                capacite[x[0]][x[1]] += rames[li-1]
        elif li == NbrL :
            for x in ens_lignes[li-1]:
                capacite[x[0]][x[1]] += rames[li-1]
        else :
            for x in ens_lignes[li-1]:
                capacite[x[0]][x[1]] += rames[li-1]
                capacite[x[1]][x[0]] += rames[li-1]
    return capacite


def creer_arcs(ens_lignes,rames,type_stations,M_adja, temps_trajet,evenement,S_eve,H_eve,NbrS,NbrL,NbrIt,dt):
    print("debut arc")
    capacite = creer_capacite(ens_lignes,rames,NbrS,NbrL)
    poid_stations = [[3,2,5,10],[4,9,6,1]] #influence de chaque station le matin et l'après_midi sur les stations voisines
    poid_sorties = [[4,1,6,9],[4,10,5,1]]
    echanges  = [[[0 for v in range(NbrS)] for u in range(NbrS)] for it in range(NbrIt)]
    arcs = [[[] for u in range(NbrS)] for it in range(NbrIt)]
    it_eve = int(H_eve*60/dt) #heure qui correspond à l'IT
    for it in range(NbrIt):
        hc = (it*dt)//60
        for u in range(NbrS):
            if u == S_eve :
                if it_eve-int(1*60/dt) <= it < it_eve :
                    Smatin = poid_sorties[0][type_stations[u]]+10
                    Saprem = poid_sorties[1][type_stations[u]]+10
                    #poid de sortie de l'evenement
            else :
                Smatin = poid_sorties[0][type_stations[u]]
                Saprem = poid_sorties[1][type_stations[u]]
            for v in range(NbrS):
                if M_adja[u][v]:
                    if hc < 7 : #matin(<12h)
                        echanges[it][u][v] = poid_stations[0][type_stations[v]] + evenement[it][u][0]
                        Smatin += echanges[it][u][v]
                    else :
                        echanges[it][u][v] = poid_stations[1][type_stations[v]] + evenement[it][u][0]
                        Saprem += echanges[it][u][v]
            for v in range(NbrS):
                if M_adja[u][v]:
                    if hc < 7 :
                        echanges[it][u][v] = echanges[it][u][v]/Smatin
                    else :
                        echanges[it][u][v] = echanges[it][u][v]/Saprem
                    #pour respecter somme probabilité = 1
                    arcs[it][u] += [[capacite[u][v],echanges[it][u][v],v,0,temps_trajet[u][v],0]]
        print("arcs",it)
    return arcs

def iteration(stations_itm1,entree,arcs,NbrS,NbrIt,it):
    print("iteration: ", it)
    stations_it = np.zeros(NbrS,int)
    NbrSortie = np.array(stations_itm1)
    for u in range(NbrS):
        reste = 0
        for N in range(len(arcs[it][u])):
            x = arcs[it-1][u][N]
            transfert = int(x[1]*stations_itm1[u])
            #print(x)
            #print(transfert)
            NbrSortie[u] = NbrSortie[u]-transfert #ceux qui sortent sont ceux qui ne vont pas vers une autre station
            if int(x[3]) != 0 :
                if int(x[4]) == int(x[-1]) :
                    #on ajoutes aux stations voisines la partie qui partent vers ces stations
                    stations_it[x[2]] += x[3]
                else :
                    #temp de trajet +1min
                    arcs[it][u][N][-1] = x[-1] + 1
                    arcs[it][u][N][3] += x[3]
                reste += transfert
            else :
                if transfert > x[0] :
                    arcs[it][u][N][3] += x[0]
                    reste += transfert-x[0]
                    arcs[it][u][N][-1] = x[-1]-1 #malus pour metro plein
                else :
                    arcs[it][u][N][3] += x[3] + transfert
                    if transfert > 0.9*arcs[it][u][N][0]:
                        arcs[it][u][N][-1] = x[-1]-1 #malus pour metro presque plein
            if it == NbrIt-1 :
                NbrSortie[x[2]] += arcs[it][u][N][3] #comptabilise les passagers encore présent apres la fin du service
            #print(reste)
        stations_it[u] += entree[it][u] + reste #on ajoute les entrees et ceux qui ne partent pas
        #print("station", u, stations_it[u])
    return stations_it,NbrSortie

def creer_aire_stations(correspondances,aire_quai_ligne,NbrS,NbrL):
    aire_stations = np.zeros(NbrS)
    for u in range(NbrS):
        for li in correspondances[u]:
            for i in range(NbrL):
                if str(li) == str(aire_quai_ligne[0][i]) :
                    aire_stations[u] += 2*aire_quai_ligne[1][i]
                    #print("aire",aire_stations[u])
    return aire_stations

def creer_densite_it(stations_it,aire_stations,NbrS):
    densite_it = np.zeros(NbrS)
    for u in range(NbrS):
        densite_it[u] = stations_it[u]/(aire_stations[u]*0.59)
        # 0.59 = aire moyenne d'une personne en m2
    return densite_it

def codage_couleur_it(densite_it,NbrS):
    couleur_it = np.zeros(NbrS,str)
    for u in range(NbrS):
        if densite_it[u] < 0.5 :
            couleur_it[u] = "green"
        elif densite_it[u] < 0.9 :
            couleur_it[u] = "yellow"
        elif densite_it[u] < 1 :
            couleur_it[u] = "red"
        else :
            couleur_it[u] = "black"
    return couleur_it


def afficher_graphe(stations,densite,name_stations,couleur,NbrS,dt):
    print("début graphe")
    i=0
    while i == 0 :
        print("itération:")
        it = int(input())
        print("1 pour densite, 0 pour stations")
        k = int(input())
        if k == 1 :
            graphe_densite(densite[it],name_stations,couleur[it],NbrS)
        else :
            graphe_personne(stations[it],name_stations,couleur[it],NbrS)
        print("fin?")
        i = int(input())
    return


def graphe_personne(stations_it,name_stations,couleur_it,NbrS):
    x = [i for i in range(NbrS)]
    height = stations_it
    heightFixe = [1000 for u in range(NbrS)]
    labelFixe = ["" for u in range(NbrS)]
    plt.bar(x,heightFixe,color = 'None',tick_label = labelFixe)
    plt.bar(x,height,color = couleur_it,tick_label = name_stations)
    plt.show()
    plt.close()
    return

def graphe_densite(densite_it,name_stations,couleur_it,NbrS):
    x = [i for i in range(NbrS)]
    heightFixe = [1 for u in range(NbrS)]
    labelFixe = ["" for u in range(NbrS)]
    plt.bar(x,heightFixe,color = 'None',tick_label = labelFixe)
    plt.bar(x,densite_it,color = couleur_it,tick_label = name_stations)
    plt.show()
    plt.close()
    return


donnees,couleur,name_stations = start_simulation(1)

DF1 = pd.DataFrame(donnees)
DF2 = pd.DataFrame(couleur)
DF3 = pd.DataFrame(name_stations)

with pd.ExcelWriter("D:\\TIPE\\resultats2024_1.xlsx") as writer:
    DF1.to_excel(writer, sheet_name='données',header=["itération","id_stations","trafic","densité","entrée","sortie"],index=False)
    DF2.to_excel(writer, sheet_name='couleur',header=name_stations,index=[i for i in range(1200)])
    DF3.to_excel(writer, sheet_name='nom',header="nom stations",index=[i for i in range(1200)])