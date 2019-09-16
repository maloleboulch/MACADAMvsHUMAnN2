#### Ce script permet de relier chaque pathway à une espèce.
import re


#Load all coverage data in a dict
ListoflinesCoverage=[]
with open("./SRR2081071_1_pathcoverage.tsv","r") as inputfile:
    for line in inputfile.readlines()[1:]:
        line=line.replace("\n","")
        line=line.split("\t")
        ListoflinesCoverage.append(line)

dPWYtoSpecies={}
lNopathway=["UNINTEGRATED","UNMAPPED"]


for item in ListoflinesCoverage:
    print (item)
    if float(item[1])>=0:
        if "|" in item[0]:
            item[0]=re.split(':|\|',item[0])
            if (item[0][0]) in lNopathway:
                if item[0][0] in dPWYtoSpecies:
                    dPWYtoSpecies[item[0][0]].append([item[0][1],item[1]])
                else:
                    dPWYtoSpecies[item[0][0]]=[item[0][1],item[1]]
            else:
                if (item[0][0],item[0][1]) in dPWYtoSpecies:
                    dPWYtoSpecies[(item[0][0],item[0][1])].append([item[0][-1],item[1]])
                else:
                    dPWYtoSpecies[(item[0][0],item[0][1])]=[[item[0][-1],item[1]]]
        else:
            item[0]=item[0].split(":")
            if item[0][0] in lNopathway:
                if item[0][0] in dPWYtoSpecies:
                    dPWYtoSpecies[item[0][0]].append(["All",item[1]])
                else:
                    dPWYtoSpecies[item[0][0]]=[["All",item[1]]]
            else:
                if (item[0][0],item[0][1]) in dPWYtoSpecies:
                    dPWYtoSpecies[(item[0][0],item[0][1])].append(["All",item[1]])
                else:
                    dPWYtoSpecies[(item[0][0],item[0][1])]=[["All",item[1]]]

#Load all abundance data in a dict

ListoflinesAbundance=[]
with open("./SRR2081071_1_pathabundance.tsv","r") as inputfile:
    for line in inputfile.readlines()[1:]:
        line=line.replace("\n","")
        line=line.split("\t")
        ListoflinesAbundance.append(line)

dPWYtoSpeciesAbundance={}
lNopathway=["UNINTEGRATED","UNMAPPED"]


for item in ListoflinesAbundance:
    if float(item[1])>=0:
        if "|" in item[0]:
            item[0]=re.split(':|\|',item[0])
            if (item[0][0]) in lNopathway:
                if item[0][0] in dPWYtoSpeciesAbundance:
                    dPWYtoSpeciesAbundance[item[0][0]].append([item[0][1],item[1]])
                else:
                    dPWYtoSpeciesAbundance[item[0][0]]=[item[0][1],item[1]]
            else:
                if (item[0][0],item[0][1]) in dPWYtoSpeciesAbundance:
                    dPWYtoSpeciesAbundance[(item[0][0],item[0][1])].append([item[0][-1],item[1]])
                else:
                    dPWYtoSpeciesAbundance[(item[0][0],item[0][1])]=[[item[0][-1],item[1]]]
        else:
            item[0]=item[0].split(":")
            if item[0][0] in lNopathway:
                if item[0][0] in dPWYtoSpeciesAbundance:
                    dPWYtoSpeciesAbundance[item[0][0]].append(["All",item[1]])
                else:
                    dPWYtoSpeciesAbundance[item[0][0]]=[["All",item[1]]]
            else:
                if (item[0][0],item[0][1]) in dPWYtoSpeciesAbundance:
                    dPWYtoSpeciesAbundance[(item[0][0],item[0][1])].append(["All",item[1]])
                else:
                    dPWYtoSpeciesAbundance[(item[0][0],item[0][1])]=[["All",item[1]]]

#Mix the dict
dTemp={}
for item in dPWYtoSpecies:
    dTemp[item]={}
    if item in dPWYtoSpeciesAbundance:
        for element in dPWYtoSpecies[item]:
            dTemp[item][element[0]]=[element[1]]
        for element in dPWYtoSpeciesAbundance[item]:
            if element[0] in dTemp[item]:
                dTemp[item][element[0]].append(element[1])

dPWYtoSpecies=dTemp

#### Print all existing pathway (ID+Name)
with open("./PresentPathway.tsv","w") as outputfile:
    for item in dPWYtoSpecies:
        if type(item) is tuple:
            outputfile.write(item[0]+"\t"+item[1].strip()+"\n")

setofPresentSpecies=set()
for item in dPWYtoSpecies:
    for element in dPWYtoSpecies[item]:
        setofPresentSpecies.add(element)

with open("./PresentSpecies.tsv","w") as outputfile:
    outputfile.write("##Present species: \n")
    for item in setofPresentSpecies:
        outputfile.write(item+"\n")

with open("./SpeciesperPathway.tsv","w") as outputfile:
    outputfile.write("#Pathway-ID\tPathway-Name\tIn-Species(Coverage,Abundance)\n")
    for item in dPWYtoSpecies:
        if type(item) is tuple:
            lListofSpecies=[]
            for tax in dPWYtoSpecies[item]:
                species=tax.split("s__")[-1]
                lListofSpecies.append(species.replace("_"," ")+" ("+dPWYtoSpecies[item][tax][0]+"; "+dPWYtoSpecies[item][tax][1]+")")
            outputfile.write(item[0]+"\t"+item[1].strip()+"\t")
            outputfile.write(", ".join(lListofSpecies)+"\n")
