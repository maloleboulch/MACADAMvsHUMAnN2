dMACADAMPWYtoInfo={}
with open("./MACADAM/GroupResults.tsv","r") as inputfile:
    lListofLine=inputfile.readlines()[1:]
    for line in lListofLine:
        line=line.replace("\n","")
        line=line.split("\t")
        SpeciesScore=line[2].split(", ")
        dTemp={}
        for element in SpeciesScore:
            element=element.split(" (")
            element[1]=element[1].replace(")","")
            element[1]=element[1].split("; ")
            dTemp[element[0]]=element[1]
        dMACADAMPWYtoInfo[line[0]]=[line[1],dTemp]


dHUMANNPWYtoInfo={}
with open("Humann/SpeciesperPathway.tsv","r") as inputfile:
    lListofLine=inputfile.readlines()[1:]
    for line in lListofLine:
        line=line.replace("\n","")
        line=line.split("\t")
        SpeciesScore=line[2].split(", ")
        dTemp={}
        for element in SpeciesScore:
            element=element.split(" (")
            element[1]=element[1].replace(")","")
            element[1]=element[1].split("; ")
            dTemp[element[0]]=element[1]
        dHUMANNPWYtoInfo[line[0]]=[line[1],dTemp]


setofMACADAMOnly=set()
setofHumannOnly=set()
setofSharedPathway=set()

for item in dMACADAMPWYtoInfo:
    if item in dHUMANNPWYtoInfo:
        setofSharedPathway.add(item)
    else:
        setofMACADAMOnly.add(item)

for item in dHUMANNPWYtoInfo:
    if item in dMACADAMPWYtoInfo:
        pass
    else:
        setofHumannOnly.add(item)

with open("./HumannOnlyPathway.tsv","w") as outputfile:
    outputfile.write("#PWY-ID\tName\tNumber of Species\tIn-Species(Pathway Coverage; Pathway Abundance)\n")
    for item in setofHumannOnly:
        all=False
        # print (dHUMANNPWYtoInfo[item][1])
        lineofspecies=[]
        for element in dHUMANNPWYtoInfo[item][1]:
            lineofspecies.append(element+" ("+dHUMANNPWYtoInfo[item][1][element][0]+"; "+dHUMANNPWYtoInfo[item][1][element][1]+")")
        for element in lineofspecies:
            if "All (" in element or "unclassified (" in element:
                all=True
        if all:
            outputfile.write(item+"\t"+dHUMANNPWYtoInfo[item][0]+"\t"+str(len(lineofspecies)-1)+"\t"+", ".join(lineofspecies)+"\n")
        else:
            outputfile.write(item+"\t"+dHUMANNPWYtoInfo[item][0]+"\t"+str(len(lineofspecies))+"\t"+", ".join(lineofspecies)+"\n")

with open("./MACADAMOnlyPathway.tsv","w") as outputfile:
    outputfile.write("#PWY-ID\tName\tNumber of Species\tIn-Species(Present in, PS, PFS)\n")
    for item in setofMACADAMOnly:
        lineofspecies=[]
        for element in dMACADAMPWYtoInfo[item][1]:
            lineofspecies.append(element+" ("+dMACADAMPWYtoInfo[item][1][element][0]+"; "+dMACADAMPWYtoInfo[item][1][element][1]+"; "+dMACADAMPWYtoInfo[item][1][element][2]+")")
        outputfile.write(item+"\t"+dMACADAMPWYtoInfo[item][0]+"\t"+str(len(lineofspecies))+"\t"+", ".join(lineofspecies)+"\n")

with open("./SharedPathway.tsv","w") as outputfile:
    outputfile.write("#PWY-ID\tName\tShared?\tNumber of Shared Species\tSpeciesShared\tNumber of HUMANN species\tIn-Species-HUMANN\tNumber of MACADAM species\tIn-Species-MACADAM\n")
    for item in setofSharedPathway:
        Shared=False
        all=False
        SpeciesInHuman=dHUMANNPWYtoInfo[item][1]
        SpeciesInMACADAM=dMACADAMPWYtoInfo[item][1]
        SharedSpecies=[]
        for element in SpeciesInHuman:
            if element in SpeciesInMACADAM:
                Shared=True
                SharedSpecies.append(element)
        lofSharedSpeciesInfos=[]
        if SharedSpecies!=[]:
            for element in SharedSpecies:
                lofSharedSpeciesInfos.append(element+" ("+" ; ".join(dHUMANNPWYtoInfo[item][1][element])+" ; "+" ; ".join(dMACADAMPWYtoInfo[item][1][element])+")")
        lofHumannSpeciesInfos=[]
        for element in dHUMANNPWYtoInfo[item][1]:
            if element not in SharedSpecies:
                lofHumannSpeciesInfos.append(element+" ("+" ; ".join(dHUMANNPWYtoInfo[item][1][element])+")")
        lofMACADAMSpeciesInfos=[]
        for element in dMACADAMPWYtoInfo[item][1]:
            if element not in SharedSpecies:
                lofMACADAMSpeciesInfos.append(element+" ("+" ; ".join(dMACADAMPWYtoInfo[item][1][element])+")")
        for element in lofHumannSpeciesInfos:
            if "All (" in element:
                all=True
        if Shared:
            if all:
                outputfile.write(item+"\t"+dMACADAMPWYtoInfo[item][0]+"\tYes\t"+str(len(lofSharedSpeciesInfos))+"\t"+", ".join(lofSharedSpeciesInfos)+"\t"+str(len(lofHumannSpeciesInfos)-1)+"\t"+", ".join(lofHumannSpeciesInfos)+"\t"+str(len(lofMACADAMSpeciesInfos))+"\t"+", ".join(lofMACADAMSpeciesInfos)+"\n")
            else:
                outputfile.write(item+"\t"+dMACADAMPWYtoInfo[item][0]+"\tYes\t"+str(len(lofSharedSpeciesInfos))+"\t"+", ".join(lofSharedSpeciesInfos)+"\t"+str(len(lofHumannSpeciesInfos))+"\t"+", ".join(lofHumannSpeciesInfos)+"\t"+str(len(lofMACADAMSpeciesInfos))+"\t"+", ".join(lofMACADAMSpeciesInfos)+"\n")
        else:
            if all:
                outputfile.write(item+"\t"+dMACADAMPWYtoInfo[item][0]+"\tNo\t0\tNA\t"+str(len(lofHumannSpeciesInfos)-1)+"\t"+", ".join(lofHumannSpeciesInfos)+"\t"+str(len(lofMACADAMSpeciesInfos))+"\t"+", ".join(lofMACADAMSpeciesInfos)+"\n")
            else:
                outputfile.write(item+"\t"+dMACADAMPWYtoInfo[item][0]+"\tNo\t0\tNA\t"+str(len(lofHumannSpeciesInfos))+"\t"+", ".join(lofHumannSpeciesInfos)+"\t"+str(len(lofMACADAMSpeciesInfos))+"\t"+", ".join(lofMACADAMSpeciesInfos)+"\n")
