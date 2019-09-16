###### Open all MACADAM results file and regroup all in a file ####
###### Filter on the 3 scores (Present in/PS/PFS) ####
import re
import os

dPWYtoSpecies={}

for file in os.listdir("./"):
    if file.endswith(".compact.tsv"):
        print (file)
        with open("./"+file,"r") as inputfile:
            lListofLines=inputfile.readlines()
            lListofIDNamePresentPSPFS=[]
            for line in lListofLines:
                if ", Taxonomy: " in line:
                    line=re.split(', Taxonomy: |\)',line)
                    sCurrentTax=line[1]
            print (sCurrentTax)
            for line in lListofLines:
                line=line.split("\t")
                if len(line)==12:
                    if line[0]!="Metabolic Pathway":
                        temp=[line[11].split("=")[-1].replace("\n",""),line[0],line[1],line[2],line[3]]
                        print (temp)
                        division=temp[2].split("/")
                        if (int(division[0])/int(division[1]))>=0: # Present in
                            if float(temp[3])>=0: # PS
                                if float(temp[4])>=0: # PFS
                                    if (temp[0],temp[1]) in dPWYtoSpecies:
                                        dPWYtoSpecies[(temp[0],temp[1])][sCurrentTax]=line[1],line[2],line[3]
                                    else:
                                        dPWYtoSpecies[(temp[0],temp[1])]={}
                                        dPWYtoSpecies[(temp[0],temp[1])][sCurrentTax]=line[1],line[2],line[3]

with open("./PresentPathwayMACADAM.tsv","w") as outputfile:
    for element in dPWYtoSpecies:
        outputfile.write(element[0]+"\t"+element[1]+"\n")

with open("./GroupResults.tsv","w") as outputfile:
    outputfile.write("#Pathway-ID\tPathway-Name\tIn-Species\n")
    for item in dPWYtoSpecies:
        if type(item) is tuple:
            lListofSpecies=[]
            for tax in dPWYtoSpecies[item]:
                print
                lListofSpecies.append(tax+" ("+dPWYtoSpecies[item][tax][0]+"; "+dPWYtoSpecies[item][tax][1]+"; "+dPWYtoSpecies[item][tax][2]+")")
            outputfile.write(item[0]+"\t"+item[1].strip()+"\t")
            outputfile.write(", ".join(lListofSpecies)+"\n")
