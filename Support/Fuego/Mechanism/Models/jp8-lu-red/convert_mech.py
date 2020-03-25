import sys

# Skeletal chemkin mechanism ready for Fuego conversion, containing TRANS
infile = sys.argv[1]

# Output file for reduced mechanis, chemkin format but containing no reactions
outfile = sys.argv[2] # File for r

# QSS species removed from Skeletal in Chemkin
QSS = ['CH2', 'CH2*', 'HCO', 'CH3O', 'C2H3', 'C2H5', 'HCCO', 'CH2CHO', 'C6H5CO']

endline = '\n'

# If the inpiut string starts with any of the strings in the list, return that string
def startswithany(astring, stringlist):
    outstring = None
    for startstring in stringlist:
        if len(astring.split()) > 0:
            if astring.split()[0] == startstring:
                outstring = startstring
    return outstring

with open(infile,'r') as indat, open(outfile,'w') as outdat:

    sections = ["ELEMENTS","SPECIES","TRANS","REACTIONS"]
    insection = None
    seclines = dict(zip(sections, [0] * len(sections) ))
    specieslist = []
    
    for line in indat:
        if line.startswith("END"):
            # Get out of section and write the entire section of lines
            for secline in seclines[insection]:
                outdat.write(secline)
            outdat.write(line)
            insection = None
            
        elif insection == "SPECIES":
            # Get rid of unwanted species
            shortlist = line.split()
            for removespec in QSS:
                if removespec in shortlist:
                    shortlist.remove(removespec)
            specieslist += shortlist
            seclines[insection] =[]
            newline = ''
            for ii,spec in enumerate(specieslist):
                newline += "{:18s}".format(spec)
                if ii % 4 == 3 or ii == len(specieslist)-1:
                    seclines[insection].append(newline + endline)
                    newline = ''
                
        elif insection == "TRANS":
            # Get rid of unwanted species
            if startswithany(line,QSS) is None:
                seclines[insection].append(line)
            
        elif insection == "REACTIONS":
            # Do not save anything
            seclines[insection] = []
            
        elif insection == "ELEMENTS":
            # Just write the line
            seclines[insection].append(line)
            
        elif startswithany(line,sections) is not None:
            # Choose to be in a section and write
            insection = startswithany(line,sections)
            seclines[insection] = []
            print (' In ' + insection)
            outdat.write(line)

        else :
            # Not in any section, just print the lines
            outdat.write(line)

                
                
                
            
