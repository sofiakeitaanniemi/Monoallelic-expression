#Made by Sofia Randelin, 11/07/2017, randelin.sofia.s@student.uta.fi

#usage: python CountReads_fromPileup.py
#example: python CountReads_fromPileup.py
#input: Samtools pileup file
#output explanation:
#tab delimited file (BED)
#CHR START END REF NUM A T C G N %-A %-T %-C %-G %-N


import re

from collections import Counter

def remove_parts(reads):

    remove = re.findall('\^\w', reads)
    for j in remove:
        reads = reads.replace(j, "")

    # Remove insertions and deletions from string
    remove2 = []
    numbers = re.findall('[0-9]', reads)
    k = 0
    while k != len(reads):
        if reads[k] in numbers:
            delete = reads[k]
            l = int(reads[k])
            apu = 1
            while l != 0:
                delete += reads[k + apu]
                apu += 1
                l -= 1
            remove2.append(delete)
        k += 1

    # insertions and deletions are replaced by an empty string.
    for i in remove2:
        reads = reads.replace(i, "")

    # Remove start of segment marks and mapping quality marks which can be mixed with matches
    reads = reads.replace("^,", "")
    reads = reads.replace("^.", "")
    reads = reads.replace("^*", "")

    return reads



def count_match(counter, match):

    # Count matches
    match += counter['.']
    match += counter[',']
    match += counter['*']

    return match



def count_other(counter, nucleotide_counts):

    # All the nucleotides are counted
    for nucl in nucleotide_counts:
        nucleotide_counts[nucl] += counter[nucl]
        nucleotide_counts[nucl] += counter[nucl.lower()]

    # Return dict
    return nucleotide_counts



def count_percentages(nucleotides, num):

    pct_A = "{:.4f}".format((nucleotides["A"] / int(num)) * 100)
    pct_T = "{:.4f}".format((nucleotides["T"] / int(num)) * 100)
    pct_C = "{:.4f}".format((nucleotides["C"] / int(num)) * 100)
    pct_G = "{:.4f}".format((nucleotides["G"] / int(num)) * 100)
    pct_N = "{:.4f}".format((nucleotides["N"] / int(num)) * 100)

    return [pct_A, pct_T, pct_C, pct_G, pct_N]




def read_data(row, treshold):

    if len(row.split()) == 6:
        chr, pos, ref, num, reads, qualities = row.split()

        # If the coverage is smaller than given treshold,
        # the row will be skipped
        if int(num) <= int(treshold):
            return []

        # Useless parts for counting, like insertions and deletions, are removed
        reads = remove_parts(reads)

        # counter is created
        counter = Counter(reads)

        # Skipping the N CIGAR operator rows
        if counter['>'] > 0 or counter['<'] > 0:
            return []

        match = count_match(counter, 0)

        # Other nucleotides are counted in count_other
        nucleotides = count_other(counter, {"A": 0, "T": 0, "C": 0, "G": 0, "N": 0})

        # The reference nucleotide count is updated to match count
        for x in nucleotides:
            if x == ref:
                nucleotides[x] = match

        # Percentages for nucleotide reads are counted in count_percentages
        percentages = count_percentages(nucleotides, num)

        if sum(nucleotides.values()) != int(num):
            print("Error in counting, position:", chr, "   ", pos)
            print("Number of reads:", num, "    A:", nucleotides["A"], "    T:",
                  nucleotides["T"], "   C:", nucleotides["C"], "    G:", nucleotides["G"],
                  nucleotides["N"])
            quit()

        start = str(int(pos) - 1)

        # All the data is added to list called information which is added to list called counted
        information = [chr, start, pos, ref, num, nucleotides["A"], nucleotides["T"], nucleotides["C"],
                       nucleotides["G"], nucleotides["N"]]

        for percent in percentages:
            information.append(str(percent))

        return information

    else:
        return []

def main():

    path = input("The program will count the number and persentage of reference reads from pileup file. \n"
                 "BED files will be created with extension '_countedPileup.bed'. \n \n"
                 "Start by entering the path of your Samtools pileup file: ")

    file_name = input("Enter name for the new parsed pileup file: ")

    treshold = 0

    #Setting the debth treshold:
    treshold_true = False

    while treshold_true != True:
        treshold = input("Enter depth threshold: ")
        try:
            if int(treshold) > 0:
                treshold_true = True
        except ValueError:
            print("Input must be integer.")
            treshold_true = False


    counted = []


    # Opening file
    try:
        with open(path, 'r') as pileup:

            print("File is opened for reading.")
            print("Counting...")

            for row in pileup:
                row.rstrip()
                if len(row.split()) == 6:

                    information = read_data(row, treshold)

                    if information != []:
                        counted.append(information)

                else:
                    continue

            # Counted pileup data is arraged according to chromosomes
            counted_order = sorted(counted, key=lambda x: x[0])

            new_name = file_name + "_countedPileup.bed"
            file = open(new_name, "w")
            head = "#CHR" + "\t" + "START" + "\t" + "END" + "\t" + "REF" + "\t" + "NUM" + "\t" +\
                     "A"+ "\t" + "T" + "\t" + "C" + "\t" + "G" + "\t" + "N" + "\t" + "%-A" +"\t" + \
                     "\t" + "%-T" + "\t" + "%-C" + "\t" + "%-G" + "\t" + "%-N" + "\n"
            file.write(head)
            for info in counted_order:
                text = ""
                for value in info:
                    text = text + str(value) + "\t"
                text += "\n"
                file.write(text)

            file.close()
            print("File created:", new_name)

    except FileNotFoundError:
        print("File not found.")

    except PermissionError:
        print("Error! Enter the path with file name.")

main()