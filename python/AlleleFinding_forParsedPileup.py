#Made by Sofia Randelin, 12/07/2017, randelin.sofia.s@student.uta.fi

#usage: python AlleleFinding_forParsedPileup.py
#example: python AlleleFinding_forParsedPileup.py

#input explanation:
# 1. Parsed pileup file: tab delimited BED file created with CountReads_fromPileup.py
#"CHR START END REF NUM A T C G N %-A %-T %-C %-G %-N"
# 2. SNP Alleles: tab delimited BED file (created from CEL files with R package Bioconductor Crlmm, for example)
#"CHR START END ALLELE"
#3. Select which alleles are processed (1, 2, 3 or all). 1 = AA. 2 = AB, 3 = BB
#4. Enter percentage treshold after allele selection. Only SNPs which have reads over treshold for a
# singular nuckleotide. Select zero (0) if you want to analyze all the SNPs.
# 5. Enter name for the new BED file

#output explanation:
#tab delimited file (BED)
#file created: entered file name + '_alleleCoverage.bed'
#CHR START END REF NUM %-A %-T %-C %-G %-N ALLELE


# Parsed pileup file is read and saved in a dict, counted_coverage
def read_pileupcoverage(counted_path, percent):
    try:
        with open(counted_path, 'r') as parsed_pileup:
            print("Opening file from: ", counted_path)

            counted_coverage = {}

            #HEAD
            parsed_pileup.readline()

            #data rows
            for row in parsed_pileup:
                row.rstrip()
                chromo, start, end, ref, num, A, T, C, G, N, pctA, pctT, pctC, pctG, pctN = row.split()
                info = [start, end, ref, num, pctA, pctT, pctC, pctG, pctN]

                percentages = [float(pctA), float(pctT), float(pctC), float(pctG), float(pctN)]

                # All the SNPs are added to dict
                if float(percent) == 0:
                    if chromo not in counted_coverage:
                        counted_coverage[chromo] = []

                    counted_coverage[chromo].append(info)
                else:
                    help = 0
                    # The data is saved, if there is a nuckleotide with percentage over the treshold.
                    for i in percentages:
                        if i <= (100-float(percent)):   # if under (100-treshold), the checking continues
                            help += 1
                        elif i >= float(percent):   # coverage over the treshold -> the data is saved in dict

                            if chromo not in counted_coverage:
                                counted_coverage[chromo] = []

                            counted_coverage[chromo].append(info)
                        else:   # if any percentage is between (100-treshold) and treshold, the row is skipped
                            continue

            parsed_pileup.close()

            print("Parsed pileup file has been read.")

            # Lists in dict values are arranged according to SNP end position
            for chromo in counted_coverage:
                counted_coverage[chromo] = sorted(counted_coverage[chromo], key=lambda x: int(x[1]))

            return counted_coverage

    except:
        print("Error in opening parsed pileup BED file.")
        return False


# Writing new file
def write_file(name, data):

    # File is created with given name and extension "_alleleCoverage.bed"
    new_file = name + "_alleleCoverage.bed"
    allele_percentages = open(new_file, "w")
    text = "CHR" + "\t" + "START" + "\t" + "END" + "\t" + "REF" + "\t" + "NUM" + "\t" + "%-A" + "\t" + \
           "\t" + "%-T" + "\t" + "%-C" + "\t" + "%-G" + "\t" + "%-N" + "\t" + "ALLELE" + "\n"
    allele_percentages.write(text)

    for x in data:
        text = ""
        for info in x:
            text = text + info + "\t"
        text += "\n"
        allele_percentages.write(text)

    allele_percentages.close()
    print(new_file, "created.")


# Binary search SNP position
def binary_search(order, position, allele, chrom):

    start = 0
    end = len(order)

    # Binary search
    while start != end:
        split = (start + end) // 2
        value = int(order[split][1])
        if value == position:
            # New list is created and data saved to it
            new_row = []
            new_row.append(chrom)
            for y in order[split]:
                new_row.append(y)
            # Allele information is added to new list
            new_row.append(allele)


            return new_row


        # Changing the search to first half of the list
        elif value > position:
            end = split

        # Searching the second half of the list
        elif value < position:
            start = split + 1


    return []



# Reading SNP allele file and searching same position from parsed pileup data
def read_and_find_alleles(path_SNP, parsed_reads, selected_allele):


    try:
        with open(path_SNP, 'r') as SNP_alleles:
            print("Opening file from:", path_SNP)

            information = []

            print("Processing...")
            for row in SNP_alleles:
                row.rstrip()
                if len(row.split()) == 4:
                    chrom, start, end, allele = row.split()

                    position = int(end)

                    if selected_allele != ("all" or "ALL" or "All"):
                        if selected_allele != allele:
                            continue

                    # Finding the same position from parese_reads (counted pileup)
                    if chrom in parsed_reads:
                        # Position is searched from parsed pileup data and saved in a SNPdata list
                        SNPdata = binary_search(parsed_reads[chrom], position, allele, chrom)

                        # If the SNP position wasn't found, the data won't be saved and we move on to the next SNP
                        if SNPdata == []:
                            continue

                        # SNPdata list is saved in information list
                        else:
                            information.append(SNPdata)


        SNP_alleles.close()


        return information


    except:
        print("Error in opening file.")
        return False




def main():

    # Path for parsed pileup file
    path = input("The program will combine info from parsed pileups and SNP alleles \n"
                 "(chr, start, end, allele). Both input files must be tab/space delimited files. \n"
                 "New BED files will be created with extension '_alleleCoverage.bed'. \n \n"
                 "Start by entering the path of your parsed pileup file (BED): ")

    # Path for SNP allele file
    path_SNP = input("Enter path of the SNP alleles (BED): ")

    allele_selection = input("Select the alleles which you want to save in the \n"
                             "processed file (1/2/3/all): ")


    treshold_percent = input("Enter SNP allele coverage treshold percentage \n"
                             "without %-sign (for example 80.5): ")



    # New file name
    file_name = input("Enter name for the new file: ")

    print()

    # Parsed pileup file is read in read_pileupcoverage.
    parsed_reads = read_pileupcoverage(path, treshold_percent)

    if parsed_reads != False:

        # Information of parsed pileup file is combined with SNP allele data and the combined data returned back to main.
        data = read_and_find_alleles(path_SNP, parsed_reads, allele_selection)

        if data != False:

            # New, combined data is written in a file.
            write_file(file_name, data)

main()