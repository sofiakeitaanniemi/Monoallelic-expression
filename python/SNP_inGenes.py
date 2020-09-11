#Made by Sofia Randelin, 12/07/2017, randelin.sofia.s@student.uta.fi

#usage: python SNP_inGenes.py
#example: python AlleleFinding_forParsedPileup.py

#input explanation:
# 1. Allele coverage file: tab delimited BED file created with AlleleFinding_forParsedPileup.py
#"CHR START END REF NUM %-A %-T %-C %-G %-N ALLELE"
# 2. Gene reference file: tab delimited BED file (USCS Gene Table, RefSeq genes, for example)
#"REFSEQ_ID CHR START END GENE_NAME"
# 3. Enter name for the new CSV file. CSV file consists of gene names and the number of SNPs in each gene.
#"GENE_NAME,NUMBER_OF_SNPs"
# 4. Enter name for the new BED position file. BED file consist of Gene names and positions for each SNP
# "CHR START END Gene_name"

#output explanation:
# 1. comma delimited file (CSV) with extension '_SNPinGenes.csv'
#"GENE_NAME,NUMBER_OF_SNPs_IN_GENE"
# 2. tab delimited file (BED) with extension '_SNPpositionsInGenes.bed'
#"CHR START END GENE_NAME"




def read_genes(path):
    ''' read_genes opens and reads a tab delimited file. The file must consist of 5 columns:
        Gene ID, chromosome, start position, end position, Gene name/symbol.
        The rows are added in a dict with Gene name as key and list of gene position data as value.
        The dict is returned to main.
        :param: path (str) 
        :return: gene_coordinates (dict)
    '''

    try:
        with open(path, 'r') as gene_file:

            print("File opened from:", path)

            gene_coordinates = {}


            #HEAD
            gene_file.readline()

            #data rows
            for row in gene_file:
                row.rstrip()
                info = row.split()
                # info = [gene ID,	chrom,	txStart,	txEnd,	gene name]

                chromosome = info[1]
                start = info[2]
                end = info[3]
                gene = info[4]


                #data is added to dict gene_coordinates
                if gene not in gene_coordinates:
                    gene_coordinates[gene] = []
                    gene_coordinates[gene].append([chromosome, start, end])

                else:

                    for i in gene_coordinates[gene]:
                        chrom = i[0]

                        if chrom == chromosome:
                            #The beginning and end positions which are saved in the dict
                            BEGIN = int(i[1])
                            END = int(i[2])
                            s = int(start)
                            e = int(end)

                            #The longest window of gene is saved

                            # Start position changes (smaller than before) while end position remains
                            if chrom == chromosome and s < BEGIN and BEGIN <= e and e <= END:
                                i[1] = start
                                break
                            # End position changes (larger than before) while start position remains
                            elif chrom == chromosome and e > END and BEGIN <= s and s <= END:
                                i[2] = end
                                break

                            # Both, start and end positions change because the window between them is bigger than
                            # originally in dict
                            elif chrom == chromosome and s < BEGIN and e > END:
                                i[1] = start
                                i[2] = end
                                break

                            #If the window is smaller than the one already existing in dict, the row data is not saved
                            elif e <= END and e >= BEGIN and s <= END and s >= BEGIN and chrom == chromosome:
                                break

                            # Gene location is in another chromosome or both the start and end coordinates
                            # were not between the previous coordinates
                            else:
                                gene_coordinates[gene].append([chromosome, start, end])


            gene_file.close()

            return gene_coordinates

    except:
        print("Error in reading gene BED file.")
        return False


def write_new_file(name, data, delimiter, extension, head):
    ''' write_new_file creates new CSV or BED file depending on the parameters. 
        The data from lists in a list (number_of_SNP or SNP_positions) is written row
        by row, comma as the delimiter (gene name,number of SNPs).
        :param: name (str), number_of_SNP (list) 
        :return: None
    '''
    # New file is written (gene, number of SNPs)
    new_name = name + extension
    new_file = open(new_name, "w")
    text = head + "\n"
    new_file.write(text)

    for info in data:
        text = ""
        for x in info:
            text = text + x + delimiter
        if delimiter == ",":
            text = text[:-1]
        text += "\n"
        new_file.write(text)

    new_file.close()
    print(new_name, "created.")


def countSNPs(genes, SNPlist, minimum):
    ''' SNPs are counted for each gene in countSNPs function. If there are 2 or more SNPs in a gene, the data is added
        as a list to a list called number_of_SNP which is returned to main.
    :param: genes (dict), SNPlist (list) 
    :return: number_of_SNP (list)
    '''
    number_of_SNP = []

    print("Counting SNPs in genes...")

    for gene in genes:
        num = 0  # Number of SNPs in a singular gene
        for x in genes[gene]:
            for SNP in SNPlist:
                if SNP[0] == x[0]:  # chromosome matches
                    if int(SNP[2]) > int(x[1]) and int(SNP[2]) < int(x[2]):  # SNP is between start and end of gene
                        num += 1  # +1 SNP inside gene

        # If there are more than minimum number of SNPs in a gene, the data is saved
        if num >= int(minimum):
            lista = [gene, str(num)]
            number_of_SNP.append(lista)

    return number_of_SNP



def find_SNP_positions(genes, SNPlist, minimum):
    ''' SNPs are counted for each gene in countSNPs function. If there are 2 or more SNPs in a gene, the data is added
        as a list to a list called number_of_SNP which is returned to main.
    :param: genes (dict), SNPlist (list) 
    :return: number_of_SNP (list)
    '''

    SNP_positions = []

    print("Finding SNPs positions in genes...")

    for gene in genes:
        help = []
        num = 0  # Number of SNPs in a singular gene
        for x in genes[gene]:
            for SNP in SNPlist:
                if SNP[0] == x[0]:  # chromosome matches
                    if int(SNP[2]) > int(x[1]) and int(SNP[2]) < int(x[2]):  # SNP is between start and end of gene
                        list = [SNP[0], SNP[1], SNP[2], SNP[3], gene]
                        help.append(list)

        # If there are more than minimum number of SNPs in a gene, the data is saved
        if len(help) >= int(minimum):
            for SNP in help:
                SNP_positions.append(SNP)

    return SNP_positions



def read_SNP_allele_coverage(path_bed):
    ''' read_SNP_allele_coverage opens and reads the BED file created in AlleleFinding_forParsedPileup.py.
        The file must consist of 10 columns:
        "chr, alku, loppu, ref, num, pctA, pctT, pctC, pctG, pctN"
        The data is added to SNPs list as a list and SNPs is returned to main.
        :param: path_bed (str) 
        :return: SNPs (list)
    '''
    try:
        with open(path_bed, 'r') as coverage_file:
            print("File opened from:", path_bed)

            SNPs = []

            #HEAD
            coverage_file.readline()

            #Data
            for row in coverage_file:
                row.rstrip()
                info = row.split()
                # info = [chr, alku, loppu, ref, num, pctA, pctT, pctC, pctG, pctN]
                if len(info) == 10:
                    SNPs.append(info)

            coverage_file.close()

            return SNPs

    except:
        print("Error in reading SNP allele coverage BED file.")
        return False


def main():


    path_bed = input("The program will search for SNPs in genes.  \n"
                     "CSV file will be created with extension '_SNPlist.csv'. \n"
                     "BED file will be created with extension '_SNPpositionsInGenes.bed' \n"
                     "CSV file consists of gene names and number of SNPs in genes. BED file has the positions"
                     "for these SNPs. \n \n"
                       "Start by entering the path of your allele coverage BED file: ")
    
    path_genes = input("Enter path of the tab delimited gene reference file (RefSeq ID, chr, start, end, gene name): ")

    file_name_CSV = input("Enter name for the new CSV file: ")

    file_name_positions = input("Enter name for the new SNP position file (BED): ")
    minimum_SNP_in_gene = input("Enter minimum number of SNPs in gene: ")

    # Tab delimited gene file (BED) is read and saved in a dict, genes
    genes = read_genes(path_genes)

    if genes != False:

        # SNP allele coverage file (BED) is read and saved in a list, SNPlist
        SNPlist = read_SNP_allele_coverage(path_bed)

        if SNPlist != False:

            # Number of SNPs are counted in each gene and saved in a list, number_of_SNP
            number_of_SNP = countSNPs(genes, SNPlist, minimum_SNP_in_gene)

            # SNP positions are saved in a list
            SNP_positions = find_SNP_positions(genes, SNPlist, minimum_SNP_in_gene)

            # The data from number_of_SNP is written in a CSV file (gene name,number of SNPs in gene)
            write_new_file(file_name_CSV, number_of_SNP, ",", "_SNPinGenes.csv", "gene,number_of_SNPs")

            # The data from SNP_positions is written in a BED file (chr start end gene_name)
            header = "#CHR" + "\t" + "START" + "\t" + "END" + "\t" + "GENE_NAME"
            write_new_file(file_name_positions, SNP_positions, "\t", "_SNPpositionsInGenes.bed", header)

main()


