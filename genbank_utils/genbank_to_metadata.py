__author__ = 'jayna'

from Bio import SeqIO
from datetime import datetime, date
import os

month_to_num = {'Jan': 1, 'Feb': 2, 'Mar': 3, 'Apr': 4, 'May': 5, 'Jun': 6, 'Jul': 7, 'Aug': 8, 'Sep': 9, 'Oct': 10,
                'Nov': 11, 'Dec': 12}


def process_gb(gb_file, working_dir, virus_family):

    ### GENBANK FILE ###
    # create a variable to read your genbank file (i.e. gb_file). Note "r", which indicates the file is readable.
    input_handle = open(gb_file, "r")
    count = 0

    ### METADATA FILE ###
    # out_tablefile is a string variable that corresponds to the filename of the metadata file. Here I am using str_replace so that the filename is similar to the output fasta filena,e

    out_tablefile = working_dir + "/" + virus_family + "_metadata.csv"

    # Similar to output_handle, create a variable so that you can write to the metadata file.
    out_tablefile_handle = open(out_tablefile, "w")

    # Write out the column names to the metadata file.
    out_tablefile_handle.write(
        "Sequence Name,Accession_no,Strain,Organism,Species/Genus,Host,Country,Collection_date,Year,Sequence_length,Pubmed_id\n")


    ### PROCESS GENBANK FILE ###

    for index, record in enumerate(SeqIO.parse(input_handle, "genbank")):

        # extract genus of virus sequence
        taxonomy = record.annotations["taxonomy"]

        if virus_family == "Coronaviridae" and len(taxonomy) > 9:

            species = taxonomy[8] + ":" + taxonomy[9]

        elif virus_family == "Phenuiviridae" and len(taxonomy) > 8:

            species = taxonomy[7] + ":" + taxonomy[8]
            print(species)

        elif virus_family == "Paramyxoviridae" and len(taxonomy) > 9:
            species = taxonomy[8] + ":" + taxonomy[9]

        elif virus_family == "Paramyxoviridae" and len(taxonomy) > 8:
            species = (taxonomy[7] +
                       ":" + taxonomy[8])

        elif virus_family == "Paramyxoviridae" and len(taxonomy) < 8:
            species = taxonomy[6] + ":" + taxonomy[7]

        elif virus_family == "Rhabdoviridae" and len(taxonomy) > 9:

            species = taxonomy[8] + ":" + taxonomy[9]

        elif len(taxonomy) >= 9:

            species = taxonomy[7] + ":" + taxonomy[8]

        else:
            species = taxonomy[len(taxonomy) - 2] + ":" + taxonomy[len(taxonomy) - 1]
            print(species)


        # For each record in the genbank file, extracts the different subfields (or qualifiers) in the FEATURES field
        features_qualifiers = record.features[0].qualifiers

        # strain_name = record.annotations['accessions'][0]

        # Metadata / fields we are interested in extracting per sequence
        d = ""
        host = "NA"
        isolate = "NA"
        country = "NA"
        accession_no = record.annotations['accessions'][0]
        organism = features_qualifiers['organism'][0]

        # Get the Pubmed id of the publication associated with the record/sequence
        reference = record.annotations['references'][0]
        pubmed_id = reference.pubmed_id

        # Extracting which host the sequence was obtained from. This could be stored in /host= field or /isolation source
        if 'host' in features_qualifiers.keys():

            host = features_qualifiers['host'][0]

            host = str(host).replace(' ', ':')

            # [optional] print out host to screen to see how the above line changes the output

        elif 'isolation_source' in features_qualifiers.keys():
            host = features_qualifiers['isolation_source'][0]

            host = str(host).replace(' ', ':')

        if 'Homo' in host:
            continue

        if "strain" in host or "NA" in host:
            continue

        # Extract the country information of where the sequence was sampled from
        if 'geo_loc_name' in features_qualifiers.keys():
            country = features_qualifiers['geo_loc_name'][0]
            country_parts = country.split(':')
            country = country_parts[0]
            country = country.replace(" ", ":")

        # Extract the date of collection of the sequence
        if 'collection_date' in features_qualifiers.keys():
            d = features_qualifiers['collection_date'][0]

        # Extract isolate name
        if 'strain' in features_qualifiers.keys():
            isolate = features_qualifiers['strain'][0]

        if 'note' in features_qualifiers.keys():
            notes = features_qualifiers['note'][0]

            if 'experiment' in notes or 'vaccine' in notes:
                continue

        if 'lab_host' in features_qualifiers.keys():
            continue

        # The Date field may need further processing.
        t_date, year = process_date(d)

        name = accession_no + "_" + country + "_" + host + "|" + str(t_date) + "|" + str.replace(species, " ",
                                                                                                 "..") + "|" + str.replace(
            organism, " ", ":")

        count += 1

        # This line writes out all the metadata collected for each record to the metadata output file (out_tablefile)
        out_tablefile_handle.write(
            name + "," +
            accession_no + "," +
            isolate + "," +
            organism + "," +
            species + "," +
            host + "," +
            country.replace(":", " ") + "," +
            str(t_date) + "," +
            str(year) + "," +
            str(len(record.seq)) + "," +
            str(pubmed_id) + "\n")

        out_tablefile_handle.flush()


    # close all the files once finished processing the genbank file.

    input_handle.close()
    out_tablefile_handle.close()


def process_gb_seg(gb_file, working_dir, virus_family):

    ### GENBANK FILE ###
    # create a variable to read your genbank file (i.e. gb_file). Note "r", which indicates the file is readable.
    input_handle = open(gb_file, "r")
    count = 0

    ### METADATA FILE ###
    # out_tablefile is a string variable that corresponds to the filename of the metadata file. Here I am using str_replace so that the filename is similar to the output fasta filena,e
    out_tablefile = working_dir + "/" + virus_family + "_metadata.csv"

    try:
        os.makedirs(working_dir + "/" + virus_family)

    except FileExistsError:
        # directory already exists
        pass

    # Similar to output_handle, create a variable so that you can write to the metadata file.
    out_tablefile_handle = open(out_tablefile, "w")

    # Write out the column names to the metadata file.
    out_tablefile_handle.write(
        "Sequence Name,Accession_no,Strain,Organism,Species/Genus,Host,Country,Collection_date,Year,Sequence_length,Pubmed_id\n")


    ### PROCESS GENBANK FILE ###

    for index, record in enumerate(SeqIO.parse(input_handle, "genbank")):

        # extract genus of virus sequence
        taxonomy = record.annotations["taxonomy"]


        # print(len(taxonomy), taxonomy)

        species_or_genus = "NA"
        taxonomy_length = len(taxonomy)
        if taxonomy_length == 8:
            species_or_genus = "unclassified " + taxonomy[7]

        elif virus_family == "Phenuiviridae" and len(taxonomy) > 8:
            species_or_genus = taxonomy[7] + ":" + taxonomy[8]

        elif taxonomy_length == 9:

            if "Mammarenavirus" in record.features[0].qualifiers['organism'][0]:
                species_or_genus = record.features[0].qualifiers['organism'][0]
            else:
                species_or_genus = taxonomy[8]

        elif taxonomy_length == 10:

            if "Orthohantavirus" in record.features[0].qualifiers['organism'][0]:
                species_or_genus = record.features[0].qualifiers['organism'][0]
            else:
                species_or_genus = taxonomy[9]

        elif taxonomy_length == 11:
            species_or_genus = taxonomy[10]

        # For each record in the genbank file, extracts the different subfields (or qualifiers) in the FEATURES field
        features_qualifiers = record.features[0].qualifiers

        # strain_name = record.annotations['accessions'][0]

        # Metadata / fields we are interested in extracting per sequence
        d = ""
        host = "NA"
        isolate = "NA"
        country = "NA"
        accession_no = record.annotations['accessions'][0]
        organism = features_qualifiers['organism'][0].replace(" ", "..")

        # Get the Pubmed id of the publication associated with the record/sequence
        reference = record.annotations['references'][0]
        pubmed_id = reference.pubmed_id

        # Extracting which host the sequence was obtained from. This could be stored in /host= field or /isolation source
        if 'host' in features_qualifiers.keys():

            host = features_qualifiers['host'][0]

            host = str(host).replace(' ', ':')

            # [optional] print out host to screen to see how the above line changes the output

        elif 'isolation_source' in features_qualifiers.keys():
            host = features_qualifiers['isolation_source'][0]

            host = str(host).replace(' ', ':')


        # Extract the country information of where the sequence was sampled from
        if 'geo_loc_name' in features_qualifiers.keys():
            country = features_qualifiers['geo_loc_name'][0]
            country_parts = country.split(':')
            country = country_parts[0]
            country = country.replace(" ", ":")

        # Extract the date of collection of the sequence
        if 'collection_date' in features_qualifiers.keys():
            d = features_qualifiers['collection_date'][0]

        # Extract isolate name
        if 'strain' in features_qualifiers.keys():
            isolate = features_qualifiers['strain'][0]

        if 'note' in features_qualifiers.keys():
            notes = features_qualifiers['note'][0]

            if 'experiment' in notes or 'vaccine' in notes:
                continue

        if 'lab_host' in features_qualifiers.keys():
            continue

        if 'cell culture' in host:
            continue

        # The Date field may need further processing.
        t_date, year = process_date(d)

        # species_or_genus = get_virus_taxonomic_rank(record, taxonomic_level)

        name = accession_no + "_" + country + "_" + host + "|" + str(t_date) + "|" + str.replace(species_or_genus, " ",
                                                                                                 "..") + "|" + str.replace(
            organism, " ", ":")

        count += 1



        # This line writes out all the metadata collected for each record to the metadata output file (out_tablefile)
        out_tablefile_handle.write(
            name + "," +
            accession_no + "," +
            isolate + "," +
            organism + "," +
            species_or_genus + "," +
            host + "," +
            country.replace(":", " ") + "," +
            str(t_date) + "," +
            str(year) + "," +
            str(len(record.seq)) + "," +
            str(pubmed_id) + "\n")

        out_tablefile_handle.flush()


    # close all the files once finished processing the genbank file.

    input_handle.close()
    out_tablefile_handle.close()


def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


def process_date(d):
    t_date = "NA"
    year = "NA"

    # If the date string is longer than 4 characters than in addition to year, likely contains information about month and day of sampling
    if len(d) > 4:

        # splits date string by "-"
        parts = d.split('-')

        # I would suggest printing out the parts variable above to better understand what the above line does

        # This section undertakes further processing to identify which components correspond to month and year parts of the date,
        # and use this information to convert the string date e.g. 2002-01-03 or Dec-2001 into decimalised date e.g. 2002.013
        if len(parts) == 2:

            # print time.strptime(d, "%b %y")
            # print record.features[2].qualifiers['isolate'][0]

            if not parts[0].isdigit():

                # parts[0] must contain a string month.
                # parts[1] must contain a integer year.
                # no day information so assume date is in middle of the month.

                year = int(parts[1])
                d_date = datetime(int(parts[1]), month_to_num[parts[0]], 15)
                days_in_year = datetime.strftime(d_date, '%j')
                total_datetime = datetime(int(parts[1]), 12, 31)
                total_days_in_year = date.strftime(total_datetime, '%j')
                decimal_date = float(parts[1]) + float(days_in_year) / float(total_days_in_year)

            elif parts[0].isdigit():
                # parts[0] must contain a string year.
                year = int(parts[0])
                d_date = datetime(int(parts[0]), int(parts[1]), 15)
                days_in_year = datetime.strftime(d_date, '%j')
                total_datetime = datetime(int(parts[0]), 12, 31)
                total_days_in_year = date.strftime(total_datetime, '%j')
                decimal_date = float(parts[0]) + float(days_in_year) / float(total_days_in_year)

            else:
                # parts[1] must contain a string year.
                year = int(parts[1])
                d_date = datetime(int(parts[1]), month_to_num[parts[0]], 15)
                days_in_year = datetime.strftime(d_date, '%j')
                total_datetime = datetime(int(parts[1]), 12, 31)
                total_days_in_year = date.strftime(total_datetime, '%j')
                decimal_date = float(parts[1]) + float(days_in_year) / float(total_days_in_year)

            t_date = str('%.4f' % decimal_date)
        elif len(parts) == 1:
            year = parts[0]
        else:
            # check if month is in a number or string format
            if parts[1].isdigit():
                # parts[0] must contain a integer year.
                # parts[1] must contain a integer month.
                # parts[2] must contain a integer day.

                year = int(parts[0])

                d_date = datetime(int(parts[0]), int(parts[1]), int(parts[2]))
                days_in_year = datetime.strftime(d_date, '%j')
                total_datetime = datetime(int(parts[0]), 12, 31)
                total_days_in_year = date.strftime(total_datetime, '%j')
                decimal_date = float(parts[0]) + float(days_in_year) / float(total_days_in_year)

                t_date = str('%.4f' % decimal_date)

            else:
                # parts[0] must contain a integer day.
                # parts[1] must contain a string month.
                # parts[2] must contain a integer year.

                if parts[1] in month_to_num and is_number(parts[2]) == True:
                    year = int(parts[2])
                    d_date = datetime(int(parts[2]), month_to_num[parts[1]], int(parts[0]))
                    days_in_year = datetime.strftime(d_date, '%j')
                    total_datetime = datetime(int(parts[2]), 12, 31)
                    total_days_in_year = date.strftime(total_datetime, '%j')
                    decimal_date = float(parts[2]) + float(days_in_year) / float(total_days_in_year)

                    t_date = str('%.4f' % decimal_date)

                # print t_date
                # print ">"+name
                # print record.seq

    # if the date is already a number then it must correspond to year
    elif is_number(d):
        t_date = float(d)
        year = t_date

    return t_date, year
