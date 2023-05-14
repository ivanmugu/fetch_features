"""Helper functions to fetch, parse and save results.

This file is part of fetch_features
BSD 3-Clause License
Copyright (c) 2023, Ivan Munoz Gutierrez
"""
import os
from datetime import datetime
import sqlite3
import pandas as pd
from pathlib import Path
from typing import IO

from Bio import Entrez
from Bio import SeqIO


def make_uid_list_from_txt(infile: Path) -> list:
    """Make list of UIDs from a txt file.

    Read a txt file that contains a list of UIDs, as accession or BioSample
    numbers. Then, create a python list with these UIDs.

    Parameters
    ----------
    infile : Path
        Path to txt file that contains UIDs as accession or BioSample
        numbers. Every UID must be separated by a new line character
        (enter). The list can be created in excell by saving the file as
        txt or using a text editor.

    Returns
    -------
    list_accessions: list
        List of accession numbers. For example, a returned list from a txt
        file with nine accession number would look like the following list:
        ['CP049609.1', 'CP028704.1', 'CP043542.1', 'CP040107.1',
        'CP041747.1', 'CP042638.1', 'CP015023.1', 'CP049163.1',
        'CP051714.1']
    """
    # Open infile.txt.
    with open(infile, 'r') as reader:
        # Create a list of unique identifiers.
        uid_list = reader.readlines()

    # Remove the '\n' character.
    for i, uid in enumerate(uid_list):
        uid_list[i] = uid.replace('\n', '')

    return uid_list


def make_uid_list_from_xlsx(infile: Path) -> list:
    """Make list of UIDs from a xlsx file.

    Read an excel file that contains a list of UIDs, as accession or
    BioSample numbers. Then, create a python list with these UIDs.
    The list must be in the first column of the first spreadsheet tab.

    Parameters
    ----------
    infile : Path
        Path to xlsx file that contains UIDs as accession or BioSample
        numbers. The list must be in the first column of the first
        spreadsheet tab.

    Returns
    -------
    list_accessions: list
        List of accession numbers. For example, a returned list from a txt
        file with nine accession number would look like the following list:
        ['CP049609.1', 'CP028704.1', 'CP043542.1', 'CP040107.1',
        'CP041747.1', 'CP042638.1', 'CP015023.1', 'CP049163.1',
        'CP051714.1']
    """
    # Open excel file and read first sheet.
    df = pd.read_excel(infile, sheet_name=0)
    # Convert the first column into a list.
    uid_list = df.iloc[:, 0].tolist()
    return uid_list


def make_uid_list(infile: Path) -> list:
    """Make a list of unique identifiers from an input file."""
    infile_extention = infile.name.split(".")[-1]
    if infile_extention == 'txt':
        uid_list = make_uid_list_from_txt(infile)
    else:
        uid_list = make_uid_list_from_xlsx(infile)
    return uid_list


def make_uid_batch_list(uid_list: list, batch_size: int = 200) -> list:
    """Make batches of commma-delimited UIDs as accession or Biosample numbers.

    Parameters
    ----------
    uid_list : list
        List of UIDs to be proccessed
    batch_size : int
        Size of batches. NCBI suggests to request a maximum of 500
        accession numbers.

    Returns
    -------
    submission_list : list
        List of strings containg UIDs separated by commas. An example of a
        submission_list created with accession numbers and a batch size of
        three look like the following:
        ['CP049609.1,CP028704.1,CP043542.1',
        'CP040107.1,CP041747.1,CP042638.1',
        'CP015023.1,CP049163.1,CP051714.1']
    """
    # Number of UIDs to process.
    number_seq = len(uid_list)
    # Declare the list of batches of comma-delimited UIDs to return after
    # processing.
    submission_list = []
    # Counter to access the list_accessions.
    counter_accessions = 0

    # Loop to create the list of UIDs by batches
    for start in range(0, number_seq, batch_size):
        end = min(number_seq, start + batch_size)
        # This list is going to save temporarily the batch of UIDs that are
        # going to be converted into a string of comma-separated UIDs
        set_list = []
        # Make batches.
        for _ in range(start, end):
            set_list.append(uid_list[counter_accessions])
            counter_accessions += 1
        # Convert the list into string.
        set_list = ','.join(set_list)
        submission_list.append(set_list)

    return submission_list


def get_biosample_numbers(submission_list: list, email_address: str) -> list:
    """Retrieve BioSample numbers from a list of batches of accession numbers.

    This function uses Entrez.epost to upload batches of accession numbers to
    the Entrez History server. To avoid problems with large batches of
    accession numbers, limit the number of accession number per batch to 200.

    Parameters
    ----------
    submission_list : list
        List of strings containg batches of accession numbers separated by
        commas. The following is an example of a submission_list containing
        batches of three accession numbers per item in the list:
        ["CP049609.1,CP028704.1,CP043542.1",
        "CP040107.1,CP041747.1,CP042638.1",
        "CP015023.1,CP049163.1,CP051714.1"]
    email_address : string
        User's email address is requested by NCBI

    Returns
    -------
    biosample_numbers : list
        List of BioSample numbers retrieved from Entrez by using the provided
        list of accession numbers. The following is a return list containing
        the corresponding BioSample numbers of the example provided in the
        Parameters section:
        ["SAMN07169263", "SAMN08875353", "SAMEA104140560",
        "SAMN05360217", "SAMN12302771", "SAMN12500846",
        "SAMN04202539", "SAMN14133047", "SAMN14609782"]

    Notes
    -----
    For more information about Entrez.epost read chapter 9.4 of Biopython
    Tutorial and Cookbook (Biopython 1.76) and visit:
    https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.EPost
    """
    # Create list to save BioSample numbers.
    biosample_numbers = []

    # Provide email address to GeneBank.
    Entrez.email = email_address

    # Initialize last accession number.
    end = 0

    # Fetch BioSample numbers from GenBank by batches.
    for _, submission in enumerate(submission_list):
        start = end
        # submission_list is a list of accession numbers separated by
        # commas. Therefore, the number of commas indicate the number of
        # accession numbers.
        batch_size = submission.count(',') + 1
        end = end + batch_size

        # Print download batch record.
        print(f"Retrieving BioSample numbers from record {start + 1} to {end}")

        # Post the submission_list.
        # Because we are requesting information from a huge list of acc
        # numbers, we have to use the ".epost" function which uploads a
        # list of UIs (acc numbers) for use in subsequent searches.
        # From .epost we can get the QueryKey and the WebEnv which define
        # our history session and can be used to performe searches of data.
        posting = Entrez.epost('nuccore', id=submission)
        search_results = Entrez.read(posting)

        # Copy cookie "WebEnv" and query "QueryKey" from our history
        # session to keep track of our batch fetching.
        # WevEnv -> Web environment string returned from a previous
        # ESearch, EPost or ELink call; QueryKey -> Integer query key
        # returned by a previous ESearch, EPost or ELink call
        webenv = search_results["WebEnv"]
        query_key = search_results["QueryKey"]

        # Get the batch information
        # db -> database, nuccore -> nuleotide, rettype -> retrieval type
        # retmode -> determines the format of the return output
        # retstart -> sequential index of the first UID in the retrieved
        # set_number to be shown in the XML output
        # retmax -> total number of UIDs from the retrieved set_number to be
        # shown in the XML output
        # idtype-> specifies the type of identifier to return for sequence
        # databases, acc -> accesion number
        fetch_handle = Entrez.efetch(
            db="nuccore",
            rettype="gb",
            retmode="text",
            retstart=0,
            retmax=batch_size,
            webenv=webenv,
            query_key=query_key,
            idtype="acc"
        )

        # Pars throw the fetched information
        for seq_record in SeqIO.parse(fetch_handle, "gb"):
            # Loop over database cross-references (dbxrefs)
            for xref in enumerate(seq_record.dbxrefs):
                list_dbxrefs = xref[1].split(':')
                # Get the BioSample number
                if 'BioSample' in list_dbxrefs:
                    biosample_numbers.append(list_dbxrefs[1])
                    break

    return biosample_numbers


def parser(fetch_handle: IO, set_number: int) -> tuple:
    """Parse information fetched from nuccore NCBI.

    Parameters
    ----------
    fetch_handle : IO Entre.efetch object
        Pointer to a network connection that will allow to download and parse
        the requested sequences from the internet.
    set_number : integer
        Number of the batch of sequences that is being analyzed. This number
        will help to keep track of the downloaded sequences.

    Returns
    -------
    records : list
        List of dictionaries. Every dictionary corresponds to features of one
        accession number.
    ref_records: list
        List of dictionaries. Every dictionary carries information of the
        references of an accession number. One accession number can have one or
        more references. Therefore, this list may be longer than `records`.
    """
    # Create a list to save records in dictionaries.
    records = []
    # Create a list to save the references related to each accession number.
    ref_records = []

    # Parse the fetched information.
    for seq_record in SeqIO.parse(fetch_handle, "gb"):
        # Dictionary to store fetched information.
        info = {}
        # Keep track of the set_number (or batch) analyzed.
        info['set_batch'] = set_number
        # Extract `description`.
        info['description'] = seq_record.description
        # Extract sequence id, i.e. `accession` number.
        info['accession'] = seq_record.id
        # Extract `size` from `.seq`.
        # `.seq` is an object with the sequence itself.
        info['size'] = len(seq_record.seq)

        #######################################################################
        # '.annotations' is a dictionary of aditional information about the
        # sequence as last modification date, topology, sequence_version,
        # organims, references, etc.
        #######################################################################
        # Extract `date` and give format.
        if 'date' in seq_record.annotations:
            mod_date = seq_record.annotations['date']
            mod_date = datetime.strptime(mod_date, '%d-%b-%Y')
            mod_date = datetime.strftime(mod_date, "%Y-%m-%d")
            info['mod_date'] = mod_date
        else:
            info['mod_date'] = 'missing'
        # Extract `topology`.
        if "topology" in seq_record.annotations:
            info['topology'] = seq_record.annotations["topology"]
        else:
            info['topology'] = "missing"

        # Check whether it is `chromosome` or `plasmid`. The 3,000,000 length
        # used in `size` is an arbitrary number that considers chormosomes of
        # that length.
        # Check if label `chromosome` is in the `description`.
        if 'chromosome' in info['description']:
            info['molecule'] = 'chromosome'
        # If not, consider a chromosome the at least 3,000,000 in legth.
        elif info['size'] >= 3_000_000 and info['topology'] == 'circular':
            info['molecule'] = 'chromosome'
        # If it is not a chromosome, the check if it is a `plasmid`.
        elif 'plasmid' in info['description']:
            info['molecule'] = 'plasmid'
        else:
            info['molecule'] = 'missing'

        # '.features' is a list of SeqFeatures objects with more structured
        # information about the features on a sequence.
        feature = seq_record.features

        # Loop over list feature.
        for index in feature:
            # '.type' is only a description of the type of feature
            # that could be source, CDS, gene, etc.
            # In source we can find organism, strain, host, country, etc.
            if index.type == "source":
                # Create a dictionary of the qualifiers from source.
                dictionary = dict(index.qualifiers)
                # '.get' gives a list
                # Extract molecule type (`mol_type`).
                if "mol_type" in dictionary:
                    info['mol_type'] = dictionary.get("mol_type")[0]
                else:
                    info['mol_type'] = 'missing'
                # Extract `organism`.
                if 'organism' in dictionary:
                    info['organism'] = dictionary.get('organism')[0]
                else:
                    info['organism'] = 'missing'
                # Extract `strain`.
                if 'strain' in dictionary:
                    info['strain'] = dictionary.get('strain')[0]
                else:
                    info['strain'] = 'missing'
                # Extract `isolation_source`.
                if 'isolation_source' in dictionary:
                    info['isolation_source'] =\
                        dictionary.get('isolation_source')[0]
                else:
                    info['isolation_source'] = 'missing'
                # Extract `host`.
                if 'host' in dictionary:
                    info['host'] = dictionary.get('host')[0]
                else:
                    info['host'] = 'missing'
                # Extract `plasmid`.
                if 'plasmid' in dictionary:
                    info['plasmid'] = dictionary.get('plasmid')[0]
                else:
                    info['plasmid'] = 'missing'
                # Extract `country`.
                if 'country' in dictionary:
                    info['country'] = dictionary.get('country')[0]
                else:
                    info['country'] = 'missing'
                # Extract `coordinates`.
                if "lat_lon" in dictionary:
                    info['lat_lon'] = dictionary.get("lat_lon")[0]
                else:
                    info['lat_lon'] = 'missing'
                # Extract `collection_date`.
                if "collection_date" in dictionary:
                    info['collection_date'] =\
                        dictionary.get("collection_date")[0]
                else:
                    info['collection_date'] = 'missing'
                # Extract `note`.
                if "note" in dictionary:
                    info['note'] = dictionary.get("note")[0]
                else:
                    info['note'] = 'missing'
                # Extract `serovar`.
                if "serovar" in dictionary:
                    info['serovar'] = dictionary.get("serovar")[0]
                else:
                    info['serovar'] = 'missing'
                # Extract `collected_by`.
                if "collected_by" in dictionary:
                    info['collected_by'] = dictionary.get("collected_by")[0]
                else:
                    info['collected_by'] = 'missing'
                # Extract `genotype`.
                if "genotype" in dictionary:
                    info['genotype'] = dictionary.get("genotype")[0]
                else:
                    info['genotype'] = 'missing'

                break

        # '.dbxrefs' is a list populated from any PROJECT or DBLINK
        # Check if .dbxrefs has any information.
        if len(seq_record.dbxrefs) == 0:
            info['BioProject'] = 'missing'
            info['BioSample'] = 'missing'

        # Convert the list dbxrefs into dictionary.
        dictionary_dbxrefs = {}
        for i in range(len(seq_record.dbxrefs)):
            s = seq_record.dbxrefs[i].split(":")
            dictionary_dbxrefs[s[0]] = s[1]

        # Extract `BioProject` and `BioSample`.
        if "BioProject" in dictionary_dbxrefs:
            info['BioProject'] = dictionary_dbxrefs.get("BioProject")
        else:
            info['BioProject'] = 'missing'
        if "BioSample" in dictionary_dbxrefs:
            info['BioSample'] = dictionary_dbxrefs.get("BioSample")
        else:
            info['BioSample'] = 'missing'

        # #####################################################################
        # Mine `structured_comment`
        # #####################################################################
        # Sometimes the data needed is organized in `Genome-Assembly-Data` or
        # in `Assembly-Data`. Check if the sequence has `structured_comment`
        # and `Genome-Assembly-Data`.
        if "structured_comment" in seq_record.annotations and (
            "Genome-Assembly-Data" in (
                seq_record.annotations["structured_comment"])):
            # Extract `Assembly Method`.
            if ("Assembly Method" in (
                seq_record.annotations["structured_comment"]
                                      ["Genome-Assembly-Data"])):
                info['Assem_Method'] = (
                    seq_record.annotations["structured_comment"]
                                          ["Genome-Assembly-Data"]
                                          ["Assembly Method"])
            else:
                info['Assem_Method'] = "missing"
            # Extract `Genome Representation`.
            if "Genome Representation" in (
                seq_record.annotations["structured_comment"]
                                      ["Genome-Assembly-Data"]):
                info['Gen_Represent'] = (
                    seq_record.annotations["structured_comment"]
                                          ["Genome-Assembly-Data"]
                                          ["Genome Representation"])
            else:
                info['Gen_Represent'] = "missing"
            # Extract `Expected Final Version`.
            if "Expected Final Version" in (
                seq_record.annotations["structured_comment"]
                                      ["Genome-Assembly-Data"]
            ):
                info['Exp_Final_Ver'] = (
                    seq_record.annotations["structured_comment"]
                                          ["Genome-Assembly-Data"]
                                          ["Expected Final Version"]
                )
            else:
                info['Exp_Final_Ver'] = "missing"
            # Extract `Genome Coverage`.
            if "Genome Coverage" in (
                seq_record.annotations["structured_comment"]
                                      ["Genome-Assembly-Data"]):
                info['Gen_Coverage'] = (
                    seq_record.annotations["structured_comment"]
                                          ["Genome-Assembly-Data"]
                                          ["Genome Coverage"])
            else:
                info['Gen_Coverage'] = "missing"
            # Exteract `Sequencing Technology`.
            if "Sequencing Technology" in (
                seq_record.annotations["structured_comment"]
                                      ["Genome-Assembly-Data"]):
                info['Seq_Technol'] = (
                    seq_record.annotations["structured_comment"]
                                          ["Genome-Assembly-Data"]
                                          ["Sequencing Technology"])
            else:
                info['Seq_Technol'] = "missing"
        # Check if the sequence has `structured_comment` and `Assembly-Data`.
        elif "structured_comment" in seq_record.annotations and (
            "Assembly-Data" in (
                seq_record.annotations["structured_comment"])):
            # Extract `Assembly Method`.
            if "Assembly Method" in (
                seq_record.annotations["structured_comment"]
                                      ["Assembly-Data"]):
                info['Assem_Method'] = (
                    seq_record.annotations["structured_comment"]
                                          ["Assembly-Data"]
                                          ["Assembly Method"])
            else:
                info['Assem_Method'] = "missing"
            # Extract `Genome Representation`.
            if "Genome Representation" in (
                seq_record.annotations["structured_comment"]
                                      ["Genome-Assembly-Data"]):
                info['Gen_Represent'] = (
                    seq_record.annotations["structured_comment"]
                                          ["Genome-Assembly-Data"]
                                          ["Genome Representation"])
            else:
                info['Gen_Represent'] = "missing"
            # Extract `Expected Final Version`.
            if "Expected Final Version" in (
                seq_record.annotations["structured_comment"]
                                      ["Genome-Assembly-Data"]
            ):
                info['Exp_Final_Ver'] = (
                    seq_record.annotations["structured_comment"]
                                          ["Genome-Assembly-Data"]
                                          ["Expected Final Version"]
                )
            else:
                info['Exp_Final_Ver'] = "missing"
            # Extract `Genome Coverage`.
            if "Genome Coverage" in (
                seq_record.annotations["structured_comment"]
                                      ["Assembly-Data"]):
                info['Gen_Coverage'] = (
                    seq_record.annotations["structured_comment"]
                                          ["Assembly-Data"]
                                          ["Genome Coverage"])
            else:
                info['Gen_Coverage'] = "missing"
            # Extract `Sequencing Technology`.
            if "Sequencing Technology" in (
                seq_record.annotations["structured_comment"]
                                      ["Assembly-Data"]):
                info['Seq_Technol'] = (
                    seq_record.annotations["structured_comment"]
                                          ["Assembly-Data"]
                                          ["Sequencing Technology"])
            else:
                info['Seq_Technol'] = "missing"
        # Otherwise, the info is missing.
        else:
            info['Assem_Method'] = "missing"
            info['Gen_Coverage'] = "missing"
            info['Seq_Technol'] = "missing"
            info['Gen_Represent'] = "missing"
            info['Exp_Final_Ver'] = "missing"

        # Append info dictionary to records list.
        records.append(info)

        # ################################################
        # Get references and make an independent database.
        # ################################################
        if 'references' in seq_record.annotations:
            ref_records.extend(get_references(seq_record))

    return (records, ref_records)


def get_references(record: IO) -> list:
    """Get all references of accession number.

    Parameters
    ----------
    record : Biopython SeqIO.parser object
        TextIO object to parse a GenBank record.

    Returns
    -------
    references = list
        List of dictionaries with all the references of the GenBank record.
    """
    references = []
    for i, reference in enumerate(record.annotations['references']):
        info = {}
        info['accession'] = record.id
        info['reference_num'] = str(i)
        if len(reference.location) == 0:
            info['location'] = "N/A"
        else:
            start = reference.location[0].start
            end = reference.location[0].end
            info['location'] = f"bases {start} to {end}"
        info['authors'] = reference.authors
        info['title'] = reference.title
        info['journal'] = reference.journal
        info['medline_id'] = reference.medline_id
        info['pubmed_id'] = reference.pubmed_id
        info['comment'] = reference.comment
        references.append(info)
    return references


def clean_references(updated_features: list, reference_records: list) -> list:
    """Remove references from outdated accession numbers.

    Parameters
    ----------
    updated_features : list
        List with features from updated accession numbers.
    reference_records : list
        List with references from all accession numbers.

    Returns
    -------
    updated_references : list
        List with references from the updated accession numbers.
    """
    # Initialized list to store updated references.
    updated_references = []
    # Make a list of updated accession numbers.
    updated_accession = [feature['accession'] for feature in updated_features]
    # Iterate over `reference_recods` and extract the references of the updated
    # accession numbers.
    for reference in reference_records:
        if reference['accession'] in updated_accession:
            updated_references.append(reference)

    return updated_references


def clean_features(raw_data: list) -> list:
    """Clean up the features obtained from a BioSample number to get the most
    updated information

    Some BioSample numbers are associated to accession numbers that does not
    have any relevant information as contigs of few hundreds of nucleotides.
    Also, in addition to updated accession numbers, some BioSample numbers have
    outdated accession numbers. This function cleans the fetched information
    from GenBank using SQL.

    Parameters
    ----------
    raw_data : list
        List of dictionaries with fetched data from BioSample numbers.

    Returns
    -------
    recods : list
        List of dictionaries with the most updated information from BioSample
        numbers.
    """
    # Create empty features.db.
    open('features.db', 'w').close()

    # Open features.db for sqlite3.
    conn = sqlite3.connect('features.db')

    # Create a cursor.
    c = conn.cursor()

    # Create SQL table for the features of the raw data fetched from GenBank.
    c.execute("""CREATE TABLE features_raw (
                set_batch int,
                description text,
                accession text,
                size integer,
                molecule text,
                mod_date text,
                topology text,
                mol_type text,
                organism text,
                strain text,
                isolation_source text,
                host text,
                plasmid text,
                country text,
                lat_lon text,
                collection_date text,
                note text,
                serovar text,
                collected_by text,
                genotype text,
                BioProject text,
                BioSample text,
                Assem_Method text,
                Gen_Coverage text,
                Seq_Technol text,
                Gen_Represent text,
                Exp_Final_Ver
                )""")
    conn.commit()

    # Transfer results obtained from GenBank into the table features_raw
    # created in SQL.
    for _, raw_result in enumerate(raw_data):
        c.execute(
            """INSERT INTO features_raw VALUES(
                :set_batch, :description, :accession, :size, :molecule,
                :mod_date, :topology, :mol_type, :organism, :strain,
                :isolation_source, :host, :plasmid, :country, :lat_lon,
                :collection_date, :note, :serovar, :collected_by, :genotype,
                :BioProject, :BioSample, :Assem_Method, :Gen_Coverage,
                :Seq_Technol, :Gen_Represent, :Exp_Final_Ver
            )""", {
                'set_batch': int(raw_result['set_batch']),
                'description': raw_result['description'],
                'accession': raw_result['accession'],
                'size': int(raw_result['size']),
                'molecule': raw_result['molecule'],
                'mod_date': raw_result['mod_date'],
                'topology': raw_result['topology'],
                'mol_type': raw_result['mol_type'],
                'organism': raw_result['organism'],
                'strain': raw_result['strain'],
                'isolation_source': raw_result['isolation_source'],
                'host': raw_result['host'],
                'plasmid': raw_result['plasmid'],
                'country': raw_result['country'],
                'lat_lon': raw_result['lat_lon'],
                'collection_date': raw_result['collection_date'],
                'note': raw_result['note'],
                'serovar': raw_result['serovar'],
                'collected_by': raw_result['collected_by'],
                'genotype': raw_result['genotype'],
                'BioProject': raw_result['BioProject'],
                'BioSample': raw_result['BioSample'],
                'Assem_Method': raw_result['Assem_Method'],
                'Gen_Coverage': raw_result['Gen_Coverage'],
                'Seq_Technol': raw_result['Seq_Technol'],
                'Gen_Represent': raw_result['Gen_Represent'],
                'Exp_Final_Ver': raw_result['Exp_Final_Ver']
            }
        )
    conn.commit()

    # Get the most updated files and order by size.
    c.execute(
        """
        SELECT *
        FROM (
            SELECT table_one.*
                FROM features_raw table_one
                WHERE table_one.mod_date = (
                    SELECT MAX(table_two.mod_date)
                        FROM features_raw table_two
                    WHERE table_two.size = table_one.size
                )
        )
        WHERE molecule != "missing"
        ORDER BY size DESC
        """)

    # Extract the updated results from SQL. The fetchall will create a list
    # of tuples.
    results = c.fetchall()

    # Convert results (list of tuples) into a list of dictionaries.
    records = []
    info = {}
    for _, updated_result in enumerate(results):
        info['set_batch'] = int(updated_result[0])
        info['description'] = updated_result[1]
        info['accession'] = updated_result[2]
        info['size'] = int(updated_result[3])
        info['molecule'] = updated_result[4]
        info['mod_date'] = updated_result[5]
        info['topology'] = updated_result[6]
        info['mol_type'] = updated_result[7]
        info['organism'] = updated_result[8]
        info['strain'] = updated_result[9]
        info['isolation_source'] = updated_result[10]
        info['host'] = updated_result[11]
        info['plasmid'] = updated_result[12]
        info['country'] = updated_result[13]
        info['lat_lon'] = updated_result[14]
        info['collection_date'] = updated_result[15]
        info['note'] = updated_result[16]
        info['serovar'] = updated_result[17]
        info['collected_by'] = updated_result[18]
        info['genotype'] = updated_result[19]
        info['BioProject'] = updated_result[20]
        info['BioSample'] = updated_result[21]
        info['Assem_Method'] = updated_result[22]
        info['Gen_Coverage'] = updated_result[23]
        info['Seq_Technol'] = updated_result[24]
        info['Gen_Represent'] = updated_result[25]
        info['Exp_Final_Ver'] = updated_result[26]
        records.append(info.copy())
    conn.close()

    return records


def save_results_as(
        results: Path, ref_results: Path, save_as: str = 'csv'
        ) -> None:
    """Save results in the requested format.

    The files are saved by default in csv before this function is called. This
    function helps to make excel files if user requested. Also, by defaults,
    all files are named `results.csv`, `ref_results.csv`, or the equivalent
    names with the .xlsx extention.

    Parameters
    ----------
    results : Path
        Path to the saved results.
    ref_results : Path
        Path to saved the references results.
    save_as : str
        Format to save results. Option: `csv`, `excel`, or `csv-excel`.
    """
    # If requested `csv`, don't do anything, the files are ready.
    if save_as == 'csv':
        return
    # If requested `excel` or `csv-excel`, make excel files.
    if save_as == 'excel' or save_as == 'csv-excel':
        # Open csv files.
        df_results = pd.read_csv(results)
        df_ref_results = pd.read_csv(ref_results)
        # Export csv files as excel.
        new_results_path = results.parent / "results.xlsx"
        new_ref_results_path = ref_results.parent / "ref_results.xlsx"
        df_results.to_excel(new_results_path, index=False)
        df_ref_results.to_excel(new_ref_results_path, index=False)
    # If requested only `excel` delete the csv files.
    if save_as == 'excel':
        # Delete csv files
        os.remove(results)
        os.remove(ref_results)

    return
