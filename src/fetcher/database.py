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

    Parameters
    ----------
    infile : Path
        Path to txt file that contains UIDs as accession or BioSample
        numbers. Every UID must be separated by a new line character
        (enter). The list can be created in Excel by saving the file as
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
    # Open infile.txt and make a list of UIDs.
    with open(infile, 'r') as reader:
        uid_list = reader.readlines()
    # Remove the `\n` character of every uid.
    uid_list = [uid.rstrip() for uid in uid_list]

    return uid_list


def make_uid_list_from_xlsx(infile: Path) -> list:
    """Make list of UIDs from a xlsx file.

    Read an Excel file that contains a list of UIDs, as accession or
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


def make_uid_batches(uid_list: list, batch_size: int = 100) -> list:
    """Make list of batches of UIDs."""
    batches = []
    for i, uid in enumerate(uid_list):
        if (i % batch_size == 0) and (i == 0):
            batch = [uid]
        elif (i % batch_size == 0):
            batches.append(batch)
            batch = [uid]
        else:
            batch.append(uid)
    batches.append(batch)
    return batches


class Info:
    """"Class to store data retrieved from the nuccore database.

    All the attributes are initiated to `missing`. Hence, if a fetched
    attribute from the GenBank file is missing we reduce the if-else statements
    in the code.
    """

    def __init__(
        self, set_batch='missing', description='missing', accession='missing',
        size='missing', mod_date='missing', molecule='missing',
        topology='missing', mol_type='missing', organism='missing',
        strain='missing', isolation_source='missing', host='missing',
        plasmid='missing', country='missing', lat_lon='missing',
        collection_date='missing', note='missing', serovar='missing',
        collected_by='missing', genotype='missing', bioproject='missing',
        biosample='missing', assem_method='missing', gen_coverage='missing',
        seq_technol='missing', gen_represent='missing', exp_final_ver='missing'
    ):
        self.set_batch = set_batch
        self.description = description
        self.accession = accession
        self.size = size
        self.molecule = molecule
        self.mod_date = mod_date
        self.topology = topology
        self.mol_type = mol_type
        self.organism = organism
        self.strain = strain
        self.isolation_source = isolation_source
        self.host = host
        self.plasmid = plasmid
        self.country = country
        self.lat_lon = lat_lon
        self.collection_date = collection_date
        self.note = note
        self.serovar = serovar
        self.collected_by = collected_by
        self.genotype = genotype
        self.bioproject = bioproject
        self.biosample = biosample
        self.assem_method = assem_method
        self.gen_coverage = gen_coverage
        self.seq_technol = seq_technol
        self.gen_represent = gen_represent
        self.exp_final_ver = exp_final_ver


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
    # Create a list to save Info classes with records.
    records = []
    # Create a list to save the references related to each accession number.
    ref_records = []

    # Parse the fetched information.
    for seq_record in SeqIO.parse(fetch_handle, "gb"):
        # Initiate an Info class to store collected data.
        info = Info()
        # Keep track of the set_number (or batch) analyzed.
        info.set_batch = set_number
        # Extract `description`.
        info.description = seq_record.description
        # Extract sequence id, i.e. `accession` number.
        info.accession = seq_record.id
        # Extract `size` from `.seq`.
        # `.seq` is an object with the sequence itself.
        info.size = len(seq_record.seq)

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
            info.mod_date = mod_date
        # Extract `topology`.
        if "topology" in seq_record.annotations:
            info.topology = seq_record.annotations["topology"]

        # Check whether it is `chromosome` or `plasmid`. The 3,000,000 length
        # used in `size` is an arbitrary number that considers chromosomes of
        # that length.
        # Check if label `chromosome` is in the `description`.
        if 'chromosome' in info.description:
            info.molecule = 'chromosome'
        # If not, consider a chromosome the at least 3,000,000 in legth.
        elif info.size >= 3_000_000 and info.topology == 'circular':
            info.molecule = 'chromosome'
        # If it is not a chromosome, then check if it is a `plasmid`.
        elif 'plasmid' in info.description:
            info.molecule = 'plasmid'

        # `.features` is a list of SeqFeatures objects with more structured
        # information about the features on a sequence.
        # Loop over list of features to find the type source.
        for feature in seq_record.features:
            # `.type` is only a description of the type of feature that could
            # be source, CDS, gene, etc. The feature that we need is `source`.
            # In type `source` we can find organism, strain, host, country, etc.
            if feature.type == "source":
                # `qualifiers` is a Python dictionary of additional information
                # about the feature.
                source = feature.qualifiers
                break
        # Extract all the info from the type `source`.
        # `source` is a dictionary where the key-values are lists. Therefore,
        # the `.get` method retrieves a list.
        # Extract molecule type (`mol_type`).
        if "mol_type" in source:
            info.mol_type = source.get("mol_type")[0]
        # Extract `organism`.
        if 'organism' in source:
            info.organism = source.get('organism')[0]
        # Extract `strain`.
        if 'strain' in source:
            info.strain = source.get('strain')[0]
        # Extract `isolation_source`.
        if 'isolation_source' in source:
            info.isolation_source = source.get('isolation_source')[0]
        # Extract `host`.
        if 'host' in source:
            info.host = source.get('host')[0]
        # Extract `plasmid`.
        if 'plasmid' in source:
            info.plasmid = source.get('plasmid')[0]
        # Extract `country`.
        if 'country' in source:
            info.country = source.get('country')[0]
        # Extract `coordinates`.
        if "lat_lon" in source:
            info.lat_lon = source.get("lat_lon")[0]
        # Extract `collection_date`.
        if "collection_date" in source:
            info.collection_date = source.get("collection_date")[0]
        # Extract `note`.
        if "note" in source:
            info.note = source.get("note")[0]
        # Extract `serovar`.
        if "serovar" in source:
            info.serovar = source.get("serovar")[0]
        # Extract `collected_by`.
        if "collected_by" in source:
            info.collected_by = source.get("collected_by")[0]
        # Extract `genotype`.
        if "genotype" in source:
            info.genotype = source.get("genotype")[0]

        # `.dbxrefs` is a list populated from any PROJECT or DBLINK.
        # If `.dbxrefs` has information, extract `BioProject` and `BioSample`.
        if len(seq_record.dbxrefs) != 0:
            for dbxref in seq_record.dbxrefs:
                # Extract `BioProject` and `BioSample`.
                if "BioProject" in dbxref:
                    info.bioproject = dbxref.split(':')[1]
                if "BioSample" in dbxref:
                    info.biosample = dbxref.split(':')[1]

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
                info.assem_method = (
                    seq_record.annotations["structured_comment"]
                                          ["Genome-Assembly-Data"]
                                          ["Assembly Method"])
            # Extract `Genome Representation`.
            if "Genome Representation" in (
                seq_record.annotations["structured_comment"]
                                      ["Genome-Assembly-Data"]):
                info.gen_represent = (
                    seq_record.annotations["structured_comment"]
                                          ["Genome-Assembly-Data"]
                                          ["Genome Representation"])
            # Extract `Expected Final Version`.
            if "Expected Final Version" in (
                seq_record.annotations["structured_comment"]
                                      ["Genome-Assembly-Data"]
            ):
                info.exp_final_ver = (
                    seq_record.annotations["structured_comment"]
                                          ["Genome-Assembly-Data"]
                                          ["Expected Final Version"]
                )
            # Extract `Genome Coverage`.
            if "Genome Coverage" in (
                seq_record.annotations["structured_comment"]
                                      ["Genome-Assembly-Data"]):
                info.gen_coverage = (
                    seq_record.annotations["structured_comment"]
                                          ["Genome-Assembly-Data"]
                                          ["Genome Coverage"])
            # Extract `Sequencing Technology`.
            if "Sequencing Technology" in (
                seq_record.annotations["structured_comment"]
                                      ["Genome-Assembly-Data"]):
                info.seq_technol = (
                    seq_record.annotations["structured_comment"]
                                          ["Genome-Assembly-Data"]
                                          ["Sequencing Technology"])
        # Check if the sequence has `structured_comment` and `Assembly-Data`.
        elif "structured_comment" in seq_record.annotations and (
            "Assembly-Data" in (
                seq_record.annotations["structured_comment"])):
            # Extract `Assembly Method`.
            if "Assembly Method" in (
                seq_record.annotations["structured_comment"]
                                      ["Assembly-Data"]):
                info.assem_method = (
                    seq_record.annotations["structured_comment"]
                                          ["Assembly-Data"]
                                          ["Assembly Method"])
            # Extract `Genome Representation`.
            if "Genome Representation" in (
                seq_record.annotations["structured_comment"]
                                      ["Genome-Assembly-Data"]):
                info.gen_represent = (
                    seq_record.annotations["structured_comment"]
                                          ["Genome-Assembly-Data"]
                                          ["Genome Representation"])
            # Extract `Expected Final Version`.
            if "Expected Final Version" in (
                seq_record.annotations["structured_comment"]
                                      ["Genome-Assembly-Data"]
            ):
                info.exp_final_ver = (
                    seq_record.annotations["structured_comment"]
                                          ["Genome-Assembly-Data"]
                                          ["Expected Final Version"]
                )
            # Extract `Genome Coverage`.
            if "Genome Coverage" in (
                seq_record.annotations["structured_comment"]
                                      ["Assembly-Data"]):
                info.gen_coverage = (
                    seq_record.annotations["structured_comment"]
                                          ["Assembly-Data"]
                                          ["Genome Coverage"])
            # Extract `Sequencing Technology`.
            if "Sequencing Technology" in (
                seq_record.annotations["structured_comment"]
                                      ["Assembly-Data"]):
                info.seq_technol = (
                    seq_record.annotations["structured_comment"]
                                          ["Assembly-Data"]
                                          ["Sequencing Technology"])

        # Convert info class into a dictionary and append to records list.
        records.append(vars(info))

        # ################################################
        # Get references and make an independent database.
        # ################################################
        if 'references' in seq_record.annotations:
            ref_records.extend(get_references(seq_record))

    return (records, ref_records)


class References:
    """Class to store references from a GenBank record."""

    def __init__(
        self, accession='missing', reference_num='missing', location='missing',
        authors='missing', title='missing', journal='missing',
        medline_id='missing', pubmed_id='missing', comment='missing'
    ):
        self.accession = accession
        self.reference_num = reference_num
        self.location = location
        self.authors = authors
        self.title = title
        self.journal = journal
        self.medline_id = medline_id
        self.pubmed_id = pubmed_id
        self.comment = comment


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
        info = References()
        info.accession = record.id
        info.reference_num = str(i + 1)
        if reference.location:
            start = reference.location[0].start
            end = reference.location[0].end
            info.location = f"bases {start} to {end}"
        if reference.authors:
            info.authors = reference.authors
        if reference.title:
            info.title = reference.title
        if reference.journal:
            info.journal = reference.journal
        if reference.medline_id:
            info.medline_id = reference.medline_id
        if reference.pubmed_id:
            info.pubmed_id = reference.pubmed_id
        if reference.comment:
            info.comment = reference.comment
        references.append(vars(info))
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
                bioproject text,
                biosample text,
                assem_method text,
                gen_coverage text,
                seq_technol text,
                gen_represent text,
                exp_final_ver
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
                :bioproject, :biosample, :assem_method, :gen_coverage,
                :seq_technol, :gen_represent, :exp_final_ver
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
                'bioproject': raw_result['bioproject'],
                'biosample': raw_result['biosample'],
                'assem_method': raw_result['assem_method'],
                'gen_coverage': raw_result['gen_coverage'],
                'seq_technol': raw_result['seq_technol'],
                'gen_represent': raw_result['gen_represent'],
                'exp_final_ver': raw_result['exp_final_ver']
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
        info['bioproject'] = updated_result[20]
        info['biosample'] = updated_result[21]
        info['assem_method'] = updated_result[22]
        info['gen_coverage'] = updated_result[23]
        info['seq_technol'] = updated_result[24]
        info['gen_represent'] = updated_result[25]
        info['exp_final_ver'] = updated_result[26]
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
