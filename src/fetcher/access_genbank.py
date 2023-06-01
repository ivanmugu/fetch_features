"""Functions to fetch features of GenBank accession or BioSample numbers.

BSD 3-Clause License
Copyright (c) 2023, Ivan Munoz Gutierrez
"""
import csv
from pathlib import Path
from typing import TextIO, IO
import urllib
import sys

from Bio import Entrez

from fetcher import database


# TODO: I need to fix the fetch_from_biosample function to work similar as the
# fetch_from_accession function.


def fetch_from_accession(
        results_writer: TextIO, references_writer: TextIO,
        submission_list: list
) -> None:
    """Fetch features from a list of accession numbers.

    Parameters
    ----------
    results_file : TextIO file object
        Opened file as csv.DictWriter to save the fetched features of UIDs.
    references_file : TextIO file object
        Opened file as csv.DictWriter to save the references of UIDs.
    submission_list : list
        List of batches of UIDs. An example of a submission_list with a batch
        size of three looks like the following:
        [['CP049609.1', 'CP028704.1', 'CP043542.1'],
         ['CP040107.1', 'CP041747.1', 'CP042638.1'],
         ['CP015023.1', 'CP049163.1', 'CP051714.1']]
    """
    # Iterate over list of batches of accession numbers to fetch gb sequences.
    end = 0
    for set_num, submission in enumerate(submission_list):
        start = end
        batch_size = len(submission)
        end = end + batch_size
        print(f"Going to download record {start + 1} to {end}\n")

        try:
            # db -> database, nuccore -> nuleotide, rettype -> retrieval type,
            # retmode -> determines the format of the return output.
            fetch_handle = Entrez.efetch(
                db="nuccore", rettype="gb", retmode="text", id=submission,
            )
        except urllib.error.HTTPError as err:
            if err.code == 400:
                print(f"{submission} has invalid nuccore UIDs")
                print(
                    "If you provided a single accession number, your " +
                    "results file is going to be empty."
                )
            else:
                print(f"An error ocurred with the internet:\n{err}")
            continue
        except urllib.error.URLError as err:
            if err.readon.errno == 8:
                print(f"An error occured:\n{err}")
                sys.exit(f"Maybe your internet is disconnected!")
            else:
                sys.exit(f"An error occured with the internet:\n{err}")

        # Parse the data fetched from NCBI.
        records, ref_records = database.parser(fetch_handle, set_num + 1)
        # Save the retrived data in the csv file.
        for record in records:
            results_writer.writerow(record)
        for ref_record in ref_records:
            references_writer.writerow(ref_record)

        # Close fetch_handle
        fetch_handle.close()


def fetch_from_biosample(
        results_file: TextIO, results_fields: list, ref_results_file: TextIO,
        ref_results_fields: list, submission_list: list
) -> None:
    """Fetch features from a list of BioSample numbers.

    Read a list of BioSample numbers to retrieve the features of all accession
    numbers associated to every BioSample number. A BioSample number can have
    one or more associated accession numbers. For example, SAMN07169263 has
    three accession numbers: CP049611.1, CP049610.1, and CP049609.1. The first
    two accession numbers belong to plasmids and the last one to a chromosome.
    In the case of SAMN07169263 this function will fetch features of the three
    above mentioned accession number.

    Parameters
    ----------
    results_file : TextIO file object
        Opened output file to save the fetched features of the accession
        numbers.
    results_fields : list
        Headers for the `results` table. This list is the information to fetch
        for every accession number associated to the BioSample number.
    ref_results_file : TextIO file object
        Opened output file to save the references of the accession numbers.
    ref_results_fields : list
        Headers for the `references_results` table. This list is the
        information to fetch for every accession number associated to the
        BioSample number.
    submission_list : list
        List of BioSample numbers requested by the user.

    Notes
    -----
    Some BioSample numbers are associated to accession numbers that does not
    have any relevant information as contigs of few hundreds of nucleotides.
    Also, in addition to updated accession numbers, some BioSample numbers have
    outdated accession numbers. To clean all the fetched information from NCBI
    this program uses the `clean_features` function provided in database.py.
    """
    # Create DictWriter
    results_writer = csv.DictWriter(results_file, results_fields)
    ref_results_writer = csv.DictWriter(ref_results_file, ref_results_fields)
    lenght_acc_list = len(submission_list)

    # Iterate over the list of BioSample numbers (submission_list).
    for query, submission in enumerate(submission_list):
        # Number to keep track set_number of sequences (query), it is
        # important in case the connection to NCBI is interrupted so we can
        # know where to continue downloading
        set_number = query + 1

        # Print download record.
        print(f"Going to download record {query + 1} of {lenght_acc_list}\n")
        print(f"Accessing BioSample number: {submission}")

        # Search for the BioSample accession number. We need usehistory
        # to get the QueryKey and the WebEnv which define our history
        # session and can be used to performe searches of data.
        search_handle = Entrez.esearch(
            db="nuccore", term=submission, usehistory="y"
        )

        # Copy information in computer memory.
        # search_results = Entrez.read(search_handle)
        search_results = use_entrez_read(search_handle)

        # Close the handle.
        search_handle.close()

        # Count the number of results (number of sequences).
        count = int(search_results["Count"])
        print(f"Number of requested sequences from BioSample: {count}")

        # Copy cookie "WebEnv" and query "QueryKey" from our history session.
        # WevEnv -> Web environment string returned from a previous ESearch,
        # EPost or ELink call; QueryKey -> Integer query key returned by a
        # previous ESearch, EPost or ELink call.
        webenv = search_results["WebEnv"]
        query_key = search_results["QueryKey"]

        # Number of sequences to be requested by batch.
        # A batch of 500 is the max that we can request.
        batch_size = 500

        # biosample_records will hold all the features of accession numbers
        # associated to the BioSample number that is being analyzed.
        # Additionaly, biosample_ref_recors will hold all the references
        # associated to each accession number.
        biosample_records = []
        biosample_ref_records = []

        # Fetch information from GenBank by batches.
        for start in range(0, count, batch_size):
            end = min(count, start + batch_size)
            # Print download batch record
            print(
                f"Going to download record {start + 1} to {end} " +
                f"from set_number {set_number}"
            )
            # Get information.
            # db -> database, nuccore -> nuleotide, rettype -> retrieval
            # type, retmode -> determines the format of the return output,
            # retstart -> sequential index of the first UID in the
            # retrieved set_number to be shown in the XML output, retmax ->
            # total number of UIDs from the retrieved set_number to be
            # shown in the XML output, idtype-> specifies the type of
            # identifier to return for sequence databases, acc -> accesion
            # number
            fetch_handle = Entrez.efetch(
                db="nuccore",
                rettype="gb",
                retmode="text",
                retstart=start,
                retmax=end,
                webenv=webenv,
                query_key=query_key,
                idtype="acc"
            )

            # Parse the data fetched from NCBI.
            records, ref_records = database.parser(fetch_handle, set_number)

            # Extend the records obtained in biosample_records for future
            # cleaning of the data, i.e to get the most updated data.
            biosample_records.extend(records)
            biosample_ref_records.extend(ref_records)

        # Use the clean_features function to clean data and obtain the
        # most updated information.
        # clean_features() returns a list of dictionaries
        updated_features = database.clean_features(biosample_records)

        # Get references from updated accession number.
        updated_references = database.clean_references(
            updated_features, biosample_ref_records
        )

        # Save the updated retrived data in the csv files.
        for updated_feature in updated_features:
            results_writer.writerow(updated_feature)
        for updated_ref in updated_references:
            ref_results_writer.writerow(updated_ref)

        print(
            "Number of sequences saved after processing: " +
            f"{len(updated_features)}\n"
        )

        # Close fetch_handle.
        fetch_handle.close()


def use_entrez_read(handle: IO) -> None:
    try:
        search_result = Entrez.read(handle)
        return search_result
    except RuntimeError:
        print('Your list has invalid UIDs.')


def fetch_features_manager(
        infile: Path, email_address: str, type_list: str = 'accession',
        output_folder: Path = Path('.'),
        access_biosample_from_accession: bool = False, save_as: str = 'csv',
        batch_size: int = 100
) -> None:
    """Fetch features from a list of accession or BioSample numbers.

    Read a txt or xlsx file with a list of UIDs. Then, convert the UIDs into a
    list of comma-separated UIDs to retrieve the information via Entrez.

    Parameters
    ----------
    infile : Path
        Path to input infile.
    email_address : string
        User's email address as requested by NCBI.
    type_list : string
        Type of UIDs in the input file, i.e. "accession" or "biosample".
    output_folder : Path
        Path to save result files.
    access_biosample_from_accession : bool
        Access the BioSample number link to an accession number and get the
        features of all the associated accession numbers.
    save_as : str
        Output file format. Options: 'csv', 'excel', or 'csv-excel'.
    batch_size : int
        Number of UIDs to be requested by batch. According to the documentation
        of BioPython, it is not recommended to make batches of more than 500
        UIDs.
    """
    # Provide email address to NCBI
    Entrez.email = email_address
    # Read infile and create a list of accession or BioSample numbers.
    list_accessions = database.make_uid_list(infile)
    print(
        f"\nRequested accession or BioSample numbers: {len(list_accessions)}\n"
    )

    # Open result files to write the fetched data in csv format.
    results_path = output_folder / "results.csv"
    ref_results_path = output_folder / "ref_results.csv"
    results = open(results_path, "w")
    ref_results = open(ref_results_path, "w")

    # results' field names for the csv table.
    fields = [
        "set_batch", "description", "accession", "size", "molecule",
        "mod_date", "topology", "mol_type", "organism", "strain",
        "isolation_source", "host", "plasmid", "country", "lat_lon",
        "collection_date", "note", "serovar", "collected_by", "genotype",
        "bioproject", "biosample", "assem_method", "gen_coverage",
        "seq_technol", "gen_represent", "exp_final_ver"
    ]
    # references' field names for the csv table.
    ref_fields = [
        "accession", "reference_num", "location", "authors", "title",
        "journal", "medline_id", "pubmed_id", "comment"
    ]
    # Create DictWriter for results.
    results_writer = csv.DictWriter(results, fields)
    references_writer = csv.DictWriter(ref_results, ref_fields)
    # Write headers into csv files.
    results_writer.writeheader()
    references_writer.writeheader()

    # ======================================================================= #
    # If user requested features from a list of accession numbers use the
    # `fetch_from_accession` function.
    # ======================================================================= #
    if type_list == 'accession' and (not access_biosample_from_accession):
        # Make batches of accession numbers.
        submission_list = database.make_uid_batches(
            uid_list=list_accessions, batch_size=batch_size
        )
        # Fetch features.
        print("\nFetching features\n")
        fetch_from_accession(
            results_writer=results_writer,
            # results_fields=fields,
            references_writer=references_writer,
            # ref_results_fields=ref_fields,
            submission_list=submission_list,
        )

    # ======================================================================= #
    # If user requested to access the BioSample numbers linked to accession
    # numbers, first use the `get_biosample_numbers` function to retrieve the
    # linked BioSample numbers. Then use the `fetch_from_biosample` function.
    # ======================================================================= #
    if type_list == 'accession' and access_biosample_from_accession:
        # Make batches of accession numbers.
        submission_list = database.make_uid_batches(
            uid_list=list_accessions, batch_size=batch_size
        )
        # Get the BioSample numbers associated to accession numbers.
        print("Retrieving BioSample numbers from list of accession numbers")
        list_biosamples = database.get_biosample_numbers(
            submission_list, email_address
        )
        # Fetch features.
        print("\nFetching features\n")
        fetch_from_biosample(
            results_file=results,
            results_fields=fields,
            ref_results_file=ref_results,
            ref_results_fields=ref_fields,
            submission_list=list_biosamples,
        )

    # ======================================================================= #
    # If user requested features from a BioSample list use the
    # `fetch_from_biosample` function.
    # ======================================================================= #
    if type_list == 'biosample':
        # Fetch features.
        print("\nFetching features\n")
        fetch_from_biosample(
            results_file=results,
            results_fields=fields,
            ref_results_file=ref_results,
            ref_results_fields=ref_fields,
            submission_list=list_accessions,
        )

    # Close result files.
    results.close()
    ref_results.close()

    # Save files in the requested format.
    database.save_results_as(
        results=results_path,
        ref_results=ref_results_path,
        save_as=save_as
    )
