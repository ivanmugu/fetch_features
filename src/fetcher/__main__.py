"""Fetch_features main function.

This file is part of fetch_features
BSD 3-Clause License
Copyright (c) 2023, Ivan Munoz Gutierrez
"""
import sys
import textwrap

from fetcher.access_genbank import fetch_features_manager
from fetcher.user_input import parse_comm_line
from fetcher import gui


def main():
    # Get user input
    user_input = parse_comm_line()
    # Check if user want to activate the GUI.
    if user_input.use_gui is True:
        app = gui.App()
        app.mainloop()
        sys.exit(0)
    # Fetch information.
    fetch_features_manager(
        infile=user_input.infile,
        email_address=user_input.email,
        type_list=user_input.uid_type,
        output_folder=user_input.output_folder,
        access_biosample_from_accession=(
            user_input.access_biosample_from_accession
        ),
        save_as=user_input.save_as
    )
    # print OK message
    print(textwrap.dedent(f"""
        \nDone!
        You should have a results.csv and ref_results.csv files in the
        path: `{user_input.output_folder}`.
        If you requested `--save-as` `excel`, you will find the result files
        with `.xlsx` extention.
        """))

    sys.exit(0)


if __name__ == "__main__":
    main()
