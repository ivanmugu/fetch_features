"""File with a class and functions to interact with the command line.

This file is part of fetch_features
BSD 3-Clause License
Copyright (c) 2023, Ivan Munoz Gutierrez
"""
import sys
import argparse
from argparse import Namespace
import textwrap
from pathlib import Path
import re
from typing import Union
import pkg_resources

# Description message.
description_msg = r"""
  __      _       _        __            _
 / _| ___| |_ ___| |__    / _| ___  __ _| |_ _   _ _ __ ___  ___
| |_ / _ \ __/ __| '_ \  | |_ / _ \/ _` | __| | | | '__/ _ \/ __|
|  _|  __/ || (__| | | | |  _|  __/ (_| | |_| |_| | | |  __/\__ \
|_|  \___|\__\___|_| |_| |_|  \___|\__,_|\__|\__,_|_|  \___||___/

Fetch features from a list of accession or BioSample numbers.
"""

# Epilog message for command line help.
epilog_msg = r"""
Usage examples:
1. The simplest command. The output is in the current working directory.
$ fetcher -i path/to/list.txt -t accession -e email@address.com

2. In this example, the output is in your Documents.
$ fetcher -i path/to/list.txt -t accession -e email@address.com -o ~/Documents

3. If you prefer the GUI version.
$ fetcher --gui
"""


class UserInput:
    """Class to save user input via the commmand line."""
    # TODO: include a `type` option for the type of identifiers in the list.

    def __init__(
            self,
            infile: Union[Path, None] = None,
            extention_infile: Union[str, None] = None,
            uid_type: Union[str, None] = None,
            email: Union[str, None] = None,
            output_folder: Path = Path('.'),
            access_biosample_from_accession: bool = False,
            save_as: str = 'csv',
            use_gui: bool = False
    ):
        self.infile = infile
        self.extention_infile = extention_infile
        self.uid_type = uid_type
        self.email = email
        self.output_folder = output_folder
        self.access_biosample_from_accession = access_biosample_from_accession
        self.save_as = save_as
        self.use_gui = use_gui


def parse_comm_line() -> UserInput:
    """Parse and check command line arguments provided by user.

    Returns
    -------
    user_input : UserInput object
    """
    # Parse arguments and provide help.
    parser = argparse.ArgumentParser(
        add_help=False,
        prog='fetch_features',
        # usage=usage_msg,
        formatter_class=argparse.RawTextHelpFormatter,
        description=description_msg,
        epilog=textwrap.dedent(epilog_msg)
    )
    # Make argument groups.
    helper = parser.add_argument_group('Help')
    required = parser.add_argument_group('Required')
    optional = parser.add_argument_group('Optional')
    # ############# #
    # Help argument #
    # ############# #
    # Make help argument.
    helper.add_argument(
        '-h', '--help', action='help',
        help='Show this help message and exit.',
    )
    prog_version = pkg_resources.get_distribution('fetch_features').version
    helper.add_argument(
        '-v', '--version', action='version',
        version=f'%(prog)s {prog_version}',
        help="Show program's version number and exit"
    )

    # ################## #
    # Required arguments #
    # ################## #
    # Make input argument.
    required.add_argument(
        '-i', '--input',
        help=(
            'Path to input file with list of unique identifiers (UIDs).\n' +
            'The user should provide the list of UIDs in a txt or xlsx file.'
        ),
    )
    # Make unique identifier type argument.
    required.add_argument(
        '-t', '--type',
        help=(
            'Type of unique identifier: `accession` or `biosample`.\n' +
            'The `biosample` option fetches the information of the most \n' +
            'updated accession numbers associated with the BioSample number.'
        ),
    )
    # Make argument to provide email.
    required.add_argument(
        '-e', '--email',
        help='Provide your email address to the NCBI.',
    )
    # ################## #
    # Optional arguments #
    # ################## #
    # Make argument to provide output path.
    optional.add_argument(
        '-o', '--output',
        help=(
            "Path to output folder (default current working directory)."
        ),
    )
    # Make argument to specify how to save the results.
    optional.add_argument(
        '-s', '--save-as',
        help=(
            "Save results as: `csv`, `excel`, or `csv-excel` (default `csv`)."
        ),
    )
    # Make argument to check if user wants to access BioSample from accession.
    optional.add_argument(
        '--access-biosample-from-accession',
        action='store_true',
        help=(
            "If you provide a list of accession numbers, get features of\n" +
            "all related accession numbers that belong to the same BioSample."
        )
    )
    # Make argument to activate gui.
    optional.add_argument(
        '--gui',
        action='store_true',
        help="Activate GUI."
    )
    # ################################################################ #
    # Parse command line arguments and store info in a UserInput class #
    # ################################################################ #
    args = parser.parse_args()
    # print(args.__dict__)
    user_input = get_command_line_input(args)
    check_required_and_gui_arguments(user_input)

    return user_input


def get_command_line_input(command_line_input: Namespace) -> UserInput:
    """Get user input from command line."""
    # Initiate UserInput.
    user_input = UserInput()
    # Get input file.
    user_input.infile = get_infile(command_line_input)
    # Check input file extention.
    if user_input.infile:
        user_input.extention_infile = check_infile_extention(
            user_input.infile
        )
    # Get type of unique identifier.
    user_input.uid_type = get_uid_type(command_line_input)
    # Get user email.
    user_input.email = get_email(command_line_input)
    # Get output folder.
    user_input.output_folder = get_output_folder(command_line_input)
    # Check if user want to access BioSample numbers from accession numbers.
    user_input.access_biosample_from_accession = (
        check_access_biosample_from_accession(command_line_input)
    )
    # Get the format of output results.
    user_input.save_as = check_format_output_file(command_line_input)
    # Check if user wants to activate gui.
    user_input.use_gui = activate_gui(command_line_input)

    return user_input


def check_required_and_gui_arguments(user_input: UserInput) -> None:
    # If `use_gui` is selected, no mandatory arguments are needed.
    if (
        (user_input.use_gui and user_input.infile) or
        (user_input.use_gui and user_input.uid_type) or
        (user_input.use_gui and user_input.email)
    ):
        sys.exit(
            'Error: too many arguments. If you want to activate the ' +
            'GUI, only provide the `--gui` flag.'
        )
    # If `use_gui` is not selected, all mandatory arguments are needed.
    if (
        (not user_input.use_gui and not user_input.infile) or
        (not user_input.use_gui and not user_input.uid_type) or
        (not user_input.use_gui and not user_input.email)
    ):
        sys.exit(
            'Error: the `--input`, `--type`, and `--email` arguments ' +
            'are required.'
        )
    # If `--access-biosample-from-accession`, check that `--type` is accession
    if (
        (user_input.uid_type != 'accession') and
        (user_input.access_biosample_from_accession)
    ):
        sys.exit(
            'Error: when requesting `--access-biosample-from-accession` ' +
            'the `--type` argument must be accession.'
        )


def get_infile(comm_line_input: Namespace) -> Union[Path, None]:
    infile = comm_line_input.input
    # If user provided input file, check if valid.
    if infile:
        infile = Path(infile)
        if not infile.exists():
            sys.exit(f"Error: {infile} does not exist")
        if not infile.is_file():
            sys.exit(f"Error: {infile} is not a file")
        return infile
    else:
        return None


def check_infile_extention(infile: Path) -> None:
    # Get file's name
    name = infile.name
    # Get file's extention
    extention = name.split(".")[-1]
    # Check if file has correct extention.
    if extention == 'txt' or extention == 'xlsx':
        return None
    else:
        sys.exit(
            f"Error: `{name}` is not a valid input file.\nProvide a `txt` or" +
            " `xlsx` file."
        )


def get_uid_type(comm_line_input: Namespace) -> Union[str, None]:
    uid_type = comm_line_input.type
    if uid_type:
        # Check correct input
        if uid_type == "accession":
            return "accession"
        elif uid_type == "biosample":
            return "biosample"
        else:
            sys.exit(
                f"Error: `{uid_type}` is not valid. Enter `accession` or " +
                "`biosample`"
            )
    else:
        return None


def get_email(comm_line_input: Namespace) -> Union[str, None]:
    email = comm_line_input.email
    # Check if user provided email.
    if not email:
        return None
    # Check if email is valid.
    if is_valid_email(email):
        return email
    else:
        sys.exit(f'Error: invalid email `{email}`')


def is_valid_email(email: str) -> bool:
    regex = r'^[a-zA-Z0-9._%+-]+@[a-zA-Z0-9.-]+\.[a-zA-Z]{2,}$'
    return re.match(regex, email) is not None


def get_output_folder(comm_line_input: Namespace) -> Path:
    # Check if user provided output folder
    output_folder = comm_line_input.output
    if not output_folder:
        return Path.cwd()
    else:
        output_folder = Path(output_folder)
    # Check if output folder is valid.
    if not output_folder.exists():
        sys.exit(f"Error: {output_folder} does not exist")
    elif not output_folder.is_dir():
        sys.exit(f"Error: {output_folder} does not exist")
    else:
        return output_folder


def check_access_biosample_from_accession(comm_line_input: Namespace) -> bool:
    access_biosample = comm_line_input.access_biosample_from_accession
    if access_biosample:
        return True
    else:
        return False


def check_format_output_file(comm_line_input: Namespace) -> str:
    save_as = comm_line_input.save_as
    # Check correct input.
    if save_as is None:
        return 'csv'
    elif save_as == 'csv':
        return save_as
    elif save_as == 'excel':
        return save_as
    elif save_as == 'csv-excel':
        return save_as
    else:
        sys.exit(
            f"Error: {save_as} is not valid. Enter `csv`, `excel` or " +
            "`csv-excel`"
        )


def activate_gui(comm_line_input: Namespace) -> bool:
    gui = comm_line_input.gui
    if gui:
        return True
    else:
        return False
