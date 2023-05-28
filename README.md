# fetch features

Fetch information from a list of accession of BioSample numbers.

## Installation

Create a virtual environment and install `fetch features` using pip as follows:

```bash
pip install fetch-features
```

## Usage and options

To view all the options run:

```bash
fetch_features --help
```

Output:

```
usage: fetch_features [-h] [-v] [-i INPUT] [-t TYPE] [-e EMAIL] [-o OUTPUT] [-s SAVE_AS]
                      [--access-biosample-from-accession] [--gui]

  __      _       _        __            _
 / _| ___| |_ ___| |__    / _| ___  __ _| |_ _   _ _ __ ___  ___
| |_ / _ \ __/ __| '_ \  | |_ / _ \/ _` | __| | | | '__/ _ \/ __|
|  _|  __/ || (__| | | | |  _|  __/ (_| | |_| |_| | | |  __/\__ \
|_|  \___|\__\___|_| |_| |_|  \___|\__,_|\__|\__,_|_|  \___||___/

Fetch features from a list of accession or BioSample numbers.

Help:
  -h, --help            Show this help message and exit.
  -v, --version         Show program's version number and exit

Required:
  -i INPUT, --input INPUT
                        Path to input file with list of unique identifiers (UIDs).
                        The user should provide the list of UIDs in a txt or xlsx file.
  -t TYPE, --type TYPE  Type of unique identifier: `accession` or `biosample`.
                        The `biosample` option fetches the information of the most
                        updated accession numbers associated with the BioSample number.
  -e EMAIL, --email EMAIL
                        Provide your email address to the NCBI.

Optional:
  -o OUTPUT, --output OUTPUT
                        Path to output folder (default current working directory).
  -s SAVE_AS, --save-as SAVE_AS
                        Save results as: `csv`, `excel`, or `csv-excel` (default `csv`).
  --access-biosample-from-accession
                        If you provide a list of accession numbers, get features of
                        all related accession numbers that belong to the same BioSample.
  --gui                 Activate GUI.
```

## Usage examples

1. The simplest command. The output is in the current working directoy.

```bash
fetch_features -i path/to/list.txt -t accession -e email@address.com
```

2. In this example, the output is in your Documents.

```bash
fetch_features -i path/to/list.txt -t accession -e email@address.com -o ~/Documents
```

3. If you prefer the GUI version.

```bash
fetch_features --gui
```

<p align='center'>
  <img src=https://github.com/ivanmugu/fetch_features/blob/main/images/fetch_features_gui.png />
</p>
