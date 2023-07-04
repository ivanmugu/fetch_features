"""Class to make the graphic user interface of fetch_features.py.

This file is part of fetch_features
BSD 3-Clause License
Copyright (c) 2023, Ivan Munoz Gutierrez
"""
from pathlib import Path
from tkinter import filedialog
from io import StringIO
import sys
import threading

import customtkinter
from customtkinter import ThemeManager

from fetcher.user_input import is_valid_email
from fetcher.access_genbank import fetch_features_manager

# TODO: make a cancel button to kill downloading and restore all parameters to
# default.
# TODO: improve the reset button to clean the display textbox.


class App(customtkinter.CTk):
    """fetch_features GUI."""

    def __init__(self):
        super().__init__()
        # Parameters to run fetch_features
        self.input_path = None
        self.output_path = Path('.')
        self.type_of_ui = 'accession'
        self.access_biosample_from_accession = False
        self.save_as = 'csv'
        self.user_email = None
        self.done = False
        # Set layout parameters
        self.geometry('650x400')
        self.title('fetch_features.py')
        self.grid_columnconfigure(1, weight=1)
        self.grid_rowconfigure(0, weight=10)
        # Parameters to have correct theme
        self._placeholder_text_color = (
            ThemeManager.theme['CTkEntry']['placeholder_text_color']
        )
        self._text_color = (
            ThemeManager.theme['CTkEntry']['text_color']
        )

        # ########## #
        # Logo frame #
        # ########## #
        self.logo_frame = customtkinter.CTkFrame(self, corner_radius=0)
        self.logo_frame.grid(row=0, column=0, sticky='nsew')
        self.logo_frame.grid_rowconfigure(
            0, weight=1)
        # Fetch features label
        self.fetch_features_label = customtkinter.CTkLabel(
            self.logo_frame, text='Fetch\nfeatures',
            font=customtkinter.CTkFont(size=20, weight='bold')
        )
        self.fetch_features_label.grid(row=0, column=0, padx=10, pady=(10, 5))

        # ############# #
        # Display frame #
        # ############# #
        self.display_frame = customtkinter.CTkFrame(
            self, corner_radius=0,
            fg_color='transparent'
        )
        self.display_frame.grid(row=0, column=1, sticky='nsew')
        self.display_frame.grid_rowconfigure(0, weight=1)
        self.display_frame.grid_columnconfigure(0, weight=1)
        self.display_window = customtkinter.CTkTextbox(
            self.display_frame, wrap='word',
        )
        self.display_window.grid(
            row=0, column=0, padx=5, pady=5, sticky='nsew')
        self.display_window.focus_set()
        self.display_window.insert(
            'end', 'Welcome to fetch_features.py!\n'
        )
        self.display_window.configure(state='disabled')

        # ################ #
        # Navigation frame #
        # ################ #
        self.navigation_frame = customtkinter.CTkFrame(self, corner_radius=0)
        self.navigation_frame.grid(
            row=1, column=0, columnspan=2, sticky='nsew')
        self.navigation_frame.grid_columnconfigure(
            (0, 1, 2, 3), weight=1)
        # Input file label
        self.input_file_label = customtkinter.CTkLabel(
            self.navigation_frame, text='Input file:'
        )
        self.input_file_label.grid(row=0, column=0, padx=(10, 5), pady=(10, 5))
        # Input file entry
        self.input_file_entry = customtkinter.CTkEntry(
            self.navigation_frame,
            placeholder_text='Browse your input file',
        )
        self.input_file_entry.grid(
            row=0, column=1, columnspan=2, padx=5, pady=(10, 5), sticky='nsew'
        )
        self.input_file_entry.configure(state='disabled')
        # Input file button
        self.input_file_button = customtkinter.CTkButton(
            self.navigation_frame,
            text='Browse',
            command=self.get_list
        )
        self.input_file_button.grid(
            row=0, column=3, padx=(5, 10), pady=(10, 5))
        # Output folder label
        self.output_folder_label = customtkinter.CTkLabel(
            self.navigation_frame, text='Output folder:'
        )
        self.output_folder_label.grid(
            row=1, column=0, padx=(10, 5), pady=5
        )
        # Output folder entry
        self.output_folder_entry = customtkinter.CTkEntry(
            self.navigation_frame,
            placeholder_text=(
                'Browse your output folder (default current folder)'
            )
        )
        self.output_folder_entry.grid(
            row=1, column=1, columnspan=2, padx=5, pady=5, sticky='nsew'
        )
        self.output_folder_entry.configure(state='disabled')
        # Select folder button
        self.select_folder_button = customtkinter.CTkButton(
            self.navigation_frame,
            text='Browse',
            command=self.get_output_path
        )
        self.select_folder_button.grid(row=1, column=3, padx=(5, 10), pady=5)
        # Email Label
        self.email_label = customtkinter.CTkLabel(
            self.navigation_frame, text='Email:', anchor='w'
        )
        self.email_label.grid(row=2, column=0, padx=(10, 5), pady=5)
        # Email entry
        self.email_entry = customtkinter.CTkEntry(
            self.navigation_frame,
            placeholder_text='Provide your email to the NCBI',
        )
        self.email_entry.grid(
            row=2, column=1, columnspan=2, padx=5, pady=5, sticky='nsew',
        )
        # Type of unique identifiers message
        self.type_identifier_label = customtkinter.CTkLabel(
            self.navigation_frame, text='Type of UIs:'
        )
        self.type_identifier_label.grid(
            row=3, column=0, padx=(10, 5), pady=5)
        # Variable to store selected type of UI
        self.type_identifier_var = customtkinter.StringVar(self, 'Accession')
        # Type of unique identifiers menu
        self.type_identifier_menu = customtkinter.CTkOptionMenu(
            self.navigation_frame,
            values=['Accession', 'BioSample'],
            variable=self.type_identifier_var,
            command=lambda _: self.update_type_of_ui()
        )
        self.type_identifier_menu.grid(
            row=4, column=0, padx=(10, 5), pady=(5, 10))
        # Access BioSample from Accession number label
        self.access_bio_from_acc_label = customtkinter.CTkLabel(
            self.navigation_frame,
            text='Access BioSample from\naccession number:'
        )
        self.access_bio_from_acc_label.grid(row=3, column=1, padx=5, pady=5)
        # Access BioSample from Accession number variable
        self.access_bio_from_acc_var = customtkinter.StringVar(self, 'No')
        # Access BioSample from Accession number menu
        self.access_bio_from_acc_menu = customtkinter.CTkOptionMenu(
            self.navigation_frame,
            values=['No', 'Yes'],
            variable=self.access_bio_from_acc_var,
            command=lambda _: self.update_access_biosample_from_accession()
        )
        self.access_bio_from_acc_menu.grid(
            row=4, column=1, padx=5, pady=(5, 10)
        )
        # Save a copy of results as excel label
        self.save_as_label = customtkinter.CTkLabel(
            self.navigation_frame, text='Save\nresults as:'
        )
        self.save_as_label.grid(row=3, column=2, padx=5, pady=5)
        # Variable to store save excel menu
        self.save_as_var = customtkinter.StringVar(self, 'csv')
        # Save a copy of results as excel
        self.save_as_menu = customtkinter.CTkOptionMenu(
            self.navigation_frame,
            values=['csv', 'excel', 'csv and excel'],
            variable=self.save_as_var,
            command=lambda _: self.update_save_as()
        )
        self.save_as_menu.grid(row=4, column=2, padx=5, pady=(5, 10))
        # Reset button
        self.reset_button = customtkinter.CTkButton(
            self.navigation_frame,
            text='Reset',
            command=self.reset_parameters
        )
        self.reset_button.grid(row=3, column=3, padx=(5, 10), pady=5)
        # Run button. Use threading to run fetch features
        self.run_button = customtkinter.CTkButton(
            self.navigation_frame,
            text='Run',
            command=lambda: (
                threading.Thread(target=self.run_fetch_features).start()
            ),
            fg_color='#b300b3',
            hover_color='#800080'
        )
        self.run_button.grid(row=4, column=3, padx=(5, 10), pady=(5, 10))

        # Parameters to redirect stdout and update of the display window with
        # results from mining the NCBI database.
        self.output = StringIO()
        sys.stdout = self.output
        self.error = StringIO()
        sys.stderr = self.error
        self.update_textbox()

    # ======================================================================= #
    #                          FUNCTIONALITY                                  #
    # ======================================================================= #
    def get_list(self) -> None:
        """Get path to list."""
        path_list = filedialog.askopenfilename(
            initialdir='.',
            filetypes=(
                ('Text', '.txt'),
                ('Excel', '.xlsx')
            )
        )
        self.input_path = Path(path_list)
        # Get name of input file
        file_name = self.input_path.name
        # Show name of input file in entry box
        self.input_file_entry.configure(state='normal')
        self.input_file_entry.delete(0, 'end')
        self.input_file_entry.configure(
            placeholder_text=file_name,
            placeholder_text_color=self._text_color)
        self.input_file_entry.configure(state='disabled')

    def check_input_file(self) -> bool:
        """Check if user provided input file."""
        if self.input_path is None:
            self.display_window.configure(state='normal')
            self.display_window.insert(
                'end',
                "Select a list of accession or BioSample numbers"
            )
            self.display_window.configure(state='disabled')
            return False
        else:
            return True

    def get_output_path(self) -> None:
        """Get path to output directory"""
        path_output = filedialog.askdirectory(
            initialdir='.',
        )
        self.output_path = Path(path_output)
        # Get name of directory
        directory = f".../{self.output_path.name}"
        # Show path in entry box
        self.output_folder_entry.configure(state='normal')
        self.output_folder_entry.delete(0, 'end')
        self.output_folder_entry.configure(
            placeholder_text=directory,
            placeholder_text_color=self._text_color
        )
        self.output_folder_entry.configure(state='disabled')

    def get_email(self) -> bool:
        """Get user's email."""
        email = self.email_entry.get()
        # Activate display window
        self.display_window.configure(state='normal')
        if len(email) == 0:
            self.display_window.insert('end', 'Provide your email.\n')
            self.display_window.configure(state='disabled')
            return False
        elif not is_valid_email(email):
            self.display_window.insert('end', 'Invalid email.\n')
            self.display_window.configure(state='disabled')
            return False
        else:
            self.user_email = email
        # Disable dispaly window
        self.display_window.configure(state='disabled')
        return True

    def textbox_has_text(self, textbox) -> bool:
        """Check if textbox has text."""
        # `end-1c` means the end excluding newline character.
        text = textbox.get('1.0', 'end-1c')
        # The text string is converted into a bool. An string with no character
        # is False.
        return bool(text.strip())

    def has_welcome_msg(self, textbox) -> bool:
        """Check if textbox has welcome message."""
        text = textbox.get('1.0', 'end')
        if 'Welcome' in text:
            return True
        else:
            return False

    def update_type_of_ui(self) -> None:
        if self.type_identifier_var.get() == 'Accession':
            self.type_of_ui = 'accession'
        else:
            self.type_of_ui = 'biosample'

    def update_access_biosample_from_accession(self) -> None:
        if self.access_bio_from_acc_var.get() == 'No':
            self.access_biosample_from_accession = False
        else:
            self.access_biosample_from_accession = True
        print("value of access_bio_from_acc_var =",
              self.access_biosample_from_accession)

    def update_save_as(self) -> None:
        if self.save_as_var.get() == 'csv':
            self.save_as = 'csv'
        elif self.save_as_var.get() == 'excel':
            self.save_as = 'excel'
        else:
            self.save_as = 'csv-excel'

    def reset_parameters(self) -> None:
        """Reset stored parameters to run fetch_features.py"""
        # Clean display
        self.display_window.configure(state='normal')
        self.display_window.delete('0.0', 'end')
        self.display_window.configure(state='disabled')
        # Reset input path
        self.input_path = None
        self.input_file_entry.delete(0, 'end')
        self.input_file_entry.configure(
            state='normal',
            placeholder_text='Browse your input file',
            placeholder_text_color=self._placeholder_text_color
        )
        self.input_file_entry._activate_placeholder()
        self.input_file_entry.configure(state='disabled')
        # Reset output folder entry
        self.output_path = Path('.')
        self.output_folder_entry.delete(0, 'end')
        self.output_folder_entry.configure(
            state='normal',
            placeholder_text='Browse your output folder',
            placeholder_text_color=self._placeholder_text_color
        )
        self.output_folder_entry._activate_placeholder()
        self.output_folder_entry.configure(state='disabled')
        # Reset user email and entry
        self.user_email = None
        self.email_entry.delete(0, 'end')
        # Reset type of unique identifier
        self.type_identifier_var.set('Accession')
        self.type_identifier_menu.configure(variable=self.type_identifier_var)
        # Reset save copy excel
        self.save_as_var.set('csv')
        self.save_as_menu.configure(variable=self.save_as_var)
        self.save_as = 'csv'

    def run_fetch_features(self) -> bool:
        """Mine GenBank database."""
        # If textbox has welcome message, clean it.
        if self.has_welcome_msg(self.display_window):
            self.display_window.configure(state='normal')
            self.display_window.delete('1.0', 'end')
        # Check email
        if not self.get_email():
            return False
        # Check input file
        if not self.check_input_file():
            return False
        # Clean display
        self.display_window.delete('1.0', 'end')
        # Disable display window
        self.display_window.configure(state='normal')
        # run fetch_features
        fetch_features_manager(
            infile=self.input_path,
            email_address=self.user_email,
            type_list=self.type_of_ui,
            output_folder=self.output_path,
            access_biosample_from_accession=(
                self.access_biosample_from_accession
            ),
            save_as=self.save_as
        )
        # print OK message
        self.display_window.insert(
            'end', (
                "Done!\nYou should have a results.csv and ref_results.csv " +
                f"files in the path: `{self.output_path}`.\nIf you " +
                "requested `--save-as-excel`, you will find the files with" +
                "`.xlsx` extention.")
        )
        self.display_window.see('end')
        self.display_window.configure(state='disabled')
        self.done = True
        return True

    def update_textbox(self) -> None:
        """Redirect stdout to update display window with results."""
        if self.done is True:
            return None
        output = self.output.getvalue()
        error = self.error.getvalue()
        if error:
            self.display_window.insert(
                'end', 'Something went wrong!\nCheck the error messages.\n\n'
            )
            self.display_window.see('end')
            self.display_window.insert('end', error)
            self.display_window.see('end')
            self.display_window.configure(state='disabled')
            return None
        self.display_window.configure(state='normal')
        self.display_window.insert('end', output)
        self.display_window.see('end')
        self.output.seek(0)
        self.output.truncate(0)
        # Call updtate constantly to check stdout and stderr
        self.after(500, self.update_textbox)

    def print_values_to_run_fetch_features(self) -> None:
        """Print values needed to run fetch features."""
        print('values for running fetch_features.py')
        print('Input paht:', self.input_path)
        print('output path:', self.output_path)
        print('type of ui:', self.type_of_ui)
        print('access bio from acc:', self.access_bio_from_acc)
        print('save results as:', self.save_as)
        print('user email:', self.user_email)


if __name__ == '__main__':
    app = App()
    app.mainloop()
