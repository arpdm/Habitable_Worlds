#!/usr/bin/env python3

from stellar_data_analyzer import Stellar_Data_Analyzer

if __name__ == "__main__":

    analyzer = Stellar_Data_Analyzer()
    star_analyzed = False

    while True:
        name, parallax, max_wavelength, flux, metallicity = input(
            "Enter star name, Parallax, Peak Wavelength, Flux, Metallicity: "
        ).split()

        analyzer.analyze_stellar_data(name, float(parallax), int(max_wavelength), float(flux), float(metallicity))
        main_sequence = input("Is it Main Sequence star? (Yes/No)")
        if main_sequence == "Yes":
            analyzer.calculate_main_sequence_mass_radius_lifetime()
        analyzer.return_stellar_data()
