#!/usr/bin/env python3

from planetary_model import Planetary_Model
from stellar_model import Stellar_Model

if __name__ == "__main__":

    star = Stellar_Model()
    planet = None
    stellar_models = []
    planets = []

    while True:

        # First ask the user for basic stellar info that they want to analyze
        name, parallax, max_wavelength, flux, metallicity = input(
            "Enter star name, Parallax, Peak Wavelength, Flux, Metallicity: "
        ).split()

        star.analyze_stellar_data(name, float(parallax), int(max_wavelength), float(flux), float(metallicity))
        #stellar_properties = star.return_stellar_data()

        # Then ask the user if this is a main sequence star. We ask this becuase some of the parameters to be calculated apply only to main sequence stars.
        # main_sequence = input("Is it Main Sequence star? (Y/N) ")
        # if main_sequence == "y":
        star.calculate_main_sequence_mass_radius_lifetime()

        # Dump all the stellar properties that was calculate on terminal and also return it as a list
        stellar_properties = star.return_stellar_data()

        # Ask the user if they want to save the stellar model
        # create_model = input("Do you want to create a stellar model for this analysis data? (Y/N) ")
        # if create_model == "y":
        stellar_models.append(star)

        # # Ask the user if they want to add a planet to the stellar model
        # add_planet = input("Add Planet? (Y/N) ")
        # if add_planet == "y":
        planet = Planetary_Model(star)

        # Enter Transit information to calculate planet properties
        transit_period = input("Enter transit period (days): ")
        planet.calculate_orbital_radius(int(transit_period))

        # Calculate Planetary Mass
        doppler_shift = input("Doppler Shift (nm): ")
        planet.calculate_mass(float(doppler_shift))

        # Calculate Planetary Radius
        brightness_drop = input("Brigthness Drop (%): ")
        planet.calculate_radius(float(brightness_drop))
        planet.calculate_density()

        # Calculate Equilibrium Temperature
        albedo = input("Planets's Albedo: ")
        planet.calculate_equilibrium_temperature(float(albedo))

        planet.get_planetary_data()
