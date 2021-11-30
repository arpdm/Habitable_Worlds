""""
    Stellar Model representing a star.
    Notes:
        - Flux Unit - W/m^2
"""

import numpy as np
from enum import Enum
from planetary_model import Planet

LY_IN_PARSEC = 3.26
METERS_IN_LY = 9.461e15
SOLAR_LUMINOSITY = 3.827e26
WIENS_DISPLACEMENT_NM_K = 2897768.5
MASS_CONSTANT = 0.2857
RADIUS_CONSTANT = 5800
LIFETIME_CONSTANT = 1e10
MILION_YEAR = 1e6

# Giant luminosity class identifier limits
GIANT_MIN_TEMP = 3000
GIANT_MAX_TEMP = 4000
GIANT_MIN_LUMINOSITY_SOLAR = 50
GIANT_MAX_LUMINOSITY_SOLAR = 1000

SUPER_GIANT_MIN_TEMP = 3000
SUPER_GIANT_MAX_TEMP = 8000
SUPER_GIANT_MIN_LUMINOSITY_SOLAR = 1e4
SUPER_GIANT_MAX_LUMINOSITY_SOLAR = 1e6

# Spectral Classes
O_CLASS_MAX_TEMP = int(6e4)
O_CLASS_MIN_TEMP = int(3e4)
B_CLASS_MAX_TEMP = int(3e4)
B_CLASS_MIN_TEMP = int(1e4)
A_CLASS_MAX_TEMP = int(1e4)
A_CLASS_MIN_TEMP = int(7500)
F_CLASS_MAX_TEMP = int(7500)
F_CLASS_MIN_TEMP = int(6000)
G_CLASS_MAX_TEMP = int(6000)
G_CLASS_MIN_TEMP = int(5000)
K_CLASS_MAX_TEMP = int(5000)
K_CLASS_MIN_TEMP = int(3500)
M_CLASS_MAX_TEMP = int(3500)
M_CLASS_MIN_TEMP = int(2500)
L_CLASS_MAX_TEMP = int(2500)
L_CLASS_MIN_TEMP = int(1200)
T_CLASS_MAX_TEMP = int(1200)
T_CLASS_MIN_TEMP = int(500)
Y_CLASS_MAX_TEMP = int(500)
Y_CLASS_MIN_TEMP = int(0)


class Spectral_Class(Enum):
    O = 0
    B = 1
    A = 2
    F = 3
    G = 4
    K = 5
    M = 6
    L = 7
    T = 8
    Y = 9


class Luminosity_Class(Enum):
    SUPER_GIANT = 0
    BRIGHT_GIANT = 1
    GIANT = 2
    MAIN_SEQUENCE = 3
    SUB_DWARF = 4
    WHITE_DWARF = 5
    UNKNOWN = 6


class EM_Wave(Enum):
    UV = 0
    VIOLET = 1
    BLUE = 2
    CYAN = 3
    GREEN = 4
    YELLOW = 5
    ORANGE = 6
    RED = 7
    IR = 8


class Stellar_Data_Analyzer:
    def analyze_stellar_data(self, name, parallax_arcseconds, peak_wavelength_nm, flux, metallicity):
        """Analyzes the stellar data and makes calculations about the star based on provided user info"""
        self.name = name
        self.parallax = parallax_arcseconds
        self.peak_weavelength_nm = peak_wavelength_nm
        self.flux = flux
        self.metallicity = metallicity
        self.spectral_class = None
        self.luminosity_class = None
        self.mass_solar = None
        self.radius_solar = None
        self.lifetime_myr = None
        self.em_wave = None
        self.habitablity_check = False

        # Calculated parameters
        (self.distance_ly, self.distance_m) = self.calculate_distance()
        (self.luminasity_w, self.luminosity_solar) = self.calculate_luminosity()
        self.temperature_k = self.calculate_temperature()

        # Spectral analysis
        self.identify_luminosity_classes()
        self.identify_EM_wave()
        self.identify_spectral_classes()
        self.is_star_worth_checking_for_habitability()

    def calculate_distance(self):
        """Distance in Light Years is calculated by dividing 3.26 by parallax in arc-seconds"""

        ly = round(LY_IN_PARSEC / self.parallax, 3)
        meters = ly * METERS_IN_LY
        return (ly, meters)

    def calculate_luminosity(self):
        """Luminosity is calculated using Flux and distance. L = Flux x 4 x PI x Distance^2"""

        l = self.flux * 4 * np.pi * pow(self.distance_m, 2)
        l_solar = round(l / SOLAR_LUMINOSITY, 3)
        return (l, l_solar)

    def calculate_temperature(self):
        """Calculate the stellar temperature using Wien's Law Lambda = b / T"""

        temperature_k = WIENS_DISPLACEMENT_NM_K / self.peak_weavelength_nm
        return int(temperature_k)

    def identify_EM_wave(self):
        """Identify the EM Wave type based on its length in nm"""

        if self.peak_weavelength_nm < 400:
            self.em_wave = EM_Wave.UV
        elif self.peak_weavelength_nm >= 400 and self.peak_weavelength_nm < 450:
            self.em_wave = EM_Wave.VIOLET
        elif self.peak_weavelength_nm >= 450 and self.peak_weavelength_nm < 475:
            self.em_wave = EM_Wave.BLUE
        elif self.peak_weavelength_nm >= 475 and self.peak_weavelength_nm < 495:
            self.em_wave = EM_Wave.CYAN
        elif self.peak_weavelength_nm >= 495 and self.peak_weavelength_nm < 570:
            self.em_wave = EM_Wave.GREEN
        elif self.peak_weavelength_nm >= 570 and self.peak_weavelength_nm < 590:
            self.em_wave = EM_Wave.YELLOW
        elif self.peak_weavelength_nm >= 590 and self.peak_weavelength_nm < 620:
            self.em_wave = EM_Wave.ORANGE
        elif self.peak_weavelength_nm >= 620 and self.peak_weavelength_nm < 750:
            self.em_wave = EM_Wave.RED
        elif self.peak_weavelength_nm >= 750:
            self.em_wave = EM_Wave.IR

    def identify_luminosity_classes(self):
        """Identify Luminosity class."""

        if (
            self.temperature_k > GIANT_MIN_TEMP
            and self.temperature_k < GIANT_MAX_TEMP
            and self.luminosity_solar > GIANT_MIN_LUMINOSITY_SOLAR
            and self.luminosity_solar < GIANT_MAX_LUMINOSITY_SOLAR
        ):
            self.luminosity_class = Luminosity_Class.GIANT
        elif (
            self.temperature_k > SUPER_GIANT_MIN_TEMP
            and self.temperature_k < SUPER_GIANT_MAX_TEMP
            and self.luminosity_solar > SUPER_GIANT_MIN_LUMINOSITY_SOLAR
            and self.luminosity_solar < SUPER_GIANT_MAX_LUMINOSITY_SOLAR
        ):
            self.luminosity_class = Luminosity_Class.SUPER_GIANT
        else:
            self.luminosity_class = Luminosity_Class.UNKNOWN

        theoretical_value = 0.0002 * self.temperature_k - 1.047
        print(f"Theoretical Value:{theoretical_value}")

    def identify_spectral_classes(self):
        """Identify Spectral Class"""
        print(self.temperature_k, F_CLASS_MIN_TEMP, F_CLASS_MAX_TEMP)
        if self.temperature_k in range(O_CLASS_MIN_TEMP, O_CLASS_MAX_TEMP):
            self.spectral_class = Spectral_Class.O
        elif self.temperature_k in range(B_CLASS_MIN_TEMP, B_CLASS_MAX_TEMP):
            self.spectral_class = Spectral_Class.B
        elif self.temperature_k in range(A_CLASS_MIN_TEMP, A_CLASS_MAX_TEMP):
            self.spectral_class = Spectral_Class.A
        elif self.temperature_k in range(F_CLASS_MIN_TEMP, F_CLASS_MAX_TEMP):
            self.spectral_class = Spectral_Class.F
        elif self.temperature_k in range(G_CLASS_MIN_TEMP, G_CLASS_MAX_TEMP):
            self.spectral_class = Spectral_Class.B
        elif self.temperature_k in range(K_CLASS_MIN_TEMP, K_CLASS_MAX_TEMP):
            self.spectral_class = Spectral_Class.K
        elif self.temperature_k in range(M_CLASS_MIN_TEMP, M_CLASS_MAX_TEMP):
            self.spectral_class = Spectral_Class.M
        elif self.temperature_k in range(L_CLASS_MIN_TEMP, L_CLASS_MAX_TEMP):
            self.spectral_class = Spectral_Class.L
        elif self.temperature_k in range(T_CLASS_MIN_TEMP, T_CLASS_MAX_TEMP):
            self.spectral_class = Spectral_Class.T
        elif self.temperature_k in range(Y_CLASS_MIN_TEMP, Y_CLASS_MAX_TEMP):
            self.spectral_class = Spectral_Class.Y

    def calculate_main_sequence_mass_radius_lifetime(self):
        """Calculate mass, radius, and lifetime of star if it is a main sequence star"""

        # Calculate stellar mass in solar masses
        self.mass_solar = pow(self.luminosity_solar, MASS_CONSTANT)

        # Calculate stellar radius in solar radius
        self.radius_solar = (pow(self.luminosity_solar, 0.5)) / pow((self.temperature_k / RADIUS_CONSTANT), 2)

        # Calculate stellar lifetime in million years (Ma)
        self.lifetime_myr = (LIFETIME_CONSTANT * pow(self.mass_solar, -2.5)) / MILION_YEAR

        return (self.mass_solar, self.radius_solar, self.lifetime_myr)

    def is_star_worth_checking_for_habitability(self):
        """Check if it is worth to check for habitable planets around the star"""

        if self.luminosity_class == Luminosity_Class.MAIN_SEQUENCE:
            if (
                self.spectral_class is Spectral_Class.A
                or self.spectral_class is Spectral_Class.F
                or self.spectral_class is Spectral_Class.G
            ):
                self.habitablity_check = True

    def return_stellar_data(self):
        """Dump all stellar model information based on calculations"""

        print(f"=========== {self.name}===========")
        print(f"Distance_LY: {self.distance_ly}")
        print(f"L_Solar: {self.luminosity_solar}")
        print(f"T_K: {self.temperature_k}")
        print(f"M_Solar: {self.mass_solar}")
        print(f"R_Solar: {self.radius_solar}")
        print(f"Lifetime_Ma: {self.lifetime_myr}")
        print(f"Luminosity Class: {self.luminosity_class.name}")
        print(f"Spectral Class: {self.spectral_class.name}")
        print(f"EM Wave Type: {self.em_wave.name}")
        print(f"Check for Habitable Worlds: {self.habitablity_check}")
