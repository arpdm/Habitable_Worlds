from enum import Enum
import numpy as np

from numpy import sqrt
from stellar_model import Stellar_Model

SECONDS_IN_A_DAY = 86400
SOLAR_MASS_G = 1.989e33
CM_IN_AU = 1.496e13
DAYS_IN_A_YEAR = 365
SPEED_OF_LIGHT = 3e8
REST_WAVELENGTH_NM = 656.3
EARTH_RADIUS_CM = 637000000
EARTH_MASS_G = 5.97e027
BOLTZMANN_CONSTANT = 5.67e-8
METERS_IN_LY = 9.461e15
SOLAR_LUMINOSITY = 3.827e26
METERS_IN_AU = 1.496e11


class Planet_Type(Enum):
    ICE_GIANT = 1
    GAS_GIANT = 2
    TERRESTRIAL = 3


class Planetary_Model:
    def __init__(self, star: Stellar_Model):
        """Initialize a new planet"""

        self.orbital_radius_au = None
        self.radius_earth = None
        self.mass_earth = None
        self.density_g_cm3 = None
        self.stellar_model = star
        self.radial_velocity_m_sec = None
        self.planet_type = None
        self.equilibrium_temperature = None

    def calculate_radius(self, brightness_drop_percent):
        """Calculate Radius of the Planet in Eart Radius units"""

        self.radius_earth = pow((brightness_drop_percent / 100), 1 / 2) * self.stellar_model.radius_solar * 109

    def calculate_mass(self, doppler_shift_nm):
        """Calculate Mass of the Planet in Earth Mass units"""

        self.radial_velocity_m_sec = (doppler_shift_nm / REST_WAVELENGTH_NM) * SPEED_OF_LIGHT
        temp_val = pow(self.radial_velocity_m_sec, 2) * self.orbital_radius_au * self.stellar_model.mass_solar
        self.mass_earth = 11.177 * pow(temp_val, 1 / 2)

    def calculate_density(self):
        """Calculate Density of the Planet in g/cm^3"""

        planet_radius_cm = self.radius_earth * EARTH_RADIUS_CM
        planetary_volume = (4 / 3) * np.pi * pow(planet_radius_cm, 3)
        planetary_mass_g = self.mass_earth * EARTH_MASS_G
        self.density_g_cm3 = planetary_mass_g / planetary_volume

        self.determine_planetary_type()

    def determine_planetary_type(self):
        """Identify tyeh type of the planet based on density"""

        if self.density_g_cm3 < 1:
            self.planet_type = Planet_Type.GAS_GIANT
        elif 2.5 > self.density_g_cm3 >= 1:
            self.planet_type = Planet_Type.ICE_GIANT
        elif self.density_g_cm3 >= 2.5:
            self.planet_type = Planet_Type.TERRESTRIAL

    def calculate_orbital_radius(self, transit_period_days):
        """Calculate Orbital Radius (AU)"""

        transit_period_yr = transit_period_days / DAYS_IN_A_YEAR
        temp_val = pow(transit_period_yr, 2) * self.stellar_model.mass_solar
        self.orbital_radius_au = pow(temp_val, 1 / 3)

    def calculate_equilibrium_temperature(self, albedo):
        """Calculate equilibrium temperature wehre no atmosphere is counted"""

        orbital_radius_m = self.orbital_radius_au * METERS_IN_AU
        luminosity = self.stellar_model.luminosity_solar * SOLAR_LUMINOSITY
        flux = luminosity / (4 * np.pi * pow(orbital_radius_m, 2))

        print(f"R: {orbital_radius_m} L: {luminosity} F: {flux}")
        
        temp_value = (flux * (1 - albedo)) / (BOLTZMANN_CONSTANT * 4)
        self.equilibrium_temperature = pow(temp_value, 0.25)

    def get_planetary_data(self):

        print(f"=========== PLANET X ===========")
        print(f"Orbital Radius (AU): {self.orbital_radius_au}")
        print(f"Planetary Mass (M_E): {self.mass_earth}")
        print(f"Planetary Radius (R_E): {self.radius_earth}")
        print(f"Planetary Density (g/cm^3): {self.density_g_cm3}")
        print(f"Planetary Type : {self.planet_type}")
        print(f"Equilibrium Temperature (K): {self.equilibrium_temperature}")
