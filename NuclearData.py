import pandas as pd
import matplotlib.pyplot as plt
import requests
import io
import math

class NuclearData:
    # constructor
    def __init__(self, accurate=False):
        # class members
        self.accurate = accurate

        # define constants 
        if self.accurate:
            self.AMU = 1.660539040E-27
            self.ELECTRONMASS = 5.48579909070E-4
            self.PROTONMASS = 1.00727646693
            self.HYDROGENMASS = 1.00782503224
            self.NEUTRONMASS = 1.00866491582
            self.SPEEDOFLIGHT = 2.99792458E8
            self.ELECTRONCHARGE = 1.0000000983 * 1.602176634E-19
            self.AMUTOKEV = 931494.0038
        else:
            self.AMU = 1.661e-27
            self.ELECTRONMASS = 5.485e-4
            self.PROTONMASS = 1.0072
            self.HYDROGENMASS = self.ELECTRONMASS + self.PROTONMASS
            self.NEUTRONMASS = 1.0084
            self.SPEEDOFLIGHT = 3.00e8
            self.ELECTRONCHARGE = 1.60e-19
            self.AMUTOKEV = self.AMU * self.SPEEDOFLIGHT**2 / self.ELECTRONCHARGE * 1e-3

        # Download the data from the server address specified by url
        url = "https://www-nds.iaea.org/amdc/ame2016/mass16.txt"
        # Make a request
        req = requests.get(url)
        if (req.status_code):
            # if the request is successful
            data = req.content.decode('ISO-8859-1')
        else:
            # terminate
            exit(0)

        # Read the experimental data into a Pandas DataFrame.
        self.df = pd.read_fwf(io.StringIO(data),
            skiprows=41,
            header=None,
            widths=(1,3,5,5,5,1,3,4,1,13,11,11,9,1,2,11,9,1,16,11),
            usecols=(2,3,4,6,9,10,11,12,18,19),
            names=('N', 'Z', 'A', 'El', 'MassExcess', 'MassExcessUnc', 'BindingEnergy', 'BindingEnergyUnc', 'AtomicMass', 'AtomicMassUnc'),
            )

        # Extrapolated values are indicated by '#' in place of the decimal place, so
        # the columns won't be numeric. Globally replace the '#' with '.'
        self.df = self.df.replace({'#': '.'}, regex=True)

        # the masses have a space in their value and are un micro-u
        self.df['AtomicMass'] = self.df['AtomicMass'].str.replace(" ", "") 

        # ensure that columns are in numeric format after removing '#'
        for column in self.df.columns[4:8]:
            self.df[column] = pd.to_numeric(self.df[column], errors='coerce')

        # convert mass into u
        for column in self.df.columns[8:]:
            self.df[column] = pd.to_numeric(self.df[column], errors='coerce') * 1e-6

        # calc nuclear mass in amu
        self.df['NuclearMassAMU'] = self.df['AtomicMass'] - self.df['Z'] * self.ELECTRONMASS

        # calc mass defect in amu
        self.df['MassDefectAMU'] = (self.df['Z'] * self.HYDROGENMASS + self.df['N'] * self.NEUTRONMASS) - self.df['AtomicMass']
        #self.df['MassDefectAMU'] = (self.df['Z'] * PROTONMASS + self.df['N'] * NEUTRONMASS) - self.df['NuclearMass']

    # return the element name
    def getElementName(self, A, Z):
        result = self.df.query('A ==' + str(A) + ' & Z ==' + str(Z))
        return result['El'].values[0]

    # return the element proton number
    def getProtonNumber(self, A, Z):
        result = self.df.query('A ==' + str(A) + ' & Z ==' + str(Z))
        return result['Z'].values[0]

    # return the element neutron number
    def getNeutronNumber(self, A, Z):
        result = self.df.query('A ==' + str(A) + ' & Z ==' + str(Z))
        return result['N'].values[0]

    # return the element nucleon number
    def getNucleonNumber(self, A, Z):
        result = self.df.query('A ==' + str(A) + ' & Z ==' + str(Z))
        return result['A'].values[0]

    # return the element atomic mass / u
    def getAtomicMassAMU(self, A, Z):
        result = self.df.query('A ==' + str(A) + ' & Z ==' + str(Z))
        return result['AtomicMass'].values[0]

    # return the element atomic mass / kg
    def getAtomicMassKG(self, A, Z):
        return self.getAtomicMassAMU(A, Z) * self.AMU

    # return the element nuclear mass / u
    def getNuclearMassAMU(self, A, Z):
        result = self.df.query('A ==' + str(A) + ' & Z ==' + str(Z))
        return result['NuclearMassAMU'].values[0]

    # return the element nuclear mass / kg
    def getNuclearMassKG(self, A, Z):
        return self.getNuclearMassAMU(A, Z) * self.AMU

    # return the element nuclear mass defect / u
    def getMassDefectAMU(self, A, Z):
        result = self.df.query('A ==' + str(A) + ' & Z ==' + str(Z))
        return result['MassDefectAMU'].values[0]

    # return the element nuclear mass defect / kg
    def getMassDefectKG(self, A, Z):
        return self.getMassDefectAMU(A, Z) * self.AMU

    # return the element binding energy / keV
    def getBindingEnergyMEV(self, A, Z):
        return self.getMassDefectAMU(A, Z) * self.AMUTOKEV / 1e3

    # return the element binding energy / J
    def getBindingEnergyJ(self, A, Z):
        return self.getMassDefectAMU(A, Z) * self.AMU * self.SPEEDOFLIGHT**2

    # return the element binding energy per nucleon / keV
    def getBindingEnergyPerNucleonMEV(self, A, Z):
        return self.getBindingEnergyMEV(A, Z) / A

    # return the element binding energy per nucleon / J
    def getBindingEnergyPerNucleonJ(self, A, Z):
        return self.getBindingEnergyJ(A, Z) / A

    # print a formatted table of data in reduced units
    def printTable(self, A, Z):
        dash = '-' * 52
        # Print the table header
        print(dash)
        print('{:<36s}{:<14s}'.format('Quantity', 'Value'))
        print(dash)
        # Now the data
        print('{:<36s}{:<14s}'.format('Element', self.getElementName(A, Z)))
        print('{:<36s}{:<14d}'.format('Proton Number', self.getProtonNumber(A, Z)))
        print('{:<36s}{:<14d}'.format('Nucleon Number', self.getNucleonNumber(A, Z)))
        print('{:<36s}{:<14d}'.format('Neutron Number', self.getNeutronNumber(A, Z)))
        print('{:<36s}{:<14.8f}'.format('Atomic Mass / u', self.getAtomicMassAMU(A, Z)))
        print('{:<36s}{:<14.8f}'.format('Nuclear Mass / u', self.getNuclearMassAMU(A, Z)))
        print('{:<36s}{:<14.8f}'.format('Mass Defect / u', self.getMassDefectAMU(A, Z)))
        print('{:<36s}{:<14.8f}'.format('Binding Energy / MeV', self.getBindingEnergyMEV(A, Z)))
        print('{:<36s}{:<14.8f}'.format('Binding Energy per Nucleon / MeV', self.getBindingEnergyPerNucleonMEV(A, Z)))
        print(dash)

    # print a formatted table of data in SI units
    def printTableSI(self, A, Z):
        dash = '-' * 52
        # Print the table header
        print(dash)
        print('{:<36s}{:<14s}'.format('Quantity', 'Value'))
        print(dash)
        # Now the data
        print('{:<36s}{:<14s}'.format('Element', self.getElementName(A, Z)))
        print('{:<36s}{:<14d}'.format('Proton Number', self.getProtonNumber(A, Z)))
        print('{:<36s}{:<14d}'.format('Nucleon Number', self.getNucleonNumber(A, Z)))
        print('{:<36s}{:<14d}'.format('Neutron Number', self.getNeutronNumber(A, Z)))
        print('{:<36s}{:<14.8e}'.format('Atomic Mass / kg', self.getAtomicMassKG(A, Z)))
        print('{:<36s}{:<14.8e}'.format('Nuclear Mass / kg', self.getNuclearMassKG(A, Z)))
        print('{:<36s}{:<14.8e}'.format('Mass Defect / kg', self.getMassDefectKG(A, Z)))
        print('{:<36s}{:<14.8e}'.format('Binding Energy / J', self.getBindingEnergyJ(A, Z)))
        print('{:<36s}{:<14.8e}'.format('Binding Energy per Nucleon / J', self.getBindingEnergyPerNucleonJ(A, Z)))
        print(dash)

    # print a formatted table of data in both units
    def printTableAll(self, A, Z):
        dash = '-' * 52
        # Print the table header
        print(dash)
        print('{:<36s}{:<14s}'.format('Quantity', 'Value'))
        print(dash)
        # Now the data
        print('{:<36s}{:<14s}'.format('Element', self.getElementName(A, Z)))
        print('{:<36s}{:<14d}'.format('Proton Number', self.getProtonNumber(A, Z)))
        print('{:<36s}{:<14d}'.format('Nucleon Number', self.getNucleonNumber(A, Z)))
        print('{:<36s}{:<14d}'.format('Neutron Number', self.getNeutronNumber(A, Z)))
        print('{:<36s}{:<14.8f}'.format('Atomic Mass / u', self.getAtomicMassAMU(A, Z)))
        print('{:<36s}{:<14.8e}'.format('Atomic Mass / kg', self.getAtomicMassKG(A, Z)))
        print('{:<36s}{:<14.8f}'.format('Nuclear Mass / u', self.getNuclearMassAMU(A, Z)))
        print('{:<36s}{:<14.8e}'.format('Nuclear Mass / kg', self.getNuclearMassKG(A, Z)))
        print('{:<36s}{:<14.8f}'.format('Mass Defect / u', self.getMassDefectAMU(A, Z)))
        print('{:<36s}{:<14.8e}'.format('Mass Defect / kg', self.getMassDefectKG(A, Z)))
        print('{:<36s}{:<14.8f}'.format('Binding Energy / MeV', self.getBindingEnergyMEV(A, Z)))
        print('{:<36s}{:<14.8e}'.format('Binding Energy / J', self.getBindingEnergyJ(A, Z)))
        print('{:<36s}{:<14.8f}'.format('Binding Energy per Nucleon / MeV', self.getBindingEnergyPerNucleonMEV(A, Z)))
        print('{:<36s}{:<14.8e}'.format('Binding Energy per Nucleon / J', self.getBindingEnergyPerNucleonJ(A, Z)))
        print(dash)

    # print a formatted table of data in both units
    def printConstants(self):
        dash = '-' * 52
        # Print the table header
        print(dash)
        print('{:<36s}{:<14s}'.format('Quantity', 'Value'))
        print(dash)
        # Now the data
        print('{:<36s}{:<14.3e}'.format('AMU', self.AMU))
        print('{:<36s}{:<14.3e}'.format('mass of electron /u', self.ELECTRONMASS))
        print('{:<36s}{:<14.4f}'.format('mass of proton / u', self.PROTONMASS))
        print('{:<36s}{:<14.4f}'.format('mass of neutron / u', self.NEUTRONMASS))
        print('{:<36s}{:<14.2e}'.format('speed of light / m/s', self.SPEEDOFLIGHT))
        print('{:<36s}{:<14.2e}'.format('electron charge / C', self.ELECTRONCHARGE))
        print(dash)

    # print a graph of binding energy per nucleon / MeV
    def printGraph(self):
        # group the data by nucleon number, A 
        gdf = self.df.groupby('A')
        # Find the rows of the grouped DataFrame with the maximum binding energy.
        maxBindingEnergy = gdf.apply(lambda t: t[t.BindingEnergy==t.BindingEnergy.max()])
        # convert binding energy to MeV
        maxBindingEnergy['BindingEnergy'] /= 1e3
        # add datacolumn for Total binding energy 
        maxBindingEnergy['TotalBindingEnergy'] = maxBindingEnergy['BindingEnergy'] * maxBindingEnergy['A']

        # Generate a plot of the experimental values.
        fig, ax1 = plt.subplots()
        bindingenergy = ax1.plot(maxBindingEnergy['A'], maxBindingEnergy['BindingEnergy'], 'o',
            alpha=0.7, ls='-', lw=1, c='b', ms=2, label='Binding Energy per Nucleon')
        ax1.set_xlabel(r'$A = N + Z$')
        ax1.set_ylabel(r'$Binding \ Energy \ per \ Nucleon \,/\mathrm{MeV}$')
        
        # instantiate a second axes that shares the same x-axis
        ax2 = ax1.twinx()
        total = ax2.plot(maxBindingEnergy['A'], maxBindingEnergy['TotalBindingEnergy'], 'o',
            alpha=0.7, ls='-', lw=1, c='m', ms=2, label='Total Binding Energy')
        ax2.set_ylabel(r'$Total \ Binding \ Energy\,/\mathrm{MeV}$')

        # added the two lines to the legend
        lns = bindingenergy + total
        labs = [l.get_label() for l in lns]
        ax1.legend(lns, labs, loc="lower right")

        # show the plot
        plt.show()

    def calcAlphaEnergyReleaseMEV(self, A, Z):
        # can work in atomic masses here...
        # get the atomic mass of an alpha
        alphaMass = self.getAtomicMassAMU(4, 2)

        # get the mass of the parent
        parentMass =  self.getAtomicMassAMU(A, Z)

        # get the mass of the daughter
        daughterMass =  self.getAtomicMassAMU(A - 4, Z - 2)

        # mass change
        deltaMass = (daughterMass + alphaMass) - parentMass

        # total energy release
        energyRelease = deltaMass * self.AMUTOKEV / 1000

        # return the value
        return energyRelease

    def calcAlphaEnergyReleaseJ(self, A, Z):
        # get energy release in MeV
        energyRelease = self.calcAlphaEnergyReleaseMEV(A, Z)

        # convert to amu
        energyRelease *= 1000 / self.AMUTOKEV 

        # convert to J
        energyRelease *= self.AMU * self.SPEEDOFLIGHT**2

        # return the value
        return energyRelease

    def calcAlphaEnergyMEV(self, A, Z):
        # get energy release in MeV
        energyRelease = abs(self.calcAlphaEnergyReleaseMEV(A, Z))

        # fraction going to alpha particle - see http://www.personal.soton.ac.uk/ab1u06/teaching/phys3002/course/07_alpha.pdf
        alphaEnergy = energyRelease * (A - 4) / A

        # return the value
        return alphaEnergy

    def calcAlphaEnergyJ(self, A, Z):
        # get alpha energy in MeV
        alphaEnergy = self.calcAlphaEnergyMEV(A, Z)

        # convert to amu
        alphaEnergy *= 1000 / self.AMUTOKEV 

        # convert to J
        alphaEnergy *= self.AMU * self.SPEEDOFLIGHT**2

        # return the value
        return alphaEnergy

    def calcAlphaSpeed(self, A, Z):
        # get alpha energy in J
        alphaEnergy = self.calcAlphaEnergyJ(A, Z)

        # get alpha mass in kg
        alphaMass = self.getNuclearMassKG(4, 2)

        # get alpha speed in m/s
        alphaSpeed = math.sqrt(2 * alphaEnergy / alphaMass)

        # return the value
        return alphaSpeed

    def calcBetaEnergyReleaseMEV(self, A, Z):
        # get the mass of the parent nucleus
        parentMass = self.getNuclearMassAMU(A, Z)

        # get the mass of the daughter nucleus
        daughterMass = self.getNuclearMassAMU(A, Z + 1)

        # mass change
        deltaMass = (daughterMass + self.ELECTRONMASS) - parentMass

        # total energy release
        energyRelease = deltaMass * self.AMUTOKEV / 1000

        # return the value
        return energyRelease

    def calcBetaEnergyReleaseJ(self, A, Z):
        # get energy release in MeV
        energyRelease = self.calcBetaEnergyReleaseMEV(A, Z)

        # convert to amu
        energyRelease *= 1000 / self.AMUTOKEV 

        # convert to J
        energyRelease *= self.AMU * self.SPEEDOFLIGHT**2

        # return the value
        return energyRelease

    def calcBetaEnergyMEV(self, A, Z):
        # get energy release in MeV
        energyRelease = abs(self.calcBetaEnergyReleaseMEV(A, Z))

        # fraction going to beta particle - see http://www.personal.soton.ac.uk/ab1u06/teaching/phys3002/course/07_alpha.pdf
        betaEnergy = energyRelease * (A - self.ELECTRONMASS) / A

        # return the value
        return betaEnergy

    def calcBetaEnergyJ(self, A, Z):
        # get beta energy in MeV
        betaEnergy = self.calcBetaEnergyMEV(A, Z)

        # convert to amu
        betaEnergy *= 1000 / self.AMUTOKEV 

        # convert to J
        betaEnergy *= self.AMU * self.SPEEDOFLIGHT**2

        # return the value
        return betaEnergy

    def calcBetaSpeed(self, A, Z):
        # get beta energy in J
        betaEnergy = self.calcBetaEnergyJ(A, Z)

        # get alpha mass in kg
        betaMass = self.ELECTRONMASS * self.AMU

        # get beta speed in m/s
        betaSpeed = math.sqrt(2 * betaEnergy / betaMass)

        # relativistic beta speed
        mc2 = betaMass * self.SPEEDOFLIGHT**2 
        betaSpeed = math.sqrt(self.SPEEDOFLIGHT**2 *(1 - (mc2/(betaEnergy + mc2))**2))

        # return the value
        return betaSpeed

    def printAlphaSpeed(self, A, Z):
        # print the value
        print('{:<20s}{:<14.4e}'.format('Alpha Speed (m/s) =', self.calcAlphaSpeed(A, Z)))

    def printAlphaEnergyMeV(self, A, Z):
        # print the value
        print('{:<17s}{:<14.4f}'.format('Alpha KE (MeV) =', self.calcAlphaEnergyMEV(A, Z)))

    def printAlphaEnergyJ(self, A, Z):
        # print the value
        print('{:<15s}{:<14.4e}'.format('Alpha KE (J) =', self.calcAlphaEnergyJ(A, Z)))

    def printAlphaEnergyReleaseMeV(self, A, Z):
        # print the value
        print('{:<23s}{:<14.4f}'.format('Energy Release (MeV) =', self.calcAlphaEnergyReleaseMEV(A, Z)))

    def printAlphaEnergyReleaseJ(self, A, Z):
        # print the value
        print('{:<21s}{:<14.4e}'.format('Energy Release (J) =', self.calcAlphaEnergyReleaseJ(A, Z)))

    def printBetaSpeed(self, A, Z):
        # print the value
        print('{:<20s}{:<14.4e}'.format('Beta Speed (m/s) =', self.calcBetaSpeed(A, Z)))

    def printBetaEnergyMeV(self, A, Z):
        # print the value
        print('{:<17s}{:<14.4f}'.format('Beta KE (MeV) =', self.calcBetaEnergyMEV(A, Z)))

    def printBetaEnergyJ(self, A, Z):
        # print the value
        print('{:<15s}{:<14.4e}'.format('Beta KE (J) =', self.calcBetaEnergyJ(A, Z)))

    def printBetaEnergyReleaseMeV(self, A, Z):
        # print the value
        print('{:<23s}{:<14.4f}'.format('Energy Release (MeV) =', self.calcBetaEnergyReleaseMEV(A, Z)))

    def printBetaEnergyReleaseJ(self, A, Z):
        # print the value
        print('{:<21s}{:<14.4e}'.format('Energy Release (J) =', self.calcBetaEnergyReleaseJ(A, Z)))

# for debugging
if __name__ == '__main__':
    data = NuclearData(accurate=False)
    A = 228
    Z = 90
    data.printAlphaSpeed(A, Z)
    data.printAlphaEnergyJ(A, Z)
    data.printAlphaEnergyMeV(A, Z)
    data.printAlphaEnergyReleaseMeV(A, Z)
    data.printAlphaEnergyReleaseJ(A, Z)

    A = 60
    Z = 27
    data.printBetaSpeed(A, Z)
    data.printBetaEnergyJ(A, Z)
    data.printBetaEnergyMeV(A, Z)
    data.printBetaEnergyReleaseMeV(A, Z)
    data.printBetaEnergyReleaseJ(A, Z)

    data.printConstants()
