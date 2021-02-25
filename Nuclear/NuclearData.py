import pandas as pd
import matplotlib.pyplot as plt
import requests
import io
import math
from collections import namedtuple

class NuclearData:
    '''
    Class to perform various nuclear energy calculations either using A-level constants or accurate constants.
    Uses the AME2016 dataset from "https://www-nds.iaea.org/amdc/ame2016/mass16.txt".
    Keyword arguments:
    accurate -- boolean type to switch between sets of fundamental constants with differing accuracy, default = False
    '''
    # define a Nucleus Type
    Nucleus = namedtuple('Nucleus', ['A', 'Z'])

    # constructor
    def __init__(self, accurate=False):

        # class members
        self.accurate = accurate

        # Download the data from the server address specified by url
        url = "https://www-nds.iaea.org/amdc/ame2016/mass16.txt"
        # Make a request
        req = requests.get(url)
        if (req.status_code):
            # if the request is successful
            data = req.content.decode('ISO-8859-1')
        else:
            # terminate the app
            print("Could not download data...")
            exit(1)

        # Read the experimental data into a Pandas DataFrame.
        self.df = pd.read_fwf(io.StringIO(data),
            skiprows=39,
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

        # Need to replace values in the first two rows to agree with constants defined above
        # define physical constants 
        if self.accurate:
            self.AMU = 1.660539040E-27
            self.SPEEDOFLIGHT = 2.99792458E8
            self.ELECTRONCHARGE = 1.0000000983 * 1.602176634E-19
            self.AMUTOKEV = 931494.0038
            self.ELECTRONMASS = 5.48579909070E-4
            self.HYDROGENMASS = self.getAtomicMassAMU(self.Nucleus(1,1))
            # Proton mass needs correction for binding energy of electron - not actually used though
            #self.PROTONMASS = self.HYDROGENMASS - self.ELECTRONMASS + 13.59844e-3 / self.AMUTOKEV 
            self.PROTONMASS = self.HYDROGENMASS - self.ELECTRONMASS 
            self.NEUTRONMASS = self.getAtomicMassAMU(self.Nucleus(1,0))

        else:
            self.AMU = 1.661e-27
            self.SPEEDOFLIGHT = 3.00e8
            self.ELECTRONCHARGE = 1.60e-19
            self.AMUTOKEV = self.AMU * self.SPEEDOFLIGHT**2 / self.ELECTRONCHARGE * 1e-3
            self.ELECTRONMASS = 5.485e-4
            self.PROTONMASS = 1.0072
            self.HYDROGENMASS = self.ELECTRONMASS + self.PROTONMASS
            self.NEUTRONMASS = 1.0084
            # change proton and neutron masses in the dataframe since we are using simple constants
            self.df.loc[0, 'AtomicMass'] = self.NEUTRONMASS
            self.df.loc[1, 'AtomicMass'] = self.HYDROGENMASS

        # conversion to MeV
        self.AMUTOMEV = self.AMUTOKEV / 1000

        # calc nuclear mass in amu
        self.df['NuclearMassAMU'] = self.df['AtomicMass'] - self.df['Z'] * self.ELECTRONMASS

        # calc mass defect in amu - as defined in the AME2016 Paper
        self.df['MassDefectAMU'] = (self.df['Z'] * self.HYDROGENMASS + self.df['N'] * self.NEUTRONMASS) - self.df['AtomicMass']
        #self.df['MassDefectAMU'] = (self.df['Z'] * PROTONMASS + self.df['N'] * NEUTRONMASS) - self.df['NuclearMass']

    # find the specified nucleus
    def findNucleus(self, nucleus):
        '''
        Internal function that returns a dataframe object of the specified nucleus.
        Keyword arguments:
        nucleus -- the nucleus of interest defined by the type: Nucleus(A, Z)
        '''
        # perform search of dataframe
        result = self.df.query('A ==' + str(nucleus.A) + ' & Z ==' + str(nucleus.Z))
        if result.empty:    # nucleus not found, terminate the calculation
            print("The nucleus (A = " + str(nucleus.A) + ", Z = " + str(nucleus.Z) + ") could not be found.")
            exit(1)
        return result
    
    def getElementName(self, nucleus):
        '''
        Returns the element name of the specified nucleus.
        Keyword arguments:
        nucleus -- the nucleus of interest defined by the type: Nucleus(A, Z)
        '''
        # find the nucleus
        result = self.findNucleus(nucleus)
        return result['El'].values[0]

    # return the element proton number
    def getProtonNumber(self, nucleus):
        result = self.findNucleus(nucleus)
        return result['Z'].values[0]

    # return the element neutron number
    def getNeutronNumber(self, nucleus):
        result = self.findNucleus(nucleus)
        return result['N'].values[0]

    # return the element nucleon number
    def getNucleonNumber(self, nucleus):
        result = self.findNucleus(nucleus)
        return result['A'].values[0]

    # return the element atomic mass / u
    def getAtomicMassAMU(self, nucleus):
        result = self.findNucleus(nucleus)
        return result['AtomicMass'].values[0]

    # return the element atomic mass / kg
    def getAtomicMassKG(self, nucleus):
        return self.getAtomicMassAMU(nucleus) * self.AMU

    # return the element nuclear mass / u
    def getNuclearMassAMU(self, nucleus):
        result = self.findNucleus(nucleus)
        return result['NuclearMassAMU'].values[0]

    # return the element nuclear mass / kg
    def getNuclearMassKG(self, nucleus):
        return self.getNuclearMassAMU(nucleus) * self.AMU

    # return the element nuclear mass defect / u
    def getMassDefectAMU(self, nucleus):
        result = self.findNucleus(nucleus)
        return result['MassDefectAMU'].values[0]

    # return the element nuclear mass defect / kg
    def getMassDefectKG(self, nucleus):
        return self.getMassDefectAMU(nucleus) * self.AMU

    # return the element binding energy / keV
    def getBindingEnergyMEV(self, nucleus):
        return self.getMassDefectAMU(nucleus) * self.AMUTOKEV / 1e3

    # return the element binding energy / J
    def getBindingEnergyJ(self, nucleus):
        return self.getMassDefectAMU(nucleus) * self.AMU * self.SPEEDOFLIGHT**2

    # return the element binding energy per nucleon / keV
    def getBindingEnergyPerNucleonMEV(self, nucleus):
        return self.getBindingEnergyMEV(nucleus) / nucleus.A

    # return the element binding energy per nucleon / J
    def getBindingEnergyPerNucleonJ(self, nucleus):
        return self.getBindingEnergyJ(nucleus) / nucleus.A

    # print a formatted table of data in reduced units
    def printTable(self, nucleus):
        dash = '-' * 52
        # Print the table header
        print(dash)
        print('{:<36s}{:<14s}'.format('Quantity', 'Value'))
        print(dash)
        # Now the data
        print('{:<36s}{:<14s}'.format('Element', self.getElementName(nucleus)))
        print('{:<36s}{:<14d}'.format('Proton Number', self.getProtonNumber(nucleus)))
        print('{:<36s}{:<14d}'.format('Nucleon Number', self.getNucleonNumber(nucleus)))
        print('{:<36s}{:<14d}'.format('Neutron Number', self.getNeutronNumber(nucleus)))
        print('{:<36s}{:<14.8f}'.format('Atomic Mass / u', self.getAtomicMassAMU(nucleus)))
        print('{:<36s}{:<14.8f}'.format('Nuclear Mass / u', self.getNuclearMassAMU(nucleus)))
        print('{:<36s}{:<14.8f}'.format('Mass Defect / u', self.getMassDefectAMU(nucleus)))
        print('{:<36s}{:<14.8f}'.format('Binding Energy / MeV', self.getBindingEnergyMEV(nucleus)))
        print('{:<36s}{:<14.8f}'.format('Binding Energy per Nucleon / MeV', self.getBindingEnergyPerNucleonMEV(nucleus)))
        print(dash)

    # print a formatted table of data in SI units
    def printTableSI(self, nucleus):
        dash = '-' * 52
        # Print the table header
        print(dash)
        print('{:<36s}{:<14s}'.format('Quantity', 'Value'))
        print(dash)
        # Now the data
        print('{:<36s}{:<14s}'.format('Element', self.getElementName(nucleus)))
        print('{:<36s}{:<14d}'.format('Proton Number', self.getProtonNumber(nucleus)))
        print('{:<36s}{:<14d}'.format('Nucleon Number', self.getNucleonNumber(nucleus)))
        print('{:<36s}{:<14d}'.format('Neutron Number', self.getNeutronNumber(nucleus)))
        print('{:<36s}{:<14.8e}'.format('Atomic Mass / kg', self.getAtomicMassKG(nucleus)))
        print('{:<36s}{:<14.8e}'.format('Nuclear Mass / kg', self.getNuclearMassKG(nucleus)))
        print('{:<36s}{:<14.8e}'.format('Mass Defect / kg', self.getMassDefectKG(nucleus)))
        print('{:<36s}{:<14.8e}'.format('Binding Energy / J', self.getBindingEnergyJ(nucleus)))
        print('{:<36s}{:<14.8e}'.format('Binding Energy per Nucleon / J', self.getBindingEnergyPerNucleonJ(nucleus)))
        print(dash)

    # print a formatted table of data in both units
    def printTableAll(self, nucleus):
        dash = '-' * 52
        # Print the table header
        print(dash)
        print('{:<36s}{:<14s}'.format('Quantity', 'Value'))
        print(dash)
        # Now the data
        print('{:<36s}{:<14s}'.format('Element', self.getElementName(nucleus)))
        print('{:<36s}{:<14d}'.format('Proton Number', self.getProtonNumber(nucleus)))
        print('{:<36s}{:<14d}'.format('Nucleon Number', self.getNucleonNumber(nucleus)))
        print('{:<36s}{:<14d}'.format('Neutron Number', self.getNeutronNumber(nucleus)))
        print('{:<36s}{:<14.8f}'.format('Atomic Mass / u', self.getAtomicMassAMU(nucleus)))
        print('{:<36s}{:<14.8e}'.format('Atomic Mass / kg', self.getAtomicMassKG(nucleus)))
        print('{:<36s}{:<14.8f}'.format('Nuclear Mass / u', self.getNuclearMassAMU(nucleus)))
        print('{:<36s}{:<14.8e}'.format('Nuclear Mass / kg', self.getNuclearMassKG(nucleus)))
        print('{:<36s}{:<14.8f}'.format('Mass Defect / u', self.getMassDefectAMU(nucleus)))
        print('{:<36s}{:<14.8e}'.format('Mass Defect / kg', self.getMassDefectKG(nucleus)))
        print('{:<36s}{:<14.8f}'.format('Binding Energy / MeV', self.getBindingEnergyMEV(nucleus)))
        print('{:<36s}{:<14.8e}'.format('Binding Energy / J', self.getBindingEnergyJ(nucleus)))
        print('{:<36s}{:<14.8f}'.format('Binding Energy per Nucleon / MeV', self.getBindingEnergyPerNucleonMEV(nucleus)))
        print('{:<36s}{:<14.8e}'.format('Binding Energy per Nucleon / J', self.getBindingEnergyPerNucleonJ(nucleus)))
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
        print('{:<36s}{:<14.3f}'.format('amu to MeV', self.AMUTOKEV / 1000))
        print(dash)

    # print a graph of binding energy per nucleon / MeV
    # to do: add options for range/number of graphs
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

    # calc energy release in alpha decay
    def calcAlphaEnergyReleaseMEV(self, nucleus):
        # can work in atomic masses here...
        # get the atomic mass of an alpha
        alphaMass = self.getAtomicMassAMU(self.Nucleus(4, 2))

        # get the mass of the parent
        parentMass = self.getAtomicMassAMU(nucleus)

        # get the mass of the daughter
        daughterMass = self.getAtomicMassAMU(self.Nucleus(nucleus.A - 4, nucleus.Z - 2))

        # mass change
        deltaMass = (daughterMass + alphaMass) - parentMass

        # total energy release
        energyRelease = deltaMass * self.AMUTOKEV / 1000

        # return the value
        return energyRelease

    def calcAlphaEnergyReleaseJ(self, nucleus):
        # get energy release in MeV
        energyRelease = self.calcAlphaEnergyReleaseMEV(nucleus)

        # convert to amu
        energyRelease *= 1000 / self.AMUTOKEV 

        # convert to J
        energyRelease *= self.AMU * self.SPEEDOFLIGHT**2

        # return the value
        return energyRelease

    def calcAlphaEnergyMEV(self, nucleus):
        # get energy release in MeV
        energyRelease = abs(self.calcAlphaEnergyReleaseMEV(nucleus))

        # fraction going to alpha particle - see http://www.personal.soton.ac.uk/ab1u06/teaching/phys3002/course/07_alpha.pdf
        alphaEnergy = energyRelease * (nucleus.A - 4) / nucleus.A

        # return the value
        return alphaEnergy

    def calcAlphaEnergyJ(self, nucleus):
        # get alpha energy in MeV
        alphaEnergy = self.calcAlphaEnergyMEV(nucleus)

        # convert to amu
        alphaEnergy *= 1000 / self.AMUTOKEV 

        # convert to J
        alphaEnergy *= self.AMU * self.SPEEDOFLIGHT**2

        # return the value
        return alphaEnergy

    def calcAlphaSpeed(self, nucleus):
        # get alpha energy in J
        alphaEnergy = self.calcAlphaEnergyJ(nucleus)

        # get alpha mass in kg
        alphaMass = self.getNuclearMassKG(self.Nucleus(4, 2))

        # get alpha speed in m/s
        alphaSpeed = math.sqrt(2 * alphaEnergy / alphaMass)

        # return the value
        return alphaSpeed

    def calcBetaEnergyReleaseMEV(self, nucleus):
        # get the mass of the parent nucleus
        parentMass = self.getNuclearMassAMU(nucleus)

        # get the mass of the daughter nucleus
        daughterMass = self.getNuclearMassAMU(self.Nucleus(nucleus.A, nucleus.Z + 1))

        # mass change
        deltaMass = (daughterMass + self.ELECTRONMASS) - parentMass

        # total energy release
        energyRelease = deltaMass * self.AMUTOKEV / 1000

        # return the value
        return energyRelease

    def calcBetaEnergyReleaseJ(self, nucleus):
        # get energy release in MeV
        energyRelease = self.calcBetaEnergyReleaseMEV(nucleus)

        # convert to amu
        energyRelease *= 1000 / self.AMUTOKEV 

        # convert to J
        energyRelease *= self.AMU * self.SPEEDOFLIGHT**2

        # return the value
        return energyRelease

    def calcBetaEnergyMEV(self, nucleus):
        # get energy release in MeV
        energyRelease = abs(self.calcBetaEnergyReleaseMEV(nucleus))

        # fraction going to beta particle - see http://www.personal.soton.ac.uk/ab1u06/teaching/phys3002/course/07_alpha.pdf
        betaEnergy = energyRelease * (nucleus.A - self.ELECTRONMASS) / nucleus.A

        # return the value
        return betaEnergy

    def calcBetaEnergyJ(self, nucleus):
        # get beta energy in MeV
        betaEnergy = self.calcBetaEnergyMEV(nucleus)

        # convert to amu
        betaEnergy *= 1000 / self.AMUTOKEV 

        # convert to J
        betaEnergy *= self.AMU * self.SPEEDOFLIGHT**2

        # return the value
        return betaEnergy

    def calcBetaSpeed(self, nucleus):
        # get beta energy in J
        betaEnergy = self.calcBetaEnergyJ(nucleus)

        # get alpha mass in kg
        betaMass = self.ELECTRONMASS * self.AMU

        # get beta speed in m/s
        betaSpeed = math.sqrt(2 * betaEnergy / betaMass)

        # relativistic beta speed
        mc2 = betaMass * self.SPEEDOFLIGHT**2 
        betaSpeed = math.sqrt(self.SPEEDOFLIGHT**2 *(1 - (mc2/(betaEnergy + mc2))**2))

        # return the value
        return betaSpeed

    def calcFissionEnergyReleaseMEV(self, parentNuc, daughterNuc1, daughterNuc2):
        '''
        Returns the energy released in MeV for the specified fission process. The
        calculation is performed using the difference in binding energies between
        the parent and dauaghter nuclei. This yields the same result as calculating
        mass difference in products and reactants.
        Keyword arguments:
        parentNuc -- the parent nucleus defined by the type: Nucleus(A, Z)
        daughterNuc1 -- the daughter nucleus defined by the type: Nucleus(A, Z)
        daughterNuc2 -- the partner daughter nucleus defined by the type: Nucleus(A, Z)
        '''
        # perform sanity check to make sure that Z is conserved
        deltaZ = parentNuc.Z - (daughterNuc1.Z + daughterNuc2.Z)
        if deltaZ != 0:
            print("The proton number is not conserved in this process, exiting...")
            exit(1)
        energyRelease = self.getBindingEnergyMEV(parentNuc) - (self.getBindingEnergyMEV(daughterNuc1) + self.getBindingEnergyMEV(daughterNuc2))
        return energyRelease

    def calcFissionEnergyReleaseJ(self, parentNuc, daughterNuc1, daughterNuc2):
        '''
        Returns the energy released in J for the specified fission process.
        Keyword arguments:
        parentNuc -- the parent nucleus defined by the type: Nucleus(A, Z)
        daughterNuc1 -- the daughter nucleus defined by the type: Nucleus(A, Z)
        daughterNuc2 -- the partner daughter nucleus defined by the type: Nucleus(A, Z)
        '''
        energyRelease = self.calcFissionEnergyReleaseMEV(parentNuc, daughterNuc1, daughterNuc2) * self.AMU * self.SPEEDOFLIGHT**2 / self.AMUTOMEV
        return energyRelease

    def calcFusionEnergyReleaseMEV(self, reactants, products):
        '''
        Returns the energy released in MeV for the specified fusion process.
        Keyword arguments:
        reactants -- list of the reactant nuclei defined by the type: Nucleus(A, Z)
        products -- list of the product nuclei defined by the type: Nucleus(A, Z)
        '''
        # perform sanity check to make sure that A, Z are conserved
        reactantA = 0
        reactantZ = 0
        for reactant in reactants:
            reactantA += reactant.A
            reactantZ += reactant.Z

        productA = 0
        productZ = 0
        for product in products:
            productA += product.A
            productZ += product.Z

        if (reactantZ - productZ):
            print("The proton number is not conserved in this process, exiting...")
            exit(1)

        if (reactantA - productA):
            print("The nucleon number is not conserved in this process, exiting...")
            exit(1)

        # total up the masses of reactants
        reactantMass = 0
        for reactant in reactants:
            if (reactant == self.Nucleus(0,-1) or reactant == self.Nucleus(0,1)):  # electron/positron detected
                reactantMass += self.ELECTRONMASS
            else:
                reactantMass += self.getNuclearMassAMU(reactant)

        # total up the masses of products
        productMass = 0
        for product in products:
            if (product == self.Nucleus(0,-1) or product == self.Nucleus(0,1)):  # electron/positron detected
                productMass += self.ELECTRONMASS
            else:
                productMass += self.getNuclearMassAMU(product)

        # calculate mass change
        deltaMass = productMass - reactantMass
        #return energy release in MeV
        energyRelease = deltaMass * self.AMUTOMEV
        return energyRelease

    def calcFusionEnergyReleaseJ(self, reactants, products):
        '''
        Returns the energy released in J for the specified fusion process.
        Keyword arguments:
        reactants -- list of the reactant nuclei defined by the type: Nucleus(A, Z)
        products -- list of the product nuclei defined by the type: Nucleus(A, Z)
        '''
        energyRelease = self.calcFusionEnergyReleaseMEV(reactants, products) * self.AMU * self.SPEEDOFLIGHT**2 / self.AMUTOMEV
        return energyRelease

    def printAlphaSpeed(self, nucleus):
        # print the value
        print('{:<20s}{:<14.4e}'.format('Alpha Speed (m/s) =', self.calcAlphaSpeed(nucleus)))

    def printAlphaEnergyMeV(self, nucleus):
        # print the value
        print('{:<17s}{:<14.4f}'.format('Alpha KE (MeV) =', self.calcAlphaEnergyMEV(nucleus)))

    def printAlphaEnergyJ(self, nucleus):
        # print the value
        print('{:<15s}{:<14.4e}'.format('Alpha KE (J) =', self.calcAlphaEnergyJ(nucleus)))

    def printAlphaEnergyReleaseMeV(self, nucleus):
        # print the value
        print('{:<23s}{:<14.4f}'.format('Energy Release (MeV) =', self.calcAlphaEnergyReleaseMEV(nucleus)))

    def printAlphaEnergyReleaseJ(self, nucleus):
        # print the value
        print('{:<21s}{:<14.4e}'.format('Energy Release (J) =', self.calcAlphaEnergyReleaseJ(nucleus)))

    def printBetaSpeed(self, nucleus):
        # print the value
        print('{:<20s}{:<14.4e}'.format('Beta Speed (m/s) =', self.calcBetaSpeed(nucleus)))

    def printBetaEnergyMeV(self, nucleus):
        # print the value
        print('{:<17s}{:<14.4f}'.format('Beta KE (MeV) =', self.calcBetaEnergyMEV(nucleus)))

    def printBetaEnergyJ(self, nucleus):
        # print the value
        print('{:<15s}{:<14.4e}'.format('Beta KE (J) =', self.calcBetaEnergyJ(nucleus)))

    def printBetaEnergyReleaseMeV(self, nucleus):
        # print the value
        print('{:<23s}{:<14.4f}'.format('Energy Release (MeV) =', self.calcBetaEnergyReleaseMEV(nucleus)))

    def printBetaEnergyReleaseJ(self, nucleus):
        # print the value
        print('{:<21s}{:<14.4e}'.format('Energy Release (J) =', self.calcBetaEnergyReleaseJ(nucleus)))

    def printFissionEnergyReleaseJ(self, parentNuc, daughterNuc1, daughterNuc2):
        # print the value
        print('{:<21s}{:<14.4e}'.format('Energy Release (J) =', self.calcFissionEnergyReleaseJ(parentNuc, daughterNuc1, daughterNuc2)))

    def printFissionEnergyReleaseMeV(self, parentNuc, daughterNuc1, daughterNuc2):
        # print the value
        print('{:<23s}{:<14.4f}'.format('Energy Release (MeV) =', self.calcFissionEnergyReleaseMEV(parentNuc, daughterNuc1, daughterNuc2)))

    def printFusionEnergyReleaseJ(self, reactants, products):
        # print the value
        print('{:<21s}{:<14.4e}'.format('Energy Release (J) =', self.calcFusionEnergyReleaseMEV(reactants, products)))

    def printFusionEnergyReleaseMeV(self, reactants, products):
        # print the value
        print('{:<23s}{:<14.4f}'.format('Energy Release (MeV) =', self.calcFusionEnergyReleaseMEV(reactants, products)))

    def performTests(self):
        # perform unit tests
        if (data.accurate):
            # alpha energy
            # create a nucleus
            nucleus = Nucleus(A=228,Z=90)

            test = data.calcAlphaEnergyReleaseMEV(nucleus) 
            assert test == -5.520152576655016, "result should be -5.520152576655016"

            test = data.calcAlphaEnergyReleaseJ(nucleus) 
            assert test == -8.84426026909684e-13, "result should be -8.84426026909684e-13"

            test = data.calcAlphaEnergyMEV(nucleus) 
            assert test == 5.423307794608436, "result should be 5.423307794608436"

            test = data.calcAlphaEnergyJ(nucleus)
            assert test == 8.689097808235492e-13, "result should be 8.689097808235492e-13"

            test = data.calcAlphaSpeed(nucleus) 
            assert test == 16172086.447361581, "result should be 16172086.447361581"

            # beta energy
            # create a nucleus
            nucleus = Nucleus(A=60,Z=27)

            test = data.calcBetaEnergyReleaseMEV(nucleus) 
            assert test == -2.822809675553812, "result should be -2.822809675553812"

            test = data.calcBetaEnergyReleaseJ(nucleus) 
            assert test == -4.5226401107649964e-13, "result should be -4.5226401107649964e-13"

            test = data.calcBetaEnergyMEV(nucleus) 
            assert test == 2.8227838666092264, "result should be 2.8227838666092264"

            test = data.calcBetaEnergyJ(nucleus)
            assert test == 4.522598760273318e-13, "result should be 4.522598760273318e-13"

            test = data.calcBetaSpeed(nucleus)
            assert test == 296249796.25722736, "result should be 296249796.25722736"

            # fusion
            # create nuclei
            reactants = [Nucleus(A=1,Z=1),Nucleus(A=1,Z=1),Nucleus(A=1,Z=1), Nucleus(A=1,Z=1), Nucleus(A=0,Z=-1), Nucleus(A=0,Z=-1)]
            products = [Nucleus(A=4, Z=2)]
            test = data.calcFusionEnergyReleaseMEV(reactants, products)
            assert test == -26.730966831944606, "result should be -26.730966831944606"

            # fission
            # create nuclei
            parentNuc = Nucleus(A=238,Z=92)
            daughterNuc1 = Nucleus(A=95,Z=38)
            daughterNuc2 = Nucleus(A=140,Z=54)
            test = data.calcFissionEnergyReleaseMEV(parentNuc, daughterNuc1, daughterNuc2)
            assert test == -171.19983397432816, "result should be -171.19983397432816"

            print("All tests passed successfully")

# for debugging
if __name__ == '__main__':
    # create object
    data = NuclearData(accurate=True)

    # define nucleus type
    Nucleus = namedtuple('Nucleus', ['A', 'Z'])

#    data.printTable(Nucleus(A=8,Z=2))
#    data.printTable(Nucleus(A=228,Z=90))
#    data.printConstants()
#    help(NuclearData)
#    parentNuc = Nucleus(A=235,Z=92)
#    daughterNuc1 = Nucleus(A=144,Z=56)
#    daughterNuc2 = Nucleus(A=90,Z=36)
#    print(data.calcFissionEnergyReleaseMEV(parentNuc, daughterNuc1, daughterNuc2))
#    reactants = [Nucleus(A=2,Z=1), Nucleus(A=2,Z=1)]
#    products = [Nucleus(A=4,Z=2)]
#    print(data.calcFusionEnergyReleaseMEV(reactants, products))

    parentNuc = Nucleus(A=235,Z=92)
    daughterNuc1 = Nucleus(A=144,Z=55)
    daughterNuc2 = Nucleus(A=90,Z=37)
    data.printFissionEnergyReleaseMeV(parentNuc, daughterNuc1, daughterNuc2)
    print(data.calcFissionEnergyReleaseMEV(parentNuc, daughterNuc1, daughterNuc2))

    deltaMass = (2 * data.getAtomicMassAMU(Nucleus(1,0)) + data.getAtomicMassAMU(Nucleus(144,55)) + data.getAtomicMassAMU(Nucleus(90,37)))
    deltaMass -= data.getAtomicMassAMU(Nucleus(1,0)) + data.getAtomicMassAMU(Nucleus(235,92))
    deltaMass *= data.AMUTOMEV
    print(deltaMass)
