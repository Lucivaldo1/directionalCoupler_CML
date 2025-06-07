import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.constants import c

'''
@author: Lucivaldo Aguiar 
References:
https://optics.ansys.com/hc/en-us/articles/360042819293-Coupled-ring-resonator-filters
Integrated ring resonator, Dominik G. Rabus
Synthesis of a parallel-coupled ring-resonator filter, Andrea Melloni
'''
def FSR2graph (FSR):
    #this graphic is valid only and if only m1 - m2 = 1
    nm = 1e-9
    um = 1e-6
    FSR1 = np.linspace(0.1*nm,15*nm, 100)
    FSR2 = FSR/(np.abs(1-(FSR/FSR1)))
    plt.plot(FSR1/nm,FSR2/nm)
    plt.ylim(min(FSR1/nm), max(FSR2/nm))
    plt.xlim(min(FSR1/nm), max(FSR1/nm))
    plt.xlabel('FSR1 [nm]')
    plt.ylabel('FSR2 [nm]')
    plt.title(rf'$FSR_2$($FSR ={FSR/nm}, FSR_1$)')
    plt.grid()
    plt.show()
    
def ringData(x, N,ng, FSR):
    """
       Computes key parameters for a ring resonator based on input parameters.

       Inputs:
           x   : Relation FSR/B (Free Spectral Range to Bandwidth ratio)
           N   : Number of rings
           ng  : Group index of the waveguide
           FSR : Free Spectral Range

       Output:
           A dictionary containing calculated values for FSR, K, Q, g, B, and L.
       """
    dictionary ={
        'FSR' : [],
        'K' : [],
        'Q' : [],
        'g' : [],
        'B' : [],
        'L' : []
        }
    #determining FSRs and ring total length
    if N > 1:
        m1 = 3
        m2 = 2
        FSR1 = FSR/m1
        FSR2 = (m1/m2)*FSR1
        dictionary['FSR'].extend([FSR1,FSR2])
        #ring total length
        l1 = (1550e-9**2)/(ng*FSR1)
        l2 = (1550e-9**2)/(ng*FSR2)
        dictionary['L'].extend([l1,l2])
    else:
        dictionary['FSR'].append(FSR)
        l = (1550e-9**2)/(ng*FSR)
        dictionary['L'].append(l)
        
    for n in range(1,3):
        #Computing parameters based on Melloni's method
        g = np.sqrt(2) * np.sin((2*n-1)/(2*N) * np.pi)
        B = FSR / x
        Q = FSR/(B*g)
        K = ((np.pi**2)/(2*Q**2))*(np.sqrt(1 + (4*Q**2)/(np.pi**2)-1))

        dictionary['g'].append(g)
        dictionary['Q'].append(FSR/(B*g))    
        dictionary['K'].append(K)
        dictionary['B'].append(B)
    #Computing coupling coefficient for the middle ring
    dictionary['K'].append(np.sqrt(0.25)*(dictionary['K'][0])**2)
    return dictionary

def placeOna(icApi,name,x,y, inputs, nmbrofpoints, startFrequency, stopFrequency):
    '''
    places an Ona at x,y coordinates
    :param name:
    :param inputs:
    :param nmbrofpoints:
    :return:
    '''
    icApi.switchtolayout()
    icApi.addelement('Optical Network Analyzer')
    icApi.set('name', name)
    icApi.setposition(name, x, y)
    icApi.set('number of input ports', inputs)
    icApi.set('input parameter', 2)
    icApi.set('start frequency',startFrequency)
    icApi.set('stop frequency', stopFrequency)
    icApi.set('number of points', nmbrofpoints)

    return 0

def getDeltaL(FSR, wavelength_center, ng):
    '''
    returns Delta L from ng and wavelength center
    :param FSR: float
    :param wavelength_center: float
    :param ng: float
    :return:
    '''
    L = wavelength_center**2 / (FSR * ng)

    return L

def getLpi(wavelength_center, neff):
    '''
    returns Lpi from neff and wavelength center
    :param FSR: float
    :param wavelength_center: float
    :param ng: float
    :return:
    '''
    Lpi = (wavelength_center) / (2*neff)

    return Lpi


def getKappa(icApi,nmbrofpoints, startFrequency, stopFrequency, Lstart, Lstop, step, elementName, toggle):
    '''
    this funciont is used to calculate the coupling coefficient of siepic pdk device
    via interconnect utilizing optical power meters
    :param icApi: interconnect api object
    units of each parameter:
    [wg_width] = um
    [gap] = um
    [radius] = um
    toggle = 1 for pdk, 2 for ideal and self made device
    '''
    icApi.switchtolayout()
    icApi.addelement(elementName)
    icApi.set('name','dc')
    icApi.setposition('dc', 25, 200)
    ##ona##
    icApi.addelement('Optical Network Analyzer')
    icApi.set('name', 'benchmark ona')
    icApi.setposition('benchmark ona', 0, 0)
    icApi.set('number of input ports', 2)
    icApi.set('input parameter', 2)
    icApi.set('start frequency', startFrequency)
    icApi.set('stop frequency', stopFrequency)
    icApi.set('number of points', nmbrofpoints)
    ##

    #adding the optical power meters
    for i in range(3):
        icApi.addelement('Optical Power Meter')
        name = f'optical power meter {i+1}'
        icApi.set('name', name)
        icApi.setposition(name, 300, 200*i)
    if toggle == 1:
        #connections
        icApi.connect('benchmark ona', 'output', 'dc', 'opt_1')
        icApi.connect('benchmark ona', 'output', 'optical power meter 1', 'input')

        icApi.connect('dc', 'opt_3', 'optical power meter 2', 'input')
        icApi.connect('dc', 'opt_3', 'benchmark ona', 'input 1')

        icApi.connect('dc', 'opt_4', 'optical power meter 3', 'input')
        icApi.connect('dc', 'opt_4', 'benchmark ona', 'input 2')
    ## self made device
    elif toggle == 2:
        # connections
        icApi.connect('benchmark ona', 'output', 'dc', 'port 3')
        icApi.connect('benchmark ona', 'output', 'optical power meter 1', 'input')

        icApi.connect('dc', 'port 2', 'optical power meter 2', 'input')
        icApi.connect('dc', 'port 2', 'benchmark ona', 'input 1')

        icApi.connect('dc', 'port 4', 'optical power meter 3', 'input')
        icApi.connect('dc', 'port 4', 'benchmark ona', 'input 2')
    outputsOna = []
    kappa = []
    LC_values = np.arange(Lstart,Lstop,step)
    for i in range(len(LC_values)):
        icApi.switchtolayout()
        icApi.select('dc')
        icApi.set('coupling_length', LC_values[i])
        icApi.run()
        xdB = icApi.getresult('optical power meter 3', 'sum/power')
        kappa.append(10**(xdB/10))
        outputsOna.append(icApi.getresult('benchmark ona', 'input 2/mode 1/gain'))
    ##dataframe to show the results

    df = pd.DataFrame(
        {
            'Coupling Length (um)': LC_values/1e-6,
            'Coupling Coefficient (dB)': kappa,
        }
    )
    ## deleting
    icApi.switchtolayout()
    icApi.select('benchmark ona')
    for i in range(3):
        icApi.shiftselect(f'optical power meter {i+1}')
    icApi.shiftselect('dc')
    icApi.delete()

    return df, outputsOna

def MZILatticefilter(icApi, neff, ng, L, delayLengths, k,LC,filenames, name, nLattice, toggle):

    #toggle: 1 for ideal device, 2 for pdk and 3 for s-parameters device
    if (toggle == 1 or toggle == 2):
        if (toggle == 1):
            dc_name = 'Waveguide Coupler'
            wg_name = 'Straight Waveguide'
            portName = 'port '
        elif (toggle == 2):
            dc_name = 'ebeam_dc_te1550'
            wg_name = 'ebeam_wg_integral_1550'
            portName = 'opt_'
        icApi.switchtolayout()
        # Adding and positioning the dcs and wgs
        for i in range(nLattice):
            icApi.addelement(dc_name)
            name_dc = f'dc{i + 1}'
            icApi.set('name', name_dc)
            icApi.setposition(name_dc, 400 + 400 * i, 100)
            if (toggle ==1):
                icApi.set('coupling coefficient 1', k[i])
            elif (toggle == 2):
                icApi.set('coupling_length', LC[i])
    if (toggle == 3):
        wg_name = 'ebeam_wg_integral_1550'
        portName = 'port '
    #S parameters device
    if (toggle == 3):
        for i in range(nLattice):
            icApi.addelement('Optical N Port S-Parameter')
            icApi.set('load from file',1)
            icApi.set('s parameters filename', filenames[i])
            name_dc = f'dc{i + 1}'
            icApi.set('name', name_dc)
            icApi.setposition(name_dc, 400 + 400 * i, 100)
    lengths = []

    for i in range(0, len(delayLengths), 2):
        lengths.extend([
            L,
            delayLengths[i],
            delayLengths[i + 1],
        ])
    if (toggle == 1 or toggle == 2 or toggle == 3):
        for i in range((nLattice-1)*2):
            icApi.addelement(wg_name)
            name_wg = f'wg{i + 1}'
            icApi.set('name', name_wg)
            if (toggle == 1):
                icApi.set('effective index 1', neff)
                icApi.set('group index 1', ng)
                icApi.set('length', lengths[i])
            elif (toggle == 2 or toggle == 3):
                icApi.set('wg_length', lengths[i])
            grupo = i // 2
            x = 600 + 400 * grupo
            y = 60 if i % 2 == 0 else 160
            icApi.setposition(name_wg, x, y)

    if (toggle == 1 or toggle == 2):
        # connecting em all together
        for i in range(nLattice-1):
            # dcs to wgs
            icApi.connect(f'dc{i + 1}', f'{portName}3', f'wg{2 * i + 1}', 'port 1')
            icApi.connect(f'dc{i + 1}', f'{portName}4', f'wg{2 * i + 2}', 'port 1')
            # wgs to dcs
            icApi.connect(f'wg{2 * i + 1}', 'port 2', f'dc{i + 2}', f'{portName}1', )
            icApi.connect(f'wg{2 * i + 2}', 'port 2', f'dc{i + 2}', f'{portName}2', )
    elif (toggle == 3):
        # connecting em all together
        for i in range(nLattice-1):
            # dcs to wgs
            icApi.connect(f'dc{i + 1}', f'{portName}2', f'wg{2 * i + 1}', 'port 1')
            icApi.connect(f'dc{i + 1}', f'{portName}4', f'wg{2 * i + 2}', 'port 1')
            # wgs to dcs
            icApi.connect(f'wg{2 * i + 1}', 'port 2', f'dc{i + 2}', f'{portName}1', )
            icApi.connect(f'wg{2 * i + 2}', 'port 2', f'dc{i + 2}', f'{portName}3', )
    # criando compound element
    icApi.select('dc1')
    for i in range(1, nLattice+1, 1):
        icApi.shiftselect(f'dc{i}')
    for j in range(0, (nLattice-1)*2, 1):
        icApi.shiftselect(f'wg{j + 1}')

    icApi.createcompound()
    icApi.select('COMPOUND_1')
    icApi.set('name', name)
    icApi.addport(name, 'port 1', 'Bidirectional', 'Optical Signal', 'Left', 0.25)
    icApi.addport(name, 'port 2', 'Bidirectional', 'Optical Signal', 'Left', 0.75)

    icApi.addport(name, 'port 3', 'Bidirectional', 'Optical Signal', 'Right', 0.25)
    icApi.addport(name, 'port 4', 'Bidirectional', 'Optical Signal', 'Right', 0.75)

    icApi.groupscope(name)
    if (toggle == 1 or toggle == 2):
        icApi.connect('RELAY_1', 'port', 'dc1', f'{portName}1')
        icApi.connect('RELAY_2', 'port', 'dc1', f'{portName}2')
        icApi.connect('RELAY_3', 'port', f'dc{nLattice}', f'{portName}3')
        icApi.connect('RELAY_4', 'port', f'dc{nLattice}', f'{portName}4')
    elif(toggle == 3):
        icApi.connect('RELAY_1', 'port', 'dc1', f'{portName}1')
        icApi.connect('RELAY_2', 'port', 'dc1', f'{portName}3')
        icApi.connect('RELAY_3', 'port', f'dc{nLattice}', f'{portName}2')
        icApi.connect('RELAY_4', 'port', f'dc{nLattice}', f'{portName}4')

    icApi.refresh()
    return 0

def MZILatticefilterPDK(icApi, L, delayLengths, LC, name, nLattice):

    icApi.switchtolayout()
    # Adding and positioning the dcs and wgs
    for i in range(nLattice):
        icApi.addelement('ebeam_dc_te1550')
        name_dc = f'dc{i + 1}'
        icApi.set('name', name_dc)
        icApi.setposition(name_dc, 400 + 400 * i, 100)
        icApi.set('coupling_length', LC[i])

    lengths = []

    for i in range(0, len(delayLengths), 2):
        lengths.extend([
            L,
            delayLengths[i],
            delayLengths[i + 1],
        ])

    for i in range((nLattice-1)*2):
        icApi.addelement('ebeam_wg_integral_1550')
        name_wg = f'wg{i + 1}'
        icApi.set('name', name_wg)
        icApi.set('wg_length', lengths[i])

        grupo = i // 2
        x = 600 + 400 * grupo
        y = 60 if i % 2 == 0 else 160
        icApi.setposition(name_wg, x, y)

    # connecting em all together
    for i in range(nLattice-1):
        # dcs to wgs
        icApi.connect(f'dc{i + 1}', 'opt_3', f'wg{2 * i + 1}', 'port 1')
        icApi.connect(f'dc{i + 1}', 'opt_4', f'wg{2 * i + 2}', 'port 1')
        # wgs to dcs
        icApi.connect(f'wg{2 * i + 1}', 'port 2', f'dc{i + 2}', 'opt_1', )
        icApi.connect(f'wg{2 * i + 2}', 'port 2', f'dc{i + 2}', 'opt_2', )

    # criando compound element
    icApi.select('dc1')
    for i in range(1, nLattice+1, 1):
        icApi.shiftselect(f'dc{i}')
    for j in range(0, (nLattice-1)*2, 1):
        icApi.shiftselect(f'wg{j + 1}')

    icApi.createcompound()
    icApi.select('COMPOUND_1')
    icApi.set('name', name)
    icApi.addport(name, 'port 1', 'Bidirectional', 'Optical Signal', 'Left', 0.25)
    icApi.addport(name, 'port 2', 'Bidirectional', 'Optical Signal', 'Left', 0.75)

    icApi.addport(name, 'port 3', 'Bidirectional', 'Optical Signal', 'Right', 0.25)
    icApi.addport(name, 'port 4', 'Bidirectional', 'Optical Signal', 'Right', 0.75)

    icApi.groupscope(name)
    icApi.connect('RELAY_1', 'port', 'dc1', 'opt_1')
    icApi.connect('RELAY_2', 'port', 'dc1', 'opt_2')
    icApi.connect('RELAY_3', 'port', f'dc{nLattice}', 'opt_3')
    icApi.connect('RELAY_4', 'port', f'dc{nLattice}', 'opt_4')

    icApi.refresh()
    return 0