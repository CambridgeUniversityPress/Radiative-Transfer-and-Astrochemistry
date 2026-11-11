#
# Script 9-1 for the textbook "Foundations of Radiative Transfer & Astrochemistry"
# Shows how much of a layered planar slab of hydrogen is photoionised by a beam of
# ionising photons, assuming photoionisation equilibrium (see section 9.2 and
# Quick Exercise 9.2)
# Includes sliders that change the number of layers, the density and thickness of
# each layer, and the incident ionising photon flux
#
#

#import libraries
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from matplotlib.patches import Rectangle

#case B recombination coefficient in m^3/s
alphaB = 2.6e-19

#maximum number of layers
maxLayers = 4

#initial values. The first three layers are those of Quick Exercise 9.2
logFlux0 = 10.0                       #log10 ionising flux, m^-2 s^-1
nLayers0 = 3
logDens0 = [6.0, 8.0, 10.0, 4.0]      #log10 densities, m^-3
width0 = [0.3, 0.8, 1.0, 0.5]         #thicknesses, pc

#colours for ionised and neutral gas
ionColour = "orange"
neutralColour = "blue"

#parsec in m
pc = 3.086e16


def ionisedThickness(phi, n, D):
    #Return the ionised thickness of a layer (in pc) and the ionising flux
    #that passes through it, given incident flux phi (m^-2 s^-1), density n (m^-3)
    #and thickness D (pc)
    
    #photons per unit area per second needed to keep the whole layer ionised
    needed = n**2 * alphaB * D * pc
    if phi >= needed:
        #the whole layer is ionised and the leftover photons carry on
        return D, phi - needed
    else:
        #only part of the layer is ionised and no photons get through
        return phi / (n**2 * alphaB) / pc, 0.0


#set up the figure: schematic at the top, sliders underneath
fig, ax = plt.subplots(figsize=(10, 8))
plt.subplots_adjust(left=0.1, right=0.95, top=0.8, bottom=0.57)

axcolor = 'lightgoldenrodyellow'
sliderLeft = 0.3
sliderWidth = 0.55

#flux and number of layers sliders
axflux = plt.axes([sliderLeft, 0.46, sliderWidth, 0.025], facecolor=axcolor)
axnlay = plt.axes([sliderLeft, 0.42, sliderWidth, 0.025], facecolor=axcolor)
sflux = Slider(axflux, r'log$_{10}$ ionising flux, m$^{-2}$s$^{-1}$', 6.0, 18.0,
               valinit=logFlux0, valstep=0.1)
snlay = Slider(axnlay, 'Number of layers', 1, maxLayers, valinit=nLayers0, valstep=1)

#density and thickness sliders for each layer
sdens = []
swidth = []
for i in range(maxLayers):
    y = 0.36 - i * 0.08
    axd = plt.axes([sliderLeft, y, sliderWidth, 0.025], facecolor=axcolor)
    axw = plt.axes([sliderLeft, y - 0.035, sliderWidth, 0.025], facecolor=axcolor)
    sdens.append(Slider(axd, r'Layer %d: log$_{10}$ density, m$^{-3}$' % (i + 1),
                        2.0, 12.0, valinit=logDens0[i], valstep=0.1))
    swidth.append(Slider(axw, 'Layer %d: thickness, pc' % (i + 1),
                         0.05, 2.0, valinit=width0[i], valstep=0.05))


#draw the schematic for the current slider values
def update(val):
    ax.cla()
    nLayers = int(snlay.val)
    phi = 10**sflux.val

    #grey out the sliders of layers that are not in use
    for i in range(maxLayers):
        alpha = 1.0 if i < nLayers else 0.3
        for s in (sdens[i], swidth[i]):
            s.ax.set_alpha(alpha)
            s.label.set_alpha(alpha)
            s.valtext.set_alpha(alpha)
            s.poly.set_alpha(alpha)

    x = 0.0
    for i in range(nLayers):
        n = 10**sdens[i].val
        D = swidth[i].val
        Dion, phi = ionisedThickness(phi, n, D)

        #ionised part of the layer, then the neutral remainder
        ax.add_patch(Rectangle((x, 0), Dion, 1, facecolor=ionColour, edgecolor="none"))
        ax.add_patch(Rectangle((x + Dion, 0), D - Dion, 1, facecolor=neutralColour,
                               edgecolor="none"))
        #layer outline
        ax.add_patch(Rectangle((x, 0), D, 1, facecolor="none", edgecolor="black", lw=1.5))

        #label each layer with its density and ionised fraction
        label = "Layer %d\nn = $10^{%.1f}$ m$^{-3}$\n%.3g%% ionised" % (
            i + 1, sdens[i].val, 100.0 * Dion / D)
        ax.text(x + D / 2, 1.05, label, ha="center", va="bottom", fontsize=10)
        x = x + D

    #arrow showing the incident beam
    ax.annotate("", xy=(0, 0.5), xytext=(-0.12 * x, 0.5),
                arrowprops=dict(arrowstyle="-|>", lw=2, color="black"))
    ax.text(-0.075 * x, 0.58, "ionising\nphotons", ha="center", va="bottom", fontsize=10)

    #how many ionising photons make it all the way through
    ax.set_title("Ionising flux escaping the far side: %.3g m$^{-2}$s$^{-1}$" % phi,
                 fontsize=12, pad=55)

    ax.set_xlim(-0.15 * x, 1.02 * x)
    ax.set_ylim(0, 1.0)
    ax.set_yticks([])
    ax.set_xlabel("Depth into the slab, pc", fontsize=12)

    #simple legend
    ax.add_patch(Rectangle((0, 0), 0, 0, facecolor=ionColour, label="ionised"))
    ax.add_patch(Rectangle((0, 0), 0, 0, facecolor=neutralColour, label="neutral"))
    ax.legend(loc="lower left", fontsize=9, frameon=False)
    fig.canvas.draw_idle()


#call update whenever a slider is moved
sflux.on_changed(update)
snlay.on_changed(update)
for s in sdens + swidth:
    s.on_changed(update)

update(None)
plt.show()
