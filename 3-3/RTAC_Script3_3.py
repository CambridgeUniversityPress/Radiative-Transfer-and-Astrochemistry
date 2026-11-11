#
# Script 3-3 for the textbook "Foundations of Radiative Transfer & Astrochemistry"
# A Monte Carlo scattered light code for a uniform spherical scattering medium that is
# illuminated from behind by a planar background radiation field
# (A solution for Quick Exercise 3.5).
#

#import libraries
import numpy as np
import matplotlib.pyplot as plt


#set up a photon packet at a random position on the back (z<0) of the sphere,
#travelling in the +z direction towards the observer. The positions are chosen to be
#uniform over the projected disc of the sphere, as for a planar radiation field
def initialisePhoton(rExtent):
    rDisc = rExtent * np.sqrt(np.random.uniform(0.0, 1.0))
    phi   = np.random.uniform(0.0, 2.0 * np.pi)

    x = rDisc * np.cos(phi)
    y = rDisc * np.sin(phi)
    z = -np.sqrt(rExtent**2 - rDisc**2)

    return x, y, z, 0.0, 0.0, 1.0

#A scattering event has occurred, choose a new random unit propagation vector for a given g.
def randomScatteringDirection(xhat, yhat, zhat, g):

    if abs(g) < 1e-6:
        # Isotropic: cos(theta) is uniform on [-1, 1]
        cos_theta = np.random.uniform(-1.0, 1.0)
    else:
        # Henyey-Greenstein scattering phase function
        xi = np.random.uniform(0.0, 1.0)
        cos_theta = (1.0 / (2.0 * g)) * (
            1.0 + g**2 - ((1.0 - g**2) / (1.0 - g + 2.0 * g * xi))**2
        )
        cos_theta = np.clip(cos_theta, -1.0, 1.0)

    sin_theta = np.sqrt(1.0 - cos_theta**2)

    # Azimuthal angle around the current direction is always uniform
    phi = np.random.uniform(0.0, 2.0 * np.pi)

    if abs(xhat) <= abs(yhat) and abs(xhat) <= abs(zhat):
        reference = np.array([1.0, 0.0, 0.0])
    elif abs(yhat) <= abs(xhat) and abs(yhat) <= abs(zhat):
        reference = np.array([0.0, 1.0, 0.0])
    else:
        reference = np.array([0.0, 0.0, 1.0])

    current_dir = np.array([xhat, yhat, zhat])

    # First perpendicular vector (normalised)
    perp1 = np.cross(current_dir, reference)
    perp1 = perp1 / np.linalg.norm(perp1)

    # Second perpendicular vector (normalised), orthogonal to both
    perp2 = np.cross(current_dir, perp1)
    perp2 = perp2 / np.linalg.norm(perp2)

    # The new photon direction
    new_dir = (sin_theta * np.cos(phi) * perp1
             + sin_theta * np.sin(phi) * perp2
             + cos_theta             * current_dir)

    # Normalise to correct any floating point drift
    new_dir = new_dir / np.linalg.norm(new_dir)

    return new_dir[0], new_dir[1], new_dir[2]

#travel some random path length along the current unit vector to the next scattering
def randomPathLength(x, y, z, xhat, yhat, zhat, alphaScat):
    xi = np.random.uniform(0.0, 1.0)
    xi = max(xi, 1e-10)

    # Path length to next scattering event
    s = -np.log(xi) / alphaScat

    # Advance the photon along its direction
    x = x + s * xhat
    y = y + s * yhat
    z = z + s * zhat

    return x, y, z


#extent of the scattering medium
rExtent=1.0

#number of photon packets to send
#more are needed than in script 3-1 because most photons do not reach the observer
#(the peel-off technique would be useful for this example if you wanted another challenge!)
nPhotons=50000

#scattering parameter
g=0.8 #g=0 is isotropic scattering

#optical depth from the centre to the edge of the medium
tau=10.0

#the observer is in the +z direction. Only photons escaping within this cone around
#the +z axis are counted in the image (cosThetaMin=1 would be a perfect detector,
#but almost no scattered photons would ever reach it)
cosThetaMin=0.9

#scattering coefficient
#note that this is not in any kind of physical units
#recall tau = alpha * L, so alpha is set from the optical depth and the medium radius
alphaScat=tau/rExtent


nBins = 50

intensityProfile = np.zeros(nBins)

binEdges=np.linspace(0.0, rExtent, nBins + 1)
binCentres=0.5 * (binEdges[:-1] + binEdges[1:])


#the main execution loop for each of the photons in the Monte Carlo calculation
for i in range(0, nPhotons):
    x, y, z, xhat, yhat, zhat = initialisePhoton(rExtent)
    inMedium = True
    #follow the photon packet until it leaves the medium
    while (inMedium):
        x, y, z = randomPathLength(x, y, z, xhat, yhat, zhat, alphaScat)
        bigR = (x**2+y**2+z**2)**0.5
        if(bigR > rExtent):
            inMedium=False
            #only count photons heading towards the observer. There is no longer
            #spherical symmetry, so photons escaping in other directions are not used
            if zhat >= cosThetaMin:
                #impact parameter of the escaping ray, i.e. where it lands on the image
                b = np.linalg.norm(np.cross([x, y, z], [xhat, yhat, zhat]))
                bin_index = int(b / rExtent * nBins)
                if bin_index < nBins:
                    intensityProfile[bin_index] += 1
        else:
            #scattering event
            xhat, yhat, zhat = randomScatteringDirection(xhat, yhat, zhat, g)



annulusArea = np.pi * (binEdges[1:]**2 - binEdges[:-1]**2)
annulusArea = np.where(annulusArea > 0, annulusArea, 1.0)
intensityProfile_norm = intensityProfile / (nPhotons * annulusArea)
fig, axes = plt.subplots(1, 2, figsize=(13, 5))
ax = axes[0]

#First plot the 1D intensity profile
ax.plot(binCentres, intensityProfile_norm, 'o-', color='blue', markersize=4)
ax.set_xlabel('Impact parameter', fontsize=12)
ax.set_ylabel('Intensity', fontsize=12)
ax.set_title(f'1D intensity profile  (g={g}, tau={tau})', fontsize=13)
ax.grid(True, alpha=0.3)

# Use the spherical symmetry to paint the 1D profile onto a 2D grid.
ax2 = axes[1]
nGrid = 200
xs = np.linspace(-rExtent, rExtent, nGrid)
ys = np.linspace(-rExtent, rExtent, nGrid)
XX, YY = np.meshgrid(xs, ys)
BB = np.sqrt(XX**2 + YY**2)   # impact parameter at each pixel

# Map each pixel's b onto a bin index and look up the intensity
image = np.zeros((nGrid, nGrid))
binIndices = np.floor(BB / rExtent * nBins).astype(int)
binIndices = np.clip(binIndices, 0, nBins - 1)
image = intensityProfile_norm[binIndices]

# Mask pixels outside the sphere
image[BB > rExtent] = np.nan

im = ax2.imshow(image, origin='lower',
                extent=[-rExtent, rExtent, -rExtent, rExtent],
                cmap='inferno')
plt.colorbar(im, ax=ax2, label='Intensity')
ax2.set_xlabel('x', fontsize=12)
ax2.set_ylabel('y', fontsize=12)
ax2.set_title(f'2D projected image', fontsize=13)

plt.tight_layout()
plt.savefig('MCScat_backlit_output.pdf')
plt.show()


