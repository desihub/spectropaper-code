#!/usr/bin/env python

from math import *
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import os
import sys

########################################################################
primary_mirror_diameter = 3797. # mm
half_angle_fov = 1.6 # deg
distance_from_primary_mirror_to_screen = 16382. #mm
distance_from_upper_ring_to_screen = 6332. #mm
distance_from_upper_ring_to_cage_top = 2259. # mm
distance_from_upper_ring_to_cage_bottom = 2000. # mm MISSING QUANTITY HERE ??? need length of cage, not found in TDR
occultation_diameter = 1800. # from desimodel/data/desi.yaml and TDR Table 2.5 : cage outer diameter
radius_of_cage = occultation_diameter/2. # mm
permaflect_reflectance = 0.935 #at 5000A
upper_ring_diameter = 5400. # mm, TDR p174 The upper ring will be a monolithic weldment 5.4 m in diameter.
radius_of_lamp_to_axis = upper_ring_diameter/2. # mm

# TDR Sec. 2.1.4 Focal Surface and Fig. 2.5
# plate scale = 68. - 75. um/arcsec
# mean plate scale = 254.8 mm/deg (Table 2.2) = 70.77 um/arcsec
# fiber diameter = 107 um
# fiber diameter = 1.57 to 1.42 arcsec , average = 1.51 arcsec
# warning : change of fiber size with angular distance !!
fiber_diameter_um = 107. # um
fiber_angular_diameter_arcsec = 1.51 # arcsec
############
d2r=pi/180. # degree to radian
foot2mm = 304.8 # foot to mm
########################################################################
test_for_sky=False # to get the sky flatfield (non uniform because of platescale and vignetting)
debug_for_sky=False # to check sky is flat if no plate scale nor vignetting

def plot_new_screen() :
    """
Replace 8 rectangles with trapezoid panels (another Dick Joyce Idea)
Reuses the 16 - 2 x 4 foot sheets
Adds 8 - 24 x 48 x 31 inch trapezoid sheets
Adds 4 - 4 x 4 foot right triangle sheets
It has a useful diameter of 5173 mm
FOV up to 2.4deg
Or a lateral placement error of +- 230 mm
This screen will fit on the current frame
    """

    print("TO DO")
    #sys.exit(12)

def plot_screen() :
    panel_length=4*foot2mm # mm
    color="gray"
    for y in panel_length/2.*np.arange(-4,4+1) :
        plt.plot([-panel_length,+panel_length],[y,y],"-",color=color)
    for y in panel_length*np.arange(-1,1+1) :
        plt.plot([-2*panel_length,+2*panel_length],[y,y],"-",color=color)
    for x in panel_length*np.arange(-1,2) :
        plt.plot([x,x],[-2*panel_length,2*panel_length],"-",color=color)
    for x in panel_length*np.array([-2,-1.5,1.5,2]) :
        plt.plot([x,x],[-panel_length,panel_length],"-",color=color)
    for i in [-1,1] :
        for j in [-1,1] :
            plt.plot([i*2*panel_length,i*panel_length],[j*panel_length,j*2*panel_length],"-",color=color)
    plt.grid()


radius_of_sensitive_area_on_screen = primary_mirror_diameter/2.+distance_from_primary_mirror_to_screen*tan(half_angle_fov*d2r)

print("radius_of_sensitive_area_on_screen   = %f mm"%radius_of_sensitive_area_on_screen)
print("diameter_of_sensitive_area_on_screen = %f mm"%(2*radius_of_sensitive_area_on_screen))



diameter_of_screen = 4*4*foot2mm # mm
half_angle_screen_fov = atan((diameter_of_screen/2.-primary_mirror_diameter/2.)/distance_from_primary_mirror_to_screen)/d2r # deg
print("diameter_of_screen = %f mm -> fov = %f deg -> displacement = %f mm"%(diameter_of_screen,half_angle_screen_fov,(diameter_of_screen/2-radius_of_sensitive_area_on_screen)))

# intensity of input isotropic beam on screen per unit area of screen is
# I = Itot/(4*pi*d**2)*cos(theta) where theta is indicence angle to normal
# this is the definition of the lambertian reflectance
#
# d = d0/cos(theta)
# I = Itot/(4*pi*d0**2)cos3(theta)
#
if test_for_sky :
    print("TEST_FOR_SKY")
    diameter_of_screen *= 1.1
    distance_from_primary_mirror_to_screen += 50.

# this is a coordinate grid on the screen in mm
nx=1000 # accuracy=2e-4 for nx=10000, <1e-4 for nx=2000
x=np.tile(np.linspace(-diameter_of_screen/2,diameter_of_screen/2,nx),(nx,1))
y=x.copy().T
# the area of one pixel on the screen
darea=(x[0,1]-x[0,0])**2 #mm

# intensities in the screen pixels
intensity=np.zeros(x.shape)
intensity_wo_vignetting=np.zeros(x.shape)

# coordinates on lamps in mm, on the plane defined by the upper ring
lamps_coords=[[radius_of_lamp_to_axis,0],[-radius_of_lamp_to_axis,0],[0,radius_of_lamp_to_axis,0],[0,-radius_of_lamp_to_axis]]
lamps=[]

cmap=plt.cm.RdYlBu


if not test_for_sky :
    # computation of the intensity of light hitting the screen pixels
    number_of_lamps=0
    lamps=[0,1,2,3]
    for i in lamps :

        number_of_lamps += 1
        lamp=lamps_coords[i]

        # indicence angle to normal of screen
        tan2_theta = ((x-lamp[0])**2+(y-lamp[1])**2)/distance_from_upper_ring_to_screen**2
        cos_theta  = 1/np.sqrt(1.+tan2_theta)
        print("lamp #%d x=%f,y=%f"%(i,lamp[0],lamp[1]))

        # coord of rays on plane of cage top
        x_cam_top = (x-lamp[0])*distance_from_upper_ring_to_cage_top/distance_from_upper_ring_to_screen+lamp[0]
        y_cam_top = (y-lamp[1])*distance_from_upper_ring_to_cage_top/distance_from_upper_ring_to_screen+lamp[1]
        occultation = (x_cam_top**2+y_cam_top**2)>radius_of_cage**2

        lamp_relative_brightness = 1. # to play

        intensity += lamp_relative_brightness*occultation/(4.*pi*distance_from_upper_ring_to_screen**2)*cos_theta**3*darea
        intensity_wo_vignetting += lamp_relative_brightness/(4.*pi*distance_from_upper_ring_to_screen**2)*cos_theta**3*darea




    # this is sub-set of the pixels that are actually on the screen
    screenmask = (np.abs(x)+np.abs(y))<diameter_of_screen*3./4.
    relative_intensity=intensity/np.max(intensity)
    min_val=np.min(relative_intensity[screenmask][relative_intensity[screenmask]>0.75])
    plt.figure("light on screen")
    plt.imshow(screenmask*relative_intensity,origin="lower",interpolation="nearest",extent=(-diameter_of_screen/2,diameter_of_screen/2,-diameter_of_screen/2,diameter_of_screen/2),vmin=min_val,cmap=cmap)
    plot_screen()
    plt.xlabel("x on screen (mm)")
    plt.ylabel("y on screen (mm)")
    plt.colorbar()

    min_val=np.min(intensity[screenmask])
    max_val=np.max(intensity[screenmask])
    print("min rel. intensity on screen (with vignetting) = %f"%(min_val/max_val))
    min_val=np.min(intensity_wo_vignetting[screenmask])
    max_val=np.max(intensity_wo_vignetting[screenmask])
    print("min rel. intensity on screen (without vignetting) = %f"%(min_val/max_val))


#plt.show()
#sys.exit(12)

if test_for_sky :
    print("WARNING, THIS IS A TEST TO EMULATE SKY")
    intensity *= 0.
    intensity += 1.



# now deal with fibers
arcsec2rad=pi/(180.*3600)
fiber_solid_angle_rad2_approx = pi*(fiber_angular_diameter_arcsec/2.*arcsec2rad)**2 # rad**2
#print "fiber_solid_angle_rad2 =",fiber_solid_angle_rad2

# fiber position on focal plane in angles

tmp=np.loadtxt("%s/data/focalplane/fiberpos.txt"%os.environ["DESIMODEL"]).T
fiber_id=tmp[0].astype(int)
fiber_positionner_id=tmp[1].astype(int)
fiber_spectro_id=tmp[2].astype(int)
fiber_x=tmp[3] # mm
fiber_y=tmp[4] # mm
tmp=np.loadtxt("%s/data/focalplane/platescale.txt"%os.environ["DESIMODEL"]).T
fp_r=tmp[0] # mm
fp_theta=tmp[1] # deg
fp_mps=tmp[6] # um/arcsec
fp_sps=tmp[7] # um/arcsec
# recompute with refined plate scale

fiber_r     = np.sqrt(fiber_x**2+fiber_y**2)
fiber_theta = np.interp(fiber_r,fp_r,fp_theta)*np.pi/180. # rad
fiber_phi   = np.arctan2(fiber_y,fiber_x)
fiber_mps   = np.interp(fiber_r,fp_r,fp_mps)
fiber_sps   = np.interp(fiber_r,fp_r,fp_sps)


fiber_solid_angle_rad2 = np.pi*(fiber_diameter_um/2.)**2/(fiber_mps*fiber_sps)*arcsec2rad**2 # rad**2

if debug_for_sky : # constant solid angle to check the flat is flat for sky
    fiber_solid_angle_rad2 *= 0.
    fiber_solid_angle_rad2 += fiber_solid_angle_rad2_approx


if False :
    plt.figure()
    plt.plot(fp_theta,fiber_solid_angle_rad2,"o")
    plt.plot(fp_theta,fiber_solid_angle_rad2_approx*np.ones(fp_theta.shape),"-")
    print("fiber_solid_angle_rad2_approx=",fiber_solid_angle_rad2_approx)
    print("mean(fiber_solid_angle_rad2)=",np.mean(fiber_solid_angle_rad2))
    plt.show()




if False :
    plt.figure()
    spectros=np.unique(fiber_spectro_id)
    for spectro in spectros :
        ok=fiber_spectro_id==spectro
        plt.plot(fiber_x[ok],fiber_y[ok],"o")

if False :
    plt.figure()
    ax=fiber_theta*np.cos(fiber_phi)
    ay=fiber_theta*np.sin(fiber_phi)
    spectros=np.unique(fiber_spectro_id)
    for spectro in spectros :
        ok=fiber_spectro_id==spectro
        plt.plot(ax[ok]*180./np.pi,ay[ok]*180./np.pi,"o")
#plt.show()



fibers=np.arange(fiber_id.size).astype(int)
fibers=[fiber_id.size-1]
#nf=200 ; fibers=(np.arange(nf)*(fiber_id.size/nf)).astype(int)

fraction_of_lamp_photons_in_fibers = np.zeros(fiber_id.size)

print(fiber_theta.shape)
print(fiber_solid_angle_rad2.shape)


for fiber in fibers :



    # intensity of the screen pixels with light collected by this fiber
    intensity_in_fiber = intensity.copy()

    # compute fraction of light on screen reflected at the fiber solid angle
    # assume an ideal diffuse reflectance
    #
    # if phi is angle to normal of reflected light
    # d2I/darea/dsolidangle = dI/darea * alpha * cos(phi)
    # this means that the total intensity seen from an inclined surface per unit solid angle ("luminance") is independent of phi
    # because the area seen in the solid angle is darea propto 1/cos(phi)
    # dI/dsolidangle = dI/darea * alpha * cos(phi) *darea = cst
    #
    #
    # CORRECT :
    # dI/darea = int_0 ^pi/2 domega(phi) * cos(phi) * alpha
    #          = int_0 ^pi/2 2 pi * sin(phi) * cos(phi) * alpha
    #          = 2*pi*alpha [ 1/2 sin^2(phi) ]_0^(pi/2)
    #          = 2*pi*alpha * 1/2
    #          = pi*alpha
    # so alpha = 1/pi
    #
    # AND :
    # dI/darea = int_0 ^pi/2 domega(phi) * cos(phi) / pi
    #
    #
    if not test_for_sky :
        intensity_in_fiber *= np.cos(fiber_theta[fiber])/np.pi*fiber_solid_angle_rad2[fiber]
    else :
        # no cos effect for sky gives uniform flat if no plate scale variation and no occultation
        # accuracy < 1e-4 for nx=2000
        intensity_in_fiber *= fiber_solid_angle_rad2[fiber]

    # now compute occultation

    tan_theta = np.tan(fiber_theta[fiber])
    cos_phi   = np.cos(fiber_phi[fiber])
    sin_phi   = np.sin(fiber_phi[fiber])

    if not debug_for_sky :

        distance_from_screen_to_cage_top = distance_from_upper_ring_to_screen-distance_from_upper_ring_to_cage_top
        distance_from_screen_to_cage_bottom = distance_from_upper_ring_to_screen+distance_from_upper_ring_to_cage_bottom

        # occultation on top of cage
        x_cam_top = x+tan_theta*cos_phi*distance_from_screen_to_cage_top
        y_cam_top = y+tan_theta*sin_phi*distance_from_screen_to_cage_top
        intensity_in_fiber *= ((x_cam_top**2+y_cam_top**2)>radius_of_cage**2)

        # occultation on bottom of cage
        x_cam_bottom = x+tan_theta*cos_phi*distance_from_screen_to_cage_bottom
        y_cam_bottom = y+tan_theta*sin_phi*distance_from_screen_to_cage_bottom

        intensity_in_fiber *= ((x_cam_bottom**2+y_cam_bottom**2)>radius_of_cage**2)

    # hitting the mirror
    x_mirror = x+tan_theta*cos_phi*distance_from_primary_mirror_to_screen
    y_mirror = y+tan_theta*sin_phi*distance_from_primary_mirror_to_screen
    r_in_mirror = np.sqrt(x_mirror**2+y_mirror**2)
    intensity_in_fiber *= (r_in_mirror<primary_mirror_diameter/2.)*(r_in_mirror>occultation_diameter/2.)

    # No need for fiber angular acceptance
    # "TDR p63 : the ability of the fibers limiting numerical aperture (NA) to fully accept the input f/# from the corrector."

    fraction_of_lamp_photons_in_fibers[fiber] = np.sum(intensity_in_fiber)
    print("fiber=%d thetha=%f deg phi=%f deg, frac of light in fiber = %g"%(fiber,fiber_theta[fiber]*180/np.pi,fiber_phi[fiber]*180/np.pi,fraction_of_lamp_photons_in_fibers[fiber]))

    if True :
        plt.figure("light on screen seen by fiber at %2.1fdeg"%(fiber_theta[fiber]*180/np.pi))
        norme=1./np.max(intensity_in_fiber)
        plt.imshow(norme*intensity_in_fiber,origin="lower",interpolation="nearest",extent=(-diameter_of_screen/2,diameter_of_screen/2,-diameter_of_screen/2,diameter_of_screen/2),vmin=norme*np.min(intensity_in_fiber[intensity_in_fiber>0])*0.98,cmap=cmap)

        # fiber numerical aperture NA = 0.22 = sin(theta) -> theta = 12.7 deg
        #fiber_acceptance_angle=7. # deg
        #radius=tan(fiber_acceptance_angle*d2r)


        plot_screen()
        plt.xlabel("x on screen (mm)")
        plt.ylabel("y on screen (mm)")
        plt.colorbar()

        #plt.show()





# obtained with test run, due to occultation
#sky_non_uniformity=np.interp(fiber_theta*180./np.pi,[0.,1.6],[1.,0.929])



#relative_to_sky = True

#if relative_to_sky :
    # correction with sky
#    fraction_of_lamp_photons_in_fibers /= sky_non_uniformity

if False and ( len(lamps)==1 or test_for_sky ) and ( not debug_for_sky ) :
    if test_for_sky :
        ofilename="fiber_illumination_fraction_sky.txt"
    else :
        ofilename="fiber_illumination_fraction_lamp_%d.txt"%(lamps[0])
    print("write file %s"%ofilename)
    file=open(ofilename,"w")
    file.write("# fiber x(mm) y(mm) theta(deg) phi(deg) illumination_fraction\n")
    for fiber in fibers :
        file.write("%d %f %f %f %f %g\n"%(fiber,fiber_x[fiber],fiber_y[fiber],fiber_theta[fiber]*180/np.pi,fiber_phi[fiber]*180/np.pi,fraction_of_lamp_photons_in_fibers[fiber]))
    file.close()

norme=1./np.max(fraction_of_lamp_photons_in_fibers)

#plt.plot(fiber_fov_angles,norme*fraction_of_lamp_photons_in_fibers,"o-")

#if not relative_to_sky :
#    plt.plot(fiber_fov_angles,sky_non_uniformity,"--",color="gray",label="SKY")
#plt.ylim(np.min(fraction_of_lamp_photons_in_fibers)*norme-0.01,1.01)
#plt.xlabel("Radial angle (deg)")
#if relative_to_sky :
#    plt.ylabel("Fiber to fiber uniformity (screen/sky)")
#else :
#    plt.ylabel("Fiber to fiber uniformity")
#plt.grid()


plt.figure("Fiber to fiber uniformity")
min_frac=np.min(fraction_of_lamp_photons_in_fibers[fibers])
max_frac=np.max(fraction_of_lamp_photons_in_fibers[fibers])

print("range of variation=",min_frac/max_frac-1)

cmap=plt.cm.RdYlBu
plt.subplot(1,2,1)
plt.scatter(fiber_x[fibers],fiber_y[fibers],c=(fraction_of_lamp_photons_in_fibers[fibers]-min_frac)/(max_frac-min_frac),cmap=cmap,edgecolors="face",s=100)
a=plt.subplot(1,2,2)
norm=mpl.colors.Normalize(vmin=min_frac/max_frac,vmax=1.)
cb=mpl.colorbar.ColorbarBase(a, cmap=cmap,norm=norm,orientation='vertical')


if not test_for_sky :
    mean_fraction_of_lamp_photons_in_fibers = np.mean(fraction_of_lamp_photons_in_fibers)
    print("max non-uniformity = ",np.min(fraction_of_lamp_photons_in_fibers)/np.max(fraction_of_lamp_photons_in_fibers))
    print("mean_fraction_of_lamp_photons_in_fibers (per lamp) =",mean_fraction_of_lamp_photons_in_fibers/number_of_lamps)
    print("(before instrument throughput and screen reflectance)")


plt.show()
