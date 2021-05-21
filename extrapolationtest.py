import math


def extrapolateTrack(pt, eta, phi, charge, x0, y0, z0):
    #R = pt / (0.3 * 3.8)
    R = pt / (0.3 * 3.8 )
    xC = x0/100.0 + R* math.cos(phi + charge * 3.14159265359/2.0)
    yC = y0/100.0 + R* math.sin(phi + charge * 3.14159265359/2.0)


    # calculate x,y intersection of track
    a = (-1* yC) / xC
    rb = 129.0 / 100.0
    RC2 = xC**2 + yC**2
    b = (RC2 - R**2 + rb**2) / (2*xC)

    qa = a**2 + 1
    qb = 2*a*b
    qc = b**2 - rb**2
    disc = qb**2 - 4*qa*qc

    y,x,y_other,x_other = 0,0,0,0
    if( disc > 0):
        # barrel can be hit, solution exists
        y1 = (-qb + math.sqrt(disc)) / (2*qa)
        y2 = (-qb - math.sqrt(disc)) / (2*qa)
        x1 = b + y1*a
        x2 = b + y2*a

        if(phi > 0):
            y = y1
            x = x1
            y_other = y2
            x_other = x2
        else:
            y = y2
            x = x2
            y_other = y1
            x_other = x1

        #Z_ECAL = z0/100.0 + R * math.sinh(eta) *(math.acos(1 - 1.29*1.29 / (2*R*R)));
        Z_ECAL = z0/100.0 - (math.asin((y-yC)/R) - phi - charge*3.14159265359/2.0) * (R*math.sinh(eta) / charge)

        if(Z_ECAL > 2.935):
            Z_ECAL = 2.935
        if(Z_ECAL < -2.935):
            Z_ECAL = -2.935
    else:
        # Barrel cannot be hit, endcap is hit (electron spirals out)
        if(eta > 0):
            Z_ECAL = 2.935
        else:
            Z_ECAL = -2.935
# -charge valtoztatva +ra
    X_ECAL = xC + R * math.cos(charge * (Z_ECAL-z0/100.0)/(R * math.sinh(eta)) + phi + charge * 3.14159265359/2.0)
    Y_ECAL = yC + R * math.sin(charge * (Z_ECAL-z0/100.0)/(R * math.sinh(eta)) + phi + charge * 3.14159265359/2.0)
  
    D_ECAL = math.sqrt(X_ECAL*X_ECAL+Y_ECAL*Y_ECAL)

    etaSC = math.asinh(Z_ECAL/D_ECAL)
    if(Y_ECAL > 0):
        phiSC = math.cos(X_ECAL/D_ECAL)
    else:
        phiSC = -1.0 * math.cos(X_ECAL/D_ECAL)


    return etaSC, phiSC, x, y, x_other, y_other, xC, yC

def phiToXY(phi):
    rb = 129.0 / 100.0
    x = rb * math.cos(phi)
    y = rb * math.sin(phi)
    return x,y

def XYToPhi(x,y):
    return math.atan2(y,x)



eta = 0.000000001
phi = 3.1415 / 2
pt = 50
x0,y0,z0=0,0,0
charge = +1


etaSC, phiSC, x, y, x_other, y_other, xC, yC = extrapolateTrack(pt, eta, phi, charge, x0, y0, z0)
print phi, phiSC, XYToPhi(x,y), xC, yC , charge
print phi, phiSC, XYToPhi(x_other,y_other), xC, yC , charge

print ""
charge = -1
etaSC, phiSC, x, y, x_other, y_other, xC, yC  = extrapolateTrack(pt, eta, phi, charge, x0, y0, z0)
print phi, phiSC, XYToPhi(x,y), xC, yC , charge
print phi, phiSC, XYToPhi(x_other,y_other), xC, yC , charge