import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm, colors
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import copy

global radPerDeg, fov, beta
radPerDeg = 2.*np.pi/360.
fov = 24*radPerDeg
beta = fov/2.

class line(object):
    ''' Line segment between two input points.
        Attributes: x,y,z give appropriate points'''
    n = 101 # Number of resolved points per line
    def __init__(self, pt1, pt2):
        self.pt1, self.pt2 = pt1, pt2
        t = np.linspace(0, 1, line.n)
        linePts = np.zeros([len(t),3])
        for i in range(len(linePts)):
            linePts[i] = self.pt1 + t[i]*(self.pt2-self.pt1)
        self.x, self.y, self.z = linePts[:,0], linePts[:,1], linePts[:,2]

class latLines(object):
    phi = np.linspace(0, 2*np.pi, 100)
    def __init__(self, theta_const):
        self.theta_const = theta_const
        self.x = np.sin(theta_const)*np.cos(latLines.phi)
        self.y = np.sin(theta_const)*np.sin(latLines.phi)
        self.z = np.cos(theta_const)
        
class longLines(object):
    theta = np.linspace(0, np.pi, 100)
    def __init__(self, phi_const):
        self.phi_const = phi_const
        self.x = np.sin(longLines.theta)*np.cos(phi_const)
        self.y = np.sin(longLines.theta)*np.sin(phi_const)
        self.z = np.cos(longLines.theta)        

class plane(object):
    '''Given three points in a plane (in quadrant 1), this object creates the plane, and
    then masks points not bound by the unique in-plane triangle defined by the points. '''
    n = line.n
    def __init__(self, pt1, pt2, pt3):
        self.pt1, self.pt2, self.pt3 = pt1, pt2, pt3
        self.normal = np.cross(self.pt2-self.pt1, self.pt3-self.pt1)  # normal vector to plane        
        
        # Plane eqn: a*x+b*y+c*z+d=0
        # [a,b,c] is the normal. Thus, we have to calculate d and we're set
        d = np.dot(self.normal, self.pt1) # any point is fine

        # calculate corresponding z on xG, yG
        self.z = (d - self.normal[0]*xG - self.normal[1]*yG) / self.normal[2]

def longLatLines(polarProj, ax):
    '''Draw longitude and latitude lines, given boolean value of 
    whether or not you visualize as polar projection'''
    phi_c = np.linspace(0., 360., 13)*radPerDeg # long lines, 30 deg sep
    if polarProj == False:
        theta_c = np.linspace(0., 180, 7)*radPerDeg # lat lines, 30 deg sep both hemi
    if polarProj == True:
        theta_c = np.linspace(0., 90., 4)*radPerDeg
    for i in range(len(theta_c)):
        latLine = latLines(theta_c[i])
        colorStr, alph, ls = 'black', 0.2, '-'
        if latLine.theta_const/radPerDeg == 90:
            colorStr, alph, ls = 'black', 0.8, '--'
        ax.plot(latLine.x, latLine.y, latLine.z, c=colorStr, alpha=alph, linestyle=ls)
    for i in range(len(phi_c)):
        longLine = longLines(phi_c[i])
        ax.plot(longLine.x, longLine.y, longLine.z, c='black', alpha=0.3)    

def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians (Euler-Rodrigues formula).
    """
    axis = np.asarray(axis)
    theta = np.asarray(theta)
    axis = axis/np.sqrt(np.dot(axis, axis))
    a = np.cos(theta/2)
    b, c, d = -axis*np.sin(theta/2)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])

def phiRotationEarly(startingVecs, phiRot):
    """
    Return input (x,y,z) points (passed as numpy array of shape 3*n*n) 
    by phiRot radians)
    """
    if np.shape(startingVecs)[0] != 3:
        print "Check input vector array to phiRotation."
    zAxis = [0,0,1]
    rotatedVecs = np.zeros(np.shape(startingVecs))
    for i in range(np.shape(startingVecs)[1]):
        for j in range(np.shape(startingVecs)[2]):
                rotatedVecs[:,i,j] = np.dot(rotation_matrix(zAxis,phiRot), startingVecs[:,i,j])
    return rotatedVecs

## EVERYTHING BELOW HERE IS POST-FIRST MAIN CHECK

class patch(object):
    numQuad = 4 # number of quadrants
    """Take [4, 3, 101, 101] (or modified nx, ny) singleCamField, make one thing"""
    def __init__(self, singleCam):
        self.singleCam = singleCam
        self.nx = np.shape(singleCam)[2]
        self.ny = np.shape(singleCam)[3]
        self.x = np.zeros([patch.numQuad*self.nx, self.ny])
        self.y = np.zeros([self.nx, patch.numQuad*self.ny])
        #hacky, but it's only once
        self.x = np.concatenate((self.singleCam[0,0,:,:],self.singleCam[1,0,:,:],self.singleCam[2,0,:,:],self.singleCam[3,0,:,:]), axis=0)
        self.y = np.concatenate((self.singleCam[0,1,:,:],self.singleCam[1,1,:,:],self.singleCam[2,1,:,:],self.singleCam[3,1,:,:]), axis=0)
        self.z = np.concatenate((self.singleCam[0,2,:,:],self.singleCam[1,2,:,:],self.singleCam[2,2,:,:],self.singleCam[3,2,:,:]), axis=0)
        
class fourCamView(object):
    """Take first camera full patch, make four cameras single object"""
    camAngleSep = fov # presumes fixed angle between cameras is 24.2deg
    def __init__(self, cam1):
        self.cam1 = cam1
        axis = [1,0,0] # will rotate about x-axis
        
        #theta = np.array([0, fourCamView.camAngleSep, 2.*fourCamView.camAngleSep, 3.*fourCamView.camAngleSep]) 
        #theta = np.array([-1.5*fourCamView.camAngleSep, -0.5*fourCamView.camAngleSep, 0.5*fourCamView.camAngleSep, 1.5*fourCamView.camAngleSep]) 
        #theta = np.array([0+180.+6., fourCamView.camAngleSep+180.+6., 2.*fourCamView.camAngleSep+180.+6., 3.*fourCamView.camAngleSep+180.+6.])
        theta = np.array([0+12.*radPerDeg, fourCamView.camAngleSep+12.*radPerDeg, 2.*fourCamView.camAngleSep+12.*radPerDeg, \
                          3.*fourCamView.camAngleSep+12.*radPerDeg]) 

        
        self.cam1pos = np.array([self.cam1.x, self.cam1.y, self.cam1.z]) # shape: (3, 404, 101)
        #again hacky, but tryna make things happen
        self.cam2pos = np.zeros([3, patch.numQuad*self.cam1.nx, self.cam1.nx])
        self.cam3pos = np.zeros([3, patch.numQuad*self.cam1.nx, self.cam1.nx])
        self.cam4pos = np.zeros([3, patch.numQuad*self.cam1.nx, self.cam1.nx])
        self.cam5pos = np.zeros([3, patch.numQuad*self.cam1.nx, self.cam1.nx]) # hax
        self.allCams = np.zeros([3, 4*patch.numQuad*self.cam1.nx, self.cam1.nx]) # concatenate along long axis, always
        for j in range(np.shape(self.cam1pos)[1]):
            for k in range(np.shape(self.cam1pos)[2]):
                self.cam5pos[:,j,k] = np.dot(rotation_matrix(axis,theta[0]), self.cam1pos[:,j,k]) # RECENT ADDON
                self.cam2pos[:,j,k] = np.dot(rotation_matrix(axis,theta[1]), self.cam1pos[:,j,k])
                self.cam3pos[:,j,k] = np.dot(rotation_matrix(axis,theta[2]), self.cam1pos[:,j,k])
                self.cam4pos[:,j,k] = np.dot(rotation_matrix(axis,theta[3]), self.cam1pos[:,j,k])
        
        self.allCams = np.concatenate((self.cam5pos, self.cam2pos, self.cam3pos, self.cam4pos),axis=1)
        self.x = self.allCams[0,:,:]
        self.y = self.allCams[1,:,:]
        self.z = self.allCams[2,:,:]

def phiRotation(view1, numPhiRot):
    """
    Relevant for yamaka plots.
    Return viewing points, rotated in phi numPhiRot times (presumes
    single year, one hemisphere)
    """
    phiRot = 360./13.*radPerDeg*numPhiRot # 27.69 deg rotations (13/hemi; watch for int division!)
    zAxis = [0,0,1]
    view = copy.deepcopy(view1)
    for j in range(np.shape(view1.allCams)[1]):
        for k in range(np.shape(view1.allCams)[2]):
            view.allCams[:,j,k] = np.dot(rotation_matrix(zAxis,phiRot), view1.allCams[:,j,k])
    view.x = view.allCams[0,:,:]
    view.y = view.allCams[1,:,:]
    view.z = view.allCams[2,:,:]
    return view

def fieldRotation(view1, axis, angle):
    """
    Return viewing points, rotated in phi by 'angle' radians about axis
    """
    view = copy.deepcopy(view1)
    for j in range(np.shape(view1.allCams)[1]):
        for k in range(np.shape(view1.allCams)[2]):
            view.allCams[:,j,k] = np.dot(rotation_matrix(axis,angle), view1.allCams[:,j,k])
    view.x = view.allCams[0,:,:]
    view.y = view.allCams[1,:,:]
    view.z = view.allCams[2,:,:]
    return view

def save_object(obj, filename):
    with open(filename, 'wb') as output:
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)


if __name__ == '__main__':

    xMax, xMin = np.sin(fov/2), 0
    yMax, yMin = np.sin(fov/2), 0
    x, y = np.linspace(xMin, xMax, 101), np.linspace(yMin, yMax, 101)
    xG, yG = np.meshgrid(x, y) # only grid defined is on 1st quadrant, 100x100-> 40k pts per patch at end.
    zSphOnG = np.sqrt(1 - xG**2 - yG**2)

    pts = np.array([[0,0,0], [0,0,np.cos(beta)], [0,np.sin(beta),np.cos(beta)], 
                    [np.sin(beta),np.sin(beta),np.cos(beta)], [np.sin(beta),0,np.cos(beta)]])

    lines = [line(pts[0],pts[1]), line(pts[0],pts[2]), 
             line(pts[0],pts[3]), line(pts[0],pts[4])]

    pl1 = plane(pts[0],pts[3],pts[4])
    pl2 = plane(pts[0],pts[3],pts[2])

    # actual view:
    zViewOnG = np.copy(zSphOnG)
    zViewOnG = np.where( (zViewOnG>pl1.z)&(zViewOnG>pl2.z), zViewOnG, np.nan)

    rot0 = np.array([xG, yG, zViewOnG])
    rot1 = phiRotationEarly(rot0, np.pi/2.)
    rot2 = phiRotationEarly(rot0, np.pi)
    rot3 = phiRotationEarly(rot0, 3.*np.pi/2.)
    singleCamField = np.array([rot0, rot1, rot2, rot3])


    cam1 = patch(singleCamField)
    view_ini = fourCamView(cam1)

    v = []
    v.append(view_ini)
    for i in range(1,13):
        view_i = phiRotation(view_ini, i)
        v.append(view_i)

    # MAKE PLOT
    plt.close('all')
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111, projection='3d')

    kwargs = dict(linewidth=0, rstride=1, cstride=1, color='blue', alpha=0.1)
    '''
    for i in range(len(singleCamField)):
        ax.plot_surface(singleCamField[i,0,:,:], singleCamField[i,1,:,:], \
                        singleCamField[i,2,:,:], **kwargs)
    ''' 
    for i in range(len(v)):
        ax.plot_surface(v[i].x, v[i].y, -v[i].z, **kwargs)
    longLatLines(False, ax) # keeps full hemisph

    ax.view_init(elev=20, azim=-30)
    ax.set(xlim=[-1,1],ylim=[-1,1],zlim=[-1,1],aspect='equal')
    ax.set_axis_off()
    plt.tight_layout()
    #plt.show()
    plt.savefig('shemiAvoid.png',dpi=256)
    plt.savefig('shemiAvoid.pdf')
