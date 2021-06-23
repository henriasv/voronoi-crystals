import pyvoro
import numpy as np
import json
import itertools
import os

_ROOT = os.path.abspath(os.path.dirname(__file__))
def get_jsondata(path):
    return os.path.join(_ROOT, 'data_json', path)


class Polycrystal:
	def __init__(self, N=2, type="regular", substance="copper", box = np.asarray([[0, 200], [0, 200], [0, 200]])):
		self.setBox(box)
		if type is "regular":
			self.regularVoronoiCenters(N, N, N)
			N = 2 * N**3
		else:
			self.randomVoroCenters(N)
		
		self.dispersion = max(self.box[:,1]-self.box[:,0])
		self.voronoi = None
		self.particles = {}
		self.element_data = json.load(open(get_jsondata("masses_and_radii.json")))
		self.substance=substance
		self.computeVoronoi()
		self.loadCrystalFromFile(self.substance+".json")
		for i in range(N):
			self.addParticlesWithReplicatedCrystal(i, 5)
			print("Added crystal %d of %d" % (i, N))

	def setBox(self, box):
		self.box = box
		self.boxlengths = self.box[:,1]-self.box[:,0]

	def randomVoroCenters(self, N):
		x = np.random.random(N)*(self.box[0,1]-self.box[0,0])+self.box[0,0]
		y = np.random.random(N)*(self.box[1,1]-self.box[1,0])+self.box[1,0]
		z = np.random.random(N)*(self.box[2,1]-self.box[2,0])+self.box[2,0]
		self.voronoiCenters = np.asarray([[x[i], y[i], z[i]] for i in range(N)])

	def boxLengths(self):
		return self.boxlengths

	def regularVoronoiCenters(self, Nx, Ny, Nz):
		x = np.linspace(self.box[0,0], self.box[0,1], Nx, endpoint=False)+0.5*self.boxLengths()[0]/Nx
		y = np.linspace(self.box[1,0], self.box[1,1], Ny, endpoint=False)+0.5*self.boxLengths()[1]/Ny
		z = np.linspace(self.box[2,0], self.box[2,1], Nz, endpoint=False)+0.5*self.boxLengths()[2]/Nz
		products = list(itertools.product(x, y, z))
		for i in range(len(products)):
			products.append([products[i][0]+0.5*self.boxLengths()[0]/Nx, products[i][1]+0.5*self.boxLengths()[1]/Ny, products[i][2]+0.5*self.boxLengths()[2]/Nz])
		self.voronoiCenters = np.asarray(products)

	def plotPoints(self):
		x = self.voronoiCenters[:,0]
		y = self.voronoiCenters[:,1]
		z = self.voronoiCenters[:,2]
		self.ax.plot(x, y, z, "o")

	def plotVoronoi(self):
		if self.voronoi is None:
			self.computeVoronoi()
		for cell in self.voronoi:
			#print cell
			vertices = cell['vertices']
			for face in cell['faces']:
				for i in range(len(face['vertices'])):
					vertex1 = vertices[face['vertices'][i]]
					vertex2 = vertices[face['vertices'][i-1]]
					self.ax.plot([vertex1[0], vertex2[0]], [vertex1[1], vertex2[1]], [vertex1[2], vertex2[2]])		

	def periodicSqDistance(self, point1, point2):
		diff = point1-point2
		#return np.sum(diff*diff, axis=1)
		diffplus = np.abs(diff + self.boxlengths)
		diffminus = np.abs(diff - self.boxlengths)
		diff = np.abs(diff)
		diff = np.minimum(diff, diffplus)
		diff = np.minimum(diff, diffminus)
		retval = np.sum(diff*diff, axis=1)
		return retval

	def insideSystemBoxAroundPoint(self, cellnum, points):
		isInside = np.ones(len(points), dtype=bool)
		cellCenter = self.voronoiCenters[cellnum]
		xmin = cellCenter[0]-self.boxlengths[0]/2
		xmax = cellCenter[0]+self.boxlengths[0]/2
		ymin = cellCenter[1]-self.boxlengths[1]/2
		ymax = cellCenter[1]+self.boxlengths[1]/2
		zmin = cellCenter[2]-self.boxlengths[2]/2
		zmax = cellCenter[2]+self.boxlengths[2]/2
		isInside = np.logical_and(isInside, points[:,0]>xmin)
		isInside = np.logical_and(isInside, points[:,0]<xmax)
		isInside = np.logical_and(isInside, points[:,1]>ymin)
		isInside = np.logical_and(isInside, points[:,1]<ymax)
		isInside = np.logical_and(isInside, points[:,2]>zmin)
		isInside = np.logical_and(isInside, points[:,2]<zmax)
		return isInside

	def insideCell(self, cellNum, points):
		p0 = self.voronoiCenters[cellNum]
		p1 = points;
		faces = self.voronoi[cellNum]['faces']
		adjacentCells = [face['adjacent_cell'] for face in faces]
		selfDistance = self.periodicSqDistance(p1, p0)
		isInsideCell = self.insideSystemBoxAroundPoint(cellNum, points)
		for adjacentCell in adjacentCells:
			otherDistance = self.periodicSqDistance(p1, self.voronoiCenters[adjacentCell])
			isInsideCell = np.logical_and(isInsideCell, (np.sqrt(otherDistance)-np.sqrt(selfDistance))>1)
		return isInsideCell


	def randomInBox(self):
		positions = np.random.random(3)
		positions *= self.boxLengths()
		positions += self.box[:,0]
		return positions

	def addRandomParticlesToCell(self, cellNum, N=100):
		for i in range(N):
			position = self.randomInBox()
			if self.insideCell(cellNum, position):
				self.particles.append(position)

	def computeVoronoi(self):
		self.voronoi = pyvoro.compute_voronoi(self.voronoiCenters, self.box, self.dispersion, periodic=[True, True, True])
		return self.voronoi
	
	def loadCrystalFromFile(self, filename):
		import json
		with open(get_jsondata(filename), 'r') as ifile:
			self.crystal = json.load(ifile)

	def rotatePoints(self, points, phi, theta):
		Rx = np.asarray([[1, 0, 0], [0, np.cos(theta), -np.sin(theta)], [0, np.sin(theta), np.cos(theta)]])
		Ry = np.asarray([[np.cos(phi), 0, -np.sin(phi)], [0, 1, 0], [np.sin(phi), 0, np.cos(phi)]])
		#print Ry
		#print points[0]
		#print np.asarray([np.dot(np.dot(Rx,Ry),points[i]).transpose() for i in range(len(points))])
		return np.dot(points, np.dot(Rx, Ry))
		#return points
		#return np.asarray([np.dot(np.dot(Rx,Ry),points[i]) for i in range(len(points))])

	def addParticlesWithReplicatedCrystal(self, cellnum, N):
		x = self.crystal['box']['x'][1]-self.crystal['box']['x'][0]
		y = self.crystal['box']['y'][1]-self.crystal['box']['y'][0]
		z = self.crystal['box']['z'][1]-self.crystal['box']['z'][0]
		xy = self.crystal['box']['xy']
		xz = self.crystal['box']['xz']
		yz = self.crystal['box']['yz']

		preproduct = list(itertools.product(range(-N, N+1), range(-N, N+1), range(-N, N+1)))
		myproduct = []
		for element in preproduct:
			i, j, k = element 
			myproduct.append([i*x+j*xy+k*xz, j*y+k*yz, k*z])
		theta = np.random.random()*np.pi
		phi = np.random.random()*np.pi
		for key, value in self.crystal['particles'].items():
			if not key in self.particles.keys():
				self.particles[key] = []
			positions = value
			tmpPositions = []
			for element in myproduct:
				tmpPositions.extend(np.asarray(positions)+np.asarray(element))
			tmpPositions = self.rotatePoints(np.asarray(tmpPositions), phi, theta)+self.voronoiCenters[cellnum]
			insideCellParticles = self.insideCell(cellnum, tmpPositions)
			moleculePad = np.ones(len(tmpPositions))*cellnum
			moleculePad = moleculePad.reshape((len(moleculePad), 1))
			tmpPositions = np.concatenate((tmpPositions, moleculePad), axis=1)
			self.particles[key].extend(tmpPositions[insideCellParticles])

	def dumpPositionsLammpsData(self, filename):
		with open(filename, "w") as ofile: 
			atomTypes = len(self.particles.keys())
			nAtoms = sum([len(element) for element in self.particles.values()])
			ofile.write("LAMMPS Description\n\n")
			ofile.write("\t\t%d atoms\n" % nAtoms)
			ofile.write("\t\t%d atom types\n\n" % atomTypes)

			ofile.write("\t%f %f xlo xhi\n" % (self.box[0,0], self.box[0,1]))
			ofile.write("\t%f %f ylo yhi\n" % (self.box[1,0], self.box[1,1]))
			ofile.write("\t%f %f zlo zhi\n\n" % (self.box[2,0], self.box[2,1]))

			keys = self.particles.keys()
			ofile.write("Masses\n\n")
			for key in keys:
				if key == "water":
					ofile.write("1\t18.0154\n")
				elif key =="methane":
					ofile.write("2\t16.04\n\n")
				elif key == "silicon":
					ofile.write("1\t28.08\n")
				elif key == "oxygen":
					ofile.write("2\t15.9994\n\n")
				elif key == "copper":
					ofile.write("1\t63.546\n\n")

			counter = 1
			ofile.write("Atoms\n\n")
			for key in keys:
				for value in self.particles[key]:
					if key == "water":
						ofile.write("%d %d 1 %f %f %f\n" % (counter, value[3], value[0], value[1], value[2]))
					elif key == "methane":
						ofile.write("%d %d 2 %f %f %f\n" % (counter, value[3], value[0], value[1], value[2]))
					elif key == "silicon":
						ofile.write("%d %d 1 %f %f %f\n" % (counter, value[3], value[0], value[1], value[2]))
					elif key == "oxygen":
						ofile.write("%d %d 2 %f %f %f\n" % (counter, value[3], value[0], value[1], value[2]))
					elif key == "copper":
						ofile.write("%d %d 1 %f %f %f\n" % (counter, value[3], value[0], value[1], value[2]))	
					counter += 1
			ofile.write("\n")



def getBindDataFromVoronoiCrystals(crystal):
	myDtype = np.dtype([('position', np.float32, (3,)), ('color', np.float32, (3,)), ('radius', np.float32)])
	length = len(crystal.voronoiCenters)+sum([len(crystal.particles[key]) for key in crystal.particles.keys()])
	myData = np.empty(length, dtype=myDtype)
	counter = 0

	myData['position'][0:len(crystal.voronoiCenters)] = np.asarray(crystal.voronoiCenters)
	myData['color'][0:len(crystal.voronoiCenters)] = (0.5, 0.5, 0.5)
	myData['radius'][0:len(crystal.voronoiCenters)] = 5
	counter += len(crystal.voronoiCenters)

	mean_position = np.mean(myData['position'][0:len(crystal.voronoiCenters)], axis=0)

	for key, value in crystal.particles.items():
		myData['position'][counter:counter+len(value)] = np.asarray(value)[:,0:3]
		myData['color'][counter:counter+len(value)] = crystal.element_data[key]["color"]
		myData['radius'][counter:counter+len(value)] = crystal.element_data[key]["radius"]    
		counter += len(value)

	myData['position']-=mean_position
	return myData

def glumpyVisualize(voronoiCrystal):
	from glumpy import app, gloo, gl, data
	from glumpy.transforms import Position, Trackball
	from glumpy.graphics.filter import Filter
	vertex = """
uniform vec3 light_position;
attribute vec3 position;
attribute vec3 color;
attribute float radius;
varying float v_size;
varying vec3 v_color;
varying float v_radius;
varying vec4 v_eye_position;
varying vec3 v_light_direction;
void main (void)
{
    v_color = color;
    v_radius = radius;
    v_eye_position = <transform.trackball_view> *
                     <transform.trackball_model> *
                     vec4(position,1.0);
    v_light_direction = normalize(light_position);
    gl_Position = <transform(position)>;
    // stackoverflow.com/questions/8608844/...
    //  ... resizing-point-sprites-based-on-distance-from-the-camera
    vec4 p = <transform.trackball_projection> *
             vec4(radius, radius, v_eye_position.z, v_eye_position.w);
    v_size = 512.0 * p.x / p.w;
    gl_PointSize = v_size + 5.0;
}
"""

	fragment = """
#include "antialias/outline.glsl"
varying float v_size;
varying vec3 v_color;
varying float v_radius;
varying vec4 v_eye_position;
varying vec3 v_light_direction;
void main()
{
    vec2 P = gl_PointCoord.xy - vec2(0.5,0.5);
    float point_size = v_size  + 5.0;
    float distance = length(P*point_size) - v_size/2;
    vec2 texcoord = gl_PointCoord* 2.0 - vec2(1.0);
    float x = texcoord.x;
    float y = texcoord.y;
    float d = 1.0 - x*x - y*y;
    if (d <= 0.0) discard;
    float z = sqrt(d);
    vec4 pos = v_eye_position;
    pos.z += v_radius*z;
    vec3 pos2 = pos.xyz;
    pos = <transform.trackball_projection> * pos;
    gl_FragDepth = 0.5*(pos.z / pos.w)+0.5;
    vec3 normal = vec3(x,y,z);
    float diffuse = clamp(dot(normal, v_light_direction), 0.0, 1.0);
    vec4 color = vec4((0.5 + 0.5*diffuse)*v_color, 1.0);
    gl_FragColor = outline(distance, 1.0, 1.0, vec4(0,0,0,1), color);
    // gl_FragColor = color;
}
	"""
	window = app.Window(width=800, height=800, color=(1,1,1,1))
	particles = gloo.Program(vertex, fragment)
	particles['light_position'] = 0., 0., 2.
	particles["transform"] = Trackball(Position(), distance=500)
	particles.bind(getBindDataFromVoronoiCrystals(voronoiCrystal).view(gloo.VertexBuffer))


	@window.event
	def on_draw(dt):
	    window.clear()
	    particles.draw(gl.GL_POINTS)

	@window.event
	def on_init():
	    gl.glEnable(gl.GL_DEPTH_TEST)

	window.attach(particles["transform"])
	app.run()

def createDatafiles(Ls=[50, 100, 150, 200, 250, 300], Nvoronoi=3):
	for L in Ls:
		myVoronoiCrystals = VoronoiCrystals(N=Nvoronoi, type="regular", box=np.asarray([[0, L], [0, L], [0, L]]))
		myVoronoiCrystals.computeVoronoi()
		myVoronoiCrystals.loadCrystalFromFile("positions.json")
		N = len(myVoronoiCrystals.voronoiCenters)
		for i in range(N):
			Nreplications = int(float(L*2)/float(12)/Nvoronoi)
			print("Nreplications: ", Nreplications)
			myVoronoiCrystals.addParticlesWithReplicatedCrystal(i, Nreplications)
			print("Added crystal %d of %d" % (i, N))
		myVoronoiCrystals.dumpPositionsLammpsData("output/polycrystal_L%d_N%d.data" % (L, Nvoronoi))
		glumpyVisualize(myVoronoiCrystals)


def run_example(substance="copper"):
	L = 100
	myVoronoiCrystals = Polycrystal(N=2, type="regular", box = np.asarray([[0, L], [0, L], [0, L]]))
	myVoronoiCrystals.dumpPositionsLammpsData("test_voronoi_crystals.data")
	glumpyVisualize(myVoronoiCrystals)

if __name__=="__main__":
	run_example()
	

