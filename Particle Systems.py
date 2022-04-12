import math

class System:
    def __init__(self, particlesIn = []):
        self.particles = particlesIn
        self.eo = 8.854 * (10 ** -12)
        self.k = 9 * (10 ** 9) #1/(4 * math.pi * eo)

    def addParticle(self, particleIn):
        self.particles.append(particleIn)

    def calculateNetForces(self):
        netForces = [None] * len(self.particles)
        for i in range(0,len(self.particles)):
            vNet = Vector(xMagIn = 0, yMagIn = 0)
            for j in range(0,len(self.particles)):
                if(i == j):
                    continue
                vNet = vNet.vectorSum(self.vectorR12(j,i).normalize().scalarMultiplication(self.forceMagnitude(j,i)))
            netForces[i] = vNet
        return netForces

    def electricFieldAtPoint(self, pointIn):
        vNet = Vector(xMagIn = 0, yMagIn = 0)
        for i in range(0,len(self.particles)):
            #print(self.vectorR12P(i,pointIn).vectorSummary())
            vNet = vNet.vectorSum(self.vectorR12P(i,pointIn).normalize().scalarMultiplication(self.particles[i].electricFieldMagnitudeFunction()(self.distanceP(i,pointIn))))
        return vNet
        

    def forceMagnitude(self, particle1Ind, particle2Ind):
        ch1 = self.particles[particle1Ind].charge
        ch2 = self.particles[particle2Ind].charge
        return (self.k * ch1 * ch2) / (self.distance(particle1Ind, particle2Ind) ** 2)

    def distance(self, particle1Ind, particle2Ind):
        x1 = self.particles[particle1Ind].x
        y1 = self.particles[particle1Ind].y
        x2 = self.particles[particle2Ind].x
        y2 = self.particles[particle2Ind].y
        return math.sqrt( ((y2-y1) ** 2) + ((x2-x1) ** 2) )

    def distanceP(self, particle1Ind, pointIn):
        x1 = self.particles[particle1Ind].x
        y1 = self.particles[particle1Ind].y
        x2 = pointIn.x
        y2 = pointIn.y
        return math.sqrt( ((y2-y1) ** 2) + ((x2-x1) ** 2) )

    def angle(self, particle1Ind, particle2Ind):
        x1 = self.particles[particle1Ind].x
        y1 = self.particles[particle1Ind].y
        x2 = self.particles[particle2Ind].x
        y2 = self.particles[particle2Ind].y
        return math.atan2(y2-y1, x2-x1)

    def vectorR12(self,particle1Ind, particle2Ind):
        x1 = self.particles[particle1Ind].x
        y1 = self.particles[particle1Ind].y
        x2 = self.particles[particle2Ind].x
        y2 = self.particles[particle2Ind].y
        return Vector(xMagIn = x2-x1, yMagIn = y2-y1)

    def vectorR12P(self, particle1Ind, pointIn):
        x1 = self.particles[particle1Ind].x
        y1 = self.particles[particle1Ind].y
        x2 = pointIn.x
        y2 = pointIn.y
        return Vector(xMagIn = x2-x1, yMagIn = y2-y1)

class Point:
    def __init__(self,vectorIn):
        self.vector = vectorIn
        self.x = vectorIn.xMag
        self.y = vectorIn.yMag
        self.r = vectorIn.magnitude
        self.angle = vectorIn.direction

class Particle:
    def __init__(self, chargeIn, xIn = None, yIn = None, rIn = None, angleIn = None):
        self.charge = chargeIn
        self.k = 9 * (10 ** 9)
        if(xIn == None and yIn == None):
            self.x = rIn * math.cos(angleIn) #radians
            self.y = rIn * math.sin(angleIn)
            self.r = rIn
            self.angle = angleIn
        else:
            self.r = math.sqrt( ((yIn) ** 2) + ((xIn) ** 2) )
            self.angle = math.atan2(yIn, xIn)
            self.x = xIn
            self.y = yIn
            
    def electricFieldMagnitudeFunction(self):
        # A funcion of distance
        def fieldMagnitude(distanceIn):
            return ((self.k * self.charge) / (distanceIn ** 2))
        return fieldMagnitude
        

class Vector:
    def __init__(self, xMagIn = None, yMagIn = None, magnitudeIn = None, directionIn = None):
        if(xMagIn == None and yMagIn == None):
            self.xMag = magnitudeIn * math.cos(directionIn) #radians
            self.yMag = magnitudeIn * math.sin(directionIn)
            self.magnitude = magnitudeIn
            self.direction = directionIn
        else:
            self.magnitude = math.sqrt( ((yMagIn) ** 2) + ((xMagIn) ** 2) )
            self.direction = math.atan2(yMagIn, xMagIn)
            self.xMag = xMagIn
            self.yMag = yMagIn

    def vectorSum(self, vectorIn):#Return sum of self and param1
        return Vector(xMagIn = (self.xMag + vectorIn.xMag), yMagIn = (self.yMag + vectorIn.yMag))

    def normalize(self):
        return Vector(xMagIn = (self.xMag / self.magnitude), yMagIn = (self.yMag / self.magnitude))

    def vectorSummary(self):
        return f"{self.xMag} i, {self.yMag} j -- Mag: {self.magnitude}, Angle: {self.direction * 180 / math.pi}"

    def scalarMultiplication(self, scalar):
        return Vector(magnitudeIn = self.magnitude * scalar, directionIn = self.direction)

#print(p2.electricFieldFunction()(4))
#print(p1.electricFieldFunction()(5))    

e = -1.60217662  * (10 ** -19)
   
p1 = Particle(6   * (10 ** -9), 0.6, 0)
p2 = Particle(-4  * (10 ** -9), 0.6, 0.8)
p3 = Particle(1, 0,0)
sys = System(particlesIn = [p1,p2,p3])
sys2= System(particlesIn = [p1,p2])

#print(sys.calculateNetForces()[2].vectorSummary())
print(p1.electricFieldMagnitudeFunction()(0.23175))
print(p2.electricFieldMagnitudeFunction()(1.16-0.23175))
print(sys2.electricFieldAtPoint(Point(Vector(xMagIn = 0, yMagIn = 0))).vectorSummary())
#print(sys2.electricFieldAtPoint(Point(Vector(xMagIn = -0.2, yMagIn = 0))).scalarMultiplication(e).vectorSummary())

''' 
mag23 = sys.forceMagnitude(1,2)
print(mag23)
print(sys.vectorR12(1,2).normalize().vectorSummary())
v23 = sys.vectorR12(1,2).normalize().scalarMultiplication(mag23)
print(v23.vectorSummary())
print("")

mag13 = sys.forceMagnitude(0,2)
print(mag13)
print(sys.vectorR12(0,2).normalize().vectorSummary())
v13 = sys.vectorR12(0,2).normalize().scalarMultiplication(mag13)
print(v13.vectorSummary())

print("")
nv3 = v23.vectorSum(v13)
print(nv3.vectorSummary())
#'''


