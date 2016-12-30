"""
Created on Tue Dec 27 19:10:15 2016

@author: A.Vasilevich
"""


import math as mth
import numpy as np
import random as rndm
import os, errno

def silentremove(filename):
    try:
        os.remove(filename)
    except OSError as e: # this would be "except OSError, e:" before Python 2.6
        if e.errno != errno.ENOENT: # errno.ENOENT = no such file or directory
            raise # re-raise exception if a different error occured

silentremove("FirstTopoChipWallsTestNew.cif")
silentremove("FirstTopoChipFeaturesTestNew.cif")  
silentremove("GeneralParameterValues.txt")  
silentremove("ParameterValuesOfFeatures.txt")  
silentremove("PrimitivesComposingFeatures.txt")  
           
            
walls = open ("FirstTopoChipWallsTestNew.cif", 'w+')
feat = open ("FirstTopoChipFeaturesTestNew.cif", 'w+')
gen = open ("GeneralParameterValues.txt", 'w+')
para = open ("ParameterValuesOfFeatures.txt", 'w+')
prim = open ("PrimitivesComposingFeatures.txt", 'w+')

lastdef=[] # Here go the statements that have to go in the last definition statement of the CIF file;
#int ldindex =0; # Used to keep track of the current position in the lastdef array


"""This script generates the design of the Topochip with the first generation of features. 
Each TopoUnit contains a feature repeated multiple times. The feature is built from primitives with the shape of a circle, 
a isosceles triangle or a line (stretched rectangle). Multiple parameter values can be set to influence: 
    the number of primitives used; the ratio between circles, lines, and triangles; the size of the primitives; 
    the likeliness of primitives being placed aligned in the same direction 
"""
def IntToStr(x): # Returns the string representation of integer x bacup:http://stackoverflow.com/questions/2267362/how-to-convert-an-integer-in-any-base-to-a-string
    x2=int(x)
    s=str(x2)
    return s

def floorInt(x):# Returns the largest integer smaller or equal to x
    x2=float(x)
    xf=mth.floor(x2)+0.01
    return int(xf)

def round(x):# Rounds a double value to the nearest integer
    x2=float(x)
    xf=mth.floor(x2+0.5)
    return int(xf)

def rad(degrees):# Converts an angle from degrees to radians
    return mth.radians(degrees)

def degree(radians):# Converts an angle from radians to degrees
    return mth.degrees(radians)

    
""" Definitions of functions for generating random values"""

#check this https://studywolf.wordpress.com/2012/10/02/a-common-random-number-generator-for-c-and-python/
def RandInt(x, y): # Returns a random integer between x and y inclusive c
    x=int(x)
    y=int(y)
    if(y<x): # Assume x<=y
        temp = x
        x = y
        y = temp
        return rndm.randint(x, y)
    else:
        return rndm.randint(x, y)

 
def RandDouble(x, y): # Returns a random double between x and y inclusive
    x=float(x)
    y=float(y)
    if(y<x): # Assume x<=y
        temp = x
        x = y
        y = temp
        return rndm.uniform(x, y)
    else:
        return rndm.uniform(x, y)


""" The class for generating normally distributed values"""


class NormalDensity: # Polar Box-Mueller method
    def __init__( self, mu, sigma):        
        self.mu=mu
        self.sigma=sigma
        
    #Constructor
    def	NormalDensity(self, m, s):
        mu = float(m)
        sigma = float(s)
    # Probability distribution
    def p(self,x):
        return 1.0/self.sigma/mth.sqrt(2.0 * mth.pi)*mth.exp(-0.5*((float(x) - self.mu)/self.sigma) * ((float(x) - self.mu)/self.sigma))
    # Auxiliary function creates 2 N(mu,sigma) deviates at a time
    def n_(self):
    # Deviates are generated in pairs 
        isset = int(0) # Indicates whether a deviate is stored
        gset=float()# Used to store a deviate	
        if(isset==0):
            v1=v2=rsq=float()
            while(rsq >= 1.0 or rsq == 0.0):
                v1 = RandDouble(-1,1)
                v2 = RandDouble(-1,1)
                rsq = v1 * v1 + v2 * v2
            sfac = float(mth.sqrt(-2.0 * mth.log(rsq)/rsq))
            gset = v1 * sfac * self.sigma + self.mu
            isset = 1
            return v2 * sfac * self.sigma + self.mu
        else:
            isset = 0
            return gset

    def Sample(self, ns):
        ns=int(ns)
        retvals=[]
        idx = None
        i=0
        while (i<ns):
            retvals.append(self.n_())
            i+=1
        return retvals 
        
def main():
    """ Define fixed parameters (these will be the same for all features of this type, 
    lengths are in 0.01 microns)
    """
    print("Defining TopoChip Parameters")
    NumberOfTopoChipsHorizontally=int(3) # Number of TopoChips horizontally on one mask
    NumberOfTopoChipsVertically = int(2) # Number of TopoChips vertically on one mask
    HorizontalDistanceBetweenTopoChips = int(800000) # Horizontal distance between two Topochips on a mask
    VerticalDistanceBetweenTopoChips = int(800000) # Vertical distance between two Topochips on a mask
    
    NumberOfUnitsHorizontally = int(66) # Number of TopoUnits in one row of the TopoChip (the script assumes this number is even)
    gen.write("NumberOfUnitsHorizontally %d\n" % NumberOfUnitsHorizontally)
    NumberOfUnitsVertically = int(66) # Number of TopoUnits in one column of the TopoChip (the script assumes this number is even)
    
    gen.write("NumberOfUnitsVertically %d\n" % NumberOfUnitsVertically)
    UnitSideLength = int(29000) # Distance between two opposite walls of a unit
    gen.write("UnitSideLength %d\n" % UnitSideLength)
    WallThickness = int(1000) # Thickness of the walls surrounding the units
    gen.write("WallThickness %d\n" % WallThickness)
    EmptySpaceNextToWall = int(500) # The width of the space next to the unit wall in which no elements are allowed
    gen.write("EmptySpaceNextToWall %d\n" % EmptySpaceNextToWall)
    TopAngleOfTriangle = float(36.0) # Since the triangle is an isosceles triangle, the other two angles are equal to (180-TopAngleOfTriangle)/2
    gen.write("TopAngleOfTriangle %lf\n" % TopAngleOfTriangle)
    LineThickness = int(300) # The thickness of a line primitive
    gen.write("LineThickness %d\n" % LineThickness)
    BaseRotation = float(0.0) # The base rotation in degrees to which the axis of a primitive is aligned. For each primitive its individual rotation is added to the base rotation before it is placed. The direction of 0.0 degrees points to the right.
    gen.write("BaseRotation %lf\n" % BaseRotation)
    
    CircleApproximationPoints = int(12) # The number of points used to approximate a circle per micrometer of circle diameter
    gen.write("CircleApproximationPoints %d\n" % CircleApproximationPoints)
    SMALL = int(0) # Used as an index for vectors in which values are stored that are different for different feature sizes
    MEDIUM = int(1) # Used as an index for vectors in which values are stored that are different for different feature sizes
    LARGE = int(2) # Used as an index for vectors in which values are stored that are different for different feature sizes
    
    NumberOfFeatureSizes = int(3) # The number of different feature sizes
    FeatureSizes=np.zeros(NumberOfFeatureSizes,int) # The vector containing the different feature sizes (It is assumed in this script that feature sizes are chosen such that features fit an integer amount of times in the part of the TopoUnit where placement of features is allowed)
    FeatureSizes[SMALL] = 1000 # The length of one side of a small feature
    FeatureSizes[MEDIUM] = 2000 # The length of one side of a medium feature
    FeatureSizes[LARGE] = 2800 # The length of one side of a large feature
	
    MinimumNumberOfPrimitives=np.zeros(NumberOfFeatureSizes,int) # The vector containing the minimum number of primitives used
    MinimumNumberOfPrimitives[SMALL] = 3 # The minimum number of primitives in a small feature
    MinimumNumberOfPrimitives[MEDIUM] = 3 # The minimum number of primitives in a medium feature
    MinimumNumberOfPrimitives[LARGE] = 3 # The minimum number of primitives in a large feature
    
    MaximumNumberOfPrimitives=np.zeros(NumberOfFeatureSizes,int) # The vector containing the minimum number of primitives used
    MaximumNumberOfPrimitives[SMALL] = 5 # The maximum number of primitives in a small feature
    MaximumNumberOfPrimitives[MEDIUM] = 12 # The maximum number of primitives in a medium feature
    MaximumNumberOfPrimitives[LARGE] = 16 # The maximum number of primitives in a large feature
    
    MinimumLengthOfLinePrimitives=np.zeros(NumberOfFeatureSizes,int) # The vector containing the minimum length of a line primitive
    MinimumLengthOfLinePrimitives[SMALL] = 300 # The minimum length of a line primitive in a small feature 
    MinimumLengthOfLinePrimitives[MEDIUM] = 300 # The minimum length of a line primitive in a medium feature
    MinimumLengthOfLinePrimitives[LARGE] = 300 # The minimum length of a line primitive in a large feature
    
    MaximumLengthOfLinePrimitives=np.zeros(NumberOfFeatureSizes,int) # The vector containing the maximum length of a line primitive
    MaximumLengthOfLinePrimitives[SMALL] = 800 # The maximum length of a line primitive in a small feature 
    MaximumLengthOfLinePrimitives[MEDIUM] = 1600 # The maximum length of a line primitive in a medium feature
    MaximumLengthOfLinePrimitives[LARGE] = 2300 # The maximum length of a line primitive in a large feature
    
    MinimumCircleDiameter=np.zeros(NumberOfFeatureSizes,int) # The vector containing the minimum diameter of a circle primitive
    MinimumCircleDiameter[SMALL] = 300 # The minimum diameter of a circle primitive in a small feature 
    MinimumCircleDiameter[MEDIUM] = 300 # The minimum diameter of a circle primitive in a medium feature
    MinimumCircleDiameter[LARGE] = 300 # The minimum diameter of a circle primitive in a large feature
    
    MaximumCircleDiameter=np.zeros(NumberOfFeatureSizes,int) # The vector containing the maximum diameter of a circle primitive
    MaximumCircleDiameter[SMALL] = 400 # The maximum diameter of a circle primitive in a small feature 
    MaximumCircleDiameter[MEDIUM] = 1000 # The maximum diameter of a circle primitive in a medium feature
    MaximumCircleDiameter[LARGE] = 1000 # The maximum diameter of a circle primitive in a large feature
    
    MinimumTriangleShortestSideLength=np.zeros(NumberOfFeatureSizes,int) # The vector containing the minimum length of the shortest side of a circle primitive
    MinimumTriangleShortestSideLength[SMALL] = 300 # The minimum length of the shortest side of a triangle primitive in a small feature
    MinimumTriangleShortestSideLength[MEDIUM] = 300 # The minimum length of the shortest side of a triangle primitive in a medium feature
    MinimumTriangleShortestSideLength[LARGE] = 300 # The minimum length of the shortest side of a triangle primitive in a large feature
    
    MaximumTriangleShortestSideLength=np.zeros(NumberOfFeatureSizes,int) # The vector containing the maximum length of the shortest side of a circle primitive
    MaximumTriangleShortestSideLength[SMALL] = 400 # The maximum length of the shortest side of a triangle primitive in a small feature
    MaximumTriangleShortestSideLength[MEDIUM] = 800 # The maximum length of the shortest side of a triangle primitive in a medium feature
    MaximumTriangleShortestSideLength[LARGE] = 1000 # The maximum length of the shortest side of a triangle primitive in a large feature
    
    MinimumSpreadRot = float(0.0) # The minimum spread of the rotation of primitives in degrees. The standard deviation of the normal distribution used to determine the rotation of a primitive is equal to the spread parameter.
    MaximumSpreadRot = float(180.0) # The maximum spread of the rotation of primitives in degrees. The standard deviation of the normal distribution used to determine the rotation of a primitive is equal to the spread parameter.

    # Define variable parameters (these will be different for all features of this type)
    print("Defining Features Parameters")
    FeatureSize=int() # The maximum size of a feature 
    NumberOfPrimitives=int() # The number of primitives placed in the feature
    NumberOfCircles=int(0) # The number of primitives that is a circle. 
    NumberOfTriangles=int(0) # The number of primitives that is a triangle
    NumberOfLines=int(0) # The number of primitives that is a line
    CircleDiameter=float() # The diameter of the circles
    TriangleShortestSideLength=float() # The length of the shortest side of the triangles
    LineLength=float() # The length of the lines
    SpreadRot=float() # The standard deviation of the normal distribution that is used to determine how much the rotation of a primitive differs from the base rotation. The lower this number, the more aligned all primitives are. The base axis of a line is along the line and the base axis of a triangle is from the middle of the short side to the top corner.                   

    # Define derived (linked) and auxiliary parameters
    print("Defining derived and auxiliary Parameters")
    FractionOfCircles=float() # The fraction of primitives that is a circle. (Derived from the number of primitives and the number of circles) 
    FractionOfTriangles=float() # The fraction of primitives that is a triangle (Derived from the number of primitives and the number of triangles)
    FractionOfLines=float() # The fraction of primitives that is a line (Derived from the number of primitives and the number of lines)
    FeatureRepetition=int() # The amount of features that are in one row (or one column) of a TopoUnit (Derived from the size of a feature, the size of a TopoUnit and the empty space next to the unit walls)
    FeatureSizeSelection=int() # Used to select one of the possible feature sizes (Linked to the feature size)
    i=j=k=m=n=int() # Used as counters
    DefinitionNumber=int(1) # Used to keep track of the current definition used to write the CIF file
    FeatureNumber=int(1)
    TotalNumberOfCircles=int(0)
    TotalNumberOfTriangles=int(0)
    TotalNumberOfLines=int(0)		
	

    #Set seed for random generator
    #srand((unsigned)time(0));
    rndm.seed(14121980);

    # Define length of one unit in the CIF file
    walls.write("(1 unit = 0.01 micron);\n")
    feat.write("(1 unit = 0.01 micron);\n")
	
    # Define parameters useful for determining where TopoChips are placed on the mask
    print("Defining TopoChip Positions on the mask")
	
    HorizontalChipSize = int(WallThickness + NumberOfUnitsHorizontally * (WallThickness + UnitSideLength))
    VerticalChipSize = int(WallThickness + NumberOfUnitsVertically * (WallThickness + UnitSideLength))
    TotalHorizontalLength = int(NumberOfTopoChipsHorizontally * HorizontalChipSize + (NumberOfTopoChipsHorizontally - 1) * HorizontalDistanceBetweenTopoChips)
    TotalVerticalLength = int(NumberOfTopoChipsVertically * VerticalChipSize + (NumberOfTopoChipsVertically - 1) * VerticalDistanceBetweenTopoChips)
    HorizontalStartPoint = int(-TotalHorizontalLength/2)
    VerticalStartPoint = int(-TotalVerticalLength/2)
	
    # Define the inside of all TopoUnits (this is used to make the walls for the TopoUnit
	
    # Make a definition of the inside of one TopoUnit
    print("Building walls")
	
    WallCoordinates=np.zeros((4,2),int) # First index is corner number, starting bottom left and counterclockwise, second is coordinate (x/y)
	
    WallCoordinates[0][0] = WallThickness
    WallCoordinates[0][1] = WallThickness
    WallCoordinates[1][0] = WallThickness + UnitSideLength
    WallCoordinates[1][1] = WallThickness
    WallCoordinates[2][0] = WallThickness + UnitSideLength
    WallCoordinates[2][1] = WallThickness + UnitSideLength
    WallCoordinates[3][0] = WallThickness
    WallCoordinates[3][1] = WallThickness + UnitSideLength
    walls.write("DS %d 1 1;\n" % DefinitionNumber)
    walls.write("L L1;\n")
    walls.write("P")
    i=0
    while (i<4):
        j=0
        while (j<2):
            walls.write(" %d" % WallCoordinates[i][j])
            j+=1
            print '.',
        i+=1
	
    walls.write(";\n")
    walls.write("DF;\n")
    DefinitionNumber+=1
	
    # Make a definition of a row of insides of TopoUnits
	
    walls.write("DS %d 1 1;\n" % DefinitionNumber)
    WallTransX = int(WallThickness + UnitSideLength)
    i=0
    while (i<NumberOfUnitsHorizontally):
        walls.write("C %d T %d 0;\n" % (DefinitionNumber - 1, i * WallTransX))
        i+=1
    walls.write("DF;\n")
    DefinitionNumber+=1
	
    # Make a definition of TopoChip of insides of TopoUnits
	
    walls.write("DS %d 1 1;\n" % DefinitionNumber)
    WallTransY = int(WallThickness + UnitSideLength)
    i=0
    while (i<NumberOfUnitsVertically):
        walls.write( "C %d T 0 %d;\n" % (DefinitionNumber - 1, i * WallTransY))
        i+=1
    walls.write("DF;\n")
    DefinitionNumber+=1
	
    # Make a definition of insides of TopoUnits of multiple TopoChips
    print("Multiply TopoChip on the mask")
	
    walls.write("DS %d 1 1;\n" % DefinitionNumber)
    m=0
    while (m<NumberOfTopoChipsHorizontally):
        n=0 
        while (n<NumberOfTopoChipsVertically):
            walls.write( "C %d T %d %d;\n" % (DefinitionNumber - 1, HorizontalStartPoint + m * (HorizontalChipSize + HorizontalDistanceBetweenTopoChips), VerticalStartPoint + n * (VerticalChipSize + VerticalDistanceBetweenTopoChips)))
            n+=1
            print '.',
        m+=1
    walls.write("DF;\n")
	
    walls.write("C %d;\n" % DefinitionNumber)
    walls.write("E\n")
	
	
    DefinitionNumber+=1
	
	
    """ Define horizontal walls
	
    int WallCoordinates[4][2]; # First index is corner number, starting bottom left and counterclockwise, second is coordinate (x/y)
	
    WallCoordinates[0][0] = 0;
    WallCoordinates[0][1] = 0;
    WallCoordinates[1][0] = WallThickness + NumberOfUnitsHorizontally * (WallThickness + UnitSideLength);
    WallCoordinates[1][1] = 0;
    WallCoordinates[2][0] = WallThickness + NumberOfUnitsHorizontally * (WallThickness + UnitSideLength);
    WallCoordinates[2][1] = WallThickness;
    WallCoordinates[3][0] = 0;
    WallCoordinates[3][1] = WallThickness;
	
    fprintf(walls,"DS %d 1 1;\n", DefinitionNumber);
    fprintf(walls,"L L1;\n");
    fprintf(walls,"P");
    for(i=0;i<4;i++)
    {
         for(j=0;j<2;j++)
         {
              fprintf(walls," %d", WallCoordinates[i][j]);
              }
    }
	
    fprintf(walls,";\n");
    fprintf(walls,"DF;\n");
	
    for(i=0;i<=NumberOfUnitsVertically;i++)
    {
         int TranslationY = i * (WallThickness + UnitSideLength);
         fprintf(walls,"C %d T 0 %d;\n", DefinitionNumber, TranslationY);
    }
    DefinitionNumber++;    
	
    # Define vertical walls
    WallCoordinates[0][0] = 0;
    WallCoordinates[0][1] = 0;
    WallCoordinates[1][0] = WallThickness;
    WallCoordinates[1][1] = 0;
    WallCoordinates[2][0] = WallThickness;
    WallCoordinates[2][1] = WallThickness + NumberOfUnitsVertically * (WallThickness + UnitSideLength);
    WallCoordinates[3][0] = 0;
    WallCoordinates[3][1] = WallThickness + NumberOfUnitsVertically * (WallThickness + UnitSideLength);
	
    fprintf(walls,"DS %d 1 1;\n", DefinitionNumber);
    fprintf(walls,"L L1;\n");
    fprintf(walls,"P");
    for(i=0;i<4;i++)
    {
         for(j=0;j<2;j++)
         {
          fprintf(walls," %d", WallCoordinates[i][j]);
          }
        }
	
    fprintf(walls,";\n");
    fprintf(walls,"DF;\n");
    for(i=0;i<=NumberOfUnitsHorizontally;i++)
    {
     int TranslationX = i * (WallThickness + UnitSideLength);
     fprintf(walls,"C %d T %d 0;\n", DefinitionNumber, TranslationX);
    }
    DefinitionNumber++;
    fprintf(walls,"E");
    """
	    

    """ Creating the features
    TopoUnits are filled with a repeating feature and each TopoUnit is copied once to the quadrant diagonally opposite on the TopoChip
    """
    print("Create Fetures unit by unit")
    row = 0
    while (row < NumberOfUnitsVertically):
        print ("Row Number %d" % row)
        col = 0
        while (col < NumberOfUnitsHorizontally/2):
            print ("Col Number %d" % col)
            # Features are generated for the left half of the TopoChip and each TopoUnit is copied one time to the right half
            
            if((row == 0 and col == 0) or (row == NumberOfUnitsVertically/4 and col == NumberOfUnitsHorizontally/4)): # There are four empty squares on the TopoChip
                col+=1
                continue
            
            para.write("FeatureNumber %d\n" % FeatureNumber)
            prim.write("FeatureNumber %d\n" % FeatureNumber)
            FeatureNumber+=1
            para.write("row %d\n" % row)
            para.write("col %d\n" % col)
            prim.write("row %d\n" % row)
            prim.write("col %d\n" % col)
	
            # Test Output
            print("Creating TopoUnit (%d,%d)\n" % (row, col)) 
	
            # Generate parameter values for this feature
            FeatureSizeSelection = floorInt(NumberOfFeatureSizes*RandDouble(0.0,1.0)) # One of the feature sizes is selected
            FeatureSize = FeatureSizes[FeatureSizeSelection] # The corresponding feature size is retrieved from the vector with feature sizes
            para.write("FeatureSize %d\n" % FeatureSize)
            prim.write("FeatureSize %d\n" % FeatureSize)
            FeatureRepetition = (UnitSideLength - 2 * EmptySpaceNextToWall)/FeatureSize # Compute the amount of times that a feature is repeated in one row (or column) of a TopoUnit
            para.write("FeatureRepetition %d\n" % FeatureRepetition)
            
            NumberOfPrimitives = floorInt(MinimumNumberOfPrimitives[FeatureSizeSelection] + (MaximumNumberOfPrimitives[FeatureSizeSelection] - MinimumNumberOfPrimitives[FeatureSizeSelection] + 1) * RandDouble(0,1)) # The number of primitives placed in the feature is chosen
            para.write("NumberOfPrimitives %d\n" % NumberOfPrimitives)
            prim.write("NumberOfPrimitives %d\n" % NumberOfPrimitives)
	
            # Determine whether each of the primitive types is represented
	
            UseCircles=int()
            UseTriangles=int()
            UseLines=int()
            NumberOfDifferentPrimitiveTypes=int()
            while(NumberOfDifferentPrimitiveTypes == 0):
                UseCircles = floorInt(2.0*RandDouble(0.0,1.0))
                UseTriangles = floorInt(2.0*RandDouble(0.0,1.0))
                UseLines = floorInt(2.0*RandDouble(0.0,1.0))
                NumberOfDifferentPrimitiveTypes = UseCircles + UseTriangles + UseLines
            
            if(NumberOfDifferentPrimitiveTypes >= 3): # All three types of primitives can be used
            # The division of primitives over circles, triangles and lines is chosen in such a way that each possible division is as likely. At least one primitive is chosen for each type
	        
                NumberOfPossibleDivisions = int((NumberOfPrimitives * (NumberOfPrimitives - 3) + 2)/2) # Number of ways to divide the primitives over the types, while having at least one primitive of each type
                SelectDivision = int(floorInt(NumberOfPossibleDivisions * RandDouble(0,1))) # Used to select one of the possible divisions
                NumberOfCircles = 1
                while(SelectDivision > (NumberOfPrimitives - 3) - (NumberOfCircles - 1)): # Choose the number of circles
                    SelectDivision -= (NumberOfPrimitives - 3) - (NumberOfCircles - 1) + 1
                    NumberOfCircles+=1
                NumberOfTriangles = 1 + floorInt((NumberOfPrimitives - NumberOfCircles - 1) * RandDouble(0,1)) # Choose the number of triangles (at least 1, at most the remaining number of primitives - 1 
                NumberOfLines = NumberOfPrimitives - NumberOfCircles - NumberOfTriangles # Compute the number of lines  
            
            if(NumberOfDifferentPrimitiveTypes == 2): # Two out of three primitive types can be used
                if(UseCircles > 0):
                    NumberOfCircles = 1 + floorInt((NumberOfPrimitives - 1) * RandDouble(0,1)) # Choose the number of circles (at least 1, at most the remaining number of primitives - 1
                else: NumberOfCircles = 0
	        
                if(UseTriangles > 0):
                    if(UseCircles > 0): NumberOfTriangles = NumberOfPrimitives - NumberOfCircles # Compute the number of triangles
                    else: NumberOfTriangles = 1 + floorInt((NumberOfPrimitives - 1) * RandDouble(0,1)) # Choose the number of triangles (at least 1, at most the remaining number of primitives - 1
                else: NumberOfTriangles = 0
                if(UseLines > 0):
                    NumberOfLines = NumberOfPrimitives - NumberOfCircles - NumberOfTriangles
                else: NumberOfLines = 0
                
            if(NumberOfDifferentPrimitiveTypes == 1): # Only one primitive type is used
                if(UseCircles > 0): NumberOfCircles = NumberOfPrimitives
                else: NumberOfCircles = 0	        	
                if(UseTriangles > 0): NumberOfTriangles = NumberOfPrimitives
                else: NumberOfTriangles = 0
                if(UseLines > 0): NumberOfLines = NumberOfPrimitives
                else: NumberOfLines = 0
	
            TotalNumberOfCircles += NumberOfCircles
            TotalNumberOfTriangles += NumberOfTriangles
            TotalNumberOfLines += NumberOfLines
	
            para.write("NumberOfCircles %d\n" % NumberOfCircles)
            para.write("NumberOfTriangles %d\n" % NumberOfTriangles)
            para.write("NumberOfLines %d\n" % NumberOfLines)
	
            FractionOfCircles = (NumberOfCircles*1.0)/(NumberOfPrimitives*1.0) # Compute the fraction of circles
            FractionOfTriangles = (NumberOfTriangles*1.0)/(NumberOfPrimitives*1.0) # Compute the fraction of triangles
            FractionOfLines = (NumberOfLines*1.0)/(NumberOfPrimitives*1.0) # Compute the fraction of lines
	
            para.write("FractionOfCircles %lf\n" % FractionOfCircles)
            para.write("FractionOfTriangles %lf\n" % FractionOfTriangles)
            para.write("FractionOfLines %lf\n" % FractionOfLines) 
	
            CircleDiameter = MinimumCircleDiameter[FeatureSizeSelection] + (MaximumCircleDiameter[FeatureSizeSelection] - MinimumCircleDiameter[FeatureSizeSelection]) * RandDouble(0,1)  # The diameter of the circle primitives is chosen
            para.write("CircleDiameter %lf\n" % CircleDiameter)
	
            TriangleShortestSideLength = MinimumTriangleShortestSideLength[FeatureSizeSelection] + (MaximumTriangleShortestSideLength[FeatureSizeSelection] - MinimumTriangleShortestSideLength[FeatureSizeSelection]) * RandDouble(0,1) # The length of the shortest side of the triangle primitives is chosen
            para.write("TriangleShortestSideLength %lf\n" % TriangleShortestSideLength)
	
            LineLength = MinimumLengthOfLinePrimitives[FeatureSizeSelection] + (MaximumLengthOfLinePrimitives[FeatureSizeSelection] - MinimumLengthOfLinePrimitives[FeatureSizeSelection]) * RandDouble(0,1) # The length of the line primitives is chosen
            para.write("LineLength %lf\n" % LineLength)
	
            SpreadRot = MinimumSpreadRot + (MaximumSpreadRot - MinimumSpreadRot) * RandDouble(0,1) # Choose the the spread in rotation
            para.write("SpreadRot %lf\n" % SpreadRot)
	
            # Create samples from the normal distribution for the spread in rotation
	
            NDR=NormalDensity(BaseRotation,SpreadRot)
            ValuesR = NDR.Sample(NumberOfPrimitives)
            VectorPosR = int(0) # Used to keep track of which values for the rotation have been used already
            
            # Start a new definition for the current feature and print the layer number
	
            feat.write("DS %d 1 1;\n" % DefinitionNumber)
            feat.write("L L2;\n")
            	
            # Test Output
            print("Generating circles\n");
            i=0
            while(i<NumberOfCircles): # Place all circles
                print '.',
                CircleRadius = float(CircleDiameter/2)
                centerXPos=centerYPos=float()
                # Determine the extreme values where the center of a circle can just be placed
                
                minXPos = float(CircleRadius)
                maxXPos = float(FeatureSize - CircleRadius)
                minYPos = float(CircleRadius)
                maxYPos = float(FeatureSize - CircleRadius)
            	        
                # Choose X coordinate
            	        
                centerXPos = minXPos + (maxXPos - minXPos) * RandDouble(0,1)
            	        	        
                # Choose Y coordinate
            	        
                centerYPos = minYPos + (maxYPos - minYPos) * RandDouble(0,1)
            	        
                prim.write("C %lf %lf %lf\n" % (centerXPos, centerYPos, CircleRadius))
            	        
                # Determine the position of the circle in the lower left feature of the TopoUnit
            	        
                CircleBaseXPosition = int(WallThickness + EmptySpaceNextToWall + col * (UnitSideLength + WallThickness)) # The X position of the bottom left corner of the TopoUnit where placement of features is allowed for the base TopoUnit
                CircleBaseYPosition = int(WallThickness + EmptySpaceNextToWall + row * (UnitSideLength + WallThickness)) # The Y position of the bottom left corner of the TopoUnit where placement of features is allowed for the base TopoUnit
                CircleXPos = float(CircleBaseXPosition + centerXPos)
                CircleYPos = float(CircleBaseYPosition + centerYPos)
            	        
                # Write to the CIF file
            	        
                feat.write("P")
            		
                NumberOfPoints = int(round(CircleApproximationPoints * CircleDiameter * 1.0/ 100))
                #printf("Number of points = %d\n", NumberOfPoints);
                print("Calculating Distances")
                m=0
                while(m < NumberOfPoints):
                    print '.',
                    XCurr = round(CircleXPos + (CircleRadius * mth.cos (m * 2 * mth.pi/ NumberOfPoints)))
                    YCurr = round(CircleYPos + (CircleRadius * mth.sin (m * 2 * mth.pi/ NumberOfPoints)))
                    feat.write(" %d %d" % (XCurr, YCurr))
                    m+=1
                feat.write(";\n")
                i+=1
                
            # Test Output
            print("Placing all triangles\n")
            #printf("Number of triangles = %d\n", NumberOfTriangles);
                
            i=0
            while(i < NumberOfTriangles): # Place all triangles
                print '.',
                TriangleHeight = float((TriangleShortestSideLength * 1.0/2)/mth.tan(rad(TopAngleOfTriangle/2)))
                CGXPos=CGYPos=minCGXPos=maxCGXPos=minCGYPos=maxCGYPos=rot=float()
            	        
                # Choose a rotation
            	        
                rot = ValuesR[VectorPosR]
                VectorPosR+=1
                if(VectorPosR >= len(ValuesR)): # Generate new deviates when all deviates are used
                    # Generate new values
            	               
                    ValuesR = NDR.Sample(NumberOfPrimitives)
                    VectorPosR = 0
            
            	        
                    # Determine positions of triangle corners relative to CIA
            	                        
                    xPosA = TriangleHeight*mth.cos(rad(rot))
                    yPosA = TriangleHeight*mth.sin(rad(rot))
                    xPosB = (TriangleShortestSideLength*1.0/2)*-mth.sin(rad(rot))
                    yPosB = (TriangleShortestSideLength*1.0/2)*mth.cos(rad(rot))
                    xPosC = (TriangleShortestSideLength*1.0/2)*mth.sin(rad(rot))
                    yPosC = (TriangleShortestSideLength*1.0/2)*-mth.cos(rad(rot))
            	                
                    # Determine positions closest to cell wall
            	                
                    lowX = min(xPosA,min(xPosB,xPosC))
                    highX = max(xPosA,max(xPosB,xPosC))
                    lowY = min(yPosA,min(yPosB,yPosC))
                    highY = max(yPosA,max(yPosB,yPosC))
            	                
                    # Test output
                    #printf("rot = ", rot, "\nlowX = ",lowX, "\nhighX = ",highX, "\nlowY = ", lowY, "\nhighY = ", highY,"\n");
            	                                        
                    # Determine extreme positions where the CIA of the triangle can just be placed
            	                        
                    minCIAXPos = - lowX
                    maxCIAXPos = FeatureSize - highX
                    minCIAYPos = - lowY
                    maxCIAYPos = FeatureSize - highY
            	                        
                    # Make corrections because CG is used to position triangles instead of CIA
            	                        
                    minCGXPos = minCIAXPos + (TriangleHeight*1.0/3)*mth.cos(rad(rot))
                    maxCGXPos = maxCIAXPos + (TriangleHeight*1.0/3)*mth.cos(rad(rot))
                    minCGYPos = minCIAYPos + (TriangleHeight*1.0/3)*mth.sin(rad(rot))
                    maxCGYPos = maxCIAYPos + (TriangleHeight*1.0/3)*mth.sin(rad(rot))
            	                        
                    # Test output
                    #printf("minCIAXPos = ", minCIAXPos, "\nmaxCIAXPos = ",maxCIAXPos, "\nminCIAYPos = ",minCIAYPos, "\nmaxCIAYPos = ", maxCIAYPos,"\n");
                    #printf("minCGXPos = ", minCGXPos, "\nmaxCGXPos = ",maxCGXPos, "\nminCGYPos = ",minCGYPos, "\nmaxCGYPos = ", maxCGYPos,"\n");        
            	        
                    # Choose X coordinate
            	        
                    CGXPos = minCGXPos + (maxCGXPos - minCGXPos) * RandDouble(0,1)
            	        
                    # Choose Y coordinate
            	        
                    CGYPos = minCGYPos + (maxCGYPos - minCGYPos) * RandDouble(0,1)
            	        
                    # Compute the position of each of the triangle corners. A is the top corner, B and C follow in counterclockwise order
                    BCRadius = mth.sqrt((TriangleHeight*1.0/3) * (TriangleHeight*1.0/3) + (TriangleShortestSideLength*1.0/2) * (TriangleShortestSideLength*1.0/2)) # The distance from the triangle center to corners B and C
                    BInitialAngle = np.arctan2((TriangleShortestSideLength*1.0/2), -(TriangleHeight*1.0/3))
                    CInitialAngle = np.arctan2(-(TriangleShortestSideLength*1.0/2), -(TriangleHeight*1.0/3))
                    #printf("TriangleHeight*2/3 = %lf\nBCRadius = %lf\nBInitialAngle= %lf\nCInitialAngle = %lf\n", TriangleHeight*2.0/3, BCRadius, BInitialAngle, CInitialAngle);
                    TriangleBaseXPosition = WallThickness + EmptySpaceNextToWall + col * (UnitSideLength + WallThickness) # The X position of the bottom left corner of the TopoUnit where placement of features is allowed for the base TopoUnit
                    TriangleBaseYPosition = WallThickness + EmptySpaceNextToWall + row * (UnitSideLength + WallThickness) # The Y position of the bottom left corner of the TopoUnit where placement of features is allowed for the base TopoUnit
            	        
                    AXPos = round(CGXPos + (TriangleHeight*2.0/3) * mth.cos(rad(rot))) + TriangleBaseXPosition
                    AYPos = round(CGYPos + (TriangleHeight*2.0/3) * mth.sin(rad(rot))) + TriangleBaseYPosition
                    BXPos = round(CGXPos + BCRadius * mth.cos(rad(rot) + BInitialAngle)) + TriangleBaseXPosition
                    BYPos = round(CGYPos + BCRadius * mth.sin(rad(rot) + BInitialAngle)) + TriangleBaseYPosition
                    CXPos = round(CGXPos + BCRadius * mth.cos(rad(rot) + CInitialAngle)) + TriangleBaseXPosition
                    CYPos = round(CGYPos + BCRadius * mth.sin(rad(rot) + CInitialAngle)) + TriangleBaseYPosition
            		
                    prim.write("T %d %d %d %d %d %d\n" % (AXPos - TriangleBaseXPosition, AYPos - TriangleBaseYPosition, BXPos - TriangleBaseXPosition, BYPos - TriangleBaseYPosition, CXPos - TriangleBaseXPosition, CYPos - TriangleBaseYPosition))
            		
                    # Write to the CIF file
            	        
                    feat.write( "P %d %d %d %d %d %d;\n"% ( AXPos, AYPos, BXPos, BYPos, CXPos, CYPos))
                i+=1
            			  
	
            # Test Output
            print("Placing all lines\n")
	
            i=0
            while(i < NumberOfLines): # Place all lines 
                print '.',
                rotL=float()
                DistCenterToCorner = mth.sqrt((LineThickness*1.0/2)*(LineThickness*1.0/2)+(LineLength*1.0/2)*(LineLength*1.0/2))
                InitRotTL = degree(mth.atan((LineThickness*1.0/2)/(LineLength*1.0/2)))
                InitRotBL = 180-InitRotTL
                InitRotBR = 180+InitRotTL
                InitRotTR = 360-InitRotTL
	
                LXPos=LYPos=LminCXPos=LmaxCXPos=LminCYPos=LmaxCYPos=float()
	        
                # Choose rotation
	        
                rotL = ValuesR[VectorPosR]
                VectorPosR+=1
                if(VectorPosR >= len(ValuesR)): # Generate new deviates when all deviates have been used

                    # Generate new values
	
                    ValuesR = NDR.Sample(NumberOfPrimitives)
                    VectorPosR = 0

	        
                    # Determine positions of line corners relative to the center of the line
	                        
                xPosTL = DistCenterToCorner*mth.cos(rad(rotL+InitRotTL))
                yPosTL = DistCenterToCorner*mth.sin(rad(rotL+InitRotTL))
                xPosBL = DistCenterToCorner*mth.cos(rad(rotL+InitRotBL))
                yPosBL = DistCenterToCorner*mth.sin(rad(rotL+InitRotBL))
                xPosBR = DistCenterToCorner*mth.cos(rad(rotL+InitRotBR))
                yPosBR = DistCenterToCorner*mth.sin(rad(rotL+InitRotBR))
                xPosTR = DistCenterToCorner*mth.cos(rad(rotL+InitRotTR))
                yPosTR = DistCenterToCorner*mth.sin(rad(rotL+InitRotTR))
	                                            
                # Determine positions closest to cell wall
	                
                LlowX = min(xPosTL,min(xPosTR,min(xPosBL,xPosBR)))
                LhighX = max(xPosTL,max(xPosTR,max(xPosBL,xPosBR)))
                LlowY = min(yPosTL,min(yPosTR,min(yPosBL,yPosBR)))
                LhighY = max(yPosTL,max(yPosTR,max(yPosBL,yPosBR)))
	                
                # Test output                        
                #printf("\nxPosTL = ", xPosTL, "\nyPosTL = ", yPosTL, "\nxPosTR = ", xPosTR, "\nyPosTR = ", yPosTR, "\nxPosBL = ", xPosBL, "\nyPosBL = ", yPosBL, "\nxPosBR = ", xPosBR, "\nyPosBR = ", yPosBR, "\nLlowX = ", LlowX, "\nLhighX = ", LhighX, "\nLlowY = ", LlowY, "\nLhighY = ", LhighY);  
                #printf("rot = ", rot, "\nlowX = ",lowX, "\nhighX = ",highX, "\nlowY = ", lowY, "\nhighY = ", highY,"\n");
	                
                # Determine extreme positions where the line can just be placed                  
	                
                LminCXPos = - LlowX
                LmaxCXPos = FeatureSize - LhighX
                LminCYPos = - LlowY
                LmaxCYPos = FeatureSize - LhighY
	                
                # Test output                        
                #printf("\nLminCXPos = ", LminCXPos, "\nLmaxCXPos = ",LmaxCXPos, "\nLminCYPos = ",LminCYPos, "\nLmaxCYPos = ", LmaxCYPos,"\n");
                #printf("minCGXPos = ", minCGXPos, "\nmaxCGXPos = ",maxCGXPos, "\nminCGYPos = ",minCGYPos, "\nmaxCGYPos = ", maxCGYPos,"\n");
	        
                # Choose X coordinate
	        
                LXPos = LminCXPos + (LmaxCXPos - LminCXPos) * RandDouble(0,1)
	        
                # Choose Y coordinate
	        
                LYPos = LminCYPos + (LmaxCYPos - LminCYPos) * RandDouble(0,1)
	        	        	        
                LineBaseXPosition = WallThickness + EmptySpaceNextToWall + col * (UnitSideLength + WallThickness) # The X position of the bottom left corner of the TopoUnit where placement of features is allowed for the base TopoUnit
                LineBaseYPosition = WallThickness + EmptySpaceNextToWall + row * (UnitSideLength + WallThickness) # The Y position of the bottom left corner of the TopoUnit where placement of features is allowed for the base TopoUnit
	        
                TLXPos = round(LXPos + xPosTL) + LineBaseXPosition
                TLYPos = round(LYPos + yPosTL) + LineBaseYPosition
                BLXPos = round(LXPos + xPosBL) + LineBaseXPosition
                BLYPos = round(LYPos + yPosBL) + LineBaseYPosition
                BRXPos = round(LXPos + xPosBR) + LineBaseXPosition
                BRYPos = round(LYPos + yPosBR) + LineBaseYPosition
                TRXPos = round(LXPos + xPosTR) + LineBaseXPosition
                TRYPos = round(LYPos + yPosTR) + LineBaseYPosition
	        
                prim.write("L %d %d %d %d %d %d %d %d\n" % (TLXPos - LineBaseXPosition, TLYPos - LineBaseYPosition, BLXPos - LineBaseXPosition, BLYPos - LineBaseYPosition, BRXPos - LineBaseXPosition, BRYPos - LineBaseYPosition, TRXPos - LineBaseXPosition, TRYPos - LineBaseYPosition))
	        		
                # Write to the CIF file
	        
                feat.write("P %d %d %d %d %d %d %d %d;\n" % (TLXPos, TLYPos, BLXPos, BLYPos, BRXPos, BRYPos, TRXPos, TRYPos))
                i+=1
            
            feat.write("DF;\n")
            
            para.write("SymbolNumber %d\n" % DefinitionNumber)
            prim.write("SymbolNumber %d\n" % DefinitionNumber)
	
            print("Compute the translation from the base TopoUnit to the copied TopoUnit")
            """ Compute the translation from the base TopoUnit to the copied TopoUnit""" 
	
            UnitTransX = (NumberOfUnitsHorizontally/2) * (UnitSideLength + WallThickness)
            UnitTransY = (NumberOfUnitsVertically/2) * (UnitSideLength + WallThickness)
            if(row >= NumberOfUnitsVertically/2): UnitTransY -= NumberOfUnitsVertically * (UnitSideLength + WallThickness) # Make a correction if the row was in the top half of the TopoChip
	
            # Make a definition with a row of features
            DefinitionNumber+=1
            feat.write("DS %d 1 1;\n" % DefinitionNumber)
            m=0
            while(m<FeatureRepetition):
                print '.',
                TransX = m * FeatureSize
                feat.write("C %d T %d 0;\n" % (DefinitionNumber-1, TransX))
                m+=1

            feat.write("DF;\n")
	
            # Make a definition with a unit of features
            DefinitionNumber+=1
            feat.write("DS %d 1 1;\n" % DefinitionNumber)
            m=0
            while(m<FeatureRepetition):

                TransY = m * FeatureSize
                feat.write("C %d T 0 %d;\n" % (DefinitionNumber-1, TransY))
                m+=1
            feat.write("DF;\n")
	
            # Call the previously created definition and write it to lastcall
	
            #fprintf(feat, "C %d T 0 0;\n", DefinitionNumber);
            #fprintf(feat, "C %d T %d %d;\n", DefinitionNumber, UnitTransX, UnitTransY);
	
            temps="".join(["C ",IntToStr(DefinitionNumber),";\n"]) # the first call statement
            lastdef.append(temps)
	
            
            temps ="".join([ "C ", IntToStr(DefinitionNumber), " T ", IntToStr(UnitTransX), " ", IntToStr(UnitTransY),";\n"])   # the translated call statement
            lastdef.append(temps)
            #lastdef[ldindex++] = "C ".append(IntToStr(DefinitionNumber).c_str()).append(" T 0 0;\n");
            #lastdef[ldindex++] = "C " + IntToStr(DefinitionNumber).c_str() + " T 0 0;\n";
            #lastdef[ldindex++] = "C " + IntToStr(DefinitionNumber).c_str() + " T " + IntToStr(UnitTransX).c_str() + " " + IntToStr(UnitTransY).c_str() + ";\n";
	
            DefinitionNumber+=1
		
            col+=1	
        row+=1
	
    feat.write("DS %d 1 1;\n" % DefinitionNumber)
    m=0
    while(m<len(lastdef)):
        feat.write( "%s"  % str(lastdef[m]))
        #print lastdef[m]
        m+=1


    feat.write("DF;\n")
    #fprintf(feat, "C %d;\n", DefinitionNumber);	
    DefinitionNumber+=1
	
    # Make a new definition calling the TopoChip definition multiple times
	
    feat.write("DS %d 1 1;\n" % DefinitionNumber)
    print("Finishing")
    m=0
    while(m<NumberOfTopoChipsHorizontally):
        n=0
        while(n<NumberOfTopoChipsVertically):
            feat.write("C %d T %d %d;\n" % (DefinitionNumber - 1, HorizontalStartPoint + m * (HorizontalChipSize + HorizontalDistanceBetweenTopoChips), VerticalStartPoint + n * (VerticalChipSize + VerticalDistanceBetweenTopoChips)))
            n+=1
            print '.',
        m+=1
    feat.write("DF;\n")
	
    feat.write("C %d;\n" % DefinitionNumber)
	
	
    feat.write("E\n")
    print("TopoChip is Ready!")
	
    #printf("TotalNumberOfCircles = %d\n", TotalNumberOfCircles);
    #printf("TotalNumberOfTriangles = %d\n", TotalNumberOfTriangles);
    #printf("TotalNumberOfLines = %d\n", TotalNumberOfLines);
main()
walls.close()
feat.close()
gen.close()
para.close()
prim.close()
