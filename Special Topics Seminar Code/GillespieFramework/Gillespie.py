import random
import math

class StructuralFormula:
    def __init__(self,chemicalNumber,coefficient):
        self.chemicalNumber = chemicalNumber
        self.coefficient = coefficient
    def getNumber(self):
        return self.chemicalNumber
    def getCoefficient(self):
        return self.coefficient

class Reaction:
    def __init__(self,reactants,products, rate):
        self.reactants = reactants
        self.products = products
        self.rate = rate
    def getPropensity(self,chemicalList):
        a = self.rate
        for chemical in self.reactants:
            for i in range(chemical.getCoefficient()):
                a *= chemicalList[chemical.getNumber()]
        return a

    def performReaction(self,chemicalList):
        # repeated code... an opportunity for refactoring
        for chemical in self.reactants:
            chemicalList[chemical.getNumber()]-= chemical.getCoefficient()
        for chemical in self.products:
            chemicalList[chemical.getNumber()]+= chemical.getCoefficient()

class GillespieObserver:
    def __init__(self,fileName):
        self.file  = open(fileName,'w')
    def __del__(self):
        self.file.close()
    def takeMeasurement(self,gillespie):
        self.file.write(str(gillespie.getT())+" ")
        self.file.write(str(gillespie.getChemicals())+"\n")
    
class Gillespie:
    def __init__(self,initialValues,reactions):
        self.initialValues = initialValues
        self.chemicalList = initialValues
        self.reactions = reactions
        self.propensities = [0.0] * len(reactions) 
        self.totalP = 0
        self.T = 0
        self.observer = None
    def setObserver(self,obs):
        self.observer = obs 
    def reset(self):
        self.T = 0
        self.chemicalList = self.initialValues
    
    def getT(self):
        return self.T
    
    def getChemicals(self):
        return self.chemicalList

    def run(self,maxT):
        while self.T<maxT:
            self.step()

    def step(self):
        self.calculatePropensities()
        self.updateTime()
        RxNum = self.chooseReaction()
        self.reactions[RxNum].performReaction(self.chemicalList)
        if self.observer is not None:
            self.observer.takeMeasurement(self)
        
    def calculatePropensities(self):
        self.totalP = 0
        for i in range(len(self.reactions)):
            self.propensities[i] = self.reactions[i].getPropensity(self.chemicalList)
            self.totalP += self.propensities[i]
    
    def updateTime(self):
        dt = 1/self.totalP*math.log(1/random.random())
        self.T += dt
    
    def chooseReaction(self):
        r = random.random()
        RxNum = 0
        tally = 0
        for a in self.propensities:
            tally += a
            if tally/self.totalP>r:
                break
            RxNum+=1
        return RxNum

if __name__=="__main__":

    rx0 = Reaction([StructuralFormula(0,1),StructuralFormula(1,1)],[StructuralFormula(2,1)],1)
    rx1 = Reaction([StructuralFormula(2,1)],[StructuralFormula(0,1),StructuralFormula(1,1)],1)

    myGillespie = Gillespie([50,50,50],[rx0,rx1])

    myObs = GillespieObserver("data.txt")

    myGillespie.setObserver(myObs)

    myGillespie.run(10)




        



