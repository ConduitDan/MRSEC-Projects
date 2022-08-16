#from dataclasses import Field
import re
import SPO
import fileinput
# TODO:
# [x] Config File Parser
# [ ] Ensemble Writer
# [x] Log File Parser



#Parser should be simple, put complexity in the file spec
class SPOEnsembleLogParser:
    def __init__(self,fileName):
        self.file = fileinput.input(fileName,inplace=True)

    def __del__(self):
        self.file.close()
    
    def parseEnsembleLog(self,ensembleID):
        finished = False
        ensembleSize = -1
        finishedCounter = 0
        for line in self.file:
            #look for lines that end in "Running " 
            runningJob = re.match("Run (\d+) of \d+ Running",line)
            if runningJob:
                #we found a job that is still running
                if int(runningJob.group(1)) == ensembleID:
                    # we are the responsible of that job,
                    # update that is is finished
                    line = line.replace("Running","Finished")
                    finishedCounter +=1
            finishedJob = re.match("Run \d+ of (\d+) Finished",line)
            if finishedJob:
                finishedCounter +=1
                ensembleSize = int(finishedJob.group(1))
            if not (finishedJob or runningJob):
                raise Exception("Could not read line %s in ensemble file"%line)

            print(line,end='')

        return finishedCounter == ensembleSize

#the parser works by filling out a file spec
# a file Spec is made out of fieldSpecs.
# field specs either are header and a rule to parse the lines after
# or 

# two options here, inherite from fileSpec every time I need a thing
# or make a 
class FileSpecFactory:
    def __init__(self):
        self.floatMatch = "([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)"
        
    def logFileSpec(self):
        logSpec = FileSpec()



        # post processing for reading from a log file
        def parseParams(paramStr):
            paramList = []
            paramPairs = paramStr.split(",")
            for param in paramPairs:
                (name,value) = param.split(':')
                paramList.append([name,float(value)])
            return paramList

        # how to use the post processing for most of the log
        def logFieldPP(self):
            for i in range(len(self.value)):
                self.value[i] = [parseParams(self.value[i][0]),float(self.value[i][1])]
        
        #how to use the post processing for the last part of the log
        def lastLogFieldPP(self):
            self.value = [parseParams(self.value[1]),int(self.value[0])]
                

        logField = FieldSpec(pattern="Step\s+\d+\s+{(.*)}\s+Residue:"+self.floatMatch,\
        postProcessing = logFieldPP)
        logSpec.addFieldSpec(logField)
        
        lastLogField = FieldSpec(pattern="Step\s+(\d+)\s+{(.*)}",multiLine=False,\
            postProcessing = lastLogFieldPP)
        logSpec.addFieldSpec(lastLogField)


        return logSpec
        
    def configFileSpec(self):
        #create the empty spec
        configSpec = HeaderFieldFileSpec()
        # the name field take everything on the next line 
        # so we can leave it blank
        nameField = HeaderFieldSpec("Name:")
        configSpec.addFieldSpec(nameField)

        simulationField = HeaderFieldSpec("Simulation:")
        configSpec.addFieldSpec(simulationField)

        partitionField = HeaderFieldSpec("Partition:")
        configSpec.addFieldSpec(partitionField)
        accountField = HeaderFieldSpec("Account:")
        configSpec.addFieldSpec(accountField)
        
        extraScriptCommandsField = HeaderFieldSpec("Extra Commands:",multiLine = True)
        configSpec.addFieldSpec(extraScriptCommandsField)

        methodField = HeaderFieldSpec("Method:")
        configSpec.addFieldSpec(methodField)
        # dataFileField matches file names
        # look for some not white space then a period then some word characters
        fileNameRegex = "(\S*\.\w+)" 
        dataFileField = HeaderFieldSpec("Data:",\
            pattern=fileNameRegex + "\s+" + fileNameRegex,multiLine=True)
        configSpec.addFieldSpec(dataFileField)
        # parameters look for alphanumerics and underscores followed by a colon
        # then a float
        def paramPostProcess(self):
            for i in range(len(self.value)):
                self.value[i] = [self.value[i][0],float(self.value[i][1])]

        paramField = HeaderFieldSpec("Parameters:",\
            pattern = "(.*):\s*"+self.floatMatch,multiLine = True,postProcessing=paramPostProcess)
        configSpec.addFieldSpec(paramField)

        # now the more complicated fields
        # runs on need to make sure we have Either Desktop or HPCC and might 
        # have a number after it. So we define some post processing to
        # check these things
        def runsOnPostProcess(self):
            # check if we match "Desktop"
            runsOnMatch = re.match("Desktop",self.value)
            if runsOnMatch:
                self.value = [SPO.SPORunsOn.Desktop,1]
            else:
            # check if we match HPCC
                line = self.value
                runsOnMatch = re.match("HPCC",line)
                if runsOnMatch:
                    self.value = [SPO.SPORunsOn.HPCC]
                    #check if we have a number after
                    runsOnSplit = line.split(' ')
                    if len(runsOnSplit)>1:
                        # if so append it
                        self.value.append(int(runsOnSplit[1]))
                    else:
                        # if not append 1
                        self.value.append(1)
                else:
                    raise Exception("Field 'Runs On' is formated incorrectly.\
                        Expected 'Desktop' or 'HPCC' possibly followed by a number")
        
        runsOnField = HeaderFieldSpec("Runs On:",postProcessing=runsOnPostProcess)
        configSpec.addFieldSpec(runsOnField)

        return configSpec

class FileSpec:
    # a Plane old file that tries to read a single field
    def __init__(self):
        self.fields = []
    def addFieldSpec(self,field):
        self.fields.append(field)
    def readFile(self,fileName):
        #for each line try each field
        parser = SPOFileParser(fileName)
        # we'll rely on the parser raising a stop iteration when its done
        # with the file
        for line in parser:
            self.findAndFillSpec(line,parser)

        return self.makeValues()

    def findAndFillSpec(self,line,parser):
        # loop though all the fields until we get one that can read the line
        for field in self.fields:
            if field.read(line,parser):
                break

    def makeValues(self):
        return list(field.value for field in self.fields)
    

class HeaderFieldFileSpec(FileSpec):
    def __init__(self):
        self.fields = {}
    def addFieldSpec(self,field):
        self.fields[field.header]=field

    def findAndFillSpec(self,line,parser):
        # Find the field with a header that matches and read it
        if line in self.fields.keys():
            nextLine = next(parser)
            if not self.fields[line].read(nextLine,parser):
                raise Exception("Field %s is formatted incorrectly could not\
                                match regex %s with line %s"%(self.fields[line].header,self.fields[line].pattern, nextLine))
        else:
            raise Exception("Unrecognized header %s in file %s"%(line,parser.fileName))

    def makeValues(self):
        #now return a dictionary of file headers and their values
        # for each field in the list of fields, make a dictionary of {header:value}
        return dict(list((self.fields[field].header,self.fields[field].value) for field in self.fields))


class FieldSpec(object):
    def __init__(self,pattern = "(.*)",multiLine = True,postProcessing=None):
        self.pattern = pattern
        self.multiLine = multiLine
        if multiLine:
            self.value = []
        else:
            self.value = None
        self.postProcessing = postProcessing
    def safeGetNextLine(self,parser):
        try:
            line = next(parser)
            return line
        except StopIteration:
            return None


    def matchLine(self,line):
        # assume what we've been given is a regular expression
        # find the matches and report all of them in a list

        # make sure line is not None
        if line:
            
            matches = re.match(self.pattern,line)
            if matches:
                # return the tuple of matches
                retval = matches.groups()
                # if we only have one value unpack the tuple
                if len(retval) == 1:
                    retval = retval[0]
                return retval
        
        return None


    def read(self, line, parser):
        # Tries to match the current line to the pattern return true if we have 
        # a match and false if we don't
        # check if our pattern is a callable if so call it
        if callable(self.pattern):
            self.pattern(parser)
        elif isinstance(self.pattern,str):
            # always read the first line
            lineValue = self.matchLine(line)
            # make sure we could match the line
            if not lineValue:
                return False


            if self.multiLine:
                # for multi line field value should be 
                # a list of value on the lines

                self.value.append(lineValue)
                line = self.safeGetNextLine(parser) 
                nextVal = self.matchLine(line)
                while nextVal:
                    self.value.append(nextVal)
                    line = self.safeGetNextLine(parser)
                    nextVal = self.matchLine(line)
                # This kind of read will always go one past the end so we need to rewind
                parser.rewind()

            else:
                self.value = lineValue


                
                
        # call postprocessing if we need to
        if self.postProcessing is not None:
            self.postProcessing(self)
        return True
        

class HeaderFieldSpec(FieldSpec):
    def __init__(self,header,pattern = "(.*)",multiLine = False,postProcessing=None):
        self.header = header
        super(HeaderFieldSpec,self).__init__(pattern = pattern,multiLine = multiLine,postProcessing =postProcessing)


class SPOFileParser:
    def __init__(self,fileName, mode = None):
        self.file = open(fileName,"r+")
        self.fileName = fileName
        #keep a list of line positions so we can rewind
        self.lastLinePos = []
        self.lastLinePos.append(self.file.tell())
    def __del__(self):
        self.file.close()

    def __next__(self):
        #grabs the next line that isn't just a new line or comment ('#')
        lastLine = self.file.tell()
        self.lastLinePos.append(lastLine)
        
        line = ""
        while len(line) == 0:
            line = self.file.readline()
            if self.file.tell() == lastLine:
                raise StopIteration
            line = self._trimLine(line)
            lastLine = self.file.tell()

        return line
        
    def __iter__(self):
        return self
    def _trimLine(self, line):
        #look for the comment symbol
        comment = re.search("(.*)#",line)
        # if its found
        if comment:
            line = comment.group(1) # grab the first group as the line
        line = line.strip() #remove any leading or trailing whitespace
        return line
    def rewind(self):
        self.file.seek(self.lastLinePos.pop())
